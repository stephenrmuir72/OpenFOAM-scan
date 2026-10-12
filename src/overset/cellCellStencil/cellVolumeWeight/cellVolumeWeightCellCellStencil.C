/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2025 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify i
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "cellVolumeWeightCellCellStencil.H"
#include "addToRunTimeSelectionTable.H"
#include "OBJstream.H"
#include "Time.H"
#include "meshToMesh.H"
#include "cellVolumeWeightMethod.H"
#include "fvMeshSubset.H"
#include "regionSplit.H"
#include "globalIndex.H"
#include "oversetFvPatch.H"
#include "zeroGradientFvPatchFields.H"
#include "syncTools.H"
#include "oversetFvMeshBase.H"
#include "oversetFvPatchField.H"
#include "SVD.H"
#include "clockTime.H"
#include "tetOverlapVolume.H"

#include <cmath>
#include <limits>
#include <algorithm>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cellCellStencils
{
    defineTypeNameAndDebug(cellVolumeWeight, 0);
    addToRunTimeSelectionTable(cellCellStencil, cellVolumeWeight, mesh);
}
}

Foam::scalar
Foam::cellCellStencils::cellVolumeWeight::defaultOverlapTolerance_ = 1e-6;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::cellCellStencils::cellVolumeWeight::normaliseVolumeWeights
(
    scalarList& weights
)
{
    scalar scale = 0;
    for (const scalar weight : weights)
    {
        if (!std::isfinite(weight) || weight < 0)
        {
            return false;
        }
        scale = max(scale, weight);
    }
    if (scale <= 0)
    {
        return false;
    }
    long double sum = 0, correction = 0;
    for (const scalar weight : weights)
    {
        const long double adjusted =
            static_cast<long double>(weight)/scale - correction;
        const long double next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
    }
    for (scalar& weight : weights)
    {
        weight = scalar((static_cast<long double>(weight)/scale)/sum);
    }
    return checkWeights(weights) == VALID_WEIGHTS;
}


bool Foam::cellCellStencils::cellVolumeWeight::affineWeights
(
    const point& sample,
    const pointList& centres,
    scalarList& weights
) const
{
    if (weights.size() != centres.size() || checkWeights(weights) != VALID_WEIGHTS)
    {
        ++weightFallbackCounts_[1];
        return false;
    }
    vector scale(Zero);
    vector coordinateScale(Zero);
    for (const point& centre : centres)
    {
        const vector d = centre - sample;
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            if (!std::isfinite(d[cmpt]))
            {
                ++weightFallbackCounts_[1];
                return false;
            }
            scale[cmpt] = max(scale[cmpt], mag(d[cmpt]));
            coordinateScale[cmpt] = max
            (
                coordinateScale[cmpt],
                max(mag(centre[cmpt]), mag(sample[cmpt]))
            );
        }
    }
    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        if (scale[cmpt] <= 128*std::numeric_limits<scalar>::epsilon()*coordinateScale[cmpt])
        {
            scale[cmpt] = 1;
        }
    }

    // Zero padding keeps the SVD tall even for one/two-donor planar stencils.
    RectangularMatrix<scalar> A(centres.size() + 4, 4, Zero);
    FixedList<long double, 4> residual({1, 0, 0, 0});
    forAll(centres, j)
    {
        A[j][0] = 1;
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            A[j][cmpt+1] = (centres[j][cmpt]-sample[cmpt])/scale[cmpt];
        }
        for (label k = 0; k < 4; ++k)
        {
            residual[k] -= static_cast<long double>(weights[j])*A[j][k];
        }
    }
    const scalar epsilon = std::numeric_limits<scalar>::epsilon();
    SVD svd(A, 100*epsilon);
    if (!svd.converged())
    {
        ++weightFallbackCounts_[0];
        return false;
    }
    const auto inverse = svd.VSinvUt();
    scalarList corrected(weights);
    forAll(corrected, j)
    {
        long double delta = 0;
        for (label k = 0; k < 4; ++k)
        {
            delta += static_cast<long double>(inverse[k][j])*residual[k];
        }
        corrected[j] += scalar(delta);
        if (!std::isfinite(corrected[j]))
        {
            ++weightFallbackCounts_[1];
            return false;
        }
    }
    if (!allowNegativeWeights_)
    {
        const scalar tolerance = weightTolerance(corrected);
        for (scalar& weight : corrected)
        {
            if (weight < -tolerance)
            {
                ++weightFallbackCounts_[2];
                return false;
            }
            weight = max(weight, scalar(0));
        }
    }
    FixedList<long double, 4> moment({-1, 0, 0, 0});
    long double sumAbs = 0;
    forAll(corrected, j)
    {
        sumAbs += mag(corrected[j]);
        for (label k = 0; k < 4; ++k)
        {
            moment[k] += static_cast<long double>(corrected[j])*A[j][k];
        }
    }
    const scalar tolerance = max(scalar(1e-10), 100*epsilon*centres.size());
    for (label k = 0; k < 4; ++k)
    {
        if
        (
            sumAbs > 1/std::sqrt(epsilon)
         || std::abs(moment[k]) > tolerance
        )
        {
            ++weightFallbackCounts_[0];
            return false;
        }
    }
    if (checkWeights(corrected) != VALID_WEIGHTS)
    {
        ++weightFallbackCounts_[1];
        return false;
    }
    weights.transfer(corrected);
    return true;
}


void Foam::cellCellStencils::cellVolumeWeight::reportWeightFallbacks
(
    label nAffineExact,
    label nAffineFallback
) const
{
    labelList counts(weightFallbackCounts_);
    Pstream::listReduce(counts, sumOp<label>());
    weightFallbackCounts_ = Zero;
    if (sum(counts))
    {
        if (!weightFallbackReported_)
        {
            WarningInFunction
                << "cellVolumeWeight affine correction fell back to validated "
                << "volume/IDW weights. Subsequent occurrences are summarized."
                << nl << endl;
            weightFallbackReported_ = true;
        }
        Info<< "cellVolumeWeight affine fallback calls: "
            << counts[0] << " infeasible/ill-conditioned, "
            << counts[1] << " invalid residual/weights, "
            << counts[2] << " negative weights." << endl;
    }
    reduce(nAffineExact, sumOp<label>());
    reduce(nAffineFallback, sumOp<label>());
    Info<< "cellVolumeWeight affine receivers: " << nAffineExact
        << " exact, " << nAffineFallback << " volume/IDW fallback." << endl;
}


Foam::label Foam::cellCellStencils::cellVolumeWeight::decodeRestartType
(
    const scalar value,
    bool& rounded,
    bool& clamped
)
{
    if (!std::isfinite(value))
    {
        FatalErrorInFunction
            << "Invalid restart cellTypes value " << value
            << ": classification must be finite." << exit(FatalError);
    }
    const scalar nearest = std::round(value);
    rounded = mag(value-nearest) > 100*SMALL;
    clamped = nearest < scalar(CALCULATED) || nearest > scalar(POROUS);
    return label(max(scalar(CALCULATED), min(scalar(POROUS), nearest)));
}


void Foam::cellCellStencils::cellVolumeWeight::refreshGeometry()
{
    geometricVolumes_.setSize(mesh_.nCells());
    geometricRoundoff_.setSize(mesh_.nCells());
    label mismatches = 0;
    scalar maxMismatch = 0;
    forAll(geometricVolumes_, celli)
    {
        const auto reference = tetOverlapVolume::cellVolumeMinDecomp(mesh_, celli);
        geometricVolumes_[celli] = reference.first();
        geometricRoundoff_[celli] = reference.second();
        const scalar mismatch = mag(reference.first()/mesh_.V()[celli]-1);
        mismatches += mismatch > reference.second();
        maxMismatch = max(maxMismatch, mismatch);
    }
    reduce(mismatches, sumOp<label>());
    reduce(maxMismatch, maxOp<scalar>());
    if (mismatches)
    {
        Info<< "cellVolumeWeight geometric/FV volume mismatch: " << mismatches
            << " owned cells, maximum relative mismatch " << maxMismatch
            << ". Using overlap-tessellation volumes; FV volumes unchanged."
            << endl;
    }
}


bool Foam::cellCellStencils::cellVolumeWeight::candidateCoverage
(
    const label celli,
    const scalarList& original,
    const scalarList& retained,
    scalar& coverage
) const
{
    long double originalVolume = 0, retainedVolume = 0;
    for (const scalar volume : original)
    {
        if (!std::isfinite(volume) || volume < 0) return false;
        originalVolume += volume;
    }
    for (const scalar volume : retained)
    {
        if (!std::isfinite(volume) || volume < 0) return false;
        retainedVolume += volume;
    }
    if (originalVolume <= 0 || retainedVolume <= 0) return false;
    const long double geometricCoverage =
        retainedVolume/geometricVolumes_[celli];
    const scalar allowance = geometricRoundoff_[celli];
    if
    (
        !std::isfinite(geometricCoverage)
     || geometricCoverage > 1 + allowance
     || geometricCoverage + allowance < minCoverage_
     || retainedVolume/originalVolume + allowance < minRetainedFraction_
    )
    {
        return false;
    }
    coverage = scalar(geometricCoverage);
    return true;
}


bool Foam::cellCellStencils::cellVolumeWeight::betterCandidate
(
    const label celli,
    const label zone,
    const scalar coverage,
    const label bestZone,
    const scalar bestCoverage
) const
{
    if (bestZone < 0) return true;
    const scalar tie = 100*SMALL*max(scalar(1), coverage);
    if (donorSelection_ == "coverage")
    {
        if (coverage > bestCoverage+tie) return true;
        if (coverage < bestCoverage-tie) return false;
    }
    const label acceptorZone = zoneID()[celli];
    const label diff = mag(acceptorZone-zone);
    const label bestDiff = mag(acceptorZone-bestZone);
    return diff < bestDiff || (diff == bestDiff && zone > bestZone);
}


Foam::scalar Foam::cellCellStencils::cellVolumeWeight::weightedDonorVolume
(
    const pointList& centres,
    const scalarList& volumes,
    scalarList weights
)
{
    if (centres.size() != volumes.size() || weights.size() != volumes.size())
    {
        FatalErrorInFunction
            << "Misaligned weighted donor geometry" << exit(FatalError);
    }
    if (!normaliseVolumeWeights(weights)) return 0;
    forAll(volumes, i)
    {
        if
        (
            !std::isfinite(volumes[i]) || volumes[i] <= 0
         || !std::isfinite(centres[i].x())
         || !std::isfinite(centres[i].y())
         || !std::isfinite(centres[i].z())
        )
        {
            FatalErrorInFunction
                << "Invalid weighted donor geometry" << exit(FatalError);
        }
    }
    labelList order(Foam::identity(weights.size()));
    std::sort(order.begin(), order.end(), [&](const label a, const label b)
    {
        for (direction c = 0; c < 3; ++c)
        {
            if (centres[a][c] != centres[b][c])
            {
                return centres[a][c] < centres[b][c];
            }
        }
        if (volumes[a] != volumes[b]) return volumes[a] < volumes[b];
        return weights[a] < weights[b];
    });
    static bool ambiguityReported = false;
    for (label i = 1; i < order.size() && !ambiguityReported; ++i)
    {
        const label a = order[i-1], b = order[i];
        if (centres[a] == centres[b])
        {
            WarningInFunction
                << "Coincident donor centres in weighted "
                << "layer scoring: physical donor identities are ambiguous; "
                << "global labels cannot resolve this reproducibly."
                << nl << endl;
            ambiguityReported = true;
        }
    }
    long double volume = 0, correction = 0;
    for (const label i : order)
    {
        const long double term =
            static_cast<long double>(weights[i])*volumes[i]-correction;
        const long double next = volume+term;
        correction = (next-volume)-term;
        volume = next;
    }
    return scalar(volume);
}


void Foam::cellCellStencils::cellVolumeWeight::finaliseDonors
(
    const globalIndex& globalCells,
    const boolList& previouslyBlanked,
    const boolList& sameZoneFallback,
    const List<DynamicList<donorCandidate>>& candidates,
    labelList& types,
    labelListList& stencil,
    scalarListList& weights,
    labelList& donorZones,
    const boolList* donorExclusions
) const
{
    labelListList joined(stencil.size());
    forAll(types, celli)
    {
        if (types[celli] == INTERPOLATED)
        {
            joined[celli] = stencil[celli];
        }
    }
    forAll(candidates, celli)
    {
        if (types[celli] != INTERPOLATED)
        {
            continue;
        }
        DynamicList<label> donors(joined[celli]);
        for (const auto& candidate : candidates[celli])
        {
            donors.append(candidate.donors);
        }
        joined[celli] = labelList(donors);
    }
    List<Map<label>> compactMap;
    const clockTime mapBuildTimer;
    mapDistribute map(globalCells, joined, compactMap);
    scalar mapBuildSeconds = mapBuildTimer.elapsedTime();
    const clockTime mapDistributionTimer;
    labelList globalIDs(Foam::identity(globalCells.range()));
    map.distribute(globalIDs);
    Map<label> slots(globalIDs.size()-mesh_.nCells());
    for (label slot = mesh_.nCells(); slot < globalIDs.size(); ++slot)
    {
        slots.insert(globalIDs[slot], slot);
    }
    labelList oldHoles(previouslyBlanked.size(), Zero);
    forAll(oldHoles, celli)
    {
        oldHoles[celli] = previouslyBlanked[celli]
            || (donorExclusions && (*donorExclusions)[celli]);
    }
    map.distribute(oldHoles);
    scalar mapDistributionSeconds = mapDistributionTimer.elapsedTime();
    labelList reasons(3, Zero);
    label rounds = 0;
    for (;;)
    {
        ++rounds;
        labelList compactTypes(types);
        const clockTime distributionTimer;
        map.distribute(compactTypes);
        mapDistributionSeconds += distributionTimer.elapsedTime();
        bool demoted = false;
        forAll(types, celli)
        {
            if (types[celli] != INTERPOLATED)
            {
                stencil[celli].clear();
                weights[celli].clear();
                continue;
            }
            auto retain = [&](const labelList& donors, const scalarList& raw,
                              labelList& kept, scalarList& keptWeights) -> bool
            {
                if (donors.size() != raw.size())
                {
                    ++reasons[0];
                    return false;
                }
                DynamicList<label> retained;
                DynamicList<scalar> retainedWeights;
                forAll(donors, j)
                {
                    label slot = -1;
                    if (globalCells.isLocal(donors[j]))
                    {
                        slot = globalCells.toLocal(donors[j]);
                    }
                    else
                    {
                        const auto iter = slots.cfind(donors[j]);
                        if (iter.good())
                        {
                            slot = iter();
                        }
                    }
                    if (slot < 0)
                    {
                        ++reasons[0];
                        return false;
                    }
                    const label type = compactTypes[slot];
                    if
                    (
                        oldHoles[slot] || type == HOLE || type == POROUS
                     || (!allowInterpolatedDonors_ && type == INTERPOLATED)
                    )
                    {
                        ++reasons[1];
                    }
                    else
                    {
                        retained.append(donors[j]);
                        retainedWeights.append(raw[j]);
                    }
                }
                kept = labelList(retained);
                keptWeights = scalarList(retainedWeights);
                return true;
            };
            labelList kept;
            scalarList keptWeights;
            bool valid = false;
            if (!candidates.empty() && !sameZoneFallback[celli])
            {
                scalar bestCoverage = -1;
                label bestZone = -1;
                for (const auto& candidate : candidates[celli])
                {
                    labelList trial;
                    scalarList trialWeights;
                    if (!retain(candidate.donors, candidate.volumes, trial, trialWeights))
                    {
                        continue;
                    }
                    scalar coverage = 0;
                    if
                    (
                        !candidateCoverage
                        (
                            celli, candidate.volumes, trialWeights, coverage
                        )
                     || !normaliseVolumeWeights(trialWeights)
                    )
                    {
                        ++reasons[2];
                        continue;
                    }
                    if (betterCandidate(celli, candidate.zone, coverage, bestZone, bestCoverage))
                    {
                        kept.transfer(trial);
                        keptWeights.transfer(trialWeights);
                        bestCoverage = coverage;
                        bestZone = candidate.zone;
                        valid = true;
                    }
                }
                if (valid)
                {
                    donorZones[celli] = bestZone;
                }
            }
            else
            {
                valid = retain(stencil[celli], weights[celli], kept, keptWeights)
                     && normaliseVolumeWeights(keptWeights);
                if (!valid)
                {
                    ++reasons[2];
                }
            }
            if (valid)
            {
                stencil[celli].transfer(kept);
                weights[celli].transfer(keptWeights);
            }
            else
            {
                types[celli] = previouslyBlanked[celli] ? HOLE : POROUS;
                stencil[celli].clear();
                weights[celli].clear();
                demoted = true;
            }
        }
        if (!returnReduceOr(demoted))
        {
            break;
        }
    }
    Pstream::listReduce(reasons, sumOp<label>());
    if (sum(reasons))
    {
        if (!invalidStencilReported_)
        {
            WarningInFunction
                << "cellVolumeWeight discarded invalid/unready donors; "
                << "subsequent occurrences are summarized." << nl << endl;
            invalidStencilReported_ = true;
        }
        Info<< "cellVolumeWeight donor validation: " << reasons[0]
            << " misaligned/addressing events, " << reasons[1]
            << " unready donor events, " << reasons[2]
            << " rejected weight/coverage events." << endl;
    }
    if (debug)
    {
        reduce(mapBuildSeconds, maxOp<scalar>());
        reduce(mapDistributionSeconds, maxOp<scalar>());
        const label halo = returnReduce
        (
            map.constructSize()-mesh_.nCells(), sumOp<label>()
        );
        Info<< "cellVolumeWeight donor validation rounds: " << rounds << endl;
        Info<< "cellVolumeWeight donor map profile: build " << mapBuildSeconds
            << " s, distribution " << mapDistributionSeconds
            << " s, summed halo slots " << halo << endl;
    }
}


void Foam::cellCellStencils::cellVolumeWeight::findHoles
(
    const globalIndex& globalCells,
    const fvMesh& mesh,
    const labelList& zoneID,
    const labelListList& stencil,
    labelList& cellTypes
) const
{
    const fvBoundaryMesh& pbm = mesh.boundary();
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    // so we start a walk from our patches and any cell we cannot reach
    // (because we walk is stopped by other-mesh patch) is a hole.


    boolList isBlockedFace(mesh.nFaces(), false);
    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        label ownType = cellTypes[own[faceI]];
        label neiType = cellTypes[nei[faceI]];
        if
        (
             (ownType == HOLE && neiType != HOLE)
          || (ownType != HOLE && neiType == HOLE)
        )
        {
            isBlockedFace[faceI] = true;
        }
    }

    labelList nbrCellTypes;
    syncTools::swapBoundaryCellList(mesh, cellTypes, nbrCellTypes);

    for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); faceI++)
    {
        label ownType = cellTypes[own[faceI]];
        label neiType = nbrCellTypes[faceI-mesh.nInternalFaces()];

        if
        (
             (ownType == HOLE && neiType != HOLE)
          || (ownType != HOLE && neiType == HOLE)
        )
        {
            isBlockedFace[faceI] = true;
        }
    }

    regionSplit cellRegion(mesh, isBlockedFace);

    Info<< typeName << " : detected " << cellRegion.nRegions()
        << " mesh regions after overset" << nl << endl;



    // Now we'll have a mesh split according to where there are cells
    // covered by the other-side patches. See what we can reach from our
    // real patches

    //  0 : region not yet determined
    //  1 : borders blockage so is not ok (but can be overridden by real
    //      patch)
    //  2 : has real patch in it so is reachable
    labelList regionType(cellRegion.nRegions(), Zero);


    // See if any regions borders blockage. Note: isBlockedFace is already
    // parallel synchronised.
    {
        for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            if (isBlockedFace[faceI])
            {
                label ownRegion = cellRegion[own[faceI]];

                if (cellTypes[own[faceI]] != HOLE)
                {
                    if (regionType[ownRegion] == 0)
                    {
                        //Pout<< "Mark region:" << ownRegion
                        //    << " on zone:" << zoneID[own[faceI]]
                        //    << " as next to blockage at:"
                        //    << mesh.faceCentres()[faceI] << endl;
                        regionType[ownRegion] = 1;
                    }
                }

                label neiRegion = cellRegion[nei[faceI]];

                if (cellTypes[nei[faceI]] != HOLE)
                {
                    if (regionType[neiRegion] == 0)
                    {
                        //Pout<< "Mark region:" << neiRegion
                        //    << " on zone:" << zoneID[nei[faceI]]
                        //    << " as next to blockage at:"
                        //    << mesh.faceCentres()[faceI] << endl;
                        regionType[neiRegion] = 1;
                    }
                }
            }
        }
        for
        (
            label faceI = mesh.nInternalFaces();
            faceI < mesh.nFaces();
            faceI++
        )
        {
            if (isBlockedFace[faceI])
            {
                label ownRegion = cellRegion[own[faceI]];

                if (regionType[ownRegion] == 0)
                {
                    //Pout<< "Mark region:" << ownRegion
                    //    << " on zone:" << zoneID[own[faceI]]
                    //    << " as next to blockage at:"
                    //    << mesh.faceCentres()[faceI] << endl;
                    regionType[ownRegion] = 1;
                }
            }
        }
    }


    // Override with real patches
    forAll(pbm, patchI)
    {
        const fvPatch& fvp = pbm[patchI];

        if (isA<oversetFvPatch>(fvp))
        {}
        else if (!fvPatch::constraintType(fvp.type()))
        {
            //Pout<< "Proper patch " << fvp.name() << " of type " << fvp.type()
            //    << endl;

            for (const label celli : fvp.faceCells())
            {
                label regionI = cellRegion[celli];

                if (cellTypes[celli] != HOLE && regionType[regionI] != 2)
                {
                    //Pout<< "reachable region : " << regionI
                    //    << " at cell " << mesh.cellCentres()[celli]
                    //    << " on zone " << zoneID[celli] << endl;
                    regionType[regionI] = 2;
                }
            }
        }
    }

    // Now we've handled
    // - cells next to blocked cells
    // - coupled boundaries
    // Only thing to handle is the interpolation between regions


    labelListList compactStencil(stencil);
    List<Map<label>> compactMap;
    mapDistribute map(globalCells, compactStencil, compactMap);

    while (true)
    {
        // Synchronise region status on processors
        // (could instead swap status through processor patches)
        Pstream::listReduce(regionType, maxOp<label>());

        // Communicate region status through interpolative cells
        labelList cellRegionType(labelUIndList(regionType, cellRegion));
        map.distribute(cellRegionType);


        label nChanged = 0;
        forAll(pbm, patchI)
        {
            const fvPatch& fvp = pbm[patchI];

            if (isA<oversetFvPatch>(fvp))
            {
                const labelUList& fc = fvp.faceCells();
                forAll(fc, i)
                {
                    label cellI = fc[i];
                    label regionI = cellRegion[cellI];

                    if (regionType[regionI] != 2)
                    {
                        const labelList& slots = compactStencil[cellI];
                        forAll(slots, i)
                        {
                            label otherType = cellRegionType[slots[i]];

                            if (otherType == 2)
                            {
                                //Pout<< "Reachable through interpolation : "
                                //    << regionI << " at cell "
                                //    << mesh.cellCentres()[cellI] << endl;
                                regionType[regionI] = 2;
                                nChanged++;
                                break;
                            }
                        }
                    }
                }
            }
        }


        if (!returnReduceOr(nChanged))
        {
            break;
        }
    }


    // See which regions have not been visited (regionType == 1)
    forAll(cellRegion, cellI)
    {
        label type = regionType[cellRegion[cellI]];
        if (type == 1 && cellTypes[cellI] != HOLE)
        {
            cellTypes[cellI] = HOLE;
        }
    }
}


void Foam::cellCellStencils::cellVolumeWeight::markPatchCells
(
    const fvMesh& mesh,
    const labelList& cellMap,
    labelList& patchCellTypes
) const
{
    const fvBoundaryMesh& pbm = mesh.boundary();

    forAll(pbm, patchI)
    {
        const fvPatch& fvp = pbm[patchI];
        const labelUList& fc = fvp.faceCells();

        if (isA<oversetFvPatch>(fvp))
        {
            //Pout<< "Marking cells on overset patch " << fvp.name() << endl;
            forAll(fc, i)
            {
                label cellI = fc[i];
                patchCellTypes[cellMap[cellI]] = OVERSET;
            }
        }
        else if (!fvPatch::constraintType(fvp.type()))
        {
            //Pout<< "Marking cells on proper patch " << fvp.name()
            //    << " with type " << fvp.type() << endl;
            forAll(fc, i)
            {
                label cellI = fc[i];
                if (patchCellTypes[cellMap[cellI]] != OVERSET)
                {
                    patchCellTypes[cellMap[cellI]] = PATCH;
                }
            }
        }
    }
}


void Foam::cellCellStencils::cellVolumeWeight::interpolatePatchTypes
(
    const labelListList& addressing,
    const labelList& patchTypes,
    labelList& result
) const
{
    forAll(result, cellI)
    {
        const labelList& slots = addressing[cellI];
        forAll(slots, i)
        {
            label type = patchTypes[slots[i]];

            if (type == OVERSET)
            {
                // 'overset' overrides anything
                result[cellI] = OVERSET;
                break;
            }
            else if (type == PATCH)
            {
                // 'patch' overrides -1 and 'other'
                result[cellI] = PATCH;
            }
            else if (result[cellI] == -1)
            {
                // 'other' overrides -1 only
                result[cellI] = OTHER;
            }
        }
    }
}


void Foam::cellCellStencils::cellVolumeWeight::interpolatePatchTypes
(
    const autoPtr<mapDistribute>& mapPtr,
    const labelListList& addressing,
    const labelList& patchTypes,
    labelList& result
) const
{
    if (result.size() != addressing.size())
    {
        FatalErrorInFunction << "result:" << result.size()
            << " addressing:" << addressing.size() << exit(FatalError);
    }


    // Initialise to not-mapped
    result = -1;

    if (mapPtr)
    {
        // Pull remote data into order of addressing
        labelList work(patchTypes);
        mapPtr().distribute(work);

        interpolatePatchTypes(addressing, work, result);
    }
    else
    {
        interpolatePatchTypes(addressing, patchTypes, result);
    }
}


void Foam::cellCellStencils::cellVolumeWeight::combineCellTypes
(
    const label subZoneID,
    const fvMesh& subMesh,
    const labelList& subCellMap,

    const label donorZoneID,
    const labelListList& addressing,
    const List<scalarList>& weights,
    const labelList& otherCells,
    const labelList& interpolatedOtherPatchTypes,

    labelListList& allStencil,
    scalarListList& allWeights,
    labelList& allCellTypes,
    labelList& allDonorID,
    List<DynamicList<donorCandidate>>* candidates
) const
{
    forAll(subCellMap, subCellI)
    {
        label cellI = subCellMap[subCellI];

        bool validDonors = true;
        switch (interpolatedOtherPatchTypes[subCellI])
        {
            case -1:
            {
                validDonors = false;
            }
            break;

            case OTHER:
            {
                // No patch interaction so keep valid
            }
            break;

            case PATCH:
            {
                // Patch-patch interaction... For now disable always
                allCellTypes[cellI] = HOLE;
                validDonors = false;
                // Alternative is to look at the amount of overlap but this
                // is not very robust
//                 if (allCellTypes[cellI] != HOLE)
//                 {
//                    scalar overlapVol = sum(weights[subCellI]);
//                    scalar v = mesh_.V()[cellI];
//                    if (overlapVol < (1.0-overlapTolerance_)*v)
//                    {
//                        //Pout<< "** Patch overlap:" << cellI
//                        //    << " at:" << mesh_.cellCentres()[cellI] << endl;
//                        allCellTypes[cellI] = HOLE;
//                        validDonors = false;
//                    }
//                 }
            }
            break;

            case OVERSET:
            {
                validDonors = true;
            }
            break;
        }


        if (validDonors)
        {
            if (candidates)
            {
                donorCandidate candidate;
                candidate.zone = donorZoneID;
                candidate.donors =
                    labelUIndList(otherCells, addressing[subCellI]);
                candidate.volumes = weights[subCellI];
                (*candidates)[cellI].append(std::move(candidate));
            }
            // There are a few possible choices how to choose between multiple
            // donor candidates:
            // 1 highest overlap volume. However this is generally already
            //   99.9% so you're just measuring truncation error.
            // 2 smallest donors cells or most donor cells. This is quite
            //   often done but can cause switching of donor zone from one
            //   time step to the other if the donor meshes are non-uniform
            //   and the acceptor cells just happens to be sweeping through
            //   some small donor cells.
            // 3 nearest zoneID. So zone 0 preferentially interpolates from
            //   zone 1, zone 1 preferentially from zone 2 etc.

            //- Option 1:
            //scalar currentVol = sum(allWeights[cellI]);
            //if (overlapVol[subCellI] > currentVol)

            //- Option 3:
            label currentDiff = mag(subZoneID-allDonorID[cellI]);
            label thisDiff = mag(subZoneID-donorZoneID);

            if
            (
                allDonorID[cellI] == -1
             || (thisDiff < currentDiff)
             || (thisDiff == currentDiff && donorZoneID > allDonorID[cellI])
            )
            {
                allWeights[cellI] = weights[subCellI];
                allStencil[cellI] =
                    labelUIndList(otherCells, addressing[subCellI]);
                allDonorID[cellI] = donorZoneID;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellCellStencils::cellVolumeWeight::cellVolumeWeight
(
    const fvMesh& mesh,
    const dictionary& dict,
    const bool doUpdate
)
:
    cellCellStencil(mesh),
    dict_(dict),
    overlapTolerance_(defaultOverlapTolerance_),
    cellTypes_(labelList(mesh.nCells(), CALCULATED)),
    interpolationCells_(0),
    cellInterpolationMap_(),
    cellStencil_(0),
    cellInterpolationWeights_(0),
    cellInterpolationWeight_
    (
        IOobject
        (
            "cellInterpolationWeight",
            mesh_.facesInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    allowInterpolatedDonors_
    (
        dict.getOrDefault("allowInterpolatedDonors", true)
    ),
    donorSelection_(dict.getOrDefault<word>("donorSelection", "legacy")),
    donorFallback_(dict.getOrDefault<word>("donorFallback", "legacy")),
    minCoverage_(dict.getOrDefault<scalar>("minCoverage", 0)),
    minRetainedFraction_(dict.getOrDefault<scalar>("minRetainedFraction", 0)),
    weightingScheme_(dict.getOrDefault<word>("weightingScheme", "volume")),
    allowNegativeWeights_(dict.getOrDefault("allowNegativeWeights", false)),
    stencilExpansionRings_(dict.getOrDefault<label>("stencilExpansionRings", 1)),
    layerDonorMetric_(dict.getOrDefault<word>("layerDonorMetric", "legacy")),
    processorMapMethod_(dict.getOrDefault<word>("processorMapMethod", "AABB")),
    geometricVolumes_(mesh.nCells()),
    geometricRoundoff_(mesh.nCells()),
    weightFallbackCounts_(3, Zero),
    weightFallbackReported_(false)
{
    // Add motion-solver fields to non-interpolated fields
    suppressMotionFields();

    // Read zoneID
    this->zoneID();

    overlapTolerance_ =
        dict_.getOrDefault("overlapTolerance", defaultOverlapTolerance_);

    if
    (
        !std::isfinite(overlapTolerance_)
     || overlapTolerance_ < 0 || overlapTolerance_ > 1
     || !std::isfinite(minCoverage_)
     || minCoverage_ < 0 || minCoverage_ > 1
     || !std::isfinite(minRetainedFraction_)
     || minRetainedFraction_ < 0 || minRetainedFraction_ > 1
     || (donorSelection_ != "legacy" && donorSelection_ != "coverage")
     || (donorFallback_ != "legacy" && donorFallback_ != "alternateZones"
         && donorFallback_ != "sameZone")
     || (weightingScheme_ != "volume" && weightingScheme_ != "linearExact")
     || (layerDonorMetric_ != "legacy" && layerDonorMetric_ != "overlapWeighted")
     || (processorMapMethod_ != "AABB" && processorMapMethod_ != "LOD")
     || stencilExpansionRings_ < 0 || stencilExpansionRings_ > 1
    )
    {
        FatalIOErrorInFunction(dict_)
            << "Invalid cellVolumeWeight configuration: tolerances must be "
            << "finite fractions in [0,1]; donorSelection is legacy|coverage, "
            << "donorFallback is legacy|alternateZones|sameZone, "
            << "weightingScheme is volume|linearExact, layerDonorMetric is "
            << "legacy|overlapWeighted, processorMapMethod is AABB|LOD, and "
            << "stencilExpansionRings is 0 or 1." << exit(FatalIOError);
    }

    // Read old-time cellTypes
    IOobject io
    (
        "cellTypes",
        mesh_.time().timeName(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );
    const bool haveRestart = io.typeHeaderOk<volScalarField>(true);
    const bool anyRestart = returnReduceOr(haveRestart);
    const label missing = returnReduce
    (
        label(!haveRestart && mesh_.nCells()), sumOp<label>()
    );
    if (anyRestart && missing)
    {
        FatalErrorInFunction
            << "Incomplete restart cellTypes: missing on " << missing
            << " non-empty ranks." << exit(FatalError);
    }
    labelList repairs(3, Zero);
    if (haveRestart)
    {
        if (debug)
        {
            Pout<< "Reading cellTypes from time " << mesh_.time().timeName()
                << endl;
        }

        const volScalarField volCellTypes(io, mesh_);
        forAll(volCellTypes, celli)
        {
            bool rounded = false, clamped = false;
            cellTypes_[celli] = decodeRestartType
            (
                volCellTypes[celli], rounded, clamped
            );
            repairs[0] += rounded || clamped;
            repairs[1] += rounded;
            repairs[2] += clamped;
        }
    }
    Pstream::listReduce(repairs, sumOp<label>());
    if (repairs[0])
    {
        WarningInFunction
            << "cellVolumeWeight restart repair: " << repairs[0]
            << " affected cells, " << repairs[1] << " rounded, "
            << repairs[2] << " clamped. cellTypes is categorical; arithmetic "
            << "mapping does not reconstruct physical holes." << nl << endl;
    }
    if (doUpdate)
    {
        update();
    }
    else
    {
        refreshGeometry();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellCellStencils::cellVolumeWeight::~cellVolumeWeight()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cellCellStencils::cellVolumeWeight::update()
{
    const clockTime updateTimer;
    const clockTime referenceTimer;
    refreshGeometry();
    const scalar referenceSeconds = referenceTimer.elapsedTime();
    boolList previouslyBlanked(mesh_.nCells(), false);
    forAll(previouslyBlanked, celli)
    {
        previouslyBlanked[celli] = cellTypes_[celli] == HOLE;
    }

    scalar layerRelax(dict_.getOrDefault("layerRelax", 1.0));
    const labelIOList& zoneID = this->zoneID();

    label nZones = gMax(zoneID)+1;
    labelList nCellsPerZone(nZones, Zero);
    forAll(zoneID, cellI)
    {
        nCellsPerZone[zoneID[cellI]]++;
    }
    Pstream::listReduce(nCellsPerZone, sumOp<label>());

    Info<< typeName << " : detected " << nZones
        << " mesh regions" << nl << endl;


    PtrList<fvMeshSubset> meshParts(nZones);

    Info<< incrIndent;
    forAll(meshParts, zonei)
    {
        Info<< indent<< "zone:" << zonei << " nCells:"
            << nCellsPerZone[zonei] << nl;

        meshParts.set
        (
            zonei,
            new fvMeshSubset(mesh_, zonei, zoneID)
        );
    }
    Info<< decrIndent;


    // Current best guess for cells. Includes best stencil. Weights should
    // add up to volume.
    labelList allCellTypes(mesh_.nCells(), CALCULATED);
    labelList allPatchTypes(mesh_.nCells(), OTHER);
    labelListList allStencil(mesh_.nCells());
    scalarListList allWeights(mesh_.nCells());
    // zoneID of donor
    labelList allDonorID(mesh_.nCells(), -1);
    const bool qualitySelection =
        donorSelection_ == "coverage"
     || minCoverage_ > 0 || minRetainedFraction_ > 0;
    const bool retainCandidates = qualitySelection || donorFallback_ != "legacy";
    List<DynamicList<donorCandidate>> candidates
    (
        retainCandidates ? mesh_.nCells() : 0
    );


    // Marking patch cells
    forAll(meshParts, partI)
    {
        const fvMesh& partMesh = meshParts[partI].subMesh();
        const labelList& partCellMap = meshParts[partI].cellMap();

        // Mark cells with
        // - overset boundary
        // - other, proper boundary
        // - other cells
        Info<< "Marking patch-cells on zone " << partI << endl;
        markPatchCells(partMesh, partCellMap, allPatchTypes);
    }

    if ((debug & 2) && mesh_.time().writeTime())
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allPatchTypes", allPatchTypes)
        );
        tfld().write();
    }


    labelList nCells(count(3, allPatchTypes));
    Info<< nl
        << "After patch analysis : nCells : "
        << returnReduce(allPatchTypes.size(), sumOp<label>()) << nl
        << incrIndent
        << indent << "other  : " << nCells[OTHER] << nl
        << indent << "patch  : " << nCells[PATCH] << nl
        << indent << "overset: " << nCells[OVERSET] << nl
        << decrIndent << endl;

    globalIndex globalCells(mesh_.nCells());

    List<boundBox> zoneBounds(nZones);
    forAll(meshParts, zonei)
    {
        if (nCellsPerZone[zonei])
        {
            zoneBounds[zonei] =
                boundBox(meshParts[zonei].subMesh().points(), true);
        }
    }
    const clockTime overlapTimer;
    for (label srcI = 0; srcI < meshParts.size()-1; srcI++)
    {
        const fvMesh& srcMesh = meshParts[srcI].subMesh();
        const labelList& srcCellMap = meshParts[srcI].cellMap();

        for (label tgtI = srcI+1; tgtI < meshParts.size(); tgtI++)
        {
            if
            (
                !nCellsPerZone[srcI] || !nCellsPerZone[tgtI]
             || !zoneBounds[srcI].overlaps(zoneBounds[tgtI])
            )
            {
                continue;
            }
            const fvMesh& tgtMesh = meshParts[tgtI].subMesh();
            const labelList& tgtCellMap = meshParts[tgtI].cellMap();

            meshToMesh mapper
            (
                srcMesh,
                tgtMesh,
                meshToMesh::interpolationMethod::imCellVolumeWeight,
                HashTable<word>(0),     // patchMap,
                wordList(0),            // cuttingPatches
                processorMapMethod_ == "LOD"
              ? meshToMesh::procMapMethod::pmLOD
              : meshToMesh::procMapMethod::pmAABB,
                false                   // do not normalise
            );

            {
                // Get tgt patch types on src mesh
                labelList interpolatedTgtPatchTypes(srcMesh.nCells(), -1);
                interpolatePatchTypes
                (
                    mapper.tgtMap(),            // How to get remote data local
                    mapper.srcToTgtCellAddr(),
                    labelList(labelUIndList(allPatchTypes, tgtCellMap)),
                    interpolatedTgtPatchTypes
                );

                // Get target cell labels in global cell indexing (on overall
                // mesh)
                labelList tgtGlobalCells(tgtMesh.nCells());
                {
                    forAll(tgtCellMap, tgtCellI)
                    {
                        label cellI = tgtCellMap[tgtCellI];
                        tgtGlobalCells[tgtCellI] = globalCells.toGlobal(cellI);
                    }
                    if (mapper.tgtMap())
                    {
                        mapper.tgtMap()->distribute(tgtGlobalCells);
                    }
                }


                combineCellTypes
                (
                    srcI,
                    srcMesh,
                    srcCellMap,

                    tgtI,
                    mapper.srcToTgtCellAddr(),
                    mapper.srcToTgtCellWght(),
                    tgtGlobalCells,
                    interpolatedTgtPatchTypes,

                    // Overall mesh data
                    allStencil,
                    allWeights,
                    allCellTypes,
                    allDonorID,
                    candidates.empty() ? nullptr : &candidates
                );
            }

            {
                // Get src patch types on tgt mesh
                labelList interpolatedSrcPatchTypes(tgtMesh.nCells(), -1);
                interpolatePatchTypes
                (
                    mapper.srcMap(),            // How to get remote data local
                    mapper.tgtToSrcCellAddr(),
                    labelList(labelUIndList(allPatchTypes, srcCellMap)),
                    interpolatedSrcPatchTypes
                );

                labelList srcGlobalCells(srcMesh.nCells());
                {
                    forAll(srcCellMap, srcCellI)
                    {
                        label cellI = srcCellMap[srcCellI];
                        srcGlobalCells[srcCellI] = globalCells.toGlobal(cellI);
                    }
                    if (mapper.srcMap())
                    {
                        mapper.srcMap()->distribute(srcGlobalCells);
                    }
                }

                combineCellTypes
                (
                    tgtI,
                    tgtMesh,
                    tgtCellMap,

                    srcI,
                    mapper.tgtToSrcCellAddr(),
                    mapper.tgtToSrcCellWght(),
                    srcGlobalCells,
                    interpolatedSrcPatchTypes,

                    // Overall mesh data
                    allStencil,
                    allWeights,
                    allCellTypes,
                    allDonorID,
                    candidates.empty() ? nullptr : &candidates
                );
            }
        }
    }


    const scalar overlapSeconds = overlapTimer.elapsedTime();
    if ((debug & 2) && mesh_.time().writeTime())
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allCellTypes", allCellTypes)
        );
        tfld().write();

        tmp<volScalarField> tdonors
        (
            createField(mesh_, "allDonorID", allDonorID)
        );
        tdonors().write();

    }


    labelList originalDonorID(retainCandidates ? allDonorID : labelList());
    if (retainCandidates)
    {
        forAll(allCellTypes, celli)
        {
            if
            (
                allCellTypes[celli] == HOLE
             || sum(allWeights[celli]) > overlapTolerance_*geometricVolumes_[celli]
            )
            {
                continue;
            }
            scalar bestCoverage = -1;
            label bestZone = -1;
            for (const auto& candidate : candidates[celli])
            {
                const scalar coverage =
                    sum(candidate.volumes)/geometricVolumes_[celli];
                if
                (
                    coverage > overlapTolerance_
                 && betterCandidate(celli, candidate.zone, coverage, bestZone, bestCoverage)
                )
                {
                    allStencil[celli] = candidate.donors;
                    allWeights[celli] = candidate.volumes;
                    allDonorID[celli] = candidate.zone;
                    bestZone = candidate.zone;
                    bestCoverage = coverage;
                }
            }
        }
    }

    // Use the patch types and weights to decide what to do
    forAll(allPatchTypes, cellI)
    {
        if (allCellTypes[cellI] != HOLE)
        {
            switch (allPatchTypes[cellI])
            {
                case OVERSET:
                {
                    // Interpolate. Check if enough overlap
                    scalar v = geometricVolumes_[cellI];
                    scalar overlapVol = sum(allWeights[cellI]);
                    if (overlapVol > overlapTolerance_*v)
                    {
                        allCellTypes[cellI] = INTERPOLATED;
                    }
                    else
                    {
                        allCellTypes[cellI] = HOLE;
                        allWeights[cellI].clear();
                        allStencil[cellI].clear();
                    }
                    break;
                }
            }
        }
    }


    if ((debug & 2) && mesh_.time().writeTime())
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allCellTypes_patch", allCellTypes)
        );
        tfld().write();

        labelList stencilSize(mesh_.nCells());
        forAll(allStencil, celli)
        {
            stencilSize[celli] = allStencil[celli].size();
        }
        tmp<volScalarField> tfldStencil
        (
            createField(mesh_, "allStencil_patch", stencilSize)
        );
        tfldStencil().write();
    }


    // Mark unreachable bits
    labelListList holeStencil;
    if (retainCandidates)
    {
        holeStencil.setSize(mesh_.nCells());
        forAll(holeStencil, celli)
        {
            DynamicList<label> donors;
            for (const auto& candidate : candidates[celli])
            {
                donors.append(candidate.donors);
            }
            holeStencil[celli] = labelList(donors);
        }
    }
    const clockTime regionTimer;
    findHoles
    (
        globalCells, mesh_, zoneID,
        retainCandidates ? holeStencil : allStencil, allCellTypes
    );
    const scalar regionSeconds = regionTimer.elapsedTime();

    if ((debug & 2) && mesh_.time().writeTime())
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allCellTypes_hole", allCellTypes)
        );
        tfld().write();

        labelList stencilSize(mesh_.nCells());
        forAll(allStencil, celli)
        {
            stencilSize[celli] = allStencil[celli].size();
        }
        tmp<volScalarField> tfldStencil
        (
            createField(mesh_, "allStencil_hole", stencilSize)
        );
        tfldStencil().write();
    }

    const clockTime policyTimer;
    scalar frontMapBuildSeconds = 0, frontMapDistributionSeconds = 0;
    label maxFrontHalo = 0;
    boolList reactivationCandidates(mesh_.nCells(), false);
    boolList sameZoneFallback(mesh_.nCells(), false);
    boolList protectedHoles(mesh_.nCells(), false);
    const bool qualifiedSelection = retainCandidates;
    boolList ordinaryCandidates
    (
        donorFallback_ == "sameZone" ? mesh_.nCells() : 0, false
    );
    boolList donorExclusions
    (
        donorFallback_ == "sameZone" ? previouslyBlanked : boolList()
    );
    boolList geometricEligible
    (
        qualifiedSelection ? mesh_.nCells() : 0, true
    );
    if (qualifiedSelection)
    {
        forAll(allCellTypes, celli)
        {
            if (allCellTypes[celli] == HOLE) continue;
            scalar bestCoverage = -1;
            label bestZone = -1;
            for (const auto& candidate : candidates[celli])
            {
                scalar coverage = 0;
                if
                (
                    candidateCoverage(celli, candidate.volumes, candidate.volumes, coverage)
                 && sum(candidate.volumes) > overlapTolerance_*geometricVolumes_[celli]
                 && betterCandidate(celli, candidate.zone, coverage, bestZone, bestCoverage)
                )
                {
                    allStencil[celli] = candidate.donors;
                    allWeights[celli] = candidate.volumes;
                    allDonorID[celli] = candidate.zone;
                    bestCoverage = coverage;
                    bestZone = candidate.zone;
                }
            }
            geometricEligible[celli] = bestZone >= 0;
        }
    }

    auto walkInterpolationFront =
        [&](const labelListList& frontStencil, labelList& types,
            scalarField& blend, boolList& protectedCells)
    {
        labelListList compact(frontStencil);
        List<Map<label>> compactMap;
        const clockTime buildTimer;
        mapDistribute frontMap(globalCells, compact, compactMap);
        frontMapBuildSeconds += buildTimer.elapsedTime();
        maxFrontHalo = max(maxFrontHalo, frontMap.constructSize()-mesh_.nCells());
        const clockTime distributionTimer;
        scalarList donorVolumes(mesh_.V());
        frontMap.distribute(donorVolumes);
        pointField donorCentres;
        if (layerDonorMetric_ == "overlapWeighted")
        {
            donorCentres = mesh_.cellCentres();
            frontMap.distribute(donorCentres);
        }
        frontMapDistributionSeconds += distributionTimer.elapsedTime();
        scalarList volumeRatio(mesh_.nCells(), Zero);
        if (layerDonorMetric_ == "overlapWeighted")
        {
            forAll(frontStencil, celli)
            {
                if (frontStencil[celli].empty()) continue;
                const scalar volume = weightedDonorVolume
                (
                    pointList(UIndirectList<point>(donorCentres, compact[celli])),
                    scalarList(UIndirectList<scalar>(donorVolumes, compact[celli])),
                    allWeights[celli]
                );
                if (volume > 0) volumeRatio[celli] = mesh_.V()[celli]/volume;
            }
        }
        walkFront
        (
            globalCells, layerRelax, frontStencil, types, blend,
            donorVolumes, compact, zoneID,
            dict_.getOrDefault("holeLayers", 1),
            dict_.getOrDefault("useLayer", -1), &protectedCells,
            layerDonorMetric_ == "overlapWeighted" ? &volumeRatio : nullptr
        );
    };

    if (donorFallback_ == "sameZone")
    {
        // Discover requested fringe cells without replacing hole geometry.
        labelList requestedTypes(allCellTypes);
        scalarField requestedBlend(mesh_.nCells(), Zero);
        boolList requestedProtected(mesh_.nCells(), false);
        walkInterpolationFront
        (
            allStencil, requestedTypes, requestedBlend, requestedProtected
        );
        forAll(ordinaryCandidates, celli)
        {
            ordinaryCandidates[celli] =
                !previouslyBlanked[celli] && !requestedProtected[celli]
             && allCellTypes[celli] != HOLE && !geometricEligible[celli]
             && (requestedTypes[celli] == INTERPOLATED
                 || requestedTypes[celli] == HOLE);
            donorExclusions[celli] =
                donorExclusions[celli] || ordinaryCandidates[celli];
            reactivationCandidates[celli] = ordinaryCandidates[celli];
        }
    }
    labelListList reactivationStencil(mesh_.nCells());
    List<scalarList> reactivationWeights(mesh_.nCells());
    DynamicList<label> reactivationCells;

    if ((debug & 2) && mesh_.time().writeTime())
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allCellTypes_old", cellTypes_)
        );
        tfld().write();
    }

    forAll(cellTypes_, celli)
    {
        if
        (
            previouslyBlanked[celli]
         &&
            (
                allCellTypes[celli] == CALCULATED
             || allCellTypes[celli] == INTERPOLATED
            )
        )
        {
            reactivationCandidates[celli] = true;
            reactivationCells.append(celli);
        }
    }

    DynamicList<label> provisionalCells(reactivationCells);
    forAll(ordinaryCandidates, celli)
    {
        if (ordinaryCandidates[celli]) provisionalCells.append(celli);
    }

    const boolList validSameZoneDonor
    (
        recoveryDonorMask
        (
            allCellTypes, previouslyBlanked, protectedHoles,
            reactivationCandidates, allowInterpolatedDonors_
        )
    );

    labelListList sameZoneCellCells;
    pointListList sameZoneCellCentres;
    globalSameZoneCellCells
    (
        globalCells,
        mesh_,
        zoneID,
        validSameZoneDonor,
        provisionalCells,
        sameZoneCellCells,
        sameZoneCellCentres,
        allCellTypes,
        allPatchTypes
    );

    forAll(reactivationCells, i)
    {
        const label celli = reactivationCells[i];
        if (allStencil[celli].empty())
        {
            sameZoneFallback[celli] = true;
            if (!sameZoneCellCells[celli].empty())
            {
                if (layerDonorMetric_ == "overlapWeighted")
                {
                    allStencil[celli] = sameZoneCellCells[celli];
                    distanceWeights
                    (
                        mesh_.cellCentres()[celli],
                        sameZoneCellCentres[celli],
                        allWeights[celli]
                    );
                }
                else
                {
                    allStencil[celli].setSize(1);
                    allStencil[celli][0] = sameZoneCellCells[celli][0];
                    allWeights[celli].setSize(1);
                    allWeights[celli][0] = 1.0;
                }
            }
        }
        reactivationStencil[celli] = allStencil[celli];
        reactivationWeights[celli] = allWeights[celli];
    }

    forAll(ordinaryCandidates, celli)
    {
        if (ordinaryCandidates[celli] && !sameZoneCellCells[celli].empty())
        {
            allStencil[celli] = sameZoneCellCells[celli];
            distanceWeights
            (
                mesh_.cellCentres()[celli], sameZoneCellCentres[celli],
                allWeights[celli]
            );
            sameZoneFallback[celli] = true;
        }
    }

    // Add buffer interpolation layer around holes
    scalarField allWeight(mesh_.nCells(), Zero);
    labelListList qualifiedStencil;
    if (qualifiedSelection)
    {
        qualifiedStencil = allStencil;
        forAll(qualifiedStencil, celli)
        {
            if (!geometricEligible[celli] && !sameZoneFallback[celli])
            {
                qualifiedStencil[celli].clear();
            }
        }
    }
    walkInterpolationFront
    (
        qualifiedSelection ? qualifiedStencil : allStencil,
        allCellTypes,
        allWeight,
        protectedHoles
    );
    forAll(allCellTypes, celli)
    {
        if
        (
            qualifiedSelection && !geometricEligible[celli]
         && !sameZoneFallback[celli] && !protectedHoles[celli]
         && allCellTypes[celli] == HOLE
        )
        {
            allCellTypes[celli] = POROUS;
        }
    }


    if ((debug & 2) && mesh_.time().writeTime())
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allCellTypes_front", allCellTypes)
        );
        tfld().write();

        labelList stencilSize(mesh_.nCells());
        forAll(allStencil, celli)
        {
            stencilSize[celli] = allStencil[celli].size();
        }
        tmp<volScalarField> tfldStencil
        (
            createField(mesh_, "allStencil_front", stencilSize)
        );
        tfldStencil().write();
    }

    forAll(reactivationCells, i)
    {
        const label celli = reactivationCells[i];
        if (protectedHoles[celli])
        {
            allCellTypes[celli] = HOLE;
            allStencil[celli].clear();
            allWeights[celli].clear();
            allWeight[celli] = 0.0;
            continue;
        }

        if (!reactivationStencil[celli].empty())
        {
            allStencil[celli] = reactivationStencil[celli];
            allWeights[celli] = reactivationWeights[celli];
            allCellTypes[celli] = INTERPOLATED;
        }
    }


    boolList ordinaryRequested
    (
        donorFallback_ == "sameZone" ? mesh_.nCells() : 0, false
    );
    if (donorFallback_ == "sameZone")
    {
        forAll(ordinaryRequested, celli)
        {
            ordinaryRequested[celli] = !previouslyBlanked[celli]
                && !protectedHoles[celli]
                &&
                (
                    (allCellTypes[celli] == INTERPOLATED && allWeight[celli] > 0)
                 || (ordinaryCandidates[celli] && allCellTypes[celli] == POROUS)
                );
        }
    }
    finaliseDonors
    (
        globalCells, previouslyBlanked, sameZoneFallback, candidates,
        allCellTypes, allStencil, allWeights, allDonorID,
        donorFallback_ == "sameZone" ? &donorExclusions : nullptr
    );

    DynamicList<label> fallbackCells;
    label nRecovered = 0;
    forAll(reactivationCells, i)
    {
        const label celli = reactivationCells[i];
        if (protectedHoles[celli])
        {
            allCellTypes[celli] = HOLE;
            allStencil[celli].clear();
            allWeights[celli].clear();
            allWeight[celli] = 0.0;
            continue;
        }

        const bool hasValidStencil =
            !sameZoneFallback[celli]
         && allCellTypes[celli] == INTERPOLATED
         && !allStencil[celli].empty()
         && allWeights[celli].size() == allStencil[celli].size()
         && sum(allWeights[celli]) > SMALL;

        if (hasValidStencil)
        {
            allWeight[celli] = 1.0;
            ++nRecovered;
        }
        else
        {
            fallbackCells.append(celli);
        }
    }

    forAll(ordinaryRequested, celli)
    {
        if
        (
            ordinaryRequested[celli]
         && (sameZoneFallback[celli] || allCellTypes[celli] == POROUS)
        )
        {
            fallbackCells.append(celli);
            ordinaryCandidates[celli] = true;
            reactivationCandidates[celli] = true;
            donorExclusions[celli] = true;
        }
    }

    const boolList fallbackDonorMask
    (
        recoveryDonorMask
        (
            allCellTypes, previouslyBlanked, protectedHoles,
            reactivationCandidates, allowInterpolatedDonors_
        )
    );

    labelListList fallbackCellCells;
    pointListList fallbackCellCentres;
    globalSameZoneCellCells
    (
        globalCells,
        mesh_,
        zoneID,
        fallbackDonorMask,
        fallbackCells,
        fallbackCellCells,
        fallbackCellCentres,
        allCellTypes,
        allPatchTypes
    );

    label nSameZoneFallback = 0;
    label nUnrecoverable = 0;
    forAll(fallbackCells, i)
    {
        const label celli = fallbackCells[i];
        const labelList& donors = fallbackCellCells[celli];

        if (donors.empty() || (!previouslyBlanked[celli] && allWeight[celli] <= 0))
        {
            allCellTypes[celli] = previouslyBlanked[celli] ? HOLE : POROUS;
            allStencil[celli].clear();
            allWeights[celli].clear();
            allWeight[celli] = 0.0;
            ++nUnrecoverable;
        }
        else
        {
            allStencil[celli] = donors;
            distanceWeights
            (
                mesh_.cellCentres()[celli],
                fallbackCellCentres[celli],
                allWeights[celli]
            );
            allCellTypes[celli] = INTERPOLATED;
            if (previouslyBlanked[celli]) allWeight[celli] = 1.0;
            sameZoneFallback[celli] = true;
            ++nRecovered;
            ++nSameZoneFallback;
        }
    }
    if (returnReduceOr(!fallbackCells.empty()))
    {
        finaliseDonors
        (
            globalCells, previouslyBlanked, sameZoneFallback, candidates,
            allCellTypes, allStencil, allWeights, allDonorID,
            donorFallback_ == "sameZone" ? &donorExclusions : nullptr
        );
    }
    nRecovered = nSameZoneFallback = nUnrecoverable = 0;
    for (const label celli : reactivationCells)
    {
        if (protectedHoles[celli])
        {
            continue;
        }
        if (allCellTypes[celli] == INTERPOLATED)
        {
            ++nRecovered;
            nSameZoneFallback += sameZoneFallback[celli];
        }
        else
        {
            ++nUnrecoverable;
        }
    }
    reduce(nRecovered, sumOp<label>());
    reduce(nSameZoneFallback, sumOp<label>());
    reduce(nUnrecoverable, sumOp<label>());
    if (debug && nRecovered)
    {
        DebugInfo<< "Previously blanked-cell recovery: "
            << nRecovered << " reactivated, "
            << nSameZoneFallback << " using same-zone fallback."
            << endl;
    }
    if (nUnrecoverable)
    {
        WarningInFunction
            << "Unable to reconstruct " << nUnrecoverable
            << " previously blanked cells from ready same-zone neighbours; "
            << "retaining them as holes." << nl << endl;
    }
    if (qualifiedSelection)
    {
        label retries = 0, ordinaryRecovered = 0, ordinaryFailed = 0;
        forAll(allCellTypes, celli)
        {
            retries += allCellTypes[celli] == INTERPOLATED
                && !sameZoneFallback[celli]
                && allDonorID[celli] != originalDonorID[celli];
            if
            (
                donorFallback_ == "sameZone"
             && ordinaryRequested[celli] && ordinaryCandidates[celli]
            )
            {
                ordinaryRecovered += allCellTypes[celli] == INTERPOLATED;
                ordinaryFailed += allCellTypes[celli] != INTERPOLATED;
            }
        }
        reduce(retries, sumOp<label>());
        reduce(ordinaryRecovered, sumOp<label>());
        reduce(ordinaryFailed, sumOp<label>());
        Info<< "cellVolumeWeight optional donor fallback: " << retries
            << " alternate-zone retries, " << ordinaryRecovered
            << " ordinary same-zone recoveries, " << ordinaryFailed
            << " unrecoverable ordinary receivers." << endl;
    }

    forAll(allCellTypes, celli)
    {
        if (allCellTypes[celli] != INTERPOLATED)
        {
            allWeight[celli] = 0.0;
        }
    }

    const scalar policySeconds = policyTimer.elapsedTime();
    const clockTime affineTimer;
    if (weightingScheme_ == "linearExact")
    {
        label nAffineReceivers = 0;
        label nAffineExact = 0;
        labelListList compact(allStencil);
        List<Map<label>> compactMap;
        mapDistribute affineMap(globalCells, compact, compactMap);
        pointField centres(mesh_.cellCentres());
        affineMap.distribute(centres);

        labelListList neighbours;
        pointListList neighbourCentres;
        labelListList neighbourZones;
        labelListList compactNeighbours;
        pointField expansionCentres;
        if (stencilExpansionRings_)
        {
            labelList needed(affineMap.constructSize(), Zero);
            for (const auto& stencil : compact)
            {
                for (const label slot : stencil)
                {
                    needed[slot] = 1;
                }
            }
            mapDistributeBase::distribute
            (
                UPstream::commsTypes::nonBlocking,
                UList<labelPair>::null(),
                mesh_.nCells(),
                affineMap.constructMap(),
                false,
                affineMap.subMap(),
                false,
                needed,
                label(0),
                maxEqOp<label>(),
                flipOp(),
                UPstream::msgType(),
                affineMap.comm()
            );
            const boolList ready
            (
                recoveryDonorMask
                (
                    allCellTypes, previouslyBlanked, protectedHoles,
                    reactivationCandidates, allowInterpolatedDonors_
                )
            );
            DynamicList<label> selected;
            forAll(needed, celli)
            {
                if (needed[celli])
                {
                    selected.append(celli);
                }
            }
            globalCellCells
            (
                globalCells, mesh_, ready, selected,
                neighbours, neighbourCentres, &zoneID, &neighbourZones
            );
            affineMap.distribute(neighbours);
            affineMap.distribute(neighbourZones);
            // Use canonical centres for the physical IDs sampled by the
            // stencil, not coupled-patch transformed neighbour images.
            compactNeighbours = neighbours;
            List<Map<label>> expansionCompactMap;
            mapDistribute expansionMap
            (
                globalCells, compactNeighbours, expansionCompactMap
            );
            expansionCentres = mesh_.cellCentres();
            expansionMap.distribute(expansionCentres);
        }
        forAll(allCellTypes, celli)
        {
            if (allCellTypes[celli] != INTERPOLATED)
            {
                continue;
            }
            ++nAffineReceivers;
            pointList donors(UIndirectList<point>(centres, compact[celli]));
            const labelList savedCounts(weightFallbackCounts_);
            if (affineWeights(mesh_.cellCentres()[celli], donors, allWeights[celli]))
            {
                ++nAffineExact;
                continue;
            }
            if (!stencilExpansionRings_ || sameZoneFallback[celli])
            {
                continue;
            }
            DynamicList<label> expanded(allStencil[celli]);
            labelHashSet seen(allStencil[celli]);
            DynamicList<point> expandedCentres(donors);
            for (const label slot : compact[celli])
            {
                forAll(neighbours[slot], j)
                {
                    if
                    (
                        neighbourZones[slot][j] == allDonorID[celli]
                     && seen.insert(neighbours[slot][j])
                    )
                    {
                        expanded.append(neighbours[slot][j]);
                        expandedCentres.append
                        (
                            expansionCentres[compactNeighbours[slot][j]]
                        );
                    }
                }
            }
            if (expanded.size() != allStencil[celli].size())
            {
                scalarList trial(allWeights[celli]);
                trial.setSize(expanded.size(), 0);
                weightFallbackCounts_ = savedCounts;
                const pointList trialCentres(expandedCentres);
                if (affineWeights(mesh_.cellCentres()[celli], trialCentres, trial))
                {
                    ++nAffineExact;
                    allStencil[celli] = labelList(expanded);
                    allWeights[celli].transfer(trial);
                }
            }
        }
        reportWeightFallbacks(nAffineExact, nAffineReceivers-nAffineExact);
    }
    const scalar affineSeconds = affineTimer.elapsedTime();

    // Write to volField for debugging
    if ((debug & 2) && mesh_.time().writeTime())
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allCellTypes_final", allCellTypes)
        );
        tfld().write();
    }


    cellTypes_.transfer(allCellTypes);
    cellStencil_.transfer(allStencil);
    cellInterpolationWeights_.transfer(allWeights);
    cellInterpolationWeight_.transfer(allWeight);
    //cellInterpolationWeight_.correctBoundaryConditions();
    oversetFvMeshBase::correctBoundaryConditions
    <
        volScalarField,
        oversetFvPatchField<scalar>,
        false
    >(cellInterpolationWeight_.boundaryFieldRef());

    DynamicList<label> interpolationCells;
    forAll(cellStencil_, cellI)
    {
        if (cellStencil_[cellI].size())
        {
            interpolationCells.append(cellI);
        }
    }
    interpolationCells_.transfer(interpolationCells);


    List<Map<label>> compactMap;
    cellInterpolationMap_.reset
    (
        new mapDistribute(globalCells, cellStencil_, compactMap)
    );

    // Dump interpolation stencil
    if ((debug & 2) && mesh_.time().writeTime())
    {
        // Dump weight
        cellInterpolationWeight_.instance() = mesh_.time().timeName();
        cellInterpolationWeight_.write();


        mkDir(mesh_.time().timePath());
        OBJstream str(mesh_.time().timePath()/"stencil2.obj");
        Info<< typeName << " : dumping to " << str.name() << endl;
        pointField cc(mesh_.cellCentres());
        cellInterpolationMap().distribute(cc);

        forAll(interpolationCells_, compactI)
        {
            label cellI = interpolationCells_[compactI];
            const labelList& slots = cellStencil_[cellI];

            Pout<< "cellI:" << cellI << " at:"
                << mesh_.cellCentres()[cellI]
                << " calculated from slots:" << slots
                << " cc:" << UIndirectList<point>(cc, slots)
                << " weights:" << cellInterpolationWeights_[cellI]
                << endl;

            forAll(slots, i)
            {
                if (cellInterpolationWeights_[cellI][i] > 0)
                {
                    const point& donorCc = cc[slots[i]];
                    const point& accCc = mesh_.cellCentres()[cellI];
                    str.writeLine(accCc, 0.1*accCc+0.9*donorCc);
                }
            }
        }
    }

    if (debug)
    {
        scalar overlapMax = overlapSeconds;
        scalar updateMax = updateTimer.elapsedTime();
        reduce(overlapMax, maxOp<scalar>());
        reduce(updateMax, maxOp<scalar>());
        Info<< "cellVolumeWeight maximum-rank times: overlap " << overlapMax
            << " s, update " << updateMax << " s" << endl;
        scalarField phases
        ({
            referenceSeconds, regionSeconds, policySeconds, affineSeconds,
            frontMapBuildSeconds, frontMapDistributionSeconds
        });
        Pstream::listReduce(phases, maxOp<scalar>());
        reduce(maxFrontHalo, sumOp<label>());
        Info<< "cellVolumeWeight phase maxima: reference " << phases[0]
            << " s, regions " << phases[1] << " s, donor policy " << phases[2]
            << " s, affine " << phases[3] << " s, front map build "
            << phases[4] << " s, front distribution " << phases[5]
            << " s, summed peak front halo slots " << maxFrontHalo << endl;
    }
    Info << this->info();

    return true;
}


void Foam::cellCellStencils::cellVolumeWeight::stencilWeights
(
    const point& sample,
    const pointList& donorCcs,
    scalarList& weights
) const
{
    distanceWeights(sample, donorCcs, weights);
    if (weightingScheme_ == "linearExact")
    {
        affineWeights(sample, donorCcs, weights);
    }
}


void Foam::cellCellStencils::cellVolumeWeight::distanceWeights
(
    const point& sample,
    const pointList& donorCcs,
    scalarList& weights
)
{
    weights.setSize(donorCcs.size());
    scalar nearest = std::numeric_limits<scalar>::infinity();
    forAll(donorCcs, i)
    {
        const vector delta = sample-donorCcs[i];
        const scalar d =
            std::hypot(std::hypot(delta.x(), delta.y()), delta.z());
        if (!std::isfinite(d))
        {
            FatalErrorInFunction
                << "Nonfinite recovery donor distance at " << sample
                << exit(FatalError);
        }
        if (d == 0)
        {
            weights = 0.0;
            weights[i] = 1.0;
            return;
        }
        weights[i] = d;
        nearest = min(nearest, d);
    }
    for (scalar& weight : weights)
    {
        weight = nearest/weight;
    }
    if (!normaliseVolumeWeights(weights))
    {
        FatalErrorInFunction
            << "No valid same-zone recovery weights at " << sample
            << exit(FatalError);
    }
}


// ************************************************************************* //
