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
#include "points0MotionSolver.H"
#include "pointIOField.H"

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
    if (affineScaling_ == "isotropic")
    {
        // One length scale for all axes keeps the minimum-norm correction
        // rotation invariant
        const scalar isoScale = cmptMax(scale);
        scale = vector::uniform
        (
            isoScale
         <= 128*std::numeric_limits<scalar>::epsilon()*cmptMax(coordinateScale)
          ? scalar(1) : isoScale
        );
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
    forAll(centres, j)
    {
        A[j][0] = 1;
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            A[j][cmpt+1] = (centres[j][cmpt]-sample[cmpt])/scale[cmpt];
        }
    }
    const scalar epsilon = std::numeric_limits<scalar>::epsilon();

    // Correction delta with A^T (w + delta) = e_0 minimising
    // sum delta_j^2/s_j^2 (s_j = 1: minimum norm; s_j^2 = w0_j + eps:
    // weighted). Rows with s_j = 0 are frozen.
    const bool weighted = affineCorrection_ != "minNorm";
    scalarList rowScale(centres.size(), scalar(1));
    if (weighted)
    {
        const scalar floor = 1e-6*max(weights);
        forAll(rowScale, j)
        {
            rowScale[j] = std::sqrt(max(weights[j], scalar(0)) + floor);
        }
    }
    auto project = [&](const scalarList& base, scalarList& result) -> bool
    {
        FixedList<long double, 4> residual({1, 0, 0, 0});
        RectangularMatrix<scalar> B(A);
        forAll(base, j)
        {
            for (label k = 0; k < 4; ++k)
            {
                residual[k] -= static_cast<long double>(base[j])*A[j][k];
                B[j][k] *= rowScale[j];
            }
        }
        SVD svd(B, affineSvdTolerance_);
        if (!svd.converged())
        {
            ++weightFallbackCounts_[0];
            return false;
        }
        const auto inverse = svd.VSinvUt();
        result = base;
        forAll(result, j)
        {
            long double delta = 0;
            for (label k = 0; k < 4; ++k)
            {
                delta += static_cast<long double>(inverse[k][j])*residual[k];
            }
            result[j] += scalar(delta*rowScale[j]);
            if (!std::isfinite(result[j]))
            {
                ++weightFallbackCounts_[1];
                return false;
            }
        }
        return true;
    };

    scalarList corrected;
    if (!project(weights, corrected))
    {
        return false;
    }
    if (affineCorrection_ == "positive")
    {
        // Clip-and-reproject: freeze negative weights at zero and project
        // the remaining moment error onto the free donors
        scalarList base(weights);
        for (label iter = 0; iter < centres.size(); ++iter)
        {
            const scalar tolerance = weightTolerance(corrected);
            label nClipped = 0;
            forAll(corrected, j)
            {
                if (rowScale[j] > 0 && corrected[j] < -tolerance)
                {
                    rowScale[j] = 0;
                    ++nClipped;
                }
            }
            if (!nClipped)
            {
                break;
            }
            forAll(base, j)
            {
                base[j] = rowScale[j] > 0 ? max(corrected[j], scalar(0)) : 0;
            }
            if (!project(base, corrected))
            {
                return false;
            }
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
            sumAbs > maxWeightAmplification_
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
    if (geometryMismatchReported_ && !debug)
    {
        // Lazy evaluation: only cells involved in overset coupling are
        // decomposed (see geometricVolume)
        geometricVolumes_ = -1;
        geometricRoundoff_ = -1;
        return;
    }
    // First update (or debug): evaluate and check all cells once
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
    geometryMismatchReported_ = true;
}


void Foam::cellCellStencils::cellVolumeWeight::evaluateGeometry
(
    const label celli
) const
{
    const auto reference = tetOverlapVolume::cellVolumeMinDecomp(mesh_, celli);
    geometricVolumes_[celli] = reference.first();
    geometricRoundoff_[celli] = reference.second();
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
        retainedVolume/geometricVolume(celli);
    const scalar allowance = geometricRoundoff(celli);
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
    if (donorSwitchMargin_ > 0 && zone != bestZone)
    {
        // Hysteresis: keep the previous donor zone unless another zone is
        // better by more than the margin
        const label previous = previousDonorZone(celli);
        if (zone == previous)
        {
            return coverage + donorSwitchMargin_ >= bestCoverage - tie;
        }
        if (bestZone == previous)
        {
            return coverage > bestCoverage + donorSwitchMargin_ + tie;
        }
    }
    if (donorSelection_ == "coverage")
    {
        if (coverage > bestCoverage+tie) return true;
        if (coverage < bestCoverage-tie) return false;
    }
    return preferZone(zoneID()[celli], zone, bestZone);
}


bool Foam::cellCellStencils::cellVolumeWeight::preferZone
(
    const label acceptorZone,
    const label zone,
    const label currentZone
) const
{
    if (donorPriority_.size())
    {
        // Lower position = higher priority; unlisted zones last
        const label rank = donorPriority_.find(zone);
        const label currentRank = donorPriority_.find(currentZone);
        if (rank != currentRank)
        {
            if (rank < 0) return false;
            if (currentRank < 0) return true;
            return rank < currentRank;
        }
    }
    const label diff = mag(acceptorZone-zone);
    const label currentDiff = mag(acceptorZone-currentZone);
    return diff < currentDiff || (diff == currentDiff && zone > currentZone);
}


Foam::scalar Foam::cellCellStencils::cellVolumeWeight::weightedDonorVolume
(
    const pointList& centres,
    const scalarList& volumes,
    scalarList weights
) const
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
    for (label i = 1; i < order.size() && !ambiguityReported_; ++i)
    {
        const label a = order[i-1], b = order[i];
        if (centres[a] == centres[b])
        {
            WarningInFunction
                << "Coincident donor centres in weighted "
                << "layer scoring: physical donor identities are ambiguous; "
                << "global labels cannot resolve this reproducibly."
                << nl << endl;
            ambiguityReported_ = true;
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
        // Report the events of the converged (final) round only
        reasons = Zero;
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
                    else if (donors[j] >= 0 && donors[j] < globalCells.totalSize())
                    {
                        // Remote slot directly from the map construction
                        const auto iter = compactMap
                        [
                            globalCells.whichProcID(donors[j])
                        ].cfind(donors[j]);
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


Foam::boolList Foam::cellCellStencils::cellVolumeWeight::anchoredCells
(
    const globalIndex& globalCells,
    const labelList& types,
    const labelListList& stencil,
    const scalarListList& weights,
    const scalarField& blend
) const
{
    // For fully interpolated rows (blend 1) the overset system contains
    // I - W_SS on any acceptor set S. If the stencils of S point only into S
    // (W_SS row sums 1) this block is singular. Rows with a positive-weight
    // path to a solved cell form a weakly chained diagonally dominant
    // M-matrix and are nonsingular.
    boolList anchored(mesh_.nCells(), true);
    labelListList compact(mesh_.nCells());
    forAll(types, celli)
    {
        if (types[celli] == INTERPOLATED && blend[celli] >= 1)
        {
            anchored[celli] = false;
            compact[celli] = stencil[celli];
        }
    }
    List<Map<label>> compactMap;
    mapDistribute map(globalCells, compact, compactMap);
    for (;;)
    {
        boolList compactAnchored(anchored);
        map.distribute(compactAnchored);
        label nNew = 0;
        forAll(anchored, celli)
        {
            if (anchored[celli])
            {
                continue;
            }
            const labelList& slots = compact[celli];
            const scalarList& w = weights[celli];
            const scalar tolerance = weightTolerance(w);
            for (label j = 0; j < min(slots.size(), w.size()); ++j)
            {
                if (mag(w[j]) > tolerance && compactAnchored[slots[j]])
                {
                    anchored[celli] = true;
                    ++nNew;
                    break;
                }
            }
        }
        if (!returnReduceOr(nNew))
        {
            break;
        }
    }
    return anchored;
}


bool Foam::cellCellStencils::cellVolumeWeight::repairUnanchored
(
    const globalIndex& globalCells,
    const boolList& anchored,
    const boolList& previouslyBlanked,
    const boolList& sameZoneFallback,
    const List<DynamicList<donorCandidate>>& candidates,
    const bool forceDemote,
    labelList& types,
    labelListList& stencil,
    scalarListList& weights,
    scalarField& blend,
    labelList& donorZones,
    labelList& counts
) const
{
    // Addressing: stencil donors first, then all candidate donors
    labelListList compact(mesh_.nCells());
    forAll(types, celli)
    {
        if (anchored[celli] || types[celli] != INTERPOLATED)
        {
            continue;
        }
        DynamicList<label> donors(stencil[celli]);
        if (!candidates.empty() && !sameZoneFallback[celli])
        {
            for (const auto& candidate : candidates[celli])
            {
                donors.append(candidate.donors);
            }
        }
        compact[celli] = labelList(donors);
    }
    List<Map<label>> compactMap;
    mapDistribute map(globalCells, compact, compactMap);
    boolList compactAnchored(anchored);
    map.distribute(compactAnchored);
    boolList compactReady(mesh_.nCells());
    forAll(compactReady, celli)
    {
        compactReady[celli] =
            types[celli] == CALCULATED && !previouslyBlanked[celli];
    }
    map.distribute(compactReady);

    bool demoted = false;
    forAll(types, celli)
    {
        if (anchored[celli] || types[celli] != INTERPOLATED)
        {
            continue;
        }
        const labelList& slots = compact[celli];
        const label nStencil = stencil[celli].size();
        bool repaired = false;

        if (!forceDemote && weights[celli].size() == nStencil)
        {
            // 1. Drop unanchored donors from the current (validated) stencil
            DynamicList<label> kept;
            DynamicList<scalar> keptWeights;
            for (label j = 0; j < nStencil; ++j)
            {
                if (compactAnchored[slots[j]] && weights[celli][j] > 0)
                {
                    kept.append(stencil[celli][j]);
                    keptWeights.append(weights[celli][j]);
                }
            }
            scalarList w(keptWeights);
            if (kept.size() && normaliseVolumeWeights(w))
            {
                stencil[celli] = labelList(kept);
                weights[celli].transfer(w);
                ++counts[0];
                repaired = true;
            }
        }

        if
        (
            !repaired && !forceDemote
         && !candidates.empty() && !sameZoneFallback[celli]
        )
        {
            // 2. Best candidate restricted to calculated (anchored) donors
            label offset = nStencil;
            scalar bestCoverage = -1;
            label bestZone = -1;
            labelList bestDonors;
            scalarList bestWeights;
            for (const auto& candidate : candidates[celli])
            {
                DynamicList<label> trial;
                DynamicList<scalar> trialVolumes;
                forAll(candidate.donors, j)
                {
                    const label slot = slots[offset+j];
                    if (compactReady[slot] && j < candidate.volumes.size())
                    {
                        trial.append(candidate.donors[j]);
                        trialVolumes.append(candidate.volumes[j]);
                    }
                }
                offset += candidate.donors.size();
                scalarList trialWeights(trialVolumes);
                scalar coverage = 0;
                if
                (
                    trial.size()
                 && candidateCoverage
                    (
                        celli, candidate.volumes, trialWeights, coverage
                    )
                 && normaliseVolumeWeights(trialWeights)
                 && betterCandidate
                    (
                        celli, candidate.zone, coverage, bestZone, bestCoverage
                    )
                )
                {
                    bestDonors = labelList(trial);
                    bestWeights.transfer(trialWeights);
                    bestCoverage = coverage;
                    bestZone = candidate.zone;
                }
            }
            if (bestZone >= 0)
            {
                stencil[celli].transfer(bestDonors);
                weights[celli].transfer(bestWeights);
                donorZones[celli] = bestZone;
                ++counts[1];
                repaired = true;
            }
        }

        if (!repaired)
        {
            // 3. Demote (solved cell or hole)
            types[celli] = previouslyBlanked[celli] ? HOLE : POROUS;
            stencil[celli].clear();
            weights[celli].clear();
            blend[celli] = 0;
            ++counts[2];
            demoted = true;
        }
    }
    return demoted;
}


Foam::label Foam::cellCellStencils::cellVolumeWeight::isolatedPockets
(
    const labelList& types,
    const scalarField& blend,
    boolList& isolated
) const
{
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    boolList isBlockedFace(mesh_.nFaces(), false);
    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        isBlockedFace[facei] =
            (types[own[facei]] == HOLE) != (types[nei[facei]] == HOLE);
    }
    labelList nbrTypes;
    syncTools::swapBoundaryCellList(mesh_, types, nbrTypes);
    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); ++facei)
    {
        isBlockedFace[facei] =
            (types[own[facei]] == HOLE)
         != (nbrTypes[facei-mesh_.nInternalFaces()] == HOLE);
    }
    const regionSplit cellRegion(mesh_, isBlockedFace);

    // 0: no solved cells, 1: solved cells only, 2: anchored
    labelList regionStatus(cellRegion.nRegions(), Zero);
    forAll(types, celli)
    {
        if (types[celli] == HOLE)
        {
            continue;
        }
        label& status = regionStatus[cellRegion[celli]];
        status = max
        (
            status,
            (types[celli] == INTERPOLATED && blend[celli] > 0) ? 2 : 1
        );
    }
    for (const fvPatch& fvp : mesh_.boundary())
    {
        if (isA<oversetFvPatch>(fvp) || fvPatch::constraintType(fvp.type()))
        {
            continue;
        }
        for (const label celli : fvp.faceCells())
        {
            if (types[celli] != HOLE)
            {
                regionStatus[cellRegion[celli]] = 2;
            }
        }
    }
    Pstream::listReduce(regionStatus, maxOp<label>());

    isolated.setSize(mesh_.nCells());
    label nIsolated = 0;
    forAll(types, celli)
    {
        isolated[celli] =
            types[celli] != HOLE && regionStatus[cellRegion[celli]] == 1;
        nIsolated += isolated[celli];
    }
    return returnReduce(nIsolated, sumOp<label>());
}


void Foam::cellCellStencils::cellVolumeWeight::findHoles
(
    const globalIndex& globalCells,
    const fvMesh& mesh,
    const labelList& zoneID,
    const labelListList& stencil,
    labelList& cellTypes,
    const boolList& holeSeeds
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

    // User hole points: 3 = forced hole (overrides real patches, e.g. for
    // bodies crossing the domain boundary, and interpolation reachability)
    forAll(holeSeeds, celli)
    {
        if (holeSeeds[celli] && cellTypes[celli] != HOLE)
        {
            regionType[cellRegion[celli]] = 3;
        }
    }

    // Now we've handled
    // - cells next to blocked cells
    // - coupled boundaries
    // Only thing to handle is the interpolation between regions


    labelListList compactStencil(stencil);
    List<Map<label>> compactMap;
    mapDistribute map(globalCells, compactStencil, compactMap);

    // Synchronise region status on processors once; afterwards only the
    // regions that became reachable (monotone change to 2) are exchanged.
    Pstream::listReduce(regionType, maxOp<label>());
    while (true)
    {
        DynamicList<label> changedRegions;

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

                    if (regionType[regionI] < 2)
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
                                changedRegions.append(regionI);
                                nChanged++;
                                break;
                            }
                        }
                    }
                }
            }
        }


        List<labelList> allChanged(Pstream::nProcs());
        allChanged[Pstream::myProcNo()] = std::move(changedRegions);
        Pstream::allGatherList(allChanged);
        label nTotalChanged = 0;
        for (const labelList& changed : allChanged)
        {
            for (const label regionI : changed)
            {
                regionType[regionI] = 2;
            }
            nTotalChanged += changed.size();
        }

        if (!nTotalChanged)
        {
            break;
        }
    }


    // See which regions have not been visited (regionType == 1) or are
    // forced holes (regionType == 3)
    label nForced = 0;
    forAll(cellRegion, cellI)
    {
        label type = regionType[cellRegion[cellI]];
        if ((type == 1 || type == 3) && cellTypes[cellI] != HOLE)
        {
            cellTypes[cellI] = HOLE;
            nForced += (type == 3);
        }
    }
    if (holeSeeds.size())
    {
        Info<< typeName << " : holePoints blanked "
            << returnReduce(nForced, sumOp<label>()) << " cells" << endl;
    }
}


Foam::label Foam::cellCellStencils::cellVolumeWeight::walkFrontLocal
(
    const globalIndex& globalCells,
    const scalar layerRelax,
    const labelListList& allStencil,
    labelList& allCellTypes,
    scalarField& allWeight,
    const scalarList& compactCellVol,
    const labelListList& compactStencil,
    const labelList& zoneID,
    const label holeLayers,
    const label useLayer,
    boolList* protectedHoles,
    const scalarList* donorVolumeRatio
) const
{
    if (!validWalkFrontLayerSettings(holeLayers, useLayer))
    {
        if (holeLayers < 1)
        {
            FatalErrorInFunction
                << "holeLayers: " << holeLayers
                << " must be greater than zero."
                << abort(FatalError);
        }

        FatalErrorInFunction
            << "useLayer: " << useLayer
            << " must be -1 (automatic) or in the range [1, "
            << holeLayers << "]."
            << abort(FatalError);
    }

    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();
    const cellList& cells = mesh_.cells();

    bitSet isFront(mesh_.nFaces());
    DynamicList<labelList> dataSet;
    DynamicList<scalarField> dAlllWeight;
    DynamicList<boolList> dProtectedHoles;
    DynamicList<scalar> daverageVolRatio;
    DynamicList<label> dNoDonor;
    bool haveVolRatioSamples = false;
    DynamicList<label> dHoles;

    const scalarField& V = mesh_.V();
    labelField nLayer(allCellTypes.size(), 0);

    for (label currLayer = 1; currLayer < holeLayers+1; ++currLayer)
    {
        isFront = false;
        setUpFront(allCellTypes, isFront);

        labelList allCellTypesWork(allCellTypes);

        bitSet isFrontWork(isFront);
        label nCurrLayer = currLayer;

        while (nCurrLayer > 1 && returnReduceOr(isFrontWork.any()))
        {
            bitSet newIsFront(mesh_.nFaces());
            for (const label facei : isFrontWork)
            {
                const label own = faceOwner[facei];
                if (allCellTypesWork[own] != HOLE)
                {
                    allCellTypesWork[own] = HOLE;
                    newIsFront.set(cells[own]);
                }
                if (mesh_.isInternalFace(facei))
                {
                    const label nei = faceNeighbour[facei];
                    if (allCellTypesWork[nei] != HOLE)
                    {
                        allCellTypesWork[nei] = HOLE;
                        newIsFront.set(cells[nei]);
                    }
                }
            }
            syncTools::syncFaceList(mesh_, newIsFront, orEqOp<unsigned int>());
            isFrontWork.transfer(newIsFront);
            nCurrLayer--;
        }

        boolList protectedHoleWork(allCellTypesWork.size(), false);
        forAll(protectedHoleWork, celli)
        {
            protectedHoleWork[celli] = allCellTypesWork[celli] == HOLE;
        }

        if ((debug & 2) && mesh_.time().writeTime())
        {
            tmp<volScalarField> tfld
            (
                createField
                (
                    mesh_,
                    "allCellTypesWork_Holes" + name(currLayer),
                    allCellTypesWork
                )
            );
            tfld().write();
        }

        if (currLayer == 1)
        {
            setUpFrontOnOversetPatch(allCellTypes, isFront);
        }
        if (currLayer > 1)
        {
            isFront = false;
            setUpFrontOnOversetPatch(allCellTypesWork, isFront);
            setUpFront(allCellTypesWork, isFront);
        }

        scalarField fraction(mesh_.nFaces(), Zero);
        scalarField volRatio(allCellTypes.size(), Zero);
        for (const label facei : isFront)
        {
            fraction[facei] = 1.0;
        }

        scalarField allWeightWork(allCellTypes.size(), Zero);
        bitSet nHoles(allCellTypes.size());
        bitSet porousWork(allCellTypes.size());
        label nNoDonor = 0;

        // Interpolate (or blank) a cell reached by the front
        auto visit = [&]
        (
            const label celli,
            const scalar faceFraction,
            bitSet& newIsFront,
            scalarField& newFraction
        )
        {
            if
            (
                allCellTypesWork[celli] == HOLE
             || porousWork.test(celli)
             || allWeightWork[celli] >= faceFraction
            )
            {
                return;
            }
            if (allStencil[celli].size())
            {
                nLayer[celli] = currLayer;
                allWeightWork[celli] = faceFraction;
                allCellTypesWork[celli] = INTERPOLATED;
                const label donorId = compactStencil[celli][0];
                volRatio[celli] = donorVolumeRatio
                    ? (*donorVolumeRatio)[celli]
                    : V[celli]/compactCellVol[donorId];
                seedCell
                (
                    celli, faceFraction-layerRelax, newIsFront, newFraction
                );
            }
            else if (stopIneligibleHoleGrowth_)
            {
                // No donors: keep solving the cell, do not grow the hole
                allWeightWork[celli] = 0.0;
                allCellTypesWork[celli] = POROUS;
                porousWork.set(celli);
                ++nNoDonor;
            }
            else
            {
                allWeightWork[celli] = 0.0;
                allCellTypesWork[celli] = HOLE;
                nHoles[celli] = true;
                ++nNoDonor;
                seedCell(celli, 1.0, newIsFront, newFraction);
            }
        };

        while (returnReduceOr(isFront.any()))
        {
            bitSet newIsFront(mesh_.nFaces());
            scalarField newFraction(fraction);

            for (const label facei : isFront)
            {
                visit(faceOwner[facei], fraction[facei], newIsFront, newFraction);
                if (mesh_.isInternalFace(facei))
                {
                    visit
                    (
                        faceNeighbour[facei], fraction[facei],
                        newIsFront, newFraction
                    );
                }
            }

            syncTools::syncFaceList(mesh_, newIsFront, orEqOp<unsigned int>());
            syncTools::syncFaceList(mesh_, newFraction, maxEqOp<scalar>());

            isFront.transfer(newIsFront);
            fraction.transfer(newFraction);
        }

        if ((debug & 2) && mesh_.time().writeTime())
        {
            tmp<volScalarField> tfld
            (
                createField
                (
                    mesh_,
                    "allCellTypesWork_Layers" + name(currLayer),
                    allCellTypesWork
                )
            );
            tfld().write();
        }

        label count = 0;
        forAll(volRatio, i)
        {
            if (volRatio[i] > 0)
            {
                count++;
            }
        }
        label nCount(count);
        reduce(nCount, sumOp<label>());

        const scalar aveVol = mag(gSum(volRatio));

        label nTotalHoles(nHoles.count());
        reduce(nTotalHoles, sumOp<label>());

        if
        (
            useLayer == -1
         && exceedsHoleGrowthLimit(currLayer, nTotalHoles, dHoles)
        )
        {
            break;
        }

        dAlllWeight.append(allWeightWork);
        dataSet.append(allCellTypesWork);
        dProtectedHoles.append(protectedHoleWork);
        dHoles.append(nTotalHoles);
        dNoDonor.append(returnReduce(nNoDonor, sumOp<label>()));

        if (nCount)
        {
            daverageVolRatio.append(aveVol/nCount);
            haveVolRatioSamples = true;
        }
        else
        {
            daverageVolRatio.append(GREAT);
        }
    }

    if ((debug & 2) && mesh_.time().writeTime())
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "walkFront_layers", nLayer)
        );
        tfld().write();
    }

    scalarList averageVolRatio;
    averageVolRatio.transfer(daverageVolRatio);
    label minVolId =
        selectLayerByVolumeRatio(averageVolRatio, haveVolRatioSamples);

    if (useLayer > -1)
    {
        Info<< nl << " Number of layers : " << averageVolRatio.size() << nl
            << " Average volumetric ratio : " << averageVolRatio << nl
            << " Number of holes cells : " << dHoles << nl
            << " Using layer : " << useLayer <<  nl << endl;
    }

    if (useLayer != -1)
    {
        minVolId = useLayer - 1;
    }

    allCellTypes.transfer(dataSet[minVolId]);
    allWeight.transfer(dAlllWeight[minVolId]);
    if (protectedHoles)
    {
        *protectedHoles = dProtectedHoles[minVolId];
    }
    return dNoDonor[minVolId];
}


void Foam::cellCellStencils::cellVolumeWeight::holePointCells
(
    const PtrList<fvMeshSubset>& meshParts,
    const labelList& zoneID,
    boolList& seeds
)
{
    seeds.setSize(mesh_.nCells());
    seeds = false;

    const pointField& points = mesh_.points();
    if (!holePoints0_ || holePoints0_().size() != points.size())
    {
        // (Re)read reference points, e.g. after a topology change
        holePoints0_.reset
        (
            new pointField
            (
                pointIOField(points0MotionSolver::points0IO(mesh_))
            )
        );
        if (holePoints0_().size() != points.size())
        {
            WarningInFunction
                << "points0 size " << holePoints0_().size()
                << " differs from mesh size " << points.size()
                << "; holePoints are used untransformed." << nl << endl;
            holePoints0_.reset(new pointField(points));
        }
    }
    const pointField& points0 = holePoints0_();

    label nNotFound = 0;
    for (const auto& hp : holePoints_)
    {
        const label zonei = hp.first();
        if (zonei < 0 || zonei >= meshParts.size())
        {
            continue;
        }
        const labelList& pointMap = meshParts[zonei].pointMap();

        // Rigid least-squares fit (Kabsch) of zone points0 -> points.
        // Processor-shared points are counted once per rank on both sides,
        // which leaves the fit of a rigid motion exact.
        vector c0(Zero), c(Zero);
        bool moved = false;
        for (const label pointi : pointMap)
        {
            c0 += points0[pointi];
            c += points[pointi];
            moved = moved || points0[pointi] != points[pointi];
        }
        const label n = returnReduce(pointMap.size(), sumOp<label>());
        reduce(c0, sumOp<vector>());
        reduce(c, sumOp<vector>());
        point x = hp.second();
        if (returnReduceOr(moved) && n)
        {
            c0 /= n;
            c /= n;
            tensor H(Zero);
            for (const label pointi : pointMap)
            {
                H += (points0[pointi]-c0)*(points[pointi]-c);
            }
            reduce(H, sumOp<tensor>());

            scalarRectangularMatrix M(3, 3);
            for (direction i = 0; i < 3; ++i)
            {
                for (direction j = 0; j < 3; ++j)
                {
                    M(i, j) = H[3*i+j];
                }
            }
            const SVD svd(M);
            tensor U(Zero), V(Zero);
            for (direction i = 0; i < 3; ++i)
            {
                for (direction j = 0; j < 3; ++j)
                {
                    U[3*i+j] = svd.U()(i, j);
                    V[3*i+j] = svd.V()(i, j);
                }
            }
            // R = V diag(1,1,d) U^T, d corrects reflections
            const scalar d = det(V & U.T()) < 0 ? -1 : 1;
            label nSmallest = 0;
            for (label i = 1; i < 3; ++i)
            {
                if (svd.S()[i] < svd.S()[nSmallest]) nSmallest = i;
            }
            tensor D(I);
            D[4*nSmallest] = d;
            const tensor R(V & D & U.T());
            x = (R & (x - c0)) + c;
        }

        label nFound = 0;
        forAll(meshParts, otherZonei)
        {
            if (otherZonei == zonei)
            {
                continue;
            }
            const label subCelli =
                meshParts[otherZonei].subMesh().findCell(x);
            if (subCelli >= 0)
            {
                const label celli = meshParts[otherZonei].cellMap()[subCelli];
                if (zoneID[celli] == otherZonei)
                {
                    seeds[celli] = true;
                    ++nFound;
                }
            }
        }
        nNotFound += !returnReduceOr(nFound);
    }
    if (nNotFound)
    {
        WarningInFunction
            << nNotFound << " holePoints are not inside any other zone."
            << nl << endl;
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


void Foam::cellCellStencils::cellVolumeWeight::patchOverlapVolumes
(
    const autoPtr<mapDistribute>& mapPtr,
    const labelListList& addressing,
    const scalarListList& weights,
    const labelList& patchTypes,
    scalarList& result
) const
{
    labelList work(patchTypes);
    if (mapPtr)
    {
        mapPtr().distribute(work);
    }
    result.setSize(addressing.size());
    forAll(addressing, celli)
    {
        const labelList& slots = addressing[celli];
        const scalarList& w = weights[celli];
        scalar v = 0;
        for (label i = 0; i < min(slots.size(), w.size()); ++i)
        {
            if (work[slots[i]] == PATCH)
            {
                v += w[i];
            }
        }
        result[celli] = v;
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
    List<DynamicList<donorCandidate>>* candidates,
    const scalarList& patchOverlap
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
                // Patch-patch interaction. Legacy: any overlap with a donor
                // cell next to a physical patch cuts the acceptor. This
                // removes all acceptors in a narrow gap between two bodies.
                bool cut = true;
                if (patchCutting_ != "legacy")
                {
                    const scalar v = geometricVolume(cellI);
                    const scalar allowance = geometricRoundoff(cellI);
                    // Hysteresis: former holes need a smaller violation
                    const scalar h =
                        cellTypes_[cellI] == HOLE ? patchCutHysteresis_ : 1;
                    if (patchCutting_ == "overlapFraction")
                    {
                        cut =
                            patchOverlap.size() > subCellI
                         && patchOverlap[subCellI]
                          > (h*patchOverlapFraction_ + allowance)*v;
                    }
                    else
                    {
                        // coverage: cut only cells extending beyond the
                        // donor domain (i.e. into the donor body)
                        cut =
                            sum(weights[subCellI])
                          < (1 - h*overlapTolerance_ - allowance)*v;
                    }
                }
                if (cut)
                {
                    allCellTypes[cellI] = HOLE;
                    validDonors = false;
                }
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

            //- Option 3 (or explicit donorPriority), with optional
            //  hysteresis towards the previous donor zone
            bool replace = allDonorID[cellI] == -1;
            if (!replace)
            {
                const label previous = previousDonorZone(cellI);
                if
                (
                    donorSwitchMargin_ > 0
                 && (previous == donorZoneID || previous == allDonorID[cellI])
                 && donorZoneID != allDonorID[cellI]
                )
                {
                    const scalar v = geometricVolume(cellI);
                    const scalar thisCoverage = sum(weights[subCellI])/v;
                    const scalar currentCoverage = sum(allWeights[cellI])/v;
                    replace =
                        previous == donorZoneID
                      ? thisCoverage + donorSwitchMargin_ >= currentCoverage
                      : thisCoverage > currentCoverage + donorSwitchMargin_;
                }
                else
                {
                    replace =
                        preferZone(subZoneID, donorZoneID, allDonorID[cellI]);
                }
            }

            if (replace)
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
    checkAnchoring_(dict.getOrDefault("checkAnchoring", true)),
    checkPockets_(dict.getOrDefault<word>("checkPockets", "report")),
    patchCutting_(dict.getOrDefault<word>("patchCutting", "legacy")),
    patchOverlapFraction_
    (
        dict.getOrDefault<scalar>("patchOverlapFraction", 0)
    ),
    patchCutHysteresis_(dict.getOrDefault<scalar>("patchCutHysteresis", 1)),
    holePoints_
    (
        dict.getOrDefault<List<Tuple2<label, point>>>
        (
            "holePoints", List<Tuple2<label, point>>()
        )
    ),
    holePoints0_(nullptr),
    stopIneligibleHoleGrowth_
    (
        dict.getOrDefault("stopIneligibleHoleGrowth", false)
    ),
    donorPriority_(dict.getOrDefault<labelList>("donorPriority", labelList())),
    donorSwitchMargin_(dict.getOrDefault<scalar>("donorSwitchMargin", 0)),
    donorZone_(),
    affineCorrection_
    (
        dict.getOrDefault<word>("affineCorrection", "minNorm")
    ),
    maxWeightAmplification_
    (
        dict.getOrDefault<scalar>
        (
            "maxWeightAmplification",
            1/std::sqrt(std::numeric_limits<scalar>::epsilon())
        )
    ),
    affineSvdTolerance_
    (
        dict.getOrDefault<scalar>
        (
            "affineSvdTolerance",
            100*std::numeric_limits<scalar>::epsilon()
        )
    ),
    affineScaling_(dict.getOrDefault<word>("affineScaling", "perAxis")),
    geometricVolumes_(mesh.nCells()),
    geometricRoundoff_(mesh.nCells()),
    weightFallbackCounts_(3, Zero),
    weightFallbackReported_(false),
    ambiguityReported_(false),
    geometryMismatchReported_(false),
    thinMeshReported_(false),
    cacheStaticPairs_(dict.getOrDefault("cacheStaticPairs", false)),
    skipStaticUpdate_(dict.getOrDefault("skipStaticUpdate", false)),
    pairOrientation_
    (
        dict.getOrDefault<word>("pairOrientation", "zoneOrder")
    ),
    overlapShortcut_(dict.getOrDefault("overlapShortcut", false)),
    pairCache_(),
    pointZone_(),
    lastPoints_(),
    haveStencil_(false)
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
     || (checkPockets_ != "none" && checkPockets_ != "report"
         && checkPockets_ != "fix")
     || (patchCutting_ != "legacy" && patchCutting_ != "overlapFraction"
         && patchCutting_ != "coverage")
     || !std::isfinite(patchOverlapFraction_)
     || patchOverlapFraction_ < 0 || patchOverlapFraction_ > 1
     || !std::isfinite(patchCutHysteresis_)
     || patchCutHysteresis_ <= 0 || patchCutHysteresis_ > 1
     || !std::isfinite(donorSwitchMargin_) || donorSwitchMargin_ < 0
     || (affineCorrection_ != "minNorm" && affineCorrection_ != "weighted"
         && affineCorrection_ != "positive")
     || !std::isfinite(maxWeightAmplification_)
     || maxWeightAmplification_ <= 0
     || !std::isfinite(affineSvdTolerance_)
     || affineSvdTolerance_ <= 0 || affineSvdTolerance_ >= 1
     || (affineScaling_ != "perAxis" && affineScaling_ != "isotropic")
     || (pairOrientation_ != "zoneOrder"
         && pairOrientation_ != "smallerTarget")
    )
    {
        FatalIOErrorInFunction(dict_)
            << "Invalid cellVolumeWeight configuration: tolerances must be "
            << "finite fractions in [0,1]; donorSelection is legacy|coverage, "
            << "donorFallback is legacy|alternateZones|sameZone, "
            << "weightingScheme is volume|linearExact, layerDonorMetric is "
            << "legacy|overlapWeighted, processorMapMethod is AABB|LOD, "
            << "stencilExpansionRings is 0 or 1, checkPockets is "
            << "none|report|fix, patchCutting is "
            << "legacy|overlapFraction|coverage, patchOverlapFraction is in "
            << "[0,1], patchCutHysteresis is in (0,1], donorSwitchMargin >= 0, "
            << "affineCorrection is minNorm|weighted|positive, "
            << "maxWeightAmplification > 0, affineSvdTolerance is in (0,1), "
            << "affineScaling is perAxis|isotropic and pairOrientation is "
            << "zoneOrder|smallerTarget." << exit(FatalIOError);
    }
    {
        // donorPriority: distinct non-negative zones (range checked in
        // update once the number of zones is known)
        labelHashSet seen;
        for (const label zonei : donorPriority_)
        {
            if (zonei < 0 || !seen.insert(zonei))
            {
                FatalIOErrorInFunction(dict_)
                    << "donorPriority " << donorPriority_
                    << " must list distinct non-negative zones."
                    << exit(FatalIOError);
            }
        }
        for (const auto& hp : holePoints_)
        {
            if
            (
                hp.first() < 0
             || !std::isfinite(hp.second().x())
             || !std::isfinite(hp.second().y())
             || !std::isfinite(hp.second().z())
            )
            {
                FatalIOErrorInFunction(dict_)
                    << "Invalid holePoints entry " << hp
                    << ": expected (zone (x y z))." << exit(FatalIOError);
            }
        }
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

    // Zone motion detection (bitwise point comparison per zone)
    boolList zoneMoved(nZones, true);
    if (cacheStaticPairs_ || skipStaticUpdate_)
    {
        const pointField& points = mesh_.points();
        const bool rebuild = returnReduceOr
        (
            pointZone_.size() != points.size()
         || lastPoints_.size() != points.size()
         || mesh_.topoChanging()
        );
        if (rebuild)
        {
            pointZone_.setSize(points.size());
            pointZone_ = -1;
            const labelListList& cellPoints = mesh_.cellPoints();
            forAll(cellPoints, celli)
            {
                for (const label pointi : cellPoints[celli])
                {
                    label& pz = pointZone_[pointi];
                    pz = (pz == -1 || pz == zoneID[celli]) ? zoneID[celli] : -2;
                }
            }
            pairCache_.clear();
            haveStencil_ = false;
        }
        else
        {
            labelList moved(nZones, Zero);
            forAll(points, pointi)
            {
                if (points[pointi] != lastPoints_[pointi])
                {
                    const label pz = pointZone_[pointi];
                    if (pz >= 0)
                    {
                        moved[pz] = 1;
                    }
                    else if (pz == -2)
                    {
                        moved = 1;
                    }
                }
            }
            Pstream::listReduce(moved, maxOp<label>());
            forAll(zoneMoved, zonei)
            {
                zoneMoved[zonei] = moved[zonei];
            }
        }
        lastPoints_ = points;

        if
        (
            skipStaticUpdate_ && haveStencil_ && !rebuild
         && !zoneMoved.found(true)
         && cellTypes_.size() == mesh_.nCells()
        )
        {
            Info<< typeName << " : no zone moved; keeping the previous "
                << "stencil" << endl;
            return false;
        }
    }

    Info<< typeName << " : detected " << nZones
        << " mesh regions" << nl << endl;

    if (donorPriority_.size())
    {
        boolList listed(nZones, false);
        for (const label zonei : donorPriority_)
        {
            if (zonei < 0 || zonei >= nZones || listed[zonei])
            {
                FatalIOErrorInFunction(dict_)
                    << "donorPriority " << donorPriority_
                    << " must list distinct zones in [0," << nZones-1 << "]."
                    << exit(FatalIOError);
            }
            listed[zonei] = true;
        }
    }


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


    if (!thinMeshReported_)
    {
        thinMeshReported_ = true;
        // Cells touching both an overset and a physical patch: donors for
        // such zones are next to walls (PATCH) and their acceptors get cut;
        // a zone only one cell thick has no interior donor layer.
        boolList onOverset(mesh_.nCells(), false);
        boolList onPhysical(mesh_.nCells(), false);
        for (const fvPatch& fvp : mesh_.boundary())
        {
            if (isA<oversetFvPatch>(fvp))
            {
                UIndirectList<bool>(onOverset, fvp.faceCells()) = true;
            }
            else if (!fvPatch::constraintType(fvp.type()))
            {
                UIndirectList<bool>(onPhysical, fvp.faceCells()) = true;
            }
        }
        labelList nThin(nZones, Zero);
        forAll(onOverset, celli)
        {
            nThin[zoneID[celli]] += onOverset[celli] && onPhysical[celli];
        }
        Pstream::listReduce(nThin, sumOp<label>());
        forAll(nThin, zonei)
        {
            if (nThin[zonei])
            {
                WarningInFunction
                    << "Zone " << zonei << " has " << nThin[zonei]
                    << " cells touching both an overset and a physical "
                    << "patch (thin donor mesh). Such cells are treated as "
                    << "overset cells; their donors next to walls cut "
                    << "acceptors. Add cell layers between the overset "
                    << "boundary and the walls." << nl << endl;
            }
        }
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

    // Overlap data of one side of a zone pair: acceptor addressing into the
    // (distributed) donor zone, overlap volumes, global donor cells,
    // interpolated donor patch types and PATCH-donor overlap volumes.
    auto computeSide = [&]
    (
        const autoPtr<mapDistribute>& donorMap,
        const labelListList& addressing,
        const scalarListList& weights,
        const labelList& donorCellMap,
        const label nAcceptorCells,
        pairSide& side
    )
    {
        side.addressing = addressing;
        side.weights = weights;
        side.patchTypes.setSize(nAcceptorCells);
        side.patchTypes = -1;
        interpolatePatchTypes
        (
            donorMap,
            addressing,
            labelList(labelUIndList(allPatchTypes, donorCellMap)),
            side.patchTypes
        );
        side.patchOverlap.clear();
        if (patchCutting_ == "overlapFraction")
        {
            patchOverlapVolumes
            (
                donorMap,
                addressing,
                weights,
                labelList(labelUIndList(allPatchTypes, donorCellMap)),
                side.patchOverlap
            );
        }
        // Donor cell labels in global cell indexing (on overall mesh)
        side.donors.setSize(donorCellMap.size());
        forAll(donorCellMap, i)
        {
            side.donors[i] = globalCells.toGlobal(donorCellMap[i]);
        }
        if (donorMap)
        {
            donorMap->distribute(side.donors);
        }
    };

    if (cacheStaticPairs_ && pairCache_.size() != nZones*nZones)
    {
        pairCache_.clear();
        pairCache_.setSize(nZones*nZones);
    }
    label nCachedPairs = 0;

    for (label srcI = 0; srcI < meshParts.size()-1; srcI++)
    {
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

            // Pair orientation: meshToMesh builds the search tree on the
            // target and loops over source cells
            const bool swapPair =
                pairOrientation_ == "smallerTarget"
             && nCellsPerZone[tgtI] > nCellsPerZone[srcI];
            const label aI = swapPair ? tgtI : srcI;
            const label bI = swapPair ? srcI : tgtI;

            const fvMesh& aMesh = meshParts[aI].subMesh();
            const labelList& aCellMap = meshParts[aI].cellMap();
            const fvMesh& bMesh = meshParts[bI].subMesh();
            const labelList& bCellMap = meshParts[bI].cellMap();

            const label pairI = srcI*nZones + tgtI;
            const bool usable =
                cacheStaticPairs_
             && pairCache_.set(pairI)
             && !zoneMoved[srcI] && !zoneMoved[tgtI];

            autoPtr<FixedList<pairSide, 2>> localSides;
            if (!usable)
            {
                localSides.reset(new FixedList<pairSide, 2>());
                FixedList<pairSide, 2>& sides = localSides();

                const bool oldShortcut =
                    cellVolumeWeightMethod::overlapShortcut;
                cellVolumeWeightMethod::overlapShortcut = overlapShortcut_;
                meshToMesh mapper
                (
                    aMesh,
                    bMesh,
                    meshToMesh::interpolationMethod::imCellVolumeWeight,
                    HashTable<word>(0),     // patchMap,
                    wordList(0),            // cuttingPatches
                    processorMapMethod_ == "LOD"
                  ? meshToMesh::procMapMethod::pmLOD
                  : meshToMesh::procMapMethod::pmAABB,
                    false                   // do not normalise
                );
                cellVolumeWeightMethod::overlapShortcut = oldShortcut;

                // Side 0: acceptors in a, donors in b
                computeSide
                (
                    mapper.tgtMap(),
                    mapper.srcToTgtCellAddr(),
                    mapper.srcToTgtCellWght(),
                    bCellMap,
                    aMesh.nCells(),
                    sides[0]
                );
                // Side 1: acceptors in b, donors in a
                computeSide
                (
                    mapper.srcMap(),
                    mapper.tgtToSrcCellAddr(),
                    mapper.tgtToSrcCellWght(),
                    aCellMap,
                    bMesh.nCells(),
                    sides[1]
                );
                if (cacheStaticPairs_)
                {
                    pairCache_.set(pairI, localSides.ptr());
                }
            }
            else
            {
                ++nCachedPairs;
            }
            const FixedList<pairSide, 2>& sides =
                localSides ? localSides() : pairCache_[pairI];

            combineCellTypes
            (
                aI,
                aMesh,
                aCellMap,

                bI,
                sides[0].addressing,
                sides[0].weights,
                sides[0].donors,
                sides[0].patchTypes,

                // Overall mesh data
                allStencil,
                allWeights,
                allCellTypes,
                allDonorID,
                candidates.empty() ? nullptr : &candidates,
                sides[0].patchOverlap
            );

            combineCellTypes
            (
                bI,
                bMesh,
                bCellMap,

                aI,
                sides[1].addressing,
                sides[1].weights,
                sides[1].donors,
                sides[1].patchTypes,

                // Overall mesh data
                allStencil,
                allWeights,
                allCellTypes,
                allDonorID,
                candidates.empty() ? nullptr : &candidates,
                sides[1].patchOverlap
            );
        }
    }
    if (cacheStaticPairs_)
    {
        Info<< typeName << " : reused " << nCachedPairs
            << " cached static zone pairs" << endl;
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
             || sum(allWeights[celli]) > overlapTolerance_*geometricVolume(celli)
            )
            {
                continue;
            }
            scalar bestCoverage = -1;
            label bestZone = -1;
            for (const auto& candidate : candidates[celli])
            {
                const scalar coverage =
                    sum(candidate.volumes)/geometricVolume(celli);
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
                    scalar v = geometricVolume(cellI);
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
    boolList holeSeeds;
    if (holePoints_.size())
    {
        holePointCells(meshParts, zoneID, holeSeeds);
    }
    findHoles
    (
        globalCells, mesh_, zoneID,
        retainCandidates ? holeStencil : allStencil, allCellTypes, holeSeeds
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
    label nFrontNoDonor = 0;
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
                 && sum(candidate.volumes) > overlapTolerance_*geometricVolume(celli)
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
        nFrontNoDonor = walkFrontLocal
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
                    // Legacy: first ready same-zone neighbour only (zeroth-
                    // order copy). Kept for reproducibility; it differs from
                    // the distance-weighted stencil used for overlapWeighted
                    // and for ordinary candidates below. The cell is
                    // interpolated only for this update (blend from the
                    // front walk) and re-evaluated on the next update.
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

    if (checkAnchoring_)
    {
        // Detect and repair closed interpolation cycles (singular I - W_SS)
        const label maxIter = 10;
        labelList repairs(3, Zero);
        label nUnanchored = 0;
        for (label iter = 0; iter < maxIter; ++iter)
        {
            const boolList anchored
            (
                anchoredCells
                (
                    globalCells, allCellTypes, allStencil, allWeights,
                    allWeight
                )
            );
            label n = 0;
            forAll(anchored, celli)
            {
                n += !anchored[celli];
            }
            if (!returnReduceOr(n))
            {
                break;
            }
            if (!iter)
            {
                nUnanchored = n;
            }
            const bool demoted = repairUnanchored
            (
                globalCells, anchored, previouslyBlanked, sameZoneFallback,
                candidates, iter == maxIter-1,
                allCellTypes, allStencil, allWeights, allWeight, allDonorID,
                repairs
            );
            if (returnReduceOr(demoted))
            {
                // Demoted donors invalidate dependent stencils
                finaliseDonors
                (
                    globalCells, previouslyBlanked, sameZoneFallback,
                    candidates, allCellTypes, allStencil, allWeights,
                    allDonorID,
                    donorFallback_ == "sameZone" ? &donorExclusions : nullptr
                );
            }
        }
        reduce(nUnanchored, sumOp<label>());
        Pstream::listReduce(repairs, sumOp<label>());
        if (nUnanchored)
        {
            const boolList anchored
            (
                anchoredCells
                (
                    globalCells, allCellTypes, allStencil, allWeights,
                    allWeight
                )
            );
            label remaining = 0;
            forAll(anchored, celli)
            {
                remaining += !anchored[celli];
            }
            reduce(remaining, sumOp<label>());
            WarningInFunction
                << "cellVolumeWeight interpolation cycles: " << nUnanchored
                << " fully interpolated cells without a donor path to a "
                << "solved cell. Repaired: " << repairs[0]
                << " by dropping unanchored donors, " << repairs[1]
                << " by re-selecting a candidate, " << repairs[2]
                << " demoted; " << remaining << " remaining."
                << " Consider allowInterpolatedDonors false or donorPriority."
                << nl << endl;
        }
    }

    if (checkPockets_ != "none")
    {
        // Solved regions without interpolation or real boundary have a
        // singular (pure Neumann) pressure system
        boolList isolated;
        const label nIsolated =
            isolatedPockets(allCellTypes, allWeight, isolated);
        if (nIsolated)
        {
            WarningInFunction
                << "cellVolumeWeight found " << nIsolated
                << " solved cells in isolated pockets (no interpolated cell "
                << "and no physical patch)"
                << (checkPockets_ == "fix" ? "; converting them to holes." : ".")
                << nl << endl;
        }
        if (nIsolated && checkPockets_ == "fix")
        {
            forAll(isolated, celli)
            {
                if (isolated[celli])
                {
                    allCellTypes[celli] = HOLE;
                    allStencil[celli].clear();
                    allWeights[celli].clear();
                    allWeight[celli] = 0;
                }
            }
            // Revalidate stencils that referenced the new holes
            finaliseDonors
            (
                globalCells, previouslyBlanked, sameZoneFallback, candidates,
                allCellTypes, allStencil, allWeights, allDonorID,
                donorFallback_ == "sameZone" ? &donorExclusions : nullptr
            );
            if (checkAnchoring_)
            {
                const boolList anchored
                (
                    anchoredCells
                    (
                        globalCells, allCellTypes, allStencil, allWeights,
                        allWeight
                    )
                );
                label remaining = 0;
                forAll(anchored, celli)
                {
                    remaining += !anchored[celli];
                }
                reduce(remaining, sumOp<label>());
                if (remaining)
                {
                    WarningInFunction
                        << "cellVolumeWeight: " << remaining
                        << " unanchored interpolated cells after pocket "
                        << "removal." << nl << endl;
                }
            }
        }
    }

    forAll(allCellTypes, celli)
    {
        if (allCellTypes[celli] != INTERPOLATED)
        {
            allWeight[celli] = 0.0;
        }
    }

    if (nFrontNoDonor)
    {
        Info<< typeName << " : " << nFrontNoDonor
            << " front cells without donors became "
            << (stopIneligibleHoleGrowth_ ? "porous (solved)" : "holes")
            << " in the interpolation-front walk" << endl;
    }

    const scalar policySeconds = policyTimer.elapsedTime();
    const clockTime affineTimer;
    const bool checkAffineAnchoring =
        checkAnchoring_ && weightingScheme_ == "linearExact";
    labelListList preAffineStencil(checkAffineAnchoring ? allStencil : labelListList());
    scalarListList preAffineWeights(checkAffineAnchoring ? allWeights : scalarListList());
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
    if (checkAffineAnchoring)
    {
        // Stencil expansion or signed weights may close a cycle: revert
        // such cells to their validated pre-affine weights.
        label nReverted = 0;
        for (label iter = 0; iter < 5; ++iter)
        {
            const boolList anchored
            (
                anchoredCells
                (
                    globalCells, allCellTypes, allStencil, allWeights,
                    allWeight
                )
            );
            label n = 0;
            forAll(anchored, celli)
            {
                if
                (
                    !anchored[celli]
                 && (
                        allStencil[celli] != preAffineStencil[celli]
                     || allWeights[celli] != preAffineWeights[celli]
                    )
                )
                {
                    allStencil[celli] = preAffineStencil[celli];
                    allWeights[celli] = preAffineWeights[celli];
                    ++n;
                }
            }
            nReverted += n;
            if (!returnReduceOr(n))
            {
                break;
            }
        }
        reduce(nReverted, sumOp<label>());
        if (nReverted)
        {
            Info<< "cellVolumeWeight anchoring: reverted " << nReverted
                << " affine stencils that closed an interpolation cycle."
                << endl;
        }
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


    if (donorSwitchMargin_ > 0)
    {
        donorZone_.setSize(mesh_.nCells());
        forAll(donorZone_, celli)
        {
            donorZone_[celli] =
                allCellTypes[celli] == INTERPOLATED ? allDonorID[celli] : -1;
        }
    }

    haveStencil_ = true;
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
