/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
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

#include "cellVolumeWeightMethod.H"
#include "addToRunTimeSelectionTable.H"
#include "treeDataCell.H"
#include "indexedOctree.H"
#include "tetOverlapVolume.H"
#include "clockTime.H"

#include <cmath>
#include <limits>

namespace Foam
{
    defineTypeNameAndDebug(cellVolumeWeightMethod, 0);
    addToRunTimeSelectionTable
    (
        meshToMeshMethod,
        cellVolumeWeightMethod,
        components
    );
}

bool Foam::cellVolumeWeightMethod::overlapShortcut = false;


bool Foam::cellVolumeWeightMethod::containsCell
(
    const polyMesh& outerMesh,
    const label outerCelli,
    const polyMesh& innerMesh,
    const label innerCelli
)
{
    const pointField& outerPoints = outerMesh.points();
    const pointField& innerPoints = innerMesh.points();
    const labelList& outerCellPoints = outerMesh.cellPoints()[outerCelli];
    const labelList& innerCellPoints = innerMesh.cellPoints()[innerCelli];
    const boundBox outerBox(outerMesh.cellBb(outerCelli));
    const scalar length = cmptMax(outerBox.span());
    const scalar tolerance = 1e3*SMALL*length;

    for (const label facei : outerMesh.cells()[outerCelli])
    {
        const face& f = outerMesh.faces()[facei];
        const bool flip = outerMesh.faceOwner()[facei] != outerCelli;
        const point& p0 = outerPoints[f[0]];
        for (label i = 1; i < f.size()-1; ++i)
        {
            vector n =
                (outerPoints[f[i]] - p0) ^ (outerPoints[f[i+1]] - p0);
            const scalar magN = mag(n);
            if (!(magN > 1e3*SMALL*sqr(length)))
            {
                // Degenerate fan triangle: no reliable half-space
                return false;
            }
            n /= (flip ? -magN : magN);
            // Convexity of the tessellated outer cell
            for (const label pointi : outerCellPoints)
            {
                if (((outerPoints[pointi] - p0) & n) > tolerance)
                {
                    return false;
                }
            }
            // Strict containment of the inner cell vertices
            for (const label pointi : innerCellPoints)
            {
                if (!(((innerPoints[pointi] - p0) & n) < -tolerance))
                {
                    return false;
                }
            }
        }
    }
    return true;
}


Foam::cellVolumeWeightMethod::cellVolumeWeightMethod
(
    const polyMesh& src,
    const polyMesh& tgt
)
:
    meshToMeshMethod(src, tgt),
    candidateTests_(0),
    discardedVolume_(0)
{}

Foam::cellVolumeWeightMethod::~cellVolumeWeightMethod()
{}

void Foam::cellVolumeWeightMethod::calculateAddressing
(
    labelListList& srcToTgtAddr,
    scalarListList& srcToTgtWght,
    pointListList& srcToTgtVec,
    labelListList& tgtToSrcAddr,
    scalarListList& tgtToSrcWght,
    pointListList& tgtToSrcVec,
    const bool centroids
)
{
    V_ = 0;
    candidateTests_ = 0;
    discardedVolume_ = 0;
    srcToTgtAddr = labelListList(src_.nCells());
    srcToTgtWght = scalarListList(src_.nCells());
    tgtToSrcAddr = labelListList(tgt_.nCells());
    tgtToSrcWght = scalarListList(tgt_.nCells());
    srcToTgtVec = pointListList(centroids ? src_.nCells() : 0);
    tgtToSrcVec = pointListList(centroids ? tgt_.nCells() : 0);
    if (!src_.nCells() || !tgt_.nCells())
    {
        return;
    }

    const clockTime timer;
    label roundoffVolumes = 0;
    label acceptedPairs = 0;
    label clampedVolumes = 0;
    scalar maxViolation = 0;
    // Reference volumes are evaluated only for source cells that are
    // intersected and, lazily, for target cells that are candidates
    // (-1: not evaluated)
    scalarField srcVol(src_.nCells(), -1), tgtVol(tgt_.nCells(), -1);
    scalarField srcRoundoff(src_.nCells(), -1), tgtRoundoff(tgt_.nCells(), -1);
    label mismatches = 0;
    auto evaluate = [&mismatches]
    (
        const polyMesh& mesh,
        const label celli,
        scalarField& vol,
        scalarField& roundoff
    )
    {
        if (vol[celli] < 0)
        {
            const auto reference =
                tetOverlapVolume::cellVolumeMinDecomp(mesh, celli);
            vol[celli] = reference.first();
            roundoff[celli] = reference.second();
            mismatches += mag(vol[celli]/mesh.cellVolumes()[celli]-1)
                > roundoff[celli];
        }
    };

    const auto& tree = tgt_.cellTree();
    List<DynamicList<label>> reverseAddr(tgt_.nCells());
    List<DynamicList<scalar>> reverseWeights(tgt_.nCells());
    List<DynamicList<point>> reverseVectors(centroids ? tgt_.nCells() : 0);
    tetOverlapVolume overlap;
    for (const label srcCelli : maskCells())
    {
        evaluate(src_, srcCelli, srcVol, srcRoundoff);
        // Box queries find contained and disconnected overlaps without
        // depending on vertex seeds, neighbour connectivity, or sliver tests.
        const treeBoundBox sourceBox(src_.cellBb(srcCelli));
        const labelList candidates(tree.findBox(sourceBox));
        DynamicList<label> addresses(candidates.size());
        DynamicList<scalar> weights(candidates.size());
        DynamicList<point> vectors(centroids ? candidates.size() : 0);
        for (const label shapei : candidates)
        {
            const label tgtCelli = tree.shapes().objectIndex(shapei);
            ++candidateTests_;
            evaluate(tgt_, tgtCelli, tgtVol, tgtRoundoff);
            const treeBoundBox targetBox(tgt_.cellBb(tgtCelli));
            scalar volume = 0;
            point centroid(Zero);
            if (centroids)
            {
                const auto moment = overlap.cellCellOverlapMomentMinDecomp
                (
                    src_, srcCelli, tgt_, tgtCelli, targetBox
                );
                volume = moment.first();
                if (volume > 0)
                {
                    centroid = moment.second()/volume;
                }
            }
            else if
            (
                overlapShortcut
             && targetBox.contains(sourceBox.min())
             && targetBox.contains(sourceBox.max())
             && containsCell(tgt_, tgtCelli, src_, srcCelli)
            )
            {
                // Source cell inside convex target cell
                volume = srcVol[srcCelli];
            }
            else if
            (
                overlapShortcut
             && sourceBox.contains(targetBox.min())
             && sourceBox.contains(targetBox.max())
             && containsCell(src_, srcCelli, tgt_, tgtCelli)
            )
            {
                // Target cell inside convex source cell
                volume = tgtVol[tgtCelli];
            }
            else
            {
                volume = overlap.cellCellOverlapVolumeMinDecomp
                (
                    src_, srcCelli, tgt_, tgtCelli, targetBox
                );
            }
            const scalar referenceVolume = min(srcVol[srcCelli], tgtVol[tgtCelli]);
            const scalar allowance =
                max(srcRoundoff[srcCelli], tgtRoundoff[tgtCelli])*referenceVolume;
            // Bound violations slightly beyond the a-priori roundoff
            // estimate (e.g. strongly warped faces) are clamped; only
            // non-finite or gross violations are fatal.
            const scalar grossAllowance =
                max(1e3*allowance, ROOTVSMALL*referenceVolume);
            if
            (
                !std::isfinite(volume)
             || volume < -grossAllowance
             || volume-referenceVolume > grossAllowance
            )
            {
                FatalErrorInFunction
                    << "Invalid intersection volume " << volume
                    << " for cells " << srcCelli << ", " << tgtCelli
                    << ", geometric bound " << referenceVolume
                    << ", roundoff allowance " << allowance
                    << exit(FatalError);
            }
            if (volume < -allowance || volume-referenceVolume > allowance)
            {
                ++clampedVolumes;
                maxViolation = max
                (
                    maxViolation,
                    max(-volume, volume-referenceVolume)/referenceVolume
                );
            }
            if (volume > referenceVolume)
            {
                volume = referenceVolume;
            }
            if (volume < 0)
            {
                ++roundoffVolumes;
                volume = 0;
            }
            if (volume > tolerance_*referenceVolume)
            {
                ++acceptedPairs;
                addresses.append(tgtCelli);
                weights.append(volume);
                reverseAddr[tgtCelli].append(srcCelli);
                reverseWeights[tgtCelli].append(volume);
                if (centroids)
                {
                    vectors.append(centroid - tgt_.cellCentres()[tgtCelli]);
                    reverseVectors[tgtCelli].append
                    (
                        centroid - src_.cellCentres()[srcCelli]
                    );
                }
                V_ += volume;
            }
            else
            {
                discardedVolume_ += volume;
            }
        }
        srcToTgtAddr[srcCelli].transfer(addresses);
        srcToTgtWght[srcCelli].transfer(weights);
        if (centroids)
        {
            srcToTgtVec[srcCelli].transfer(vectors);
        }
    }
    forAll(reverseAddr, celli)
    {
        tgtToSrcAddr[celli].transfer(reverseAddr[celli]);
        tgtToSrcWght[celli].transfer(reverseWeights[celli]);
        if (centroids)
        {
            tgtToSrcVec[celli].transfer(reverseVectors[celli]);
        }
    }
    if (clampedVolumes)
    {
        // Local (per-rank) report: this routine is not called collectively
        // on ranks without cells, so no reduction is possible here.
        Pout<< "Warning: cellVolumeWeight clamped " << clampedVolumes
            << " intersection volumes exceeding the roundoff allowance "
            << "(maximum relative violation " << maxViolation
            << "). Check cell quality (warped faces)." << endl;
    }
    if (debug)
    {
        Pout<< "cellVolumeWeight search: " << candidateTests_
            << " exact candidate tests, " << acceptedPairs
            << " accepted pairs, " << discardedVolume_
            << " discarded overlap volume, " << roundoffVolumes
            << " roundoff-volume contacts, " << V_ << " overlap volume, "
            << timer.elapsedTime() << " s, " << src_.nCells()
            << " source and " << tgt_.nCells() << " target cells" << endl;
        Pout<< "cellVolumeWeight local geometric/FV mismatches: "
            << mismatches << " (evaluated cells only; includes distributed "
            << "geometry copies)" << endl;
    }
}

void Foam::cellVolumeWeightMethod::calculate
(
    labelListList& srcToTgtAddr,
    scalarListList& srcToTgtWght,
    pointListList& srcToTgtVec,
    labelListList& tgtToSrcAddr,
    scalarListList& tgtToSrcWght,
    pointListList& tgtToSrcVec
)
{
    calculateAddressing
    (
        srcToTgtAddr, srcToTgtWght, srcToTgtVec,
        tgtToSrcAddr, tgtToSrcWght, tgtToSrcVec, false
    );
}

// ************************************************************************* //
