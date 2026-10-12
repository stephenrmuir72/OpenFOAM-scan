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
    scalarField srcVol(src_.nCells()), tgtVol(tgt_.nCells());
    scalarField srcRoundoff(src_.nCells()), tgtRoundoff(tgt_.nCells());
    label mismatches = 0;
    forAll(srcVol, celli)
    {
        const auto reference = tetOverlapVolume::cellVolumeMinDecomp(src_, celli);
        srcVol[celli] = reference.first();
        srcRoundoff[celli] = reference.second();
        mismatches += mag(srcVol[celli]/src_.cellVolumes()[celli]-1)
            > srcRoundoff[celli];
    }
    forAll(tgtVol, celli)
    {
        const auto reference = tetOverlapVolume::cellVolumeMinDecomp(tgt_, celli);
        tgtVol[celli] = reference.first();
        tgtRoundoff[celli] = reference.second();
        mismatches += mag(tgtVol[celli]/tgt_.cellVolumes()[celli]-1)
            > tgtRoundoff[celli];
    }

    const auto& tree = tgt_.cellTree();
    List<DynamicList<label>> reverseAddr(tgt_.nCells());
    List<DynamicList<scalar>> reverseWeights(tgt_.nCells());
    List<DynamicList<point>> reverseVectors(centroids ? tgt_.nCells() : 0);
    tetOverlapVolume overlap;
    for (const label srcCelli : maskCells())
    {
        // Box queries find contained and disconnected overlaps without
        // depending on vertex seeds, neighbour connectivity, or sliver tests.
        const labelList candidates
        (
            tree.findBox(treeBoundBox(src_.cellBb(srcCelli)))
        );
        DynamicList<label> addresses(candidates.size());
        DynamicList<scalar> weights(candidates.size());
        DynamicList<point> vectors(centroids ? candidates.size() : 0);
        for (const label shapei : candidates)
        {
            const label tgtCelli = tree.shapes().objectIndex(shapei);
            ++candidateTests_;
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
            if
            (
                !std::isfinite(volume)
             || volume < -allowance
             || volume-referenceVolume > allowance
            )
            {
                FatalErrorInFunction
                    << "Invalid intersection volume " << volume
                    << " for cells " << srcCelli << ", " << tgtCelli
                    << ", geometric bound " << referenceVolume
                    << ", roundoff allowance " << allowance
                    << exit(FatalError);
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
            << mismatches << " (includes distributed geometry copies)" << endl;
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
