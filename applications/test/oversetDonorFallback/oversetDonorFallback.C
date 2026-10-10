/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | www.openfoam.com
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Keysight Technologies
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

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "columnFvMesh.H"
#include "globalIndex.H"
#include "Pstream.H"
#include "inverseDistanceCellCellStencil.H"
#include "leastSquaresCellCellStencil.H"
#include "trackingInverseDistanceCellCellStencil.H"

#include <cmath>
#include <iostream>

using namespace Foam;


template<class StencilType>
class donorFallbackHarness
:
    public StencilType
{
    using donorCandidate = typename StencilType::donorCandidate;

public:

    donorFallbackHarness
    (
        const fvMesh& mesh,
        const dictionary& dict
    )
    :
        StencilType(mesh, dict, false)
    {}


    bool run(const word& method)
    {
        const label nCells = this->mesh_.nCells();
        const globalIndex globalCells(nCells);

        labelList localCell;
        labelList globalCell;
        labelList cellOwner;
        if
        (
            !fixtureCellAddressing
            (
                method,
                globalCells,
                localCell,
                globalCell,
                cellOwner
            )
        )
        {
            return false;
        }

        label failure = 0;
        labelList zoneID(nCells, 0);
        labelList patchTypes(nCells, this->OTHER);
        labelList allCellTypes(nCells, this->CALCULATED);
        labelList allDonor(nCells, -1);
        labelListList allStencil(nCells);
        List<DynamicList<donorCandidate>> allCandidates(nCells);

        const label physicalZones[] = {2, 1, 2, 0};
        forAll(localCell, physicalCell)
        {
            const label celli = localCell[physicalCell];
            if (celli >= 0)
            {
                zoneID[celli] = physicalZones[physicalCell];
            }
        }

        if (localCell[0] >= 0)
        {
            allCellTypes[localCell[0]] = this->HOLE;
        }

        bool rejectedHoleDonor = false;
        if (localCell[1] >= 0)
        {
            const label acceptor = localCell[1];
            const label acceptorZone = zoneID[acceptor];

            // This zone-2 source would win the zone tie-break, but is a hole.
            rejectedHoleDonor =
                !this->donorEligible
                (
                    this->HOLE,
                    this->OTHER,
                    patchTypes[acceptor]
                );

            // Keep the invalid candidate in this synthetic ranked list to
            // verify that post-classification filtering also skips it.
            this->recordDonorCandidate
            (
                acceptor,
                acceptorZone,
                physicalZones[0],
                globalCell[0],
                allCandidates,
                allDonor
            );

            this->recordDonorCandidate
            (
                acceptor,
                acceptorZone,
                physicalZones[3],
                globalCell[3],
                allCandidates,
                allDonor
            );

        }

        if (localCell[2] >= 0)
        {
            const label acceptor = localCell[2];
            const label acceptorZone = zoneID[acceptor];

            this->recordDonorCandidate
            (
                acceptor,
                acceptorZone,
                physicalZones[1],
                globalCell[1],
                allCandidates,
                allDonor
            );
            this->recordDonorCandidate
            (
                acceptor,
                acceptorZone,
                physicalZones[3],
                globalCell[3],
                allCandidates,
                allDonor
            );
        }

        this->orderDonorCandidates
        (
            zoneID,
            allCandidates,
            allStencil,
            allDonor
        );

        // Cell 1 also needs a usable stencil so the front walk classifies it
        // as INTERPOLATED instead of HOLE.
        if (localCell[1] >= 0 && allStencil[localCell[1]].empty())
        {
            failure = 1;
        }

        boolList reactivationCandidates(nCells, false);
        boolList sameZoneFallback(nCells, false);
        labelListList reactivationStencil(nCells);
        scalarField allWeight(nCells, Zero);

        this->classifyAndWalkFront
        (
            globalCells,
            zoneID,
            patchTypes,
            allCandidates,
            1.0,
            1,
            -1,
            allCellTypes,
            allStencil,
            allWeight,
            reactivationCandidates,
            sameZoneFallback,
            reactivationStencil
        );

        if (localCell[1] >= 0 && !rejectedHoleDonor)
        {
            failure = 1;
        }

        if
        (
            localCell[1] >= 0
         &&
            (
                allCellTypes[localCell[1]] != this->INTERPOLATED
             || allStencil[localCell[1]].empty()
             || allStencil[localCell[1]][0] != globalCell[3]
            )
        )
        {
            failure = 1;
        }

        if
        (
            localCell[2] >= 0
         &&
            (
                allStencil[localCell[2]].empty()
             || allStencil[localCell[2]][0] != globalCell[3]
            )
        )
        {
            failure = 1;
        }

        if
        (
            Pstream::nProcs() > 1
         &&
            (
                cellOwner[1] == cellOwner[3]
             || cellOwner[2] == cellOwner[1]
            )
        )
        {
            failure = 1;
        }

        reduce(failure, maxOp<label>());
        if (failure)
        {
            return fail(method, "ranked donor fallback did not select the valid alternative");
        }

        if (Pstream::master())
        {
            Info<< "Passed donor fallback checks: " << method << endl;
        }
        return true;
    }


    bool runExpandedStencilEligibility
    (
        const word& method,
        const label adjacentType
    )
    {
        const label nCells = this->mesh_.nCells();
        const globalIndex globalCells(nCells);

        labelList localCell;
        labelList globalCell;
        labelList cellOwner;
        if
        (
            !fixtureCellAddressing
            (
                method,
                globalCells,
                localCell,
                globalCell,
                cellOwner
            )
        )
        {
            return false;
        }

        this->cellTypes_.setSize(nCells);
        this->cellTypes_ = this->CALCULATED;
        forAll(localCell, physicalCell)
        {
            const label celli = localCell[physicalCell];
            if (celli >= 0)
            {
                if (physicalCell == 0 || physicalCell == 1)
                {
                    this->cellTypes_[celli] = this->INTERPOLATED;
                }
                else if (physicalCell == 3)
                {
                    this->cellTypes_[celli] = adjacentType;
                }
            }
        }

        this->interpolationCells_.setSize(localCell[0] >= 0 ? 1 : 0);
        if (localCell[0] >= 0)
        {
            this->interpolationCells_[0] = localCell[0];
        }

        this->cellStencil_.setSize(nCells);
        forAll(this->cellStencil_, celli)
        {
            this->cellStencil_[celli].clear();
        }

        this->cellInterpolationWeights_.setSize(nCells);
        forAll(this->cellInterpolationWeights_, celli)
        {
            this->cellInterpolationWeights_[celli].clear();
        }

        if (localCell[0] >= 0)
        {
            this->cellStencil_[localCell[0]].setSize(1);
            this->cellStencil_[localCell[0]][0] = globalCell[2];
        }

        List<Map<label>> compactMap;
        this->cellInterpolationMap_.reset
        (
            new mapDistribute(globalCells, this->cellStencil_, compactMap)
        );

        this->createStencil(globalCells, this->allowHoleDonors_);

        labelList compactCellTypes(this->cellTypes_);
        labelList compactPhysicalCells(nCells, -1);
        forAll(localCell, physicalCell)
        {
            if (localCell[physicalCell] >= 0)
            {
                compactPhysicalCells[localCell[physicalCell]] = physicalCell;
            }
        }
        this->cellInterpolationMap().distribute(compactCellTypes);
        this->cellInterpolationMap().distribute(compactPhysicalCells);

        label failure = 0;
        if (localCell[0] >= 0)
        {
            const labelList& slots = this->cellStencil_[localCell[0]];
            const scalarList& weights =
                this->cellInterpolationWeights_[localCell[0]];

            bool foundPrimary = false;
            bool foundInterpolatedNeighbour = false;
            bool foundAdjacentNeighbour = false;
            boolList seenPhysicalCell(4, false);

            if
            (
                slots.empty()
             || slots.size() != weights.size()
             || mag(sum(weights) - 1.0) > 100*SMALL
            )
            {
                failure = 1;
            }

            forAll(slots, sloti)
            {
                const label compactI = slots[sloti];
                if
                (
                    compactI < 0
                 || compactI >= compactPhysicalCells.size()
                )
                {
                    failure = 1;
                    continue;
                }

                const label physicalCell = compactPhysicalCells[compactI];
                const label donorType = compactCellTypes[compactI];

                if (physicalCell < 0 || physicalCell >= seenPhysicalCell.size())
                {
                    failure = 1;
                    continue;
                }
                if (seenPhysicalCell[physicalCell])
                {
                    failure = 1;
                }
                seenPhysicalCell[physicalCell] = true;

                if (physicalCell == 2)
                {
                    foundPrimary = true;
                    if (donorType != this->CALCULATED)
                    {
                        failure = 1;
                    }
                }
                else if (physicalCell == 1)
                {
                    foundInterpolatedNeighbour = true;
                    if (donorType != this->INTERPOLATED)
                    {
                        failure = 1;
                    }
                }
                else if (physicalCell == 3)
                {
                    foundAdjacentNeighbour = true;
                    if (donorType != adjacentType)
                    {
                        failure = 1;
                    }
                }
                else
                {
                    failure = 1;
                }

                if
                (
                    !this->allowInterpolatedDonors_
                 && donorType == this->INTERPOLATED
                )
                {
                    failure = 1;
                }

                if
                (
                    !this->allowHoleDonors_
                 && donorType == this->HOLE
                )
                {
                    failure = 1;
                }
            }

            if
            (
                !foundPrimary
             || foundInterpolatedNeighbour != this->allowInterpolatedDonors_
             || foundAdjacentNeighbour
                !=
                (
                    adjacentType != this->HOLE
                 || this->allowHoleDonors_
                )
            )
            {
                failure = 1;
            }
        }

        if
        (
            Pstream::nProcs() > 1
         &&
            (
                cellOwner[0] == cellOwner[2]
             || cellOwner[1] == cellOwner[2]
            )
        )
        {
            failure = 1;
        }

        reduce(failure, maxOp<label>());
        if (failure)
        {
            return fail
            (
                method,
                "expanded stencil contains an ineligible or missing donor"
            );
        }
        return true;
    }


    bool runCoincidentAcceptorStencil(const word& method)
    {
        const label nCells = this->mesh_.nCells();
        const globalIndex globalCells(nCells);

        labelList localCell;
        labelList globalCell;
        labelList cellOwner;
        if
        (
            !fixtureCellAddressing
            (
                method,
                globalCells,
                localCell,
                globalCell,
                cellOwner
            )
        )
        {
            return false;
        }

        label failure = 0;
        if
        (
            Pstream::nProcs() > 1
         &&
            (
                cellOwner[0] != cellOwner[1]
             || cellOwner[0] == cellOwner[2]
            )
        )
        {
            failure = 1;
        }

        if ((localCell[0] >= 0) != (localCell[1] >= 0))
        {
            failure = 1;
        }

        this->cellTypes_.setSize(nCells);
        this->cellTypes_ = this->HOLE;
        for (label physicalCell = 0; physicalCell < 2; ++physicalCell)
        {
            if (localCell[physicalCell] >= 0)
            {
                this->cellTypes_[localCell[physicalCell]] =
                    this->INTERPOLATED;
            }
        }
        if (localCell[2] >= 0)
        {
            this->cellTypes_[localCell[2]] = this->CALCULATED;
        }

        this->interpolationCells_.setSize
        (
            (localCell[0] >= 0) + (localCell[1] >= 0)
        );
        label interpolationCellI = 0;
        for (label physicalCell = 0; physicalCell < 2; ++physicalCell)
        {
            if (localCell[physicalCell] >= 0)
            {
                this->interpolationCells_[interpolationCellI++] =
                    localCell[physicalCell];
            }
        }

        this->cellStencil_.setSize(nCells);
        forAll(this->cellStencil_, celli)
        {
            this->cellStencil_[celli].clear();
        }

        this->cellInterpolationWeights_.setSize(nCells);
        forAll(this->cellInterpolationWeights_, celli)
        {
            this->cellInterpolationWeights_[celli].clear();
        }

        for (label acceptorI = 0; acceptorI < 2; ++acceptorI)
        {
            const label celli = localCell[acceptorI];
            if (celli >= 0)
            {
                this->cellStencil_[celli].setSize(1);
                this->cellStencil_[celli][0] = globalCell[2];
            }
        }

        List<Map<label>> compactMap;
        this->cellInterpolationMap_.reset
        (
            new mapDistribute(globalCells, this->cellStencil_, compactMap)
        );

        // Emulate coincident centres in overlapping overset regions.
        vectorField& cellCentres =
            const_cast<vectorField&>(this->mesh_.cellCentres());
        point originalSecondAcceptorCentre(Zero);
        if (localCell[0] >= 0 && localCell[1] >= 0)
        {
            originalSecondAcceptorCentre = cellCentres[localCell[1]];
            cellCentres[localCell[1]] = cellCentres[localCell[0]];
        }

        this->createStencil(globalCells, this->allowHoleDonors_);

        if (localCell[0] >= 0 && localCell[1] >= 0)
        {
            cellCentres[localCell[1]] = originalSecondAcceptorCentre;
        }

        labelList compactGlobalCell(nCells);
        forAll(compactGlobalCell, celli)
        {
            compactGlobalCell[celli] = globalCells.toGlobal(celli);
        }
        this->cellInterpolationMap().distribute(compactGlobalCell);

        if (localCell[0] >= 0 && localCell[1] >= 0)
        {
            const label firstAcceptor = localCell[0];
            const label secondAcceptor = localCell[1];
            const labelList& firstStencil =
                this->cellStencil_[firstAcceptor];
            const labelList& secondStencil =
                this->cellStencil_[secondAcceptor];
            const scalarList& firstWeights =
                this->cellInterpolationWeights_[firstAcceptor];
            const scalarList& secondWeights =
                this->cellInterpolationWeights_[secondAcceptor];

            if
            (
                this->cellTypes_[firstAcceptor] != this->INTERPOLATED
             || this->cellTypes_[secondAcceptor] != this->INTERPOLATED
             || firstStencil.empty()
             || firstStencil.size() != secondStencil.size()
             || firstWeights.size() != firstStencil.size()
             || secondWeights.size() != secondStencil.size()
             || mag(sum(firstWeights) - 1.0) > 100*SMALL
             || mag(sum(secondWeights) - 1.0) > 100*SMALL
            )
            {
                failure = 1;
            }
            else
            {
                forAll(firstStencil, sloti)
                {
                    if
                    (
                        firstStencil[sloti] != secondStencil[sloti]
                     || mag(firstWeights[sloti] - secondWeights[sloti])
                        > 100*SMALL
                     || firstStencil[sloti] < 0
                     || firstStencil[sloti] >= compactGlobalCell.size()
                     || compactGlobalCell[firstStencil[sloti]]
                        != globalCell[2]
                    )
                    {
                        failure = 1;
                    }
                }
            }
        }

        reduce(failure, maxOp<label>());
        if (failure)
        {
            return fail
            (
                method,
                "coincident acceptors did not both receive the shared stencil"
            );
        }

        if (Pstream::master())
        {
            Info<< "Passed coincident acceptor stencil checks: "
                << method << endl;
        }
        return true;
    }


    bool runWalkFrontChecks(const word& method)
    {
        const label nCells = this->mesh_.nCells();
        const globalIndex globalCells(nCells);

        labelList localCell;
        labelList globalCell;
        labelList cellOwner;
        if
        (
            !fixtureCellAddressing
            (
                method,
                globalCells,
                localCell,
                globalCell,
                cellOwner
            )
        )
        {
            return false;
        }

        label failure = 0;

        if
        (
            !this->validWalkFrontLayerSettings(2, -1)
         || !this->validWalkFrontLayerSettings(2, 1)
         || !this->validWalkFrontLayerSettings(2, 2)
         || this->validWalkFrontLayerSettings(2, -2)
         || this->validWalkFrontLayerSettings(2, 0)
         || this->validWalkFrontLayerSettings(2, 3)
         || this->validWalkFrontLayerSettings(0, -1)
         || this->validWalkFrontLayerSettings(-1, -1)
        )
        {
            failure = 1;
        }

        // These counts include the current layer, as in walkFront().
        labelList holeCounts(2);
        holeCounts[0] = 1;
        holeCounts[1] = 3;
        if
        (
            !this->exceedsHoleGrowthLimit(2, 3, holeCounts)
         || this->exceedsHoleGrowthLimit(2, 2, holeCounts)
         || this->exceedsHoleGrowthLimit(1, 3, holeCounts)
        )
        {
            failure = 1;
        }

        scalarList emptyScores(2, GREAT);
        scalarList populatedScores(2);
        populatedScores[0] = 0.8;
        populatedScores[1] = 1.0;
        scalarList scoresWithEmptyLayer(2);
        scoresWithEmptyLayer[0] = GREAT;
        scoresWithEmptyLayer[1] = 0.8;
        if
        (
            this->selectLayerByVolumeRatio(emptyScores, false) != 0
         || this->selectLayerByVolumeRatio(populatedScores, true) != 0
         || this->selectLayerByVolumeRatio(scoresWithEmptyLayer, true) != 1
        )
        {
            failure = 1;
        }

        reduce(failure, maxOp<label>());
        if (failure)
        {
            return fail(method, "walk-front layer validation or safeguard checks failed");
        }

        labelList zoneID(nCells, 0);
        labelList allCellTypes(nCells, this->CALCULATED);
        labelListList allStencil(nCells);
        labelListList compactStencil(nCells);
        scalarField allWeight(nCells, Zero);
        scalarList compactCellVol(4, Zero);

        const scalar referenceVolume = this->mesh_.V()[0];
        compactCellVol[0] = referenceVolume/0.8;
        compactCellVol[1] = referenceVolume;
        compactCellVol[2] = referenceVolume;
        compactCellVol[3] = referenceVolume;

        forAll(localCell, physicalCell)
        {
            const label celli = localCell[physicalCell];
            if (celli >= 0)
            {
                if (physicalCell == 1)
                {
                    allCellTypes[celli] = this->HOLE;
                }
                else
                {
                    allStencil[celli].setSize(1);
                    allStencil[celli][0] = 0;
                    compactStencil[celli].setSize(1);
                    compactStencil[celli][0] =
                        physicalCell == 3 ? 1 : 0;
                }
            }
        }

        this->walkFront
        (
            globalCells,
            1.0,
            allStencil,
            allCellTypes,
            allWeight,
            compactCellVol,
            compactStencil,
            zoneID,
            2,
            -1
        );

        failure = 0;
        forAll(localCell, physicalCell)
        {
            const label celli = localCell[physicalCell];
            if (celli >= 0)
            {
                const label expectedType =
                    physicalCell == 1
                  ? this->HOLE
                  : physicalCell == 3
                  ? this->CALCULATED
                  : this->INTERPOLATED;

                if (allCellTypes[celli] != expectedType)
                {
                    failure = 1;
                }
            }
        }
        reduce(failure, maxOp<label>());
        if (failure)
        {
            return fail
            (
                method,
                "automatic layer selection did not choose the lower-ratio layer"
            );
        }

        // Empty stencils on both sides of the initial hole exercise owner and
        // neighbour hole accounting. The run log should report two new holes.
        allCellTypes = this->CALCULATED;
        allWeight = Zero;
        forAll(allStencil, celli)
        {
            allStencil[celli].clear();
            compactStencil[celli].clear();
        }
        forAll(localCell, physicalCell)
        {
            const label celli = localCell[physicalCell];
            if (celli >= 0)
            {
                if (physicalCell == 1)
                {
                    allCellTypes[celli] = this->HOLE;
                }
                else if (physicalCell == 3)
                {
                    allStencil[celli].setSize(1);
                    allStencil[celli][0] = 0;
                    compactStencil[celli].setSize(1);
                    compactStencil[celli][0] = 1;
                }
            }
        }

        this->walkFront
        (
            globalCells,
            1.0,
            allStencil,
            allCellTypes,
            allWeight,
            compactCellVol,
            compactStencil,
            zoneID,
            2,
            1
        );

        failure = 0;
        forAll(localCell, physicalCell)
        {
            const label celli = localCell[physicalCell];
            if (celli >= 0)
            {
                const label expectedType =
                    physicalCell == 3 ? this->INTERPOLATED : this->HOLE;
                if (allCellTypes[celli] != expectedType)
                {
                    failure = 1;
                }
            }
        }
        reduce(failure, maxOp<label>());
        if (failure)
        {
            return fail
            (
                method,
                "empty-stencil front cells were not retained as holes"
            );
        }

        // A missing stencil on the second layer grows the hole count from
        // zero. Explicit layer selection must still be able to use layer 3.
        allCellTypes = this->CALCULATED;
        allWeight = Zero;
        forAll(allStencil, celli)
        {
            allStencil[celli].clear();
            compactStencil[celli].clear();
        }
        forAll(localCell, physicalCell)
        {
            const label celli = localCell[physicalCell];
            if (celli >= 0)
            {
                if (physicalCell == 1)
                {
                    allCellTypes[celli] = this->HOLE;
                }
                else if (physicalCell == 0 || physicalCell == 2)
                {
                    allStencil[celli].setSize(1);
                    allStencil[celli][0] = 0;
                    compactStencil[celli].setSize(1);
                    compactStencil[celli][0] = 0;
                }
            }
        }

        this->walkFront
        (
            globalCells,
            1.0,
            allStencil,
            allCellTypes,
            allWeight,
            compactCellVol,
            compactStencil,
            zoneID,
            3,
            3
        );

        failure = 0;
        forAll(allCellTypes, celli)
        {
            if (allCellTypes[celli] != this->HOLE)
            {
                failure = 1;
            }
        }
        reduce(failure, maxOp<label>());
        if (failure)
        {
            return fail
            (
                method,
                "explicit layer selection was truncated by the leak guard"
            );
        }

        // No fronts means no score samples; automatic selection must safely
        // retain the first, all-hole layer.
        allCellTypes = this->HOLE;
        allWeight = Zero;
        forAll(allStencil, celli)
        {
            allStencil[celli].clear();
            compactStencil[celli].clear();
        }

        this->walkFront
        (
            globalCells,
            1.0,
            allStencil,
            allCellTypes,
            allWeight,
            compactCellVol,
            compactStencil,
            zoneID,
            2,
            -1
        );

        failure = 0;
        forAll(allCellTypes, celli)
        {
            if
            (
                allCellTypes[celli] != this->HOLE
             || mag(allWeight[celli]) > SMALL
            )
            {
                failure = 1;
            }
        }
        reduce(failure, maxOp<label>());
        if (failure)
        {
            return fail(method, "empty-layer automatic fallback failed");
        }

        if (Pstream::master())
        {
            Info<< "Passed walk-front checks: " << method << endl;
        }
        return true;
    }


private:

    bool fixtureCellAddressing
    (
        const word& method,
        const globalIndex& globalCells,
        labelList& localCell,
        labelList& globalCell,
        labelList& cellOwner
    ) const
    {
        List<labelList> procCellData(Pstream::nProcs());
        procCellData[Pstream::myProcNo()].setSize(8, -1);
        localCell.setSize(4, -1);

        forAll(this->mesh_.cellCentres(), celli)
        {
            const label physicalCell =
                label(this->mesh_.cellCentres()[celli].x());

            if (physicalCell >= 0 && physicalCell < 4)
            {
                localCell[physicalCell] = celli;
                procCellData[Pstream::myProcNo()][physicalCell] =
                    globalCells.toGlobal(celli);
                procCellData[Pstream::myProcNo()][4 + physicalCell] =
                    Pstream::myProcNo();
            }
        }

        Pstream::allGatherList(procCellData);

        globalCell.setSize(4);
        globalCell = -1;
        cellOwner.setSize(4);
        cellOwner = -1;
        label failure = 0;

        forAll(procCellData, procI)
        {
            for (label physicalCell = 0; physicalCell < 4; ++physicalCell)
            {
                const label globalI = procCellData[procI][physicalCell];
                if (globalI >= 0)
                {
                    if (globalCell[physicalCell] >= 0)
                    {
                        failure = 1;
                    }
                    else
                    {
                        globalCell[physicalCell] = globalI;
                        cellOwner[physicalCell] =
                            procCellData[procI][4 + physicalCell];
                    }
                }
            }
        }

        forAll(globalCell, physicalCell)
        {
            if (globalCell[physicalCell] < 0)
            {
                failure = 1;
            }
        }

        reduce(failure, maxOp<label>());
        if (failure)
        {
            return fail(method, "could not identify all fixture cells");
        }
        return true;
    }


    bool fail(const word& method, const char* reason) const
    {
        if (Pstream::master())
        {
            std::cerr << "Donor fallback regression failed for " << method
                << ": " << reason << '\n';
        }
        return false;
    }
};


template<class StencilType>
bool runWeightingChecks(const fvMesh& mesh, const word& method)
{
    label failure = 0;

    dictionary defaultDict;
    donorFallbackHarness<StencilType> defaultWeights(mesh, defaultDict);

    pointList idwDonors(3);
    idwDonors[0] = point(-1, 0, 0);
    idwDonors[1] = point(0, 0, 0);
    idwDonors[2] = point(1, 0, 0);
    scalarList weights;
    defaultWeights.stencilWeights(point(0.25, 0, 0), idwDonors, weights);

    scalar idwValue = 0;
    forAll(weights, i)
    {
        idwValue += weights[i]*idwDonors[i].x();
    }
    const scalar idwExpected = 0.0869565;
    const scalar idwTolerance = 1e-6*max(scalar(1), mag(idwExpected));
    if
    (
        weights.size() != idwDonors.size()
     || mag(sum(weights) - 1) > idwTolerance
     || mag(idwValue - idwExpected) > idwTolerance
    )
    {
        failure = 1;
    }

    pointList donors(4);
    donors[0] = point(0, 0, 0);
    donors[1] = point(2, 0, 0);
    donors[2] = point(0, 3, 0);
    donors[3] = point(0, 0, 4);
    const point sample(0.4, 0.6, 0.8);
    const scalar scalarConstant = 2.75;
    scalarList scalarField(4, scalarConstant);
    scalarList affineScalar(4);
    List<vector> affineVector(4);
    forAll(donors, i)
    {
        affineScalar[i] =
            1.5 + 2.0*donors[i].x() - 0.75*donors[i].y()
          + 0.5*donors[i].z();
        affineVector[i] = vector
        (
            donors[i].x() + 2.0*donors[i].y() - donors[i].z(),
            2.0*donors[i].x() - donors[i].y() + 3.0*donors[i].z(),
            -donors[i].x() + 0.5*donors[i].y() + 2.0*donors[i].z()
        );
    }

    const scalar expectedScalar =
        1.5 + 2.0*sample.x() - 0.75*sample.y() + 0.5*sample.z();
    const vector expectedVector
    (
        sample.x() + 2.0*sample.y() - sample.z(),
        2.0*sample.x() - sample.y() + 3.0*sample.z(),
        -sample.x() + 0.5*sample.y() + 2.0*sample.z()
    );

    const word schemes[] = {"leastSquares", "linearExact"};
    for (const word& scheme : schemes)
    {
        dictionary linearDict;
        linearDict.add("weightingScheme", scheme);
        donorFallbackHarness<StencilType> linearWeights(mesh, linearDict);
        linearWeights.stencilWeights(sample, donors, weights);

        scalar actualConstant = 0;
        scalar actualScalar = 0;
        vector actualVector(Zero);
        forAll(weights, i)
        {
            actualConstant += weights[i]*scalarField[i];
            actualScalar += weights[i]*affineScalar[i];
            actualVector += weights[i]*affineVector[i];
        }

        const scalar tolerance =
            1e-6*max(scalar(1), mag(expectedScalar));
        const scalar vectorTolerance =
            1e-6*max(scalar(1), mag(expectedVector));
        if
        (
            weights.size() != donors.size()
         || mag(sum(weights) - 1) > tolerance
         || mag(actualConstant - scalarConstant) > tolerance
         || mag(actualScalar - expectedScalar) > tolerance
         || mag(actualVector - expectedVector) > vectorTolerance
        )
        {
            failure = 1;
        }
    }

    pointList exactDonors(4);
    exactDonors[0] = point(-2, 1, 0);
    exactDonors[1] = point(3, -1, 2);
    exactDonors[2] = point(0.5, -0.25, 1.25);
    exactDonors[3] = point(0, 2, -1);
    dictionary exactDict;
    exactDict.add("weightingScheme", "linearExact");
    donorFallbackHarness<StencilType> exactWeights(mesh, exactDict);
    exactWeights.stencilWeights(exactDonors[2], exactDonors, weights);
    forAll(weights, i)
    {
        if (mag(weights[i] - (i == 2 ? scalar(1) : scalar(0))) > 1e-12)
        {
            failure = 1;
        }
    }

    pointList degenerateDonors(4);
    degenerateDonors[0] = point(-1, 0, 0);
    degenerateDonors[1] = point(0, 0, 0);
    degenerateDonors[2] = point(1, 0, 0);
    degenerateDonors[3] = point(2, 0, 0);
    exactWeights.stencilWeights(point(0.25, 1, 0), degenerateDonors, weights);
    scalar weightSum = 0;
    forAll(weights, i)
    {
        if (!std::isfinite(weights[i]))
        {
            failure = 1;
        }
        weightSum += weights[i];
    }
    if
    (
        weights.size() != degenerateDonors.size()
     || !std::isfinite(weightSum)
     || mag(weightSum - 1) > 1e-6
    )
    {
        failure = 1;
    }

    reduce(failure, maxOp<label>());
    if (failure)
    {
        if (Pstream::master())
        {
            std::cerr << "Weighting regression failed for " << method << '\n';
        }
        return false;
    }

    if (Pstream::master())
    {
        Info<< "Passed weighting scheme checks: " << method << endl;
    }
    return true;
}


bool runLeastSquaresMethodChecks(const fvMesh& mesh)
{
    dictionary dict;
    donorFallbackHarness<cellCellStencils::leastSquares> weightsMethod
    (
        mesh,
        dict
    );

    const pointList donors
    ({
        point(0, 0, 0),
        point(2, 0, 0),
        point(0, 3, 0),
        point(0, 0, 4)
    });
    const point sample(0.4, 0.6, 0.8);
    scalarList weights;
    weightsMethod.stencilWeights(sample, donors, weights);

    scalar actual = 0;
    forAll(weights, i)
    {
        actual += weights[i]
            *(1.5 + 2.0*donors[i].x() - 0.75*donors[i].y()
              + 0.5*donors[i].z());
    }
    const scalar expected =
        1.5 + 2.0*sample.x() - 0.75*sample.y() + 0.5*sample.z();
    label failure =
        weights.size() != donors.size()
     || mag(sum(weights) - 1) > 1e-6
     || mag(actual - expected) > 1e-6*max(scalar(1), mag(expected));

    const pointList exactDonors
    ({
        point(-2, 1, 0),
        point(3, -1, 2),
        point(0.5, -0.25, 1.25),
        point(0, 2, -1)
    });
    weightsMethod.stencilWeights(exactDonors[2], exactDonors, weights);
    forAll(weights, i)
    {
        if (mag(weights[i] - (i == 2 ? scalar(1) : scalar(0))) > 1e-12)
        {
            failure = 1;
        }
    }

    reduce(failure, maxOp<label>());
    if (failure)
    {
        if (Pstream::master())
        {
            std::cerr << "Weighting regression failed for leastSquares\n";
        }
        return false;
    }

    if (Pstream::master())
    {
        Info<< "Passed weighting scheme checks: leastSquares" << endl;
    }
    return true;
}


template<class StencilType>
bool runMethodChecks(const fvMesh& mesh, const word& method)
{
    dictionary disallowInterpolated;
    disallowInterpolated.add("allowInterpolatedDonors", false);
    donorFallbackHarness<StencilType> noInterpolated(mesh, disallowInterpolated);

    if
    (
        !noInterpolated.run(method)
     || !noInterpolated.runCoincidentAcceptorStencil(method)
     || !noInterpolated.runWalkFrontChecks(method)
     || !noInterpolated.runExpandedStencilEligibility
        (
            method,
            cellCellStencil::POROUS
        )
     || !noInterpolated.runExpandedStencilEligibility
        (
            method,
            cellCellStencil::HOLE
        )
    )
    {
        return false;
    }

    dictionary allowInterpolated;
    allowInterpolated.add("allowInterpolatedDonors", true);
    donorFallbackHarness<StencilType> withInterpolated(mesh, allowInterpolated);
    if
    (
        !withInterpolated.runExpandedStencilEligibility
        (
            method,
            cellCellStencil::POROUS
        )
    )
    {
        return false;
    }

    dictionary allowHole;
    allowHole.add("allowInterpolatedDonors", false);
    allowHole.add("allowHoleDonors", true);
    donorFallbackHarness<StencilType> withHole(mesh, allowHole);
    if
    (
        !withHole.runExpandedStencilEligibility
        (
            method,
            cellCellStencil::HOLE
        )
    )
    {
        return false;
    }

    if (Pstream::master())
    {
        Info<< "Passed expanded stencil eligibility checks: " << method << endl;
    }
    return true;
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if
    (
        !runMethodChecks<cellCellStencils::inverseDistance>
        (
            mesh,
            "inverseDistance"
        )
    )
    {
        return 1;
    }

    if
    (
        !runMethodChecks<cellCellStencils::trackingInverseDistance>
        (
            mesh,
            "trackingInverseDistance"
        )
    )
    {
        return 1;
    }

    bool weightingChecksPassed =
        runWeightingChecks<cellCellStencils::inverseDistance>
        (
            mesh,
            "inverseDistance"
        );
    weightingChecksPassed =
        runWeightingChecks<cellCellStencils::trackingInverseDistance>
        (
            mesh,
            "trackingInverseDistance"
        ) && weightingChecksPassed;
    weightingChecksPassed =
        runLeastSquaresMethodChecks(mesh) && weightingChecksPassed;
    if (!weightingChecksPassed)
    {
        return 1;
    }

    return 0;
}

// ************************************************************************* //
