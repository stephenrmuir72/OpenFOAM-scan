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
#include "trackingInverseDistanceCellCellStencil.H"

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
bool runMethodChecks(const fvMesh& mesh, const word& method)
{
    dictionary disallowInterpolated;
    disallowInterpolated.add("allowInterpolatedDonors", false);
    donorFallbackHarness<StencilType> noInterpolated(mesh, disallowInterpolated);

    if
    (
        !noInterpolated.run(method)
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

    return 0;
}

// ************************************************************************* //
