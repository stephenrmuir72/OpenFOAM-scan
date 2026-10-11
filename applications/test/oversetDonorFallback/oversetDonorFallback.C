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
#include "cellVolumeWeightCellCellStencil.H"
#include "patchFaceIntersection.H"
#include "Random.H"

#include <cmath>
#include <iostream>
#include <limits>

using namespace Foam;


template<class StencilType>
class donorFallbackHarness
:
    public StencilType
{
    using donorCandidate = typename StencilType::donorCandidate;

public:

    using StencilType::reportWeightFallbacks;

    donorFallbackHarness
    (
        const fvMesh& mesh,
        const dictionary& dict
    )
    :
        StencilType(mesh, dict, false)
    {}

    bool checkAffineFixture
    (
        const word& method,
        const label fixture,
        const point& sample,
        const pointList& donors,
        const bool requireNegative = false
    ) const
    {
        const label before = sum(this->weightFallbackCounts_);
        scalarList weights;
        this->stencilWeights(sample, donors, weights);
        bool valid =
            weights.size() == donors.size()
         && this->checkWeights(weights) == this->VALID_WEIGHTS
         && sum(this->weightFallbackCounts_) == before;
        long double actualScalar = 0;
        FixedList<long double, 3> actualVector(Zero);
        long double errorScale = 1;
        bool negative = false;
        for (label i = 0; i < min(weights.size(), donors.size()); ++i)
        {
            const vector d = donors[i] - sample;
            actualScalar += weights[i]*(2.75 + 2*d.x() - 0.75*d.y() + d.z());
            const vector value
            (
                1 + d.x() + 2*d.y() - d.z(),
                2 + 2*d.x() - d.y() + 3*d.z(),
                3 - d.x() + 0.5*d.y() + 2*d.z()
            );
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                actualVector[cmpt] += weights[i]*value[cmpt];
                errorScale += static_cast<long double>(mag(weights[i]))
                    *(mag(donors[i][cmpt]) + mag(sample[cmpt]) + mag(d[cmpt]));
            }
            negative |= weights[i] < 0;
        }
        const long double tolerance =
            128*std::numeric_limits<scalar>::epsilon()*errorScale;
        valid = valid && std::abs(actualScalar - 2.75L) <= tolerance;
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            valid = valid && std::abs(actualVector[cmpt] - (cmpt + 1)) <= tolerance;
        }
        valid = valid && (!requireNegative || negative);
        if (!valid && Pstream::master())
        {
            const auto precision = Info().precision();
            Info().precision(std::numeric_limits<scalar>::max_digits10);
            Info<< "Affine fixture failed: " << method << " " << fixture
                << " fallback calls=" << sum(this->weightFallbackCounts_) - before
                << " sample=" << sample << " donors=" << donors << endl;
            Info().precision(precision);
        }
        return valid;
    }

    bool checkOffPlaneFixture
    (
        const point& sample,
        const pointList& donors
    ) const
    {
        const label before = sum(this->weightFallbackCounts_);
        scalarList weights;
        this->stencilWeights(sample, donors, weights);
        label nearest = 0;
        forAll(donors, i)
        {
            if (magSqr(donors[i] - sample) < magSqr(donors[nearest] - sample))
            {
                nearest = i;
            }
        }
        bool valid =
            sum(this->weightFallbackCounts_) == before + 1
         && weights.size() == donors.size();
        forAll(weights, i)
        {
            valid = valid && weights[i] == (i == nearest ? scalar(1) : scalar(0));
        }
        if (!valid && Pstream::master())
        {
            Info<< "Off-plane fixture failed: sample=" << sample
                << " fallback calls=" << sum(this->weightFallbackCounts_) - before
                << " weights=" << weights << endl;
        }
        return valid;
    }

    bool runDiagnosticReportingChecks(const word& method, const label instance)
    {
        const label nCells = this->mesh_.nCells();
        const globalIndex globalCells(nCells);
        labelList localCell, globalCell, cellOwner;
        if (!fixtureCellAddressing
            (method, globalCells, localCell, globalCell, cellOwner))
        {
            return false;
        }
        label failure = this->invalidStencilReported_ || this->weightFallbackReported_;
        for (label update = 0; update < 3; ++update)
        {
            if (Pstream::master())
            {
                Info<< "Diagnostic fixture begin: " << method
                    << " instance=" << instance << " update=" << update << endl;
            }
            this->cellTypes_ = this->CALCULATED;
            this->cellStencil_.setSize(nCells);
            this->cellInterpolationWeights_.setSize(nCells);
            forAll(this->cellStencil_, celli)
            {
                this->cellStencil_[celli].clear();
                this->cellInterpolationWeights_[celli].clear();
            }
            this->cellInterpolationWeight_.primitiveFieldRef() = Zero;
            if (update < 2 && localCell[0] >= 0)
            {
                const label celli = localCell[0];
                this->cellTypes_[celli] = this->INTERPOLATED;
                this->cellStencil_[celli] = labelList({globalCell[2]});
                this->cellInterpolationWeights_[celli] = scalarList({0.5});
                this->cellInterpolationWeight_[celli] = 1;
            }
            List<Map<label>> compactMap;
            this->cellInterpolationMap_.reset
            (
                new mapDistribute(globalCells, this->cellStencil_, compactMap)
            );
            const boolList mask(nCells, false);
            const labelList patches(nCells, this->OTHER);
            this->recoverSameZoneInterpolations
            (
                globalCells, this->zoneID(), mask, mask, mask, mask, patches,
                false, false, this->cellTypes_, this->cellStencil_,
                this->cellInterpolationWeights_,
                this->cellInterpolationWeight_.primitiveFieldRef(),
                this->interpolationCells_, this->cellInterpolationMap_
            );
            if (localCell[0] >= 0)
            {
                failure |= this->cellTypes_[localCell[0]]
                    != (update < 2 ? this->HOLE : this->CALCULATED);
            }
            if (update < 2 && Pstream::master())
            {
                scalarList weights;
                this->nearestWeights
                (
                    point(0.25, 1, 0),
                    pointList({point(-1, 0, 0), point(0, 0, 0)}),
                    weights, 1
                );
            }
            this->reportWeightFallbacks();
            failure |= sum(this->weightFallbackCounts_) != 0;
            failure |= !this->invalidStencilReported_
                || !this->weightFallbackReported_;
            if (Pstream::master())
            {
                Info<< "Diagnostic fixture end" << endl;
            }
        }
        reduce(failure, maxOp<label>());
        return !failure;
    }


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

            // Insert the two tied zone-0 candidates before the higher-ranked
            // zone-2 candidate. Moving the latter must preserve their order.
            this->recordDonorCandidate
            (
                acceptor,
                acceptorZone,
                physicalZones[3],
                globalCell[3],
                allCandidates,
                allDonor
            );

            this->recordDonorCandidate
            (
                acceptor,
                acceptorZone,
                physicalZones[3],
                globalCell[2],
                allCandidates,
                allDonor
            );

            // This hole donor wins the zone tie-break, but must be skipped.
            this->recordDonorCandidate
            (
                acceptor,
                acceptorZone,
                physicalZones[0],
                globalCell[0],
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
        boolList previouslyBlanked(nCells, false);
        boolList protectedHoles(nCells, false);
        labelListList reactivationStencil(nCells);
        scalarField allWeight(nCells, Zero);

        this->classifyAndWalkFront
        (
            globalCells,
            zoneID,
            patchTypes,
            previouslyBlanked,
            allCandidates,
            1.0,
            1,
            -1,
            allCellTypes,
            allStencil,
            allWeight,
            reactivationCandidates,
            sameZoneFallback,
            reactivationStencil,
            protectedHoles
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


    bool runOversetExpandedStencilEligibility(const word& method)
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
        reduce(failure, maxOp<label>());
        if (failure)
        {
            return fail
            (
                method,
                "fringe-policy fixture does not keep the coincident acceptors together"
            );
        }

        for (label scenario = 0; scenario < 2; ++scenario)
        {
            this->cellTypes_.setSize(nCells);
            this->cellTypes_ = this->CALCULATED;
            this->cellStencil_.setSize(nCells);
            this->cellInterpolationWeights_.setSize(nCells);
            forAll(this->cellStencil_, celli)
            {
                this->cellStencil_[celli].clear();
                this->cellInterpolationWeights_[celli].clear();
            }

            labelList patchTypes(nCells, this->OTHER);
            for (label physicalCell = 0; physicalCell < 2; ++physicalCell)
            {
                if (localCell[physicalCell] >= 0)
                {
                    this->cellTypes_[localCell[physicalCell]] =
                        this->INTERPOLATED;
                    patchTypes[localCell[physicalCell]] =
                        physicalCell == 0 ? this->OVERSET : this->OTHER;
                }
            }
            if (localCell[2] >= 0)
            {
                this->cellTypes_[localCell[2]] = this->CALCULATED;
            }
            if (localCell[3] >= 0)
            {
                this->cellTypes_[localCell[3]] =
                    scenario == 0 ? this->INTERPOLATED : this->CALCULATED;
                patchTypes[localCell[3]] =
                    scenario == 0 ? this->OTHER : this->OVERSET;
            }

            if (localCell[0] >= 0)
            {
                this->cellStencil_[localCell[0]].setSize(2);
                this->cellStencil_[localCell[0]][0] = globalCell[2];
                this->cellStencil_[localCell[0]][1] = globalCell[3];
                this->cellInterpolationWeights_[localCell[0]].setSize(2);
                this->cellInterpolationWeights_[localCell[0]][0] = 0.4;
                this->cellInterpolationWeights_[localCell[0]][1] = 0.6;
            }
            if (localCell[1] >= 0)
            {
                this->cellStencil_[localCell[1]].setSize(2);
                this->cellStencil_[localCell[1]][0] = globalCell[2];
                this->cellStencil_[localCell[1]][1] = globalCell[3];
                this->cellInterpolationWeights_[localCell[1]].setSize(2);
                this->cellInterpolationWeights_[localCell[1]][0] = 0.4;
                this->cellInterpolationWeights_[localCell[1]][1] = 0.6;
            }

            if (scenario == 0 && localCell[3] >= 0)
            {
                this->cellStencil_[localCell[3]].setSize(1);
                this->cellStencil_[localCell[3]][0] = globalCell[2];
                this->cellInterpolationWeights_[localCell[3]].setSize(1);
                this->cellInterpolationWeights_[localCell[3]][0] = 1.0;
            }

            List<Map<label>> compactMap;
            this->cellInterpolationMap_.reset
            (
                new mapDistribute(globalCells, this->cellStencil_, compactMap)
            );
            DynamicList<label> interpolationCells;
            if (localCell[0] >= 0)
            {
                interpolationCells.append(localCell[0]);
            }
            if (localCell[1] >= 0)
            {
                interpolationCells.append(localCell[1]);
            }
            if (scenario == 0 && localCell[3] >= 0)
            {
                interpolationCells.append(localCell[3]);
            }
            this->interpolationCells_.transfer(interpolationCells);
            this->cellInterpolationWeight_.primitiveFieldRef() = Zero;

            boolList noCandidates(nCells, false);
            boolList noFallback(nCells, false);
            boolList noOldHoles(nCells, false);
            boolList noProtectedHoles(nCells, false);
            vector originalSecondCentre(Zero);
            vectorField& cellCentres =
                const_cast<vectorField&>(this->mesh_.cellCentres());
            if (localCell[0] >= 0 && localCell[1] >= 0)
            {
                originalSecondCentre = cellCentres[localCell[1]];
                cellCentres[localCell[1]] = cellCentres[localCell[0]];
            }
            this->recoverSameZoneInterpolations
            (
                globalCells,
                this->zoneID(),
                noCandidates,
                noFallback,
                noOldHoles,
                noProtectedHoles,
                patchTypes,
                this->allowHoleDonors_,
                this->allowInterpolatedDonors_,
                this->cellTypes_,
                this->cellStencil_,
                this->cellInterpolationWeights_,
                this->cellInterpolationWeight_.primitiveFieldRef(),
                this->interpolationCells_,
                this->cellInterpolationMap_
            );
            if (localCell[0] >= 0 && localCell[1] >= 0)
            {
                cellCentres[localCell[1]] = originalSecondCentre;
            }

            labelList compactGlobalCell(nCells, -1);
            forAll(compactGlobalCell, celli)
            {
                compactGlobalCell[celli] = globalCells.toGlobal(celli);
            }
            this->cellInterpolationMap().distribute(compactGlobalCell);

            if (localCell[0] >= 0)
            {
                const labelList& slots = this->cellStencil_[localCell[0]];
                bool hasPrimary = false;
                bool hasForbiddenFringe = false;
                forAll(slots, sloti)
                {
                    if
                    (
                        slots[sloti] < 0
                     || slots[sloti] >= compactGlobalCell.size()
                    )
                    {
                        failure = 1;
                        continue;
                    }
                    const label globalDonor = compactGlobalCell[slots[sloti]];
                    hasPrimary |= globalDonor == globalCell[2];
                    hasForbiddenFringe |= globalDonor == globalCell[3];
                }
                if
                (
                    this->cellTypes_[localCell[0]] != this->INTERPOLATED
                 || !hasPrimary
                 || hasForbiddenFringe
                 || this->cellInterpolationWeights_[localCell[0]].size()
                    != slots.size()
                 || mag(sum(this->cellInterpolationWeights_[localCell[0]]) - 1)
                    > 1e-12
                )
                {
                    failure = 1;
                }
            }

            if (localCell[1] >= 0)
            {
                const labelList& slots = this->cellStencil_[localCell[1]];
                bool hasPrimary = false;
                bool hasAllowedNeighbour = false;
                forAll(slots, sloti)
                {
                    if
                    (
                        slots[sloti] < 0
                     || slots[sloti] >= compactGlobalCell.size()
                    )
                    {
                        failure = 1;
                        continue;
                    }
                    const label globalDonor = compactGlobalCell[slots[sloti]];
                    hasPrimary |= globalDonor == globalCell[2];
                    hasAllowedNeighbour |= globalDonor == globalCell[3];
                }
                if
                (
                    this->cellTypes_[localCell[1]] != this->INTERPOLATED
                 || !hasPrimary
                 || !hasAllowedNeighbour
                 || this->cellInterpolationWeights_[localCell[1]].size()
                    != slots.size()
                 || mag(sum(this->cellInterpolationWeights_[localCell[1]]) - 1)
                    > 1e-12
                )
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
                "coincident acceptors did not receive independently filtered fringe stencils"
            );
        }
        if (Pstream::master())
        {
            Info<< "Passed receiver-specific fringe checks: " << method << endl;
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

        // The helper compares the current count to preceding accepted layers.
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
        boolList selectedProtectedHoles(nCells, false);

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
            -1,
            &selectedProtectedHoles
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
                if
                (
                    (physicalCell == 1) != selectedProtectedHoles[celli]
                )
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

        // Layer 1 has no missing-stencil holes and a score of 1. Layer 2
        // grows one unfillable hole, then reaches a cell with score 0.05.
        // Automatic selection must reject that layer and its protected mask.
        allCellTypes = this->CALCULATED;
        allWeight = Zero;
        forAll(allStencil, celli)
        {
            allStencil[celli].clear();
            compactStencil[celli].clear();
        }
        compactCellVol = referenceVolume;
        compactCellVol[0] = referenceVolume/0.05;
        forAll(localCell, physicalCell)
        {
            const label celli = localCell[physicalCell];
            if (celli >= 0)
            {
                if (physicalCell == 0)
                {
                    allCellTypes[celli] = this->HOLE;
                }
                else if (physicalCell == 1 || physicalCell == 3)
                {
                    allStencil[celli].setSize(1);
                    allStencil[celli][0] = globalCell[0];
                    compactStencil[celli].setSize(1);
                    compactStencil[celli][0] =
                        physicalCell == 1 ? 1 : 0;
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
            -1,
            &selectedProtectedHoles
        );

        failure = 0;
        forAll(localCell, physicalCell)
        {
            const label celli = localCell[physicalCell];
            if (celli >= 0)
            {
                const label expectedType =
                    physicalCell == 0
                  ? this->HOLE
                  : physicalCell == 1
                  ? this->INTERPOLATED
                  : this->CALCULATED;
                if
                (
                    allCellTypes[celli] != expectedType
                 || selectedProtectedHoles[celli] != (physicalCell == 0)
                )
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
                "automatic selection retained a leaking low-score layer"
            );
        }

        // Explicit selection is not subject to automatic rejection. The
        // unfillable front cell is a new hole, not a protected expansion.
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
                if (physicalCell == 0)
                {
                    allCellTypes[celli] = this->HOLE;
                }
                else if (physicalCell == 1 || physicalCell == 3)
                {
                    allStencil[celli].setSize(1);
                    allStencil[celli][0] = globalCell[0];
                    compactStencil[celli].setSize(1);
                    compactStencil[celli][0] =
                        physicalCell == 1 ? 1 : 0;
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
            2,
            &selectedProtectedHoles
        );

        failure = 0;
        forAll(localCell, physicalCell)
        {
            const label celli = localCell[physicalCell];
            if (celli >= 0)
            {
                const bool protectedCell =
                    physicalCell == 0 || physicalCell == 1;
                const label expectedType =
                    physicalCell < 3 ? this->HOLE : this->INTERPOLATED;
                if
                (
                    allCellTypes[celli] != expectedType
                 || selectedProtectedHoles[celli] != protectedCell
                )
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
                "explicit leaking layer or protected-hole provenance was lost"
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


    bool runRecoveryReadinessChecks(const word& method)
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

        auto resetState = [&]()
        {
            this->cellTypes_.setSize(nCells);
            this->cellTypes_ = this->CALCULATED;
            this->cellStencil_.setSize(nCells);
            this->cellInterpolationWeights_.setSize(nCells);
            forAll(this->cellStencil_, celli)
            {
                this->cellStencil_[celli].clear();
                this->cellInterpolationWeights_[celli].clear();
            }
            this->interpolationCells_.clear();
            this->cellInterpolationWeight_.primitiveFieldRef() = Zero;
        };

        auto buildMap = [&]()
        {
            List<Map<label>> compactMap;
            this->cellInterpolationMap_.reset
            (
                new mapDistribute(globalCells, this->cellStencil_, compactMap)
            );
        };

        auto recover = [&]
        (
            const boolList& candidates,
            const boolList& sameZoneFallback,
            const boolList& previouslyBlanked,
            const boolList& protectedHoles,
            const labelList& patchTypes
        )
        {
            this->recoverSameZoneInterpolations
            (
                globalCells,
                this->zoneID(),
                candidates,
                sameZoneFallback,
                previouslyBlanked,
                protectedHoles,
                patchTypes,
                this->allowHoleDonors_,
                this->allowInterpolatedDonors_,
                this->cellTypes_,
                this->cellStencil_,
                this->cellInterpolationWeights_,
                this->cellInterpolationWeight_.primitiveFieldRef(),
                this->interpolationCells_,
                this->cellInterpolationMap_
            );
        };

        // Two recovering cells share a ready donor and an old hole that has
        // already been relabelled INTERPOLATED. Its poison value must not
        // affect either scalar or vector interpolation.
        resetState();
        boolList candidates(nCells, false);
        boolList sameZoneFallback(nCells, false);
        boolList previouslyBlanked(nCells, false);
        boolList protectedHoles(nCells, false);
        labelList patchTypes(nCells, this->OTHER);
        DynamicList<label> receivers;
        for (label physicalCell = 0; physicalCell < 2; ++physicalCell)
        {
            const label celli = localCell[physicalCell];
            if (celli >= 0)
            {
                this->cellTypes_[celli] = this->INTERPOLATED;
                this->cellStencil_[celli].setSize(2);
                this->cellStencil_[celli][0] = globalCell[2];
                this->cellStencil_[celli][1] = globalCell[3];
                this->cellInterpolationWeights_[celli].setSize(2);
                this->cellInterpolationWeights_[celli][0] = 0.25;
                this->cellInterpolationWeights_[celli][1] = 0.75;
                this->cellInterpolationWeight_[celli] = 1.0;
                candidates[celli] = true;
                previouslyBlanked[celli] = true;
                receivers.append(celli);
            }
        }
        if (localCell[2] >= 0)
        {
            this->cellTypes_[localCell[2]] = this->CALCULATED;
        }
        if (localCell[3] >= 0)
        {
            this->cellTypes_[localCell[3]] = this->INTERPOLATED;
            previouslyBlanked[localCell[3]] = true;
            this->cellStencil_[localCell[3]] = labelList(1, globalCell[2]);
            this->cellInterpolationWeights_[localCell[3]] = scalarList(1, 1);
        }
        this->interpolationCells_.transfer(receivers);
        buildMap();
        recover
        (
            candidates,
            sameZoneFallback,
            previouslyBlanked,
            protectedHoles,
            patchTypes
        );

        labelList compactGlobalCell(nCells, -1);
        forAll(compactGlobalCell, celli)
        {
            compactGlobalCell[celli] = globalCells.toGlobal(celli);
        }
        this->cellInterpolationMap().distribute(compactGlobalCell);

        label failure = 0;
        for (label physicalCell = 0; physicalCell < 2; ++physicalCell)
        {
            const label celli = localCell[physicalCell];
            if (celli >= 0)
            {
                const labelList& slots = this->cellStencil_[celli];
                const scalarList& weights =
                    this->cellInterpolationWeights_[celli];
                if
                (
                    slots.size() != 1
                 || slots[0] < 0
                 || slots[0] >= compactGlobalCell.size()
                 || weights.size() != 1
                )
                {
                    failure = 1;
                    continue;
                }

                if
                (
                    this->cellTypes_[celli] != this->INTERPOLATED
                 || mag(weights[0] - 1.0) > 1e-12
                 || compactGlobalCell[slots[0]] != globalCell[2]
                 || mag(this->cellInterpolationWeight_[celli] - 1.0)
                    > 1e-12
                )
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
                "recovery retained a previously blanked donor"
            );
        }

        scalarField scalarPsi(nCells, Zero);
        vectorField vectorPsi(nCells, Zero);
        if (localCell[2] >= 0)
        {
            scalarPsi[localCell[2]] = 1.0;
            vectorPsi[localCell[2]] = vector(1, 2, 3);
        }
        if (localCell[3] >= 0)
        {
            scalarPsi[localCell[3]] = 999.0;
            vectorPsi[localCell[3]] = vector(999, -999, 999);
        }
        for (label physicalCell = 0; physicalCell < 2; ++physicalCell)
        {
            if (localCell[physicalCell] >= 0)
            {
                scalarPsi[localCell[physicalCell]] = -50.0;
                vectorPsi[localCell[physicalCell]] = vector(-50, -50, -50);
            }
        }
        this->interpolate(this->mesh_, scalarPsi);
        this->interpolate(this->mesh_, vectorPsi);

        for (label physicalCell = 0; physicalCell < 2; ++physicalCell)
        {
            const label celli = localCell[physicalCell];
            if
            (
                celli >= 0
             &&
                (
                    mag(scalarPsi[celli] - 1.0) > 1e-12
                 || mag(vectorPsi[celli] - vector(1, 2, 3)) > 1e-12
                )
            )
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
                "recovered field values depended on stale donor poison"
            );
        }

        if (localCell[3] >= 0)
        {
            scalarPsi[localCell[3]] = -1e6;
            vectorPsi[localCell[3]] = vector(-1e6, 2e6, -3e6);
        }
        for (label physicalCell = 0; physicalCell < 2; ++physicalCell)
        {
            if (localCell[physicalCell] >= 0)
            {
                scalarPsi[localCell[physicalCell]] = 50;
                vectorPsi[localCell[physicalCell]] = vector(50, 50, 50);
            }
        }
        this->interpolate(this->mesh_, scalarPsi);
        this->interpolate(this->mesh_, vectorPsi);
        for (label physicalCell = 0; physicalCell < 2; ++physicalCell)
        {
            const label celli = localCell[physicalCell];
            if
            (
                celli >= 0
             &&
                (
                    mag(scalarPsi[celli] - 1.0) > 1e-12
                 || mag(vectorPsi[celli] - vector(1, 2, 3)) > 1e-12
                )
            )
            {
                failure = 1;
            }
        }
        reduce(failure, maxOp<label>());
        if (failure)
        {
            return fail(method, "changing the old-hole poison changed recovery");
        }

        // With no ready donor, recovery must remain a hole. An explicitly
        // protected front hole must also stay a hole despite ready neighbours.
        resetState();
        this->cellTypes_ = this->HOLE;
        candidates = false;
        sameZoneFallback = false;
        previouslyBlanked = false;
        protectedHoles = false;
        if (localCell[0] >= 0)
        {
            const label celli = localCell[0];
            this->cellTypes_[celli] = this->INTERPOLATED;
            this->cellStencil_[celli].setSize(1);
            this->cellStencil_[celli][0] = globalCell[3];
            this->cellInterpolationWeights_[celli].setSize(1);
            this->cellInterpolationWeights_[celli][0] = 1.0;
            this->cellInterpolationWeight_[celli] = 1.0;
            candidates[celli] = true;
            previouslyBlanked[celli] = true;
        }
        if (localCell[3] >= 0)
        {
            this->cellTypes_[localCell[3]] = this->INTERPOLATED;
            previouslyBlanked[localCell[3]] = true;
        }
        buildMap();
        this->interpolationCells_.setSize(localCell[0] >= 0 ? 1 : 0);
        if (localCell[0] >= 0)
        {
            this->interpolationCells_[0] = localCell[0];
        }
        recover
        (
            candidates,
            sameZoneFallback,
            previouslyBlanked,
            protectedHoles,
            patchTypes
        );

        if
        (
            localCell[0] >= 0
         &&
            (
                this->cellTypes_[localCell[0]] != this->HOLE
             || !this->cellStencil_[localCell[0]].empty()
             || !this->cellInterpolationWeights_[localCell[0]].empty()
             || mag(this->cellInterpolationWeight_[localCell[0]]) > SMALL
            )
        )
        {
            failure = 1;
        }
        reduce(failure, maxOp<label>());
        if (failure)
        {
            return fail(method, "unready or protected recovery was reactivated");
        }

        resetState();
        candidates = false;
        sameZoneFallback = false;
        previouslyBlanked = false;
        protectedHoles = false;
        if (localCell[0] >= 0)
        {
            const label celli = localCell[0];
            this->cellTypes_[celli] = this->INTERPOLATED;
            this->cellStencil_[celli].setSize(1);
            this->cellStencil_[celli][0] = globalCell[2];
            this->cellInterpolationWeights_[celli].setSize(1);
            this->cellInterpolationWeights_[celli][0] = 1.0;
            this->cellInterpolationWeight_[celli] = 1.0;
            candidates[celli] = true;
            previouslyBlanked[celli] = true;
            protectedHoles[celli] = true;
        }
        if (localCell[2] >= 0)
        {
            this->cellTypes_[localCell[2]] = this->CALCULATED;
        }
        buildMap();
        this->interpolationCells_.setSize(localCell[0] >= 0 ? 1 : 0);
        if (localCell[0] >= 0)
        {
            this->interpolationCells_[0] = localCell[0];
        }
        recover
        (
            candidates,
            sameZoneFallback,
            previouslyBlanked,
            protectedHoles,
            patchTypes
        );

        if
        (
            localCell[0] >= 0
         &&
            (
                this->cellTypes_[localCell[0]] != this->HOLE
             || !this->cellStencil_[localCell[0]].empty()
             || !this->cellInterpolationWeights_[localCell[0]].empty()
             || mag(this->cellInterpolationWeight_[localCell[0]]) > SMALL
            )
        )
        {
            failure = 1;
        }
        reduce(failure, maxOp<label>());
        if (failure)
        {
            return fail(method, "front-selected hole was reactivated");
        }

        if (Pstream::master())
        {
            Info<< "Passed recovery readiness checks: " << method << endl;
        }
        return true;
    }

    bool runReviewValidationChecks(const word& method)
    {
        const label nCells = this->mesh_.nCells();
        const globalIndex globalCells(nCells);
        labelList localCell, globalCell, cellOwner;
        if (!fixtureCellAddressing
            (method, globalCells, localCell, globalCell, cellOwner))
        {
            return false;
        }

        label failure = 0;
        const scalar epsilon = std::numeric_limits<scalar>::epsilon();
        failure |= this->checkWeights(scalarList({1 + 4*epsilon}))
                != this->VALID_WEIGHTS;
        failure |= this->checkWeights(scalarList({1 + 16*epsilon}))
                != this->UNNORMALIZED_WEIGHTS;
        const scalar maxWeight = std::numeric_limits<scalar>::max();
        failure |= this->checkWeights(scalarList({maxWeight, maxWeight}))
                != this->UNNORMALIZED_WEIGHTS;
        for (label scenario = 0; scenario < 7; ++scenario)
        {
            this->cellTypes_ = this->CALCULATED;
            this->cellStencil_.setSize(nCells);
            this->cellInterpolationWeights_.setSize(nCells);
            forAll(this->cellStencil_, celli)
            {
                this->cellStencil_[celli].clear();
                this->cellInterpolationWeights_[celli].clear();
            }
            this->cellInterpolationWeight_.primitiveFieldRef() = Zero;
            boolList candidates(nCells, false), oldHoles(nCells, false);
            boolList fallback(nCells, false), protectedHoles(nCells, false);
            labelList patches(nCells, this->OTHER);
            if (localCell[3] >= 0 && (scenario == 0 || scenario == 6))
            {
                this->cellTypes_[localCell[3]] =
                    scenario == 0 ? this->HOLE : this->INTERPOLATED;
                oldHoles[localCell[3]] = true;
            }
            if (localCell[0] >= 0)
            {
                const label celli = localCell[0];
                this->cellTypes_[celli] = this->INTERPOLATED;
                this->cellInterpolationWeight_[celli] = 1;
                this->cellStencil_[celli] = labelList
                    ({globalCell[2], globalCell[3]});
                const scalar roundOff =
                    std::sqrt(std::numeric_limits<scalar>::epsilon());
                const scalar scale = 1/roundOff;
                this->cellInterpolationWeights_[celli] =
                    scenario == 0 || scenario == 6 ? scalarList({0.5, 0.5})
                  : scenario == 1 ? scalarList({scale, -scale + 1 + 2*roundOff})
                  : scalarList({0.5, 0.6});
                if (scenario == 3)
                {
                    this->cellInterpolationWeights_[celli][0] =
                        std::numeric_limits<scalar>::quiet_NaN();
                }
                if (scenario == 4)
                {
                    this->cellInterpolationWeights_[celli][0] =
                        std::numeric_limits<scalar>::infinity();
                }
                if (scenario == 5)
                {
                    this->cellInterpolationWeights_[celli].setSize(1);
                }
            }
            for (label update = 0; update < 2; ++update)
            {
                List<Map<label>> compactMap;
                this->cellInterpolationMap_.reset
                (
                    new mapDistribute
                        (globalCells, this->cellStencil_, compactMap)
                );
                this->recoverSameZoneInterpolations
                (
                    globalCells, this->zoneID(), candidates, fallback,
                    oldHoles, protectedHoles, patches, this->allowHoleDonors_,
                    this->allowInterpolatedDonors_, this->cellTypes_,
                    this->cellStencil_, this->cellInterpolationWeights_,
                    this->cellInterpolationWeight_.primitiveFieldRef(),
                    this->interpolationCells_, this->cellInterpolationMap_
                );
                if (localCell[0] >= 0)
                {
                    const label celli = localCell[0];
                    if (scenario >= 2 && scenario != 6)
                    {
                        failure |= this->cellTypes_[celli] != this->HOLE
                            || !this->cellStencil_[celli].empty()
                            || this->cellInterpolationWeight_[celli] != 0;
                    }
                    else
                    {
                        const bool excludesOldHole =
                            (scenario == 0 && !this->allowHoleDonors_)
                         || (scenario == 6
                            && (!this->allowHoleDonors_
                             || !this->allowInterpolatedDonors_));
                        failure |= this->cellTypes_[celli] != this->INTERPOLATED
                            || this->cellStencil_[celli].size()
                            != (excludesOldHole ? 1 : 2);
                    }
                }
                if (update == 0)
                {
                    labelList ids(nCells);
                    forAll(ids, celli)
                    {
                        ids[celli] = globalCells.toGlobal(celli);
                    }
                    this->cellInterpolationMap_->distribute(ids);
                    for (labelList& stencil : this->cellStencil_)
                    {
                        for (label& slot : stencil)
                        {
                            slot = ids[slot];
                        }
                    }
                }
            }
        }
        reduce(failure, maxOp<label>());
        if (failure)
        {
            return fail(method, "ordinary old-hole policy or weight validation");
        }
        if (Pstream::master())
        {
            Info<< "Passed review validation checks: " << method << endl;
        }
        return true;
    }

    bool runRecoveryFallbackPolicies(const word& method)
    {
        const label nCells = this->mesh_.nCells();
        const globalIndex globalCells(nCells);
        labelList localCell, globalCell, cellOwner;
        if
        (
            !fixtureCellAddressing
            (
                method, globalCells, localCell, globalCell, cellOwner
            )
        )
        {
            return false;
        }

        label failure = 0;
        for (label scenario = 0; scenario < 5; ++scenario)
        {
            this->cellTypes_ = labelList(nCells, this->CALCULATED);
            this->cellStencil_ = labelListList(nCells);
            this->cellInterpolationWeights_ = scalarListList(nCells);
            this->cellInterpolationWeight_.primitiveFieldRef() = Zero;
            boolList candidates(nCells, false), fallback(nCells, false);
            boolList oldHoles(nCells, false), protectedHoles(nCells, false);
            labelList patches(nCells, this->OTHER);
            const bool oldInterpolated = scenario == 4;
            const bool interpolatedNeighbour = scenario % 2 == 0;
            if (localCell[0] >= 0)
            {
                this->cellTypes_[localCell[0]] = this->HOLE;
                oldHoles[localCell[0]] = true;
            }
            if (localCell[1] >= 0)
            {
                const label receiver = localCell[1];
                this->cellTypes_[receiver] = this->INTERPOLATED;
                oldHoles[receiver] = true;
                candidates[receiver] = fallback[receiver] = !oldInterpolated;
                patches[receiver] =
                    scenario < 2 ? this->OVERSET : this->OTHER;
                if (oldInterpolated)
                {
                    this->cellStencil_[receiver] =
                        labelList(1, globalCell[0]);
                    this->cellInterpolationWeights_[receiver] =
                        scalarList(1, 1);
                }
            }
            if (localCell[2] >= 0)
            {
                const label donor = localCell[2];
                this->cellTypes_[donor] =
                    interpolatedNeighbour && !oldInterpolated
                  ? this->INTERPOLATED : this->CALCULATED;
                patches[donor] =
                    scenario % 2 ? this->OVERSET : this->OTHER;
                if (this->cellTypes_[donor] == this->INTERPOLATED)
                {
                    this->cellStencil_[donor] = labelList(1, globalCell[3]);
                    this->cellInterpolationWeights_[donor] = scalarList(1, 1);
                }
            }
            List<Map<label>> compactMap;
            this->cellInterpolationMap_.reset
            (
                new mapDistribute(globalCells, this->cellStencil_, compactMap)
            );
            this->recoverSameZoneInterpolations
            (
                globalCells, labelList(nCells, 0), candidates, fallback,
                oldHoles, protectedHoles, patches, this->allowHoleDonors_,
                this->allowInterpolatedDonors_, this->cellTypes_,
                this->cellStencil_, this->cellInterpolationWeights_,
                this->cellInterpolationWeight_.primitiveFieldRef(),
                this->interpolationCells_, this->cellInterpolationMap_
            );
            labelList compactGlobals(nCells);
            forAll(compactGlobals, celli)
            {
                compactGlobals[celli] = globalCells.toGlobal(celli);
            }
            this->cellInterpolationMap().distribute(compactGlobals);
            if (localCell[1] >= 0)
            {
                const label receiver = localCell[1];
                const labelList& slots = this->cellStencil_[receiver];
                const bool shouldRecover =
                    oldInterpolated
                 || (scenario >= 2
                  && (!interpolatedNeighbour || this->allowInterpolatedDonors_));
                if (shouldRecover)
                {
                    if
                    (
                        this->cellTypes_[receiver] != this->INTERPOLATED
                     || slots.size() != 1
                     || slots[0] < 0
                     || slots[0] >= compactGlobals.size()
                     || compactGlobals[slots[0]] != globalCell[2]
                     || this->cellInterpolationWeights_[receiver].size() != 1
                     || !std::isfinite(this->cellInterpolationWeights_[receiver][0])
                     || mag(this->cellInterpolationWeights_[receiver][0] - 1)
                        > 1e-12
                     || this->cellInterpolationWeight_[receiver] != 1
                    )
                    {
                        failure = 1;
                    }
                }
                else if
                (
                    this->cellTypes_[receiver] != this->HOLE
                 || !slots.empty()
                 || !this->cellInterpolationWeights_[receiver].empty()
                 || this->cellInterpolationWeight_[receiver] != 0
                )
                {
                    failure = 1;
                }
            }
        }
        reduce(failure, maxOp<label>());
        if (failure)
        {
            return fail(method, "same-zone fallback bypassed receiver policy or missed an old interpolated hole");
        }
        if (Pstream::master())
        {
            Info<< "Passed recovery fallback policy checks: " << method << endl;
        }
        return true;
    }

    bool runReachabilityRetry(const word& method)
    {
        const label nCells = this->mesh_.nCells();
        const globalIndex globalCells(nCells);
        labelList localCell, globalCell, cellOwner;
        if
        (
            !fixtureCellAddressing
            (
                method, globalCells, localCell, globalCell, cellOwner
            )
        )
        {
            return false;
        }
        label failure = 0;
        for (label scenario = 0; scenario < 4; ++scenario)
        {
            this->cellTypes_ = labelList(nCells, this->CALCULATED);
            labelList types(nCells, this->CALCULATED);
            labelList patches(nCells, this->OTHER);
            List<DynamicList<donorCandidate>> choices(nCells);
            boolList previouslyBlanked(nCells, false);
            if (localCell[0] >= 0)
            {
                const label receiver = localCell[0];
                types[receiver] = this->INTERPOLATED;
                patches[receiver] = this->OVERSET;
                if (scenario != 1 && scenario != 2)
                {
                    choices[receiver].append(donorCandidate(1, globalCell[1]));
                }
                choices[receiver].append(donorCandidate(1, globalCell[3]));
                if (scenario == 3)
                {
                    this->cellTypes_[receiver] = this->HOLE;
                    previouslyBlanked[receiver] = true;
                }
            }
            if (scenario == 2 && localCell[1] >= 0)
            {
                choices[localCell[1]].append(donorCandidate(1, globalCell[3]));
                choices[localCell[1]].append(donorCandidate(1, globalCell[0]));
            }
            if (localCell[2] >= 0)
            {
                types[localCell[2]] = this->HOLE;
            }
            labelListList stencil(nCells), saved(nCells);
            scalarField weights(nCells, Zero);
            boolList candidates(nCells, false), fallback(nCells, false);
            boolList protectedHoles(nCells, false);
            this->classifyAndWalkFront
            (
                globalCells, labelList(nCells, 0), patches,
                previouslyBlanked, choices, 1, scenario == 2 ? 2 : 1,
                scenario == 2 ? 2 : -1, types, stencil, weights,
                candidates, fallback, saved, protectedHoles
            );
            if (scenario != 2 && localCell[0] >= 0)
            {
                const label receiver = localCell[0];
                if
                (
                    types[receiver] != this->INTERPOLATED
                 || stencil[receiver].size() != 1
                 || stencil[receiver][0] != globalCell[3]
                 || protectedHoles[receiver]
                 || (scenario == 3 && !candidates[receiver])
                )
                {
                    Pout<< "Retry scenario " << scenario
                        << " receiver type=" << types[receiver]
                        << " stencil=" << stencil[receiver]
                        << " protected=" << protectedHoles[receiver]
                        << " recovering=" << candidates[receiver] << endl;
                    failure = 1;
                }
            }
            if (scenario == 2 && localCell[1] >= 0)
            {
                if
                (
                    types[localCell[1]] != this->HOLE
                 || !protectedHoles[localCell[1]]
                 || stencil[localCell[1]].size() != 1
                 || stencil[localCell[1]][0] != globalCell[3]
                )
                {
                    Pout<< "Selected-layer receiver type="
                        << types[localCell[1]]
                        << " stencil=" << stencil[localCell[1]]
                        << " protected=" << protectedHoles[localCell[1]]
                        << endl;
                    failure = 1;
                }
            }
        }
        reduce(failure, maxOp<label>());
        if (failure)
        {
            return fail(method, "reachability holes prevented donor retry or selected-layer protection was lost");
        }
        if (Pstream::master())
        {
            Info<< "Passed reachability retry checks: " << method << endl;
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


class volumeFallbackHarness : public cellCellStencils::cellVolumeWeight
{
public:
    volumeFallbackHarness(const fvMesh& mesh, const dictionary& dict)
    :
        cellVolumeWeight(mesh, dict, false)
    {}

    bool run()
    {
        const globalIndex globalCells(mesh_.nCells());
        const labelList selected(identity(mesh_.nCells()));
        const labelList zones(mesh_.nCells(), Zero);
        label failure = 0;
        for (label scenario = 0; scenario < 4; ++scenario)
        {
            const labelList types
                (mesh_.nCells(), scenario == 1 ? INTERPOLATED : CALCULATED);
            const labelList patches
                (mesh_.nCells(), scenario == 0 ? OVERSET : OTHER);
            const boolList oldHoles(mesh_.nCells(), scenario == 3);
            const boolList protectedHoles(mesh_.nCells(), false);
            const boolList candidates(mesh_.nCells(), false);
            const boolList mask
            (
                recoveryDonorMask
                (
                    types, oldHoles, protectedHoles, candidates,
                    allowInterpolatedDonors_
                )
            );
            labelListList neighbours;
            pointListList centres;
            globalSameZoneCellCells
            (
                globalCells, mesh_, zones, mask, selected, neighbours,
                centres, types, patches
            );
            const bool expected =
                scenario == 2 || (scenario == 1 && allowInterpolatedDonors_);
            for (const label celli : selected)
            {
                failure |= neighbours[celli].empty() == expected
                    || neighbours[celli].size() != centres[celli].size();
            }
        }
        reduce(failure, maxOp<label>());
        if (failure)
        {
            std::cerr << "Volume-weight fallback donor policy failed\n";
            return false;
        }
        if (Pstream::master())
        {
            Info<< "Passed volume-weight fallback policy checks" << endl;
        }
        return true;
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
    if (weights[1] != 1)
    {
        failure = 1;
    }
    std::swap(degenerateDonors[0], degenerateDonors[1]);
    exactWeights.stencilWeights(point(0.25, 1, 0), degenerateDonors, weights);
    if (weights[0] != 1)
    {
        failure = 1;
    }
    const scalar distant = 2*std::sqrt(std::numeric_limits<scalar>::max());
    const pointList distantDonors
    ({
        point(2*distant, 0, 0), point(distant, 0, 0), point(3*distant, 0, 0)
    });
    exactWeights.stencilWeights(Zero, distantDonors, weights);
    if (weights[1] != 1 || sum(weights) != 1)
    {
        failure = 1;
    }
    exactWeights.reportWeightFallbacks();

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
bool runAffineRoundoffChecks(const fvMesh& mesh, const word& method)
{
    label failure = 0;
    for (const word& scheme : {word("leastSquares"), word("linearExact")})
    {
        dictionary dict;
        dict.add("weightingScheme", scheme);
        donorFallbackHarness<StencilType> harness(mesh, dict);
        const pointList planar
        ({
            point(-1, -1, 0), point(1, -1, 0), point(1, 1, 0),
            point(-1, 1, 0), point(0.3, 0.2, 0)
        });
        const point sample(0.15, -0.1, 0);
        const tensor rotation
        (
            0.6, -0.64, 0.48,
            0.8, 0.48, -0.36,
            0, 0.6, 0.8
        );
        for (label fixture = 0; fixture < 5; ++fixture)
        {
            const vector offset = fixture == 1 ? vector(0, 0, 0.005)
                : fixture == 4 ? vector(1e4, -2e4, 3e4)
                : fixture >= 2 ? vector(12.3, -4.5, 7.8) : vector(Zero);
            const scalar scale = fixture
                ? max
                (
                    scalar(2e-4),
                    1024*std::numeric_limits<scalar>::epsilon()*mag(offset)
                )
                : 1;
            pointList donors(planar);
            const bool rotate = fixture >= 3;
            for (point& donor : donors)
            {
                donor = offset + scale*(rotate ? rotation & donor : donor);
            }
            point receiver
            (
                offset + scale*(rotate ? rotation & sample : sample)
            );
            if (fixture == 1)
            {
                // Independent centre calculations differ by a few coordinate ulps.
                for (label i = 0; i < 3; ++i)
                {
                    receiver.z() = std::nextafter
                    (
                        receiver.z(), std::numeric_limits<scalar>::infinity()
                    );
                }
            }
            failure |= !harness.checkAffineFixture(method, fixture, receiver, donors);
            const vector normal = rotate ? rotation & vector(0, 0, 1) : vector(0, 0, 1);
            failure |= !harness.checkOffPlaneFixture(receiver + 0.1*scale*normal, donors);
        }
        failure |= !harness.checkAffineFixture
        (
            method, 5, point(2, 0.1, 0), planar, true
        );
        const pointList anisotropic
        ({
            point(0.021640926264394362, 0.021511702597657972, 8.6632908161664963e-05),
            point(0.038426219370354189, 0.0057644683332752364, 0.010146476509941304),
            point(0.039741760485973845, 0.0023795135054964132, 0.020783161356641953),
            point(0.043555766077231377, 1.0526642102772572, 0.0002269181873430377),
            point(0.065767225094314194, 1.0549197846147056, 0.010675898748506121),
            point(0.0036175745358908683, 1.0488304179368257, 0.020298837374149361),
            point(0.009708889660957766, 2.0803663283284775, 0.00021658178137174883),
            point(0.067100049928227462, 2.0469296387200786, 0.010789302297544623),
            point(0.09640847035564569, 2.038076369835832, 0.020137150161493249),
            point(1.0542982286033229, 0.057155057489893059, 7.1612606203204677e-05),
            point(1.0382898844107249, 0.070730138914978521, 0.010237522188081567),
            point(1.0471058350869924, 0.097170411903055587, 0.02038586596686303),
            point(1.099946676038553, 1.0340079403101583, 6.462774683505314e-05),
            point(1.0861523737705443, 1.0836760352129204, 0.010048951485245856),
            point(1.0652357127235637, 1.0804109878799952, 0.020930352264845072),
            point(1.0291957726650496, 2.0327279498622395, 0.00030037568349984735),
            point(1.0302359488172954, 2.0237450867124078, 0.010483688120411134),
            point(1.0021320389108961, 2.0081275965097065, 0.020795888464374034),
            point(2.0384307033989413, 0.034385872708802154, 0.00012002078520514389),
            point(2.0768920950129455, 0.089658658850108153, 0.010201523606199161),
            point(2.0767672036265936, 0.0046384233198845694, 0.020027509602223981),
            point(2.03583508888295, 1.0937449475427476, 0.00051050246830427169),
            point(2.0840295072707802, 1.0018278030583221, 0.010958246633192129),
            point(2.0816771735033601, 1.0034933350597648, 0.020355447301400136),
            point(2.0246092495363155, 2.0474331440263276, 0.00061345553573458533),
            point(2.0568626902547398, 2.0805096986525649, 0.010516363409707272),
            point(2.0318393630095009, 2.031516356537169, 0.020262465000714353)
        });
        failure |= !harness.checkAffineFixture
        (
            method, 6,
            point(1.0762600959182054, 1.0696263659380569, 0.0037560394615102377),
            anisotropic
        );
        const pointList nearHit
        ({
            point(0.01035714286, 0.0016428571424999996, 0.0049999999999999992),
            point(0.010500000000000001, 0.0016428571425000003, 0.0050000000000000001),
            point(0.01035714286, 0.0017857142855, 0.0049999999999999992),
            point(0.01035714286, 0.0015000000000000005, 0.0050000000000000001),
            point(0.010214285715, 0.0016428571425000003, 0.0050000000000000001)
        });
        failure |= !harness.checkAffineFixture
        (
            method, 7,
            point(0.010358045023469278, 0.0016417646055956233, 0.0050000000000000001),
            nearHit
        );
        Random random(7);
        for (label fixture = 0; fixture < 1000; ++fixture)
        {
            pointList donors(27);
            label donorI = 0;
            for (label i = 0; i < 3; ++i)
            {
                for (label j = 0; j < 3; ++j)
                {
                    for (label k = 0; k < 3; ++k)
                    {
                        donors[donorI++] = point
                        (
                            i + 0.1*random.sample01<scalar>(),
                            j + 0.1*random.sample01<scalar>(),
                            (k + 0.1*random.sample01<scalar>())/100
                        );
                    }
                }
            }
            const point receiver
            (
                0.25 + random.sample01<scalar>(),
                0.25 + random.sample01<scalar>(),
                (0.25 + random.sample01<scalar>())/100
            );
            failure |= !harness.checkAffineFixture
            (
                method, fixture + 8, receiver, donors
            );
        }
    }
    reduce(failure, maxOp<label>());
    if (!failure && Pstream::master())
    {
        Info<< "Passed affine roundoff checks: " << method << endl;
    }
    return !failure;
}

bool runAffineRoundoffChecks(const fvMesh& mesh)
{
    bool passed = runAffineRoundoffChecks<cellCellStencils::inverseDistance>
        (mesh, "inverseDistance");
    passed = runAffineRoundoffChecks<cellCellStencils::trackingInverseDistance>
        (mesh, "trackingInverseDistance") && passed;
    passed = runAffineRoundoffChecks<cellCellStencils::leastSquares>
        (mesh, "leastSquares") && passed;
    return passed;
}

template<class StencilType>
bool runDiagnosticReportingChecks(const fvMesh& mesh, const word& method)
{
    for (label instance = 0; instance < 2; ++instance)
    {
        dictionary dict;
        donorFallbackHarness<StencilType> harness(mesh, dict);
        if (!harness.runDiagnosticReportingChecks(method, instance))
        {
            return false;
        }
    }
    return true;
}


class trackingGeometryToleranceHarness
:
    public cellCellStencils::trackingInverseDistance
{
public:

    trackingGeometryToleranceHarness
    (
        const fvMesh& mesh,
        const dictionary& dict
    )
    :
        cellCellStencils::trackingInverseDistance(mesh, dict, false)
    {}

    const vector& geometryTolerance() const
    {
        return smallVec_;
    }
};


bool runTrackingGeometryToleranceChecks(fvMesh& mesh)
{
    const List<labelVector> searchBoxDivisions
    (
        1,
        labelVector(4, 2, 2)
    );

    dictionary zeroToleranceDict;
    zeroToleranceDict.add("searchBoxDivisions", searchBoxDivisions);
    zeroToleranceDict.add("tolerance", scalar(0));

    trackingGeometryToleranceHarness zeroTolerance
    (
        mesh,
        zeroToleranceDict
    );
    zeroTolerance.update();

    label failure = mag(zeroTolerance.geometryTolerance()) > SMALL;

    const scalar tolerance = 1e-6;
    dictionary scaledToleranceDict;
    scaledToleranceDict.add("searchBoxDivisions", searchBoxDivisions);
    scaledToleranceDict.add("tolerance", tolerance);

    trackingGeometryToleranceHarness scaledTolerance
    (
        mesh,
        scaledToleranceDict
    );
    scaledTolerance.update();

    const vector initialExpected(mesh.bounds().span()*tolerance);
    if
    (
        mag(scaledTolerance.geometryTolerance() - initialExpected)
     > 1e-12*max(scalar(1), mag(initialExpected))
    )
    {
        failure = 1;
    }

    pointField movedPoints(mesh.points());
    for (point& p : movedPoints)
    {
        p *= 2.0;
    }
    mesh.movePoints(movedPoints);
    scaledTolerance.update();

    const vector movedExpected(mesh.bounds().span()*tolerance);
    if
    (
        mag(scaledTolerance.geometryTolerance() - movedExpected)
     > 1e-12*max(scalar(1), mag(movedExpected))
     || mag(movedExpected - initialExpected) <= SMALL
    )
    {
        failure = 1;
    }

    reduce(failure, maxOp<label>());
    if (failure)
    {
        if (Pstream::master())
        {
            std::cerr << "Tracking geometry tolerance regression failed\n";
        }
        return false;
    }

    if (Pstream::master())
    {
        Info<< "Passed tracking geometry tolerance checks" << endl;
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
     || !noInterpolated.runReviewValidationChecks(method)
     || !noInterpolated.runCoincidentAcceptorStencil(method)
     || !noInterpolated.runRecoveryReadinessChecks(method)
     || !noInterpolated.runRecoveryFallbackPolicies(method)
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
        !withInterpolated.runReviewValidationChecks(method)
     || !withInterpolated.runRecoveryReadinessChecks(method)
     || !withInterpolated.runRecoveryFallbackPolicies(method)
     || !withInterpolated.runOversetExpandedStencilEligibility(method)
     ||
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
        !withHole.runReviewValidationChecks(method)
     || !withHole.runRecoveryReadinessChecks(method)
     || !withHole.runExpandedStencilEligibility
        (
            method,
            cellCellStencil::HOLE
        )
    )
    {
        return false;
    }

    dictionary holeAndInterpolated;
    holeAndInterpolated.add("allowInterpolatedDonors", true);
    holeAndInterpolated.add("allowHoleDonors", true);
    donorFallbackHarness<StencilType> allDonors(mesh, holeAndInterpolated);
    if
    (
        !allDonors.runReviewValidationChecks(method)
     || !allDonors.runRecoveryReadinessChecks(method)
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
    argList::addBoolOption("retry", "Test donor-graph reachability retries");
    argList::addBoolOption("affineValidation", "Test affine rounding and rank controls");
    argList::addBoolOption("diagnosticValidation", "Test repeated global reporting");
    argList::addOption("voxelBudget", "count", "Test configured patch voxel budget");
    argList::addOption("weightingScheme", "name", "Test weighting dictionary validation");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (args.found("affineValidation"))
    {
        return runAffineRoundoffChecks(mesh) ? 0 : 1;
    }
    if (args.found("diagnosticValidation"))
    {
        return
            runDiagnosticReportingChecks<cellCellStencils::inverseDistance>
                (mesh, "inverseDistance")
         && runDiagnosticReportingChecks<cellCellStencils::trackingInverseDistance>
                (mesh, "trackingInverseDistance")
         && runDiagnosticReportingChecks<cellCellStencils::leastSquares>
                (mesh, "leastSquares") ? 0 : 1;
    }
    if (args.found("weightingScheme"))
    {
        dictionary dict;
        dict.add("weightingScheme", args.get<word>("weightingScheme"));
        donorFallbackHarness<cellCellStencils::inverseDistance> stencil(mesh, dict);
        return 0;
    }
    if (args.found("voxelBudget"))
    {
        const label budget = args.get<label>("voxelBudget");
        dictionary dict;
        dict.add("maxPatchVoxelCount", budget);
        dict.add
        (
            "searchBoxDivisions",
            List<labelVector>(gMax(cellCellStencil::zoneID(mesh)) + 1,
                labelVector(2, 2, 1))
        );
        donorFallbackHarness<cellCellStencils::trackingInverseDistance> stencil(mesh, dict);
        stencil.update();
        label count = 0;
        if (budget >= 27000000
         && (!patchFaceIntersection::checkedVoxelCount
                (labelVector(300, 300, 300), budget, count)
          || count != 27000000))
        {
            return 1;
        }
        Info<< "Passed configured voxel budget checks" << endl;
        return 0;
    }
    if (args.found("retry"))
    {
        dictionary dict;
        dict.add("allowInterpolatedDonors", true);
        donorFallbackHarness<cellCellStencils::inverseDistance> inverse(mesh, dict);
        donorFallbackHarness<cellCellStencils::trackingInverseDistance> tracking(mesh, dict);
        return
            inverse.runReachabilityRetry("inverseDistance")
         && tracking.runReachabilityRetry("trackingInverseDistance") ? 0 : 1;
    }

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

    if (!runAffineRoundoffChecks(mesh))
    {
        return 1;
    }

    if (!runTrackingGeometryToleranceChecks(mesh))
    {
        return 1;
    }
    for (const bool allowInterpolated : {false, true})
    {
        dictionary dict;
        dict.add("allowInterpolatedDonors", allowInterpolated);
        volumeFallbackHarness volume(mesh, dict);
        if (!volume.run())
        {
            return 1;
        }
    }

    return 0;
}

// ************************************************************************* //
