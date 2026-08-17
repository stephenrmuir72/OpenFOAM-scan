/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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

#include "perFaceAggregation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(perFaceAggregation, 0);
addToRunTimeSelectionTable
(
    separationAggregationModel,
    perFaceAggregation,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

perFaceAggregation::perFaceAggregation
(
    const liquidFilmBase& film,
    const dictionary& dict
)
:
    separationAggregationModel(type(), film, dict)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void perFaceAggregation::aggregateFaces
(
    scalarField& massAccumulator,
    vectorField& momentumAccumulator,
    scalarField& diameterAccumulator,
    scalarField& accumulationTime,
    const scalarField& eventTime,
    scalarField& massToInject,
    scalarField& diameterToInject,
    vectorField& velocityToInject
)
{
    const label nFaces = massAccumulator.size();

    if
    (
        momentumAccumulator.size() != nFaces
     || diameterAccumulator.size() != nFaces
     || accumulationTime.size() != nFaces
     || eventTime.size() != nFaces
     || massToInject.size() != nFaces
     || diameterToInject.size() != nFaces
     || velocityToInject.size() != nFaces
    )
    {
        FatalErrorInFunction
            << "All per-face aggregation fields must have the same size"
            << abort(FatalError);
    }

    massToInject = Zero;
    diameterToInject = Zero;
    velocityToInject = Zero;

    forAll(massAccumulator, facei)
    {
        const scalar mass = massAccumulator[facei];

        if (mass > VSMALL && eventTime[facei] <= VSMALL)
        {
            FatalErrorInFunction
                << "Active per-face aggregation face " << facei
                << " has a non-positive physical event time "
                << eventTime[facei] << abort(FatalError);
        }

        if
        (
            mass > VSMALL
         && accumulationTime[facei] >= eventTime[facei]
        )
        {
            massToInject[facei] = mass;
            diameterToInject[facei] = diameterAccumulator[facei]/mass;
            velocityToInject[facei] = momentumAccumulator[facei]/mass;

            massAccumulator[facei] = 0;
            momentumAccumulator[facei] = Zero;
            diameterAccumulator[facei] = 0;
            accumulationTime[facei] = 0;
        }
    }
}


void perFaceAggregation::aggregate
(
    scalarField& massAccumulator,
    vectorField& momentumAccumulator,
    scalarField& diameterAccumulator,
    scalarField& accumulationTime,
    const scalarField& eventTime,
    const scalarField& coherenceLength,
    const labelField& morphology,
    scalarField& massToInject,
    scalarField& diameterToInject,
    vectorField& velocityToInject
) const
{
    (void)coherenceLength;
    (void)morphology;

    aggregateFaces
    (
        massAccumulator,
        momentumAccumulator,
        diameterAccumulator,
        accumulationTime,
        eventTime,
        massToInject,
        diameterToInject,
        velocityToInject
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
