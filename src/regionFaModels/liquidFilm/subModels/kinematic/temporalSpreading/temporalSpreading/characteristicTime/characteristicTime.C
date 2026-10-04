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

#include "characteristicTime.H"
#include "addToRunTimeSelectionTable.H"
#include "timeScaleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(characteristicTime, 0);
addToRunTimeSelectionTable
(
    temporalSpreading,
    characteristicTime,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

characteristicTime::characteristicTime
(
    const faMesh& mesh,
    const dictionary& dict
)
:
    temporalSpreading(typeName, mesh, dict),
    timeScalePtr_(timeScaleModel::New(dict))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

characteristicTime::~characteristicTime()
{}  // timeScaleModel was forward declared


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void characteristicTime::calculateDuration
(
    areaScalarField& duration,
    const areaScalarField& mass,
    const areaVectorField& momentum,
    const areaScalarField& normalImpulse,
    const areaScalarField& diameter,
    const areaScalarField& rho,
    const areaScalarField& sigma
) const
{
    timeScalePtr_->calculate
    (
        duration,
        mass,
        momentum,
        normalImpulse,
        diameter,
        rho,
        sigma
    );
}

scalar characteristicTime::calculateDurationAt
(
    const label face,
    const scalar mass,
    const vector& momentum,
    const scalar normalImpulse,
    const scalar diameter,
    const scalar rho,
    const scalar sigma
) const
{
    return timeScalePtr_->calculateAt
    (
        face, mass, momentum, normalImpulse, diameter, rho, sigma
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
