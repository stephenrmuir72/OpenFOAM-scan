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

#include "advective.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(advective, 0);
addToRunTimeSelectionTable(timeScaleModel, advective, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

advective::advective(const dictionary& dict)
:
    timeScaleModel(typeName, dict),
    Ca_(dict.getOrDefault<scalar>("Ca", 1))
{
    if (!(Ca_ > 0) || !std::isfinite(Ca_))
    {
        FatalIOErrorInFunction(dict)
            << "Ca must be finite and > 0 for the advective time scale, "
            << "received " << Ca_
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void advective::calculate
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
    duration.primitiveFieldRef() = Zero;

    forAll(duration, facei)
    {
        const bool active =
            mag(mass[facei]) > VSMALL
         || mag(momentum[facei]) > VSMALL
         || mag(normalImpulse[facei]) > VSMALL;

        if (!active)
        {
            continue;
        }

        duration[facei] = calculateAt
        (
            facei, mass[facei], momentum[facei], normalImpulse[facei],
            diameter[facei], rho[facei], sigma[facei]
        );
    }

    duration.correctBoundaryConditions();
}

scalar advective::calculateAt
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
    const scalar sourceMass = mag(mass);
    if
    (
        !(sourceMass > VSMALL)
     || !(diameter > VSMALL)
     || !std::isfinite(sourceMass)
     || !std::isfinite(diameter)
    )
    {
        FatalErrorInFunction
            << "Cannot calculate an advective duration for face " << face
            << ": mass=" << mass << ", diameter=" << diameter
            << abort(FatalError);
    }

    const scalar Un = mag(normalImpulse)/sourceMass;
    if (!(Un > VSMALL) || !std::isfinite(Un))
    {
        FatalErrorInFunction
            << "Invalid normal impact velocity " << Un << " on face " << face
            << abort(FatalError);
    }

    const scalar tau = Ca_*diameter/Un;
    if (!(tau > VSMALL) || !std::isfinite(tau))
    {
        FatalErrorInFunction
            << "Invalid advective duration " << tau << " on face " << face
            << abort(FatalError);
    }
    return tau;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
