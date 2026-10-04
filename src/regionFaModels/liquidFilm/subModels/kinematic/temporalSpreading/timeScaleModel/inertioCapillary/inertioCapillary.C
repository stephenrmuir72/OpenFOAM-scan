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

#include "inertioCapillary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(inertioCapillary, 0);
addToRunTimeSelectionTable(timeScaleModel, inertioCapillary, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inertioCapillary::inertioCapillary(const dictionary& dict)
:
    timeScaleModel(typeName, dict),
    Csigma_(dict.getOrDefault<scalar>("Csigma", 1))
{
    if (!(Csigma_ > 0) || !std::isfinite(Csigma_))
    {
        FatalIOErrorInFunction(dict)
            << "Csigma must be finite and > 0 for the inertio-capillary "
            << "time scale, received " << Csigma_
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void inertioCapillary::calculate
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

scalar inertioCapillary::calculateAt
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
    if
    (
        !(mag(mass) > VSMALL)
     || !(diameter > VSMALL)
     || !(rho > VSMALL)
     || !(sigma > VSMALL)
     || !std::isfinite(mass)
     || !std::isfinite(diameter)
     || !std::isfinite(rho)
     || !std::isfinite(sigma)
    )
    {
        FatalErrorInFunction
            << "Cannot calculate an inertio-capillary duration for face "
            << face << ": mass=" << mass << ", diameter=" << diameter
            << ", rho=" << rho << ", sigma=" << sigma
            << abort(FatalError);
    }

    const scalar tau = Csigma_*Foam::sqrt(rho*pow3(diameter)/sigma);
    if (!(tau > VSMALL) || !std::isfinite(tau))
    {
        FatalErrorInFunction
            << "Invalid inertio-capillary duration " << tau
            << " on face " << face << abort(FatalError);
    }
    return tau;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
