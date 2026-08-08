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

        const scalar sourceMass = mag(mass[facei]);
        const scalar d = diameter[facei];

        if (!(sourceMass > VSMALL) || !(d > VSMALL))
        {
            FatalErrorInFunction
                << "Cannot calculate an advective temporal-spreading duration "
                << "for active face " << facei
                << ": mass=" << mass[facei] << ", diameter=" << d << "."
                << " Active characteristic-time sources require positive "
                << "mass and diameter."
                << abort(FatalError);
        }

        const scalar Un = mag(normalImpulse[facei])/sourceMass;

        if (!(Un > VSMALL) || !std::isfinite(Un))
        {
            FatalErrorInFunction
                << "Cannot calculate an advective temporal-spreading duration "
                << "for active face " << facei << ": |Un|=" << Un
                << ". The advective time scale requires a finite positive "
                << "normal impact velocity."
                << abort(FatalError);
        }

        const scalar tau = Ca_*d/Un;

        if (!(tau > VSMALL) || !std::isfinite(tau))
        {
            FatalErrorInFunction
                << "Invalid advective temporal-spreading duration " << tau
                << " on active face " << facei
                << abort(FatalError);
        }

        duration[facei] = tau;
    }

    duration.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
