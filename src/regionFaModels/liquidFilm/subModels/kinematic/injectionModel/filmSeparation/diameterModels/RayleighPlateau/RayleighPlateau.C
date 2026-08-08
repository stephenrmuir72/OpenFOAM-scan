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

#include "RayleighPlateau.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(RayleighPlateau, 0);
addToRunTimeSelectionTable(diameterModel, RayleighPlateau, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

RayleighPlateau::RayleighPlateau
(
    const liquidFilmBase& film,
    const dictionary& dict
)
:
    diameterModel(type(), film, dict),
    Cligament_(dict.getOrDefault<scalar>("Cligament", 1.89)),
    viscousCorrection_(dict.getOrDefault<bool>("viscousCorrection", false))
{
    if (Cligament_ <= VSMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Ligament-to-drop ratio Cligament must be positive, received "
            << Cligament_ << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar RayleighPlateau::dParcel
(
    const scalar h,
    const scalar rho,
    const scalar mu,
    const scalar sigma,
    const scalar Cligament,
    const bool viscousCorrection
)
{
    const scalar hSafe = Foam::max(h, scalar(0));

    if (hSafe <= VSMALL)
    {
        return 0;
    }

    if (viscousCorrection && rho > VSMALL && sigma > ROOTVSMALL)
    {
        const scalar denom = Foam::sqrt(rho*sigma*hSafe);
        const scalar Oh = mu/Foam::max(denom, ROOTVSMALL);

        return Cligament*hSafe*Foam::pow(1.0 + 3.0*Oh, 1.0/6.0);
    }

    return Cligament*hSafe;
}


tmp<scalarField> RayleighPlateau::diameter() const
{
    const scalarField& h = film().h().primitiveField();
    const scalarField& rho = film().rho().primitiveField();
    const scalarField& mu = film().mu().primitiveField();
    const scalarField& sigma = film().sigma().primitiveField();

    auto td = tmp<scalarField>::New(h.size(), Zero);
    auto& d = td.ref();

    forAll(d, facei)
    {
        d[facei] =
            dParcel
            (
                h[facei],
                rho[facei],
                mu[facei],
                sigma[facei],
                Cligament_,
                viscousCorrection_
            );
    }

    clip(d);
    return td;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
