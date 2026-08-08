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

#include "aerodynamicDiameter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(aerodynamicDiameter, 0);
addToRunTimeSelectionTable(diameterModel, aerodynamicDiameter, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

aerodynamicDiameter::aerodynamicDiameter
(
    const liquidFilmBase& film,
    const dictionary& dict
)
:
    diameterModel(type(), film, dict),
    rhop_(dict.get<scalar>("rhop")),
    WeCrit_(dict.getOrDefault<scalar>("WeCrit", 12))
{
    if (rhop_ <= VSMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Primary-phase density rhop must be positive, received "
            << rhop_ << exit(FatalIOError);
    }

    if (WeCrit_ <= VSMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Critical Weber number WeCrit must be positive, received "
            << WeCrit_ << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar aerodynamicDiameter::dParcel
(
    const scalar sigma,
    const scalar magSqrUrel,
    const scalar rhop,
    const scalar WeCrit,
    const scalar hFallback
)
{
    const scalar denom = rhop*magSqrUrel;

    if (sigma > ROOTVSMALL && denom > ROOTVSMALL)
    {
        return WeCrit*sigma/denom;
    }

    return Foam::max(hFallback, scalar(0));
}


tmp<scalarField> aerodynamicDiameter::diameter() const
{
    const scalarField& h = film().h().primitiveField();
    const scalarField& sigma = film().sigma().primitiveField();

    const areaVectorField Up(film().Up());
    const vectorField& Upi = Up.primitiveField();
    const vectorField& Ufi = film().Uf().primitiveField();

    auto td = tmp<scalarField>::New(h.size(), Zero);
    auto& d = td.ref();

    forAll(d, facei)
    {
        d[facei] = 
            dParcel
            (
                sigma[facei],
                magSqr(Upi[facei] - Ufi[facei]),
                rhop_,
                WeCrit_,
                h[facei]
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
