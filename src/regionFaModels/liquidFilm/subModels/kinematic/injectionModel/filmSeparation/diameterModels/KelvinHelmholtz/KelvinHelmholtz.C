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

#include "KelvinHelmholtz.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(KelvinHelmholtz, 0);
addToRunTimeSelectionTable(diameterModel, KelvinHelmholtz, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

KelvinHelmholtz::KelvinHelmholtz
(
    const liquidFilmBase& film,
    const dictionary& dict
)
:
    diameterModel(type(), film, dict),
    rhop_(dict.get<scalar>("rhop")),
    Cdrop_(dict.getOrDefault<scalar>("Cdrop", 1.88))
{
    if (rhop_ <= VSMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Primary-phase density rhop must be positive, received "
            << rhop_ << exit(FatalIOError);
    }

    if (Cdrop_ <= VSMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Ligament-to-drop ratio Cdrop must be positive, received "
            << Cdrop_ << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar KelvinHelmholtz::dParcel
(
    const scalar h,
    const scalar rho,
    const scalar mu,
    const scalar sigma,
    const scalar magSqrUrel,
    const scalar rhop,
    const scalar Cdrop,
    const scalar hFallback
)
{
    const scalar hSafe = Foam::max(h, scalar(0));
    const scalar denom = rhop*magSqrUrel;

    if (hSafe > VSMALL && sigma > ROOTVSMALL && denom > ROOTVSMALL)
    {
        const scalar dLigament = Foam::sqrt(16.0*sigma*hSafe/denom);

        if (dLigament > ROOTVSMALL && rho > VSMALL)
        {
            const scalar ohDenom = Foam::sqrt(rho*sigma*dLigament);
            const scalar Oh = mu/Foam::max(ohDenom, ROOTVSMALL);

            return Cdrop*dLigament*Foam::pow(1.0 + 3.0*Oh, 1.0/6.0);
        }
    }

    return Foam::max(hFallback, scalar(0));
}


tmp<scalarField> KelvinHelmholtz::diameter() const
{
    const scalarField& h = film().h().primitiveField();
    const scalarField& rho = film().rho().primitiveField();
    const scalarField& mu = film().mu().primitiveField();
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
                h[facei],
                rho[facei],
                mu[facei],
                sigma[facei],
                magSqr(Upi[facei] - Ufi[facei]),
                rhop_,
                Cdrop_,
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
