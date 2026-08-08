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

#include "regimeMapMorphology.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(regimeMapMorphology, 0);
addToRunTimeSelectionTable(morphologyModel, regimeMapMorphology, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

regimeMapMorphology::regimeMapMorphology
(
    const liquidFilmBase& film,
    const dictionary& dict
)
:
    morphologyModel(type(), film, dict),
    rhop_(0),
    WeDripToLigament_(0),
    WeLigamentToSheet_(0),
    lowWeMorphology_(DRIP)
{
    const dictionary& coeffs = dict.subDict("regimeMapCoeffs");

    rhop_ = coeffs.get<scalar>("rhop");
    WeDripToLigament_ = coeffs.get<scalar>("WeDripToLigament");
    WeLigamentToSheet_ = coeffs.get<scalar>("WeLigamentToSheet");
    lowWeMorphology_ =
        morphologyTypeNames.getOrDefault("lowWeMorphology", coeffs, DRIP);

    if (!(rhop_ > VSMALL) || !std::isfinite(rhop_))
    {
        FatalIOErrorInFunction(coeffs)
            << "rhop must be finite and positive, received " << rhop_
            << exit(FatalIOError);
    }

    if
    (
        !(WeDripToLigament_ >= 0)
     || !(WeLigamentToSheet_ > WeDripToLigament_)
     || !std::isfinite(WeDripToLigament_)
     || !std::isfinite(WeLigamentToSheet_)
    )
    {
        FatalIOErrorInFunction(coeffs)
            << "Require finite thresholds 0 <= WeDripToLigament < "
            << "WeLigamentToSheet. Received "
            << WeDripToLigament_ << " and " << WeLigamentToSheet_
            << exit(FatalIOError);
    }

    if
    (
        lowWeMorphology_ != DIRECT_DROPLET
     && lowWeMorphology_ != DRIP
    )
    {
        FatalIOErrorInFunction(coeffs)
            << "lowWeMorphology must be directDroplet or drip, received "
            << morphologyTypeNames[lowWeMorphology_]
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<labelField> regimeMapMorphology::morphology() const
{
    const scalarField& h = film_.h().primitiveField();
    const scalarField& sigma = film_.sigma().primitiveField();

    const areaVectorField Up(film_.Up());
    const vectorField& Upi = Up.primitiveField();
    const vectorField& Ufi = film_.Uf().primitiveField();

    auto tTypes = tmp<labelField>::New(h.size(), label(lowWeMorphology_));
    auto& types = tTypes.ref();

    forAll(types, facei)
    {
        if (h[facei] <= VSMALL || sigma[facei] <= ROOTVSMALL)
        {
            continue;
        }

        const scalar WeGas =
            rhop_*magSqr(Upi[facei] - Ufi[facei])*h[facei]/sigma[facei];

        if (WeGas >= WeLigamentToSheet_)
        {
            types[facei] = SHEET;
        }
        else if (WeGas >= WeDripToLigament_)
        {
            types[facei] = LIGAMENT;
        }
    }

    return tTypes;
}


boolList regimeMapMorphology::possibleTypes() const
{
    boolList result(nMorphologyTypes, false);
    result.set(label(lowWeMorphology_));
    result.set(label(LIGAMENT));
    result.set(label(SHEET));
    return result;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
