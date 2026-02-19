/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Keysight Ltd.
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

#include "GEKO.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void GEKO<BasicTurbulenceModel>::correctNut()
{

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
GEKO<BasicTurbulenceModel>::GEKO
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Cmu_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    Comega2_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Comega2",
            this->coeffDict_,
            0.083
        )
    ),
    sigmakTilde_
    (
        dimensionedScalar::getOrAddToDict
        (
            "sigmakTilde",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaOmegaTilde_
    (
        dimensionedScalar::getOrAddToDict
        (
            "sigmaOmegaTilde",
            this->coeffDict_,
            1.17
        )
    ),
    CRealize_
    (
        dimensionedScalar::getOrAddToDict
        (
            "CRealize",
            this->coeffDict_,
            scalar(1.0)/std::sqrt(scalar(3.0))  // (MM:Eq. 17)
        )
    ),
    kappa_
    (
        dimensionedScalar::getOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    CFbTurb_
    (
        dimensionedScalar::getOrAddToDict
        (
            "CFbTurb",
            this->coeffDict_,
            2.0
        )
    ),
    CFbLam_
    (
        dimensionedScalar::getOrAddToDict
        (
            "CFbLam",
            this->coeffDict_,
            1.0  // increase to 25 in case of transition simulations (MM:p. 4)
        )
    ),
    Aplus_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Aplus",
            this->coeffDict_,
            15.0
        )
    ),
    Cc1_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cc1",
            this->coeffDict_,
            1.7
        )
    ),
    Cc2_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cc2",
            this->coeffDict_,
            1.4
        )
    ),
    CPKlim_
    (
        dimensionedScalar::getOrAddToDict
        (
            "CPKlim",
            this->coeffDict_,
            10.0
        )
    ),

{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool GEKO<BasicTurbulenceModel>::read()
{
    if (!eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        return false;
    }

    Cmu_.readIfPresent(this->coeffDict());

    return true;
}


template<class BasicTurbulenceModel>
void GEKO<BasicTurbulenceModel>::correct()
{

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
