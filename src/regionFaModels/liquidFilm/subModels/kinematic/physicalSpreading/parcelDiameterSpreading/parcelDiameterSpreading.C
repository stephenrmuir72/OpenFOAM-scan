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

#include "parcelDiameterSpreading.H"
#include "Function1.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(parcelDiameterSpreading, 0);
addToRunTimeSelectionTable
(
    physicalSpreading,
    parcelDiameterSpreading,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

parcelDiameterSpreading::parcelDiameterSpreading
(
    const faMesh& mesh,
    const dictionary& dict
)
:
    physicalSpreading(type(), mesh, dict),
    betaPtr_
    (
        Function1<scalar>::NewIfPresent
        (
            "beta",
            dict,
            word::null,
            &mesh.thisDb()
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

parcelDiameterSpreading::~parcelDiameterSpreading()
{}  // Function1 was forward declared


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void parcelDiameterSpreading::spread
(
    areaScalarField& mass,
    areaVectorField& momentum,
    areaScalarField& normalImpulse,
    const areaScalarField& diameter
) const
{
    const scalar beta =
        betaPtr_
      ? betaPtr_->value(mesh_.time().value())
      : scalar(1);

    if (!(beta > 0))
    {
        FatalIOErrorInFunction(dict_)
            << "beta must evaluate to a value > 0 at simulation time "
            << mesh_.time().value() << ", received " << beta
            << exit(FatalIOError);
    }

    spreadSources(mass, momentum, normalImpulse, diameter, beta);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
