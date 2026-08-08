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

#include "diameterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(diameterModel, 0);
defineRunTimeSelectionTable(diameterModel, dictionary);


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void diameterModel::clip
(
    scalarField& d,
    const scalar dMin,
    const scalar dMax
)
{
    for (auto& di : d)
    {
        di = Foam::max(dMin, Foam::min(dMax, di));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

diameterModel::diameterModel
(
    const word&,
    const liquidFilmBase& film,
    const dictionary& dict
)
:
    film_(film),
    mesh_(film.regionMesh()),
    dMin_(dict.getOrDefault<scalar>("dMin", 0)),
    dMax_(dict.getOrDefault<scalar>("dMax", GREAT))
{
    if (dMin_ < 0 || dMax_ <= dMin_)
    {
        FatalIOErrorInFunction(dict)
            << "Invalid diameter clip bounds: dMin = " << dMin_
            << ", dMax = " << dMax_ << nl
            << "Require 0 <= dMin < dMax."
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
