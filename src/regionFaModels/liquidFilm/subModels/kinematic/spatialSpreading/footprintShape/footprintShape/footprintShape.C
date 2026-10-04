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

#include "footprintShape.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(footprintShape, 0);
defineRunTimeSelectionTable(footprintShape, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

footprintShape::footprintShape
(
    const word& modelType,
    const dictionary& dict
)
:
    dict_(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void footprintShape::planarOffset
(
    const footprintGeometry& fp,
    const point& target,
    scalar& x,
    scalar& y
) const
{
    const vector delta(target - fp.centre);
    x = delta & fp.majorAxis;
    y = delta & fp.minorAxis;
}


scalar footprintShape::boundaryRadius
(
    const scalar x,
    const scalar y,
    const scalar majorRadius,
    const scalar minorRadius
) const
{
    const scalar radius = Foam::sqrt(Foam::sqr(x) + Foam::sqr(y));
    if (radius <= ROOTVSMALL)
    {
        return majorRadius;
    }

    const scalar cosTheta = x/radius;
    const scalar sinTheta = y/radius;

    return 1.0/Foam::sqrt
    (
        Foam::sqr(cosTheta/majorRadius)
      + Foam::sqr(sinTheta/minorRadius)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
