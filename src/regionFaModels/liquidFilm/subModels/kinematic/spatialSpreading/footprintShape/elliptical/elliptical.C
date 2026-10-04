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

#include "elliptical.H"
#include "addToRunTimeSelectionTable.H"
#include "spatialSpreadingData.H"

#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(elliptical, 0);
addToRunTimeSelectionTable(footprintShape, elliptical, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

elliptical::elliptical(const dictionary& dict)
:
    footprintShape(typeName, dict),
    aspectRatio_(dict.get<scalar>("aspectRatio")),
    direction_(dict.get<vector>("direction"))
{
    if (!(aspectRatio_ >= 1) || !std::isfinite(aspectRatio_))
    {
        FatalIOErrorInFunction(dict)
            << "aspectRatio must be finite and >= 1 for an elliptical "
            << "footprint, "
            << "received " << aspectRatio_
            << exit(FatalIOError);
    }

    const scalar directionMagnitude = Foam::mag(direction_);
    if
    (
        !std::isfinite(directionMagnitude)
     || directionMagnitude <= ROOTVSMALL
    )
    {
        FatalIOErrorInFunction(dict)
            << "direction must be finite and non-zero for an elliptical "
            << "footprint"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

footprintShape::footprintGeometry elliptical::footprintAt
(
    const spatialSpreadingData& data,
    const label sourceFace,
    const scalar minorRadius
) const
{
    footprintGeometry fp;
    fp.centre = data.centres()[sourceFace];
    fp.normal = data.normals()[sourceFace];

    fp.majorAxis = direction_ - fp.normal*(direction_ & fp.normal);

    const scalar majorAxisMagnitude = Foam::mag(fp.majorAxis);
    if
    (
        !std::isfinite(majorAxisMagnitude)
     || majorAxisMagnitude <= ROOTVSMALL
    )
    {
        FatalIOErrorInFunction(dict_)
            << "direction " << direction_
            << " is parallel to the normal of active source face "
            << sourceFace
            << "; its tangent-plane projection is degenerate"
            << exit(FatalIOError);
    }

    fp.majorAxis /= majorAxisMagnitude;
    fp.minorAxis = fp.normal^fp.majorAxis;
    const scalar minorAxisMagnitude = Foam::mag(fp.minorAxis);
    if (!std::isfinite(minorAxisMagnitude) || minorAxisMagnitude <= ROOTVSMALL)
    {
        FatalIOErrorInFunction(dict_)
            << "Cannot construct a finite tangent-plane basis for source "
            << "face " << sourceFace
            << exit(FatalIOError);
    }
    fp.minorAxis /= minorAxisMagnitude;

    fp.majorRadius = aspectRatio_*minorRadius;
    fp.minorRadius = minorRadius;

    return fp;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
