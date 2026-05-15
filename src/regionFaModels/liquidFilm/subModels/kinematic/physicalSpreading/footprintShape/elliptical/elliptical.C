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
#include "physicalSpreadingGeometry.H"

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
    if (aspectRatio_ < 1)
    {
        FatalIOErrorInFunction(dict)
            << "aspectRatio must be >= 1 for an elliptical footprint, "
            << "received " << aspectRatio_
            << exit(FatalIOError);
    }

    if (mag(direction_) <= ROOTVSMALL)
    {
        FatalIOErrorInFunction(dict)
            << "direction must be non-zero for an elliptical footprint"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

footprintShape::footprintGeometry elliptical::footprintAt
(
    const physicalSpreadingGeometry& geom,
    const label sourceGlobalFace,
    const scalar minorRadius
) const
{
    footprintGeometry fp;
    fp.centre = geom.centres()[sourceGlobalFace];
    fp.normal = geom.normals()[sourceGlobalFace];

    fp.majorAxis = direction_ - fp.normal*(direction_ & fp.normal);

    if (mag(fp.majorAxis) <= ROOTVSMALL)
    {
        FatalIOErrorInFunction(dict_)
            << "direction " << direction_
            << " is parallel to the normal of active global face "
            << sourceGlobalFace
            << "; its tangent-plane projection is degenerate"
            << exit(FatalIOError);
    }

    fp.majorAxis /= mag(fp.majorAxis);
    fp.minorAxis = fp.normal ^ fp.majorAxis;
    fp.minorAxis /= mag(fp.minorAxis);

    fp.majorRadius = aspectRatio_*minorRadius;
    fp.minorRadius = minorRadius;

    return fp;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
