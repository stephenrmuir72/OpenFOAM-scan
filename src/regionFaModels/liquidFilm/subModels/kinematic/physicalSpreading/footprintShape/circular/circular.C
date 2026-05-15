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

#include "circular.H"
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

defineTypeNameAndDebug(circular, 0);
addToRunTimeSelectionTable(footprintShape, circular, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

circular::circular(const dictionary& dict)
:
    footprintShape(typeName, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

footprintShape::footprintGeometry circular::footprintAt
(
    const physicalSpreadingGeometry& geom,
    const label sourceGlobalFace,
    const scalar minorRadius
) const
{
    footprintGeometry fp;
    fp.centre = geom.centres()[sourceGlobalFace];
    fp.normal = geom.normals()[sourceGlobalFace];

    vector reference(1, 0, 0);
    if (mag(fp.normal.y()) < mag(fp.normal.x()))
    {
        reference = vector(0, 1, 0);
    }
    if
    (
        mag(fp.normal.z()) < mag(fp.normal.x())
     && mag(fp.normal.z()) < mag(fp.normal.y())
    )
    {
        reference = vector(0, 0, 1);
    }

    fp.majorAxis = reference - fp.normal*(reference & fp.normal);
    fp.majorAxis /= mag(fp.majorAxis);
    fp.minorAxis = fp.normal ^ fp.majorAxis;
    fp.minorAxis /= mag(fp.minorAxis);

    fp.majorRadius = minorRadius;
    fp.minorRadius = minorRadius;

    return fp;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
