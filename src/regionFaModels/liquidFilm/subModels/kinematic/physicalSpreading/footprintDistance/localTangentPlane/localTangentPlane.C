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

#include "localTangentPlane.H"
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

defineTypeNameAndDebug(localTangentPlane, 0);
addToRunTimeSelectionTable(footprintDistance, localTangentPlane, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

localTangentPlane::localTangentPlane(const dictionary& dict)
:
    footprintDistance(typeName, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void localTangentPlane::findCandidates
(
    const physicalSpreadingGeometry& geom,
    const footprintShape& shape,
    const footprintShape::footprintGeometry& fp,
    const label sourceGlobalFace,
    const scalar searchRadius,
    DynamicList<label>& candidateFaces,
    DynamicList<scalar>& selectedDistances
) const
{
    const labelList candidates
    (
        geom.candidatesInSphere(fp.centre, searchRadius)
    );

    candidateFaces.clear();
    selectedDistances.clear();
    candidateFaces.reserve(candidates.size());
    selectedDistances.reserve(candidates.size());

    for (const label targetGlobalFace : candidates)
    {
        scalar x = 0;
        scalar y = 0;
        shape.planarOffset(fp, geom.centres()[targetGlobalFace], x, y);

        candidateFaces.append(targetGlobalFace);
        selectedDistances.append(Foam::sqrt(Foam::sqr(x) + Foam::sqr(y)));
    }
}


scalar localTangentPlane::kernelDistance
(
    const scalar selectedDistance,
    const scalar planarDistance,
    const scalar centroidPlanarDistance
) const
{
    return centroidPlanarDistance;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
