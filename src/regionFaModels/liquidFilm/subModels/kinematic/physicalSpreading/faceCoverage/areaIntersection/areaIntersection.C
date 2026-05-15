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

#include "areaIntersection.H"
#include "addToRunTimeSelectionTable.H"
#include "faceAreaIntersect.H"
#include "mathematicalConstants.H"
#include "physicalSpreadingGeometry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(areaIntersection, 0);
addToRunTimeSelectionTable(faceCoverage, areaIntersection, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar areaIntersection::overlapArea
(
    const physicalSpreadingGeometry& geom,
    const footprintShape::footprintGeometry& fp,
    const label targetGlobalFace,
    point& overlapCentre
) const
{
    static constexpr label nFootprintPoints = 64;

    pointField footprintPoints(nFootprintPoints);
    face footprintFace(nFootprintPoints);

    for (label pointi = 0; pointi < nFootprintPoints; ++pointi)
    {
        const scalar theta =
            constant::mathematical::twoPi
           *scalar(pointi)/scalar(nFootprintPoints);

        footprintPoints[pointi] =
            fp.centre
          + fp.majorRadius*Foam::cos(theta)*fp.majorAxis
          + fp.minorRadius*Foam::sin(theta)*fp.minorAxis;
        footprintFace[pointi] = pointi;
    }

    const pointField& targetFacePoints = geom.facePoints()[targetGlobalFace];

    // Fewer than three points, cannot define an area
    if (targetFacePoints.size() < 3)
    {
        overlapCentre = fp.centre;
        return 0;
    }

    pointField projectedTargetPoints(targetFacePoints.size());
    face targetFace(targetFacePoints.size());

    // Project target vertex p onto the footprint plane
    forAll(targetFacePoints, pointi)
    {
        const point& p = targetFacePoints[pointi];
        projectedTargetPoints[pointi] =
            p - fp.normal*((p - fp.centre) & fp.normal);
        targetFace[pointi] = pointi;
    }

    faceList footprintTriList;
    faceList targetTriList;
    faceAreaIntersect::triangulate
    (
        footprintFace,
        footprintPoints,
        faceAreaIntersect::tmFan,
        footprintTriList
    );
    faceAreaIntersect::triangulate
    (
        targetFace,
        projectedTargetPoints,
        faceAreaIntersect::tmMesh,
        targetTriList
    );

    DynamicList<face> footprintTris(footprintTriList.size());
    DynamicList<face> targetTris(targetTriList.size());

    for (const face& f : footprintTriList)
    {
        footprintTris.append(f);
    }
    for (const face& f : targetTriList)
    {
        targetTris.append(f);
    }

    scalar area = 0;
    vector centroid(Zero);
    faceAreaIntersect intersection
    (
        footprintPoints,
        projectedTargetPoints,
        footprintTris,
        targetTris
    );
    intersection.calc
    (
        footprintFace,
        targetFace,
        fp.normal,
        area,
        centroid
    );

    overlapCentre = centroid;
    return Foam::max(area, scalar(0));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

areaIntersection::areaIntersection(const dictionary& dict)
:
    faceCoverage(typeName, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

faceCoverage::coverageResult areaIntersection::evaluate
(
    const physicalSpreadingGeometry& geom,
    const footprintShape::footprintGeometry& fp,
    const label targetGlobalFace,
    const scalar selectedDistance,
    const scalar boundaryRadius
) const
{
    coverageResult result;

    if
    (
        selectedDistance
      > fp.enclosingRadius() + geom.faceRadii()[targetGlobalFace]
    )
    {
        return result;
    }

    point overlapCentre(geom.centres()[targetGlobalFace]);
    const scalar area = overlapArea(geom, fp, targetGlobalFace, overlapCentre);

    if (area <= VSMALL)
    {
        return result;
    }

    result.included = true;
    result.measure = area;
    result.centre = overlapCentre;

    return result;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
