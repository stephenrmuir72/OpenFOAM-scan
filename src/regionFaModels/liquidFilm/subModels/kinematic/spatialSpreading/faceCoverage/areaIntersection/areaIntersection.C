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
#include "spatialSpreadingData.H"

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

areaIntersection::overlapResult areaIntersection::calcOverlap
(
    const spatialSpreadingData& data,
    const footprintShape::footprintGeometry& fp,
    const label targetFacei
) const
{
    overlapResult result;

    // Create a 64-corner polygon for the parcel footprint
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

    // Get corners of the target face
    const pointField& targetFacePoints = data.facePoints()[targetFacei];

    // Fewer than three points, cannot define an area
    if (targetFacePoints.size() < 3)
    {
        result.centre = fp.centre;
        return result;
    }

    pointField projectedTargetPoints(targetFacePoints.size());
    face targetPolygon(targetFacePoints.size());
    vector projectedArea(Zero);

    // Project target vertex p onto the footprint plane
    forAll(targetFacePoints, pointi)
    {
        const point& p = targetFacePoints[pointi];
        projectedTargetPoints[pointi] =
            p - fp.normal*((p - fp.centre) & fp.normal);
        targetPolygon[pointi] = pointi;
    }

    forAll(projectedTargetPoints, pointi)
    {
        const label nextPoint =
            (pointi + 1) % projectedTargetPoints.size();

        projectedArea +=
            (projectedTargetPoints[pointi] - fp.centre)
          ^ (projectedTargetPoints[nextPoint] - fp.centre);
    }

    // Triangulate footprint and target polygons for 'faceAreaIntersect'
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
        targetPolygon,
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

    const bool reverseB = ((projectedArea & fp.normal) < 0);
    scalar area = 0;
    vector centroid(Zero);
    faceAreaIntersect intersection
    (
        footprintPoints,
        projectedTargetPoints,
        footprintTris,
        targetTris,
        reverseB
    );
    intersection.calc
    (
        footprintFace,
        targetPolygon,
        fp.normal,
        area,
        centroid  // area-weighted centre of the overlap
    );

    result.centre = centroid;
    result.area = Foam::max(area, scalar(0));

    return result;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

areaIntersection::areaIntersection(const dictionary& dict)
:
    faceCoverage(typeName, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

faceCoverage::coverageResult areaIntersection::evaluate
(
    const spatialSpreadingData& data,
    const footprintShape::footprintGeometry& fp,
    const label targetFacei,
    const scalar selectedDistance,
    const scalar boundaryRadius
) const
{
    coverageResult result;

    if
    (
        selectedDistance
      > (fp.enclosingRadius() + data.faceRadii()[targetFacei])
    )
    {
        return result;
    }

    const overlapResult overlap = calcOverlap(data, fp, targetFacei);

    if (overlap.area <= VSMALL)
    {
        return result;
    }

    result.included = true;
    result.measure = overlap.area;
    result.centre = overlap.centre;

    return result;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
