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

faceCoverage::coverageResult areaIntersection::calcOverlap
(
    const spatialSpreadingData& data,
    const footprintShape::footprintGeometry& fp,
    const label targetFacei
) const
{
    // Build a 64-point outline of the footprint, flatten the
    // target face onto the footprint's plane, triangulate both shapes,
    // compute their intersection, and return its area and centre.

    coverageResult result;

    // Create a 64-corner polygon for the parcel footprint
    static constexpr label nFootprintPoints = 64;

    pointField footprintPoints(nFootprintPoints);
    // Indices connecting footprintPoints into one polygon
    face footprintFace(nFootprintPoints);

    // Place each coordinate around the footprint
    for (label pointi = 0; pointi < nFootprintPoints; ++pointi)
    {
        const scalar theta =
            constant::mathematical::twoPi
           *scalar(pointi)/scalar(nFootprintPoints);

        // Standard ellipse parameterisation
        footprintPoints[pointi] =
            fp.centre
          + fp.majorRadius*Foam::cos(theta)*fp.majorAxis
          + fp.minorRadius*Foam::sin(theta)*fp.minorAxis;

        footprintFace[pointi] = pointi;
    }

    // Get corner coords of the target face
    const pointField& targetFacePoints = data.facePoints()[targetFacei];

    // Fewer than three points, polygon cannot define an area
    if (targetFacePoints.size() < 3)
    {
        result.centre = fp.centre;
        return result;  // with zero area
    }

    // Allocate storage for target face's projected points and polygon
    pointField projectedTargetPoints(targetFacePoints.size());
    face targetPolygon(targetFacePoints.size());
    vector projectedArea(zero{});

    // Project target vertex p onto the footprint plane
    forAll(targetFacePoints, pointi)
    {
        const point& p = targetFacePoints[pointi];

        // Remove the component in the fp.normal direction
        projectedTargetPoints[pointi] =
            p - fp.normal*((p - fp.centre).inner(fp.normal));

        targetPolygon[pointi] = pointi;
    }

    // Determine projected polygon's rotational-order direction, CW vs CCW
    forAll(projectedTargetPoints, pointi)
    {
        const label nextPoint = (pointi + 1) % projectedTargetPoints.size();

        // Only the sign of the plane normal is needed later
        projectedArea +=
            (projectedTargetPoints[pointi] - fp.centre)
            .cross(projectedTargetPoints[nextPoint] - fp.centre);
    }

    // Triangulate footprint and target polygons to use 'faceAreaIntersect'
    faceList targetTriList;
    DynamicList<face> footprintTris(nFootprintPoints - 2);
    faceAreaIntersect::triangleFan(footprintFace, footprintTris);
    faceAreaIntersect::triangulate
    (
        targetPolygon,
        projectedTargetPoints,
        faceAreaIntersect::tmMesh,  // triangulation mode
        targetTriList
    );

    // Copy target triangulation into input type expected by 'faceAreaIntersect'
    DynamicList<face> targetTris(targetTriList.size());

    for (const face& f : targetTriList)
    {
        targetTris.append(f);
    }

    // reverseB is true if projected target polygon's face normal is opposite
    // to fp.normal; if so, reverse faces of B which is the target polygon
    const bool reverseB = (projectedArea.inner(fp.normal) < 0);

    scalar area = 0;
    vector centroid(zero{});

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
        area,     // output
        centroid  // output: area-weighted centre of the overlap
    );

    result.measure = Foam::max(area, scalar(0));
    result.centre = centroid;

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
    const scalar,
    const scalar
) const
{
    coverageResult result;

    scalar x = 0;
    scalar y = 0;
    footprintShape::planarOffset
    (
        fp,
        data.centres()[targetFacei],
        x,
        y
    );
    const scalar candidateRadius =
        fp.enclosingRadius() + data.faceRadii()[targetFacei];

    if (Foam::hypot(x, y) > candidateRadius)
    {
        return result;
    }

    coverageResult overlap(calcOverlap(data, fp, targetFacei));

    if (overlap.measure <= VSMALL)
    {
        return result;  // with zero area
    }

    overlap.included = true;

    return overlap;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
