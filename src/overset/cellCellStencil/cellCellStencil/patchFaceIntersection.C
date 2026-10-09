/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | www.openfoam.com
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Keysight Technologies
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "patchFaceIntersection.H"
#include "fvMesh.H"
#include "fvPatch.H"
#include "oversetFvPatch.H"
#include "polyMesh.H"
#include "triangle.H"

namespace Foam
{
namespace patchFaceIntersection
{

namespace
{

inline boundBox triangleBounds(const triPoints& tri)
{
    return boundBox(tri.box());
}


inline bool pointNearTriangle
(
    const point& p,
    const triPoints& tri,
    const scalar tolerance
)
{
    const triPointRef triRef(tri.a(), tri.b(), tri.c());
    const pointHit nearest(triRef.nearestPoint(p));
    return magSqr(nearest.point() - p) <= sqr(tolerance);
}


inline bool segmentIntersectsTriangle
(
    const point& p0,
    const point& p1,
    const triPoints& tri,
    const scalar tolerance
)
{
    if
    (
        pointNearTriangle(p0, tri, tolerance)
     || pointNearTriangle(p1, tri, tolerance)
    )
    {
        return true;
    }

    const vector edge(p1 - p0);
    const scalar edgeLength(mag(edge));

    if (edgeLength <= VSMALL)
    {
        return false;
    }

    const triPointRef triRef(tri.a(), tri.b(), tri.c());
    const scalar relativeTolerance(tolerance/edgeLength);
    const pointHit hit
    (
        triRef.intersection
        (
            p0,
            edge,
            intersection::FULL_RAY,
            relativeTolerance
        )
    );

    return
        hit.hit()
     && hit.distance() >= -relativeTolerance
     && hit.distance() <= 1 + relativeTolerance;
}


scalar segmentSegmentDistanceSqr
(
    const point& p1,
    const point& q1,
    const point& p2,
    const point& q2
)
{
    const vector d1(q1 - p1);
    const vector d2(q2 - p2);
    const vector r(p1 - p2);
    const scalar a(d1 & d1);
    const scalar e(d2 & d2);
    const scalar f(d2 & r);

    scalar s = 0;
    scalar t = 0;

    if (a <= VSMALL && e <= VSMALL)
    {
        return magSqr(r);
    }

    if (a <= VSMALL)
    {
        t = max(scalar(0), min(scalar(1), f/e));
    }
    else
    {
        const scalar c(d1 & r);

        if (e <= VSMALL)
        {
            s = max(scalar(0), min(scalar(1), -c/a));
        }
        else
        {
            const scalar b(d1 & d2);
            const scalar denom(a*e - sqr(b));

            if (denom > VSMALL)
            {
                s = max
                (
                    scalar(0),
                    min(scalar(1), (b*f - c*e)/denom)
                );
            }

            t = (b*s + f)/e;

            if (t < 0)
            {
                t = 0;
                s = max(scalar(0), min(scalar(1), -c/a));
            }
            else if (t > 1)
            {
                t = 1;
                s = max(scalar(0), min(scalar(1), (b - c)/a));
            }
        }
    }

    return magSqr((p1 + s*d1) - (p2 + t*d2));
}


bool voxelRange
(
    const boundBox& searchBb,
    const labelVector& nDivs,
    const boundBox& queryBb,
    labelVector& minIds,
    labelVector& maxIds
)
{
    if (!searchBb.overlaps(queryBb))
    {
        return false;
    }

    const vector span(searchBb.span());

    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        if (nDivs[cmpt] <= 0)
        {
            FatalErrorInFunction
                << "Invalid voxel divisions " << nDivs
                << exit(FatalError);
        }

        if (span[cmpt] <= VSMALL)
        {
            minIds[cmpt] = 0;
            maxIds[cmpt] = nDivs[cmpt] - 1;
        }
        else
        {
            minIds[cmpt] = label
            (
                floor
                (
                    (queryBb.min()[cmpt] - searchBb.min()[cmpt])
                   /span[cmpt]*nDivs[cmpt]
                )
            );
            maxIds[cmpt] = label
            (
                floor
                (
                    (queryBb.max()[cmpt] - searchBb.min()[cmpt])
                   /span[cmpt]*nDivs[cmpt]
                )
            );

            minIds[cmpt] = max(label(0), minIds[cmpt]);
            maxIds[cmpt] = min(nDivs[cmpt] - 1, maxIds[cmpt]);
        }

        if (minIds[cmpt] > maxIds[cmpt])
        {
            return false;
        }
    }

    return true;
}


triPoints faceTriangle
(
    const face& f,
    const UList<point>& points
)
{
    return triPoints(points, f[0], f[1], f[2]);
}


void forFaceTriangles
(
    const face& f,
    const UList<point>& points,
    const boundBox& searchBb,
    const labelVector& nDivs,
    const labelVector& strides,
    const vector& tolerance,
    const PackedList<2>* voxels,
    PackedList<2>* writeVoxels,
    const unsigned int value,
    bool& found
)
{
    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        if (nDivs[cmpt] <= 0)
        {
            FatalErrorInFunction
                << "Invalid voxel divisions " << nDivs
                << exit(FatalError);
        }
    }

    faceList triangleFaces(f.nTriangles());
    label triI = 0;
    f.triangles(points, triI, triangleFaces);

    const vector width
    (
        searchBb.span()[0]/nDivs[0],
        searchBb.span()[1]/nDivs[1],
        searchBb.span()[2]/nDivs[2]
    );

    forAll(triangleFaces, triFacei)
    {
        const triPoints tri(faceTriangle(triangleFaces[triFacei], points));
        boundBox triBb(triangleBounds(tri));
        triBb.grow(tolerance);

        labelVector minIds;
        labelVector maxIds;
        if (!voxelRange(searchBb, nDivs, triBb, minIds, maxIds))
        {
            continue;
        }

        for (label i = minIds[0]; i <= maxIds[0]; ++i)
        {
            for (label j = minIds[1]; j <= maxIds[1]; ++j)
            {
                for (label k = minIds[2]; k <= maxIds[2]; ++k)
                {
                    const point voxelMin
                    (
                        searchBb.min()[0] + i*width[0],
                        searchBb.min()[1] + j*width[1],
                        searchBb.min()[2] + k*width[2]
                    );
                    const boundBox voxelBb(voxelMin, voxelMin + width);

                    if (!triangleIntersectsBox(tri, voxelBb, tolerance))
                    {
                        continue;
                    }

                    const label voxeli
                    (
                        i*strides[0] + j*strides[1] + k*strides[2]
                    );

                    if (writeVoxels)
                    {
                        (*writeVoxels)[voxeli] = value;
                    }
                    else if (voxels && (*voxels)[voxeli] == value)
                    {
                        found = true;
                        return;
                    }
                }
            }
        }
    }
}

} // End anonymous namespace


bool triangleIntersectsBox
(
    const triPoints& tri,
    const boundBox& box,
    const vector& tolerance
)
{
    boundBox expandedBox(box);
    expandedBox.grow(tolerance);

    if (!expandedBox.overlaps(triangleBounds(tri)))
    {
        return false;
    }

    const point centre(0.5*(expandedBox.min() + expandedBox.max()));
    const vector half(0.5*expandedBox.span());

    const vector edges[3] =
    {
        tri.b() - tri.a(),
        tri.c() - tri.b(),
        tri.a() - tri.c()
    };
    const vector boxAxes[3] =
    {
        vector(1, 0, 0),
        vector(0, 1, 0),
        vector(0, 0, 1)
    };

    const auto separated =
    [&tri, &centre, &half](const vector& axis) -> bool
    {
        if (magSqr(axis) <= VSMALL)
        {
            return false;
        }

        scalar triMin(axis & tri.a());
        scalar triMax(triMin);
        for (const point& p : {tri.b(), tri.c()})
        {
            const scalar proj(axis & p);
            triMin = min(triMin, proj);
            triMax = max(triMax, proj);
        }

        const scalar boxCentre(axis & centre);
        const scalar boxRadius
        (
            half[0]*mag(axis[0])
          + half[1]*mag(axis[1])
          + half[2]*mag(axis[2])
        );

        return
            triMin > boxCentre + boxRadius
         || triMax < boxCentre - boxRadius;
    };

    for (const vector& axis : boxAxes)
    {
        if (separated(axis))
        {
            return false;
        }
    }

    if (separated(edges[0] ^ edges[1]))
    {
        return false;
    }

    for (const vector& edge : edges)
    {
        for (const vector& axis : boxAxes)
        {
            if (separated(edge ^ axis))
            {
                return false;
            }
        }
    }

    return true;
}


bool trianglesIntersect
(
    const triPoints& triA,
    const triPoints& triB,
    const scalar tolerance
)
{
    boundBox bbA(triangleBounds(triA));
    boundBox bbB(triangleBounds(triB));
    const vector toleranceVec(tolerance, tolerance, tolerance);
    bbA.grow(toleranceVec);

    if (!bbA.overlaps(bbB))
    {
        return false;
    }

    for (const point& p : triA)
    {
        if (pointNearTriangle(p, triB, tolerance))
        {
            return true;
        }
    }
    for (const point& p : triB)
    {
        if (pointNearTriangle(p, triA, tolerance))
        {
            return true;
        }
    }

    const point triAVertices[3] = {triA.a(), triA.b(), triA.c()};
    const point triBVertices[3] = {triB.a(), triB.b(), triB.c()};

    for (direction edgeA = 0; edgeA < 3; ++edgeA)
    {
        const point& a0 = triAVertices[edgeA];
        const point& a1 = triAVertices[(edgeA + 1) % 3];

        if (segmentIntersectsTriangle(a0, a1, triB, tolerance))
        {
            return true;
        }

        for (direction edgeB = 0; edgeB < 3; ++edgeB)
        {
            if
            (
                segmentSegmentDistanceSqr
                (
                    a0,
                    a1,
                    triBVertices[edgeB],
                    triBVertices[(edgeB + 1) % 3]
                )
             <= sqr(tolerance)
            )
            {
                return true;
            }
        }
    }

    for (direction edgeB = 0; edgeB < 3; ++edgeB)
    {
        if
        (
            segmentIntersectsTriangle
            (
                triBVertices[edgeB],
                triBVertices[(edgeB + 1) % 3],
                triA,
                tolerance
            )
        )
        {
            return true;
        }
    }

    return false;
}


void fillFaceVoxels
(
    PackedList<2>& voxels,
    const boundBox& searchBb,
    const labelVector& nDivs,
    const labelVector& strides,
    const face& f,
    const UList<point>& points,
    const vector& tolerance,
    const unsigned int value
)
{
    bool found = false;
    forFaceTriangles
    (
        f,
        points,
        searchBb,
        nDivs,
        strides,
        tolerance,
        nullptr,
        &voxels,
        value,
        found
    );
}


bool overlapsFaceVoxels
(
    const boundBox& searchBb,
    const labelVector& nDivs,
    const labelVector& strides,
    const face& f,
    const UList<point>& points,
    const vector& tolerance,
    const PackedList<2>& voxels,
    const unsigned int value
)
{
    bool found = false;
    forFaceTriangles
    (
        f,
        points,
        searchBb,
        nDivs,
        strides,
        tolerance,
        &voxels,
        nullptr,
        value,
        found
    );
    return found;
}


pointField patchTriangles(const fvMesh& mesh)
{
    const fvBoundaryMesh& boundary = mesh.boundary();
    label nTriangles = 0;

    forAll(boundary, patchi)
    {
        const fvPatch& fvp = boundary[patchi];
        if
        (
            fvPatch::constraintType(fvp.type())
         || isA<oversetFvPatch>(fvp)
        )
        {
            continue;
        }

        const polyPatch& pp = fvp.patch();
        forAll(pp, facei)
        {
            nTriangles += pp[facei].nTriangles();
        }
    }

    pointField result(3*nTriangles);
    label trianglei = 0;

    forAll(boundary, patchi)
    {
        const fvPatch& fvp = boundary[patchi];
        if
        (
            fvPatch::constraintType(fvp.type())
         || isA<oversetFvPatch>(fvp)
        )
        {
            continue;
        }

        const polyPatch& pp = fvp.patch();
        forAll(pp, facei)
        {
            const face& f = pp[facei];
            faceList triangleFaces(f.nTriangles());
            label triI = 0;
            f.triangles(pp.points(), triI, triangleFaces);

            forAll(triangleFaces, triFacei)
            {
                const face& triFace = triangleFaces[triFacei];
                for (direction pointi = 0; pointi < 3; ++pointi)
                {
                    result[3*trianglei + pointi] =
                        pp.points()[triFace[pointi]];
                }
                ++trianglei;
            }
        }
    }

    return result;
}


pointField trianglesOverlapping
(
    const pointField& trianglePoints,
    const boundBox& targetBb,
    const vector& tolerance
)
{
    if (trianglePoints.size() % 3)
    {
        FatalErrorInFunction
            << "Triangle point count " << trianglePoints.size()
            << " is not divisible by three"
            << exit(FatalError);
    }

    label nPoints = 0;
    for (label pointi = 0; pointi < trianglePoints.size(); pointi += 3)
    {
        const triPoints tri
        (
            trianglePoints[pointi],
            trianglePoints[pointi + 1],
            trianglePoints[pointi + 2]
        );
        boundBox triBb(triangleBounds(tri));
        triBb.grow(tolerance);
        if (triBb.overlaps(targetBb))
        {
            nPoints += 3;
        }
    }

    pointField result(nPoints);
    label outPointi = 0;

    for (label pointi = 0; pointi < trianglePoints.size(); pointi += 3)
    {
        const triPoints tri
        (
            trianglePoints[pointi],
            trianglePoints[pointi + 1],
            trianglePoints[pointi + 2]
        );
        boundBox triBb(triangleBounds(tri));
        triBb.grow(tolerance);
        if (triBb.overlaps(targetBb))
        {
            for (direction pointi2 = 0; pointi2 < 3; ++pointi2)
            {
                result[outPointi++] = trianglePoints[pointi + pointi2];
            }
        }
    }

    return result;
}


bool intersectsCell
(
    const pointField& trianglePoints,
    const fvMesh& mesh,
    const label celli,
    const scalar tolerance
)
{
    if (trianglePoints.size() % 3)
    {
        FatalErrorInFunction
            << "Triangle point count " << trianglePoints.size()
            << " is not divisible by three"
            << exit(FatalError);
    }

    if (trianglePoints.empty())
    {
        return false;
    }

    const vector toleranceVec(tolerance, tolerance, tolerance);
    boundBox cellBb(mesh.cellBb(celli));
    cellBb.grow(toleranceVec);

    for (label pointi = 0; pointi < trianglePoints.size(); pointi += 3)
    {
        const triPoints tri
        (
            trianglePoints[pointi],
            trianglePoints[pointi + 1],
            trianglePoints[pointi + 2]
        );
        boundBox triBb(triangleBounds(tri));
        triBb.grow(toleranceVec);

        if (!cellBb.overlaps(triBb))
        {
            continue;
        }

        for (const point& p : tri)
        {
            if (mesh.pointInCell(p, celli, polyMesh::CELL_TETS))
            {
                return true;
            }
        }
    }

    label nCellTriangles = 0;
    const cell& c = mesh.cells()[celli];
    const faceList& faces = mesh.faces();

    forAll(c, i)
    {
        nCellTriangles += faces[c[i]].nTriangles();
    }

    List<triPoints> cellTriangles(nCellTriangles);
    label cellTrianglei = 0;
    const pointField& meshPoints = mesh.points();

    forAll(c, i)
    {
        const face& f = faces[c[i]];
        faceList triangleFaces(f.nTriangles());
        label triI = 0;
        f.triangles(meshPoints, triI, triangleFaces);

        forAll(triangleFaces, triFacei)
        {
            const face& triFace = triangleFaces[triFacei];
            cellTriangles[cellTrianglei++] =
                triPoints
                (
                    meshPoints,
                    triFace[0],
                    triFace[1],
                    triFace[2]
                );
        }
    }

    for (label pointi = 0; pointi < trianglePoints.size(); pointi += 3)
    {
        const triPoints tri
        (
            trianglePoints[pointi],
            trianglePoints[pointi + 1],
            trianglePoints[pointi + 2]
        );
        boundBox triBb(triangleBounds(tri));
        triBb.grow(toleranceVec);

        if (!cellBb.overlaps(triBb))
        {
            continue;
        }

        for (const triPoints& cellTri : cellTriangles)
        {
            if (trianglesIntersect(tri, cellTri, tolerance))
            {
                return true;
            }
        }
    }

    return false;
}


bool checkedVoxelCount
(
    const labelVector& nDivs,
    const label maxVoxels,
    label& nVoxels
)
{
    nVoxels = 0;

    if (maxVoxels <= 0)
    {
        return false;
    }

    label candidateCount = 1;

    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        if
        (
            nDivs[cmpt] <= 0
         || candidateCount > maxVoxels/nDivs[cmpt]
        )
        {
            return false;
        }

        candidateCount *= nDivs[cmpt];
    }

    nVoxels = candidateCount;
    return true;
}


bool refineVoxelDivisions
(
    const labelVector& nDivs,
    const labelVector& geometricD,
    const label maxVoxels,
    labelVector& refinedDivs,
    label& refinedVoxelCount
)
{
    label currentVoxelCount = 0;
    if (!checkedVoxelCount(nDivs, maxVoxels, currentVoxelCount))
    {
        refinedDivs = nDivs;
        refinedVoxelCount = 0;
        return false;
    }

    refinedDivs = nDivs;

    forAll(geometricD, cmpt)
    {
        if (geometricD[cmpt] != -1)
        {
            refinedDivs[cmpt] *= 2;
        }
    }

    if (refinedDivs == nDivs)
    {
        refinedVoxelCount = currentVoxelCount;
        return false;
    }

    return checkedVoxelCount(refinedDivs, maxVoxels, refinedVoxelCount);
}

} // End namespace patchFaceIntersection
} // End namespace Foam

// ************************************************************************* //
