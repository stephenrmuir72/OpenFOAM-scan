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
#include "tetIndices.H"
#include "triangle.H"

#include <cmath>
#include <algorithm>

namespace Foam
{
namespace patchFaceIntersection
{

class triangleData
{
    const pointField& points_;

public:

    static constexpr int nDim() noexcept { return 3; }

    explicit triangleData(const pointField& points)
    :
        points_(points)
    {}

    bool empty() const { return points_.empty(); }
    label size() const { return points_.size()/3; }
    label objectIndex(const label index) const noexcept { return index; }

    const point& centre(const label index) const
    {
        return points_[3*index];
    }

    treeBoundBox bounds(const labelUList& indices) const
    {
        treeBoundBox result;
        for (const label index : indices)
        {
            result.min() = min(result.min(), bounds(index).min());
            result.max() = max(result.max(), bounds(index).max());
        }
        return result;
    }

    bool overlaps(const label index, const treeBoundBox& box) const
    {
        return bounds(index).overlaps(box);
    }

    bool overlaps
    (
        const label index,
        const point& centrePoint,
        const scalar radiusSqr
    ) const
    {
        return bounds(index).overlaps(centrePoint, radiusSqr);
    }

private:

    treeBoundBox bounds(const label index) const
    {
        treeBoundBox bb
        (
            points_[3*index],
            points_[3*index]
        );
        bb.min() = min(bb.min(), points_[3*index + 1]);
        bb.max() = max(bb.max(), points_[3*index + 1]);
        bb.min() = min(bb.min(), points_[3*index + 2]);
        bb.max() = max(bb.max(), points_[3*index + 2]);
        return bb;
    }
};

namespace
{

inline void validateTolerance(const scalar tolerance)
{
    if (!std::isfinite(tolerance) || tolerance < 0)
    {
        FatalErrorInFunction
            << "Geometry tolerance must be finite and non-negative, got "
            << tolerance << exit(FatalError);
    }
}

inline void validateTolerance(const vector& tolerance)
{
    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        if (!std::isfinite(tolerance[cmpt]) || tolerance[cmpt] < 0)
        {
            FatalErrorInFunction
                << "Geometry tolerance must be finite and non-negative, got "
                << tolerance << exit(FatalError);
        }
    }
}

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
    validateTolerance(tolerance);
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
    validateTolerance(tolerance);
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
    const pointHit hit
    (
        triRef.intersection
        (
            p0,
            edge,
            intersection::FULL_RAY,
            0
        )
    );

    return
        hit.hit()
     && hit.distance() >= 0
     && hit.distance() <= 1;
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


} // End anonymous namespace


bool closedVoxelRange
(
    const boundBox& searchBb,
    const labelVector& nDivs,
    const boundBox& queryBb,
    labelVector& minIds,
    labelVector& maxIds
)
{
    boundBox closedQueryBb(queryBb);
    const scalar closureTol
    (
        10*SMALL*max(searchBb.mag(), scalar(1))
    );
    closedQueryBb.grow(vector(closureTol, closureTol, closureTol));

    if (!searchBb.overlaps(closedQueryBb))
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
            const scalar lower
            (
                (closedQueryBb.min()[cmpt] - searchBb.min()[cmpt])
               /span[cmpt]*nDivs[cmpt]
            );
            const scalar upper
            (
                (closedQueryBb.max()[cmpt] - searchBb.min()[cmpt])
               /span[cmpt]*nDivs[cmpt]
            );
            const scalar clippedLower
            (
                max(scalar(0), min(lower, scalar(nDivs[cmpt])))
            );
            const scalar clippedUpper
            (
                max(scalar(0), min(upper, scalar(nDivs[cmpt])))
            );
            minIds[cmpt] = label(floor(clippedLower));
            maxIds[cmpt] = label(floor(clippedUpper));

            const scalar lowerNearest(round(clippedLower));
            if
            (
                clippedLower > 0
             && clippedLower < nDivs[cmpt]
             && mag(clippedLower - lowerNearest) <= 10*SMALL
            )
            {
                --minIds[cmpt];
            }

            minIds[cmpt] = max(label(0), min(nDivs[cmpt] - 1, minIds[cmpt]));
            maxIds[cmpt] = max(label(0), min(nDivs[cmpt] - 1, maxIds[cmpt]));
        }

        if (minIds[cmpt] > maxIds[cmpt])
        {
            return false;
        }
    }

    return true;
}

namespace
{

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
    validateTolerance(tolerance);

    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        if (nDivs[cmpt] <= 0)
        {
            FatalErrorInFunction
                << "Invalid voxel divisions " << nDivs
                << exit(FatalError);
        }
    }

    const scalar toleranceMax
    (
        max(tolerance.x(), max(tolerance.y(), tolerance.z()))
    );
    const vector rasterTolerance
    (
        toleranceMax,
        toleranceMax,
        toleranceMax
    );

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
        triBb.grow(rasterTolerance);

        labelVector minIds;
        labelVector maxIds;
        if (!closedVoxelRange(searchBb, nDivs, triBb, minIds, maxIds))
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

                    if
                    (
                        !triangleIntersectsBox
                        (
                            tri,
                            voxelBb,
                            rasterTolerance
                        )
                    )
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
    validateTolerance(tolerance);
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
    validateTolerance(tolerance);
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
    validateTolerance(tolerance);
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
    validateTolerance(tolerance);
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


triangleSearch::triangleSearch
(
    const pointField& trianglePoints,
    const vector& tolerance
)
:
    trianglePoints_(trianglePoints),
    treePtr_(nullptr),
    tolerance_(tolerance)
{
    validateTolerance(tolerance);

    if (trianglePoints_.size() % 3)
    {
        FatalErrorInFunction
            << "Triangle point count " << trianglePoints_.size()
            << " is not divisible by three"
            << exit(FatalError);
    }

    if (trianglePoints_.empty())
    {
        return;
    }

    treeBoundBox bb
    (
        trianglePoints_[0],
        trianglePoints_[0]
    );
    for (const point& p : trianglePoints_)
    {
        bb.min() = min(bb.min(), p);
        bb.max() = max(bb.max(), p);
    }
    const vector span(bb.span());
    if (cmptMin(span) <= 1e-8*max(cmptMax(span), VSMALL))
    {
        // Flat geometry on an octree split plane is duplicated into both
        // sides and can stop subdivision. Median splits in the varying axes
        // store each triangle once, with full 3-D bounds for conservative queries.
        const label nTriangles = size();
        triangleBounds_.setSize(nTriangles);
        triangleOrder_.setSize(nTriangles);
        for (label trianglei = 0; trianglei < nTriangles; ++trianglei)
        {
            const triPoints tri
            (
                trianglePoints_[3*trianglei],
                trianglePoints_[3*trianglei + 1],
                trianglePoints_[3*trianglei + 2]
            );
            triangleBounds_[trianglei] = triangleBounds(tri);
            triangleOrder_[trianglei] = trianglei;
        }
        boundsNodes_.reserve(nTriangles/2 + 1);
        buildBoundsTree(0, nTriangles);
        return;
    }
    const scalar toleranceMax
    (
        max(tolerance_.x(), max(tolerance_.y(), tolerance_.z()))
    );
    bb.grow(toleranceMax);
    const scalar extent
    (
        max(toleranceMax, 10*SMALL*max(bb.mag(), scalar(1)))
    );
    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        if (bb.span(cmpt) <= VSMALL)
        {
            bb.min()[cmpt] -= extent;
            bb.max()[cmpt] += extent;
        }
    }

    treePtr_.reset
    (
        new indexedOctree<triangleData>
        (
            triangleData(trianglePoints_),
            bb,
            10,
            10,
            2.0
        )
    );
}


triangleSearch::triangleSearch
(
    const pointField& trianglePoints,
    const scalar tolerance
)
:
    triangleSearch
    (
        trianglePoints,
        vector(tolerance, tolerance, tolerance)
    )
{}


triangleSearch::~triangleSearch() = default;


label triangleSearch::buildBoundsTree(const label begin, const label end)
{
    boundsNode node;
    node.begin = begin;
    node.end = end;
    node.bb = triangleBounds_[triangleOrder_[begin]];
    boundBox centres
    (
        node.bb.centre(),
        node.bb.centre()
    );
    for (label i = begin + 1; i < end; ++i)
    {
        const boundBox& bb = triangleBounds_[triangleOrder_[i]];
        node.bb.min() = min(node.bb.min(), bb.min());
        node.bb.max() = max(node.bb.max(), bb.max());
        centres.min() = min(centres.min(), bb.centre());
        centres.max() = max(centres.max(), bb.centre());
    }
    const label nodei = boundsNodes_.size();
    boundsNodes_.append(node);
    if (end - begin > 10)
    {
        const vector span(centres.span());
        direction axis = 0;
        for (direction cmpt = 1; cmpt < 3; ++cmpt)
        {
            if (span[cmpt] > span[axis])
            {
                axis = cmpt;
            }
        }
        const label middle = begin + (end - begin)/2;
        std::nth_element
        (
            triangleOrder_.begin() + begin,
            triangleOrder_.begin() + middle,
            triangleOrder_.begin() + end,
            [&](const label a, const label b)
            {
                const scalar ca = triangleBounds_[a].centre()[axis];
                const scalar cb = triangleBounds_[b].centre()[axis];
                return ca < cb || (ca == cb && a < b);
            }
        );
        const label left = buildBoundsTree(begin, middle);
        const label right = buildBoundsTree(middle, end);
        boundsNodes_[nodei].left = left;
        boundsNodes_[nodei].right = right;
    }
    return nodei;
}


void triangleSearch::queryBoundsTree
(
    const label nodei,
    const boundBox& bb,
    DynamicList<label>& candidates,
    queryStats& stats
) const
{
    ++stats.nodes;
    const boundsNode& node = boundsNodes_[nodei];
    if (!node.bb.overlaps(bb))
    {
        return;
    }
    if (node.left >= 0)
    {
        queryBoundsTree(node.left, bb, candidates, stats);
        queryBoundsTree(node.right, bb, candidates, stats);
    }
    else
    {
        for (label i = node.begin; i < node.end; ++i)
        {
            ++stats.triangleTests;
            const label trianglei = triangleOrder_[i];
            if (triangleBounds_[trianglei].overlaps(bb))
            {
                candidates.append(trianglei);
            }
        }
    }
}


void triangleSearch::queryOctree
(
    const label nodei,
    const treeBoundBox& bb,
    labelHashSet& candidates,
    queryStats& stats
) const
{
    using treeType = indexedOctree<triangleData>;
    ++stats.nodes;
    const auto& tree = treePtr_();
    const auto& node = tree.nodes()[nodei];
    for (direction octant = 0; octant < 8; ++octant)
    {
        const labelBits index = node.subNodes_[octant];
        if (treeType::isNode(index))
        {
            const label child = treeType::getNode(index);
            if (tree.nodes()[child].bb_.overlaps(bb))
            {
                queryOctree(child, bb, candidates, stats);
            }
        }
        else if
        (
            treeType::isContent(index)
         && node.bb_.subOverlaps(octant, bb)
        )
        {
            for (const label trianglei : tree.getIndices(index))
            {
                ++stats.triangleTests;
                if (tree.shapes().overlaps(trianglei, bb))
                {
                    candidates.insert(trianglei);
                }
            }
        }
    }
}


labelList triangleSearch::query
(
    const boundBox& targetBb,
    queryStats* stats
) const
{
    queryStats work;
    boundBox queryBb(targetBb);
    queryBb.grow(cmptMax(tolerance_));
    labelList result;
    if (!boundsNodes_.empty())
    {
        DynamicList<label> candidates;
        queryBoundsTree(0, queryBb, candidates, work);
        result.transfer(candidates);
    }
    else if (treePtr_.valid() && treePtr_().bb().overlaps(queryBb))
    {
        // Unlike indexedOctree::findBox(), reserve storage only for results,
        // not in proportion to the entire source geometry.
        labelHashSet candidates;
        queryOctree(0, treeBoundBox(queryBb), candidates, work);
        result = candidates.toc();
    }
    if (stats)
    {
        *stats = work;
    }
    return result;
}


label triangleSearch::nodeCount() const
{
    return treePtr_.valid() ? treePtr_().nodes().size() : boundsNodes_.size();
}


label triangleSearch::referenceCount() const
{
    label count = triangleOrder_.size();
    if (treePtr_.valid())
    {
        for (const auto& node : treePtr_().nodes())
        {
            for (const labelBits index : node.subNodes_)
            {
                count += treePtr_().getIndices(index).size();
            }
        }
    }
    return count;
}


std::size_t triangleSearch::storageBytes() const
{
    std::size_t bytes =
        trianglePoints_.size()*sizeof(point)
      + triangleBounds_.size()*sizeof(boundBox)
      + triangleOrder_.size()*sizeof(label)
      + boundsNodes_.capacity()*sizeof(boundsNode);
    if (treePtr_.valid())
    {
        bytes += treePtr_().nodes().size()
                *sizeof(indexedOctree<triangleData>::node);
        for (const auto& node : treePtr_().nodes())
        {
            for (const labelBits index : node.subNodes_)
            {
                if (indexedOctree<triangleData>::isContent(index))
                {
                    bytes += sizeof(labelList);
                }
            }
        }
        bytes += referenceCount()*sizeof(label);
    }
    return bytes;
}


pointField triangleSearch::trianglesOverlapping
(
    const boundBox& targetBb
) const
{
    const labelList candidates(query(targetBb));
    pointField result(3*candidates.size());

    forAll(candidates, candidatei)
    {
        const label pointi = 3*candidates[candidatei];
        for (direction pointi2 = 0; pointi2 < 3; ++pointi2)
        {
            result[3*candidatei + pointi2] =
                trianglePoints_[pointi + pointi2];
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
    validateTolerance(tolerance);
    const triangleSearch search(trianglePoints, tolerance);
    return search.trianglesOverlapping(targetBb);
}


bool triangleSearch::intersectsCell
(
    const fvMesh& mesh,
    const label celli
) const
{
    const scalar tolerance
    (
        max(tolerance_.x(), max(tolerance_.y(), tolerance_.z()))
    );
    if (trianglePoints_.size() % 3)
    {
        FatalErrorInFunction
            << "Triangle point count " << trianglePoints_.size()
            << " is not divisible by three"
            << exit(FatalError);
    }

    if (trianglePoints_.empty())
    {
        return false;
    }

    const vector toleranceVec(tolerance, tolerance, tolerance);
    boundBox cellBb(mesh.cellBb(celli));
    cellBb.grow(toleranceVec);
    const labelList candidates(query(mesh.cellBb(celli)));
    if (candidates.empty())
    {
        return false;
    }

    forAll(candidates, candidatei)
    {
        const label pointi = 3*candidates[candidatei];
        const triPoints tri
        (
            trianglePoints_[pointi],
            trianglePoints_[pointi + 1],
            trianglePoints_[pointi + 2]
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

    const cell& c = mesh.cells()[celli];
    const faceList& faces = mesh.faces();
    DynamicList<triPoints> cellTriangles;
    forAll(c, i)
    {
        const label facei = c[i];
        const face& f = faces[facei];
        if (mesh.tetBasePtIs()[facei] < 0 || f.size() < 3)
        {
            FatalErrorInFunction
                << "Cell " << celli << " face " << facei
                << " has no valid polyMesh::CELL_TETS base point"
                << exit(FatalError);
        }

        for (label tetPointi = 1; tetPointi < f.size() - 1; ++tetPointi)
        {
            const triPointRef faceTri
            (
                tetIndices(celli, facei, tetPointi).faceTri(mesh)
            );
            cellTriangles.append
            (
                triPoints(faceTri.a(), faceTri.b(), faceTri.c())
            );
        }
    }

    forAll(candidates, candidatei)
    {
        const label pointi = 3*candidates[candidatei];
        const triPoints tri
        (
            trianglePoints_[pointi],
            trianglePoints_[pointi + 1],
            trianglePoints_[pointi + 2]
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


bool intersectsCell
(
    const pointField& trianglePoints,
    const fvMesh& mesh,
    const label celli,
    const scalar tolerance
)
{
    const triangleSearch search(trianglePoints, tolerance);
    return search.intersectsCell(mesh, celli);
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
