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

#ifndef Foam_patchFaceIntersection_H
#define Foam_patchFaceIntersection_H

#include "boundBox.H"
#include "face.H"
#include "labelVector.H"
#include "PackedList.H"
#include "pointField.H"
#include "triangle.H"
#include "indexedOctree.H"
#include "treeBoundBox.H"
#include "DynamicList.H"
#include "HashSet.H"

#include <cstddef>

namespace Foam
{

class fvMesh;

namespace patchFaceIntersection
{

class triangleData;

//- Update-local triangle index. Axis-flat/near-flat sets use a median-split
//  bounds tree; ordinary 3-D sets use indexedOctree. Queries expand by the
//  maximum tolerance component and return unique candidates.
class triangleSearch
{
public:

    struct queryStats
    {
        label nodes = 0;
        label triangleTests = 0;
    };

private:

    struct boundsNode
    {
        boundBox bb;
        label begin;
        label end;
        label left = -1;
        label right = -1;
    };

    pointField trianglePoints_;
    autoPtr<indexedOctree<triangleData>> treePtr_;
    vector tolerance_;
    List<boundBox> triangleBounds_;
    labelList triangleOrder_;
    DynamicList<boundsNode> boundsNodes_;

    label buildBoundsTree(const label begin, const label end);

    void queryBoundsTree
    (
        const label nodei,
        const boundBox& bb,
        DynamicList<label>& candidates,
        queryStats& stats
    ) const;

    void queryOctree
    (
        const label nodei,
        const treeBoundBox& bb,
        labelHashSet& candidates,
        queryStats& stats
    ) const;

public:

    triangleSearch
    (
        const pointField& trianglePoints,
        const vector& tolerance
    );

    triangleSearch
    (
        const pointField& trianglePoints,
        const scalar tolerance
    );

    ~triangleSearch();

    labelList query
    (
        const boundBox& targetBb,
        queryStats* stats = nullptr
    ) const;

    label nodeCount() const;
    label referenceCount() const;

    //- Estimated owned storage, excluding allocator overhead
    std::size_t storageBytes() const;

    pointField trianglesOverlapping(const boundBox& targetBb) const;

    bool intersectsCell(const fvMesh& mesh, const label celli) const;

    label size() const noexcept
    {
        return trianglePoints_.size()/3;
    }
};

//- Exact triangle/box overlap, with an axis-wise geometric tolerance
bool triangleIntersectsBox
(
    const triPoints& tri,
    const boundBox& box,
    const vector& tolerance
);

//- Exact triangle/triangle overlap, with a geometric tolerance
bool trianglesIntersect
(
    const triPoints& triA,
    const triPoints& triB,
    const scalar tolerance
);

//- Conservative voxel index range for a closed query box
bool closedVoxelRange
(
    const boundBox& searchBb,
    const labelVector& nDivs,
    const boundBox& queryBb,
    labelVector& minIds,
    labelVector& maxIds
);

//- Fill voxels intersected by a triangulated face
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
);

//- Test if a triangulated face intersects a voxel set to value
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
);

//- Collect non-coupled, non-overset boundary faces as flattened triangle points
pointField patchTriangles(const fvMesh& mesh);

//- Keep only triangles whose tolerance-expanded bounds overlap the target box
pointField trianglesOverlapping
(
    const pointField& trianglePoints,
    const boundBox& targetBb,
    const vector& tolerance
);

//- Test if any flattened triangle intersects a mesh cell
bool intersectsCell
(
    const pointField& trianglePoints,
    const fvMesh& mesh,
    const label celli,
    const scalar tolerance
);

//- Calculate the voxel count without overflow and subject to a budget
bool checkedVoxelCount
(
    const labelVector& nDivs,
    const label maxVoxels,
    label& nVoxels
);

//- Double active divisions if the resulting grid remains within its budget
bool refineVoxelDivisions
(
    const labelVector& nDivs,
    const labelVector& geometricD,
    const label maxVoxels,
    labelVector& refinedDivs,
    label& refinedVoxelCount
);

} // End namespace patchFaceIntersection
} // End namespace Foam

#endif

// ************************************************************************* //
