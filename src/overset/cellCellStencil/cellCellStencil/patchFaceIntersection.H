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

namespace Foam
{

class fvMesh;

namespace patchFaceIntersection
{

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
