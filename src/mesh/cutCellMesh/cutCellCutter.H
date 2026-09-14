/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
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

Description
    Per-cell cutting of a Cartesian background cell by surface triangles.

    This is the S1 cutting core (07 section 6.6). Every distinct plane of the
    surface inside the cell splits it, in a fixed order, so a cell crossed by
    several facets, which happens where a surface edge or corner passes through
    it, comes out as the arrangement of those planes: one convex piece per
    combination of sides.

    A piece is not a region. A facet can cover only part of the cross-section,
    so two pieces can lie in the same region and one region can be made of
    several pieces. The driver classifies each piece against the surface, which
    is why it, and not the cutter, needs the surface.

    A background cell refined by hexRef8 carries faces that are subdivided or
    have extra points inserted, so the pieces are general convex polyhedra.

SourceFiles
    cutCellCutter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cutCellCutter_H
#define Foam_cutCellCutter_H

#include "cutCellPolyhedron.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Result of attempting to cut one cell
enum class cutState
{
    untouched,   //- no triangle crosses the cell
    cut,         //- the cell was split, and at least two pieces came out
    refused      //- triangles cross it but produced no usable pieces
};


//- One piece of a cut cell
struct cutCellPiece
{
    //- Geometry, in the coordinates of the cell that was cut
    cutCellPolyhedron poly;

    //- Side of each of the cell's planes, +1 or -1. The pieces are ordered by
    //- this vector, negative first, so a cell always cuts the same way and the
    //- first piece is the one inside a single plane.
    labelList side;
};


//- A face shared by two pieces of a cut cell
struct cutCellPieceFace
{
    label pieceA = -1;
    label faceInA = -1;
    label pieceB = -1;
    label faceInB = -1;
};


//- Result of cutting one cell
struct cutCellResult
{
    //- Distinct supporting planes of the surface inside the cell, outward
    //- oriented and in a fixed order
    List<plane> planes;

    //- The pieces, ordered negative first against the planes
    List<cutCellPiece> pieces;

    //- The faces shared by two pieces. These are the interfaces where the two
    //- pieces are in different regions, and the internal faces the cell is
    //- divided by where they are not.
    List<cutCellPieceFace> sharedFaces;

    //- Number of triangles that clipped to a polygon of at least three points,
    //- that is the number that actually cross the cell
    label nPolygons = 0;
};


//- Clip a planar polygon by a plane, keeping the part on the negative side.
//  The polygon must be in cyclic order and the result is empty if fewer than
//  three points survive.
pointField clipPolygon
(
    const pointField& poly,
    const plane& pln,
    const scalar tol
);


//- The part of a triangle inside an axis-aligned box, empty if there is none.
//  This is REF-1 section 2.2 A: the triangle is clipped against the six planes
//  of the cell box, which both tests the intersection and yields the per-cell
//  polygon.
pointField triangleBoxPolygon
(
    const pointField& tri,
    const boundBox& bb,
    const scalar tol
);


//- Area of a planar polygon
scalar polygonArea(const pointField& poly);

//- True if a planar polygon is degenerate rather than merely thin: its area is
//- no more than the tolerance times its longest edge, which is a width of the
//- order of the tolerance. Such a polygon cannot be a face of a mesh, so a
//- split that would make one is better refused than built.
bool polygonDegenerate(const pointField& poly, const scalar tol);


//- Group the pieces of one cut cell by the faces they share, so that the pieces
//  of one region which touch each other come out as one group and become one
//  cell. 'pieceZone' is the region of each piece.
//
//  Returns the group index of each piece. Groups are numbered from zero in the
//  order their first piece appears, so the group of the first piece is zero,
//  which is the group that keeps the label of the background cell.
//
//  Only pieces that share a face are joined. A cell has to be face connected, so
//  pieces meeting at an edge or a point stay separate, and only pieces in the
//  same region are joined, since joining across a region would destroy the
//  interface between them.
labelList groupPieces
(
    const label nPieces,
    const List<cutCellPieceFace>& sharedFaces,
    const labelList& pieceZone
);


//- Cut a cell by the surface triangles that cross it, into the arrangement of
//  their distinct planes inside the cell.
//
//  The plane normals are those of the surface polygons, so for an outward
//  oriented surface the positive side of a plane is the outside of the
//  surface. The pieces come out ordered by their side against the planes,
//  negative first, so the first piece of a cell cut by one plane is the piece
//  inside the surface.
//
//  'result' is left empty unless the state is cut.
cutState cutCellByTriangles
(
    const cutCellPolyhedron& cell,
    const List<pointField>& tris,
    const scalar tol,
    cutCellResult& result
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
