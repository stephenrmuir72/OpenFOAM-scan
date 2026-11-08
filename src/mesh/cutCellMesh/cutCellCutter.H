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


//- Why a cut was refused. The sites are different defects: one is a face the
//- surface crosses in a closed loop, whose far side cannot be a polyMesh face,
//- and the other is a piece whose faces cannot be oriented consistently. They
//- are counted apart because a fix for one is not a fix for another, and which
//- one dominates decides what to build.
enum class cutRefusal
{
    none,                  //- not refused
    faceWithHole,          //- a face is crossed by a closed loop, so the far
                           //- side of that face would be an annulus, which a
                           //- polyMesh face cannot be
    inconsistentPieceFaces //- the faces of a piece cannot be oriented
                           //- consistently
};


//- One piece of a cut cell
struct cutCellPiece
{
    //- Geometry, in the coordinates of the cell that was cut
    cutCellPolyhedron poly;
};


//- A face shared by two pieces of a cut cell
struct cutCellPieceFace
{
    label pieceA = -1;
    label faceInA = -1;
    label pieceB = -1;
    label faceInB = -1;

    //- The surface that produced this face
    label surfacei = -1;

    //- The triangle whose plane produced this face. Internal: it is the
    //- provenance the four-argument form can record, and the five-argument
    //- form resolves it to a region.
    label tri = -1;

    //- The surface region this face belongs to, inherited from that triangle.
    //- -1 until a caller supplies the regions, so a caller that has no surface
    //- is given no claim rather than a plausible-looking index.
    label patch = -1;
};


//- Result of cutting one cell
struct cutCellResult
{
    //- Why the cut was refused, when it was. Set at the site that refused it.
    cutRefusal refusal = cutRefusal::none;

    //- Distinct supporting planes of the surface inside the cell, outward
    //- oriented and in a fixed order
    List<plane> planes;

    //- The triangle each plane came from, parallel to planes. The grouping
    //- merges facets that share a plane, so without this the provenance of a
    //- plane is lost and a face cannot say which surface region it belongs to.
    labelList planeTri;

    //- The pieces, ordered negative first against the planes
    List<cutCellPiece> pieces;

    //- The faces shared by two pieces. These are the interfaces where the two
    //- pieces are in different regions, and the internal faces the cell is
    //- divided by where they are not.
    List<cutCellPieceFace> sharedFaces;

    //- The regions of each of the cell's faces, and the piece each belongs to.
    //- The assembly needs these rather than a point query: the regions tile the
    //- face by construction, so which piece a sub-face is in is a lookup in a
    //- set we built, not a question asked of a polyhedron.
    List<List<pointField>> faceRegions;
    List<List<label>> faceRegionPiece;

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


//- The parts of the convex polygons qs that lie outside the convex polygon r,
//  as disjoint convex polygons whose union is that outside part. All polygons
//  are in the plane of normal n, with their vertices in order about n.
List<pointField> outsideConvex
(
    const List<pointField>& qs,
    const pointField& r,
    const vector& n
);


//- Cut a cell by the surface facets by way of the closed loop where the surface
//- meets the cell's boundary, rather than by arranging the facet planes. The
//- facets are still clipped to the cell and grouped by plane - that is what an
//- interface face is - but the pieces come from the loop: the cell's boundary is
//- split along it, the regions are stitched into pieces by shared edges, and
//- every piece sees the interface faces. Cost is linear in the facets crossing
//- the cell, where the arrangement of their planes is cubic and its regions are
//- thin wedges for a doubly-curved surface, which is where it refuses.
//- Turn the loop's per-cell-edge report on or off. It prints, for each edge of
//- the cell, the sub-edges each of its two faces contributes and whether they
//- agree, which is where regions fail to stitch into pieces.
void cutCellLoopDebug(const bool on);


//- True when the mesher is dumping this one cell (dumpCell in the dictionary).
//  It lifts the print budget of the part and nearest-region-edge probes for that
//  cell alone, so that a part the enumeration drops can be read whole without
//  turning on cutCellLoopDebug for the entire run. Added 26-Sep-2026 to dump one
//  part out of IFACE's `0` bucket, per 06 B.3 item 0.
cutState cutCellByLoop
(
    const cutCellPolyhedron& cell,
    const List<pointField>& tris,
    const labelList& patchOf,
    const scalar tol,
    cutCellResult& result,
    const bool dump = false
);


//- The part of a triangle inside a cell, empty if there is none. The triangle
//  is clipped against the half-space of each cell face, so the polygon follows
//  the actual faces of the cell rather than its bounding box. The cell must be
//  convex for this to be exact; for a mildly non-convex snapped cell it is a
//  conservative superset.
pointField triangleCellPolygon
(
    const pointField& tri,
    const cutCellPolyhedron& cell,
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


//- True when two pieces have the same zone membership across all surfaces
bool samePieceZones(const labelList& a, const labelList& b);


//- Group the pieces of one cut cell by the faces they share, so that pieces
//  with the same per-surface zone membership become one group.
labelList groupPieces
(
    const label nPieces,
    const List<cutCellPieceFace>& sharedFaces,
    const List<labelList>& pieceZones
);


//- Single-zone compatibility form of groupPieces()
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
