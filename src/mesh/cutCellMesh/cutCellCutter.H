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

    This is the first slice of the S1 cutting core (07 section 6.6). It handles
    cells that one plane of the surface passes through. Cells crossed by more
    than one plane, which happens where a surface edge or corner passes through
    a cell, and cells whose plane does not span the whole cross-section, are
    reported as refused rather than cut incorrectly.

    The result of a cut is two pieces and the interface polygon between them.
    The pieces are general polyhedra, since a background cell refined by hexRef8
    carries faces that are subdivided or have extra points inserted.

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
    cut,         //- the cell was split; the outputs are set
    refused      //- more than one plane crosses it, or the plane is partial
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


//- Cut a cell by the surface triangles that cross it.
//
//  'minus' and 'plus' are the pieces on either side of the plane, and
//  'interface' is the polygon between them. The interface appears as a face of
//  both pieces, oppositely oriented, as a coupled interface pair requires.
cutState cutCellByTriangles
(
    const cutCellPolyhedron& cell,
    const List<pointField>& tris,
    const scalar tol,
    cutCellPolyhedron& minus,
    cutCellPolyhedron& plus,
    pointField& interface
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
