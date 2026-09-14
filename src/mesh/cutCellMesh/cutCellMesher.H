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

Class
    Foam::cutCellMesher

Description
    Cuts a mesh by a surface, replacing every cell the surface passes through
    with one sub-cell on each side, in a single topology change.

    This is the first slice of the S1 cutting core (07 section 6.6). It handles
    cells that one plane of the surface passes through. Cells crossed by more
    than one plane, which is what happens where a surface edge or corner passes
    through a cell, and cells whose plane does not span the whole cross-section,
    are reported with their index and left untouched rather than cut
    incorrectly.

    The result is a conformal polyMesh: a face the surface crosses is split
    once, and both of the cells it borders see the same two sub-faces, so there
    are no hanging nodes and no duplicated points. Three things keep it that
    way.

    - Positions where the surface crosses a mesh edge are made once, through a
      cache keyed by the mesh edge, so the cells sharing that edge all refer to
      one point. A cut cell computes its own plane, and two cells sharing an
      edge may differ in the last bits, so the cache resolves them by position.
    - A crossed face is split from its own point list rather than from the two
      cells' overlapping views of it, and the two sub-polygons are shared.
    - A cell's faces are implied by their owner and neighbour, so splitting a
      face updates both of its cells without either being touched explicitly.

    Each sub-cell carries the label of the background cell it came from, which
    is the state S4 needs to restore a cell the surface has left (07 section
    6.3). One sub-cell keeps the label of its background cell, so labels of
    untouched cells, and the parent map, are both affordable to maintain.

    Output is a cellZone per region and a faceZone per interface, oriented from
    the inside region to the outside, with no baffles (07 section 6.4).

SourceFiles
    cutCellMesher.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cutCellMesher_H
#define Foam_cutCellMesher_H

#include "cutCellCutter.H"
#include "cutCellSurface.H"
#include "polyMesh.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyTopoChange;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class cutCellMesher Declaration
\*---------------------------------------------------------------------------*/

class cutCellMesher
{
public:

    //- Report of one cutting pass
    struct report
    {
        //- Cells in the mesh that was cut
        label nCells = 0;

        //- Cells actually crossed by the surface
        label nCandidate = 0;

        //- Cells replaced by two sub-cells
        label nCut = 0;

        //- Cells reported and left untouched
        label nUnsupported = 0;

        //- Their indices, in mesh order
        labelList unsupported;

        //- Why each of them was unsupported
        wordList unsupportedReason;

        //- Points added to the mesh
        label nNewPoints = 0;

        //- Total area of the interface
        scalar interfaceArea = 0;
    };


private:

    // Private Data

        //- Mesh, cut in place by cut()
        polyMesh& mesh_;

        //- Surfaces to cut with
        const PtrList<cutCellSurface>& surfaces_;

        //- Geometric tolerance
        scalar tol_;

        //- Report of the last cut()
        report info_;

        //- Parent background cell of every cell of the cut mesh
        labelList backgroundCell_;

        //- Region cellZone of every background cell
        labelList cellRegionZone_;


    // Private Data Types

        //- One cell to cut, and the two regions it is split into
        struct cellCut
        {
            //- Background cell
            label celli;

            //- Supporting plane
            vector normal;
            point refPoint;

            //- Region on the side the normal points away from, and on the other
            label minusZone;
            label plusZone;

            //- Sub-cell added on the positive side. The negative side keeps
            //- the label of the background cell.
            label plusCelli;
        };


    // Private Member Functions

        //- Find the cells the surface crosses, with the single plane through
        //- each, and report the cells that cannot be cut in this slice
        void findCutCells(List<cellCut>& cuts);

        //- Polyhedron of one mesh cell, faces outward oriented
        static cutCellPolyhedron cellPolyhedron
        (
            const polyMesh& mesh,
            const label celli
        );

        //- Cut the mesh, in one topology change
        autoPtr<mapPolyMesh> assemble(const List<cellCut>& cuts);

        //- Zone ID of a named cellZone, creating it empty if it does not exist
        label cellZoneID(const word& name);

        //- Zone ID of a named faceZone, creating it empty if it does not exist
        label faceZoneID(const word& name);


public:

    // Constructors

        //- Construct from a mesh and the surfaces to cut it with
        cutCellMesher
        (
            polyMesh& mesh,
            const PtrList<cutCellSurface>& surfaces
        );


    // Member Functions

        //- Cut the mesh. Returns the topology map, which is empty when the
        //- surface crosses nothing.
        autoPtr<mapPolyMesh> cut();

        //- Report of the last cut()
        const report& info() const noexcept
        {
            return info_;
        }

        //- Parent background cell of every cell of the cut mesh. Only
        //- meaningful after cut().
        const labelList& backgroundCell() const noexcept
        {
            return backgroundCell_;
        }

        //- Write backgroundCell next to the cut mesh
        void writeBackgroundCell() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
