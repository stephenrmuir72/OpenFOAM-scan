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

    - A cell crossed by several facets is cut into the arrangement of their
      planes, one convex piece per combination of sides. A piece is not a
      region, since a facet can cover part of the cross-section only, so every
      piece is classified against the surface instead.
    - Positions where the surface crosses a mesh edge are made once, through a
      cache keyed by the mesh edge, so the cells sharing that edge refer to one
      point. Two cells sharing an edge may compute it from different planes and
      differ in the last bits, so the cache resolves them by position.
    - A crossed face is subdivided once, by the planes of every cell that meets
      it, so the cells on both sides see the same sub-faces. This is what keeps
      the mesh conformal when the two cells are cut by different planes.
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

        //- Cells whose cut was refused as degenerate and left whole. A grazing
        //- surface can make two pieces meet along a line rather than over an
        //- area, and a face between them would have no area, so the cell is
        //- left as one cell and reported. This is a geometric outcome rather
        //- than an unsupported case, so it is counted and does not fail the run.
        //- The most pieces any one cell was cut into
        label maxPieces = 0;

        //- Cells cut into more than two pieces
        label nManyPieces = 0;

        //- Their labels, to the tenth, so a failing cell can be named
        labelList manyPieces;

        label nDegenerate = 0;

        //- Of those refusals, how many were each site. See cutRefusal.
        label nRefusedMalformed = 0;
        label nRefusedDegenerateFace = 0;
        label nRefusedFaceWithHole = 0;

        //- Their indices, in mesh order
        labelList degenerate;

        //- Cells reported and left untouched
        label nUnsupported = 0;

        //- Their indices, in mesh order
        labelList unsupported;

        //- Why each of them was unsupported
        wordList unsupportedReason;

        //- Sub-cells the cut cells were replaced by, before merging
        label nPieces = 0;

        //- Cells the merged cut cells were replaced by
        label nGroups = 0;

        //- Sub-cells absorbed by the merge
        label nMerged = 0;

        //- Faces dropped because the two pieces they separated are one cell
        label nMergedFaces = 0;

        //- Cut cells by the number of planes crossing them, byPlanes[n] being
        //- the number cut by n planes
        labelList byPlanes;

        //- Cut cells whose pieces all fall in one region, which means the
        //- surface passes through without separating the cell
        label nSingleRegion = 0;

        //- Cells absorbed into a neighbour because they were too small
        label nAbsorbed = 0;

        //- Mesh points the collapse moved, and the largest move. A cut
        //- crossing within vertexSnapFraction of an edge endpoint is collapsed
        //- onto that endpoint by moving the point onto the plane of the facet,
        //- so the cut passes through it and no near-vertex sliver is made.
        label nSnapped = 0;
        scalar maxSnapDistance = 0;

        //- Collapses refused because the resolved move was larger than
        //- maxSnapMoveFraction of the shortest edge at the point
        label nSnapRejected = 0;

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

        //- How close, as a fraction of a cell edge's length, a cut crossing has
        //- to be to the edge's endpoint to be collapsed onto it. 05 section
        //- 7.1.2 gives 0.05 as the default and 0.10 as the limit. Zero disables
        //- the collapse.
        scalar vertexSnapFraction_ = 0.05;

        //- The largest move a collapse may make, as a fraction of the shortest
        //- edge at the point. The trigger bounds a single pull - the crossing it
        //- collapses onto lies within vertexSnapFraction of the edge's endpoint,
        //- so the move is at most that fraction of that edge - but the solve for
        //- a point pulled by several planes is unbounded: two nearly parallel
        //- planes meet along a line that runs away as their dihedral angle
        //- closes, and a coarsely tessellated curved surface is full of such
        //- pairs. Zero disables the guard.
        scalar maxSnapMoveFraction_ = 0.05;

        //- Which cut to use: "plane" arranges the facet planes, "loop" traces
        //- where the surface meets the cell's boundary and takes the pieces as
        //- the regions that separates. Both produce the same cutCellResult; the
        //- loop is linear in the facets a cell meets where the arrangement is
        //- cubic in them and its regions degenerate for a doubly-curved surface.
        word cutMethod_ = "plane";

        //- True when the loop is cutting, in which case a piece is identified
        //- by containment rather than by its sign against each plane
        bool loopCut_ = false;

        //- Dump this one cell and the triangles crossing it, for replay
        label dumpCell_ = -1;

        //- Displacement of each mesh point by the collapse, zero for the ones
        //- that did not move. Sized by snapPoints().
        pointField snapDelta_;

        //- Parent background cell of every cell of the cut mesh
        labelList backgroundCell_;

        //- Background cells every cell of the cut mesh was built from. Entry
        //- zero is backgroundCell_[celli]; a merged cell lists all of them.
        List<labelList> backgroundParents_;

        //- Cell volumes of the mesh before it was cut, by background cell
        scalarField backgroundVolumes_;

        //- A cell below this fraction of its parents' volume is a sliver and is
        //- absorbed into a neighbour, until the merged group reaches
        //- minVolumeRatio_ of it. Either being zero disables the merge.
        //- Off by default: removeFaces still refuses some configurations, see
        //- 07 section 4.1, and a mesh that aborts is worse than a sliver.
        scalar sliverVolumeRatio_ = 0.0;
        scalar minVolumeRatio_ = 0.30;

        //- Region cellZone of every background cell
        labelList cellRegionZone_;


    // Private Data Types

        //- One cell to cut, and the pieces it is split into
        struct cellCut
        {
            //- The regions of each cell face and their pieces, on the loop path
            List<List<pointField>> faceRegions;
            List<List<label>> faceRegionPiece;

            //- Background cell
            label celli;

            //- Distinct planes of the surface inside it, in the cutter's order
            List<plane> planes;

            //- The pieces, ordered negative first against those planes
            List<cutCellPiece> pieces;

            //- Faces shared by two pieces, which are the interfaces and the
            //- internal faces the cell is divided by
            List<cutCellPieceFace> sharedFaces;

            //- Mesh cell of each piece. The pieces of one region that touch
            //- each other share one cell, and the group of the first piece
            //- keeps the label of the background cell.
            labelList pieceCell;

            //- Group of each piece, from groupPieces()
            labelList pieceGroup;

            //- Region cellZone of each piece
            labelList pieceZone;
        };


    // Private Member Functions

        //- Move every mesh point that a cut crossing falls too close to onto
        //- the plane of the facet that crosses there, so that the cut passes
        //- through the point. Must run before findCutCells() and before the
        //- mesh geometry is otherwise used.
        void snapPoints();

        //- Absorb every cell below sliverVolumeRatio_ of its parents' volume
        //- into the neighbour the merge scores best, in one topology change
        void mergeSlivers();

        //- Find the cells the surface crosses, cut them into pieces by the
        //- planes of every facet that crosses them, and classify each piece
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

        //- Write the background cells each cell of the cut mesh came from
        void writeBackgroundParents() const;

        //- A cell below this fraction of its parents' volume is absorbed into
        //- a neighbour of the same region. Zero disables the merge.
        void sliverVolumeRatio(const scalar f) { sliverVolumeRatio_ = f; }

        //- The volume a merged group is grown to, as a fraction of its
        //- parents' volume. Zero disables the merge.
        void minVolumeRatio(const scalar f) { minVolumeRatio_ = f; }

        //- Set the collapse fraction, 0 to 0.10, 0 disabling it
        void vertexSnapFraction(const scalar f)
        {
            vertexSnapFraction_ = f;
        }

        //- Set the largest move a collapse may make, as a fraction of the
        //- shortest edge at the point. Zero disables the guard.
        void maxSnapMoveFraction(const scalar f)
        {
            maxSnapMoveFraction_ = f;
        }

        //- Set which cut to use, "plane" or "loop"
        void cutMethod(const word& m)
        {
            cutMethod_ = m;
            loopCut_ = (m == "loop");
        }

        //- Set the cell to dump for replay, or -1 for none
        void dumpCell(const label celli)
        {
            dumpCell_ = celli;
        }

        //- Displacement of each mesh point by the collapse, zero where it did
        //- not move. Only meaningful after cut().
        const pointField& snapDisplacement() const noexcept
        {
            return snapDelta_;
        }

        //- Write snapDisplacement next to the cut mesh
        void writeSnapDisplacement() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
