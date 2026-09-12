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
    Foam::cutCellPolyhedron

Description
    A closed convex polyhedron, given as points and outward-oriented faces, that
    can be split by a plane.

    Used by the cut-cell mesh generator to split a background cell that a
    surface passes through into one piece on each side of the surface. The
    pieces are general polyhedra: a background cell refined by hexRef8 is a
    convex box whose faces may be subdivided, a coarse cell carrying four
    coplanar faces where its refined neighbour carries one, and the split
    inherits that structure.

    Convexity is assumed. The interface polygon is ordered by angle about the
    plane normal, which is only correct for a convex cross-section.

    A split returns the interface polygon as well as the two pieces. That
    polygon appears as a face of both pieces with opposite orientation, which is
    what a coupled interface pair requires.

SourceFiles
    cutCellPolyhedron.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cutCellPolyhedron_H
#define Foam_cutCellPolyhedron_H

#include "pointField.H"
#include "faceList.H"
#include "plane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class cutCellPolyhedron Declaration
\*---------------------------------------------------------------------------*/

class cutCellPolyhedron
{
    // Private Data

        //- Points, in the numbering used by faces_
        pointField points_;

        //- Faces, outward oriented, as lists of local point labels
        faceList faces_;


public:

    // Constructors

        //- Construct null
        cutCellPolyhedron() = default;

        //- Construct from points and faces
        cutCellPolyhedron(const pointField& pts, const faceList& fcs);


    // Member Functions

        //- Points
        const pointField& points() const noexcept
        {
            return points_;
        }

        //- Faces
        const List<face>& faces() const noexcept
        {
            return faces_;
        }

        //- True if there are no points
        bool empty() const noexcept
        {
            return points_.empty();
        }

        //- Number of faces
        label nFaces() const noexcept
        {
            return faces_.size();
        }

        //- Signed volume by the divergence theorem, positive when the faces
        //  are outward oriented
        scalar volume() const;

        //- Volume centroid
        point centre() const;

        //- Number of edges not shared by exactly two faces
        label nBadEdges() const;

        //- True if every edge is shared by exactly two faces
        bool closed() const;

        //- Split by a plane.
        //
        //  Returns false if the plane does not properly divide the polyhedron,
        //  that is if it misses it, touches it, or contains one of its faces.
        //
        //  On success 'negative' holds the piece on the side the plane normal
        //  points away from and 'positive' the other piece. 'cutFace' holds the
        //  interface polygon, ordered counter-clockwise about the plane normal,
        //  and the same polygon appears in both pieces with opposite
        //  orientation.
        bool split
        (
            const plane& pln,
            cutCellPolyhedron& negative,
            cutCellPolyhedron& positive,
            pointField& cutFace
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
