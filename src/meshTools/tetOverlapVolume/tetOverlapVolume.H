/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2014 OpenFOAM Foundation
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
    Foam::tetOverlapVolume

Description
    Calculates the overlap volume of two cells using tetrahedral decomposition

SourceFiles
    tetOverlapVolume.C
    tetOverlapVolumeTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_tetOverlapVolume_H
#define Foam_tetOverlapVolume_H

#include "FixedList.H"
#include "labelList.H"
#include "treeBoundBox.H"
#include "Tuple2.H"
#include "tetrahedron.H"

namespace Foam
{

// Forward Declarations
class primitiveMesh;
class polyMesh;

/*---------------------------------------------------------------------------*\
                      Class tetOverlapVolume Declaration
\*---------------------------------------------------------------------------*/

class tetOverlapVolume
{
    // Private Classes

        //- tetPoints handling : sum resulting volumes
        class sumMomentOp
        {
        public:
            Tuple2<scalar, point> vol_;

            inline sumMomentOp()
            :
                vol_(0.0, Zero)
            {}

            inline void operator()(const tetPoints& tet)
            {
                const tetPointRef t(tet.tet());
                scalar tetVol = t.mag();
                vol_.first() += tetVol;
                vol_.second() += (tetVol*t.centre());
            }
        };

        //- tetPoints combining : check for overlap
        class hasOverlapOp
        {
        public:

            const scalar threshold_;
            tetPointRef::sumVolOp iop_;
            bool ok_;

            inline hasOverlapOp(const scalar threshold)
            :
                threshold_(threshold),
                iop_(),
                ok_(false)
            {}

            //- Overlap two positively oriented tets. Early exit is only
            //  valid for all-positive fans (see cellCellOverlapMinDecomp).
            inline bool operator()
            (
                const tetPoints& A,
                const tetPoints& B,
                const scalar sign
            )
            {
                if (sign > 0)
                {
                    // Same accumulation order as the unsigned engine
                    tetTetOverlap<tetPointRef::sumVolOp>(A, B, iop_);
                }
                else
                {
                    tetPointRef::sumVolOp piece;
                    tetTetOverlap<tetPointRef::sumVolOp>(A, B, piece);
                    iop_.vol_ -= piece.vol_;
                }
                ok_ = (iop_.vol_ > threshold_);
                return ok_;
            }
        };

        //- tetPoints combining : sum overlap volume
        class sumOverlapOp
        {
        public:

            tetPointRef::sumVolOp iop_;

            inline sumOverlapOp()
            :
                iop_()
            {}

            //- Overlap two positively oriented tets; accumulate with the
            //  product of the fan orientations
            inline bool operator()
            (
                const tetPoints& A,
                const tetPoints& B,
                const scalar sign
            )
            {
                if (sign > 0)
                {
                    // Same accumulation order as the unsigned engine
                    tetTetOverlap<tetPointRef::sumVolOp>(A, B, iop_);
                }
                else
                {
                    tetPointRef::sumVolOp piece;
                    tetTetOverlap<tetPointRef::sumVolOp>(A, B, piece);
                    iop_.vol_ -= piece.vol_;
                }
                return false;
            }
        };

        //- tetPoints combining : sum overlap volume
        class sumOverlapMomentOp
        {
        public:

            sumMomentOp iop_;

            inline sumOverlapMomentOp()
            :
                iop_()
            {}

            //- Overlap two positively oriented tets; accumulate with the
            //  product of the fan orientations
            inline bool operator()
            (
                const tetPoints& A,
                const tetPoints& B,
                const scalar sign
            )
            {
                if (sign > 0)
                {
                    tetTetOverlap<sumMomentOp>(A, B, iop_);
                }
                else
                {
                    sumMomentOp piece;
                    tetTetOverlap<sumMomentOp>(A, B, piece);
                    iop_.vol_.first() -= piece.vol_.first();
                    iop_.vol_.second() -= piece.vol_.second();
                }
                return false;
            }
        };


    // Private member functions

        //- Tet overlap calculation
        template<class tetPointsOp>
        static void tetTetOverlap
        (
            const tetPoints& tetA,
            const tetPoints& tetB,
            tetPointsOp& insideOp
        );

        //- Cell overlap calculation
        template<class tetsOp>
        static void cellCellOverlapMinDecomp
        (
            const primitiveMesh& meshA,
            const label cellAI,
            const primitiveMesh& meshB,
            const label cellBI,
            const treeBoundBox& cellBbB,
            tetsOp& combineTetsOp
        );

        //- True if every non-degenerate first-vertex fan tet of the cell
        //  is positively oriented about the cell centre
        static bool positiveFan(const primitiveMesh& mesh, const label celli);

        //- Return a const treeBoundBox
        static treeBoundBox pyrBb
        (
            const pointField& points,
            const face& f,
            const point& fc
        );


public:

    //- Runtime type information
    ClassName("tetOverlapVolume");


    // Constructors

        //- Default construct
        tetOverlapVolume() = default;


    // Public members

        //- Return a list of cells in meshA which overlaps with cellBI in
        //- meshB
        labelList overlappingCells
        (
            const polyMesh& meshA,
            const polyMesh& meshB,
            const label cellBI
        ) const;

        //- Return true if overlap volume is greater than threshold
        bool cellCellOverlapMinDecomp
        (
            const primitiveMesh& meshA,
            const label cellAI,
            const primitiveMesh& meshB,
            const label cellBI,
            const treeBoundBox& cellBbB,
            const scalar threshold = 0.0
        ) const;

        //- Signed first-vertex fan volume and relative roundoff allowance.
        //  Inverted fan tets (concave or warped faces) are accepted with a
        //  negative sign; the signed sum equals the volume enclosed by the
        //  fan-triangulated cell surface. Overlap volumes use the same
        //  signed decomposition, so both remain consistent.
        static Tuple2<scalar, scalar> cellVolumeMinDecomp
        (
            const primitiveMesh& mesh,
            const label celli
        );

        //- Calculates the overlap volume
        scalar cellCellOverlapVolumeMinDecomp
        (
            const primitiveMesh& meshA,
            const label cellAI,

            const primitiveMesh& meshB,
            const label cellBI,
            const treeBoundBox& cellBbB
        ) const;

        //- Calculates the overlap volume and moment
        Tuple2<scalar, point> cellCellOverlapMomentMinDecomp
        (
            const primitiveMesh& meshA,
            const label cellAI,

            const primitiveMesh& meshB,
            const label cellBI,
            const treeBoundBox& cellBbB
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "tetOverlapVolumeTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
