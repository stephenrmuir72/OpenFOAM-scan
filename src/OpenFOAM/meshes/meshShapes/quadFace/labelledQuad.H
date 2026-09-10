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
    Foam::labelledQuad

Description
    A quadFace with additional (region) index.

SourceFiles
    labelledQuadI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_labelledQuad_H
#define Foam_labelledQuad_H

#include "quadFace.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class labelledQuad;
inline Istream& operator>>(Istream&, labelledQuad&);
inline Ostream& operator<<(Ostream&, const labelledQuad&);

/*---------------------------------------------------------------------------*\
                         Class labelledQuad Declaration
\*---------------------------------------------------------------------------*/

class labelledQuad
:
    public quadFace
{
    // Private Data

        //- The object index (region)
        label index_;


    // Private Member Functions

        //- Assign from 4 or 5 labels. The 5th label (optional) is the region.
        template<class ListType>
        inline void assignList(const ListType& list);


public:

    // Constructors

        //- Default construct, with invalid point labels and region/index (-1).
        inline labelledQuad();

        //- Construct from quadFace
        //- and optional region/index (0 if unspecified)
        inline labelledQuad(const quadFace& f, label idx = 0);

        //- Construct from four point labels
        //- and optional region/index (0 if unspecified)
        inline labelledQuad
        (
            label p0, label p1, label p2, label p3,
            label idx = 0
        );

        //- Construct from a list of 4 or 5 labels. Default region/index is 0.
        inline explicit labelledQuad(const labelUList& list);

        //- Construct from a list of 4 or 5 labels. Default region/index is 0.
        inline explicit labelledQuad(std::initializer_list<label>);

        //- Copy construct from a subset of point labels
        //- and optional region/index (0 if unspecified)
        inline labelledQuad
        (
            const labelUList& list,
            const FixedList<label, 4>& indices,
            label idx = 0
        );

        //- Copy construct from a subset of point labels
        //- with region index from input labelledQuad
        inline labelledQuad
        (
            const labelUList& list,
            const labelledQuad& indices
        );

        //- Construct from Istream
        inline explicit labelledQuad(Istream& is);


    // Member Functions

        //- Return the index (eg, the region)
        label index() const noexcept { return index_; }

        //- Non-const access to the index (eg, the region)
        label& index() noexcept { return index_; }

        //- Set the index (eg, the region)
        void setIndex(label i) noexcept { index_ = i; }

        //- The region index
        label region() const noexcept { return index_; }

        //- Non-const access to the region index
        label& region() noexcept { return index_; }

        //- Set the region/index
        void region(label i) noexcept { index_ = i; }


    // Splitting

        //- Inherit all triangles() methods
        using quadFace::triangles;

        //- Split into triangles using given points.
        inline int triangles
        (
            const UList<point>& points,
            label& triI,
            UList<labelledTri>& triFaces,
            //! Virtual flipping before the split?
            const flipType flipping = flipType::NONE
        ) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, labelledQuad&);
        friend Ostream& operator<<(Ostream&, const labelledQuad&);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for labelledQuad
template<> struct is_contiguous<labelledQuad> : std::true_type {};

//- Contiguous label data for labelledQuad
template<> struct is_contiguous_label<labelledQuad> : std::true_type {};


//- Specialization to offset faces, used in ListListOps::combineOffset
template<>
struct offsetOp<labelledQuad>
{
    labelledQuad operator()(const labelledQuad& x, const label offset) const
    {
        labelledQuad f(x);
        f += offset;
        return f;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "labelledQuadI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
