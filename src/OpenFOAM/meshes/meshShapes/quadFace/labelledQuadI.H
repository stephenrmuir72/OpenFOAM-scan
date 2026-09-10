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

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ListType>
inline void Foam::labelledQuad::assignList(const ListType& list)
{
    const label len = list.size();

    // checkSize
    if (len < quadFace::size() || len > quadFace::size()+1)
    {
         FatalErrorInFunction
            << "size " << len << " != ("
            << quadFace::size() << " or " << (quadFace::size()+1) << ")"
            << Foam::abort(Foam::FatalError);
    }

    auto iter = list.begin();

    for (label i = 0; i < quadFace::size(); ++i)
    {
        (*this)[i] = *iter;
        ++iter;
    }

    if (len > quadFace::size()) { index_ = *iter; }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::labelledQuad::labelledQuad()
:
    quadFace(),
    index_(-1)
{}


inline Foam::labelledQuad::labelledQuad
(
    const quadFace& f,
    label idx
)
:
    quadFace(f),
    index_(idx)
{}


inline Foam::labelledQuad::labelledQuad
(
    label p0,
    label p1,
    label p2,
    label p3,
    label idx
)
:
    quadFace(p0, p1, p2, p3),
    index_(idx)
{}


inline Foam::labelledQuad::labelledQuad(const labelUList& list)
:
    quadFace(),
    index_(0)
{
    assignList(list);
}


inline Foam::labelledQuad::labelledQuad(std::initializer_list<label> list)
:
    quadFace(),
    index_(0)
{
    assignList(list);
}


inline Foam::labelledQuad::labelledQuad
(
    const labelUList& list,
    const FixedList<label, 4>& indices,
    label idx
)
:
    quadFace(list, indices),
    index_(idx)
{}


inline Foam::labelledQuad::labelledQuad
(
    const labelUList& list,
    const labelledQuad& indices
)
:
    quadFace(list, indices),
    index_(indices.index())
{}


inline Foam::labelledQuad::labelledQuad(Istream& is)
{
    operator>>(is, *this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline int Foam::labelledQuad::triangles
(
    const UList<point>& points,
    label& triI,
    UList<labelledTri>& tris,
    const flipType flipping
) const
{
    return quadFace::triangles(points, triI, tris, this->index(), flipping);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Foam::Istream& Foam::operator>>(Istream& is, labelledQuad& f)
{
    if (is.format() == IOstreamOption::BINARY)
    {
        Detail::readContiguous<labelledQuad>
        (
            is,
            reinterpret_cast<char*>(&f),
            sizeof(labelledQuad)
        );
    }
    else
    {
        is.readBegin("labelledQuad");
        is >> static_cast<quadFace&>(f) >> f.index();
        is.readEnd("labelledQuad");
    }

    is.check(FUNCTION_NAME);
    return is;
}


inline Foam::Ostream& Foam::operator<<(Ostream& os, const labelledQuad& f)
{
    if (os.format() == IOstreamOption::BINARY)
    {
        os.write
        (
            reinterpret_cast<const char*>(&f),
            sizeof(labelledQuad)
        );
    }
    else
    {
        os  << token::BEGIN_LIST
            << static_cast<const quadFace&>(f) << token::SPACE << f.index()
            << token::END_LIST;
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
