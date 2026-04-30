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
#include "face.H"
#include "triangle.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline int Foam::quadFace::compare(const quadFace& a, const quadFace& b)
{
    // Check all 4 forward rotations
    if
    (
        (a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3])
     || (a[0] == b[1] && a[1] == b[2] && a[2] == b[3] && a[3] == b[0])
     || (a[0] == b[2] && a[1] == b[3] && a[2] == b[0] && a[3] == b[1])
     || (a[0] == b[3] && a[1] == b[0] && a[2] == b[1] && a[3] == b[2])
    )
    {
        // Identical
        return 1;
    }
    else if
    (
        // Check all 4 reverse rotations (b reversed: [0],[3],[2],[1])
        (a[0] == b[0] && a[1] == b[3] && a[2] == b[2] && a[3] == b[1])
     || (a[0] == b[3] && a[1] == b[2] && a[2] == b[1] && a[3] == b[0])
     || (a[0] == b[2] && a[1] == b[1] && a[2] == b[0] && a[3] == b[3])
     || (a[0] == b[1] && a[1] == b[0] && a[2] == b[3] && a[3] == b[2])
    )
    {
        // Same face, but reversed orientation
        return -1;
    }

    return 0;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::vector
Foam::quadFace::approxNormal(const UList<point>& pts) const
{
    // Use first three vertices to define the face plane
    const auto& p0 = pts[a()];
    const auto& p1 = pts[b()];
    const auto& p2 = pts[c()];

    return (p1 - p0) ^ (p2 - p0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::quadFace::quadFace()
:
    FixedList<label, 4>(-1)
{}


inline Foam::quadFace::quadFace(label p0, label p1, label p2, label p3) noexcept
{
    a() = p0;
    b() = p1;
    c() = p2;
    d() = p3;
}


inline Foam::quadFace::quadFace(std::initializer_list<label> list)
:
    FixedList<label, 4>(list)
{}


inline Foam::quadFace::quadFace(const labelUList& list)
:
    FixedList<label, 4>(list)
{}


inline Foam::quadFace::quadFace
(
    const labelUList& list,
    const FixedList<label, 4>& indices
)
:
    FixedList<label, 4>(list, indices)
{}


inline Foam::quadFace::quadFace(Istream& is)
:
    FixedList<label, 4>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::quadFace::good() const noexcept
{
    return
    (
        a() >= 0 && a() != b() && a() != c() && a() != d()
     && b() >= 0 && b() != c() && b() != d()
     && c() >= 0 && c() != d()
     && d() >= 0
    );
}


inline Foam::label Foam::quadFace::collapse()
{
    // Cannot resize FixedList, so mark duplicates with '-1'
    // (the lower vertex is retained).
    // Catch any existing '-1' (eg, if called multiple times).

    label n = 4;

    // Forward consecutive duplicates
    for (label i = 1; i < 4; ++i)
    {
        if (operator[](i) == operator[](i-1) || operator[](i) == -1)
        {
            operator[](i) = -1;
            --n;
        }
    }

    // Wrap-around duplicate: [3] == [0]
    if (operator[](3) != -1 && operator[](3) == operator[](0))
    {
        operator[](3) = -1;
        --n;
    }

    return n;
}


inline void Foam::quadFace::flip()
{
    // Reverse in-place; starting vertex [0] is unchanged.
    // [0,1,2,3] -> [0,3,2,1]
    std::swap(get<1>(), get<3>());
}


inline Foam::pointField Foam::quadFace::points(const UList<point>& pts) const
{
    pointField p(4);

    p[0] = pts[a()];
    p[1] = pts[b()];
    p[2] = pts[c()];
    p[3] = pts[d()];

    return p;
}


inline Foam::face Foam::quadFace::face() const
{
    return Foam::face(*this);
}


inline Foam::point Foam::quadFace::average(const UList<point>& pts) const
{
    // Average of the four vertices.
    // Unlike face::average, there are few points so don't yet bother
    // with enforcing double precision
    return 0.25*(pts[a()] + pts[b()] + pts[c()] + pts[d()]);
}


inline Foam::vector
Foam::quadFace::unitNormal(const UList<point>& pts) const
{
    const vector n(this->areaNormal(pts));
    const scalar s(::Foam::mag(n));
    return s < ROOTVSMALL ? Zero : n/s;
}


inline Foam::scalar Foam::quadFace::mag(const UList<point>& pts) const
{
    return ::Foam::mag(areaNormal(pts));
}


inline Foam::scalar Foam::quadFace::magSqr(const UList<point>& pts) const
{
    return ::Foam::magSqr(areaNormal(pts));
}


inline Foam::Pair<Foam::point>
Foam::quadFace::box(const UList<point>& pts) const
{
    Pair<point> bb(pts[get<0>()], pts[get<0>()]);

    for (label i = 1; i < quadFace::size(); ++i)
    {
        bb.first()  = Foam::min(bb.first(),  pts[(*this)[i]]);
        bb.second() = Foam::max(bb.second(), pts[(*this)[i]]);
    }
    return bb;
}


inline Foam::quadFace Foam::quadFace::reverseFace() const
{
    // The starting vertex of the original and reverse face are identical.
    return quadFace(a(), d(), c(), b());
}


inline Foam::label Foam::quadFace::which(const label vertex) const
{
    return FixedList<label, 4>::find(vertex);
}


inline Foam::vector Foam::quadFace::edge
(
    label edgei,
    const UList<point>& pts
) const
{
    return vector(pts[nextLabel(edgei)] - pts[thisLabel(edgei)]);
}


inline Foam::edge Foam::quadFace::rcEdge(label edgei) const
{
    // Edge 0 (forward and reverse) always starts at [0]
    // for consistency with face flipping.
    // rcEdge maps: 0->0, 1->3, 2->2, 3->1
    const label pointi = edgei ? (quadFace::size() - edgei) : 0;
    return Foam::edge(thisLabel(pointi), prevLabel(pointi));
}


inline Foam::vector Foam::quadFace::rcEdge
(
    label edgei,
    const UList<point>& pts
) const
{
    const label pointi = edgei ? (quadFace::size() - edgei) : 0;
    return vector(pts[prevLabel(pointi)] - pts[thisLabel(pointi)]);
}


inline Foam::edgeList Foam::quadFace::edges() const
{
    edgeList theEdges(4);

    theEdges[0].first()  = a();
    theEdges[0].second() = b();

    theEdges[1].first()  = b();
    theEdges[1].second() = c();

    theEdges[2].first()  = c();
    theEdges[2].second() = d();

    theEdges[3].first()  = d();
    theEdges[3].second() = a();

    return theEdges;
}


inline Foam::edgeList Foam::quadFace::rcEdges() const
{
    edgeList theEdges(4);

    theEdges[0].first()  = a();
    theEdges[0].second() = d();

    theEdges[1].first()  = d();
    theEdges[1].second() = c();

    theEdges[2].first()  = c();
    theEdges[2].second() = b();

    theEdges[3].first()  = b();
    theEdges[3].second() = a();

    return theEdges;
}


inline int Foam::quadFace::edgeDirection(const Foam::edge& e) const
{
    for (label i = 0; i < quadFace::size(); ++i)
    {
        if (e.first() == thisLabel(i))
        {
            if (e.second() == nextLabel(i)) return  1;  // forward
            if (e.second() == prevLabel(i)) return -1;  // reverse
        }
    }

    return 0;  // Not found
}


inline Foam::label Foam::quadFace::find(const Foam::edge& e) const
{
    for (label i = 0; i < quadFace::size(); ++i)
    {
        if
        (
            (e.first() == thisLabel(i) && e.second() == nextLabel(i))
         || (e.first() == nextLabel(i) && e.second() == thisLabel(i))
        )
        {
            return i;
        }
    }

    return -1;  // Not found
}


inline bool Foam::quadFace::contains(const Foam::edge& e) const
{
    return (edgeDirection(e) != 0);
}


inline Foam::scalar Foam::quadFace::sweptVol
(
    const UList<point>& oldPoints,
    const UList<point>& newPoints
) const
{
    // Central fan decomposition into 4 triangles.
    // Note: centre point always comes first for best accuracy (cf. face.C).

    const point oldCenter(this->centre(oldPoints));
    const point newCenter(this->centre(newPoints));

    scalar sv = 0;

    for (label pointi = 0; pointi < quadFace::size(); ++pointi)
    {
        // Note: for best accuracy, centre point comes first
        const auto prev = prevLabel(pointi);
        const auto curr = thisLabel(pointi);

        sv += triPointRef
        (
            oldCenter, oldPoints[prev], oldPoints[curr]
        ).sweptVol
        (
            triPointRef(newCenter, newPoints[prev], newPoints[curr])
        );
    }

    return sv;
}


inline Foam::tensor Foam::quadFace::inertia
(
    const UList<point>& pts,
    const point& refPt,
    scalar density
) const
{
    // Central fan decomposition into 4 triangles
    const point ctr(this->centre(pts));

    return
    (
        triPointRef(pts[a()], pts[b()], ctr).inertia(refPt, density)
      + triPointRef(pts[b()], pts[c()], ctr).inertia(refPt, density)
      + triPointRef(pts[c()], pts[d()], ctr).inertia(refPt, density)
      + triPointRef(pts[d()], pts[a()], ctr).inertia(refPt, density)
    );
}


inline int Foam::quadFace::sign
(
    const point& p,
    const UList<point>& pts,
    const scalar tol
) const
{
    // Use first three vertices to define the face plane
    return triPointRef(pts, a(), b(), c()).sign(p, tol);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::quadFace::operator+=(label vertexOffset)
{
    if (vertexOffset)
    {
        a() += vertexOffset;
        b() += vertexOffset;
        c() += vertexOffset;
        d() += vertexOffset;
    }
}


// ************************************************************************* //
