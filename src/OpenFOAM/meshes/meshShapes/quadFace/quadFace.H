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
    Foam::quadFace

Description
    A quadrilateral face using a FixedList of labels corresponding to mesh
    vertices.

See also
    Foam::face, Foam::triFace

SourceFiles
    quadFaceI.H
    quadFace.cxx
    quadFace.txx
    quadFaceIntersection.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_quadFace_H
#define Foam_quadFace_H

#include "FixedList.H"
#include "edgeList.H"
#include "pointHit.H"
#include "intersection.H"
#include "pointField.H"
#include "triangle.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class face;
class quadFace;
class triFace;
class labelledTri;
template<class T, int SizeMin> class DynamicList;

typedef List<quadFace> quadFaceList;    //!< List of quadFace

/*---------------------------------------------------------------------------*\
                          Class quadFace Declaration
\*---------------------------------------------------------------------------*/

class quadFace
:
    public FixedList<label, 4>
{
public:

    // Public Data Types

        //- Face flipping
        enum class flipType { NONE, FLIP };


private:

    // Typedefs

        //- The data content type
        typedef FixedList<label, 4> tuple_type;


    // Private Member Functions

        //- Approximate normal (for directional calculations).
        //  Normal vector obtained from three points on the face.
        //  Used when only the basic direction of the face is needed.
        //  For speed, the vector is not normalized and its magnitude has
        //  no particular meaning.
        inline vector approxNormal(const UList<point>& points) const;

        //- Find index of largest internal angle
        label mostConcaveAngle(const UList<point>& points) const;

        //- Split quad into triangles at the flattest internal angle.
        //  \Return number of faces created (normally two).
        template<class TriFaceT>
        int split
        (
            //! Reference points for the face geometry
            const UList<point>& points,
            //! [in,out] Index counter for triangles
            label& triI,
            //! [out] Output of triangle results
            UList<TriFaceT>& triFaces,
            //! Virtual flipping before the split?
            const flipType flipping
        ) const;

public:

    // Generated Methods

        //- Copy construct
        quadFace(const quadFace&) = default;

        //- Move construct
        quadFace(quadFace&&) = default;

        //- Copy assignment
        quadFace& operator=(const quadFace&) = default;

        //- Move assignment
        quadFace& operator=(quadFace&&) = default;

        //- The front() accessor (from FixedList) has no purpose
        void front() = delete;

        //- The back() accessor (from FixedList) has no purpose
        void back() = delete;


    // Constructors

        //- Default construct, with invalid vertex labels (-1)
        inline quadFace();

        //- Construct from four vertex labels
        inline quadFace(label p0, label p1, label p2, label p3) noexcept;

        //- Construct from an initializer list of four vertex labels
        inline explicit quadFace(std::initializer_list<label> list);

        //- Copy construct from a list of four vertex labels
        inline explicit quadFace(const labelUList& list);

        //- Copy construct from a subset of vertex labels
        inline quadFace
        (
            const labelUList& list,
            const FixedList<label, 4>& indices
        );

        //- Construct from Istream
        inline quadFace(Istream& is);


    // Member Functions

    // Access

        //- The first vertex
        constexpr label a() const noexcept { return get<0>(); }

        //- The second vertex
        constexpr label b() const noexcept { return get<1>(); }

        //- The third vertex
        constexpr label c() const noexcept { return get<2>(); }

        //- The fourth vertex
        constexpr label d() const noexcept { return get<3>(); }

        //- The first vertex
        constexpr label& a() noexcept { return get<0>(); }

        //- The second vertex
        constexpr label& b() noexcept { return get<1>(); }

        //- The third vertex
        constexpr label& c() noexcept { return get<2>(); }

        //- The fourth vertex
        constexpr label& d() noexcept { return get<3>(); }

        //- Next vertex on face (circular order)
        label nextLabel(label i) const { return (*this)[fcIndex(i)]; }

        //- Previous vertex on face (circular order)
        label prevLabel(label i) const { return (*this)[rcIndex(i)]; }

        //- The vertex on face (naming as per nextLabel, prevLabel),
        //- identical to operator[]
        label thisLabel(label i) const { return (*this)[i]; }


    // Queries

        //- True if vertices are unique and non-negative.
        inline bool good() const noexcept;


    // Other

        //- 'Collapse' face by marking duplicate vertex labels.
        //  Duplicates vertex labels are marked with '-1'
        //  (the lower vertex is retained).
        //  Return the collapsed size.
        inline label collapse();

        //- Flip the face in-place.
        //  The starting vertex of the original and flipped face are identical.
        inline void flip();

        //- Return the points corresponding to this face
        inline pointField points(const UList<point>& pts) const;

        //- Return as a Foam::face
        inline Foam::face face() const;

        //- Return centre (centroid)
        point centre(const UList<point>& points) const;

        //- Average of face points: a quick estimate of the face centre
        inline point average(const UList<point>& points) const;

        //- Calculate average value at centroid of face
        template<class Type>
        Type average(const UList<point>& unused, const Field<Type>& fld) const;

        //- The area normal - with magnitude equal to area of face.
        //- Calculated as the sum of triangulations to the centre point
        vector areaNormal(const UList<point>& points) const;

        //- The unit normal
        inline vector unitNormal(const UList<point>& points) const;

        //- Magnitude of face area.
        //- Calculated as the sum of triangulations to the centre point
        inline scalar mag(const UList<point>& points) const;

        //- Magnitude squared of face area
        inline scalar magSqr(const UList<point>& points) const;

        //- The enclosing (bounding) box for the face
        inline Pair<point> box(const UList<point>& points) const;

        //- Number of triangles after splitting == 2
        static constexpr label nTriangles() noexcept { return 2; }

        //- Return face with reverse direction.
        //  The starting vertex of the original and reverse face are identical.
        inline quadFace reverseFace() const;

        //- Find local vertex on face for the vertex label, same as find()
        //  \return position in face (0,1,2,3) or -1 if not found.
        inline label which(const label vertex) const;

        //- Return swept-volume from old-points to new-points
        inline scalar sweptVol
        (
            const UList<point>& oldPoints,
            const UList<point>& newPoints
        ) const;

        //- Return the inertia tensor, with optional reference
        //- point and density specification
        inline tensor inertia
        (
            const UList<point>& points,
            const point& refPt = vector::zero,
            scalar density = 1.0
        ) const;

        //- Return point intersection with a ray starting at p, in direction q.
        pointHit ray
        (
            const point& p,
            const vector& q,
            const UList<point>& points,
            const intersection::algorithm = intersection::FULL_RAY,
            const intersection::direction dir = intersection::VECTOR
        ) const;

        //- Fast intersection with a ray.
        //  Does face-centre fan decomposition and returns triangle intersection
        //  point closest to p.
        pointHit intersection
        (
            const point& p,
            const vector& q,
            const UList<point>& points,
            const intersection::algorithm alg,
            const scalar tol = 0.0
        ) const;

        pointHit intersection
        (
            const point& p,
            const vector& q,
            const point& ctr,
            const UList<point>& points,
            const intersection::algorithm alg,
            const scalar tol = 0.0
        ) const;

        //- Return nearest point to face
        pointHit nearestPoint
        (
            const point& p,
            const UList<point>& points
        ) const;

        //- Return nearest point to face and classify it:
        //  + near point (nearType=POINT, nearLabel=0, 1, 2, 3)
        //  + near edge  (nearType=EDGE,  nearLabel=0, 1, 2, 3)
        pointHit nearestPointClassify
        (
            const point& p,
            const UList<point>& points,
            label& nearType,
            label& nearLabel
        ) const;

        //- The sign for which side of the face plane the point is on.
        //  Uses the supplied tolerance for rounding around zero.
        //  \return
        //  -  0: on plane
        //  - +1: front-side
        //  - -1: back-side
        inline int sign
        (
            const point& p,
            const UList<point>& points,
            const scalar tol = SMALL
        ) const;

        //- The number of edges == 4
        static constexpr label nEdges() noexcept { return 4; }

        //- Return i-th face edge (forward walk order).
        //  The edge 0 is from [0] to [1]
        //  and edge 3 is from [3] to [0]
        inline Foam::edge edge(label edgei) const;

        //- Return vector of i-th face edge (forward walk order).
        //  The edge 0 is from [0] to [1]
        //  and edge 3 is from [3] to [0]
        inline vector edge(label edgei, const UList<point>& pts) const;

        //- Return i-th face edge in reverse walk order.
        //  The rcEdge 0 is from [0] to [3]
        //  and rcEdge 3 is from [1] to [0]
        inline Foam::edge rcEdge(label edgei) const;

        //- Return vector of i-th face edge in reverse walk order.
        //  The rcEdge 0 is from [0] to [3]
        //  and rcEdge 3 is from [1] to [0]
        inline vector rcEdge(label edgei, const UList<point>& pts) const;

        //- Return list of edges in forward walk order.
        //  The edge 0 is from [0] to [1]
        //  and edge 3 is from [3] to [0]
        inline edgeList edges() const;

        //- Return list of edges in reverse walk order.
        //  The rcEdge 0 is from [0] to [3]
        //  and rcEdge 3 is from [1] to [0]
        inline edgeList rcEdges() const;


    // Searching

        //- Regular contains(pointi) tests
        using tuple_type::contains;

        //- True if face contains(edge)
        inline bool contains(const Foam::edge& e) const;

        //- Regular find pointi within face
        using tuple_type::find;

        //- Find the edge within the face.
        //  \return the edge index or -1 if not found
        inline label find(const Foam::edge& e) const;

        //- Test the edge direction on the face
        //  \return
        //  - +1: forward (counter-clockwise) on the face
        //  - -1: reverse (clockwise) on the face
        //  -  0: edge not found on the face
        inline int edgeDirection(const Foam::edge& e) const;

        //- Compare quadFaces
        //  \return:
        //  -  0: different
        //  - +1: identical
        //  - -1: same face, but different orientation
        static inline int compare(const quadFace& a, const quadFace& b);


    // Splitting

        //- Split into triangles using given points.
        //  Splits intelligently to maximize triangle quality.
        //  \Return number of faces created.
        int triangles
        (
            const UList<point>& points,
            //! [in,out] write position into tris
            label& triI,
            //! [out] The triangle results, written starting at triI
            UList<Foam::face>& tris,
            //! Virtual flipping before the split?
            const flipType flipping = flipType::NONE
        ) const;

        //- Split into triangles using given points.
        //  \Return number of faces created.
        template<int SizeMin>
        int triangles
        (
            const UList<point>& points,
            //! [out] The triangle results (appended to list)
            DynamicList<Foam::face, SizeMin>& tris,
            //! Virtual flipping before the split?
            const flipType flipping = flipType::NONE
        ) const;

        //- Split into triangles using given points.
        //  \Return number of faces created.
        int triangles
        (
            const UList<point>& points,
            //! [in,out] write position into tris
            label& triI,
            //! [out] The triangle results, written starting at triI
            UList<triFace>& tris,
            //! Virtual flipping before the split?
            const flipType flipping = flipType::NONE
        ) const;

        //- Split into triangles using given points.
        //  \Return number of faces created.
        template<int SizeMin>
        int triangles
        (
            const UList<point>& points,
            //! [out] The triangle results (appended to list)
            DynamicList<triFace, SizeMin>& tris,
            //! Virtual flipping before the split?
            const flipType flipping = flipType::NONE
        ) const;

        //- Split into triangles using given points.
        //  \Return number of faces created.
        int triangles
        (
            const UList<point>& points,
            //! [in,out] write position into tris
            label& triI,
            //! [out] The triangle results, written starting at triI
            UList<labelledTri>& tris,
            const label region = 0,
            //! Virtual flipping before the split?
            const flipType flipping = flipType::NONE
        ) const;

        //- Split into triangles using given points.
        //  \Return number of faces created.
        template<int SizeMin>
        int triangles
        (
            const UList<point>& points,
            //! [out] The triangle results (appended to list)
            DynamicList<labelledTri, SizeMin>& tris,
            const label region = 0,
            //! Virtual flipping before the split?
            const flipType flipping = flipType::NONE
        ) const;


    // Member Operators

        //- Inherit assignment
        using tuple_type::operator=;

        //- Increment (offset) vertices by given amount
        inline void operator+=(label vertexOffset);

        //- Equality: quadFaces are equal if they share the same vertices
        //- in the same or opposite order (orientation-agnostic)
        bool operator==(const quadFace& other) const noexcept
        {
            return (quadFace::compare(*this, other) != 0);
        }

        //- Inequality
        bool operator!=(const quadFace& other) const noexcept
        {
            return (quadFace::compare(*this, other) == 0);
        }


    // Hashing

        //- The (commutative) hash value for quadFace, hashes lower value first
        inline unsigned hash_code(unsigned seed=0) const noexcept
        {
            // Sort smallest-first
            label p = a(), q = b(), r = c(), s = d();
            if (p > q) std::swap(p, q);
            if (r > s) std::swap(r, s);
            if (p > r) std::swap(p, r);
            if (q > s) std::swap(q, s);
            if (q > r) std::swap(q, r);

            Foam::Hash<label> op;
            return op(s, op(r, op(q, op(p, seed))));
        }

        //- Hashing functor for quadFace (commutative)
        struct hasher
        {
            unsigned operator()(const quadFace& obj, unsigned seed=0) const
            {
                return obj.hash_code(seed);
            }
        };


    // Housekeeping

        //- Legacy name for areaNormal()
        //  \deprecated(2018-06) Deprecated for new use
        FOAM_DEPRECATED_FOR(2018-12, "areaNormal() or unitNormal()")
        vector normal(const UList<point>& points) const
        {
            return areaNormal(points);  // Legacy definition
        }
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for quadFace
template<> struct is_contiguous<quadFace> : std::true_type {};

//- Contiguous label data for quadFace
template<> struct is_contiguous_label<quadFace> : std::true_type {};

//- Hashing for quadFace uses commutative (incremental) hash
template<> struct Hash<quadFace> : quadFace::hasher {};


//- Specialization to offset faces, used in ListListOps::combineOffset
template<>
struct offsetOp<quadFace>
{
    quadFace operator()(const quadFace& x, const label offset) const
    {
        quadFace f(x);
        f += offset;
        return f;
    }
};


// * * * * * * * * * * * * * * * * * ADL Lookups  * * * * * * * * * * * * * //

template<std::size_t I>
constexpr label& get(quadFace& t) noexcept
{
    static_assert(I < 4, "Index out-of-bounds for quadFace");
    return Foam::get<I>(static_cast<FixedList<label, 4>&>(t));
}

template<std::size_t I>
constexpr const label& get(const quadFace& t) noexcept
{
    static_assert(I < 4, "Index out-of-bounds for quadFace");
    return Foam::get<I>(static_cast<const FixedList<label, 4>&>(t));
}

template<std::size_t I>
constexpr label&& get(quadFace&& t) noexcept
{
    static_assert(I < 4, "Index out-of-bounds for quadFace");
    return Foam::get<I>(static_cast<FixedList<label, 4>&&>(std::move(t)));
}

template<std::size_t I>
constexpr const label&& get(const quadFace&& t) noexcept
{
    static_assert(I < 4, "Index out-of-bounds for quadFace");
    return Foam::get<I>(static_cast<const FixedList<label, 4>&&>(std::move(t)));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specializations (structured binding)
namespace std
{

//- The tuple_size of quadFace == 4
template<>
struct tuple_size<Foam::quadFace> : std::integral_constant<std::size_t, 4> {};

//- Identical type (label) for all elements of quadFace
template<std::size_t I>
struct tuple_element<I, Foam::quadFace> { using type = Foam::label; };

} // End namespace std

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "quadFaceI.H"

#ifdef NoRepository
    #include "quadFace.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
