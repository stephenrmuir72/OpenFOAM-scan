/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2025 OpenCFD Ltd.
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

#include "face.H"
#include "triangle.H"
#include "triFace.H"
#include "quadFace.H"
#include "labelledTri.H"
#include "mathematicalConstants.H"
#include "Circulator.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::face::typeName = "face";


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Calculate unit edge vectors for a polygon
void calcEdgeVectors
(
    Foam::UList<Foam::vector>& edgeVecs,
    const Foam::UList<Foam::point>& points,
    const Foam::UList<Foam::label>& verts
)
{
    using namespace Foam;

    const label nEdges = std::min(edgeVecs.size(), verts.size());

    for (label i = 0; i < nEdges; ++i)
    {
        edgeVecs[i] = Foam::normalised
        (
            points[verts[(i == nEdges-1) ? 0 : i+1]]  // fcIndex
          - points[verts[i]]
        );
    }
}


// Find index of largest internal between face edges
// \returns the face point index of the most concave angle
Foam::label mostConcaveAngle
(
    //! Normalized edge vectors
    const Foam::UList<Foam::vector>& edges,
    //! Reference face normal
    const Foam::vector& n,
    //! [out] the max internal angle (radians)
    Foam::scalar& maxAngle
)
{
    using namespace Foam;

    label index = 0;
    maxAngle = -GREAT;

    const auto nEdges = edges.size();

    forAll(edges, i)
    {
        const vector& leftEdge = edges[(i ? i-1 : nEdges-1)];  // rcIndex
        const vector& rightEdge = edges[i];

        // NOTE: is -ve angle since left edge pointing in other direction
        scalar edgeCos = (leftEdge & rightEdge);
        scalar angle = Foam::acos(clamp(edgeCos, -1, 1));

        vector edgeNormal = (rightEdge ^ leftEdge);
        if ((edgeNormal & n) > 0)
        {
            // Concave angle.
            angle = constant::mathematical::pi + angle;
        }
        else
        {
            // Convex angle. Note '-' to take into account that rightEdge
            // and leftEdge are head-to-tail connected.
            angle = constant::mathematical::pi - angle;
        }

        if (maxAngle < angle)
        {
            maxAngle = angle;
            index = i;
        }
    }

    return index;
}

}  // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<Foam::face::splitMode Mode, class TriFaceT, class QuadFaceT>
Foam::label Foam::face::split
(
    [[maybe_unused]] const UList<point>& points,
    label& triI,
    label& quadI,
    [[maybe_unused]] UList<TriFaceT>& triFaces,
    [[maybe_unused]] UList<QuadFaceT>& quadFaces,
    [[maybe_unused]] const flipType flipping
) const
{
    // QuadFaceT is allowed to be an integral (non-face type), which signals
    // that tris/quads go into a single output. For this case,
    // TriFaceT must be resizable (ie, not derived from FixedList)

    constexpr bool singleOutput = std::is_integral_v<QuadFaceT>;

    static_assert
    (
        !(
            singleOutput &&
            (Mode == splitMode::SPLIT_QUAD) &&
            std::is_base_of_v<Foam::triFace, TriFaceT>
        ),
        "Single-output required a resizable tri-face type!"
    );

    const label oldIndices = (triI + quadI);

    if (size() < 3)
    {
        FatalErrorInFunction
            << "Serious problem: asked to split a face with < 3 vertices"
            << abort(FatalError);
    }

    // Counting triangles is trivial and independent of any internal splitting
    if constexpr (Mode == splitMode::COUNT_TRIANGLE)
    {
        triI += nTriangles();
        return nTriangles();  // Early return
    }
    else if constexpr (Mode == splitMode::SPLIT_FAN)
    {
        // Fan split around vertex 0 (trivial)
        if (flipping == flipType::FLIP)
        {
            for (label i = size()-1; i >= 2; --i)
            {
                auto& f = triFaces[triI++];
                f.resize(3);  // (no-op for FixedList)
                f[0] = (*this)[0];
                f[1] = (*this)[i];
                f[2] = (*this)[i-1];
            }
        }
        else
        {
            for (label i = 2; i < size(); ++i)
            {
                auto& f = triFaces[triI++];
                f.resize(3);  // (no-op for FixedList)
                f[0] = (*this)[0];
                f[1] = (*this)[i-1];
                f[2] = (*this)[i];
            }
        }

        return nTriangles();  // Early return
    }
    else if (size() == 3)
    {
        // Triangle. Just copy.

        if constexpr
        (
            Mode == splitMode::COUNT_TRIANGLE  // already handled above
         || Mode == splitMode::COUNT_QUAD
        )
        {
            ++triI;
        }
        else  // All split modes
        {
            auto& f = triFaces[triI++];
            f.resize(3);  // (no-op for FixedList)
            f[0] = (*this)[0];
            f[1] = (*this)[1];
            f[2] = (*this)[2];

            if (flipping == flipType::FLIP)
            {
                // Handle flipped triangle
                std::swap(f[1], f[2]);
            }
        }
    }
    else if (size() == 4)
    {
        if constexpr (Mode == splitMode::COUNT_TRIANGLE)
        {
            // already handled above
            triI += 2;
        }
        else if constexpr (Mode == splitMode::COUNT_QUAD)
        {
            ++quadI;
        }
        else if constexpr (Mode == splitMode::SPLIT_TRIANGLE)
        {
            // Start at point with largest internal angle.
            const vector fNormal(approxNormal(points));
            FixedList<vector, 4> workEdges;
            UList<vector> edges(workEdges.data(), workEdges.size());
            calcEdgeVectors(edges, points, *this);

            scalar minAngle(0);
            const auto idx0 = mostConcaveAngle(edges, fNormal, minAngle);
            const auto idx1 = fcIndex(idx0);
            const auto idx2 = fcIndex(idx1);
            const auto idx3 = fcIndex(idx2);

            // Triangle 1
            auto& f1 = triFaces[triI++];
            {
                f1.resize(3);  // (no-op for FixedList)
                f1[0] = (*this)[idx0];
                f1[1] = (*this)[idx1];  // <- middle vertex
                f1[2] = (*this)[idx2];
            }

            // Triangle 2
            auto& f2 = triFaces[triI++];
            {
                f2.resize(3);  // (no-op for FixedList)
                f2[0] = (*this)[idx2];
                f2[1] = (*this)[idx3];  // <- middle vertex
                f2[2] = (*this)[idx0];
            }

            if (flipping == flipType::FLIP)
            {
                // Handle flipped quad face (see note in quadFace)
                std::swap(f1[1], f2[1]);  // <- middle vertices
            }
        }
        else if constexpr (Mode == splitMode::SPLIT_QUAD)
        {
            if constexpr (singleOutput)
            {
                // Write quads into the same output as tris
                auto& f = triFaces[triI++];
                f.resize(4);  // (no-op for FixedList)
                f[0] = (*this)[0];
                f[1] = (*this)[1];
                f[2] = (*this)[2];
                f[3] = (*this)[3];

                if (flipping == flipType::FLIP)
                {
                    // With flipped face orientation
                    std::swap(f[1], f[3]);
                }
            }
            else
            {
                auto& f = quadFaces[quadI++];
                f.resize(4);  // (no-op for FixedList)
                f[0] = (*this)[0];
                f[1] = (*this)[1];
                f[2] = (*this)[2];
                f[3] = (*this)[3];

                if (flipping == flipType::FLIP)
                {
                    // With flipped face orientation
                    std::swap(f[1], f[3]);
                }
            }
        }
    }
    else
    {
        // General polygon.
        // Use iterative approach with work stack to avoid temporary face
        // allocations and recursion.
        //
        // Stack stores explicit vertices:
        // [v0, v1, v2, ..., N_1, v0, v1, v2, ..., N_2, ...]
        // where N_1, N_2 are the counts for the first, second sub-face etc

        vector fNormal(approxNormal(points));
        List<vector> workEdges(size());
        List<label> workVerts(size());
        DynamicList<label> vertexStack;

        if (flipping == flipType::FLIP)
        {
            // Initialize stack with face vertices [flip order]
            fNormal *= -1;
            vertexStack.push_back((*this)[0]);
            for (label i = size()-1; i > 0; --i)
            {
                vertexStack.push_back((*this)[i]);
            }
        }
        else
        {
            // Initialize stack with face vertices [forward order]
            vertexStack.push_back(*this);
        }
        // ... and push the vertex count
        vertexStack.push_back(size());

        while (!vertexStack.empty())
        {
            // Pop count
            const label nVerts = vertexStack.back();
            vertexStack.pop_back();

            // The vertices (and edges) for the sub-face
            auto edges = workEdges.slice(0, nVerts);
            auto verts = workVerts.slice(0, nVerts);

            verts = vertexStack.slice(vertexStack.size()-nVerts);
            vertexStack.pop_back(nVerts);

            if (nVerts == 3)
            {
                // Triangle. Just copy.

                if constexpr
                (
                    Mode == splitMode::COUNT_TRIANGLE  // already handled above
                 || Mode == splitMode::COUNT_QUAD
                )
                {
                    ++triI;
                }
                else  // (SPLIT_TRIANGLE | SPLIT_QUAD)
                {
                    auto& f = triFaces[triI++];
                    f.resize(3);  // (no-op for FixedList)
                    f[0] = verts[0];
                    f[1] = verts[1];
                    f[2] = verts[2];
                }
            }
            else if (nVerts == 4)
            {
                if constexpr (Mode == splitMode::COUNT_TRIANGLE)
                {
                    // Note: already handled above
                    triI += 2;
                }
                else if constexpr (Mode == splitMode::COUNT_QUAD)
                {
                    ++quadI;
                }
                else if constexpr (Mode == splitMode::SPLIT_TRIANGLE)
                {
                    // Split quad directly into triangles

                    // Calculate edges for sub-face
                    calcEdgeVectors(edges, points, verts);

                    scalar angle(0);
                    const auto idx0 = mostConcaveAngle(edges, fNormal, angle);
                    const auto idx1 = ((idx0+1) % nVerts);
                    const auto idx2 = ((idx1+1) % nVerts);
                    const auto idx3 = ((idx2+1) % nVerts);

                    // Triangle 1
                    auto& f1 = triFaces[triI++];
                    {
                        f1.resize(3);  // (no-op for FixedList)
                        f1[0] = verts[idx0];
                        f1[1] = verts[idx1];
                        f1[2] = verts[idx2];
                    }

                    // Triangle 2
                    auto& f2 = triFaces[triI++];
                    {
                        f2.resize(3);  // (no-op for FixedList)
                        f2[0] = verts[idx2];
                        f2[1] = verts[idx3];
                        f2[2] = verts[idx0];
                    }
                }
                else if constexpr (Mode == splitMode::SPLIT_QUAD)
                {
                    if constexpr (singleOutput)
                    {
                        // Write quads into the same output as tris
                        auto& f = triFaces[triI++];
                        f.resize(4);
                        f[0] = verts[0];
                        f[1] = verts[1];
                        f[2] = verts[2];
                        f[3] = verts[3];
                    }
                    else
                    {
                        auto& f = quadFaces[quadI++];
                        f.resize(4);  // (no-op for FixedList)
                        f[0] = verts[0];
                        f[1] = verts[1];
                        f[2] = verts[2];
                        f[3] = verts[3];
                    }
                }
            }
            else
            {
                // Polygon: find split point for this sub-face
                // (including the new diagonal edges)

                // Calculate edges for sub-face
                calcEdgeVectors(edges, points, verts);

                // Find most concave angle
                scalar angle(0);
                const label startIdx = mostConcaveAngle(edges, fNormal, angle);
                const scalar bisectAngle = (angle / 2);

                // Edge at startIdx for angle comparisons
                const vector& startEdge = edges[startIdx];

                // Find best split point (at least two points away from start)
                label minIndex = (startIdx + 2) % nVerts;
                scalar minDiff = constant::mathematical::pi;

                for (label i = 0; i < nVerts - 3; ++i)
                {
                    const label idx = (startIdx + 2 + i) % nVerts;
                    vector splitEdge = Foam::normalised
                    (
                        points[verts[idx]]
                      - points[verts[startIdx]]
                    );

                    const scalar splitCos = splitEdge & startEdge;
                    const scalar splitAngle = acos(clamp(splitCos, -1, 1));
                    const scalar angleDiff = fabs(splitAngle - bisectAngle);

                    if (minDiff > angleDiff)
                    {
                        minDiff = angleDiff;
                        minIndex = idx;
                    }
                }

                // Split into two subshapes.
                //     face1: startIdx to minIndex
                //     face2: minIndex to startIdx

                // Get sizes of the two subshapes
                label diff = 0;
                if (minIndex > startIdx)
                {
                    diff = (minIndex - startIdx);
                }
                else
                {
                    // Folded around
                    diff = (minIndex - startIdx + nVerts);
                }

                const label nPoints1 = (diff + 1);
                const label nPoints2 = (nVerts - diff + 1);

                // Push parts in reverse order (LIFO, so pop_back first)
                // - only push if the sub-face has at least 3 points

                if (nPoints2 >= 3)
                {
                    // Sub-face 2: from minIdx to startIdx
                    for (label i = minIndex; i != startIdx; i = (i+1) % nVerts)
                    {
                        vertexStack.push_back(verts[i]);
                    }
                    vertexStack.push_back(verts[startIdx]);
                    vertexStack.push_back(nPoints2);
                }

                if (nPoints1 >= 3)
                {
                    // Sub-face 1: from startIdx to minIndex
                    for (label i = startIdx; i != minIndex; i = (i+1) % nVerts)
                    {
                        vertexStack.push_back(verts[i]);
                    }
                    vertexStack.push_back(verts[minIndex]);
                    vertexStack.push_back(nPoints1);
                }
            }
        }
    }

    return (triI + quadI - oldIndices);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::face::face(const triFace& f)
:
    labelList(f)
{}


Foam::face::face(const quadFace& f)
:
    labelList(f)
{}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

int Foam::face::compare(const face& a, const face& b)
{
    // Basic rule: we assume that the sequence of labels in each list
    // will be circular in the same order (but not necessarily in the
    // same direction or from the same starting point).

    const label sizeA = a.size();
    const label sizeB = b.size();

    if (sizeA != sizeB)
    {
        // Trivial reject: faces have different sizes
        return 0;
    }
    else if (sizeA == 0)
    {
        // Both faces with zero vertices. Always identical
        return 1;
    }
    else if (sizeA == 1)
    {
        // Both faces with a single vertex. Simple check
        return (a[0] == b[0] ? 1 : 0);
    }

    ConstCirculator<face> aCirc(a);
    ConstCirculator<face> bCirc(b);

    // Rotate face b until its element matches the starting element of face a.
    do
    {
        if (aCirc() == bCirc())
        {
            // Set bCirc fulcrum to its iterator and increment the iterators
            bCirc.setFulcrumToIterator();
            ++aCirc;
            ++bCirc;

            break;
        }
    } while (bCirc.circulate(CirculatorBase::CLOCKWISE));

    // If the circulator has stopped then faces a and b do not share a matching
    // point. Doesn't work on matching, single element face.
    if (!bCirc.circulate())
    {
        return 0;
    }

    // Look forwards around the faces for a match
    do
    {
        if (*aCirc != *bCirc)
        {
            break;
        }
    }
    while
    (
        aCirc.circulate(CirculatorBase::CLOCKWISE),
        bCirc.circulate(CirculatorBase::CLOCKWISE)
    );

    // If the circulator has stopped then faces a and b matched.
    if (!aCirc.circulate())
    {
        return 1;
    }
    else
    {
        // Reset the circulators back to their fulcrum
        aCirc.setIteratorToFulcrum();
        bCirc.setIteratorToFulcrum();
        ++aCirc;
        --bCirc;
    }

    // Look backwards around the faces for a match
    do
    {
        if (*aCirc != *bCirc)
        {
            break;
        }
    }
    while
    (
        aCirc.circulate(CirculatorBase::CLOCKWISE),
        bCirc.circulate(CirculatorBase::ANTICLOCKWISE)
    );

    // If the circulator has stopped then faces a and b matched.
    if (!aCirc.circulate())
    {
        return -1;
    }

    return 0;
}


bool Foam::face::sameVertices(const face& a, const face& b)
{
    const label sizeA = a.size();
    const label sizeB = b.size();

    // Trivial reject: faces are different size
    if (sizeA != sizeB)
    {
        return false;
    }
    // Trivial: face with a single vertex
    else if (sizeA == 1)
    {
        return (a[0] == b[0]);
    }

    for (const auto pointi : a)
    {
        // Number of occurrences within 'a'
        label diff = 0;
        for (auto val : a)
        {
            if (pointi == val) ++diff;
        }

        // Versus number of occurrences within 'b'
        for (auto val : b)
        {
            if (pointi == val) --diff;
        }

        // Same number of occurrences in 'a' and 'b'?
        if (diff != 0) return false;
    }

    return true;
}


unsigned Foam::face::symmhash_code(const UList<label>& f, unsigned seed)
{
    Foam::Hash<label> op;

    label len = f.size();

    if (!len)
    {
        // Trivial: zero-sized
        return seed;
    }
    else if (len == 1)
    {
        // Trivial: single vertex
        return op(f[0], seed);
    }

    // Find location of the min vertex
    label pivot = 0;
    for (label i = 1; i < len; ++i)
    {
        if (f[pivot] > f[i])
        {
            pivot = i;
        }
    }

    // Use next lowest value for deciding direction to circulate
    if (f.fcValue(pivot) < f.rcValue(pivot))
    {
        // Forward circulate
        while (len--)
        {
            seed = op(f[pivot], seed);
            pivot = f.fcIndex(pivot);
        }
    }
    else
    {
        // Reverse circulate
        while (len--)
        {
            seed = op(f[pivot], seed);
            pivot = f.rcIndex(pivot);
        }
    }

    return seed;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::face::collapse()
{
    if (size() > 1)
    {
        label ci = 0;
        for (label i=1; i<size(); i++)
        {
            if (operator[](ci) != operator[](i))
            {
                operator[](++ci) = operator[](i);
            }
        }

        if (operator[](ci) != operator[](0))
        {
            ++ci;
        }

        labelList::resize(ci);
    }

    return size();
}


void Foam::face::flip()
{
    const label n = size();

    if (n > 2)
    {
        for (label i=1; i < (n+1)/2; ++i)
        {
            std::swap(operator[](i), operator[](n-i));
        }
    }
}


Foam::point Foam::face::centre(const UList<point>& points) const
{
    // Calculate the centre by breaking the face into triangles and
    // area-weighted averaging their centres

    const label nPoints = size();

    // If the face is a triangle, do a direct calculation
    if (nPoints == 3)
    {
        return triPointRef::centre
        (
            points[(*this)[0]],
            points[(*this)[1]],
            points[(*this)[2]]
        );
    }
    // else if (nPoints == 0)
    // {
    //     // Should not happen
    //     return Zero;
    // }


    // Note: use double precision to avoid overflows when summing

    // Estimated centre by averaging the face points
    const point centrePoint(this->average(points));

    doubleScalar sumA(0);
    doubleVector sumAc(Zero);

    for (label pointi = 0; pointi < nPoints; ++pointi)
    {
        const auto& thisPoint = points[thisLabel(pointi)];
        const auto& nextPoint = points[nextLabel(pointi)];

        // Calculate 3*triangle centre
        const vector ttc(thisPoint + nextPoint + centrePoint);

        // Calculate 2*triangle area
        const scalar ta = Foam::mag
        (
            (thisPoint - centrePoint)
          ^ (nextPoint - centrePoint)
        );

        sumA += ta;
        sumAc += ta*ttc;
    }

    if (sumA > VSMALL)
    {
        return sumAc/(3.0*sumA);
    }
    else
    {
        return centrePoint;
    }
}


Foam::vector Foam::face::areaNormal(const UList<point>& points) const
{
    const label nPoints = size();

    // Calculate the area normal by summing the face triangle area normals.
    // Changed to deal with small concavity by using a central decomposition

    // If the face is a triangle, do a direct calculation
    if (nPoints == 3)
    {
        return triPointRef::areaNormal
        (
            points[(*this)[0]],
            points[(*this)[1]],
            points[(*this)[2]]
        );
    }
    // else if (nPoints == 0)
    // {
    //     // Should not happen
    //     return Zero;
    // }

    // Estimated centre by averaging the face points
    const point centrePoint(this->average(points));

    vector norm(Zero);

    for (label pointi = 0; pointi < nPoints; ++pointi)
    {
        const auto& thisPoint = points[thisLabel(pointi)];
        const auto& nextPoint = points[nextLabel(pointi)];

        // Calculate 2*triangle area normal
        norm +=
        (
            (thisPoint - centrePoint)
          ^ (nextPoint - centrePoint)
        );
    }

    return (0.5*norm);
}


Foam::face Foam::face::reverseFace() const
{
    // Reverse the label list and return
    // The starting points of the original and reverse face are identical.

    const labelUList& origFace = *this;
    const label len = origFace.size();

    face newFace(len);
    if (len)
    {
        newFace[0] = origFace[0];
        for (label i=1; i < len; ++i)
        {
            newFace[i] = origFace[len - i];
        }
    }

    return newFace;
}


Foam::scalar Foam::face::sweptVol
(
    const UList<point>& oldPoints,
    const UList<point>& newPoints
) const
{
    // Calculate the swept volume by breaking the face into triangles and
    // summing their swept volumes.
    // Changed to deal with small concavity by using a central decomposition

    const point centreOldPoint(this->centre(oldPoints));
    const point centreNewPoint(this->centre(newPoints));

    scalar sv = 0;
    const label nPoints = size();

    for (label pointi = 0; pointi < nPoints; ++pointi)
    {
        // Note: for best accuracy, centre point comes first
        const auto prev = prevLabel(pointi);
        const auto curr = thisLabel(pointi);

        sv += triPointRef
        (
            centreOldPoint, oldPoints[prev], oldPoints[curr]
        ).sweptVol
        (
            triPointRef(centreNewPoint, newPoints[prev], newPoints[curr])
        );
    }

    return sv;
}


Foam::tensor Foam::face::inertia
(
    const UList<point>& p,
    const point& refPt,
    scalar density
) const
{
    // If the face is a triangle, do a direct calculation
    if (size() == 3)
    {
        return triPointRef
        (
            p, (*this)[0], (*this)[1], (*this)[2]
        ).inertia(refPt, density);
    }

    const point ctr = centre(p);

    tensor J = Zero;

    forAll(*this, i)
    {
        J +=
            triPointRef
            (
                p[thisLabel(i)], p[nextLabel(i)], ctr
            ).inertia(refPt, density);
    }

    return J;
}


Foam::edgeList Foam::face::edges() const
{
    const labelList& verts = *this;
    const label nVerts = verts.size();

    edgeList theEdges(nVerts);

    // Last edge closes the polygon
    theEdges.back().first() = verts.back();
    theEdges.back().second() = verts.front();

    for (label verti = 0; verti < nVerts - 1; ++verti)
    {
        theEdges[verti].first() = verts[verti];
        theEdges[verti].second() = verts[verti + 1];
    }

    return theEdges;
}


Foam::edgeList Foam::face::rcEdges() const
{
    const labelList& verts = *this;
    const label nVerts = verts.size();

    edgeList theEdges(nVerts);

    // First edge closes the polygon
    theEdges.front().first() = verts.front();
    theEdges.front().second() = verts.back();

    for (label verti = 1; verti < nVerts; ++verti)
    {
        theEdges[verti].first() = verts[nVerts - verti];
        theEdges[verti].second() = verts[nVerts - verti - 1];
    }

    return theEdges;
}


Foam::label Foam::face::find(const Foam::edge& e) const
{
    // NB: the point index is simultaneously the edge index.
    // ie, face point 0 starts edge 0, point 1 starts edge 1, ...

    if (const label idx = labelList::find(e.first()); idx >= 0)
    {
        if (e.second() == nextLabel(idx))
        {
            // Current edge matches forward
            return idx;
        }
        if (e.second() == prevLabel(idx))
        {
            // Previous edge matches reverse
            return rcIndex(idx);
        }
    }

    return -1;  // Not found
}


int Foam::face::edgeDirection(const Foam::edge& e) const
{
    if (const label idx = labelList::find(e.first()); idx >= 0)
    {
        if (e.second() == nextLabel(idx))
        {
            // Current edge matches forward. Could encode: (idx+1)
            return 1;
        }
        if (e.second() == prevLabel(idx))
        {
            // Previous edge matches reverse. Could encode: -(rcIndex(idx)+1)
            return -1;
        }
    }

    return 0;  // Not found
}


Foam::label Foam::face::nTriangles(const UList<point>&) const
{
    return nTriangles();
}


Foam::label Foam::face::triangles
(
    const UList<point>& points,
    label& triI,
    UList<face>& tris,
    const flipType flipping
) const
{
    label quadI(0);    // placeholder
    UList<int> quads;  // placeholder

    return split<SPLIT_TRIANGLE>(points, triI, quadI, tris, quads, flipping);
}


Foam::label Foam::face::triangles
(
    const UList<point>& points,
    label& triI,
    UList<triFace>& tris,
    const flipType flipping
) const
{
    label quadI(0);    // placeholder
    UList<int> quads;  // placeholder

    return split<SPLIT_TRIANGLE>(points, triI, quadI, tris, quads, flipping);
}


Foam::label Foam::face::triangles
(
    const UList<point>& points,
    label& triI,
    UList<labelledTri>& tris,
    const label region,
    const flipType flipping
) const
{
    label quadI(0);    // placeholder
    UList<int> quads;  // placeholder

    const label start = triI;

    auto total =
        split<SPLIT_TRIANGLE>(points, triI, quadI, tris, quads, flipping);

    // Set region/index information
    for (label i = start; i < triI; ++i)
    {
        tris[i].setIndex(region);
    }

    return total;
}


Foam::label Foam::face::trianglesFan
(
    const UList<point>& points,
    label& triI,
    UList<face>& tris,
    const flipType flipping
) const
{
    label quadI(0);    // placeholder
    UList<int> quads;  // placeholder

    return split<SPLIT_FAN>(points, triI, quadI, tris, quads, flipping);
}


Foam::label Foam::face::trianglesFan
(
    const UList<point>& points,
    label& triI,
    UList<triFace>& tris,
    const flipType flipping
) const
{
    label quadI(0);    // placeholder
    UList<int> quads;  // placeholder

    return split<SPLIT_FAN>(points, triI, quadI, tris, quads, flipping);
}


Foam::label Foam::face::trianglesFan
(
    const UList<point>& points,
    label& triI,
    UList<labelledTri>& tris,
    const label region,
    const flipType flipping
) const
{
    label quadI(0);    // placeholder
    UList<int> quads;  // placeholder

    const label start = triI;

    auto total = split<SPLIT_FAN>(points, triI, quadI, tris, quads, flipping);

    // Set region/index information
    for (label i = start; i < triI; ++i)
    {
        tris[i].setIndex(region);
    }

    return total;
}


Foam::label Foam::face::nTrianglesQuads
(
    const UList<point>& points,
    label& triI,
    label& quadI,
    const flipType flipping
) const
{
    UList<triFace> tris;    // unused
    UList<quadFace> quads;  // unused

    return split<COUNT_QUAD>(points, triI, quadI, tris, quads, flipping);
}


Foam::label Foam::face::trianglesQuads
(
    const UList<point>& points,
    label& triI,
    label& quadI,
    UList<face>& tris,
    UList<face>& quads,
    const flipType flipping
) const
{
    return split<SPLIT_QUAD>(points, triI, quadI, tris, quads, flipping);
}


Foam::label Foam::face::trianglesQuads
(
    const UList<point>& points,
    label& faceI,
    UList<face>& splitFaces,
    const flipType flipping
) const
{
    label quadI(0);    // placeholder
    UList<int> quads;  // integral type == use single output

    return split<SPLIT_QUAD>(points, faceI, quadI, splitFaces, quads, flipping);
}


Foam::label Foam::face::trianglesQuads
(
    const UList<point>& points,
    label& triI,
    label& quadI,
    UList<triFace>& tris,
    UList<quadFace>& quads,
    const flipType flipping
) const
{
    return split<SPLIT_QUAD>(points, triI, quadI, tris, quads, flipping);
}


Foam::label Foam::face::longestEdge(const UList<point>& pts) const
{
    const labelList& verts = *this;
    const label nVerts = verts.size();

    // Last edge closes the polygon. Use it to initialize loop
    label longest = nVerts - 1;
    scalar longestLen = pts[verts.front()].distSqr(pts[verts.back()]);

    // Examine other edges
    for (label edgei = 0; edgei < nVerts - 1; ++edgei)
    {
        scalar edgeLen = pts[verts[edgei]].distSqr(pts[verts[edgei+1]]);

        if (longestLen < edgeLen)
        {
            longestLen = edgeLen;
            longest = edgei;
        }
    }

    return longest;
}


// ************************************************************************* //
