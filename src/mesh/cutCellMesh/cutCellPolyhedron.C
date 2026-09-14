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

\*---------------------------------------------------------------------------*/

#include "cutCellPolyhedron.H"
#include "DynamicList.H"
#include "Map.H"

#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Scalar triple product a . (b x c), the signed volume of the tetrahedron
//  formed with the origin
static inline scalar triple
(
    const point& a,
    const point& b,
    const point& c
)
{
    return a & (b ^ c);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

cutCellPolyhedron::cutCellPolyhedron
(
    const pointField& pts,
    const faceList& fcs
)
:
    points_(pts),
    faces_(fcs)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar cutCellPolyhedron::volume() const
{
    //- Reference the tetrahedra to a point of the polyhedron. Taking the
    //- origin as the apex instead cancels terms of order the coordinate cubed
    //- against the volume, which loses most of the significant digits for a
    //- thin cell far from the origin.
    const point& ref = points_[0];

    scalar sumT = 0;

    forAll(faces_, fi)
    {
        const face& f = faces_[fi];

        const vector p0(points_[f[0]] - ref);

        for (label i = 1; i < f.size() - 1; i++)
        {
            sumT += triple(p0, points_[f[i]] - ref, points_[f[i + 1]] - ref);
        }
    }

    return sumT/6.0;
}


point cutCellPolyhedron::centre() const
{
    //- Referenced for the same reason as the volume, and the reference is added
    //- back to the volume-weighted result at the end
    const point& ref = points_[0];

    scalar sumT = 0;
    vector sumTc = Zero;

    forAll(faces_, fi)
    {
        const face& f = faces_[fi];

        const vector p0(points_[f[0]] - ref);

        for (label i = 1; i < f.size() - 1; i++)
        {
            const vector p1(points_[f[i]] - ref);
            const vector p2(points_[f[i + 1]] - ref);

            const scalar t = triple(p0, p1, p2);

            sumT += t;
            sumTc += t*(p0 + p1 + p2);
        }
    }

    if (mag(sumT) < VSMALL)
    {
        return Zero;
    }

    return ref + sumTc/(4.0*sumT);
}


label cutCellPolyhedron::nBadEdges() const
{
    Map<label> edgeCount;

    forAll(faces_, fi)
    {
        const face& f = faces_[fi];

        if (f.size() < 3)
        {
            return 1;
        }

        forAll(f, i)
        {
            const label a = f[i];
            const label b = f[f.fcIndex(i)];

            //- Pack the edge so that its direction does not matter
            const label key =
                min(a, b)*max(label(1), label(points_.size())) + max(a, b);

            Map<label>::iterator iter = edgeCount.find(key);

            if (iter == edgeCount.end())
            {
                edgeCount.insert(key, 1);
            }
            else
            {
                iter()++;
            }
        }
    }

    label nBad = 0;

    forAllConstIters(edgeCount, iter)
    {
        if (iter() != 2)
        {
            nBad++;
        }
    }

    return nBad;
}


bool cutCellPolyhedron::closed() const
{
    return nBadEdges() == 0;
}


bool cutCellPolyhedron::split
(
    const plane& pln,
    cutCellPolyhedron& negative,
    cutCellPolyhedron& positive,
    pointField& cutFace
) const
{
    cutCellDetail detail;

    return split(pln, negative, positive, cutFace, detail);
}


bool cutCellPolyhedron::split
(
    const plane& pln,
    cutCellPolyhedron& negative,
    cutCellPolyhedron& positive,
    pointField& cutFace,
    cutCellDetail& detail
) const
{
    negative = cutCellPolyhedron();
    positive = cutCellPolyhedron();
    cutFace.clear();
    detail.negativeSource.clear();
    detail.positiveSource.clear();
    detail.interfaceSource.clear();

    if (points_.empty() || faces_.empty())
    {
        return false;
    }

    const vector& n = pln.normal();

    //- Tolerance for deciding whether a point lies on the plane. The scale is
    //- the smallest extent of the polyhedron rather than the largest, so that a
    //- thin cell keeps a tolerance small compared with its thickness, with a
    //- floor from the coordinate magnitude so that a vertex that is exactly on
    //- the plane is still recognised when the cell is far from the origin.
    scalar minExtent = GREAT;
    scalar maxDist = 0;

    for (direction axis = 0; axis < 3; axis++)
    {
        scalar lo = GREAT;
        scalar hi = -GREAT;

        forAll(points_, pi)
        {
            lo = min(lo, points_[pi][axis]);
            hi = max(hi, points_[pi][axis]);
        }

        minExtent = min(minExtent, hi - lo);
    }

    forAll(points_, pi)
    {
        maxDist = max(maxDist, mag(points_[pi]));
    }

    const scalar tol = 1e-9*max(minExtent, SMALL) + 1e-12*maxDist;

    //- Signed distance of every point, and whether both sides are occupied
    scalarField dist(points_.size());

    bool hasNeg = false;
    bool hasPos = false;

    forAll(points_, pi)
    {
        dist[pi] = pln.signedDistance(points_[pi]);

        if (dist[pi] > tol)
        {
            hasPos = true;
        }
        else if (dist[pi] < -tol)
        {
            hasNeg = true;
        }
    }

    if (!hasNeg || !hasPos)
    {
        //- The plane misses the polyhedron, or only touches it
        return false;
    }


    //- A piece contains every point on its own side, plus every point on the
    //- plane, which is shared
    labelList toNeg(points_.size(), -1);
    labelList toPos(points_.size(), -1);

    DynamicList<point> negPts;
    DynamicList<point> posPts;

    //- Where each of those points came from, in the same order
    DynamicList<cutPointSource> negSrc;
    DynamicList<cutPointSource> posSrc;

    forAll(points_, pi)
    {
        if (dist[pi] <= tol)
        {
            toNeg[pi] = negPts.size();
            negPts.append(points_[pi]);
            negSrc.append(cutPointSource(pi, -1));
        }

        if (dist[pi] >= -tol)
        {
            toPos[pi] = posPts.size();
            posPts.append(points_[pi]);
            posSrc.append(cutPointSource(pi, -1));
        }
    }


    //- Intersection points, one per cut edge, keyed by the edge so that the
    //- two faces sharing it agree and the pieces stay closed. Each piece has
    //- its own numbering.
    Map<label> negEdgePts;
    Map<label> posEdgePts;

    auto edgePoint = [&](const label a, const label b)
    {
        const scalar w = dist[a]/(dist[a] - dist[b]);

        return points_[a] + w*(points_[b] - points_[a]);
    };

    auto cutEdgePoint = [&]
    (
        Map<label>& edgePts,
        DynamicList<point>& pts,
        DynamicList<cutPointSource>& src,
        const label a,
        const label b
    )
    {
        const label key =
            min(a, b)*max(label(1), label(points_.size())) + max(a, b);

        Map<label>::const_iterator iter = edgePts.find(key);

        if (iter != edgePts.end())
        {
            return iter();
        }

        const label newPtI = pts.size();

        pts.append(edgePoint(a, b));
        src.append(cutPointSource(a, b));
        edgePts.insert(key, newPtI);

        return newPtI;
    };


    //- The interface polygon, as the labels it has in each piece
    DynamicList<label> cutNegL;
    DynamicList<label> cutPosL;
    DynamicList<point> cutPts;

    //- Where each interface point came from, in the same order
    DynamicList<cutPointSource> cutSrc;


    auto toFace = [](const DynamicList<label>& lst)
    {
        face f(lst.size());

        forAll(lst, i)
        {
            f[i] = lst[i];
        }

        return f;
    };


    faceList negFaces;
    faceList posFaces;

    forAll(faces_, fi)
    {
        const face& f = faces_[fi];

        DynamicList<label> negPoly;
        DynamicList<label> posPoly;

        forAll(f, i)
        {
            const label a = f[i];
            const label b = f[f.fcIndex(i)];

            if (dist[a] <= tol)
            {
                negPoly.append(toNeg[a]);
            }

            if (dist[a] >= -tol)
            {
                posPoly.append(toPos[a]);
            }

            if
            (
                (dist[a] > tol && dist[b] < -tol)
             || (dist[a] < -tol && dist[b] > tol)
            )
            {
                //- An edge is shared by two faces, so its crossing is met
                //- twice. It must appear once in the interface polygon, or the
                //- polygon picks up duplicate points.
                const label key =
                    min(a, b)*max(label(1), label(points_.size())) + max(a, b);

                const bool isNew = !negEdgePts.found(key);

                negPoly.append(cutEdgePoint(negEdgePts, negPts, negSrc, a, b));
                posPoly.append(cutEdgePoint(posEdgePts, posPts, posSrc, a, b));

                if (isNew)
                {
                    cutNegL.append(negPoly.last());
                    cutPosL.append(posPoly.last());
                    cutPts.append(edgePoint(a, b));
                    cutSrc.append(cutPointSource(a, b));
                }
            }
        }

        if (negPoly.size() >= 3)
        {
            negFaces.append(toFace(negPoly));
        }

        if (posPoly.size() >= 3)
        {
            posFaces.append(toFace(posPoly));
        }
    }


    //- Points of the polyhedron lying on the plane are also vertices of the
    //- interface
    forAll(points_, pi)
    {
        if (mag(dist[pi]) <= tol)
        {
            cutNegL.append(toNeg[pi]);
            cutPosL.append(toPos[pi]);
            cutPts.append(points_[pi]);
            cutSrc.append(cutPointSource(pi, -1));
        }
    }


    if (cutPts.size() < 3)
    {
        return false;
    }


    //- Order the interface counter-clockwise about the plane normal. A centre
    //- and an in-plane basis are enough because the cross-section is convex.
    point c = Zero;

    forAll(cutPts, pi)
    {
        c += cutPts[pi];
    }
    c /= scalar(cutPts.size());

    vector e1 = Zero;
    scalar maxR2 = 0;

    forAll(cutPts, pi)
    {
        const vector d = cutPts[pi] - c;

        if (magSqr(d) > maxR2)
        {
            maxR2 = magSqr(d);
            e1 = d;
        }
    }

    if (maxR2 < SMALL)
    {
        //- Degenerate cross-section
        return false;
    }

    e1 /= mag(e1);

    const vector e2 = n ^ e1;

    scalarField angle(cutPts.size());

    forAll(cutPts, pi)
    {
        const vector d = cutPts[pi] - c;

        angle[pi] = std::atan2(d & e2, d & e1);
    }

    //- Sort by angle, carrying the per-piece labels along
    for (label i = 0; i < cutPts.size() - 1; i++)
    {
        for (label j = i + 1; j < cutPts.size(); j++)
        {
            if (angle[j] < angle[i])
            {
                Swap(angle[i], angle[j]);
                Swap(cutPts[i], cutPts[j]);
                Swap(cutNegL[i], cutNegL[j]);
                Swap(cutPosL[i], cutPosL[j]);
                Swap(cutSrc[i], cutSrc[j]);
            }
        }
    }


    //- The interface is outward for each piece, so its orientation is reversed
    //- on the positive side
    DynamicList<label> cutNegOrd;
    DynamicList<label> cutPosOrd;

    forAll(cutPts, i)
    {
        cutNegOrd.append(cutNegL[i]);
    }

    forAllReverse(cutPts, i)
    {
        cutPosOrd.append(cutPosL[i]);
    }

    negFaces.append(toFace(cutNegOrd));
    posFaces.append(toFace(cutPosOrd));


    pointField negPoints(negPts.size());

    forAll(negPts, i)
    {
        negPoints[i] = negPts[i];
    }

    pointField posPoints(posPts.size());

    forAll(posPts, i)
    {
        posPoints[i] = posPts[i];
    }

    cutFace.setSize(cutPts.size());

    forAll(cutPts, i)
    {
        cutFace[i] = cutPts[i];
    }

    detail.negativeSource.setSize(negPts.size());

    forAll(negPts, i)
    {
        detail.negativeSource[i] = negSrc[i];
    }

    detail.positiveSource.setSize(posPts.size());

    forAll(posPts, i)
    {
        detail.positiveSource[i] = posSrc[i];
    }

    detail.interfaceSource.setSize(cutPts.size());

    forAll(cutPts, i)
    {
        detail.interfaceSource[i] = cutSrc[i];
    }

    negative = cutCellPolyhedron(negPoints, negFaces);
    positive = cutCellPolyhedron(posPoints, posFaces);

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
