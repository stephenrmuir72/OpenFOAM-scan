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

#include "cutCellCutter.H"
#include "DynamicList.H"

#include <algorithm>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool polygonDegenerate(const pointField& poly, const scalar tol)
{
    if (poly.size() < 3)
    {
        return true;
    }

    scalar longest = 0;

    forAll(poly, i)
    {
        longest = max(longest, mag(poly[poly.fcIndex(i)] - poly[i]));
    }

    return polygonArea(poly) <= max(tol, VSMALL)*longest;
}


scalar polygonArea(const pointField& poly)
{
    if (poly.size() < 3)
    {
        return 0;
    }

    point c = Zero;

    forAll(poly, i)
    {
        c += poly[i];
    }
    c /= scalar(poly.size());

    vector A = Zero;

    forAll(poly, i)
    {
        A += (poly[i] - c) ^ (poly[(i + 1) % poly.size()] - c);
    }

    return 0.5*mag(A);
}


pointField clipPolygon
(
    const pointField& poly,
    const plane& pln,
    const scalar tol
)
{
    const label n = poly.size();

    if (n < 3)
    {
        return pointField();
    }

    scalarField dist(n);

    forAll(poly, i)
    {
        dist[i] = pln.signedDistance(poly[i]);
    }

    DynamicList<point> out(poly.size() + 1);

    forAll(poly, i)
    {
        const label j = (i + 1) % n;

        //- Keep the vertex when it is on or below the plane
        if (dist[i] <= tol)
        {
            out.append(poly[i]);
        }

        //- Where the edge crosses, add the crossing point
        if
        (
            (dist[i] < -tol && dist[j] > tol)
         || (dist[i] > tol && dist[j] < -tol)
        )
        {
            const scalar w = dist[i]/(dist[i] - dist[j]);

            out.append(poly[i] + w*(poly[j] - poly[i]));
        }
    }

    if (out.size() < 3)
    {
        return pointField();
    }

    pointField result(out.size());

    forAll(out, i)
    {
        result[i] = out[i];
    }

    return result;
}


List<pointField> outsideConvex
(
    const List<pointField>& qs,
    const pointField& r,
    const vector& n
)
{
    List<pointField> out;

    forAll(qs, qi)
    {
        pointField rest(qs[qi]);

        for (label k = 0; k < r.size() && rest.size() >= 3; ++k)
        {
            const point& p0 = r[k];
            const vector m(n ^ (r[(k + 1) % r.size()] - p0));
            const scalar magM = mag(m);

            if (magM <= VSMALL)
            {
                continue;
            }

            const vector mh = m/magM;

            const pointField outside(clipPolygon(rest, plane(p0, mh), 0));

            if (outside.size() >= 3)
            {
                out.append(outside);
            }

            rest = clipPolygon(rest, plane(p0, -mh), 0);
        }
    }

    return out;
}


pointField triangleCellPolygon
(
    const pointField& tri,
    const cutCellPolyhedron& cell,
    const scalar tol
)
{
    pointField poly(tri);

    const pointField& cPts = cell.points();

    forAll(cell.faces(), fi)
    {
        const face& f = cell.faces()[fi];

        const vector fn(f.normal(cPts));

        if (mag(fn) <= VSMALL)
        {
            continue;
        }

        const vector n(fn/mag(fn));
        const point fc(f.centre(cPts));

        //- The outward deviation of the cell beyond this face's plane. A snapped
        //- cell can be slightly non-planar, so a face bent by more than 10 tol
        //- is pushed out to contain every vertex. Smaller deviations are noise,
        //- and the face stays exact.
        scalar dev = 0;

        forAll(cPts, pi)
        {
            dev = max(dev, (cPts[pi] - fc) & n);
        }

        const scalar shift = dev > 10*tol ? dev : 0;

        //- Keep the part on the inner side of the face's plane
        poly = clipPolygon(poly, plane(fc + shift*n, n), tol);

        if (poly.size() < 3)
        {
            return pointField();
        }
    }

    return poly;
}


pointField triangleBoxPolygon
(
    const pointField& tri,
    const boundBox& bb,
    const scalar tol
)
{
    pointField poly(tri);

    for (direction axis = 0; axis < 3; axis++)
    {
        //- Keep the part with coordinate above the lower box plane. The normal
        //- points along -axis, so the inside is the negative side.
        vector nLo(Zero);
        nLo[axis] = -1;

        point pLo(Zero);
        pLo[axis] = bb.min()[axis];

        poly = clipPolygon(poly, plane(pLo, nLo), tol);

        if (poly.size() < 3)
        {
            return pointField();
        }

        //- Keep the part below the upper box plane
        vector nHi(Zero);
        nHi[axis] = 1;

        point pHi(Zero);
        pHi[axis] = bb.max()[axis];

        poly = clipPolygon(poly, plane(pHi, nHi), tol);

        if (poly.size() < 3)
        {
            return pointField();
        }
    }

    return poly;
}

bool samePieceZones(const labelList& a, const labelList& b)
{
    if (a.size() != b.size())
    {
        return false;
    }

    forAll(a, i)
    {
        if (a[i] != b[i])
        {
            return false;
        }
    }

    return true;
}


labelList groupPieces
(
    const label nPieces,
    const List<cutCellPieceFace>& sharedFaces,
    const List<labelList>& pieceZones
)
{
    //- Union-find over the pieces, with the group index of each piece kept as
    //- the parent, so no separate array is needed
    labelList group(nPieces);

    forAll(group, pi)
    {
        group[pi] = pi;
    }

    auto root = [&](label pi) -> label
    {
        while (group[pi] != pi)
        {
            group[pi] = group[group[pi]];
            pi = group[pi];
        }

        return pi;
    };

    forAll(sharedFaces, si)
    {
        const cutCellPieceFace& shared = sharedFaces[si];

        //- A face two pieces share is a face between two regions, which is the
        //- interface, or a face inside one region, which is what merges
        if
        (
            !samePieceZones
            (
                pieceZones[shared.pieceA],
                pieceZones[shared.pieceB]
            )
        )
        {
            continue;
        }

        const label a = root(shared.pieceA);
        const label b = root(shared.pieceB);

        if (a != b)
        {
            //- Keep the lower root, so the group of the earlier piece wins
            group[max(a, b)] = min(a, b);
        }
    }

    //- Number the groups from zero in the order their first piece appears
    labelList groupOf(nPieces, -1);
    labelList index(nPieces, -1);
    label nGroups = 0;

    forAll(group, pi)
    {
        const label r = root(pi);

        if (index[r] < 0)
        {
            index[r] = nGroups++;
        }

        groupOf[pi] = index[r];
    }

    return groupOf;
}


labelList groupPieces
(
    const label nPieces,
    const List<cutCellPieceFace>& sharedFaces,
    const labelList& pieceZone
)
{
    List<labelList> pieceZones(nPieces);

    forAll(pieceZones, pi)
    {
        pieceZones[pi] = labelList(1, pieceZone[pi]);
    }

    return groupPieces(nPieces, sharedFaces, pieceZones);
}


} // End namespace Foam

// ************************************************************************* //
