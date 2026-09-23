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

namespace
{

//- The side of a convex piece against a plane. The sign of the largest distance
//- over its points is taken rather than the sign at its centre, which a sliver
//- would make meaningless.
label pieceSide(const cutCellPolyhedron& piece, const plane& pln)
{
    scalar best = 0;

    forAll(piece.points(), pi)
    {
        const scalar d = pln.signedDistance(piece.points()[pi]);

        if (mag(d) > mag(best))
        {
            best = d;
        }
    }

    return (best > 0 ? 1 : -1);
}


//- Order the pieces by their side against the planes, negative first
bool lowerSide(const labelList& a, const labelList& b)
{
    forAll(a, i)
    {
        if (a[i] != b[i])
        {
            return a[i] < b[i];
        }
    }

    return false;
}


//- Do two faces of two pieces describe the same polygon? The clipper gives both
//- copies of a cut face the same point coordinates, one reversed, so matching
//- the point sets with a tolerance finds the faces the pieces share.
bool sameFace
(
    const face& a,
    const pointField& aPoints,
    const face& b,
    const pointField& bPoints,
    const scalar tol
)
{
    if (a.size() != b.size())
    {
        return false;
    }

    forAll(a, i)
    {
        bool found = false;

        forAll(b, j)
        {
            if (mag(aPoints[a[i]] - bPoints[b[j]]) <= tol)
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            return false;
        }
    }

    return true;
}

}   // End anonymous namespace


labelList groupPieces
(
    const label nPieces,
    const List<cutCellPieceFace>& sharedFaces,
    const labelList& pieceZone
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
        if (pieceZone[shared.pieceA] != pieceZone[shared.pieceB])
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


cutState cutCellByTriangles
(
    const cutCellPolyhedron& cell,
    const List<pointField>& tris,
    const scalar tol,
    cutCellResult& result
)
{
    result = cutCellResult();

    if (cell.empty())
    {
        return cutState::untouched;
    }

    const boundBox bb(cell.points(), false);

    //- The distinct planes of the parts of the triangles inside the cell.
    //- Facets of a closed surface that lie on one plane share its orientation,
    //- so the orientation of the plane is the orientation of the surface there.
    DynamicList<plane> planes(4);

    //- The triangle each plane came from, kept in step with planes
    DynamicList<label> planeTri(4);

    forAll(tris, ti)
    {
        const pointField poly(triangleBoxPolygon(tris[ti], bb, tol));

        if (poly.size() < 3)
        {
            continue;
        }

        result.nPolygons++;

        const plane pln(poly[0], poly[1], poly[2]);

        label found = -1;

        forAll(planes, pi)
        {
            if
            (
                (planes[pi].normal() & pln.normal()) > 1 - tol
             && mag
                (
                    planes[pi].normal()
                  & (pln.refPoint() - planes[pi].refPoint())
                ) <= tol
            )
            {
                found = pi;
                break;
            }
        }

        if (found < 0)
        {
            planes.append(pln);

            //- A new plane, so its provenance is this triangle. A plane
            //- already found keeps the triangle that created it: facets on one
            //- plane are one region in practice, and the first is as good as
            //- any if they are not.
            planeTri.append(ti);
        }
    }

    if (planes.empty())
    {
        return cutState::untouched;
    }

    //- A fixed order, so that a cell always cuts the same way. An index list is
    //- sorted rather than the planes themselves, so the triangle list takes the
    //- same permutation: sorting in place would leave the triangles in their old
    //- order and every patch silently wrong.
    labelList order(planes.size());

    forAll(order, pi)
    {
        order[pi] = pi;
    }

    std::sort
    (
        order.begin(),
        order.end(),
        [&planes](const label a, const label b)
        {
            const vector na(planes[a].normal());
            const vector nb(planes[b].normal());

            for (direction i = 0; i < 3; i++)
            {
                if (mag(na[i] - nb[i]) > SMALL)
                {
                    return na[i] < nb[i];
                }
            }

            return (planes[a].normal() & planes[a].refPoint())
                 < (planes[b].normal() & planes[b].refPoint());
        }
    );

    result.planes.setSize(planes.size());
    result.planeTri.setSize(planeTri.size());

    forAll(order, pi)
    {
        result.planes[pi] = planes[order[pi]];
        result.planeTri[pi] = planeTri[order[pi]];
    }

    //- The arrangement of the planes inside the cell: one plane at a time, on
    //- every piece, carrying on the pieces a plane does not divide
    DynamicList<cutCellPolyhedron> frontier(2);
    frontier.append(cell);

    forAll(result.planes, pi)
    {
        const plane& pln = result.planes[pi];

        DynamicList<cutCellPolyhedron> next(frontier.size() + 1);

        forAll(frontier, fi)
        {
            cutCellPolyhedron negative, positive;
            pointField cutFace;

            if (frontier[fi].split(pln, negative, positive, cutFace))
            {
                next.append(negative);
                next.append(positive);
            }
            else
            {
                //- The plane touches this piece but does not divide it
                next.append(frontier[fi]);
            }
        }

        frontier.transfer(next);
    }

    if (frontier.size() < 2)
    {
        //- The surface touches the cell without dividing it
        return cutState::untouched;
    }

    //- The pieces, with their side against every plane, ordered negative first
    result.pieces.setSize(frontier.size());

    forAll(frontier, fi)
    {
        cutCellPiece& piece = result.pieces[fi];

        piece.poly = frontier[fi];
        //- A polygon with two vertices within the tolerance is not a polygon. The
    //- driver resolves a position to one mesh point at this tolerance, so a face
    //- built from such a polygon would have a repeated point, no area, and would
    //- leave the cells either side open. The clipping itself works to a tighter
    //- tolerance, which is exactly why this has to be checked here: a polygon it
    //- considers distinct can still be one point to the driver.
    forAll(result.pieces, ppi)
    {
        const cutCellPolyhedron& poly = result.pieces[ppi].poly;

        forAll(poly.faces(), fi)
        {
            const face& f = poly.faces()[fi];

            forAll(f, i)
            {
                for (label j = i + 1; j < f.size(); j++)
                {
                    if
                    (
                        mag(poly.points()[f[i]] - poly.points()[f[j]]) <= tol
                    )
                    {
                        result.pieces.clear();
                        result.sharedFaces.clear();
                        result.refusal = cutRefusal::malformedPiece;

                        return cutState::refused;
                    }
                }
            }
        }
    }


    piece.side.setSize(result.planes.size());

        forAll(result.planes, pi)
        {
            piece.side[pi] = pieceSide(piece.poly, result.planes[pi]);
        }
    }

    std::sort
    (
        result.pieces.begin(),
        result.pieces.end(),
        [](const cutCellPiece& a, const cutCellPiece& b)
        {
            return lowerSide(a.side, b.side);
        }
    );

    //- The faces two pieces share. A cell has only a handful of pieces, so
    //- matching the polygons directly is affordable and needs no bookkeeping
    //- through the splits.
    forAll(result.pieces, a)
    {
        const cutCellPiece& pieceA = result.pieces[a];

        for (label b = a + 1; b < result.pieces.size(); b++)
        {
            const cutCellPiece& pieceB = result.pieces[b];

            forAll(pieceA.poly.faces(), facei)
            {
                forAll(pieceB.poly.faces(), facej)
                {
                    if
                    (
                        sameFace
                        (
                            pieceA.poly.faces()[facei], pieceA.poly.points(),
                            pieceB.poly.faces()[facej], pieceB.poly.points(),
                            tol
                        )
                    )
                    {
                        //- Two pieces that meet along a line rather than over an
                        //- area share a degenerate polygon, which cannot be a
                        //- mesh face: the cells either side would touch along an
                        //- edge with nothing between them. Refuse the whole cut
                        //- of this cell instead, so it stays one cell and the
                        //- driver reports it. The test is scale free, the area
                        //- against the squared perimeter, because a collinear
                        //- polygon's area comes out as cancellation noise and
                        //- is not reliably zero. The margin is wide: the
                        //- thinnest real sliver in the oblique case is around
                        //- 1e-8 of a face, which is 1e-8 of the squared
                        //- perimeter for a face about 0.5 across, against a
                        //- threshold of 1e-12.
                        const face& sharedFace = pieceA.poly.faces()[facei];
                        const pointField& sharedPts = pieceA.poly.points();

                        pointField ring(sharedFace.size());

                        forAll(sharedFace, k)
                        {
                            ring[k] = sharedPts[sharedFace[k]];
                        }

                        scalar perimeter = 0;

                        forAll(ring, k)
                        {
                            perimeter += mag(ring[(k + 1) % ring.size()] - ring[k]);
                        }

                        const scalar ringArea = polygonArea(ring);

                        if
                        (
                            ringArea <= tol*tol
                         || ringArea <= 1e-12*perimeter*perimeter
                        )
                        {
                            result.pieces.clear();
                            result.sharedFaces.clear();
                            result.refusal =
                                cutRefusal::degenerateSharedFace;

                            return cutState::refused;
                        }

                        cutCellPieceFace shared;
                        shared.pieceA = a;
                        shared.faceInA = facei;
                        shared.pieceB = b;
                        shared.faceInB = facej;

                        //- The patch is inherited from the triangle whose
                        //- plane produced this face. This loop is over pieces
                        //- and faces and does not know the plane, so the face's
                        //- own points are matched against the result's planes,
                        //- the comparison the grouping itself performs. A face
                        //- that matches none keeps -1 rather than a guess.
                        {
                            const face& sf = pieceA.poly.faces()[facei];
                            const pointField& spts = pieceA.poly.points();

                            if (sf.size() >= 3)
                            {
                                const plane fpln
                                (
                                    spts[sf[0]], spts[sf[1]], spts[sf[2]]
                                );

                                forAll(result.planes, pi)
                                {
                                    if
                                    (
                                        (result.planes[pi].normal()
                                      & fpln.normal()) > 1 - tol
                                     && mag
                                        (
                                            result.planes[pi].normal()
                                          & (fpln.refPoint()
                                          - result.planes[pi].refPoint())
                                        ) <= tol
                                    )
                                    {
                                        shared.tri = result.planeTri[pi];
                                        break;
                                    }
                                }
                            }
                        }

                        result.sharedFaces.append(shared);
                    }
                }
            }
        }
    }

    return cutState::cut;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

cutState cutCellByTriangles
(
    const cutCellPolyhedron& cell,
    const List<pointField>& tris,
    const labelList& patchOf,
    const scalar tol,
    cutCellResult& result
)
{
    const cutState state = cutCellByTriangles(cell, tris, tol, result);

    if (state != cutState::cut)
    {
        return state;
    }

    //- A face inherits the region of the triangle that produced its plane.
    forAll(result.sharedFaces, i)
    {
        cutCellPieceFace& sf = result.sharedFaces[i];

        if (sf.tri >= 0 && sf.tri < patchOf.size())
        {
            sf.patch = patchOf[sf.tri];
        }
        else
        {
            sf.patch = -1;
        }
    }

    return state;
}
} // End namespace Foam

// ************************************************************************* //


