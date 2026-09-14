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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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


cutState cutCellByTriangles
(
    const cutCellPolyhedron& cell,
    const List<pointField>& tris,
    const scalar tol,
    cutCellPolyhedron& minus,
    cutCellPolyhedron& plus,
    pointField& interface
)
{
    minus = cutCellPolyhedron();
    plus = cutCellPolyhedron();
    interface.clear();

    if (cell.empty())
    {
        return cutState::untouched;
    }

    //- The cell box, for the triangle intersection test and for clipping the
    //- triangles down to the cell
    const boundBox bb(cell.points(), false);

    //- The distinct planes of the triangle parts inside the cell, and the total
    //- area they cover on each
    DynamicList<plane> planes;
    DynamicList<scalar> planeArea;

    forAll(tris, ti)
    {
        const pointField poly(triangleBoxPolygon(tris[ti], bb, tol));

        if (poly.size() < 3)
        {
            continue;
        }

        const plane pln(poly[0], poly[1], poly[2]);

        label found = -1;

        forAll(planes, pi)
        {
            //- Same plane, allowing for the opposite orientation
            if
            (
                mag(planes[pi].normal() & pln.normal()) > 1 - tol
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
            planeArea.append(0);
            found = planes.size() - 1;
        }

        planeArea[found] += polygonArea(poly);
    }

    if (planes.empty())
    {
        return cutState::untouched;
    }

    if (planes.size() > 1)
    {
        //- A surface edge or corner passes through the cell. Splitting by every
        //- plane here would create pieces in the same region, which needs the
        //- region classification and merging of the next slice.
        return cutState::refused;
    }

    if (!cell.split(planes[0], minus, plus, interface))
    {
        //- The plane touches the cell but does not divide it
        return cutState::untouched;
    }

    //- The plane must span the whole cross-section of the cell. Where it only
    //- covers part of it, the surface enters and leaves the cell, and the two
    //- pieces can lie in the same region.
    //-
    //- The tolerance is relative to the area of the plane's own polygons, with
    //- no absolute floor: a floor tied to a unit area would, at the scale of a
    //- thin feature, amount to accepting any coverage at all.
    if
    (
        mag(polygonArea(interface) - planeArea[0])
      > 1e-6*max(SMALL, planeArea[0])
    )
    {
        minus = cutCellPolyhedron();
        plus = cutCellPolyhedron();
        interface.clear();

        return cutState::refused;
    }

    return cutState::cut;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
