/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.com
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * //

//- Set by cutCellLoopDebug, reported by cutCellByLoop
static bool loopDebug = false;


void cutCellLoopDebug(const bool on)
{
    loopDebug = on;
}


namespace
{

//- True if two planes are the same plane, allowing for the tolerance
bool samePlane(const plane& a, const plane& b, const scalar tol)
{
    return
    (
        (a.normal() & b.normal()) > 1 - tol
     && mag(a.normal() & (a.refPoint() - b.refPoint())) <= tol
    );
}


//- The polygon where a plane cuts a convex polyhedron. Empty when it misses.
pointField planeSection
(
    const cutCellPolyhedron& cell,
    const plane& pln,
    const scalar tol
)
{
    const pointField& pts = cell.points();

    DynamicList<point> section(8);

    forAll(cell.faces(), fi)
    {
        const face& f = cell.faces()[fi];

        forAll(f, i)
        {
            const label a = f[i];
            const label b = f[f.fcIndex(i)];

            const scalar da = pln.signedDistance(pts[a]);
            const scalar db = pln.signedDistance(pts[b]);

            //- A vertex lying on the plane is part of the section. A strict
            //- crossing test alone loses it, and with it the sliver between that
            //- vertex and the crossings either side of it. At the sphere's pole
            //- the facet's plane passes through the cell's corner, and the section
            //- then came out 4.68e-05 SMALLER than the patch clipped to the same
            //- cell - which the patch cannot be, and the enumeration dropped that
            //- much surface in each of the six pole cells (26-Sep-2026, 07
            //- section 22).
            const bool onA = mag(da) <= tol;
            const bool onB = mag(db) <= tol;

            DynamicList<point> add(2);

            if (onA)
            {
                add.append(pts[a]);
            }

            if (onB)
            {
                add.append(pts[b]);
            }

            if (!onA && !onB && da*db < 0)
            {
                add.append(pts[a] + (da/(da - db))*(pts[b] - pts[a]));
            }

            //- The cell has no edge list, so its edges are walked once per face
            //- and each crossing is met twice. Kept once, or the section polygon
            //- carries coincident points and its edges cannot pair up.
            forAll(add, ai)
            {
                const point& x = add[ai];

                bool known = false;

                forAll(section, k)
                {
                    if (mag(section[k] - x) <= tol)
                    {
                        known = true;
                        break;
                    }
                }

                if (!known)
                {
                    section.append(x);
                }
            }
        }
    }

    if (section.size() < 3)
    {
        return pointField();
    }

    //- The cell is convex, so the section is a convex polygon and its points
    //- can be ordered about their own centre by angle in the plane.
    point c = Zero;

    forAll(section, i)
    {
        c += section[i];
    }

    c /= scalar(section.size());

    const vector n(pln.normal());
    const vector e1 = normalised(section[0] - c);
    const vector e2 = n ^ e1;

    pointField result(section.size());

    forAll(section, i)
    {
        const vector d(section[i] - c);

        result[i] = section[i];

        //- Stored with the angle, then sorted below
        section[i] = point(std::atan2(d & e2, d & e1), 0, 0);
    }

    labelList order(result.size());

    forAll(order, i)
    {
        order[i] = i;
    }

    std::sort
    (
        order.begin(),
        order.end(),
        [&](const label a, const label b) { return section[a].x() < section[b].x(); }
    );

    pointField sorted(result.size());

    forAll(order, i)
    {
        sorted[i] = result[order[i]];
    }

    return sorted;
}


//- True if a point lies on a polygonal face, within the tolerance
bool onFace
(
    const face& f,
    const pointField& pts,
    const point& x,
    const scalar tol
)
{
    if (f.size() < 3)
    {
        return false;
    }

    //- The face's own normal, and a point known to be inside it.
    //-
    //- From the face's area normal, not from its first three points. The
    //- conformality step adds points along a face's edges, so a face can
    //- begin with three collinear points - face 11 of 05_cube_cut cell 716
    //- is (0,-0.0625,-0.25) (0,-0.0625,-0.21875) (0,-0.0625,-0.1875) - and
    //- then the cross product is zero, the normal is degenerate, and the
    //- in-plane test below passes for EVERY point. Every chord was then read
    //- as lying on that face, the chain assembled for it was the whole closed
    //- loop, and a closed chain cannot split a region: two faces were left
    //- unsplit, one region straddled the cutting plane, the join bridged the
    //- two halves and the cell came out as one piece - untouched, with its
    //- surface and volume lost. The area normal sums over the polygon, so
    //- collinear leading points cannot degenerate it.
    const vector area(f.areaNormal(pts));

    if (mag(area) <= VSMALL)
    {
        return false;
    }

    const vector n(area/mag(area));

    point c = Zero;

    forAll(f, i)
    {
        c += pts[f[i]];
    }

    c /= scalar(f.size());

    if (mag(n & (x - pts[f[0]])) > tol)
    {
        return false;
    }

    //- Inside a convex polygon: every edge has x on the same side as the
    //- centre does. The sign is taken from the centre rather than assumed, so
    //- a face wound either way is handled.
    forAll(f, i)
    {
        const point& pa = pts[f[i]];
        const point& pb = pts[f[f.fcIndex(i)]];

        const scalar ac = ((pb - pa) ^ (c - pa)) & n;
        const scalar ax = ((pb - pa) ^ (x - pa)) & n;

        if (ac*ax < 0)
        {
            return false;
        }
    }

    return true;
}


//- True if a point lies on one of a convex polyhedron's edges
bool onCellEdge
(
    const cutCellPolyhedron& cell,
    const point& x,
    const scalar tol
)
{
    const pointField& pts = cell.points();

    forAll(cell.faces(), fi)
    {
        const face& f = cell.faces()[fi];

        forAll(f, i)
        {
            const point& a = pts[f[i]];
            const point& b = pts[f[f.fcIndex(i)]];

            const vector ab(b - a);
            const scalar t = ((x - a) & ab)/(ab & ab);

            if
            (
                t >= -SMALL && t <= 1 + SMALL
             && mag(x - (a + t*ab)) <= tol
            )
            {
                return true;
            }
        }
    }

    return false;
}


//- True if a point lies on the boundary of a convex polyhedron
bool onBoundary
(
    const cutCellPolyhedron& cell,
    const point& x,
    const scalar tol
)
{
    forAll(cell.faces(), fi)
    {
        if (onFace(cell.faces()[fi], cell.points(), x, tol))
        {
            return true;
        }
    }

    return false;
}


//- The distance within which two points are the same mesh point. The cell's
//- geometry produces a vertex on more than one polygon - a facet's clipping, a
//- region's walk, an interface part - and each carries its own arithmetic
//- error, so identity needs a margin over the tolerance itself. Ten times is
//- the margin the rest of this file was already written with; naming it here
//- keeps the builder and every matcher that reads what it built in step.
inline scalar samePointTol(const scalar tol)
{
    return 10*tol;
}


//- Build a polyhedron from a set of polygons that share their edges, keeping
//- one point per distinct location and orienting every face outwards
//-
//- The merge has to admit that margin: measured on 06_oblique_box, whose
//- oblique intersections carry more error than an axis-aligned case, 146
//- vertices missed a merge at 1.07 to 1.20 times tol, and one pair at 2.2. A
//- vertex that misses is a second label for the same point, the faces either
//- side of it do not pair up, and the polyhedron is then not a closed solid -
//- its volume and centre, integrals over the whole boundary, mean nothing. On
//- 07_sphere_cut the same probe found no near-miss at all, so this is not a
//- general widening, it is the measured margin. The distinct vertices of a
//- piece are a cell edge or more apart, on the order of a thousandth of the
//- cell, which ten times a tolerance that is a billionth of the domain cannot
//- merge.
cutCellPolyhedron polyhedronFrom
(
    const List<pointField>& polys,
    const point& inside,
    const scalar tol
)
{
    const scalar mergeTol = samePointTol(tol);

    pointField pts;
    faceList fcs(polys.size());

    forAll(polys, pi)
    {
        const pointField& p = polys[pi];

        face f(p.size());

        forAll(p, i)
        {
            label found = -1;

            forAll(pts, j)
            {
                if (mag(pts[j] - p[i]) <= mergeTol)
                {
                    found = j;
                    break;
                }
            }

            if (found < 0)
            {
                found = pts.size();
                pts.append(p[i]);
            }

            f[i] = found;
        }

        //- Outward orientation, tested against a point known to be inside
        const point& p0 = pts[f[0]];

        vector sum = Zero;

        for (label i = 1; i < f.size() - 1; i++)
        {
            sum += (pts[f[i]] - p0) ^ (pts[f[i + 1]] - p0);
        }

        if ((sum & (p0 - inside)) < 0)
        {
            reverse(f);
        }

        fcs[pi] = f;
    }

    return cutCellPolyhedron(pts, fcs);
}


//- How many of the polyhedron's edges are not shared by exactly two of its
//- faces. Not used as a guard: a piece's own faces legitimately meet at
//- T-junctions, because a facet's clipped polygon and a region's walk do not
//- subdivide the chain alike, and the mesher resolves those when it emits the
//- faces. Kept for a dump, where the count says how far from closed a piece is.
label nOpenEdges(const cutCellPolyhedron& poly)
{
    const faceList& fcs = poly.faces();

    DynamicList<Pair<label>> edges;
    edges.reserve(4*fcs.size());

    forAll(fcs, fi)
    {
        const face& f = fcs[fi];

        forAll(f, i)
        {
            edges.append
            (
                Pair<label>(min(f[i], f.fcIndex(i)), max(f[i], f.fcIndex(i)))
            );
        }
    }

    label nOpen = 0;

    forAll(edges, ei)
    {
        label n = 0;

        forAll(edges, ej)
        {
            if (edges[ej] == edges[ei])
            {
                n++;
            }
        }

        if (n != 2)
        {
            nOpen++;
        }
    }

    return nOpen;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

cutState cutCellByLoop
(
    const cutCellPolyhedron& cell,
    const List<pointField>& tris,
    const labelList& patchOf,
    const scalar tol,
    cutCellResult& result,
    const bool dump
)
{
    const boundBox bb(cell.points(), false);
    const pointField& cPts = cell.points();

    //- 1. Clip every facet to the cell, keeping each one, and group them by
    //- their supporting plane. The group is what becomes one interface face,
    //- so a run of coplanar facets of a developable surface stays one face.
    List<pointField> facet;
    labelList facetTri;

    List<plane> planes;
    labelList planeTri;
    List<labelList> planeFacets;

    forAll(tris, ti)
    {
        pointField p(triangleBoxPolygon(tris[ti], bb, tol));

        if (p.size() < 3)
        {
            continue;
        }

        result.nPolygons++;

        const label i = facet.size();
        facet.append(p);
        facetTri.append(ti);

        const plane pln(p[0], p[1], p[2]);

        label found = -1;

        forAll(planes, pi)
        {
            if (samePlane(planes[pi], pln, tol))
            {
                found = pi;
                break;
            }
        }

        if (found < 0)
        {
            found = planes.size();
            planes.append(pln);
            planeTri.append(ti);
            planeFacets.append(labelList());
        }

        planeFacets[found].append(i);
    }

    if (planes.empty())
    {
        return cutState::untouched;
    }

    result.planes = planes;
    result.planeTri = planeTri;

    //- 2. The loop: the clipped facets' edges that lie on the cell's boundary.
    //- Those are the curves that separate the cell's boundary into regions.
    //- Edges interior to the cell, coplanar ones included, are not part of it.
    List<pointField> chords;

    //- The facet each chord came from. Without it the chords cannot be linked,
    //- and linking them by coordinates alone merges curves that only touch: two
    //- crossing sheets meet at a point where four chords end, and the adjacency
    //- of the facets says which two belong to each sheet.
    labelList chordFacet;

    forAll(facet, fi)
    {
        const pointField& p = facet[fi];

        forAll(p, i)
        {
            const point& a = p[i];
            const point& b = p[p.fcIndex(i)];

            if
            (
                mag(a - b) > tol
             && onBoundary(cell, 0.5*(a + b), tol)
            )
            {
                chords.append(pointField({a, b}));
                chordFacet.append(fi);
            }
        }
    }

    //- True when two facets meet along an edge that ends at a point, which is
    //- what makes two chords consecutive on the same surface. A structural test
    //- of shared edges, not a coincidence of coordinates: this is REF-1 2.6.1's
    //- twin pointers, reached by comparing the edges rather than storing them.
    auto adjacentAt = [&](const label fa, const pointField& Ca,
                          const label fb, const pointField& Cb,
                          const point& Q) -> bool
    {
        //- An edge of a facet's polygon matching a chord, ending at Q. The
        //- chords ARE the facets' edges, so this asks whether the facet has
        //- that edge at all and whether Q is one of its ends.
        auto edgeAt = [&](const label fi, const pointField& C) -> bool
        {
            const pointField& P = facet[fi];

            forAll(P, i)
            {
                const point& p0 = P[i];
                const point& p1 = P[P.fcIndex(i)];

                if (mag(p0 - Q) > 10*tol && mag(p1 - Q) > 10*tol)
                {
                    continue;
                }

                if
                (
                    (mag(p0 - C[0]) <= 10*tol && mag(p1 - C[1]) <= 10*tol)
                 || (mag(p0 - C[1]) <= 10*tol && mag(p1 - C[0]) <= 10*tol)
                )
                {
                    return true;
                }
            }

            return false;
        };

        if (fa == fb)
        {
            //- Consecutive edges of one facet's own polygon, which means both
            //- of them end at Q. Returning true for any two edges of the same
            //- facet was wrong: a facet clipped so that it meets the cell's
            //- boundary in three edges could have two chords linked that are
            //- not neighbours, leaving the real neighbour unlinked and the walk
            //- with a wrong turn in it.
            return edgeAt(fa, Ca) && edgeAt(fb, Cb);
        }

        const pointField& A = facet[fa];
        const pointField& B = facet[fb];

        forAll(A, i)
        {
            const point& a0 = A[i];
            const point& a1 = A[A.fcIndex(i)];

            if (mag(a0 - Q) > 10*tol && mag(a1 - Q) > 10*tol)
            {
                continue;
            }

            forAll(B, j)
            {
                const point& b0 = B[j];
                const point& b1 = B[B.fcIndex(j)];

                if
                (
                    (mag(a0 - b0) <= 10*tol && mag(a1 - b1) <= 10*tol)
                 || (mag(a0 - b1) <= 10*tol && mag(a1 - b0) <= 10*tol)
                )
                {
                    return true;
                }
            }
        }

        return false;
    };

    //- 3. Split each cell face along the CHAINS the loop makes on it. The chords
    //- on a face are consecutive segments of one curve, meeting at interior
    //- vertices, so cutting by each chord's LINE cuts far past the curve.
    //- Measured on one cell of the sphere: a face carrying four chords in a
    //- chain makes TWO regions, and the line cuts made EIGHT. Everything
    //- downstream followed from it - regions with no cell edge to stitch to,
    //- pieces that do not close, eight pieces where the geometry owes two.
    //- Cutting along the chain gives the regions as plain boundary walks, so no
    //- cut ever extends beyond the loop.
    List<List<pointField>> faceRegions(cell.nFaces());

    label nChains = 0;
    label nChainsSkipped = 0;

    //- Set when a chain closed on itself and split no region of its face. Such a
    //- chain cannot split that face at all - its far side would be the face minus
    //- the loop's interior, an annulus - so the cell is refused rather than cut
    //- with the face left whole. See the refusal at the end of this function.
    bool closedChainSplitNothing = false;

    //- Where a point sits on a polygon's boundary, as an edge and a parameter
    auto onBoundaryEdge = [&](const pointField& poly, const point& x,
                              label& edgei, scalar& t) -> bool
    {
        forAll(poly, i)
        {
            const point& a = poly[i];
            const point& b = poly[poly.fcIndex(i)];

            const vector ab(b - a);
            const scalar len2 = ab & ab;

            if (len2 <= VSMALL)
            {
                continue;
            }

            const scalar u = ((x - a) & ab)/len2;

            if
            (
                u >= -SMALL && u <= 1 + SMALL
             && mag(x - (a + u*ab)) <= 10*tol
            )
            {
                edgei = i;
                t = u;
                return true;
            }
        }

        return false;
    };

    //- A polygon's boundary walked forward from one position to another, the
    //- walk starting and ending at the positions themselves
    auto walk = [&](const pointField& poly, const label ea, const scalar ta,
                    const label eb, const scalar tb, pointField& out) -> bool
    {
        const label n = poly.size();

        const point P = poly[ea] + ta*(poly[poly.fcIndex(ea)] - poly[ea]);
        const point Q = poly[eb] + tb*(poly[poly.fcIndex(eb)] - poly[eb]);

        out.clear();
        out.append(P);

        if (ea == eb)
        {
            //- Both positions on one edge, and this walk is the part between
            //- them, which exists only when the first comes first
            if (ta < tb - SMALL)
            {
                out.append(Q);
                return true;
            }

            return false;
        }

        label i = poly.fcIndex(ea);
        label nSteps = 0;

        while (i != eb)
        {
            out.append(poly[i]);
            i = poly.fcIndex(i);
            nSteps++;

            if (nSteps > n)
            {
                return false;
            }
        }

        out.append(poly[eb]);
        out.append(Q);

        return true;
    };

    //- Split a region along a chain whose two ends lie on its boundary. The two
    //- regions are the two boundary arcs, each closed by the chain.
    auto splitByChain = [&](const pointField& region, const pointField& chain,
                            pointField& A, pointField& B) -> bool
    {
        if (region.size() < 3 || chain.size() < 2)
        {
            return false;
        }

        label ea = -1;
        label eb = -1;
        scalar ta = 0;
        scalar tb = 0;

        if (!onBoundaryEdge(region, chain.first(), ea, ta)) return false;
        if (!onBoundaryEdge(region, chain.last(), eb, tb)) return false;

        //- A chain that begins and ends at one point closes on itself inside the
        //- face, and its two sides cannot both be faces of a polyMesh. It is
        //- counted and left alone rather than guessed at.
        if (ea == eb && mag(ta - tb) <= 1e-9)
        {
            return false;
        }

        pointField arcA;
        pointField arcB;

        if (!walk(region, ea, ta, eb, tb, arcA))
        {
            return false;
        }

        if (!walk(region, eb, tb, ea, ta, arcB))
        {
            if (ea != eb)
            {
                return false;
            }

            //- Both positions on one edge: the second arc is the rest of the
            //- boundary, from the far position all the way round to the near one
            arcB.clear();
            arcB.append(arcA.last());

            label i = region.fcIndex(eb);

            forAll(region, k)
            {
                arcB.append(region[i]);
                i = region.fcIndex(i);
            }

            arcB.append(arcA.first());
        }

        A = arcA;

        for (label k = chain.size() - 2; k >= 1; k--)
        {
            A.append(chain[k]);
        }

        B = arcB;

        for (label k = 1; k <= chain.size() - 2; k++)
        {
            B.append(chain[k]);
        }

        return A.size() >= 3 && B.size() >= 3;
    };

    //- The first part of a chain that lies inside a region with both of its ends
    //- on that region's boundary. A chain that crosses a cut made by an earlier
    //- chain has its ends in different regions, so it must be divided at the
    //- crossing before either part can be a cut.
    auto chainInRegion = [&](const pointField& chain, const pointField& region,
                             pointField& sub) -> bool
    {
        if (region.size() < 3 || chain.size() < 2)
        {
            return false;
        }

        //- The region's inward normals, in the face's plane, each sign taken
        //- from the region's own centroid
        point rc = Zero;

        forAll(region, k)
        {
            rc += region[k];
        }

        rc /= scalar(region.size());

        //- The region's own normal, from the area normal of its polygon. Its
        //- first three points can be collinear once the conformality step has
        //- added points along an edge - face 11's region of 05_cube_cut cell
        //- 716 begins (0,-0.0625,-0.25) (0,-0.0625,-0.21875)
        //- (0,-0.0625,-0.1875) - and then the cross product is zero, every
        //- inward normal below is degenerate, and the clipping accepts nothing:
        //- the chain that crosses the region is rejected, the region is never
        //- split, and the cell comes out as one piece.
        vector areaSum = Zero;

        forAll(region, k)
        {
            areaSum += region[k] ^ region[region.fcIndex(k)];
        }

        areaSum *= 0.5;

        if (mag(areaSum) <= VSMALL)
        {
            return false;
        }

        const vector rn(areaSum/mag(areaSum));

        List<vector> inward(region.size());

        forAll(region, k)
        {
            const point& ra = region[k];
            const point& rb = region[region.fcIndex(k)];

            const vector e = (rb - ra) ^ rn;

            inward[k] = (((rc - ra) & e) >= 0 ? e : -e);
        }

        //- The segments of the chain, clipped to the region, merged into runs
        DynamicList<point> run(chain.size() + 2);
        bool haveRun = false;

        for (label i = 1; i < chain.size(); i++)
        {
            const point& p0 = chain[i - 1];
            const point& p1 = chain[i];

            scalar t0 = 0;
            scalar t1 = 1;

            forAll(region, k)
            {
                const point& ra = region[k];

                const scalar da = inward[k] & (p0 - ra);
                const scalar db = inward[k] & (p1 - ra);

                if (da < -10*tol && db < -10*tol)
                {
                    t0 = 1;
                    t1 = 0;
                    break;
                }

                if (da >= -10*tol && db >= -10*tol)
                {
                    continue;
                }

                const scalar t = da/(da - db);

                if (da < 0)
                {
                    t0 = max(t0, t);
                }
                else
                {
                    t1 = min(t1, t);
                }
            }

            if (t0 >= t1)
            {
                if (haveRun)
                {
                    break;
                }

                continue;
            }

            const point A = p0 + t0*(p1 - p0);
            const point B = p0 + t1*(p1 - p0);

            if (!haveRun)
            {
                run.clear();
                run.append(A);
                haveRun = true;
            }
            else if (mag(run.last() - A) > 10*tol)
            {
                //- The chain left the region and came back
                break;
            }

            run.append(B);
        }

        if (run.size() < 2)
        {
            return false;
        }

        sub = run;

        //- Both ends must be on the region's boundary to be a cut of it
        label ea = -1;
        label eb = -1;
        scalar ta = 0;
        scalar tb = 0;

        if
        (
            !onBoundaryEdge(region, sub.first(), ea, ta)
         || !onBoundaryEdge(region, sub.last(), eb, tb)
         || (ea == eb && mag(ta - tb) <= 1e-9)
        )
        {
            return false;
        }

        //- And the chain must run through the region, not along its boundary.
        //- Once a region has been cut by this chain, the chain lies on the new
        //- region's own boundary and satisfies the two ends test again, so
        //- without this the retry cuts the same chain over and over: a hex cut
        //- by one plane came out as 38 regions and 30 pieces.
        {
            //- A segment midpoint strictly inside the region. Not the mean of
            //- the chain's points, which for a bent chain need not be inside
            //- the region at all and rejected every curved patch in the sphere.
            bool through = false;

            for (label k = 0; k < sub.size() - 1 && !through; k++)
            {
                const point mid(0.5*(sub[k] + sub[k + 1]));

                through = true;

                forAll(region, j)
                {
                    if ((inward[j] & (mid - region[j])) <= 10*tol)
                    {
                        through = false;
                        break;
                    }
                }
            }

            if (!through)
            {
                return false;
            }
        }

        return true;
    };

    forAll(cell.faces(), ffi)
    {
        const face& f = cell.faces()[ffi];

        pointField fp(f.size());

        forAll(f, i)
        {
            fp[i] = cPts[f[i]];
        }

        faceRegions[ffi].append(fp);

        //- The chords on this face
        labelList onThisFace;

        forAll(chords, ci)
        {
            if (onFace(f, cPts, 0.5*(chords[ci][0] + chords[ci][1]), tol))
            {
                onThisFace.append(ci);
            }
        }

        if (onThisFace.empty())
        {
            continue;
        }

        //- How many chords of this face end at a point. Where only two meet
        //- there is nothing to disambiguate and they must continue each other;
        //- adjacency is needed only where more than two meet, which is where
        //- two curves cross. Requiring adjacency everywhere left 192 of the
        //- sphere's cells whole: a loop that turns at a point where two facets
        //- meet without sharing an edge was refused.
        auto degreeAt = [&](const point& Q) -> label
        {
            label n = 0;

            forAll(onThisFace, ci)
            {
                if
                (
                    mag(chords[onThisFace[ci]][0] - Q) <= tol
                 || mag(chords[onThisFace[ci]][1] - Q) <= tol
                )
                {
                    n++;
                }
            }

            return n;
        };

        //- Link them into chains by their ends. Each end meets one other end, so
        //- the linking is unambiguous, and a chain runs from where the loop
        //- crosses one cell edge to where it crosses another.
        List<pointField> chains;

        {
            labelList used(onThisFace.size(), 0);

            forAll(onThisFace, c0)
            {
                if (used[c0])
                {
                    continue;
                }

                DynamicList<point> forward(4);
                DynamicList<point> backward(4);

                forward.append(chords[onThisFace[c0]][0]);
                forward.append(chords[onThisFace[c0]][1]);
                used[c0] = 1;

                label fwdFacet = chordFacet[onThisFace[c0]];
                label bwdFacet = fwdFacet;

                bool grew = true;

                while (grew)
                {
                    grew = false;

                    forAll(onThisFace, cj)
                    {
                        if (used[cj])
                        {
                            continue;
                        }

                        //- Only a chord whose facet meets this one along an
                        //- edge ending here continues the same curve
                        if
                        (
                            degreeAt(forward.last()) > 2
                         && !adjacentAt
                            (
                                fwdFacet,
                                pointField({forward[forward.size() - 2], forward.last()}),
                                chordFacet[onThisFace[cj]],
                                chords[onThisFace[cj]],
                                forward.last()
                            )
                        )
                        {
                            continue;
                        }

                        if (mag(forward.last() - chords[onThisFace[cj]][0]) <= tol)
                        {
                            forward.append(chords[onThisFace[cj]][1]);
                        }
                        else if (mag(forward.last() - chords[onThisFace[cj]][1]) <= tol)
                        {
                            forward.append(chords[onThisFace[cj]][0]);
                        }
                        else
                        {
                            continue;
                        }

                        used[cj] = 1;
                        fwdFacet = chordFacet[onThisFace[cj]];
                        grew = true;
                        break;
                    }
                }

                grew = true;

                while (grew)
                {
                    grew = false;

                    forAll(onThisFace, cj)
                    {
                        if (used[cj])
                        {
                            continue;
                        }

                        //- The head moves as the chain grows backwards: it is
                        //- the last point appended, not the original first
                        const point head =
                            (backward.empty() ? forward.first() : backward.last());

                        if
                        (
                            degreeAt(head) > 2
                         && !adjacentAt
                            (
                                bwdFacet,
                                pointField
                                (
                                    backward.size() >= 2
                                  ? pointField({backward[backward.size() - 2], head})
                                  : pointField({forward[0], forward[1]})
                                ),
                                chordFacet[onThisFace[cj]],
                                chords[onThisFace[cj]],
                                head
                            )
                        )
                        {
                            continue;
                        }

                        if (mag(head - chords[onThisFace[cj]][1]) <= tol)
                        {
                            backward.append(chords[onThisFace[cj]][0]);
                        }
                        else if (mag(head - chords[onThisFace[cj]][0]) <= tol)
                        {
                            backward.append(chords[onThisFace[cj]][1]);
                        }
                        else
                        {
                            continue;
                        }

                        used[cj] = 1;
                        bwdFacet = chordFacet[onThisFace[cj]];
                        grew = true;
                        break;
                    }
                }

                pointField chain(backward.size() + forward.size());

                label n = 0;

                for (label k = backward.size() - 1; k >= 0; k--)
                {
                    chain[n++] = backward[k];
                }

                forAll(forward, k)
                {
                    chain[n++] = forward[k];
                }

                chains.append(chain);
            }
        }

        if (loopDebug)
        {
            //- How far apart are the ends that should have linked?
            forAll(onThisFace, ca)
            {
                forAll(onThisFace, cb)
                {
                    if (ca == cb) continue;

                    forAll(labelList({0, 1}), e1)
                    {
                        forAll(labelList({0, 1}), e2)
                        {
                            const point pa = chords[onThisFace[ca]][e1];
                            const point pb = chords[onThisFace[cb]][e2];

                            const scalar d = mag(pa - pb);

                            if (d > tol && d < 1e-5)
                            {
                                Info<< "NEAREND " << ffi << " " << d
                                    << " " << pa << " " << pb << endl;
                            }
                        }
                    }
                }
            }

            //- Pairs whose ends coincide and which the adjacency rejected:
            //- the chords a chain should have continued through
            forAll(onThisFace, ca)
            {
                forAll(onThisFace, cb)
                {
                    if (ca == cb)
                    {
                        continue;
                    }

                    forAll(labelList({0, 1}), e1)
                    {
                        forAll(labelList({0, 1}), e2)
                        {
                            const point Q = chords[onThisFace[ca]][e1];

                            if (mag(Q - chords[onThisFace[cb]][e2]) > tol)
                            {
                                continue;
                            }

                            if
                            (
                                !adjacentAt
                                (
                                    chordFacet[onThisFace[ca]],
                                    chords[onThisFace[ca]],
                                    chordFacet[onThisFace[cb]],
                                    chords[onThisFace[cb]],
                                    Q
                                )
                            )
                            {
                                Info<< "UNLINKED face " << ffi << " facets "
                                    << chordFacet[onThisFace[ca]] << " and "
                                    << chordFacet[onThisFace[cb]]
                                    << " at " << Q << endl;

                                //- What edges each facet actually has there
                                forAll(labelList({0, 1}), side)
                                {
                                    const label fi = (side == 0
                                      ? chordFacet[onThisFace[ca]]
                                      : chordFacet[onThisFace[cb]]);

                                    Info<< "   facet " << fi << " edges near Q:";

                                    forAll(facet[fi], k)
                                    {
                                        const point& e0 = facet[fi][k];
                                        const point& e1 = facet[fi][facet[fi].fcIndex(k)];

                                        if
                                        (
                                            mag(e0 - Q) <= 10*tol
                                         || mag(e1 - Q) <= 10*tol
                                        )
                                        {
                                            Info<< nl << "     " << e0
                                                << " -> " << e1;
                                        }
                                    }

                                    Info<< endl;
                                }
                            }
                        }
                    }
                }
            }

            forAll(chains, chi)
            {
                Info<< "CHAIN face " << ffi << " " << chi
                    << " size " << chains[chi].size()
                    << " from " << chains[chi].first()
                    << " to " << chains[chi].last() << endl;

                forAll(chains[chi], ck)
                {
                    Info<< "      pt " << ck << ' ' << chains[chi][ck] << endl;
                }

                Info<< "   face " << ffi << " regions "
                    << faceRegions[ffi].size() << endl;

                forAll(faceRegions[ffi], rk)
                {
                    Info<< "      region " << rk << " size "
                        << faceRegions[ffi][rk].size() << " :";

                    forAll(faceRegions[ffi][rk], pk)
                    {
                        Info<< ' ' << faceRegions[ffi][rk][pk];
                    }

                    Info<< endl;
                }
            }
        }

        forAll(chains, chi)
        {
            nChains++;

            bool done = false;
            bool progress = true;
            label nPasses = 0;

            while (progress && nPasses < 8)
            {
                progress = false;
                nPasses++;

                forAll(faceRegions[ffi], ri)
                {
                    pointField sub;

                    if (!chainInRegion(chains[chi], faceRegions[ffi][ri], sub))
                    {
                        continue;
                    }

                    pointField A;
                    pointField B;

                    if (splitByChain(faceRegions[ffi][ri], sub, A, B))
                    {
                        faceRegions[ffi][ri] = A;
                        faceRegions[ffi].append(B);
                        done = true;
                        progress = true;
                        break;
                    }
                }
            }

            if (!done)
            {
                nChainsSkipped++;

                //- A chain whose two ends are the same point closed on itself. If
                //- it split nothing, that face cannot be split by it.
                if
                (
                    chains[chi].size() >= 3
                 && mag(chains[chi].first() - chains[chi].last())
                    <= samePointTol(tol)
                )
                {
                    closedChainSplitNothing = true;
                }

                if (loopDebug)
                {
                    Info<< "   SKIPPED chain " << chi << " on face " << ffi
                        << " size " << chains[chi].size()
                        << " regions " << faceRegions[ffi].size() << endl;
                }
            }
            else if (loopDebug)
            {
                Info<< "   SPLIT chain " << chi << " on face " << ffi
                    << " -> regions " << faceRegions[ffi].size() << endl;
            }
        }
    }

    //- 3b. Subdivide every cell edge ONCE, and cut every region's edge by that
    //- one list. Refining each region from the points it happens to know was
    //- tried and does not work: two faces put different cuts on the same cell
    //- edge - measured, on one edge, at 0.62632, 0.668335, 0.700654 on one face
    //- and 0.640531, 0.675121, 0.700654 on the other - and only the one point
    //- they share can be stitched. Taking the union of everything that lands on
    //- the edge and imposing it on both sides makes them agree by construction,
    //- and the extra points become vertices of the mesh's faces, which is what
    //- a shared edge requires of them.
    List<Pair<label>> cellEdges;

    forAll(cell.faces(), ffi)
    {
        const face& f = cell.faces()[ffi];

        forAll(f, i)
        {
            const label a = f[i];
            const label b = f[f.fcIndex(i)];

            const Pair<label> e(min(a, b), max(a, b));

            bool known = false;

            forAll(cellEdges, k)
            {
                if (cellEdges[k] == e)
                {
                    known = true;
                    break;
                }
            }

            if (!known)
            {
                cellEdges.append(e);
            }
        }
    }

    //- The cuts on each cell edge, as a parameter along it
    List<DynamicList<scalar>> edgeCuts(cellEdges.size());

    auto edgeOn = [&](const point& a, const point& b) -> label
    {
        forAll(cellEdges, ei)
        {
            const point& e0 = cPts[cellEdges[ei].first()];
            const vector ed(cPts[cellEdges[ei].second()] - e0);
            const scalar len2 = ed & ed;

            if (len2 <= VSMALL)
            {
                continue;
            }

            const scalar ta = ((a - e0) & ed)/len2;
            const scalar tb = ((b - e0) & ed)/len2;

            if
            (
                ta >= -SMALL && ta <= 1 + SMALL
             && tb >= -SMALL && tb <= 1 + SMALL
             && mag(a - (e0 + ta*ed)) <= 10*tol
             && mag(b - (e0 + tb*ed)) <= 10*tol
            )
            {
                return ei;
            }
        }

        return -1;
    };

    auto addCut = [&](const label ei, const scalar t)
    {
        if (ei < 0 || t <= SMALL || t >= 1 - SMALL)
        {
            return;
        }

        forAll(edgeCuts[ei], k)
        {
            if (mag(edgeCuts[ei][k] - t) <= 1e-9)
            {
                return;
            }
        }

        edgeCuts[ei].append(t);
    };

    //- Every chord end that lands on a cell edge
    forAll(chords, ci)
    {
        const label ei = edgeOn(chords[ci][0], chords[ci][1]);

        if (ei < 0)
        {
            continue;
        }

        const point& e0 = cPts[cellEdges[ei].first()];
        const vector ed(cPts[cellEdges[ei].second()] - e0);
        const scalar len2 = ed & ed;

        addCut(ei, ((chords[ci][0] - e0) & ed)/len2);
        addCut(ei, ((chords[ci][1] - e0) & ed)/len2);
    }

    //- And every region vertex that lands on one, so a cut a region already
    //- carries is imposed on the face opposite it
    forAll(faceRegions, ffi)
    {
        forAll(faceRegions[ffi], ri)
        {
            const pointField& r = faceRegions[ffi][ri];

            forAll(r, i)
            {
                forAll(cellEdges, ei)
                {
                    const point& e0 = cPts[cellEdges[ei].first()];
                    const vector ed(cPts[cellEdges[ei].second()] - e0);
                    const scalar len2 = ed & ed;

                    if (len2 <= VSMALL)
                    {
                        continue;
                    }

                    const scalar t = ((r[i] - e0) & ed)/len2;

                    if
                    (
                        t >= -SMALL && t <= 1 + SMALL
                     && mag(r[i] - (e0 + t*ed)) <= 10*tol
                    )
                    {
                        addCut(ei, t);
                    }
                }
            }
        }
    }

    //- Cut every region's edge that lies on a cell edge by that edge's list
    forAll(faceRegions, ffi)
    {
        forAll(faceRegions[ffi], ri)
        {
            const pointField r(faceRegions[ffi][ri]);

            DynamicList<point> refined(2*r.size());

            forAll(r, i)
            {
                const point& a = r[i];
                const point& b = r[r.fcIndex(i)];

                const vector ab(b - a);
                const scalar len2 = ab & ab;

                refined.append(a);

                if (len2 <= VSMALL)
                {
                    continue;
                }

                const label ei = edgeOn(a, b);

                if (ei < 0)
                {
                    continue;
                }

                const scalar ta = 0;
                const scalar tb = 1;

                //- The cell edge's own parameterisation, not the region edge's,
                //- so both faces place the same cut in the same place
                const point& e0 = cPts[cellEdges[ei].first()];
                const vector ed(cPts[cellEdges[ei].second()] - e0);
                const scalar elen2 = ed & ed;

                const scalar t0 = ((a - e0) & ed)/elen2;
                const scalar t1 = ((b - e0) & ed)/elen2;

                List<Tuple2<scalar, point>> inside(edgeCuts[ei].size());
                label n = 0;

                forAll(edgeCuts[ei], k)
                {
                    const scalar t = edgeCuts[ei][k];

                    const bool between =
                        (t0 < t1 ? (t > t0 && t < t1) : (t > t1 && t < t0));

                    if (between)
                    {
                        inside[n++] = Tuple2<scalar, point>
                        (
                            mag(t - t0),
                            e0 + t*ed
                        );
                    }
                }

                inside.setSize(n);

                std::sort
                (
                    inside.begin(),
                    inside.end(),
                    [](const Tuple2<scalar, point>& x, const Tuple2<scalar, point>& y)
                    {
                        return x.first() < y.first();
                    }
                );

                forAll(inside, k)
                {
                    refined.append(inside[k].second());
                }

                (void)ta;
                (void)tb;
            }

            faceRegions[ffi][ri].transfer(refined);
        }
    }

    //- Every chord, and which cell face the splitting thinks it lies on. A
    //- segment can only be on one face unless it lies along a cell edge, so a
    //- chord claiming two faces is not a crossing of the loop but a segment
    //- running along an edge, and it splits both faces when it should not.
    if (loopDebug)
    {
        forAll(chords, ci)
        {
            Info<< "CHORD " << ci << " " << chords[ci][0]
                << " -> " << chords[ci][1] << " on faces";

            forAll(cell.faces(), ffi)
            {
                if (onFace(cell.faces()[ffi], cPts, 0.5*(chords[ci][0] + chords[ci][1]), tol))
                {
                    Info<< ' ' << ffi;
                }
            }

            Info<< " len " << mag(chords[ci][1] - chords[ci][0]) << endl;
        }
    }

    //- Is the loop a set of clean closed curves? Every chord end should be met
    //- by exactly two chords, so that the chords pair up into loops that
    //- separate the cell's boundary. An end met by one or by three means the
    //- curve touches or crosses itself, and the regions either side of a chord
    //- cannot then be told apart.
    if (loopDebug)
    {
        List<point> ends;
        labelList endCount;
        label nOne = 0;
        label nTwo = 0;
        label nMore = 0;

        forAll(chords, ci)
        {
            forAll(chords[ci], k)
            {
                label found = -1;

                forAll(ends, e)
                {
                    if (mag(chords[ci][k] - ends[e]) <= tol)
                    {
                        found = e;
                        break;
                    }
                }

                if (found < 0)
                {
                    ends.append(chords[ci][k]);
                    endCount.append(0);
                    found = ends.size() - 1;
                }

                endCount[found]++;
            }
        }

        forAll(ends, e)
        {
            if (endCount[e] == 1) nOne++;
            else if (endCount[e] == 2) nTwo++;
            else nMore++;
        }

        Info<< "LOOP curve: chords " << chords.size()
            << " ends " << ends.size()
            << " met by 1 " << nOne
            << ", by 2 " << nTwo
            << ", by 3 or more " << nMore << endl;
    }

    //- Every region, so their structure can be seen at once
    if (loopDebug)
    {
        label n = 0;

        forAll(faceRegions, ffi)
        {
            forAll(faceRegions[ffi], ri)
            {
                const pointField& r = faceRegions[ffi][ri];

                point rc = Zero;

                forAll(r, k)
                {
                    rc += r[k];
                }

                rc /= scalar(r.size());

                scalar zlo = GREAT;
                scalar zhi = -GREAT;

                forAll(r, k)
                {
                    zlo = min(zlo, r[k].z());
                    zhi = max(zhi, r[k].z());
                }

                Info<< "REG " << n++ << " face " << ffi
                    << " size " << r.size()
                    << " area " << polygonArea(r)
                    << " z [" << zlo << " " << zhi << "]"
                    << " centre " << rc << endl;
            }
        }
    }

    //- The report, when asked for: for each edge of the cell, the cuts each of
    //- its two faces contributes. Regions cannot stitch across an edge their
    //- faces cut differently.
    if (loopDebug)
    {
        label nRegions = 0;

        forAll(faceRegions, ffi)
        {
            nRegions += faceRegions[ffi].size();
        }

        Info<< "LOOP report: cell edges " << cellEdges.size()
            << " regions " << nRegions << endl;

        //- Do the regions tile each face? Everything the assembly does with
        //- them assumes it, and a point found outside every piece would be
        //- evidence of a gap rather than of a tolerance.
        forAll(cell.faces(), ffi)
        {
            const face& f = cell.faces()[ffi];

            pointField fp(f.size());

            forAll(f, i)
            {
                fp[i] = cPts[f[i]];
            }

            scalar sum = 0;

            forAll(faceRegions[ffi], ri)
            {
                sum += polygonArea(faceRegions[ffi][ri]);
            }

            const scalar faceArea = polygonArea(fp);

            if (mag(sum - faceArea) > 1e-6*max(faceArea, SMALL))
            {
                Info<< "GAP face " << ffi << " regions " << faceRegions[ffi].size()
                    << " area " << sum << " of " << faceArea
                    << " difference " << (sum - faceArea) << endl;
            }
        }

        forAll(cellEdges, ei)
        {
            const point& e0 = cPts[cellEdges[ei].first()];
            const vector ed(cPts[cellEdges[ei].second()] - e0);
            const scalar len2 = ed & ed;

            forAll(cell.faces(), ffi)
            {
                DynamicList<scalar> ts;

                forAll(faceRegions[ffi], ri)
                {
                    const pointField& r = faceRegions[ffi][ri];

                    forAll(r, k)
                    {
                        const point& a = r[k];
                        const point& b = r[r.fcIndex(k)];

                        const scalar ta = ((a - e0) & ed)/len2;
                        const scalar tb = ((b - e0) & ed)/len2;

                        if
                        (
                            ta >= -SMALL && ta <= 1 + SMALL
                         && tb >= -SMALL && tb <= 1 + SMALL
                         && mag(a - (e0 + ta*ed)) <= 10*tol
                         && mag(b - (e0 + tb*ed)) <= 10*tol
                        )
                        {
                            bool seenA = false;
                            bool seenB = false;

                            forAll(ts, m)
                            {
                                if (mag(ts[m] - ta) <= 1e-9) seenA = true;
                                if (mag(ts[m] - tb) <= 1e-9) seenB = true;
                            }

                            if (!seenA) ts.append(ta);
                            if (!seenB) ts.append(tb);
                        }
                    }
                }

                if (ts.size() < 2)
                {
                    continue;
                }

                std::sort(ts.begin(), ts.end());

                Info<< "  edge " << ei << " face " << ffi << " cuts";

                forAll(ts, k)
                {
                    Info<< ' ' << ts[k];
                }

                Info<< endl;
            }
        }
    }

    //- 4. Group the regions into the pieces, by sharing an edge that lies on a
    //- CELL edge. Two regions either side of a chord lie on the same face and
    //- share that chord, but the chord is the loop: they are on opposite sides of
    //- it and belong to different pieces, so a join across one is never made.
    List<pointField> allRegions;

    forAll(faceRegions, ffi)
    {
        forAll(faceRegions[ffi], ri)
        {
            allRegions.append(faceRegions[ffi][ri]);
        }
    }

    labelList group(allRegions.size());

    forAll(group, i)
    {
        group[i] = i;
    }

    label nOnEdge = 0;
    static label nShown = 0;
    label nJoins = 0;
    label nCandidates = 0;
    label nUnmatched = 0;

    auto root = [&](label i) -> label
    {
        while (group[i] != i)
        {
            group[i] = group[group[i]];
            i = group[i];
        }

        return i;
    };

    forAll(allRegions, i)
    {
        const pointField& a = allRegions[i];

        forAll(a, ai)
        {
            const point a0 = a[ai];
            const point a1 = a[a.fcIndex(ai)];

            const label edgei = edgeOn(a0, a1);

            if (edgei < 0)
            {
                continue;
            }

            nOnEdge++;

            static label nShownEdge = 0;

            if (nShownEdge < 2)
            {
                nShownEdge++;

                const point mid(0.5*(a0 + a1));

                Info<< "LOOP onEdge " << i << " " << a0 << " -> " << a1 << nl;

                forAll(allRegions, ri)
                {
                    const pointField& rr = allRegions[ri];

                    forAll(rr, e)
                    {
                        const point r0 = rr[e];
                        const point r1 = rr[rr.fcIndex(e)];

                        //- Collinear with the outer edge and overlapping it
                        const vector d1(a1 - a0);
                        const vector d2(r1 - r0);

                        //- Parallel, by the cosine rather than a cross product
                        if
                        (
                            mag(d1) > VSMALL && mag(d2) > VSMALL
                         && mag((d1 & d2)/(mag(d1)*mag(d2))) >= 1 - 1e-9
                        )
                        {
                            Info<< "    region " << ri << " "
                                << r0 << " -> " << r1 << nl;
                        }
                    }
                }
            }

            for (label j = i + 1; j < allRegions.size(); j++)
            {
                const pointField& b = allRegions[j];

                bool joined = false;

                point bb0 = Zero;
                point bb1 = Zero;

                forAll(b, bi)
                {
                    const point b0 = b[bi];
                    const point b1 = b[b.fcIndex(bi)];

                    bb0 = b0;
                    bb1 = b1;

                    //- Overlap along the cell edge the outer edge lies on, not
                    //- equality of ends. The cuts are canonical now, so both
                    //- regions are unions of the same sub-edges, but not of the
                    //- same ones: a region spanning 0.62632 to 0.700654 faces
                    //- two that split it, and all three are one piece. Requiring
                    //- the ends to match left six of this cell's twenty-four
                    //- regions alone in a group of their own.
                    if (edgei < 0 || edgeOn(b0, b1) != edgei)
                    {
                        continue;
                    }

                    const point& ce0 = cPts[cellEdges[edgei].first()];
                    const vector ced(cPts[cellEdges[edgei].second()] - ce0);
                    const scalar clen2 = ced & ced;

                    const scalar at0 = ((a0 - ce0) & ced)/clen2;
                    const scalar at1 = ((a1 - ce0) & ced)/clen2;
                    const scalar bt0 = ((b0 - ce0) & ced)/clen2;
                    const scalar bt1 = ((b1 - ce0) & ced)/clen2;

                    const scalar lo = max(min(at0, at1), min(bt0, bt1));
                    const scalar hi = min(max(at0, at1), max(bt0, bt1));

                    if (hi - lo > 1e-14)
                    {
                        joined = true;
                        break;
                    }
                }

                nCandidates++;

                if (joined)
                {
                    nJoins++;

                    static label nJoinPrint = 0;

                    if (loopDebug && nJoinPrint < 40)
                    {
                        nJoinPrint++;

                        Info<< "JOIN " << i << " " << j
                            << " on edge " << edgei
                            << " a " << a0 << " -> " << a1
                            << " b " << bb0 << " -> " << bb1 << endl;
                    }

                    const label ra = root(i);
                    const label rb = root(j);

                    if (ra != rb)
                    {
                        group[rb] = ra;
                    }
                }
                else
                {
                    nUnmatched++;

                    if (nUnmatched == 1 && nShown < 2)
                    {
                        nShown++;
                    }
                }
            }
        }
    }

    labelList pieceOf(allRegions.size(), -1);

    DynamicList<label> roots;

    forAll(group, i)
    {
        const label r = root(i);

        if (pieceOf[r] < 0)
        {
            pieceOf[r] = roots.size();
            roots.append(r);
        }

        pieceOf[i] = pieceOf[r];
    }

    const label nPieces = roots.size();

    {
        static label nReported = 0;

        if (nReported < 4)
        {
            nReported++;

            label nCuts = 0;

            forAll(edgeCuts, ei)
            {
                nCuts += edgeCuts[ei].size();
            }

            Info<< "LOOP edgeCuts " << nCuts
                << " chords " << chords.size() << endl;

            Info<< "LOOP groups:";

            forAll(pieceOf, i)
            {
                Info<< ' ' << i << "->" << pieceOf[i];
            }

            Info<< endl;

            //- For a region alone in its group, which cell edge its edges on
            //- cell edges are attributed to, and who else is on those edges
            labelList groupSize(pieceOf.size(), 0);

            forAll(pieceOf, i)
            {
                groupSize[pieceOf[i]]++;
            }

            label nShownIso = 0;

            forAll(pieceOf, i)
            {
                if (groupSize[pieceOf[i]] != 1 || nShownIso >= 2)
                {
                    continue;
                }

                nShownIso++;

                point rc = Zero;

                forAll(allRegions[i], k)
                {
                    rc += allRegions[i][k];
                }

                rc /= scalar(allRegions[i].size());

                Info<< "ISO region " << i << " group " << pieceOf[i]
                    << " size " << allRegions[i].size()
                    << " area " << polygonArea(allRegions[i])
                    << " centre " << rc << nl;

                forAll(allRegions[i], e)
                {
                    const point& r0 = allRegions[i][e];
                    const point& r1 = allRegions[i][allRegions[i].fcIndex(e)];

                    forAll(cellEdges, ei)
                    {
                        const point& e0 = cPts[cellEdges[ei].first()];
                        const vector ed(cPts[cellEdges[ei].second()] - e0);
                        const scalar len2 = ed & ed;

                        const scalar ta = ((r0 - e0) & ed)/len2;
                        const scalar tb = ((r1 - e0) & ed)/len2;

                        if
                        (
                            ta >= -SMALL && ta <= 1 + SMALL
                         && tb >= -SMALL && tb <= 1 + SMALL
                         && mag(r0 - (e0 + ta*ed)) <= 10*tol
                         && mag(r1 - (e0 + tb*ed)) <= 10*tol
                         && mag(r1 - r0) > VSMALL
                        )
                        {
                            Info<< "   edge " << ei << " [" << ta << " "
                                << tb << "]";

                            //- Who else has an edge on this cell edge
                            label nOther = 0;
                            labelList others;

                            forAll(allRegions, j)
                            {
                                if (j == i)
                                {
                                    continue;
                                }

                                forAll(allRegions[j], f)
                                {
                                    const point& s0 = allRegions[j][f];
                                    const point& s1 = allRegions[j][allRegions[j].fcIndex(f)];

                                    const scalar ua = ((s0 - e0) & ed)/len2;
                                    const scalar ub = ((s1 - e0) & ed)/len2;

                                    if
                                    (
                                        ua >= -SMALL && ua <= 1 + SMALL
                                     && ub >= -SMALL && ub <= 1 + SMALL
                                     && mag(s0 - (e0 + ua*ed)) <= 10*tol
                                     && mag(s1 - (e0 + ub*ed)) <= 10*tol
                                     && mag(s1 - s0) > VSMALL
                                    )
                                    {
                                        if (findIndex(others, j) < 0)
                                        {
                                            others.append(j);
                                            nOther++;

                                            Info<< " r" << j << "[" << ua << " "
                                                << ub << "]";
                                        }
                                    }
                                }
                            }

                            Info<< " others " << nOther << nl;
                        }
                    }
                }
            }

            Info<< "LOOP chains " << nChains
                << " skipped " << nChainsSkipped << endl;

            Info<< "LOOP cell: regions " << allRegions.size()
                << " candidates " << nCandidates
                << " joins " << nJoins
                << " unmatched " << nUnmatched
                << " edgesOnCellEdge " << nOnEdge
                << " pieces " << nPieces << endl;
        }
    }

    if (nPieces < 2)
    {
        return cutState::untouched;
    }

    //- 5. Build the pieces from their regions alone. The interface faces are
    //- found afterwards, by looking at the piece on each side of every part of
    //- every facet plane, which is what an interface is. Building them by
    //- clipping a plane's section by the other planes' inside half-spaces was
    //- tried and is wrong: it gives one side of the plane rather than both, so
    //- a cell whose surface is two crossing sheets came out with one interface
    //- where it has two.
    List<labelList> pieceRegions(nPieces);

    forAll(pieceOf, i)
    {
        pieceRegions[pieceOf[i]].append(i);
    }

    result.pieces.setSize(nPieces);
    result.sharedFaces.clear();
    result.planeTri = planeTri;

    forAll(result.pieces, piecei)
    {
        List<pointField> polys(pieceRegions[piecei].size());

        forAll(pieceRegions[piecei], i)
        {
            polys[i] = allRegions[pieceRegions[piecei][i]];
        }

        //- A point inside the piece, from the mean of its regions' centres. A
        //- piece bounds the cell on the regions' side and the surface on the
        //- other, so this lies on the piece's surface, which is enough to
        //- orient its faces.
        point inside = Zero;
        label nR = 0;

        forAll(pieceRegions[piecei], i)
        {
            const pointField& r = allRegions[pieceRegions[piecei][i]];

            point c = Zero;

            forAll(r, k)
            {
                c += r[k];
            }

            inside += c/scalar(r.size());
            nR++;
        }

        inside /= scalar(max(nR, 1));

        result.pieces[piecei].poly = polyhedronFrom(polys, inside, tol);
        result.pieces[piecei].side.setSize(planes.size(), 0);
    }

    //- 6. The interface faces. Every part of a facet plane inside the cell
    //- separates the piece on one side of it from the piece on the other, so the
    //- parts are enumerated by splitting the plane's section by the other planes
    //- and then asking which piece contains a point just off each side.
    //- Do two collinear segments share more than a point?
    //-
    //- The meeting test below used to compare whole-edge MIDPOINTS, which holds
    //- only when the two edges are the same segment. It fails whenever one is a
    //- sub-segment of the other, and that is the normal case: the conformality
    //- step subdivides the regions' edges while a part's edge is not subdivided
    //- the same way. Measured on 09_pipe_cut cell 493, the cap plane's parts
    //- were 0.011 to 0.055 from the nearest region edge midpoint against a tol
    //- of 1.7e-9, so every edge missed, sides came out empty, and the whole cap
    //- interface was discarded - 0.0986 of the pipe's 2.0627 surface, all of it
    //- in eight cells. Step 4's join test was given this same treatment
    //- earlier, for the same reason.
    auto edgesOverlap = [&](const point& a0, const point& a1,
                            const point& b0, const point& b1) -> bool
    {
        const vector da(a1 - a0);

        if ((da & da) <= VSMALL)
        {
            return false;
        }

        const vector db(b1 - b0);

        if ((db & db) <= VSMALL)
        {
            return false;
        }

        //- Collinear: each end of one lies on the line of the other
        if
        (
            mag(da ^ (b0 - a0)) > 10*tol*mag(da)
         || mag(da ^ (b1 - a0)) > 10*tol*mag(da)
         || mag(db ^ (a0 - b0)) > 10*tol*mag(db)
         || mag(db ^ (a1 - b0)) > 10*tol*mag(db)
        )
        {
            return false;
        }

        //- And they overlap along a by more than the tolerance
        const scalar t0 = ((b0 - a0) & da)/(da & da);
        const scalar t1 = ((b1 - a0) & da)/(da & da);

        const scalar lo = max(min(t0, t1), 0.0);
        const scalar hi = min(max(t0, t1), 1.0);

        return (hi - lo)*mag(da) > 10*tol;
    };

    //- The pieces this cell was grouped into. Two is the generic answer - one
    //- cap crossing a cell - and it is the fallback pair for a patch part whose
    //- own boundary reaches no region at all; see the decision in the parts loop.
    labelList cellPieces;

    forAll(pieceOf, i)
    {
        if (findIndex(cellPieces, pieceOf[i]) < 0)
        {
            cellPieces.append(pieceOf[i]);
        }
    }

    //- Probe: the cell's own face planes, built only for the dumped cell. The
    //- regions are pieces of these faces, so every region edge lies in one of
    //- them - which is checked below rather than assumed - and a part's boundary
    //- can then be tested against the same planes the regions live in. Added
    //- 26-Sep-2026 to settle why one plane's patch is dropped where nine others
    //- are accepted, per 06 B.3 item (0).
    DynamicList<plane> facePlanes;

    if (dump)
    {
        const pointField& cps = cell.points();

        forAll(cell.faces(), fi)
        {
            const face& f = cell.faces()[fi];

            facePlanes.append(plane(cps[f[0]], cps[f[1]], cps[f[2]]));
        }

        //- A number that must be zero: regions carrying a vertex off the plane
        //- of the cell face they are a piece of. If it is not zero then the
        //- regions are not pieces of those faces, and every conclusion drawn
        //- from that - including the whole of this probe - has to be re-read.
        label nOffPlane = 0;

        forAll(faceRegions, ffi)
        {
            forAll(faceRegions[ffi], ri)
            {
                forAll(faceRegions[ffi][ri], pk)
                {
                    if
                    (
                        mag
                        (
                            facePlanes[ffi].signedDistance
                            (
                                faceRegions[ffi][ri][pk]
                            )
                        ) > 10*tol
                    )
                    {
                        nOffPlane++;
                        break;
                    }
                }
            }
        }

        Info<< "REGIONPLANE cellfaces " << cell.nFaces()
            << " planes " << facePlanes.size()
            << " regionsOffTheirPlane " << nOffPlane
            << " pieces " << cellPieces.size() << " " << cellPieces
            << nl << "            pieceOf " << pieceOf << nl;
    }

    //- Which pieces meet over a part: the ones owning the regions along its
    //- boundary edges. Every edge of the loop lies on the cell's boundary, so it
    //- is shared with a region of exactly one piece, and a part of the surface
    //- separates two pieces. That holds for the loop; a facet wholly inside the
    //- cell has no edge on that boundary, which is why the caller keeps a
    //- fallback pair (26-Sep-2026, 06 B.3 item (0)).
    //-
    //- Asking the pieces directly with a containment test was tried and does not
    //- work: a piece is built from its regions alone at this point, so it is
    //- open, and a point outside it tests as inside.
    auto piecesOf = [&](const pointField& part) -> labelList
    {
        labelList sides;

        forAll(part, m)
        {
            const point& pa = part[m];
            const point& pb = part[part.fcIndex(m)];

            forAll(allRegions, ri)
            {
                const pointField& r = allRegions[ri];

                bool matches = false;

                forAll(r, e)
                {
                    if (edgesOverlap(pa, pb, r[e], r[r.fcIndex(e)]))
                    {
                        matches = true;
                        break;
                    }
                }

                if (matches)
                {
                    const label gp = pieceOf[ri];

                    if (findIndex(sides, gp) < 0)
                    {
                        sides.append(gp);
                    }
                }
            }
        }

        return sides;
    };

    List<pointField> ifacePoly;
    List<label> ifacePlane;

    //- Probe: what the enumeration discards, how much of it, and where the pair
    //- for a part whose boundary reached no region came from.
    scalar droppedArea = 0;
    label nPairedFromPatch = 0;
    label nPairedFromCell = 0;
    label nPairless = 0;
    label nDropped = 0;
    label nDropped3 = 0;
    label nDropped1 = 0;
    label nDropped0 = 0;


    forAll(planes, pi)
    {
        const pointField sec(planeSection(cell, planes[pi], tol));

        if (sec.size() < 3)
        {
            continue;
        }

        //- The surface actually in this plane, as the facets of the plane
        //- clipped to the cell. The parts below are pieces of the plane's WHOLE
        //- cross-section, and that is the surface only where a facet fills the
        //- cross-section. At a corner it is not: the cross-section runs on past
        //- the surface, and the parts beyond it still have a piece on each side,
        //- because both pieces exist there - they are merely the same two pieces
        //- - so the meeting test passes them and they are counted as interface.
        scalar patchArea = 0;

        forAll(planeFacets[pi], fj)
        {
            patchArea += polygonArea(facet[planeFacets[pi][fj]]);
        }

        //- The part of a polygon that lies in the patch, as one polygon per
        //- facet it meets. Both shapes are convex - a facet is a clipped
        //- triangle, a part is a section of a hex cut by planes - so each of a
        //- facet's edges is a half-space, `clipPolygon` does the work, and the
        //- facet's own centre says which side is inside.
        auto clipToPatch = [&](const pointField& poly) -> DynamicList<pointField>
        {
            DynamicList<pointField> pieces;

            forAll(planeFacets[pi], fj)
            {
                const pointField& f = facet[planeFacets[pi][fj]];

                if (f.size() < 3)
                {
                    continue;
                }

                point fc = Zero;

                forAll(f, k)
                {
                    fc += f[k];
                }

                fc /= scalar(f.size());

                pointField out(poly);

                forAll(f, i)
                {
                    const point& a = f[i];
                    const vector n = (f[f.fcIndex(i)] - a) ^ planes[pi].normal();

                    //- The facet's centre must be on the side that is KEPT, and
                    //- `clipPolygon` keeps the negative one, so the plane's
                    //- normal points away from it.
                    out = clipPolygon(out, plane(a, ((fc - a) & n) <= 0 ? n : -n), tol);

                    if (out.size() < 3)
                    {
                        break;
                    }
                }

                if (out.size() >= 3 && polygonArea(out) > tol*tol)
                {
                    pieces.append(out);
                }
            }

            return pieces;
        };

        DynamicList<pointField> parts(4);
        parts.append(sec);

        forAll(planes, pj)
        {
            if (pj == pi)
            {
                continue;
            }

            DynamicList<pointField> next(parts.size() + 2);

            forAll(parts, k)
            {
                const pointField& part = parts[k];

                scalar minSide = GREAT;
                scalar maxSide = -GREAT;

                forAll(part, m)
                {
                    const scalar d = planes[pj].signedDistance(part[m]);

                    minSide = min(minSide, d);
                    maxSide = max(maxSide, d);
                }

                if (minSide > -tol || maxSide < tol)
                {
                    next.append(part);
                    continue;
                }

                const pointField neg(clipPolygon(part, planes[pj], tol));

                const plane other(planes[pj].refPoint(), -planes[pj].normal());

                const pointField pos(clipPolygon(part, other, tol));

                if (neg.size() >= 3 && pos.size() >= 3)
                {
                    //- A part narrower than the tolerance is the arrangement's
                    //- own degeneracy appearing in two dimensions. It is
                    //- dropped rather than made into a face of no area.
                    const scalar w =
                        max(mag(planes[pj].signedDistance(neg[0])),
                            mag(planes[pj].signedDistance(pos[0])));

                    if (minSide < -100*tol && maxSide > 100*tol)
                    {
                        next.append(neg);
                        next.append(pos);
                        continue;
                    }

                    (void)w;
                }

                next.append(part);
            }

            parts.transfer(next);
        }

        //- The pair this plane's patch separates, taken from any part of it
        //- whose boundary does reach the regions, and every part's own sides -
        //- matched once here so that the decision below reads rather than
        //- re-matches, and so that a pair found later in the plane still counts
        //- for a part earlier in it.
        List<labelList> sidesOf(parts.size());
        labelList patchPair;

        forAll(parts, k)
        {
            if (parts[k].size() < 3)
            {
                continue;
            }

            sidesOf[k] = piecesOf(parts[k]);

            if (patchPair.empty() && sidesOf[k].size() == 2)
            {
                patchPair = sidesOf[k];
            }
        }

        forAll(parts, k)
        {
            const pointField& part = parts[k];

            if (part.size() < 3)
            {
                continue;
            }

            point c = Zero;

            forAll(part, m)
            {
                c += part[m];
            }

            c /= scalar(part.size());


            //- This part's own sides, matched once for the whole plane above.
            const labelList& sides = sidesOf[k];

            //- Probe: what each interface candidate found
            {
                static label nPartPrint = 0;

                if (loopDebug && nPartPrint < 30)
                {
                    nPartPrint++;

                    Info<< "PART pi " << pi << " k " << k
                        << " size " << part.size()
                        << " sides " << sides.size()
                        << " nRegions " << allRegions.size()
                        << " nPieces " << pieceOf.size()
                        << " area " << polygonArea(part)
                        << " secArea " << polygonArea(sec)
                        << " patchArea " << patchArea
                        << " nFacets " << planeFacets[pi].size() << endl;
                }
            }

            //- The dump of this cell's parts, one line each, with the polygon
            //- when nothing was found to separate. Printed for the one cell the
            //- mesher was asked to dump, and for no other, so a dropped part can
            //- be read whole without loopDebug's budget cutting it off.
            if (dump)
            {
                //- How much surface this part would carry if it were accepted:
                //- the area of the part that lies in the patch. For a part the
                //- edge test drops, it is the surface the drop costs, so the
                //- dropped parts' clipped areas must sum to the cell's IFACE
                //- `lost` - and if they do not, the loss is in the clip and not
                //- in the drop decision.
                scalar clippedArea = 0;

                {
                    const DynamicList<pointField> cp(clipToPatch(part));

                    forAll(cp, m)
                    {
                        clippedArea += polygonArea(cp[m]);
                    }
                }

                Info<< "PARTDUMP pi " << pi << " k " << k
                    << " size " << part.size()
                    << " sides " << sides.size()
                    << " area " << polygonArea(part)
                    << " secArea " << polygonArea(sec)
                    << " patchArea " << patchArea
                    << " clipped " << clippedArea
                    << " nFacets " << planeFacets[pi].size()
                    << " nRegions " << allRegions.size()
                    << " nPieces " << pieceOf.size() << nl;

                Info<< "         part " << part << nl;

                //- Probe: how this part's boundary sits relative to the planes
                //- the regions live in. `touch` counts its vertices lying on a
                //- cell-face plane, `onFace` its edges lying in one, and `match`
                //- the edges the meeting test reached. A dropped part with
                //- onFace 0 is one the edge test *could not* have matched, and
                //- one with onFace greater than 0 is a tolerance or predicate
                //- defect instead - which is the whole question.
                label touch = 0;
                label onFace = 0;
                label match = 0;

                forAll(part, m)
                {
                    forAll(facePlanes, fi)
                    {
                        if
                        (
                            mag(facePlanes[fi].signedDistance(part[m]))
                          < 10*tol
                        )
                        {
                            touch++;
                            break;
                        }
                    }
                }

                forAll(part, m)
                {
                    const point& pa = part[m];
                    const point& pb = part[part.fcIndex(m)];

                    forAll(facePlanes, fi)
                    {
                        if
                        (
                            mag(facePlanes[fi].signedDistance(pa)) < 10*tol
                         && mag(facePlanes[fi].signedDistance(pb)) < 10*tol
                        )
                        {
                            onFace++;
                            break;
                        }
                    }

                    forAll(allRegions, ri)
                    {
                        const pointField& r = allRegions[ri];

                        bool hit = false;

                        forAll(r, e)
                        {
                            if (edgesOverlap(pa, pb, r[e], r[r.fcIndex(e)]))
                            {
                                hit = true;
                                break;
                            }
                        }

                        if (hit)
                        {
                            match++;
                            break;
                        }
                    }
                }

                Info<< "EDGES onFace " << onFace
                    << " of " << part.size()
                    << " touch " << touch
                    << " match " << match
                    << " sides " << sides.size() << " " << sides << nl;
            }

            //- The interface is the surface, not the plane. The part is a piece of
            //- the plane's WHOLE cross-section, so clip it to the patch: that is
            //- the identity wherever the part is already inside - every cell of
            //- every case that was already exact, so nothing in those can move -
            //- and a cut at a corner, which is where the surplus is.
            const DynamicList<pointField> patchPiece(clipToPatch(part));

            scalar patchPieceArea = 0;

            forAll(patchPiece, m)
            {
                patchPieceArea += polygonArea(patchPiece[m]);
            }

            //- The pair this part separates. Its own boundary is the primary
            //- source; a part of the patch whose boundary reaches no region has
            //- none, and then the pair is the patch's own - every accepted part
            //- of a cell was measured naming the same pair (26-Sep-2026, cell
            //- 2952: nine of nine are (0 1)) - or the cell's two pieces when the
            //- whole patch found nothing. Only a part that carries surface is
            //- paired this way: without that, it lies outside the patch and its
            //- drop below is unchanged.
            labelList pair(sides);

            if (pair.size() != 2 && patchPieceArea > 0)
            {
                if (!patchPair.empty())
                {
                    pair = patchPair;
                    nPairedFromPatch++;
                }
                else if (cellPieces.size() == 2)
                {
                    pair = cellPieces;
                    nPairedFromCell++;
                }
            }

            if (pair.size() != 2)
            {
                //- Counted, so that the part the policy still cannot place is a
                //- number in the log and not a silence.
                if (patchPieceArea > 0)
                {
                    nPairless++;
                }

                droppedArea += polygonArea(part);
                nDropped++;

                if (sides.size() > 2) nDropped3++;
                else if (sides.size() == 1) nDropped1++;
                else nDropped0++;

                //- Probe: how far each of this part's edge midpoints is from the
                //- nearest edge midpoint of any region.
                {
                    static label nFar = 0;

                    if ((loopDebug && nFar < 4) || dump)
                    {
                        nFar++;

                        forAll(part, m)
                        {
                            const point mm(0.5*(part[m] + part[part.fcIndex(m)]));
                            scalar best = GREAT;

                            forAll(allRegions, ri)
                            {
                                const pointField& r = allRegions[ri];

                                forAll(r, e)
                                {
                                    best = min
                                    (
                                        best,
                                        mag(mm - 0.5*(r[e] + r[r.fcIndex(e)]))
                                    );
                                }
                            }

                            Info<< "   FAR pi " << pi << " edge " << m
                                << " mid " << mm
                                << " nearestRegEdge " << best << endl;
                        }
                    }
                }

                continue;
            }

            if (patchPiece.empty())
            {
                //- The part lies entirely outside the surface. It separates the
                //- two pieces, but there is no surface along it, so it is not an
                //- interface: the plane extended past the patch is not a face.
                droppedArea += polygonArea(part);
                nDropped++;
                nDropped0++;
                continue;
            }

            //- Kept as one polygon where the part lies inside the patch, which is
            //- what makes a run of coplanar facets of a developable surface ONE
            //- interface face rather than one per facet; the pieces otherwise.
            DynamicList<pointField> emit(1);

            if (mag(patchPieceArea - polygonArea(part)) <= 1e-6*polygonArea(part))
            {
                emit.append(part);
            }
            else
            {
                forAll(patchPiece, m)
                {
                    emit.append(patchPiece[m]);
                }
            }

            //- Probe: is this part the same polygon as one already appended for
            //- this cell? Matching is order-free, because the same polygon can
            //- start at a different vertex. The interface faces arrived
            //- duplicated in the mesh, so the enumeration is counted here.
            if (loopDebug)
            forAll(ifacePoly, qi)
            {
                if (ifacePoly[qi].size() != part.size())
                {
                    continue;
                }

                boolList used(ifacePoly[qi].size(), false);
                bool samePart = true;

                forAll(part, m)
                {
                    bool hit = false;

                    forAll(ifacePoly[qi], qm)
                    {
                        if (!used[qm] && mag(ifacePoly[qi][qm] - part[m]) <= tol)
                        {
                            used[qm] = true;
                            hit = true;
                            break;
                        }
                    }

                    if (!hit)
                    {
                        samePart = false;
                        break;
                    }
                }

                if (samePart)
                {
                    static label nDupPart = 0;
                    nDupPart++;

                    if (nDupPart <= 8)
                    {
                        Info<< "DUPPART plane " << pi << " repeats earlier "
                            << "plane " << ifacePlane[qi] << " size "
                            << part.size() << " area " << polygonArea(part)
                            << " pieces " << pair[0] << "/" << pair[1] << endl;
                    }

                    break;
                }
            }

            if (loopDebug)
            {
                static label nAppend = 0;
                nAppend++;

                if (nAppend <= 10)
                {
                    Info<< "APPEND n " << nAppend
                        << " pi " << pi << " size " << part.size()
                        << " area " << polygonArea(part)
                        << " emit " << emit.size()
                        << " pieceArea " << patchPieceArea
                        << " patchArea " << patchArea
                        << " p0 " << part[0] << endl;
                }
            }

            forAll(emit, ei)
            {
                cutCellPieceFace shared;
                shared.pieceA = pair[0];
                shared.pieceB = pair[1];
                shared.tri = planeTri[pi];
                shared.patch = (shared.tri >= 0 && shared.tri < patchOf.size())
                  ? patchOf[shared.tri] : -1;
                shared.faceInA = -1;
                shared.faceInB = -1;

                result.sharedFaces.append(shared);
                ifacePoly.append(emit[ei]);
                ifacePlane.append(pi);
            }
        }
    }

    if (loopDebug)
    {
        //- The interface the cell should carry is the clipped facets
        //- themselves: they are the surface inside the cell. Comparing the two
        //- sums says whether step 6 loses surface, and how much.
        scalar ifaceArea = 0;

        forAll(ifacePoly, ii)
        {
            ifaceArea += polygonArea(ifacePoly[ii]);
        }

        scalar facetArea = 0;

        forAll(facet, fi)
        {
            facetArea += polygonArea(facet[fi]);
        }

        point ccen = Zero;

        forAll(cell.points(), k)
        {
            ccen += cell.points()[k];
        }

        ccen /= scalar(cell.points().size());

        Info<< "PAIR fromPatch " << nPairedFromPatch
            << " fromCell " << nPairedFromCell
            << " pairless " << nPairless << nl;

        Info<< "IFACE cellc " << ccen
            << " parts " << ifacePoly.size()
            << " area " << ifaceArea << " facets " << facet.size()
            << " area " << facetArea
            << " lost " << (facetArea - ifaceArea)
            << " (" << (facetArea > SMALL
                ? 100*(facetArea - ifaceArea)/facetArea : 0) << "%)"
            << " | dropped " << nDropped << " area " << droppedArea
            << " (>2 " << nDropped3 << ", 1 " << nDropped1
            << ", 0 " << nDropped0 << ")" << endl;
    }

    //- 7. Put each interface face into both of its pieces. The pieces are
    //- rebuilt with it, and the polyhedron builder orients every face away from
    //- the point known to be inside the piece, so the two copies of a shared
    //- face come out opposed as they must.
    //-
    //- Each piece gets ALL of its interface faces in ONE rebuild. Rebuilding it
    //- once per shared face shifted the faces already added, and every
    //- faceInA/faceInB was recorded as index 0 - so after the second rebuild
    //- that index held the LAST interface face, and a cell with n shared faces
    //- emitted that one face n times: 2, 4 and 8 copies of one interface in the
    //- mesh, with the other interfaces lost.
    {
        List<List<pointField>> ifaceOf(result.pieces.size());

        forAll(result.sharedFaces, si)
        {
            const cutCellPieceFace& shared = result.sharedFaces[si];

            ifaceOf[shared.pieceA].append(ifacePoly[si]);
            ifaceOf[shared.pieceB].append(ifacePoly[si]);
        }

        forAll(result.pieces, piecei)
        {
            //- The point inside this piece, from the mean of its regions'
            //- centres, which lies on the piece's surface and is enough to
            //- orient every face of it
            point inside = Zero;
            label nR = 0;

            forAll(pieceRegions[piecei], i)
            {
                const pointField& r = allRegions[pieceRegions[piecei][i]];

                point rc = Zero;

                forAll(r, k)
                {
                    rc += r[k];
                }

                inside += rc/scalar(r.size());
                nR++;
            }

            inside /= scalar(max(nR, 1));

            const pointField& pts = result.pieces[piecei].poly.points();
            const faceList& fcs = result.pieces[piecei].poly.faces();

            List<pointField> add(ifaceOf[piecei].size() + fcs.size());

            label n = 0;

            forAll(ifaceOf[piecei], i)
            {
                add[n++] = ifaceOf[piecei][i];
            }

            forAll(fcs, fj)
            {
                pointField fp(fcs[fj].size());

                forAll(fcs[fj], k)
                {
                    fp[k] = pts[fcs[fj][k]];
                }

                add[n++] = fp;
            }

            result.pieces[piecei].poly = polyhedronFrom(add, inside, tol);
        }

        //- Where each interface face landed, matched by its own points: the
        //- builder need not keep the order of the list it was given, so the
        //- position cannot be assumed.
        auto faceOf = [&](const cutCellPolyhedron& poly, const pointField& want)
        -> label
        {
            forAll(poly.faces(), fi)
            {
                const face& f = poly.faces()[fi];

                if (f.size() != want.size())
                {
                    continue;
                }

                boolList used(want.size(), false);
                bool same = true;

                forAll(f, k)
                {
                    bool hit = false;

                    forAll(want, m)
                    {
                        if
                        (
                            !used[m]
                         && mag(poly.points()[f[k]] - want[m])
                            <= samePointTol(tol)
                        )
                        {
                            used[m] = true;
                            hit = true;
                            break;
                        }
                    }

                    if (!hit)
                    {
                        same = false;
                        break;
                    }
                }

                if (same)
                {
                    return fi;
                }
            }

            return -1;
        };

        forAll(result.sharedFaces, si)
        {
            cutCellPieceFace& shared = result.sharedFaces[si];

            shared.faceInA =
                faceOf(result.pieces[shared.pieceA].poly, ifacePoly[si]);
            shared.faceInB =
                faceOf(result.pieces[shared.pieceB].poly, ifacePoly[si]);

            if (shared.faceInA < 0 || shared.faceInB < 0)
            {
                FatalErrorInFunction
                    << "the interface face " << si << " of this cell is not a "
                    << "face of the pieces it separates (A " << shared.faceInA
                    << ", B " << shared.faceInB << ")" << exit(FatalError);
            }
        }
    }

    //- Where a piece fails to close. An edge used once is either a boundary
    //- the piece is missing or a face that does not reach it; the coordinates
    //- say which: on a cell edge, on a line where two facet planes meet, or
    //- neither.
    if (loopDebug)
    {
        forAll(result.pieces, pi)
        {
            const cutCellPolyhedron& p = result.pieces[pi].poly;

            Map<label> count;
            labelList once;

            forAll(p.faces(), fi)
            {
                const face& f = p.faces()[fi];

                forAll(f, i)
                {
                    const label a = f[i];
                    const label b = f.fcIndex(i) < f.size()
                      ? f[f.fcIndex(i)] : f[0];

                    const label key =
                        min(a, b)*max(label(1), p.points().size()) + max(a, b);

                    if (count.found(key))
                    {
                        count[key]++;
                    }
                    else
                    {
                        count.insert(key, 1);
                    }
                }
            }

            Info<< "  piece " << pi << " unclosed edges:";

            label nOnce = 0;

            forAllConstIters(count, iter)
            {
                if (iter() != 2)
                {
                    nOnce++;

                    if (nOnce <= 8)
                    {
                        //- Decode the key back to its two point labels
                        const label lo = iter.key()/max(label(1), p.points().size());
                        const label hi = iter.key()%max(label(1), p.points().size());

                        const point mid(0.5*(p.points()[lo] + p.points()[hi]));

                        Info<< nl << "    " << iter() << "x "
                            << p.points()[lo] << " -> " << p.points()[hi]
                            << "  mid " << mid
                            << "  onCellEdge "
                            << (onCellEdge(cell, mid, 100*tol) ? "yes" : "no");
                    }
                }
            }

            Info<< " (" << nOnce << " of " << count.size() << ")" << endl;
        }
    }

    //- Hand the regions of each face to the caller with their pieces. Each
    //- region is one of the cell's faces cut by the loop's chains, and the
    //- pieces grouped from them are what it was cut into.
    result.faceRegions.setSize(cell.nFaces());
    result.faceRegionPiece.setSize(cell.nFaces());

    forAll(faceRegions, ffi)
    {
        result.faceRegions[ffi] = faceRegions[ffi];
        result.faceRegionPiece[ffi].setSize(faceRegions[ffi].size());

        label n = 0;

        forAll(faceRegions, fj)
        {
            forAll(faceRegions[fj], ri)
            {
                if (fj == ffi)
                {
                    result.faceRegionPiece[ffi][ri] = pieceOf[n];
                }

                n++;
            }
        }
    }

    if (result.sharedFaces.empty())
    {
        result.pieces.clear();
        return cutState::untouched;
    }

    //- A chain that closed on itself and split nothing means a face is crossed
    //- by a closed loop. That face cannot be split: its far side would be the
    //- face minus the loop's interior, which is an annulus, and a polyMesh face
    //- cannot have a hole. Left whole, it is given to one piece and the other
    //- piece is open across it - and an open piece is not a solid. Its volume
    //- comes out larger than the cell it came from and its centre can land in
    //- any region, which is how such a cell came out with every piece in one
    //- region and its surface represented nowhere. Refused, and named.
    //-
    //- Measured on 06_oblique_box: the box's corner sits 0.0038 below one cell's
    //- face, so the loop round it covers 0.3 per cent of that face and the rest
    //- of the face is the other side of the surface. Refusing is honest - the
    //- method cannot represent the cell - where the pieces it was emitting were
    //- not.
    if (closedChainSplitNothing)
    {
        //- The pieces this cell produced are what the refusal is about: the one
        //- that should have taken the rest of the looped face is open across it.
        if (loopDebug)
        {
            forAll(result.pieces, piecei)
            {
                Info<< "   HOLECHAIN piece " << piecei << " openEdges "
                    << nOpenEdges(result.pieces[piecei].poly)
                    << " nfaces " << result.pieces[piecei].poly.faces().size()
                    << " vol " << result.pieces[piecei].poly.volume() << endl;
            }
        }

        result.pieces.clear();
        result.sharedFaces.clear();
        result.refusal = cutRefusal::faceWithHole;
        return cutState::refused;
    }

    return cutState::cut;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
