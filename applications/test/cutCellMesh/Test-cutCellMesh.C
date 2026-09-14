/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    Test-cutCellMesh

Description
    Checks cutCellPolyhedron::split against analytic values.

    Four polyhedra are used:

    - a unit hexahedron, to check the basic split and the analytic area of an
      axis-aligned and of an oblique cross-section;
    - a transition cell, matching the 9-face structure found in the S0 cube mesh:
      the refined-neighbour side is four coplanar quads and the four adjacent
      sides are pentagons carrying the subdivision vertex, to check that a split
      preserves the subdivision and stays closed;
    - degenerate cases, where the plane misses the polyhedron, touches it, or
      passes through a single vertex, all of which must be refused.

    Every successful split must satisfy the volume invariant
    V(negative) + V(positive) = V(original).

\*---------------------------------------------------------------------------*/

#include "cutCellPolyhedron.H"
#include "cutCellCutter.H"
#include "plane.H"
#include "DynamicList.H"
#include "ListOps.H"
#include "IOstreams.H"

#include <cmath>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label nChecks = 0;
label nFail = 0;


void checkTrue(const bool ok, const char* what)
{
    nChecks++;

    Info<< (ok ? "  ok    " : "  FAIL  ") << what << endl;

    if (!ok)
    {
        nFail++;
    }
}


void checkClose(const scalar got, const scalar want, const char* what)
{
    const bool ok = mag(got - want) <= 1e-9*max(1.0, mag(want));

    nChecks++;

    Info<< (ok ? "  ok    " : "  FAIL  ") << what
        << "  got " << got << " want " << want << endl;

    if (!ok)
    {
        nFail++;
    }
}


void checkLabel(const label got, const label want, const char* what)
{
    const bool ok = (got == want);

    nChecks++;

    Info<< (ok ? "  ok    " : "  FAIL  ") << what
        << "  got " << got << " want " << want << endl;

    if (!ok)
    {
        nFail++;
    }
}


point P(const scalar x, const scalar y, const scalar z)
{
    return point(x, y, z);
}


//- Build a polyhedron from polygons given as coordinates. Points are
//- deduplicated, and each face is oriented so that its normal points away from
//- 'inside'.
cutCellPolyhedron makePolyhedron(const List<List<point>>& polys, const point& inside)
{
    DynamicList<point> pts;
    faceList faces(polys.size());

    forAll(polys, qi)
    {
        face f(polys[qi].size());

        forAll(polys[qi], i)
        {
            const point& fp = polys[qi][i];

            label found = -1;

            forAll(pts, pi)
            {
                if (mag(pts[pi] - fp) < 1e-12)
                {
                    found = pi;
                    break;
                }
            }

            if (found < 0)
            {
                found = pts.size();
                pts.append(fp);
            }

            f[i] = found;
        }

        faces[qi] = f;
    }

    pointField pts0(pts.size());

    forAll(pts, i)
    {
        pts0[i] = pts[i];
    }

    forAll(faces, fi)
    {
        face& f = faces[fi];

        vector A = Zero;

        forAll(f, i)
        {
            A += pts0[f[i]] ^ pts0[f[(i + 1) % f.size()]];
        }

        point fc = Zero;

        forAll(f, i)
        {
            fc += pts0[f[i]];
        }
        fc /= scalar(f.size());

        if ((A & (fc - inside)) < 0)
        {
            reverse(f);
        }
    }

    return cutCellPolyhedron(pts0, faces);
}


//- Split and check the volume invariant
bool splitAndCheck
(
    const cutCellPolyhedron& poly,
    const plane& pln,
    cutCellPolyhedron& neg,
    cutCellPolyhedron& pos,
    pointField& cut
)
{
    const scalar vol = poly.volume();

    const bool ok = poly.split(pln, neg, pos, cut);

    if (!ok)
    {
        return false;
    }

    checkClose(neg.volume() + pos.volume(), vol, "volume invariant");

    return true;
}


//- Total area of the faces two pieces of a cut share
scalar sharedArea(const cutCellResult& result)
{
    scalar area = 0;

    forAll(result.sharedFaces, si)
    {
        const cutCellPieceFace& shared = result.sharedFaces[si];

        const cutCellPolyhedron& piece = result.pieces[shared.pieceA].poly;
        const face& f = piece.faces()[shared.faceInA];

        pointField poly(f.size());

        forAll(f, i)
        {
            poly[i] = piece.points()[f[i]];
        }

        area += polygonArea(poly);
    }

    return area;
}


//- Axis-aligned box polyhedron from two opposite corners. Lets the tests below
//- choose their own size and position, which matters for aspect ratio and for
//- the distance from the origin.
cutCellPolyhedron boxPolyhedron(const point& lo, const point& hi)
{
    const scalar x0 = lo.x();
    const scalar y0 = lo.y();
    const scalar z0 = lo.z();
    const scalar x1 = hi.x();
    const scalar y1 = hi.y();
    const scalar z1 = hi.z();

    List<List<point>> faces
    ({
        {P(x0,y0,z0), P(x1,y0,z0), P(x1,y1,z0), P(x0,y1,z0)},
        {P(x0,y0,z1), P(x1,y0,z1), P(x1,y1,z1), P(x0,y1,z1)},
        {P(x0,y0,z0), P(x1,y0,z0), P(x1,y0,z1), P(x0,y0,z1)},
        {P(x0,y1,z0), P(x1,y1,z0), P(x1,y1,z1), P(x0,y1,z1)},
        {P(x0,y0,z0), P(x0,y1,z0), P(x0,y1,z1), P(x0,y0,z1)},
        {P(x1,y0,z0), P(x1,y1,z0), P(x1,y1,z1), P(x1,y0,z1)}
    });

    return makePolyhedron(faces, 0.5*(lo + hi));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main()
{
    cutCellPolyhedron neg, pos;
    pointField cut;


    //  Unit hexahedron
    // ~~~~~~~~~~~~~~~

    Info<< nl << "unit hexahedron" << endl;

    List<List<point>> hexQuads
    ({
        {P(0,0,0), P(1,0,0), P(1,1,0), P(0,1,0)},
        {P(0,0,1), P(1,0,1), P(1,1,1), P(0,1,1)},
        {P(0,0,0), P(1,0,0), P(1,0,1), P(0,0,1)},
        {P(0,1,0), P(1,1,0), P(1,1,1), P(0,1,1)},
        {P(0,0,0), P(0,1,0), P(0,1,1), P(0,0,1)},
        {P(1,0,0), P(1,1,0), P(1,1,1), P(1,0,1)}
    });

    const cutCellPolyhedron hex
    (
        makePolyhedron(hexQuads, P(0.5, 0.5, 0.5))
    );

    checkLabel(hex.nFaces(), 6, "hex face count");
    checkTrue(hex.closed(), "hex closed");
    checkClose(hex.volume(), 1.0, "hex volume");
    checkClose(hex.centre().x(), 0.5, "hex centroid x");


    //  Axis-aligned split at the centre
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Info<< nl << "axis-aligned split" << endl;

    const plane xMid(P(0.5, 0, 0), vector(1, 0, 0));

    checkTrue(splitAndCheck(hex, xMid, neg, pos, cut), "split returns true");

    checkClose(neg.volume(), 0.5, "negative volume");
    checkClose(pos.volume(), 0.5, "positive volume");
    checkTrue(neg.closed(), "negative closed");
    checkTrue(pos.closed(), "positive closed");
    checkLabel(neg.nFaces(), 6, "negative face count");
    checkLabel(pos.nFaces(), 6, "positive face count");
    checkLabel(cut.size(), 4, "cut polygon point count");
    checkClose(polygonArea(cut), 1.0, "cut polygon area");


    //  Oblique split through the centre
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Info<< nl << "oblique split" << endl;

    const plane mid111(P(0.5, 0.5, 0.5), vector(1, 1, 1));

    checkTrue(splitAndCheck(hex, mid111, neg, pos, cut), "split returns true");

    checkClose(neg.volume(), 0.5, "negative volume");
    checkClose(pos.volume(), 0.5, "positive volume");
    checkTrue(neg.closed(), "negative closed");
    checkTrue(pos.closed(), "positive closed");
    checkLabel(cut.size(), 6, "cut polygon point count");

    //- A regular hexagon through the six edge midpoints of the unit hex,
    //- side sqrt(1/2), area (3 sqrt(3) / 2) side^2
    checkClose(polygonArea(cut), 0.75*std::sqrt(3.0), "cut polygon area");


    //  Transition cell: a box whose +x side carries four coplanar faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Info<< nl << "transition cell" << endl;

    //- A transition cell as hexRef8 leaves it. A check on the S0 cube mesh
    //- showed one with 9 faces: the refined-neighbour side is four coplanar
    //- quads, and the four adjacent sides are pentagons, each carrying the
    //- extra vertex where the subdivision meets its edge. The sides are not
    //- split, which is what keeps the mesh free of hanging nodes.
    List<List<point>> transFaces
    ({
        //- z = 0: the refined-neighbour side, four coplanar quads
        {P(0,0,0), P(1,0,0), P(1,1,0), P(0,1,0)},
        {P(1,0,0), P(2,0,0), P(2,1,0), P(1,1,0)},
        {P(0,1,0), P(1,1,0), P(1,2,0), P(0,2,0)},
        {P(1,1,0), P(2,1,0), P(2,2,0), P(1,2,0)},

        //- z = 2: the opposite side, whole
        {P(0,0,2), P(2,0,2), P(2,2,2), P(0,2,2)},

        //- the four adjacent sides, pentagons
        {P(0,0,0), P(0,1,0), P(0,2,0), P(0,2,2), P(0,0,2)},
        {P(2,0,0), P(2,1,0), P(2,2,0), P(2,2,2), P(2,0,2)},
        {P(0,0,0), P(1,0,0), P(2,0,0), P(2,0,2), P(0,0,2)},
        {P(0,2,0), P(1,2,0), P(2,2,0), P(2,2,2), P(0,2,2)}
    });

    const cutCellPolyhedron trans
    (
        makePolyhedron(transFaces, P(1, 1, 1))
    );

    checkLabel(trans.nFaces(), 9, "transition cell face count");
    checkTrue(trans.closed(), "transition cell closed");
    checkClose(trans.volume(), 8.0, "transition cell volume");

    //- Cut across the subdivision
    const plane xOne(P(1, 0, 0), vector(1, 0, 0));

    checkTrue
    (
        splitAndCheck(trans, xOne, neg, pos, cut),
        "cross-subdivision split returns true"
    );

    checkClose(neg.volume(), 4.0, "negative volume");
    checkClose(pos.volume(), 4.0, "positive volume");
    checkTrue(neg.closed(), "negative closed");
    checkTrue(pos.closed(), "positive closed");
    checkLabel(neg.nFaces(), 7, "negative face count");
    checkLabel(pos.nFaces(), 7, "positive face count");

    //- The interface inherits the subdivision: four corners plus the extra
    //- vertex where the refined side meets it
    checkLabel(cut.size(), 5, "cut polygon point count");
    checkClose(polygonArea(cut), 4.0, "cut polygon area");

    //- Cut between the subdivision lines
    const plane yOne(P(0, 1, 0), vector(0, 1, 0));

    checkTrue
    (
        splitAndCheck(trans, yOne, neg, pos, cut),
        "in-plane split returns true"
    );

    checkClose(neg.volume(), 4.0, "negative volume");
    checkClose(pos.volume(), 4.0, "positive volume");
    checkTrue(neg.closed(), "negative closed");
    checkTrue(pos.closed(), "positive closed");
    checkLabel(cut.size(), 5, "cut polygon point count");
    checkClose(polygonArea(cut), 4.0, "cut polygon area");


    //  Degenerate cases, all of which must be refused
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Info<< nl << "degenerate cases" << endl;

    checkTrue
    (
        !hex.split(plane(P(-1,0,0), vector(1,0,0)), neg, pos, cut),
        "plane missing the polyhedron is refused"
    );

    checkTrue
    (
        !hex.split(plane(P(0,0,0), vector(1,0,0)), neg, pos, cut),
        "plane containing a face is refused"
    );

    checkTrue
    (
        !hex.split(plane(P(0,0,0), vector(1,1,1)), neg, pos, cut),
        "plane through a single vertex is refused"
    );


    //  Per-cell cutting
    // ~~~~~~~~~~~~~~~~~

    Info<< nl << "per-cell cutting" << endl;

    //- A triangle in the plane x = 1/2, wide enough to span the whole
    //- cross-section of the unit hex
    const pointField spanTri({P(0.5,-1,-1), P(0.5,3,-1), P(0.5,-1,3)});

    cutCellResult result;

    cutState state = cutCellByTriangles
    (
        hex,
        List<pointField>(1, spanTri),
        1e-9,
        result
    );

    checkTrue(state == cutState::cut, "spanning triangle cuts the hex");
    checkLabel(result.planes.size(), 1, "one plane");
    checkLabel(result.pieces.size(), 2, "two pieces");
    checkClose(result.pieces[0].poly.volume(), 0.5, "minus volume");
    checkClose(result.pieces[1].poly.volume(), 0.5, "plus volume");
    checkTrue(result.pieces[0].poly.closed(), "minus closed");
    checkTrue(result.pieces[1].poly.closed(), "plus closed");
    checkClose(sharedArea(result), 1.0, "interface area");
    checkClose(mag(result.planes[0].normal()), 1.0, "plane normal is a unit vector");

    //- A triangle in the same plane covering only part of the cross-section. The
    //- cut is the same, because the plane is the same: cutting is the arrangement
    //- of the planes, and a piece is not a region, so what the interface is, is
    //- decided by classifying the pieces, not here. A partial facet is no longer
    //- refused, which is what lets a surface edge cross a cell.
    const pointField partTri({P(0.5,0,0), P(0.5,0.5,0), P(0.5,0,0.5)});

    state = cutCellByTriangles
    (
        hex,
        List<pointField>(1, partTri),
        1e-9,
        result
    );

    checkTrue(state == cutState::cut, "partial triangle cuts the plane it lies on");
    checkLabel(result.pieces.size(), 2, "two pieces, whatever the triangle covers");
    checkClose(sharedArea(result), 1.0, "the shared face is the whole cross-section");

    //- Two triangles in perpendicular planes, as where a surface edge crosses a
    //- cell. Their arrangement inside the cell is four pieces.
    List<pointField> twoPlanes(2);
    twoPlanes[0] = spanTri;
    twoPlanes[1] = pointField({P(-1,0.5,-1), P(3,0.5,-1), P(-1,0.5,3)});

    state = cutCellByTriangles(hex, twoPlanes, 1e-9, result);

    checkTrue(state == cutState::cut, "two planes cut the hex");
    checkLabel(result.planes.size(), 2, "two planes");
    checkLabel(result.pieces.size(), 4, "four pieces");
    checkClose(result.pieces[0].poly.volume(), 0.25, "first quarter volume");
    checkClose(result.pieces[3].poly.volume(), 0.25, "last quarter volume");

    scalar quarterTotal = 0;

    forAll(result.pieces, pi)
    {
        quarterTotal += result.pieces[pi].poly.volume();
    }

    checkClose(quarterTotal, 1.0, "quarters sum to the cell volume");
    checkClose(sharedArea(result), 2.0, "both plane cross-sections are shared");

    //- A triangle that misses the cell
    const pointField farTri({P(5,0,0), P(5,1,0), P(5,0,1)});

    state = cutCellByTriangles
    (
        hex,
        List<pointField>(1, farTri),
        1e-9,
        result
    );

    checkTrue(state == cutState::untouched, "distant triangle leaves the cell");
    checkLabel(result.nPolygons, 0, "no polygon for a distant triangle");

    //- The transition cell spans [0,2]^3, so its cross-section is 2 by 2 and the
    //- triangle has to reach y + z = 4 to cover it
    const pointField transTri({P(1,-1,-1), P(1,5,-1), P(1,-1,5)});

    state = cutCellByTriangles
    (
        trans,
        List<pointField>(1, transTri),
        1e-9,
        result
    );

    checkTrue(state == cutState::cut, "spanning triangle cuts the transition cell");
    checkLabel(result.pieces.size(), 2, "two pieces");
    checkClose(result.pieces[0].poly.volume(), 4.0, "minus volume");
    checkClose(result.pieces[1].poly.volume(), 4.0, "plus volume");
    checkTrue(result.pieces[0].poly.closed(), "minus closed");
    checkTrue(result.pieces[1].poly.closed(), "plus closed");
    checkClose(sharedArea(result), 4.0, "interface area");


    //  Thin cells, and cells far from the origin
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Info<< nl << "thin cells and cells far from the origin" << endl;

    //- High aspect ratio: 1 x 1 x 1e-4
    const cutCellPolyhedron thin(boxPolyhedron(P(0, 0, 0), P(1, 1, 1e-4)));

    checkClose(thin.volume(), 1e-4, "thin cell volume");

    checkTrue
    (
        splitAndCheck(thin, plane(P(0.5, 0, 0), vector(1, 0, 0)), neg, pos, cut),
        "thin cell splits"
    );

    checkClose(neg.volume(), 5e-5, "thin minus volume");
    checkClose(pos.volume(), 5e-5, "thin plus volume");
    checkTrue(neg.closed(), "thin minus closed");
    checkTrue(pos.closed(), "thin plus closed");
    checkClose(polygonArea(cut), 1e-4, "thin interface area");

    //- The same proportional shape, placed far from the origin. The tetrahedron
    //- formula must not lose the answer to cancellation between terms of order
    //- the coordinate cubed.
    const point farLo(1000, 1000, 1000);
    const point farHi(farLo + vector(0.125, 0.125, 0.001));
    const scalar farVol = 0.125*0.125*0.001;

    const cutCellPolyhedron farThin(boxPolyhedron(farLo, farHi));

    checkClose(farThin.volume(), farVol, "far thin cell volume");
    checkClose(farThin.centre().z(), 1000.0005, "far thin cell centre z");

    checkTrue
    (
        splitAndCheck
        (
            farThin,
            plane(farLo + vector(0.0625, 0, 0), vector(1, 0, 0)),
            neg,
            pos,
            cut
        ),
        "far thin cell splits"
    );

    checkClose(neg.volume(), 0.5*farVol, "far thin minus volume");
    checkClose(pos.volume(), 0.5*farVol, "far thin plus volume");
    checkTrue(neg.closed(), "far thin minus closed");
    checkTrue(pos.closed(), "far thin plus closed");

    //- A small cell, 1e-3 across, so its cross-section is 1e-6. The coverage
    //- guard has to stay relative at this scale instead of falling back to an
    //- absolute tolerance that would accept almost anything.
    const cutCellPolyhedron small(boxPolyhedron(P(0, 0, 0), P(1e-3, 1e-3, 1e-3)));

    checkClose(small.volume(), 1e-9, "small cell volume");

    //- A triangle in the mid plane covering a quarter of the cross-section
    const pointField smallPart
    ({
        P(5e-4, 0, 0), P(5e-4, 5e-4, 0), P(5e-4, 0, 5e-4)
    });

    cutState smallState = cutCellByTriangles
    (
        small,
        List<pointField>(1, smallPart),
        1e-9,
        result
    );

    checkTrue
    (
        smallState == cutState::cut,
        "a partial triangle on a small cell still cuts the plane"
    );
    checkLabel(result.pieces.size(), 2, "two pieces on a small cell");

    //- The same cell cut properly, as a positive control
    const pointField smallSpan
    ({
        P(5e-4, -1, -1), P(5e-4, 3, -1), P(5e-4, -1, 3)
    });

    smallState = cutCellByTriangles
    (
        small,
        List<pointField>(1, smallSpan),
        1e-9,
        result
    );

    checkTrue(smallState == cutState::cut, "spanning triangle cuts a small cell");
    checkTrue(result.pieces[0].poly.closed(), "small minus closed");
    checkTrue(result.pieces[1].poly.closed(), "small plus closed");
    checkClose(sharedArea(result), 1e-6, "small interface area");


    //  Where the points of a cut came from
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Info<< nl << "source of the cut points" << endl;

    //- The hexahedron split at x = 0.5. Every interface point lies on an edge
    //- of the polyhedron, and the driver needs that edge so that the cells
    //- sharing it can be given one shared point.
    {
        cutCellPolyhedron m, p;
        pointField iface;
        cutCellDetail detail;

        checkTrue
        (
            hex.split(xMid, m, p, iface, detail),
            "split reports detail"
        );

        checkLabel(detail.interfaceSource.size(), iface.size(), "interface source count");
        checkLabel(detail.negativeSource.size(), m.points().size(), "minus source count");
        checkLabel(detail.positiveSource.size(), p.points().size(), "plus source count");

        //- Every interface point of this cut is an edge crossing
        label nEdge = 0;

        forAll(detail.interfaceSource, i)
        {
            const cutPointSource& src = detail.interfaceSource[i];

            if (src.second() >= 0)
            {
                nEdge++;
            }
        }

        checkLabel(nEdge, iface.size(), "every interface point is an edge crossing");

        //- A point reported on an edge lies on the segment between the two
        //- points of that edge, which is what lets two cells agree on it
        bool onEdges = true;

        forAll(detail.negativeSource, i)
        {
            const cutPointSource& src = detail.negativeSource[i];

            if (src.second() < 0)
            {
                onEdges = onEdges
                  && mag(m.points()[i] - hex.points()[src.first()]) < 1e-12;
            }
            else
            {
                const point& a = hex.points()[src.first()];
                const point& b = hex.points()[src.second()];

                const vector ab(b - a);
                const scalar t = ((m.points()[i] - a) & ab)/(ab & ab);

                onEdges = onEdges
                  && t >= -1e-12 && t <= 1 + 1e-12
                  && mag(m.points()[i] - (a + t*ab)) < 1e-12;
            }
        }

        checkTrue(onEdges, "minus points sit on their reported edge or point");

        //- Each interface point is shared with the pieces, and the negative
        //- piece reports the same source for it. That is what lets the driver
        //- give the cells sharing an edge one single point, which is what keeps
        //- the assembled mesh conformal.
        bool sourcesMatch = true;

        forAll(detail.interfaceSource, i)
        {
            const cutPointSource& src = detail.interfaceSource[i];

            label found = -1;

            forAll(m.points(), j)
            {
                if (mag(m.points()[j] - iface[i]) < 1e-12)
                {
                    found = j;
                    break;
                }
            }

            if (found < 0)
            {
                sourcesMatch = false;
            }
            else
            {
                sourcesMatch = sourcesMatch
                  && detail.negativeSource[found].first() == src.first()
                  && detail.negativeSource[found].second() == src.second();
            }
        }

        checkTrue(sourcesMatch, "interface points are shared with the minus piece");
    }


    //- The cutter reports the plane it accepted, and how many triangles really
    //- cross the cell, so the driver can tell a cell nothing crosses from one
    //- that a plane fails to divide
    {
        const pointField spanning
        ({
            P(0.5, -1, -1), P(0.5, 3, -1), P(0.5, -1, 3)
        });

        cutCellResult result;

        const cutState state = cutCellByTriangles
        (
            hex, List<pointField>(1, spanning), 1e-9, result
        );

        checkTrue(state == cutState::cut, "the result reports a cut");
        checkLabel(result.nPolygons, 1, "one clipped polygon");
        checkLabel(result.planes.size(), 1, "one plane");
        checkClose(mag(result.planes[0].normal()), 1.0, "reported plane normal is a unit vector");
        checkClose
        (
            result.planes[0].normal() & (P(0.5, 0, 0) - result.planes[0].refPoint()),
            0.0,
            "reported plane passes through the cut"
        );

        //- A triangle that misses the cell gives no polygon at all
        const pointField missing
        ({
            P(2, 0, 0), P(3, 0, 0), P(3, 1, 0)
        });

        const cutState missState = cutCellByTriangles
        (
            hex, List<pointField>(1, missing), 1e-9, result
        );

        checkTrue(missState == cutState::untouched, "missing triangle is untouched");
        checkLabel(result.nPolygons, 0, "no clipped polygon for a missing triangle");
    }


    //  Arrangement of three planes
    //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Info<< nl << "arrangement of three planes" << endl;

    //- Three planes meeting inside the cell, as at a corner of the cube of
    //- 05_cube_cut. The arrangement of three planes is eight pieces, and the
    //- areas in between them are the three plane cross-sections.
    {
        List<pointField> threePlanes(3);
        threePlanes[0] = pointField({P(0.25,-1,-1), P(0.25,3,-1), P(0.25,-1,3)});
        threePlanes[1] = pointField({P(-1,0.75,-1), P(3,0.75,-1), P(-1,0.75,3)});
        threePlanes[2] = pointField({P(-1,-1,0.5), P(3,-1,0.5), P(-1,3,0.5)});

        cutCellResult result;

        const cutState state = cutCellByTriangles
        (
            hex, threePlanes, 1e-9, result
        );

        checkTrue(state == cutState::cut, "three planes cut the hex");
        checkLabel(result.planes.size(), 3, "three planes");
        checkLabel(result.pieces.size(), 8, "eight pieces");

        scalar total = 0;

        forAll(result.pieces, pi)
        {
            total += result.pieces[pi].poly.volume();
        }

        checkClose(total, 1.0, "the eight pieces sum to the cell volume");

        //- Ordered negative first against the planes, so the first piece is the
        //- one on the negative side of every plane and the last on the positive
        //- side of every one
        checkLabel(result.pieces.first().side[0], -1, "first piece sides negative");
        checkLabel(result.pieces.last().side[0], 1, "last piece sides positive");
        checkLabel(result.pieces.first().side[2], -1, "first piece sides negative in z");
        checkLabel(result.pieces.last().side[2], 1, "last piece sides positive in z");

        //- Every internal face is between two pieces, and the three plane
        //- cross-sections add up to three unit squares
        checkClose(sharedArea(result), 3.0, "internal faces of the arrangement");
    }

    //- A plane coinciding with a face of the cell divides nothing, so the cell
    //- is untouched rather than cut into degenerate pieces
    {
        const pointField onFace({P(0,-1,-1), P(0,3,-1), P(0,-1,3)});

        cutCellResult result;

        const cutState state = cutCellByTriangles
        (
            hex, List<pointField>(1, onFace), 1e-9, result
        );

        checkTrue
        (
            state == cutState::untouched,
            "a plane on a face leaves the cell whole"
        );
        checkLabel(result.nPolygons, 1, "the triangle still crosses the cell");
    }

    //- The same plane shifted off the edge by a small but real distance does
    //- divide the cell properly, into two pieces meeting over a thin rectangle.
    //- It is a sliver, which is S3's business, not a degeneracy, so the
    //- threshold must not catch it.
    {
        const pointField thinCut
        ({
            P(1.5,-1,-0.499999), P(1.5,3,-0.499999), P(0.5,-1,0.500001)
        });

        cutCellResult result;

        const cutState state = cutCellByTriangles
        (
            hex, List<pointField>(1, thinCut), 1e-9, result
        );

        checkTrue
        (
            state == cutState::cut,
            "a plane just off the edge still cuts, giving a sliver"
        );
        checkLabel(result.pieces.size(), 2, "into two pieces");
        checkLabel(result.sharedFaces.size(), 1, "sharing one thin face");
    }


    //  Grouping the pieces of a cell
    //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Info<< nl << "grouping the pieces of a cell" << endl;

    //- One piece, one group
    {
        const labelList zone(1, 0);

        const labelList group = groupPieces(1, List<cutCellPieceFace>(), zone);

        checkLabel(group.size(), 1, "one piece gives one group");
        checkLabel(group[0], 0, "the first piece is group zero");
    }

    //- Two pieces sharing a face, one region each: not merged
    {
        cutCellPieceFace shared;
        shared.pieceA = 0;
        shared.faceInA = 0;
        shared.pieceB = 1;
        shared.faceInB = 0;

        const labelList zone({0, 1});

        const labelList group = groupPieces
        (
            2, List<cutCellPieceFace>(1, shared), zone
        );

        checkLabel(max(group) + 1, 2, "two regions stay two groups");
        checkLabel(group[0], 0, "the inside piece keeps group zero");
        checkLabel(group[1], 1, "the outside piece is group one");
    }

    //- The same pair in one region: merged into one cell
    {
        cutCellPieceFace shared;
        shared.pieceA = 0;
        shared.faceInA = 0;
        shared.pieceB = 1;
        shared.faceInB = 0;

        const labelList zone({0, 0});

        const labelList group = groupPieces
        (
            2, List<cutCellPieceFace>(1, shared), zone
        );

        checkLabel(max(group) + 1, 1, "one region sharing a face is one group");
        checkLabel(group[0], 0, "the first piece is group zero");
        checkLabel(group[1], 0, "the second piece joins it");
    }

    //- Two pieces of one region that share no face stay apart, since a cell has
    //- to be face connected
    {
        const labelList zone({0, 0});

        const labelList group = groupPieces(2, List<cutCellPieceFace>(), zone);

        checkLabel(max(group) + 1, 2, "no shared face means no merge");
    }

    //- A chain of three in one region collapses to one group
    {
        cutCellPieceFace first;
        first.pieceA = 0;
        first.faceInA = 0;
        first.pieceB = 1;
        first.faceInB = 0;

        cutCellPieceFace second;
        second.pieceA = 1;
        second.faceInA = 1;
        second.pieceB = 2;
        second.faceInB = 0;

        List<cutCellPieceFace> faces(2);
        faces[0] = first;
        faces[1] = second;

        const labelList zone({0, 0, 0});

        const labelList group = groupPieces(3, faces, zone);

        checkLabel(max(group) + 1, 1, "a chain of one region is one group");
        checkLabel(group[2], 0, "the end of the chain joins it");
    }

    //- On a real cut: a hex cut by three planes, grouped by the side of the
    //- first plane, which is the surface there
    {
        List<pointField> threePlanes(3);
        threePlanes[0] = pointField({P(0.25,-1,-1), P(0.25,3,-1), P(0.25,-1,3)});
        threePlanes[1] = pointField({P(-1,0.75,-1), P(3,0.75,-1), P(-1,0.75,3)});
        threePlanes[2] = pointField({P(-1,-1,0.5), P(3,-1,0.5), P(-1,3,0.5)});

        cutCellResult result;

        cutCellByTriangles(hex, threePlanes, 1e-9, result);

        labelList zone(result.pieces.size());

        forAll(result.pieces, pi)
        {
            const scalar d = result.planes[0].signedDistance
            (
                result.pieces[pi].poly.centre()
            );

            zone[pi] = (d < 0 ? 0 : 1);
        }

        const labelList group = groupPieces
        (
            result.pieces.size(), result.sharedFaces, zone
        );

        checkLabel(result.pieces.size(), 8, "eight pieces to group");
        checkLabel(max(group) + 1, 2, "one group per side of the first plane");
        checkLabel(group[0], 0, "the first piece is group zero");

        //- A group never mixes the two sides, which is what keeps the interface
        bool sameSide = true;

        forAll(group, gi)
        {
            forAll(group, gj)
            {
                if (group[gi] == group[gj] && zone[gi] != zone[gj])
                {
                    sameSide = false;
                }
            }
        }

        checkTrue(sameSide, "a group never mixes regions");
    }


    Info<< nl << nChecks << " checks, " << nFail << " failed" << nl << endl;

    return nFail ? 1 : 0;
}


// ************************************************************************* //
