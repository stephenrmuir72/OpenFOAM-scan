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

    cutState state = cutCellByTriangles
    (
        hex,
        List<pointField>(1, spanTri),
        1e-9,
        neg,
        pos,
        cut
    );

    checkTrue(state == cutState::cut, "spanning triangle cuts the hex");
    checkClose(neg.volume(), 0.5, "minus volume");
    checkClose(pos.volume(), 0.5, "plus volume");
    checkTrue(neg.closed(), "minus closed");
    checkTrue(pos.closed(), "plus closed");
    checkClose(polygonArea(cut), 1.0, "interface area");

    //- A triangle in the same plane that covers only part of the cross-section:
    //- the surface enters and leaves the cell, so both pieces can lie in the
    //- same region and the cut is refused rather than guessed
    const pointField partTri({P(0.5,0,0), P(0.5,0.5,0), P(0.5,0,0.5)});

    state = cutCellByTriangles
    (
        hex,
        List<pointField>(1, partTri),
        1e-9,
        neg,
        pos,
        cut
    );

    checkTrue(state == cutState::refused, "partial triangle is refused");

    //- Two triangles in different planes, as where a surface edge crosses the
    //- cell: refused in this slice, to be handled by the arrangement
    List<pointField> twoPlanes(2);
    twoPlanes[0] = spanTri;
    twoPlanes[1] = pointField({P(-1,0.5,-1), P(3,0.5,-1), P(-1,0.5,3)});

    state = cutCellByTriangles(hex, twoPlanes, 1e-9, neg, pos, cut);

    checkTrue
    (
        state == cutState::refused,
        "two planes crossing the cell are refused"
    );

    //- A triangle that misses the cell
    const pointField farTri({P(5,0,0), P(5,1,0), P(5,0,1)});

    state = cutCellByTriangles
    (
        hex,
        List<pointField>(1, farTri),
        1e-9,
        neg,
        pos,
        cut
    );

    checkTrue(state == cutState::untouched, "distant triangle leaves the cell");

    //- The transition cell spans [0,2]^3, so its cross-section is 2 by 2 and the
    //- triangle has to reach y + z = 4 to cover it
    const pointField transTri({P(1,-1,-1), P(1,5,-1), P(1,-1,5)});

    state = cutCellByTriangles
    (
        trans,
        List<pointField>(1, transTri),
        1e-9,
        neg,
        pos,
        cut
    );

    checkTrue(state == cutState::cut, "spanning triangle cuts the transition cell");
    checkClose(neg.volume(), 4.0, "minus volume");
    checkClose(pos.volume(), 4.0, "plus volume");
    checkTrue(neg.closed(), "minus closed");
    checkTrue(pos.closed(), "plus closed");
    checkClose(polygonArea(cut), 4.0, "interface area");


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
        neg,
        pos,
        cut
    );

    checkTrue
    (
        smallState == cutState::refused,
        "partial triangle on a small cell is refused"
    );

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
        neg,
        pos,
        cut
    );

    checkTrue(smallState == cutState::cut, "spanning triangle cuts a small cell");
    checkTrue(neg.closed(), "small minus closed");
    checkTrue(pos.closed(), "small plus closed");
    checkClose(polygonArea(cut), 1e-6, "small interface area");


    Info<< nl << nChecks << " checks, " << nFail << " failed" << nl << endl;

    return nFail ? 1 : 0;
}


// ************************************************************************* //
