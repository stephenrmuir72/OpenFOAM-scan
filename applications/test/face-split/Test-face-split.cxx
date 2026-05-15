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

Application
    Test-face-split

Description
    Test face splitting behaviour

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "labelledTri.H"
#include "edge.H"
#include "faceList.H"
#include "pointList.H"
#include "ListOps.H"

#include "Random.H"
#include "vector.H"
#include "Ostream.H"
#include "fileName.H"
#include "OFstream.H"
#include "foamVtkSurfaceWriter.H"
#include "mathematicalConstants.H"

using namespace Foam;

namespace
{

enum class refSplitMode
{
    TRIANGLES,
    TRIANGLES_QUADS_SINGLE,
    TRIANGLES_QUADS_SEPARATE
};


void calcEdgeVectors
(
    List<vector>& edgeVecs,
    const UList<point>& points,
    const face& f
)
{
    edgeVecs.resize(f.size());

    forAll(edgeVecs, i)
    {
        edgeVecs[i] = Foam::normalised
        (
            points[f.nextLabel(i)]
          - points[f.thisLabel(i)]
        );
    }
}


label mostConcaveAngle
(
    const UList<point>& points,
    const face& f,
    const UList<vector>& edges,
    scalar& maxAngle
)
{
    const vector n(f.areaNormal(points));

    label index = 0;
    maxAngle = -GREAT;

    forAll(edges, i)
    {
        const vector& leftEdge = edges[f.rcIndex(i)];
        const vector& rightEdge = edges[i];

        const scalar edgeCos = (leftEdge & rightEdge);
        const scalar edgeAngle = Foam::acos(clamp(edgeCos, -1, 1));

        scalar angle;

        if (((rightEdge ^ leftEdge) & n) > 0)
        {
            angle = constant::mathematical::pi + edgeAngle;
        }
        else
        {
            angle = constant::mathematical::pi - edgeAngle;
        }

        if (maxAngle < angle)
        {
            maxAngle = angle;
            index = i;
        }
    }

    return index;
}


void referenceSplit
(
    const face& f,
    const UList<point>& points,
    DynamicList<face>& primary,
    DynamicList<face>& quads,
    const refSplitMode mode
)
{
    if (f.size() < 3)
    {
        FatalErrorInFunction
            << "Cannot split a face with fewer than three vertices: " << f
            << abort(FatalError);
    }
    else if (f.size() == 3)
    {
        primary.push_back(f);
    }
    else if (f.size() == 4)
    {
        if (mode == refSplitMode::TRIANGLES)
        {
            List<vector> edges;
            calcEdgeVectors(edges, points, f);

            scalar angle;
            const label idx0 = mostConcaveAngle(points, f, edges, angle);
            const label idx1 = f.fcIndex(idx0);
            const label idx2 = f.fcIndex(idx1);
            const label idx3 = f.fcIndex(idx2);

            primary.push_back(face{f[idx0], f[idx1], f[idx2]});
            primary.push_back(face{f[idx2], f[idx3], f[idx0]});
        }
        else if (mode == refSplitMode::TRIANGLES_QUADS_SINGLE)
        {
            primary.push_back(f);
        }
        else
        {
            quads.push_back(f);
        }
    }
    else
    {
        List<vector> edges;
        calcEdgeVectors(edges, points, f);

        scalar angle;
        const label startIndex = mostConcaveAngle(points, f, edges, angle);
        const scalar bisectAngle = angle/2;
        const vector& rightEdge = edges[startIndex];

        label index = f.fcIndex(f.fcIndex(startIndex));

        label minIndex = index;
        scalar minDiff = constant::mathematical::pi;

        for (label i = 0; i < f.size() - 3; ++i)
        {
            const vector splitEdge = Foam::normalised
            (
                points[f[index]]
              - points[f[startIndex]]
            );

            const scalar splitCos = splitEdge & rightEdge;
            const scalar splitAngle = Foam::acos(clamp(splitCos, -1, 1));
            const scalar angleDiff = Foam::mag(splitAngle - bisectAngle);

            if (angleDiff < minDiff)
            {
                minDiff = angleDiff;
                minIndex = index;
            }

            index = f.fcIndex(index);
        }

        label diff = 0;
        if (minIndex > startIndex)
        {
            diff = minIndex - startIndex;
        }
        else
        {
            diff = minIndex + f.size() - startIndex;
        }

        const label nPoints1 = diff + 1;
        const label nPoints2 = f.size() - diff + 1;

        face face1(nPoints1);
        index = startIndex;
        for (label i = 0; i < nPoints1; ++i)
        {
            face1[i] = f[index];
            index = f.fcIndex(index);
        }

        face face2(nPoints2);
        index = minIndex;
        for (label i = 0; i < nPoints2; ++i)
        {
            face2[i] = f[index];
            index = f.fcIndex(index);
        }

        referenceSplit(face1, points, primary, quads, mode);
        referenceSplit(face2, points, primary, quads, mode);
    }
}


bool sameFace(const face& a, const face& b)
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


void checkFaces
(
    const word& name,
    const word& context,
    const UList<face>& expected,
    const UList<face>& actual
)
{
    if (expected.size() != actual.size())
    {
        FatalErrorInFunction
            << name << " " << context << " count mismatch: expected "
            << expected.size() << " got " << actual.size() << nl
            << "expected: " << expected << nl
            << "actual: " << actual << nl
            << abort(FatalError);
    }

    forAll(expected, i)
    {
        if (!sameFace(expected[i], actual[i]))
        {
            FatalErrorInFunction
                << name << " " << context << " face " << i
                << " mismatch: expected " << expected[i]
                << " got " << actual[i] << nl
                << "expected: " << expected << nl
                << "actual: " << actual << nl
                << abort(FatalError);
        }
    }
}


face flippedFace(const face& f, const face::flipType flipping)
{
    face work(f);

    if (flipping == face::flipType::FLIP)
    {
        work.flip();
    }

    return work;
}


word flipName(const face::flipType flipping)
{
    return
    (
        flipping == face::flipType::FLIP
      ? word("flip")
      : word("none")
    );
}


void checkTriangles
(
    const word& name,
    const UList<point>& points,
    const face& f,
    const face::flipType flipping
)
{
    DynamicList<face> expected;
    DynamicList<face> unusedQuads;
    referenceSplit
    (
        flippedFace(f, flipping),
        points,
        expected,
        unusedQuads,
        refSplitMode::TRIANGLES
    );

    DynamicList<face> actual;
    f.triangles(points, actual, flipping);

    checkFaces(name, word("triangles_") + flipName(flipping), expected, actual);
}


void checkTrianglesQuads
(
    const word& name,
    const UList<point>& points,
    const face& f,
    const face::flipType flipping
)
{
    const face refFace(flippedFace(f, flipping));

    DynamicList<face> expectedSingle;
    DynamicList<face> unusedQuads;
    referenceSplit
    (
        refFace,
        points,
        expectedSingle,
        unusedQuads,
        refSplitMode::TRIANGLES_QUADS_SINGLE
    );

    DynamicList<face> actualSingle;
    f.trianglesQuads(points, actualSingle, flipping);

    checkFaces
    (
        name,
        word("trianglesQuadsSingle_") + flipName(flipping),
        expectedSingle,
        actualSingle
    );

    DynamicList<face> expectedTris;
    DynamicList<face> expectedQuads;
    referenceSplit
    (
        refFace,
        points,
        expectedTris,
        expectedQuads,
        refSplitMode::TRIANGLES_QUADS_SEPARATE
    );

    DynamicList<face> actualTris;
    DynamicList<face> actualQuads;
    f.trianglesQuads(points, actualTris, actualQuads, flipping);

    checkFaces
    (
        name,
        word("trianglesQuadsTris_") + flipName(flipping),
        expectedTris,
        actualTris
    );
    checkFaces
    (
        name,
        word("trianglesQuadsQuads_") + flipName(flipping),
        expectedQuads,
        actualQuads
    );

    label triI = 7;
    label quadI = 11;
    const label total = f.nTrianglesQuads(points, triI, quadI, flipping);

    if
    (
        triI != 7 + expectedTris.size()
     || quadI != 11 + expectedQuads.size()
     || total != expectedTris.size() + expectedQuads.size()
    )
    {
        FatalErrorInFunction
            << name << " nTrianglesQuads_" << flipName(flipping)
            << " count mismatch: triI=" << triI
            << " quadI=" << quadI
            << " total=" << total
            << " expected triI=" << 7 + expectedTris.size()
            << " expected quadI=" << 11 + expectedQuads.size()
            << " expected total="
            << expectedTris.size() + expectedQuads.size()
            << abort(FatalError);
    }
}


void checkCase
(
    const word& name,
    const UList<point>& points,
    const face& f,
    const bool triangulate
)
{
    for
    (
        const face::flipType flipping
      : {face::flipType::NONE, face::flipType::FLIP}
    )
    {
        if (triangulate)
        {
            checkTriangles(name, points, f, flipping);
        }
        else
        {
            checkTrianglesQuads(name, points, f, flipping);
        }
    }
}


void runReferenceChecks(const bool triangulate)
{
    Info<< "Running deterministic face split checks"
        << (triangulate ? " (triangles)" : " (triangles/quads)")
        << nl;

    checkCase
    (
        "triangle",
        pointField
        ({
            {0, 0, 0},
            {1, 0, 0},
            {0, 1, 0}
        }),
        face{0, 1, 2},
        triangulate
    );

    checkCase
    (
        "convexQuad",
        pointField
        ({
            {0, 0, 0},
            {1, 0, 0.05},
            {1, 1, -0.04},
            {0, 1, 0.02}
        }),
        face{0, 1, 2, 3},
        triangulate
    );

    checkCase
    (
        "concaveQuad",
        pointField
        ({
            {0, 0, 0},
            {1, 0, 0},
            {0.25, 0.25, 0},
            {0, 1, 0}
        }),
        face{0, 1, 2, 3},
        triangulate
    );

    checkCase
    (
        "convexPentagon",
        pointField
        ({
            {0, 0, 0},
            {1, 0, 0.03},
            {1.4, 0.8, -0.02},
            {0.45, 1.35, 0.04},
            {-0.35, 0.65, -0.03}
        }),
        face{0, 1, 2, 3, 4},
        triangulate
    );

    checkCase
    (
        "concavePentagon",
        pointField
        ({
            {0, 0, 0},
            {1.1, 0, 0},
            {0.45, 0.35, 0.02},
            {1.0, 1.1, -0.02},
            {-0.1, 0.9, 0}
        }),
        face{0, 1, 2, 3, 4},
        triangulate
    );

    checkCase
    (
        "warpedHexagon",
        pointField
        ({
            {0, 0, 0},
            {1.0, -0.1, 0.08},
            {1.6, 0.55, -0.05},
            {1.0, 1.35, 0.07},
            {0.05, 1.15, -0.04},
            {-0.45, 0.4, 0.05}
        }),
        face{0, 1, 2, 3, 4, 5},
        triangulate
    );

    Info<< "Face split checks passed" << nl << nl;
}

}  // End anonymous namespace


// Create a planar regular polygon with edge length variance
List<point> createPolygon
(
    label nSides,
    scalar baseRadius,
    scalar edgeLengthVariance,
    Random& rnd,
    const point& centre
)
{
    List<point> pts(nSides);

    // Two orthogonal vectors in the z-plane
    const vector e1(1, 0, 0);
    const vector e2(0, 1, 0);

    for (label i = 0; i < nSides; ++i)
    {
        scalar angle = (i * constant::mathematical::twoPi) / nSides;

        // Base position on circle
        vector radial
        (
            baseRadius * Foam::cos(angle) * e1
          + Foam::sin(angle) * e2
        );

        // Edge length variance
        scalar variance =
        (
            (2.0 * rnd.sample01<scalar>() - 1.0)
          * edgeLengthVariance
        );

        pts[i] = centre + (radial * (1.0 + variance));
    }

    return pts;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addBoolOption("check", "Run deterministic regression checks");
    argList::addBoolOption("tri", "Triangulate");
    argList::addBoolOption("variance", "Edge-length variance (percent)");
    argList::addVerboseOption("Report face information of split");

    argList args(argc, argv);

    const int optVerbose = args.verbose();
    const int optTriangulate = args.found("tri");
    const bool optCheck = args.found("check");

    if (optCheck)
    {
        runReferenceChecks(optTriangulate);
        return 0;
    }

    scalar edgeLengthVariance = 0.20;  // 20% by default

    if (args.readIfPresent("variance", edgeLengthVariance))
    {
        edgeLengthVariance *= 0.01;
    }

    Info<< "Creating test faces with edge length variance: "
        << (edgeLengthVariance * 100) << "%" << nl << nl;


    Random rnd(12345);  // Some reproducible seed

    DynamicList<point> dynPoints;
    DynamicList<face> allFaces;

    // Create faces with 5-10 edges
    for (label nEdges = 5; nEdges <= 10; ++nEdges)
    {
        // Create multiple variants of the same polygon
        for (label variant = 0; variant < 2; ++variant)
        {
            // Different base radii and positions for variety
            scalar baseRadius = 1.0 + 0.5 * variant;
            point centre(3.0 * (nEdges - 5), 2.0 * variant, 0);

            List<point> pts = createPolygon
            (
                nEdges,
                baseRadius,
                edgeLengthVariance,
                rnd,
                centre
            );

            // Face connectivity
            label startIdx = dynPoints.size();
            allFaces.emplace_back(identity(pts.size(), startIdx));
            dynPoints.push_back(std::move(pts));
        }
    }

    // Need a plain pointField
    pointField allPoints(std::move(dynPoints));

    // Output results
    Info<< "Generated " << allFaces.size() << " test faces with "
        << allPoints.size() << " total points" << nl << nl;

    {
        vtk::surfaceWriter writer
        (
            allPoints,
            allFaces,
            vtk::formatType::INLINE_ASCII,
            "input-faces",
            false  // serial only
        );

        writer.writeGeometry();

        Info<< "Writing VTK to "
            << args.relativePath(writer.output()) << nl;
    }

    DynamicList<face> splitFaces;
    {
        // Not particular efficient...

        // Worst-case sizing
        {
            label maxFaces = 0;
            for (const auto& f : allFaces)
            {
                maxFaces += f.nTriangles();
            }
            splitFaces.reserve(maxFaces);
        }

        for (const auto& f : allFaces)
        {
            auto first = splitFaces.size();

            if (optTriangulate)
            {
                f.triangles(allPoints, splitFaces);
            }
            else
            {
                f.trianglesQuads(allPoints, splitFaces);
            }
            auto last = splitFaces.size();

            if (optVerbose)
            {
                Info<< f << " :";
                while (first != last)
                {
                    Info<< " " << splitFaces[first];
                    ++first;
                }
                Info<< nl;
            }
        }

        vtk::surfaceWriter writer
        (
            allPoints,
            splitFaces,
            vtk::formatType::INLINE_ASCII,
            "split-faces",
            false  // serial only
        );

        writer.writeGeometry();

        Info<< "Writing VTK to "
            << args.relativePath(writer.output()) << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
