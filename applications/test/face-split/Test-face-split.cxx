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

using namespace Foam;

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
    argList::addBoolOption("tri", "Triangulate");
    argList::addBoolOption("variance", "Edge-length variance (percent)");
    argList::addVerboseOption("Report face information of split");

    argList args(argc, argv);

    const int optVerbose = args.verbose();
    const int optTriangulate = args.found("tri");

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
