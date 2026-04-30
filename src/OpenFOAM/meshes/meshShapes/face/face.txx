/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
#include "triFace.H"
#include "quadFace.H"
#include "labelledTri.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<int SizeMin>
Foam::label Foam::face::triangles
(
    const UList<point>& points,
    DynamicList<face, SizeMin>& tris,
    const flipType flipping
) const
{
    label triI = tris.size();

    tris.resize(triI + nTriangles());
    auto total = triangles(points, triI, tris, flipping);
    tris.resize(triI);  // safety
    return total;
}


template<int SizeMin>
Foam::label Foam::face::triangles
(
    const UList<point>& points,
    DynamicList<triFace, SizeMin>& tris,
    const flipType flipping
) const
{
    label triI = tris.size();

    tris.resize(triI + nTriangles());
    auto total = triangles(points, triI, tris, flipping);
    tris.resize(triI);  // safety
    return total;
}


template<int SizeMin>
Foam::label Foam::face::triangles
(
    const UList<point>& points,
    DynamicList<labelledTri, SizeMin>& tris,
    const label region,
    const flipType flipping
) const
{
    label triI = tris.size();

    tris.resize(triI + nTriangles());
    auto total = triangles(points, triI, tris, region, flipping);
    tris.resize(triI);  // safety
    return total;
}


template<int SizeMin>
Foam::label Foam::face::trianglesFan
(
    const UList<point>& points,
    DynamicList<face, SizeMin>& tris,
    const flipType flipping
) const
{
    label triI = tris.size();

    tris.resize(triI + nTriangles());
    auto total = trianglesFan(points, triI, tris, flipping);
    tris.resize(triI);  // safety
    return total;
}


template<int SizeMin>
Foam::label Foam::face::trianglesFan
(
    const UList<point>& points,
    DynamicList<triFace, SizeMin>& tris,
    const flipType flipping
) const
{
    label triI = tris.size();

    tris.resize(triI + nTriangles());
    auto total = trianglesFan(points, triI, tris, flipping);
    tris.resize(triI);  // safety
    return total;
}


template<int SizeMin>
Foam::label Foam::face::trianglesFan
(
    const UList<point>& points,
    DynamicList<labelledTri, SizeMin>& tris,
    const label region,
    const flipType flipping
) const
{
    label triI = tris.size();

    tris.resize(triI + nTriangles());
    auto total = trianglesFan(points, triI, tris, region, flipping);
    tris.resize(triI);  // safety
    return total;
}


template<int SizeMin1, int SizeMin2>
Foam::label Foam::face::trianglesQuads
(
    const UList<point>& points,
    DynamicList<face, SizeMin1>& tris,
    DynamicList<face, SizeMin2>& quads,
    const flipType flipping
) const
{
    label triI = tris.size();
    label quadI = quads.size();

    // Upper bound on number of possible split faces
    tris.resize(triI + nTriangles());
    quads.resize(quadI + nTriangles());

    auto total = trianglesQuads(points, triI, quadI, tris, quads, flipping);

    // The actual size used
    tris.resize(triI);
    quads.resize(quadI);

    return total;
}


template<int SizeMin>
Foam::label Foam::face::trianglesQuads
(
    const UList<point>& points,
    DynamicList<face, SizeMin>& splitFaces,
    const flipType flipping
) const
{
    label faceI = splitFaces.size();

    // Upper bound on number of possible split faces
    splitFaces.resize(faceI + nTriangles());

    auto total = trianglesQuads(points, faceI, splitFaces, flipping);

    // The actual size used
    splitFaces.resize(faceI);

    return total;
}


template<int SizeMin1, int SizeMin2>
Foam::label Foam::face::trianglesQuads
(
    const UList<point>& points,
    DynamicList<triFace, SizeMin1>& tris,
    DynamicList<quadFace, SizeMin2>& quads,
    const flipType flipping
) const
{
    label triI = tris.size();
    label quadI = quads.size();

    // Upper bound on number of possible split faces
    tris.resize(triI + nTriangles());
    quads.resize(quadI + nTriangles());

    auto total = trianglesQuads(points, triI, quadI, tris, quads, flipping);

    // The actual size used
    tris.resize(triI);
    quads.resize(quadI);

    return total;
}


template<class Type>
Type Foam::face::average
(
    const UList<point>& meshPoints,
    const Field<Type>& fld
) const
{
    // Calculate the average by breaking the face into triangles and
    // area-weighted averaging their averages

    const label nPoints = size();

    // If the face is a triangle, do a direct calculation
    if (nPoints == 3)
    {
        return
        (
            (1.0/3.0) *
            (fld[(*this)[0]] + fld[(*this)[1]] + fld[(*this)[2]])
        );
    }

    Type cf = Zero;
    for (label pointi = 0; pointi < nPoints; ++pointi)
    {
        cf += fld[(*this)[pointi]];
    }
    cf /= nPoints;

    const point centrePoint(this->average(meshPoints));

    scalar sumA = 0;
    Type sumAf = Zero;

    for (label pointi = 0; pointi < nPoints; ++pointi)
    {
        const auto curr = thisLabel(pointi);
        const auto next = nextLabel(pointi);

        // Calculate 3*triangle centre field value
        Type ttcf = (fld[curr] + fld[next] + cf);

        // Calculate 2*triangle area
        const scalar ta = Foam::mag
        (
            (meshPoints[curr] - centrePoint)
          ^ (meshPoints[next] - centrePoint)
        );

        sumA += ta;
        sumAf += ta*ttcf;
    }

    if (sumA > VSMALL)
    {
        return sumAf/(3*sumA);
    }
    else
    {
        return cf;
    }
}


// ************************************************************************* //
