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

Note
    When splitting a flipped quad into triangles, the winding order of each
    triangle is reversed and the triangle output order is also reversed.

    Normal split yields:
      * f1 = [idx0, idx1, idx2]
      * f2 = [idx2, idx3, idx0]

    Splitting a flipped quad yields:
      * f1_flip = [idx0, idx3, idx2]  (is f2 with reversed winding)
      * f2_flip = [idx2, idx1, idx0]  (is f1 with reversed winding)

    Thus swapping the middle vertices (the idx1, idx3 entries)
    will do both:
       * Change winding in each triangle
       * Reverses the output order via cross-diagonal

\*---------------------------------------------------------------------------*/

#include "quadFace.H"
#include "triFace.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class TriFaceT>
int Foam::quadFace::split
(
    const UList<point>& points,
    label& triI,
    UList<TriFaceT>& triFaces,
    const flipType flipping
) const
{
    // Split quad directly into triangles
    // Start at point with largest internal angle.
    const auto idx0 = mostConcaveAngle(points);
    const auto idx1 = fcIndex(idx0);
    const auto idx2 = fcIndex(idx1);
    const auto idx3 = fcIndex(idx2);

    // Triangle 1
    auto& f1 = triFaces[triI++];
    {
        f1.resize(3);  // (no-op for FixedList)
        f1[0] = (*this)[idx0];
        f1[1] = (*this)[idx1];  // <- middle vertex
        f1[2] = (*this)[idx2];
    }

    // Triangle 2
    auto& f2 = triFaces[triI++];
    {
        f2.resize(3);  // (no-op for FixedList)
        f2[0] = (*this)[idx2];
        f2[1] = (*this)[idx3];  // <- middle vertex
        f2[2] = (*this)[idx0];
    }

    if (flipping == flipType::FLIP)
    {
        // Handle flipped quad face - see note above
        std::swap(f1[1], f2[1]);  // <- middle vertices
    }

    return 2;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::quadFace::average
(
    const UList<point>& meshPoints,
    const Field<Type>& fld
) const
{
    Type cf = 0.25*(fld[a()] + fld[b()] + fld[c()] + fld[d()]);

    const point centrePoint(this->average(meshPoints));

    scalar sumA = 0;
    Type sumAf = Zero;

    for (label pointi = 0; pointi < quadFace::size(); ++pointi)
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


template<int SizeMin>
int Foam::quadFace::triangles
(
    const UList<point>& points,
    DynamicList<Foam::face, SizeMin>& tris,
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
int Foam::quadFace::triangles
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
int Foam::quadFace::triangles
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


// ************************************************************************* //
