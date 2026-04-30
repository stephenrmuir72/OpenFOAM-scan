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

\*---------------------------------------------------------------------------*/

#include "quadFace.H"
#include "face.H"
#include "triFace.H"
#include "labelledTri.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label
Foam::quadFace::mostConcaveAngle(const UList<point>& points) const
{
    const vector n(approxNormal(points));

    // Unit edge vectors
    FixedList<vector, quadFace::size()> edges;
    for (label i = 0; i < quadFace::size(); ++i)
    {
        edges[i] = Foam::normalised
        (
            points[nextLabel(i)]
          - points[thisLabel(i)]
        );
    }

    label index(0);
    scalar maxAngle(-10);

    for (label i = 0; i < quadFace::size(); ++i)
    {
        const vector& leftEdge = edges[rcIndex(i)];
        const vector& rightEdge = edges[i];

        // NOTE: is -ve angle since left edge pointing in other direction
        scalar edgeCos = (leftEdge & rightEdge);
        scalar angle = Foam::acos(clamp(edgeCos, -1, 1));

        vector edgeNormal = (rightEdge ^ leftEdge);
        if ((edgeNormal & n) > 0)
        {
            // Concave angle.
            angle = constant::mathematical::pi + angle;
        }
        else
        {
            // Convex angle. Note '-' to take into account that rightEdge
            // and leftEdge are head-to-tail connected.
            angle = constant::mathematical::pi - angle;
        }

        if (maxAngle < angle)
        {
            maxAngle = angle;
            index = i;
        }
    }

    return index;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::quadFace::centre(const UList<point>& points) const
{
    // Calculate the centre by breaking the face into triangles and
    // area-weighted averaging their centres

    // Note: use double precision to avoid overflows when summing

    // Estimated centre by averaging the face points
    const point centrePoint(this->average(points));

    doubleScalar sumA(0);
    doubleVector sumAc(Zero);

    for (label pointi = 0; pointi < size(); ++pointi)
    {
        const auto& thisPoint = points[thisLabel(pointi)];
        const auto& nextPoint = points[nextLabel(pointi)];

        // Calculate 3*triangle centre
        const vector ttc(thisPoint + nextPoint + centrePoint);

        // Calculate 2*triangle area
        const scalar ta = Foam::mag
        (
            (thisPoint - centrePoint)
          ^ (nextPoint - centrePoint)
        );

        sumA += ta;
        sumAc += ta*ttc;
    }

    if (sumA > VSMALL)
    {
        return sumAc/(3.0*sumA);
    }
    else
    {
        return centrePoint;
    }
}


Foam::vector Foam::quadFace::areaNormal(const UList<point>& points) const
{
    // Decompose quad into four triangles from centre

    // The hub : the estimated centre
    const point hub(this->average(points));

    // The spokes : centre to vertex triangulation vectors
    FixedList<vector, quadFace::size()> spokes;
    for (label i = 0; i < quadFace::size(); ++i)
    {
        spokes[i] = (points[(*this)[i]] - hub);
    }

    return
    (
        0.5 *
        // Sum of (2*triangle) area normals
        (
            (spokes[0] ^ spokes[1])
          + (spokes[1] ^ spokes[2])
          + (spokes[2] ^ spokes[3])
          + (spokes[3] ^ spokes[0])
        )
    );
}


int Foam::quadFace::triangles
(
    const UList<point>& points,
    label& triI,
    UList<Foam::face>& tris,
    const flipType flipping
) const
{
    return split(points, triI, tris, flipping);
}


int Foam::quadFace::triangles
(
    const UList<point>& points,
    label& triI,
    UList<triFace>& tris,
    const flipType flipping
) const
{
    return split(points, triI, tris, flipping);
}


int Foam::quadFace::triangles
(
    const UList<point>& points,
    label& triI,
    UList<labelledTri>& tris,
    const label region,
    const flipType flipping
) const
{
    const label start = triI;

    auto total = split(points, triI, tris, flipping);

    // Set region/index information
    for (label i = start; i < triI; ++i)
    {
        tris[i].setIndex(region);
    }

    return total;
}


// ************************************************************************* //
