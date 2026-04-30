/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
#include "pointHit.H"
#include "triangle.H"
#include "line.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::pointHit Foam::face::ray
(
    const point& p,
    const vector& n,
    const UList<point>& meshPoints,
    const intersection::algorithm alg,
    const intersection::direction dir
) const
{
    // Return potential intersection with face with a ray starting
    // at p, direction n (does not need to be normalized)
    // Does face-center decomposition and returns triangle intersection
    // point closest to p.

    // In case of miss the point is the nearest point intersection of the
    // face plane and the ray and the distance is the distance between the
    // intersection point and the nearest point on the face

    const face& f = *this;

    // If the face is a triangle, do a direct calculation
    if (f.size() == 3)
    {
        return triPointRef
        (
            meshPoints, f[0], f[1], f[2]
        ).ray(p, n, alg, dir);
    }

    const point ctr = Foam::average(points(meshPoints));

    scalar nearestHitDist = GREAT;
    scalar nearestMissDist = GREAT;
    bool eligible = false;

    // Initialize to miss, distance = GREAT
    pointHit nearest(p);

    const label nPoints = f.size();

    for (label pointi = 0; pointi < nPoints; ++pointi)
    {
        // Note: for best accuracy, centre point always comes last
        pointHit curHit = triPointRef
        (
            meshPoints[f.thisLabel(pointi)],
            meshPoints[f.nextLabel(pointi)],
            ctr
        ).ray(p, n, alg, dir);

        if (curHit.hit())
        {
            if (Foam::mag(curHit.distance()) < Foam::mag(nearestHitDist))
            {
                nearestHitDist = curHit.distance();
                nearest.hitPoint(curHit.point());
            }
        }
        else if (!nearest.hit())
        {
            // Miss and no hit yet. Update miss statistics.
            if (curHit.eligibleMiss())
            {
                eligible = true;

                // Miss distance is the distance between the plane intersection
                // point and the nearest point of the triangle
                scalar missDist = curHit.point().dist(p + curHit.distance()*n);

                if (missDist < nearestMissDist)
                {
                    nearestMissDist = missDist;
                    nearest.setDistance(curHit.distance());
                    nearest.setPoint(curHit.point());
                }
            }
        }
    }

    if (nearest.hit())
    {
        nearest.setDistance(nearestHitDist);
    }
    else
    {
        // Haven't hit a single face triangle
        nearest.setMiss(eligible);
    }

    return nearest;
}


Foam::pointHit Foam::face::intersection
(
    const point& p,
    const vector& q,
    const point& ctr,
    const UList<point>& meshPoints,
    const intersection::algorithm alg,
    const scalar tol
) const
{
    const face& f = *this;

    // If the face is a triangle, do a direct calculation
    if (f.size() == 3)
    {
        return triPointRef
        (
            meshPoints, f[0], f[1], f[2]
        ).intersection(p, q, alg, tol);
    }

    scalar nearestHitDist = VGREAT;

    // Initialize to miss, distance = GREAT
    pointHit nearest(p);

    const label nPoints = f.size();

    for (label pointi = 0; pointi < nPoints; ++pointi)
    {
        // Note: for best accuracy, centre point always comes last
        pointHit curHit = triPointRef
        (
            meshPoints[f.thisLabel(pointi)],
            meshPoints[f.nextLabel(pointi)],
            ctr
        ).intersection(p, q, alg, tol);

        if (curHit.hit())
        {
            if (Foam::mag(curHit.distance()) < Foam::mag(nearestHitDist))
            {
                nearestHitDist = curHit.distance();
                nearest.hitPoint(curHit.point());
            }
        }
    }

    if (nearest.hit())
    {
        nearest.setDistance(nearestHitDist);
    }

    return nearest;
}


Foam::pointHit Foam::face::nearestPoint
(
    const point& p,
    const UList<point>& meshPoints
) const
{
    // Dummy labels
    label nearType(-1);
    label nearLabel(-1);

    return nearestPointClassify(p, meshPoints, nearType, nearLabel);
}


Foam::pointHit Foam::face::nearestPointClassify
(
    const point& p,
    const UList<point>& meshPoints,
    label& nearType,
    label& nearLabel
) const
{
    const face& f = *this;

    // If the face is a triangle, do a direct calculation
    if (f.size() == 3)
    {
        return triPointRef
        (
            meshPoints, f[0], f[1], f[2]
        ).nearestPointClassify(p, nearType, nearLabel);
    }

    const point ctr = this->centre(meshPoints);

    // Initialize to miss, distance=GREAT
    pointHit nearest(p);

    nearType = -1;
    nearLabel = -1;

    const label nPoints = f.size();

    for (label pointi = 0; pointi < nPoints; ++pointi)
    {
        // Note: for best accuracy, centre point always comes last
        triPointRef tri
        (
            meshPoints[f.thisLabel(pointi)],
            meshPoints[f.nextLabel(pointi)],
            ctr
        );

        label tmpNearType(-1);
        label tmpNearLabel(-1);

        pointHit curHit = tri.nearestPointClassify
        (
            p,
            tmpNearType,
            tmpNearLabel
        );

        if (Foam::mag(curHit.distance()) < Foam::mag(nearest.distance()))
        {
            nearest.setDistance(curHit.distance());

            // Assume at first that the near type is NONE on the
            // triangle (i.e. on the face of the triangle) then it is
            // therefore also for the face.

            nearType = proxType::NONE;

            if (tmpNearType == triPointRef::EDGE && tmpNearLabel == 0)
            {
                // If the triangle edge label is 0, then this is also
                // an edge of the face, if not, it is on the face

                nearType = proxType::EDGE;
                nearLabel = pointi;
            }
            else if (tmpNearType == triPointRef::POINT && tmpNearLabel < 2)
            {
                // If the triangle point label is 0 or 1, then this is
                // also a point of the face, if not, it is on the face

                nearType = proxType::POINT;
                nearLabel = pointi + tmpNearLabel;
            }

            if (curHit.hit())
            {
                nearest.hitPoint(curHit.point());
            }
            else
            {
                // In nearest point, miss is always eligible
                nearest.setMiss(true);
                nearest.setPoint(curHit.point());
            }
        }
    }

    return nearest;
}


int Foam::face::sign
(
    const point& p,
    const UList<point>& points,
    const scalar tol
) const
{
    // Take three points [0, 1/3, 2/3] from the face
    // - assumes the face is not severely warped

    return triPointRef
    (
        points, (*this)[0], (*this)[size()/3], (*this)[(2*size())/3]
    ).sign(p, tol);
}


// ************************************************************************* //
