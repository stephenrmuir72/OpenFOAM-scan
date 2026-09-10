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
#include "pointHit.H"
#include "triangle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::pointHit Foam::quadFace::ray
(
    const point& p,
    const vector& n,
    const UList<point>& meshPoints,
    const intersection::algorithm alg,
    const intersection::direction dir
) const
{
    // Fan decomposition from face centre into 4 triangles.
    // Returns triangle intersection point closest to p.
    // Centre point comes last in each triangle for best accuracy.

    const point ctr = centre(meshPoints);

    scalar nearestHitDist = GREAT;
    scalar nearestMissDist = GREAT;
    bool eligible = false;

    // Initialize to miss
    pointHit nearest(p);

    for (label i = 0; i < quadFace::size(); ++i)
    {
        const pointHit curHit = triPointRef
        (
            meshPoints[thisLabel(i)],
            meshPoints[nextLabel(i)],
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
            if (curHit.eligibleMiss())
            {
                eligible = true;

                // Miss distance: between the plane intersection point
                // and the nearest point on the triangle
                const scalar missDist =
                    curHit.point().dist(p + curHit.distance()*n);

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
        nearest.setMiss(eligible);
    }

    return nearest;
}


Foam::pointHit Foam::quadFace::intersection
(
    const point& p,
    const vector& q,
    const UList<point>& meshPoints,
    const intersection::algorithm alg,
    const scalar tol
) const
{
    return intersection(p, q, centre(meshPoints), meshPoints, alg, tol);
}


Foam::pointHit Foam::quadFace::intersection
(
    const point& p,
    const vector& q,
    const point& ctr,
    const UList<point>& meshPoints,
    const intersection::algorithm alg,
    const scalar tol
) const
{
    // Fan decomposition from face centre into 4 triangles.
    // Centre point comes last in each triangle for best accuracy.

    scalar nearestHitDist = VGREAT;

    // Initialize to miss
    pointHit nearest(p);

    for (label i = 0; i < quadFace::size(); ++i)
    {
        const pointHit curHit = triPointRef
        (
            meshPoints[thisLabel(i)],
            meshPoints[nextLabel(i)],
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


Foam::pointHit Foam::quadFace::nearestPoint
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


Foam::pointHit Foam::quadFace::nearestPointClassify
(
    const point& p,
    const UList<point>& meshPoints,
    label& nearType,
    label& nearLabel
) const
{
    // Decompose from face centre into 4 triangles
    const point ctr(this->centre(meshPoints));

    // Initialize to miss, distance = GREAT
    pointHit nearest(p);

    nearType  = -1;
    nearLabel = -1;

    for (label i = 0; i < quadFace::size(); ++i)
    {
        label tmpNearType  = -1;
        label tmpNearLabel = -1;

        // Centre point comes last for best accuracy
        const triPointRef tri
        (
            meshPoints[thisLabel(i)],
            meshPoints[nextLabel(i)],
            ctr
        );

        const pointHit curHit =
            tri.nearestPointClassify(p, tmpNearType, tmpNearLabel);

        if (Foam::mag(curHit.distance()) < Foam::mag(nearest.distance()))
        {
            nearest.setDistance(curHit.distance());

            // Map triangle proximity back to face proximity.
            // Default: on the face interior (NONE).
            nearType = face::NONE;

            if (tmpNearType == triPointRef::EDGE && tmpNearLabel == 0)
            {
                // Triangle edge 0 is face edge i
                nearType  = face::EDGE;
                nearLabel = i;
            }
            else if (tmpNearType == triPointRef::POINT && tmpNearLabel < 2)
            {
                // Triangle point 0 -> face point i
                // Triangle point 1 -> face point at nextLabel(i)
                nearType  = face::POINT;
                nearLabel = (tmpNearLabel == 0) ? i : (i+1) % quadFace::size();
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


// ************************************************************************* //
