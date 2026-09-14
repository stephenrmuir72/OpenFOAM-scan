/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 OpenFOAM Foundation
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

#include "cutCellMerge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

//- The area vector of a planar polygon, from a triangle fan about its first
//- point. This is the measure the mesh will use for the face.
vector polygonAreaVector(const pointField& poly)
{
    vector sum = Zero;

    if (poly.size() < 3)
    {
        return sum;
    }

    const point& p0 = poly[0];

    for (label k = 1; k < poly.size() - 1; k++)
    {
        sum += (poly[k] - p0)^(poly[k + 1] - p0);
    }

    return 0.5*sum;
}


//- The centre of a polygon, as the mean of its points. Only ever used to
//- orient a face and to test concavity, never for a volume.
point polygonCentre(const pointField& poly)
{
    point c = Zero;

    forAll(poly, i)
    {
        c += poly[i];
    }

    return c/max(scalar(poly.size()), scalar(1));
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

mergeQuality evaluateMerge
(
    const List<pointField>& boundary,
    const point& reference
)
{
    mergeQuality quality{false, 0, GREAT, 0};

    scalar sumPyramid = 0;
    scalar minPyramid = GREAT;

    vectorField normals(boundary.size());
    pointField centres(boundary.size());

    forAll(boundary, fi)
    {
        const pointField& poly = boundary[fi];

        vector n = polygonAreaVector(poly);

        centres[fi] = polygonCentre(poly);

        //- Orient it away from the reference point. The branch did the same:
        //- the point is inside the cell, so the face's outward normal has a
        //- positive component along the direction to its centre.
        if ((n & (centres[fi] - reference)) < 0)
        {
            n = -n;
        }

        normals[fi] = n;

        const scalar pyramid = (n & (centres[fi] - reference))/3.0;

        sumPyramid += pyramid;

        minPyramid = min(minPyramid, pyramid);
    }

    //- A face is concave when its outward normal points at another boundary
    //- face's centre, which for a convex cell cannot happen
    label concave = 0;

    forAll(boundary, fi)
    {
        const vector unitNormal = normals[fi]/max(mag(normals[fi]), VSMALL);

        forAll(boundary, fj)
        {
            if (fi == fj)
            {
                continue;
            }

            const vector dir =
                (centres[fj] - centres[fi])
               /max(mag(centres[fj] - centres[fi]), VSMALL);

            if ((dir & unitNormal) > 1e-6)
            {
                ++concave;
                break;
            }
        }
    }

    quality.volume = sumPyramid;
    quality.minPyramidRatio =
        minPyramid/max(mag(sumPyramid), VSMALL);
    quality.concaveFaces = concave;
    quality.valid = (sumPyramid > VSMALL);

    return quality;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
