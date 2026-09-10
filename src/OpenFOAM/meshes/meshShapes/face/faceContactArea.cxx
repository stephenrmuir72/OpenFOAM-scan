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
#include "scalarField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::face::areaInContact
(
    const UList<point>& meshPoints,
    const scalarField& v
) const
{
    // Calculate area in contact given displacement of vertices relative to
    // the face plane. Positive displacement is above the face (no contact);
    // negative is in contact

    // Sample the vertex values
    const List<scalar> vertexValue(v, static_cast<const labelUList&>(*this));

    // Loop through vertexValue.
    // - if all greater than 0 : return 0 (no contact);
    // - if all less than zero : return 1
    // - all zeros is assumed to be in contact.

    bool allPositive = true;
    bool allNegative = true;

    for (auto val : vertexValue)
    {
        if (val > 0)
        {
            allNegative = false;
        }
        else
        {
            allPositive = false;
        }
    }

    if (allPositive)
    {
        return 0;
    }
    if (allNegative)
    {
        return 1;
    }

    // There is a partial contact.
    // Algorithm:
    // Go through all edges. if both vertex values for the edge are
    // positive, discard. If one is positive and one is negative,
    // create a point and start the edge with it. If both are
    // negative, add the edge into the new face.  When finished,
    // calculate area of new face and return relative area (0<x<1)

    // Dimension new point list to max possible size
    const labelUList& faceLabels = *this;

    const label nPoints = size();

    DynamicList<point> newFacePoints(2*nPoints);

    for (label idx0 = 0; idx0 < nPoints; ++idx0)
    {
        label idx1 = fcIndex(idx0);
        const auto& pt0 = meshPoints[faceLabels[idx0]];
        const auto& pt1 = meshPoints[faceLabels[idx1]];

        if (vertexValue[idx0] <= 0)
        {
            // This is a point in contact
            newFacePoints.push_back(pt0);
        }

        if
        (
            (vertexValue[idx0] > 0 && vertexValue[idx1] < 0)
         || (vertexValue[idx0] < 0 && vertexValue[idx1] > 0)
        )
        {
            // Edge intersection. Calculate intersection point and add to list
            point intersection =
            (
                pt0
              + vertexValue[idx0]/(vertexValue[idx1] - vertexValue[idx0])
              * (pt0 - pt1)
            );

            newFacePoints.push_back(intersection);
        }
    }

    // Make the sub-face (points are ordered!)
    face subf(Foam::identity(newFacePoints.size()));

    // Calculate relative area
    return subf.mag(newFacePoints)/(this->mag(meshPoints) + VSMALL);
}


// ************************************************************************* //
