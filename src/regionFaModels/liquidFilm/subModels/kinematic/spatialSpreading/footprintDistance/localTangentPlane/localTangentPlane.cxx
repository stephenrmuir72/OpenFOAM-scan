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

#include "localTangentPlane.H"
#include "addToRunTimeSelectionTable.H"
#include "spatialSpreadingData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(localTangentPlane, 0);
addToRunTimeSelectionTable(footprintDistance, localTangentPlane, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

localTangentPlane::localTangentPlane(const dictionary& dict)
:
    footprintDistance(typeName, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void localTangentPlane::findCandidates
(
    const spatialSpreadingData& data,
    const footprintShape::footprintGeometry& fp,
    const label sourceProc,
    const scalar searchRadius,
    labelList& candidateFaces,
    scalarField& candidateDistances
) const
{
    const labelList candidates
    (
        data.candidatesInTangentPlane
        (
            fp.centre,
            fp.normal,
            fp.majorAxis,
            fp.minorAxis,
            searchRadius
        )
    );

    candidateFaces.resize(candidates.size());
    candidateDistances.resize(candidates.size());
    label nSelected = 0;

    // Identify face/distance pairs with radial dist on tangent plane
    for (const label targetFacei : candidates)
    {
        scalar x = 0;
        scalar y = 0;

        footprintShape::planarOffset
        (
            fp,
            data.centres()[targetFacei],
            x,
            y
        );

        const scalar dist = Foam::hypot(x, y);

        if (dist <= searchRadius)
        {
            candidateFaces[nSelected] = targetFacei;
            candidateDistances[nSelected] = dist;
            ++nSelected;
        }
    }

    candidateFaces.resize(nSelected);
    candidateDistances.resize(nSelected);
}


scalar localTangentPlane::kernelDistance
(
    const scalar candidateDistance,
    const scalar planarDistance,
    const scalar centroidPlanarDistance
) const
{
    return centroidPlanarDistance;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
