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

#include "geodesic.H"
#include "addToRunTimeSelectionTable.H"
#include "physicalSpreadingGeometry.H"

#include <functional>
#include <queue>
#include <utility>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(geodesic, 0);
addToRunTimeSelectionTable(footprintDistance, geodesic, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

geodesic::geodesic(const dictionary& dict)
:
    footprintDistance(typeName, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void geodesic::findCandidates
(
    const physicalSpreadingGeometry& geom,
    const footprintShape& shape,
    const footprintShape::footprintGeometry& fp,
    const label sourceGlobalFace,
    const scalar searchRadius,
    DynamicList<label>& candidateFaces,
    DynamicList<scalar>& selectedDistances
) const
{
    const vectorField& centres = geom.centres();
    const labelListList& adjacency = geom.adjacency();

    scalarField distance(centres.size(), GREAT);
    distance[sourceGlobalFace] = 0;

    using queueEntry = std::pair<scalar, label>;
    std::priority_queue
    <
        queueEntry,
        std::vector<queueEntry>,
        std::greater<queueEntry>
    > queue;
    queue.emplace(0, sourceGlobalFace);

    candidateFaces.clear();

    while (!queue.empty())
    {
        const scalar d = queue.top().first;
        const label facei = queue.top().second;
        queue.pop();

        if (d > distance[facei] || d > searchRadius)
        {
            continue;
        }

        candidateFaces.append(facei);

        for (const label nbr : adjacency[facei])
        {
            const scalar nextDistance =
                d + mag(centres[nbr] - centres[facei]);

            if (nextDistance < distance[nbr] && nextDistance <= searchRadius)
            {
                distance[nbr] = nextDistance;
                queue.emplace(nextDistance, nbr);
            }
        }
    }

    Foam::sort(candidateFaces);

    selectedDistances.clear();
    selectedDistances.reserve(candidateFaces.size());
    for (const label facei : candidateFaces)
    {
        selectedDistances.append(distance[facei]);
    }
}


scalar geodesic::kernelDistance
(
    const scalar selectedDistance,
    const scalar planarDistance,
    const scalar centroidPlanarDistance
) const
{
    return max
    (
        scalar(0),
        selectedDistance + centroidPlanarDistance - planarDistance
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
