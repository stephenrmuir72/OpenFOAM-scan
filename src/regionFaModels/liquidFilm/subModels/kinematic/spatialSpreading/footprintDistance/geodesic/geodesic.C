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
#include "spatialSpreadingData.H"

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
    const spatialSpreadingData& data,
    const footprintShape& shape,
    const footprintShape::footprintGeometry& fp,
    const label sourceProc,
    const label sourceFace,
    const scalar searchRadius,
    DynamicList<label>& candidateFaces,
    DynamicList<scalar>& selectedDistances
) const
{
    data.geodesicCandidates
    (
        sourceProc,
        sourceFace,
        searchRadius,
        candidateFaces,
        selectedDistances
    );
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
