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

Class
    Foam::regionModels::areaSurfaceFilmModels::geodesic

Description
    Geodesic footprint distance.

    Finds candidates with a distributed weighted-shortest-path (Dijkstra)
    traversal of the partitioned finite-area face graph, using face-centre
    edge lengths and processor-boundary frontiers, bounded by the footprint
    search radius. The traversal cannot cross into a different connected
    component and is decomposition independent.

Usage
    Minimal example:
    \verbatim
    spatialSpreading
    {
        // Optional entries
        footprintDistance  geodesic;
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      footprintDistance | Type name: geodesic           | word | yes  | -
    \endtable

SourceFiles
    geodesic.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_geodesic_H
#define regionFaModels_geodesic_H

#include "footprintDistance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                          Class geodesic Declaration
\*---------------------------------------------------------------------------*/

class geodesic
:
    public footprintDistance
{
public:

    //- Runtime type information
    TypeName("geodesic");


    // Constructors

        //- Construct from dictionary
        explicit geodesic(const dictionary& dict);


    //- Destructor
    virtual ~geodesic() = default;


    // Member Functions

        //- Find candidates with a distributed bounded Dijkstra traversal
        virtual void findCandidates
        (
            const spatialSpreadingData& data,
            const footprintShape& shape,
            const footprintShape::footprintGeometry& fp,
            const label sourceProc,
            const label sourceFace,
            const scalar searchRadius,
            DynamicList<label>& candidateFaces,
            DynamicList<scalar>& selectedDistances
        ) const;

        //- Offset the geodesic distance by the tangent-plane centroid
        //- correction relative to the target-centre planar distance
        virtual scalar kernelDistance
        (
            const scalar selectedDistance,
            const scalar planarDistance,
            const scalar centroidPlanarDistance
        ) const;

        virtual bool usesGeodesicSearch() const noexcept { return true; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
