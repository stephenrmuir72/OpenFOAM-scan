/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O  peration     |
    \\  /    A  nd           |
     \\/     M anipulation   |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Keysight Technologies
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::regionModels::areaSurfaceFilmModels::spatialSpreadingData

Description
    Local finite-area geometry and distributed graph services used by the
    spatial-spreading models. Local centre bounds are cached with geometry;
    geodesic frontiers are communicated in batches. No replicated finite-area
    face data is stored.

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_spatialSpreadingData_H
#define regionFaModels_spatialSpreadingData_H

#include "faMesh.H"
#include "indexedOctree.H"
#include "treeDataPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class spatialSpreadingData Declaration
\*---------------------------------------------------------------------------*/

class spatialSpreadingData
{
    // Private Data

        //- Reference to the finite-area mesh
        const faMesh& mesh_;

        //- Local topology and search data
        mutable labelListList adjacency_;
        mutable labelListList boundaryEdges_;
        mutable bool boundaryEdgesValid_;
        mutable labelList localComponents_;
        mutable autoPtr<indexedOctree<treeDataPoint>> centreTreePtr_;

        //- Axis-aligned bounds of local face centres
        mutable point centreBoundsMin_;
        mutable point centreBoundsMax_;
        mutable bool centreBoundsValid_;

        //- Optional local polygon data
        mutable List<pointField> facePoints_;
        mutable scalarField faceRadii_;
        mutable scalar maxFaceRadius_;

        //- Cache state
        mutable bool topologyValid_;
        mutable bool componentsValid_;
        mutable bool facePointsValid_;
        mutable label nFaces_;
        mutable label nPoints_;
        //- Last mesh point event used to build the geometry caches
        mutable label pointEventNo_;
        mutable fileName pointsInstance_;
        mutable fileName facesInstance_;

        //- Search diagnostics
        mutable label tangentBoundCornerChecks_;
        mutable label centreBoundsBuilds_;
        mutable label centreBoundsFaceVisits_;
        mutable label geodesicExchangePhases_;


    // Private Member Functions

        //- Rebuild local topology and invalidate dependent data
        void updateTopology() const;

        //- Build local processor-boundary edge addressing
        void updateBoundaryEdges() const;

        //- Rebuild the local centre bounds
        void updateCentreBounds() const;

        //- Invalidate geometry-dependent search data after point motion
        void updateGeometry() const;

        //- Rebuild distributed connected-component labels
        void updateComponents() const;

        //- Rebuild optional local face polygons and radii
        void updateFacePoints() const;

        //- Return whether the local topology cache is current
        bool topologyCurrent() const;

        //- Return whether the local geometry cache is current
        bool current() const;

        //- Return the point-field event number, including same-time motion
        label currentPointEventNo() const;


public:

    // Constructors

        //- Construct from the finite-area mesh
        explicit spatialSpreadingData(const faMesh& mesh);


    // Generated Methods

        //- No copy construct
        spatialSpreadingData(const spatialSpreadingData&) = delete;

        //- No copy assignment
        void operator=(const spatialSpreadingData&) = delete;


    // Member Functions

        //- Refresh local geometry and optional topology services
        void update(const bool needFacePoints = false) const;

        //- Return finite-area mesh
        const faMesh& mesh() const noexcept { return mesh_; }

        //- Return local face centres
        const vectorField& centres() const;

        //- Return the local face-centre bounds, if any faces are present
        bool centreBounds(point& min, point& max) const;

        //- Return local unit face normals
        const vectorField& normals() const;

        //- Return local face areas
        const scalarField& areas() const;

        //- Return local face polygons
        const List<pointField>& facePoints() const;

        //- Return local face radii
        const scalarField& faceRadii() const;

        //- Return maximum local face radius
        scalar maxFaceRadius() const;

        //- Return local face adjacency
        const labelListList& adjacency() const;

        //- Return distributed component labels for local faces
        const labelList& components() const;

        //- Return local faces whose centres lie within a sphere
        labelList candidatesInSphere
        (
            const point& centre,
            const scalar radius
        ) const;

        //- Number of projected bounding-box corners tested by tangent-plane
        //- searches
        label tangentBoundCornerChecks() const noexcept
        {
            return tangentBoundCornerChecks_;
        }

        //- Number of geometry-bound rebuilds
        label centreBoundsBuilds() const noexcept
        {
            return centreBoundsBuilds_;
        }

        //- Number of local face centres visited while rebuilding bounds
        label centreBoundsFaceVisits() const noexcept
        {
            return centreBoundsFaceVisits_;
        }

        //- Number of batched geodesic communication rounds
        label geodesicExchangePhases() const noexcept
        {
            return geodesicExchangePhases_;
        }

        //- Return local faces within a conservative tangent-plane cylinder
        labelList candidatesInTangentPlane
        (
            const point& centre,
            const vector& normal,
            const vector& majorAxis,
            const vector& minorAxis,
            const scalar radius
        ) const;

        //- Find local candidates using a distributed geodesic search.
        void geodesicCandidates
        (
            const label sourceProc,
            const label sourceFace,
            const scalar searchRadius,
            DynamicList<label>& candidateFaces,
            DynamicList<scalar>& selectedDistances
        ) const;

        //- Find candidates for multiple geodesic sources with one frontier
        //- exchange per propagation round
        void geodesicCandidatesBatch
        (
            const labelUList& sourceProcs,
            const labelUList& sourceFaces,
            const scalarField& searchRadii,
            List<labelList>& candidateFaces,
            List<scalarField>& selectedDistances
        ) const;

        //- Distinct sources may use the same geometric seed face
        void geodesicCandidatesBatch
        (
            const labelUList& sourceProcs,
            const labelUList& sourceIds,
            const labelUList& sourceFaces,
            const scalarField& searchRadii,
            List<labelList>& candidateFaces,
            List<scalarField>& selectedDistances
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
