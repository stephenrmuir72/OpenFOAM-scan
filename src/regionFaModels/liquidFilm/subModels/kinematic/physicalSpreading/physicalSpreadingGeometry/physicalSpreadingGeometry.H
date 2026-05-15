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
    Foam::regionModels::areaSurfaceFilmModels::physicalSpreadingGeometry

Description
    Shared read-only global finite-area geometry and adjacency context for
    the \c physicalSpreading footprint component models.

    The cache is rebuilt only when the mesh topology/points change or when
    the caller newly requires face polygon points (needed for
    area-intersection coverage), and is otherwise reused across calls.

SourceFiles
    physicalSpreadingGeometry.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_physicalSpreadingGeometry_H
#define regionFaModels_physicalSpreadingGeometry_H

#include "areaFields.H"
#include "autoPtr.H"
#include "faMesh.H"
#include "globalIndex.H"
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
                 Class physicalSpreadingGeometry Declaration
\*---------------------------------------------------------------------------*/

class physicalSpreadingGeometry
{
    // Private Data

        //- Const reference to the finite-area mesh
        const faMesh& mesh_;

        //- Local face count at the last geometry assembly
        mutable label nFaces_;

        //- Local point count at the last geometry assembly
        mutable label nPoints_;

        //- Point geometry instance at the last geometry assembly
        mutable fileName pointsInstance_;

        //- Face topology instance at the last geometry assembly
        mutable fileName facesInstance_;

        //- Whether the last assembly included global face polygons
        mutable bool hasFacePoints_;

        //- Global area-face numbering
        mutable autoPtr<globalIndex> globalFacesPtr_;

        //- Search tree for global face centres
        mutable autoPtr<indexedOctree<treeDataPoint>> centreTreePtr_;

        //- Global face centres
        mutable vectorField centres_;

        //- Global unit face normals
        mutable vectorField normals_;

        //- Global face areas
        mutable scalarField areas_;

        //- Global face polygons (only assembled when required)
        mutable List<pointField> facePoints_;

        //- Global finite-area face adjacency
        mutable labelListList adjacency_;

        //- Connected-component identifier per global face
        mutable labelList components_;

        //- Maximum centre-to-vertex distance per global face
        mutable scalarField faceRadii_;

        //- Maximum centre-to-vertex distance over the global mesh
        mutable scalar maxFaceRadius_;


public:

    // Constructors

        //- Construct from the finite-area mesh
        explicit physicalSpreadingGeometry(const faMesh& mesh);


    // Generated Methods

        //- No copy construct
        physicalSpreadingGeometry(const physicalSpreadingGeometry&) = delete;

        //- No copy assignment
        void operator=(const physicalSpreadingGeometry&) = delete;


    // Member Functions

        //- Assemble or refresh the shared global geometry and adjacency.
        //- Rebuilds when the mesh moves/changes topology or when
        //- face polygons are newly required.
        void update(const bool needFacePoints) const;


    // Access

        //- Return the finite-area mesh
        const faMesh& mesh() const noexcept { return mesh_; }

        //- Return the global area-face numbering
        const globalIndex& globalFaces() const { return *globalFacesPtr_; }

        //- Return the global face centres
        const vectorField& centres() const noexcept { return centres_; }

        //- Return the global unit face normals
        const vectorField& normals() const noexcept { return normals_; }

        //- Return the global face areas
        const scalarField& areas() const noexcept { return areas_; }

        //- Return the global face polygons
        const List<pointField>& facePoints() const noexcept
        {
            return facePoints_;
        }

        //- Return the global finite-area face adjacency
        const labelListList& adjacency() const noexcept { return adjacency_; }

        //- Return the connected-component identifier per global face
        const labelList& components() const noexcept { return components_; }

        //- Return the maximum centre-to-vertex distance per global face
        const scalarField& faceRadii() const noexcept { return faceRadii_; }

        //- Return the maximum centre-to-vertex distance over the mesh
        scalar maxFaceRadius() const noexcept { return maxFaceRadius_; }

        //- Global faces whose centre lies within the given sphere, sorted
        //- by global face index. Falls back to all global faces, in order,
        //- when no search tree has been assembled.
        labelList candidatesInSphere
        (
            const point& centre,
            const scalar radius
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
