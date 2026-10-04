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
    Local finite-area geometry and distributed component data used by the
    spatial-spreading models. Geometry-dependent data is prepared once for a
    fixed mesh. No replicated finite-area face data is stored.

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_spatialSpreadingData_H
#define regionFaModels_spatialSpreadingData_H

#include "faMesh.H"
#include "indexedOctree.H"
#include "treeBoundBox.H"
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

        // Local topology and search data

        //- One group label for each local face
        //  Faces that can reach one another by walking across connected edges
        //  end up with the same group label, including across proc boundaries
        mutable labelList localComponents_;

        //- Search tree for the face centres
        mutable autoPtr<indexedOctree<treeDataPoint>> centreTreePtr_;

        // Axis-aligned bounds of local face centres

        //- Bounds of this proc's face centres (empty on ranks without faces)
        mutable treeBoundBox centreBounds_;

        // Optional local polygon data

        //- List of corner points for each local face
        mutable List<pointField> facePoints_;

        //- The largest distance from face centre to any of that face's corners
        mutable scalarField faceRadii_;

        //- The largest of the face radii
        mutable scalar maxFaceRadius_;

        // Cache state

        //- Flag to mark local search data as built, including empty ranks
        mutable bool centreSearchBuilt_;

        //- Flag to mark distributed component labels as built
        mutable bool componentsValid_;

        //- Flag to mark optional face corners and radii as built
        mutable bool facePointsValid_;


    // Private Member Functions

        //- Build local centre bounds and search tree
        void buildCentreSearch() const;

        //- Build distributed connected-component labels
        void buildComponents() const;

        //- Build optional local face polygons and radii
        void buildFacePoints() const;


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

        //- Prepare fixed-mesh data on all ranks before accessing its services
        void prepare(const bool needFacePoints) const;

        //- Return local face centres
        const vectorField& centres() const;

        //- Return the local face-centre bounds (empty if no faces are present)
        const treeBoundBox& centreBounds() const;

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

        //- Return distributed component labels for local faces
        const labelList& components() const;

        //- Return local faces within a conservative tangent-plane cylinder
        labelList candidatesInTangentPlane
        (
            const point& centre,
            const vector& normal,
            const vector& majorAxis,
            const vector& minorAxis,
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
