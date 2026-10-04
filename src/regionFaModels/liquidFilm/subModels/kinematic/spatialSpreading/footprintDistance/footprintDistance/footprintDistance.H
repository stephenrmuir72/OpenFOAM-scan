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
    Foam::regionModels::areaSurfaceFilmModels::footprintDistance

Description
    Abstract base class for the runtime-selectable footprint-distance
    submodel of the \c spatialSpreading model.

    A footprint-distance model performs a candidate-face search around a
    source footprint and outputs the selected distance metric value for every
    retained candidate, together with the distance used for subsequent kernel
    evaluation when the coverage measure supplies a projected centroid.

Usage
    Minimal example:
    spatialSpreading
    {
        // Optional entries
        footprintDistance    <word>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property          | Description     | Type | Reqd  | Deflt
      footprintDistance | Type name       | word | no    | localTangentPlane
    \endtable

    Options for the \c footprintDistance entry:
    \verbatim
      localTangentPlane    | Local tangent-plane footprint distance (default)
    \endverbatim

SourceFiles
    footprintDistance.cxx
    footprintDistanceNew.cxx

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_footprintDistance_H
#define regionFaModels_footprintDistance_H

#include "footprintShape.H"
#include "scalarField.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// Forward Declarations
class spatialSpreadingData;

/*---------------------------------------------------------------------------*\
                     Class footprintDistance Declaration
\*---------------------------------------------------------------------------*/

class footprintDistance
{
protected:

    // Protected Data

        //- Const reference to model dictionary
        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("footprintDistance");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            footprintDistance,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Constructors

        //- Construct from model type and dictionary
        footprintDistance(const word& modelType, const dictionary& dict);


    // Generated Methods

        //- No copy construct
        footprintDistance(const footprintDistance&) = delete;

        //- No copy assignment
        void operator=(const footprintDistance&) = delete;


    // Selectors

        //- Select a footprint-distance model
        static autoPtr<footprintDistance> New(const dictionary& dict);


    //- Destructor
    virtual ~footprintDistance();


    // Member Functions

        //- Find local candidate faces around a source footprint and calculate
        //- the selected-distance metric value for each retained candidate
        virtual void findCandidates
        (
            const spatialSpreadingData& data,
            const footprintShape::footprintGeometry& fp,
            const label sourceProc,
            const scalar searchRadius,
            labelList& candidateFaces,
            scalarField& candidateDistances
        ) const = 0;

        //- Kernel-evaluation distance used when the coverage measure
        //- supplies a projected centroid different from the target
        //- face centre, i.e. area-intersection coverage
        virtual scalar kernelDistance
        (
            const scalar candidateDistance,
            const scalar planarDistance,
            const scalar centroidPlanarDistance
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
