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
    Foam::regionModels::areaSurfaceFilmModels::faceCoverage

Description
    Abstract base class for the runtime-selectable face-coverage submodel
    of the \c spatialSpreading model.

    A face-coverage model decides whether a candidate face is included in
    the source footprint, and outputs the raw coverage measure and centroid
    for the subsequent kernel-distance evaluation.

Usage
    Minimal example:
    spatialSpreading
    {
        // Optional entries
        faceCoverage    <word>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                 | Type | Reqd  | Deflt
      faceCoverage | Type name                   | word | no    | faceCentre
    \endtable

    Options for the \c faceCoverage entry:
    \verbatim
      faceCentre        | Return whole face area (default)
      areaIntersection  | Return partial intersection face area
    \endverbatim

SourceFiles
    faceCoverage.cxx
    faceCoverageNew.cxx

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_faceCoverage_H
#define regionFaModels_faceCoverage_H

#include "footprintShape.H"
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
                       Class faceCoverage Declaration
\*---------------------------------------------------------------------------*/

class faceCoverage
{
protected:

    // Protected Data

        //- Const reference to model dictionary
        const dictionary& dict_;


public:

    // Public Data Types

        //- Coverage data of a candidate face
        struct coverageResult
        {
            //- Whether the candidate face is included in the footprint
            bool included = false;

            //- Raw area coverage measure before kernel weighting
            scalar measure = 0;

            //- Coverage centroid used for kernel-distance evaluation
            point centre = zero{};
        };


    //- Runtime type information
    TypeName("faceCoverage");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            faceCoverage,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Constructors

        //- Construct from model type and dictionary
        faceCoverage(const word& modelType, const dictionary& dict);


    // Generated Methods

        //- No copy construct
        faceCoverage(const faceCoverage&) = delete;

        //- No copy assignment
        void operator=(const faceCoverage&) = delete;


    // Selectors

        //- Select a face-coverage model
        static autoPtr<faceCoverage> New(const dictionary& dict);


    //- Destructor
    virtual ~faceCoverage();


    // Member Functions

        //- Flag to check whether corner points of target faces are needed
        //  e.g. 'areaIntersection' needs face polygons while 'faceCentre'
        //  can work without them
        virtual bool needsFacePoints() const noexcept { return false; }

        //- Flag to check whether the coverage centre is the target face centre
        virtual bool usesProjectedCentroid() const noexcept { return false; }

        //- Calculate face-coverage data for a candidate face
        virtual coverageResult evaluate
        (
            const spatialSpreadingData& data,  // target-face geometry
            const footprintShape::footprintGeometry& fp,
            const label targetFacei,
            const scalar candidateDistance,
            const scalar footprintRadius
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
