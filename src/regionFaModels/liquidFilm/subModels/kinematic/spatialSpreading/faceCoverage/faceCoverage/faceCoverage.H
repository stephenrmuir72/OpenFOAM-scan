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
    a source footprint and outputs the raw coverage measure before kernel
    weighting and the coverage centroid used for kernel-distance
    evaluation.

Usage
    Minimal example:
    spatialSpreading
    {
        // Optional entries
        faceCoverage    <word>;
    }
    \endverbatim

    Options for the \c faceCoverage entry:
    \verbatim
      faceCentre        | Use whole face area (default)
      areaIntersection  | Use intersection face area
    \endverbatim

SourceFiles
    faceCoverage.C
    faceCoverageNew.C

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

        //- Model dictionary for run-time input diagnostics
        const dictionary dict_;


public:

    // Public Data Types

        //- Result of evaluating coverage for a candidate face
        struct coverageResult
        {
            //- Whether the candidate face is included in the footprint
            bool included = false;

            //- Raw area coverage measure before kernel weighting
            scalar measure = 0;

            //- Coverage centroid used for kernel-distance evaluation
            point centre = Zero;
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

        //- Whether this coverage measure needs local target face polygons
        //- to be assembled
        virtual bool needsFacePoints() const noexcept { return false; }

        //- Whether the coverage centroid can differ from the target
        //- face centre, requiring a kernel-distance projection
        virtual bool usesProjectedCentroid() const noexcept { return false; }

        //- Evaluate inclusion and the raw coverage measure for a
        //- candidate face
        virtual coverageResult evaluate
        (
            const spatialSpreadingData& data,
            const footprintShape::footprintGeometry& fp,
            const label targetFace,
            const scalar selectedDistance,
            const scalar boundaryRadius
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
