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
    Foam::regionModels::areaSurfaceFilmModels::spatialSpreading

Description
    Abstract base class for the optional conservative spatial spreading of
    Lagrangian-to-film mass and normal-impulse sources. Deposited momentum
    is projected onto each target film face; its normal component is a wall
    reaction and is not added to the film momentum or parcel-impact impulse.
    Source descriptors are exchanged only with spatially relevant ranks;
    normalization and geodesic frontier communication are batched.

Usage
    // Minimal example
    \verbatim
    spatialSpreading
    {
        // Mandatory entries
        type                <word>;

        // Optional entries
        footprintShape      <word>;
        footprintDistance   <word>;
        footprintKernel     <word>;
        faceCoverage        <word>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property        | Description            | Type | Reqd | Deflt
      footprintShape  | Footprint-shape model  | word | no   | circular
      footprintDistance | Footprint-distance mdl | word | no | localTangentPlane
      footprintKernel | Footprint-kernel model | word | no | uniform
      faceCoverage    | Face-coverage model    | word | no | faceCentre
    \endtable

    The inherited entries are elaborated in:
      - \link footprintShape.H \endlink
      - \link footprintDistance.H \endlink
      - \link footprintKernel.H \endlink
      - \link faceCoverage.H \endlink

SourceFiles
    spatialSpreading.C
    spatialSpreadingNew.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_spatialSpreading_H
#define regionFaModels_spatialSpreading_H

#include "areaFields.H"
#include "runTimeSelectionTables.H"
#include "spatialSpreadingData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// Forward Declarations
class footprintShape;
class footprintDistance;
class footprintKernel;
class faceCoverage;

/*---------------------------------------------------------------------------*\
                      Class spatialSpreading Declaration
\*---------------------------------------------------------------------------*/

class spatialSpreading
{
public:

    //- Integrated transfer from one impact on a local finite-area face
    struct sourceInput
    {
        label face = -1;
        scalar mass = 0;
        vector momentum = Zero;
        scalar normalImpulse = 0;
        scalar diameter = 0;
        scalar duration = 0;
    };

    //- A locally deposited, independently timed source-to-target transfer
    struct targetContribution
    {
        label sourceProc = -1;
        label sourceId = -1;
        label face = -1;
        scalar mass = 0;
        vector momentum = Zero;
        scalar normalImpulse = 0;
        scalar duration = 0;
    };


protected:

    // Protected Data

        //- Const reference to the finite-area mesh
        const faMesh& mesh_;

        //- Model dictionary for run-time input diagnostics
        const dictionary dict_;

        //- Local finite-area geometry and distributed graph services
        mutable spatialSpreadingData data_;

        //- Selected footprint-shape component model
        autoPtr<footprintShape> shapePtr_;

        //- Selected footprint-distance component model
        autoPtr<footprintDistance> distancePtr_;

        //- Selected distribution-kernel component model
        autoPtr<footprintKernel> kernelPtr_;

        //- Selected face-coverage component model
        autoPtr<faceCoverage> coveragePtr_;

        //- Search diagnostics used by regression tests and profiling
        mutable label sourceExchangePhases_;
        mutable label weightSumExchangePhases_;
        mutable label weightFactorExchangePhases_;
        mutable label sourceDescriptorsSent_;
        mutable label sourceDescriptorsReceived_;
        mutable label maxSourcesRetained_;
        mutable label candidateFacesVisited_;


    // Protected Member Functions

        //- Spread all source channels using one shared set of weights,
        //- delegating shape/distance/coverage/kernel evaluation to the
        //- selected component models
        void spreadSources
        (
            areaScalarField& mass,
            areaVectorField& momentum,
            areaScalarField& normalImpulse,
            const areaScalarField& diameter,
            const scalar beta,
            areaScalarField* temporalWeight = nullptr,
            areaScalarField* temporalMoment = nullptr,
            const UList<sourceInput>* impacts = nullptr,
            const areaScalarField* legacyDuration = nullptr,
            DynamicList<targetContribution>* transfers = nullptr
        ) const;


public:

    //- Runtime type information
    TypeName("spatialSpreading");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            spatialSpreading,
            dictionary,
            (
                const faMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Constructors

        //- Construct from model type, mesh and dictionary
        spatialSpreading
        (
            const word& modelType,
            const faMesh& mesh,
            const dictionary& dict
        );


    // Selectors

        //- Select a physical-spreading model
        static autoPtr<spatialSpreading> New
        (
            const faMesh& mesh,
            const dictionary& dict
        );

        //- Select the optional spatialSpreading sub-dictionary
        static autoPtr<spatialSpreading> NewIfPresent
        (
            const faMesh& mesh,
            const dictionary& filmDict
        );


    // Generated Methods

        //- No copy construct
        spatialSpreading(const spatialSpreading&) = delete;

        //- No copy assignment
        void operator=(const spatialSpreading&) = delete;


    //- Destructor
    virtual ~spatialSpreading();


    // Member Functions

        //- Spread mass, momentum and normal impulse together
        virtual void spread
        (
            areaScalarField& mass,
            areaVectorField& momentum,
            areaScalarField& normalImpulse,
            const areaScalarField& diameter
        ) const = 0;

        //- Spread mapped legacy sources and individual parcel impacts
        virtual void spread
        (
            areaScalarField& mass,
            areaVectorField& momentum,
            areaScalarField& normalImpulse,
            const areaScalarField& diameter,
            const UList<sourceInput>& impacts
        ) const;

        //- Spread each timed source into individually timed target records
        virtual void spread
        (
            areaScalarField& mass,
            areaVectorField& momentum,
            areaScalarField& normalImpulse,
            const areaScalarField& diameter,
            const areaScalarField& legacyDuration,
            const UList<sourceInput>& impacts,
            DynamicList<targetContribution>& transfers
        ) const;

        //- Spread transfer channels and temporal duration moments together
        virtual void spread
        (
            areaScalarField& mass,
            areaVectorField& momentum,
            areaScalarField& normalImpulse,
            areaScalarField& temporalWeight,
            areaScalarField& temporalMoment,
            const areaScalarField& diameter
        ) const;

        //- Flag indicating whether parcel diameter mapping is required
        virtual bool needsDiameter() const noexcept { return false; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
