/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    Foam::regionModels::areaSurfaceFilmModels::filmSeparation

Description
    The \c filmSeparation is a collection of curvature thin-film separation
    models designed to predict the onset of film separation and mass separation
    in geometries featuring sharp and/or rounded corners.

    When \c separationAggregation is not specified, separated mass is
    transferred immediately on a per-face basis. When it is specified, the
    selected model receives buffered mass, momentum, diameter moments,
    accumulation ages, morphology, and physical event/coherence scales and
    controls the release of aggregate parcels.

    The optional \c morphology model classifies the physical liquid structure
    before a morphology-specific \c diameterModel calculates the resulting
    droplet diameter.  The classifier and diameter model are independent of
    computational parcelisation; all emitted particles remain ordinary
    spherical Lagrangian parcels.

    The \c perFace separation-aggregation model keeps this buffered transfer
    local to each finite-area face and therefore avoids the processor-wide
    topology communication used by \c connectedComponent.

    Aggregation requires a runtime-selected \c scaleModel dictionary.  The
    scale model uses the existing morphology, physical diameter, film
    properties, and primary-phase velocity to calculate the event time and
    coherence reach.  The former dimensional threshold entries
    \c minParcelMass, \c maxAccumulationTime, \c targetParcelMass, and
    \c aggregationLength are not supported.

Usage
    Minimal example by using boundary-condition files:
    \verbatim
    injectionModels
    {
        // Mandatory entries
        filmSeparation
    }

    filmSeparationCoeffs
    {
        model           <word>;

        // Conditional entries

            // Option-1: when model == OwenRyley

            // Option-2: when model == Friedrich

        // Optional separation aggregation
        separationAggregation
        {
            type            perFace | connectedComponent;

            scaleModel
            {
                type morphologyBased;
                ...
            }
        }

        // Optional physical morphology and diameter routing
        morphology
        {
            type            <word>;
            ...
        }

        diameterModels
        {
            directDroplet { type filmThickness; }
            drip          { type filmThickness; }
            ligament      { type RayleighPlateau; }
            sheet         { type KelvinHelmholtz; }
        }

        // Optional aggregate velocity perturbation
        velocityJitter      <scalar>;
        seed                <label>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      model        | Name of the filmSeparation model   | word | yes  | -
      separationAggregation | Optional mass/momentum aggregation
                  | dict | no | -
      scaleModel   | Physical event/coherence scale model when aggregation
                  | dict | conditional | -
      morphology   | Optional physical morphology classifier | dict | no | -
      diameterModels | Per-morphology physical diameter models | dict | no | -
      velocityJitter | Aggregate velocity jitter fraction | scalar | no | 0
      seed         | Random number generator seed      | label | no | 0
    \endtable

    Options for the \c model entry:
    \verbatim
      OwenRyley    | Model proposed by Owen-Ryley (1985)
      Friedrich    | Model proposed by Friedrich et al. (2008)
    \endverbatim

    Morphology options are \c fixed and \c regimeMap.  The initial physical
    states are \c directDroplet, \c drip, \c ligament, and \c sheet. The
    default mapping is film thickness for direct droplets and drips,
    Rayleigh-Plateau breakup for ligaments, and Kelvin-Helmholtz breakup for
    sheets.      The \c fixedDiameter option of \c connectedComponent aggregation remains
    an aggregate-level diameter choice and cannot be combined with an
    explicit morphology model.

    With \c DebugSwitches { filmSeparation 1; }, \c morphology,
    \c diameterToInject, \c velocityToInject, and \c massToInject are written
    on the finite-area mesh. \c morphologyMassTotal is a restartable
    four-entry field containing cumulative separated mass for
    directDroplet, drip, ligament, and sheet in that order.

    The inherited entries are elaborated in:
    - \link injectionModel.H \endlink

SourceFiles
    filmSeparation.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_filmSeparation_H
#define Foam_filmSeparation_H

#include "injectionModel.H"
#include "PtrList.H"
#include "Random.H"
#include "areaFields.H"
#include "IOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class filmSeparationModel;

namespace regionModels
{
namespace areaSurfaceFilmModels
{

// Forward Declarations
class separationAggregationModel;
class aggregationScaleModel;
class morphologyModel;
class diameterModel;

/*---------------------------------------------------------------------------*\
                        Class filmSeparation Declaration
\*---------------------------------------------------------------------------*/

class filmSeparation
:
    public injectionModel
{
    // Private Data

        //- Film-separation model
        autoPtr<filmSeparationModel> filmSeparationModelPtr_;

        //- Separation-aggregation model
        autoPtr<separationAggregationModel> separationAggregationModelPtr_;

        //- Physical event/coherence scale model
        autoPtr<aggregationScaleModel> aggregationScaleModelPtr_;

        //- Physical liquid-structure classifier
        autoPtr<morphologyModel> morphologyModelPtr_;

        //- Physical droplet-diameter models, indexed by morphology state
        PtrList<diameterModel> diameterModels_;


        //- Per-face accumulated mass [kg] when aggregation is enabled
        IOField<scalar> massAccumulator_;

        //- Per-face accumulated momentum [kg.m/s] when aggregation is enabled
        IOField<vector> momentumAccumulator_;

        //- Per-face accumulated mass-diameter moment [kg.m]
        IOField<scalar> diameterAccumulator_;

        //- Per-face accumulation time [s]
        IOField<scalar> accumulationTime_;

        //- Restartable global cumulative separated mass by morphology [kg]
        IOField<scalar> morphologyMassTotal_;

        //- Locally accumulated separated mass by morphology [kg]
        scalarField morphologySeparatedMass_;

        //- Locally accumulated mass-transfer closure residual [kg]
        scalar massConservationResidual_;

        //- Locally accumulated momentum-transfer closure residual [kg.m/s]
        vector momentumConservationResidual_;


        // Velocity jitter

            //- Velocity jitter fraction (0 = none, 0.1 = 10% of |U|)
            scalar velocityJitter_;

            //- Random number generator
            Random rndGen_;


        // Velocity communication

            //- Registered injection velocity field for Lagrangian reading
            autoPtr<areaVectorField> UfInjectionPtr_;


    // Private Member Functions

        //- Return an isotropic random unit vector (uniform on the unit
        //- sphere) for the optional velocity jitter
        vector randomDirection();


public:

    //- Runtime type information
    TypeName("filmSeparation");


    // Constructors

        //- Construct from base film model and dictionary
        filmSeparation
        (
            liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~filmSeparation();


    // Member Functions

        //- Correct film properties due to the film separation
        virtual void correct
        (
            scalarField& availableMass,
            scalarField& massToInject,
            scalarField& diameterToInject
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
