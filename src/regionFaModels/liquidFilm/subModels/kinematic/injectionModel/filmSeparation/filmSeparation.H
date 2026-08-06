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
    selected model receives buffered mass, momentum, diameter moments, and
    accumulation times and controls the release of aggregate parcels.

    The optional \c morphology model classifies the physical liquid structure
    before a morphology-specific \c diameterModel calculates the resulting
    droplet diameter.  The classifier and diameter model are independent of
    computational parcelisation; all emitted particles remain ordinary
    spherical Lagrangian parcels.

    The \c perFace separation-aggregation model keeps this buffered transfer
    local to each finite-area face and therefore avoids the processor-wide
    topology communication used by \c connectedComponent.

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
            ...
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
      velocityJitter | Aggregate velocity jitter fraction | scalar | no | 0
      seed         | Random number generator seed      | label | no   | 0
    \endtable

    Options for the \c model entry:
    \verbatim
      OwenRyley    | Model proposed by Owen-Ryley (1985)
      Friedrich    | Model proposed by Friedrich et al. (2008)
    \endverbatim

    The inherited entries are elaborated in:
    - \link injectionModel.H \endlink

SourceFiles
    filmSeparation.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_filmSeparation_H
#define Foam_filmSeparation_H

#include "injectionModel.H"
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

        //- Per-face accumulated mass [kg] when aggregation is enabled
        IOField<scalar> massAccumulator_;

        //- Per-face accumulated momentum [kg.m/s] when aggregation is enabled
        IOField<vector> momentumAccumulator_;

        //- Per-face accumulated mass-diameter moment [kg.m]
        IOField<scalar> diameterAccumulator_;

        //- Per-face accumulation time [s]
        IOField<scalar> accumulationTime_;


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
