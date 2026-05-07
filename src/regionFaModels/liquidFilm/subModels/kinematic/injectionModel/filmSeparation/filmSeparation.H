/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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

    When the optional \c accumulation mode is enabled, the model decouples
    the film stripping rate from the Lagrangian parcel injection frequency
    by buffering separated mass (and momentum) per face. A Lagrangian parcel
    is only injected when the accumulated mass on a given face exceeds the
    \c minParcelMass threshold or when the \c maxAccumulationTime is reached.
    This prevents the "particle explosion" problem where
    \f$N_{\text{parcels}} \propto N_{\text{faces}} \times N_{\text{steps}}\f$.

    The injected parcel diameter can be set independently of the accumulated
    mass via the \c diameterModel entry. The injected velocity is the
    mass-weighted average of the film velocity during accumulation, with an
    optional stochastic jitter for physical spray dispersion.

    Reference:
    \verbatim
        Segment-based Eulerian-Lagrangian transition (tag:BCD):
            Bao, J., Cleary, M. J., & Day, M. S. (2024).
            Segment-based Eulerian-Lagrangian transition method
            for flat nozzle spray atomization simulation.
            Engineering Applications of Computational Fluid Mechanics,
            18(1), 2391448.
            DOI:10.1080/19942060.2024.2391448
    \endverbatim

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

        // Optional entries for mass accumulation
        accumulation        <bool>;
        minParcelMass       <scalar>;
        maxAccumulationTime <scalar>;
        diameterModel       <word>;
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
      accumulation | Enable mass accumulation buffer  | bool | no   | false
      minParcelMass  | Min accumulated mass to inject  | scalar | no | 1e-8
      maxAccumulationTime | Max time before forced injection <!--
                   -->                                | scalar | no | GREAT
      diameterModel | Diameter model for injected parcels <!--
                   -->                                 | word | no | filmThickness
      velocityJitter | Velocity jitter fraction         | scalar | no | 0
      seed         | Random number generator seed      | label | no   | 0
    \endtable

    Options for the \c model entry:
    \verbatim
      OwenRyley    | Model proposed by Owen-Ryley (1985)
      Friedrich    | Model proposed by Friedrich et al. (2008)
    \endverbatim

    Options for the \c diameterModel entry:
    \verbatim
      filmThickness  | Use local film thickness as parcel diameter
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
#include "Enum.H"
#include "areaFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class filmSeparationModel;

namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                        Class filmSeparation Declaration
\*---------------------------------------------------------------------------*/

class filmSeparation
:
    public injectionModel
{
public:

    // Public Enumerations

        //- Options for the diameter model
        enum diameterModelType : char
        {
            FILM_THICKNESS = 0     //!< "Use local film thickness"
        };

        //- Names for diameterModelType
        static const Enum<diameterModelType> diameterModelTypeNames;


private:

    // Private Data

        //- Film-separation model
        autoPtr<filmSeparationModel> filmSeparationModelPtr_;


        // Mass accumulation

            //- Flag to enable mass accumulation buffer
            bool accumulationEnabled_;

            //- Per-face accumulated mass [kg]
            scalarField massAccumulator_;

            //- Per-face accumulated momentum [kg.m/s]
            Field<vector> momentumAccumulator_;

            //- Per-face accumulation time [s]
            scalarField accumulationTime_;


        // Threshold parameters

            //- Minimum accumulated mass to trigger injection [kg]
            scalar minParcelMass_;

            //- Maximum accumulation time before forced injection [s]
            scalar maxAccumulationTime_;


        // Diameter model

            //- Diameter model type for injected parcels
            diameterModelType diameterModel_;


        // Velocity jitter

            //- Velocity jitter fraction (0 = none, 0.1 = 10% of |U|)
            scalar velocityJitter_;

            //- Random number generator
            Random rndGen_;


        // Velocity communication

            //- Registered injection velocity field for Lagrangian reading
            autoPtr<areaVectorField> UfInjectionPtr_;


    // Private Member Functions

        //- Compute the physical diameter for the injected parcel
        scalar computeDiameter(const label faceI) const;


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
