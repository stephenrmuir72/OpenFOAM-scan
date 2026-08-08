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
    Foam::regionModels::areaSurfaceFilmModels::temporalSpreading

Description
    Abstract, conservative temporal distribution of integrated
    Lagrangian-to-film transfers.

    The \c temporalSpreading framework distributes an already integrated
    Lagrangian-to-film transfer over a finite physical duration. Instead of
    applying all incoming mass, momentum, and normal impulse immediately, it
    stores each source batch as a pending event and releases the exact fraction
    that belongs to the current time interval.

    The framework has two independent runtime-selected choices:

    - A duration model that determines the physical duration for each active
    source face.
    - \c temporalKernel: A temporal kernel that determines how the transfer is
    released during that duration.

Usage
    Minimal example:
    \verbatim
    temporalSpreading
    {
        // Mandatory entries
        type                 <word>;

        // Conditional mandatory entries

            // when 'type' is 'characteristicTime'
            timeScaleModel   <word>;

        // Optional entries
        temporalKernel       <word>;
        minStepsPerDuration  <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property              | Description                | Type | Reqd | Deflt
      type                  | Duration model             | word | yes  | -
      timeScaleModel        | Characteristic-time scale model   | word | no | -
      temporalKernel        | Temporal kernel model      | word | no   | topHat
      minStepsPerDuration   | Min steps per duration   | scalar | no   | -1
    \endtable

SourceFiles
    temporalSpreading.C
    temporalSpreadingNew.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_temporalSpreading_H
#define regionFaModels_temporalSpreading_H

#include "PtrList.H"
#include "areaFields.H"
#include "autoPtr.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// Forward Declarations
class temporalKernel;

/*---------------------------------------------------------------------------*\
                  Class temporalSpreading Declaration
\*---------------------------------------------------------------------------*/

class temporalSpreading
{
    // Private Data

        //- Const reference to the finite-area mesh
        const faMesh& mesh_;

        //- Model dictionary for input diagnostics
        const dictionary dict_;

        //- Selected temporal-spreading type
        const word modelType_;

        //- Normalized temporal release kernel
        autoPtr<temporalKernel> kernelPtr_;

        //- Minimum number of time steps resolving an event duration
        scalar minStepsPerDuration_;

        //- Restartable, processor-local event metadata
        localIOdictionary state_;

        //- Next unique local event identifier
        label nextEventId_;


    // Private Member Classes

        //- Pending field state for one source batch
        class event
        {
            // Private Data

                //- Unique local event identifier
                label id_;

                //- Physical start time of the event
                scalar startTime_;

                //- Pending mass field
                autoPtr<areaScalarField> massPtr_;

                //- Pending momentum field
                autoPtr<areaVectorField> momentumPtr_;

                //- Pending normal impulse field
                autoPtr<areaScalarField> normalImpulsePtr_;

                //- Pending physical-duration field
                autoPtr<areaScalarField> durationPtr_;

        public:

            // Constructors

                //- Construct a pending event and optionally restore its fields
                event
                (
                    const faMesh& mesh,
                    const label id,
                    const scalar startTime,
                    const bool readExisting
                );


            // Generated Methods

                //- No copy construct
                event(const event&) = delete;

                //- No copy assignment
                void operator=(const event&) = delete;

                //- Default move construct
                event(event&&) = default;

                //- Default move assignment
                event& operator=(event&&) = default;


            // Member Functions

                //- Return the event identifier
                label id() const noexcept { return id_; }

                //- Return the event start time
                scalar startTime() const noexcept { return startTime_; }

                //- Return the pending mass field
                areaScalarField& mass() noexcept { return massPtr_(); }

                //- Return the pending momentum field
                areaVectorField& momentum() noexcept { return momentumPtr_(); }

                //- Return the pending normal impulse field
                areaScalarField& normalImpulse() noexcept
                {
                    return normalImpulsePtr_();
                }

                //- Return the pending physical-duration field
                areaScalarField& duration() noexcept { return durationPtr_(); }
        };

        //- Pending events. Completed entries are cleared in place.
        PtrList<event> events_;


    // Private Member Functions

        //- Return the restart-state IO object
        static IOobject stateIO(const faMesh& mesh);

        //- Return the IO object for a pending event field
        static IOobject eventIO
        (
            const word& fieldName,
            const faMesh& mesh,
            const label id,
            const bool readExisting
        );

        //- Recreate pending fields from restart metadata
        void restoreState();

        //- Synchronize auto-written metadata with pending event fields
        void updateState();

        //- Issue the optional under-resolution warning for a new source batch
        void checkResolution
        (
            const areaScalarField& duration,
            const areaScalarField& mass,
            const areaVectorField& momentum,
            const areaScalarField& normalImpulse
        ) const;

        //- Calculate model-specific duration at source faces
        virtual void calculateDuration
        (
            areaScalarField& duration,
            const areaScalarField& mass,
            const areaVectorField& momentum,
            const areaScalarField& normalImpulse,
            const areaScalarField& diameter,
            const areaScalarField& rho,
            const areaScalarField& sigma
        ) const = 0;


public:

    //- Runtime type information
    TypeName("temporalSpreading");


    // Declare runtime constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        temporalSpreading,
        dictionary,
        (
            const faMesh& mesh,
            const dictionary& dict
        ),
        (mesh, dict)
    );


    // Constructors

        //- Construct from model type, mesh and dictionary
        temporalSpreading
        (
            const word& modelType,
            const faMesh& mesh,
            const dictionary& dict
        );


    // Generated Methods

        //- No copy construct
        temporalSpreading(const temporalSpreading&) = delete;

        //- No copy assignment
        void operator=(const temporalSpreading&) = delete;


    // Selectors

        //- Select a temporal-spreading model
        static autoPtr<temporalSpreading> New
        (
            const faMesh& mesh,
            const dictionary& dict
        );

        //- Select the optional temporalSpreading sub-dictionary
        static autoPtr<temporalSpreading> NewIfPresent
        (
            const faMesh& mesh,
            const dictionary& filmDict
        );


    //- Destructor
    virtual ~temporalSpreading();


    // Member Functions

        //- True when the model needs mass-weighted parcel diameter mapping
        virtual bool needsDiameter() const noexcept { return false; }

        //- Calculate physical source-face durations and check their resolution
        void duration
        (
            areaScalarField& duration,
            const areaScalarField& mass,
            const areaVectorField& momentum,
            const areaScalarField& normalImpulse,
            const areaScalarField& diameter,
            const areaScalarField& rho,
            const areaScalarField& sigma
        ) const;

        //- Store an already spatially distributed transfer batch
        void enqueue
        (
            const areaScalarField& mass,
            const areaVectorField& momentum,
            const areaScalarField& normalImpulse,
            const areaScalarField& duration,
            const scalar startTime
        );

        //- Release exact temporal fractions over [beginTime, endTime]
        void release
        (
            areaScalarField& mass,
            areaVectorField& momentum,
            areaScalarField& normalImpulse,
            const scalar beginTime,
            const scalar endTime
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
