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
    Foam::regionModels::areaSurfaceFilmModels::separationAggregationModel

Description
    Abstract base class for aggregating finite-area film-separation mass
    into physical event-scale Lagrangian parcel injection events.

    The model operates on buffered separated mass and momentum.  Derived
    classes select representative finite-area faces for injection while
    conserving the accumulated mass and momentum and resetting all member
    faces that contributed to emitted parcels.  Release and grouping decisions
    consume the per-face physical event time and coherence reach supplied by
    \c filmSeparation.

Usage
    Minimal example:
    \verbatim
    separationAggregation
    {
        type    <word>;

        scaleModel
        {
            type <word>;
            ...
        }
    }
    \endverbatim

    The \c scaleModel dictionary is mandatory whenever aggregation is enabled.
    The removed dimensional threshold entries \c minParcelMass,
    \c maxAccumulationTime, \c targetParcelMass, and \c aggregationLength are
    rejected explicitly.

    Options for the \c type entry:
    \verbatim
      perFace            | Aggregate each finite-area face independently
      connectedComponent | Aggregate mass and momentum collectively
    \endverbatim

SourceFiles
    separationAggregationModel.C
    separationAggregationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_separationAggregationModel_H
#define regionFaModels_separationAggregationModel_H

#include "liquidFilmBase.H"
#include "scalarField.H"
#include "vectorField.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "labelField.H"
#include "scalarField.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                 Class separationAggregationModel Declaration
\*---------------------------------------------------------------------------*/

class separationAggregationModel
{
protected:

    // Protected Data

        //- Const reference to the liquid film model
        const liquidFilmBase& film_;

        //- Const reference to the finite-area mesh
        const faMesh& mesh_;


public:

    //- Runtime type information
    TypeName("separationAggregationModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            separationAggregationModel,
            dictionary,
            (
                const liquidFilmBase& film,
                const dictionary& dict
            ),
            (film, dict)
        );


    // Constructors

        //- Construct from model type, film, and dictionary
        separationAggregationModel
        (
            const word& modelType,
            const liquidFilmBase& film,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected model
        static autoPtr<separationAggregationModel> New
        (
            const liquidFilmBase& film,
            const dictionary& dict
        );


    // Generated Methods

        //- No copy construct
        separationAggregationModel(const separationAggregationModel&) = delete;

        //- No copy assignment
        void operator=(const separationAggregationModel&) = delete;


    //- Destructor
    virtual ~separationAggregationModel() = default;


    // Member Functions

        //- True when this legacy aggregation model replaces the supplied
        //  physical droplet diameter with an aggregate-level value
        virtual bool overridesDiameter() const noexcept
        {
            return false;
        }

        //- Aggregate buffered separated mass into injection fields using
        //  physical event and coherence scales. All fields are finite-area
        //  face fields. The diameter accumulator contains the mass-weighted
        //  physical droplet diameter selected before computational
        //  parcelisation.
        virtual void aggregate
        (
            scalarField& massAccumulator,
            vectorField& momentumAccumulator,
            scalarField& diameterAccumulator,
            scalarField& accumulationTime,
            const scalarField& eventTime,
            const scalarField& coherenceLength,
            const labelField& morphology,
            scalarField& massToInject,
            scalarField& diameterToInject,
            vectorField& velocityToInject
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
