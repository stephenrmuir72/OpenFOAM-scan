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
    Foam::regionModels::areaSurfaceFilmModels::aggregationScaleModel

Description
    Runtime-selectable physical event-time and coherence-length model for
    film-separation aggregation.

    The model consumes the morphology and physical diameter already selected
    by \c filmSeparation.  It returns per-face physical scales used by both
    separation-aggregation policies.  No dimensional parcel-mass, time, or
    grouping-length threshold is part of this interface.

Usage
    \verbatim
    separationAggregation
    {
        type connectedComponent;

        scaleModel
        {
            type morphologyBased;
            ...
        }
    }
    \endverbatim

    Available scale-model types:
    \verbatim
      morphologyBased | Morphology-attached capillary/breakup scales
    \endverbatim

SourceFiles
    aggregationScaleModel.C
    aggregationScaleModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_aggregationScaleModel_H
#define regionFaModels_aggregationScaleModel_H

#include "autoPtr.H"
#include "boolList.H"
#include "dictionary.H"
#include "labelField.H"
#include "liquidFilmBase.H"
#include "runTimeSelectionTables.H"
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
                    Class aggregationScaleModel Declaration
\*---------------------------------------------------------------------------*/

class aggregationScaleModel
{
protected:

    // Protected Data

        //- Const reference to the liquid film model
        const liquidFilmBase& film_;

        //- Const reference to the finite-area mesh
        const faMesh& mesh_;


public:

    //- Runtime type information
    TypeName("aggregationScaleModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            aggregationScaleModel,
            dictionary,
            (
                const liquidFilmBase& film,
                const dictionary& dict
            ),
            (film, dict)
        );


    // Constructors

        //- Construct from model type, film, and dictionary
        aggregationScaleModel
        (
            const word& modelType,
            const liquidFilmBase& film,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected physical scale model
        static autoPtr<aggregationScaleModel> New
        (
            const liquidFilmBase& film,
            const dictionary& dict
        );


    // Generated Methods

        //- No copy construct
        aggregationScaleModel(const aggregationScaleModel&) = delete;

        //- No copy assignment
        void operator=(const aggregationScaleModel&) = delete;


    //- Destructor
    virtual ~aggregationScaleModel() = default;


    // Member Functions

        //- Calculate the per-face physical event time [s] and coherence
        //  reach [m].  Inactive faces may be left unconfigured by the model.
        virtual void calculate
        (
            const labelField& morphology,
            const scalarField& physicalDiameter,
            const vectorField& Up,
            const boolList& activeFace,
            scalarField& eventTime,
            scalarField& coherenceLength
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
