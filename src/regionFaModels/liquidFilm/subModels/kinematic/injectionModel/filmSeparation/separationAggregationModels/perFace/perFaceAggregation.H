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
    Foam::regionModels::areaSurfaceFilmModels::perFaceAggregation

Description
    Conservative separation aggregation performed independently for each
    finite-area face.

    Unlike connected-component aggregation, this model never communicates
    face topology or buffered values between processors.  A face emits its
    own buffered mass when its physical accumulation age reaches the local
    event time supplied by the selected \c aggregationScaleModel.  The
    emitted diameter and velocity are the mass-weighted diameter moment and
    momentum of that face.

    This policy is face-local and can therefore depend on finite-area
    resolution when one physical structure spans multiple faces.

Usage
    \verbatim
    separationAggregation
    {
        type                   perFace;

        scaleModel
        {
            type morphologyBased;
            ...
        }
    }
    \endverbatim

SourceFiles
    perFaceAggregation.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_perFaceAggregation_H
#define regionFaModels_perFaceAggregation_H

#include "separationAggregationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                     Class perFaceAggregation Declaration
\*---------------------------------------------------------------------------*/

class perFaceAggregation
:
    public separationAggregationModel
{
public:

    //- Runtime type information
    TypeName("perFace");


    // Constructors

        //- Construct from film and dictionary
        perFaceAggregation
        (
            const liquidFilmBase& film,
            const dictionary& dict
        );


    // Generated Methods

        //- No copy construct
        perFaceAggregation(const perFaceAggregation&) = delete;

        //- No copy assignment
        void operator=(const perFaceAggregation&) = delete;


    //- Destructor
    virtual ~perFaceAggregation() = default;


    // Member Functions

        //- Apply per-face aggregation to buffered fields.
        static void aggregateFaces
        (
            scalarField& massAccumulator,
            vectorField& momentumAccumulator,
            scalarField& diameterAccumulator,
            scalarField& accumulationTime,
            const scalarField& eventTime,
            scalarField& massToInject,
            scalarField& diameterToInject,
            vectorField& velocityToInject
        );

        //- Aggregate buffered separated mass into injection fields.
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
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
