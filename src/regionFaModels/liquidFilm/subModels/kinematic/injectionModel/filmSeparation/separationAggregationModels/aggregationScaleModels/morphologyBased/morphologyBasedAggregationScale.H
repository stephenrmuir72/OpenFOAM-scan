/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
      \\/     M anipulation  |
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
    Foam::regionModels::areaSurfaceFilmModels::
    morphologyBasedAggregationScale

Description
    Morphology-attached physical event and coherence scales for
    film-separation aggregation.

    The event time is an inertio-capillary time for direct droplets and drips,
    a Rayleigh-Plateau ligament-breakup time for ligaments, and an advective
    Kelvin-Helmholtz time for sheets.  The coherence reach is obtained from
    the selected physical diameter and the distance transported during the
    same event time.

Usage
    \verbatim
    scaleModel
    {
        type morphologyBased;

        // Dimensionless, literature-backed scale coefficients
        Ccapillary  1;
        Cligament  2.915;
        Csheet     1;
        Ccoherence 1;
    }
    \endverbatim

    \c Cligament defaults to the inverse of the fastest-growing inviscid
    Rayleigh-Plateau growth coefficient (approximately 0.343).  The
    Ohnesorge correction uses the same \f$(1+3Oh)^{1/6}\f$ form as the
    existing Rayleigh-Plateau diameter model.

SourceFiles
    morphologyBasedAggregationScale.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_morphologyBasedAggregationScale_H
#define regionFaModels_morphologyBasedAggregationScale_H

#include "aggregationScaleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
             Class morphologyBasedAggregationScale Declaration
\*---------------------------------------------------------------------------*/

class morphologyBasedAggregationScale
:
    public aggregationScaleModel
{
    // Private Data

        //- Inertio-capillary event-time coefficient
        scalar Ccapillary_;

        //- Rayleigh-Plateau ligament event-time coefficient
        scalar Cligament_;

        //- Kelvin-Helmholtz sheet event-time coefficient
        scalar Csheet_;

        //- Coherence-reach coefficient
        scalar Ccoherence_;


public:

    //- Runtime type information
    TypeName("morphologyBased");


    // Static Member Functions

        //- Inertio-capillary time [s]
        static scalar inertioCapillaryTime
        (
            const scalar rho,
            const scalar sigma,
            const scalar diameter,
            const scalar coefficient
        );

        //- Viscous Rayleigh-Plateau breakup time [s]
        static scalar ligamentBreakupTime
        (
            const scalar rho,
            const scalar mu,
            const scalar sigma,
            const scalar diameter,
            const scalar coefficient
        );

        //- Advective sheet-breakup time [s]
        static scalar sheetBreakupTime
        (
            const scalar diameter,
            const scalar relativeSpeed,
            const scalar coefficient
        );


    // Constructors

        //- Construct from the film and model dictionary
        morphologyBasedAggregationScale
        (
            const liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~morphologyBasedAggregationScale() = default;


    // Member Functions

        //- Calculate the per-face physical event time and coherence reach
        virtual void calculate
        (
            const labelField& morphology,
            const scalarField& physicalDiameter,
            const vectorField& Up,
            const boolList& activeFace,
            scalarField& eventTime,
            scalarField& coherenceLength
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
