/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
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
    Foam::regionModels::areaSurfaceFilmModels::custom

Description
    Temporal kernel defined by an integrable Function1<scalar> release rate
    over normalized event age.

SourceFiles
    custom.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_temporalKernel_custom_H
#define regionFaModels_temporalKernel_custom_H

#include "temporalKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                           Class custom Declaration
\*---------------------------------------------------------------------------*/

class custom
:
    public temporalKernel
{
    // Private Data

        //- Release rate as a function of normalized event age
        autoPtr<Function1<scalar>> releaseRatePtr_;

        //- Integral of the release rate over normalized event age [0, 1]
        scalar totalIntegral_;


public:

    //- Runtime type information
    TypeName("custom");


    // Constructors

        //- Construct from dictionary
        explicit custom(const dictionary& dict);


    //- Destructor
    virtual ~custom() = default;


    // Member Functions

        //- Return the normalized cumulative release fraction
        virtual scalar cumulative
        (
            const scalar age,
            const scalar duration
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
