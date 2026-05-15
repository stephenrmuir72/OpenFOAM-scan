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
    Foam::regionModels::areaSurfaceFilmModels::gaussian

Description
    Gaussian distribution kernel.

    Multiplies the coverage measure by
    \f$ \exp(-0.5 (d/\sigma)^2) \f$
    where \c d is the kernel-evaluation distance and \c sigma is an
    absolute standard deviation.

Usage
    Minimal example:
    \verbatim
    physicalSpreading
    {
        // Optional entries
        footprintKernel     gaussian;

        // Conditional mandatory entries
        sigma               <scalar>;
    }
    \endverbatim

    \table
      Property | Description                           | Type | Reqd | Deflt
      footprintKernel | Type name: gaussian            | word | yes  | -
      sigma    | Absolute Gaussian standard deviation  | scalar | yes | -
    \endtable

SourceFiles
    gaussian.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_footprintKernel_gaussian_H
#define regionFaModels_footprintKernel_gaussian_H

#include "footprintKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                          Class gaussian Declaration
\*---------------------------------------------------------------------------*/

class gaussian
:
    public footprintKernel
{
    // Private Data

        //- Absolute Gaussian standard deviation
        const scalar sigma_;


public:

    //- Runtime type information
    TypeName("gaussian");


    // Constructors

        //- Construct from dictionary
        explicit gaussian(const dictionary& dict);


    //- Destructor
    virtual ~gaussian() = default;


    // Member Functions

        //- Return the coverage measure scaled by the Gaussian kernel
        virtual scalar weight
        (
            const scalar coverageMeasure,
            const scalar kernelDistance
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
