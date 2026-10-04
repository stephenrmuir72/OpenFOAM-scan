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
    Foam::regionModels::areaSurfaceFilmModels::temporalKernel

Description
    Abstract runtime-selectable interface for the normalized cumulative release
    function

    The kernel receives 'event age' and 'physical duration' and returns the
    cumulative fraction in [0,1] released by that age. The \c temporalSpreading
    model uses differences in this cumulative fraction to release the event
    over each time interval.

Usage
    Minimal example:
    \verbatim
    temporalSpreading
    {
        // Inherited entries
        ...

        // Optional entries
        temporalKernel    <word>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property       | Description                  | Type | Reqd | Deflt
      temporalKernel | Cumulative release kernel    | word | no   | topHat
    \endtable

    Options for the \c temporalKernel entry:
    \verbatim
      topHat        | Constant-rate temporal kernel
      custom        | User-defined release-rate function
    \endverbatim

    The \c custom kernel uses a \c Function1<scalar> named \c releaseRate,
    evaluated as a function of normalized age \f$x = a/\tau\f$. Its cumulative
    release fraction is the integral of the release rate, normalized by its
    integral over \f$x \in [0,1]\f$. The release rate must be non-negative and
    its selected Function1 type must implement integration.

    The former symmetric triangular kernel is represented by this piecewise
    linear release-rate table:
    \verbatim
    temporalSpreading
    {
        temporalKernel  custom;

        releaseRate
        {
            type    table;
            values
            (
                (0    0)
                (0.5  2)
                (1    0)
            );
        }
    }
    \endverbatim

SourceFiles
    temporalKernel.C
    temporalKernelNew.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_temporalKernel_H
#define regionFaModels_temporalKernel_H

#include "autoPtr.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type> class Function1;

namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                        Class temporalKernel Declaration
\*---------------------------------------------------------------------------*/

class temporalKernel
{
protected:

    // Protected Data

        //- Model dictionary for input diagnostics
        const dictionary dict_;


public:

    //- Runtime type information
    TypeName("temporalKernel");


    // Declare runtime constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        temporalKernel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from model type and dictionary
        temporalKernel(const word& modelType, const dictionary& dict);


    // Generated Methods

        //- No copy construct
        temporalKernel(const temporalKernel&) = delete;

        //- No copy assignment
        void operator=(const temporalKernel&) = delete;


    // Selectors

        //- Select a temporal kernel
        static autoPtr<temporalKernel> New(const dictionary& dict);


    //- Destructor
    virtual ~temporalKernel() = default;


    // Member Functions

        //- Return the normalized cumulative release at the supplied age
        virtual scalar cumulative
        (
            const scalar age,
            const scalar duration
        ) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
