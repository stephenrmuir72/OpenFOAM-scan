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
    Foam::regionModels::areaSurfaceFilmModels::footprintKernel

Description
    Abstract base class for the runtime-selectable distribution-kernel
    submodel of the \c spatialSpreading model.

    A footprint kernel scales a raw face-coverage measure, target face area
    or overlap area, by a distance-dependent weighting factor before
    per-source normalization.

Usage
    Minimal example:
    spatialSpreading
    {
        // Optional entries
        footprintKernel    <word>;
    }
    \endverbatim

    Options for the \c faceCoverage entry:
    \verbatim
      uniform   | Use uniform distribution (default)
      gaussian  | Use the standard Gaussian distribution
    \endverbatim

SourceFiles
    footprintKernel.C
    footprintKernelNew.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_footprintKernel_H
#define regionFaModels_footprintKernel_H

#include "autoPtr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class footprintKernel Declaration
\*---------------------------------------------------------------------------*/

class footprintKernel
{
protected:

    // Protected Data

        //- Model dictionary for run-time input diagnostics
        const dictionary dict_;


public:

    //- Runtime type information
    TypeName("footprintKernel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            footprintKernel,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Constructors

        //- Construct from model type and dictionary
        footprintKernel(const word& modelType, const dictionary& dict);


    // Generated Methods

        //- No copy construct
        footprintKernel(const footprintKernel&) = delete;

        //- No copy assignment
        void operator=(const footprintKernel&) = delete;


    // Selectors

        //- Select a footprint-kernel model
        static autoPtr<footprintKernel> New(const dictionary& dict);


    //- Destructor
    virtual ~footprintKernel() = default;


    // Member Functions

        //- Weight a raw coverage measure using the kernel evaluated at the
        //- given distance from the source footprint
        virtual scalar weight
        (
            const scalar coverageMeasure,
            const scalar kernelDistance
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
