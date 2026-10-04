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
    Foam::regionModels::areaSurfaceFilmModels::uniform

Description
    Uniform distribution kernel.

    Returns the coverage measure unchanged, i.e. every retained candidate
    face is weighted only by its coverage measure, target face area or
    overlap area for area-intersection coverage.

Usage
    Minimal example:
    \verbatim
    spatialSpreading
    {
        // Optional entries
        footprintKernel     uniform;
    }
    \endverbatim

    \table
      Property | Description                           | Type | Reqd | Deflt
      footprintKernel | Type name: uniform             | word | yes  | -
    \endtable

SourceFiles
    uniform.cxx

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_footprintKernel_uniform_H
#define regionFaModels_footprintKernel_uniform_H

#include "footprintKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                           Class uniform Declaration
\*---------------------------------------------------------------------------*/

class uniform
:
    public footprintKernel
{
public:

    //- Runtime type information
    TypeName("uniform");


    // Constructors

        //- Construct from dictionary
        explicit uniform(const dictionary& dict);


    //- Destructor
    virtual ~uniform() = default;


    // Member Functions

        //- Return the coverage measure unchanged
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
