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
    Foam::regionModels::areaSurfaceFilmModels::triangular

Description
    Symmetric triangular temporal kernel with a release rate that rises to
    the midpoint and then falls. Its cumulative fraction is:

    \f[
        F(x) = 2\,x^2 \, \, \text{if } x \leq 0.5 \\
               1 - 2(1 - x)^2 \, \, \text{else}
    \f]

    with

    \f[
        x = \frac{a}{\tau}
    \f]

    where:
    \vartable
      F(x)     | Cumulative fraction [-]
      x        | Normalised age [-]
      a        | Age [s]
      \tau     | Physical duration [s]
    \endvartable

Usage
    Minimal example:
    \verbatim
    temporalSpreading
    {
        // Optional entries
        temporalKernel    triangular;
    }
    \endverbatim

    \table
      Property       | Description                   | Type | Reqd | Deflt
      temporalKernel | Type name: triangular         | word | no   | topHat
    \endtable

SourceFiles
    triangular.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_temporalKernel_triangular_H
#define regionFaModels_temporalKernel_triangular_H

#include "temporalKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                         Class triangular Declaration
\*---------------------------------------------------------------------------*/

class triangular
:
    public temporalKernel
{
public:

    //- Runtime type information
    TypeName("triangular");


    // Constructors

        //- Construct from dictionary
        explicit triangular(const dictionary& dict);


    //- Destructor
    virtual ~triangular() = default;


    // Member Functions

        //- Return the symmetric triangular cumulative release fraction
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
