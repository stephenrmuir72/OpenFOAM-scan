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
    Foam::regionModels::areaSurfaceFilmModels::topHat

Description
    Constant-rate temporal kernel. Its cumulative fraction is:

    \f[
        F(a) = \min(1, \max(0, \frac{a}{\tau}))
    \f]

    where:
    \vartable
      \tau      | Duration [s]
      a         | Age [s]
      d         | Mapped parcel diameter [m]
      F(a)      | Cumulative fraction [-]
    \endvartable

    No Lagrangian-to-film transfer is released before the event starts,
    release is uniform during the physical duration, and the cumulative
    fraction is one at the event end.

Usage
    Minimal example:
    \verbatim
    temporalSpreading
    {
        // Inherited entries
        ...

        // Optional entries
        temporalKernel      topHat;
    }
    \endverbatim

    \table
      Property       | Description                | Type | Reqd | Deflt
      temporalKernel | Type name: topHat          | word | no   | topHat
    \endtable

SourceFiles
    topHat.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_temporalKernel_topHat_H
#define regionFaModels_temporalKernel_topHat_H

#include "temporalKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                            Class topHat Declaration
\*---------------------------------------------------------------------------*/

class topHat
:
    public temporalKernel
{
public:

    //- Runtime type information
    TypeName("topHat");


    // Constructors

        //- Construct from dictionary
        explicit topHat(const dictionary& dict);


    //- Destructor
    virtual ~topHat() = default;


    // Member Functions

        //- Return the constant-rate cumulative release fraction
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
