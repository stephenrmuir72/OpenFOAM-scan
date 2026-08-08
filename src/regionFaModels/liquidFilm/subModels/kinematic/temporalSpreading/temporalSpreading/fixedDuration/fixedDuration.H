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
    Foam::regionModels::areaSurfaceFilmModels::fixedDuration

Description
    Temporal spreading model with a configured \c duration [s].

Usage
    Minimal example:
    \verbatim
    temporalSpreading
    {
        // Mandatory entries
        type            fixedDuration;
        duration        <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    \table
      Property | Description                 | Type       | Reqd | Deflt
      type     | Type name: fixedDuration    | word       | yes  | -
      duration | Fixed physical duration [s] | scalar     | yes  | -
    \endtable

SourceFiles
    fixedDuration.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_fixedDuration_H
#define regionFaModels_fixedDuration_H

#include "temporalSpreading.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                        Class fixedDuration Declaration
\*---------------------------------------------------------------------------*/

class fixedDuration
:
    public temporalSpreading
{
    // Private Data

        //- Physical impact duration [s]
        const scalar duration_;


protected:

    // Protected Member Functions

        //- Set the configured physical duration on every face
        virtual void calculateDuration
        (
            areaScalarField& duration,
            const areaScalarField& mass,
            const areaVectorField& momentum,
            const areaScalarField& normalImpulse,
            const areaScalarField& diameter,
            const areaScalarField& rho,
            const areaScalarField& sigma
        ) const;


public:

    //- Runtime type information
    TypeName("fixedDuration");


    // Constructors

        //- Construct from mesh and dictionary
        fixedDuration(const faMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~fixedDuration() = default;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
