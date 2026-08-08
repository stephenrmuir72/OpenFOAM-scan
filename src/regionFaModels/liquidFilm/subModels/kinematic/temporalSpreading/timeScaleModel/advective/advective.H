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
    Foam::regionModels::areaSurfaceFilmModels::advective

Description
    Advective impact \c timeScaleModel, which is calculated as:

    \f[
        \tau = \frac{C_a \, d}{|U_n|}
    \f]

    with

    \f[
        |U_n|=\frac{|I_n|}{|m|}
    \f]

    where:
    \vartable
      \tau      | Temporal spreading duration [s]
      C_a       | User-defined duration multiplier [-]
      d         | Mapped parcel diameter [m]
      U_n       | Normal velocity [m/s]
      I_n       | Normal impulse [kg m/s]
      m         | Parcel mass [kg]
    \endvartable

Usage
    Minimal example:
    \verbatim
    temporalSpreading
    {
        timeScaleModel  advective;

        // Optional entries
        Ca              <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                   | Type   | Reqd | Deflt
      timeScaleModel | Type name: advective    | word   | yes  | -
      Ca       | Advective duration multiplier | scalar | no   | 1
    \endtable

SourceFiles
    advective.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_timeScaleModel_advective_H
#define regionFaModels_timeScaleModel_advective_H

#include "timeScaleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                         Class advective Declaration
\*---------------------------------------------------------------------------*/

class advective
:
    public timeScaleModel
{
    // Private Data

        //- Duration multiplier [-]
        const scalar Ca_;


public:

    //- Runtime type information
    TypeName("advective");


    // Constructors

        //- Construct from dictionary
        explicit advective(const dictionary& dict);


    //- Destructor
    virtual ~advective() = default;


    // Member Functions

        //- Calculate advective physical durations for active source faces
        virtual void calculate
        (
            areaScalarField& duration,
            const areaScalarField& mass,
            const areaVectorField& momentum,
            const areaScalarField& normalImpulse,
            const areaScalarField& diameter,
            const areaScalarField& rho,
            const areaScalarField& sigma
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
