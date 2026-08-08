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
    Foam::regionModels::areaSurfaceFilmModels::inertioCapillary

Description
    Inertio-capillary impact \c timeScaleModel, which is calculated as:

    \f[
        \tau = C_\sigma \sqrt{\rho \frac{d^3}{\sigma}}
    \f]

    where:
    \vartable
      \tau      | Temporal spreading duration [s]
      C_\sigma  | User-defined duration multiplier [-]
      d         | Mapped parcel diameter [m]
      \rho      | Parcel density [kg/m^3]
      \sigma    | Surface tension [kg/s^2]
    \endvartable

Usage
    Minimal example:
    \verbatim
    temporalSpreading
    {
        timeScaleModel  inertioCapillary;

        // Optional entries
        Csigma          <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                           | Type   | Reqd | Deflt
      timeScaleModel | Type name: inertioCapillary     | word   | yes  | -
      Csigma   | Inertio-capillary duration multiplier | scalar | no   | 1
    \endtable

SourceFiles
    inertioCapillary.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_timeScaleModel_inertioCapillary_H
#define regionFaModels_timeScaleModel_inertioCapillary_H

#include "timeScaleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class inertioCapillary Declaration
\*---------------------------------------------------------------------------*/

class inertioCapillary
:
    public timeScaleModel
{
    // Private Data

        //- Duration multiplier [-]
        const scalar Csigma_;


public:

    //- Runtime type information
    TypeName("inertioCapillary");


    // Constructors

        //- Construct from dictionary
        explicit inertioCapillary(const dictionary& dict);


    //- Destructor
    virtual ~inertioCapillary() = default;


    // Member Functions

        //- Calculate inertio-capillary physical durations for active faces
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
