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
    Foam::regionModels::areaSurfaceFilmModels::timeScaleModel

Description
    Abstract runtime-selectable interface for the \c characteristicTime
    temporal-spreading model.

    Its \c calculate() writes a finite positive duration for each active
    source face.

Usage
    Minimal example:
    \verbatim
    temporalSpreading
    {
        type            characteristicTime;

        // Mandatory entries
        timeScaleModel  <word>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property       | Description             | Type | Reqd | Deflt
      timeScaleModel | Physical duration model | word | yes  | -
    \endtable

    The inherited entries are elaborated in:
      - \link advective.H \endlink
      - \link inertioCapillary.H \endlink

    Custom characteristic-time models can override \c calculateAt to
    support independent durations for parcel impacts. The field-based
    \c calculate API remains available for aggregated source fields.

SourceFiles
    timeScaleModel.C
    timeScaleModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_timeScaleModel_H
#define regionFaModels_timeScaleModel_H

#include "areaFields.H"
#include "autoPtr.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                        Class timeScaleModel Declaration
\*---------------------------------------------------------------------------*/

class timeScaleModel
{
protected:

    // Protected Data

        //- Model dictionary for input diagnostics
        const dictionary dict_;


public:

    //- Runtime type information
    TypeName("timeScaleModel");


    // Declare runtime constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        timeScaleModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from model type and dictionary
        timeScaleModel(const word& modelType, const dictionary& dict);


    // Generated Methods

        //- No copy construct
        timeScaleModel(const timeScaleModel&) = delete;

        //- No copy assignment
        void operator=(const timeScaleModel&) = delete;


    // Selectors

        //- Select a physical duration model
        static autoPtr<timeScaleModel> New(const dictionary& dict);


    //- Destructor
    virtual ~timeScaleModel() = default;


    // Member Functions

        //- Set physical durations for active incoming source faces
        virtual void calculate
        (
            areaScalarField& duration,
            const areaScalarField& mass,
            const areaVectorField& momentum,
            const areaScalarField& normalImpulse,
            const areaScalarField& diameter,
            const areaScalarField& rho,
            const areaScalarField& sigma
        ) const = 0;

        //- Calculate a single independently timed source impact
        virtual scalar calculateAt
        (
            const label face,
            const scalar mass,
            const vector& momentum,
            const scalar normalImpulse,
            const scalar diameter,
            const scalar rho,
            const scalar sigma
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
