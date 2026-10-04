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
    Foam::regionModels::areaSurfaceFilmModels::characteristicTime

Description
    Temporal spreading model whose duration is selected by \c timeScaleModel.

Usage
    Minimal example:
    \verbatim
    temporalSpreading
    {
        // Mandatory entries
        type            characteristicTime;
        timeScaleModel  <word>;

        // Inherited entries
        ...
    }
    \endverbatim

    \table
      Property       | Description                     | Type | Reqd | Deflt
      type           | Type name: characteristicTime   | word | yes  | -
      timeScaleModel | Physical duration model         | word | yes  | -
    \endtable

SourceFiles
    characteristicTime.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_characteristicTime_H
#define regionFaModels_characteristicTime_H

#include "temporalSpreading.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// Forward Declarations
class timeScaleModel;

/*---------------------------------------------------------------------------*\
                   Class characteristicTime Declaration
\*---------------------------------------------------------------------------*/

class characteristicTime
:
    public temporalSpreading
{
    // Private Data

        //- Selected physical impact-time model
        autoPtr<timeScaleModel> timeScalePtr_;


protected:

    // Protected Member Functions

        //- Calculate durations using the selected physical time-scale model
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

        virtual scalar calculateDurationAt
        (
            const label face,
            const scalar mass,
            const vector& momentum,
            const scalar normalImpulse,
            const scalar diameter,
            const scalar rho,
            const scalar sigma
        ) const;


public:

    //- Runtime type information
    TypeName("characteristicTime");


    // Constructors

        //- Construct from mesh and dictionary
        characteristicTime(const faMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~characteristicTime();


    // Member Functions

        //- Characteristic-time models require mapped parcel diameter
        virtual bool needsDiameter() const noexcept { return true; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
