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
    Foam::regionModels::areaSurfaceFilmModels::parcelDiameterSpreading

Description
    Physical spreading based on the mapped mass-weighted parcel diameter.

Usage
    // Minimal example
    \verbatim
    physicalSpreading
    {
        // Mandatory entries
        type              parcelDiameterSpreading;

        // Optional entries
        beta                <Function1<scalar>>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property        | Description                | Type   | Reqd | Deflt
      type            | Type name: parcelDiameterSpreading  | word | yes | -
      beta  | Positive diameter multiplier     | Function1<scalar> | no  | 1
    \endtable

    The inherited entries are elaborated in:
      - \link physicalSpreading.H \endlink

Note
  - The circular radius is \c 0.5*beta*d. For an ellipse this is the
    minor radius and the major radius is \c aspectRatio times larger.

SourceFiles
    parcelDiameterSpreading.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_parcelDiameterSpreading_H
#define regionFaModels_parcelDiameterSpreading_H

#include "physicalSpreading.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class Function1;

namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                  Class parcelDiameterSpreading Declaration
\*---------------------------------------------------------------------------*/

class parcelDiameterSpreading
:
    public physicalSpreading
{
    // Private Data

        //- Time-dependent diameter multiplier
        autoPtr<Function1<scalar>> betaPtr_;


    // Private Member Functions

        //- Evaluate and validate the active footprint multiplier
        scalar beta() const;


public:

    //- Runtime type information
    TypeName("parcelDiameterSpreading");


    // Constructors

        //- Construct from mesh and dictionary
        parcelDiameterSpreading(const faMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~parcelDiameterSpreading();


    // Member Functions

        //- Spread all source channels using parcel-diameter footprints
        virtual void spread
        (
            areaScalarField& mass,
            areaVectorField& momentum,
            areaScalarField& normalImpulse,
            const areaScalarField& diameter
        ) const;

        //- Spread transfer channels and temporal duration moments together
        virtual void spread
        (
            areaScalarField& mass,
            areaVectorField& momentum,
            areaScalarField& normalImpulse,
            areaScalarField& temporalWeight,
            areaScalarField& temporalMoment,
            const areaScalarField& diameter
        ) const;

        //- Parcel-diameter spreading requires mapped parcel diameter
        virtual bool needsDiameter() const noexcept { return true; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
