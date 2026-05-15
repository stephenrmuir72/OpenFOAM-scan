/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 OpenCFD Ltd.
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
    Foam::regionModels::areaSurfaceFilmModels::sourceSpreadingModel

Description
    Abstract base class for spatial spreading of Lagrangian-to-film source
    terms (mass, momentum, normal pressure impulse) over neighbouring
    finite-area faces.

    When a Lagrangian parcel is absorbed into the film, its source terms
    are initially mapped one-to-one onto a single finite-area face.  For
    large super-parcels this can create impulsive local sources that
    drive pressure gradients and destabilise the simulation.  Derived
    classes redistribute (spread) these sources to a neighbourhood of
    faces while preserving conservation.

    The spreading is applied to the raw impulse fields (mass [kg],
    momentum [kg m/s], normal impulse [kg m/s]) after
    \c volSurfaceMapping::mapToSurface() and before the rate-per-area
    conversion in \c liquidFilmModel::preEvolveRegion().

Usage
    Configured via a \c sourceSpreading sub-dictionary in the film model
    dictionary:
    \verbatim
    sourceSpreading
    {
        type        none;   // none | parcelDiameter | gaussian
    }
    \endverbatim

SourceFiles
    sourceSpreadingModel.C
    sourceSpreadingModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_sourceSpreadingModel_H
#define regionFaModels_sourceSpreadingModel_H

#include "faMesh.H"
#include "scalarField.H"
#include "vectorField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class sourceSpreadingModel Declaration
\*---------------------------------------------------------------------------*/

class sourceSpreadingModel
{
    // Private Member Functions

        //- No copy construct
        sourceSpreadingModel(const sourceSpreadingModel&) = delete;

        //- No copy assignment
        void operator=(const sourceSpreadingModel&) = delete;


protected:

    // Protected Data

        //- Reference to the finite-area mesh
        const faMesh& mesh_;


public:

    //- Runtime type information
    TypeName("sourceSpreadingModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             sourceSpreadingModel,
             dictionary,
             (
                const faMesh& mesh,
                const dictionary& dict
             ),
             (mesh, dict)
         );


    // Constructors

        //- Construct from mesh and dictionary
        sourceSpreadingModel
        (
            const word& modelType,
            const faMesh& mesh,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected model
        static autoPtr<sourceSpreadingModel> New
        (
            const faMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sourceSpreadingModel() = default;


    // Member Functions

        //- Spread a scalar source field in place.
        //  The diameter field provides the mass-weighted average parcel
        //  diameter at each source face (used by some models).
        virtual void spread
        (
            scalarField& field,
            const scalarField& diameter
        ) const = 0;

        //- Spread a vector source field in place.
        //  The diameter field provides the mass-weighted average parcel
        //  diameter at each source face (used by some models).
        virtual void spread
        (
            vectorField& field,
            const scalarField& diameter
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
