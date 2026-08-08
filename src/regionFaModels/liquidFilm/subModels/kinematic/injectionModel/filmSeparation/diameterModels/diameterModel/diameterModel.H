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
    Foam::regionModels::areaSurfaceFilmModels::diameterModel

Description
    Run-time-selectable physical droplet-diameter model for separated
    finite-area film.

Usage
    Minimal example:
    \verbatim
    filmSeparationCoeffs
    {
        // Optional entries
        diameterModel  <word>;
        dMin           <scalar>;
        dMax           <scalar>;

        // Inherited entries
        ...
    }

SourceFiles
    diameterModel.C
    diameterModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_diameterModel_H
#define regionFaModels_diameterModel_H

#include "autoPtr.H"
#include "dictionary.H"
#include "liquidFilmBase.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                        Class diameterModel Declaration
\*---------------------------------------------------------------------------*/

class diameterModel
{
protected:

    // Protected Data

        //- Const reference to the liquid film model
        const liquidFilmBase& film_;

        //- Const reference to the finite-area mesh
        const faMesh& mesh_;

        //- Minimum admissible physical droplet diameter [m]
        scalar dMin_;

        //- Maximum admissible physical droplet diameter [m]
        scalar dMax_;


    // Protected Member Functions

        //- Clip a per-face diameter field to [dMin_, dMax_]
        void clip(scalarField& d) const
        {
            clip(d, dMin_, dMax_);
        }


public:

    //- Runtime type information
    TypeName("diameterModel");


    // Static Member Functions

        //- Clip a per-face diameter field to [dMin, dMax]
        static void clip
        (
            scalarField& d,
            const scalar dMin,
            const scalar dMax
        );


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            diameterModel,
            dictionary,
            (
                const liquidFilmBase& film,
                const dictionary& dict
            ),
            (film, dict)
        );


    // Constructors

        //- Construct from model type, film, and dictionary
        diameterModel
        (
            const word& modelType,
            const liquidFilmBase& film,
            const dictionary& dict
        );


    // Generated Methods

        //- No copy construct
        diameterModel(const diameterModel&) = delete;

        //- No copy assignment
        void operator=(const diameterModel&) = delete;


    // Selectors

        //- Select using the legacy diameterModel entry
        static autoPtr<diameterModel> New
        (
            const liquidFilmBase& film,
            const dictionary& dict
        );

        //- Select an explicitly named model
        static autoPtr<diameterModel> New
        (
            const liquidFilmBase& film,
            const dictionary& dict,
            const word& modelType
        );


    //- Destructor
    virtual ~diameterModel() = default;


    // Member Functions

        //- Return const access to the film properties
        const liquidFilmBase& film() const noexcept { return film_; }

        //- Return const access to the finite-area mesh
        const faMesh& mesh() const noexcept { return mesh_; }

        //- Return the per-face physical droplet diameter [m]
        virtual tmp<scalarField> diameter() const = 0;

        //- Preserve legacy direct-transfer behaviour where requested
        virtual bool scalesWithSeparationRatio() const noexcept
        {
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
