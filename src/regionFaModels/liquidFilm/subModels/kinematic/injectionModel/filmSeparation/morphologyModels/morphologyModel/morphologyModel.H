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
    Foam::regionModels::areaSurfaceFilmModels::morphologyModel

Description
    Physical liquid-structure classifier for separated finite-area film.
    Morphology determines the downstream physical droplet-diameter model.
    The state IDs are \c directDroplet=0, \c drip=1, \c ligament=2, and
    \c sheet=3.

Usage
    \verbatim
    filmSeparationCoeffs
    {
        // Optional entries
        morphology
        {
            type  <word>;
            ...
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                   | Type | Reqd | Deflt
      morphology   | Morphology settings           | dict | no   | -
    \endtable

SourceFiles
    morphologyModel.C
    morphologyModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_morphologyModel_H
#define regionFaModels_morphologyModel_H

#include "autoPtr.H"
#include "boolList.H"
#include "dictionary.H"
#include "Enum.H"
#include "labelField.H"
#include "liquidFilmBase.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                       Class morphologyModel Declaration
\*---------------------------------------------------------------------------*/

class morphologyModel
{
public:

    // Public Enumerations

        //- Initial liquid structures produced by film separation
        enum morphologyType : char
        {
            DIRECT_DROPLET = 0,
            DRIP,
            LIGAMENT,
            SHEET,
            nMorphologyTypes
        };

        //- Names for morphologyType
        static const Enum<morphologyType> morphologyTypeNames;


protected:

    // Protected Data

        //- Const reference to the liquid film model
        const liquidFilmBase& film_;

        //- Const reference to the finite-area mesh
        const faMesh& mesh_;


public:

    //- Runtime type information
    TypeName("morphologyModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            morphologyModel,
            dictionary,
            (
                const liquidFilmBase& film,
                const dictionary& dict
            ),
            (film, dict)
        );


    // Constructors

        //- Construct from model type, film, and dictionary
        morphologyModel
        (
            const word& modelType,
            const liquidFilmBase& film,
            const dictionary& dict
        );


    // Generated Methods

        //- No copy construct
        morphologyModel(const morphologyModel&) = delete;

        //- No copy assignment
        void operator=(const morphologyModel&) = delete;


    // Selectors

        //- Select a morphology model from a morphology dictionary
        static autoPtr<morphologyModel> New
        (
            const liquidFilmBase& film,
            const dictionary& dict
        );

        //- Select a morphology model if the film dictionary supplies one
        static autoPtr<morphologyModel> NewIfPresent
        (
            const liquidFilmBase& film,
            const dictionary& filmDict
        );


    //- Destructor
    virtual ~morphologyModel() = default;


    // Member Functions

        //- Return the local morphology ID for every finite-area face
        virtual tmp<labelField> morphology() const = 0;

        //- Return the morphology IDs that this model can produce
        virtual boolList possibleTypes() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
