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
    Foam::regionModels::areaSurfaceFilmModels::fixedMorphology

Description
    Deterministic morphology modeld.

Usage
    \verbatim
    morphology
    {
        // Mandatory entries
        type            fixed;

        fixedCoeffs
        {
            morphology ligament;
        }
    }
    \endverbatim

    The selected morphology is returned for every finite-area face.  The
    model can produce only the configured morphology type.

SourceFiles
    fixedMorphology.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_fixedMorphology_H
#define regionFaModels_fixedMorphology_H

#include "morphologyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                        Class fixedMorphology Declaration
\*---------------------------------------------------------------------------*/

class fixedMorphology
:
    public morphologyModel
{
    // Private Data

        //- Selected morphology
        morphologyType morphology_;


public:

    //- Runtime type information
    TypeName("fixed");


    // Constructors

        //- Construct from the film and model dictionary
        fixedMorphology
        (
            const liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fixedMorphology() = default;


    // Member Functions

        //- Return the configured morphology for every film face
        virtual tmp<labelField> morphology() const;

        //- Return the configured morphology as the only possible type
        virtual boolList possibleTypes() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
