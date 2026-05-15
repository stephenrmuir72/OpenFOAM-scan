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
    Foam::regionModels::areaSurfaceFilmModels::noSpreading

Description
    No source spreading.  Source terms are left unchanged on the single
    finite-area face to which they were mapped from the Lagrangian parcel.

    This is the default model and reproduces the original (pre-spreading)
    behaviour.

Usage
    \verbatim
    sourceSpreading
    {
        type    none;
    }
    \endverbatim

SourceFiles
    noSpreading.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_noSpreading_H
#define regionFaModels_noSpreading_H

#include "sourceSpreadingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                         Class noSpreading Declaration
\*---------------------------------------------------------------------------*/

class noSpreading
:
    public sourceSpreadingModel
{
    // Private Member Functions

        //- No copy construct
        noSpreading(const noSpreading&) = delete;

        //- No copy assignment
        void operator=(const noSpreading&) = delete;


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from mesh and dictionary
        noSpreading(const faMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~noSpreading() = default;


    // Member Functions

        //- No-op: leave scalar field unchanged
        virtual void spread
        (
            scalarField& field,
            const scalarField& diameter
        ) const;

        //- No-op: leave vector field unchanged
        virtual void spread
        (
            vectorField& field,
            const scalarField& diameter
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
