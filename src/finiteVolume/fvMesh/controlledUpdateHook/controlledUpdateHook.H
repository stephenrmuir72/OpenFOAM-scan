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
    Foam::controlledUpdateHook

Description
    Runtime-selectable hook executed from dynamicFvMesh::controlledUpdate().

SourceFiles
    controlledUpdateHook.C
    controlledUpdateHookNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_controlledUpdateHook_H
#define Foam_controlledUpdateHook_H

#include "autoPtr.H"
#include "className.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                       Class controlledUpdateHook Declaration
\*---------------------------------------------------------------------------*/

class controlledUpdateHook
{
    // Private Data

        //- Hook entry name
        word name_;

        //- Finite-volume mesh
        fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("controlledUpdateHook");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            controlledUpdateHook,
            dictionary,
            (
                const word& name,
                fvMesh& mesh,
                const dictionary& dict
            ),
            (name, mesh, dict)
        );


    // Selectors

        //- Return a reference to the selected controlled-update hook
        static autoPtr<controlledUpdateHook> New
        (
            const word& name,
            fvMesh& mesh,
            const dictionary& dict
        );


    // Constructors

        //- Construct from hook name, mesh and dictionary
        controlledUpdateHook
        (
            const word& name,
            fvMesh& mesh,
            const dictionary& dict
        );


    // Generated Methods

        //- No copy construct
        controlledUpdateHook(const controlledUpdateHook&) = delete;

        //- No copy assignment
        void operator=(const controlledUpdateHook&) = delete;


    //- Destructor
    virtual ~controlledUpdateHook() = default;


    // Member Functions

    // Access

        //- Return hook name
        const word& name() const noexcept
        {
            return name_;
        }

        //- Return mesh
        fvMesh& mesh() const noexcept
        {
            return mesh_;
        }


    // Evaluation

        //- Execute the hook and return mesh-changed status
        virtual bool execute() = 0;


    // I/O

        //- Read settings from a dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
