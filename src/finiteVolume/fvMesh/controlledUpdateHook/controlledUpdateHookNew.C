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

\*---------------------------------------------------------------------------*/

#include "controlledUpdateHook.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::controlledUpdateHook> Foam::controlledUpdateHook::New
(
    const word& name,
    fvMesh& mesh,
    const dictionary& dict
)
{
    const word hookType(dict.getWord("type"));

    Info<< "Selecting controlled-update hook '" << name
        << "' of type " << hookType << endl;

    if (!dictionaryConstructorTablePtr_)
    {
        FatalErrorInFunction
            << "controlledUpdateHook table is empty" << exit(FatalError);
    }

    auto* ctorPtr = dictionaryConstructorTable(hookType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "type",
            hookType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<controlledUpdateHook>(ctorPtr(name, mesh, dict));
}


// ************************************************************************* //
