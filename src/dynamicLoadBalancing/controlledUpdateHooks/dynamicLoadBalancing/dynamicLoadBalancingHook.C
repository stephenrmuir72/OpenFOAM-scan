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

#include "dynamicLoadBalancingHook.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace controlledUpdateHooks
{
    defineTypeNameAndDebug(dynamicLoadBalancing, 0);
    addToRunTimeSelectionTable(controlledUpdateHook, dynamicLoadBalancing, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::controlledUpdateHooks::dynamicLoadBalancing::dynamicLoadBalancing
(
    const word& name,
    fvMesh& mesh,
    const dictionary& dict
)
:
    controlledUpdateHook(name, mesh, dict),
    loadBalancer_(mesh, dict)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::controlledUpdateHooks::dynamicLoadBalancing::execute()
{
    mesh().clearMeshPhi();
    const bool balanced = loadBalancer_.balance();
    mesh().createMeshPhi();

    return balanced;
}


bool Foam::controlledUpdateHooks::dynamicLoadBalancing::read
(
    const dictionary& dict
)
{
    return controlledUpdateHook::read(dict) && loadBalancer_.read(dict);
}


// ************************************************************************* //
