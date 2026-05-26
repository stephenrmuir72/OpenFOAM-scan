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
    Foam::controlledUpdateHooks::dynamicLoadBalancing

Description
    Controlled-update hook that executes Foam::dynamicLoadBalancer.

SourceFiles
    dynamicLoadBalancingHook.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_controlledUpdateHooks_dynamicLoadBalancing_H
#define Foam_controlledUpdateHooks_dynamicLoadBalancing_H

#include "controlledUpdateHook.H"
#include "dynamicLoadBalancer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace controlledUpdateHooks
{

/*---------------------------------------------------------------------------*\
                    Class dynamicLoadBalancing Declaration
\*---------------------------------------------------------------------------*/

class dynamicLoadBalancing
:
    public controlledUpdateHook
{
    // Private Data

        //- Dynamic load-balancing orchestrator
        dynamicLoadBalancer loadBalancer_;


public:

    //- Runtime type information
    TypeName("dynamicLoadBalancing");


    // Constructors

        //- Construct from hook name, mesh and dictionary
        dynamicLoadBalancing
        (
            const word& name,
            fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~dynamicLoadBalancing() = default;


    // Member Functions

        //- Execute the hook and return mesh-changed status
        virtual bool execute();

        //- Read settings from a dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace controlledUpdateHooks
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
