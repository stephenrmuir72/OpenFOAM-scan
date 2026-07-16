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
    Foam::selectivePairGAMGAgglomeration

Description
    Pair agglomeration restricted to strong internal matrix connections.

    The non-negative connection weight is the magnitude of the upper
    coefficient for symmetric matrices, or the maximum magnitude of the upper
    and lower coefficients for asymmetric matrices. An internal face is strong
    when its positive weight is at least \c strengthThreshold times the maximum
    positive incident weight of either endpoint. The default threshold is 0.25.

    This remains aggregation with piecewise-constant transfer operators; it is
    not a C/F splitting or interpolation scheme. Processor and coupled-interface
    edges do not participate in strength selection. Their existing restriction
    and optional processor-agglomeration handling are unchanged.

Usage
    \table
        Property             | Description                    | Required | Default
        agglomerator         | Agglomerator type selectivePair | yes      |
        strengthThreshold    | Strong-connection fraction [0,1]| no       | 0.25
    \endtable

SourceFiles
    selectivePairGAMGAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef selectivePairGAMGAgglomeration_H
#define selectivePairGAMGAgglomeration_H

#include "pairGAMGAgglomeration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class selectivePairGAMGAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class selectivePairGAMGAgglomeration
:
    public pairGAMGAgglomeration
{

public:

    //- Runtime type information
    TypeName("selectivePair");


    // Constructors

        //- Construct from matrix and controls
        selectivePairGAMGAgglomeration
        (
            const lduMatrix& matrix,
            const dictionary& controlDict
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
