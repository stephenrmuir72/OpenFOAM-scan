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

#include "selectivePairGAMGAgglomeration.H"
#include "lduMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(selectivePairGAMGAgglomeration, 0);

    addToRunTimeSelectionTable
    (
        GAMGAgglomeration,
        selectivePairGAMGAgglomeration,
        lduMatrix
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::selectivePairGAMGAgglomeration::selectivePairGAMGAgglomeration
(
    const lduMatrix& matrix,
    const dictionary& controlDict
)
:
    pairGAMGAgglomeration(matrix.mesh(), controlDict)
{
    const scalar strengthThreshold
    (
        controlDict.getOrDefault<scalar>("strengthThreshold", 0.25)
    );

    if (!(strengthThreshold >= 0 && strengthThreshold <= 1))
    {
        FatalIOErrorInFunction(controlDict)
            << "strengthThreshold " << strengthThreshold
            << " is outside the range [0,1]"
            << exit(FatalIOError);
    }

    if (matrix.hasLower())
    {
        agglomerateSelective
        (
            nCellsInCoarsestLevel_,
            0,
            max(mag(matrix.upper()), mag(matrix.lower())),
            strengthThreshold,
            true
        );
    }
    else
    {
        agglomerateSelective
        (
            nCellsInCoarsestLevel_,
            0,
            mag(matrix.upper()),
            strengthThreshold,
            true
        );
    }
}


// ************************************************************************* //
