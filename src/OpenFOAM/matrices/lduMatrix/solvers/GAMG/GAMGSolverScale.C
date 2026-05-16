/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "GAMGSolver.H"
#include "FixedList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::GAMGSolver::scale
(
    solveScalarField& field,
    solveScalarField& Acf,
    const lduMatrix& A,
    const FieldField<Field, scalar>& interfaceLevelBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaceLevel,
    const solveScalarField& source,
    const direction cmpt
) const
{
    A.Amul
    (
        Acf,
        field,
        interfaceLevelBouCoeffs,
        interfaceLevel,
        cmpt
    );


    const label nCells = field.size();
    solveScalar* __restrict__ fieldPtr = field.begin();
    const solveScalar* const __restrict__ sourcePtr = source.begin();
    const solveScalar* const __restrict__ AcfPtr = Acf.begin();


    FixedList<solveScalar, 2> scalingFactor(Zero);

    const labelRange cellRange(nCells);

    #ifdef FOAM_OFFLOAD
    auto transform_op = [&](label i) -> FixedList<solveScalar, 2>
    {
        solveScalar v = fieldPtr[i];
        return FixedList<solveScalar, 2>
        {
            v * sourcePtr[i],
            v * AcfPtr[i]
        };
    };

    auto reduce_op = [](const FixedList<solveScalar, 2>& a,
        const FixedList<solveScalar, 2>& b) -> FixedList<solveScalar, 2> 
    {
        return FixedList<solveScalar, 2>
        {
            a[0] + b[0],
            a[1] + b[1]
        };
    };

    apply_policy(cellRange.size(), [&](auto policy)
    {
        scalingFactor = std::transform_reduce
        (
            policy,
            cellRange.begin(),
            cellRange.end(),
            FixedList<solveScalar, 2>(Zero),
            reduce_op,
            transform_op
        );
    });
    #else
    for (label i=0; i<nCells; i++)
    {
        scalingFactor[0] += fieldPtr[i]*sourcePtr[i];
        scalingFactor[1] += fieldPtr[i]*AcfPtr[i];
    }
    #endif

    A.mesh().reduce(scalingFactor, sumOp<solveScalar>());

    const solveScalar sf =
    (
        scalingFactor[0]
      / stabilise(scalingFactor[1], pTraits<solveScalar>::vsmall)
    );

    if (debug >= 2)
    {
        Pout<< sf << " ";
    }

    const scalarField& D = A.diag();
    const scalar* const __restrict__ DPtr = D.begin();

    apply_policy(cellRange.size(), [&](auto policy)
    {
        std::for_each
        (
            policy,
            cellRange.begin(),
            cellRange.end(),
            [&](const auto& i)
            {
                fieldPtr[i] =
                    sf*fieldPtr[i] + (sourcePtr[i] - sf*AcfPtr[i])/DPtr[i];
            }
        );
    });
}


// ************************************************************************* //
