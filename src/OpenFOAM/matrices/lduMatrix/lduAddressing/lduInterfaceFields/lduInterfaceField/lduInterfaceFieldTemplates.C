/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::lduInterfaceField::addToInternalField
(
    Field<Type>& result,
    const bool add,
    const labelUList& faceCells,
    const scalarField& coeffs,
    const Field<Type>& vals
) const
{
    if (add)
    {
        forAll(faceCells, elemi)
        {
            result[faceCells[elemi]] += coeffs[elemi]*vals[elemi];
        }
    }
    else
    {
        forAll(faceCells, elemi)
        {
            result[faceCells[elemi]] -= coeffs[elemi]*vals[elemi];
        }
    }
}

template<class Type>
void Foam::lduInterfaceField::addToInternalField
(
    Field<Type>& result,
    const bool add,
    const lduAddressing& addr,
    const label patchId,
    const scalarField& coeffs,
    const Field<Type>& vals
) const
{
    const labelUList& patchCells = addr.patchSortCells(patchId);
    const labelUList& sort       = addr.patchSortAddr(patchId);
    const labelUList& sortStart  = addr.patchSortStartAddr(patchId);
    
    const scalar sign = add ? scalar(1) : scalar(-1);
    const labelRange cellRange(patchCells.size());

    auto neiFun =
        [sign, c = coeffs.cdata(), v = vals.cdata()]
        (const label /*cellIdx*/, const label facei) -> Type
        {
            return sign * c[facei] * v[facei];
        };

    // Accumulation into result cell
    auto neiOp =
        [](Type& res, const Type& val)
        {
            res += val;
        };

    apply_policy(cellRange.size(), [&](auto policy)
    {
        std::for_each
        (
            policy,
            cellRange.begin(),
            cellRange.end(),
            patchCellReduceFunctor
            <
                Type,
                decltype(neiFun),
                decltype(neiOp)
            >
            (
                result,
                sortStart,
                sort,
                patchCells,
                neiFun,
                neiOp
            )
        );
    });
}

// ************************************************************************* //
