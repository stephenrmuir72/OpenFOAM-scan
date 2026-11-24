/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2009 TU Delft
    Copyright (C) 2009 Dubravko Matijasevic, FSB Zagreb
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

#include "RBFInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, class CtrlField>
Foam::tmp<Foam::Field<Type>>
Foam::RBFInterpolation::interpolate_impl(const CtrlField& ctrlField) const
{
    #ifdef FULLDEBUG
    if (ctrlField.size() != controlPoints_.size())
    {
        FatalErrorInFunction
            << "Incorrect size of source field."
            << " Size = " << ctrlField.size()
            << " nControlPoints = " << controlPoints_.size() << endl
            << abort(FatalError);
    }
    #endif

    auto tresult = tmp<Field<Type>>::New(interpPoints_.size(), Foam::zero{});
    auto& result = tresult.ref();

    // Frank Bos (2008-12-21) :
    // 1) Calculate alpha and beta coefficients using the Inverse
    // 2) Calculate displacements of internal nodes using RBF values,
    //    alpha's and beta's
    // 3) Return displacements using tresult

    const label nControlPoints = controlPoints_.size();

    // The interpolation matrix: size == nControlPoints + (0 or 4)
    const auto& mat = this->interpMatrix();

    // Number of polynomial coefficients (0 or 4)
    const int nPolyCoeffs = static_cast<int>(mat.n() - nControlPoints);

    if (FOAM_UNLIKELY(nPolyCoeffs && nPolyCoeffs != 4))
    {
        FatalErrorInFunction
            << "Programming issue, matrix dimensioned with "
            << nPolyCoeffs << " polynomials instead of 4" << endl
            << abort(FatalError);
    }

    // Determine interpolation coefficients
    Field<Type> alpha(nControlPoints, Foam::zero{});
    FixedList<Type, 4> beta(Foam::zero{});

    constexpr bool packedInterpolation = true;

    if constexpr (!packedInterpolation)
    {
        // Simple element-wise coefficient evaluation
        // instead of the packed-triangle sweep

        for (label row = 0; row < nControlPoints; ++row)
        {
            for (label col = 0; col < nControlPoints; ++col)
            {
                alpha[row] += mat(row, col)*ctrlField[col];
            }
        }

        for (int i = 0; i < nPolyCoeffs; ++i)
        {
            const label row = (nControlPoints + i);

            for (label col = 0; col < nControlPoints; col++)
            {
                beta[i] += mat(row, col)*ctrlField[col];
            }
        }
    }
    else
    {
        // Determine interpolation coefficients (alpha, beta) as the product
        // of the symmetric matrix with the zero-padded control field.
        // Single sweep over the packed upper triangle: each stored
        // off-diagonal value contributes to both (row,col) and (col,row).
        // The polynomial columns (if any) trail each row and only affect beta.

        // Rows >= nControlPoints only see zero-padded entries
        // (the control field is zero there), so they can be skipped.
        for (label row = 0; row < nControlPoints; ++row)
        {
            const Type& xrow = ctrlField[row];

            // Trailing portion of the row: columns (row..n-1)
            // [diagonal, off-diagonals, polynomial columns]
            const scalar* __restrict__ vptr = mat.rowData(row);

            // Diagonal
            Type sum = (*vptr)*xrow;
            ++vptr;

            // Off-diagonals (control point columns)
            for (label col = row + 1; col < nControlPoints; ++col, ++vptr)
            {
                const scalar wt = *vptr;

                sum += wt*ctrlField[col];
                alpha[col] += wt*xrow;
            }

            alpha[row] += sum;

            // Polynomial columns
            for (int i = 0; i < nPolyCoeffs; ++i, ++vptr)
            {
                beta[i] += (*vptr)*xrow;
            }
        }
    }


    // Work array for weights
    scalarField weights(nControlPoints);

    const auto& RBFfunc = RBFptr_();

    forAll(interpPoints_, pointi)
    {
        const auto& refPoint = interpPoints_[pointi];

        // Cut-off function to justify neglecting outer boundary points

        const scalar t =
        (
            (refPoint.dist(focalPoint_) - innerDist_)
          / (outerDist_ - innerDist_)
        );

        if (t < 1)
        {
            // Full calculation of weights
            RBFfunc.weights(weights, controlPoints_, refPoint);

            forAll(controlPoints_, i)
            {
                result[pointi] += weights[i]*alpha[i];
            }

            if (nPolyCoeffs == 4)
            {
                result[pointi] +=
                (
                    beta[0]
                  + beta[1]*refPoint.x()
                  + beta[2]*refPoint.y()
                  + beta[3]*refPoint.z()
                );
            }

            scalar w = 1.0;

            if (t > VSMALL)
            {
                // Cubic smooth step
                w = 1 - Foam::sqr(t)*(3-2*t);
            }

            result[pointi] = w*result[pointi];
        }
    }

    return tresult;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RBFInterpolation::interpolate
(
    const Field<Type>& ctrlField
) const
{
    if (ctrlField.size() != controlPoints_.size())
    {
        FatalErrorInFunction
            << "Incorrect size of source field."
            << " Size = " << ctrlField.size()
            << " nControlPoints = " << controlPoints_.size() << endl
            << abort(FatalError);
    }

    return this->interpolate_impl<Type>(ctrlField);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RBFInterpolation::interpolate
(
    const UIndirectList<Type>& ctrlField
) const
{
    if (ctrlField.size() != controlPoints_.size())
    {
        FatalErrorInFunction
            << "Incorrect size of source field."
            << " Size = " << ctrlField.size()
            << " nControlPoints = " << controlPoints_.size() << endl
            << abort(FatalError);
    }

    return this->interpolate_impl<Type>(ctrlField);
}


// ************************************************************************* //
