/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2009 Frank Bos, TU Delft
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
#include "dictionary.H"
#include "scalarMatrices.H"
#include "SubField.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Invert the (symmetric) matrix using LU decomposition with pivoting,
// directly into compact storage.
// Cholesky is not applicable: the matrix is indefinite with polynomials.
//
// On input, only the upper triangle of 'mat' needs to be populated.
// The matrix is used as the LU workspace, its contents are destroyed.
static RBFInterpolation::InterpolationMatrix
LUinvert
(
    SquareMatrix<scalar>&& A
)
{
    // Modify the original input matrix for the substitutions
    SquareMatrix<scalar> mat(std::move(A));

    const label n = mat.n();

    // Mirror upper triangle for the pivoting LU
    for (label i = 0; i < n; ++i)
    {
        for (label j = i + 1; j < n; ++j)
        {
            mat(j, i) = mat(i, j);
        }
    }

    scalarList column(n);
    labelList pivotIndices(n);

    Foam::LUDecompose(mat, pivotIndices);

    // The inverse is symmetric: only the upper triangle is stored
    RBFInterpolation::InterpolationMatrix result(n);

    for (label j = 0; j < n; ++j)
    {
        column = Foam::zero{};
        column[j] = scalar(1);

        Foam::LUBacksubstitute(mat, pivotIndices, column);

        for (label i = 0; i <= j; ++i)
        {
            result(i, j) = column[i];
        }
    }

    return result;
}

}  // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::RBFInterpolation::calcInterpolationMatrix() const
{
    if (matrixPtr_)
    {
        FatalErrorInFunction
            << "The interpolation matrix has already been calculated"
            << exit(FatalError);
    }

    // Determine inverse of boundary connectivity matrix

    const label nControlPoints = controlPoints_.size();

    // Fill Nb x Nb matrix
    SquareMatrix<scalar> A
    (
        nControlPoints + (withPolynomials_ ? 4 : 0),
        Foam::zero{}
    );

    const auto& RBFfunc = RBFptr_();

    // Build the RBF coefficient matrix, A(i,j) = phi(|x_i - x_j|),
    // where phi is the selected radial basis function.
    // A is symmetric, so only the upper triangle (j >= i) is filled.
    // Row r therefore only needs the weights between control point r
    // and control points r..N-1 (the diagonal included).

    for (label row = 0; row < nControlPoints; ++row)
    {
        const auto& refPoint = controlPoints_[row];
        const label subLen = (nControlPoints - row);

        // The upper-triangle points
        const auto subPoints = controlPoints_.slice(row, subLen);

        // The weights are the corresponding upper-triangle matrix coeffs
        UList<scalar> subWeights(A.rowData(row) + row, subLen);

        // Fill in the [0-1] weights
        RBFfunc.weights(subWeights, subPoints, refPoint);
    }

    if (withPolynomials_)
    {
        // Polynomial columns: 1, x, y, z
        for (label col = 0; col < nControlPoints; ++col)
        {
            const point& pt = controlPoints_[col];

            A(col, nControlPoints)     = 1.0;
            A(col, nControlPoints + 1) = pt.x();
            A(col, nControlPoints + 2) = pt.y();
            A(col, nControlPoints + 3) = pt.z();
        }

        // This is not needed (already zeroed in initializer)
        //
        // 4x4 zero block (upper triangle)
        // for (label row = nControlPoints; row < (nControlPoints + 4); ++row)
        // {
        //     for (label col = row; col < (nControlPoints + 4); ++col)
        //     {
        //         A(row, col) = 0.0;
        //     }
        // }
    }

    // Debugging
    if constexpr (false)
    {
        // For debugging, also need to mirror the upper-triangle
        const label n = A.n();

        // Mirror upper triangle
        for (label i = 0; i < n; ++i)
        {
            for (label j = i + 1; j < n; ++j)
            {
                A(j, i) = A(i, j);
            }
        }

        Info<< "Control Points: "; controlPoints_.writeList(Info) << nl;
        Info<< "A[" << A.m() << "][" << A.n() << "]" << nl;
        A.writeMatrix(Info);
    }

    Info<< "Inverting RBF motion matrix: size=" << A.m() << endl;

    // Invert (A is consumed) directly into compact storage
    matrixPtr_ = std::make_unique<InterpolationMatrix>
    (
        LUinvert(std::move(A))
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBFInterpolation::RBFInterpolation
(
    const dictionary& dict,
    const pointField& controlPoints,
    const pointField& interpPoints
)
:
    controlPoints_(controlPoints),
    interpPoints_(interpPoints),
    RBFptr_
    (
        RBFFunction::New(dict.get<word>("RBF"), dict)
    ),
    matrixPtr_(nullptr),
    focalPoint_(dict.get<point>("focalPoint")),
    innerDist_(dict.get<scalar>("innerDistance")),
    outerDist_(dict.get<scalar>("outerDistance")),
    withPolynomials_(dict.getOrDefault<bool>("polynomials", false))
{}


Foam::RBFInterpolation::RBFInterpolation
(
    const RBFInterpolation& interp
)
:
    controlPoints_(interp.controlPoints_),
    interpPoints_(interp.interpPoints_),
    RBFptr_(interp.RBFptr_.clone()),
    matrixPtr_(nullptr),
    focalPoint_(interp.focalPoint_),
    innerDist_(interp.innerDist_),
    outerDist_(interp.outerDist_),
    withPolynomials_(interp.withPolynomials_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// Non-default in header (incomplete types)
Foam::RBFInterpolation::~RBFInterpolation() = default;


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RBFInterpolation::clearOut()
{
    matrixPtr_.reset(nullptr);
}


// ************************************************************************* //
