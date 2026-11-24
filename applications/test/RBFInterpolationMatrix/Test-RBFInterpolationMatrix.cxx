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
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Application
    Test-RBFInterpolationMatrix

Description
    Basic checks of RBFInterpolation::InterpolationMatrix
    (packed upper-triangle symmetric matrix).

\*---------------------------------------------------------------------------*/

#include "RBFInterpolation.H"
#include "SquareMatrix.H"
#include "IOmanip.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    typedef RBFInterpolation::InterpolationMatrix matrixType;

    matrixType empty;
    Info<< "default: n=" << empty.n() << " size=" << empty.size()
        << " empty=" << empty.empty() << nl;

    const label n = 4;
    SquareMatrix<scalar> mat(n, Zero);
    scalar val = 0;
    for (label i = 0; i < n; ++i)
    {
        for (label j = 0; j < n; ++j)
        {
            mat(i, j) = val++;  // deliberately non-symmetric
        }
    }

    const matrixType sym(mat);
    Info<< "n=" << sym.n() << " stored=" << sym.size()
        << " (expect " << (n*(n+1))/2 << ")" << nl;

    Info<< nl << "Source matrix (non-symmetric):" << nl << mat << nl;

    Info<< nl << "Packed upper-triangle storage:" << nl << "(";
    for (label k = 0; k < sym.size(); ++k)
    {
        Info<< (k ? " " : "") << sym.cdata()[k];
    }
    Info<< ")" << nl;

    Info<< nl << "Expanded via operator(i,j):" << nl;
    for (label i = 0; i < n; ++i)
    {
        for (label j = 0; j < n; ++j)
        {
            Info<< ' ' << Foam::setw(4) << sym(i, j);
        }
        Info<< nl;
    }

    {
        matrixType copy(mat);
        copy(3, 0) = -1;  // writes through to the (0,3) storage
        Info<< nl << "After copy(3,0) = -1: copy(0,3) = " << copy(0, 3)
            << nl;
    }

    Info<< nl;

    label nFail = 0;
    for (label i = 0; i < n; ++i)
    {
        for (label j = 0; j < n; ++j)
        {
            const scalar expected = (i <= j) ? mat(i, j) : mat(j, i);
            if (sym(i, j) != expected)
            {
                ++nFail;
                Info<< "FAIL (" << i << ',' << j << "): "
                    << sym(i, j) << " != " << expected << nl;
            }
        }
    }

    if (sym.size() != (n*(n+1))/2 || sym.empty() || !sym.cdata())
    {
        ++nFail;
    }

    Info<< (nFail ? "FAILED" : "PASSED") << nl;
    return nFail ? 1 : 0;
}

// ************************************************************************* //
