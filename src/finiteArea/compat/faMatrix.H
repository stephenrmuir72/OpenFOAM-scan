/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

Description
    finiteArea backward-compatibility forwarding header.

    Legacy finiteArea shipped a dedicated \c faMatrix.H / \c faScalarMatrix.H;
    in finiteArea \c faMatrix<Type> is an alias for \c fvMatrix<Type> and all
    the matrix typedefs are collected in \c faMatrices.H. This header simply
    forwards to \c faMatrices.H so that legacy \c \#include "faMatrix.H" keeps
    resolving (a permanent, retained compatibility forwarder).

\*---------------------------------------------------------------------------*/

#ifndef Foam_compat_faMatrix_H
#define Foam_compat_faMatrix_H

#include "faMatrices.H"

#endif

// ************************************************************************* //
