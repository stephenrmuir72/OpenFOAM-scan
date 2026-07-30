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

    In finiteArea area patch fields ARE finite-volume patch fields
    (\c faPatchField<Type> = \c fvPatchField<Type>), so the concrete boundary
    condition set is provided by \c fvPatchFields.H and selected at run time by
    the same type names. This header pulls in the finite-volume patch-field set
    plus the \c fa* aliases so that legacy \c \#include "faPatchFields.H" keeps
    resolving (a permanent, retained compatibility forwarder).

    Note: legacy per-type C++ class names such as \c zeroGradientFaPatchField
    are NOT provided; C++ code naming a concrete faPatchField subclass must use
    the \c fv* equivalent. Dictionary-driven run-time BC selection is unchanged.

\*---------------------------------------------------------------------------*/

#ifndef Foam_compat_faPatchFields_H
#define Foam_compat_faPatchFields_H

#include "fvPatchFields.H"
#include "areaFieldsFwd.H"

#endif

// ************************************************************************* //
