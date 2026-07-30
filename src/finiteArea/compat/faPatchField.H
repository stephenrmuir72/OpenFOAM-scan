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

    In finiteArea \c faPatchField<Type> is an alias for \c fvPatchField<Type>
    (declared in \c areaFieldsFwd.H). This header pulls in the full
    \c fvPatchField definition plus that alias so that legacy
    \c \#include "faPatchField.H" keeps resolving (a permanent, retained compatibility forwarder).

\*---------------------------------------------------------------------------*/

#ifndef Foam_compat_faPatchField_H
#define Foam_compat_faPatchField_H

#include "fvPatchField.H"
#include "areaFieldsFwd.H"   // template<class> using faPatchField = fvPatchField

#endif

// ************************************************************************* //
