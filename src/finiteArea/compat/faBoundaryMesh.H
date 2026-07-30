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

    In finiteArea \c faMesh is a thin \c fvMesh subclass, so its boundary is a
    \c fvBoundaryMesh (\c faMesh::boundary() returns \c const fvBoundaryMesh&).
    This header provides the legacy \c faBoundaryMesh name as an alias for
    \c fvBoundaryMesh so that legacy \c \#include "faBoundaryMesh.H" and uses of
    the \c faBoundaryMesh type keep resolving (now and after the finiteArea ->
    finiteArea rename).

    Note: the legacy \c faBoundaryMesh was an edge-addressed \c faPatch
    container; the finiteArea alias is the reduced-order \c fvBoundaryMesh
    (faces == edges). The public accessor surface (size(), operator[],
    findPatchID(), names(), types(), ...) matches.

\*---------------------------------------------------------------------------*/

#ifndef Foam_compat_faBoundaryMesh_H
#define Foam_compat_faBoundaryMesh_H

#include "fvBoundaryMesh.H"

namespace Foam
{
    //- Backward-compatibility alias (finiteArea): the finite-area boundary
    //  mesh is the reduced-order finite-volume boundary mesh.
    typedef fvBoundaryMesh faBoundaryMesh;
}

#endif

// ************************************************************************* //
