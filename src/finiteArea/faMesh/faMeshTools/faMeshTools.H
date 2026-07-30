/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 OpenCFD Ltd.
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

Class
    Foam::faMeshTools

Description
    A collection of tools for operating on a finiteArea faMesh (a
    reduced-order 2D fvMesh: FA "faces" are fvMesh cells, FA "edges" are
    fvMesh faces).

    Currently exposes the mesh-check reporting used by the makeFaMesh
    utility. The legacy finiteArea decompose/reconstruct helpers are not
    ported: decomposition of finiteArea meshes is handled by
    decomposePar / redistributePar directly.

SourceFiles
    faMeshToolsChecks.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faMeshTools_H
#define Foam_faMeshTools_H

#include "faMesh.H"
#include "edgeFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class faMeshTools Declaration
\*---------------------------------------------------------------------------*/

class faMeshTools
{
public:

    //- Report mesh information (counts, patches, geometry ranges).
    //  Human-facing labels match the legacy finiteArea output; the
    //  underlying quantities are read from the reduced-order fvMesh.
    static void printMeshChecks
    (
        const faMesh& mesh,
        const int verbose = 1
    );

    //- Flatten an edge field (internal and boundary) into a single linear
    //  field addressed by reduced-order face (edge) index.
    //  In finiteArea an edge field IS a surfaceField, so this concatenates
    //  the internal-face values with the boundary-face values in mesh-face
    //  ordering (length == mesh.nFaces()).
    template<class Type>
    static void flattenEdgeField
    (
        UList<Type>& result,
        const GeometricField<Type, faePatchField, edgeMesh>& fld
    );

    //- Flatten an edge field (internal and boundary) into linear addressing
    template<class Type>
    static tmp<Field<Type>> flattenEdgeField
    (
        const GeometricField<Type, faePatchField, edgeMesh>& fld
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "faMeshTools.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
