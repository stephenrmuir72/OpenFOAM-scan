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

Class
    Foam::fvMesh2DTools

Description
    A collection of tools for operating on a 2D fvMesh.

SourceFiles
    fvMesh2DTools.C
    fvMesh2DToolsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvMesh2DTools_H
#define Foam_fvMesh2DTools_H

#include "fvMesh.H"
#include "EdgeMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyTopoChange;

/*---------------------------------------------------------------------------*\
                        Class fvMesh2DTools Declaration
\*---------------------------------------------------------------------------*/

class fvMesh2DTools
{
protected:

    //- Insert changes to extract the faces as 2D cells
    template<class PatchType>
    static void patchToMesh
    (
        const polyMesh& mesh,
        const PatchType& pp,
        const labelList& addressing,
        const EdgeMap<label>& edgeToPatch,
        polyTopoChange& meshMod,
        labelList& pointMap     // From point to point on parent mesh
                                // (cells are in order of the faceZone)
    );


public:

    //- Calculate per fz edge the patch. Currently does not handle side-ness of
    //- faceZone.
    template<class PatchType>
    static EdgeMap<label> edgePatch
    (
        const polyMesh& mesh,
        const label defaultPatch,
        const PatchType& pp,
        const labelList& addressing
    );

    //- Create 2D mesh from faceZone
    template<class PatchType>
    static autoPtr<fvMesh> createMesh
    (
        const IOobject& io,     // name, registry
        const fvMesh& mesh,     // parent mesh
        const PatchType& pp,
        const labelList& addressing,
        const EdgeMap<label>& edgeToPatch,
        labelList& pointMap     // From point to point on parent mesh
                                // (cells are in order of the faceZone)
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMesh2DToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
