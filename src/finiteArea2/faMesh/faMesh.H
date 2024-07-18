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
    Foam::faMesh

Description
    Finite area mesh (used for 2-D non-Euclidian finite area method)
    defined using a \em patch of faces on a polyMesh
    (ie, uindirectPrimitivePatch).

    The ordering of faces and points on the faMesh corresponds to
    the localFaces and localPoints as per Foam::PrimitivePatch but
    the edge addressing is handled slightly differently.
    The internal edges of the faMesh will generally correspond identically
    to the internalEdges of the PrimitivePatch (may change in the future)
    but the boundaryEdges will be reordered compared to the PrimitivePatch
    to allow edge boundary slices to be obtained.

SourceFiles
    faMesh.C
    faMeshDemandDrivenData.C
    faMeshPatches.C
    faMeshTopology.C
    faMeshUpdate.C

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#ifndef Foam_faMesh_H
#define Foam_faMesh_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef fvMesh faMesh;
//class faMesh
//:
//public fvMesh
//{
//    labelPairList whichPatchFaces_;
//
//public:
//
//    //- Runtime type information
//    TypeName("faMesh");
//
//
//    // Constructors
//
//        //- Construct from IOobject
//        explicit faMesh
//        (
//            const IOobject& io,
//            const labelUList& meshFaces,
//            const bool doInit=true
//        )
//        :
//            fvMesh(io, doInit)
//        {
//            const refPtr<fvMesh> tmesh = fvMesh::parentMesh();
//            if (tmesh)
//            {
//                const auto& mesh = tmesh();
//                whichPatchFaces_.resize_nocopy(meshFaces.size());
//                forAll(meshFaces, i)
//                {
//                    whichPatchFaces_[i] =
//                        mesh.boundaryMesh().whichPatchFace(meshFaces[i]);
//                }
//            }
//        }
//
//
//    // Access
//
//            //- The polyPatch/local-face for each faceLabels()
//            inline const List<labelPair>& whichPatchFaces() const
//            {
//                return whichPatchFaces_;
//            }
//}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
