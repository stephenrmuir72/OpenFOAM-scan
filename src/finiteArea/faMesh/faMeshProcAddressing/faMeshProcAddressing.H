/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | www.openfoam.com
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Class
    Foam::faMeshProcAddressing

Description
    Thin derivation of a finiteArea region's processor addressing from the
    parent (finiteVolume) decomposition, for the reduced-order faMesh model
    (area-face = cell, edge = face).

    For a coupled faMesh the FA decomposition is slaved to the parent: an
    area-face is a parent boundary face and an FA point is a parent point. Given
    a serial reference faMesh (built from the undecomposed parent) this class
    maps a per-processor faMesh onto it, producing standard fvMesh-style
    \c {cell,face,point,boundary}ProcAddressing:

      - \b cellProcAddressing : coupling faceLabels mapped through the parent
        \c faceProcAddressing, ranked into the serial (sorted-parent-face) order.
      - \b pointProcAddressing : FA point -> parent point (faMesh::pointLabels)
        mapped through the parent \c pointProcAddressing.
      - \b faceProcAddressing : each edge identified by its global FA point pair
        and looked up in the serial FA edge->face map; orientation encoded as the
        signed \c +/-(serialFace+1) convention. Processor-boundary edges resolve
        through the shared points with no extra communication.
      - \b boundaryProcAddressing : identity, \c -1 for processor patches.

    This replaces the collective addressing core of the legacy
    faMeshReconstructor: no serial-mesh rebuild and no cross-processor half-edge
    matching are required.

SourceFiles
    faMeshProcAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faMeshProcAddressing_H
#define Foam_faMeshProcAddressing_H

#include "faMesh.H"
#include "HashTable.H"
#include "EdgeMap.H"
#include "UPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class faMeshProcAddressing Declaration
\*---------------------------------------------------------------------------*/

class faMeshProcAddressing
{
    // Private Data

        //- Serial area-cell -> parent face (== serialFa.faceLabels())
        labelList serialFaceLabels_;

        //- Serial FA point -> parent point (== serialFa.pointLabels())
        labelList serialPointLabels_;

        //- Serial FA faces (reduced-order, 2-point) in serial-point ids
        faceList serialFaces_;

        //- Serial FA boundary patch names
        wordList serialPatchNames_;

        //- Global parent face -> serial area-cell (derived)
        Map<label> serialCellOfParentFace_;

        //- Global parent point -> serial FA point (derived)
        Map<label> serialPointOfParentPoint_;

        //- Serial FA point pair -> serial FA face (edge) (derived)
        EdgeMap<label> serialFaceOfEdge_;


    // Private Member Functions

        //- (Re)build the lookup maps from the serial reference lists. Called by
        //- the constructor and after broadcast().
        void buildMaps();


public:

    // Constructors

        //- Construct empty (for subsequent broadcast() from master)
        faMeshProcAddressing() = default;

        //- Construct from a serial reference area mesh. The serial topology is
        //- captured up front, so the source mesh need not outlive this object.
        explicit faMeshProcAddressing(const faMesh& serialFa);


    // Member Functions

        //- Broadcast the serial reference topology from the master rank to all
        //- ranks and rebuild the lookup maps. Enables SPMD use (e.g.
        //- redistributePar) where only the master holds the complete faMesh.
        void broadcast(const label comm = UPstream::worldComm);

        //- Derive the FA region processor addressing for a per-processor
        //- faMesh from the parent decomposition addressing.
        void calc
        (
            const faMesh& procFa,
            const labelUList& parentFaceProcAddr,
            const labelUList& parentPointProcAddr,
            labelList& cellProcAddr,
            labelList& faceProcAddr,
            labelList& pointProcAddr,
            labelList& boundaryProcAddr
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
