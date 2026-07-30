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

\*---------------------------------------------------------------------------*/

#include "faMeshProcAddressing.H"
#include "processorPolyPatch.H"
#include "polyBoundaryMesh.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::faMeshProcAddressing::buildMaps()
{
    serialCellOfParentFace_.clear();
    serialPointOfParentPoint_.clear();
    serialFaceOfEdge_.clear();

    serialCellOfParentFace_.resize(2*serialFaceLabels_.size());
    serialPointOfParentPoint_.resize(2*serialPointLabels_.size());
    serialFaceOfEdge_.resize(2*serialFaces_.size());

    // Serial area-cell -> parent face (global, since serial)
    forAll(serialFaceLabels_, celli)
    {
        serialCellOfParentFace_.insert(serialFaceLabels_[celli], celli);
    }

    // Serial FA point -> parent point (global, since serial)
    forAll(serialPointLabels_, pointi)
    {
        serialPointOfParentPoint_.insert(serialPointLabels_[pointi], pointi);
    }

    // Serial FA edge (point pair) -> serial FA face
    forAll(serialFaces_, facei)
    {
        const face& f = serialFaces_[facei];   // 2-point reduced-order face
        serialFaceOfEdge_.insert(edge(f[0], f[1]), facei);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMeshProcAddressing::faMeshProcAddressing(const faMesh& serialFa)
:
    serialFaceLabels_(serialFa.faceLabels()),
    serialPointLabels_(serialFa.pointLabels()),
    serialFaces_(serialFa.faces()),
    serialPatchNames_(serialFa.boundaryMesh().names())
{
    buildMaps();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faMeshProcAddressing::broadcast(const label comm)
{
    Pstream::broadcasts
    (
        comm,
        serialFaceLabels_,
        serialPointLabels_,
        serialFaces_,
        serialPatchNames_
    );

    buildMaps();
}

void Foam::faMeshProcAddressing::calc
(
    const faMesh& procFa,
    const labelUList& parentFaceProcAddr,
    const labelUList& parentPointProcAddr,
    labelList& cellProcAddr,
    labelList& faceProcAddr,
    labelList& pointProcAddr,
    labelList& boundaryProcAddr
) const
{
    // cellProcAddressing: area-face = cell
    const labelList& faceLabels = procFa.faceLabels();
    cellProcAddr.resize_nocopy(faceLabels.size());
    forAll(faceLabels, celli)
    {
        const label globalParentFace = mag(parentFaceProcAddr[faceLabels[celli]]) - 1;
        cellProcAddr[celli] = serialCellOfParentFace_[globalParentFace];
    }

    // pointProcAddressing: FA point -> parent point -> serial FA point
    const labelList& pointLabels = procFa.pointLabels();
    pointProcAddr.resize_nocopy(pointLabels.size());
    forAll(pointLabels, pointi)
    {
        const label globalParentPoint = parentPointProcAddr[pointLabels[pointi]];
        pointProcAddr[pointi] = serialPointOfParentPoint_[globalParentPoint];
    }

    // faceProcAddressing: edge identified by its global FA point pair, signed by
    // orientation relative to the serial FA face.
    const faceList& faces = procFa.faces();
    faceProcAddr.resize_nocopy(faces.size());
    forAll(faces, facei)
    {
        const face& lf = faces[facei];
        const label sp0 = pointProcAddr[lf[0]];
        const label sp1 = pointProcAddr[lf[1]];

        const label serialFacei = serialFaceOfEdge_[edge(sp0, sp1)];
        const face& sfc = serialFaces_[serialFacei];

        faceProcAddr[facei] =
        (
            (sfc[0] == sp0 && sfc[1] == sp1)
          ? (serialFacei + 1)      // same orientation
          : -(serialFacei + 1)     // flipped
        );
    }

    // boundaryProcAddressing: map processor patch -> serial patch by name
    // (patch ordering/count differ between the per-proc and serial meshes, e.g.
    // appended processor patches, so index identity is not valid). Processor
    // patches (absent in serial) map to -1.
    const polyBoundaryMesh& pbm = procFa.boundaryMesh();
    boundaryProcAddr.resize_nocopy(pbm.size());
    forAll(pbm, patchi)
    {
        boundaryProcAddr[patchi] =
        (
            isA<processorPolyPatch>(pbm[patchi])
          ? -1
          : serialPatchNames_.find(pbm[patchi].name())
        );
    }
}


// ************************************************************************* //
