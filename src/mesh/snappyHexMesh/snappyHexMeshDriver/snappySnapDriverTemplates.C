/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "snappySnapDriver.H"
#include "PatchTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceList>
Foam::labelList Foam::snappySnapDriver::getFacePoints
(
    const indirectPrimitivePatch& pp,
    const FaceList& faces
)
{
    // Could use PrimitivePatch & localFaces to extract points but might just
    // as well do it ourselves.

    boolList pointOnZone(pp.nPoints(), false);

    forAll(faces, i)
    {
        const face& f = faces[i];

        forAll(f, fp)
        {
            label meshPointi = f[fp];

            const auto iter = pp.meshPointMap().cfind(meshPointi);

            if (iter.good())
            {
                const label pointi = iter.val();
                pointOnZone[pointi] = true;
            }
        }
    }

    return findIndices(pointOnZone, true);
}


// Gather and merge patch data onto master processor, returning the
// merged patch and maps:
// - globalFaces : globalIndex from merged patch faces (on master proc only)
//   to the input patch. Send information from merged to local with
//   globalFaces.scatter
// - points : since points can be duplicated is better done through a
//            mapDistribute:
// const pointField testFld(mergedPatch.localPoints());
// pointMapPtr().distribute(testFld);
// forAll(ppToMerged, patchPointi)
// {
//    Pout<< "local:" << pp.localPoints()[patchPointi]
//        << " merged:" << testFld[ppToMerged[patchPointi]]
// }
template<class PatchType>
Foam::autoPtr<Foam::PrimitivePatch<Foam::faceList, Foam::pointField>>
Foam::snappySnapDriver::gatherAndMerge
(
    const polyMesh& mesh,
    const PatchType& pp,

    autoPtr<globalIndex>& globalFacesPtr,
    autoPtr<mapDistribute>& pointMapPtr,
    labelList& ppToMerged
)
{
    // Collect merge information
    labelList pointToGlobal;
    labelList uniqueMeshPointLabels;
    autoPtr<globalIndex> globalPointsPtr;
    List<face> mergedFaces;
    pointField mergedPoints;
    PatchTools::gatherAndMerge
    (
        mesh,
        pp.localFaces(),
        pp.meshPoints(),
        pp.meshPointMap(),

        pointToGlobal,
        uniqueMeshPointLabels,
        globalPointsPtr,
        globalFacesPtr,
        mergedFaces,
        mergedPoints
    );

    // Create merged patch
    autoPtr<PrimitivePatch<faceList, pointField>> mergedPatchPtr
    (
        new PrimitivePatch<faceList, pointField>
        (
            std::move(mergedFaces),
            std::move(mergedPoints)
        )
    );
    const auto& mergedPatch = mergedPatchPtr();


    // Addressing
    // ~~~~~~~~~~

    // Scatter the merged point indices (on master) back to the
    // senders in global mergedPoints indices. Note that since the
    // mergedPatch has faces only on proc0 we could avoid the 'toGlobal'...
    const globalIndex globalMerged(mergedPatch.nPoints());
    const labelList masterPts
    (
        globalPointsPtr->scatter
        (
            globalMerged.toGlobal(0, identity(mergedPatch.nPoints()))
        )
    );

    // Re-work on the local processor back to my pp.localPoint
    ppToMerged.resize_nocopy(pp.nPoints());
    ppToMerged = -1;
    forAll(uniqueMeshPointLabels, i)
    {
        const label meshPointi = uniqueMeshPointLabels[i];
        const auto fnd = pp.meshPointMap().cfind(meshPointi);
        if (fnd.good())
        {
            const label patchPointi = fnd.val();
            ppToMerged[patchPointi] = masterPts[i];
        }
    }

    // Now ppToMerged is still -1 for slave coupled points. Use the
    // coupled points to sync for max
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        ppToMerged,
        maxEqOp<label>(),
        -1
    );

    // Rework the global indices into a communication pattern
    List<Map<label>> compactMap;
    pointMapPtr.reset
    (
        new mapDistribute
        (
            globalMerged,
            ppToMerged,
            compactMap
        )
    );
    if (UPstream::master())
    {
        // Fix map. See issue 3509. TBD.
        pointMapPtr().subMap()[0] = identity(pp.nPoints());
        pointMapPtr().constructMap()[0] = ppToMerged;
    }

    return mergedPatchPtr;
}


// ************************************************************************* //
