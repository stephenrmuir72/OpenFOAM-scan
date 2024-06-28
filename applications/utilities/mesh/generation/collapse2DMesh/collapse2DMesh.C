/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019, 2024 OpenCFD Ltd.
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

Application
    collapse2DMesh

Group
    grpMeshGenerationUtilities

Description
    Takes 3D mesh with empty front and back and generate 2D version.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "topoSet.H"
#include "processorMeshes.H"
#include "emptyPolyPatch.H"
#include "removePoints.H"
#include "bitSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Converts 2D polyMesh into polyMesh2D"
    );

    argList::addArgument("patch", "The (empty) patch giving the geometry");

    #include "addOverwriteOption.H"
    #include "addRegionOption.H"

    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    #include "createTime.H"

    #include "createNamedMesh.H"

    const word oldInstance = mesh.pointsInstance();

    const wordRe patchName(args[1]);
    const bool overwrite = args.found("overwrite");

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    const label keepPatchi = pbm.findPatchID(patchName);
    if (keepPatchi == -1 || !isA<emptyPolyPatch>(pbm[keepPatchi]))
    {
        FatalErrorInFunction
            << "Patch " << patchName << " does not exist or is not empty"
            << exit(FatalError);
    }


    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    {
        // Constructing a 2D mesh
        polyTopoChange meshMod(mesh, true, 2);

        const polyPatch& keepPatch = pbm[keepPatchi];

        // Get local points on the patch. Note: no need to parallel sync
        const bitSet isPatchPoint(mesh.nPoints(), keepPatch.meshPoints());

        labelList cellToPatchFace(mesh.nCells(), -1);
        const labelUList& fc = keepPatch.faceCells();
        forAll(fc, patchFacei)
        {
            cellToPatchFace[fc[patchFacei]] = patchFacei;
        }


        face newFace(2);
        DynamicList<label> newVerts(2);

        for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
        {
            const face& oldFace = mesh.faces()[facei];

            newVerts.clear();
            for (const auto& pointi : oldFace)
            {
                if (isPatchPoint[pointi])
                {
                    newVerts.append(pointi);
                }
            }

            if (newVerts.size() != 2)
            {
                FatalErrorInFunction << "Problem: face "
                    << facei << " at " << mesh.faceCentres()[facei]
                    << " has more than 2 vertices on the patch"
                    << exit(FatalError);
            }
            forAll(newFace, i)
            {
                newFace[i] = newVerts[i];
            }

            // Check that patch face ordering is opposite of owner side
            const face& patchFace = keepPatch[cellToPatchFace[own[facei]]];
            const label index = patchFace.find(newFace[0]);
            if (patchFace[patchFace.fcIndex(index)] != newFace[1])
            {
                Swap(newFace[0], newFace[1]);
            }

            meshMod.modifyFace
            (
                newFace,
                facei,
                own[facei],
                nei[facei],
                false,
                -1,
                -1,
                false
            );
        }
        forAll(pbm, patchi)
        {
            //if (patchi != keepPatchi)
            const polyPatch& pp = pbm[patchi];
            if (!isA<emptyPolyPatch>(pp))   // || patchi == keepPatchi)
            {
                forAll(pp, patchFacei)
                {
                    const label meshFacei = pp.start()+patchFacei;
                    const face& oldFace = mesh.faces()[meshFacei];

                    newVerts.clear();
                    for (const auto& pointi : oldFace)
                    {
                        if (isPatchPoint[pointi])
                        {
                            newVerts.append(pointi);
                        }
                    }

                    if (newVerts.size() != 2)
                    {
                        FatalErrorInFunction << "Problem: face "
                            << meshFacei << " on patch " << patchi
                            << " at " << mesh.faceCentres()[meshFacei]
                            << " has more than 2 vertices on the patch"
                            << exit(FatalError);
                    }
                    forAll(newFace, i)
                    {
                        newFace[i] = newVerts[i];
                    }

                    // Check that patch face ordering is opposite of owner side
                    const label keepFacei = cellToPatchFace[own[meshFacei]];
                    const face& patchFace = keepPatch[keepFacei];
                    const label index = patchFace.find(newFace[0]);
                    if (patchFace[patchFace.fcIndex(index)] != newFace[1])
                    {
                        Swap(newFace[0], newFace[1]);
                    }

                    meshMod.modifyFace
                    (
                        newFace,
                        meshFacei,
                        own[meshFacei],
                        -1,
                        false,
                        patchi,
                        -1,
                        false
                    );
                }
            }
            else
            {
                forAll(pp, patchFacei)
                {
                    meshMod.removeFace(pp.start()+patchFacei, -1);
                }
            }
        }

        // Create a mesh from topo changes.
        autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

        mesh.updateMesh(morphMap);
    }

    // Remove excess points
    {
        removePoints pointRemover(mesh, false);
        boolList pointCanBeDeleted;
        pointRemover.countPointUsage(-GREAT, pointCanBeDeleted);

        polyTopoChange meshMod(mesh, true, 2);
        pointRemover.setRefinement(pointCanBeDeleted, meshMod);

        autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

        mesh.updateMesh(morphMap);
        pointRemover.updateMesh(morphMap);
    }


    if (!overwrite)
    {
        runTime++;
    }
    else
    {
        mesh.setInstance(oldInstance);
    }

    Info<< "\nWriting collapsed mesh to time = " << runTime.timeName()
        << nl << endl;

    mesh.write();
    //topoSet::removeFiles(mesh);   // only (boundary)faceSets are wrong
    processorMeshes::removeFiles(mesh);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
