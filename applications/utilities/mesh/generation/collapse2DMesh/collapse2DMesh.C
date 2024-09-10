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
    Takes 3D mesh and generates 2D slice

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "topoSet.H"
#include "processorMeshes.H"
#include "emptyPolyPatch.H"
#include "bitSet.H"
#include "fvMesh2DTools.H"
#include "IOmapPolyMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Converts patch of polyMesh into polyMesh2D"
    );

    argList::addArgument("patch", "The patch or faceZone giving the geometry");

    //#include "addOverwriteOption.H"
    #include "addRegionOption.H"

    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const word oldInstance = mesh.pointsInstance();

    const wordRe patchName(args[1]);

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    const label keepPatchi = pbm.findPatchID(patchName);
    if (keepPatchi == -1)   // || !isA<emptyPolyPatch>(pbm[keepPatchi]))
    {
        FatalErrorInFunction
            << "Patch " << patchName << " does not exist. Valid patches are "
            << pbm.names() << exit(FatalError);
    }
    const label defaultPatchID = -1;

    const polyPatch& pp = pbm[keepPatchi];
//    const edgeList& edges = pp.edges();
//    const labelListList& edgeFaces = pp.edgeFaces();
//    const labelList& mp = pp.meshPoints();
//
//    forAll(edgeFaces, edgei)
//    {
//        if (edgeFaces[edgei].size() > 2)
//        {
//            FatalErrorInFunction << "Patch is non-manifold at edge "
//                << edges[edgei].line(pp.localPoints())
//                << exit(FatalError);
//        }
//    }
//
//    // Mark on pp outside the connected patches. Note: maybe also have
//    // orientation so if extruding faceZone we can pick the 'correct' patch
//
//    EdgeMap<label> edgeToPatch(pp.nEdges()-pp.nInternalEdges());
//
//    {
//        // Insert (with -1) all wanted edges
//        for (label edgei = pp.nInternalEdges(); edgei < pp.nEdges(); edgei++)
//        {
//            const edge meshE(mp, edges[edgei]);
//            edgeToPatch.insert(meshE, -1);
//        }
//
//        // Find on other patches
//        forAll(pbm, patchi)
//        {
//            if (patchi != keepPatchi)
//            {
//                for (const edge& patchE : pbm[patchi].edges())
//                {
//                    const edge meshE(pbm[patchi].meshPoints(), patchE);
//                    auto fnd = edgeToPatch.find(meshE);
//                    if (fnd)
//                    {
//                        if (fnd() == -1)
//                        {
//                            //Pout<< "For mesh edge:" << meshE
//                            //    << " points:" << meshE.line(mesh.points())
//                            //    << " found patch " << patchi
//                            //    << " name:" << pbm[patchi].name()
//                            //    << endl;
//                            fnd() = patchi;
//                        }
//                    }
//                }
//            }
//        }
//    }


    const labelList addressing(pp.range().labels());

    const EdgeMap<label> edgeToPatch
    (
        fvMesh2DTools::edgePatch
        (
            mesh,
            defaultPatchID,
            pp,
            addressing
        )
    );



//    // Constructing a 2D mesh
//    polyTopoChange meshMod(pbm.size(), true, 2);
//    // Pre-size with expected number of cells, faces, points
//    meshMod.setCapacity(mp.size(), edges.size(), pp.size());
//
//    labelList faceToCell(pp.size());
//    forAll(pp, facei)
//    {
//        faceToCell[facei] = meshMod.addCell
//        (
//            -1,
//            -1,
//            pp.start()+facei,   // original mesh face label
//            -1,
//            -1
//        );
//    }
//
//    labelList pointToPoint(pp.nPoints());
//    forAll(pp.localPoints(), pointi)
//    {
//        pointToPoint[pointi] = meshMod.addPoint
//        (
//            pp.localPoints()[pointi],
//            mp[pointi],     // original mesh point label
//            -1,
//            true
//        );
//    }
//
//    labelList edgeToFace(edges.size(), -1);
//    face newFace(2);
//    forAll(pp, facei)
//    {
//        const auto& f = pp.localFaces()[facei];
//        const auto& fEdges = pp.faceEdges()[facei];
//
//        forAll(f, fp)
//        {
//            const label edgei = fEdges[fp];
//
//            if (edgeToFace[edgei] != -1)
//            {
//                // Edge already done
//                continue;
//            }
//
//            // Get neighbour
//            const auto& eFaces = edgeFaces[edgei];
//
//            // Get edge in same order as face
//            newFace[0] = pointToPoint[f[fp]];
//            newFace[1] = pointToPoint[f[f.fcIndex(fp)]];
//
//            if (eFaces.size() == 2)
//            {
//                label ownCell = faceToCell[facei];
//                label nbrCell = -1;
//                if (eFaces[0] == facei)
//                {
//                    nbrCell = faceToCell[eFaces[1]];
//                }
//                else
//                {
//                    nbrCell = faceToCell[eFaces[0]];
//                }
//
//                if (nbrCell < ownCell)
//                {
//                    std::swap(newFace[0], newFace[1]);
//                    std::swap(ownCell, nbrCell);
//                }
//
//                edgeToFace[edgei] = meshMod.addFace
//                (
//                    newFace,
//                    ownCell,
//                    nbrCell,
//                    -1,             //masterPointID,
//                    -1,             //masterEdgeID,
//                    -1,             //masterFaceID,
//                    false,          // flipFaceFlux
//                    -1,             // patchID
//                    -1,             // zoneID
//                    false           // zoneFlip
//                );
//            }
//            else
//            {
//                const label ownCell = faceToCell[facei];
//                const edge meshE(mp, edges[edgei]);
//                const label patchi = edgeToPatch[meshE];
//
//                edgeToFace[edgei] = meshMod.addFace
//                (
//                    newFace,
//                    ownCell,
//                    -1,
//                    -1,             //masterPointID,
//                    -1,             //masterEdgeID,
//                    -1,             //masterFaceID,
//                    false,          // flipFaceFlux
//                    patchi,         // patchID
//                    -1,             // zoneID
//                    false           // zoneFlip
//                );
//            }
//        }
//    }
//
//
//    // Create dummy fvSchemes, fvSolution
//    fvMeshTools::createDummyFvMeshFiles
//    (
//        mesh,
//        patchName,
//        true
//    );
//
//    // Create a mesh from topo changes.
//    //autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);
//    //mesh.updateMesh(morphMap);
//    autoPtr<fvMesh> subMeshPtr;
//    autoPtr<mapPolyMesh> morphMap = meshMod.makeMesh
//    (
//        subMeshPtr,
//        IOobject
//        (
//            patchName,
//            mesh.pointsInstance(),
//            //"finite-area",
//            mesh.thisDb(),
//            IOobjectOption::NO_READ
//        ),
//        mesh
//    );

    labelList pointMap;

    autoPtr<fvMesh> subMeshPtr
    (
        fvMesh2DTools::createMesh
        (
            IOobject
            (
                patchName,
                mesh.pointsInstance(),
                //"finite-area",
                mesh.thisDb(),
                IOobjectOption::NO_READ
            ),
            mesh,
            pp,
            addressing,
            edgeToPatch,
    
            pointMap         // From point to parent mesh point
        )
    );

    // TBD:
    // - remove excess processor patches
    // - write as IOmapPolyMesh
    // - check where to write. Check e.g. subsetMesh? decomposePar? Is map
    //   back to original mesh.
    //      - pointMap back to patch-point
    //      - faceMap back to original mesh faces (only if structured since
    //        face actually created from original mesh edge)
    //      - cellMap back to original mesh boundary faces
    // - subset fields? Is a bit like subsetMesh.
    // - extract from faceZone. Have defaultFaces patch if zone is fully inside?

    subMeshPtr().setInstance(oldInstance);

    Info<< "\nWriting collapsed mesh to time = " << runTime.timeName()
        << nl << endl;

    subMeshPtr().write();

    //- Does not work yet. TBD.
    //IOmapPolyMesh ioMorphMap
    //(
    //    IOobject
    //    (
    //        "map",
    //        subMeshPtr().facesInstance(),
    //        fvMesh::meshSubDir,
    //        subMeshPtr().thisDb(),
    //        IOobjectOption::NO_READ
    //    ),
    //    morphMap()
    //);
    //ioMorphMap.write();
    topoSet::removeFiles(subMeshPtr());
    processorMeshes::removeFiles(subMeshPtr());

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
