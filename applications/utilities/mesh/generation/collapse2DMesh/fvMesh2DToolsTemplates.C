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

\*---------------------------------------------------------------------------*/

#include "polyTopoChange.H"
#include "syncTools.H"
#include "dummyTransform.H"
#include "fvMeshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class PatchType>
void Foam::fvMesh2DTools::patchToMesh
(
    const polyMesh& mesh,
    const PatchType& pp,
    const labelList& addressing,

    const EdgeMap<label>& edgeToPatch,
    polyTopoChange& meshMod,
    labelList& pointMap
)
{
    const edgeList& edges = pp.edges();
    const labelListList& edgeFaces = pp.edgeFaces();
    const labelListList& faceEdges = pp.faceEdges();
    const labelList& mp = pp.meshPoints();
    const faceList& localFaces = pp.localFaces();

    // Check
    forAll(edgeFaces, edgei)
    {
        if (edgeFaces[edgei].size() > 2)
        {
            FatalErrorInFunction
                << "Cannot create mesh from faceZone on mesh " << mesh.name()
                << " since edge " << pp.edges()[edgei].line(pp.localPoints())
                << " is used by more than 2 faces" << exit(FatalError);
        }
    }

    // From sub mesh to parent mesh points
    pointMap = mp;

    // Pre-size with expected number of points, faces, cells
    meshMod.setCapacity(mp.size(), edges.size(), pp.size());

    labelList faceToCell(pp.size());
    forAll(pp, facei)
    {
        faceToCell[facei] = meshMod.addCell
        (
            -1,
            -1,
            addressing[facei],   // original mesh face label
            -1,
            -1
        );
    }

    labelList pointToPoint(pp.nPoints());
    forAll(pp.localPoints(), pointi)
    {
        pointToPoint[pointi] = meshMod.addPoint
        (
            pp.localPoints()[pointi],
            mp[pointi],     // original mesh point label
            -1,
            true
        );
    }

    labelList edgeToFace(edges.size(), -1);
    face newFace(2);
    forAll(localFaces, facei)
    {
        const auto& f = localFaces[facei];
        const auto& fEdges = faceEdges[facei];

        forAll(f, fp)
        {
            const label edgei = fEdges[fp];

            if (edgeToFace[edgei] != -1)
            {
                // Edge already done
                continue;
            }

            // Get neighbour
            const auto& eFaces = edgeFaces[edgei];

            // Get edge in same order as face (already flipped if needed)
            newFace[0] = pointToPoint[f[fp]];
            newFace[1] = pointToPoint[f[f.fcIndex(fp)]];

            if (eFaces.size() == 2)
            {
                // Edge becomes internal face

                label ownCell = faceToCell[facei];
                label nbrCell = -1;
                if (eFaces[0] == facei)
                {
                    nbrCell = faceToCell[eFaces[1]];
                }
                else
                {
                    nbrCell = faceToCell[eFaces[0]];
                }

                if (nbrCell < ownCell)
                {
                    std::swap(newFace[0], newFace[1]);
                    std::swap(ownCell, nbrCell);
                }

                edgeToFace[edgei] = meshMod.addFace
                (
                    newFace,
                    ownCell,
                    nbrCell,
                    -1,             //masterPointID,
                    -1,             //masterEdgeID,
                    -1,             //masterFaceID,
                    false,          // flipFaceFlux
                    -1,             // patchID
                    -1,             // zoneID
                    false           // zoneFlip
                );
            }
            else
            {
                const label ownCell = faceToCell[facei];
                const edge meshE(mp, edges[edgei]);
                const label patchi = edgeToPatch[meshE];

                edgeToFace[edgei] = meshMod.addFace
                (
                    newFace,
                    ownCell,
                    -1,
                    -1,             //masterPointID,
                    -1,             //masterEdgeID,
                    -1,             //masterFaceID,
                    false,          // flipFaceFlux
                    patchi,         // patchID
                    -1,             // zoneID
                    false           // zoneFlip
                );
            }
        }
    }
}


template<class PatchType>
Foam::EdgeMap<Foam::label> Foam::fvMesh2DTools::edgePatch
(
    const polyMesh& mesh,
    const label defaultPatch,
    const PatchType& pp,
    const labelList& addressing
)
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    const edgeList& edges = pp.edges();
    const labelList& mp = pp.meshPoints();

    EdgeMap<label> edgeToPatch(pp.nEdges()-pp.nInternalEdges());

    // Mark (with -1) all wanted edges
    for (label edgei = pp.nInternalEdges(); edgei < pp.nEdges(); edgei++)
    {
        edgeToPatch.insert(edge(mp, edges[edgei]), -1);
    }

    //  Push marker to all
    syncTools::syncEdgeMap
    (
        mesh,
        edgeToPatch,
        maxEqOp<label>(),
        dummyTransform()
    );


    // Loop over all non-processor patches and collect patches connected
    // to edge

    const bitSet isZoneFace(addressing);

    for (label patchi = 0; patchi < pbm.nNonProcessor(); patchi++)
    {
        const auto& pp = pbm[patchi];
        forAll(pp, i)
        {
            const label facei = pp.start()+i;
            if (!isZoneFace[facei])
            {
                const face& f = mesh.faces()[facei];
                forAll(f, fp)
                {
                    const edge meshE(f.edge(fp));
                    auto fnd = edgeToPatch.find(meshE);
                    if (fnd)
                    {
                        if (fnd() == -1)
                        {
                            //Pout<< "For mesh edge:" << meshE
                            //    << " points:" << meshE.line(mesh.points())
                            //    << " found patch " << patchi
                            //    << " name:" << pbm[patchi].name()
                            //    << endl;
                            fnd() = patchi;
                        }
                        else
                        {
                            // Take max
                            maxEqOp<label>()(fnd(), patchi);
                        }
                    }
                }
            }
        }
    }

    // Max patch lists
    syncTools::syncEdgeMap
    (
        mesh,
        edgeToPatch,
        maxEqOp<label>(),
        dummyTransform()
    );

    // Override any open edge with default edge
    for (auto& patchi : edgeToPatch)
    {
        if (patchi == -1)
        {
            patchi = defaultPatch;
        }
    }

    return edgeToPatch;
}


template<class PatchType>
Foam::autoPtr<Foam::fvMesh> Foam::fvMesh2DTools::createMesh
(
    const IOobject& io,
    const fvMesh& mesh,
    const PatchType& pp,
    const labelList& addressing,
    const EdgeMap<label>& edgeToPatch,
    
    labelList& pointMap         // From point to parent mesh point
    //edgeList& faceMap           // From face to parent mesh edge. Edge given
                                // as two points to avoid unecessary addressing.
)
{
    // Constructing a 2D mesh
    polyTopoChange meshMod(mesh.boundaryMesh().size(), true, 2);

    // Insert mesh changes
    patchToMesh
    (
        mesh,
        pp,
        addressing,
        edgeToPatch,
        meshMod,
        pointMap
    );

    // Create dummy fvSchemes, fvSolution
    fvMeshTools::createDummyFvMeshFiles
    (
        mesh,
        io.name(),
        true
    );

    autoPtr<fvMesh> subMeshPtr;
    //autoPtr<mapPolyMesh> morphMap =
    meshMod.makeMesh
    (
        subMeshPtr,
        io,
        mesh
    );
    return subMeshPtr;
}


// ************************************************************************* //
