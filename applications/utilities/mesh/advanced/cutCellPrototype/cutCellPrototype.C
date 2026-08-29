/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Description
        Split one hexahedral cell using the intersection of a closed
        triangulated surface with the cell edges.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "closedTriSurfaceMesh.H"
#include "cellCuts.H"
#include "meshCutter.H"
#include "removeFaces.H"
#include "edgeVertex.H"
#include "DynamicList.H"
#include "IOobject.H"
#include "OFstream.H"

using namespace Foam;

namespace
{

bool faceContainsEdge
(
    const polyMesh& mesh,
    const label facei,
    const label edgei
)
{
    return mesh.faceEdges()[facei].found(edgei);
}


label nextCutEdge
(
    const polyMesh& mesh,
    const labelList& cellFaces,
    const label currentEdge,
    const label previousEdge,
    const labelList& cutEdges
)
{
    forAll(cellFaces, i)
    {
        const label facei = cellFaces[i];

        if (!faceContainsEdge(mesh, facei, currentEdge))
        {
            continue;
        }

        label other = -1;
        label nCut = 0;

        forAll(cutEdges, j)
        {
            if (faceContainsEdge(mesh, facei, cutEdges[j]))
            {
                ++nCut;

                if (cutEdges[j] != currentEdge)
                {
                    other = cutEdges[j];
                }
            }
        }

        if (nCut == 2 && other != previousEdge)
        {
            return other;
        }
    }

    return -1;
}


bool buildLoop
(
    const polyMesh& mesh,
    const label celli,
    const labelList& cutEdges,
    labelList& loopEdges
)
{
    const labelList& cellFaces = mesh.cells()[celli];

    if (cutEdges.size() < 3)
    {
        return false;
    }

    forAll(cellFaces, i)
    {
        label nCut = 0;

        forAll(cutEdges, j)
        {
            nCut += faceContainsEdge(mesh, cellFaces[i], cutEdges[j]);
        }

        if (nCut > 2)
        {
            WarningInFunction
                << "Cell " << celli << " has " << nCut
                << " intersections on face " << cellFaces[i]
                << "; only one cut loop with two points per face is supported."
                << endl;
            return false;
        }
    }

    DynamicList<label> ordered(cutEdges.size());
    label previous = -1;
    label current = cutEdges[0];

    for (label step = 0; step < cutEdges.size(); ++step)
    {
        if (ordered.found(current))
        {
            break;
        }

        ordered.append(current);

        const label next =
            nextCutEdge(mesh, cellFaces, current, previous, cutEdges);

        if (next < 0)
        {
            return false;
        }

        previous = current;
        current = next;
    }

    if (current != ordered[0] || ordered.size() != cutEdges.size())
    {
        WarningInFunction
            << "Cell " << celli
            << " does not produce one closed intersection loop." << endl;
        return false;
    }

    loopEdges.transfer(ordered);
    return true;
}


label commonFace
(
    const polyMesh& mesh,
    const label celli,
    const label neighbour
)
{
    const labelList& cellFaces = mesh.cells()[celli];
    const labelList& neighbourFaces = mesh.cells()[neighbour];

    forAll(cellFaces, i)
    {
        if (neighbourFaces.found(cellFaces[i]))
        {
            return cellFaces[i];
        }
    }

    return -1;
}


bool isCutFace
(
    const Map<label>& addedFaces,
    const label facei
)
{
    forAllConstIters(addedFaces, iter)
    {
        if (iter.val() == facei)
        {
            return true;
        }
    }

    return false;
}


struct mergeQuality
{
    bool valid;
    scalar volume;
    scalar minPyramidRatio;
    label concaveFaces;
};


mergeQuality evaluateMergedCell
(
    const polyMesh& mesh,
    const label celli,
    const label neighbour,
    const label removedFace
)
{
    mergeQuality quality{false, -GREAT, -GREAT, labelMax};
    DynamicList<label> boundaryFaces;

    forAll(mesh.cells()[celli], i)
    {
        const label facei = mesh.cells()[celli][i];

        if (facei != removedFace)
        {
            boundaryFaces.append(facei);
        }
    }

    forAll(mesh.cells()[neighbour], i)
    {
        const label facei = mesh.cells()[neighbour][i];

        if (facei != removedFace && !boundaryFaces.found(facei))
        {
            boundaryFaces.append(facei);
        }
    }

    scalar mergedVolume = mesh.cellVolumes()[celli]
        + mesh.cellVolumes()[neighbour];
    point mergedCentre =
    (
        mesh.cellCentres()[celli]*mesh.cellVolumes()[celli]
      + mesh.cellCentres()[neighbour]*mesh.cellVolumes()[neighbour]
    )/max(mergedVolume, VSMALL);

    const pointField& faceCentres = mesh.faceCentres();
    const vectorField& faceAreas = mesh.faceAreas();
    const labelList& faceOwner = mesh.faceOwner();
    scalar signedVolume = 0;
    scalar minPyramidVolume = GREAT;
    label concaveFaces = 0;

    forAll(boundaryFaces, i)
    {
        const label facei = boundaryFaces[i];
        vector faceNormal = faceAreas[facei];

        if
        (
            faceOwner[facei] != celli
         && faceOwner[facei] != neighbour
        )
        {
            FatalErrorInFunction
                << "Face " << facei << " is not owned by either cell "
                << celli << " or " << neighbour << exit(FatalError);
        }

        if (faceOwner[facei] != celli)
        {
            faceNormal *= -1;
        }

        if
        (
            (faceNormal & (faceCentres[facei] - mergedCentre)) < 0
        )
        {
            faceNormal *= -1;
        }

        const vector unitNormal = faceNormal
            /max(mag(faceNormal), VSMALL);
        const scalar pyramidVolume =
            (faceNormal & (faceCentres[facei] - mergedCentre))/3;

        signedVolume += pyramidVolume;
        minPyramidVolume = min(minPyramidVolume, pyramidVolume);

        forAll(boundaryFaces, j)
        {
            if (i == j)
            {
                continue;
            }

            vector centreDirection = faceCentres[boundaryFaces[j]]
                - faceCentres[facei];
            centreDirection /= max(mag(centreDirection), VSMALL);

            if ((centreDirection & unitNormal) > 1e-6)
            {
                ++concaveFaces;
                break;
            }
        }
    }

    quality.volume = signedVolume;
    quality.minPyramidRatio =
        minPyramidVolume/max(mag(signedVolume), VSMALL);
    quality.concaveFaces = concaveFaces;
    quality.valid = signedVolume > VSMALL;

    return quality;
}


bool findMergeCandidate
(
    const polyMesh& mesh,
    const Map<label>& addedFaces,
    const label sliverCell,
    label& mergeCell,
    label& mergeFace
)
{
    const label sliverZone = mesh.cellZones().whichZone(sliverCell);
    const labelList& neighbours = mesh.cellCells()[sliverCell];
    label bestConcave = labelMax;
    scalar bestPyramidRatio = -GREAT;
    scalar bestArea = -1;
    label bestNeighbour = -1;
    label bestFace = -1;
    mergeQuality bestQuality{false, -GREAT, -GREAT, labelMax};

    forAll(neighbours, i)
    {
        const label neighbour = neighbours[i];
        const label facei = commonFace(mesh, sliverCell, neighbour);
        const bool eligible =
            mesh.isInternalFace(facei)
         && sliverZone == mesh.cellZones().whichZone(neighbour)
         && !isCutFace(addedFaces, facei)
         && mesh.faceZones().whichZone(facei) < 0;

        if (!eligible)
        {
            continue;
        }

        const mergeQuality quality =
            evaluateMergedCell(mesh, sliverCell, neighbour, facei);

        if (!quality.valid)
        {
            Info<< "Rejected merge candidate " << neighbour
                << " across face " << facei
                << ": merged volume " << quality.volume
                << ", minimum face-pyramid ratio "
                << quality.minPyramidRatio << endl;
            continue;
        }

        const label concave = quality.concaveFaces;
        const scalar area = mag(mesh.faceAreas()[facei]);

        if
        (
            concave < bestConcave
         || (
                concave == bestConcave
             && (
                    quality.minPyramidRatio > bestPyramidRatio
                 || (
                        quality.minPyramidRatio == bestPyramidRatio
                     && (
                    area > bestArea
                 || (area == bestArea && neighbour < bestNeighbour)
                        )
                    )
                )
            )
        )
        {
            bestConcave = concave;
            bestPyramidRatio = quality.minPyramidRatio;
            bestArea = area;
            bestNeighbour = neighbour;
            bestFace = facei;
            bestQuality = quality;
        }
    }

    if (bestNeighbour < 0)
    {
        return false;
    }

    mergeCell = bestNeighbour;
    mergeFace = bestFace;

    Info<< "Selected merge candidate " << mergeCell
        << " across face " << mergeFace
        << " with merged-cell concavity score " << bestConcave
    << ", minimum face-pyramid ratio "
    << bestQuality.minPyramidRatio
    << ", volume " << bestQuality.volume
    << " and shared-face area " << bestArea << endl;

    return true;
}


void reportMergeCandidates
(
    const polyMesh& mesh,
    const Map<label>& addedFaces,
    const label sliverCell
)
{
    const label sliverZone = mesh.cellZones().whichZone(sliverCell);
    const labelList& neighbours = mesh.cellCells()[sliverCell];

    Info<< "Sliver cell " << sliverCell
        << " is in cellZone " << sliverZone
        << "; merge candidates:" << endl;

    if (neighbours.empty())
    {
        Info<< "  none (cell has no neighbours)" << endl;
        return;
    }

    forAll(neighbours, i)
    {
        const label neighbour = neighbours[i];
        const label facei = commonFace(mesh, sliverCell, neighbour);
        const label neighbourZone = mesh.cellZones().whichZone(neighbour);
        const label patchi = mesh.boundaryMesh().whichPatch(facei);
        const label faceZonei = mesh.faceZones().whichZone(facei);
        const bool cutFace = isCutFace(addedFaces, facei);
        const bool sameZone = sliverZone == neighbourZone;
        const bool internalFace = mesh.isInternalFace(facei);
        const bool eligible =
            internalFace
         && sameZone
         && !cutFace
         && faceZonei < 0;

        Info<< "  neighbour " << neighbour
            << ", face " << facei
            << ", neighbourZone " << neighbourZone
            << ", internal " << internalFace
            << ", cutInterface " << cutFace
            << ", patch " << patchi
            << ", faceZone " << faceZonei
            << ", eligible " << eligible
            << endl;
    }
}


void verifyCutCellZones
(
    const polyMesh& mesh,
    const meshCutter& cutter,
    const labelList& originalCutCells,
    const mapPolyMesh& morphMap
)
{
    const Map<label>& addedCells = cutter.addedCells();

    forAll(originalCutCells, i)
    {
        const label oldCelli = originalCutCells[i];
        const label celli = morphMap.reverseCellMap()[oldCelli];
        const label addedCelli = addedCells[celli];
        const label cellZone = mesh.cellZones().whichZone(celli);
        const label addedZone = mesh.cellZones().whichZone(addedCelli);

        if (cellZone != addedZone)
        {
            FatalErrorInFunction
                << "Cell-zone inheritance failed for cut cell " << celli
                << " and added cell " << addedCelli
                << ": zones are " << cellZone << " and " << addedZone
                << exit(FatalError);
        }

        Info<< "Cut cell " << celli << " and added cell " << addedCelli
            << " inherit cellZone " << cellZone << endl;
    }
}


void applyZoneOverrides
(
    polyMesh& mesh,
    const triSurface& surface,
    const Time& runTime,
    const meshCutter& cutter,
    const labelList& originalCutCells,
    const mapPolyMesh& morphMap,
    const label insideZone,
    const label outsideZone
)
{
    closedTriSurfaceMesh volumeSurface
    (
        IOobject
        (
            "cutCellVolumeSurface",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        surface
    );

    const Map<label>& addedCells = cutter.addedCells();
    labelList overrideCells(2*originalCutCells.size(), -1);
    labelList overrideZones(2*originalCutCells.size(), -1);
    label nOverrides = 0;

    forAll(originalCutCells, i)
    {
        const label celli = morphMap.reverseCellMap()[originalCutCells[i]];
        const label addedCelli = addedCells[celli];
        pointField samples(2);
        samples[0] = mesh.cellCentres()[celli];
        samples[1] = mesh.cellCentres()[addedCelli];

        List<volumeType> types;
        volumeSurface.getVolumeType(samples, types);

        forAll(samples, samplei)
        {
            label targetZone = -1;

            if (types[samplei] == volumeType::INSIDE)
            {
                targetZone = insideZone;
            }
            else if (types[samplei] == volumeType::OUTSIDE)
            {
                targetZone = outsideZone;
            }
            else
            {
                FatalErrorInFunction
                    << "Could not classify cut-cell representative point "
                    << samplei << " for cell " << celli
                    << "; volume type is " << types[samplei]
                    << exit(FatalError);
            }

            overrideCells[nOverrides] =
                samplei == 0 ? celli : addedCelli;
            overrideZones[nOverrides] = targetZone;
            ++nOverrides;
        }
    }

    overrideCells.setSize(nOverrides);
    overrideZones.setSize(nOverrides);

    cellZoneMesh& zones = mesh.cellZones();
    List<DynamicList<label>> members(zones.size());
    boolList overridden(mesh.nCells(), false);

    forAll(overrideCells, i)
    {
        overridden[overrideCells[i]] = true;
    }

    forAll(zones, zonei)
    {
        for (const label celli : zones[zonei])
        {
            if (!overridden[celli])
            {
                members[zonei].append(celli);
            }
        }
    }

    forAll(overrideCells, i)
    {
        members[overrideZones[i]].append(overrideCells[i]);
    }

    forAll(zones, zonei)
    {
        members[zonei].shrink();
        zones[zonei].resetAddressing(labelList(members[zonei]));
    }

    forAll(overrideCells, i)
    {
        Info<< "Geometry override assigned cell " << overrideCells[i]
            << " to cellZone " << overrideZones[i] << endl;
    }
}

}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Prototype: split intersected hexahedral cells using a triangulated-surface cut"
    );
    argList::addOption
    (
        "cell",
        "label",
        "Cell label to cut (default: all intersected cells)"
    );
    argList::addOption("surface", "file", "STL/OBJ surface file");
    argList::addOption
    (
        "tol",
        "scalar",
        "Endpoint tolerance as a fraction of edge length (default 1e-8)"
    );
    argList::addOption
    (
        "minVolumeRatio",
        "scalar",
        "Reject a cut if either side is below this parent-volume ratio "
        "(default 0.3; use 0 to disable)"
    );
    argList::addBoolOption
    (
        "mergeSlivers",
        "Merge a rejected sliver into a safe adjacent cell"
    );
    argList::addOption
    (
        "insideZone",
        "zone",
        "Existing cellZone to assign to cut descendants classified inside"
    );
    argList::addOption
    (
        "outsideZone",
        "zone",
        "Existing cellZone to assign to cut descendants classified outside"
    );
    argList::noParallel();
    argList::noFunctionObjects();

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    const fileName surfaceName = args.get<fileName>("surface");
    const scalar endpointTol = args.getOrDefault<scalar>("tol", 1e-8);
    const scalar minVolumeRatio =
        args.getOrDefault<scalar>("minVolumeRatio", 0.3);
    const bool mergeSlivers = args.found("mergeSlivers");
    const bool hasInsideZone = args.optionFound("insideZone");
    const bool hasOutsideZone = args.optionFound("outsideZone");

    if (hasInsideZone != hasOutsideZone)
    {
        FatalErrorInFunction
            << "-insideZone and -outsideZone must be supplied together."
            << exit(FatalError);
    }

    label selectedCell = -1;
    if (args.optionFound("cell"))
    {
        selectedCell = args.get<label>("cell");

        if (selectedCell < 0 || selectedCell >= mesh.nCells())
        {
            FatalErrorInFunction
                << "Cell " << selectedCell << " is outside [0, "
                << mesh.nCells() << ')' << exit(FatalError);
        }
    }

    if (endpointTol <= 0 || endpointTol >= 0.5)
    {
        FatalErrorInFunction
            << "Endpoint tolerance must be in (0, 0.5), got "
            << endpointTol << exit(FatalError);
    }

    if (minVolumeRatio < 0 || minVolumeRatio >= 1.0)
    {
        FatalErrorInFunction
            << "Minimum volume ratio must be in [0, 1), got "
            << minVolumeRatio << exit(FatalError);
    }

    Info<< "Reading surface " << surfaceName << endl;
    triSurface surface(surfaceName);
    triSurfaceSearch surfaceSearch(surface);
    const scalarField oldCellVolumes(mesh.cellVolumes());
    const edgeList& edges = mesh.edges();
    const pointField& points = mesh.points();

    boolList edgeIsCut(mesh.nEdges(), false);
    scalarField edgeWeights(mesh.nEdges(), -GREAT);
    boolList edgeToCheck(mesh.nEdges(), selectedCell < 0);

    if (selectedCell >= 0)
    {
        const labelList& selectedCellEdges = mesh.cellEdges()[selectedCell];

        forAll(selectedCellEdges, i)
        {
            edgeToCheck[selectedCellEdges[i]] = true;
        }
    }

    forAll(edges, edgei)
    {
        if (!edgeToCheck[edgei])
        {
            continue;
        }

        const edge& e = edges[edgei];
        const pointField starts(1, points[e.start()]);
        const pointField ends(1, points[e.end()]);
        List<List<pointIndexHit>> hits;

        surfaceSearch.findLineAll(starts, ends, hits);

        if (hits.size() != 1 || hits[0].empty())
        {
            continue;
        }

        if (hits[0].size() != 1)
        {
            FatalErrorInFunction
                << "Edge " << edgei << " has " << hits[0].size()
                << " surface intersections; multiple hits "
                << "on one edge are unsupported." << exit(FatalError);
        }

        const point& hit = hits[0][0].hitPoint();
        const vector edgeVec = points[e.end()] - points[e.start()];
        const scalar edgeLengthSqr = magSqr(edgeVec);
        const scalar weight = ((hit - points[e.start()]) & edgeVec)
            / edgeLengthSqr;

        if (weight <= endpointTol || weight >= 1 - endpointTol)
        {
            FatalErrorInFunction
                << "Edge " << edgei << " intersects too close to an existing "
                << "mesh vertex. "
                << "Vertex cuts are not supported by this prototype."
                << exit(FatalError);
        }

        edgeIsCut[edgei] = true;
        edgeWeights[edgei] = weight;
    }

    edgeVertex cutsEncoding(mesh);
    DynamicList<label> cutCells;
    DynamicList<labelList> loops;
    DynamicList<scalarField> loopWeights;

    for (label celli = 0; celli < mesh.nCells(); ++celli)
    {
        if (selectedCell >= 0 && celli != selectedCell)
        {
            continue;
        }

        const labelList& cellEdges = mesh.cellEdges()[celli];
        DynamicList<label> cellCutEdges(cellEdges.size());

        forAll(cellEdges, i)
        {
            if (edgeIsCut[cellEdges[i]])
            {
                cellCutEdges.append(cellEdges[i]);
            }
        }

        if (cellCutEdges.empty())
        {
            continue;
        }

        labelList orderedEdges;
        if (!buildLoop(mesh, celli, cellCutEdges, orderedEdges))
        {
            FatalErrorInFunction
                << "Could not construct one supported intersection loop for "
                << "cell " << celli << exit(FatalError);
        }

        labelList loop(orderedEdges.size());
        scalarField weights(orderedEdges.size(), -GREAT);

        forAll(orderedEdges, i)
        {
            const label edgei = orderedEdges[i];
            loop[i] = cutsEncoding.edgeToEVert(edgei);
            weights[i] = edgeWeights[edgei];
        }

        cutCells.append(celli);
        loops.append(loop);
        loopWeights.append(weights);

        Info<< "Constructed one loop with " << orderedEdges.size()
            << " edge intersections for cell " << celli << endl;
    }

    cutCells.shrink();
    loops.shrink();
    loopWeights.shrink();

    if (cutCells.empty())
    {
        FatalErrorInFunction
            << "Surface does not cross any selected mesh cell"
            << exit(FatalError);
    }

    labelListList cellLoops(cutCells.size());
    List<scalarField> allLoopWeights(cutCells.size());

    forAll(cutCells, i)
    {
        cellLoops[i] = loops[i];
        allLoopWeights[i] = loopWeights[i];
    }

    Info<< "Constructed " << cutCells.size() << " supported cut loop"
        << (cutCells.size() == 1 ? "" : "s") << endl;

    cellCuts cuts(mesh, cutCells, cellLoops, allLoopWeights);

    if (cuts.nLoops() != cutCells.size())
    {
        FatalErrorInFunction
            << "OpenFOAM accepted " << cuts.nLoops() << " of "
            << cutCells.size() << " supplied cut descriptions"
            << exit(FatalError);
    }

    polyTopoChange meshMod(mesh);
    meshCutter cutter(mesh);
    cutter.setRefinement(cuts, meshMod);

    autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

    if (morphMap().hasMotionPoints())
    {
        mesh.movePoints(morphMap().preMotionPoints());
    }

    cutter.updateMesh(morphMap());

    verifyCutCellZones(mesh, cutter, cutCells, morphMap());

    if (hasInsideZone)
    {
        const word insideZoneName = args.get<word>("insideZone");
        const word outsideZoneName = args.get<word>("outsideZone");
        const label insideZone = mesh.cellZones().findZoneID(insideZoneName);
        const label outsideZone = mesh.cellZones().findZoneID(outsideZoneName);

        if (insideZone < 0 || outsideZone < 0)
        {
            FatalErrorInFunction
                << "Geometry override zones must already exist. Requested "
                << insideZoneName << " and " << outsideZoneName
                << exit(FatalError);
        }

        applyZoneOverrides
        (
            mesh,
            surface,
            runTime,
            cutter,
            cutCells,
            morphMap(),
            insideZone,
            outsideZone
        );
    }

    if (minVolumeRatio > 0)
    {
        const Map<label>& addedCells = cutter.addedCells();
        const Map<label>& addedFaces = cutter.addedFaces();

        forAll(cutCells, i)
        {
            const label celli = cutCells[i];
            const label addedCelli = addedCells[celli];
            const scalar parentVolume = oldCellVolumes[celli];
            const scalar originalRatio =
                mesh.cellVolumes()[celli]/parentVolume;
            const scalar addedRatio =
                mesh.cellVolumes()[addedCelli]/parentVolume;
            const scalar minimumRatio = min(originalRatio, addedRatio);

            if (minimumRatio < minVolumeRatio)
            {
                const label sliverCell =
                    originalRatio < addedRatio ? celli : addedCelli;

                reportMergeCandidates(mesh, addedFaces, sliverCell);

                if (mergeSlivers)
                {
                    label mergeCell = -1;
                    label mergeFace = -1;

                    if
                    (
                        !findMergeCandidate
                        (
                            mesh,
                            addedFaces,
                            sliverCell,
                            mergeCell,
                            mergeFace
                        )
                    )
                    {
                        FatalErrorInFunction
                            << "No safe same-zone merge candidate exists for "
                            << "sliver cell " << sliverCell
                            << exit(FatalError);
                    }

                    removeFaces faceRemover(mesh, 1.0);
                    labelList piercedFaces(1, mergeFace);
                    labelList cellRegion;
                    labelList cellRegionMaster;
                    labelList allPiercedFaces;

                    const label nRegions =
                        faceRemover.compatibleRemoves
                        (
                            piercedFaces,
                            cellRegion,
                            cellRegionMaster,
                            allPiercedFaces
                        );

                    if (nRegions != 1)
                    {
                        FatalErrorInFunction
                            << "Expected one merge region for face "
                            << mergeFace << ", found " << nRegions
                            << exit(FatalError);
                    }

                    polyTopoChange mergeMod(mesh);
                    faceRemover.setRefinement
                    (
                        allPiercedFaces,
                        cellRegion,
                        cellRegionMaster,
                        mergeMod
                    );

                    autoPtr<mapPolyMesh> mergeMap =
                        mergeMod.changeMesh(mesh, false);

                    if (mergeMap().hasMotionPoints())
                    {
                        mesh.movePoints(mergeMap().preMotionPoints());
                    }

                    const label masterCell =
                        cellRegionMaster[cellRegion[sliverCell]];

                    Info<< "Merged sliver cell " << sliverCell
                        << " with cell " << mergeCell
                        << " across face " << mergeFace
                        << "; master cell is " << masterCell << endl;

                    ++runTime;
                    mesh.write();

                    Info<< "Wrote repaired cut mesh with "
                        << mesh.nCells() << " cells at time "
                        << runTime.timeName() << endl;

                    return 0;
                }

                FatalErrorInFunction
                    << "Cut cell " << celli << " creates a sliver: "
                    << "side volume ratios are " << originalRatio << " and "
                    << addedRatio << ", minimum " << minimumRatio
                    << " is below the configured threshold "
                    << minVolumeRatio
                    << ". Enable -mergeSlivers to attempt a quality-ranked "
                    << "same-zone repair."
                    << exit(FatalError);
            }
        }
    }

    ++runTime;
    mesh.write();

    Info<< "Wrote cut mesh with " << mesh.nCells()
        << " cells at time " << runTime.timeName() << endl;

    return 0;
}

// ************************************************************************* //
