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
#include "removeCells.H"
#include "polyPatch.H"
#include "wallPolyPatch.H"
#include "edgeVertex.H"
#include "DynamicList.H"
#include "IOobject.H"
#include "OFstream.H"
#include "cuboidClipper.H"

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


bool isCutInterfaceFace
(
    const Map<label>& interfaceFaces,
    const label facei
)
{
    return interfaceFaces.found(facei);
}


labelList uniqueRegions
(
    const labelList& loopEdges,
    const labelList& edgeRegions
)
{
    DynamicList<label> regions(loopEdges.size());

    forAll(loopEdges, i)
    {
        const label regioni = edgeRegions[loopEdges[i]];

        if (regioni >= 0 && !regions.found(regioni))
        {
            regions.append(regioni);
        }
    }

    regions.shrink();
    return labelList(regions);
}


wordList regionNamesFor
(
    const labelList& regions,
    const wordList& regionNames
)
{
    wordList names(regions.size());

    forAll(regions, i)
    {
        if (regions[i] >= 0 && regions[i] < regionNames.size())
        {
            names[i] = regionNames[regions[i]];
        }
        else
        {
            names[i] = word("region" + Foam::name(regions[i]));
        }
    }

    return names;
}


void validateRegionNames(const wordList& regionNames)
{
    forAll(regionNames, regioni)
    {
        if (regionNames[regioni].empty())
        {
            FatalErrorInFunction
                << "STL region " << regioni << " has an empty name."
                << exit(FatalError);
        }

        for (label other = 0; other < regioni; ++other)
        {
            if (regionNames[other] == regionNames[regioni])
            {
                FatalErrorInFunction
                    << "STL regions " << other << " and " << regioni
                    << " have the duplicate patch name "
                    << regionNames[regioni] << "."
                    << exit(FatalError);
            }
        }
    }
}


label requireSingleRegion
(
    const label celli,
    const labelList& regions,
    const wordList& regionNames
)
{
    if (regions.empty())
    {
        FatalErrorInFunction
            << "Cut cell " << celli
            << " has no STL-region provenance for its generated cut face."
            << exit(FatalError);
    }

    if (regions.size() != 1)
    {
        FatalErrorInFunction
            << "Cut cell " << celli
            << " intersects multiple STL regions "
            << regionNamesFor(regions, regionNames) << nl
            << "Retained-inside extraction requires exact subdivision of "
            << "edge/corner cut surfaces into region-specific boundary faces. "
            << "That subdivision is not implemented safely yet, so this "
            << "utility refuses to assign the whole exposed face to one "
            << "arbitrary region."
            << exit(FatalError);
    }

    return regions[0];
}


label ensureRegionPatch
(
    polyMesh& mesh,
    const word& patchName
)
{
    const polyBoundaryMesh& boundary = mesh.boundaryMesh();
    const label existingPatchi = boundary.findPatchID(patchName);

    if (existingPatchi >= 0)
    {
        if (boundary[existingPatchi].coupled())
        {
            FatalErrorInFunction
                << "Cannot use coupled patch " << patchName
                << " for retained cut-surface extraction."
                << exit(FatalError);
        }

        return existingPatchi;
    }

    mesh.clearOut();

    polyBoundaryMesh& patches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    const label patchi = patches.size();

    dictionary patchDict;
    patchDict.set("type", wallPolyPatch::typeName);
    patchDict.set("nFaces", 0);
    patchDict.set("startFace", mesh.nFaces());

    patches.push_back
    (
        polyPatch::New
        (
            patchName,
            patchDict,
            patchi,
            patches
        )
    );

    Info<< "Added retained extraction patch " << patchName
        << " at index " << patchi << endl;

    return patchi;
}


void assignCellsToZone
(
    polyMesh& mesh,
    const labelList& cells,
    const label targetZone
)
{
    if (cells.empty())
    {
        return;
    }

    cellZoneMesh& zones = mesh.cellZones();

    if (targetZone < 0 || targetZone >= zones.size())
    {
        FatalErrorInFunction
            << "Invalid target cellZone " << targetZone
            << exit(FatalError);
    }

    boolList selected(mesh.nCells(), false);

    forAll(cells, i)
    {
        const label celli = cells[i];

        if (celli < 0 || celli >= mesh.nCells())
        {
            FatalErrorInFunction
                << "Cannot assign invalid cell " << celli
                << " to cellZone " << targetZone
                << exit(FatalError);
        }

        selected[celli] = true;
    }

    List<DynamicList<label>> members(zones.size());

    forAll(zones, zonei)
    {
        for (const label celli : zones[zonei])
        {
            if (!selected[celli])
            {
                members[zonei].append(celli);
            }
        }
    }

    forAll(cells, i)
    {
        members[targetZone].append(cells[i]);
    }

    forAll(zones, zonei)
    {
        members[zonei].shrink();
        zones[zonei].resetAddressing(labelList(members[zonei]));
    }
}


void retainedInsideExtraction
(
    polyMesh& mesh,
    const triSurface& surface,
    const Time& runTime,
    const meshCutter& cutter,
    const labelList& originalCutCells,
    const List<labelList>& cutCellRegions,
    const Map<label>& interfaceFaces,
    const Map<label>& clippedRetainedCells,
    const List<volumeType>& preClipCellTypes,
    const wordList& regionNames,
    const label keepZone
)
{
    List<volumeType> cellTypes(preClipCellTypes);

    const Map<label>& addedCells = cutter.addedCells();

    forAll(originalCutCells, i)
    {
        if (cutCellRegions[i].size() <= 1)
        {
            continue;
        }

        const label celli = originalCutCells[i];
        const label addedCelli = addedCells[celli];

        const bool originalRetained =
            clippedRetainedCells.found(celli);
        const bool addedRetained =
            clippedRetainedCells.found(addedCelli);

        if (originalRetained == addedRetained)
        {
            FatalErrorInFunction
                << "Missing retained-polyhedron classification for "
                << "multi-region cut cell " << celli
                << exit(FatalError);
        }

        cellTypes[celli] =
            originalRetained ? volumeType::INSIDE : volumeType::OUTSIDE;
        cellTypes[addedCelli] =
            addedRetained ? volumeType::INSIDE : volumeType::OUTSIDE;
    }

    DynamicList<label> cellsToRemove(mesh.nCells()/2);
    DynamicList<label> cellsToKeep(mesh.nCells()/2);
    label nInside = 0;
    label nOutside = 0;

    forAll(cellTypes, celli)
    {
        if (cellTypes[celli] == volumeType::INSIDE)
        {
            ++nInside;
            cellsToKeep.append(celli);
        }
        else if (cellTypes[celli] == volumeType::OUTSIDE)
        {
            ++nOutside;
            cellsToRemove.append(celli);
        }
        else
        {
            FatalErrorInFunction
                << "Could not classify cell " << celli
                << " for retained-inside extraction; volume type is "
                << cellTypes[celli]
                << exit(FatalError);
        }
    }

    cellsToKeep.shrink();
    cellsToRemove.shrink();

    Info<< "Retained-inside classification: " << nInside
        << " cells inside, " << nOutside << " cells outside" << endl;

    label nInsideOriginalDescendants = 0;
    label nInsideAddedDescendants = 0;

    forAll(originalCutCells, i)
    {
        const label celli = originalCutCells[i];

        if (!addedCells.found(celli))
        {
            FatalErrorInFunction
                << "Missing cutter descendant/provenance for cut cell "
                << celli << exit(FatalError);
        }

        const label addedCelli = addedCells[celli];
        if (cellTypes[celli] == volumeType::INSIDE)
        {
            ++nInsideOriginalDescendants;
        }
        if (cellTypes[addedCelli] == volumeType::INSIDE)
        {
            ++nInsideAddedDescendants;
        }

        if (cellTypes[celli] == cellTypes[addedCelli])
        {
            FatalErrorInFunction
                << "Cut descendants " << celli << " and " << addedCelli
                << " are both classified " << cellTypes[celli]
                << ". Expected one inside and one outside descendant."
                << exit(FatalError);
        }
    }

    Info<< "Cut-descendant classification: "
        << nInsideOriginalDescendants << " original descendants inside, "
        << nInsideAddedDescendants << " added descendants inside" << endl;

    assignCellsToZone(mesh, cellsToKeep, keepZone);

    labelList regionPatchIDs(regionNames.size(), -1);
    boolList expectedRegions(regionNames.size(), false);
    bool allowEmptyRegions = false;

    forAll(regionNames, regioni)
    {
        regionPatchIDs[regioni] = ensureRegionPatch(mesh, regionNames[regioni]);
    }

    forAll(cutCellRegions, i)
    {
        allowEmptyRegions =
            allowEmptyRegions || cutCellRegions[i].size() > 1;

        forAll(cutCellRegions[i], j)
        {
            expectedRegions[cutCellRegions[i][j]] = true;
        }
    }

    removeCells cellRemover(mesh);
    const labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
    labelList exposedPatchIDs(exposedFaces.size(), -1);
    labelList exposedRegionCounts(regionNames.size(), 0);

    forAll(exposedFaces, i)
    {
        const label facei = exposedFaces[i];

        if (!interfaceFaces.found(facei))
        {
            FatalErrorInFunction
                << "Retained-inside extraction would expose face " << facei
                << " without exact STL-region provenance. Only generated "
                << "single-region cut faces can become final STL patches in "
                << "this phase."
                << exit(FatalError);
        }

        const label regioni = interfaceFaces[facei];

        if (regioni < 0 || regioni >= regionPatchIDs.size())
        {
            FatalErrorInFunction
                << "Invalid STL region " << regioni
                << " for exposed face " << facei
                << exit(FatalError);
        }

        exposedPatchIDs[i] = regionPatchIDs[regioni];
        ++exposedRegionCounts[regioni];
    }

    Info<< "Retained cut-surface patch assignment:" << endl;
    forAll(exposedRegionCounts, regioni)
    {
        Info<< "  " << regionNames[regioni] << ": "
            << exposedRegionCounts[regioni] << " exposed faces" << endl;
    }

    polyTopoChange removeMod(mesh);
    cellRemover.setRefinement
    (
        cellsToRemove,
        exposedFaces,
        exposedPatchIDs,
        removeMod
    );

    autoPtr<mapPolyMesh> removeMap = removeMod.changeMesh(mesh, false);

    if (removeMap().hasMotionPoints())
    {
        mesh.movePoints(removeMap().preMotionPoints());
    }

    forAll(regionPatchIDs, regioni)
    {
        const label patchi = regionPatchIDs[regioni];

        if
        (
            mesh.boundaryMesh()[patchi].empty()
         && (!allowEmptyRegions || expectedRegions[regioni])
        )
        {
            FatalErrorInFunction
                << "Retained-inside extraction produced no faces for STL "
                << "region " << regionNames[regioni]
                << ". Exact final region-patch reconstruction failed."
                << exit(FatalError);
        }

        if (expectedRegions[regioni])
        {
            Info<< "Final STL-region patch " << regionNames[regioni]
                << " (patch " << patchi << ") contains "
                << mesh.boundaryMesh()[patchi].size() << " faces" << endl;
        }
    }
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
    const Map<label>& interfaceFaces,
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
         && !isCutInterfaceFace(interfaceFaces, facei)
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
    const Map<label>& interfaceFaces,
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
        const bool cutFace = isCutInterfaceFace(interfaceFaces, facei);
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
    argList::addBoolOption
    (
        "audit",
        "Report cut-cell and STL-region provenance without changing the mesh"
    );
    argList::addBoolOption
    (
        "extractInside",
        "Cut the surface, retain cells classified inside in -keepZone, "
        "remove outside cells, and patch exposed cut faces by STL region"
    );
    argList::addOption
    (
        "keepZone",
        "zone",
        "CellZone expected to contain the retained side during an audit or "
        "retained-inside extraction"
    );
    argList::noParallel();
    argList::noFunctionObjects();

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    const fileName surfaceName = args.get<fileName>("surface");
    const scalar endpointTol = args.getOrDefault<scalar>("tol", 1e-8);
    const bool extractInside = args.found("extractInside");
    const scalar minVolumeRatio =
        args.getOrDefault<scalar>("minVolumeRatio", extractInside ? 0.0 : 0.3);
    const bool mergeSlivers = args.found("mergeSlivers");
    const bool hasInsideZone = args.found("insideZone");
    const bool hasOutsideZone = args.found("outsideZone");
    const bool audit = args.found("audit");

    if (hasInsideZone != hasOutsideZone)
    {
        FatalErrorInFunction
            << "-insideZone and -outsideZone must be supplied together."
            << exit(FatalError);
    }

    if (extractInside && audit)
    {
        FatalErrorInFunction
            << "-extractInside changes the mesh and cannot be combined with "
            << "-audit." << exit(FatalError);
    }

    if (extractInside && hasInsideZone)
    {
        FatalErrorInFunction
            << "-extractInside assigns retained cells to -keepZone; do not "
            << "combine it with -insideZone/-outsideZone."
            << exit(FatalError);
    }

    if (extractInside && mergeSlivers)
    {
        FatalErrorInFunction
            << "-extractInside cannot currently be combined with "
            << "-mergeSlivers because sliver repair changes the cut-cell "
            << "topology before retained-side extraction."
            << exit(FatalError);
    }

    label selectedCell = -1;
    if (args.found("cell"))
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

    label keepZone = -1;
    if (args.found("keepZone"))
    {
        const word keepZoneName = args.get<word>("keepZone");
        keepZone = mesh.cellZones().findZoneID(keepZoneName);

        if (keepZone < 0)
        {
            FatalErrorInFunction
                << "Audit keepZone " << keepZoneName
                << " does not exist in the mesh" << exit(FatalError);
        }
    }
    else if (extractInside)
    {
        FatalErrorInFunction
            << "-extractInside requires -keepZone <existingCellZone>."
            << exit(FatalError);
    }

    Info<< "Reading surface " << surfaceName << endl;
    triSurface surface(surfaceName);
    triSurfaceSearch surfaceSearch(surface);
    wordList regionNames(surface.patches().size());
    forAll(regionNames, regioni)
    {
        regionNames[regioni] = surface.patches()[regioni].name();
    }
    validateRegionNames(regionNames);

    const scalarField oldCellVolumes(mesh.cellVolumes());
    const edgeList& edges = mesh.edges();
    const pointField& points = mesh.points();

    boolList edgeIsCut(mesh.nEdges(), false);
    scalarField edgeWeights(mesh.nEdges(), -GREAT);
    labelList edgeRegions(mesh.nEdges(), -1);
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
        const vector edgeVec = points[e.end()] - points[e.start()];
        List<List<pointIndexHit>> hits;

        surfaceSearch.findLineAll(starts, ends, hits);

        if (hits.size() != 1 || hits[0].empty())
        {
            continue;
        }

        DynamicList<pointIndexHit> uniqueHits(hits[0].size());
        const scalar hitTolerance =
            max(endpointTol*mag(edgeVec), SMALL);

        forAll(hits[0], hitI)
        {
            const pointIndexHit& hit = hits[0][hitI];
            bool duplicate = false;

            forAll(uniqueHits, uniqueI)
            {
                if
                (
                    mag(hit.hitPoint() - uniqueHits[uniqueI].hitPoint())
                  <= hitTolerance
                )
                {
                    duplicate = true;
                    break;
                }
            }

            if (!duplicate)
            {
                uniqueHits.append(hit);
            }
        }

        if (uniqueHits.size() != 1)
        {
            FatalErrorInFunction
                << "Edge " << edgei << " has " << hits[0].size()
                << " reported surface intersections and "
                << uniqueHits.size() << " distinct intersections; "
                << "multiple distinct hits "
                << "on one edge are unsupported." << exit(FatalError);
        }

        if (uniqueHits.size() != hits[0].size())
        {
            Info<< "Edge " << edgei << " has "
                << hits[0].size()
                << " coincident triangle hits; treating them as one "
                << "geometric intersection." << endl;
        }

        const pointIndexHit& uniqueHit = uniqueHits[0];
        const point& hit = uniqueHit.hitPoint();
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
        edgeRegions[edgei] = surface[uniqueHit.index()].region();
    }

    edgeVertex cutsEncoding(mesh);
    DynamicList<label> cutCells;
    DynamicList<labelList> loops;
    DynamicList<labelList> surfaceLoops;
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
        surfaceLoops.append(orderedEdges);
        loopWeights.append(weights);

        Info<< "Constructed one loop with " << orderedEdges.size()
            << " edge intersections for cell " << celli << endl;
    }

    cutCells.shrink();
    loops.shrink();
    surfaceLoops.shrink();
    loopWeights.shrink();

    List<labelList> cutCellRegions(cutCells.size());

    forAll(cutCells, i)
    {
        cutCellRegions[i] = uniqueRegions(surfaceLoops[i], edgeRegions);
    }

    if (cutCells.empty())
    {
        FatalErrorInFunction
            << "Surface does not cross any selected mesh cell"
            << exit(FatalError);
    }

    if (audit)
    {
        label nKeepCells = 0;
        label nBackgroundCells = 0;
        label nOtherZoneCells = 0;
        labelList regionCellCounts(regionNames.size(), 0);

        forAll(mesh.cellZones(), zonei)
        {
            const label nCells = mesh.cellZones()[zonei].size();

            Info<< "CellZone " << mesh.cellZones()[zonei].name()
                << " contains " << nCells << " cells" << endl;
        }

        for (label celli = 0; celli < mesh.nCells(); ++celli)
        {
            const label zonei = mesh.cellZones().whichZone(celli);

            if (zonei == keepZone)
            {
                ++nKeepCells;
            }
            else if (zonei < 0)
            {
                ++nBackgroundCells;
            }
            else
            {
                ++nOtherZoneCells;
            }
        }

        Info<< nl
            << "Audit summary" << nl
            << "-------------" << nl
            << "Total cells             : " << mesh.nCells() << nl
            << "Intersected cells       : " << cutCells.size() << nl
            << "Retained-zone cells     : " << nKeepCells << nl
            << "Background cells        : " << nBackgroundCells << nl
            << "Other cell-zone cells   : " << nOtherZoneCells << nl
            << "STL regions             : " << regionNames << nl
            << endl;

        forAll(cutCells, i)
        {
            const label celli = cutCells[i];
            const label zonei = mesh.cellZones().whichZone(celli);
            labelList cellRegions;

            cellRegions = cutCellRegions[i];

            Info<< "Cut cell " << celli
                << " zone " << zonei
                << " STL regions:";

            if (cellRegions.empty())
            {
                Info<< " none";
            }
            else
            {
                forAll(cellRegions, j)
                {
                    Info<< ' ' << regionNames[cellRegions[j]];
                    ++regionCellCounts[cellRegions[j]];
                }
            }

            Info<< endl;
        }

        Info<< nl << "STL-region cut-cell coverage" << nl
            << "----------------------------" << endl;

        forAll(regionCellCounts, regioni)
        {
            Info<< regionNames[regioni] << ": "
                << regionCellCounts[regioni] << " cut cells" << endl;
        }

        Info<< nl
            << "Audit completed without modifying the mesh." << endl;

        return 0;
    }

    if (extractInside)
    {
        labelList nCutCellsPerRegion(regionNames.size(), 0);
        bool hasMultiRegionCutCell = false;

        forAll(cutCells, i)
        {
            const labelList& regions = cutCellRegions[i];

            if (regions.size() == 1)
            {
                ++nCutCellsPerRegion
                [
                    requireSingleRegion(cutCells[i], regions, regionNames)
                ];
            }
            else if (regions.size() > 1)
            {
                hasMultiRegionCutCell = true;
                Info<< "Retained-inside preflight defers multi-region cell "
                    << cutCells[i] << " to the cuboid clipping component: "
                    << regionNamesFor(regions, regionNames) << endl;
            }
        }

        Info<< "Retained-inside extraction preflight accepted "
            << cutCells.size()
            << " cut cells. Region provenance:" << endl;

        forAll(nCutCellsPerRegion, regioni)
        {
            Info<< "  " << regionNames[regioni] << ": "
                << nCutCellsPerRegion[regioni] << " cut cells" << endl;
        }

        if (hasMultiRegionCutCell)
        {
            cuboidClipper::applyDirect
            (
                mesh,
                surface,
                max(endpointTol*mag(mesh.bounds().span()), SMALL),
                keepZone
            );

            ++runTime;
            mesh.write();

            Info<< "Wrote retained-inside cuboid mesh with "
                << mesh.nCells() << " cells at time "
                << runTime.timeName() << endl;

            return 0;
        }
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

    List<labelList> originalCellFaces(cutCells.size());
    List<faceList> originalCellFaceGeometry(cutCells.size());

    forAll(cutCells, i)
    {
        originalCellFaces[i] = mesh.cells()[cutCells[i]];
        originalCellFaceGeometry[i].setSize(originalCellFaces[i].size());

        forAll(originalCellFaces[i], facei)
        {
            originalCellFaceGeometry[i][facei] =
                mesh.faces()[originalCellFaces[i][facei]];
        }
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

    List<volumeType> preClipCellTypes;

    if (extractInside)
    {
        closedTriSurfaceMesh volumeSurface
        (
            IOobject
            (
                "cutCellPreClipVolumeSurface",
                runTime.constant(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            surface
        );

        volumeSurface.getVolumeType(mesh.cellCentres(), preClipCellTypes);
    }

    labelList mappedCutCells(cutCells.size(), -1);

    forAll(cutCells, i)
    {
        mappedCutCells[i] = morphMap().reverseCellMap()[cutCells[i]];

        if (mappedCutCells[i] < 0)
        {
            FatalErrorInFunction
                << "Cut cell " << cutCells[i]
                << " was not mapped after cutting."
                << exit(FatalError);
        }
    }

    Map<label> interfaceFaces;
    Map<label> clippedRetainedCells;

    forAll(cutCells, i)
    {
        if (cutCellRegions[i].size() == 1)
        {
            interfaceFaces.insert
            (
                cutter.addedFaces()[mappedCutCells[i]],
                cutCellRegions[i][0]
            );
        }
    }

    cuboidClipper::apply
    (
        mesh,
        surface,
        cutter,
        cutCells,
        mappedCutCells,
        originalCellFaces,
        originalCellFaceGeometry,
        morphMap(),
        cutCellRegions,
        max(endpointTol*mag(mesh.bounds().span()), SMALL),
        extractInside,
        interfaceFaces,
        clippedRetainedCells
    );

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

        forAll(cutCells, i)
        {
            const label oldCelli = cutCells[i];
            const label celli = mappedCutCells[i];
            const label addedCelli = addedCells[celli];
            const scalar parentVolume = oldCellVolumes[oldCelli];
            const scalar originalRatio =
                mesh.cellVolumes()[celli]/parentVolume;
            const scalar addedRatio =
                mesh.cellVolumes()[addedCelli]/parentVolume;
            const scalar minimumRatio = min(originalRatio, addedRatio);

            if (minimumRatio < minVolumeRatio)
            {
                const label sliverCell =
                    originalRatio < addedRatio ? celli : addedCelli;

                reportMergeCandidates(mesh, interfaceFaces, sliverCell);

                if (mergeSlivers)
                {
                    label mergeCell = -1;
                    label mergeFace = -1;

                    if
                    (
                        !findMergeCandidate
                        (
                            mesh,
                            interfaceFaces,
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

    if (extractInside)
    {
        retainedInsideExtraction
        (
            mesh,
            surface,
            runTime,
            cutter,
            mappedCutCells,
            cutCellRegions,
            interfaceFaces,
            clippedRetainedCells,
            preClipCellTypes,
            regionNames,
            keepZone
        );
    }

    ++runTime;
    mesh.write();

    Info<< "Wrote " << (extractInside ? "retained-inside " : "")
        << "cut mesh with " << mesh.nCells()
        << " cells at time " << runTime.timeName() << endl;

    return 0;
}

// ************************************************************************* //
