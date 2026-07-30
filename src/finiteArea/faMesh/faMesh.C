/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024-2026 OpenCFD Ltd.
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

#include "faMesh.H"
#include "faMeshesRegistry.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "uindirectPrimitivePatch.H"
#include "polyBoundaryMesh.H"
#include "polyPatch.H"
#include "processorPolyPatch.H"
#include "fvBoundaryMesh.H"
#include "fvPatch.H"
#include "surfaceFields.H"
#include "ListOps.H"
#include "dictionary.H"
#include "IOdictionary.H"
#include "wordRes.H"
#include "EdgeMap.H"
#include "syncTools.H"
#include "dummyTransform.H"
#include "bitSet.H"
#include "labelIOList.H"
#include "vectorIOField.H"
#include "HashSet.H"
#include "DynamicList.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::faMesh::prefix_("finite-area");

Foam::word Foam::faMesh::meshSubDir("faMesh");


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Calculate per boundary-edge the parent polyPatch (adapted from
// fvMesh2DTools::edgePatch). Currently does not handle side-ness of faceZone.
static EdgeMap<label> edgePatchFn
(
    const polyMesh& mesh,
    const label defaultPatch,
    const uindirectPrimitivePatch& pp,
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

    // Push marker to all
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
        const auto& patch = pbm[patchi];
        forAll(patch, i)
        {
            const label facei = patch.start()+i;
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
                            fnd() = patchi;
                        }
                        else
                        {
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


// Per (parent-local) edge, the neighbour rank of the parent processor patch it
// lies on. Used to turn the reduced-order region's inter-domain boundary edges
// (whose adjacent area-face is on another processor) into genuine processor
// patches. Only edges carried by a parent processorPolyPatch appear; if an edge
// lies on more than one (corner), the higher neighbour rank wins (deterministic).
static EdgeMap<label> edgeToNeighbProcFn(const polyMesh& mesh)
{
    // Scans the parent boundary for processorPolyPatches. This works both in a
    // live parallel run (the rank's own decomposition) and in a serial
    // decomposePar, where each processor mesh is loaded in turn and carries
    // topological processor patches read from disk.
    EdgeMap<label> edgeToProc;

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    const faceList& faces = mesh.faces();

    forAll(pbm, patchi)
    {
        const auto* ppp = isA<processorPolyPatch>(pbm[patchi]);
        if (!ppp)
        {
            continue;
        }
        const label nbrProc = ppp->neighbProcNo();

        const label start = pbm[patchi].start();
        const label sz = pbm[patchi].size();
        for (label facei = start; facei < start + sz; ++facei)
        {
            const face& f = faces[facei];
            forAll(f, fp)
            {
                const edge meshE(f.edge(fp));
                auto iter = edgeToProc.find(meshE);
                if (iter.good())
                {
                    maxEqOp<label>()(iter.val(), nbrProc);
                }
                else
                {
                    edgeToProc.insert(meshE, nbrProc);
                }
            }
        }
    }

    return edgeToProc;
}


// This processor's rank as recorded on the parent mesh. In a live parallel run
// this equals UPstream::myProcNo(); in a serial decomposePar it is the rank of
// the decomposed processor mesh currently being processed (read from the parent
// processor patches, which decomposePar constructs with explicit myProcNo).
static label parentMyProcFn(const polyMesh& mesh)
{
    for (const polyPatch& pp : mesh.boundaryMesh())
    {
        const auto* ppp = isA<processorPolyPatch>(pp);
        if (ppp)
        {
            return ppp->myProcNo();
        }
    }
    return UPstream::myProcNo();
}


// Insert changes to extract the patch faces as reduced-order (2D) cells
// (adapted from fvMesh2DTools::patchToMesh).
static void patchToMeshFn
(
    const polyMesh& mesh,
    const uindirectPrimitivePatch& pp,
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

    forAll(edgeFaces, edgei)
    {
        if (edgeFaces[edgei].size() > 2)
        {
            FatalErrorInFunction
                << "Cannot create area mesh from patch on mesh " << mesh.name()
                << " since edge " << pp.edges()[edgei].line(pp.localPoints())
                << " is used by more than 2 faces" << exit(FatalError);
        }
    }

    // From sub mesh to parent mesh points
    pointMap = mp;

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

    // Collected boundary face, sorted before being added so that faces on the
    // same (processor) patch appear in a globally-consistent order across ranks.
    // The sort key is the edge's parent point coordinates, which both ranks
    // agree on bit-for-bit for a shared edge.
    struct boundaryFace
    {
        label patchi;
        label ownCell;
        face f;
        point pLo;   // ordered endpoint coordinates (lexicographic)
        point pHi;
    };

    const pointField& parentPts = mesh.points();

    const auto lexLess = [](const point& a, const point& b) noexcept
    {
        if (a.x() != b.x()) return a.x() < b.x();
        if (a.y() != b.y()) return a.y() < b.y();
        return a.z() < b.z();
    };

    DynamicList<boundaryFace> bfaces(edges.size() - pp.nInternalEdges());

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
                continue;
            }
            edgeToFace[edgei] = 1;  // mark done (id unused within this function)

            const auto& eFaces = edgeFaces[edgei];

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

                meshMod.addFace
                (
                    newFace,
                    ownCell,
                    nbrCell,
                    -1, -1, -1,
                    false,
                    -1,
                    -1,
                    false
                );
            }
            else
            {
                const edge meshE(mp, edges[edgei]);
                const label patchi = edgeToPatch[meshE];

                point pa = parentPts[meshE[0]];
                point pb = parentPts[meshE[1]];
                if (!lexLess(pa, pb))
                {
                    std::swap(pa, pb);
                }

                bfaces.append
                (
                    boundaryFace{patchi, faceToCell[facei], newFace, pa, pb}
                );
            }
        }
    }

    Foam::sort
    (
        bfaces,
        [&](const boundaryFace& a, const boundaryFace& b) noexcept
        {
            if (a.patchi != b.patchi) return a.patchi < b.patchi;
            if (a.pLo != b.pLo) return lexLess(a.pLo, b.pLo);
            return lexLess(a.pHi, b.pHi);
        }
    );

    for (const boundaryFace& bf : bfaces)
    {
        meshMod.addFace
        (
            bf.f,
            bf.ownCell,
            -1,
            -1, -1, -1,
            false,
            bf.patchi,
            -1,
            false
        );
    }
}


// Lightweight finiteArea boundary-patch definition, used to map each
// reduced-order boundary face (parent edge) to a named/typed patch.
// The owner is the parent polyPatch of the area face; the neighbour is the
// adjacent parent polyPatch across the edge (from edgePatchFn).
struct faPatchDef
{
    word name_;
    word type_;
    label ownerPolyPatchId_ = -1;
    label neighPolyPatchId_ = -1;

    //- For a "processor" patch: the neighbour rank across the inter-domain
    //  edges (otherwise -1)
    label procNeighbour_ = -1;

    //- Ordered match with (owner,neighbour) parent-patch ids (cf.
    //  faPatchData::matchPatchPair): 0=none, 2=neighbour, 3=owner+neighbour
    int matchPatchPair(const label ownerP, const label neighP) const noexcept
    {
        int ret = 0;
        if (ownerP >= 0 && ownerP == ownerPolyPatchId_) ret |= 1;
        if (neighP >= 0 && neighP == neighPolyPatchId_) ret |= 2;
        return ret;
    }
};

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::faMesh::readFaceLabels
(
    const word& areaName,
    const polyMesh& pMesh
)
{
    const word regionSub =
    (
        (areaName.empty() || areaName == polyMesh::defaultRegion)
      ? word::null : areaName
    );

    const fileName localDir =
    (
        regionSub.empty()
      ? (faMesh::prefix()/meshSubDir)
      : (faMesh::prefix()/regionSub/meshSubDir)
    );

    const word inst(pMesh.time().findInstance(localDir, "faceLabels"));

    labelIOList labels
    (
        IOobject
        (
            "faceLabels",
            inst,
            localDir,
            pMesh.thisDb(),
            IOobjectOption::MUST_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        )
    );

    return labelList(std::move(static_cast<labelList&>(labels)));
}


void Foam::faMesh::readCouplingAddressing()
{
    // Only meaningful for a parent-coupled region that has not already had its
    // addressing populated (e.g. the component/definition constructors).
    if (!meshFaces_.empty() || !hasParent())
    {
        return;
    }

    const refPtr<fvMesh> pRef(parentMesh());
    const polyMesh& pMesh = pRef();

    // The single (default) area region persists its coupling addressing
    // directly under <prefix>/faMesh (no <region> subdir). The solvers name
    // that region faMesh::prefix() ("finite-area") via createFaMesh.H, whereas
    // makeFaMesh uses polyMesh::defaultRegion; treat both as the default so
    // the read path matches the on-disk location written by writeObject().
    const word regionSub =
    (
        (
            areaName_.empty()
         || areaName_ == polyMesh::defaultRegion
         || areaName_ == faMesh::prefix()
        )
      ? word::null : areaName_
    );

    const fileName localDir =
    (
        regionSub.empty()
      ? (faMesh::prefix()/meshSubDir)
      : (faMesh::prefix()/regionSub/meshSubDir)
    );

    // faceLabels: written by writeObject() alongside the reduced-order mesh.
    // Required for a parent-coupled read; missing => leave empty (standalone).
    {
        IOobject io
        (
            "faceLabels",
            pMesh.time().findInstance
            (
                localDir, "faceLabels", IOobjectOption::READ_IF_PRESENT
            ),
            localDir,
            pMesh.thisDb(),
            IOobjectOption::READ_IF_PRESENT,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        );

        if (io.typeHeaderOk<labelIOList>(false))
        {
            meshFaces_ = labelIOList(io);
        }
    }

    // pointLabels: persisted by writeObject() since the reduced-order rework.
    // READ_IF_PRESENT so meshes written before this was added still load.
    {
        IOobject io
        (
            "pointLabels",
            pMesh.time().findInstance
            (
                localDir, "pointLabels", IOobjectOption::READ_IF_PRESENT
            ),
            localDir,
            pMesh.thisDb(),
            IOobjectOption::READ_IF_PRESENT,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        );

        if (io.typeHeaderOk<labelIOList>(false))
        {
            meshPointLabels_ = labelIOList(io);
        }
    }
}


Foam::autoPtr<Foam::fvMesh> Foam::faMesh::buildRegionMesh
(
    const word& areaName,
    const polyMesh& pMesh,
    const labelList& faceLabels,
    labelList& pointLabels
)
{
    const fvMesh& fvm = refCast<const fvMesh>(pMesh);

    // Indirect patch over the selected parent boundary faces
    uindirectPrimitivePatch pp
    (
        UIndirectList<face>(pMesh.faces(), faceLabels),
        pMesh.points()
    );

    const labelList addressing(faceLabels);
    const label defaultPatchID = -1;

    const EdgeMap<label> edgeToPatch
    (
        edgePatchFn(pMesh, defaultPatchID, pp, addressing)
    );

    // Constructing a reduced-order (2D) mesh
    polyTopoChange meshMod(pMesh.boundaryMesh().size(), true, 2);

    labelList pointMap;
    patchToMeshFn(pMesh, pp, addressing, edgeToPatch, meshMod, pointMap);

    autoPtr<fvMesh> regionPtr;
    autoPtr<mapPolyMesh> map = meshMod.makeMesh
    (
        regionPtr,
        IOobject
        (
            (areaName.empty() ? word(polyMesh::defaultRegion) : areaName),
            fvm.pointsInstance(),
            pMesh.thisDb(),
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        ),
        fvm
    );

    // new-point -> parent (old) point: the coupling point addressing
    pointLabels = map().pointMap();

    return regionPtr;
}


Foam::faMesh::meshData Foam::faMesh::buildMeshData
(
    const word& areaName,
    const polyMesh& pMesh
)
{
    labelList faceLabels(readFaceLabels(areaName, pMesh));
    labelList pointLabels;
    autoPtr<fvMesh> region
    (
        buildRegionMesh(areaName, pMesh, faceLabels, pointLabels)
    );

    meshData data;
    data.region = std::move(region);
    data.faceLabels = std::move(faceLabels);
    data.pointLabels = std::move(pointLabels);
    return data;
}


Foam::labelList Foam::faMesh::selectPatchFaces
(
    const polyBoundaryMesh& pbm,
    const wordRes& patchNames
)
{
    const labelList patchIDs
    (
        pbm.indices(patchNames, true)  // useGroups
    );

    if (patchIDs.empty())
    {
        FatalErrorInFunction
            << "No polyPatches selected for finite-area mesh from "
            << patchNames << nl
            << "Known patch names: " << flatOutput(pbm.names()) << nl
            << exit(FatalError);
    }

    label nFaces = 0;
    for (const label patchi : patchIDs)
    {
        nFaces += pbm[patchi].size();
    }

    labelList faceLabels(nFaces);
    nFaces = 0;
    for (const label patchi : patchIDs)
    {
        const polyPatch& pp = pbm[patchi];
        forAll(pp, i)
        {
            faceLabels[nFaces++] = pp.start() + i;
        }
    }

    return faceLabels;
}


Foam::autoPtr<Foam::fvMesh> Foam::faMesh::buildRegionMeshFromDefinition
(
    const word& areaName,
    const polyMesh& pMesh,
    const labelList& faceLabels,
    const dictionary& faMeshDefinition,
    wordList& patchNames,
    wordList& patchTypes,
    labelList& patchNeighbProcs,
    labelList& pointLabels
)
{
    const fvMesh& fvm = refCast<const fvMesh>(pMesh);
    const polyBoundaryMesh& pbm = pMesh.boundaryMesh();

    // Indirect patch over the selected parent boundary faces
    uindirectPrimitivePatch pp
    (
        UIndirectList<face>(pMesh.faces(), faceLabels),
        pMesh.points()
    );

    const labelList& addressing = faceLabels;

    // Neighbour parent polyPatch per boundary edge (open edges -> -1)
    const EdgeMap<label> edgeToNeigh
    (
        edgePatchFn(pMesh, -1, pp, addressing)
    );

    // Parallel coupling of the area region across processor boundaries.
    // edgeToProc: neighbour rank of the parent processor patch each edge lies
    // on. faEdgeCount: how many processors carry an area-face on the edge
    // (>= 2 => the edge is interior to the area region but cut by the domain
    // decomposition, so it must become a genuine processor patch rather than a
    // default/undefined boundary).
    const EdgeMap<label> edgeToProc(edgeToNeighbProcFn(pMesh));
    const label myProc = parentMyProcFn(pMesh);

    // In a live parallel run we can confirm an edge is interior to the area
    // region (an area-face on >= 2 ranks) by synchronising a per-edge count
    // across the parent processor patches. In a serial decomposePar there is no
    // live neighbour to communicate with, so we fall back to treating every
    // unmatched area boundary edge that lies on a parent processor patch as an
    // inter-domain edge (the decomposition of a full-surface area region always
    // continues the region across an internal cut).
    const bool haveEdgeCount = UPstream::parRun();

    EdgeMap<label> faEdgeCount(2*(pp.nEdges() - pp.nInternalEdges()));
    if (haveEdgeCount)
    {
        const edgeList& ppEdges = pp.edges();
        const labelList& ppMp = pp.meshPoints();
        for (label edgei = pp.nInternalEdges(); edgei < pp.nEdges(); ++edgei)
        {
            faEdgeCount.insert(edge(ppMp, ppEdges[edgei]), 1);
        }
        syncTools::syncEdgeMap
        (
            pMesh,
            faEdgeCount,
            plusEqOp<label>(),
            dummyTransform()
        );
    }

    // Owner parent polyPatch per area face (which polyMeshPatch it came from)
    labelList faceOwnerPatch(pp.size());
    forAll(faceLabels, facei)
    {
        faceOwnerPatch[facei] = pbm.whichPatch(faceLabels[facei]);
    }

    // Transcribe the "boundary" sub-dictionary into patch definitions
    const dictionary& bndDict = faMeshDefinition.subDict("boundary");

    DynamicList<faPatchDef> defs(bndDict.size() + 4);
    for (const entry& dEntry : bndDict)
    {
        if (!dEntry.isDict())
        {
            continue;
        }
        const dictionary& patchDict = dEntry.dict();

        auto& def = defs.emplace_back();
        def.name_ = dEntry.keyword();
        def.type_ = patchDict.get<word>("type");

        word patchName;
        if (patchDict.readIfPresent("ownerPolyPatch", patchName))
        {
            def.ownerPolyPatchId_ = pbm.findPatchID(patchName);
            if (def.ownerPolyPatchId_ < 0)
            {
                FatalErrorInFunction
                    << "ownerPolyPatch " << patchName << " not found" << nl
                    << exit(FatalError);
            }
        }
        if (patchDict.readIfPresent("neighbourPolyPatch", patchName))
        {
            def.neighPolyPatchId_ = pbm.findPatchID(patchName);
            if (def.neighPolyPatchId_ < 0)
            {
                FatalErrorInFunction
                    << "neighbourPolyPatch " << patchName << " not found" << nl
                    << exit(FatalError);
            }
        }
    }

    const label nUserDefs = defs.size();

    // Optional 'empty' placeholder patch (always emitted, typically size 0)
    label emptyDefIndex = -1;
    {
        word emptyName;
        if (faMeshDefinition.readIfPresent("emptyPatch", emptyName))
        {
            emptyDefIndex = defs.size();
            auto& def = defs.emplace_back();
            def.name_ = emptyName;
            def.type_ = "empty";
        }
    }

    // Catch-all for otherwise undefined edges (dropped later if unused)
    const label undefDefIndex = defs.size();
    {
        auto& def = defs.emplace_back();
        def.name_ = "undefined";
        def.type_ = "patch";
        if (const dictionary* dPtr = faMeshDefinition.findDict("defaultPatch"))
        {
            dPtr->readIfPresent("name", def.name_);
            dPtr->readIfPresent("type", def.type_);
        }
    }

    // Assign each boundary edge to a definition patch
    const edgeList& edges = pp.edges();
    const labelList& mp = pp.meshPoints();
    const labelListList& edgeFaces = pp.edgeFaces();

    EdgeMap<label> edgeToDef(2*(pp.nEdges() - pp.nInternalEdges()));
    bitSet usedDef(defs.size());

    // Lazily-created processor patch definitions, keyed by neighbour rank
    Map<label> procDefIndex;

    for (label edgei = pp.nInternalEdges(); edgei < pp.nEdges(); ++edgei)
    {
        const edge meshE(mp, edges[edgei]);
        const label neighP = edgeToNeigh[meshE];
        const label ownerP = faceOwnerPatch[edgeFaces[edgei][0]];

        label best = -1;
        for (label defi = 0; defi < nUserDefs; ++defi)
        {
            const int m = defs[defi].matchPatchPair(ownerP, neighP);
            if (m == 3)
            {
                best = defi;
                break;
            }
            else if (m == 2 && best < 0)
            {
                best = defi;
            }
        }

        if (best < 0)
        {
            // Inter-domain edge: interior to the area region (>= 2 procs carry
            // an area-face on it) but cut by the parent decomposition. Route it
            // onto a per-neighbour processor patch so the reduced-order region
            // couples across ranks exactly like any other fvMesh.
            const label nbrProc = edgeToProc.lookup(meshE, -1);
            const bool interDomain =
            (
                nbrProc >= 0
             && (!haveEdgeCount || faEdgeCount.lookup(meshE, 0) >= 2)
            );
            if (interDomain)
            {
                auto iter = procDefIndex.find(nbrProc);
                if (iter.good())
                {
                    best = iter.val();
                }
                else
                {
                    best = defs.size();
                    procDefIndex.insert(nbrProc, best);
                    auto& def = defs.emplace_back();
                    def.type_ = processorPolyPatch::typeName;
                    def.procNeighbour_ = nbrProc;
                    def.name_ = processorPolyPatch::newName(myProc, nbrProc);
                }
            }
            else
            {
                best = (neighP < 0 && emptyDefIndex >= 0)
                     ? emptyDefIndex
                     : undefDefIndex;
            }
        }

        edgeToDef.insert(meshE, best);
        usedDef.set(best);
    }

    // Finalise the emitted patch list: all user defs + empty placeholder are
    // always emitted; the undefined catch-all only when actually used. Build
    // the (possibly sparse) old->new remap so the emitted patches are dense.
    labelList defRemap(defs.size(), -1);
    label nEmit = 0;

    // Non-processor patches first, in definition order
    for (label defi = 0; defi < defs.size(); ++defi)
    {
        if (defs[defi].procNeighbour_ >= 0)
        {
            continue;  // processor patches emitted last
        }

        const bool keep =
        (
            defi < nUserDefs
         || defi == emptyDefIndex
         || (defi == undefDefIndex && usedDef.test(defi))
        );
        if (keep)
        {
            defRemap[defi] = nEmit++;
        }
    }

    // Processor patches last, ordered by neighbour rank so the emitted order is
    // deterministic and independent of edge-traversal order
    {
        labelList procDefs(procDefIndex.size());
        label n = 0;
        forAllConstIters(procDefIndex, it)
        {
            procDefs[n++] = it.val();
        }
        Foam::sort
        (
            procDefs,
            [&](const label a, const label b) noexcept
            {
                return defs[a].procNeighbour_ < defs[b].procNeighbour_;
            }
        );
        for (const label defi : procDefs)
        {
            defRemap[defi] = nEmit++;
        }
    }

    if (usedDef.test(undefDefIndex) && defRemap[undefDefIndex] >= 0)
    {
        WarningInFunction
            << "finite-area mesh '" << areaName
            << "' has boundary edges not matched by any patch definition;"
            << " added to '" << defs[undefDefIndex].name_ << "'" << nl;
    }

    // Remap edge assignments to the dense emitted indices
    forAllIters(edgeToDef, iter)
    {
        iter.val() = defRemap[iter.val()];
    }

    // Record the emitted patch names/types/neighbour-procs (dense, by region
    // patch index)
    patchNames.resize(nEmit);
    patchTypes.resize(nEmit);
    patchNeighbProcs.resize(nEmit, -1);
    for (label defi = 0; defi < defs.size(); ++defi)
    {
        const label emit = defRemap[defi];
        if (emit >= 0)
        {
            patchNames[emit] = defs[defi].name_;
            patchTypes[emit] = defs[defi].type_;
            patchNeighbProcs[emit] = defs[defi].procNeighbour_;
        }
    }

    // Construct the reduced-order (2D) mesh, grouping boundary faces by the
    // emitted definition-patch index
    polyTopoChange meshMod(nEmit, true, 2);

    labelList pointMap;
    patchToMeshFn(pMesh, pp, addressing, edgeToDef, meshMod, pointMap);

    // Template patch clone source for makeMesh (names/types are overridden by
    // the faMesh constructor afterwards). Map each reduced-order patch onto a
    // parent template of the SAME KIND: a processor patch to rank k onto the
    // parent processorPolyPatch to rank k (guaranteed to exist -- the
    // inter-domain edges were derived from exactly that parent patch), a
    // non-processor patch onto a distinct parent non-processor patch (avoids
    // transient duplicate-name warnings). The old positional map instead cloned
    // arbitrary parent patches, so on a CLOSED manifold -- where the region has
    // no non-processor boundary at all -- the parent's leading non-processor
    // patch shifted every processor patch onto a template with the wrong (or
    // no) neighbour. That mismatch fed reorderCoupledFaces an asymmetric comm
    // graph and segfaulted (see the syncParallel note below).
    Map<label> parentProcPatch(2*pbm.size());
    DynamicList<label> parentNonProc(pbm.size());
    forAll(pbm, patchi)
    {
        const auto* ppp = isA<processorPolyPatch>(pbm[patchi]);
        if (ppp)
        {
            parentProcPatch.insert(ppp->neighbProcNo(), patchi);
        }
        else
        {
            parentNonProc.append(patchi);
        }
    }
    const label fallbackTemplate =
        (parentNonProc.empty() ? 0 : parentNonProc.front());

    labelList patchMap(nEmit);
    label nonProcSeen = 0;
    forAll(patchMap, i)
    {
        const label nbrProc = patchNeighbProcs[i];
        if (nbrProc >= 0)
        {
            patchMap[i] = parentProcPatch[nbrProc];
        }
        else
        {
            patchMap[i] =
            (
                nonProcSeen < parentNonProc.size()
              ? parentNonProc[nonProcSeen]
              : fallbackTemplate
            );
            ++nonProcSeen;
        }
    }

    // makeMesh with syncParallel=false: do NOT run the inherited cross-processor
    // reorderCoupledFaces on this reduced-order mesh. That routine matches
    // coupled faces via processorPolyPatch::order, which is built for 3D
    // polygonal faces -- it compares face centres/normals/anchor points. Here a
    // "processor face" is a single reduced-order boundary EDGE (a 2-point face
    // of zero area and undefined normal), so that geometric matcher is
    // degenerate: it rotates the 2-point faces, reversing an edge relative to
    // its owner area-cell and breaking the cell's edge loop (basicFvGeometry
    // scheme's primitiveMesh::loop) as well as the parallel face correspondence.
    //
    // The coupled-face correspondence the reorder would provide is instead
    // supplied by construction: patchToMeshFn emits every boundary face sorted
    // by its parent edge-endpoint coordinates (a global, decomposition-invariant
    // key), so the two sides of a processor patch already list the same shared
    // edges in the same order with opposite (outward-from-own-cell) orientation
    // -- exactly the processor convention. Disabling the reorder therefore keeps
    // that valid ordering rather than corrupting it. syncParallel gates ONLY
    // reorderCoupledFaces inside compactAndReorder (no other collective), so this
    // is also collective-free; the boundary is finalised uniformly on every rank
    // later by finaliseParallelBoundary().
    autoPtr<fvMesh> regionPtr;
    autoPtr<mapPolyMesh> map = meshMod.makeMesh
    (
        regionPtr,
        IOobject
        (
            (areaName.empty() ? word(polyMesh::defaultRegion) : areaName),
            fvm.pointsInstance(),
            pMesh.thisDb(),
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        ),
        fvm,
        patchMap,
        /*syncParallel=*/false
    );

    // new-point -> parent (old) point: the coupling point addressing
    pointLabels = map().pointMap();

    return regionPtr;
}


Foam::faMesh::meshData Foam::faMesh::buildMeshDataFromDefinition
(
    const word& areaName,
    const polyMesh& pMesh,
    const dictionary& faMeshDefinition
)
{
    labelList faceLabels
    (
        selectPatchFaces
        (
            pMesh.boundaryMesh(),
            faMeshDefinition.get<wordRes>("polyMeshPatches")
        )
    );

    // faceLabels are gathered patch-by-patch; keep them sorted so the
    // reduced-order cell ordering is deterministic
    Foam::sort(faceLabels);

    meshData data;
    data.region = buildRegionMeshFromDefinition
    (
        areaName,
        pMesh,
        faceLabels,
        faMeshDefinition,
        data.patchNames,
        data.patchTypes,
        data.patchNeighbProcs,
        data.pointLabels
    );
    data.faceLabels = std::move(faceLabels);
    return data;
}


Foam::IOobject Foam::faMesh::regionIO
(
    const word& areaName,
    const polyMesh& pMesh
)
{
    return IOobject
    (
        (areaName.empty() ? word(polyMesh::defaultRegion) : areaName),
        pMesh.time().timeName(),
        faMeshesRegistry::Registry(pMesh),
        IOobjectOption::NO_READ,
        IOobjectOption::AUTO_WRITE,
        IOobjectOption::REGISTER
    );
}


Foam::dictionary Foam::faMesh::readLegacyControlDict
(
    const word& faDictName,
    const objectRegistry& db
)
{
    // Legacy finiteArea cases ship these controls in system/finite-area/
    // (== faMesh::prefix()), at case level. Resolve relative to the case root
    // (Time) so the path is correct regardless of whether the passed registry
    // is the parent mesh (createFaMesh path) or the faMeshesRegistry (which
    // already carries the finite-area dbDir). Read as READ_IF_PRESENT so a
    // missing file simply yields an empty dictionary (no fallback available).
    IOobject controlIO
    (
        faDictName,
        db.time().system(),
        faMesh::prefix(),
        db.time(),
        IOobjectOption::READ_IF_PRESENT,
        IOobjectOption::NO_WRITE,
        IOobjectOption::NO_REGISTER
    );

    if (controlIO.typeHeaderOk<IOdictionary>(true))
    {
        return dictionary(IOdictionary(controlIO));
    }

    return dictionary();
}


void Foam::faMesh::calcWhichPatchFaces() const
{
    const polyBoundaryMesh& pbm = parentMesh()().boundaryMesh();

    polyPatchFacesPtr_ =
        std::make_unique<List<labelPair>>(pbm.whichPatchFace(meshFaces_));

    labelHashSet ids;
    for (const labelPair& tup : *polyPatchFacesPtr_)
    {
        ids.insert(tup.first());
    }
    ids.erase(-1);  // Without internal faces (patchi == -1)

    // parSync
    Foam::reduce
    (
        ids,
        bitOrOp<labelHashSet>(),
        UPstream::msgType(),
        parentMesh()().comm()
    );

    polyPatchIdsPtr_ = std::make_unique<labelList>(ids.sortedToc());
}


void Foam::faMesh::finaliseParallelBoundary()
{
    // Parallel-safe boundary finalisation for the reduced-order region.
    //
    // The region is an fvMesh whose geometry is computed by faGeometryScheme.
    // Evaluating that geometry (faceCentres()) pulls in two worldComm
    // collectives: globalData() (via the lduAddr() comms schedule) and the
    // point-normal syncTools swap in calcPointAreaNormals(). Additionally
    // nTopologicalD() (checked inside the scheme) is itself a worldComm reduce.
    //
    // With the stock validBoundary=true path these fire lazily from
    // processorPolyPatch::initUpdateMesh()/initGeometry() during
    // polyBoundaryMesh::updateMesh()/calcGeometry() -- i.e. ONLY on ranks whose
    // reduced-order region carries a processor patch. When the parent
    // decomposition leaves the area region absent from a subdomain, that rank
    // has no processor patch, skips those collectives, and races ahead to
    // checkDefinition()'s all-reduce: the collectives mismatch across ranks and
    // the run deadlocks.
    //
    // Fix: drive the finalisation explicitly and in the same order on every
    // rank, forcing the geometry-scheme collectives (faceCentres()) up front so
    // they fire uniformly -- the empty subdomain contributes empty data but
    // still participates. The subsequent processor-patch topology/geometry
    // exchange and checkDefinition() are then reached in lockstep.
    polyBoundaryMesh& pbm =
        const_cast<polyBoundaryMesh&>(this->boundaryMesh());

    // 0. Seed the topological dimension uniformly. nTopologicalD() performs a
    //    worldComm reduce on first evaluation; the local geometry calculators
    //    no longer trigger it, so force it here (all ranks, including an empty
    //    subdomain) before any lazy/asymmetric path can.
    (void)this->nTopologicalD();

    // 1. Geometry-scheme worldComm collectives, uniform on all ranks. Caches
    //    faceCentres/faceAreas so the later processor-patch init (which reads
    //    them) no longer triggers these collectives lazily/asymmetrically.
    (void)this->faceCentres();

    // 2. Processor-patch topology (neighbour face ordering) -- collective.
    pbm.updateMesh();

    // 3. Patch geometry incl. the processor-patch neighbour face-centre
    //    exchange -- collective; the primitive geometry is now cached, so no
    //    lazy geometry-scheme re-trigger here. movePoints() is the public
    //    equivalent of the (private) polyBoundaryMesh::calcGeometry().
    pbm.movePoints(this->points());

    // 4. Boundary sanity checks -- collective all-reduce, now in lockstep.
    pbm.checkDefinition();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMesh::faMesh
(
    const word& areaName,
    const polyMesh& pMesh,
    faMesh::meshData&& data
)
:
    // Build the reduced-order fvMesh from the region components, honouring
    // legacy finiteArea faSchemes/faSolution as the fvSchemes/fvSolution
    // fallback (the region is built from its parent and typically has no own
    // system/<region>/fvSchemes; mirrors the IOobject constructor's D2
    // backward-compat behaviour so createFaMesh- and solver-built regions
    // resolve their differencing/solver settings identically).
    fvMesh
    (
        regionIO(areaName, pMesh),
        pointField(data.region->points()),
        faceList(data.region->faces()),
        labelList(data.region->faceOwner()),
        labelList(data.region->faceNeighbour()),
        readLegacyControlDict("faSchemes", pMesh.thisDb()),
        readLegacyControlDict("faSolution", pMesh.thisDb())
    ),
    areaName_(areaName),
    meshFaces_(std::move(data.faceLabels)),
    meshPointLabels_(std::move(data.pointLabels))
{
    const polyBoundaryMesh& srcbm = data.region->boundaryMesh();

    if (data.patchNames.empty())
    {
        // Clone the boundary patches from the temporary region mesh
        // (parent-derived patch names/types)
        polyPatchList newPatches(srcbm.size());
        forAll(srcbm, patchi)
        {
            newPatches.set(patchi, srcbm[patchi].clone(this->boundaryMesh()));
        }

        // Insert the patches WITHOUT the collective boundary finalisation
        // (validBoundary=false); it is driven uniformly on every rank by
        // finaliseParallelBoundary() below.
        addFvPatches(newPatches, /*validBoundary=*/false);
    }
    else
    {
        // Build named/typed boundary patches from the faMeshDefinition.
        // The region boundary faces are already grouped contiguously by
        // definition-patch index (see buildRegionMeshFromDefinition), and
        // the def list has been trimmed to exactly the emitted patches, so
        // each region patch maps 1:1 onto a named patch.
        const label myProc = parentMyProcFn(pMesh);

        polyPatchList newPatches(srcbm.size());
        forAll(srcbm, patchi)
        {
            const label nbrProc =
            (
                patchi < data.patchNeighbProcs.size()
              ? data.patchNeighbProcs[patchi]
              : -1
            );

            if (nbrProc >= 0)
            {
                // Inter-domain area edges: a genuine fvMesh processor patch.
                // Because the reduced-order region is an fvMesh, the standard
                // processorFvPatchField coupling handles the halo exchange with
                // no finite-area-specific machinery.
                newPatches.set
                (
                    patchi,
                    new processorPolyPatch
                    (
                        srcbm[patchi].size(),
                        srcbm[patchi].start(),
                        patchi,
                        this->boundaryMesh(),
                        myProc,
                        nbrProc
                    )
                );
            }
            else
            {
                newPatches.set
                (
                    patchi,
                    polyPatch::New
                    (
                        data.patchTypes[patchi],
                        data.patchNames[patchi],
                        srcbm[patchi].size(),
                        srcbm[patchi].start(),
                        patchi,
                        this->boundaryMesh()
                    )
                );
            }
        }

        addFvPatches(newPatches, /*validBoundary=*/false);
    }

    // Publish the parent-point map BEFORE the boundary finalisation triggers
    // the geometry build: finaliseParallelBoundary() -> faceCentres() runs the
    // reduced-order geometry scheme's calcPointAreaNormals, which consults this
    // map to sum point-area-normals through the parent's global point coupling.
    registerParentPointMap();

    // Finalise the boundary topology/geometry collectively and uniformly on
    // every rank (parallel-safe; avoids the empty-subdomain deadlock).
    finaliseParallelBoundary();

    // Release the temporary region mesh
    data.region.reset(nullptr);
}


Foam::faMesh::faMesh
(
    const word& areaName,
    const polyMesh& pMesh,
    const bool doInit
)
:
    faMesh(areaName, pMesh, buildMeshData(areaName, pMesh))
{}


Foam::faMesh::faMesh
(
    const word& areaName,
    const polyMesh& pMesh,
    const dictionary& faMeshDefinition,
    const bool doInit
)
:
    faMesh
    (
        areaName,
        pMesh,
        buildMeshDataFromDefinition(areaName, pMesh, faMeshDefinition)
    )
{}


Foam::faMesh::faMesh(const IOobject& io, const bool doInit)
:
    // D2 backward-compat: read the region fvSchemes/fvSolution if present,
    // otherwise fall back to legacy system/finite-area/faSchemes|faSolution.
    fvMesh
    (
        io,
        doInit,
        readLegacyControlDict("faSchemes", io.db()),
        readLegacyControlDict("faSolution", io.db())
    ),
    areaName_(io.name())
{
    // The fvMesh base only read the reduced-order polyMesh; repopulate the
    // parent coupling addressing (meshFaces_/meshPointLabels_) from disk so
    // faceLabels()/pointLabels()/patch() work for createFaMesh-style reads.
    readCouplingAddressing();

    // Publish the parent-point map so the reduced-order geometry scheme sums
    // point-area-normals through the parent's global point coupling (see
    // registerParentPointMap). The geometry build is demand-driven here, so
    // registering now (before first geometry access) is sufficient.
    registerParentPointMap();

    // Notify (once) when the legacy faSchemes fallback is actually in force,
    // i.e. no region fvSchemes exists but a legacy faSchemes does.
    if (io.isReadRequired() || io.isReadOptional())
    {
        const bool haveFvSchemes = IOobject
        (
            "fvSchemes",
            time().system(),
            *this,
            IOobjectOption::NO_READ
        ).typeHeaderOk<IOdictionary>(false);

        const bool haveFaSchemes = IOobject
        (
            "faSchemes",
            time().system(),
            faMesh::prefix(),
            time(),
            IOobjectOption::NO_READ
        ).typeHeaderOk<IOdictionary>(false);

        if (!haveFvSchemes && haveFaSchemes)
        {
            Info<< "faMesh: region '" << areaName_
                << "' has no system/" << areaName_
                << "/fvSchemes; using legacy system/" << faMesh::prefix()
                << "/faSchemes (deprecated - please migrate to fvSchemes)"
                << endl;
        }
    }
}


// * * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * //

const Foam::objectRegistry* Foam::faMesh::registry(const polyMesh& pMesh)
{
    return pMesh.cfindObject<objectRegistry>(faMesh::prefix());
}


const Foam::objectRegistry& Foam::faMesh::Registry(const polyMesh& pMesh)
{
    return faMeshesRegistry::Registry(pMesh);
}


const Foam::faMesh& Foam::faMesh::mesh
(
    const polyMesh& pMesh,
    const word& areaRegion
)
{
    const objectRegistry* obr = faMesh::registry(pMesh);

    if (!obr)
    {
        obr = &(pMesh.thisDb());
    }

    if (areaRegion.empty())
    {
        return obr->lookupObject<faMesh>(polyMesh::defaultRegion);
    }

    return obr->lookupObject<faMesh>(areaRegion);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::List<Foam::labelPair>& Foam::faMesh::whichPatchFaces() const
{
    if (!polyPatchFacesPtr_ || !polyPatchIdsPtr_)
    {
        calcWhichPatchFaces();
    }
    return *polyPatchFacesPtr_;
}


const Foam::labelList& Foam::faMesh::whichPolyPatches() const
{
    if (!polyPatchFacesPtr_ || !polyPatchIdsPtr_)
    {
        calcWhichPatchFaces();
    }
    return *polyPatchIdsPtr_;
}


const Foam::uindirectPrimitivePatch& Foam::faMesh::patch() const
{
    if (!patchPtr_)
    {
        const refPtr<fvMesh> pRef(parentMesh());
        const polyMesh& pMesh = pRef();

        patchPtr_ = std::make_unique<uindirectPrimitivePatch>
        (
            UIndirectList<face>(pMesh.faces(), meshFaces_),
            pMesh.points()
        );
    }
    return *patchPtr_;
}


void Foam::faMesh::calcNgbPolyPatchFaces() const
{
    const label nBndFaces = max(nFaces() - nInternalFaces(), label(0));

    ngbPolyPatchFacesPtr_ = std::make_unique<labelList>(nBndFaces, label(-1));
    auto& ngbFace = *ngbPolyPatchFacesPtr_;

    // Standalone area region (no parent), or parent point addressing
    // unavailable/inconsistent: no coupling to derive.
    if (!hasParent() || meshPointLabels_.size() != nPoints())
    {
        return;
    }

    const refPtr<fvMesh> pRef(parentMesh());
    const polyMesh& pMesh = pRef();
    const polyBoundaryMesh& pbm = pMesh.boundaryMesh();
    const faceList& pFaces = pMesh.faces();

    // Parent faces that make up this area region: excluded as neighbours
    bitSet isZoneFace(pMesh.nFaces());
    isZoneFace.set(meshFaces_);

    // The reduced-order boundary faces ARE the finite-area boundary edges
    // (2-point faces in reduced-point numbering). Map each to its parent edge.
    const label nIntF = nInternalFaces();
    const faceList& rFaces = faces();

    EdgeMap<label> edgeToBndFace(2*nBndFaces);
    for (label facei = nIntF; facei < nFaces(); ++facei)
    {
        const face& f = rFaces[facei];
        const edge pe(meshPointLabels_[f[0]], meshPointLabels_[f[1]]);
        edgeToBndFace.insert(pe, facei - nIntF);
    }

    // For each parent non-processor, non-region boundary face, tag the rim
    // edges it shares with the region. The neighbour is local, so no halo
    // exchange is required (maxEqOp keeps corner cases deterministic).
    for (label patchi = 0; patchi < pbm.nNonProcessor(); ++patchi)
    {
        const polyPatch& pp = pbm[patchi];
        forAll(pp, i)
        {
            const label facei = pp.start() + i;
            if (isZoneFace[facei])
            {
                continue;
            }

            const face& f = pFaces[facei];
            forAll(f, fp)
            {
                const auto iter = edgeToBndFace.cfind(f.edge(fp));
                if (iter.good())
                {
                    maxEqOp<label>()(ngbFace[iter.val()], facei);
                }
            }
        }
    }
}


const Foam::labelList& Foam::faMesh::ngbPolyPatchFaces() const
{
    if (!ngbPolyPatchFacesPtr_)
    {
        calcNgbPolyPatchFaces();
    }
    return *ngbPolyPatchFacesPtr_;
}


void Foam::faMesh::calcPatchPointMap() const
{
    const uindirectPrimitivePatch& pp = patch();

    patchPointMapPtr_ = std::make_unique<labelList>(nPoints(), label(-1));
    auto& map = *patchPointMapPtr_;

    // Standalone area region (no parent) or parent point addressing
    // unavailable/inconsistent: the reduced-order and patch numberings then
    // coincide by construction, so fall back to identity.
    if (!hasParent() || meshPointLabels_.size() != nPoints())
    {
        forAll(map, i)
        {
            map[i] = i;
        }
        return;
    }

    // The reduced-order mesh and patch() share the same free-surface points,
    // both keyed on parent-mesh point labels, but the two enumerations may
    // differ. Route reduced-order -> parent -> patch-local.
    const Map<label>& parentToLocal = pp.meshPointMap();
    forAll(map, rp)
    {
        const auto iter = parentToLocal.cfind(meshPointLabels_[rp]);
        map[rp] = iter.good() ? iter.val() : -1;
    }
}


const Foam::labelList& Foam::faMesh::patchPointMap() const
{
    if (!patchPointMapPtr_)
    {
        calcPatchPointMap();
    }
    return *patchPointMapPtr_;
}


void Foam::faMesh::calcPatchEdgeToBoundaryFace() const
{
    const uindirectPrimitivePatch& pp = patch();

    patchEdgeToBndFacePtr_ =
        std::make_unique<labelList>(pp.nEdges(), label(-1));
    auto& map = *patchEdgeToBndFacePtr_;

    // Each reduced-order boundary face is a 2-point face whose points are the
    // reduced-order point labels of a finite-area boundary edge. Translate them
    // to patch-local numbering and locate the coincident parent-patch edge.
    const labelList& ppMap = patchPointMap();
    const faceList& rFaces = faces();
    for (label facei = nInternalFaces(); facei < nFaces(); ++facei)
    {
        const face& f = rFaces[facei];
        const label la = ppMap[f[0]];
        const label lb = ppMap[f[1]];
        if (la == -1 || lb == -1)
        {
            continue;
        }

        const label pe = pp.findEdge(edge(la, lb));
        if (pe != -1)
        {
            map[pe] = facei;
        }
    }
}


const Foam::labelList& Foam::faMesh::patchEdgeToBoundaryFace() const
{
    if (!patchEdgeToBndFacePtr_)
    {
        calcPatchEdgeToBoundaryFace();
    }
    return *patchEdgeToBndFacePtr_;
}


Foam::label Foam::faMesh::ngbPolyPatchIndex(const label faPatchi) const
{
    const refPtr<fvMesh> pRef(parentMesh());
    if (!pRef.good())
    {
        return -1;
    }

    const labelList& ngbFace = ngbPolyPatchFaces();
    const polyBoundaryMesh& pbm = pRef().boundaryMesh();

    // Use the polyPatch view: a constraint fvPatch (e.g. emptyFvPatch) reports
    // size() == 0 even though its reduced-order boundary faces are real, so the
    // fvPatch face count/addressing would hide the rim edges entirely.
    const polyPatch& pp = boundary()[faPatchi].patch();
    const label base = pp.start() - nInternalFaces();

    label ngbPatch = -1;
    for (label i = 0; i < pp.size(); ++i)
    {
        const label pf = ngbFace[base + i];
        if (pf != -1)
        {
            maxEqOp<label>()(ngbPatch, pbm.whichPatch(pf));
        }
    }

    return ngbPatch;
}


Foam::tmp<Foam::vectorField>
Foam::faMesh::ngbPolyPatchFaceNormals(const label faPatchi) const
{
    // polyPatch view: emptyFvPatch (and other constraint fvPatches) report
    // size() == 0, so one normal per real reduced-order boundary face is only
    // available through the underlying polyPatch.
    const polyPatch& pp = boundary()[faPatchi].patch();

    auto tfn = tmp<vectorField>::New(pp.size(), Zero);
    auto& fn = tfn.ref();

    const refPtr<fvMesh> pRef(parentMesh());
    if (!pRef.good())
    {
        return tfn;
    }

    const labelList& ngbFace = ngbPolyPatchFaces();
    const polyMesh& pMesh = pRef();
    const faceList& pFaces = pMesh.faces();
    const pointField& pPts = pMesh.points();
    const label base = pp.start() - nInternalFaces();

    forAll(fn, i)
    {
        const label pf = ngbFace[base + i];
        if (pf != -1)
        {
            fn[i] = pFaces[pf].unitNormal(pPts);
        }
    }

    return tfn;
}


Foam::tmp<Foam::vectorField>
Foam::faMesh::ngbPolyPatchPointNormals(const label faPatchi) const
{
    // Neighbour parent-face normal per boundary edge (patch face)
    const vectorField fn(ngbPolyPatchFaceNormals(faPatchi));

    // The reduced-order boundary patch faces ARE the finite-area edges, so
    // pointFaces() gives the rim edges incident to each patch point.
    const polyPatch& pp = boundary()[faPatchi].patch();
    const labelListList& pointFaces = pp.pointFaces();

    auto tpn = tmp<vectorField>::New(pp.nPoints(), Zero);
    auto& pn = tpn.ref();

    forAll(pn, pointi)
    {
        vector& n = pn[pointi];
        for (const label facei : pointFaces[pointi])
        {
            n += fn[facei];
        }
        n.normalise();
    }

    return tpn;
}


Foam::tmp<Foam::scalarField> Foam::faMesh::lPN() const
{
    // Straight-line P-N distance per internal finite-area edge
    // (reduced-order internal face centre-to-centre distance).
    const tmp<surfaceVectorField> tdelta(delta());
    return tmp<scalarField>::New(mag(tdelta().primitiveField()));
}


Foam::labelList Foam::faMesh::boundaryPoints() const
{
    labelList bp(patch().boundaryPoints());
    Foam::sort(bp);
    return bp;
}


Foam::labelList Foam::faMesh::internalPoints() const
{
    bitSet markPoints(patch().nPoints());
    markPoints.set(patch().boundaryPoints());
    markPoints.flip();
    return markPoints.sortedToc();
}


void Foam::faMesh::updatePointNormalsBoundaryCorrection() const
{
    // Consumed by basicFvGeometryScheme::calcPointAreaNormals (finiteVolume):
    // publish, per reduced-order point, the accumulated neighbour (wall)
    // polyPatch face normal whose component must be removed from the point-area
    // normal at a contact line, so the reduced-order geometry scheme enforces
    // the finite-area/neighbour-polyPatch perpendicularity that curvature (the
    // surface-tension restoring term) depends on. Mirrors the legacy
    // faMesh::calcPointAreaNormals haloFaceNormals removeCollinear correction,
    // restricted to patches opted in via correctPatchPointNormals()
    // (interfaceTrackingFvMesh's pointNormalsCorrectionPatches).
    //
    // Computed from the (already moved) parent points, keyed on the stable
    // reduced-order point labels, so it is valid to build before the base
    // fvMesh::movePoints rebuilds the reduced-order geometry.
    if (!correctPatchPointNormalsPtr_ || !hasParent())
    {
        return;
    }

    vectorField corr(nPoints(), Zero);
    const faceList& rFaces = faces();

    forAll(boundary(), patchi)
    {
        const fvPatch& fap = boundary()[patchi];

        // Legacy applies the correction only on opted-in, non-coupled patches.
        if (!correctPatchPointNormals(patchi) || fap.coupled())
        {
            continue;
        }

        // polyPatch view: a constraint fvPatch (e.g. emptyFvPatch) reports
        // size() == 0, so the neighbour normals (one per reduced-order boundary
        // face = finite-area edge) are only available through the polyPatch.
        const polyPatch& pp = fap.patch();
        const tmp<vectorField> tfn(ngbPolyPatchFaceNormals(patchi));
        const vectorField& fn = tfn();

        // The reduced-order boundary faces of this patch span
        // [pp.start(), pp.start()+pp.size()) in reduced-order face numbering;
        // fn is indexed in the same patch-face order. Each 2-point face is a
        // finite-area boundary edge, so accumulate its neighbour normal onto
        // both reduced-order end points.
        for (label i = 0; i < pp.size(); ++i)
        {
            const face& f = rFaces[pp.start() + i];
            corr[f[0]] += fn[i];
            corr[f[1]] += fn[i];
        }
    }

    // (Re)publish in this mesh's own registry. Name kept in sync with
    // basicFvGeometryScheme::calcPointAreaNormals.
    auto* existing =
        getObjectPtr<vectorIOField>("pointAreaNormalsBoundaryCorrection");

    if (existing && existing->size() == corr.size())
    {
        static_cast<vectorField&>(*existing) = corr;
    }
    else
    {
        regIOobject::store
        (
            new vectorIOField
            (
                IOobject
                (
                    "pointAreaNormalsBoundaryCorrection",
                    pointsInstance(),
                    *this,
                    IOobjectOption::NO_READ,
                    IOobjectOption::NO_WRITE,
                    IOobjectOption::REGISTER
                ),
                std::move(corr)
            )
        );
    }
}


bool Foam::faMesh::movePoints()
{
    // Follow the parent-boundary motion: the reduced-order point i coincides
    // with parent-mesh point meshPointLabels_[i], so re-read those positions
    // from the (already moved) parent and hand them to fvMesh::movePoints,
    // which clears/rebuilds the reduced-order geometry. Without a usable
    // parent point map (standalone region, or a mesh read before pointLabels
    // was persisted) there is nothing to track.
    if (!hasParent() || meshPointLabels_.size() != nPoints())
    {
        return false;
    }

    // Refresh the contact-line point-normal correction from the (already moved)
    // parent BEFORE the geometry rebuild, so the reduced-order pointAreaNormals
    // scheme applies it during fvMesh::movePoints and it propagates through the
    // eagerly-derived edge normals / Le into faceCurvatures.
    updatePointNormalsBoundaryCorrection();

    const pointField& parentPts = parentMesh()().points();

    pointField newPoints(nPoints());
    forAll(newPoints, pointi)
    {
        newPoints[pointi] = parentPts[meshPointLabels_[pointi]];
    }

    fvMesh::movePoints(newPoints);

    return true;
}


Foam::boolList& Foam::faMesh::correctPatchPointNormals() const
{
    if (!correctPatchPointNormalsPtr_)
    {
        correctPatchPointNormalsPtr_ =
            std::make_unique<boolList>(boundary().size(), false);
    }
    return *correctPatchPointNormalsPtr_;
}


bool Foam::faMesh::correctPatchPointNormals(const label patchID) const
{
    if (correctPatchPointNormalsPtr_)
    {
        return (*correctPatchPointNormalsPtr_)[patchID];
    }
    return false;
}


void Foam::faMesh::registerParentPointMap()
{
    // Nothing to publish for a standalone region, when the parent point
    // addressing is unavailable (a mesh written before pointLabels was
    // persisted), or when it is inconsistent with the reduced-order mesh.
    if
    (
        meshPointLabels_.empty()
     || meshPointLabels_.size() != nPoints()
     || !hasParent()
    )
    {
        return;
    }

    // Publish, in this mesh's own registry, the parent-mesh point label of each
    // reduced-order point. The reduced-order (2D) geometry scheme
    // (basicFvGeometryScheme::calcPointAreaNormals, finiteVolume) looks this up
    // by name and, together with parentMesh(), sums the locally-accumulated
    // point-area-normals through the PARENT mesh's correct global point
    // coupling. The parent's globalPointSlaves resolve BOTH adjacent-rank and
    // multiply-shared "global" points (a sphere pole at a 2x2x2 domain corner
    // is shared by four non-adjacent ranks) that the reduced-order mesh's own
    // collapsed coupling cannot -- the parent-global analogue of the legacy
    // faMesh::calcPointAreaNormals processor exchange + sharedPoints
    // combineReduce.
    //
    // A registry entry (rather than a virtual or data member on fvMesh) is used
    // deliberately: it keeps the fvMesh ABI untouched (fvMesh has many
    // subclasses compiled in other libraries) while still letting the generic
    // finiteVolume scheme reach faMesh-specific addressing. The parent mesh's
    // globalData is independent of this reduced-order mesh, so the sum is safe
    // to evaluate from the construction-time geometry build
    // (finaliseParallelBoundary -> faceCentres()), which runs uniformly on
    // every rank; no lifecycle deferral is required.
    //
    // Name kept in sync with basicFvGeometryScheme::calcPointAreaNormals.
    if (foundObject<labelIOList>("pointAreaNormalsParentMap"))
    {
        return;
    }

    regIOobject::store
    (
        new labelIOList
        (
            IOobject
            (
                "pointAreaNormalsParentMap",
                pointsInstance(),
                *this,
                IOobjectOption::NO_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::REGISTER
            ),
            meshPointLabels_
        )
    );
}


bool Foam::faMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    bool ok = fvMesh::writeObject(streamOpt, writeOnProc);

    // Persist the parent-mesh coupling addressing (faceLabels) alongside the
    // reduced-order mesh so it can be read back by faMesh(word, polyMesh).
    // Written under the parent-mesh registry to
    // <parentInstance>/finite-area[/<region>]/faMesh/faceLabels, matching the
    // location searched by readFaceLabels().
    if (!meshFaces_.empty() && hasParent())
    {
        const refPtr<fvMesh> pRef(parentMesh());
        const polyMesh& pMesh = pRef();

        // Match readCouplingAddressing(): the default single region (named
        // either polyMesh::defaultRegion or faMesh::prefix() by the solvers)
        // persists directly under <prefix>/faMesh with no <region> subdir.
        const word regionSub =
        (
            (
                areaName_.empty()
             || areaName_ == polyMesh::defaultRegion
             || areaName_ == faMesh::prefix()
            )
          ? word::null : areaName_
        );

        const fileName localDir =
        (
            regionSub.empty()
          ? (faMesh::prefix()/meshSubDir)
          : (faMesh::prefix()/regionSub/meshSubDir)
        );

        labelIOList faceLabelsIO
        (
            IOobject
            (
                "faceLabels",
                pMesh.facesInstance(),
                localDir,
                pMesh.thisDb(),
                IOobjectOption::NO_READ,
                IOobjectOption::AUTO_WRITE,
                IOobjectOption::NO_REGISTER
            ),
            meshFaces_
        );

        ok = faceLabelsIO.writeObject(streamOpt, writeOnProc) && ok;

        // Also persist the area-point -> parent-point addressing so that a
        // subsequent read (readCouplingAddressing) can rebuild pointLabels()
        // without recomputing the region point map. Written next to faceLabels.
        if (!meshPointLabels_.empty())
        {
            labelIOList pointLabelsIO
            (
                IOobject
                (
                    "pointLabels",
                    pMesh.facesInstance(),
                    localDir,
                    pMesh.thisDb(),
                    IOobjectOption::NO_READ,
                    IOobjectOption::AUTO_WRITE,
                    IOobjectOption::NO_REGISTER
                ),
                meshPointLabels_
            );

            ok = pointLabelsIO.writeObject(streamOpt, writeOnProc) && ok;
        }
    }

    return ok;
}


// ************************************************************************* //
