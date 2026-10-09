/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2025 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify i
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

#include "trackingInverseDistanceCellCellStencil.H"
#include "addToRunTimeSelectionTable.H"
#include "OBJstream.H"
#include "Time.H"
#include "fvMeshSubset.H"
#include "globalIndex.H"
#include "oversetFvPatch.H"
#include "zeroGradientFvPatchFields.H"
#include "syncTools.H"
#include "voxelMeshSearch.H"
#include "dynamicOversetFvMesh.H"
#include "waveMethod.H"
#include "patchFaceIntersection.H"

namespace
{
constexpr Foam::label maxPatchVoxelCount = 16*1024*1024;
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cellCellStencils
{
    defineTypeNameAndDebug(trackingInverseDistance, 0);
    addToRunTimeSelectionTable(cellCellStencil, trackingInverseDistance, mesh);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::cellCellStencils::trackingInverseDistance::markBoundaries
(
    const label zonei,
    const fvMesh& mesh,
    const vector& smallVec,

    const boundBox& bb,
    labelVector& nDivs,
    PackedList<2>& patchTypes,

    const labelList& cellMap,
    labelList& patchCellTypes
)
{
    // Mark voxels intersected by triangulated patch faces

    static bool hasWarned = false;

    const fvBoundaryMesh& pbm = mesh.boundary();

    patchTypes = patchCellType::OTHER;
    const labelVector strides(voxelMeshSearch::offset(nDivs));

    // Mark non-coupled patch boundaries
    forAll(pbm, patchi)
    {
        const fvPatch& fvp = pbm[patchi];
        const labelUList& fc = fvp.faceCells();

        if (!fvPatch::constraintType(fvp.type()))
        {
            //Info<< "Marking cells on proper patch " << fvp.name()
            //    << " with type " << fvp.type() << endl;
            const polyPatch& pp = fvp.patch();
            forAll(pp, i)
            {
                // Mark in overall patch types
                patchCellTypes[cellMap[fc[i]]] = patchCellType::PATCH;

                // Mark in voxel mesh
                boundBox faceBb(pp.points(), pp[i], false);
                faceBb.grow(smallVec);

                if (bb.overlaps(faceBb))
                {
                    patchFaceIntersection::fillFaceVoxels
                    (
                        patchTypes,
                        bb,
                        nDivs,
                        strides,
                        pp[i],
                        pp.points(),
                        smallVec,
                        patchCellType::PATCH
                    );
                }
            }
        }
    }

    // Override with overset boundaries
    forAll(pbm, patchi)
    {
        const fvPatch& fvp = pbm[patchi];
        const labelUList& fc = fvp.faceCells();

        if (isA<oversetFvPatch>(fvp))
        {
            //Info<< "Marking cells on overset patch " << fvp.name() << endl;
            const polyPatch& pp = fvp.patch();
            const vectorField::subField faceCentres(pp.faceCentres());
            forAll(pp, i)
            {
                // Mark in overall patch types
                patchCellTypes[cellMap[fc[i]]] = patchCellType::OVERSET;

                // Mark in voxel mesh
                boundBox faceBb(pp.points(), pp[i], false);
                faceBb.grow(smallVec);

                if (!bb.contains(faceCentres[i]))
                {
                    if (!hasWarned)
                    {
                        hasWarned = true;
                        WarningInFunction << "Patch " << pp.name()
                            << " : face at " << faceCentres[i]
                            << " is not inside search bounding box of"
                            << " voxel mesh " << bb << endl
                            << "    Is your 'searchBox' specification correct?"
                            << endl;
                    }
                }
                else
                {
                    // Test for having voxels already marked as patch
                    // -> voxel mesh is too coarse
                    if
                    (
                        patchFaceIntersection::overlapsFaceVoxels
                        (
                            bb,
                            nDivs,
                            strides,
                            pp[i],
                            pp.points(),
                            smallVec,
                            patchTypes,
                            static_cast<unsigned int>(patchCellType::PATCH)
                        )
                    )
                    {
                        // Determine number of voxels from number of cells
                        // in mesh
                        const labelVector& dim = mesh.geometricD();
                        labelVector refinedDivs;
                        label refinedVoxelCount = 0;
                        if
                        (
                            !patchFaceIntersection::refineVoxelDivisions
                            (
                                nDivs,
                                dim,
                                maxPatchVoxelCount,
                                refinedDivs,
                                refinedVoxelCount
                            )
                        )
                        {
                            FatalErrorInFunction
                                << "trackingInverseDistance cannot refine its "
                                << "voxel grid within the configured budget"
                                << " for zone " << zonei
                                << " in search box " << bb
                                << ". Current divisions: " << nDivs
                                << ", requested divisions: " << refinedDivs
                                << ", maximum voxel count: "
                                << maxPatchVoxelCount
                                << ". Overset patch: " << pp.name()
                                << ", face: " << i
                                << ", face bounding box: " << faceBb
                                << exit(FatalError);
                        }
                        nDivs = refinedDivs;

                        Pout<< "Voxel mesh too coarse. Bounding box "
                            << faceBb
                            << " contains both non-overset and overset patches"
                            << ". Refining voxel mesh to " << nDivs << endl;

                        return false;
                    }

                    patchFaceIntersection::fillFaceVoxels
                    (
                        patchTypes,
                        bb,
                        nDivs,
                        strides,
                        pp[i],
                        pp.points(),
                        smallVec,
                        patchCellType::OVERSET
                    );
                }
            }
        }
    }

    return true;
}


void Foam::cellCellStencils::trackingInverseDistance::markPatchesAsHoles
(
    PstreamBuffers& pBufs,
    const List<treeBoundBoxList>& patchBb,
    const List<labelVector>& patchDivisions,
    const PtrList<PackedList<2>>& patchParts,

    const label srcI,
    const label tgtI,
    labelList& allCellTypes
) const
{
    const treeBoundBoxList& srcPatchBbs = patchBb[srcI];
    const treeBoundBoxList& tgtPatchBbs = patchBb[tgtI];
    const labelList& tgtCellMap = meshParts_[tgtI].cellMap();
    const pointField srcTriangles
    (
        patchFaceIntersection::patchTriangles
        (
            meshParts_[srcI].subMesh()
        )
    );
    const scalar geometricTolerance
    (
        max(max(smallVec_[0], smallVec_[1]), smallVec_[2])
    );

    // 1. do processor-local src-tgt patch overlap
    {
        const treeBoundBox& srcPatchBb = srcPatchBbs[Pstream::myProcNo()];
        const treeBoundBox& tgtPatchBb = tgtPatchBbs[Pstream::myProcNo()];

        if (srcPatchBb.overlaps(tgtPatchBb))
        {
            const PackedList<2>& srcPatchTypes = patchParts[srcI];
            const labelVector& srcDivs = patchDivisions[srcI];
            const pointField localTriangles
            (
                patchFaceIntersection::trianglesOverlapping
                (
                    srcTriangles,
                    tgtPatchBb,
                    smallVec_
                )
            );

            forAll(tgtCellMap, tgtCelli)
            {
                label celli = tgtCellMap[tgtCelli];
                boundBox cBb(mesh_.cellBb(celli));
                cBb.grow(smallVec_);

                if
                (
                    voxelMeshSearch::overlaps
                    (
                        srcPatchBb,
                        srcDivs,
                        cBb,
                        srcPatchTypes,
                        static_cast<unsigned int>(patchCellType::PATCH)
                    )
                )
                {
                    if
                    (
                        patchFaceIntersection::intersectsCell
                        (
                            localTriangles,
                            mesh_,
                            celli,
                            geometricTolerance
                        )
                    )
                    {
                        allCellTypes[celli] = HOLE;
                    }
                }
            }
        }
    }


    // 2. Send over srcMesh bits that overlap tgt and do calculation
    pBufs.clear();
    for (const int procI : Pstream::allProcs())
    {
        if (procI != Pstream::myProcNo())
        {
            const treeBoundBox& srcPatchBb = srcPatchBbs[Pstream::myProcNo()];
            const treeBoundBox& tgtPatchBb = tgtPatchBbs[procI];

            if (srcPatchBb.overlaps(tgtPatchBb))
            {
                const pointField procTriangles
                (
                    patchFaceIntersection::trianglesOverlapping
                    (
                        srcTriangles,
                        tgtPatchBb,
                        smallVec_
                    )
                );

                // Send over complete patch voxel map. Tbd: could
                // subset
                UOPstream os(procI, pBufs);
                os  << srcPatchBb
                    << patchDivisions[srcI]
                    << patchParts[srcI]
                    << procTriangles;
            }
        }
    }
    pBufs.finishedSends();
    for (const int procI : Pstream::allProcs())
    {
        if (procI != Pstream::myProcNo())
        {
            const treeBoundBox& srcPatchBb = srcPatchBbs[procI];
            const treeBoundBox& tgtPatchBb = tgtPatchBbs[Pstream::myProcNo()];

            if (srcPatchBb.overlaps(tgtPatchBb))
            {
                UIPstream is(procI, pBufs);
                const treeBoundBox receivedBb(is);
                const labelVector srcDivs(is);
                const PackedList<2> srcPatchTypes(is);
                const pointField remoteTriangles(is);

                // Verify validity
                if (srcPatchBb != receivedBb)
                {
                    FatalErrorInFunction
                        << "proc:" << procI
                        << " srcPatchBb:" << srcPatchBb
                        << " receivedBb:" << receivedBb
                        << exit(FatalError);
                }

                forAll(tgtCellMap, tgtCelli)
                {
                    label celli = tgtCellMap[tgtCelli];
                    boundBox cBb(mesh_.cellBb(celli));
                    cBb.grow(smallVec_);

                    if
                    (
                        voxelMeshSearch::overlaps
                        (
                            srcPatchBb,
                            srcDivs,
                            cBb,
                            srcPatchTypes,
                            static_cast<unsigned int>(patchCellType::PATCH)
                        )
                    )
                    {
                        if
                        (
                            patchFaceIntersection::intersectsCell
                            (
                                remoteTriangles,
                                mesh_,
                                celli,
                                geometricTolerance
                            )
                        )
                        {
                            allCellTypes[celli] = HOLE;
                        }
                    }
                }
            }
        }
    }
}


void Foam::cellCellStencils::trackingInverseDistance::markDonors
(
    PstreamBuffers& pBufs,
    const List<treeBoundBoxList>& meshBb,
    const PtrList<voxelMeshSearch>& meshSearches,
    const labelList& allCellTypes,
    const labelList& allPatchTypes,

    const label srcI,
    const label tgtI,
    List<DynamicList<donorCandidate>>& allCandidates,
    labelList& allDonor
) const
{
    const treeBoundBoxList& srcBbs = meshBb[srcI];
    const treeBoundBoxList& tgtBbs = meshBb[tgtI];

    const fvMesh& srcMesh = meshParts_[srcI].subMesh();
    const labelList& srcCellMap = meshParts_[srcI].cellMap();
    const fvMesh& tgtMesh = meshParts_[tgtI].subMesh();
    const pointField& tgtCc = tgtMesh.cellCentres();
    const labelList& tgtCellMap = meshParts_[tgtI].cellMap();

    // 1. do processor-local src/tgt overlap
    {
        labelList tgtToSrcAddr;
        waveMethod::calculate(tgtMesh, srcMesh, tgtToSrcAddr);

        forAll(tgtCellMap, tgtCelli)
        {
            const label srcCelli = tgtToSrcAddr[tgtCelli];

            if
            (
                srcCelli != -1
             && allCellTypes[tgtCellMap[tgtCelli]] != HOLE
             && donorEligible
                (
                    allCellTypes[srcCellMap[srcCelli]],
                    allPatchTypes[srcCellMap[srcCelli]],
                    allPatchTypes[tgtCellMap[tgtCelli]]
                )
            )
            {
                label celli = tgtCellMap[tgtCelli];

                recordDonorCandidate
                (
                    celli,
                    tgtI,
                    srcI,
                    globalCells_.toGlobal(srcCellMap[srcCelli]),
                    allCandidates,
                    allDonor
                );
            }
        }
    }


    // 2. Send over tgtMesh bits that overlap src and do calculation on
    //    srcMesh.


    // (remote) processors where the tgt overlaps my src
    DynamicList<label> tgtOverlapProcs(Pstream::nProcs());
    // (remote) processors where the src overlaps my tgt
    DynamicList<label> srcOverlapProcs(Pstream::nProcs());
    for (const int procI : Pstream::allProcs())
    {
        if (procI != Pstream::myProcNo())
        {
            if (tgtBbs[procI].overlaps(srcBbs[Pstream::myProcNo()]))
            {
                tgtOverlapProcs.append(procI);
            }
            if (srcBbs[procI].overlaps(tgtBbs[Pstream::myProcNo()]))
            {
                srcOverlapProcs.append(procI);
            }
        }
    }



    // Indices of tgtcells to send over to each processor
    List<DynamicList<label>> tgtSendCells(Pstream::nProcs());
    forAll(srcOverlapProcs, i)
    {
        label procI = srcOverlapProcs[i];
        tgtSendCells[procI].reserve(tgtMesh.nCells()/srcOverlapProcs.size());
    }


    forAll(tgtCellMap, tgtCelli)
    {
        label celli = tgtCellMap[tgtCelli];
        if (srcOverlapProcs.size())
        {
            treeBoundBox subBb(mesh_.cellBb(celli));
            subBb.grow(smallVec_);

            forAll(srcOverlapProcs, i)
            {
                label procI = srcOverlapProcs[i];
                if (subBb.overlaps(srcBbs[procI]))
                {
                    tgtSendCells[procI].append(tgtCelli);
                }
            }
        }
    }

    // Send target cell centres to overlapping processors
    pBufs.clear();

    forAll(srcOverlapProcs, i)
    {
        label procI = srcOverlapProcs[i];
        const labelList& cellIDs = tgtSendCells[procI];

        labelList tgtPatchTypes(cellIDs.size());
        forAll(cellIDs, j)
        {
            tgtPatchTypes[j] = allPatchTypes[tgtCellMap[cellIDs[j]]];
        }

        UOPstream os(procI, pBufs);
        os << UIndirectList<point>(tgtCc, cellIDs);
        os << tgtPatchTypes;
    }
    pBufs.finishedSends();

    // Receive bits of target processors; find; send back
    (void)srcMesh.tetBasePtIs();
    forAll(tgtOverlapProcs, i)
    {
        label procI = tgtOverlapProcs[i];

        UIPstream is(procI, pBufs);
        pointList samples(is);
        labelList samplePatchTypes(is);

        labelList donors(samples.size(), -1);
        forAll(samples, sampleI)
        {
            const point& sample = samples[sampleI];
            label srcCelli = srcMesh.findCell(sample, polyMesh::CELL_TETS);
            if
            (
                srcCelli != -1
             && donorEligible
                (
                    allCellTypes[srcCellMap[srcCelli]],
                    allPatchTypes[srcCellMap[srcCelli]],
                    samplePatchTypes[sampleI]
                )
            )
            {
                donors[sampleI] = globalCells_.toGlobal(srcCellMap[srcCelli]);
            }
        }

        // Use same pStreamBuffers to send back.
        UOPstream os(procI, pBufs);
        os << donors;
    }
    pBufs.finishedSends();

    forAll(srcOverlapProcs, i)
    {
        label procI = srcOverlapProcs[i];
        const labelList& cellIDs = tgtSendCells[procI];

        UIPstream is(procI, pBufs);
        labelList donors(is);

        if (donors.size() != cellIDs.size())
        {
            FatalErrorInFunction<< "problem : cellIDs:" << cellIDs.size()
                << " donors:" << donors.size() << abort(FatalError);
        }

        forAll(donors, donorI)
        {
            label globalDonor = donors[donorI];

            if (globalDonor != -1)
            {
                label celli = tgtCellMap[cellIDs[donorI]];

                if (allCellTypes[celli] != HOLE)
                {
                    recordDonorCandidate
                    (
                        celli,
                        tgtI,
                        srcI,
                        globalDonor,
                        allCandidates,
                        allDonor
                    );
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellCellStencils::trackingInverseDistance::trackingInverseDistance
(
    const fvMesh& mesh,
    const dictionary& dict,
    const bool doUpdate
)
:
    inverseDistance(mesh, dict, false),
    globalCells_(mesh_.nCells())
{
    // Initialise donor cell
    globalDonor_.setSize(mesh_.nCells());
    globalDonor_ = -1;

    // Initialise mesh partitions
    const labelIOList& zoneID = this->zoneID();
    label nZones = gMax(zoneID)+1;

    labelList nCellsPerZone(nZones, Zero);
    forAll(zoneID, celli)
    {
        nCellsPerZone[zoneID[celli]]++;
    }
    Pstream::listReduce(nCellsPerZone, sumOp<label>());

    meshParts_.setSize(nZones);
    forAll(meshParts_, zonei)
    {
        meshParts_.set
        (
            zonei,
            new fvMeshSubset(mesh_, zonei, zoneID)
        );

        // Trigger early evaluation of mesh dimension
        // (in case there are locally zero cells in mesh)
        (void)meshParts_[zonei].subMesh().nGeometricD();
    }


    // Print a bit
    {
        Info<< typeName << " : detected " << nZones
            << " mesh regions" << endl;
        Info<< incrIndent;
        forAll(nCellsPerZone, zonei)
        {
            Info<< indent<< "zone:" << zonei
                << " nCells:" << nCellsPerZone[zonei]
                << endl;
        }
        Info<< decrIndent;
    }


    // Do geometry update
    if (doUpdate)
    {
        update();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellCellStencils::trackingInverseDistance::~trackingInverseDistance()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cellCellStencils::trackingInverseDistance::update()
{
    DebugInfo<< FUNCTION_NAME << " : Start of analysis" << endl;

    scalar layerRelax(dict_.getOrDefault("layerRelax", 1.0));
    const labelIOList& zoneID = this->zoneID();
    label nZones = meshParts_.size();

    // Update stored mesh partitions for geometry changes
    forAll(meshParts_, zonei)
    {
        pointField subPoints(mesh_.points(), meshParts_[zonei].pointMap());

        fvMesh& subMesh = meshParts_[zonei].subMesh();
        subMesh.movePoints(subPoints);
    }

    DebugInfo<< FUNCTION_NAME << " : Moved zone sub-meshes" << endl;

    // Calculate fast search structure for each zone
    List<labelVector> searchBoxDivisions(dict_.lookup("searchBoxDivisions"));
    if (searchBoxDivisions.size() != nZones)
    {
        FatalIOErrorInFunction(dict_) << "Number of searchBoxDivisions "
            << searchBoxDivisions.size()
            << " should  equal the number of zones " << nZones
            << exit(FatalIOError);
    }

    const boundBox& allBb = mesh_.bounds();

    List<treeBoundBoxList> meshBb(nZones);

    // Determine zone meshes and bounding boxes
    {
        // Per processor, per zone the bounding box
        List<treeBoundBoxList> procBb(Pstream::nProcs());
        procBb[Pstream::myProcNo()].setSize(nZones);

        forAll(meshParts_, zonei)
        {
            const fvMesh& subMesh = meshParts_[zonei].subMesh();

            if (subMesh.nPoints())
            {
                procBb[Pstream::myProcNo()][zonei] =
                    treeBoundBox(subMesh.points());
                procBb[Pstream::myProcNo()][zonei].inflate(1e-6);
            }
            else
            {
                // No part of zone on this processor. Make up bb.
                procBb[Pstream::myProcNo()][zonei] = treeBoundBox
                (
                    allBb.min() - 2*allBb.span(),
                    allBb.min() - allBb.span()
                );
                procBb[Pstream::myProcNo()][zonei].inflate(1e-6);
            }
        }

        Pstream::allGatherList(procBb);

        // Move local bounding boxes to per-mesh indexing
        forAll(meshBb, zonei)
        {
            treeBoundBoxList& bbs = meshBb[zonei];
            bbs.setSize(Pstream::nProcs());
            forAll(procBb, proci)
            {
                bbs[proci] = procBb[proci][zonei];
            }
        }
    }

    DebugInfo<< FUNCTION_NAME << " : Calculated bounding boxes" << endl;


    // Determine patch bounding boxes. These are either global and provided
    // by the user or processor-local as a copy of the mesh bounding box.

    List<treeBoundBoxList> patchBb(nZones);
    List<labelVector> patchDivisions(searchBoxDivisions);
    PtrList<PackedList<2>> patchParts(nZones);
    labelList allPatchTypes(mesh_.nCells(), OTHER);

    {
        treeBoundBox globalPatchBb;
        if (dict_.readIfPresent("searchBox", globalPatchBb))
        {
            // All processors, all zones have the same bounding box
            patchBb = treeBoundBoxList(Pstream::nProcs(), globalPatchBb);
        }
        else
        {
            // Use the meshBb (differing per zone, per processor)
            patchBb = meshBb;
        }
    }

    forAll(patchParts, zonei)
    {
        while (true)
        {
            label nVoxels = 0;
            if
            (
                !patchFaceIntersection::checkedVoxelCount
                (
                    patchDivisions[zonei],
                    maxPatchVoxelCount,
                    nVoxels
                )
            )
            {
                FatalErrorInFunction
                    << "trackingInverseDistance voxel budget exceeded"
                    << " for zone " << zonei
                    << " in search box "
                    << patchBb[zonei][Pstream::myProcNo()]
                    << ". Divisions: " << patchDivisions[zonei]
                    << ", maximum voxel count: " << maxPatchVoxelCount
                    << exit(FatalError);
            }

            patchParts.set
            (
                zonei,
                new PackedList<2>(nVoxels)
            );

            bool ok = markBoundaries
            (
                zonei,
                meshParts_[zonei].subMesh(),
                smallVec_,

                patchBb[zonei][Pstream::myProcNo()],
                patchDivisions[zonei],
                patchParts[zonei],

                meshParts_[zonei].cellMap(),
                allPatchTypes
            );

            //if (returnReduceAnd(ok))
            if (ok)
            {
                break;
            }
        }
    }
    DebugInfo<< FUNCTION_NAME << " : Calculated boundary voxel meshes" << endl;


    PtrList<voxelMeshSearch> meshSearches(meshParts_.size());
    forAll(meshParts_, zonei)
    {
        meshSearches.set
        (
            zonei,
            new voxelMeshSearch
            (
                meshParts_[zonei].subMesh(),
                patchBb[zonei][Pstream::myProcNo()],
                patchDivisions[zonei],
                true
            )
        );
    }

    DebugInfo<< FUNCTION_NAME << " : Constructed cell voxel meshes" << endl;

    PtrList<fvMesh> voxelMeshes;
    if (debug)
    {
        voxelMeshes.setSize(meshSearches.size());
        forAll(meshSearches, zonei)
        {
            IOobject meshIO
            (
                word("voxelMesh" + Foam::name(zonei)),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ
            );

            Pout<< "Constructing voxel mesh " << meshIO.objectPath() << endl;
            voxelMeshes.set(zonei, meshSearches[zonei].makeMesh(meshIO));
        }
    }


    if (debug & 2)
    {
        forAll(patchParts, zonei)
        {
            const labelList voxelTypes(patchParts[zonei].values());
            const fvMesh& vm = voxelMeshes[zonei];
            tmp<volScalarField> tfld
            (
                createField<label>
                (
                    vm,
                    "patchTypes",
                    voxelTypes
                )
            );
            tfld().write();
        }
    }

    // Current best guess for cells
    labelList allCellTypes(mesh_.nCells(), CALCULATED);
    labelListList allStencil(mesh_.nCells());
    List<DynamicList<donorCandidate>> allCandidates(mesh_.nCells());
    // zoneID of donor
    labelList allDonorID(mesh_.nCells(), -1);

    const globalIndex globalCells(mesh_.nCells());

    PstreamBuffers pBufs;

    DebugInfo<< FUNCTION_NAME << " : Allocated donor-cell structures" << endl;

    for (label srci = 0; srci < meshParts_.size()-1; srci++)
    {
        for (label tgti = srci+1; tgti < meshParts_.size(); tgti++)
        {
            markPatchesAsHoles
            (
                pBufs,

                patchBb,
                patchDivisions,
                patchParts,

                srci,
                tgti,
                allCellTypes
            );
            markPatchesAsHoles
            (
                pBufs,

                patchBb,
                patchDivisions,
                patchParts,

                tgti,
                srci,
                allCellTypes
            );
        }
    }

    for (label srci = 0; srci < meshParts_.size()-1; srci++)
    {
        for (label tgti = srci+1; tgti < meshParts_.size(); tgti++)
        {
            markDonors
            (
                pBufs,
                meshBb,
                meshSearches,
                allCellTypes,   // to exclude hole donors
                allPatchTypes,

                tgti,
                srci,
                allCandidates,
                allDonorID
            );
            markDonors
            (
                pBufs,
                meshBb,
                meshSearches,
                allCellTypes,   // to exclude hole donors
                allPatchTypes,

                srci,
                tgti,
                allCandidates,
                allDonorID
            );
        }
    }

    orderDonorCandidates(zoneID, allCandidates, allStencil, allDonorID);

    DebugInfo<< FUNCTION_NAME << " : Determined holes and donor-acceptors"
        << endl;

    if ((debug & 2) && mesh_.time().writeTime())
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allCellTypes", allCellTypes)
        );
        tfld().write();

        tmp<volScalarField> tallDonorID
        (
            createField(mesh_, "allDonorID", allDonorID)
        );
        tallDonorID().write();
    }


    // Use the patch types and weights to decide what to do
    forAll(allPatchTypes, celli)
    {
        if (allCellTypes[celli] != HOLE)
        {
            switch (allPatchTypes[celli])
            {
                case OVERSET:
                {
                    // Require interpolation. See if possible.

                    if (allStencil[celli].size())
                    {
                        allCellTypes[celli] = INTERPOLATED;
                    }
                    else
                    {
                        allCellTypes[celli] = HOLE;
                    }
                }
            }
        }
    }
    DebugInfo<< FUNCTION_NAME << " : Removed bad donors" << endl;

    if ((debug & 2) && mesh_.time().writeTime())
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allCellTypes_patch", allCellTypes)
        );
        tfld().write();

        tmp<volScalarField> tfldOld
        (
            createField(mesh_, "allCellTypes_old", cellTypes_)
        );
        tfldOld().write();
    }

    boolList reactivationCandidates(mesh_.nCells(), false);
    boolList sameZoneFallback(mesh_.nCells(), false);
    labelListList reactivationStencil(mesh_.nCells());
    scalarField allWeight(mesh_.nCells(), Zero);

    classifyAndWalkFront
    (
        globalCells,
        zoneID,
        allPatchTypes,
        allCandidates,
        layerRelax,
        dict_.getOrDefault("holeLayers", 1),
        dict_.getOrDefault("useLayer", -1),
        allCellTypes,
        allStencil,
        allWeight,
        reactivationCandidates,
        sameZoneFallback,
        reactivationStencil
    );

    DebugInfo<< FUNCTION_NAME << " : Flood-filled holes and walked front"
        << endl;

    if ((debug & 2) && mesh_.time().writeTime())
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "allCellTypes_front", allCellTypes)
        );
        tfld().write();
    }

    labelListList compactStencil(allStencil);
    List<Map<label>> compactStencilMap;
    mapDistribute map(globalCells, compactStencil, compactStencilMap);

    // Restore cells that were holes but are calculated on the new geometry.
    // A missing or no-longer-valid cross-zone donor is replaced by a valid
    // same-zone neighbour when possible.
    {
        labelList compactCellTypesForRecovery(allCellTypes);
        map.distribute(compactCellTypesForRecovery);

        DynamicList<label> fallbackCells;
        forAll(reactivationCandidates, celli)
        {
            if (!reactivationCandidates[celli])
            {
                continue;
            }

            bool hasValidDonor = false;

            if (!sameZoneFallback[celli])
            {
                const labelList& slots = compactStencil[celli];
                forAll(slots, sloti)
                {
                    const label donorType = compactCellTypesForRecovery[slots[sloti]];
                    if
                    (
                        donorType != HOLE
                     && donorType != POROUS
                     && (allowInterpolatedDonors_ || donorType != INTERPOLATED)
                    )
                    {
                        hasValidDonor = true;
                        break;
                    }
                }
            }

            if (hasValidDonor)
            {
                allStencil[celli] = reactivationStencil[celli];
                allCellTypes[celli] = INTERPOLATED;
                allWeight[celli] = 1.0;
            }
            else
            {
                sameZoneFallback[celli] = true;
                fallbackCells.append(celli);
            }
        }

        boolList fallbackDonorMask(mesh_.nCells(), false);
        forAll(allCellTypes, celli)
        {
            fallbackDonorMask[celli] =
                !reactivationCandidates[celli]
             &&
                (
                    allCellTypes[celli] == CALCULATED
                 || (allowInterpolatedDonors_ && allCellTypes[celli] == INTERPOLATED)
                );
        }

        labelListList fallbackCellCells;
        pointListList fallbackCellCentres;
        globalSameZoneCellCells
        (
            globalCells,
            mesh_,
            zoneID,
            fallbackDonorMask,
            fallbackCells,
            fallbackCellCells,
            fallbackCellCentres
        );

        forAll(fallbackCells, i)
        {
            const label celli = fallbackCells[i];
            const labelList& donors = fallbackCellCells[celli];

            if (donors.empty())
            {
                allCellTypes[celli] = HOLE;
                allStencil[celli].clear();
                allWeight[celli] = 0.0;
            }
            else
            {
                allStencil[celli].setSize(1);
                allStencil[celli][0] = donors[0];
                allCellTypes[celli] = INTERPOLATED;
                allWeight[celli] = 1.0;
            }
        }
    }

    labelListList finalCompactStencil(allStencil);
    List<Map<label>> finalCompactStencilMap;
    mapDistribute finalMap
    (
        globalCells,
        finalCompactStencil,
        finalCompactStencilMap
    );
    labelList compactCellTypes(allCellTypes);
    finalMap.distribute(compactCellTypes);
    labelList compactCellPatchTypes(allPatchTypes);
    finalMap.distribute(compactCellPatchTypes);

    // Convert cell-cell addressing to stencil in compact notation

    cellTypes_.transfer(allCellTypes);
    cellStencil_.setSize(mesh_.nCells());
    cellInterpolationWeights_.setSize(mesh_.nCells());
    DynamicList<label> interpolationCells;

    label nHoleDonors = 0;
    forAll(cellTypes_, celli)
    {
        if (cellTypes_[celli] == INTERPOLATED)
        {
            const labelList& slots = finalCompactStencil[celli];
            if (slots.size())
            {
                if
                (
                    !donorEligible
                    (
                        compactCellTypes[slots[0]],
                        compactCellPatchTypes[slots[0]],
                        allPatchTypes[celli]
                    )
                )
                {
                    cellTypes_[celli] = POROUS;
                    cellStencil_[celli].clear();
                    cellInterpolationWeights_[celli].clear();
                    nHoleDonors++;
                }
                else
                {
                    cellStencil_[celli].transfer(allStencil[celli]);
                    cellInterpolationWeights_[celli].setSize(1);
                    cellInterpolationWeights_[celli][0] = 1.0;
                    interpolationCells.append(celli);
                }
            }
            else
            {
                cellTypes_[celli] = HOLE;
                cellStencil_[celli].clear();
                cellInterpolationWeights_[celli].clear();
                allWeight[celli] = 0.0;
            }
        }
        else
        {
            cellStencil_[celli].clear();
            cellInterpolationWeights_[celli].clear();
        }
    }

    reduce(nHoleDonors, sumOp<label>());

    DebugInfo<< FUNCTION_NAME << "nHole Donors : " << nHoleDonors << endl;

    interpolationCells_.transfer(interpolationCells);

    List<Map<label>> compactMap;
    cellInterpolationMap_.reset
    (
        new mapDistribute
        (
            globalCells,
            cellStencil_,
            compactMap
        )
    );
    cellInterpolationWeight_.transfer(allWeight);
    oversetFvMeshBase::correctBoundaryConditions
    <
        volScalarField,
        oversetFvPatchField<scalar>,
        false
    >(cellInterpolationWeight_.boundaryFieldRef());


    if ((debug & 2) && mesh_.time().writeTime())
    {
        // Dump mesh
        mesh_.time().write();

        // Dump stencil
        mkDir(mesh_.time().timePath());
        OBJstream str(mesh_.time().timePath()/"injectionStencil.obj");
        Pout<< typeName << " : dumping injectionStencil to "
            << str.name() << endl;
        pointField cc(mesh_.cellCentres());
        cellInterpolationMap().distribute(cc);

        forAll(cellStencil_, celli)
        {
            const labelList& slots = cellStencil_[celli];
            if (slots.size())
            {
                const point& accCc = mesh_.cellCentres()[celli];
                forAll(slots, i)
                {
                    const point& donorCc = cc[slots[i]];
                    str.writeLine(accCc, 0.1*accCc+0.9*donorCc);
                }
            }
        }
    }

    DebugInfo<< FUNCTION_NAME << " : Transferred donor to stencil" << endl;


    // Extend stencil to get inverse distance weighted neighbours
    createStencil(globalCells, allowHoleDonors_);
    DebugInfo<< FUNCTION_NAME << " : Extended stencil" << endl;

    recoverSameZoneInterpolations
    (
        globalCells,
        zoneID,
        reactivationCandidates,
        sameZoneFallback,
        allowInterpolatedDonors_,
        cellTypes_,
        cellStencil_,
        cellInterpolationWeights_,
        cellInterpolationWeight_.primitiveFieldRef(),
        interpolationCells_,
        cellInterpolationMap_
    );
    oversetFvMeshBase::correctBoundaryConditions
    <
        volScalarField,
        oversetFvPatchField<scalar>,
        false
    >(cellInterpolationWeight_.boundaryFieldRef());

    // Optional: convert hole cells next to non-hole cells into
    // interpolate-from-neighbours (of cell type SPECIAL)
    if (allowHoleDonors_)
    {
        holeExtrapolationStencil(globalCells_);
    }


    if ((debug & 2) && mesh_.time().writeTime())
    {
        // Dump weight
        cellInterpolationWeight_.instance() = mesh_.time().timeName();
        cellInterpolationWeight_.write();

        // Dump stencil
        mkDir(mesh_.time().timePath());
        OBJstream str(mesh_.time().timePath()/"stencil.obj");
        Pout<< typeName << " : dumping to " << str.name() << endl;
        pointField cc(mesh_.cellCentres());
        cellInterpolationMap().distribute(cc);

        forAll(cellStencil_, celli)
        {
            const labelList& slots = cellStencil_[celli];
            if (slots.size())
            {
                const point& accCc = mesh_.cellCentres()[celli];
                forAll(slots, i)
                {
                    const point& donorCc = cc[slots[i]];
                    str.writeLine(accCc, 0.1*accCc+0.9*donorCc);
                }
            }
        }
    }


    // Print some stats
    Info<< this->info();

    DebugInfo<< FUNCTION_NAME << " : Finished analysis" << endl;

    // Tbd: detect if anything changed. Most likely it did!
    return true;
}


// ************************************************************************* //
