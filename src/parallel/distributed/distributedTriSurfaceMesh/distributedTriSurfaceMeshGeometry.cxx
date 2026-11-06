/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Keysight Technologies
    Copyright (C) 2026 Mark Olesen
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

Description
    Mostly surface geometry related routines
    (not searchable surface routines)

\*---------------------------------------------------------------------------*/

#include "distributedTriSurfaceMesh.H"
#include "triangle.H"
#include "triSurfaceTools.H"
#include "matchPoints.H"
#include "globalIndex.H"

#include "decompositionMethod.H"
#include "treeDecomp.H"
#include "profiling.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Calculate point-faces.
// uses "mini-format" (map from smallest vertex to face number)
void calculate_pointFaces
(
    const Foam::label numPoints,
    const Foam::UList<Foam::labelledTri>& tris,
    Foam::labelListList& pointFaces,
    const Foam::bitSet& eligiblePoints
)
{
    using namespace Foam;

    pointFaces.resize_nocopy(numPoints);

    if (eligiblePoints.none())
    {
        for (auto& sub : pointFaces)
        {
            sub.clear();
        }
        return;
    }

    // mini format:
    // only map from smallest vertex to face number

    // The output list sizes
    labelList sizes(numPoints, Foam::zero{});

    for (const auto& tri : tris)
    {
        const auto slot = tri.min();

        if (eligiblePoints.test(slot) && slot < numPoints)
        {
            ++sizes[slot];
        }
    }


    // Size the output
    forAll(sizes, i)
    {
        pointFaces[i].resize(sizes[i]);
    }

    // Fill output
    sizes = 0;
    forAll(tris, trii)
    {
        const auto& tri = tris[trii];

        const auto slot = tri.min();

        if (eligiblePoints.test(slot) && slot < numPoints)
        {
            pointFaces[slot][sizes[slot]++] = trii;
        }
    }
}


// Calculate point-faces
// "regular" - full density mapping
// "mini"    - map from smallest vertex to face number
void calculate_pointFaces
(
    const Foam::label numPoints,
    const Foam::UList<Foam::labelledTri>& tris,
    Foam::labelListList& pointFaces,
    const bool miniFormat = false
)
{
    using namespace Foam;

    if (!miniFormat)
    {
        // The dense format
        Foam::invertManyToMany(numPoints, tris, pointFaces);
    }
    else
    {
        // mini format:
        // only map from smallest vertex to face number

        // The output list sizes
        labelList sizes(numPoints, Foam::zero{});

        for (const auto& tri : tris)
        {
            const auto slot = tri.min();

            if (slot >= 0 && slot < numPoints)
            {
                ++sizes[slot];
            }
        }

        // Size output
        pointFaces.resize_nocopy(numPoints);
        forAll(sizes, i)
        {
            pointFaces[i].resize(sizes[i]);
        }

        // Fill output
        sizes = 0;
        forAll(tris, trii)
        {
            const auto& tri = tris[trii];

            const auto slot = tri.min();

            if (slot >= 0 && slot < numPoints)
            {
                pointFaces[slot][sizes[slot]++] = trii;
            }
        }
    }
}


// Calculate face-faces
void calculate_faceFaces
(
    const Foam::triSurface& surf,
    const Foam::labelListList& pointFaces,
    Foam::labelListList& faceFaces
)
{
    using namespace Foam;

    faceFaces.resize_nocopy(surf.size());

    DynamicList<label> nbrs(4);

    forAll(surf, facei)
    {
        const auto& f = surf[facei];

        nbrs.clear();

        forAll(f, fp)
        {
            const edge e = f.edge(fp);

            for (label otherFacei : pointFaces[e[0]])
            {
                if
                (
                    (otherFacei != facei)
                 && surf[otherFacei].contains(e)
                )
                {
                    nbrs.push_uniq(otherFacei);
                }
            }
        }
        faceFaces[facei] = std::move(nbrs);
    }
}


// Calculate face-faces (compact layout)
// FUTURE?
#if 0
void calculate_faceFacesCompact
(
    const Foam::triSurface& surf,
    Foam::CompactListList<Foam::label>& pointFaces,
    Foam::CompactListList<Foam::label>& faceFaces
)
{
    using namespace Foam;

    const label numFaces = surf.size();

    // Sanity checks
    if (numFaces <= 1)
    {
        // 0/1 faces => no neighbours
        faceFaces.clear();
        faceFaces.resize(numFaces, 0);
        return;
    }

    // Recover any pre-existing memory

    DynamicList<label> offsets(std::move(faceFaces.offsets()));
    offsets.clear();

    DynamicList<label> values(std::move(faceFaces.values()));
    values.clear();

    // Initial size (over)estimate :
    // each triangle has three neighbours

    offsets.reserve_nocopy(numFaces+1);
    values.reserve_nocopy(3*numFaces);

    // Local work array
    DynamicList<label> nbrs(4);

    // Start offset for face0 = zero.
    offsets.push_back(0);

    forAll(surf, facei)
    {
        const auto& f = surf[facei];

        nbrs.clear();

        forAll(f, fp)
        {
            const edge e = f.edge(fp);

            for (label otherFacei : pointFaces[e[0]])
            {
                if (otherFacei != facei && surf[otherFacei].contains(e))
                {
                    nbrs.push_uniq(otherFacei);
                }
            }
        }

        // Accumulate content, update offsets
        values.push_back(nbrs);
        offsets.push_back(values.size());
    }

    // Transfer
    faceFaces.offsets() = std::move(offsets);
    faceFaces.values() = std::move(values);
}
#endif

} // End anonymous namespace


// * * * * * * * * * * * * * Surface Exchange Helper * * * * * * * * * * * * //

namespace
{

using receiveOrder = Foam::distributedTriSurfaceMesh::receiveOrder;

//- Sparse point-to-point exchange of surfaces between all ranks.
//
//  - sendOp(proci)  -> autoPtr<triSurface>   (nullptr: nothing to send)
//  - sendGeom(proci, const triSurface&)      post non-blocking send
//  - recvGeom(proci, triSurface&)            post non-blocking receive
//  - procOp(proci, triSurface&)              consume a completed receive
//
//  The geometry send/receive is two messages (points, faces), so each
//  receive occupies a contiguous block of requests, which is completed
//  as a unit (pairwise). Does not handle the local rank.
//
//  With syncPatches, the patches of all sent surfaces are merged and
//  shared once and assigned to every received surface.

template
<
    class SendOp,
    class SendGeomOp,
    class RecvGeomOp,
    class ProcessOp,
    class PreRecvOp = Foam::identityOp
>
void exchangeSurfacesGeneric
(
    const int communicator,
    const receiveOrder recvOrder,

    //! The send or extract operation
    SendOp&& sendOp,

    SendGeomOp&& geomSendOp,
    RecvGeomOp&& geomRecvOp,

    //! Processing of the received surface
    ProcessOp&& procOp,
    const bool syncPatches = false,

    //! Called with the receive sizes (nPoints, nFaces) per rank,
    //! before any receives are posted. Eg, for pre-sizing storage.
    PreRecvOp&& preRecvOp = PreRecvOp()
)
{
    using namespace Foam;

    if (!UPstream::is_parallel(communicator)) return;

    const auto myProci = UPstream::myProcNo(communicator);
    const auto numProc = UPstream::nProcs(communicator);

    PtrList<triSurface> sendBuffers(numProc);
    PtrList<triSurface> recvBuffers(numProc);

    // Meta-data: (nPoints, nFaces), zero if nothing is sent
    List<labelPair> sendSizes(numProc, labelPair(0, 0));
    List<labelPair> recvSizes(numProc, labelPair(0, 0));

    const auto startOfRequests = UPstream::nRequests();


    // ~~~~~~~~~~~~~
    // Initial stage
    //
    // Surface extract, storage, remember sizing,
    // optionally collect surface patch information.

    for (label proci = 0; proci < numProc; ++proci)
    {
        // No exchange with self
        if (proci == myProci) continue;

        // [Callback] for extraction (does no sending)
        autoPtr<triSurface> ptr(sendOp(proci));

        if (ptr && (ptr->size() > 0))
        {
            sendSizes[proci] = { ptr->points().size(), ptr->size() };
            sendBuffers.set(proci, std::move(ptr));
        }
    }

    // [optional]
    // Common patch definitions consolidated from the send data

    List<geometricSurfacePatch> commonPatches;

    if (syncPatches)
    {
        Map<wordPair> patchMap;  // (index : (name, type))

        // Accumulate the local information
        forAll(sendBuffers, i)
        {
            if (const auto* ptr = sendBuffers.get(i))
            {
                for (const auto& p : ptr->patches())
                {
                    patchMap.emplace(p.index(), p.name(), p.geometricType());
                }
            }
        }

        // Combine is per value:
        // for a duplicate index, keep the existing entry
        Pstream::mapCombineReduce
        (
            patchMap,
            // No merging of content, we just want to collect entries
            [](wordPair&, const wordPair&) {},
            UPstream::msgType(),
            communicator
        );

        commonPatches.resize(patchMap.size());

        auto output = commonPatches.begin();

        for (const auto& iter : patchMap.csorted())
        {
            *output++ = geometricSurfacePatch(iter.key(), iter.val());
        }
    }

    // ~~~~~~~~~~~~~~~~~~
    // Intermediate stage
    //
    // Exchange the sizes (PEX approach)

    UPstream::allToAll(sendSizes, recvSizes, communicator);

    preRecvOp(recvSizes);   // [Callback] for any pre-receive operations


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Initiate receive communication
    //
    // Post all receives, recording the request block of each

    const auto beginOfRecvRequests = UPstream::nRequests();

    // Track rank and (first last) request positions [absolute]
    struct pendingRecv
    {
        int proci;

        // Absolute request positions
        int reqs_first;    // first of the requests
        int reqs_last;     // last (one-past) of the requests

        int nreqs() const noexcept { return (reqs_last - reqs_first); }

        void wait() const { UPstream::waitRequests(reqs_first, nreqs()); }
    };

    DynamicList<pendingRecv> pending(numProc);

    for (label proci = 0; proci < numProc; ++proci)
    {
        const auto [numPoints, numFaces] = recvSizes[proci];

        // No exchange with self, no receives for empty surfaces
        if (proci == myProci || numFaces <= 0)
        {
            continue;
        }

        auto& sub = recvBuffers.emplace_set(proci);
        sub.storedPoints().resize_nocopy(numPoints);
        sub.storedFaces().resize_nocopy(numFaces);

        if (syncPatches)
        {
            sub.patches() = commonPatches;
        }

        const auto beg = UPstream::nRequests();

        geomRecvOp(proci, sub);  // [Callback] for receive

        const auto end = UPstream::nRequests();

        // Note: emplace_back() doesn't work with the ad hoc struct
        pending.push_back( pendingRecv{ proci, beg, end } );
    }

    // Done posting receives
    const auto numRecvRequests(UPstream::nRequests() - beginOfRecvRequests);

    //
    // Post all sends
    //
    for (label proci = 0; proci < numProc; ++proci)
    {
        if (const auto* ptr = sendBuffers.get(proci))
        {
            geomSendOp(proci, *ptr);
        }
    }

    // ~~~~~~~~~~~~~~~~~~~~~~
    // Complete the receives.
    // We are almost there...
    //
    // sequential :
    //   - just walk through 'pending', they are already in sequential order
    //
    // any/polling :
    //   - from 'pending' determine how many requests correspond to the
    //     given proci and their position. Poll with wait-some to get
    //     as much progress as possible per call.
    //     From the finished requests, lookup and count down how many were
    //     expected. When the number of waiting is 0 for that proci,
    //     can perform the actual work on the data.

    if (recvOrder == receiveOrder::sequenced)
    {
        for (const auto& recv : pending)
        {
            recv.wait();
            procOp(recv.proci, recvBuffers[recv.proci]);
        }
    }
    else   // (receiveOrder::unsequenced) == polling
    {
        // Map request (relative to beginOfRecvRequests) to its pending entry
        // As sanity check :
        //  numRecvRequests = 2*pending.size()

        labelList requestToPending(numRecvRequests, label(-1));
        labelList numWaiting(pending.size(), label(0));

        forAll(pending, i)
        {
            const auto& recv = pending[i];

            numWaiting[i] = recv.nreqs();

            requestToPending.slice
            (
                (recv.reqs_first - beginOfRecvRequests),
                recv.nreqs()
            ) = i;
        }

        // Polling with count down and dispatch
        DynamicList<int> indices(pending.size());

        while
        (
            UPstream::waitSomeRequests
            (
                beginOfRecvRequests,
                numRecvRequests,
               &indices
            )
        )
        {
            for (const int reqi : indices)
            {
                const label i = requestToPending[reqi];

                // The block is complete when all of its requests are done
                if (i >= 0 && --numWaiting[i] == 0)
                {
                    auto proci = pending[i].proci;

                    // [Callback] for processing the received surface
                    procOp(proci, recvBuffers[proci]);
                }
            }
        }
    }

    // Sends *must* complete before the buffers go out of scope
    UPstream::waitRequests(startOfRequests);
}


// Calling patterns
// The member exchangeSurfaces() wraps exchangeSurfacesGeneric with
// the sendSurfaceGeometry/recvSurfaceGeometry:
//
// 1. Order-independent processing (cf. markMixed)
//
//      exchangeSurfaces
//      (
//          communicator,
//          receiveOrder::unsequenced,
//          baseMessageTag,
//
//          [&](label proci) -> autoPtr<triSurface>
//          {
//              if (condition)
//              {
//                  return autoPtr<triSurface>::New(s.subsetMesh(sel));
//              }
//              return nullptr;
//          },
//
//          [&](label proci, triSurface& sub) { do_something }
//      );
//
// 2. Processing in rank order (cf. merge of the received surfaces).
//    Only the order and the last lambda differ:
//
//          receiveOrder::sequenced, ...
//
//          [&](label proci, triSurface& sub) { merge( ... ); }
//
// 3. Non-default transport: call exchangeSurfacesGeneric directly with
//    custom sendGeom/recvGeom lambdas.

} // End anonymous namespace


// * * * * * * * * * * * * * Communication Handling  * * * * * * * * * * * * //

template<class SendOp, class ProcessOp, class PreRecvOp>
void Foam::distributedTriSurfaceMesh::exchangeSurfaces
(
    const int communicator,
    const receiveOrder recvOrder,
    const int baseMessageTag,
    SendOp&& sendOp,
    ProcessOp&& procOp,
    PreRecvOp&& preRecvOp
)
{
    exchangeSurfacesGeneric
    (
        communicator,
        recvOrder,
        std::forward<SendOp>(sendOp),

        [&](int proci, const triSurface& s)
        {
            sendSurfaceGeometry(proci, s, baseMessageTag, communicator);
        },

        [&](int proci, triSurface& s)
        {
            recvSurfaceGeometry(proci, s, baseMessageTag, communicator);
        },

        std::forward<ProcessOp>(procOp),

        false,  // syncPatches

        std::forward<PreRecvOp>(preRecvOp)
    );
}


// Non-blocking send of surface primitives.
//
// The primitives have different message tags to disambiguate
// the send/recv ordering.
void Foam::distributedTriSurfaceMesh::sendSurfaceGeometry
(
    int toProcNo,
    const UList<point>& points,
    const UList<labelledTri>& faces,
    int baseMessageTag,
    int communicator,
    UPstream::commsTypes commsType
)
{
    if (UPstream::is_parallel(communicator))
    {
        const auto startOfRequests = UPstream::nRequests();

        // points
        UOPstream::write
        (
            UPstream::commsTypes::nonBlocking,
            toProcNo,
            points,
            (baseMessageTag + 101),
            communicator
        );

        // faces
        UOPstream::write
        (
            UPstream::commsTypes::nonBlocking,
            toProcNo,
            faces,
            (baseMessageTag + 102),
            communicator
        );

        // Wait for completion if the caller did not specify non-blocking
        if (commsType != UPstream::commsTypes::nonBlocking)
        {
            UPstream::waitRequests(startOfRequests);
        }
    }
}


// Non-blocking receive of surface primitives.
//
// The primitives have different message tags to disambiguate
// the send/recv ordering.
void Foam::distributedTriSurfaceMesh::recvSurfaceGeometry
(
    int fromProcNo,
    List<point>& points,
    List<labelledTri>& faces,
    int baseMessageTag,
    int communicator,
    UPstream::commsTypes commsType
)
{
    if (UPstream::is_parallel(communicator))
    {
        const auto startOfRequests = UPstream::nRequests();

        // points
        UIPstream::read
        (
            UPstream::commsTypes::nonBlocking,
            fromProcNo,
            points,
            (baseMessageTag + 101),
            communicator
        );

        // faces
        UIPstream::read
        (
            UPstream::commsTypes::nonBlocking,
            fromProcNo,
            faces,
            (baseMessageTag + 102),
            communicator
        );

        // Wait for completion if the caller did not specify non-blocking
        if (commsType != UPstream::commsTypes::nonBlocking)
        {
            UPstream::waitRequests(startOfRequests);
        }
    }
}


// * * * * * * * * * * * Private Static Member Functions  * * * * * * * * * * //

Foam::distributedTriSurfaceMesh::sendStyles
Foam::distributedTriSurfaceMesh::analyzeSendMaps
(
    const triSurface& surf,
    const UList<labelList>& faceSendMap,
    int communicator
)
{
    const auto numProc = UPstream::nProcs(communicator);

    auto sendStyle = sendStyles::mode_general;

    // The first non-empty send map on this rank
    const label firstProci =
    (
        ListOps::find_if
        (
            faceSendMap,
            [](const auto& map) { return !map.empty(); }
        )
    );

    // Check if any sub-ranks are sending
    bool subranks_sending =
    (
        UPstream::is_subrank(communicator) && (firstProci >= 0)
    );

    // Only need result on the master, so MPI_Reduce instead of Allreduce
    UPstream::mpiReduceOr(&subranks_sending, 1, communicator);


    // Determine the send mode.
    // decomposing if subranks do not send
    if
    (
        UPstream::master(communicator)
     && (!subranks_sending && firstProci >= 0)
    )
    do
    {
        sendStyle = sendStyles::mode_decompose;

        // The first send target
        const auto& firstMap = faceSendMap[firstProci];

        // Test for broadcast only (ie, master sends identical content).
        // Compare the remaining targets

        bool isBroadcast = true;

        for (label proci = firstProci+1; proci < faceSendMap.size(); ++proci)
        {
            if (const auto& map = faceSendMap[proci]; !map.empty())
            {
                isBroadcast = (map == firstMap);
                if (!isBroadcast)
                {
                    break;
                }
            }
        }

        // Not broadcasting - no further analysis
        if (!isBroadcast)
        {
            break;
        }

        // Appears to be broadcasting
        // --------------------------
        sendStyle = sendStyles::mode_broadcast;

        // Check if the send content itself is simply an identity
        // relationship (eg, for NODEMASTER mode)

        if (firstMap.size() == surf.size())
        {
            bool isIdentity = true;

            forAll(firstMap, i)
            {
                isIdentity = (firstMap[i] == i);
                if (!isIdentity)
                {
                    break;
                }
            }

            // An identity broadcast (no subsetting)
            if (isIdentity)
            {
                sendStyle = sendStyles::mode_identity;
            }
        }
    } while (false);

    int value = static_cast<int>(sendStyle);

    // Let other ranks know what was analyzed
    Pstream::broadcast(value, communicator);

    labelList sendSizes;

    if (debug)
    {
        sendSizes.resize_fill(numProc, 0);

        forAll(faceSendMap, proci)
        {
            sendSizes[proci] = faceSendMap[proci].size();
        }
    }

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::analyzeSendMaps :"
            << " send style:" << int(sendStyle) << endl;

        Pout<< "send sizes: "; sendSizes.writeList(Pout) << endl;
    }

    // Extra debugging
    if (debug & 2)
    {
        labelList recvSizes(numProc, Foam::zero{});

        UPstream::allToAll(sendSizes, recvSizes, communicator);

        Pout<< "recv sizes: "; recvSizes.writeList(Pout) << endl;
    }

    // Deeper analysis of the decompose send maps
    if
    (
        (sendStyle == sendStyles::mode_decompose)
     && UPstream::master(communicator)
     && debug
    )
    {
        // Find which send maps are possibly identical to each other
        Map<labelHashSet> candidates;

        // Grouped according to length
        forAll(sendSizes, proci)
        {
            if (const auto count = sendSizes[proci]; count > 0)
            {
                candidates(count).insert(proci);
            }
        }

        // Only consider anything with two or more identical lengths
        Map<DynamicList<labelHashSet>> matches;

        forAllConstIters(candidates, iter)
        {
            const auto key = iter.key();    // size of the that map
            const auto& procs = iter.val();

            if (procs.size() < 2)
            {
                continue;
            }

            // From same sized send map, find equal contents
            const labelList procList(procs.sortedToc());
            bitSet assigned(procList.size());

            forAll(procList, i)
            {
                if (assigned.test(i))
                {
                    continue;
                }

                const auto& firstMap = faceSendMap[procList[i]];

                labelHashSet group;

                for (label j = i+1; j < procList.size(); ++j)
                {
                    const auto& currMap = faceSendMap[procList[j]];

                    if (!assigned.test(j) && firstMap == currMap)
                    {
                        assigned.set(j);
                        group.insert(procList[j]);
                    }
                }

                // Only keep groups with at least two identical maps
                if (!group.empty())
                {
                    group.insert(procList[i]);
                    matches(key).push_back(std::move(group));
                }
            }
        }

        Pout<< "identical send maps (groups of procs): " << matches << endl;
    }

    return sendStyles(value);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::decompositionMethod&
Foam::distributedTriSurfaceMesh::decomposer() const
{
    // Expected target number of domains
    const auto numProc = UPstream::nProcs(UPstream::worldComm);


    // Set numProc and master nodes (if any)
    const auto assignTreeDecomp = [this,numProc](treeDecomp* ptr) -> bool
    {
        if (ptr)
        {
            ptr->nDomains(numProc);

            if (masterProcIDs_.empty())
            {
                ptr->clearWeights();
            }
            else
            {
                bitSet selection;
                for (auto i : masterProcIDs_)
                {
                    selection.set(i);
                }
                ptr->setWeights(selection);
            }
        }
        return bool(ptr);
    };


    // Special hardcoded handling for tree decomposition,
    // which operates with a minimum of external inputs.
    if
    (
        (distType_ == distributionType::FULLTREE)
     || (distType_ == distributionType::SPARSETREE)
        // (distType_ == distributionType::NODETREE)
    )
    {
        // Discard existing decomposer type if not already treeDecomp
        if (decomposer_ && !dynamic_cast<treeDecomp*>(decomposer_.get()))
        {
            // Wrong decomposer type
            if (debug)
            {
                Pout<< "Had " << decomposer_->type()
                    << " decomposer instead of tree, replacing it."
                    << endl;
            }

            decomposer_.reset(nullptr);
            decomposeParDict_.reset(nullptr);
        }

        // Create tree decomposer
        if (!decomposer_)
        {
            decomposer_.reset(new treeDecomp(numProc));
        }

        // Set numProc and master nodes (if any)
        assignTreeDecomp(dynamic_cast<treeDecomp*>(decomposer_.get()));

        if (debug)
        {
            Pout<< "Setting " << decomposer_->type()
                << " decomposer (distributionType:"
                << distributionTypeNames_[distType_]
                << "), " << masterProcIDs_.size()
                <<  " processor weights" << endl;
        }
    }

    if (decomposer_)
    {
        // Done
        return *decomposer_;
    }


    // Check if local dictionary contains decomposition entries
    if
    (
        word methodName;
        dict_.readIfPresent("method", methodName, keyType::LITERAL)
    )
    {
        if (debug)
        {
            Pout<< "distributedTriSurfaceMesh::decomposer() :"
                << " constructing decomposer from provided dictionary "
                << dict_ << endl;
        }

        decomposeParDict_.reset(nullptr);

        if (methodName == treeDecomp::typeName)
        {
            decomposer_.reset(new treeDecomp(numProc));

            // Set numProc and master nodes (if any)
            if
            (
                assignTreeDecomp(dynamic_cast<treeDecomp*>(decomposer_.get()))
            )
            {
                if (debug)
                {
                    Pout<< "Setting tree-decomposer (from method name)" << nl;
                }
                return *decomposer_;
            }
        }
        else
        {
            // Construct from local/merged dictionary content
            decomposer_ = decompositionMethod::New(dict_);
        }
    }


    // Nothing local, get from decomposeParDict. Cannot use
    // decompositionModel here since we've only got Time and not a mesh.

    if (!decomposer_)
    {
        const auto* dictPtr =
            searchableSurface::time().findObject<IOdictionary>
            (
                // == decompositionModel::canonicalName
                "decomposeParDict"
            );

        if (dictPtr)
        {
            if (debug)
            {
                Pout<< "distributedTriSurfaceMesh::decomposer() :"
                    << " constructing decomposer from registered"
                    << " dictionary " << *dictPtr << endl;
            }
        }
        else
        {
            if (!decomposeParDict_)
            {
                decomposeParDict_.reset
                (
                    new IOdictionary
                    (
                        IOobject
                        (
                            // == decompositionModel::canonicalName
                            "decomposeParDict",
                            searchableSurface::time().system(),
                            searchableSurface::time(),
                            IOobjectOption::MUST_READ,
                            IOobjectOption::NO_WRITE
                        )
                    )
                );

                if (debug)
                {
                    Pout<< "distributedTriSurfaceMesh::decomposer() :"
                        << " reading "
                        << decomposeParDict_().objectRelPath() << endl;
                }
            }

            dictPtr = decomposeParDict_.get();
        }

        decomposer_ = decompositionMethod::New(*dictPtr);
    }


    // Check sanity
    if (UPstream::master() && decomposer_->nDomains() != numProc)
    {
        FatalIOErrorInFunction(dict_)
            << "Number of domains provided in dictionary " << dict_
            << " does not correspond to the number of domains for the mesh"
            << exit(FatalIOError);
    }

    return decomposer_();
}


// Find bounding boxes that guarantee a more or less uniform distribution
// of triangles. Decomposition in here is only used to get the bounding
// boxes, actual decomposition is done later on.
// Returns a per processor a list of bounding boxes that most accurately
// describe the shape. For now just a single bounding box per processor but
// optimisation might be to determine a better fitting shape.
void Foam::distributedTriSurfaceMesh::calcDistribution
(
    const triSurface& s,
    // per local triangle the destination processor
    labelList& distribution
) const
{
    addProfiling
    (
        calc_distribution,
        "distributedTriSurfaceMesh::calcDistribution"
    );

    // Purge old information
    distribution.clear();

    const globalIndex& triIndexer = globalTris();

    // Pure decompose mode if everything currently only resides on the master
    bool pureDecompose = (triIndexer.localSize(0) == triIndexer.totalSize());

    if (pureDecompose)
    {
        if (debug)
        {
            Pout<< "distributedTriSurfaceMesh::calcDistribution :"
                << " determining master-only decomposition for " << s.size()
                << " centroids for " << searchableSurface::name() << endl;
        }

        // Triangle centres
        tmp<pointField> tfaceCentres;
        if (s.hasFaceCentres())
        {
            tfaceCentres.cref(s.faceCentres());
        }
        else
        {
            tfaceCentres = tmp<pointField>::New();
            triSurface::centres(tfaceCentres.ref());
        }
        const auto& triCentres = tfaceCentres();

        if (decomposer().geometric())
        {
            // Purely geometric method (including treeDecomp)
            if (UPstream::master())
            {
                const bool oldParRun = UPstream::parRun(false);
                distribution = decomposer().decompose(triCentres);
                UPstream::parRun(oldParRun);  // Restore parallel state
            }
        }
        else
        {
            // Some other (non-geometric) method, eg scotch etc.
            // Use connectivity

            if (UPstream::master())
            {
                // For very large surfaces, the connectivity array may be
                // too large for compact storage

                // FUTURE: with/without compact storage
                // as CompactListList : if (s.size() >= labelMax/3)

                labelListList pointFaces;
                calculate_pointFaces(s.points().size(), s, pointFaces);

                labelListList faceFaces;
                calculate_faceFaces(s, pointFaces, faceFaces);

                const bool oldParRun = UPstream::parRun(false);
                distribution = decomposer().decompose(faceFaces, triCentres);
                UPstream::parRun(oldParRun);  // Restore parallel state
            }
        }
    }
    else if (distType_ == distributionType::DISTRIBUTED)
    {
        if (debug)
        {
            Pout<< "distributedTriSurfaceMesh::calcDistribution :"
                << " distributed decomposition" << endl;
        }

        // Fully distributed decomposition.
        // Does not filter duplicate triangles.
        if (!decomposer().parallelAware())
        {
            FatalErrorInFunction
                << "The decomposition method " << decomposer().typeName
                << " does not decompose in parallel."
                << " Please choose one that does." << exit(FatalError);
        }

        if (debug)
        {
            Pout<< "distributedTriSurfaceMesh::calcDistribution :"
                << " determining decomposition for " << s.size()
                << " centroids of surface " << searchableSurface::name()
                << endl;
        }

        // Triangle centres
        tmp<pointField> tfaceCentres;
        if (s.hasFaceCentres())
        {
            tfaceCentres.cref(s.faceCentres());
        }
        else
        {
            tfaceCentres = tmp<pointField>::New();
            triSurface::centres(tfaceCentres.ref());
        }
        const auto& triCentres = tfaceCentres();

        // Updated
        distribution = decomposer().decompose(triCentres);
    }
    else
    {
        // Master-only decomposition.
        // Filters duplicate triangles so repeatable.

        if (debug)
        {
            Pout<< "distributedTriSurfaceMesh::calcDistribution :"
                << " collecting all centroids for surface "
                << searchableSurface::name() << endl;
        }

        // Collect all triangle centres
        pointField allCentres(s.size());
        {
            if (s.hasFaceCentres())
            {
                allCentres = s.faceCentres();
            }
            else
            {
                allCentres.resize(s.size());

                forAll(s, trii)
                {
                    allCentres[trii] = s[trii].centre(s.points());
                }
            }

            globalTris().gatherInplace(allCentres);
        }

        // Determine local decomposition
        {
            labelList allDistribution;
            if (UPstream::master())
            {
                labelList allToMerged;
                const label nUnique = mergePoints
                (
                    allCentres,
                    mergeDist_,
                    false,          // verbose = false
                    allToMerged
                );
                const bool hadDuplicates = (nUnique < allCentres.size());

                if (debug)
                {
                    Pout<< "distributedTriSurfaceMesh::calcDistribution :"
                        << " surface:" << searchableSurface::name()
                        << " merged " << allCentres.size()
                        << " centroids down to " << nUnique << endl;
                }

                // The decomposition will be based on the merged points:
                pointField mergedPoints;

                if (hadDuplicates)
                {
                    // Collect merged points
                    mergedPoints =
                        UIndirectList<point>(allCentres, allToMerged);

                    allCentres.clear();  // done with the points
                }
                else
                {
                    // Can recycle the original points
                    mergedPoints = std::move(allCentres);
                }

                labelList mergedDistrib;
                {
                    // Decompose merged centres
                    const bool oldParRun = UPstream::parRun(false);
                    mergedDistrib = decomposer().decompose(mergedPoints);
                    UPstream::parRun(oldParRun);  // Restore parallel state
                }

                // Convert to all
                if (hadDuplicates)
                {
                    allDistribution =
                        UIndirectList<label>(mergedDistrib, allToMerged);
                }
                else
                {
                    allDistribution = std::move(mergedDistrib);
                }
            }

            // Scatter back to processors
            distribution.clear();
            globalTris().scatter(allDistribution, distribution);
        }
    }
}


void Foam::distributedTriSurfaceMesh::markMixedOverlap
(
    const indexedOctree<treeDataTriSurface>& tree,
    PackedList<2>& nodeTypes
) const
{
    const int communicator = UPstream::worldComm;

    if (!UPstream::is_parallel(communicator))
    {
        // not used/useful for non-parallel
        return;
    }

    const int baseMessageTag = UPstream::msgType() + 10;

    const triSurface& surf = *this;

    // At this point each triangle is on a unique processor so we can just
    // go and check what additional triangles would be sent over based on the
    // bounding boxes.
    // Note that the surface at this point is already uniquely distributed
    // (every triangle is on a single processor only) and the tree is only
    // for those.


    // polling receive and process
    if
    (
        UPstream::is_parallel(communicator)
     && (control_ & controlType::ctrlPollingOverlaps)
    )
    {
        // Sender: find/extract overlapping faces

        bitSet faceSelect(surf.size());

        // Extraction (for sending)
        const auto extractOverlap = [&](label proci) -> autoPtr<triSurface>
        {
            if
            (
                triSurfaceTools::findOverlapFaces
                (
                    surf,
                    procBb_[proci],
                    faceSelect
                )
            )
            {
                return autoPtr<triSurface>::New(surf.subsetMesh(faceSelect));
            }
            return nullptr;
        };

        // Receiver:
        // Set all tree nodes overlapping buffer triangles to mixed.
        // This state will be preserved when caching the volume type
        // (see below). This uses the actual local triangles.

        const auto processOverlap = [&](label proci, triSurface& s)
        {
            // node=0 (root node)
            markMixed(tree, 0, s, nodeTypes);
            s.clear();  // discard geometry after processing
        };

        exchangeSurfaces
        (
            communicator,
            receiveOrder::unsequenced,  // process in any order
            baseMessageTag,

            extractOverlap,   // Sender
            processOverlap    // Receiver

            // No pre-receive callback
        );
    }

    // non-polling : PstreamBuffers
    if
    (
        UPstream::is_parallel(communicator)
     && !(control_ & controlType::ctrlPollingOverlaps)
    )
    {
        // Sender: find/extract overlapping faces

        const auto myProci = UPstream::myProcNo(communicator);

        bitSet faceSelect(surf.size());

        PstreamBuffers pBufs(communicator, baseMessageTag);

        // Test & send
        forAll(procBb_, proci)
        {
            if
            (
                (proci != myProci) &&
                triSurfaceTools::findOverlapFaces
                (
                    surf,
                    procBb_[proci],
                    faceSelect
                )
            )
            {
                triSurface subSurface(surf.subsetMesh(faceSelect));

                UOPstream os(proci, pBufs);
                os << subSurface;
            }
        }

        pBufs.finishedSends();

        // Receiver:
        // Set all tree nodes overlapping buffer triangles to mixed.
        // This state will be preserved when caching the volume type
        // (see below). This uses the actual local triangles.

        for (const int proci : pBufs.allProcs())
        {
            if (proci != myProci && pBufs.recvDataCount(proci))
            {
                UIPstream is(proci, pBufs);
                triSurface subSurface(is);

                // node=0 (root node)
                markMixed(tree, 0, subSurface, nodeTypes);
            }
        }
    }
}


// Merge into allTris, allPoints
void Foam::distributedTriSurfaceMesh::merge
(
    const scalar mergeDist,
    const UList<labelledTri>& subTris,
    const UList<point>& subPoints,

    DynamicList<labelledTri>& allTris,
    DynamicList<point>& allPoints,

    labelList& faceConstructMap,
    labelList& pointConstructMap
)
{
    matchPoints
    (
        subPoints,
        allPoints,
        mergeDist,  // match distance (uniform)
        false,      // verbose = false
        pointConstructMap
    );


    // Add all unmatched points
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    const label nOldAllPoints = allPoints.size();
    const label nOldAllFaces = allTris.size();

    label allPointi = nOldAllPoints;
    forAll(pointConstructMap, pointi)
    {
        if (pointConstructMap[pointi] == -1)
        {
            pointConstructMap[pointi] = allPointi++;
        }
    }

    if (allPointi > nOldAllPoints)
    {
        allPoints.resize(allPointi);

        forAll(pointConstructMap, pointi)
        {
            if (pointConstructMap[pointi] >= nOldAllPoints)
            {
                allPoints[pointConstructMap[pointi]] = subPoints[pointi];
            }
        }
    }


    // Add all unmatched triangles
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label allTrii = allTris.size();
    allTris.resize(allTris.size() + subTris.size());

    faceConstructMap.resize_nocopy(subTris.size());

    // The original faces
    const auto oldTris = allTris.slice(0, nOldAllFaces);

    // Test if remapped subTri is a candidate for matching within old
    // allTris. Only possible when all of its points map within
    // nOldAllPoints.
    const auto eligible_oldTri = [nOldAllPoints](const triFace& tri)
    {
        return (tri.max() < nOldAllPoints);
    };


    // On the reasonable assumption that there will be relatively few
    // triangles that match the existing ones (eg, overlapping faces),
    // we pre-scan to see which of the subTris are candidates.
    //
    // If the remapped subTri has any vertex outside of the original
    // points (ie, >= nOldAllPoints) it cannot match.

    bitSet eligiblePoints(nOldAllPoints);

    forAll(subTris, trii)
    {
        // Using triFace (not labelledTri) is enough for this check
        const triFace mappedTri(pointConstructMap, subTris[trii]);

        // Candidate only if all points map to old points
        if (eligible_oldTri(mappedTri))
        {
            eligiblePoints.set(mappedTri.min());
        }
    }

    labelListList oldPointFaces;
    calculate_pointFaces
    (
        nOldAllPoints, allTris, oldPointFaces, eligiblePoints
    );

    forAll(subTris, trii)
    {
        const auto& subTri = subTris[trii];

        // Get triangle in new numbering
        const labelledTri mappedTri(pointConstructMap, subTri);

        label matchedTri = -1;

        if (eligible_oldTri(mappedTri))
        {
            // All three points map within old points.
            // See if same triangle.

            // Only check the original faces!
            matchedTri = findTriangle(oldTris, oldPointFaces, mappedTri);
        }

        if (matchedTri >= 0)
        {
            faceConstructMap[trii] = matchedTri;
        }
        else
        {
            // Add
            faceConstructMap[trii] = allTrii;
            allTris[allTrii] = mappedTri;
            allTrii++;
        }
    }

    // Resize to the current total used
    allTris.resize(allTrii);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::distributedTriSurfaceMesh::distributeGeometry
(
    labelListList&& faceSendMap,
    autoPtr<mapDistribute>& faceMap,
    autoPtr<mapDistribute>& pointMap,
    List<labelledTri>& allTris,
    pointField& allPoints
) const
{
    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);
    const auto numProc = UPstream::nProcs(communicator);

    // The base message tag for communications
    const int baseMessageTag = UPstream::msgType() + 10;

    const triSurface& surf = *this;

    // For mapDistribute
    labelListList pointSendMap(numProc);

    labelListList faceConstructMap(numProc);
    labelListList pointConstructMap(numProc);


    // Create a subset surface from faceSendMap
    // and imprint additional information onto pointSendMap
    const auto create_submesh =
        [&](label proci) -> triSurface
        {
            return surf.subSampleMesh(pointSendMap[proci], faceSendMap[proci]);
        };


    // What type of exchange is used?
    const auto sendStyle = analyzeSendMaps(surf, faceSendMap, communicator);


    // Exchange surfaces
    // ~~~~~~~~~~~~~~~~~

    // We are mostly fortunate and will be decomposing,
    // which means that direct, non-blocking transmission of the surfaces
    // is possible.

    if
    (
        (
            sendStyle == sendStyles::mode_decompose ||
            sendStyle == sendStyles::mode_broadcast ||
            sendStyle == sendStyles::mode_identity
        )
    )
    {
        // Storage for non-blocking sends
        PtrList<triSurface> sendBuffers;

        if (debug)
        {
            Info<< "distributedTriSurfaceMesh::distributeGeometry :"
                << " decomposing (style:" << int(sendStyle) << ")" << endl;
        }

        const auto startOfRequests = UPstream::nRequests();

        // The (nPoints, nFaces) per destination rank.
        // accumulated on master, scattered to sub-ranks
        List<labelPair> sendSizes;

        if (UPstream::master(communicator))
        {
            sendBuffers.resize(numProc);
            sendSizes.resize_fill(numProc, {0, 0});

            // Since "broadcast" mode sends the same content to everyone,
            // create a single subset surface once and just keep sending it
            // to different receivers.
            //
            // In theory we could further optimize when the faceMap is an
            // identity (avoid sub-setting entirely) but it may not have
            // compact points etc so don't give that special treatment for
            // now.

            // The first non-empty send map when broadcasting, otherwise -1.
            const int bcastProci =
            (
                (
                    sendStyle == sendStyles::mode_broadcast
                 || sendStyle == sendStyles::mode_identity
                )
              ? static_cast<int>
                (
                    ListOps::find_if
                    (
                        faceSendMap,
                        [](const auto& map) { return (map.size() > 0); }
                    )
                )
              : -1
            );

            if (int proci = bcastProci; proci >= 0)
            {
                auto& sub =
                    sendBuffers.emplace_set(proci, create_submesh(proci));

                // Sizes not yet needed elsewhere, but cheap to query
                // (numPoints, numFaces)
                sendSizes[proci] = { sub.points().size(), sub.size() };

                if (debug)
                {
                    Info<< "Using broadcast surface=" << proci
                        << " (numPoints numFaces) : "
                        << sendSizes[proci] << nl;
                }
            }

            for (int proci : UPstream::subProcs(communicator))
            {
                if (faceSendMap[proci].size() > 0)
                {
                    // Send a broadcasting surface if possible,
                    // or generate another sub-surface.

                    const auto* sendptr = sendBuffers.get(bcastProci);

                    if (sendptr)
                    {
                        // May need to update the send maps.
                        // - faceMap should be OK (was used to establish
                        //   the send style).
                        // - pointMap needs to be filled.

                        if (proci != bcastProci)
                        {
                            pointSendMap[proci] = pointSendMap[bcastProci];
                        }
                    }
                    else
                    {
                        sendBuffers.emplace_set(proci, create_submesh(proci));
                        sendptr = sendBuffers.get(proci);
                    }
                    const auto& sub = *sendptr;

                    // (numPoints, numFaces)
                    sendSizes[proci] = { sub.points().size(), sub.size() };

                    // Send surface primitives (non-blocking)
                    sendSurfaceGeometry
                    (
                        proci,
                        sub,
                        baseMessageTag,
                        communicator
                    );
                }
            }
        }

        // Scatter the meta-data (sizing)
        const labelPair receiving =
            UPstream::listScatterValues(sendSizes, communicator);


        // On master:
        // - mode_decompose : extract local sample
        // - mode_broadcast : wait for all sends to complete and recover
        //   the broadcast geometry from sendBuffers[0]

        if (UPstream::master(communicator))
        {
            if (faceSendMap[0].size() > 0)
            {
                // If broadcasting, then sendBuffers[0] will be populated.
                // This means we also need to wait for that data to
                // finish being sent before we can use it ourselves.

                if (auto* surfptr = sendBuffers.get(0))
                {
                    // Wait for broadcasting to complete
                    UPstream::waitRequests(startOfRequests);

                    auto& sub = *surfptr;

                    allTris = std::move(sub.storedFaces());
                    allPoints = std::move(sub.storedPoints());
                }
                else
                {
                    // No sendBuffers[0] used, so no surface to recover
                    // so just subset now.

                    triSurface sub(create_submesh(0));

                    allTris = std::move(sub.storedFaces());
                    allPoints = std::move(sub.storedPoints());
                }
            }
            else
            {
                allTris.clear();
                allPoints.clear();
            }
        }

        // On sub-ranks:
        // - use the numFaces/numPoints to setup the receives
        else
        {
            // sub-ranks
            const auto& [numPoints, numFaces] = receiving;

            if (numFaces > 0)
            {
                // Receive surface primitives (non-blocking)
                // directly into allPoints, allTris

                allTris.resize_nocopy(numFaces);
                allPoints.resize_nocopy(numPoints);

                recvSurfaceGeometry
                (
                    UPstream::masterNo(),
                    allPoints,
                    allTris,
                    baseMessageTag,
                    communicator
                );
            }
            else
            {
                allTris.clear();
                allPoints.clear();
            }
        }

        // Final wait
        UPstream::waitRequests(startOfRequests);

        // All ranks have populated allTris, allPoints (above)

        // The corresponding construct maps are simple.
        // - everyone constructs from master only, and the data
        //   has already arrived in (allTris, allPoints)

        {
            faceConstructMap[0] = Foam::identity(allTris.size());
            pointConstructMap[0] = Foam::identity(allPoints.size());
        }

    }
    else  // (sendStyles::mode_general)
    {
        // General communication
        // ~~~~~~~~~~~~~~~~~~~~~

        // Accumulated surface geometry (resized as needed)
        DynamicList<labelledTri> dynTris;
        DynamicList<point> dynPoints;

        // My own surface first
        {
            triSurface subSurface(create_submesh(myProci));

            dynTris = std::move(subSurface.storedFaces());
            dynPoints = std::move(subSurface.storedPoints());

            faceConstructMap[myProci] = identity(faceSendMap[myProci].size());
            pointConstructMap[myProci] = identity(pointSendMap[myProci].size());
        }


        if (UPstream::is_parallel(communicator))
        {
            // Sender:
            // create a send surface for an non-empty faceSendMap.
            // Uses the local create_submesh() helper, which also imprints
            // the pointSendMap.

            const auto subsetAndSend = [&](label proci) -> autoPtr<triSurface>
            {
                if (proci != myProci && faceSendMap[proci].size() > 0)
                {
                    return autoPtr<triSurface>::New(create_submesh(proci));
                }
                return nullptr;
            };

            // Receiver:
            // Merge the received surface in sequence, update the construct maps

            const auto mergeReceived = [&](label proci, triSurface& sub)
            {
                merge
                (
                    mergeDist_,
                    sub,
                    sub.points(),
                    dynTris,
                    dynPoints,
                    faceConstructMap[proci],
                    pointConstructMap[proci]
                );
                sub.clear();  // discard geometry after processing
            };


            // Allocation/resizing before receiving
            // Upper bound for the merged sizes (no de-dup)

            const auto preAllocator = [&](const UList<labelPair>& recvSizes)
            {
                auto nPointsMax = dynPoints.size();
                auto nFacesMax  = dynTris.size();

                for (const auto& [np, nf] : recvSizes)
                {
                    nPointsMax += np;
                    nFacesMax  += nf;
                }
                dynPoints.reserve_exact(nPointsMax);
                dynTris.reserve_exact(nFacesMax);
            };


            exchangeSurfaces
            (
                communicator,
                receiveOrder::sequenced,  // process in rank order
                baseMessageTag,

                subsetAndSend,  // Sender
                mergeReceived,  // Receiver

                preAllocator    // pre-allocate before receives
            );
        }

        allTris.transfer(dynTris);
        allPoints.transfer(dynPoints);

    }   // (sendStyles)


    // -------------------------------------------------
    // Done all communication - define the mapDistribute
    // -------------------------------------------------

    faceMap.emplace
    (
        allTris.size(),
        std::move(faceSendMap),
        std::move(faceConstructMap)
    );

    pointMap.emplace
    (
        allPoints.size(),
        std::move(pointSendMap),
        std::move(pointConstructMap)
    );
}


void Foam::distributedTriSurfaceMesh::distributeGeometry_buffered
(
    labelListList&& faceSendMap,
    autoPtr<mapDistribute>& faceMap,
    autoPtr<mapDistribute>& pointMap,
    List<labelledTri>& allTris,
    pointField& allPoints
) const
{
    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);
    const auto numProc = UPstream::nProcs(communicator);

    // The base message tag for communications
    const int baseMessageTag = UPstream::msgType() + 10;

    const triSurface& surf = *this;

    // For mapDistribute
    labelListList pointSendMap(numProc);

    labelListList faceConstructMap(numProc);
    labelListList pointConstructMap(numProc);


    // Create a subset surface from faceSendMap
    // and imprint additional information onto pointSendMap
    const auto create_submesh =
        [&](label proci) -> triSurface
        {
            return surf.subSampleMesh(pointSendMap[proci], faceSendMap[proci]);
        };


    // What type of exchange is used?
    const auto sendStyle = analyzeSendMaps(surf, faceSendMap, communicator);


    // Exchange surfaces
    // ~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Info<< "distributedTriSurfaceMesh::distributeGeometry_buffered :"
            << " decomposing (style:" << int(sendStyle) << ")" << endl;
    }


    // General communication
    // - with PstreamBuffers (standard for OpenFOAM-2606 and earlier)

    // My own surface first
    // ~~~~~~~~~~~~~~~~~~~~
    DynamicList<labelledTri> dynTris;
    DynamicList<point> dynPoints;

    {
        triSurface subSurface(create_submesh(myProci));

        dynTris = std::move(subSurface.storedFaces());
        dynPoints = std::move(subSurface.storedPoints());

        faceConstructMap[myProci] = identity(faceSendMap[myProci].size());
        pointConstructMap[myProci] = identity(pointSendMap[myProci].size());
    }


    // Send all
    // ~~~~~~~~

    PstreamBuffers pBufs(communicator, baseMessageTag);

    // Meta-data: (nPoints, nFaces), zero if nothing is sent
    List<labelPair> sendSizes;

    constexpr bool with_presizing = false;

    if constexpr (with_presizing)
    {
        sendSizes.resize_fill(numProc, labelPair(0, 0));
    }

    forAll(faceSendMap, proci)
    {
        if (proci != myProci && (faceSendMap[proci].size() > 0))
        {
            UOPstream os(proci, pBufs);

            triSurface sub(create_submesh(proci));

            if constexpr (with_presizing)
            {
                sendSizes[proci] = { sub.points().size(), sub.size() };
            }

            os << sub;
        }
    }

    pBufs.finishedSends();

    if constexpr (with_presizing)
    {
        List<labelPair> recvSizes(numProc, labelPair(0, 0));

        UPstream::allToAll(sendSizes, recvSizes, communicator);

        auto nPointsMax = dynPoints.size();
        auto nFacesMax  = dynTris.size();

        for (const auto& [np, nf] : recvSizes)
        {
            nPointsMax += np;
            nFacesMax  += nf;
        }
        dynPoints.reserve_exact(nPointsMax);
        dynTris.reserve_exact(nFacesMax);
    };


    // Receive and merge all
    // ~~~~~~~~~~~~~~~~~~~~~

    for (int proci : pBufs.allProcs())
    {
        if (proci != myProci && (pBufs.recvDataCount(proci) > 0))
        {
            UIPstream is(proci, pBufs);
            triSurface subSurface(is);

            // Merge the received surface in sequence
            merge
            (
                mergeDist_,
                subSurface,
                subSurface.points(),

                dynTris,
                dynPoints,
                faceConstructMap[proci],
                pointConstructMap[proci]
            );
        }
    }

    allTris.transfer(dynTris);
    allPoints.transfer(dynPoints);


    // -------------------------------------------------
    // Done all communication - define the mapDistribute
    // -------------------------------------------------

    faceMap.emplace
    (
        allTris.size(),
        std::move(faceSendMap),
        std::move(faceConstructMap)
    );

    pointMap.emplace
    (
        allPoints.size(),
        std::move(pointSendMap),
        std::move(pointConstructMap)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
