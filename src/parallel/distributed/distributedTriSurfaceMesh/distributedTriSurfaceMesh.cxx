/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2025 OpenCFD Ltd.
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

Notes
    Hidden, undocumented developer switches - may disappear or change at
    any time.

    - distributed.direct-exchange (compile-time: on) :
      Use direct binary point-to-point geometry sending.

      When off/disabled : use PstreamBuffers to serialize and exchange
      (as per 2606 and earlier).

    - distributed.overlap-nodeMaster (compile-time: on)
      Use double buffering and local buffers to overlap communication
      and calculation for NODEMASTER queries.

      When off/disabled : use globalIndex to gather all input parameters
      sample everything on master and then scatter back
      (as per 2606 and earlier).

    - distributed.overlap-polling (compile-time: on)
      Add any distributed ragged surface edges into the octree bookkeeping
      with sampling and non-blocking sends. The receiving side
      polls the waiting request and processes whenever something arrives.

      When off/disabled : use PstreamBuffers to serialize the sub-surfaces
      and exchange [waits for everybody] and then locally process
      whatever can be reserialized from PstreamBuffers
      (as per 2606 and earlier).

    Hidden, undocumented testing switches.
    Names will either disappear entirely or be renamed and promoted
    to documented user switches.

    - pss.minTriangles :
      If surface loaded on master is below this threshold,
      switch from node-master, distributed etc to "transparent mode",
      which means running with a copy of the surface on every node.
      Makes <distributedTriSurfaceMesh> behave like <triSurfaceMesh> instead.

    - pss.masterTriangles :
      If surface loaded on master is below this threshold,
      switch from distributed to node-master mode
      (as if "nodeMaster" had been specified in the dictionary).

\*---------------------------------------------------------------------------*/

#include "distributedTriSurfaceMesh.H"
#include "mapDistribute.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"
#include "triangle.H"
#include "triSurfaceTools.H"
#include "matchPoints.H"
#include "globalIndex.H"
#include "Time.H"
#include "DynamicField.H"

#include "IFstream.H"
#include "vectorList.H"
#include "bitSet.H"
#include "PatchTools.H"
#include "OBJstream.H"
#include "labelBits.H"
#include "orientedSurface.H"
#include "decompositionMethod.H"
#include "profiling.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(distributedTriSurfaceMesh, 0);
    addToRunTimeSelectionTable
    (
        searchableSurface,
        distributedTriSurfaceMesh,
        dict
    );

} // End namespace Foam


const Foam::Enum
<
    Foam::distributedTriSurfaceMesh::distributionType
>
Foam::distributedTriSurfaceMesh::distributionTypeNames_
({
    { distributionType::FROZEN, "frozen" },
    { distributionType::FOLLOW, "follow" },
    { distributionType::TRANSPARENT, "transparent" },
    { distributionType::INDEPENDENT, "independent" },
    { distributionType::DISTRIBUTED, "distributed" },
    { distributionType::NODEMASTER, "nodeMaster" },
    { distributionType::FULLTREE, "fullTree" },
    { distributionType::SPARSETREE, "sparseTree" },
    // FUTURE { distributionType::NODETREE, "nodeTree" },
});


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Quick debug control
void check_debug(const Foam::dictionary& dict)
{
    int flag = 0;

    // triSurfaceMesh::debug
    if
    (
        dict.readIfPresent("debug.triSurfaceMesh", flag)
     && (flag > 0)
    )
    {
        Foam::triSurfaceMesh::debug = flag;
    }

    if
    (
        dict.readIfPresent("debug.distributedTriSurfaceMesh", flag)
     && (flag > 0)
    )
    {
        Foam::distributedTriSurfaceMesh::debug = flag;
    }
}


// Dynamic (pre-reservable) equivalent of an input container:
// Field -> DynamicField, everything else (List, UList, ...) -> DynamicList.
//
// Simpler (more fragile?) to just check
// std::is_base_of_v<Foam::refCount, T>

template<class T>
using dynamic_input_t = std::conditional_t
<
    stdFoam::is_template_base_of_v<Foam::Field, T>,
    Foam::DynamicField<typename T::value_type>,
    Foam::DynamicList<typename T::value_type>
>;


// Node-master query: sub-ranks send their inputs to the (node-local)
// master, which receives, processes and sends back the results one
// sub-rank at a time. Receives are double-buffered so that the next
// sub-rank's data is in flight while the current one is processed,
// and the sends of results are non-blocking.
//
// The output parameter and each input parameter uses a separate
// message tag to esnure precise message matching.
//
// The query functor is called as fn(inputs..., output).

template
<
    class OutputType,
    class Func,
    class... Inputs
>
void nodeMasterQuery
(
    const int communicator,
    const int tag,  // base message tag
    Func&& fn,
    OutputType& output,
    const Inputs&... inputs
)
{
    using namespace Foam;

    // Message tag for results after the inputs
    const int tagResult = tag + 4 + int(sizeof...(Inputs));

    const auto startOfRequests = UPstream::nRequests();

    const auto numProc = UPstream::nProcs(communicator);

    static_assert
    (
        sizeof...(Inputs) >= 1,
        "nodeMasterQuery requires one or more inputs"
    );

    // Contiguous input content is needed for the send/recv
    static_assert
    (
        (is_contiguous_v<typename Inputs::value_type> && ...),
        "nodeMasterQuery requires contiguous input value types"
    );

    // The local size of the first input dictates the sizes of all other
    // input fields and the output field.

    const label localSize =
    (
        std::get<0>(std::forward_as_tuple(inputs...)).size()
    );

    const globalIndex gi
    (
        globalIndex::gatherOnly{},
        localSize,
        communicator
    );

    constexpr bool contiguous_output =
        is_contiguous_v<typename OutputType::value_type>;


    // Persistent storage for master sending to sub-ranks.
    // contiguous output or non-contiguous output
    [[maybe_unused]] PtrList<OutputType> sendBuffers;
    [[maybe_unused]] PtrList<DynamicList<char>> charBuffers;


    if (UPstream::master(communicator))
    {
        // Double buffering when receiving from sub-ranks.
        // - the dynamic_input_t helper maps inputs (Field, List,..) to their
        //   dynamic equivalent. Avoids allocation thrashing in loops.

        std::tuple<dynamic_input_t<Inputs>...> recvBuffers[2];
        // non-dynamic: std::tuple<Inputs...> recvBuffers[2];

        // The max receive size, per double-buffer slot
        label maxRecvSize[2] = {0, 0};

        // The first/last range of recv requests, per double-buffer slot
        int req_first[2] = {0, 0};
        int req_last[2] = {0, 0};

        // Sub-ranks with data to process
        DynamicList<int> procs(numProc);

        for (int proci : UPstream::subProcs(communicator))
        {
            if (label recvSize = gi.localSize(proci); recvSize > 0)
            {
                const int slot = (procs.size() % 2);  // buffer slot (0,1)
                maxRecvSize[slot] = std::max(maxRecvSize[slot], recvSize);

                procs.push_back(proci);
            }
        }

        if constexpr (contiguous_output)
        {
            sendBuffers.resize(numProc);
        }
        else
        {
            charBuffers.resize(numProc);
        }

        // Pre-reserve the receive buffers
        for (int slot = 0; slot < 2; ++slot)
        {
            std::apply
            (
                [&](auto&... buf) { (buf.reserve(maxRecvSize[slot]), ...); },
                recvBuffers[slot]
            );
        }

        // Action for posting a receive for given item in the procs sequence
        auto postReceive = [&](const label seq)
        {
            const int slot = (seq % 2);  // buffer slot (0,1)
            const int proci = procs[seq];
            auto& bufs = recvBuffers[slot];

            req_first[slot] = UPstream::nRequests();

            std::apply
            (
                [&](auto&... buf)
                {
                    int itag = tag;
                    (
                        (
                            buf.resize_nocopy(gi.localSize(proci)),
                            UIPstream::read
                            (
                                UPstream::commsTypes::nonBlocking,
                                proci,
                                buf,
                              ++itag,
                                communicator
                            )
                        ),
                        ...
                    );
                },
                bufs
            );

            req_last[slot] = UPstream::nRequests();
        };

        // Post the initial receive
        if (!procs.empty())
        {
            postReceive(0);
        }

        // Overlap own query with the first receive
        if (localSize > 0)
        {
            fn(inputs..., output);
        }
        else
        {
            output.clear();
        }

        // Process the input sequence
        for (label seq = 0; seq < procs.size(); ++seq)
        {
            const int slot = (seq % 2);  // buffer slot (0,1)
            const int proci = procs[seq];

            // Post the next receive
            if (seq+1 < procs.size())
            {
                postReceive(seq+1);
            }

            // Wait for the previous receive
            UPstream::waitRequests
            (
                req_first[slot],
                (req_last[slot] - req_first[slot])
            );

            // For contiguous, use sendBuffers
            if constexpr (contiguous_output)
            {
                auto& procOutput = sendBuffers.emplace_set(proci);

                std::apply
                (
                    [&](auto&... buf) { fn(buf..., procOutput); },
                    recvBuffers[slot]
                );

                UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    procOutput.data(),
                    procOutput.size(),
                    tagResult,
                    communicator
                );
            }
            else
            {
                // For non-contiguous, use charBuffers for storage
                // and UOPstream to format into that storage.

                auto& serialized = charBuffers.emplace_set(proci);
                {
                    OutputType procOutput;
                    std::apply
                    (
                        [&](auto&... buf) { fn(buf..., procOutput); },
                        recvBuffers[slot]
                    );

                    UOPstream os(serialized);
                    os << procOutput;
                }

                UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    serialized,
                    tagResult,
                    communicator
                );
            }
        }
    }
    else  // sub-rank
    {
        if (localSize > 0)
        {
            // send data from sub-rank to master
            int itag = tag;
            (
                UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    UPstream::masterNo(),
                    inputs,
                  ++itag,
                    communicator
                ),
                ...
            );

            // receive data from master
            if constexpr (contiguous_output)
            {
                output.resize_nocopy(localSize);

                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    UPstream::masterNo(),
                    output.data(),
                    output.size(),
                    tagResult,
                    communicator
                );
            }
            else
            {
                IPstream::recv
                (
                    output,
                    UPstream::masterNo(),
                    tagResult,
                    communicator
                );
            }
        }
        else
        {
            output.clear();
        }
    }

    // Final wait
    UPstream::waitRequests(startOfRequests);
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Find edge-connected face.
//
// Accepts "mini-format" point-faces since it does lookups
// starting with the min face vertex

static label findOtherFace
(
    const triSurface& surf,
    const UList<labelList>& pointFaces,
    label nearFacei,
    label nearLabel
)
{
    if (nearFacei < 0 || nearFacei >= surf.size())
    {
        return -1;
    }
    const auto& f = surf[nearFacei];

    // if (nearLabel < 0 || nearLabel >= f.size());
    // {
    //     return -1;
    // }

    // Search for attached faces by interrogating
    // the point-faces based on edge.
    // Use the edge min vertex in case we have compact point-faces

    const edge e = f.edge(nearLabel);

    const label minVert = e.min();

    for (const auto facei : pointFaces[minVert])
    {
        if (facei != nearFacei)
        {
            if (surf[facei].contains(e))
            {
                return facei;
            }
        }
    }
    return -1;
}

}  // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //- Combine operator for volume types
    struct volumeCombineOp
    {
        void operator()(volumeType& x, const volumeType& y) const
        {
            if (x == volumeType::MIXED || y == volumeType::MIXED)
            {
                FatalErrorInFunction << "Illegal volume type "
                    << volumeType::names[x]
                    << "," << volumeType::names[y] << exit(FatalError);
            }
            else
            {
                switch (x)
                {
                    case volumeType::UNKNOWN :
                    {
                        if (y == volumeType::INSIDE || y == volumeType::OUTSIDE)
                        {
                            x = y;
                        }
                    }
                    break;
                    case volumeType::INSIDE :
                    case volumeType::OUTSIDE :
                    {
                        if (x != y)
                        {
                            FatalErrorInFunction << "Conflicting volume types "
                                << volumeType::names[x] << ","
                                << volumeType::names[y] << exit(FatalError);
                        }
                    }
                    break;
                    case volumeType::MIXED:
                    break;
                }
            }
        }
    };


    ////- Same as above but with some more debugging
    //typedef Tuple2<Pair<point>, volumeType> NearType;
    //
    ////- Combine operator for volume types
    //struct NearTypeCombineOp
    //{
    //    //const auto& this_;
    //
    //    void operator()(NearType& nearX, const NearType& nearY) const
    //    {
    //        auto x = nearX.second();
    //        auto y = nearY.second();
    //
    //        if (x == volumeType::MIXED || y == volumeType::MIXED)
    //        {
    //            FatalErrorInFunction << "Illegal volume type "
    //                << volumeType::names[x]
    //                << "," << volumeType::names[y]
    //                << " nearX:" << nearX << " nearY:" << nearY
    //                << exit(FatalError);
    //        }
    //        else
    //        {
    //            switch (x)
    //            {
    //                case volumeType::UNKNOWN :
    //                {
    //                    if
    //                    (
    //                        y == volumeType::INSIDE
    //                     || y == volumeType::OUTSIDE
    //                    )
    //                    {
    //                        nearX = nearY;
    //                    }
    //                }
    //                break;
    //                case volumeType::INSIDE :
    //                case volumeType::OUTSIDE :
    //                {
    //                    if (x != y)
    //                    {
    //                        FatalErrorInFunction
    //                            << "Conflicting volume types "
    //                            << volumeType::names[x] << ","
    //                            << volumeType::names[y]
    //                            << " nearX:" << nearX << " nearY:" << nearY
    //                            << exit(FatalError);
    //                    }
    //                }
    //                break;
    //                case volumeType::MIXED :
    //                break;
    //            }
    //        }
    //    }
    //};


    //- Combine operator for nearest
    typedef Tuple2<pointIndexHit, scalar> nearestAndDist;

    const nearestAndDist nearestZero(pointIndexHit(), -GREAT);

    struct nearestEqOp
    {
        //- Min distance (MINLOC) operation for nearestAndDist,
        //- that only propagates through 'hit' entries
        //- with distance tracked by the second part of the tuple
        void operator()(nearestAndDist& x, const nearestAndDist& y) const
        {
            if
            (
                y.first().hit() &&
                (x.first().not_hit() || (y.second() < x.second()))
            )
            {
                x = y;
            }
        }
    };

}  // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::distributedTriSurfaceMesh::upgradeDictionary(dictionary& dict)
{
    // bounds
    if (entry* eptr = dict.findEntry("bounds", keyType::LITERAL))
    {
        // v2606 and earlier: was treeBoundBoxList
        // now: treeBoundBox
        treeBoundBox bb;

        enum class actions
        {
            no_action, remove_nonstream, remove_empty, replace_entry
        };

        actions action = actions::no_action;

        if (auto* streamPtr = eptr->streamPtr(); !streamPtr)
        {
            action = actions::remove_nonstream;
        }
        else
        {
            const bool oldThrowingError = FatalError.throwing(true);
            const bool oldThrowingIOerr = FatalIOError.throwing(true);

            try
            {
                ITstream& is = *streamPtr;
                treeBoundBoxList list(is);

                if (list.empty())
                {
                    action = actions::remove_empty;
                }
                else
                {
                    action = actions::replace_entry;
                    bb = list[0];
                }
            }
            catch (...)
            {
                // Read failure. Must be a plain boundBox ?
            }

            FatalError.throwing(oldThrowingError);
            FatalIOError.throwing(oldThrowingIOerr);
        }

        switch (action)
        {
            case actions::remove_nonstream :
            {
                Info<< "Removing non-primitive 'bounds' from dictionary."
                    << endl;
                dict.remove("bounds");
                break;
            }
            case actions::remove_empty :
            {
                Info<< "Removing empty 'bounds' list from dictionary."
                    << endl;
                dict.remove("bounds");
                break;
            }
            case actions::replace_entry :
            {
                Info<< "Replacing 'bounds' (list) with single entry."
                    << endl;
                dict.set("bounds", bb);
                break;
            }
            default : break;  // no_action
        }
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::distributedTriSurfaceMesh::findLocalInstance
(
    const IOobject& io
)
{
    // Modified findInstance which also looks in parent directory

    word instance
    (
        io.time().findInstance
        (
            io.local(),
            word::null,
            IOobjectOption::READ_IF_PRESENT,
            word::null,             // No stop instance
            !UPstream::parRun()     // Fallback to "constant" for non-parallel
        )
    );

    if (instance.size())
    {
        return instance;
    }


    // The rest of this code is only when parRun == true ...

    // Replicate findInstance operation but now on parent directory

    // Search in parent directory
    fileName parentDir =
    (
        io.rootPath()/io.globalCaseName()
       /io.instance()/io.db().dbDir()/io.local()/io.name()
    );

    if (fileHandler().isDir(parentDir))
    {
        return io.instance();
    }


    const scalar startValue = io.time().timeOutputValue();

    instantList ts = io.time().times();

    label instIndex = ts.size()-1;

    // Backward search for first time that is <= startValue
    for (; instIndex >= 0; --instIndex)
    {
        if (ts[instIndex] <= startValue)
        {
            break;
        }
    }

    // Continue searching from here
    for (; instIndex >= 0; --instIndex)
    {
        // Shortcut: if actual directory is the timeName we've already tested it
        if (ts[instIndex].name() == io.instance())
        {
            continue;
        }

        parentDir =
        (
            io.rootPath()/io.globalCaseName()
           /ts[instIndex].name()/io.db().dbDir()/io.local()/io.name()
        );

        if (fileHandler().isDir(parentDir))
        {
            return ts[instIndex].name();
        }
    }

    // times() usually already includes the constant() so would have been
    // checked above. Re-test if
    // - times() is empty. Sometimes this can happen (e.g. decomposePar with
    //   collated)
    // - times()[0] is not constant
    if (!ts.size() || ts[0].name() != io.time().constant())
    {
        // Note. This needs to be a hard-coded constant, rather than the
        // constant function of the time, because the latter points to
        // the case constant directory in parallel cases

        parentDir =
        (
            io.rootPath()/io.globalCaseName()
           /io.time().constant()/io.db().dbDir()/io.local()/io.name()
        );

        if (fileHandler().isDir(parentDir))
        {
            return io.time().constant();
        }
    }

    // FatalErrorInFunction
    //     << "Cannot find directory " << io.local() << " in times " << ts
    //     << exit(FatalError);
    //
    // return word::null;

    // Nothing found in parent?
    // - treat like regular findInstance with "constant" fallback
    return io.time().constant();
}


void Foam::distributedTriSurfaceMesh::setNodeMasters
(
    int nProcessorsPerMaster,
    bool needsCommunicator
)
{
    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);
    const auto numProc = UPstream::nProcs(communicator);

    // Free the old (NODEMASTER) communicator and provide a reasonable
    // fallback

    freeCommunicator();
    //FUTURE: comm_ = UPstream::worldComm;

    masterProcIDs_.clear();

    if (nProcessorsPerMaster < 1)
    {
        // Use the regular (shared memory) node master ids
        masterProcIDs_ = UPstream::procID(UPstream::commInterNode());

        if (needsCommunicator)
        {
            // Use shared memory node
            comm_ = UPstream::commLocalNode();
        }
    }
    else if (nProcessorsPerMaster == 1)
    {
        // Everyone is a master!
        masterProcIDs_.resize(numProc);
        Foam::identity(masterProcIDs_);

        if (needsCommunicator)
        {
            // Only talk to ourselves
            comm_ = UPstream::commSelf();
        }
    }
    else  // (nProcessorsPerMaster > 1)
    {
        masterProcIDs_.clear();
        masterProcIDs_.reserve_exact(1+numProc/nProcessorsPerMaster);

        for (label proci = 0; proci < numProc; proci += nProcessorsPerMaster)
        {
            masterProcIDs_.push_back(proci);
        }

        if (masterProcIDs_.empty())  // being overly cautious
        {
            masterProcIDs_.push_back(0);
        }

        if (needsCommunicator)
        {
            // Colour ranks according to their nearest master node.
            // - starts with colour=1 (or any other positive number)

            const int colour =
            (
                static_cast<int>
                (
                    std::upper_bound
                    (
                        masterProcIDs_.begin(),
                        masterProcIDs_.end(),
                        myProci
                    ) - masterProcIDs_.begin()
                )
            );

            comm_ = UPstream::splitCommunicator(UPstream::worldComm, colour);
        }
    }
}


// side-effects:
// - adjusts procBb_ (from dictionary or from current surface)
// - rebuilds exclusionBb_
// - sets/reset communicator
// - resets masterProcIDs_

bool Foam::distributedTriSurfaceMesh::readSettings(bool undecomposed)
{
    const int communicator = UPstream::worldComm;
    const auto numProc = UPstream::nProcs(communicator);

    // hidden, undocumented developer controls
    // - may disappear or change without warning
    {
        bool specified = false;

        if (bool on; dict_.readIfPresent("distributed.direct-exchange", on))
        {
            specified = true;
            if (on) control_ |= (controlType::ctrlDirectExchange);
            else    control_ &= (~controlType::ctrlDirectExchange);
        }
        if (bool on; dict_.readIfPresent("distributed.overlap-nodeMaster", on))
        {
            specified = true;
            if (on) control_ |= (controlType::ctrlOverlappingNodeMaster);
            else    control_ &= (~controlType::ctrlOverlappingNodeMaster);
        }
        if (bool on; dict_.readIfPresent("distributed.overlap-polling", on))
        {
            specified = true;
            if (on) control_ |= (controlType::ctrlPollingOverlaps);
            else    control_ &= (~controlType::ctrlPollingOverlaps);
        }

        // Report
        if (specified || debug)
        {
            Info<< "[controls]"
            << " surface-exchange:"
            << ((control_ & controlType::ctrlDirectExchange) ?
                "direct" : "buffered")
            << " node-master:"
            << ((control_ & controlType::ctrlOverlappingNodeMaster) ?
                "overlapping" : "buffered")
            << " mixed-surfaces:"
            << ((control_ & controlType::ctrlPollingOverlaps) ?
                "polling" : "buffered")
            << endl;
        }
    }

    // Get bb of all domains.a
    {
        List<treeBoundBox> bbs(procBb_.release_boxes());
        bbs.resize_nocopy(numProc);

        treeBoundBox localBb;

        if (!dict_.readIfPresent("bounds", localBb))
        {
            // Not found. Use local triangles' bounding box
            localBb = triSurfaceMesh::bounds();

            dict_.add("bounds", localBb);
        }

        boundBox::allGather(bbs, localBb, communicator);

        procBb_.setBoxes(std::move(bbs));
        procBb_.rebuild();
    }

    decomposeUsingBbs_ = true;
    doExclusionBb_ = true;

    dict_.readIfPresent("decomposeUsingBbs", decomposeUsingBbs_);
    dict_.readIfPresent("exclusionBb", doExclusionBb_);

    // Wanted distribution type, with most sensible fallback type
    constexpr distributionType fallbackType = distributionType::DISTRIBUTED;
    if
    (
       !distributionTypeNames_.readIfPresent
        (
            "distributionType",
            dict_,
            distType_
        )
    )
    {
        distType_ = fallbackType;
        dict_.add("distributionType", distributionTypeNames_[distType_]);
    }

    if
    (
        (distType_ == distributionType::TRANSPARENT)
     || (distType_ == distributionType::NODEMASTER)
    )
    {
        if (decomposeUsingBbs_)
        {
            IOWarningInFunction(dict_)
                << "Switching off bounding box distribution"
                << " ('decomposeUsingBbs') since not"
                << " decomposing the surface. (distributionType: "
                << distributionTypeNames_[distType_] << ")"
                << endl;
            decomposeUsingBbs_ = false;
            dict_.set("decomposeUsingBbs", false);
        }
        doExclusionBb_ = false;
        dict_.set("exclusionBb", false);
    }

    // Free the old (NODEMASTER) communicator and provide a
    // reasonable fallback.

    freeCommunicator();
    //FUTURE: comm_ = UPstream::worldComm;

    masterProcIDs_.clear();

    // Handle node-master designations and setup communicators etc.
    if
    (
        (distType_ == distributionType::NODEMASTER)
     || (distType_ == distributionType::SPARSETREE)
        //(distType_ == distributionType::NODETREE)
    )
    {
        int nProcessorsPerMaster(0);

        // Need an appropriate communicator for NODEMASTER mode
        bool needsCommunicator(distType_ == distributionType::NODEMASTER);

        if
        (
            dict_.readIfPresent
            (
                "nProcessorsPerMaster",
                nProcessorsPerMaster,
                keyType::LITERAL
            )
         && (nProcessorsPerMaster < 1)
        )
        {
            WarningInFunction
                << "Found nProcessorsPerMaster = " << nProcessorsPerMaster
                << " but it should be >= 1. Ignoring that input."
                << endl;

            // Fallback as 0 [SIC], which we permit as an internal
            // value but not as user input (since it is meaningless)
            nProcessorsPerMaster = 0;
        }

        setNodeMasters(nProcessorsPerMaster, needsCommunicator);
    }

    // Merge distance
    if (!dict_.readIfPresent("mergeDistance", mergeDist_))
    {
        mergeDist_ = SMALL;
        dict_.add("mergeDistance", mergeDist_);
    }

    // Is closed?
    if (int closed; dict_.readIfPresent("closed", closed))
    {
        surfaceClosed_ = closed;
    }
    else
    {
        closed = false;

        if (undecomposed)
        {
            closed = triSurfaceMesh::hasVolumeType();
            Pstream::broadcast(closed, communicator);
        }

        dict_.add("closed", closed);

        surfaceClosed_ = closed;
    }

    // Do we know volume inside/outside?
    if
    (
        volumeType::type vt;
        volumeType::names.readIfPresent("outsideVolumeType", dict_, vt)
    )
    {
        outsideVolType_ = vt;
    }
    else
    {
        if (undecomposed && surfaceClosed_)
        {
            if (UPstream::master(communicator))
            {
                // Determine nearest to do inside/outside testing.
                // Do linear search to avoid building tree on this large
                // undecomposed surface. See e.g. triSurfaceTools::surfaceSide

                const triSurface& surf = *this;
                const auto& points = surf.points();
                const boundBox& localBb = triSurfaceMesh::bounds();
                const point outsidePt(localBb.max()+localBb.centre());

                label nearTrii = -1;
                scalar nearest(GREAT);
                forAll(surf, trii)
                {
                    pointHit info = surf[trii].nearestPoint(outsidePt, points);

                    if (nearest > info.distance())
                    {
                        nearest = info.distance();
                        nearTrii = trii;
                    }
                }

                if (nearTrii != -1)
                {
                    auto side =
                        triSurfaceTools::surfaceSide
                        (
                            surf,
                            outsidePt,
                            nearTrii
                        );

                    outsideVolType_ = triSurfaceTools::sideToVolumeType(side);
                }
            }

            // Send the volumeType information to everyone
            {
                int volType = outsideVolType_;
                Pstream::broadcast(volType, communicator);
                outsideVolType_ = volumeType(volType);
            }
        }
        else
        {
            // Mark as to-be-done
            outsideVolType_ = volumeType::UNKNOWN;
        }

        dict_.add("outsideVolumeType", volumeType::names[outsideVolType_]);
    }

    // Have new procBb. Update (if necessary) exclusionBb
    calcExclusionBb();

    if (debug && UPstream::master(communicator))
    {
        Info<< nl
            << "processor boxes: " << procBb_.info() << nl
            << "exclusion boxes: " << exclusionBb_.info() << endl;
    }

    return true;
}


bool Foam::distributedTriSurfaceMesh::optimizeStrategy()
{
    if (!UPstream::parRun())
    {
        return false;
    }

    const int communicator = UPstream::worldComm;
    // const auto myProci = UPstream::myProcNo(communicator);
    // const auto numProc = UPstream::nProcs(communicator);

    // Check if pure decompose
    const triSurface& s = *this;
    const label numFaces = s.size();

    // Pure decompose mode if everything currently only resides on the master
    const bool pureDecompose =
    (
        returnReduceOr
        (
            UPstream::master(communicator) ? numFaces : !numFaces,
            communicator
        )
    );

    if (!pureDecompose)
    {
        Info<< "Cannot adjust decomposition strategy: already distributed"
            << nl;

        return false;
    }

    // Something that will never conflict with distributionType
    constexpr int invalid_strategy = -1000;

    // Some very Setup with some very FUTURE
    int strategy_change = invalid_strategy;

    int nProcessorsPerMaster(0);

    // Assume we might need this again. Don't need to warn this time
    dict_.readIfPresent
    (
        "nProcessorsPerMaster",
        nProcessorsPerMaster,
        keyType::LITERAL
    );

    // hidden, undocumented developer controls
    // - may disappear or change without warning

    if (UPstream::master(communicator) && (numFaces > 0)) do
    {
        label limit1;

        if (dict_.readIfPresent("pss.minTriangles", limit1))
        {
            // Optimization for small surfaces
            if (numFaces < limit1)
            {
                strategy_change = int(distributionType::TRANSPARENT);

                Info<< "Applying min-size strategy (faces:"
                    << numFaces << " < min:" << limit1
                    << ") -> transparent mode" << nl;
                break;
            }
        }

        if (dict_.readIfPresent("pss.masterTriangles", limit1))
        {
            // Optimization for medium size surfaces
            if (numFaces < limit1)
            {
                strategy_change = int(distributionType::NODEMASTER);

                Info<< "Applying node-master strategy (faces:"
                    << numFaces << " < min:" << limit1
                    << ") -> node-master mode" << nl;
                break;
            }
        }
    }
    while (false);

    UPstream::broadcast(&strategy_change, 1, communicator);

    if
    (
        (strategy_change == invalid_strategy)
     || (strategy_change == int(distributionType::FROZEN))
     || (strategy_change == int(distType_))
    )
    {
        if (debug)
        {
            if (strategy_change == int(distType_))
            {
                Info<< "... already using "
                    << distributionTypeNames_[distType_] << endl;
            }
            else
            {
                Info<< "No change in strategy" << nl;
            }
        }
        return false;
    }

    bool changed = false;

    if
    (
        distType_ == distributionType::TRANSPARENT
     || distType_ == distributionType::FOLLOW
    )
    {
        Info<< "Cannot change from "
            << distributionTypeNames_[distType_]
            << " to a different mode" << endl;
    }
    else
    {
        distributionType newDistType = distributionType(strategy_change);

        if (newDistType == distributionType::TRANSPARENT)
        {
            // Easy to change
            changed = true;
            distType_ = newDistType;

            freeCommunicator();
            masterProcIDs_.clear();
        }
        else if (newDistType == distributionType::NODEMASTER)
        {
            changed = true;
            distType_ = newDistType;
            setNodeMasters(nProcessorsPerMaster, true);
        }
        else
        {
            Info<< "Did not switch from "
                << distributionTypeNames_[distType_] << " to "
                << distributionTypeNames_[newDistType] << endl;
        }
    }

    // May need to change our approach to bounding boxes with the
    // change in mode.
    if
    (
        changed &&
        (
            distType_ == distributionType::TRANSPARENT
         || distType_ == distributionType::NODEMASTER
        )
    )
    {
        if (decomposeUsingBbs_)
        {
            decomposeUsingBbs_ = false;
            dict_.set("decomposeUsingBbs", false);
        }
        doExclusionBb_ = false;
        dict_.set("exclusionBb", false);
    }

    return changed;
}


void Foam::distributedTriSurfaceMesh::calcVertexNormals
(
    const triSurface& surf,
    List<FixedList<vector, 3>>& vn
) const
{
    // For now don't do features so just use the point normals
    vectorField pointNormals(surf.points().size(), vector::zero);
    for (const auto& f : surf)
    {
        const vector n = f.unitNormal(surf.points());
        for (label pointi : f)
        {
            pointNormals[pointi] += n;
        }
    }
    pointNormals.normalise();

    vn.resize_fill(surf.size(), FixedList<vector, 3>(Foam::zero{}));

    forAll(surf, facei)
    {
        const auto& f = surf[facei];
        forAll(f, fp)
        {
            vn[facei][fp] = pointNormals[f[fp]];
        }
    }
}


bool Foam::distributedTriSurfaceMesh::isInExclusion
(
    const point& start,
    const point& end
) const
{
    const auto myProci = UPstream::myProcNo();

    return
    (
        (myProci < exclusionBb_.size())
     && (exclusionBb_[myProci].contains_all(start, end))
    );
}


Foam::treeBoundBox Foam::distributedTriSurfaceMesh::subtract
(
    const treeBoundBox& a,
    const treeBoundBox& b
)
{
    // Subtract b from a
    // We calculate the biggest single bounding box that is exclusively
    // in a. Done by slicing with all 6 faces of bb and seeing which
    // one returns largest remaining volume.
    //
    // Since the volume is invariant, the largest remaining volume
    // corresponds to the largest remaining fraction in a given direction.

    // New approach with better handling for collapsed boxes (OCT-2026)

    enum states : uint8_t { ignore_dir, keep_lower, keep_upper };

    FixedList<uint8_t, 3> state(states::keep_upper);
    FixedList<scalar, 3> maxFraction(-1);

    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        if (b.min()[dir] >= a.max()[dir] || b.max()[dir] <= a.min()[dir])
        {
            // No overlap
            return a;
        }

        const scalar a_span = a.span(dir);

        if (a_span < VSMALL)
        {
            // Collapsed in this direction...
            state[dir] = states::ignore_dir;
            continue;
        }

        if (b.min()[dir] > a.min()[dir])
        {
            // Plane at bb.min. Keep bit to the left; shift max to left

            scalar fraction = (b.min()[dir] - a.min()[dir]) / a_span;
            if (fraction > 1.0)
            {
                fraction = 1.0;
            }

            maxFraction[dir] = fraction;
            state[dir] = states::keep_lower;
        }

        if (b.max()[dir] > a.min()[dir])
        {
            // Plane at bb.max. Keep bit to the right i.e. modify a.min

            scalar fraction = (a.max()[dir] - b.max()[dir]) / a_span;
            if (fraction > 1.0)
            {
                fraction = 1.0;
            }

            if (maxFraction[dir] < fraction)
            {
                maxFraction[dir] = fraction;
                state[dir] = states::keep_upper;
            }
        }
    }

    // Find direction that keeps max volume (or area for a collapsed box).

    direction maxIndex = vector::nComponents;
    scalar maxValue(-1);

    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        if
        (
            (state[dir] != states::ignore_dir) &&
            (
                (maxIndex == vector::nComponents)
             || (maxValue < maxFraction[dir])
            )
        )
        {
            maxValue = maxFraction[dir];
            maxIndex = dir;
        }
    }

    if (maxValue <= 0 || maxIndex == vector::nComponents)
    {
        return treeBoundBox::null();
    }
    else
    {
        treeBoundBox subBb(a);
        if (state[maxIndex] == states::keep_upper)
        {
            // keep_upper => adjust min
            subBb.min()[maxIndex] = b.max()[maxIndex];
        }
        else
        {
            // keep_lower => adjust max
            subBb.max()[maxIndex] = b.min()[maxIndex];
        }
        return subBb;
    }
}


void Foam::distributedTriSurfaceMesh::calcExclusionBb()
{
    // Update exclusionBb_. Assumes all procs have same procBb_.
    if (!doExclusionBb_)
    {
        exclusionBb_.clear();
        return;
    }
    if (debug & 4)
    {
        InfoInFunction<< endl;
    }

    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);
    const auto numProc = UPstream::nProcs(communicator);

    // Find a local bb that has none of the remote triangles
    // (or rather bb)
    treeBoundBox exclusionBb(procBb_[myProci]);

    forAll(procBb_, proci)
    {
        if (proci != myProci)
        {
            const auto& otherBb = procBb_[proci];
            if (!otherBb.empty() && exclusionBb.overlaps(otherBb))
            {
                exclusionBb = subtract(exclusionBb, otherBb);
            }
        }
    }

    if (!returnReduceOr(exclusionBb.good()))
    {
        // No valid exclusion bb on any processor. Switch off
        WarningInFunction
            << "Did not calculate a valid exclusionBb"
            << " on any processor. Switching exclusion check off"
            << endl;
        exclusionBb_.clear();
        return;
    }

    // Define the exclusion bound boxes
    {
        List<treeBoundBox> bbs(exclusionBb_.release_boxes());
        bbs.resize_nocopy(numProc);

        boundBox::allGather(bbs, exclusionBb, communicator);

        exclusionBb_.setBoxes(std::move(bbs));
        exclusionBb_.rebuild();
    }

    if (debug && UPstream::master(communicator))
    {
        forAll(exclusionBb_, proci)
        {
            Pout<< "proc:" << proci << nl
                << "\tbound-box:" << procBb_[proci] << nl
                << "\texclusion:" << exclusionBb_[proci]
                << endl;
        }

        if (debug & 2)
        {
            forAll(exclusionBb_, proci)
            {
                OBJstream os
                (
                    searchableSurface::time().path()
                   /searchableSurface::name()
                 + "_exclusion_" + Foam::name(proci) + ".obj"
                );
                os.write(exclusionBb_[proci], true);
                Pout<< "** Written bb to " << os.name()
                    << endl;
            }
        }
    }
}


//void Foam::distributedTriSurfaceMesh::splitSegment
//(
//    const label segmenti,
//    const point& start,
//    const point& end,
//    const treeBoundBox& bb,
//
//    DynamicList<segment>& allSegments,
//    DynamicList<label>& allSegmentMap,
//    DynamicList<label> sendMap
//) const
//{
//    // Work points
//    point clipPt0, clipPt1;
//
//    if (bb.contains(start))
//    {
//        // start within, trim end to bb
//        if (bb.intersects(end, start, clipPt0))
//        {
//            // segment from start to clippedStart passes
//            // through proc.
//            sendMap[proci].push_back(allSegments.size());
//            allSegmentMap.push_back(segmenti);
//            allSegments.emplace_back(start, clipPt0);
//        }
//    }
//    else if (bb.contains(end))
//    {
//        // end within, trim start to bb
//        if (bb.intersects(start, end, clipPt0))
//        {
//            sendMap[proci].push_back(allSegments.size());
//            allSegmentMap.push_back(segmenti);
//            allSegments.emplace_back(clipPt0, end);
//        }
//    }
//    else
//    {
//        // trim both
//        if
//        (
//            bb.intersects(start, end, clipPt0)
//         && bb.intersects(end, clipPt0, clipPt1)
//        )
//        {
//            // middle part of segment passes through proc.
//            sendMap[proci].push_back(allSegments.size());
//            allSegmentMap.push_back(segmenti);
//            allSegments.emplace_back(clipPt0, clipPt1);
//        }
//    }
//}


void Foam::distributedTriSurfaceMesh::distributeSegment
(
    const label segmenti,
    const point& start,
    const point& end,

    DynamicList<segment>& allSegments,
    DynamicList<label>& allSegmentMap,
    List<DynamicList<label>>& sendMap
) const
{
    if (debug & 4)
    {
        InfoInFunction<< endl;
    }

    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);
    // const auto numProc = UPstream::nProcs(communicator);

    DynamicList<label> hits;

    // All of the 'contains' targets:
    {
        int targetProci = -1;

        // Is there a region that does not overlap any remote triangles?
        // In which case this is the only processor that needs to test.
        // Note that the calling function might have truncated the end point
        // so now the segment is fully local.
        if (isInExclusion(start, end))
        {
            targetProci = myProci;
        }
        else if (decomposeUsingBbs_ && procBb_.find_contains(start, end, hits))
        {
            // Fully inside local or another processor.
            // In that case, only need to send to that one processor
            // even if it intersects another.
            //
            // Rare occurrence but cheap if it happens

            targetProci = (hits.contains(myProci) ? myProci : hits[0]);
        }
        else
        {
            // 3a. If not contained in single processor send to single exclusion
            //     processor.
            //     Note: if exclusion bb is small compared to overall this
            //     extra test might be more expensive than adding it to the map?

            targetProci = exclusionBb_.first_contains(start, end);
        }

        if (targetProci >= 0)
        {
            sendMap[targetProci].push_back(allSegments.size());
            allSegmentMap.push_back(segmenti);
            allSegments.emplace_back(start, end);
            return;
        }
    }


    // 3. If not contained in single processor send to all intersecting
    // processors.

    // Scheme a: any processor that intersects the segment gets
    // the whole segment.

    for (auto proci : procBb_.all_intersects(start, end, hits))
    {
        sendMap[proci].push_back(allSegments.size());
        allSegmentMap.push_back(segmenti);
        allSegments.emplace_back(start, end);
    }

    // Alternative: any processor only gets clipped bit of
    // segment. This gives small problems with additional
    // truncation errors.

    // forAll(procBb_, proci)
    // {
    //     const treeBoundBox& bb = procBb_[proci];
    //
    //     splitSegment
    //     (
    //         segmenti,
    //         start,
    //         end,
    //         bb,
    //
    //         allSegments,
    //         allSegmentMap,
    //         sendMap[proci]
    //     );
    // }
}


Foam::autoPtr<Foam::mapDistribute>
Foam::distributedTriSurfaceMesh::distributeSegments
(
    const pointField& start,
    const pointField& end,
    const bitSet& doTest,

    List<segment>& allSegments,
    labelList& allSegmentMap
) const
{
    if (debug & 4)
    {
        InfoInFunction<< endl;
    }

    // Determine send map
    // ~~~~~~~~~~~~~~~~~~

    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);
    const auto numProc = UPstream::nProcs(communicator);

    labelListList sendMap(numProc);

    const label nTests = doTest.count();

    {
        // Since intersection test is quite expensive compared to memory
        // allocation we use DynamicList to immediately store the segment
        // in the correct bin.

        // Segments to test
        DynamicList<segment> dynAllSegments(nTests);
        // Original index of segment
        DynamicList<label> dynAllSegmentMap(nTests);
        // Per processor indices into allSegments to send
        List<DynamicList<label>> dynSendMap(numProc);

        // Pre-size
        forAll(dynSendMap, proci)
        {
            dynSendMap[proci].reserve
            (
                (proci == myProci) ? (nTests) : (nTests/numProc)
            );
        }

        for (const label segmenti : doTest)
        {
            distributeSegment
            (
                segmenti,
                start[segmenti],
                end[segmenti],

                dynAllSegments,
                dynAllSegmentMap,
                dynSendMap
            );
        }

        // Convert dynamic to regular list
        forAll(sendMap, proci)
        {
            sendMap[proci].transfer(dynSendMap[proci]);
        }

        allSegments.transfer(dynAllSegments);
        allSegmentMap.transfer(dynAllSegmentMap);
    }

    return autoPtr<mapDistribute>::New(std::move(sendMap));
}


void Foam::distributedTriSurfaceMesh::findLine
(
    const bool nearestIntersection,
    const pointField& start,
    const pointField& initialEnd,
    List<pointIndexHit>& infos
) const
{
    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);
    // const auto numProc = UPstream::nProcs(communicator);

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::findLine :"
            << " intersecting surface " << searchableSurface::name()
            << " with " << start.size() << " rays" << endl;
    }
    addProfiling(findLine, "distributedTriSurfaceMesh::findLine");

    const indexedOctree<treeDataTriSurface>& octree = tree();

    // end-point - updated with local hits
    pointField end(initialEnd);

    // Initialise
    infos.resize_nocopy(start.size());
    for (auto& pi : infos)
    {
        pi.setMiss();
    }

    // Important:force synchronised construction of indexing
    const globalIndex& triIndexer = globalTris();


    // Do any local queries
    // ~~~~~~~~~~~~~~~~~~~~

    // label nLocal = 0;

    // Tests to be done
    bitSet doTest(start.size(), true);

    forAll(start, i)
    {
        const auto& startp = start[i];
        auto& endp = end[i];
        auto& info = infos[i];

        if (isInExclusion(startp, endp))
        {
            // Ray segment is exclusively on myProcNo!

            if (nearestIntersection)
            {
                info = octree.findLine(startp, endp);
            }
            else
            {
                info = octree.findLineAny(startp, endp);
            }

            if (auto idx = info.index(); info.hit())
            {
                info.setIndex(triIndexer.toGlobal(myProci, idx));

                // Update endpoint
                endp = info.point();
            }

            doTest.unset(i);
        }
        else if (procBb_[myProci].contains_all(startp, endp))
        {
            // Ray segment is inside procBb on myProcNo:
            //  decomposeUsingBbs_ : we are guaranteed to have all the triangles
            if (nearestIntersection)
            {
                info = octree.findLine(startp, endp);
            }
            else
            {
                info = octree.findLineAny(startp, endp);
            }

            if (auto idx = info.index(); info.hit())
            {
                info.setIndex(triIndexer.toGlobal(myProci, idx));
                // Update endpoint
                endp = info.point();
            }

            // Done all tests for this ray if:
            // - have hit & any intersection is ok
            // - decomposeUsingBbs : we have all the triangles to make decision
            if (decomposeUsingBbs_ || (!nearestIntersection && info.hit()))
            {
                doTest.unset(i);
            }
        }
    }

    if (returnReduceAnd(doTest.none()))
    {
        // Done all tests (for whatever reason)
        return;
    }


    // Not all can be resolved locally. Build segments and map,
    // send over segments, do intersections, send back and merge.


    // Construct queries (segments)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Segments to test
    List<segment> allSegments;
    // Original index of segment
    labelList allSegmentMap;

    const autoPtr<mapDistribute> mapPtr
    (
        distributeSegments
        (
            start,
            end,
            doTest,
            allSegments,
            allSegmentMap
        )
    );
    const mapDistribute& map = mapPtr();

    const label nOldAllSegments = allSegments.size();


    // Exchange the segments
    // ~~~~~~~~~~~~~~~~~~~~~

    map.distribute(allSegments);


    // Do tests i need to do
    // ~~~~~~~~~~~~~~~~~~~~~

    // Intersections
    List<pointIndexHit> intersections(allSegments.size());

    forAll(allSegments, i)
    {
        auto& allInfo = intersections[i];

        if (nearestIntersection)
        {
            allInfo = octree.findLine
            (
                allSegments[i].first(),
                allSegments[i].second()
            );
        }
        else
        {
            allInfo = octree.findLineAny
            (
                allSegments[i].first(),
                allSegments[i].second()
            );
        }

        // Convert triangle index to global numbering
        if (auto idx = allInfo.index(); allInfo.hit())
        {
            allInfo.setIndex(triIndexer.toGlobal(myProci, idx));
        }
    }


    // Exchange the intersections (opposite to segments)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Note: segments are unique so no potential overlap. Otherwise could
    //       build combineOp to do the min(magSqr) decision below into
    //       the reverse distribute.
    map.reverseDistribute(nOldAllSegments, intersections);


    // Extract the hits
    // ~~~~~~~~~~~~~~~~

    forAll(intersections, i)
    {
        const auto& allInfo = intersections[i];
        const auto segmenti = allSegmentMap[i];

        const auto& startp = start[segmenti];
        auto& hitInfo = infos[segmenti];

        if
        (
            allInfo.hit() &&
            (
                hitInfo.not_hit() ||
                (
                    // Nearest intersection
                    nearestIntersection &&
                    (
                        allInfo.point().distSqr(startp)
                      < hitInfo.point().distSqr(startp)
                    )
                )
            )
        )
        {
            hitInfo = allInfo;
        }
    }
}


void Foam::distributedTriSurfaceMesh::convertTriIndices
(
    UList<pointIndexHit>& infos
) const
{
    // Important:force synchronised construction of indexing
    const globalIndex& triIndexer = globalTris();

    for (auto& pi : infos)
    {
        if (pi.hit())
        {
            pi.setIndex(triIndexer.toGlobal(pi.index()));
        }
        else
        {
            // Probably only without this:
            // pi.setIndex(-1);
        }
    }
}


// Exchanges indices to the processor they come from.
// - calculates exchange map
// - uses map to calculate local triangle index
Foam::autoPtr<Foam::mapDistribute>
Foam::distributedTriSurfaceMesh::calcLocalQueries
(
    const UList<pointIndexHit>& info,
    labelList& triangleIndex
) const
{
    if (debug & 4)
    {
        Pout<< "distributedTriSurfaceMesh::calcLocalQueries :"
               " (processing)" << endl;
    }

    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);
    const auto numProc = UPstream::nProcs(communicator);

    // Note: does not filter duplicate queries so a triangle that gets requested
    //       from more than one processor also get local queried more than
    //       once.

    triangleIndex.resize_fill(info.size(), -1);

    const globalIndex& triIndexer = globalTris();


    // Determine send map
    // ~~~~~~~~~~~~~~~~~~

    // Since determining which processor the query should go to is
    // cheap we do a multi-pass algorithm to save some memory temporarily.

    labelListList sendMap(numProc);

    {
        // 1. Count. Also cache the target proc in triangleIndex to avoid
        // repeating whichProcID() queries.

        labelList sendCount(numProc, Foam::zero{});

        forAll(info, i)
        {
            if (info[i].hit())
            {
                auto proci = triIndexer.whichProcID(myProci, info[i].index());

                ++sendCount[proci];

                // proc index [temporary cache]
                triangleIndex[i] = proci;
            }
        }

        // 2. Size sendMap
        forAll(sendCount, proci)
        {
            sendMap[proci].resize_nocopy(sendCount[proci]);
            sendCount[proci] = 0;
        }

        // 3. Fill sendMap. The triangle proc info already filled in stage 1
        forAll(info, i)
        {
            // recover proc index from [temporary cache]
            auto proci = triangleIndex[i];
            auto globali = info[i].index();

            if (proci >= 0 && info[i].hit())
            {
                sendMap[proci][sendCount[proci]++] = i;

                triangleIndex[i] = triIndexer.toLocal(proci, globali);
            }
        }
    }

    // Pack into distribution map
    auto mapPtr = autoPtr<mapDistribute>::New(std::move(sendMap));

    // Send over queries
    mapPtr().distribute(triangleIndex);

    return mapPtr;
}


// Generate queries for parallel distance calculation
// - calculates exchange map
// - uses map to exchange points and radius
Foam::autoPtr<Foam::mapDistribute>
Foam::distributedTriSurfaceMesh::calcLocalQueries
(
    const bool includeLocalProcessor,
    const pointField& centres,
    const scalarField& radiusSqr,

    pointField& allCentres,
    scalarField& allRadiusSqr,
    labelList& allSegmentMap
) const
{
    if (debug & 4)
    {
        Pout<< "distributedTriSurfaceMesh::calcLocalQueries :"
               " ..." << endl;
    }

    // Determine queries
    // ~~~~~~~~~~~~~~~~~

    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);
    const auto numProc = UPstream::nProcs(communicator);

    labelListList sendMap(numProc);

    {
        // Queries
        DynamicList<point> dynAllCentres(centres.size());
        DynamicList<scalar> dynAllRadiusSqr(centres.size());
        // Original index of segment
        DynamicList<label> dynAllSegmentMap(centres.size());
        // Per processor indices into allSegments to send
        List<DynamicList<label>> dynSendMap(numProc);

        // Pre-size
        forAll(dynSendMap, proci)
        {
            dynSendMap[proci].reserve
            (
                (proci == myProci)
              ? centres.size()
              : centres.size()/numProc
            );
        }

        // Work array
        DynamicList<label> hits;

        forAll(centres, centrei)
        {
            const auto& cc = centres[centrei];
            const scalar radsqr = radiusSqr[centrei];

            if
            (
                label proci = exclusionBb_.first_contains(cc, radsqr);
                proci >= 0
            )
            {
                dynSendMap[proci].push_back(dynAllCentres.size());
                dynAllSegmentMap.push_back(centrei);
                dynAllCentres.push_back(cc);
                dynAllRadiusSqr.push_back(radsqr);
                // Proci contains all the triangles it can potentially
                // interact with
                continue;
            }

            for (auto proci : procBb_.all_overlaps(cc, radsqr, hits))
            {
                if (includeLocalProcessor || (proci != myProci))
                {
                    dynSendMap[proci].push_back(dynAllCentres.size());
                    dynAllSegmentMap.push_back(centrei);
                    dynAllCentres.push_back(cc);
                    dynAllRadiusSqr.push_back(radsqr);
                }
            }
        }

        // Convert dynamic to regular list
        forAll(sendMap, proci)
        {
            sendMap[proci].transfer(dynSendMap[proci]);
        }

        allCentres.transfer(dynAllCentres);
        allRadiusSqr.transfer(dynAllRadiusSqr);
        allSegmentMap.transfer(dynAllSegmentMap);
    }

    return autoPtr<mapDistribute>::New(std::move(sendMap));
}


void Foam::distributedTriSurfaceMesh::surfaceSide
(
    const pointField& samples,
    const UList<pointIndexHit>& nearestInfo,
    List<volumeType>& volType
) const
{
    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::surfaceSide :"
            << " on surface " << searchableSurface::name()
            << " finding surface side given points on surface for "
            << samples.size() << " samples"
            << " cached vertexnormals " << vertexNormals_.valid() << endl;
    }

    // Use global index to send local tri and nearest back to originating
    // processor

    labelList triangleIndex(nearestInfo.size());
    autoPtr<mapDistribute> mapPtr
    (
        calcLocalQueries
        (
            nearestInfo,
            triangleIndex
        )
    );
    const mapDistribute& map = mapPtr();

    // Send over samples
    pointField localSamples(samples);
    map.distribute(localSamples);

    // Send over nearest point since we have it anyway
    List<pointIndexHit> localNearestInfo(nearestInfo);
    map.distribute(localNearestInfo);


    // Do my tests
    // ~~~~~~~~~~~

    volType.resize_fill(triangleIndex.size(), volumeType::UNKNOWN);

    const triSurface& surf = *this;
    const auto& points = surf.points();

    if (vertexNormals_)
    {
        // Use smooth interpolation for the normal

        forAll(triangleIndex, i)
        {
            const label facei = triangleIndex[i];
            const auto& f = surf[facei];
            const auto& vn = vertexNormals_()[facei];

            const point& sample = localSamples[i];

            const point& nearest = localNearestInfo[i].point();

            const vector sampleNearestVec(sample - nearest);

            const barycentric2D w(f.tri(points).pointToBarycentric(nearest));

            const vector n = w[0]*vn[0]+w[1]*vn[1]+w[2]*vn[2];

            const scalar c = (sampleNearestVec & n);

            if (c > 0)
            {
                volType[i] = volumeType::OUTSIDE;
            }
            else
            {
                volType[i] = volumeType::INSIDE;
            }
        }
    }
    else
    {
        //const labelListList& pointFaces = surf.pointFaces();
        // Construct pointFaces. Let's hope surface has compact point
        // numbering ...

        // labelListList tpointFaces;
        // if (!surf.hasPointFaces())
        // {
        //     Foam::invertManyToMany(points.size(), surf, tpointFaces);
        // }
        // const labelListList& pointFaces =
        // (
        //     surf.hasPointFaces() ? surf.pointFaces() : tpointFaces
        // );

        labelListList pointFaces;
        invertManyToMany(points.size(), surf, pointFaces);

        EdgeMap<labelPair> edgeToFaces;

        forAll(triangleIndex, i)
        {
            const label facei = triangleIndex[i];
            const auto& f = surf[facei];
            const point& sample = localSamples[i];

            // Find where point is on face
            label nearType, nearLabel;
            pointHit pHit =
                f.nearestPointClassify(sample, points, nearType, nearLabel);

            const point& nearestPoint = pHit.point();

            if (nearType == triPointRef::NONE)
            {
                const vector sampleNearestVec = (sample - nearestPoint);

                // Nearest to face interior. Use faceNormal to determine side
                //scalar c = sampleNearestVec & surf.faceNormals()[facei];
                const scalar c = (sampleNearestVec & f.areaNormal(points));

                if (c > 0)
                {
                    volType[i] = volumeType::OUTSIDE;
                }
                else
                {
                    volType[i] = volumeType::INSIDE;
                }
            }
            else if (nearType == triPointRef::EDGE)
            {
                // Nearest to edge nearLabel. Note that this can only be a
                // knife-edge
                // situation since otherwise the nearest point could never be
                // the edge.

                // NB: edgeSide also handles an open edge (== UNKNOWN)
                auto side = triSurfaceTools::edgeSide
                (
                    surf,
                    sample,
                    nearestPoint,
                    facei,
                    findOtherFace(*this, pointFaces, facei, nearLabel)
                );

                volType[i] = triSurfaceTools::sideToVolumeType(side);
            }
            else
            {
                // Nearest to point. Could use pointNormal here but is not
                // correct.
                // Instead determine which edge using point is nearest and
                // use test above (nearType == triPointRef::EDGE).

                const label pointi = f[nearLabel];
                const auto& pFaces = pointFaces[pointi];
                const vector sampleNearestVec = (sample - nearestPoint);

                // Loop over all faces. Check if have both edge faces. Do
                // test.
                edgeToFaces.clear();

                scalar maxCosAngle = -GREAT;
                labelPair maxEdgeFaces(-1, -1);

                for (const label facei : pFaces)
                {
                    const auto& tri = surf[facei];
                    const label fp = tri.find(pointi);

                    if (fp < 0)
                    {
                        // Should never happen; just to avoid compiler
                        // warning.
                        FatalErrorInFunction << "Point not found in triangle"
                            << exit(FatalError);
                    }

                    Pair<edge> edges
                    (
                        edge(pointi, tri.fcValue(fp)),  // edge forward
                        edge(pointi, tri.rcValue(fp))   // edge backward

                    );

                    // Check forward edge and backward edge
                    // determine distance/angle to nearPoint
                    for (const edge& e : edges)
                    {
                        if (auto iter = edgeToFaces.find(e); iter.good())
                        {
                            auto& efaces = iter.val();

                            if (efaces.second() == -1)
                            {
                                // Found second face. Now we have edge+faces
                                efaces.second() = facei;

                                vector eVec(e.vec(points));
                                scalar magEVec = mag(eVec);

                                if (magEVec > VSMALL)
                                {
                                    eVec /= magEVec;

                                    // Determine edge most in direction of
                                    // sample
                                    scalar cosAngle(sampleNearestVec&eVec);
                                    if (cosAngle > maxCosAngle)
                                    {
                                        maxCosAngle = cosAngle;
                                        maxEdgeFaces = efaces;
                                    }
                                }
                            }
                            else
                            {
                                FatalErrorInFunction << "Not closed"
                                    << exit(FatalError);
                            }
                        }
                        else
                        {
                            edgeToFaces.insert(e, labelPair(facei, -1));
                        }
                    }
                }


                // Check that surface is closed
                bool closed = true;
                for (auto iter : edgeToFaces)
                {
                    if (iter[0] == -1 || iter[1] == -1)
                    {
                        closed = false;
                        break;
                    }
                }
                if (closed)
                {
                    auto side = triSurfaceTools::edgeSide
                    (
                        surf,
                        sample,
                        nearestPoint,
                        maxEdgeFaces.first(),
                        maxEdgeFaces.second()
                    );

                    volType[i] = triSurfaceTools::sideToVolumeType(side);
                }
            }
        }
    }


    // Send back results
    // ~~~~~~~~~~~~~~~~~

    // Note that we make sure that if multiple processors hold data only
    // the one with inside/outside wins. (though this should already be
    // handled by the fact we have a unique nearest triangle so we only
    // send the volume-query to a single processor)


    //{
    //- Alternative which also sends over source sample so we get better
    //- debug messages.
    //    // Pack sample, nearest and volType so we get nicer error messages
    //    typedef Tuple2<Pair<point>, volumeType> NearType;
    //    List<NearType> nearTypes(volType.size());
    //    forAll(localSamples, i)
    //    {
    //        const point& sample = localSamples[i];
    //        const point& near = localNearestInfo[i].point();
    //        nearTypes[i] = NearType(Pair<point>(sample, near), volType[i]);
    //    }
    //
    //
    //    const NearType zero(Pair<point>(Zero, Zero), volumeType::UNKNOWN);
    //    mapDistributeBase::distribute
    //    (
    //        UPstream::commsTypes::nonBlocking,
    //        UList<labelPair>::null(),
    //        nearestInfo.size(),
    //        map.constructMap(),
    //        map.constructHasFlip(),
    //        map.subMap(),
    //        map.subHasFlip(),
    //        nearTypes,
    //        zero,
    //        NearTypeCombineOp(),
    //        identityOp(),           // No flipping
    //        UPstream::msgType(),
    //        map.comm()
    //    );
    //    volType.setSize(nearTypes.size());
    //    forAll(nearTypes, i)
    //    {
    //        volType[i] = nearTypes[i].second();
    //    }
    //}


    const volumeType zero(volumeType::UNKNOWN);
    mapDistributeBase::distribute
    (
        UPstream::commsTypes::nonBlocking,
        UList<labelPair>::null(),
        nearestInfo.size(),
        map.constructMap(),
        map.constructHasFlip(),
        map.subMap(),
        map.subHasFlip(),
        volType,
        zero,
        volumeCombineOp(),
        identityOp(),           // No flipping
        UPstream::msgType(),
        map.comm()
    );

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::surfaceSide :"
            << " finished finding surface " << searchableSurface::name()
            << " given points on surface "
            << searchableSurface::name() << " for "
            << samples.size() << " samples" << endl;
    }
}


Foam::volumeType Foam::distributedTriSurfaceMesh::markMixed
(
    const indexedOctree<treeDataTriSurface>& tree,
    const label nodei,
    const triPointRef& tri,
    PackedList<2>& nodeTypes
) const
{
    using tree_base = indexedOctreeBase;
    const auto& nod = tree.nodes()[nodei];

    //Pout<< indent
    //    << "Entering node:" << nodei
    //    << " with bb:" << nod.bb_ << endl;
    //Pout<< incrIndent;

    volumeType myType = volumeType::UNKNOWN;

    for (direction octant = 0; octant < nod.subNodes_.size(); octant++)
    {
        const labelBits index = nod.subNodes_[octant];
        const treeBoundBox subBb(nod.bb_.subBbox(octant));

        if (subBb.intersects(tri))
        {
            volumeType subType;

            if (tree_base::isNode(index))
            {
                //Pout<< indent
                //    << "Testing node octant:" << octant << " with node:"
                //    << indexedOctree<treeDataTriSurface>::getNode(index)
                //    << " subBb:" << subBb << endl;

                // tree node. Recurse.
                subType = markMixed
                (
                    tree,
                    tree_base::getNode(index),
                    tri,
                    nodeTypes
                );
            }
            else if (tree_base::isContent(index))
            {
                // Content inside so we cannot pre-determine whether
                // whole node box is inside or outside

                //const auto indices = tree.getIndices(index);
                //Pout<< indent
                //    << "Testing content octant:" << octant
                //    << " subBb:" << subBb
                //    << " holding shapes:" << flatOutput(indices)
                //    << endl;

                // Store octant type
                subType = volumeType::MIXED;
            }
            else
            {
                // Nothing in it but intersects triangle. Can't cache

                //Pout<< indent
                //    << "Testing empty octant:" << octant
                //    << " subBb:" << subBb << endl;

                // Store octant type
                subType = volumeType::MIXED;
            }

            // Store octant type
            nodeTypes.set(labelBits::pack(nodei, octant), subType);

            // Combine sub node types into type for treeNode. Result is
            // 'mixed' if types differ among subnodes.
            if (myType == volumeType::UNKNOWN)
            {
                myType = subType;
            }
            else if (subType != myType)
            {
                //Pout<< indent
                //    << "for octant:" << octant
                //    << " siwtched from myType:" << myType
                //    << " and subType:" << subType
                //    << " switched to MIXED" << endl;
                myType = volumeType::MIXED;
            }
        }
    }

    //Pout<< indent
    //    << "Result for node:" << nodei << " myType:" << myType << endl;
    //Pout<< decrIndent;

    return myType;
}


void Foam::distributedTriSurfaceMesh::markMixed
(
    const indexedOctree<treeDataTriSurface>& tree,
    const label nodei,
    const triSurface& surf,
    PackedList<2>& nodeTypes
) const
{
    for (const auto& f : surf)
    {
        markMixed(tree, nodei, f.tri(surf.points()), nodeTypes);
    }
}


void Foam::distributedTriSurfaceMesh::collectLeafMids
(
    const label nodeI,
    DynamicField<point>& midPoints
) const
{
    const auto& nod = tree().nodes()[nodeI];

    for (direction octant = 0; octant < nod.subNodes_.size(); octant++)
    {
        const labelBits index = nod.subNodes_[octant];

        if (indexedOctree<treeDataTriSurface>::isNode(index))
        {
            // tree node. Recurse.
            collectLeafMids
            (
                indexedOctree<treeDataTriSurface>::getNode(index),
                midPoints
            );
        }
        else if (indexedOctree<treeDataTriSurface>::isContent(index))
        {}
        else
        {
            // No data in this octant. Set type for octant acc. to the mid
            // of its bounding box.
            const treeBoundBox subBb = nod.bb_.subBbox(octant);
            midPoints.push_back(subBb.centre());
        }
    }
}


Foam::volumeType Foam::distributedTriSurfaceMesh::calcVolumeType
(
    const UList<volumeType>& midPointTypes,
    label& midPointi,
    PackedList<2>& nodeTypes,
    const label nodeI
) const
{
    // Pre-calculates wherever possible the volume status per node/subnode.
    // Recurses to determine status of lowest level boxes. Level above is
    // combination of octants below.

    const auto& nod = tree().nodes()[nodeI];

    volumeType myType = volumeType::UNKNOWN;

    for (direction octant = 0; octant < nod.subNodes_.size(); octant++)
    {
        volumeType subType;

        const labelBits index = nod.subNodes_[octant];

        if (indexedOctree<treeDataTriSurface>::isNode(index))
        {
            // tree node. Recurse.
            subType = calcVolumeType
            (
                midPointTypes,
                midPointi,
                nodeTypes,
                indexedOctree<treeDataTriSurface>::getNode(index)
            );
        }
        else if (indexedOctree<treeDataTriSurface>::isContent(index))
        {
            // Contents. Depending on position in box might be on either
            // side.
            subType = volumeType::MIXED;
        }
        else
        {
            // No data in this octant. Set type for octant acc. to the mid
            // of its bounding box.
            //Pout<< "    for leaf at bb:" << nod.bb_.subBbox(octant)
            //    << " node:" << nodeI
            //    << " octant:" << octant
            //    << " caching volType to:" << midPointTypes[midPointi] << endl;
            subType = midPointTypes[midPointi++];
        }

        // Store octant type
        nodeTypes.set(labelBits::pack(nodeI, octant), subType);

        // Combine sub node types into type for treeNode. Result is 'mixed' if
        // types differ among subnodes.
        if (myType == volumeType::UNKNOWN)
        {
            myType = subType;
        }
        else if (subType != myType)
        {
            myType = volumeType::MIXED;
        }
    }
    return myType;
}


void Foam::distributedTriSurfaceMesh::cacheVolumeType(PackedList<2>& nt) const
{
    // Get local tree
    const indexedOctree<treeDataTriSurface>& t = tree();
    const auto& nodes = t.nodes();

    // Collect midpoints
    DynamicField<point> midPoints(label(0.5*nodes.size()));
    if (nodes.size())
    {
        collectLeafMids(0, midPoints);
    }

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::cacheVolumeType :"
            << " surface " << searchableSurface::name()
            << " triggering orientation caching for "
            << midPoints.size() << " leaf mids" << endl;
    }

    // Get volume type of mid points
    List<volumeType> midVolTypes;

    // Note: could recurse here (since now outsideVolType_ is set)
    // but this would use the cached mid point data which we're trying
    // to calculate. Instead bypass caching and do a full search
    {
        List<pointIndexHit> nearestInfo;
        findNearest
        (
            midPoints,
            scalarField(midPoints.size(), Foam::sqr(GREAT)),
            nearestInfo
        );
        surfaceSide(midPoints, nearestInfo, midVolTypes);

        if (debug & 2)
        {
            OBJstream insideStr
            (
                searchableSurface::time().path()
               /searchableSurface::name() + "_cacheVolumeType_inside.obj"
            );
            OBJstream outsideStr
            (
                searchableSurface::time().path()
               /searchableSurface::name() + "_cacheVolumeType_outside.obj"
            );
            forAll(midVolTypes, i)
            {
                const linePointRef ln(midPoints[i], nearestInfo[i].point());
                if (midVolTypes[i] == volumeType::INSIDE)
                {
                    insideStr.write(ln);
                }
                else if (midVolTypes[i] == volumeType::OUTSIDE)
                {
                    outsideStr.write(ln);
                }
            }
            Pout<< "Whilst caching " << searchableSurface::name()
                << " have inside:" << insideStr.nVertices()
                << " have outside:" << outsideStr.nVertices()
                << endl;
        }
    }

    // Cache on local tree
    if (nodes.size())
    {
        label index = 0;
        calcVolumeType
        (
            midVolTypes,
            index,
            nt,
            0               // nodeI
        );
    }
    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::cacheVolumeType :"
            << " surface " << searchableSurface::name()
            << " done orientation caching for "
            << midPoints.size() << " leaf mids" << endl;
    }
}


Foam::volumeType Foam::distributedTriSurfaceMesh::cachedVolumeType
(
    const label nodeI,
    const point& sample
) const
{
    const auto& nod = tree().nodes()[nodeI];

    direction octant = nod.bb_.subOctant(sample);

    volumeType octantType = volumeType::type
    (
        tree().nodeTypes().get(labelBits::pack(nodeI, octant))
    );

    if (octantType == volumeType::INSIDE)
    {
        return octantType;
    }
    else if (octantType == volumeType::OUTSIDE)
    {
        return octantType;
    }
    else if (octantType == volumeType::UNKNOWN)
    {
        // Can happen for e.g. non-manifold surfaces.
        return octantType;
    }
    else if (octantType == volumeType::MIXED)
    {
        labelBits index = nod.subNodes_[octant];

        if (indexedOctree<treeDataTriSurface>::isNode(index))
        {
            // Recurse
            volumeType subType = cachedVolumeType
            (
                indexedOctree<treeDataTriSurface>::getNode(index),
                sample
            );

            return subType;
        }
        else if (indexedOctree<treeDataTriSurface>::isContent(index))
        {
            // Content.
            return volumeType::UNKNOWN;
        }
        else
        {
            // Empty node. Cannot have 'mixed' as its type since not divided
            // up and has no items inside it.
            FatalErrorInFunction
                << "Sample:" << sample << " node:" << nodeI
                << " with bb:" << nod.bb_ << nl
                << "Empty subnode has invalid volume type MIXED."
                << abort(FatalError);

            return volumeType::UNKNOWN;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Sample:" << sample << " at node:" << nodeI
            << " octant:" << octant
            << " with bb:" << nod.bb_.subBbox(octant) << nl
            << "Node has invalid volume type " << octantType
            << abort(FatalError);

        return volumeType::UNKNOWN;
    }
}


Foam::label Foam::distributedTriSurfaceMesh::findTriangle
(
    const UList<labelledTri>& allFaces,  // haystack
    const UList<labelList>& allPointFaces,
    const labelledTri& f                 // needle
)
{
    // Search allFaces for a match by checking the point-faces
    // attached to the min vertex

    const label minVert = f.min();

    for (label facei : allPointFaces[minVert])
    {
        const auto& a = allFaces[facei];

        // same region, same vertices and orientation
        if (a.region() == f.region() && triFace::compare(a, f) > 0)
        {
            return facei;
        }
    }

    return -1;
}


Foam::label Foam::distributedTriSurfaceMesh::findTriangle
(
    const UList<labelledTri>& allFaces,  // haystack
    const Map<DynamicList<label>>& allPointFaces,
    const labelledTri& f                 // needle
)
{
    // Search allFaces for a match by checking the point-faces
    // attached to the min vertex

    const label minVert = f.min();

    if (const auto iter = allPointFaces.find(minVert); iter.good())
    {
        for (label facei : iter.val())
        {
            const auto& a = allFaces[facei];

            // same region, same vertices and orientation
            if (a.region() == f.region() && triFace::compare(a, f) > 0)
            {
                return facei;
            }
        }
    }

    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributedTriSurfaceMesh::distributedTriSurfaceMesh
(
    const IOobject& io,
    const triSurface& s,
    const dictionary& dict
)
:
    triSurfaceMesh(io, s),
    dict_
    (
        IOobject
        (
            searchableSurface::name() + "Dict",
            searchableSurface::instance(),
            searchableSurface::local(),
            searchableSurface::db(),
            searchableSurface::NO_READ,
            searchableSurface::writeOpt(),
            searchableSurface::registerObject()
        ),
        dict
    ),

    decomposeUsingBbs_(true),
    doExclusionBb_(true),
    distType_(distributionType::DISTRIBUTED),
    currentDistType_(distributionType::FROZEN),  // to trigger re-distribution
    mergeDist_(SMALL),
    control_(controlType::ctrlDefault),
    comm_(-1)
{
    check_debug(dict);

    // Construct from components. If empty dictionary get all the info from the
    // supplied (undecomposed) triSurface.

    const bool undecomposed = dict_.empty();

    upgradeDictionary(dict_);

    readSettings(undecomposed);

    bounds().reduce(UPstream::worldComm);

    if (debug)
    {
        InfoInFunction
            << "Constructed from triSurface "
            << searchableSurface::name() << nl;

        writeStats(Info);
        printBounds();
    }
}


Foam::distributedTriSurfaceMesh::distributedTriSurfaceMesh(const IOobject& io)
:
    triSurfaceMesh
    (
        IOobject
        (
            io.name(),
            findLocalInstance(io),  // findInstance with parent searching
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            io.registerObject()
        ),
        triSurfaceMesh::masterOnly    // allow parent searching
    ),
    dict_
    (
        IOobject
        (
            searchableSurface::name() + "Dict",
            searchableSurface::instance(),
            searchableSurface::local(),
            searchableSurface::db(),
            IOobjectOption::lazierRead(searchableSurface::readOpt()),
            searchableSurface::writeOpt(),
            searchableSurface::registerObject()
        ),
        dictionary()
    ),

    decomposeUsingBbs_(true),
    doExclusionBb_(true),
    distType_(distributionType::DISTRIBUTED),
    currentDistType_(distributionType::FROZEN),  // to trigger re-distribution
    mergeDist_(SMALL),
    control_(controlType::ctrlDefault),
    comm_(-1)
{
    // Try and find out where the triSurface was loaded from. On slave this
    // might give empty filename.
    const fileName actualFile(triSurfaceMesh::findFile(io, true));

    // Was it read undecomposed?
    const bool readFromMaster =
    (
        actualFile.empty()
     || actualFile != io.localFilePath(triSurfaceMesh::typeName)
    );

    upgradeDictionary(dict_);

    readSettings(readFromMaster);

    bounds().reduce(UPstream::worldComm);

    if (readFromMaster)
    {
        // Make sure orientation is consistent in case we
        // want to use for inside/outside testing
        DebugInFunction
            << "Forcing consistent normals since undecomposed" << endl;

        auto& surf = static_cast<triSurface&>(*this);
        orientedSurface::orientConsistent(surf);
    }

    if (readFromMaster && UPstream::parRun())
    {
        DebugInFunction
            << "Checking for adjustment to the decomposition strategy"
            << endl;

        optimizeStrategy();
    }
    else if (!readFromMaster && UPstream::parRun())
    {
        DebugInFunction
            << "No adjustment to the decomposition strategy"
            << " (not read from master)" << endl;
    }

    if
    (
        readFromMaster && !decomposeUsingBbs_
     && (distType_ != distributionType::NODEMASTER)
    )
    {
        // No fill-in so store normals
        DebugInFunction
            << "Determining vertex based normals since undecomposed" << endl;

        const triSurface& surf = *this;
        vertexNormals_ = autoPtr<List<FixedList<vector, 3>>>::New();
        calcVertexNormals(surf, vertexNormals_());
    }


    const bool canDistribute =
    (
        distType_ != distributionType::FROZEN
     && distType_ != distributionType::FOLLOW
    );

    if (readFromMaster && canDistribute)
    {
        DebugInFunction
            << "Read distributedTriSurface " << io.name()
            << " from parent path " << actualFile << endl;

        if (UPstream::parRun())
        {
            // Distribute (checks that distType != currentDistType_ so should
            // always trigger re-distribution)

            autoPtr<mapDistribute> faceMap;
            autoPtr<mapDistribute> pointMap;
            distribute
            (
                List<treeBoundBox>::null(),  // followBb
                true,       // keep unmapped triangles
                faceMap,
                pointMap
            );
        }
    }
    else
    {
        if (debug)
        {
            InfoInFunction
                << "Read distributedTriSurface " << io.name()
                << " from actual path " << actualFile << ':' << nl;

            printBounds();
        }
    }

    if (debug)
    {
        InfoInFunction
            << "Read distributedTriSurface " << io.name() << ':' << endl;
        writeStats(Info);
    }
}


Foam::distributedTriSurfaceMesh::distributedTriSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict
)
:
    triSurfaceMesh
    (
        IOobject
        (
            io.name(),
            findLocalInstance(io),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            io.registerObject()
        ),
        dict,
        triSurfaceMesh::masterOnly
    ),
    dict_
    (
        IOobject
        (
            searchableSurface::name() + "Dict",
            searchableSurface::instance(),
            searchableSurface::local(),
            searchableSurface::db(),
            IOobjectOption::lazierRead(searchableSurface::readOpt()),
            searchableSurface::writeOpt(),
            searchableSurface::registerObject()
        ),
        dictionary()
    ),

    decomposeUsingBbs_(true),
    doExclusionBb_(true),
    distType_(distributionType::DISTRIBUTED),
    currentDistType_(distributionType::FROZEN),  // to trigger re-distribution
    mergeDist_(SMALL),
    control_(controlType::ctrlDefault),
    comm_(-1)
{
    check_debug(dict);
    if (debug)
    {
        Info<< "distributedTriSurface: construct from dictionary"
            << dict << endl;
    }


    // Try and find out where the triSurface was loaded from. On slave this
    // might give empty filename.
    const fileName actualFile(triSurfaceMesh::findFile(io, dict, true));

    // Was it read undecomposed?
    const bool readFromMaster =
    (
        actualFile.empty()
     || actualFile != io.localFilePath(triSurfaceMesh::typeName)
    );

    upgradeDictionary(dict_);

    // Add supplied dictionary entries (override existing)
    dict_ <<= dict;

    readSettings(readFromMaster);

    if (readFromMaster && UPstream::parRun())
    {
        DebugInFunction
            << "Checking for adjustment to the decomposition strategy"
            << endl;

        optimizeStrategy();
    }
    else if (!readFromMaster && UPstream::parRun())
    {
        DebugInFunction
            << "No adjustment to the decomposition strategy"
            << " (not read from master)" << endl;
    }

    bounds().reduce(UPstream::worldComm);

    if
    (
        readFromMaster && !decomposeUsingBbs_
     && (distType_ != distributionType::NODEMASTER)
    )
    {
        // No fill-in so store enough normals to calculate
        DebugInFunction
            << "Determining vertex based normals since undecomposed" << endl;

        const triSurface& surf = *this;
        vertexNormals_ = autoPtr<List<FixedList<vector, 3>>>::New();
        calcVertexNormals(surf, vertexNormals_());
    }


    const bool canDistribute =
    (
        distType_ != distributionType::FROZEN
     && distType_ != distributionType::FOLLOW
    );

    if (readFromMaster && canDistribute)
    {
        DebugInFunction
            << "Read distributedTriSurface " << io.name()
            << " from parent path " << actualFile
            << " and dictionary" << endl;

        if (UPstream::parRun())
        {
            // Distribute (checks that distType != currentDistType_ so should
            // always trigger re-distribution)
            autoPtr<mapDistribute> faceMap;
            autoPtr<mapDistribute> pointMap;
            distribute
            (
                List<treeBoundBox>::null(),  // followBb
                true,       // keep unmapped triangles
                faceMap,
                pointMap
            );
        }
    }
    else
    {
        if (debug)
        {
            InfoInFunction
                << "Read distributedTriSurface " << io.name()
                << " from actual path " << actualFile
                << " and dictionary:" << nl;

            printBounds();
        }
    }

    if (debug)
    {
        InfoInFunction
            << "Read distributedTriSurface " << io.name() << ':' << endl;
        writeStats(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributedTriSurfaceMesh::~distributedTriSurfaceMesh()
{
    clearOut();
    freeCommunicator();
}


void Foam::distributedTriSurfaceMesh::clearOut()
{
    globalTris_.clear();
    triSurfaceMesh::clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::globalIndex& Foam::distributedTriSurfaceMesh::globalTris() const
{
    if (!globalTris_)
    {
        label localSize = triSurface::size();

        if
        (
            distType_ == distributionType::TRANSPARENT
         || distType_ == distributionType::NODEMASTER
        )
        {
            // NodeMaster distribution.
            // Behave as if all triangles on proc0 only.

            if (!UPstream::master())
            {
                localSize = 0;
            }
        }

        globalTris_.emplace(localSize);
    }

    return *globalTris_;
}


void Foam::distributedTriSurfaceMesh::flip()
{
    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::flip :"
            << " surface " << searchableSurface::name()
            << " with current orientation "
            << volumeType::names[outsideVolType_]
            << " vertexNormals_ " << vertexNormals_.good()
            << endl;
    }

    triSurfaceMesh::flip();

    if (vertexNormals_)
    {
        for (auto& normals : *vertexNormals_)
        {
            for (auto& n : normals)
            {
                n = -n;
            }
        }
    }
}


void Foam::distributedTriSurfaceMesh::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::boundingSpheres :" << endl;
    }

    if
    (
        (distType_ == distributionType::NODEMASTER)
     && (!runningTransparent() && !UPstream::master(comm_))
    )
    {
        // No surface on sub-ranks
        centres.clear();
        radiusSqr.clear();
        return;
    }

    triSurfaceMesh::boundingSpheres(centres, radiusSqr);
}


void Foam::distributedTriSurfaceMesh::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);
    const auto numProc = UPstream::nProcs(communicator);

    if (runningTransparent())
    {
        triSurfaceMesh::findNearest(samples, nearestDistSqr, info);
        return;
    }

    addProfiling
    (
        findNearest,
        "distributedTriSurfaceMesh::findNearest"
    );

    if (debug)
    {
        scalar maxDistSqr = 0;
        if (!samples.empty())
        {
            if (nearestDistSqr.empty())
            {
                maxDistSqr = searchableSurface::defaultMaxDistSqr();
            }
            else
            {
                maxDistSqr = Foam::max(nearestDistSqr);
            }
        }

        Pout<< "distributedTriSurfaceMesh::findNearest :"
            << " surface " << searchableSurface::name()
            << " trying to find nearest for " << samples.size()
            << " samples with max sphere " << std::sqrt(maxDistSqr) << endl;
    }

    // -----------
    // Node master
    // -----------

    if
    (
        (distType_ == distributionType::NODEMASTER)
     && (control_ & controlType::ctrlOverlappingNodeMaster)
    )
    {
        // Service all queries via the node-local master
        nodeMasterQuery
        (
            comm_,
            UPstream::msgType(),
            [this]
            (
                const pointField& s,
                const scalarField& d,
                List<pointIndexHit>& out
            )
            {
                triSurfaceMesh::findNearest(s, d, out);
            },
            info,
            samples, nearestDistSqr
        );
        return;
    }

    if
    (
        (distType_ == distributionType::NODEMASTER)
     && !(control_ & controlType::ctrlOverlappingNodeMaster)
    )
    {
        // Gather all inputs, handle on node-local master, scatter result
        constexpr auto commType = UPstream::commsTypes::nonBlocking;
        const int tag = UPstream::msgType();

        const globalIndex gi(globalIndex::gatherOnly{}, samples.size(), comm_);

        const pointField allSamples
        (
            gi.gather(samples, tag, commType, comm_)
        );
        const scalarField allDistSqr
        (
            gi.gather(nearestDistSqr, tag, commType, comm_)
        );

        List<pointIndexHit> allInfo;
        if (UPstream::master(comm_))
        {
            triSurfaceMesh::findNearest(allSamples, allDistSqr, allInfo);
        }

        info.resize_nocopy(samples.size());
        gi.scatter(allInfo, info, tag, commType, comm_);
        return;
    }


    // ----------------
    // Regular handling
    // ----------------

    const globalIndex& triIndexer = globalTris();

    // Finding nearest:
    // 1. sphere (centre+radius) completely outside all bbs. Easy.
    // 2. sphere completely inside a single exclusive bb. Only test on this bb.
    // 3. centre outside all bbs, but sphere potentially inside. Send to
    //    all interacting bbs. Do two pass.
    // 4. sphere inside potentially multiple bbs. Do two pass.

    // Two-pass searching:
    // 1. send the sample to the processor whose bb contains it. This is
    //    most likely also the one that holds the nearest triangle. (In case
    //    there is no containing processor send to nearest processors. Note
    //    that this might cause a lot of traffic if this is likely)
    //    Send the resulting nearest point back.
    // 2. with the find from 1 look at which other processors might have a
    //    better triangle. Since hopefully step 1) will have produced a tight
    //    bounding box this should limit the amount of points to be retested


    // 1. Test samples (centre only) on processor(s) that contains them
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    scalarField distSqr(samples.size());
    {
        scalar maxDistSqr = searchableSurface::defaultMaxDistSqr();

        forAll(samples, i)
        {
            if (i < nearestDistSqr.size()) { maxDistSqr = nearestDistSqr[i]; }

            distSqr[i] = maxDistSqr;
        }
    }

    autoPtr<mapDistribute> map1Ptr;


    // With exclusionBb we can do a single pass since we know we
    // have all the triangles local to make the correct decision.
    bitSet doneSample(samples.size());

    label nExclusive = 0;   // Num samples exclusively on one proc
    label nInside = 0;      // Num samples inside on one or more procs
    label nMixed = 0;       // Num samples overlapping procs (rare?)
    label nOutside = 0;     // Total fail?

    {
        List<DynamicList<label>> dynSendMap(numProc);
        // Pre-size. Assume samples are uniformly distributed
        forAll(dynSendMap, proci)
        {
            dynSendMap[proci].reserve(samples.size()/numProc);
        }

        DynamicList<label> hits;

        forAll(samples, samplei)
        {
            const auto& centre = samples[samplei];
            const auto& radiusSqr = distSqr[samplei];

            if
            (
                label proci = exclusionBb_.first_contains(centre, radiusSqr);
                proci >= 0
            )
            {
                ++nExclusive;

                dynSendMap[proci].push_back(samplei);
                doneSample.set(samplei);

                if (debug & 2)
                {
                    Pout<< " sample[" << samplei << "] exclude="
                        << proci << endl;
                }
                continue;
            }

            Tuple2<scalar, volumeType> found =
                procBb_.findNearest(centre, radiusSqr, hits);

            scalar tightest = found.first();

            switch (found.second())
            {
                case volumeType::INSIDE :
                {
                    // processor fully contains sphere
                    ++nInside;

                    // Send to all processors containing centre.
                    // This is usually (hopefully) only one.
                    // Can already adjust distSqr for max distance to
                    // search.

                    distSqr[samplei] = tightest;

                    for (auto proci : hits)
                    {
                        dynSendMap[proci].push_back(samplei);
                    }

                    if (debug & 2)
                    {
                        Pout<< " sample[" << samplei << "] contains: ";
                        hits.writeList(Pout) << endl;
                    }
                    break;
                }

                case volumeType::MIXED :
                {
                    // Sample is outside all bb.
                    ++nMixed;

                    distSqr[samplei] = tightest;

                    // Send to single, 'best', processor
                    dynSendMap[hits[0]].push_back(samplei);

                    // Alternative: send to all overlapping processors
                    //
                    // for (auto proci : hits)
                    // {
                    //     dynSendMap[proci].push_back(samplei);
                    // }

                    if (debug & 2)
                    {
                        Pout<< " sample[" << samplei << "] best=" << hits[0]
                            << " overlaps: ";
                        hits.writeList(Pout) << endl;
                    }
                    break;
                }

                default :  // OUTSIDE | UNKNOWN
                {
                    ++nOutside;
                    break;
                }
            }
        }

        // Convert dynamic to regular list
        labelListList sendMap(numProc);
        forAll(sendMap, proci)
        {
            sendMap[proci].transfer(dynSendMap[proci]);
        }
        map1Ptr.emplace(std::move(sendMap));
    }
    const mapDistribute& map1 = *map1Ptr;


    if (debug)
    {
        Pout<< searchableSurface::name() << " Pass1:"
            << " of " << samples.size() << " samples sending to" << endl;
        label nSend = 0;
        forAll(map1.subMap(), proci)
        {
            Pout<< "    " << proci << "\t" << map1.subMap()[proci].size()
                << endl;
            nSend += map1.subMap()[proci].size();
        }
        Pout<< "    sum\t" << nSend << nl
            << "    exclusive:" << nExclusive
            << " inside:" << nInside
            << " overlaps:" << nMixed
            << " outside:" << nOutside << endl;
    }


    List<nearestAndDist> nearestInfo;
    {
        // Get the points I need to test and test locally
        pointField localPoints(samples);
        map1.distribute(localPoints);

        scalarField localDistSqr(distSqr);
        map1.distribute(localDistSqr);

        List<pointIndexHit> localInfo;
        triSurfaceMesh::findNearest(localPoints, localDistSqr, localInfo);
        convertTriIndices(localInfo);

        // Pack into structure for combining information from multiple
        // processors
        nearestInfo.resize_fill
        (
            localInfo.size(),
            nearestAndDist(pointIndexHit(), Foam::sqr(GREAT))
        );

        label nHit = 0;
        label nIgnoredHit = 0;

        const treeBoundBox& myProcBb = procBb_[myProci];

        forAll(nearestInfo, i)
        {
            const pointIndexHit& info = localInfo[i];

            if (const auto& p = info.point(); info.hit())
            {
                nHit++;

                if
                (
                    decomposeUsingBbs_
                 && surfaceClosed_
                 && !myProcBb.contains(p)
                )
                {
                    // Nearest point is not on local processor so the
                    // the triangle is only there because some other bit
                    // of it is on it. Assume there is another processor that
                    // holds the full surrounding of the triangle so we can
                    // ignore this particular nearest.
                    nIgnoredHit++;
                }
                else
                {
                    // 'Normal' situation : we've found the hopefully only
                    // occurence of nearest point

                    nearestInfo[i].first()  = info;
                    nearestInfo[i].second() = p.distSqr(localPoints[i]);
                }
            }
        }

        if (debug)
        {
            scalar maxRadius2 = 0;
            if (localDistSqr.size())
            {
                maxRadius2 = std::sqrt(Foam::max(localDistSqr));
            }

            Pout<< "distributedTriSurfaceMesh::findNearest :"
                << " surface " << searchableSurface::name()
                << " searched locally for " << localPoints.size()
                << " samples with max sphere " << maxRadius2
                << " found hits:" << nHit
                << " of which outside local bb:" << nIgnoredHit
                << endl;
        }
    }

    // Send back to originating processor. Choose best if sent to multiple
    // processors. Note that afterwards all unused entries have the unique
    // value nearestZero (distance < 0). This is used later on to see if
    // the sample was sent to any processor.
    mapDistributeBase::distribute
    (
        UPstream::commsTypes::nonBlocking,
        UList<labelPair>::null(),
        samples.size(),
        map1.constructMap(),
        map1.constructHasFlip(),
        map1.subMap(),
        map1.subHasFlip(),
        nearestInfo,
        nearestZero,
        nearestEqOp(),
        identityOp(),           // No flipping
        UPstream::msgType(),
        map1.comm()
    );


    // 2. Test samples on other processor(s) that overlap
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Now we have (in nearestInfo) for every input sample the current best
    // hit (on the processor that originates the sample). See if we can
    // improve it by sending the queries to any other processors

    autoPtr<mapDistribute> map2Ptr;
    {
        DynamicList<label> hits;
        List<DynamicList<label>> dynSendMap(numProc);

        forAll(nearestInfo, samplei)
        {
            if (doneSample.test(samplei))
            {
                // Sample fully inside exclusionBb so above findNearest
                // would already have found the best possible and is already
                // in nearestInfo.
                continue;
            }

            const point& sample = samples[samplei];
            const auto& near = nearestInfo[samplei];
            const auto& info = near.first();

            label hitProci = -1;
            scalar radsqr = distSqr[samplei];

            if (info.hit())
            {
                hitProci = triIndexer.whichProcID(info.index());
                radsqr = near.second();
            }

            // Find the processors this sample+radius overlaps.

            for (auto proci : procBb_.all_overlaps(sample, radsqr, hits))
            {
                // Check this sample wasn't already handled above. This
                // could be improved since the sample might have been
                // searched on multiple processors. We now only exclude the
                // processor where the point was inside.
                if (proci != hitProci)
                {
                    dynSendMap[proci].push_back(samplei);
                }
            }
        }

        // Convert dynamic to regular list
        labelListList sendMap(numProc);
        forAll(sendMap, proci)
        {
            sendMap[proci].transfer(dynSendMap[proci]);
        }
        map2Ptr.emplace(std::move(sendMap));
    }

    const mapDistribute& map2 = map2Ptr();

    if (debug)
    {
        Pout << searchableSurface::name() << " Pass2:"
            << " of " << samples.size() << " samples sending to" << endl;
        label nSend = 0;
        forAll(map2.subMap(), proci)
        {
            Pout<< "    " << proci << "\t" << map2.subMap()[proci].size()
                << endl;
            nSend += map2.subMap()[proci].size();
        }
        Pout<< "    sum\t" << nSend << endl;
    }

    // Send samples and current best distance
    pointField localSamples(samples);
    map2.distribute(localSamples);

    scalarField localDistSqr(distSqr);
    forAll(nearestInfo, samplei)
    {
        const auto& ni = nearestInfo[samplei];
        if (ni.first().hit())
        {
            localDistSqr[samplei] = ni.second();
        }
    }
    map2.distribute(localDistSqr);

    // Do local test
    List<pointIndexHit> localInfo;
    triSurfaceMesh::findNearest(localSamples, localDistSqr, localInfo);
    convertTriIndices(localInfo);

    // Pack and send back
    List<nearestAndDist> localBest(localSamples.size());
    label nHit = 0;
    label nIgnoredHit = 0;

    const treeBoundBox& myProcBb = procBb_[myProci];

    forAll(localInfo, i)
    {
        const pointIndexHit& info = localInfo[i];

        if (const auto& p = info.point(); info.hit())
        {
            nHit++;
            if
            (
                decomposeUsingBbs_
            &&  surfaceClosed_
            && !myProcBb.contains(info.point())
            )
            {
                // See above
                nIgnoredHit++;
            }
            else
            {
                localBest[i] = { info, p.distSqr(localSamples[i]) };
            }
        }
    }

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::findNearest :"
            << " surface " << searchableSurface::name()
            << " searched locally for " << localSamples.size()
            << " samples with max sphere "
            << (localDistSqr.size() ? Foam::sqrt(max(localDistSqr)) : Zero)
            << " found hits:" << nHit
            << " of which outside local bb:" << nIgnoredHit
            << endl;
    }

    mapDistributeBase::distribute
    (
        UPstream::commsTypes::nonBlocking,
        UList<labelPair>::null(),
        samples.size(),
        map2.constructMap(),
        map2.constructHasFlip(),
        map2.subMap(),
        map2.subHasFlip(),
        localBest,
        nearestZero,
        nearestEqOp(),
        identityOp(),           // No flipping
        UPstream::msgType(),
        map2.comm()
    );

    // Combine with nearestInfo
    info.resize_nocopy(samples.size());
    forAll(samples, samplei)
    {
        nearestAndDist ni(nearestInfo[samplei]);
        nearestEqOp()(ni, localBest[samplei]);

        info[samplei] = ni.first();
    }
}


void Foam::distributedTriSurfaceMesh::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    const labelList& regionIndices,
    List<pointIndexHit>& info
) const
{
    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);

    if (runningTransparent())
    {
        triSurfaceMesh::findNearest
        (
            samples,
            nearestDistSqr,
            regionIndices,
            info
        );
        return;
    }
    else if (returnReduceAnd(regionIndices.empty(), communicator))
    {
        findNearest(samples, nearestDistSqr, info);
        return;
    }

    addProfiling
    (
        findNearestRegion,
        "distributedTriSurfaceMesh::findNearestRegion"
    );

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::findNearest :"
            << " surface " << searchableSurface::name()
            << " trying to find nearest and region for " << samples.size()
            << " samples with max sphere "
            << (samples.size() ? Foam::sqrt(max(nearestDistSqr)) : Zero)
            << endl;
    }

    // -----------
    // Node master
    // -----------

    if
    (
        (distType_ == distributionType::NODEMASTER)
     && (control_ & controlType::ctrlOverlappingNodeMaster)
    )
    {
        // Service all queries via the node-local master
        nodeMasterQuery
        (
            comm_,
            UPstream::msgType(),
            [this]
            (
                const pointField& s,
                const scalarField& d,
                const labelList& r,
                List<pointIndexHit>& out
            )
            {
                triSurfaceMesh::findNearest(s, d, r, out);
            },
            info,
            samples, nearestDistSqr, regionIndices
        );
        return;
    }

    if
    (
        (distType_ == distributionType::NODEMASTER)
     && !(control_ & controlType::ctrlOverlappingNodeMaster)
    )
    {
        // Gather all inputs, handle on node-local master, scatter result
        constexpr auto commType = UPstream::commsTypes::nonBlocking;
        const int tag = UPstream::msgType();

        const globalIndex gi
        (
            globalIndex::gatherOnly{},
            samples.size(),
            comm_
        );

        const pointField allSamples
        (
            gi.gather(samples, tag, commType, comm_)
        );
        const scalarField allDistSqr
        (
            gi.gather(nearestDistSqr, tag, commType, comm_)
        );
        const labelList allRegionIndices
        (
            gi.gather(regionIndices, tag, commType, comm_)
        );

        List<pointIndexHit> allInfo;
        if (UPstream::master(comm_))
        {
            triSurfaceMesh::findNearest
            (
                allSamples,
                allDistSqr,
                allRegionIndices,
                allInfo
            );
        }

        info.resize_nocopy(samples.size());
        gi.scatter(allInfo, info, tag, commType, comm_);
        return;
    }


    // ---------------
    // Normal handling
    // ---------------


    // Calculate queries and exchange map
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    pointField allCentres;
    scalarField allRadiusSqr;
    labelList allSegmentMap;
    autoPtr<mapDistribute> mapPtr
    (
        calcLocalQueries
        (
            true,      // also send to local processor
            samples,
            nearestDistSqr,

            allCentres,
            allRadiusSqr,
            allSegmentMap
        )
    );
    const mapDistribute& map = mapPtr();


    // swap samples to local processor
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    map.distribute(allCentres);
    map.distribute(allRadiusSqr);


    // Do my tests
    // ~~~~~~~~~~~

    List<pointIndexHit> allInfo(allCentres.size());
    triSurfaceMesh::findNearest
    (
        allCentres,
        allRadiusSqr,
        regionIndices,
        allInfo
    );
    convertTriIndices(allInfo);

    const treeBoundBox& myProcBb = procBb_[myProci];

    forAll(allInfo, i)
    {
        if (allInfo[i].hit())
        {
            if
            (
                decomposeUsingBbs_
            &&  surfaceClosed_
            && !myProcBb.contains(allInfo[i].point())
            )
            {
                // Nearest point is not on local processor so the
                // the triangle is only there because some other bit of it
                // is on it. Assume there is another processor that holds
                // the full surrounding of the triangle so we can clear
                // this particular nearest.
                allInfo[i].setMiss();
                allInfo[i].setIndex(-1);
            }
        }
    }


    // Send back results
    // ~~~~~~~~~~~~~~~~~

    map.reverseDistribute(allSegmentMap.size(), allInfo);


    // Extract information
    // ~~~~~~~~~~~~~~~~~~~

    forAll(allInfo, i)
    {
        if (allInfo[i].hit())
        {
            label pointi = allSegmentMap[i];

            if
            (
                info[pointi].not_hit() ||
                (
                    // Nearest intersection
                    samples[pointi].distSqr(allInfo[i].point())
                  < samples[pointi].distSqr(info[pointi].point())
                )
            )
            {
                info[pointi] = allInfo[i];
            }
        }
    }
}


void Foam::distributedTriSurfaceMesh::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    if (runningTransparent())
    {
        triSurfaceMesh::findLine(start, end, info);
        return;
    }

    // -----------
    // Node Master
    // -----------

    if
    (
        (distType_ == distributionType::NODEMASTER)
     && (control_ & controlType::ctrlOverlappingNodeMaster)
    )
    {
        // Service all queries via the node-local master
        nodeMasterQuery
        (
            comm_,
            UPstream::msgType(),
            [this]
            (
                const pointField& s,
                const pointField& e,
                List<pointIndexHit>& out
            )
            {
                triSurfaceMesh::findLine(s, e, out);
            },
            info,
            start, end
        );
        return;
    }

    if
    (
        (distType_ == distributionType::NODEMASTER)
     && !(control_ & controlType::ctrlOverlappingNodeMaster)
    )
    {
        // Gather all inputs, handle on node-local master, scatter result
        constexpr auto commType = UPstream::commsTypes::nonBlocking;
        const int tag = UPstream::msgType();

        const globalIndex gi(globalIndex::gatherOnly{}, start.size(), comm_);

        const pointField allStart
        (
            gi.gather(start, tag, commType, comm_)
        );
        const pointField allEnd
        (
            gi.gather(end, tag, commType, comm_)
        );

        List<pointIndexHit> allInfo;
        if (UPstream::master(comm_))
        {
            triSurfaceMesh::findLine(allStart, allEnd, allInfo);
        }

        info.resize_nocopy(start.size());
        gi.scatter(allInfo, info, tag, commType, comm_);
        return;
    }

    // ----------------
    // Regular handling
    // ----------------

    constexpr bool nearestIntersection = true;
    findLine(nearestIntersection, start, end, info);
}


void Foam::distributedTriSurfaceMesh::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    if (runningTransparent())
    {
        triSurfaceMesh::findLineAny(start, end, info);
        return;
    }

    // -----------
    // Node Master
    // -----------

    if
    (
        (distType_ == distributionType::NODEMASTER)
     && (control_ & controlType::ctrlOverlappingNodeMaster)
    )
    {
        // Service all queries via the node-local master
        nodeMasterQuery
        (
            comm_,
            UPstream::msgType(),
            [this]
            (
                const pointField& s,
                const pointField& e,
                List<pointIndexHit>& out
            )
            {
                triSurfaceMesh::findLineAny(s, e, out);
            },
            info,
            start, end
        );
        return;
    }

    if
    (
        (distType_ == distributionType::NODEMASTER)
     && !(control_ & controlType::ctrlOverlappingNodeMaster)
    )
    {
        // Gather all inputs, handle on node-local master, scatter result
        constexpr auto commType = UPstream::commsTypes::nonBlocking;
        const int tag = UPstream::msgType();

        const globalIndex gi(globalIndex::gatherOnly{}, start.size(), comm_);

        const pointField allStart
        (
            gi.gather(start, tag, commType, comm_)
        );
        const pointField allEnd
        (
            gi.gather(end, tag, commType, comm_)
        );

        List<pointIndexHit> allInfo;
        if (UPstream::master(comm_))
        {
            triSurfaceMesh::findLineAny(allStart, allEnd, allInfo);
        }

        info.resize_nocopy(start.size());
        gi.scatter(allInfo, info, tag, commType, comm_);
        return;
    }

    // ----------------
    // Regular handling
    // ----------------

    constexpr bool nearestIntersection = true;

    findLine(nearestIntersection, start, end, info);
}


void Foam::distributedTriSurfaceMesh::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    if (runningTransparent())
    {
        triSurfaceMesh::findLineAll(start, end, info);
        return;
    }

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::findLineAll :"
            << " surface " << searchableSurface::name()
            << " intersecting with "
            << start.size() << " rays" << endl;
    }

    addProfiling
    (
        findLineAll,
        "distributedTriSurfaceMesh::findLineAll"
    );

    // -----------
    // Node Master
    // -----------

    if
    (
        (distType_ == distributionType::NODEMASTER)
     && (control_ & controlType::ctrlOverlappingNodeMaster)
    )
    {
        // Service all queries via the node-local master
        nodeMasterQuery
        (
            comm_,
            UPstream::msgType(),
            [this]
            (
                const pointField& s,
                const pointField& e,
                List<List<pointIndexHit>>& out
            )
            {
                triSurfaceMesh::findLineAll(s, e, out);
            },
            info,
            start, end
        );
        return;
    }

    if
    (
        (distType_ == distributionType::NODEMASTER)
     && !(control_ & controlType::ctrlOverlappingNodeMaster)
    )
    {
        // Gather all inputs, handle on node-local master, scatter result
        constexpr auto commType = UPstream::commsTypes::nonBlocking;
        const int tag = UPstream::msgType();

        const globalIndex gi(globalIndex::gatherOnly{}, start.size(), comm_);

        const pointField allStart
        (
            gi.gather(start, tag, commType, comm_)
        );
        const pointField allEnd
        (
            gi.gather(end, tag, commType, comm_)
        );

        List<List<pointIndexHit>> allInfo;
        if (UPstream::master(comm_))
        {
            triSurfaceMesh::findLineAll(allStart, allEnd, allInfo);
        }

        info.resize_nocopy(start.size());
        gi.scatter(allInfo, info, tag, commType, comm_);
        return;
    }


    // ----------------
    // Regular handling
    // ----------------


    // Reuse fineLine. We could modify all of findLine to do multiple
    // intersections but this would mean a lot of data transferred so
    // for now we just find nearest intersection and retest from that
    // intersection onwards.

    // Work array.
    List<pointIndexHit> hits(start.size());
    {
        constexpr bool nearestIntersection = true;
        findLine(nearestIntersection, start, end, hits);
    }

    // Tolerances:
    // To find all intersections we add a small vector to the last intersection
    // This is chosen such that
    // - it is significant (SMALL is smallest representative relative tolerance;
    //   we need something bigger since we're doing calculations)
    // - if the start-end vector is zero we still progress

    // Copy to input and compact any hits

    DynamicList<label> pointMap(start.size());
    DynamicField<point> e0(start.size());
    DynamicField<point> e1(start.size());

    info.resize_nocopy(hits.size());
    forAll(hits, pointi)
    {
        info[pointi].clear();

        if (hits[pointi].hit())
        {
            info[pointi].resize_fill(1, hits[pointi]);

            const auto& starts = start[pointi];
            const auto& ends = end[pointi];

            // For progression
            const vector dirVec = (ends - starts);
            const scalar magSqrDirVec = dirVec.magSqr();
            const vector smallVec
            (
                ROOTSMALL*dirVec + vector::uniform(ROOTVSMALL)
            );

            point pt = hits[pointi].point() + smallVec;

            // The current coordinate along ray
            const scalar d = ((pt-starts) & dirVec);

            if (d <= magSqrDirVec)
            {
                e0.push_back(pt);
                e1.push_back(ends);
                pointMap.push_back(pointi);
            }
        }
    }


    label iter = 0;
    while (returnReduceOr(e0.size()))
    {
        constexpr bool nearestIntersection = true;
        findLine(nearestIntersection, e0, e1, hits);

        // Extract
        label compacti = 0;
        forAll(hits, i)
        {
            if (hits[i].hit())
            {
                label pointi = pointMap[i];

                const auto& starts = start[pointi];
                const auto& ends = end[pointi];

                // For progression
                const vector dirVec = (ends - starts);
                const scalar magSqrDirVec = dirVec.magSqr();
                const vector smallVec
                (
                    ROOTSMALL*dirVec + vector::uniform(ROOTVSMALL)
                );

                // push_back, even although it is List, not DynamicList!
                info[pointi].push_back(hits[i]);

                point pt = hits[i].point() + smallVec;

                // The current coordinate along ray
                const scalar d = ((pt-starts) & dirVec);

                // Note check for d>0. Very occasionally the octree will find
                // an intersection to the left of the ray due to tolerances.
                if (d > 0 && d <= magSqrDirVec)
                {
                    e0[compacti] = pt;
                    e1[compacti] = ends;
                    pointMap[compacti] = pointi;
                    ++compacti;
                }
            }
        }

        // Trim
        e0.resize(compacti);
        e1.resize(compacti);
        pointMap.resize(compacti);

        iter++;

        if (iter == 1000)
        {
            // pointMap holds indices into start/end. e0/e1 are already
            // compacted to match pointMap's size, so they are dumped
            // directly. (e1 is omitted: it is end[pointMap] by
            // construction.)
            if (!e0.empty())
            {
                Pout<< "distributedTriSurfaceMesh::findLineAll :"
                    << " Exiting loop due to excessive number of"
                    << " intersections along " << e0.size() << " rays"
                    << " start:" << UIndirectList<point>(start, pointMap)
                    << " end:" << UIndirectList<point>(end, pointMap)
                    << " stuckAt:" << e0
                    << endl;
            }
            break;
        }
    }

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::findLineAll :"
            << " surface " << searchableSurface::name()
            << " finished intersecting with "
            << start.size() << " rays" << endl;
    }
}


void Foam::distributedTriSurfaceMesh::getRegion
(
    const UList<pointIndexHit>& info,
    labelList& region
) const
{
    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::getRegion :"
            << " surface " << searchableSurface::name()
            << " getting region for "
            << info.size() << " triangles" << endl;
    }

    addProfiling(getRegion, "distributedTriSurfaceMesh::getRegion");

    if (runningTransparent())
    {
        triSurfaceMesh::getRegion(info, region);

        if (debug)
        {
            Pout<< "distributedTriSurfaceMesh::getRegion :"
                << " surface " << searchableSurface::name()
                << " finished getting region for "
                << info.size() << " triangles" << endl;
        }
        return;
    }


    // -----------
    // Node Master
    // -----------

    if
    (
        (distType_ == distributionType::NODEMASTER)
     && (control_ & controlType::ctrlOverlappingNodeMaster)
    )
    {
        // Service all queries via the node-local master
        // - simplify the data exchange by pre-extracting all hit indices

        labelList indices(info.size());
        std::transform
        (
            info.begin(), info.end(),
            indices.begin(),
            [](const auto& x) { return x.hit() ? x.index() : -1; }
        );

        nodeMasterQuery
        (
            comm_,
            UPstream::msgType(),
            [this](const labelUList& i, labelList& out)
            {
                triSurface::regionIDs(i, out);
            },
            region,
            indices
        );
        return;
    }

    if
    (
        (distType_ == distributionType::NODEMASTER)
     && !(control_ & controlType::ctrlOverlappingNodeMaster)
    )
    {
        // Gather all inputs, handle on node-local master, scatter results
        // - simplify data exchange by pre-extracting all hit indices

        labelList indices(info.size());
        std::transform
        (
            info.begin(), info.end(),
            indices.begin(),
            [](const auto& x) { return x.hit() ? x.index() : -1; }
        );

        const label localSize = indices.size();

        constexpr auto commType = UPstream::commsTypes::nonBlocking;
        const int tag = UPstream::msgType();

        const globalIndex gi(globalIndex::gatherOnly{}, localSize, comm_);

        const labelList allIndices(gi.gather(indices, tag, commType, comm_));

        labelList allRegion;
        if (UPstream::master(comm_))
        {
            triSurface::regionIDs(allIndices, allRegion);
        }

        region.resize_nocopy(localSize);
        gi.scatter(allRegion, region, tag, commType, comm_);
        return;
    }


    // ----------------
    // Regular handling
    // ----------------

    // Get query data (= local index of triangle)
    // ~~~~~~~~~~~~~~

    labelList triangleIndex(info.size());
    autoPtr<mapDistribute> mapPtr
    (
        localQueries
        (
            info,
            triangleIndex
        )
    );
    const mapDistribute& map = mapPtr();


    // Extract the local regions
    triSurface::regionIDs(triangleIndex, region);


    // Send back results
    // ~~~~~~~~~~~~~~~~~

    map.reverseDistribute(info.size(), region);

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::getRegion :"
            << " surface " << searchableSurface::name()
            << " finished getting region for "
            << info.size() << " triangles" << endl;
    }
}


void Foam::distributedTriSurfaceMesh::getNormal
(
    const UList<pointIndexHit>& info,
    vectorField& normal
) const
{
    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::getNormal :"
            << " surface " << searchableSurface::name()
            << " getting normal for "
            << info.size() << " triangles" << endl;
    }

    addProfiling(getNormal, "distributedTriSurfaceMesh::getNormal");

    if (runningTransparent())
    {
        triSurfaceMesh::getNormal(info, normal);
        return;
    }

    // -----------
    // Node master
    // -----------

    if
    (
        (distType_ == distributionType::NODEMASTER)
     && (control_ & controlType::ctrlOverlappingNodeMaster)
    )
    {
        // Service all queries via the node-local master
        nodeMasterQuery
        (
            comm_,
            UPstream::msgType(),
            [this](const UList<pointIndexHit>& i, vectorField& out)
            {
                triSurfaceMesh::getNormal(i, out);
            },
            normal,
            info
        );
        return;
    }

    if
    (
        (distType_ == distributionType::NODEMASTER)
     && !(control_ & controlType::ctrlOverlappingNodeMaster)
    )
    {
        // Gather all inputs, handle on node-local master, scatter results
        constexpr auto commType = UPstream::commsTypes::nonBlocking;
        const int tag = UPstream::msgType();

        const globalIndex gi(globalIndex::gatherOnly{}, info.size(), comm_);

        const List<pointIndexHit> allInfo
        (
            gi.gather(info, tag, commType, comm_)
        );

        vectorField allNormal;
        if (UPstream::master(comm_))
        {
            // Do findNearest on master
            triSurfaceMesh::getNormal(allInfo, allNormal);
            // if (vertexNormals_)
            // {
            //     // Use smooth interpolation for the normal
            //     forAll(triangleIndex, i)
            //     {
            //         const point& nearest = allInfo[i].point();
            //         const label trii = allInfo[i].index();
            //         const triPointRef tri(s[trii].tri(s.points()));
            //         const auto& vn = vertexNormals_()[trii];
            //         const barycentric2D w(tri.pointToBarycentric(nearest));
            //         normal[i] = w[0]*vn[0]+w[1]*vn[1]+w[2]*vn[2];
            //     }
            // }
            // else
            // {
            //     forAll(triangleIndex, i)
            //     {
            //         const label trii = triangleIndex[i];
            //         normal[i] = s[trii].unitNormal(s.points());
            //     }
            // }
        }

        normal.resize_nocopy(info.size());
        gi.scatter(allNormal, normal, tag, commType, comm_);
        return;
    }


    // Get query data (= local index of triangle)
    // ~~~~~~~~~~~~~~

    labelList triangleIndex(info.size());
    autoPtr<mapDistribute> mapPtr
    (
        localQueries
        (
            info,
            triangleIndex
        )
    );
    const mapDistribute& map = mapPtr();


    // Do my tests
    // ~~~~~~~~~~~

    const triSurface& s = *this;

    normal.resize_nocopy(triangleIndex.size());

    if (vertexNormals_)
    {
        // Use smooth interpolation for the normal

        // Send over nearest point since we have it anyway
        List<pointIndexHit> localNearestInfo(info);
        map.distribute(localNearestInfo);

        forAll(triangleIndex, i)
        {
            const point& nearest = localNearestInfo[i].point();

            const label trii = triangleIndex[i];
            const triPointRef tri(s[trii].tri(s.points()));
            const auto& vn = vertexNormals_()[trii];
            const barycentric2D w(tri.pointToBarycentric(nearest));
            normal[i] = w[0]*vn[0]+w[1]*vn[1]+w[2]*vn[2];
        }
    }
    else
    {
        forAll(triangleIndex, i)
        {
            const label trii = triangleIndex[i];
            normal[i] = s[trii].unitNormal(s.points());
        }
    }


    // Send back results
    // ~~~~~~~~~~~~~~~~~

    map.reverseDistribute(info.size(), normal);

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::getNormal :"
            << " surface " << searchableSurface::name()
            << " finished getting normal for "
            << info.size() << " triangles" << endl;
    }
}


Foam::label Foam::distributedTriSurfaceMesh::measure
(
    List<vector>* centres,
    List<scalar>* radiusSqr,
    List<label>* regions,
    List<vector>* normals
) const
{
    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::measure :"
            << " surface " << searchableSurface::name() << endl;
    }

    if
    (
        (distType_ == distributionType::NODEMASTER)
     && (!runningTransparent() && !UPstream::master(comm_))
    )
    {
        // No surface on sub-ranks
        // Same as searchableSurface::measure()

        if (centres) centres->clear();
        if (radiusSqr) radiusSqr->clear();
        if (regions) regions->clear();
        if (normals) normals->clear();
        return 0;
    }

    return triSurfaceMesh::measure(centres, radiusSqr, regions, normals);
}


//void Foam::distributedTriSurfaceMesh::getVolumeTypeUncached
//(
//    const pointField& samples,
//    List<volumeType>& volType
//) const
//{
//    if (runningTransparent())
//    {
//        triSurfaceMesh::getVolumeType(samples, volType);
//        return;
//    }
//
//
//    if (!hasVolumeType())
//    {
//        FatalErrorInFunction
//            << "Volume type only supported for closed distributed surfaces."
//            << exit(FatalError);
//    }
//
//    // Trigger (so parallel synchronised) construction of outside type.
//    // Normally this would get triggered from inside individual searches
//    // so would not be parallel synchronised
//    if (outsideVolType_ == volumeType::UNKNOWN)
//    {
//        // Determine nearest (in parallel)
//        const point outsidePt(bounds().max() + 0.5*bounds().span());
//        if (debug)
//        {
//            Pout<< "distributedTriSurfaceMesh::outsideVolumeType :"
//                << " triggering outsidePoint" << outsidePt
//                << " orientation" << endl;
//        }
//
//        const pointField outsidePts(1, outsidePt);
//        List<pointIndexHit> nearestInfo;
//        findNearest
//        (
//            outsidePts,
//            scalarField(1, Foam::sqr(GREAT)),
//            nearestInfo
//        );
//
//        List<volumeType> outsideVolTypes;
//        surfaceSide(outsidePts, nearestInfo, outsideVolTypes);
//        outsideVolType_ = outsideVolTypes[0];
//
//        if (debug)
//        {
//            Pout<< "distributedTriSurfaceMesh::outsideVolumeType :"
//                << " determined outsidePoint" << outsidePt
//                << " to be " << volumeType::names[outsideVolType_] << endl;
//        }
//    }
//
//    // Determine nearest (in parallel)
//    List<pointIndexHit> nearestInfo(samples.size());
//    findNearest
//    (
//        samples,
//        scalarField(samples.size(), Foam::sqr(GREAT)),
//        nearestInfo
//    );
//
//    // Determine orientation (in parallel)
//    surfaceSide(samples, nearestInfo, volType);
//}


void Foam::distributedTriSurfaceMesh::getVolumeType
(
    const pointField& samples,
    List<volumeType>& volType
) const
{
    const int communicator = UPstream::worldComm;
    const auto numProc = UPstream::nProcs(communicator);

    if (runningTransparent())
    {
        triSurfaceMesh::getVolumeType(samples, volType);
        return;
    }

    // -----------
    // Node Master
    // -----------

    if
    (
        (distType_ == distributionType::NODEMASTER)
     && (control_ & controlType::ctrlOverlappingNodeMaster)
    )
    {
        // Service all queries via the node-local master
        nodeMasterQuery
        (
            comm_,
            UPstream::msgType(),
            [this](const pointField& s, List<volumeType>& out)
            {
                triSurfaceMesh::getVolumeType(s, out);
            },
            volType,
            samples
        );
        return;
    }

    if
    (
        (distType_ == distributionType::NODEMASTER)
     && !(control_ & controlType::ctrlOverlappingNodeMaster)
    )
    {
        // Gather all inputs, handle on node-local master, scatter results
        constexpr auto commType = UPstream::commsTypes::nonBlocking;
        const int tag = UPstream::msgType();

        const globalIndex gi(globalIndex::gatherOnly{}, samples.size(), comm_);

        const pointField allSamples
        (
            gi.gather(samples, tag, commType, comm_)
        );

        List<volumeType> allVolType;
        if (UPstream::master(comm_))
        {
            triSurfaceMesh::getVolumeType(allSamples, allVolType);
        }

        volType.resize_nocopy(samples.size());
        gi.scatter(allVolType, volType, tag, commType, comm_);
        return;
    }

    // ----------------
    // Regular handling
    // ----------------

    if (!hasVolumeType())
    {
        FatalErrorInFunction
            << "Volume type only supported for closed distributed surfaces."
            << exit(FatalError);
    }

    // Trigger (so parallel synchronised) construction of outside type.
    // Normally this would get triggered from inside individual searches
    // so would not be parallel synchronised
    if (outsideVolType_ == volumeType::UNKNOWN)
    {
        addProfiling
        (
            getVolumeType,
            "distributedTriSurfaceMesh::getCachedVolumeType"
        );

        // Determine nearest (in parallel)
        const point outsidePt(bounds().max() + 0.5*bounds().span());
        if (debug)
        {
            Pout<< "distributedTriSurfaceMesh::getVolumeType :"
                << " surface " << searchableSurface::name()
                << " triggering outsidePoint" << outsidePt
                << " orientation" << endl;
        }

        const pointField outsidePts(1, outsidePt);
        List<pointIndexHit> nearestInfo;
        findNearest
        (
            outsidePts,
            scalarField(1, Foam::sqr(GREAT)),
            nearestInfo
        );

        List<volumeType> outsideVolTypes;
        surfaceSide(outsidePts, nearestInfo, outsideVolTypes);

        // All processors (that have enough surface) will return the same
        // status, all others will return UNKNOWN. Make INSIDE/OUTSIDE win.
        outsideVolType_ = volumeType
        (
            returnReduce(int(outsideVolTypes[0]), maxOp<int>())
        );

        if (debug)
        {
            Pout<< "distributedTriSurfaceMesh::getVolumeType :"
                << " surface " << searchableSurface::name()
                << " determined outsidePoint" << outsidePt
                << " to be " << volumeType::names[outsideVolType_] << endl;
        }

        if
        (
            outsideVolType_ == volumeType::INSIDE
         || outsideVolType_ == volumeType::OUTSIDE
        )
        {
            const indexedOctree<treeDataTriSurface>& t = tree();
            const auto& nodes = t.nodes();
            PackedList<2>& nt = t.nodeTypes();
            nt.resize(nodes.size());
            nt = volumeType::UNKNOWN;

            if (!decomposeUsingBbs_)
            {
                markMixedOverlap(t, nt);
            }

            cacheVolumeType(nt);
        }
    }


    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::getVolumeType :"
            << " surface " << searchableSurface::name()
            << " finding orientation for " << samples.size()
            << " samples" << endl;
    }

    addProfiling
    (
        getVolumeType,
        "distributedTriSurfaceMesh::getVolumeType"
    );


    DynamicList<label> outsideSamples;

    // Distribute samples to relevant processors
    autoPtr<mapDistribute> mapPtr;
    {
        List<DynamicList<label>> dynSendMap(numProc);
        DynamicList<label> hits;

        forAll(samples, samplei)
        {
            // Find the processors this sample overlaps.
            if (procBb_.find_contains(samples[samplei], hits))
            {
                for (auto proci : hits)
                {
                    dynSendMap[proci].push_back(samplei);
                }
            }
            else
            {
                // Special case: point is outside all bbs.
                // These would not get sent to anyone so handle locally.
                // Note that is the equivalent of the test in triSurfaceMesh
                // against the local tree bb

                outsideSamples.push_back(samplei);
            }
        }

        // From dynamic to regular storage
        labelListList sendMap(numProc);
        forAll(sendMap, proci)
        {
            sendMap[proci].transfer(dynSendMap[proci]);
        }

        mapPtr.emplace(std::move(sendMap));
    }
    const mapDistribute& map = *mapPtr;

    // Get the points I need to test
    pointField localPoints(samples);
    map.distribute(localPoints);

    volType.resize_fill(localPoints.size(), volumeType::UNKNOWN);

    // Split the local queries into those that I can look up on the tree and
    // those I need to search the nearest for
    DynamicField<point> fullSearchPoints(localPoints.size());
    DynamicList<label> fullSearchMap(localPoints.size());

    forAll(localPoints, i)
    {
        // If no-fill in:
        // - tree nodes on the outside will not include fill-in triangles
        // - so might decide the wrong state (might be 'mixed' with fill-in
        //   triangles)
        // - for now do not cache
        // - TBD. Note that the cached volume types already include remote
        //   fill-in triangles. There is still a problem there though. TBD.

        if (decomposeUsingBbs_ && tree().nodes().size())
        {
            volType[i] = cachedVolumeType(0, localPoints[i]);
        }
        if (volType[i] == volumeType::UNKNOWN)
        {
            fullSearchMap.append(i);
            fullSearchPoints.append(localPoints[i]);
        }
    }

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::getVolumeType :"
            << " surface " << searchableSurface::name()
            << " original samples:" << samples.size()
            << " resulting in local queries:"
            << localPoints.size() << " of which cached:"
            << localPoints.size()-fullSearchPoints.size()
            << endl;
    }

    // Determine nearest (in parallel)
    List<pointIndexHit> nearestInfo;
    findNearest
    (
        fullSearchPoints,
        scalarField(fullSearchPoints.size(), Foam::sqr(GREAT)),
        nearestInfo
    );

    // Determine orientation (in parallel)
    List<volumeType> fullSearchType;
    surfaceSide(fullSearchPoints, nearestInfo, fullSearchType);

    // Combine
    forAll(fullSearchMap, i)
    {
        volumeType& vt = volType[fullSearchMap[i]];
        if (vt == volumeType::UNKNOWN)
        {
            vt = fullSearchType[i];
        }
        else if (vt != fullSearchType[i])
        {
            Pout<< "At point:" << fullSearchPoints[i]
                << " or point:" << localPoints[fullSearchMap[i]]
                << " have cached status:" << vt
                << " have full-search status:" << fullSearchType[i]
                << endl;
        }
        //volType[fullSearchMap[i]] = fullSearchType[i];
    }

    //{
    //    //- Variant of below that produces better error messages
    //    typedef Tuple2<Pair<point>, volumeType> NearType;
    //    const NearType zero(Pair<point>(Zero, Zero), volumeType::UNKNOWN);
    //
    //    // Combine nearest and volType so we can give nicer error messages
    //    List<NearType> nearTypes(volType.size(), zero);
    //
    //    NearTypeCombineOp cop;
    //
    //    // For all volType we do have the originating point but not the
    //    // nearest
    //    forAll(localPoints, i)
    //    {
    //        nearTypes[i] =
    //            NearType(Pair<point>(localPoints[i], Zero), volType[i]);
    //    }
    //    // But for all full-search we also have the nearest
    //    forAll(fullSearchMap, i)
    //    {
    //        const NearType nt
    //        (
    //            Pair<point>(fullSearchPoints[i], nearestInfo[i].point()),
    //            fullSearchType[i]
    //        );
    //
    //        cop(nearTypes[fullSearchMap[i]], nt);
    //    }
    //    mapDistributeBase::distribute
    //    (
    //        UPstream::commsTypes::nonBlocking,
    //        UList<labelPair>::null(),
    //        samples.size(),
    //        map.constructMap(),
    //        map.constructHasFlip(),
    //        map.subMap(),
    //        map.subHasFlip(),
    //        nearTypes,
    //        zero,
    //        cop,                //NearTypeCombineOp(),
    //        identityOp(),           // No flipping
    //        UPstream::msgType(),
    //        map.comm()
    //    );
    //}

    // Send back to originator. In case of multiple answers choose inside or
    // outside
    const volumeType zero(volumeType::UNKNOWN);
    mapDistributeBase::distribute
    (
        UPstream::commsTypes::nonBlocking,
        UList<labelPair>::null(),
        samples.size(),
        map.constructMap(),
        map.constructHasFlip(),
        map.subMap(),
        map.subHasFlip(),
        volType,
        zero,
        volumeCombineOp(),
        identityOp(),           // No flipping
        UPstream::msgType(),
        map.comm()
    );

    // Add the points outside the bounding box
    for (label samplei : outsideSamples)
    {
        volType[samplei] = outsideVolType_;
    }

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::getVolumeType :"
            << " surface " << searchableSurface::name()
            << " finished finding orientation for " << samples.size()
            << " samples" << endl;
    }
}


Foam::label Foam::distributedTriSurfaceMesh::getField(labelList& values) const
{
    if
    (
        (distType_ == distributionType::NODEMASTER)
     && (!runningTransparent() && !UPstream::master(comm_))
    )
    {
        // This is wildly inconsistent!
        // triSurfaceMesh treats a missing field as "ignore"
        // (does not touch the values).
        // But if the field does exist and either surface or field are
        // zero length, it will clear the values.

        // So for now pretend that this sub-rank has a field
        // but that it has 0 values (and obviously also has 0 faces)

        values.clear();
        return 0;
    }

    return triSurfaceMesh::getField(values);
}


void Foam::distributedTriSurfaceMesh::getField
(
    const UList<pointIndexHit>& info,
    labelList& values
) const
{
    if (runningTransparent())
    {
        triSurfaceMesh::getField(info, values);
        return;
    }

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::getField :"
            << " surface " << searchableSurface::name()
            << " retrieving field for "
            << info.size() << " triangles" << endl;
    }

    addProfiling(getField, "distributedTriSurfaceMesh::getField");

    // -----------
    // Node master
    // -----------

    // Simple test and exit if missing
    if (distType_ == distributionType::NODEMASTER)
    {
        int hasField =
        (
            // ie, bool value
            UPstream::master(comm_) && triSurfaceMesh::hasValuesField()
        );

        Pstream::broadcast(hasField, comm_);

        if (!hasField)
        {
            return;
        }
    }

    if
    (
        (distType_ == distributionType::NODEMASTER)
     && (control_ & controlType::ctrlOverlappingNodeMaster)
    )
    {
        // Service all queries via the node-local master
        // - simplify the data exchange by pre-extracting all hit indices

        // Already checked above that hasValuesField is true on master...

        labelList indices(info.size());
        std::transform
        (
            info.begin(), info.end(),
            indices.begin(),
            [](const auto& x) { return x.hit() ? x.index() : -1; }
        );

        nodeMasterQuery
        (
            comm_,
            UPstream::msgType(),
            [this](const labelUList& i, labelList& out)
            {
                triSurfaceMesh::getField_at(i, out);
            },
            values,
            indices
        );
        return;
    }

    if
    (
        (distType_ == distributionType::NODEMASTER)
     && !(control_ & controlType::ctrlOverlappingNodeMaster)
    )
    {
        // Gather all inputs, handle on node-local master, scatter results
        // - simplify data exchange by pre-extracting all hit indices

        // Already checked above that hasValuesField is true on master...

        labelList indices(info.size());
        std::transform
        (
            info.begin(), info.end(),
            indices.begin(),
            [](const auto& x) { return x.hit() ? x.index() : -1; }
        );

        constexpr auto commType = UPstream::commsTypes::nonBlocking;
        const int tag = UPstream::msgType();

        const globalIndex gi(globalIndex::gatherOnly{}, indices.size(), comm_);

        const labelList allIndices(gi.gather(indices, tag, commType, comm_));

        labelList allValues;
        if (UPstream::master(comm_))
        {
            triSurfaceMesh::getField_at(allIndices, allValues);
        }

        values.resize(indices.size());
        gi.scatter(allValues, values, tag, commType, comm_);
        return;
    }


    // ----------------
    // Regular handling
    // ----------------

    if (triSurfaceMesh::hasValuesField())
    {
        // Get query data (= local index of triangle)
        // ~~~~~~~~~~~~~~

        labelList triangleIndex(info.size());
        autoPtr<mapDistribute> mapPtr
        (
            localQueries
            (
                info,
                triangleIndex
            )
        );
        const mapDistribute& map = mapPtr();

        // Do my tests
        // ~~~~~~~~~~~

        triSurfaceMesh::getField_at(triangleIndex, values);

        // Send back results
        // ~~~~~~~~~~~~~~~~~

        map.reverseDistribute(info.size(), values);
    }

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::getField :"
            << " surface " << searchableSurface::name()
            << " finished retrieving field for "
            << info.size() << " triangles" << endl;
    }
}


// Exchanges indices to the processor they come from.
// - calculates exchange map
// - uses map to calculate local triangle index
Foam::autoPtr<Foam::mapDistribute>
Foam::distributedTriSurfaceMesh::localQueries
(
    const UList<pointIndexHit>& info,
    labelList& triangleIndex
) const
{
    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);
    const auto numProc = UPstream::nProcs(communicator);

    triangleIndex.resize_fill(info.size(), -1);

    const globalIndex& triIndexer = globalTris();


    // Determine send map
    // ~~~~~~~~~~~~~~~~~~

    // Since determining which processor the query should go to is
    // cheap we do a multi-pass algorithm to save some memory temporarily.

    labelListList sendMap(numProc);

    {
        // 1. Count. Also cache the target proc in triangleIndex to avoid
        // repeating whichProcID() queries.

        labelList sendCount(numProc, Foam::zero{});

        forAll(info, i)
        {
            if (info[i].hit())
            {
                auto proci = triIndexer.whichProcID(myProci, info[i].index());

                ++sendCount[proci];

                // proc index [temporary cache]
                triangleIndex[i] = proci;
            }
        }

        // 2. Size sendMap
        forAll(sendCount, proci)
        {
            sendMap[proci].resize_nocopy(sendCount[proci]);
            sendCount[proci] = 0;
        }

        // 3. Fill sendMap. The triangle proc info already filled in stage 1
        forAll(info, i)
        {
            // recover proc index from [temporary cache]
            auto proci = triangleIndex[i];
            auto globali = info[i].index();

            if (proci >= 0 && info[i].hit())
            {
                sendMap[proci][sendCount[proci]++] = i;

                triangleIndex[i] = triIndexer.toLocal(proci, globali);
            }
        }
    }

    // Pack into distribution map

    // auto mapPtr = autoPtr<mapDistribute>::New
    // (
    //     std::move(sendMap),
    //     false, false,  // no flip maps
    //     communicator
    // );

    auto mapPtr = autoPtr<mapDistribute>::New(std::move(sendMap));

    // Send over queries
    mapPtr().distribute(triangleIndex);

    return mapPtr;
}


bool Foam::distributedTriSurfaceMesh::getDistribution
(
    const distributionType distType,
    const UList<treeBoundBox>& followBb,
    labelList& distribution,
    // Per processor the bounding box of all (destination) triangles
    List<treeBoundBox>& newProcBb
) const
{
    if (!UPstream::parRun())
    {
        // No meaning for non-parallel
        newProcBb.clear();
        distribution.clear();
        return false;
    }

    const int communicator = UPstream::worldComm;
    const auto numProc = UPstream::nProcs(communicator);

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::getDistribution :"
            << " surface " << searchableSurface::name()
            << " distributing surface according to method:"
            << distributionTypeNames_[distType_]
            << " fill-in:" << decomposeUsingBbs_
            << " follow bb:" << flatOutput(followBb) << endl;
    }

    addProfiling
    (
        get_distribution,
        "distributedTriSurfaceMesh::getDistribution"
    );

    // Get bbs of all domains
    // ~~~~~~~~~~~~~~~~~~~~~~

    switch (distType)
    {
        // Nothing to do
        case distributionType::FROZEN :
        {
            // No updates, no distribution applicable
            newProcBb.clear();
            distribution.clear();

            return false;
            break;
        }

        case distributionType::FOLLOW :
        {
            // Update bounding boxes
            newProcBb.resize_nocopy(numProc);

            treeBoundBox localBb;
            for (const auto& bb : followBb)
            {
                localBb.add(bb);
            }
            boundBox::allGather(newProcBb, localBb, communicator);

            // No distribution applicable
            distribution.clear();

            return false;
            break;
        }

        case distributionType::INDEPENDENT :
        case distributionType::DISTRIBUTED :
        case distributionType::FULLTREE :
        case distributionType::SPARSETREE :
        // This will not work directly with NODEMASTER ...
        {
            // The destination processor per triangle
            calcDistribution(*this, distribution);

            // Update bounding boxes.
            // Initialize with inverted and add contributions
            newProcBb.resize_fill(numProc, treeBoundBox::null());

            // Find bounding box for all triangles on new distribution.
            if (debug)
            {
                Pout<< "distributedTriSurfaceMesh::getDistribution :"
                    << " determining processor bounding boxes for "
                    << searchableSurface::name() << endl;
            }

            if (const triSurface& surf = *this; !distribution.empty())
            {
                forAll(surf, facei)
                {
                    auto proci = distribution[facei];
                    newProcBb[proci].add(surf.points(), surf[facei]);
                }

                // Inflate a bit,
                // for if triangle coincide with processor boundaries
                for (auto& bb : newProcBb)
                {
                    bb.grow(mergeDist_);
                }
            }

            // Same size on all ranks: reduce all together.
            // This is not an allGather operation.
            boundBox::reduce(newProcBb.data(), newProcBb.size(), communicator);

            return true;
            break;
        }

        case distributionType::TRANSPARENT :
        case distributionType::NODEMASTER :
        {
            // No updates, no distribution applicable
            newProcBb.clear();
            distribution.clear();

            FatalErrorInFunction
                << "distributionType <"
                <<  distributionTypeNames_[distType_]
                << "> not supported in getDistribution."
                << exit(FatalError);
            break;
        }

        default :
        {
            // No updates, no distribution applicable
            newProcBb.clear();
            distribution.clear();

            FatalErrorInFunction
                << "Unsupported distribution type:" << int(distType)
                << exit(FatalError);
            break;
        }
    }

    return false;
}


// side-effects:
// - adjusts procBb_ and updates its entry in the dictionary

void Foam::distributedTriSurfaceMesh::distribute
(
    const UList<treeBoundBox>& followBb,
    const bool keepNonLocal,
    autoPtr<mapDistribute>& faceMap,
    autoPtr<mapDistribute>& pointMap
)
{
    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);
    const auto numProc = UPstream::nProcs(communicator);

    if (!UPstream::parRun())
    {
        return;
    }

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::distribute :"
            << " surface " << searchableSurface::name()
            << " distributing surface according to method:"
            << distributionTypeNames_[distType_]
            << " fill-in:" << decomposeUsingBbs_
            << " follow bb:" << flatOutput(followBb) << endl;
    }

    addProfiling(distribute, "distributedTriSurfaceMesh::distribute");


    // Debug information
    if (debug)
    {
        const label numFaces = triSurface::size();

        labelList nTris
        (
            UPstream::listGatherValues<label>(numFaces, communicator)
        );

        if (UPstream::master(communicator))
        {
            InfoInFunction
                << "before distribution:" << endl << "\tproc\ttris" << endl;

            forAll(nTris, proci)
            {
                Info<< '\t' << proci << '\t' << nTris[proci] << endl;
            }
            Info<< endl;
        }
    }


    // Use distType_ to decide whichfaces go where
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelListList faceSendMap(numProc);

    if (distType_ == distributionType::NODEMASTER && !runningTransparent())
    {
        // Send from master to all node masters

        if (distType_ == currentDistType_)
        {
            return;
        }

        // Clear out old information
        List<treeBoundBox> bbs(procBb_.release_boxes());
        bbs.resize_fill(numProc, treeBoundBox::null());

        boundBox masterBounds;
        if (UPstream::master(communicator))
        {
            masterBounds = triSurfaceMesh::bounds();
        }
        masterBounds.broadcast(communicator);

        // All node-masters get the same bounding box
        for (int proci : masterProcIDs_)
        {
            bbs[proci] = masterBounds;
        }
        dict_.set("bounds", bbs[myProci]);

        procBb_.setBoxes(std::move(bbs));
        procBb_.rebuild();

        // Send actual surface from world-master to node-masters
        if (UPstream::master(communicator))
        {
            for (int proci : masterProcIDs_)
            {
                faceSendMap[proci] = Foam::identity(triSurface::size());
            }
        }
    }
    else if (runningTransparent())
    {
        // Send from master to all node masters

        if (distType_ == currentDistType_)
        {
            return;
        }

        // Clear out old information
        List<treeBoundBox> bbs(procBb_.release_boxes());
        bbs.resize_fill(numProc, treeBoundBox::null());

        boundBox masterBounds;
        if (UPstream::master(communicator))
        {
            masterBounds = triSurfaceMesh::bounds();
        }
        masterBounds.broadcast(communicator);

        // All nodes get the same bounding box
        bbs[myProci] = masterBounds;
        dict_.set("bounds", masterBounds);

        procBb_.setBoxes(std::move(bbs));
        procBb_.rebuild();

        // Send actual surface from world-master to all nodes
        // (including itself) - a slightly wasteful way to broadcast
        if (UPstream::master(communicator))
        {
            for (label proci = 0; proci < numProc; ++proci)
            {
                faceSendMap[proci] = Foam::identity(triSurface::size());
            }
        }
    }
    else
    {
        // independent/distributed only redistribute when the distribution
        // type itself changes. follow must keep tracking the supplied bbs
        // (usually the mesh bounds).
        if
        (
            (distType_ == currentDistType_) &&
            (
                distType_ == distributionType::INDEPENDENT
             || distType_ == distributionType::DISTRIBUTED
             || distType_ == distributionType::FULLTREE
             || distType_ == distributionType::SPARSETREE
             // distType_ == distributionType::NODETREE
            )
        )
        {
            return;
        }

        // Per triangle the destination processor
        labelList distribution;
        {
            // Per processor the bounding box of all (destination) triangles
            List<treeBoundBox> newProcBb;

            getDistribution
            (
                distType_,
                followBb,
                // Per triangle the destination processor
                // or just -1 everywhere for FOLLOW etc...
                distribution,
                // Per processor the bounding box of all (destination) triangles
                newProcBb
            );

            if (newProcBb.empty() || (newProcBb == procBb_.boxes()))
            {
                // No change
                return;
            }
            else
            {
                // Update proc bound-boxes. Trust the veracity of the new
                // boxes that getDistribution delivers

                if (newProcBb.size() != numProc)
                {
                    FatalErrorInFunction
                        << "getDistribution returned " << newProcBb.size()
                        << " bound boxes instead of " << numProc << endl;
                }

                dict_.set("bounds", newProcBb[myProci]);

                procBb_.clear();
                procBb_.setBoxes(std::move(newProcBb));
                procBb_.rebuild();
            }
        }

        const triSurface& s = *this;

        // Work array for face selection
        bitSet faceSelection;

        if (decomposeUsingBbs_)
        {
            forAll(procBb_, proci)
            {
                triSurfaceTools::findOverlapFaces
                (
                    s,
                    procBb_[proci],   //<- the updated bound box
                    faceSelection
                );

                faceSendMap[proci] = faceSelection.toc();
            }
        }
        else
        {
            // Build the subset maps according to the decomposition colouring
            for (label proci = 0; proci < numProc; ++proci)
            {
                PatchTools::subsetMapFaces
                (
                    s,
                    proci,         // colour
                    distribution,  // colouring
                    faceSendMap[proci]
                );
            }
        }

        if (keepNonLocal && decomposeUsingBbs_)
        {
            // Include in faceSendMap/pointSendMap the triangles that are
            // not mapped to any processor so they stay local.

            // Only possible/relevant when decomposeUsingBbs_ is used
            // and the bounding boxes do not arise from a decomposer.

            // Start with no faces
            faceSelection.resize_fill(s.size(), false);

            // All mapped faces
            forAll(faceSendMap, proci)
            {
                faceSelection.set(faceSendMap[proci]);
            }

            // Flip to get the unmapped (missed) faces
            faceSelection.flip();
            faceSelection.resize(s.size());  // safety

            if (faceSelection.any())
            {
                // Add unmapped into local face selection
                auto& map = faceSendMap[myProci];
                const label oldLen = map.size();
                const label count = static_cast<label>(faceSelection.count());

                map.resize(oldLen + count);

                label i = oldLen;
                for (label missed : faceSelection)
                {
                    map[i++] = missed;
                }

                // TBD: useful to have consistent monotonic order?
                // Foam::sort(map);
            }
        }
    }


    // Exchange surfaces
    // ~~~~~~~~~~~~~~~~~
    {
        // Storage for resulting surface
        List<labelledTri> allTris;
        pointField allPoints;

        if (control_ & controlType::ctrlDirectExchange)
        {
            distributeGeometry
            (
                std::move(faceSendMap),
                faceMap,
                pointMap,
                allTris,
                allPoints
            );
        }
        else
        {
            distributeGeometry_buffered
            (
                std::move(faceSendMap),
                faceMap,
                pointMap,
                allTris,
                allPoints
            );
        }

        auto& s = static_cast<triSurface&>(*this);

        // Construct triSurface. Reuse storage.
        s = triSurface(allTris, patches(), allPoints, true);
    }


    // Clear trees, preserve topological info (surfaceClosed, outsidePointType)
    clearOut();

    // boundBox totalBb(static_cast<const triSurface&>(*this).box());
    // totalBb.reduce();

    // Set the bounds() value to the boundBox of the undecomposed surface
    bounds() = boundBox(points(), true);

    // Update exclusionBb if procBb has changed
    calcExclusionBb();
    currentDistType_ = distType_;

    // Regions stays same
    // Volume type stays same.

    distributeFields<label>(faceMap());
    distributeFields<scalar>(faceMap());
    distributeFields<vector>(faceMap());
    distributeFields<sphericalTensor>(faceMap());
    distributeFields<symmTensor>(faceMap());
    distributeFields<tensor>(faceMap());

    if (vertexNormals_)
    {
        List<FixedList<vector, 3>>& vn = vertexNormals_();
        faceMap().distribute(vn);

        if (debug & 2)
        {
            OBJstream str
            (
                searchableSurface::time().path()
               /searchableSurface::name() + "_vertex_normals.obj"
            );
            const triSurface& s = *this;
            forAll(s, trii)
            {
                const auto& f = s[trii];
                const scalar l = Foam::sqrt(f.mag(s.points()));
                forAll(f, fp)
                {
                    const point& pt = s.points()[f[fp]];
                    const vector& n = vn[trii][fp];
                    str.write(linePointRef(pt, pt+l*n));
                }
            }
            Pout<< "Dumped local vertex normals to " << str.name() << endl;
        }
    }


    if
    (
        (
            outsideVolType_ == volumeType::INSIDE ||
            outsideVolType_ == volumeType::OUTSIDE
        ) &&
        (
            distType_ != distributionType::TRANSPARENT
         && distType_ != distributionType::NODEMASTER
        )
    )
    {
        // Re-build tree & inside/outside for closed surfaces

        const indexedOctree<treeDataTriSurface>& t = tree();
        const auto& nodes = t.nodes();
        PackedList<2>& nt = t.nodeTypes();
        nt.resize(nodes.size());
        nt = volumeType::UNKNOWN;

        // Send over any overlapping but not included triangles
        if (!decomposeUsingBbs_)
        {
            // We do not backfill the bounding box so the local tree might
            // not have all the triangles to cache inside/outside correctly.
            // Do a pre-pass to mark 'mixed' the tree leaves that intersect
            // the triangles of the (would be) fill-in triangles.
            markMixedOverlap(t, nt);
        }

        cacheVolumeType(nt);
    }


    if (debug)
    {
        InfoInFunction<< "after distribution:" << endl;
        printBounds();

        if (debug & 2)
        {
            OBJstream str
            (
                searchableSurface::time().path()
               /searchableSurface::name() + "_after.obj"
            );
            Info<< "Writing local bounding box to " << str.name() << endl;
            {
                str.write(procBb_[myProci], true);  // lines
            }
        }
        if (debug & 2)
        {
            OBJstream str
            (
                searchableSurface::time().path()
               /searchableSurface::name() + "_after_all.obj"
            );
            Info<< "Writing all bounding boxes to " << str.name() << endl;
            forAll(procBb_, proci)
            {
                str.write(procBb_[proci], true);  // lines
            }
        }
    }

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::distribute :"
            << " surface " << searchableSurface::name()
            << " done distributing surface according to method:"
            << distributionTypeNames_[distType_]
            << " follow bb:";
        // if (followBb.good())
        // {
        //     Pout<< followBb << endl;
        // }
        // else
        // {
        //     Pout<< "none" << endl;
        // }
    }
}


bool Foam::distributedTriSurfaceMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::writeObject :"
            << " surface " << searchableSurface::name()
            << " writing surface:" << writeOnProc << endl;
    }

    // Make sure dictionary goes to same directory as surface
    const_cast<fileName&>(dict_.instance()) = searchableSurface::instance();

    // Copy of triSurfaceMesh::writeObject except for the sorting of
    // triangles by region. This is done so we preserve region names,
    // even if locally we have zero triangles.
    {
        fileName fullPath(searchableSurface::objectPath());

        if (!Foam::mkDir(fullPath.path()))
        {
            return false;
        }

        // Important: preserve any zero-sized patches
        triSurface::write(fullPath, true);

        if (!Foam::isFile(fullPath))
        {
            return false;
        }
    }

    // Dictionary needs to be written in ascii - binary output not supported.
    streamOpt.format(IOstreamOption::ASCII);
    bool ok = dict_.writeObject(streamOpt, true);

    if (debug)
    {
        Pout<< "distributedTriSurfaceMesh::writeObject :"
            << " surface " << searchableSurface::name()
            << " done writing surface" << endl;
    }

    return ok;
}


void Foam::distributedTriSurfaceMesh::writeStats(Ostream& os) const
{
    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);

    boundBox bb;
    label numPoints = 0;
    PatchTools::calcBounds(*this, bb, numPoints);
    bb.reduce(communicator);

    label numFaces = triSurface::size();

    Foam::reduce(numFaces , sumOp<label>(), UPstream::msgType(), communicator);
    Foam::reduce(numPoints, sumOp<label>(), UPstream::msgType(), communicator);

    os  << "Triangles      : " << numFaces << nl
        << "Vertices       : " << numPoints << nl
        << "Vertex normals : " << bool(vertexNormals_) << nl
        << "Bounding Box   : " << bb << nl
        << "Fill-in:       : " << decomposeUsingBbs_ << nl;
    if (myProci < exclusionBb_.size())
    {
        os  << "Exclusion Box  : " << exclusionBb_[myProci] << nl;
    }
    os  << "Closed         : " << int(surfaceClosed_) << nl
        << "Outside type   : " << volumeType::names[outsideVolType_] << nl
        << "Distribution   : " << distributionTypeNames_[distType_] << endl;
}


void Foam::distributedTriSurfaceMesh::printBounds() const
{
    const int communicator = UPstream::worldComm;

    const labelList nTris
    (
        UPstream::listGatherValues<label>(triSurface::size(), communicator)
    );

    if (UPstream::master(communicator))
    {
        if (exclusionBb_.size())
        {
            Info<< nl << "\tproc\ttris\tbb\texclusion" << nl;
            forAll(nTris, proci)
            {
                Info<< '\t' << proci << '\t' << nTris[proci]
                    << '\t' << procBb_[proci]
                    << '\t' << exclusionBb_[proci] << nl;
            }
        }
        else
        {
            Info<< nl << "\tproc\ttris\tbb" << nl;
            forAll(nTris, proci)
            {
                Info<< '\t' << proci << '\t' << nTris[proci]
                    << '\t' << procBb_[proci] << nl;
            }
        }
        Info<< endl;
    }
}


// ************************************************************************* //
