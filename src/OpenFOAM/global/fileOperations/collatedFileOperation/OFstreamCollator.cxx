/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
    Copyright (C) 2019-2025 OpenCFD Ltd.
    Copyright (C) 2026 Keysight Technologies
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

#include "OFstreamCollator.H"
#include "OFstream.H"
#include "decomposedBlockData.H"
#include "dictionary.H"
#include "masterUncollatedFileOperation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(OFstreamCollator, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::OFstreamCollator::writeFile
(
    const int communicator,
    const word& objectType,
    const fileName& fName,
    std::string_view localData,
    const UList<int64_t>& recvSizes,
    const UList<std::string_view>& procData,  // optional proc data
    IOstreamOption streamOpt,
    IOstreamOption::atomicType atomic,
    IOstreamOption::appendType append,
    const dictionary& headerEntries
)
{
    if (debug)
    {
        Pout<< "OFstreamCollator : Writing local " << int64_t(localData.size())
            << " bytes to " << fName << " using comm " << communicator
            << " and " << procData.size() << " sub-ranks" << endl;

        forAll(procData, proci)
        {
            Pout<< "    " << proci << " size:"
                << int64_t(procData[proci].size()) << nl;
        }
    }

    autoPtr<OSstream> osPtr;
    if (UPstream::master(communicator))
    {
        Foam::mkDir(fName.path());
        osPtr.reset(new OFstream(atomic, fName, streamOpt, append));
        auto& os = *osPtr;

        if (append == IOstreamOption::NO_APPEND)
        {
            // No IOobject so cannot use IOobject::writeHeader

            // FoamFile
            decomposedBlockData::writeHeader
            (
                os,
                streamOpt,      // streamOpt for container
                objectType,
                "",             // note
                "",             // location (leave empty, otherwise inaccurate)
                fName.name(),   // object name
                headerEntries
            );
        }
    }

    // Assuming threaded writing hides any slowness so we
    // can use scheduled communication to send the data to
    // the master processor in order. However can be unstable
    // for some mpi so default is non-blocking.
    const UPstream::commsTypes myCommsType
    (
        // Blocking when buffer size is 0
        Foam::mag
        (
            fileOperations::masterUncollatedFileOperation::
            maxMasterFileBufferSize
        ) < 1
      ? UPstream::commsTypes::scheduled
      : UPstream::commsTypes::nonBlocking
    );


    List<std::streamoff> blockOffsets;  // Optional
    decomposedBlockData::writeBlocks
    (
        communicator,
        osPtr,
        blockOffsets,  // or List<std::streamoff>::null()
        localData,
        recvSizes,
        procData,
        myCommsType,
        false       // do not sync return state
    );

    if (osPtr && !osPtr->good())
    {
        FatalIOErrorInFunction(*osPtr)
            << "Failed writing to " << fName << exit(FatalIOError);
    }

    if (debug)
    {
        Pout<< "OFstreamCollator : Finished writing "
            << localData.size() << " bytes";

        if (UPstream::master(communicator))
        {
            int64_t total = 0;
            for (auto recv : recvSizes)
            {
                total += int64_t(recv);
            }
            // Use std::to_string to display long int
            Pout<< " (overall " << std::to_string(total) << ')';
        }
        Pout<< " to " << fName
            << " using comm " << communicator << endl;
    }

    return true;
}


void* Foam::OFstreamCollator::writeAll(void *threadarg)
{
    OFstreamCollator& handler = *static_cast<OFstreamCollator*>(threadarg);

    // Consume stack
    while (true)
    {
        std::unique_ptr<writeData> ptr;

        {
            std::lock_guard<std::mutex> guard(handler.mutex_);

            if (handler.objects_.size())
            {
                // FIFO
                ptr.reset(handler.objects_.front());
                handler.objects_.pop_front();
            }
        }

        if (!ptr)
        {
            break;
        }

        writeData& obj = *ptr;

        // Obtain views from storage
        List<std::string_view> procData(obj.procData_.size());
        forAll(procData, proci)
        {
            procData[proci] = obj.procView(proci);
        }

        bool ok = writeFile
        (
            obj.comm_,
            obj.objectType_,
            obj.pathName_,
            obj.localView(),
            obj.sizes_,
            procData,
            obj.streamOpt_,
            obj.atomic_,
            obj.append_,
            obj.headerEntries_
        );

        if (!ok)
        {
            FatalIOErrorInFunction(obj.pathName_)
                << "Failed writing " << obj.pathName_
                << exit(FatalIOError);
        }
        //sleep(1);
    }

    if (debug)
    {
        Pout<< "OFstreamCollator : Exiting write thread " << endl;
    }

    {
        std::lock_guard<std::mutex> guard(handler.mutex_);
        handler.threadRunning_ = false;
    }

    return nullptr;
}


void Foam::OFstreamCollator::waitForBufferSpace(const int64_t wantedSize) const
{
    while (true)
    {
        // The pending output size(s)
        int64_t totalSize = 0;

        {
            std::lock_guard<std::mutex> guard(mutex_);
            for (const writeData* ptr : objects_)
            {
                if (ptr) totalSize += int64_t(ptr->size());
            }
        }

        if
        (
            totalSize == 0
         || (wantedSize >= 0 && (totalSize+wantedSize) <= maxBufferSize_)
        )
        {
            break;
        }

        if (debug)
        {
            std::lock_guard<std::mutex> guard(mutex_);
            Pout<< "OFstreamCollator : Waiting for buffer space."
                << " Currently in use:" << totalSize
                << " limit:" << maxBufferSize_
                << " files:" << objects_.size()
                << endl;
        }

        sleep(5);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OFstreamCollator::OFstreamCollator(int64_t maxBufferSize)
:
    OFstreamCollator(maxBufferSize, UPstream::worldComm)
{}


Foam::OFstreamCollator::OFstreamCollator
(
    int64_t maxBufferSize,
    int communicator
)
:
    maxBufferSize_(maxBufferSize),
    threadRunning_(false),
    localComm_(communicator),
    threadComm_(UPstream::dupCommunicator(localComm_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::OFstreamCollator::~OFstreamCollator()
{
    if (thread_)
    {
        if (debug)
        {
            Pout<< "~OFstreamCollator : Waiting for write thread" << endl;
        }
        thread_->join();
        thread_.reset(nullptr);
    }

    UPstream::freeCommunicator(threadComm_);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// The address of the content if it is of the (native) storage type,
// otherwise nullptr.
template<class Storage, class T>
Storage* ownedStorage(T& content)
{
    if constexpr (std::is_same_v<Storage, T>)
    {
        return &content;
    }
    else
    {
        return nullptr;
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::OFstreamCollator::write
(
    const word& objectType,
    const fileName& fName,
    DynamicList<char>&& localData,
    IOstreamOption streamOpt,
    IOstreamOption::atomicType atomic,
    IOstreamOption::appendType append,
    const bool useThread,
    const dictionary& headerEntries
)
{
    // Can take ownership if it is the native storage type
    auto* owner = ownedStorage<writeData::dynamic_storage>(localData);

    return write_dispatch
    (
        objectType,
        fName,
        std::string_view(localData.data(), localData.size()),
        owner,
        streamOpt, atomic, append, useThread, headerEntries
    );
}


bool Foam::OFstreamCollator::write
(
    const word& objectType,
    const fileName& fName,
    std::vector<char>&& localData,
    IOstreamOption streamOpt,
    IOstreamOption::atomicType atomic,
    IOstreamOption::appendType append,
    const bool useThread,
    const dictionary& headerEntries
)
{
    // Can take ownership if it is the native storage type
    auto* owner = ownedStorage<writeData::dynamic_storage>(localData);

    return write_dispatch
    (
        objectType,
        fName,
        std::string_view(localData.data(), localData.size()),
        owner,
        streamOpt, atomic, append, useThread, headerEntries
    );
}


bool Foam::OFstreamCollator::write
(
    const word& objectType,
    const fileName& fName,
    std::string&& localData,
    IOstreamOption streamOpt,
    IOstreamOption::atomicType atomic,
    IOstreamOption::appendType append,
    const bool useThread,
    const dictionary& headerEntries
)
{
    // Can take ownership if it is the native storage type
    auto* owner = ownedStorage<writeData::dynamic_storage>(localData);

    return write_dispatch
    (
        objectType,
        fName,
        std::string_view(localData.data(), localData.size()),
        owner,
        streamOpt, atomic, append, useThread, headerEntries
    );
}


bool Foam::OFstreamCollator::write
(
    const word& objectType,
    const fileName& fName,
    std::string_view localData,
    IOstreamOption streamOpt,
    IOstreamOption::atomicType atomic,
    IOstreamOption::appendType append,
    const bool useThread,
    const dictionary& headerEntries
)
{
    // No ownership of a string_view is possible.
    // Use view directly, content is copied only if required (threaded)

    return write_dispatch
    (
        objectType,
        fName,
        localData,
        nullptr,  // owner
        streamOpt, atomic, append, useThread, headerEntries
    );
}


bool Foam::OFstreamCollator::write_dispatch
(
    const word& objectType,
    const fileName& fName,
    std::string_view localData,
    writeData::dynamic_storage* localDataOwner,
    IOstreamOption streamOpt,
    IOstreamOption::atomicType atomic,
    IOstreamOption::appendType append,
    const bool useThread,
    const dictionary& headerEntries
)
{
    // Determine (on master) sizes to receive. Note: do NOT use thread
    // communicator

    const List<int64_t> recvSizes
    (
        UPstream::listGatherValues(int64_t(localData.size()), localComm_)
    );

    // Determine the dimensions (maxLocalSize, totalSize)
    int64_t maxSizes[2] = { 0, 0 };

    for (auto recvSize : recvSizes)
    {
        maxSizes[0] = std::max(maxSizes[0], recvSize);
        maxSizes[1] += recvSize;
    }

    // Broadcast the information to everyone
    UPstream::broadcast(maxSizes, 2, localComm_);

    // All ranks now have the same maxLocalSize, totalSize information
    // (from the master) and can stay in lockstep later
    const int64_t maxLocalSize = maxSizes[0];
    const int64_t totalSize = maxSizes[1];


    // Determine how things will be gathered and written...

    enum class dispatchModes { DIRECT_WRITE, THREADED_WRITE, FULL_THREADED };

    dispatchModes dispatch(dispatchModes::DIRECT_WRITE);

    if (!useThread || maxBufferSize_ == 0 || maxLocalSize > maxBufferSize_)
    {
        dispatch = dispatchModes::DIRECT_WRITE;
    }
    else if (totalSize <= maxBufferSize_)
    {
        // Total size can be stored locally
        // - gather all data now and only do the writing in the thread

        dispatch = dispatchModes::THREADED_WRITE;
    }
    else
    {
        // Gather data and write in the thread

        dispatch = dispatchModes::FULL_THREADED;

        if (!UPstream::haveThreads())
        {
            WarningInFunction
                << "MPI not initialized with thread support." << nl
                << "    maxThreadFileBufferSize = 0 to disable threading" << nl
                << " or maxThreadFileBufferSize > " << totalSize
                << " to collate before threaded writing." << nl << nl;

            dispatch = dispatchModes::DIRECT_WRITE;
        }
    }


    // -----------
    // Dispatching
    // -----------

    if (dispatch == dispatchModes::DIRECT_WRITE)
    {
        if (debug)
        {
            Pout<< "OFstreamCollator : non-thread gather "
                << "(local comm: " << localComm_
                << "); non-thread write of "
                << fName << endl;
        }

        // Direct collating and writing (so master blocks until all written!)
        return writeFile
        (
            localComm_,
            objectType,
            fName,
            localData,
            recvSizes,
            UList<std::string_view>::null(),  // dummy proc data
            streamOpt,
            atomic,
            append,
            headerEntries
        );
    }
    else if (dispatch == dispatchModes::THREADED_WRITE)
    {
        if (debug)
        {
            Pout<< "OFstreamCollator : non-thread gather";

            if (UPstream::master(localComm_))
            {
                Pout<< "[" << (localDataOwner ? "copying" : "moving")
                    << " localData prior to threading]";
            }

            Pout<< " (local comm: " << localComm_
                << "); thread write of "
                << fName << endl;
        }

        if (UPstream::master(localComm_))
        {
            waitForBufferSpace(totalSize);
        }

        std::unique_ptr<writeData> fileAndDataPtr
        (
            new writeData
            (
                threadComm_,
                objectType,
                fName,
                recvSizes,
                streamOpt,
                atomic,
                append,
                headerEntries
            )
        );
        auto& procData = fileAndDataPtr->procData_;

        if (UPstream::master(localComm_))
        {
            // Move in (or copy) local data (master only!)
            if (localDataOwner)
            {
                fileAndDataPtr->transfer(*localDataOwner);
            }
            else
            {
                fileAndDataPtr->assign(localData);
            }

            // Storage for receive data
            procData.resize(UPstream::nProcs(localComm_));
        }
        else if (UPstream::is_subrank(localComm_))
        {
            // Requires sizing for decomposedBlockData::writeBlocks() logic
            procData.resize(UPstream::nProcs(localComm_));
        }


        // Gather all data onto master. Is done in local communicator since
        // not in write thread.
        const auto startOfRequests = UPstream::nRequests();

        const int messageTag = (UPstream::msgType() + 256);

        if (UPstream::master(localComm_))
        {
            // The receive procs, large messages first
            DynamicList<int> recvProcs(recvSizes.size());

            forAll(recvSizes, proci)
            {
                // Ignore empty slots and rank=0 [self]
                if (recvSizes[proci] > 0 && proci != UPstream::masterNo())
                {
                    recvProcs.push_back(proci);
                }
            }

            // Reverse ordered
            std::stable_sort
            (
                recvProcs.begin(), recvProcs.end(),
                [&](int a, int b){ return (recvSizes[b] < recvSizes[a]); }
            );

            // Start the receives (already filtered to remove empty slots)
            for (const int proci : recvProcs)
            {
                auto& slot = procData[proci];

                slot.resize(recvSizes[proci]);

                UIPstream::read
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    slot.data(),
                    slot.size(),
                    messageTag,
                    localComm_
                );
            }
        }
        else if (UPstream::is_subrank(localComm_) && !localData.empty())
        {
            // Send to content to master
            if
            (
               !UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    UPstream::masterNo(),
                    localData.data(),
                    localData.size(),
                    messageTag,
                    localComm_
                )
            )
            {
                FatalErrorInFunction
                    << "Failure to send message (size: "
                    << localData.size() << ") to master" << nl
                    << Foam::abort(FatalError);
            }
        }

        UPstream::waitRequests(startOfRequests);

        // The localData has been moved (master) or communicated
        if (localDataOwner)
        {
            localDataOwner->clear();
            localDataOwner->shrink_to_fit();
        }


        // Queue up for threading
        {
            std::lock_guard<std::mutex> guard(mutex_);

            // Add to thread buffer (as FIFO), take ownership
            objects_.push_back(fileAndDataPtr.release());

            // Start thread if not running
            if (!threadRunning_)
            {
                if (thread_)
                {
                    if (debug)
                    {
                        Pout<< "OFstreamCollator : Waiting for write thread"
                            << endl;
                    }
                    thread_->join();
                }

                if (debug)
                {
                    Pout<< "OFstreamCollator : Starting write thread"
                        << endl;
                }
                thread_.reset(new std::thread(writeAll, this));
                threadRunning_ = true;
            }
        }

        return true;
    }
    else if (dispatch == dispatchModes::FULL_THREADED)
    {
        if (debug)
        {
            Pout<< "OFstreamCollator : thread gather; thread write"
                << "[" << (localDataOwner ? "copying" : "moving")
                << " localData prior to threading]";

            Pout<< " (thread comm: " << threadComm_
                << ") of " << fName << endl;
        }

        if (UPstream::master(localComm_))
        {
            waitForBufferSpace(localData.size());
        }

        std::unique_ptr<writeData> fileAndDataPtr
        (
            new writeData
            (
                threadComm_,
                objectType,
                fName,
                recvSizes,
                streamOpt,
                atomic,
                append,
                headerEntries
            )
        );

        // Move in (or copy) local data (all procs)
        if (localDataOwner)
        {
            fileAndDataPtr->transfer(*localDataOwner);
        }
        else
        {
            fileAndDataPtr->assign(localData);
        }


        // Queue up for threading
        {
            std::lock_guard<std::mutex> guard(mutex_);

            // Add to thread buffer (as FIFO), take ownership
            objects_.push_back(fileAndDataPtr.release());

            // Note: no proc data provided
            // so it will trigger communication inside the thread!!!

            if (!threadRunning_)
            {
                if (thread_)
                {
                    if (debug)
                    {
                        Pout<< "OFstreamCollator : Waiting for write thread"
                            << endl;
                    }
                    thread_->join();
                }

                if (debug)
                {
                    Pout<< "OFstreamCollator : Starting write thread" << endl;
                }
                thread_.reset(new std::thread(writeAll, this));
                threadRunning_ = true;
            }
        }

        return true;
    }

    FatalErrorInFunction
        << "Unknown dispatch mode: " << int(dispatch)
        << " - programming error?" << abort(FatalError);

    return false;
}


void Foam::OFstreamCollator::waitAll()
{
    // Wait for all buffer space to be available
    // - ie, wait for all jobs to finish

    if (UPstream::master(localComm_))
    {
        if (debug)
        {
            Pout<< "OFstreamCollator : waiting for thread to have consumed all"
                << endl;
        }
        waitForBufferSpace(-1);
    }
}


// ************************************************************************* //
