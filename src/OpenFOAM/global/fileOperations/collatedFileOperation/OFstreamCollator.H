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

Class
    Foam::OFstreamCollator

Description
    Threaded file writer.

    Collects all data from all processors and writes as single
    'decomposedBlockData' file. The operation is determined by the
    buffer size (maxThreadFileBufferSize setting):
    - local size of data is larger than buffer: receive and write processor
    by processor (i.e. 'scheduled'). Does not use a thread, no file size
    limit.
    - total size of data is larger than buffer (but local is not):
    thread does all the collecting and writing of the processors. No file
    size limit.
    - total size of data is less than buffer:
    collecting is done locally; the thread only does the writing
    (since the data has already been collected)

SourceFiles
    OFstreamCollator.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_OFstreamCollator_H
#define Foam_OFstreamCollator_H

#include "IOstream.H"
#include "CircularBuffer.H"  // As FIFO
#include "DynamicList.H"
#include "dictionary.H"

#include <algorithm>
#include <mutex>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class OFstreamCollator Declaration
\*---------------------------------------------------------------------------*/

class OFstreamCollator
{
    // Private Class

        //- Holds data to be written
        struct writeData
        {
            //- Storage type for local data
            using dynamic_storage = std::vector<char>;
            //- Storage type for individual proc data
            using static_storage = std::vector<char>;

            const int comm_;
            const word objectType_;
            const fileName pathName_;

            dynamic_storage localData_;      // local data
            List<static_storage> procData_;  // subrank data (if any)

            const List<int64_t> sizes_;
            const IOstreamOption streamOpt_;
            IOstreamOption::atomicType atomic_;
            IOstreamOption::appendType append_;
            const dictionary headerEntries_;

            writeData() = delete;                      // No default construct
            writeData(const writeData&) = delete;      // No copy construct
            writeData(writeData&&) = delete;           // No move construct
            void operator=(const writeData&) = delete; // No copy assign
            void operator=(writeData&&) = delete;      // No move assign

            //- Construct without local or proc data
            writeData
            (
                const int communicator,
                const word& objectType,
                const fileName& pathName,
                const UList<int64_t>& sizes,
                IOstreamOption streamOpt,
                IOstreamOption::atomicType atomic,
                IOstreamOption::appendType append,
                const dictionary& headerEntries
            )
            :
                comm_(communicator),
                objectType_(objectType),
                pathName_(pathName),
                sizes_(sizes),
                streamOpt_(streamOpt),
                atomic_(atomic),
                append_(append),
                headerEntries_(headerEntries)
            {}

            //- Copy assign local data
            void assign(std::string_view sv)
            {
                localData_.resize(sv.size());
                std::copy(sv.begin(), sv.end(), localData_.begin());
            }

            //- Move reset local data
            void transfer(dynamic_storage& localData)
            {
                localData_ = std::move(localData);
            }

            //- A view of local data
            std::string_view localView() const
            {
                return std::string_view(localData_.data(), localData_.size());
            }

            //- A view of optional proc data
            std::string_view procView(int proci) const
            {
                if (proci >= 0 && proci < procData_.size())
                {
                    return std::string_view
                    (
                        procData_[proci].data(),
                        procData_[proci].size()
                    );
                }
                else
                {
                    return std::string_view();
                }
            }

            //- The (approximate) size of local + any optional proc data
            int64_t size() const
            {
                int64_t total = int64_t(localData_.size());
                for (const auto& data : procData_)
                {
                    total += int64_t(data.size());
                }
                return total;
            }
        };


    // Private Data

        //- Total amount of storage to use for object stack below
        const int64_t maxBufferSize_;

        mutable std::mutex mutex_;

        std::unique_ptr<std::thread> thread_;

        //- FIFO of files to write and their contents
        CircularBuffer<writeData*> objects_;

        //- Whether thread is running (and not exited)
        bool threadRunning_;

        //- Communicator to use for all parallel ops (in simulation thread)
        int localComm_;

        //- Communicator to use for all parallel ops (in write thread)
        int threadComm_;


    // Private Member Functions

        //- Write actual file
        static bool writeFile
        (
            const int communicator,
            const word& objectType,
            const fileName& fName,
            std::string_view localData,
            const UList<int64_t>& recvSizes,
            //! optional proc data
            const UList<std::string_view>& procData,
            IOstreamOption streamOpt,
            IOstreamOption::atomicType atomic,
            IOstreamOption::appendType append,
            const dictionary& headerEntries
        );

        //- Write implementation/dispatcher.
        //
        //  The content is always accessed via the string_view of \p localData.
        //
        //  The \p localDataOwner refers to the storage behind the view that
        //  can be moved (no copy) when queuing for the write thread.
        //  A nullptr indicates ownership cannot be transferred and
        //  thus a deep copy of the \p localData string_view will be used
        //  for the threaded writer.
        bool write_dispatch
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
        );

        //- Write all files in stack
        static void* writeAll(void *threadarg);

        //- Wait for total size of objects_ (master + optional slave data)
        //  to be wantedSize less than overall maxBufferSize.
        void waitForBufferSpace(const int64_t wantedSize) const;


public:

    // Declare name of the class and its debug switch
    TypeName("OFstreamCollator");


    // Constructors

        //- Construct from buffer size (0 = do not use thread)
        //- and with worldComm
        explicit OFstreamCollator(int64_t maxBufferSize);

        //- Construct from buffer size (0 = do not use thread)
        //- and specified communicator
        OFstreamCollator(int64_t maxBufferSize, int communicator);


    //- Destructor
    virtual ~OFstreamCollator();


    // Member Functions

        //- Write file with contents, possibly taking ownership of the
        //- content (if it is the native storage type).
        //
        //  Blocks until write-thread has space available
        //  (total file sizes < maxBufferSize)
        bool write
        (
            const word& objectType,
            const fileName& fName,
            DynamicList<char>&& localData,
            IOstreamOption streamOpt,
            IOstreamOption::atomicType atomic,
            IOstreamOption::appendType append,
            const bool useThread = true,
            const dictionary& headerEntries = dictionary::null
        );

        //- Write file with contents, possibly taking ownership of the
        //- content (if it is the native storage type).
        //
        //  Blocks until write-thread has space available
        //  (total file sizes < maxBufferSize)
        bool write
        (
            const word& objectType,
            const fileName& fName,
            std::vector<char>&& localData,
            IOstreamOption streamOpt,
            IOstreamOption::atomicType atomic,
            IOstreamOption::appendType append,
            const bool useThread = true,
            const dictionary& headerEntries = dictionary::null
        );

        //- Write file with contents, possibly taking ownership of the
        //- content (if it is the native storage type).
        //
        //  Blocks until write-thread has space available
        //  (total file sizes < maxBufferSize)
        bool write
        (
            const word& objectType,
            const fileName& fName,
            std::string&& localData,
            IOstreamOption streamOpt,
            IOstreamOption::atomicType atomic,
            IOstreamOption::appendType append,
            const bool useThread = true,
            const dictionary& headerEntries = dictionary::null
        );

        //- Write file with contents.
        //
        //  Blocks until write-thread has space available
        //  (total file sizes < maxBufferSize)
        FOAM_DEPRECATED_FOR(2023-09, "prefer write with movable content")
        bool write
        (
            const word& objectType,
            const fileName& fName,
            std::string_view localData,
            IOstreamOption streamOpt,
            IOstreamOption::atomicType atomic,
            IOstreamOption::appendType append,
            const bool useThread = true,
            const dictionary& headerEntries = dictionary::null
        );

        //- Wait for all thread actions to have finished
        void waitAll();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
