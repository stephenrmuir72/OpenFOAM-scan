/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
    Copyright (C) 2020-2025 OpenCFD Ltd.
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

#include "decomposedBlockData.H"
#include "Fstream.H"
#include "IPstream.H"
#include "OPstream.H"
#include "SpanStream.H"
#include "dictionary.H"
#include "objectRegistry.H"
#include "masterUncollatedFileOperation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decomposedBlockData, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::decomposedBlockData::isCollatedType(const word& objectType)
{
    return (objectType == decomposedBlockData::typeName);
}


bool Foam::decomposedBlockData::isCollatedType(const IOobject& io)
{
    // same: return io.isHeaderClass(decomposedBlockData::typeName);
    // same: return isCollatedType(io.headerClassName());
    return io.isHeaderClass<decomposedBlockData>();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decomposedBlockData::decomposedBlockData
(
    const int communicator,
    const IOobject& io,
    const UPstream::commsTypes commsType
)
:
    regIOobject(io),
    commsType_(commsType),
    comm_(communicator)
{
    // Temporary warning
    if (readOpt() == IOobjectOption::READ_MODIFIED)
    {
        WarningInFunction
            << "decomposedBlockData " << name()
            << " constructed with READ_MODIFIED"
            " but decomposedBlockData does not support automatic rereading."
            << endl;
    }
    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        read();
    }
}


// * * * * * * * * * * * * * * * Members Functions * * * * * * * * * * * * * //

std::streamoff Foam::decomposedBlockData::writeBlockEntry
(
    OSstream& os,
    const label blocki,
    const char* str,
    const size_t len
)
{
    // Offset to the beginning of this output
    // This should generally be OK for non-compressed streams
    // (eg, std::ofstream)

    std::streamoff blockOffset = os.stdStream().tellp();

    const word procName("processor" + Foam::name(blocki));

    // Write as primitiveEntry or commented content
    constexpr bool isDictFormat = false;

    if constexpr (isDictFormat)
    {
        // Like writeKeyword()
        os << nl << procName << nl;
    }
    else
    {
        // Human-readable comments
        os << nl << "// " << procName << nl;
    }

    // Data
    if (str && len > 0)
    {
        // Special treatment for char data (binary I/O only)
        const auto oldFmt = os.format(IOstreamOption::BINARY);

        os << int64_t(len) << nl;
        os.write(str, len);
        os << nl;

        os.format(oldFmt);
    }
    else
    {
        os << int64_t(0) << nl;
    }

    if constexpr (isDictFormat)
    {
        os.endEntry();
    }

    return blockOffset;
}


namespace
{

// Handler for readBlockEntry with
// - std::string,
// - std::vector<char>
// - List<char>
//
// The List<char> version currently could handle any of these:
// - 0.  NCHARS (...)
// - 1.  List<char> NCHARS (...)
// - 2.  processorN  List<char> NCHARS (...) ;
// - 3.  processorN  NCHARS (...) ;
//
// The std::string, std::vector<char> versions can only handle these:
// - 0.  NCHARS (...)
// - 3.  processorN  NCHARS (...) ;
//
// May wish to drop the idea of dictionary format anyhow,
// since it was never better than a simple proof-of-concept.

template<class BufferType>
bool readBlockEntry_impl(Foam::Istream& is, BufferType& charData)
{
    using namespace Foam;

    // Report the public API name instead of the template's mangled signature
    constexpr const char* funcName = "decomposedBlockData::readBlockEntry";

    // Is the buffer a List<char> type?
    constexpr bool is_list_buffer = std::is_base_of_v<List<char>, BufferType>;

    // static_assert
    // (
    //     (
    //         std::is_base_of_v<std::string, BufferType>
    //      || std::is_base_of_v<std::vector<char>, BufferType>
    //      || is_list_buffer
    //     ),
    //     "Helper is for string, vector<char>, List<char>"
    // );

    // As per List<char> version, but reading into
    // std::string or std::vector.

    is.fatalCheck(funcName);
    token tok(is);
    is.fatalCheck(funcName);

    // Dictionary format?
    // - has a keyword, which is not a compound token:

    const bool isDictFormat = (tok.isWord() && !tok.isCompound());

    if constexpr (is_list_buffer)
    {
        // Token is compound or label, so put it back for readList
        if (!isDictFormat && tok.good())
        {
            is.putBack(tok);
        }
    }
    else
    {
        // std::string or std::vector, so no compound token possible...
        if (isDictFormat)
        {
            is.read(tok);  // Skip the keyword, get another token
        }
        if (tok.isCompound())
        {
            is.read(tok);  // Skip the compound, get another token
        }
    }

    if constexpr (is_list_buffer)
    {
        charData.readList(is);
    }
    else
    {
        // Any future possibility of using dictionary format would fall apart
        // here since we wouldn't be able to steal the compoundToken
        // contents into std::string/std::vector and has label restrictions
        // etc.

        // Handle like List<char>::readList() - binary

        charData.clear();

        if (tok.is_int64())  // int64, int32, etc
        {
            const int64_t len = tok.int64Token();

            if (len > 0)
            {
                charData.resize(std::size_t(len));

                // Special treatment for char data (binary I/O only)
                const auto oldFmt = is.format(IOstreamOption::BINARY);

                // read(...) includes surrounding start/end delimiters
                is.read(charData.data(), charData.size());

                is.fatalCheck
                (
                    "decomposedBlockData::readBlockEntry : [binary block]"
                );

                is.format(oldFmt);
            }
            // - len == 0 : there is no further input, not even brackets.
            // - len < 0  : corrupt, let someone else deal with it
        }
        else
        {
            // FatalIOErrorInFunction(is)
            FatalIOErrorIn(funcName, is)
                << "Expected (int32/int64) for the block size, found "
                << tok.info()
                << Foam::exit(FatalIOError);
        }
    }

    if (isDictFormat)
    {
        is.fatalCheck(funcName);
        is >> tok;
        is.fatalCheck(funcName);

        // Swallow trailing ';'
        if (tok.good() && !tok.isPunctuation(token::END_STATEMENT))
        {
            is.putBack(tok);
        }
    }

    return true;
}

} // End anonymous namespace


bool Foam::decomposedBlockData::readBlockEntry
(
    Istream& is,
    std::string& charData
)
{
    return readBlockEntry_impl(is, charData);
}


bool Foam::decomposedBlockData::readBlockEntry
(
    Istream& is,
    std::vector<char>& charData
)
{
    return readBlockEntry_impl(is, charData);
}


bool Foam::decomposedBlockData::readBlockEntry
(
    Istream& is,
    List<char>& charData
)
{
    return readBlockEntry_impl(is, charData);
}


bool Foam::decomposedBlockData::skipBlockEntry(Istream& is)
{
    // As per readBlockEntry but seeks instead of reading.
    // Internals like charList::readList - ie, always binary

    // Handle any of these:
    // 0.  NCHARS (...)
    // 1.  List<char> NCHARS (...)
    // 2.  processorN  List<char> NCHARS (...) ;
    // 3.  processorN  NCHARS (...) ;

    if (!is.good()) return false;
    token tok(is);
    if (!is.good()) return false;

    // Dictionary format?
    // - has a keyword, which is not a compound token:
    const bool isDictFormat = (tok.isWord() && !tok.isCompound());

    if (isDictFormat)
    {
        is >> tok;
        if (!is.good()) return false;
    }


    bool handled = false;

    // Like charList::readList
    if (tok.isCompound())
    {
        handled = true;
    }
    else if (tok.is_int64())  // int64, int32, etc
    {
        // Label: could be int(..) or just a plain '0'

        const int64_t len = tok.int64Token();

        // Special treatment for char data (binary I/O only)
        const auto oldFmt = is.format(IOstreamOption::BINARY);

        if (len)
        {
            // read(...) includes surrounding start/end delimiters.

            // Note: nullptr to ignore instead of reading
            is.read(nullptr, std::streamsize(len));
        }
        is.format(oldFmt);

        handled = true;
    }
    else
    {
        // Incorrect token
        return false;
    }

    if (isDictFormat)
    {
        is.fatalCheck(FUNCTION_NAME);
        is >> tok;
        is.fatalCheck(FUNCTION_NAME);

        // Swallow trailing ';'
        if (tok.good() && !tok.isPunctuation(token::END_STATEMENT))
        {
            is.putBack(tok);
        }
    }

    return handled;
}


Foam::label Foam::decomposedBlockData::getNumBlocks
(
    Istream& is,
    const label maxNumBlocks
)
{
    label nBlocks = 0;

    // Handle OpenFOAM header if it is the first entry
    if (is.good())
    {
        token tok(is);

        if (is.good() && tok.isWord("FoamFile"))
        {
            dictionary headerDict(is);  // Read sub-dictionary content

            if (headerDict.readIfPresent("version", tok))
            {
                is.version(tok);
            }

            word formatName;
            if (headerDict.readIfPresent("format", formatName))
            {
                is.format(formatName);
            }

            //// Obtain number of blocks directly
            ///  This may not be reliable...
            //if (headerDict.readIfPresent("blocks", nBlocks))
            //{
            //    return nBlocks;
            //}
        }
        else if (tok.good())
        {
            is.putBack(tok);
        }
    }

    while (is.good() && skipBlockEntry(is))
    {
        ++nBlocks;

        if (maxNumBlocks == nBlocks)
        {
            break;
        }
    }

    return nBlocks;
}


bool Foam::decomposedBlockData::hasBlock(Istream& is, const label blockNumber)
{
    return
    (
        blockNumber >= 0
     && (blockNumber < getNumBlocks(is, blockNumber+1))
    );
}


std::streamoff Foam::decomposedBlockData::writeBlockEntry
(
    OSstream& os,
    IOstreamOption streamOptData,
    const regIOobject& io,
    const label blocki,
    const bool withLocalHeader
)
{
    // Serialize content to write
    OCharStream serialized(streamOptData);
    {
        auto& buf = serialized;
        buf.reserve(4*1024);  // Start with a slightly larger buffer

        bool ok = true;

        // Generate FoamFile header on master, without comment banner
        if (withLocalHeader)
        {
            const bool old = IOobject::bannerEnabled(false);

            ok = io.writeHeader(buf);

            IOobject::bannerEnabled(old);
        }

        // Serialize the output
        ok = ok && io.writeData(buf);

        if (!ok)
        {
            return std::streamoff(-1);
        }
    }

    return decomposedBlockData::writeBlockEntry
    (
        os,
        blocki,
        serialized.view()
    );
}


Foam::autoPtr<Foam::ISstream>
Foam::decomposedBlockData::readBlock
(
    const label blocki,
    ISstream& is,
    IOobject& headerIO
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readBlock:"
            << " stream:" << is.name() << " attempt to read block " << blocki
            << endl;
    }

    // The character input stream for the specified block
    autoPtr<ISstream> blockIsPtr;

    // Extracted header information
    IOstreamOption streamOptData;
    unsigned labelWidth = is.labelByteSize();
    unsigned scalarWidth = is.scalarByteSize();

    // Read master for header
    List<char> blockChars;
    decomposedBlockData::readBlockEntry(is, blockChars);

    if (blocki == 0)  // ie, UPstream::masterNo()
    {
        blockIsPtr.reset(new ICharStream(std::move(blockChars)));
        blockIsPtr->name() = is.name();

        {
            // Read header from first block,
            // advancing the stream position
            if (!headerIO.readHeader(*blockIsPtr))
            {
                FatalIOErrorInFunction(*blockIsPtr)
                    << "Problem while reading object header "
                    << is.relativeName() << nl
                    << exit(FatalIOError);
            }
        }
    }
    else
    {
        {
            // Read header from first block,
            // without advancing the stream position
            ISpanStream headerStream(blockChars.data(), blockChars.size());
            if (!headerIO.readHeader(headerStream))
            {
                FatalIOErrorInFunction(headerStream)
                    << "Problem while reading object header "
                    << is.relativeName() << nl
                    << exit(FatalIOError);
            }
            streamOptData = static_cast<IOstreamOption>(headerStream);
            labelWidth = headerStream.labelByteSize();
            scalarWidth = headerStream.scalarByteSize();
        }

        // Skip intermediate blocks
        for (label i = 1; i < blocki; ++i)
        {
            decomposedBlockData::skipBlockEntry(is);
        }

        // Read the block of interest
        decomposedBlockData::readBlockEntry(is, blockChars);

        blockIsPtr.reset(new ICharStream(std::move(blockChars)));
        blockIsPtr->name() = is.name();

        // Apply stream settings
        {
            auto& iss = blockIsPtr();
            iss.format(streamOptData.format());
            iss.version(streamOptData.version());
            iss.setLabelByteSize(labelWidth);
            iss.setScalarByteSize(scalarWidth);
        }
    }

    return blockIsPtr;
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Read data (on master) and scatter to individual ranks
//
// BufferType could be
// - std::string
// - std::vector<char>
// - List<char>
// - DynamicList<char>

template<class BufferType>
static bool readBlocks_scatter
(
    const int communicator,
    // [in] The input stream (only valid on master)
    autoPtr<ISstream>& isPtr,
    // [out] The processor local data
    BufferType& localData
)
{
    // Report the public API name instead of the template's mangled signature
    constexpr const char* funcName = "decomposedBlockData::readBlocks_scatter";

    if (decomposedBlockData::debug)
    {
        Pout<< "decomposedBlockData::readBlocks_scatter:"
            << " stream:" << (isPtr ? isPtr->name() : "<null>")
            << " non-blocking comm:" << communicator << endl;
    }

    // Read data on master and transmit. Always non-blocking
    bool ok = false;

    // The send buffers
    List<BufferType> procBuffers;

    // Some unique tag for this read/write grouping (as extra safety)
    const int messageTag = (UPstream::msgType() + 256);

    const auto startOfRequests = UPstream::nRequests();

    if (UPstream::master(communicator))
    {
        auto& is = isPtr();
        is.fatalCheck(FUNCTION_NAME);

        // Read master data
        decomposedBlockData::readBlockEntry(is, localData);

        // Read proc data and setup non-blocking sends
        procBuffers.resize(UPstream::nProcs(communicator));
        for (const int proci : UPstream::subProcs(communicator))
        {
            auto& slot = procBuffers[proci];

            decomposedBlockData::readBlockEntry(is, slot);

            // Send content (non-blocking)
            UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                slot.data(),
                slot.size(),
                messageTag,
                communicator
            );
        }

        ok = is.good();
    }
    else if (UPstream::is_subrank(communicator))
    {
        auto& slot = localData;

        // Probe for the message size
        auto [fromProci, num_bytes] =
            UPstream::probeMessage
            (
                UPstream::masterNo(),
                messageTag,
                communicator
            );

        if (num_bytes < 0)
        {
            FatalErrorIn(funcName)
                << "probed to receive negative bytes from rank:" << fromProci
                << Foam::exit(FatalError);
        }

        if (decomposedBlockData::debug)
        {
            Pout<< "probed to receive " << num_bytes
                << " from " << fromProci << endl;
        }

        slot.resize(num_bytes);

        // Receive content (can also be zero-sized)
        UIPstream::read
        (
            UPstream::commsTypes::nonBlocking,
            UPstream::masterNo(),
            slot.data(),
            slot.size(),
            messageTag,
            communicator
        );
    }

    UPstream::waitRequests(startOfRequests);
    procBuffers.clear();

    // Sync the status
    Pstream::broadcast(ok, communicator);

    return ok;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::ISstream>
Foam::decomposedBlockData::readBlocks
(
    const int communicator,
    const fileName& fName,
    autoPtr<ISstream>& isPtr,
    IOobject& headerIO,
    const UPstream::commsTypes  /* unused */
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readBlocks:"
            << " stream:" << (isPtr ? isPtr->name() : "<null>")
            << " non-blocking" << endl;
    }

    // Read data on master and transmit. Always non-blocking
    bool ok = false;
    List<char> localData;
    List<List<char>> procBuffers;
    autoPtr<ISstream> blockIsPtr;

    // Some unique tag for this read/write/probe grouping
    const int messageTag = (UPstream::msgType() + 256);

    const auto startOfRequests = UPstream::nRequests();

    if (UPstream::master(communicator))
    {
        auto& is = *isPtr;
        is.fatalCheck(FUNCTION_NAME);

        // Read master data
        decomposedBlockData::readBlockEntry(is, localData);

        // Move block data into a stream
        blockIsPtr.reset(new ICharStream(std::move(localData)));
        blockIsPtr->name() = fName;

        {
            // Read header from first block,
            // advancing the stream position
            if (!headerIO.readHeader(*blockIsPtr))
            {
                FatalIOErrorInFunction(*blockIsPtr)
                    << "Problem while reading object header "
                    << is.relativeName() << nl
                    << exit(FatalIOError);
            }
        }

        // Read proc data and setup non-blocking sends
        procBuffers.resize(UPstream::nProcs(communicator));
        for (const int proci : UPstream::subProcs(communicator))
        {
            auto& slot = procBuffers[proci];

            decomposedBlockData::readBlockEntry(is, slot);

            // Send content - non-blocking mode
            UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                slot.data(),
                slot.size(),
                messageTag,
                communicator
            );
        }

        ok = is.good();
    }
    else if (UPstream::is_subrank(communicator))
    {
        auto& slot = localData;

        // Probe for the message size
        auto [fromProci, num_bytes] =
            UPstream::probeMessage
            (
                UPstream::masterNo(),
                messageTag,
                communicator
            );

        if (num_bytes < 0)
        {
            FatalErrorInFunction
                << "probed to receive negative bytes from rank:" << fromProci
                << Foam::exit(FatalError);
        }

        if (debug)
        {
            Pout<< "probed to receive " << num_bytes
                << " from " << fromProci << endl;
        }

        slot.resize(num_bytes);

        // Receive content (can also be zero-sized)
        UIPstream::read
        (
            UPstream::commsTypes::nonBlocking,
            UPstream::masterNo(),
            slot.data(),
            slot.size(),
            messageTag,
            communicator
        );
    }

    UPstream::waitRequests(startOfRequests);
    procBuffers.clear();

    if (UPstream::is_subrank(communicator))
    {
        // Move block data into a stream
        blockIsPtr.reset(new ICharStream(std::move(localData)));
        blockIsPtr->name() = fName;
    }

    // Broadcast master header info,
    // set stream properties from blockIsPtr on master

    int verValue(0);
    int fmtValue(0);
    unsigned labelWidth(0);
    unsigned scalarWidth(0);
    word headerName(headerIO.name());

    // The stream characteristics
    //
    // unsigned formatSizes
    // (
    //     ((static_cast<unsigned>(iss.format()) & 0xFF) << 16)
    //   | ((iss.labelByteSize() & 0xFF) << 8)
    //   | ((iss.scalarByteSize() & 0xFF))
    // );

    if (UPstream::master(communicator))
    {
        auto& iss = blockIsPtr();
        verValue = iss.version().canonical();
        fmtValue = static_cast<int>(iss.format());
        labelWidth = iss.labelByteSize();
        scalarWidth = iss.scalarByteSize();
    }

    Pstream::broadcasts
    (
        communicator,
        verValue,
        fmtValue,
        labelWidth,
        scalarWidth,
        headerName,
        headerIO.headerClassName(),
        headerIO.note()
        // Unneeded: headerIO.instance()
        // Unneeded: headerIO.local()
    );

    if (blockIsPtr)
    {
        auto& iss = *blockIsPtr;
        iss.version(IOstreamOption::versionNumber::canonical(verValue));
        iss.format(IOstreamOption::streamFormat(fmtValue));
        iss.setLabelByteSize(labelWidth);
        iss.setScalarByteSize(scalarWidth);
    }

    headerIO.rename(headerName);

    if (debug)
    {
        Info<< "reading ok:" << ok << endl;
    }

    return blockIsPtr;
}


void Foam::decomposedBlockData::gatherProcData
(
    const int communicator,
    std::string_view localData,
    const UList<int64_t>& recvSizes,  // required on master

    const labelRange& whichProcs,     // required on all procs

    CSRAddressing<int64_t>& recvAddr,
    DynamicList<char>& recvData,
    const UPstream::commsTypes commsType
)
{
    const auto myProci = UPstream::myProcNo(communicator);
    const auto numProc = UPstream::nProcs(communicator);

    // Some unique tag for this read/write grouping (as extra safety)
    const int messageTag = (UPstream::msgType() + 256);

    recvData.clear();

    // On master, calculate sizing/offsets and resize the recv buffer.
    // Do not need sliceSizes when nonBlocking
    if (UPstream::master(communicator))
    {
        auto& offsets = recvAddr.offsets();

        // Start with zero
        offsets.resize_fill(numProc+1, 0);

        int64_t totalSize = 0;

        // The offsets for the selected range
        for (auto proci : whichProcs)
        {
            offsets[proci] = totalSize;
            totalSize += recvSizes[proci];
        }

        // Fill the remaining offsets as zero-sized slots
        for (label i = whichProcs.end_value(); i < offsets.size(); ++i)
        {
            offsets[i] = totalSize;
        }

        recvData.resize(totalSize);
    }
    else
    {
        // Not meaningful for sub-ranks
        recvAddr.clear();
    }

    const auto startOfRequests = UPstream::nRequests();

    {
        if (UPstream::master(communicator))
        {
            for (auto proci : whichProcs)
            {
                const auto localStart = recvAddr.localStart(proci);
                const auto localSize = recvAddr.localSize(proci);

                if (localSize == 0)
                {
                    continue;
                }
                else if (proci == UPstream::masterNo())
                {
                    // No self-communication, although masterNo is normally
                    // not contained in whichProcs range anyhow.
                    std::copy_n
                    (
                        localData.data(), localSize,
                        (recvData.data() + localStart)
                    );
                }
                else
                {
                    // Receive non-zero content
                    UIPstream::read
                    (
                        UPstream::commsTypes::nonBlocking,
                        proci,
                        (recvData.data() + localStart),
                        localSize,
                        messageTag,
                        communicator
                    );
                }
            }
        }
        else if (whichProcs.contains(myProci) && !localData.empty())
        {
            // Send non-zero content
            UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                UPstream::masterNo(),
                localData.data(),
                localData.size(),
                messageTag,
                communicator
            );
        }
    }


    // If not non-blocking, the caller expected us wait until completion
    if (UPstream::commsTypes::nonBlocking != commsType)
    {
        UPstream::waitRequests(startOfRequests);
    }
}


bool Foam::decomposedBlockData::writeBlocks
(
    const int communicator,
    autoPtr<OSstream>& osPtr,
    List<std::streamoff>& blockOffset,
    std::string_view localData,

    const UList<int64_t>& recvSizes,  // required on master
    const UList<std::string_view>& procData,

    const UPstream::commsTypes commsType,
    const bool syncReturnState
)
{
    const label nProcs = UPstream::nProcs(communicator);

    bool ok = true;

    // Recovery of blockOffset is optional
    if (UPstream::master(communicator) && notNull(blockOffset))
    {
        blockOffset.resize(nProcs);
    }

    // Max proc data size to be received
    const int64_t maxNonLocalSize =
    (
        (
            UPstream::master(communicator)
         && (procData.empty() && recvSizes.size() > 1)
        )
      ? *(std::max_element(recvSizes.begin()+1, recvSizes.end()))
      : int64_t(0)
    );


    if (debug)
    {
        Pout<< " stream:" << (osPtr ? osPtr->name() : "<null>")
            << " data:" << localData.size()
            << " proc-data:" << procData.size()
            << " max-size:" << maxNonLocalSize
            << " " << UPstream::commsTypeNames[commsType] << endl;
    }

    if (procData.size())
    {
        // --------
        // With pre-gathered proc data
        // --------

        if (UPstream::master(communicator))
        {
            auto& os = osPtr();

            std::streamoff currOffset =
                decomposedBlockData::writeBlockEntry
                (
                    os,
                    UPstream::masterNo(),
                    localData
                );

            if (blockOffset.size() > UPstream::masterNo())
            {
                blockOffset[UPstream::masterNo()] = currOffset;
            }

            // Write all pre-gathered proc data.
            for (label proci = 1; proci < nProcs; ++proci)
            {
                currOffset =
                    decomposedBlockData::writeBlockEntry
                    (
                        os,
                        proci,
                        procData[proci]
                    );

                if (blockOffset.size() > proci)
                {
                    blockOffset[proci] = currOffset;
                }
            }

            ok = os.good();
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        // --------
        // Gather/write each rank, one at a time.
        // Note: This is often associated with maxMasterFileBufferSize == 0
        // --------

        // Some unique tag for this read/write grouping (as extra safety)
        const int messageTag = (UPstream::msgType() + 256);

        if (UPstream::master(communicator))
        {
            auto& os = osPtr();

            std::streamoff currOffset =
                decomposedBlockData::writeBlockEntry
                (
                    os,
                    UPstream::masterNo(),
                    localData
                );

            if (blockOffset.size() > UPstream::masterNo())
            {
                blockOffset[UPstream::masterNo()] = currOffset;
            }


            DynamicList<char> recvData(maxNonLocalSize);
            for (label proci = 1; proci < nProcs; ++proci)
            {
                recvData.clear();
                recvData.resize(recvSizes[proci]);

                if (!recvData.empty())
                {
                    UIPstream::read
                    (
                        UPstream::commsTypes::scheduled,
                        proci,
                        recvData.data(),
                        recvData.size(),
                        messageTag,
                        communicator
                    );
                }

                currOffset =
                    decomposedBlockData::writeBlockEntry
                    (
                        os,
                        proci,
                        recvData
                    );

                if (blockOffset.size() > proci)
                {
                    blockOffset[proci] = currOffset;
                }
            }

            ok = os.good();
        }
        else if (UPstream::is_subrank(communicator) && !localData.empty())
        {
            UOPstream::write
            (
                UPstream::commsTypes::scheduled,
                UPstream::masterNo(),
                localData.data(),
                localData.size(),
                messageTag,
                communicator
            );
        }
    }
    else
    {
        // --------
        // Gather/write ranks, packing together several smaller gathers
        // into a single buffer space
        // --------

        DynamicList<char> recvData;
        CSRAddressing<int64_t> recvAddr;  // CSR addressing into recvData

        // Offsets of combined ranks for communication.
        // Never includes master rank (handled separately)
        labelList procOffsets(nProcs, Foam::zero{});

        // Max combined data to be received (master only)
        label maxRecvCount = 0;

        if (UPstream::master(communicator))
        {
            // Find out how many ranks can be received into
            // maxMasterFileBufferSize and the corresponding schedule

            int64_t maxBufferSize
            (
                fileOperations::masterUncollatedFileOperation::
                maxMasterFileBufferSize
            );

            // Buffer must fit the largest off-processor size
            if (maxBufferSize < maxNonLocalSize)
            {
                maxBufferSize = maxNonLocalSize;
            }

            // Max combined proc data size to be received
            int64_t maxCollected = 0;

            for (label proci = 1, nChunks = 0; proci < nProcs; /*nil*/)
            {
                procOffsets[nChunks] = proci;

                // At least one proc, regardless of maxBufferSize.
                // Also handles the corner case when the first proc has
                // size 0, but the next one is too large.

                for
                (
                    int64_t total = 0;
                    (
                        proci < nProcs
                     && (!total || (total + recvSizes[proci] < maxBufferSize))
                    );
                    ++proci
                )
                {
                    total += recvSizes[proci];

                    if (maxCollected < total)
                    {
                        maxCollected = total;
                    }
                }

                procOffsets[++nChunks] = proci;
            }

            maxRecvCount = maxCollected;
        }

        if (debug && UPstream::master(communicator))
        {
            OStringStream ranges;

            for (label nChunks = 1; nChunks < nProcs; ++nChunks)
            {
                const labelRange whichProcs
                (
                    procOffsets[nChunks-1],
                    procOffsets[nChunks]-procOffsets[nChunks-1]
                );

                if (whichProcs.start() >= nProcs || whichProcs.size() <= 0)
                {
                    break;
                }

                ranges << ' ' << whichProcs.min() << '-' << whichProcs.max();
            }

            Pout<< " write-schedule:" << ranges.str().c_str() << endl;
        }


        // Same schedule to be known by everyone
        UPstream::broadcast
        (
            procOffsets.data(), procOffsets.size(),
            communicator
        );

        // (master only)
        recvData.clear();
        recvData.resize(maxRecvCount);

        if (UPstream::master(communicator))
        {
            auto& os = osPtr();

            std::streamoff currOffset =
                decomposedBlockData::writeBlockEntry
                (
                    os,
                    UPstream::masterNo(),
                    localData
                );

            if (blockOffset.size() > UPstream::masterNo())
            {
                blockOffset[UPstream::masterNo()] = currOffset;
            }
        }

        for (label nChunks = 1; nChunks < nProcs; ++nChunks)
        {
            const labelRange whichProcs
            (
                procOffsets[nChunks-1],
                procOffsets[nChunks]-procOffsets[nChunks-1]
            );

            if (whichProcs.start() >= nProcs || whichProcs.size() <= 0)
            {
                break;
            }

            const auto startOfRequests = UPstream::nRequests();

            // Setup non-blocking send/recv, optionally waiting
            // - uses (UPstream::msgType()+256)
            gatherProcData
            (
                communicator,
                localData,
                recvSizes,

                whichProcs,

                recvAddr,
                recvData,
                commsType  // ie, blocking or non-blocking
            );

            // Wait for everyone, or sequentially?
            if (UPstream::commsTypes::nonBlocking != commsType)
            {
                UPstream::waitRequests(startOfRequests);
            }

            // For sanity checks
            // const auto endOfRequests = UPstream::nRequests();

            if (UPstream::master(communicator))
            {
                auto& os = osPtr();

                // Write received data
                auto currRequest = startOfRequests;
                for (const auto proci : whichProcs)
                {
                    const auto localStart = recvAddr.localStart(proci);
                    const auto localSize = recvAddr.localSize(proci);

                    std::string_view procSlice
                    (
                        (recvData.data() + localStart),
                        localSize
                    );

                    if
                    (
                        (UPstream::commsTypes::nonBlocking == commsType)
                     && (proci != UPstream::masterNo())
                     && !procSlice.empty()
                    )
                    {
                        UPstream::waitRequest(currRequest);
                        ++currRequest;
                    }

                    std::streamoff currOffset =
                        decomposedBlockData::writeBlockEntry
                        (
                            os,
                            proci,
                            procSlice
                        );

                    if (blockOffset.size() > proci)
                    {
                        blockOffset[proci] = currOffset;
                    }
                }

                ok = os.good();
            }

            UPstream::waitRequests(startOfRequests);
        }
    }

    if (syncReturnState)
    {
        //- Enable to get synchronised error checking.
        //  Ensures that all procs are as slow as the master
        //  (which does all the writing)
        Pstream::broadcast(ok, communicator);
    }

    return ok;
}


bool Foam::decomposedBlockData::read()
{
    autoPtr<ISstream> isPtr;
    fileName objPath(fileHandler().filePath(false, *this, word::null));
    if (UPstream::master(comm_))
    {
        isPtr.reset(new IFstream(objPath));
        IOobject::readHeader(*isPtr);
    }

    return readBlocks_scatter(comm_, isPtr, contentData_);
}


bool Foam::decomposedBlockData::writeData(Ostream& os) const
{
    IOobject io(*this);
    IOstreamOption streamOpt(os);

    int verValue(0);
    int fmtValue(0);
    // Unneeded: word masterName(name());
    fileName masterLocation(instance()/db().dbDir()/local());

    // Re-read my own data to find out the header information
    if (UPstream::master(comm_))
    {
        ISpanStream headerStream(contentView());
        io.readHeader(headerStream);

        verValue = headerStream.version().canonical();
        fmtValue = static_cast<int>(headerStream.format());
    }

    // Broadcast header information
    Pstream::broadcasts
    (
        comm_,
        verValue,
        fmtValue,
        // Unneeded: masterName
        io.headerClassName(),
        io.note(),
        // Unneeded: io.instance()
        // Unneeded: io.local()
        masterLocation
    );

    streamOpt.version(IOstreamOption::versionNumber::canonical(verValue));
    streamOpt.format(IOstreamOption::streamFormat(fmtValue));

    if (UPstream::is_subrank(comm_))
    {
        decomposedBlockData::writeHeader
        (
            os,
            streamOpt,  // streamOpt for data
            io.headerClassName(),
            io.note(),
            masterLocation,
            name(),
            dictionary()
        );
    }

    const auto sv = contentView();

    // Write the character data
    if (isA<OSstream>(os))
    {
        // Serial stream - can output characters directly
        os.writeRaw(sv.data(), sv.size());
    }
    else
    {
        // Other cases are less fortunate, and no std::string_view
        os.writeQuoted(sv.data(), sv.size(), false);
    }

    if (UPstream::is_subrank(comm_))
    {
        IOobject::writeEndDivider(os);
    }

    return os.good();
}


bool Foam::decomposedBlockData::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    autoPtr<OSstream> osPtr;
    if (UPstream::master(comm_))
    {
        // Note: always write binary. These are strings so readable anyway.
        //       They have already be tokenised on the sending side.

        osPtr.reset(new OFstream(objectPath(), IOstreamOption::BINARY));

        // Update meta-data for current state
        const_cast<regIOobject&>
        (
            static_cast<const regIOobject&>(*this)
        ).updateMetaData();

        decomposedBlockData::writeHeader
        (
            *osPtr,
            streamOpt,  // streamOpt for data
            static_cast<const IOobject&>(*this)
        );
    }

    const List<int64_t> recvSizes
    (
        UPstream::listGatherValues<int64_t>(contentData_.size(), comm_)
    );

    List<std::streamoff> blockOffsets;  // Optional
    return writeBlocks
    (
        comm_,
        osPtr,
        blockOffsets,
        contentView(),
        recvSizes,
        UList<std::string_view>::null(),  // dummy proc data
        commsType_
    );
}


// ************************************************************************* //
