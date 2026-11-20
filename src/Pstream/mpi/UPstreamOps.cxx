/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "UPstream.H"
#include "UPstreamWrapping.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * Reduction * * * * * * * * * * * * * * * * //

int Foam::UPstream::find_first(bool condition, int communicator)
{
    if (UPstream::is_parallel(communicator))
    {
        // Ensure we use exact and known data types
        typedef int32_t IntType;

        auto proci = static_cast<IntType>(UPstream::myProcNo(communicator));
        auto nproc = static_cast<IntType>(UPstream::nProcs(communicator));

        if (!condition)
        {
            // The out-of-range value (min reduction)
            proci = nproc;
        }

        PstreamDetail::allReduce
        (
            &proci,
            1,
            MPI_INT32_T,
            MPI_MIN,
            communicator
        );

        // Found: the first rank (min). Not found: -1
        return (proci < nproc ? proci : -1);
    }
    else
    {
        // Not parallel (serial), so we are rank=0
        // return 0 (found) or -1 (not found)
        return (condition ? 0 : -1);
    }
}


namespace
{

// Implementation detail
template<class IntType>
IntType findOffset_impl(IntType localSize, int communicator)
{
    using namespace Foam;

    if (!UPstream::is_parallel(communicator))
    {
        return 0;
    }

    // Backport: no constexpr or UPstreamTraits etc, so hardcoded instead:
    MPI_Datatype datatype =
    (
        (sizeof(IntType) == sizeof(int32_t)) ? MPI_INT32_T : MPI_INT64_T
    );

    // Probably safe without copying localSize -> input, but it's cheap.
    IntType input(localSize);
    IntType work(localSize);

    MPI_Exscan
    (
        &input,
        &work,
        1,
        datatype,
        MPI_SUM,
        PstreamGlobals::MPICommunicators_[communicator]
    );

    // In exclusive mode, the value for rank=0 is degenerate and the
    // standard doesn't say which value (if any) it will have.
    //
    // For Exscan_sum, setting it as 0 for rank=0 is what we do.

    if (UPstream::master(communicator))
    {
        work = 0;
    }

    return work;
}

} // End anonymous namespace


int32_t Foam::UPstream::find_offset(int32_t localSize, int communicator)
{
    return findOffset_impl(localSize, communicator);
}

int64_t Foam::UPstream::find_offset(int64_t localSize, int communicator)
{
    return findOffset_impl(localSize, communicator);
}


namespace
{

// Implementation detail
template<class Type>
void minLocation_impl
(
    Foam::Tuple2<Type, int>* tuples,
    int64_t count,
    int communicator
)
{
    using namespace Foam;

    if (!count || !tuples || !UPstream::is_parallel(communicator))
    {
        return;
    }

    // only (float|double), so basic type is same as Type
    typedef Type basic_type;

    // The C++ container
    using  cppTuple = Tuple2<Type, int>;

    // The struct for MPI
    struct mpiTuple { basic_type first; int second; };

    static_assert
    (
        (sizeof(cppTuple) >= sizeof(mpiTuple)),
        "Tuple sizeof is not >= sizeof MPI struct. Transcribing will fail!"
    );


    // Backport: no constexpr or UPstreamTraits etc, so hardcoded instead:
    MPI_Datatype datatype =
    (
        (sizeof(Type) == sizeof(float))
      ? MPI_FLOAT_INT
      : MPI_DOUBLE_INT
    );

    // Backport: no UPstreamTraits etc, so hardcoded instead:
    MPI_Op optype = MPI_MINLOC;

    if (count == 1)
    {
        // Fast(er) path
        auto& data = *tuples;

        // Copy into the C-struct for exact alignment etc.
        mpiTuple blob { basic_type(data.first()), data.second() };

        if
        (
            MPI_Allreduce
            (
                MPI_IN_PLACE,  // recv is also send
                &blob,
                1,
                datatype,
                optype,
                PstreamGlobals::MPICommunicators_[communicator]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Allreduce (min_loc) failed for count:1"
                << Foam::abort(FatalError);
        }

        // Copy back from C-struct to original locations
        data.first() = Type(blob.first);
        data.second() = blob.second;
        return;
    }

    // Regular path
    auto* blobs = reinterpret_cast<mpiTuple*>(tuples);

    // Transcribe to ensure data layout is exact, without padding
    for (int64_t i = 0; i < count; ++i)
    {
        mpiTuple t { basic_type(tuples[i].first()), tuples[i].second() };
        blobs[i] = t;
    }

    // The MPI (int) message-size limit, but since we are chunking
    // anyhow can take a smaller size too. Using half of that range.
    constexpr auto maxChunk = int64_t(INT_MAX/2+1);

    // Operate on chunks, single pass when no chunking is needed
    for (int64_t start = 0; start < count; start += maxChunk)
    {
        int span_count = int(std::min(maxChunk, count - start));

        if
        (
            MPI_Allreduce
            (
                MPI_IN_PLACE,  // recv is also send
                (blobs + start),
                span_count,
                datatype,
                optype,
                PstreamGlobals::MPICommunicators_[communicator]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Allreduce (min_loc) failed for count:" << span_count
                << Foam::abort(FatalError);
        }
    }

    // Copy back from C-struct to original locations. Reverse order!
    for (int64_t i = count-1; i >= 0; --i)
    {
        cppTuple t( Type(blobs[i].first), blobs[i].second );
        tuples[i] = t;
    }
}

} // End anonymous namespace


#undef  Pstream_defineMethods
#define Pstream_defineMethods(Type)                                           \
                                                                              \
void Foam::UPstream::minLocation                                              \
(                                                                             \
    Tuple2<Type, int>* data_tuples, int64_t count, int communicator           \
)                                                                             \
{                                                                             \
    minLocation_impl(data_tuples, count, communicator);                       \
}                                                                             \
                                                                              \
void Foam::UPstream::minLocation(UList<Tuple2<Type, int>>& tuples, int comm)  \
{                                                                             \
    minLocation_impl(tuples.data(), tuples.size(), comm);                     \
}                                                                             \
                                                                              \
int Foam::UPstream::minLocation                                               \
(                                                                             \
    Tuple2<Type, int>& tuple,                                                 \
    int communicator                                                          \
)                                                                             \
{                                                                             \
    minLocation_impl(&tuple, 1, communicator);                                \
    return tuple.second();                                                    \
}

// For backport only need (float|double) and minLocation
Pstream_defineMethods(float);
Pstream_defineMethods(double);

#undef Pstream_defineMethods


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Special reductions for bool

void Foam::UPstream::reduceAnd(bool values[], label count, int communicator)
{
    PstreamDetail::allReduce
    (
        values, int(count), MPI_C_BOOL, MPI_LAND, communicator
    );
}


void Foam::UPstream::reduceOr(bool values[], label count, int communicator)
{
    PstreamDetail::allReduce
    (
        values, int(count), MPI_C_BOOL, MPI_LOR, communicator
    );
}


// ************************************************************************* //
