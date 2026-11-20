/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "UPstream.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * Reduction * * * * * * * * * * * * * * * * //

int Foam::UPstream::find_first(bool condition, int communicator)
{
    // Not parallel (serial), so we are rank=0
    // return 0 (found) or -1 (not found)
    return (condition ? 0 : -1);
}

int32_t Foam::UPstream::find_offset(int32_t localSize, int communicator)
{
    return 0;
}

int64_t Foam::UPstream::find_offset(int64_t localSize, int communicator)
{
    return 0;
}


#undef  Pstream_defineMethods
#define Pstream_defineMethods(Type)                                           \
                                                                              \
int Foam::UPstream::minLocation                                               \
(                                                                             \
    Tuple2<Type, int>& tuple,                                                 \
    int communicator                                                          \
)                                                                             \
{                                                                             \
    return tuple.second();                                                    \
}                                                                             \
                                                                              \
void Foam::UPstream::minLocation                                              \
(                                                                             \
    Tuple2<Type, int>* tuples,                                                \
    int64_t count,                                                            \
    int communicator                                                          \
)                                                                             \
{}                                                                            \
                                                                              \
void Foam::UPstream::minLocation                                              \
(                                                                             \
    UList<Tuple2<Type, int>>& tuples,                                         \
    int communicator                                                          \
)                                                                             \
{}                                                                            \

// For backport only need (float|double) and minLocation
Pstream_defineMethods(float);
Pstream_defineMethods(double);

#undef Pstream_defineMethods


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Special reductions for bool

void Foam::UPstream::reduceAnd(bool values[], label count, int communicator)
{}


void Foam::UPstream::reduceOr(bool values[], label count, int communicator)
{}


// ************************************************************************* //
