/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "IOmapPolyMesh.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(IOmapPolyMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::IOmapPolyMesh::readIOcontents()
{
    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        readStream(typeName) >> *this;
        close();
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOmapPolyMesh::IOmapPolyMesh(const IOobject& io, const polyMesh& mesh)
:
    regIOobject(io),
    mapPolyMesh(mesh)
{
    // Warn for MUST_READ_IF_MODIFIED
    warnNoRereading<IOmapPolyMesh>();

    readIOcontents();
}


Foam::IOmapPolyMesh::IOmapPolyMesh
(
    const IOobject& io,
    const mapPolyMesh& map
)
:
    regIOobject(io),
    mapPolyMesh(map.mesh())
{
    // Warn for MUST_READ_IF_MODIFIED
    warnNoRereading<IOmapPolyMesh>();

    if (!readIOcontents())
    {
        mapPolyMesh::operator=(map);
    }
}


//Foam::IOmapPolyMesh::IOmapPolyMesh
//(
//    const IOobject& io,
//    mapPolyMesh&& map
//)
//:
//    regIOobject(io),
//    mapPolyMesh(map.mesh())
//{
//    // Warn for MUST_READ_IF_MODIFIED
//    warnNoRereading<IOmapPolyMesh>();
//
//    mapPolyMesh::transfer(map);
//
//    readIOcontents();
//}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::IOmapPolyMesh::readData(Istream& is)
{
    is >> *this;
    return is.good();
}


bool Foam::IOmapPolyMesh::writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}


// ************************************************************************* //
