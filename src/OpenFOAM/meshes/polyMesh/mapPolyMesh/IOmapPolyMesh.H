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

Class
    Foam::IOmapPolyMesh

Description
    IOmapPolyMesh is derived from mapPolyMesh and
    IOobject to give the mapPolyMesh
    automatic IO functionality via the objectRegistry.

SourceFiles
    IOmapPolyMesh.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_IOmapPolyMesh_H
#define Foam_IOmapPolyMesh_H

#include "mapPolyMesh.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class IOmapPolyMesh Declaration
\*---------------------------------------------------------------------------*/

class IOmapPolyMesh
:
    public regIOobject,
    public mapPolyMesh
{
    // Private Member Functions

        //- Read if IOobject flags set. Return true if read.
        bool readIOcontents();

public:

    //- Runtime type information
    TypeName("mapPolyMesh");


    // Constructors

        //- Construct from IOobject
        explicit IOmapPolyMesh(const IOobject& io, const polyMesh& mesh);

        //- Construct from IOobject, copying mapPolyMesh contents
        IOmapPolyMesh
        (
            const IOobject& io,
            const mapPolyMesh& map
        );

        ////- Construct from IOobject, moving mapPolyMesh contents
        //IOmapPolyMesh
        //(
        //    const IOobject& io,
        //    mapPolyMesh&& map
        //);


    //- Destructor
    virtual ~IOmapPolyMesh() = default;


    // Member Functions

        //- The readData method for regIOobject read operation
        virtual bool readData(Istream& is);

        //- The writeData method for regIOobject write operation
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
