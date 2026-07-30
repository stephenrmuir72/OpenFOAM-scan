/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024-2026 OpenCFD Ltd.
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
    Foam::faGeometryScheme

Description
    Geometry calculation scheme for reduced-order (2D) manifold meshes, as
    used by the finite-area method: the mesh "cells" are surface polygons and
    the mesh "faces" are their bounding edges (lines).

    Reuses all of the standard (3D) machinery of basicFvGeometryScheme
    (weights, delta coefficients, non-orthogonal corrections) and adds the
    2D-manifold specific geometry: point/face area normals, edge length
    vectors (Le) and the re-orientation of the face areas perpendicular to
    the local tangential direction.

    This keeps the 2D-manifold responsibility out of basicFvGeometryScheme
    (single-responsibility) and is selected automatically for meshes with
    nTopologicalD() == 2, or explicitly via the "geometry" scheme entry.

SourceFiles
    faGeometryScheme.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faGeometryScheme_H
#define Foam_faGeometryScheme_H

#include "basicFvGeometryScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class faGeometryScheme Declaration
\*---------------------------------------------------------------------------*/

class faGeometryScheme
:
    public basicFvGeometryScheme
{
    // Private Member Functions

        //- No copy construct
        faGeometryScheme(const faGeometryScheme&) = delete;

        //- No copy assignment
        void operator=(const faGeometryScheme&) = delete;


public:

    //- Runtime type information
    TypeName("fa");


    // Constructors

        //- Construct from mesh
        faGeometryScheme(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~faGeometryScheme() = default;


    // Member Functions

        //- Do what is necessary if the mesh has moved.
        //  Adds the reduced-order (2D) manifold geometry on top of the
        //  standard geometry provided by basicFvGeometryScheme.
        virtual void movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
