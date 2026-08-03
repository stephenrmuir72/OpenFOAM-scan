/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2021-2026 OpenCFD Ltd.
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

Application
    checkFaMesh (finiteArea)

Description
    Check a finite-area mesh.

    This is a thin wrapper around the finiteArea (reduced-order fvMesh)
    faMesh: read-construct the area mesh from the on-disk faceLabels and
    report the finite-area-specific mesh checks (faMeshTools::printMeshChecks).

    Topological validity of the reduced-order region is now covered by
    checkMesh -region <faRegion>, so this utility keeps only the FA-specific
    report. The legacy -geometryOrder option is dropped (finiteArea uses the
    runtime-selected faGeometryScheme).

Original Authors
    Zeljko Tukovic, FAMENA
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "argList.H"
#include "faMesh.H"
#include "faMeshTools.H"
#include "polyMesh.H"
#include "regionProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Check a finite-area mesh"
    );

    argList::addBoolOption
    (
        "write-vtk",
        "Write mesh as a vtp (vtk) file for display or debugging"
    );
    argList::addBoolOption
    (
        "write-edges-vtk",
        "Write mesh edges as obj files (one per processor)",
        true  // advanced option (debugging only)
    );

    #include "addRegionOption.H"
    #include "addAllFaRegionOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"

    // Handle area region selections
    #include "getAllFaRegionOptions.H"

    // ------------------------------------------------------------------------

    #include "createNamedMesh.H"

    if (args.found("write-vtk") || args.found("write-edges-vtk"))
    {
        Info<< "Note: VTK output is not yet ported to finiteArea"
            << " - ignoring -write-vtk" << nl << endl;
    }

    for (const word& areaName : areaRegionNames)
    {
        Info<< "Create faMesh";
        if (!polyMesh::regionName(areaName).empty())
        {
            Info<< " [" << areaName << "]";
        }
        Info<< " for time = " << runTime.timeName() << nl;

        // Presence probe (same faceLabels check the legacy TryNew used);
        // preserve the graceful skip when the region was never written.
        if (!faMesh::hasMeshFiles(areaName, mesh))
        {
            Info<< "    ...failed to create area-mesh";
            if (!polyMesh::regionName(areaName).empty())
            {
                Info<< " [" << areaName << "]";
            }
            Info<< endl;
            continue;
        }
        else
        {
            Info<< endl;
        }

        // Read-construct from the PERSISTED reduced-order mesh (as
        // createFaMesh.H / the finiteArea solvers do), so the report is
        // faithful to what a solver loads and preserves the definition
        // boundary patches instead of re-deriving them from the parent.
        // Using the finite-area registry with the loop region name maps a
        // default region (polyMesh::defaultRegion) and named regions alike
        // onto the same on-disk layout that makeFaMesh wrote.
        faMesh aMesh
        (
            IOobject
            (
                areaName,
                runTime.timeName(),
                faMesh::Registry(mesh),
                IOobject::MUST_READ
            ),
            false
        );
        aMesh.init(true);   // initialise all (lower levels and current)

        // Mesh information (verbose)
        faMeshTools::printMeshChecks(aMesh, 1);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
