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
    makeFaMesh (finiteArea)

Description
    A mesh generator for the finiteArea (reduced-order fvMesh) faMesh.

    This is a thin wrapper around the finiteArea generation constructor
    faMesh(areaName, polyMesh, faMeshDefinition): read the faMeshDefinition,
    construct the reduced-order mesh, report the mesh checks and write the
    faceLabels addressing.

    Unlike the legacy finiteArea makeFaMesh, this utility does NOT act like
    decomposePar in parallel: procAddressing creation and field decomposition
    are handled by decomposePar / redistributePar directly on the
    reduced-order mesh.

Original Authors
    Zeljko Tukovic, FAMENA
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "argList.H"
#include "faMesh.H"
#include "faMeshTools.H"
#include "IOdictionary.H"
#include "regionProperties.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "A mesh generator for the finiteArea (reduced-order) faMesh"
    );
    argList::addOption
    (
        "empty-patch",
        "name",
        "Specify name for a default empty patch",
        false  // An advanced option, but not enough to worry about that
    );
    argList::addOption("dict", "file", "Alternative faMeshDefinition");

    argList::addDryRunOption
    (
        "Create but do not write"
    );
    argList::addBoolOption
    (
        "write-vtk",
        "Write mesh as a vtp (vtk) file for display or debugging"
    );

    #include "addRegionOption.H"
    #include "addAllFaRegionOptions.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    // Handle area region selections
    #include "getAllFaRegionOptions.H"

    if (args.found("write-vtk"))
    {
        Info<< "Note: VTK output is not yet ported to finiteArea"
            << " - ignoring -write-vtk" << nl << endl;
    }

    // ------------------------------------------------------------------------

    for (const word& areaRegionName : areaRegionNames)
    {
        // Reading faMeshDefinition dictionary
        #include "findMeshDefinitionDict.H"

        if (!faMeshDefinitionPtr)
        {
            if (args.dryRun())
            {
                break;
            }
            else
            {
                FatalErrorInFunction
                    << "Did not find area-mesh definition"<< nl
                    << exit(FatalError);
            }
        }

        auto& meshDefDict = (*faMeshDefinitionPtr);


        // Inject/overwrite name for optional 'empty' patch
        if (word name; args.readIfPresent("empty-patch", name))
        {
            meshDefDict.add("emptyPatch", name, true);
        }

        Info<< "Create area-mesh";
        if (!polyMesh::regionName(areaRegionName).empty())
        {
            Info<< " [" << areaRegionName << "]";
        }
        Info<< " with polyMesh at time = " << runTime.timeName() << nl;

        // Create
        faMesh aMesh(areaRegionName, mesh, meshDefDict);

        // Mesh information (less verbose)
        faMeshTools::printMeshChecks(aMesh, 0);

        if (args.dryRun())
        {
            Info<< "\ndry-run: not writing mesh\n" << nl;
        }
        else
        {
            // More precision (for points data)
            IOstream::minPrecision(10);

            Info<< nl << "Write finite area mesh";
            if (!polyMesh::regionName(aMesh.areaName()).empty())
            {
                Info<< " [" << aMesh.areaName() << "]";
            }
            Info<< nl;

            aMesh.write();

            Info<< endl;
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
