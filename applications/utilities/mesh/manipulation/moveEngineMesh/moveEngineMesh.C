/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    moveEngineMesh

Group
    grpMeshManipulationUtilities

Description
    Solver for moving meshes for engine calculations  incl. optional, auxiliary
    3D mesh motion.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "motionSolver.H"
#include "engineTime.H"
#include "engineMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "A solver utility for moving meshes for engine calculations incl."
        " additional, auxiliary mesh motion"
    );

    #include "setRootCase.H"
    #include "createEngineTime.H"
    #include "createEngineMesh.H"

    // Check for aux. motion switch. If true, any 3D motion involved on top.
    // Possibly not the best way to access engineDict/'auxDynamicMesh' here.
    IOdictionary engineGeometryDict
    (
        IOobject
        (
            "engineGeometry",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    const bool auxDynamicMesh
               = engineGeometryDict.getOrDefault<bool>("auxDynamicMesh", false);

    // 'motionSolver' needs 'dynamicMeshDict' so tied directly to
    // 'auxDynamicMesh' presence
    autoPtr<motionSolver> motionPtr;
    if (auxDynamicMesh)
    {
        Info<< "\nSelecting auxiliary motion solver" << endl;
        motionPtr = motionSolver::New(mesh);
    }

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        // deltaTheta normally more informative than deltaT
        // More appropriate  implement next to deltaT???
        if (runTime.controlDict().get< bool >("adjustTimeStep"))
        {
            Info<< "deltaCAD = " << runTime.deltaTheta() << runTime.unit()
            << endl;
        }

        Info<< "Engine time = " << runTime.theta() << runTime.unit()
        << endl;

        // Default benaviour - no aux. mesh motion
        Info<< "Engine mesh motion (1D)" << endl;
        mesh.move();

        // Optional behaviour - also 3D motion, e.g., for canted valves
        if (auxDynamicMesh)
        {
            Info<< "Auxiliary mesh motion (3D)" << endl;
            mesh.movePoints(motionPtr->newPoints());
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
