/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    engineFoam

Description
    Transient solver for compressible, turbulent engine flow with a spray
    particle cloud incl. optional, auxiliary 3D mesh motion.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "motionSolver.H"
#include "engineTime.H"
#include "engineMesh.H"
#include "turbulentFluidThermoModel.H"
#include "basicSprayCloud.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "radiationModel.H"
#include "SLGThermo.H"
#include "pimpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for compressible, turbulent engine flow"
        " with a spray particle cloud incl. additional, auxiliary mesh motion."
    );

    #define CREATE_TIME createEngineTime.H
    #define CREATE_MESH createEngineMesh.H
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createEngineTime.H"
    #include "createEngineMesh.H"
    #include "createControl.H"
    #include "readEngineTimeControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createRhoUf.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"
    #include "initContinuityErrs.H"
    #include "startSummary.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

    while (runTime.run())
    {
        #include "readEngineTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        // 'deltaTheta' normally much informative than 'deltaT'. Possibly more
        // appropriate to implement next to 'deltaT' instead of here?
        if (runTime.controlDict().get< bool >("adjustTimeStep"))
        {
            Info<< "deltaCAD = " << runTime.deltaTheta() << runTime.unit()
            << endl;
        }

        Info<< "Engine time = " << runTime.theta() << runTime.unit()
        << endl;

        // Default behaviour - only 1D motion
        Info<< "Engine mesh motion (1D)" << endl;
        mesh.move();

        // Optional behaviour - aux. 3D motion, e.g., for usable canted valves
        if (auxDynamicMesh)
        {
            Info<< "Auxiliary mesh motion (3D)" << endl;
            mesh.movePoints(motionPtr->newPoints());
        }

        parcels.evolve();

        #include "rhoEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "YEqn.H"
            #include "EEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        #include "logSummary.H"

        rho = thermo.rho();

        if (runTime.write())
        {
            combustion->Qdot()().write();
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
