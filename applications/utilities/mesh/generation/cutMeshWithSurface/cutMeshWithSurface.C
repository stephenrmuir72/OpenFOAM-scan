/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
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
    cutMeshWithSurface

Description
    Cut a mesh by a surface, replacing every cell the surface passes through
    with one sub-cell on each side, in a single topology change.

    Stage S1 of the cut-cell Cartesian mesh generator. It reads a castellated
    background mesh from constant/polyMesh and writes the cut mesh to the time
    the case is at, so that both stages stay available and the parent mapping of
    every sub-cell can be checked against the background it came from.

    The result carries a cellZone per region and a faceZone per interface, and
    no baffles. Converting the faceZone into paired patches is left to
    createBaffles, and splitting into per-region meshes to splitMeshRegions.

    This is the first slice: it handles cells that one plane of the surface
    passes through. Cells crossed by more than one plane, which happens where a
    surface edge or corner passes through a cell, are reported with their index
    and left untouched. The mesh is still written, and the utility exits
    non-zero so that the gap cannot pass unnoticed.

    The dictionary is cutCellMeshDict, the same one cutCellCartesian reads, so
    the geometry is specified once.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "IOdictionary.H"
#include "PtrList.H"
#include "cutCellMesher.H"
#include "cutCellSurface.H"
#include "mapPolyMesh.H"
#include "error.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    argList::addNote
    (
        "Cut a mesh by a surface, one plane per cell, in one topology change"
    );

    argList::addOption("dict", "file", "Alternative cutCellMeshDict");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    //- The cut mesh goes to the time the case is at, so the background mesh in
    //- constant stays available. That is also what S4 needs, and it makes a
    //- re-run legible: reading a mesh that has already been cut and cutting it
    //- again would degrade it a little more each time (07 section 6.3).
    const word writeInstance(runTime.timeName());

    if (mesh.pointsInstance() != "constant")
    {
        FatalErrorInFunction
            << "The mesh was read from instance '" << mesh.pointsInstance()
            << "', but the background mesh must be the one in constant/." << nl
            << "  This run would write the cut mesh to '" << writeInstance
            << "', so re-running would cut a mesh that has already been cut."
            << nl
            << "  Remove " << writeInstance << "/polyMesh, or set startTime to"
            << " a time that does not exist yet." << exit(FatalError);
    }

    const word dictName("cutCellMeshDict");
    #include "setSystemMeshDictionaryIO.H"

    const IOdictionary meshDict(dictIO);

    const dictionary& geometryDict(meshDict.subDict("geometry"));

    Info<< "Geometry" << nl << "--------" << endl;

    PtrList<cutCellSurface> surfaces(geometryDict.size());

    label surfI = 0;

    for (const entry& e : geometryDict)
    {
        surfaces.set
        (
            surfI++,
            new cutCellSurface(e.keyword(), e.dict(), runTime)
        );
    }

    //- How close a cut crossing has to be to the endpoint of a cell edge, as a
    //- fraction of that edge's length, to be collapsed onto it. 05 section
    //- 7.1.2 gives 0.05 as the default and 0.10 as the limit, and zero disables
    //- the collapse.
    scalar vertexSnapFraction
    (
        meshDict.getOrDefault<scalar>("vertexSnapFraction", 0.05)
    );

    if (vertexSnapFraction < 0 || vertexSnapFraction > 0.10)
    {
        WarningInFunction
            << "vertexSnapFraction " << vertexSnapFraction << " is outside the"
            << " allowable range 0 to 0.10, and is clamped to it" << endl;

        vertexSnapFraction = min(max(vertexSnapFraction, 0.0), 0.10);
    }

    cutCellMesher mesher(mesh, surfaces);

    mesher.vertexSnapFraction(vertexSnapFraction);

    //- The sliver thresholds, either being zero disabling the merge. The sliver
    //- ratio defaults to zero so that a case which does not ask for a merge is
    //- untouched; 0.01 is the value the merge was designed around.
    const scalar sliverVolumeRatio
    (
        meshDict.getOrDefault<scalar>("sliverVolumeRatio", 0.0)
    );

    const scalar minVolumeRatio
    (
        meshDict.getOrDefault<scalar>("minVolumeRatio", 0.30)
    );

    mesher.sliverVolumeRatio(sliverVolumeRatio);
    mesher.minVolumeRatio(minVolumeRatio);

    Info<< "Sliver merge thresholds       : ratio " << sliverVolumeRatio
        << ", minVolumeRatio " << minVolumeRatio
        << (sliverVolumeRatio > 0 && minVolumeRatio > 0 ? " (on)" : " (off)")
        << nl << endl;

    autoPtr<mapPolyMesh> map = mesher.cut();

    const cutCellMesher::report& info = mesher.info();

    if (map.valid())
    {
        mesher.writeBackgroundCell();
    mesher.writeBackgroundParents();
        mesher.writeSnapDisplacement();

        mesh.setInstance(writeInstance);

        //- The zones are registered without a write option, so mesh.write()
        //- leaves them alone. They also have to follow the mesh to its new
        //- instance, or they would land in constant/polyMesh without the mesh
        //- they describe.
        mesh.cellZones().instance() = writeInstance;
        mesh.faceZones().instance() = writeInstance;
        mesh.cellZones().writeOpt() = IOobject::AUTO_WRITE;
        mesh.faceZones().writeOpt() = IOobject::AUTO_WRITE;

        Info<< nl << "Writing cut mesh to instance " << writeInstance << nl
            << endl;

        mesh.write();
    }

    if (info.nUnsupported)
    {
        WarningInFunction
            << info.nUnsupported << " cell(s) could not be cut in this slice"
            << " and were left untouched, so this is not a complete cut of the"
            << " surface. They are listed above. The mesh was written anyway,"
            << " so the state can be inspected." << endl;

        return 1;
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
