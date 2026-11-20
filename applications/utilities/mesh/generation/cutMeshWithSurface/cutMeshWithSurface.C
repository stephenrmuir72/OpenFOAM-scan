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

    Each cell the surface crosses is cut along the loop where the surface meets
    its boundary. A cell the cutter cannot divide is reported with its index and
    left untouched. The mesh is still written, and the utility exits non-zero so
    that the gap cannot pass unnoticed.

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

#include <algorithm>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    argList::addNote
    (
        "Cut a mesh by a surface along each cell's boundary loop"
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
    const dictionary& castellatedControls
    (
        meshDict.subDict("castellatedMeshControls")
    );

    if (castellatedControls.found("locationsInMesh"))
    {
        FatalErrorInFunction
            << "locationsInMesh is not supported; use one locationInMesh point"
            << exit(FatalError);
    }

    const bool hasLocationInMesh =
        castellatedControls.found("locationInMesh");
    point locationInMesh(Zero);

    if (hasLocationInMesh)
    {
        locationInMesh = castellatedControls.get<point>("locationInMesh");
    }

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

    //- The snap, and it is a mesh motion rather than a cut setting. A surface
    //- vertex within snapPlanes of a background plane makes the surface graze
    //- that plane, which is the annulus's multiply-connected face; the plane is
    //- then moved past the vertex by the vertex's own gap plus the clearance,
    //- each plane on its own, because the two vertices of an antipodal pair want
    //- opposite moves and no rigid move satisfies both (07 section 32.2). Zero
    //- disables it. A plane is moved only where the mesh's own points lie on it,
    //- so a refined background needs no assumption beyond that.
    const scalar snapPlanes
    (
        meshDict.getOrDefault<scalar>("snapPlanes", 0.0)
    );

    const scalar snapPlanesClearance
    (
        meshDict.getOrDefault<scalar>("snapPlanesClearance", 0.005)
    );

    if (snapPlanes > 0)
    {
        const pointField& pts = mesh.points();

        pointField moved(pts);

        const scalar coordTol =
            1e-9*mag(mesh.bounds().max() - mesh.bounds().min());

        label nMovedPlanes = 0;

        Info<< "Snap planes" << nl << "-----------" << nl;

        for (direction dir = 0; dir < 3; dir++)
        {
            scalarField all(pts.size());

            forAll(pts, i)
            {
                all[i] = pts[i][dir];
            }

            std::sort(all.begin(), all.end());

            //- The distinct coordinates along this axis, and its finest spacing.
            scalarField coords(all.size());

            label nC = 0;

            forAll(all, i)
            {
                if (nC == 0 || all[i] - coords[nC - 1] > coordTol)
                {
                    coords[nC++] = all[i];
                }
            }

            coords.setSize(nC);

            scalar cellSize = GREAT;

            for (label i = 1; i < nC; i++)
            {
                cellSize = min(cellSize, coords[i] - coords[i - 1]);
            }

            if (nC < 2 || cellSize <= VSMALL)
            {
                continue;
            }

            //- The move each plane needs: the largest over the vertices near it.
            scalarField move(nC, 0.0);

            forAll(surfaces, si)
            {
                const pointField& spts = surfaces[si].surface().points();

                forAll(spts, vi)
                {
                    const scalar v = spts[vi][dir];

                    label best = -1;
                    scalar bestD = GREAT;

                    for (label i = 0; i < nC; i++)
                    {
                        const scalar d = mag(coords[i] - v);

                        if (d < bestD)
                        {
                            bestD = d;
                            best = i;
                        }
                    }

                    if (best >= 0 && bestD <= snapPlanes*cellSize)
                    {
                        const scalar want
                        (
                            (v > coords[best] ? 1 : -1)
                          * (bestD + snapPlanesClearance*cellSize)
                        );

                        if (mag(want) > mag(move[best]))
                        {
                            move[best] = want;
                        }
                    }
                }
            }

            //- Apply, to every point that lies on a moved plane.
            for (label j = 0; j < nC; j++)
            {
                if (mag(move[j]) <= VSMALL)
                {
                    continue;
                }

                forAll(pts, i)
                {
                    if (mag(pts[i][dir] - coords[j]) <= coordTol)
                    {
                        moved[i][dir] = coords[j] + move[j];
                    }
                }

                nMovedPlanes++;

                Info<< "    " << ("xyz"[dir]) << " = " << coords[j]
                    << " moved by " << move[j] << " to "
                    << coords[j] + move[j] << nl;
            }
        }

        Info<< "    " << nMovedPlanes << " plane(s) moved, "
            << mag(snapPlanes)*100 << "% trigger, "
            << mag(snapPlanesClearance)*100 << "% clearance" << nl << endl;

        if (nMovedPlanes)
        {
            mesh.movePoints(moved);
        }
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

    const dictionary* materialRegions =
        meshDict.found("materialRegions")
      ? &meshDict.subDict("materialRegions")
      : nullptr;

    cutCellMesher mesher(mesh, surfaces, materialRegions);

    //- The collapse may move a point by at most this fraction of the point's
    //- shortest edge. The trigger alone does not bound it: a point pulled by
    //- several planes is resolved by a least-squares solve that runs away when
    //- the planes are nearly parallel, which is what a coarsely tessellated
    //- curved surface gives. Measured before the guard existed: a 320-triangle
    //- sphere in a 16-cubed background moved a point by a third of a cell, and
    //- the assembly aborted with no mesh written. Zero disables the guard.
    scalar maxSnapMoveFraction
    (
        meshDict.getOrDefault<scalar>("maxSnapMoveFraction", 0.05)
    );

    if (maxSnapMoveFraction < 0 || maxSnapMoveFraction > 0.5)
    {
        WarningInFunction
            << "maxSnapMoveFraction " << maxSnapMoveFraction << " is outside the"
            << " sensible range 0 to 0.5, and is clamped to it" << endl;

        maxSnapMoveFraction = min(max(maxSnapMoveFraction, 0.0), 0.5);
    }

    //- The cut is cutCellByLoop, which traces where the surface meets the
    //- cell's boundary. It is the only cut, so cutMethod is not a choice: loop
    //- is accepted and changes nothing, and plane is refused.
    const word cutMethod(meshDict.getOrDefault<word>("cutMethod", "loop"));

    if (cutMethod == "plane")
    {
        FatalErrorInFunction
            << "cutMethod plane has been removed, because the loop is the only"
            << " cut. Use cutMethod loop, or remove the keyword, in "
            << meshDict.name() << exit(FatalError);
    }

    if (cutMethod != "loop")
    {
        FatalErrorInFunction
            << "cutMethod '" << cutMethod << "' is not loop" << exit(FatalError);
    }

    if (meshDict.getOrDefault<bool>("cutCellLoopDebug", false))
    {
        cutCellLoopDebug(true);
    }

    //- One cell to dump for replay, by label. -1 is none.
    mesher.dumpCell(meshDict.getOrDefault<label>("dumpCell", -1));
    mesher.vertexSnapFraction(vertexSnapFraction);
    mesher.maxSnapMoveFraction(maxSnapMoveFraction);

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

    if (hasLocationInMesh)
    {
        if (!map.valid())
        {
            FatalErrorInFunction
                << "locationInMesh was specified, but the surface cut made no "
                << "topology change" << exit(FatalError);
        }

        mesher.retainLocation(locationInMesh);
    }

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
