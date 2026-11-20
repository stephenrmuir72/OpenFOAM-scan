/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.com
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    cutCellCartesian

Group
    grpMeshGenerationUtilities

Description
    Castellation-only background generator for the cut-cell Cartesian mesh
    generator.

    A refined, unsnapped Cartesian background is produced by reusing
    OpenFOAM's snappyHexMesh refinement library directly:

    - surface, feature-edge and volume (shell) refinement are supplied by
      refinementSurfaces, refinementFeatures and shellSurfaces;
    - splitting and 2:1 balancing are performed by hexRef8, reached through
      meshRefinement and snappyRefineDriver;
    - the snap and layer drivers are never constructed, so background points
      are not displaced from the Cartesian lattice;
    - refinement is driven by snappyRefineDriver::doRefineOnly, which performs
      the refinement stages of doRefine but not its region extraction. All
      cells are therefore retained, and no baffles are introduced at the
      surface, so the whole Cartesian background survives to be cut later.

    Castellation in snappyHexMesh normally also removes the cells that are not
    reachable from the keep points and baffles the surface. Neither is wanted
    here: the cutting stage introduces the surface, and for conjugate heat
    transfer both sides of it must remain in the mesh.

    The resulting mesh carries the per-cell cellLevel, per-point pointLevel
    and level0Edge fields, which the later cutting stages use for O(1) cell
    identification.

    The dictionary names follow snappyHexMeshDict so that the library readers
    work unmodified; only 'geometry' and 'castellatedMeshControls' are
    required.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "IOdictionary.H"
#include "IOmanip.H"

#include "searchableSurfaces.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "shellSurfaces.H"
#include "surfaceZonesInfo.H"
#include "wallPolyPatch.H"

#include "meshRefinement.H"
#include "refinementHistory.H"
#include "refinementParameters.H"
#include "snapParameters.H"
#include "snappyRefineDriver.H"

#include "decompositionMethod.H"
#include "fvMeshDistribute.H"

#include "coordSetWriter.H"
#include "surfaceWriter.H"
#include "vtkCoordSetWriter.H"
#include "vtkSurfaceWriter.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

//- Absolute merge distance from a relative merge tolerance
scalar calcMergeDistance(const polyMesh& mesh, const scalar mergeTol)
{
    const boundBox& meshBb = mesh.bounds();
    const scalar mergeDist = mergeTol*meshBb.mag();

    Info<< nl
        << "Overall mesh bounding box  : " << meshBb << nl
        << "Relative tolerance         : " << mergeTol << nl
        << "Absolute matching distance : " << mergeDist << nl
        << endl;

    return mergeDist;
}


//- Create a patch for every surface region and fill the region-to-patch maps.
//  This mirrors the equivalent step in snappyHexMesh so that surface regions
//  become named patches of the castellated background.
void addSurfaceRegionPatches
(
    meshRefinement& meshRefiner,
    const searchableSurfaces& allGeometry,
    const refinementSurfaces& surfaces,
    labelList& globalToMasterPatch,
    labelList& globalToSlavePatch,
    const bool dryRun
)
{
    globalToMasterPatch.setSize(surfaces.nRegions(), -1);
    globalToSlavePatch.setSize(surfaces.nRegions(), -1);

    const labelList& surfaceGeometry = surfaces.surfaces();
    const PtrList<dictionary>& surfacePatchInfo = surfaces.patchInfo();

    Info<< nl
        << "Adding patches for surface regions" << nl
        << "----------------------------------" << nl
        << endl;

    //- Add one patch for a region, using its patchInfo when supplied
    auto addPatch = [&](const label globalRegioni, const word& name) -> label
    {
        if (surfacePatchInfo.set(globalRegioni))
        {
            return meshRefiner.addMeshedPatch
            (
                name,
                surfacePatchInfo[globalRegioni]
            );
        }

        dictionary patchInfo;
        patchInfo.set("type", wallPolyPatch::typeName);

        return meshRefiner.addMeshedPatch(name, patchInfo);
    };

    forAll(surfaceGeometry, surfi)
    {
        const label geomi = surfaceGeometry[surfi];

        const wordList& regNames = allGeometry.regionNames()[geomi];
        const wordList& fzNames = surfaces.surfZones()[surfi].faceZoneNames();

        Info<< surfaces.names()[surfi] << ':' << nl << nl;

        if (fzNames.empty())
        {
            // Normal surface: master and slave share a single patch
            forAll(regNames, i)
            {
                const label globalRegioni = surfaces.globalRegion(surfi, i);

                const label patchi = addPatch(globalRegioni, regNames[i]);

                globalToMasterPatch[globalRegioni] = patchi;
                globalToSlavePatch[globalRegioni] = patchi;

                if (!dryRun)
                {
                    Info<< setf(ios_base::left)
                        << setw(6) << patchi
                        << setw(20)
                        << meshRefiner.mesh().boundaryMesh()[patchi].type()
                        << setw(30) << regNames[i] << nl;
                }
            }
        }
        else
        {
            // Zoned surface: separate master and slave patches plus a faceZone
            forAll(regNames, i)
            {
                const label globalRegioni = surfaces.globalRegion(surfi, i);

                const label masterPatchi = addPatch(globalRegioni, regNames[i]);

                const word slaveName = regNames[i] + "_slave";
                const label slavePatchi = addPatch(globalRegioni, slaveName);

                globalToMasterPatch[globalRegioni] = masterPatchi;
                globalToSlavePatch[globalRegioni] = slavePatchi;

                if (!dryRun)
                {
                    Info<< setf(ios_base::left)
                        << setw(6) << masterPatchi
                        << setw(20)
                        << meshRefiner.mesh().boundaryMesh()[masterPatchi].type()
                        << setw(30) << regNames[i] << nl
                        << setw(6) << slavePatchi
                        << setw(20)
                        << meshRefiner.mesh().boundaryMesh()[slavePatchi].type()
                        << setw(30) << slaveName << nl;
                }
            }

            const polyBoundaryMesh& pbm = meshRefiner.mesh().boundaryMesh();

            forAll(fzNames, fzi)
            {
                const label globalRegioni = surfaces.globalRegion(surfi, fzi);

                meshRefiner.addFaceZone
                (
                    fzNames[fzi],
                    pbm[globalToMasterPatch[globalRegioni]].name(),
                    pbm[globalToSlavePatch[globalRegioni]].name(),
                    surfaces.surfZones()[surfi].faceType()
                );
            }
        }

        Info<< nl;
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Cut-cell Cartesian background generator: castellation only, "
        "no snapping, no layer addition, no cell removal"
    );

    #include "addOverwriteOption.H"

    argList::addDryRunOption
    (
        "Check the case set-up only, without changing the mesh"
    );

    argList::addOption("dict", "file", "Alternative cutCellMeshDict");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const bool dryRun = args.dryRun();
    const bool overwrite = args.found("overwrite");

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const word dictName("cutCellMeshDict");
    #include "setSystemMeshDictionaryIO.H"

    const IOdictionary meshDict(dictIO);

    // These two entries are required.
    // Held by value rather than by reference: meshRefinement::subDict returns a
    // reference into meshDict, but gcc-13's -Wdangling-reference cannot see that
    // and flags the binding, and FOAM_NO_DANGLING_REFERENCE only supplies
    // [[gnu::no_dangling]] for gcc >= 14. Copying these small sub-dictionaries
    // keeps the build log clean without changing behaviour.
    const dictionary geometryDict
    (
        meshRefinement::subDict(meshDict, "geometry", dryRun)
    );

    const dictionary refineDict
    (
        meshRefinement::subDict(meshDict, "castellatedMeshControls", dryRun)
    );

    // The remaining entries are only needed by parts of the driver that are
    // not exercised here, so they may be omitted
    const dictionary& motionDict =
        meshDict.subOrEmptyDict("meshQualityControls");

    const dictionary& snapDict =
        meshDict.subOrEmptyDict("snapControls");

    const scalar mergeDist = calcMergeDistance
    (
        mesh,
        meshRefinement::get<scalar>(meshDict, "mergeTolerance", dryRun)
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Writers, only needed to satisfy the snappyRefineDriver constructor.
    // The vtk writers are registered in libmeshTools and libsurfMesh, both of
    // which are linked.
    autoPtr<coordSetWriter> setFormatter
    (
        coordSetWriter::New
        (
            coordSetWriters::vtkWriter::typeName,
            meshDict.subOrEmptyDict("formatOptions").optionalSubDict
            (
                coordSetWriters::vtkWriter::typeName
            )
        )
    );

    refPtr<surfaceWriter> surfFormatter;
    surfFormatter = surfaceWriter::New
    (
        surfaceWriters::vtkWriter::typeName,
        meshDict.subOrEmptyDict("formatOptions").optionalSubDict
        (
            surfaceWriters::vtkWriter::typeName
        )
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "Reading geometry." << endl;

    searchableSurfaces allGeometry
    (
        IOobject
        (
            "abc",                      // dummy name
            mesh.time().constant(),     // instance
            "triSurface",               // local
            mesh.time(),                // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        geometryDict,
        meshDict.getOrDefault("singleRegionName", true)
    );

    Info<< "Reading refinement surfaces." << endl;

    refinementSurfaces surfaces
    (
        allGeometry,
        meshRefinement::subDict(refineDict, "refinementSurfaces", dryRun),
        refineDict.getOrDefault("gapLevelIncrement", 0),
        dryRun
    );

    Info<< "Reading refinement regions." << endl;

    shellSurfaces shells
    (
        allGeometry,
        refineDict.subOrEmptyDict("refinementRegions"),
        dryRun
    );

    surfaces.setMinLevelFields(shells);

    // No limit regions: the whole Cartesian background is retained
    shellSurfaces limitShells
    (
        allGeometry,
        refineDict.subOrEmptyDict("limitRegions"),
        dryRun
    );

    Info<< "Reading features." << endl;

    autoPtr<refinementFeatures> featuresPtr;

    if (refineDict.found("features"))
    {
        featuresPtr.reset
        (
            new refinementFeatures
            (
                mesh,
                PtrList<dictionary>
                (
                    meshRefinement::lookup(refineDict, "features", dryRun)
                ),
                dryRun
            )
        );
    }
    else
    {
        featuresPtr.reset(new refinementFeatures(mesh, PtrList<dictionary>(), dryRun));
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Main refinement engine. hexRef8 lives inside and keeps cellLevel,
    // pointLevel and level0Edge in memory.
    meshRefinement meshRefiner
    (
        mesh,
        mergeDist,              // tolerance used in sorting coordinates
        overwrite,              // overwrite mesh files?
        surfaces,               // surface intersection refinement
        featuresPtr(),          // feature edge/point based refinement
        shells,                 // volume (inside/outside) refinement
        limitShells,            // limit of volume refinement
        labelList(),            // initial faces to test
        meshRefinement::CASTELLATED,
        dryRun
    );

    if (!dryRun)
    {
        // Recalculate intersections for all faces
        labelList allFaces(mesh.nFaces());

        forAll(allFaces, facei)
        {
            allFaces[facei] = facei;
        }

        meshRefiner.updateIntersections(allFaces);
    }

    meshRefiner.printMeshInfo(false, "Initial mesh", true);

    const refinementParameters refineParams(refineDict, dryRun);
    const snapParameters snapParams(snapDict, dryRun);

    // Decomposition and distribution engines. They are only exercised by the
    // load-balancing paths, which are not triggered here (the cell budgets
    // are large and balancing is not requested), but the driver requires them
    // to exist.
    dictionary decomposeDict;

    if (meshDict.found("decomposeParDict"))
    {
        decomposeDict = meshDict.subDict("decomposeParDict");
    }
    else
    {
        // Serial default, matching snappyHexMesh: 'none' is the dummy
        // decomposition method and needs no further entries. It is only a
        // placeholder, because it is exercised solely by the load-balancing
        // paths, which are not reached with the cell budgets used here.
        decomposeDict.add("method", "none");
        decomposeDict.add("numberOfSubdomains", label(1));
    }

    autoPtr<decompositionMethod> decomposerPtr
    (
        decompositionMethod::New(decomposeDict)
    );
    decompositionMethod& decomposer = *decomposerPtr;

    fvMeshDistribute distributor(mesh);

    // From global region number to mesh patch
    labelList globalToMasterPatch;
    labelList globalToSlavePatch;

    addSurfaceRegionPatches
    (
        meshRefiner,
        allGeometry,
        surfaces,
        globalToMasterPatch,
        globalToSlavePatch,
        dryRun
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    {
        snappyRefineDriver refineDriver
        (
            meshRefiner,
            decomposer,
            distributor,
            globalToMasterPatch,
            globalToSlavePatch,
            setFormatter(),
            surfFormatter,
            dryRun
        );

        if (!overwrite && !dryRun)
        {
            ++runTime;
        }

        // Refinement only: no region extraction, so all cells are retained,
        // no points are moved and no baffles are introduced. Note that
        // co-planar face merging is not performed at all, which preserves the
        // split-hex cell shape of the Cartesian background.
        refineDriver.doRefineOnly
        (
            refineDict,
            refineParams,
            snapParams,
            motionDict
        );
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    if (!dryRun)
    {
        const word writeInstance(meshRefiner.timeName());

        // Keep the level fields with the mesh they describe
        meshRefiner.setInstance(writeInstance);

        Info<< nl
            << "Writing refined background mesh to instance "
            << writeInstance << nl << endl;

        // Remove any stale refinement history before writing
        refinementHistory::removeFiles(mesh);

        // Writes the mesh, and cellLevel/pointLevel/level0Edge with it
        meshRefiner.write
        (
            meshRefinement::debugType(0),
            meshRefinement::writeType(meshRefinement::WRITEMESH),
            runTime.path()/writeInstance
        );
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
