/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2026 OpenCFD Ltd.
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
    mergeMeshes

Group
    grpMeshManipulationUtilities

Description
    Merges two meshes.

    With -fields, vol and surface fields are read from both cases and the
    combined fields are written with the mesh. The meshes are appended.
    Coincident faces are not coupled.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "mergePolyMesh.H"
#include "topoSet.H"
#include "processorMeshes.H"
#include "fvMesh.H"
#include "fvMeshAdder.H"
#include "fvMeshTools.H"
#include "faceCoupleInfo.H"
#include "IOobjectList.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "processorPolyPatch.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

// Vol and surface fields read from one mesh
struct meshFields
{
    PtrList<volScalarField> vsFlds;
    PtrList<volVectorField> vvFlds;
    PtrList<volSphericalTensorField> vstFlds;
    PtrList<volSymmTensorField> vsymtFlds;
    PtrList<volTensorField> vtFlds;

    PtrList<surfaceScalarField> ssFlds;
    PtrList<surfaceVectorField> svFlds;
    PtrList<surfaceSphericalTensorField> sstFlds;
    PtrList<surfaceSymmTensorField> ssymtFlds;
    PtrList<surfaceTensorField> stFlds;

    void read(const fvMesh& mesh)
    {
        IOobjectList objects(mesh, mesh.time().timeName());

        // Volume fields
        ReadFields(mesh, objects, vsFlds);
        ReadFields(mesh, objects, vvFlds);
        ReadFields(mesh, objects, vstFlds);
        ReadFields(mesh, objects, vsymtFlds);
        ReadFields(mesh, objects, vtFlds);

        // Surface fields
        ReadFields(mesh, objects, ssFlds);
        ReadFields(mesh, objects, svFlds);
        ReadFields(mesh, objects, sstFlds);
        ReadFields(mesh, objects, ssymtFlds);
        ReadFields(mesh, objects, stFlds);
    }

    void setInstance(const word& inst)
    {
        auto setInst = [&inst](auto& flds)
        {
            for (auto& fld : flds)
            {
                fld.instance() = inst;
            }
        };

        setInst(vsFlds);
        setInst(vvFlds);
        setInst(vstFlds);
        setInst(vsymtFlds);
        setInst(vtFlds);

        setInst(ssFlds);
        setInst(svFlds);
        setInst(sstFlds);
        setInst(ssymtFlds);
        setInst(stFlds);
    }
};


//- Skip writing a master field that has no counterpart to merge.
//  Warn about fields that exist only on the mesh being added.
template<class GeoField>
void suppressUnmapped
(
    PtrList<GeoField>& masterFlds,
    const PtrList<GeoField>& addFlds
)
{
    for (auto& fld : masterFlds)
    {
        bool found = false;

        for (const auto& addFld : addFlds)
        {
            if (addFld.name() == fld.name())
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            // fvMeshAdder leaves this field at the old size
            fld.writeOpt(IOobject::NO_WRITE);
        }
    }

    for (const auto& addFld : addFlds)
    {
        bool found = false;

        for (const auto& fld : masterFlds)
        {
            if (fld.name() == addFld.name())
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            WarningInFunction
                << "Not mapping field " << addFld.name()
                << " - not present on the master mesh" << endl;
        }
    }
}


//- Add zero-valued target fields for fields that only exist on the other mesh,
//  so fvMeshAdder can map/write the combined field.
template<class GeoField, class Type>
void createMissingMasterFields
(
    PtrList<GeoField>& masterFlds,
    const PtrList<GeoField>& addFlds,
    const fvMesh& masterMesh
)
{
    for (const auto& addFld : addFlds)
    {
        bool found = false;

        for (const auto& fld : masterFlds)
        {
            if (fld.name() == addFld.name())
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            Info<< "Creating zero field " << addFld.name()
                << " for mesh field mapping" << endl;

            const label fieldi = masterFlds.size();
            masterFlds.resize(fieldi + 1);

            masterFlds.set
            (
                fieldi,
                new GeoField
                (
                    IOobject
                    (
                        addFld.name(),
                        masterMesh.time().timeName(),
                        masterMesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    masterMesh,
                    dimensioned<Type>("zero", addFld.dimensions(), Zero)
                )
            );
        }
    }
}


void createMissingMasterFields(meshFields& master, const meshFields& added, const fvMesh& masterMesh)
{
    createMissingMasterFields<volScalarField, scalar>(master.vsFlds, added.vsFlds, masterMesh);
    createMissingMasterFields<volVectorField, vector>(master.vvFlds, added.vvFlds, masterMesh);
    createMissingMasterFields<volSphericalTensorField, sphericalTensor>(master.vstFlds, added.vstFlds, masterMesh);
    createMissingMasterFields<volSymmTensorField, symmTensor>(master.vsymtFlds, added.vsymtFlds, masterMesh);
    createMissingMasterFields<volTensorField, tensor>(master.vtFlds, added.vtFlds, masterMesh);

    createMissingMasterFields<surfaceScalarField, scalar>(master.ssFlds, added.ssFlds, masterMesh);
    createMissingMasterFields<surfaceVectorField, vector>(master.svFlds, added.svFlds, masterMesh);
    createMissingMasterFields<surfaceSphericalTensorField, sphericalTensor>(master.sstFlds, added.sstFlds, masterMesh);
    createMissingMasterFields<surfaceSymmTensorField, symmTensor>(master.ssymtFlds, added.ssymtFlds, masterMesh);
    createMissingMasterFields<surfaceTensorField, tensor>(master.stFlds, added.stFlds, masterMesh);
}


void createMissingAddedFields(const meshFields& master, meshFields& added, const fvMesh& meshToAdd)
{
    createMissingMasterFields<volScalarField, scalar>(added.vsFlds, master.vsFlds, meshToAdd);
    createMissingMasterFields<volVectorField, vector>(added.vvFlds, master.vvFlds, meshToAdd);
    createMissingMasterFields<volSphericalTensorField, sphericalTensor>(added.vstFlds, master.vstFlds, meshToAdd);
    createMissingMasterFields<volSymmTensorField, symmTensor>(added.vsymtFlds, master.vsymtFlds, meshToAdd);
    createMissingMasterFields<volTensorField, tensor>(added.vtFlds, master.vtFlds, meshToAdd);

    createMissingMasterFields<surfaceScalarField, scalar>(added.ssFlds, master.ssFlds, meshToAdd);
    createMissingMasterFields<surfaceVectorField, vector>(added.svFlds, master.svFlds, meshToAdd);
    createMissingMasterFields<surfaceSphericalTensorField, sphericalTensor>(added.sstFlds, master.sstFlds, meshToAdd);
    createMissingMasterFields<surfaceSymmTensorField, symmTensor>(added.ssymtFlds, master.ssymtFlds, meshToAdd);
    createMissingMasterFields<surfaceTensorField, tensor>(added.stFlds, master.stFlds, meshToAdd);
}


void suppressUnmapped(meshFields& master, const meshFields& added)
{
    suppressUnmapped(master.vsFlds, added.vsFlds);
    suppressUnmapped(master.vvFlds, added.vvFlds);
    suppressUnmapped(master.vstFlds, added.vstFlds);
    suppressUnmapped(master.vsymtFlds, added.vsymtFlds);
    suppressUnmapped(master.vtFlds, added.vtFlds);

    suppressUnmapped(master.ssFlds, added.ssFlds);
    suppressUnmapped(master.svFlds, added.svFlds);
    suppressUnmapped(master.sstFlds, added.sstFlds);
    suppressUnmapped(master.ssymtFlds, added.ssymtFlds);
    suppressUnmapped(master.stFlds, added.stFlds);
}


bool reorderProcessorPatchesLast(fvMesh& mesh)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const label nPatches = patches.size();

    labelList oldToNewPatch(nPatches, -1);
    labelList newToOldPatch(nPatches, -1);

    label newPatchi = 0;
    bool reordered = false;

    forAll(patches, patchi)
    {
        if (!isA<processorPolyPatch>(patches[patchi]))
        {
            oldToNewPatch[patchi] = newPatchi;
            newToOldPatch[newPatchi] = patchi;
            reordered = reordered || (newPatchi != patchi);
            ++newPatchi;
        }
    }

    forAll(patches, patchi)
    {
        if (isA<processorPolyPatch>(patches[patchi]))
        {
            oldToNewPatch[patchi] = newPatchi;
            newToOldPatch[newPatchi] = patchi;
            reordered = reordered || (newPatchi != patchi);
            ++newPatchi;
        }
    }

    if (!reordered)
    {
        return false;
    }

    const label nInternalFaces = mesh.nInternalFaces();

    labelList oldToNewFace(mesh.nFaces(), -1);
    faceList newFaces(mesh.faces().size());
    labelList newOwner(mesh.faceOwner().size());
    labelList newNeighbour(mesh.faceNeighbour().size());

    for (label facei = 0; facei < nInternalFaces; ++facei)
    {
        oldToNewFace[facei] = facei;
        newFaces[facei] = mesh.faces()[facei];
        newOwner[facei] = mesh.faceOwner()[facei];
        newNeighbour[facei] = mesh.faceNeighbour()[facei];
    }

    labelList patchSizes(nPatches, Zero);
    labelList patchStarts(nPatches, Zero);

    label newFacei = nInternalFaces;

    forAll(newToOldPatch, patchi)
    {
        const label oldPatchi = newToOldPatch[patchi];
        const polyPatch& pp = patches[oldPatchi];

        patchSizes[patchi] = pp.size();
        patchStarts[patchi] = newFacei;

        forAll(pp, patchFacei)
        {
            const label oldFacei = pp.start() + patchFacei;

            oldToNewFace[oldFacei] = newFacei;
            newFaces[newFacei] = mesh.faces()[oldFacei];
            newOwner[newFacei] = mesh.faceOwner()[oldFacei];

            ++newFacei;
        }
    }

    if (newFacei != mesh.nFaces())
    {
        FatalErrorInFunction
            << "Processor patch reorder mapped " << newFacei
            << " faces, but mesh has " << mesh.nFaces()
            << abort(FatalError);
    }

    for (faceZone& fz : mesh.faceZones())
    {
        labelList addressing(fz.size());
        boolList flipMap(fz.flipMap());

        forAll(fz, i)
        {
            addressing[i] = oldToNewFace[fz[i]];
        }

        fz.resetAddressing(addressing, flipMap);
    }

    pointField points(mesh.points());

    fvMeshTools::reorderPatches(mesh, oldToNewPatch, nPatches, false);

    mesh.resetMotion();
    mesh.resetPrimitives
    (
        autoPtr<pointField>::New(std::move(points)),
        autoPtr<faceList>::New(std::move(newFaces)),
        autoPtr<labelList>::New(std::move(newOwner)),
        autoPtr<labelList>::New(std::move(newNeighbour)),
        patchSizes,
        patchStarts,
        true
    );

    return true;
}

} // End anonymous namespace

void getRootCase(fileName& casePath)
{
    casePath.clean();  // Remove unneeded ".."

    if (casePath.empty() || casePath == ".")
    {
        // handle degenerate form and '.'
        casePath = cwd();
    }
    else if (casePath[0] != '/' && casePath.name() == "..")
    {
        // avoid relative cases ending in '..' - makes for very ugly names
        casePath = cwd()/casePath;
        casePath.clean();  // Remove unneeded ".."
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Merge two meshes, and optionally their fields"
    );

    #include "addOverwriteOption.H"

    argList::addArgument("masterCase");
    argList::addOption
    (
        "masterRegion",
        "name",
        "Specify alternative mesh region for the master mesh"
    );

    argList::addArgument("addCase");
    argList::addOption
    (
        "addRegion",
        "name",
        "Specify alternative mesh region for the additional mesh"
    );
    argList::addOption
    (
        "resultTime",
        "time",
        "Specify a time for the resulting mesh"
    );
    argList::addBoolOption
    (
        "fields",
        "Read, merge and write vol and surface fields"
    );

    argList args(argc, argv);
    if (!args.check())
    {
         FatalError.exit();
    }

    const bool overwrite = args.found("overwrite");
    const bool doFields = args.found("fields");

    auto masterCase = args.get<fileName>(1);
    auto addCase = args.get<fileName>(2);

    const word masterRegion =
        args.getOrDefault<word>("masterRegion", polyMesh::defaultRegion);

    const word addRegion =
        args.getOrDefault<word>("addRegion", polyMesh::defaultRegion);

    // Since we don't use argList processor directory detection, add it to
    // the casename ourselves so it triggers the logic inside TimePath.
    const fileName& cName = args.caseName();
    const auto pos = cName.find("processor");
    if (pos != string::npos && pos != 0)
    {
        fileName processorName = cName.substr(pos);
        masterCase += '/' + processorName;
        addCase += '/' + processorName;
    }


    getRootCase(masterCase);
    getRootCase(addCase);

    Info<< "Master:      " << masterCase << "  region " << masterRegion << nl
        << "mesh to add: " << addCase    << "  region " << addRegion << endl;

    #include "createTimes.H"

    if (doFields)
    {
        Info<< "Reading master mesh for time = "
            << runTimeMaster.timeName() << nl;

        Info<< "Create mesh\n" << endl;
        fvMesh masterMesh
        (
            IOobject
            (
                masterRegion,
                runTimeMaster.timeName(),
                runTimeMaster
            )
        );

        Info<< "Reading mesh to add for time = "
            << runTimeToAdd.timeName() << nl;

        Info<< "Create mesh\n" << endl;
        fvMesh meshToAdd
        (
            IOobject
            (
                addRegion,
                runTimeToAdd.timeName(),
                runTimeToAdd
            )
        );

        // Read before any time change so fields come from the current time
        meshFields masterFld;
        meshFields addFld;

        Info<< "Reading fields from master" << endl;
        masterFld.read(masterMesh);

        Info<< "Reading fields from mesh to add" << endl;
        addFld.read(meshToAdd);

        createMissingMasterFields(masterFld, addFld, masterMesh);
        createMissingAddedFields(masterFld, addFld, meshToAdd);
        suppressUnmapped(masterFld, addFld);

        word meshInstance = masterMesh.pointsInstance();

        const bool specifiedInstance =
        (
            !overwrite
         && args.readIfPresent("resultTime", meshInstance)
        );

        if (specifiedInstance)
        {
            runTimeMaster.setTime(instant(meshInstance), 0);
        }
        else if (!overwrite)
        {
            runTimeMaster++;
        }

        Info<< "Writing combined mesh to "
            << runTimeMaster.timeName() << endl;

        // Empty coupling: append points, faces and cells.
        // Coincident faces stay as boundary faces.
        const faceCoupleInfo couples
        (
            masterMesh,
            labelList(),
            meshToAdd,
            labelList(),
            Zero,   // tolerance
            true,   // perfectMatch
            true,   // orderedFaces
            false   // patchDivision
        );

        fvMeshAdder::add(masterMesh, meshToAdd, couples);

        if (reorderProcessorPatchesLast(masterMesh))
        {
            Info<< "Reordering processor patches last" << endl;
        }

        if (overwrite || specifiedInstance)
        {
            masterMesh.setInstance(meshInstance);
        }

        masterFld.setInstance(runTimeMaster.timeName());

        masterMesh.write();
        topoSet::removeFiles(masterMesh);
        processorMeshes::removeFiles(masterMesh);

        Info<< "End\n" << endl;

        return 0;
    }

    Info<< "Reading master mesh for time = " << runTimeMaster.timeName() << nl;

    Info<< "Create mesh\n" << endl;
    mergePolyMesh masterMesh
    (
        IOobject
        (
            masterRegion,
            runTimeMaster.timeName(),
            runTimeMaster
        )
    );

    Info<< "Reading mesh to add for time = " << runTimeToAdd.timeName() << nl;
    Info<< "Create mesh\n" << endl;
    polyMesh meshToAdd
    (
        IOobject
        (
            addRegion,
            runTimeToAdd.timeName(),
            runTimeToAdd
        )
    );

    word meshInstance = masterMesh.pointsInstance();

    const bool specifiedInstance =
    (
        !overwrite
     && args.readIfPresent("resultTime", meshInstance)
    );

    if (specifiedInstance)
    {
        runTimeMaster.setTime(instant(meshInstance), 0);
    }
    else if (!overwrite)
    {
        runTimeMaster++;
    }

    Info<< "Writing combined mesh to " << runTimeMaster.timeName() << endl;

    masterMesh.addMesh(meshToAdd);
    masterMesh.merge();

    if (overwrite || specifiedInstance)
    {
        masterMesh.setInstance(meshInstance);
    }

    masterMesh.write();
    topoSet::removeFiles(masterMesh);
    processorMeshes::removeFiles(masterMesh);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
