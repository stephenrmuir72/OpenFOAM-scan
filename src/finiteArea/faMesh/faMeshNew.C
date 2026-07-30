/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "faMesh.H"
#include "polyMesh.H"
#include "fileOperation.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

bool Foam::faMesh::hasMeshFiles(const word& areaName, const polyMesh& pMesh)
{
    // The reduced-order faMesh persists only its parent coupling addressing as
    //   <instance>/finite-area[/<region>]/faMesh/faceLabels
    // (see faMesh::readFaceLabels / faMesh::writeObject); the reduced mesh
    // itself is rebuilt in-memory from the parent + faceLabels. So the
    // existence probe is a faceLabels probe, mirroring the layout used there.
    const word regionSub =
    (
        (areaName.empty() || areaName == polyMesh::defaultRegion)
      ? word::null : areaName
    );

    const fileName localDir =
    (
        regionSub.empty()
      ? (faMesh::prefix()/meshSubDir)
      : (faMesh::prefix()/regionSub/meshSubDir)
    );

    // READ_IF_PRESENT so a missing file yields a fallback instance rather than
    // a fatal error during the probe.
    const word instance = pMesh.time().findInstance
    (
        localDir,
        "faceLabels",
        IOobjectOption::READ_IF_PRESENT
    );

    IOobject meshIOobject
    (
        "faceLabels",
        instance,
        localDir,
        pMesh.thisDb(),
        IOobjectOption::READ_IF_PRESENT,
        IOobjectOption::NO_WRITE,
        IOobjectOption::NO_REGISTER
    );

    const fileName found
    (
        Foam::fileHandler().filePath(false, meshIOobject, "labelList")
    );

    bool looksValid = !found.empty();

    // Everybody needs it, or they all fail (matches legacy TryNew semantics).
    Pstream::reduceAnd(looksValid);

    return looksValid;
}


Foam::autoPtr<Foam::faMesh> Foam::faMesh::TryNew
(
    const word& areaName,
    const polyMesh& pMesh
)
{
    if (faMesh::hasMeshFiles(areaName, pMesh))
    {
        return autoPtr<faMesh>::New(areaName, pMesh);
    }

    return nullptr;
}


Foam::autoPtr<Foam::faMesh> Foam::faMesh::TryNew(const polyMesh& pMesh)
{
    return faMesh::TryNew(word::null, pMesh);
}


// ************************************************************************* //
