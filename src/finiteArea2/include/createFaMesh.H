/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 Wikki Ltd
    Copyright (C) 2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Create a finite-area mesh for the defaultRegion

Required Variables
    - mesh [polyMesh]

Provided Variables
    - aMesh [faMesh]

\*---------------------------------------------------------------------------*/

Foam::Info << "Create finite-area mesh";
Foam::Info << Foam::nl;

//// "getRegionOption.H"
//Foam::word aRegionName
//(
//    args.getOrDefault<word>("region", Foam::polyMesh::defaultRegion)
//);

Foam::faMesh aMesh
(
    Foam::IOobject
    (
        regionName,
        runTime.timeName(),
        mesh,   //runTime,
        Foam::IOobject::MUST_READ
    ),
    false
);
aMesh.init(true);   // initialise all (lower levels and current)

// Make sure all 'fv' properties are of the correct dimensions
const_cast<surfaceVectorField&>(aMesh.Sf()).dimensions() = dimLength;
const_cast<DimensionedField<scalar, volMesh>&>(aMesh.V()).dimensions() =
    dimArea;

// ************************************************************************* //
