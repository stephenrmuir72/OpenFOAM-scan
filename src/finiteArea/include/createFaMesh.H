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
        aRegionName,
        runTime.timeName(),
        mesh,   //runTime,
        Foam::IOobject::MUST_READ
    ),
    false
);
aMesh.init(true);   // initialise all (lower levels and current)

// Note: the reduced-order (2D) mesh natively constructs Sf() with
// dimLength and V() with dimArea (see fvMesh::makeSf()/V()), so no
// dimension fix-up is required here.

// ************************************************************************* //
