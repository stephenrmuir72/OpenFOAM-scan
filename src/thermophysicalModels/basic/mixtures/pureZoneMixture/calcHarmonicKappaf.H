/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 OpenCFD Ltd.
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

Class
    Foam::calcHarmonicKappaf

Description
    Compute distance-weighted harmonic-averaged face-normal anisotropic
    thermal conductivity (kappa) for zone-based mixtures.

    For each internal face:
        kappaf = 1 / ((1-w)/k_own + w/k_nei)    // harmonic of n·(taniAlpha*cp)·n

    where w is the standard linear interpolation weight (|Cf - C_nei|/|d|),
    accounting for non-uniform cell sizes.

    At boundaries: directly project the adjacent cell tensor onto the face
    normal and multiply by cp.

Usage
    \code
    tmp<surfaceScalarField> tkappaf =
        calcHarmonicKappaf(mesh, taniAlpha, cp);
    \endcode

SourceFiles
    calcHarmonicKappaf.H (this file)

\*---------------------------------------------------------------------------*/

#ifndef calcHarmonicKappaf_H
#define calcHarmonicKappaf_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Function calcHarmonicKappaf
\*---------------------------------------------------------------------------*/

//- Compute harmonic-averaged face-normal kappa (kappaf) for anisotropic
//  zone mixtures. Uses distance-weighted harmonic mean via the standard
//  mesh interpolation weights (same as OpenFOAM's harmonic scheme).
inline tmp<surfaceScalarField> calcHarmonicKappaf
(
    const fvMesh& mesh,
    const volSymmTensorField& taniAlpha,
    const volScalarField& cp
)
{
    IOobject ioKappa
    (
        "kappaf",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    auto tKappa = taniAlpha*cp;
    const dimensionSet dimKappa = tKappa().dimensions();

    auto tkappaf = tmp<surfaceScalarField>::New
    (
        ioKappa,
        mesh,
        dimKappa
    );

    auto& kappaf = tkappaf.ref();

    // Standard linear interpolation weights from the mesh
    // w[facei] = |Cf - C_nei| / |C_nei - C_own|
    const surfaceScalarField& w = mesh.surfaceInterpolation::weights();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    // Internal faces: distance-weighted harmonic mean
    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        const vector n = Sf[facei] / magSf[facei];

        const scalar aOwn = n & taniAlpha[owner[facei]] & n;
        const scalar aNei = n & taniAlpha[neighbour[facei]] & n;
        const scalar kOwn = aOwn * cp[owner[facei]];
        const scalar kNei = aNei * cp[neighbour[facei]];

        const scalar wf = w[facei];  // weight toward neighbour

        kappaf[facei] = 1.0 / ((1.0 - wf)/kOwn + wf/kNei + VSMALL);
    }

    // Boundary faces: project adjacent cell tensor onto the face normal
    forAll(kappaf.boundaryField(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        const labelUList& faceCells = patch.faceCells();
        const vectorField n(patch.nf());

        auto& kappab = kappaf.boundaryFieldRef()[patchi];

        forAll(patch, facei)
        {
            const label celli = faceCells[facei];
            const scalar an = n[facei] & taniAlpha[celli] & n[facei];

            kappab[facei] = an * cp[celli];
        }
    }

    return tkappaf;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
