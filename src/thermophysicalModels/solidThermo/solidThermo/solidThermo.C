/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "solidThermo.H"
#include "fvMesh.H"
#include "fvm.H"
#include "fvc.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"


/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(solidThermo, 0);
    defineRunTimeSelectionTable(solidThermo, fvMesh);
    defineRunTimeSelectionTable(solidThermo, dictionary);
    defineRunTimeSelectionTable(solidThermo, fvMeshDictPhase);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace
{

//- Compute harmonic-averaged face-normal kappaf for anisotropic zone mixtures
//  kappaf = 1 / ((1-w)/k_own + w/k_nei)  where k = n·(taniAlpha·Cp)·n
tmp<surfaceScalarField> calcHarmonicKappaf
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

} // End anonymous namespace
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidThermo::solidThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    rhoThermo(mesh, phaseName)
{}


Foam::solidThermo::solidThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    rhoThermo(mesh, dict, phaseName)
{}


Foam::solidThermo::solidThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictionaryName
)
:
     rhoThermo(mesh, phaseName, dictionaryName)
{}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidThermo> Foam::solidThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<solidThermo>(mesh, phaseName);
}


Foam::autoPtr<Foam::solidThermo> Foam::solidThermo::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
{
    return basicThermo::New<solidThermo>(mesh, dict, phaseName);
}


Foam::autoPtr<Foam::solidThermo> Foam::solidThermo::New
(
     const fvMesh& mesh,
     const word& phaseName,
     const word& dictName
)
{
    return basicThermo::New<solidThermo>(mesh, phaseName, dictName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidThermo::~solidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::solidThermo::rho() const
{
    return rho_;
}


Foam::tmp<Foam::scalarField> Foam::solidThermo::rho(const label patchi) const
{
    return rho_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::solidThermo::rho()
{
    return rho_;
}


bool Foam::solidThermo::read()
{
    return regIOobject::read();
}


Foam::tmp<Foam::fvScalarMatrix> Foam::solidThermo::heatDiffusion
(
    const volScalarField& betav,
    const volScalarField& h
) const
{
    if (!isZoneMixture())
    {
        // Standard case: single implicit laplacian
        if (isotropic())
        {
            return fvm::laplacian(betav*alpha(), h, "laplacian(alpha,h)");
        }
        else
        {
            return fvm::laplacian(betav*(*Anialpha()), h, "laplacian(alpha,h)");
        }
    }
    else
    {
        // Zone mixture: add explicit correction terms.
        // At convergence, fvm::laplacian(alpha, h) and
        // fvc::laplacian(alpha, h) cancel, leaving
        // -fvc::laplacian(kappa, T), which is the correct diffusion.
        if (isotropic())
        {
            return
                fvm::laplacian(betav*alpha(), h, "laplacian(alpha,h)")
              - fvc::laplacian(betav*alpha(), h, "laplacian(alpha,h)")
              + fvc::laplacian(betav*kappa(), T(), "laplacian(kappa,T)");
        }
        else
        {
            const auto& aniAlpha = *Anialpha();

            const auto tkappaf =
                fvc::interpolate(betav)
              * calcHarmonicKappaf(h.mesh(), aniAlpha, this->Cp());

            return
                fvm::laplacian(betav*aniAlpha, h, "laplacian(alpha,h)")
              - fvc::laplacian(betav*aniAlpha, h, "laplacian(alpha,h)")
              + fvc::laplacian(tkappaf, T(), "laplacian(kappa,T)");
        }
    }
}


// ************************************************************************* //
