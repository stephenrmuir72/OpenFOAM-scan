/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "basicFvGeometryScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicFvGeometryScheme, 0);
    addToRunTimeSelectionTable(fvGeometryScheme, basicFvGeometryScheme, dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicFvGeometryScheme::basicFvGeometryScheme
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvGeometryScheme(mesh, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::basicFvGeometryScheme::movePoints()
{
    fvGeometryScheme::movePoints();

    if (debug)
    {
        Pout<< "basicFvGeometryScheme::movePoints() : "
            << "recalculating primitiveMesh centres" << endl;
    }

    // Use lower level to calculate the geometry
    const_cast<fvMesh&>(mesh_).primitiveMesh::updateGeom();

    if (mesh_.nTopologicalD() == 2)
    {
        // Now we are parallel aware:
        // - calculate cell direction1 (= normal, since is now a face)
        // - average cell direction1 onto face
        // - re-orient face area to be perpendicular to direction1

        tmp<vectorField> tcDir0(mesh_.cellDir0());
        const auto& cDir0 = tcDir0();

        autoPtr<vectorField> tfDir1
        (
            autoPtr<vectorField>::New
            (
                mesh_.nFaces(),
                Zero
            )
        );
        auto& fDir1 = tfDir1.ref();

        const labelUList& owner = mesh_.faceOwner();
        const labelUList& neighbour = mesh_.neighbour();

        forAll(owner, facei)
        {
            fDir1[facei] = cDir0[owner[facei]];
        }

        vectorField nbrCellDir0;
        syncTools::swapBoundaryCellList(mesh_, cDir0, nbrCellDir0);

        // Finish internal faces
        forAll(neighbour, facei)
        {
            fDir1[facei] += cDir0[neighbour[facei]];
        }

        // Finish coupled faces
        for (auto& fpp : mesh_.boundary())
        {
            const auto* ppPtr = isA<coupledPolyPatch>(fpp.patch());
            if (ppPtr)
            {
                vectorField s
                (
                    fpp.boundarySlice(nbrCellDir0)
                  + fpp.patchInternalField(cDir0)
                );

                if (!ppPtr->owner())
                {
                    // Orientation of face is not reversed since only two
                    // vertices and 0th one should be matched
                    s = -s;
                }

                SubList<vector>(fDir1, fpp.size(), fpp.start()) = s;
            }
        }

        fDir1.normalise();

        // Re-do face-area to be perpendicular to the tangential face-direction
        vectorField& faceAreas = const_cast<vectorField&>(mesh_.faceAreas());
        const scalarField magS(mag(faceAreas));

        faceAreas = (faceAreas ^ fDir1);
        faceAreas.normalise();
        faceAreas *= magS;

        const_cast<fvMesh&>(mesh_).primitiveMesh::resetPrimitives
        (
            autoPtr<vectorField>(),
            autoPtr<vectorField>(),
            std::move(tfDir1),
            autoPtr<vectorField>()
        );
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::basicFvGeometryScheme::weights() const
{
    if (debug)
    {
        Pout<< "basicFvGeometryScheme::weights() : "
            << "Constructing weighting factors for face interpolation"
            << endl;
    }

    auto tweights =
        tmp<surfaceScalarField>::New
        (
            IOobject
            (
                "weights",
                mesh_.pointsInstance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimless
        );

    auto& weights = tweights.ref();
    weights.setOriented();

    // Set local references to mesh data
    // Note that we should not use fvMesh sliced fields at this point yet
    // since this causes a loop when generating weighting factors in
    // coupledFvPatchField evaluation phase
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& C = mesh_.cellCentres();
    const vectorField& Sf = mesh_.faceAreas();

    // ... and reference to the internal field of the weighting factors
    scalarField& w = weights.primitiveFieldRef();

    forAll(owner, facei)
    {
        // Note: mag in the dot-product.
        // For all valid meshes, the non-orthogonality will be less than
        // 90 deg and the dot-product will be positive.  For invalid
        // meshes (d & s <= 0), this will stabilise the calculation
        // but the result will be poor.
        scalar SfdOwn = mag(Sf[facei] & (Cf[facei] - C[owner[facei]]));
        scalar SfdNei = mag(Sf[facei] & (C[neighbour[facei]] - Cf[facei]));

        if (mag(SfdOwn + SfdNei) > ROOTVSMALL)
        {
            w[facei] = SfdNei/(SfdOwn + SfdNei);
        }
        else
        {
            w[facei] = 0.5;
        }
    }

    auto& wBf = weights.boundaryFieldRef();

    forAll(mesh_.boundary(), patchi)
    {
        mesh_.boundary()[patchi].makeWeights(wBf[patchi]);
    }

    if (debug)
    {
        Pout<< "basicFvGeometryScheme::weights : "
            << "Finished constructing weighting factors for face interpolation"
            << endl;
    }
    return tweights;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::basicFvGeometryScheme::deltaCoeffs() const
{
    if (debug)
    {
        Pout<< "basicFvGeometryScheme::deltaCoeffs() : "
            << "Constructing differencing factors array for face gradient"
            << endl;
    }

    // Force the construction of the weighting factors
    // needed to make sure deltaCoeffs are calculated for parallel runs.
    (void)mesh_.weights();

    auto tdeltaCoeffs =
        tmp<surfaceScalarField>::New
        (
            IOobject
            (
                "deltaCoeffs",
                mesh_.pointsInstance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimless/dimLength
        );

    auto& deltaCoeffs = tdeltaCoeffs.ref();
    deltaCoeffs.setOriented();


    // Set local references to mesh data
    const vectorField& C = mesh_.cellCentres();
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    forAll(owner, facei)
    {
        deltaCoeffs[facei] = 1.0/mag(C[neighbour[facei]] - C[owner[facei]]);
    }

    auto& deltaCoeffsBf = deltaCoeffs.boundaryFieldRef();

    forAll(deltaCoeffsBf, patchi)
    {
        const fvPatch& p = mesh_.boundary()[patchi];
        deltaCoeffsBf[patchi] = 1.0/mag(p.delta());

        // Optionally correct
        p.makeDeltaCoeffs(deltaCoeffsBf[patchi]);
    }

    return tdeltaCoeffs;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::basicFvGeometryScheme::nonOrthDeltaCoeffs() const
{
    if (debug)
    {
        Pout<< "basicFvGeometryScheme::nonOrthDeltaCoeffs() : "
            << "Constructing differencing factors array for face gradient"
            << endl;
    }

    // Force the construction of the weighting factors
    // needed to make sure deltaCoeffs are calculated for parallel runs.
    weights();

    auto tnonOrthDeltaCoeffs =
        tmp<surfaceScalarField>::New
        (
            IOobject
            (
                "nonOrthDeltaCoeffs",
                mesh_.pointsInstance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimless/dimLength
        );

    auto& nonOrthDeltaCoeffs = tnonOrthDeltaCoeffs.ref();
    nonOrthDeltaCoeffs.setOriented();


    // Set local references to mesh data
    const volVectorField& C = mesh_.C();
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    forAll(owner, facei)
    {
        vector delta = C[neighbour[facei]] - C[owner[facei]];
        vector unitArea = Sf[facei]/magSf[facei];

        // Standard cell-centre distance form
        //NonOrthDeltaCoeffs[facei] = (unitArea & delta)/magSqr(delta);

        // Slightly under-relaxed form
        //NonOrthDeltaCoeffs[facei] = 1.0/mag(delta);

        // More under-relaxed form
        //NonOrthDeltaCoeffs[facei] = 1.0/(mag(unitArea & delta) + VSMALL);

        // Stabilised form for bad meshes
        nonOrthDeltaCoeffs[facei] = 1.0/max(unitArea & delta, 0.05*mag(delta));
    }

    auto& nonOrthDeltaCoeffsBf = nonOrthDeltaCoeffs.boundaryFieldRef();

    forAll(nonOrthDeltaCoeffsBf, patchi)
    {
        fvsPatchScalarField& patchDeltaCoeffs = nonOrthDeltaCoeffsBf[patchi];

        const fvPatch& p = patchDeltaCoeffs.patch();

        const vectorField patchDeltas(mesh_.boundary()[patchi].delta());

        forAll(p, patchFacei)
        {
            vector unitArea =
                Sf.boundaryField()[patchi][patchFacei]
               /magSf.boundaryField()[patchi][patchFacei];

            const vector& delta = patchDeltas[patchFacei];

            patchDeltaCoeffs[patchFacei] =
                1.0/max(unitArea & delta, 0.05*mag(delta));
        }

        // Optionally correct
        p.makeNonOrthoDeltaCoeffs(patchDeltaCoeffs);
    }
    return tnonOrthDeltaCoeffs;
}


Foam::tmp<Foam::surfaceVectorField>
Foam::basicFvGeometryScheme::nonOrthCorrectionVectors() const
{
    if (debug)
    {
        Pout<< "surfaceInterpolation::makeNonOrthCorrectionVectors() : "
            << "Constructing non-orthogonal correction vectors"
            << endl;
    }

    auto tnonOrthCorrectionVectors =
        tmp<surfaceVectorField>::New
        (
            IOobject
            (
                "nonOrthCorrectionVectors",
                mesh_.pointsInstance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimless
        );

    auto& corrVecs = tnonOrthCorrectionVectors.ref();
    corrVecs.setOriented();

    // Set local references to mesh data
    const volVectorField& C = mesh_.C();
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    tmp<surfaceScalarField> tNonOrthDeltaCoeffs(nonOrthDeltaCoeffs());
    const surfaceScalarField& NonOrthDeltaCoeffs = tNonOrthDeltaCoeffs();

    forAll(owner, facei)
    {
        vector unitArea(Sf[facei]/magSf[facei]);
        vector delta(C[neighbour[facei]] - C[owner[facei]]);

        corrVecs[facei] = unitArea - delta*NonOrthDeltaCoeffs[facei];
    }

    // Boundary correction vectors set to zero for boundary patches
    // and calculated consistently with internal corrections for
    // coupled patches

    auto& corrVecsBf = corrVecs.boundaryFieldRef();

    forAll(corrVecsBf, patchi)
    {
        fvsPatchVectorField& patchCorrVecs = corrVecsBf[patchi];

        const fvPatch& p = patchCorrVecs.patch();

        if (!patchCorrVecs.coupled())
        {
            patchCorrVecs = Zero;
        }
        else
        {
            const auto& patchNonOrthDeltaCoeffs =
                NonOrthDeltaCoeffs.boundaryField()[patchi];

            const vectorField patchDeltas(mesh_.boundary()[patchi].delta());

            forAll(p, patchFacei)
            {
                vector unitArea =
                    Sf.boundaryField()[patchi][patchFacei]
                   /magSf.boundaryField()[patchi][patchFacei];

                const vector& delta = patchDeltas[patchFacei];

                patchCorrVecs[patchFacei] =
                    unitArea - delta*patchNonOrthDeltaCoeffs[patchFacei];
            }
        }

        // Optionally correct
        p.makeNonOrthoCorrVectors(patchCorrVecs);
    }

    if (debug)
    {
        Pout<< "surfaceInterpolation::makeNonOrthCorrectionVectors() : "
            << "Finished constructing non-orthogonal correction vectors"
            << endl;
    }
    return tnonOrthCorrectionVectors;
}


// ************************************************************************* //
