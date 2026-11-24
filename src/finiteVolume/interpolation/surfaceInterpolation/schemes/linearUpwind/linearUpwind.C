/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "linearUpwind.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::linearUpwind<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "linearUpwind::correction(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensioned<Type>(vf.name(), vf.dimensions(), Zero)
        )
    );

    GeometricField<Type, fvsPatchField, surfaceMesh>& sfCorr = tsfCorr.ref();

    const surfaceScalarField& faceFlux = this->faceFlux_;

    const auto& owner = mesh.owner();
    const auto& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    tmp<fv::gradScheme<scalar>> gradScheme_
    (
        fv::gradScheme<scalar>::New
        (
            mesh,
            mesh.gradScheme(gradSchemeName_)
        )
    );

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        tmp<volVectorField> tgradVf =
            gradScheme_().grad(vf.component(cmpt), gradSchemeName_);

        const volVectorField& gradVf = tgradVf();

        labelRange faceRange(faceFlux.size());

        // Capture raw pointers by value — the GPU kernel must access them
        // and references to UList/Field point to host memory.
        const scalar* __restrict__ faceFluxPtr = faceFlux.cdata();
        const label* __restrict__ ownerPtr = owner.cdata();
        const label* __restrict__ neighbourPtr = neighbour.cdata();
        const vector* __restrict__ CfPtr = Cf.cdata();
        const vector* __restrict__ CPtr = C.cdata();
        const vector* __restrict__ gradVfPtr = gradVf.cdata();
        Type* __restrict__ sfCorrPtr = sfCorr.data();

        apply_policy(faceFlux.size(),
        [
            &,
            faceFluxPtr,
            ownerPtr,
            neighbourPtr,
            CfPtr,
            CPtr,
            gradVfPtr,
            cmpt,
            sfCorrPtr
        ](auto policy)
        {
            std::for_each
            (
                policy,
                faceRange.begin(),
                faceRange.end(),
                [faceFluxPtr, ownerPtr, neighbourPtr, CfPtr, CPtr, gradVfPtr, cmpt, sfCorrPtr]
                (size_t facei)
                {
                    const label celli =
                        (faceFluxPtr[facei] > 0)
                      ? ownerPtr[facei] : neighbourPtr[facei];

                    setComponent(sfCorrPtr[facei], cmpt) =
                        (CfPtr[facei] - CPtr[celli]) & gradVfPtr[celli];
                }
            );
        });

        typename GeometricField<Type, fvsPatchField, surfaceMesh>::
            Boundary& bSfCorr = sfCorr.boundaryFieldRef();

        forAll(bSfCorr, patchi)
        {
            fvsPatchField<Type>& pSfCorr = bSfCorr[patchi];

            if (pSfCorr.coupled())
            {
                const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
                const vectorField& pCf = Cf.boundaryField()[patchi];
                const scalarField& pFaceFlux = faceFlux.boundaryField()[patchi];

                const vectorField pGradVfNei
                (
                    gradVf.boundaryField()[patchi].patchNeighbourField()
                );

                // Build the d-vectors
                const vectorField pd
                (
                    Cf.boundaryField()[patchi].patch().delta()
                );

                const labelRange faceRange(pOwner.size());

                // Capture raw pointers by value
                const label* __restrict__ pOwnerPtr = pOwner.cdata();
                const vector* __restrict__ pCfPtr = pCf.cdata();
                const scalar* __restrict__ pFaceFluxPtr = pFaceFlux.cdata();
                const vector* __restrict__ pdPtr = pd.cdata();
                const vector* __restrict__ pGradVfNeiPtr = pGradVfNei.cdata();
                const vector* __restrict__ CPtr = C.cdata();
                const vector* __restrict__ gradVfPtr = gradVf.cdata();
                Type* __restrict__ pSfCorrPtr = pSfCorr.data();

                apply_policy(pOwner.size(),
                [
                    &,
                    pOwnerPtr,
                    pCfPtr,
                    pFaceFluxPtr,
                    pdPtr,
                    pGradVfNeiPtr,
                    CPtr,
                    gradVfPtr,
                    cmpt,
                    pSfCorrPtr
                ](auto policy)
                {
                    std::for_each
                    (
                        policy,
                        faceRange.begin(),
                        faceRange.end(),
                        [pOwnerPtr, pCfPtr, pFaceFluxPtr, pdPtr,
                         pGradVfNeiPtr, CPtr, gradVfPtr, cmpt, pSfCorrPtr]
                        (size_t facei)
                        {
                            label own = pOwnerPtr[facei];

                            if (pFaceFluxPtr[facei] > 0)
                            {
                                setComponent(pSfCorrPtr[facei], cmpt) =
                                    (pCfPtr[facei] - CPtr[own])
                                  & gradVfPtr[own];
                            }
                            else
                            {
                                setComponent(pSfCorrPtr[facei], cmpt) =
                                    (pCfPtr[facei] - pdPtr[facei] - CPtr[own])
                                  & pGradVfNeiPtr[facei];
                            }
                        }
                    );
                });
            }
        }
    }

    return tsfCorr;
}


template<>
Foam::tmp<Foam::surfaceVectorField>
Foam::linearUpwind<Foam::vector>::correction
(
    const volVectorField& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceVectorField> tsfCorr
    (
        new surfaceVectorField
        (
            IOobject
            (
                "linearUpwind::correction(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensioned<vector>(vf.name(), vf.dimensions(), Zero)
        )
    );

    surfaceVectorField& sfCorr = tsfCorr.ref();

    const surfaceScalarField& faceFlux = this->faceFlux_;

    const auto& owner = mesh.owner();
    const auto& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    tmp<fv::gradScheme<vector>> gradScheme_
    (
        fv::gradScheme<vector>::New
        (
            mesh,
            mesh.gradScheme(gradSchemeName_)
        )
    );

    tmp<volTensorField> tgradVf = gradScheme_().grad(vf, gradSchemeName_);
    const volTensorField& gradVf = tgradVf();

    labelRange faceRange(faceFlux.size());

    // Capture raw pointers by value
    const scalar* __restrict__ faceFluxPtr = faceFlux.cdata();
    const label* __restrict__ ownerPtr = owner.cdata();
    const label* __restrict__ neighbourPtr = neighbour.cdata();
    const vector* __restrict__ CfPtr = Cf.cdata();
    const vector* __restrict__ CPtr = C.cdata();
    const tensor* __restrict__ gradVfPtr = gradVf.cdata();
    vector* __restrict__ sfCorrPtr = sfCorr.data();

    apply_policy(faceFlux.size(),
    [
        faceFluxPtr,
        ownerPtr,
        neighbourPtr,
        CfPtr,
        CPtr,
        gradVfPtr,
        sfCorrPtr,
        faceRange
    ](auto policy)
    {
        std::for_each
        (
            policy,
            faceRange.begin(),
            faceRange.end(),
            [faceFluxPtr, ownerPtr, neighbourPtr, CfPtr, CPtr, gradVfPtr, sfCorrPtr]
            (size_t facei)
            {
                const label celli =
                    (faceFluxPtr[facei] > 0)
                  ? ownerPtr[facei] : neighbourPtr[facei];
                sfCorrPtr[facei] = (CfPtr[facei] - CPtr[celli]) & gradVfPtr[celli];
            }
        );
    });


    typename surfaceVectorField::Boundary& bSfCorr = sfCorr.boundaryFieldRef();

    forAll(bSfCorr, patchi)
    {
        fvsPatchVectorField& pSfCorr = bSfCorr[patchi];

        if (pSfCorr.coupled())
        {
            const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
            const vectorField& pCf = Cf.boundaryField()[patchi];
            const scalarField& pFaceFlux = faceFlux.boundaryField()[patchi];

            const tensorField pGradVfNei
            (
                gradVf.boundaryField()[patchi].patchNeighbourField()
            );

            // Build the d-vectors
            vectorField pd(Cf.boundaryField()[patchi].patch().delta());

            const labelRange faceRange(pOwner.size());

            // Capture raw pointers by value
            const label* __restrict__ pOwnerPtr = pOwner.cdata();
            const vector* __restrict__ pCfPtr = pCf.cdata();
            const scalar* __restrict__ pFaceFluxPtr = pFaceFlux.cdata();
            const tensor* __restrict__ pGradVfNeiPtr = pGradVfNei.cdata();
            const vector* __restrict__ pdPtr = pd.cdata();
            const vector* __restrict__ CPtr = C.cdata();
            const tensor* __restrict__ gradVfPtr = gradVf.cdata();
            vector* __restrict__ pSfCorrPtr = pSfCorr.data();

            apply_policy(pOwner.size(),
            [
                pOwnerPtr,
                pCfPtr,
                pFaceFluxPtr,
                pGradVfNeiPtr,
                pdPtr,
                CPtr,
                gradVfPtr,
                pSfCorrPtr,
                faceRange
            ](auto policy)
            {
                std::for_each
                (
                    policy,
                    faceRange.begin(),
                    faceRange.end(),
                    [pOwnerPtr, pCfPtr, pFaceFluxPtr, pGradVfNeiPtr,
                     pdPtr, CPtr, gradVfPtr, pSfCorrPtr]
                    (size_t facei)
                    {
                        label own = pOwnerPtr[facei];

                        if (pFaceFluxPtr[facei] > 0)
                        {
                            pSfCorrPtr[facei] =
                                (pCfPtr[facei] - CPtr[own])
                              & gradVfPtr[own];
                        }
                        else
                        {
                            pSfCorrPtr[facei] =
                                (pCfPtr[facei] - pdPtr[facei] - CPtr[own])
                              & pGradVfNeiPtr[facei];
                        }
                    }
                );
            });
        }
    }

    return tsfCorr;
}


namespace Foam
{
    makelimitedSurfaceInterpolationScheme(linearUpwind)
}

// ************************************************************************* //
