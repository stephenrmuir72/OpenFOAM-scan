/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2026 UKRI-STFC.
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

    Some parts here added have been inspired by the work done at the following
    repository https://github.com/SimFlowCFD/RapidCFD-dev

\*---------------------------------------------------------------------------*/

#include "fvcReconstruct.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"

#ifdef FOAM_OFFLOAD
#include "lduAddressingFunctors.H"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef FOAM_OFFLOAD
struct reconstructMagFunctor
{
    const surfaceVectorField& Sf_;
    const surfaceVectorField& Cf_;
    const volVectorField& C_;
    const surfaceScalarField& ssf_;
    const surfaceScalarField& magSf_;


    reconstructMagFunctor
    (
        const surfaceVectorField& Sf,
        const surfaceVectorField& Cf,
        const volVectorField& C,
        const surfaceScalarField& ssf,
        const surfaceScalarField& magSf
    )
    :
        Sf_(Sf),
        Cf_(Cf),
        C_(C),
        ssf_(ssf),
        magSf_(magSf)
    {}


    scalar operator()(const label celli, const label facei) const
    {
        return
        (
            Sf_[facei] & (Cf_[facei] - C_[celli])
        )*ssf_[facei]/magSf_[facei];
    }
};
#endif

tmp<volScalarField> reconstructMag(const surfaceScalarField& ssf)
{
    const fvMesh& mesh = ssf.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    tmp<volScalarField> treconField
    (
        new volScalarField
        (
            IOobject
            (
                "reconstruct("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(ssf.dimensions()/dimArea, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    scalarField& rf = treconField.ref();
    #ifdef FOAM_OFFLOAD
    reconstructMagFunctor fun
    (
        Sf,
        Cf,
        C,
        ssf,
        magSf
    );

     // v is a random-access sequence of size n, all zeros
    std::vector<int> v(owner.size(), 0);
    matrixOperation
    (
        v.begin(),
        rf,
        mesh.lduAddr(),
        fun,
        fun,
        sumOp<scalar>(),
        minusOp<scalar>()
    );
    #else
    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        rf[own] += (Sf[facei] & (Cf[facei] - C[own]))*ssf[facei]/magSf[facei];
        rf[nei] -= (Sf[facei] & (Cf[facei] - C[nei]))*ssf[facei]/magSf[facei];
    }
    #endif

    const surfaceScalarField::Boundary& bsf = ssf.boundaryField();

    forAll(bsf, patchi)
    {
        const fvsPatchScalarField& psf = bsf[patchi];

        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];
        const vectorField& pSf = Sf.boundaryField()[patchi];
        const scalarField& pMagSf = magSf.boundaryField()[patchi];

        #ifdef FOAM_OFFLOAD
        const labelList& pcells = mesh.lduAddr().patchSortCells(patchi);

        const labelList& losort = mesh.lduAddr().patchSortAddr(patchi);
        const labelList& losortStart = 
            mesh.lduAddr().patchSortStartAddr(patchi);

        labelRange index(pcells.size());
        apply_policy(pcells.size(), [&](auto policy)
        {
            std::for_each
            (
                policy,
                index.begin(),
                index.end(),
                [&](size_t i)
                {
                    label celli = pcells[i];
                    scalar out = rf[celli];
                    label nStart = losortStart[i];
                    label nSize = losortStart[i+1] - nStart;
                    for (label j = 0; j < nSize; ++j)
                    {
                        label facei = losort[nStart + j];
                        out += (Sf[facei] & (Cf[facei] - C[celli]))
                            *ssf[facei]/magSf[facei];
                    }
                    rf[celli] = out;
                }
            );
        });
        #else
        forAll(pOwner, pFacei)
        {
            label own = pOwner[pFacei];
            rf[own] +=
                (pSf[pFacei] & (pCf[pFacei] - C[own]))
               *psf[pFacei]/pMagSf[pFacei];
        }
        #endif
    }

    rf /= mesh.V();

    treconField.ref().correctBoundaryConditions();

    return treconField;
}


tmp<volScalarField> reconstructMag(const tmp<surfaceScalarField>& tssf)
{
    tmp<volScalarField> tvf
    (
        fvc::reconstructMag(tssf())
    );
    tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
