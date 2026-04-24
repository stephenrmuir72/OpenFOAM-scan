/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "fvcSurfaceIntegrate.H"
#include "fvMesh.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
template<class Type, bool integrate>
struct SurfaceIntegrateFunctor
{
    const Type zero_;

    const Field<Type>& issf_;
    const labelUList& ownStart_;
    const labelUList& neiStart_;
    const labelUList& own_;
    const labelUList& nei_;
    const labelUList& losort_;


    SurfaceIntegrateFunctor
    (
        const Field<Type>& issf,
        const labelUList& ownStart,
        const labelUList& neiStart,
        const labelUList& own,
        const labelUList& nei,
        const labelUList& losort
    )
    :
        zero_(pTraits<Type>::zero),
        issf_(issf),
        ownStart_(ownStart),
        neiStart_(neiStart),
        own_(own),
        nei_(nei),
        losort_(losort)
    {}


    Type operator()(const label celli) const
    {
        Type result(zero_);

        const label ownStart = ownStart_[celli];
        const label ownSize = ownStart_[celli + 1] - ownStart;

        for (label i = 0; i < ownSize; ++i)
        {
            const label facei = ownStart + i;
            result += issf_[facei];
        }

        const label neiStart = neiStart_[celli];
        const label neiSize = neiStart_[celli + 1] - neiStart;

        for (label i = 0; i < neiSize; ++i)
        {
            const label facei = losort_[neiStart + i];

            if (integrate)
            {
                result -= issf_[facei];
            }
            else
            {
                result += issf_[facei];
            }
        }

        return result;
    }
};


template<class Type>
struct surfaceIntegratePatchFunctor
{
    const Field<Type>& issf_;
    Field<Type>& igInterpolate_;

    const labelUList& neiStart_;
    const labelUList& losort_;
    const labelUList& patchCells_;

    surfaceIntegratePatchFunctor
    (
        const Field<Type>& issf,
        Field<Type>& igInterpolate,
        const labelUList& neiStart,
        const labelUList& losort,
        const labelUList& patchCells
    )
    :
        issf_(issf),
        igInterpolate_(igInterpolate),
        neiStart_(neiStart),
        losort_(losort),
        patchCells_(patchCells)
    {}


    bool operator()(const label id) const
    {
        const label celli = patchCells_[id];

        Type result(igInterpolate_[celli]);

        const label neiStart = neiStart_[id];
        const label neiSize = neiStart_[id + 1] - neiStart;

        for (label i = 0; i < neiSize; ++i)
        {
            const label facei = losort_[neiStart + i];
            result += issf_[facei];
        }

        igInterpolate_[celli] = result;

        return true;
    }
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void surfaceIntegrate
(
    Field<Type>& ivf,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const Field<Type>& issf = ssf;

    #ifdef FOAM_OFFLOAD
    labelRange indices(ivf.size());
    apply_policy(ivf.size(), [&](auto policy)
    {
        std::transform
        (   
            policy, 
            indices.begin(), 
            indices.end(),
            ivf.begin(),
            SurfaceIntegrateFunctor<Type,true>
            (   
                issf,
                mesh.lduAddr().ownerStartAddr(),
                mesh.lduAddr().losortStartAddr(),
                mesh.lduAddr().lowerAddr(),
                mesh.lduAddr().upperAddr(),
                mesh.lduAddr().losortAddr()
            )
        );
    });
    #else
    forAll(owner, facei)
    {
        ivf[owner[facei]] += issf[facei];
        ivf[neighbour[facei]] -= issf[facei];
    }
    #endif

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

        #ifdef FOAM_OFFLOAD
        const labelUList& pcells = mesh.lduAddr().patchSortCells(patchi);
        const labelUList& plosort = mesh.lduAddr().patchSortAddr(patchi);
        const labelUList& plosortStart = mesh.lduAddr().patchSortStartAddr(patchi);

        labelRange indices(pcells.size());
        apply_policy(pcells.size(), [&](auto policy)
        {
            std::for_each
            (
                policy,
                indices.begin(),
                indices.end(),
                surfaceIntegratePatchFunctor<Type>
                (   
                    pssf,
                    ivf,
                    plosortStart,
                    plosort,
                    pcells
                )
            );
        });
        #else
        forAll(mesh.boundary()[patchi], facei)
        {
            ivf[pFaceCells[facei]] += pssf[facei];
        }
        #endif
    }

    ivf /= mesh.Vsc();
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
surfaceIntegrate
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceIntegrate("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>(ssf.dimensions()/dimVol, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf.ref();

    surfaceIntegrate(vf.primitiveFieldRef(), ssf);
    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
surfaceIntegrate
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tssf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        fvc::surfaceIntegrate(tssf())
    );
    tssf.clear();
    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
surfaceSum
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceSum("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>(ssf.dimensions(), Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf.ref();

    #ifdef FOAM_OFFLOAD
    labelRange indices(vf.size());

    apply_policy(vf.size(), [&](auto policy)
    {
        std::transform
        (   
            policy, 
            indices.begin(), 
            indices.end(),
            vf.begin(),
            SurfaceIntegrateFunctor<Type,false>
            (   
                ssf,
                mesh.lduAddr().ownerStartAddr(),
                mesh.lduAddr().losortStartAddr(),
                mesh.lduAddr().lowerAddr(),
                mesh.lduAddr().upperAddr(),
                mesh.lduAddr().losortAddr()
            )
        );
    });
    #else
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        vf[owner[facei]] += ssf[facei];
        vf[neighbour[facei]] += ssf[facei];
    }
    #endif

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

        #ifdef FOAM_OFFLOAD
        const labelUList& pcells = mesh.lduAddr().patchSortCells(patchi);
        const labelUList& plosort = mesh.lduAddr().patchSortAddr(patchi);
        const labelUList& plosortStart = mesh.lduAddr().patchSortStartAddr(patchi);
        
        labelRange indices(pcells.size());

        apply_policy(pcells.size(),
        [&](auto policy)
        {
            std::for_each
            (
                policy,
                indices.begin(),
                indices.end(),
                surfaceIntegratePatchFunctor<Type>
                (   
                    pssf,
                    vf.primitiveFieldRef(),
                    plosortStart,
                    plosort,
                    pcells
                )
            );
        });
        #else
        forAll(mesh.boundary()[patchi], facei)
        {
            vf[pFaceCells[facei]] += pssf[facei];
        }
        #endif
    }

    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> surfaceSum
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tssf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf = surfaceSum(tssf());
    tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
