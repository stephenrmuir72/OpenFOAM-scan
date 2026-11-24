/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Keysight Technologies
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

Description
    Testing std::execution::par_unseq policy

\*---------------------------------------------------------------------------*/

#include "fvMatrix.H"
#include "UList.H"
#include "argList.H"
#include "Time.H"
#include "OSspecific.H"
#include "FieldOps.H"
#include "clockTime.H"
#include "Random.H"
#include "vectorField.H"
#include "fvCFD.H"

using namespace Foam;

// r/lib/x86_64-linux-gnu/libstdc++.so.6
// #7  0x00007fffd58abc6c in std::terminate() () from /usr/lib/x86_64-linux-gnu/libstdc++.so.6
// #8  0x00007fffd58c3901 in __cxa_throw () from /usr/lib/x86_64-linux-gnu/libstdc++.so.6
// #9  0x000055555560a720 in thrust::THRUST_200805_400500_NS::hip_rocprim::parallel_for<thrust::THRUST_200805_400500_NS::hip_rocprim::par_t, thrust::THRUST_200805_400500_NS::hip_rocprim::__fill::functor<char*, char>, long>(thrust::THRUST_200805_400500_NS::hip_rocprim::execution_policy<thrust::THRUST_200805_400500_NS::hip_rocprim::par_t>&, thrust::THRUST_200805_400500_NS::hip_rocprim::__fill::functor<char*, char>, long)::workaround::par(thrust::THRUST_200805_400500_NS::hip_rocprim::execution_policy<thrust::THRUST_200805_400500_NS::hip_rocprim::par_t>&, thrust::THRUST_200805_400500_NS::hip_rocprim::__fill::functor<char*, char>, long) ()
// #10 0x000055555560a48b in Foam::List<Foam::Vector<double> >::List(int, Foam::zero) ()
// #11 0x0000555555609655 in Foam::fvMatrix<Foam::Vector<double> >::fvMatrix(Foam::GeometricField<Foam::Vector<double>, Foam::fvPatchField, Foam::volMesh> const&, Foam::dimensionSet const&) ()
// #12 0x00007ffff7cb53c8 in Foam::fv::gaussConvectionScheme<Foam::Vector<double> >::fvmDiv(Foam::GeometricField<double, Foam::fvsPatchField, Foam::surfaceMesh> const&, Foam::GeometricField<Foam::Vector<double>, Foam::fvPatchField, Foam::volMesh> const&) const ()
// --Type <RET> for more, q to quit, c to continue without paging--
// using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >
>
gradf
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
    const word& name
)
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    const fvMesh& mesh = ssf.mesh();

    tmp<GradFieldType> tgGrad
    (
        new GradFieldType
        (
            IOobject
            (
                name,
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>(ssf.dimensions()/dimLength, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    GradFieldType& gGrad = tgGrad.ref();

    const vectorField& Sf = mesh.Sf();

    Field<GradType>& igGrad = gGrad;
    const Field<Type>& issf = ssf;

    #ifdef FOAM_OFFLOAD
    auto fun = [Sf = Sf.cdata(), issf = issf.cdata()]
    (const label, const label facei) -> GradType
    {
        return Sf[facei]*issf[facei];
    };

    auto ownOp = [](GradType& result, const GradType& value)
    {
        result += value;
    };

    auto neiOp = [](GradType& result, const GradType& value)
    {
        result -= value;
    };

    const labelRange cellRange(igGrad.size());

    apply_policy(igGrad.size(),
    [
        cellRange,
        igGrad = igGrad.begin(),
        ownStart = mesh.lduAddr().ownerStartAddr().cdata(),
        neiStart = mesh.lduAddr().losortStartAddr().cdata(),
        losort = mesh.lduAddr().losortAddr().cdata(),
        &fun,
        &ownOp,
        &neiOp
    ](auto policy)
    {
        std::transform
        (
            policy,
            cellRange.begin(),
            cellRange.end(),
            igGrad,   // d = old igGrad[celli]
            cellReduceFunctor
            <
                GradType,
                decltype(fun),
                decltype(fun),
                decltype(ownOp),
                decltype(neiOp)
            >
            (
                pTraits<GradType>::zero,
                ownStart,
                neiStart,
                losort,
                fun,
                fun,
                ownOp,
                neiOp
            )
        );
    });
    #else
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        const GradType Sfssf = Sf[facei]*issf[facei];

        igGrad[owner[facei]] += Sfssf;
        igGrad[neighbour[facei]] -= Sfssf;
    }
    #endif

    forAll(mesh.boundary(), patchi)
    {
        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];

        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

        #ifdef FOAM_OFFLOAD
        const labelUList& patchCells = mesh.lduAddr().patchSortCells(patchi);
        const labelUList& sort = mesh.lduAddr().patchSortAddr(patchi);
        const labelUList& sortStart = mesh.lduAddr().patchSortStartAddr(patchi);

        const labelRange cellRange(patchCells.size());

        auto neiFun =
        [pSf = pSf.cdata(), pssf = pssf.cdata()]
        (const label, const label facei) -> GradType
        {
            return pSf[facei]*pssf[facei];
        };

        auto neiOp = [](GradType& result, const GradType& value)
        {
            result += value;
        };

        apply_policy(cellRange.size(), [&](auto policy)
        {
            std::for_each
            (
                policy,
                cellRange.begin(),
                cellRange.end(),
                patchCellReduceFunctor
                <
                    GradType,
                    decltype(neiFun),
                    decltype(neiOp)
                >
                (
                    igGrad.data(),
                    sortStart.cdata(),
                    sort.cdata(),
                    patchCells.cdata(),
                    neiFun,
                    neiOp
                )
            );
        });
        #else
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        forAll(mesh.boundary()[patchi], facei)
        {
            igGrad[pFaceCells[facei]] += pSf[facei]*pssf[facei];
        }
        #endif
    }

    igGrad /= mesh.V();

    gGrad.correctBoundaryConditions();

    return tgGrad;
}

void residual
(
    const lduMatrix& m,
    solveScalarField& rA,
    const solveScalarField& psi,
    const scalarField& source,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
)
{
    solveScalar* __restrict__ rAPtr = rA.begin();

    const solveScalar* const __restrict__ psiPtr = psi.begin();
    const scalar* const __restrict__ diagPtr = m.diag().begin();
    const scalar* const __restrict__ sourcePtr = source.begin();

    const label* const __restrict__ uPtr = m.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = m.lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ upperPtr = m.upper().begin();
    const scalar* const __restrict__ lowerPtr = m.lower().begin();

    // Parallel boundary initialisation.
    // Note: there is a change of sign in the coupled
    // interface update.  The reason for this is that the
    // internal coefficients are all located at the l.h.s. of
    // the matrix whereas the "implicit" coefficients on the
    // coupled boundaries are all created as if the
    // coefficient contribution is of a source-kind (i.e. they
    // have a sign as if they are on the r.h.s. of the matrix.
    // To compensate for this, it is necessary to turn the
    // sign of the contribution.

    const label startRequest = UPstream::nRequests();

    // Initialise the update of interfaced interfaces
    m.initMatrixInterfaces
    (
        false,
        interfaceBouCoeffs,
        interfaces,
        psi,
        rA,
        cmpt
    );

    const label nCells = m.diag().size();
    const label nFaces = m.upper().size();

    #ifdef FOAM_OFFLOAD
    const labelRange cellIndices(nCells);
DebugVar(nCells);

    apply_policy
    (
        static_cast<std::size_t>(nCells),
        [
            =,
            ownerStartAddr = m.lduAddr().ownerStartAddr().cdata(),
            losortStartAddr = m.lduAddr().losortStartAddr().cdata(),
            upperAddr = m.lduAddr().upperAddr().cdata(),
            lowerAddr = m.lduAddr().lowerAddr().cdata(),
            upperCoeffs = m.upper().cdata(),
            lowerCoeffs = m.lower().cdata()
        ](auto policy)
    {
        std::for_each
        (
            policy,
            cellIndices.begin(),
            cellIndices.end(),
            [=](const label i)
            {
                // Start with source - diag*psi
                rAPtr[i] = sourcePtr[i] - diagPtr[i]*psiPtr[i];

                // Add contributions from owner faces (upper matrix)
                const label oStart = ownerStartAddr[i];
                const label oSize = ownerStartAddr[i+1] - oStart;
                for (label j = 0; j < oSize; ++j)
                {
                    const label facei = oStart + j;
                    rAPtr[upperAddr[facei]] -= upperCoeffs[facei]*psiPtr[i];
                }

                // Add contributions from losort faces (lower matrix)
                const label nStart = losortStartAddr[i];
                const label nSize = losortStartAddr[i+1] - nStart;
                for (label j = 0; j < nSize; ++j)
                {
                    const label facei = nStart + j;
                    rAPtr[lowerAddr[facei]] -= lowerCoeffs[facei]*psiPtr[i];
                }
            }
        );
    });
    #else
    for (label cell=0; cell<nCells; cell++)
    {
        rAPtr[cell] = sourcePtr[cell] - diagPtr[cell]*psiPtr[cell];
    }


    const label nFaces = upper().size();

    for (label face=0; face<nFaces; face++)
    {
        rAPtr[uPtr[face]] -= lowerPtr[face]*psiPtr[lPtr[face]];
        rAPtr[lPtr[face]] -= upperPtr[face]*psiPtr[uPtr[face]];
    }
    #endif

    // Update interface interfaces
    m.updateMatrixInterfaces
    (
        false,
        interfaceBouCoeffs,
        interfaces,
        psi,
        rA,
        cmpt,
        startRequest
    );
}



// Main program:
int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const label nCells = 10000000;
    const label nIntFaces = nCells-1;

    DebugVar(ListPolicy::is_aligned_type<scalar>());
    DebugVar(ListPolicy::is_aligned_type<vector>());
    DebugVar(ListPolicy::is_aligned_type<tensor>());
    DebugVar(is_contiguous_v<scalar>);
    DebugVar(is_contiguous_v<vector>);
    DebugVar(is_contiguous_v<tensor>);
    DebugVar(bool(std::is_constructible_v<scalar, int>));
    DebugVar(bool(std::is_constructible_v<vector, int>));
    DebugVar(bool(std::is_constructible_v<tensor, int>));

    DebugVar(bool(std::is_constructible_v<scalar, Foam::zero>));
    DebugVar(bool(std::is_constructible_v<vector, Foam::zero>));
    DebugVar(bool(std::is_constructible_v<tensor, Foam::zero>));

    // size_t len = 10;
    // for (label i = 0; i < 10; i++)
    // {
    //     Pout<< "Allocating fld size:" << len << endl;
    //     DebugVar(ListPolicy::use_alignment(len));
    //     DebugVar(ListPolicy::use_memory_pool(len));
    //     DebugVar(ListPolicy::use_offload(len));
    //     Field<vector> fld(len);
    //     Pout<< "Allocated fld size:" << fld.size()
    //         << " bytes:" << fld.size_bytes() << endl;
    //     fld = Zero;
    //     // fill_uniform<vector>(&fld, Zero);
    //     Pout<< "Done zeroing size:" << fld.size()
    //         << " bytes:" << fld.size_bytes() << endl;
    //     len *= 10;
    // }

{
    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
DebugVar("here");

    fvMatrix<scalar> fvm(fvm::laplacian(p));
DebugVar("here");

    solveScalarField rA(p.size());
    residual
    (
        fvm,
        rA,
        p.internalField(),
        fvm.source(),
        fvm.boundaryCoeffs(),
        p.boundaryField().scalarInterfaces(),
        0   // direction cmpt
    );

return 0;

    auto tgradp
    (
        gradf
        (
            linearInterpolate(p)(),
            "banana"
        )
    );
}
{
    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
DebugVar("here");

    fvMatrix<vector> fvm
    (
        U,
        dimless
    );
DebugVar("here");
    List<vector> junk(1368, Zero);
DebugVar("here");
}


    List<label> cellVals(nCells);

    Info<< "Testing fill of large list" << nl;
    {
        clockTime timing;
        for (label i = 0; i < 1000; ++i)
        {
            // Do some work to ensure we have a reasonable amount of time in the
            // test to allow for attaching a profiler, etc.
            std::fill_n(cellVals.begin(), cellVals.size(), 123);
        }
        Info<< "Sequential fill took " << timing.timeIncrement() << "s" << nl;
        for (label i = 0; i < 1000; ++i)
        {
            // Do some work to ensure we have a reasonable amount of time in the
            // test to allow for attaching a profiler, etc.
            std::fill_n
            (
                std::execution::par_unseq,
                cellVals.begin(),
                cellVals.size(),
                123
            );
        }
        Info<< "Parallel fill took " << timing.timeIncrement() << "s" << nl;
        for (label i = 0; i < 1000; ++i)
        {
            // Do some work to ensure we have a reasonable amount of time in the
            // test to allow for attaching a profiler, etc.
            std::fill_n(cellVals.begin(), cellVals.size(), 123);
        }
        Info<< "Sequential fill took " << timing.timeIncrement() << "s" << nl;
    }
    Info<< nl << endl;


    // Construct lduAddressing for all internal faces
    // to allow testing of interpolation with non-trivial addressing

    labelList P(nIntFaces);
    labelList N(P.size());
    forAll(P, facei)
    {
        P[facei] = facei;
        N[facei] = P[facei]+1;
    }
    // Interpolation weights for linear interpolation
    const scalarField lambda(P.size(), 0.5);
    // Faked face area vectors for testing
    const vectorField Sfi(P.size(), vector(1, 0, 0));



    Info<< "Testing dotInterpolate (cell to face) using lambda" << nl;
    {
        Field<scalar> sfi(P.size());
        Field<vector> vfi(nCells);
        Random::uniformGeneratorOp<vector> gen(123);
        FieldOps::assign(vfi, vfi, gen);

        const labelRange faceRange(P.size());

        auto dotInterpolate =
        [
            sfiPtr = sfi.data(),
            SfiPtr = Sfi.cdata(),
            lambdaPtr = lambda.cdata(),
            vfiPtr = vfi.cdata(),
            PPtr = P.cdata(),
            NPtr = N.cdata()
        ]
        (const label fi)
        {
            sfiPtr[fi] =
                SfiPtr[fi]
              & (lambdaPtr[fi]*(vfiPtr[PPtr[fi]] - vfiPtr[NPtr[fi]]) + vfiPtr[NPtr[fi]]);
        };

        clockTime timing;
        for (label i = 0; i < 1000; ++i)
        {
            std::for_each
            (
                faceRange.begin(),
                faceRange.end(),
                dotInterpolate
            );
        }
        Info<< "Sequential interpolation took " << timing.timeIncrement()
            << "s" << nl;
        for (label i = 0; i < 1000; ++i)
        {
            std::for_each
            (
                std::execution::par_unseq,
                faceRange.begin(),
                faceRange.end(),
                dotInterpolate
            );
        }
        Info<< "Parallel interpolation took " << timing.timeIncrement()
            << "s" << nl;
        for (label i = 0; i < 1000; ++i)
        {
            std::for_each
            (
                faceRange.begin(),
                faceRange.end(),
                dotInterpolate
            );
        }
        Info<< "Sequential interpolation took " << timing.timeIncrement()
            << "s" << nl;
        for (label i = 0; i < 1000; ++i)
        {
            std::for_each
            (
                std::execution::par_unseq,
                faceRange.begin(),
                faceRange.end(),
                dotInterpolate
            );
        }
        Info<< "Parallel interpolation took " << timing.timeIncrement()
            << "s" << nl;
    }

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
