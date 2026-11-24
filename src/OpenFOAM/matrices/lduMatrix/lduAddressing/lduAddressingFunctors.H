/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Description
    functors for LDU-level offloading

\*---------------------------------------------------------------------------*/

#ifdef FOAM_OFFLOAD

#ifndef lduAddressingFunctors_H
#define lduAddressingFunctors_H

#include "lduAddressing.H"
#include "ops.H"

namespace Foam
{
    template<class Type, class OwnFun, class NeiFun, class OwnOp, class NeiOp>
    struct lduAddressingFunctor
    {
        // Store raw pointers (captured by value), not references.
        // References to UList/Field point to host memory that the GPU
        // cannot access — the GPU kernel must dereference these pointers
        // and they must be valid device addresses (or mapped/pinned).
        const label* __restrict__ ownStartPtr_;
        const label* __restrict__ neiStartPtr_;
        const label* __restrict__ losortPtr_;

        OwnFun ownFun_;
        NeiFun neiFun_;
        OwnOp ownOp_;
        NeiOp neiOp_;

        lduAddressingFunctor
        (
            const label* __restrict__ ownStart,
            const label* __restrict__ neiStart,
            const label* __restrict__ losort,
            const OwnFun ownFun,
            const NeiFun neiFun,
            const OwnOp ownOp,
            const NeiOp neiOp
        )
        :
            ownStartPtr_(ownStart),
            neiStartPtr_(neiStart),
            losortPtr_(losort),
            ownFun_(ownFun),
            neiFun_(neiFun),
            ownOp_(ownOp),
            neiOp_(neiOp)
        {}

        Type operator()(const label& id, const Type& d)
        {
            Type out = d;
            label oStart = ownStartPtr_[id];
            label oSize = ownStartPtr_[id+1] - oStart;

            label nStart = neiStartPtr_[id];
            label nSize = neiStartPtr_[id+1] - nStart;

            for (label i = 0; i < oSize; ++i)
            {
                label face = oStart + i;
                out = ownOp_(out, ownFun_(id, face));
            }

            for (label i = 0; i < nSize; ++i)
            {
                label face = losortPtr_[nStart + i];
                out = neiOp_(out, neiFun_(id, face));
            }
            return out;
        }
    };


    template<class Type, class LUType, class Op>
    struct matrixCoeffsMultiplyFunctor
    {
        // Store raw pointers (captured by value), not references.
        const Type* __restrict__ psiPtr_;
        const scalar* __restrict__ coeffsPtr_;
        const label* __restrict__ addrPtr_;
        Op op_;

        matrixCoeffsMultiplyFunctor
        (
            const Type* __restrict__ psi,
            const scalar* __restrict__ coeffs,
            const label* __restrict__ addr,
            const Op& op
        )
        :
            psiPtr_(psi),
            coeffsPtr_(coeffs),
            addrPtr_(addr),
            op_(op)
        {}

        Type operator()(const label& cell, const label& face)
        {
            return op_(coeffsPtr_[face]*psiPtr_[addrPtr_[face]]);
        }
    };


    template<class Type, class Op>
    struct matrixCoeffsFunctor
    {
        // Store raw pointer (captured by value), not a reference.
        const scalar* __restrict__ coeffsPtr_;
        Op op_;

        matrixCoeffsFunctor
        (
            const scalar* __restrict__ coeffs,
            const Op& op
        )
        :
            coeffsPtr_(coeffs),
            op_(op)
        {}

        Type operator()(const label& cell, const label& face)
        {
            return op_(coeffsPtr_[face]);
        }
    };


    template<class Type, class OwnFun, class NeiFun, class OwnOp,
        class NeiOp>
    inline void matrixOperation
    (
        const List<Type>& in,
        List<Type>& out,
        const lduAddressing& addr,
        OwnFun o,
        NeiFun n,
        OwnOp oo,
        NeiOp no
    )
    {
        const labelRange indices(out.size());

        // Capture raw pointers from the lduAddressing by value.
        // The lduAddressing object outlives this lambda (it's a reference
        // passed in), so the cdata() pointers remain valid.
        apply_policy
        (
        indices.size(),
        [
            in = in.cbegin(),
            ownStart = addr.ownerStartAddr().cdata(),
            neiStart = addr.losortStartAddr().cdata(),
            losort = addr.losortAddr().cdata(),
            out = out.begin(),
            &indices,
            o,
            n,
            oo,
            no
        ](auto policy)
        {
            std::transform
            (
                policy,
                indices.begin(),
                indices.end(),
                in,
                out,
                lduAddressingFunctor<Type, OwnFun, NeiFun, OwnOp, NeiOp>
                (
                    ownStart,
                    neiStart,
                    losort,
                    o,
                    n,
                    oo,
                    no
                )
            );
        });
    }


    template<class Type, class OwnFun, class NeiFun>
    inline void matrixOperation
    (
        const List<Type>& in,
        List<Type>& out,
        const lduAddressing& addr,
        OwnFun o,
        NeiFun n
    )
    {
        matrixOperation
        (
            in,
            out,
            addr,
            o,
            n,
            sumOp<Type>(),
            sumOp<Type>()
        );
    }
}

#endif

#endif  // FOAM_OFFLOAD

// ************************************************************************* //
