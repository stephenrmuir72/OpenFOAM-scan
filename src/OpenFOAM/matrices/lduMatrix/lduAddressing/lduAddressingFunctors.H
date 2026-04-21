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
        const labelUList& ownStart;
        const labelUList& neiStart;
        const labelUList& losort;
        OwnFun oFun;
        NeiFun nFun;
        OwnOp oOp;
        NeiOp nOp;
    
        lduAddressingFunctor
        (
            const labelUList& _ownStart,
            const labelUList& _neiStart,
            const labelUList& _losort,
            const OwnFun _oFun,
            const NeiFun _nFun,
            const OwnOp _oOp,
            const NeiOp _nOp
        )
        :
            ownStart(_ownStart),
            neiStart(_neiStart),
            losort(_losort),
            oFun(_oFun),
            nFun(_nFun),
            oOp(_oOp),
            nOp(_nOp)
        {}
    
        Type operator()(const label& id,const Type& d)
        {
            Type out = d;
            label oStart = ownStart[id];
            label oSize = ownStart[id+1] - oStart;
    
            label nStart = neiStart[id];
            label nSize = neiStart[id+1] - nStart;
    
            for(label i = 0; i<oSize; i++)
            {
                label face = oStart + i;
                out = oOp(out, oFun(id, face));
            }
    
            for(label i = 0; i<nSize; i++)
            {
                label face = losort[nStart + i];
                out = nOp(out,nFun(id,face));
            }
            return out;
        }
    };
    
    template<class Type,class Fun,class Op>
    struct lduAddressingPatchFunctor
    {
        Field<Type>& field;
        const labelUList& neiStart;
        const labelUList& losort;
        const labelUList& pcells;
        Fun fun;
        Op op;
    
        lduAddressingPatchFunctor
        (
            Field<Type>& _field,
            const labelUList& _neiStart,
            const labelUList& _losort,
            const labelUList& _pcells,
            const Fun _fun,
            const Op _op
        )
        :
            field(_field),
            neiStart(_neiStart),
            losort(_losort),
            pcells(_pcells),
            fun(_fun),
            op(_op)
        {}
        void operator()(const label& id)
        {
            Type out = field[pcells[id]];
    
            label nStart = neiStart[id];
            label nSize = neiStart[id+1] - nStart;
            label cellI = pcells[id];
    
            for(label i = 0; i<nSize; i++)
            {
                label face = losort[nStart + i];
                out = op(out,fun(cellI,face));
            }

            field[pcells[id]] = out;
        }
    };
    
    template<class Type, class LUType, class Op>
    struct matrixCoeffsMultiplyFunctor
    {
        const Field<Type>& psi_;
        const scalarField& coeffs_;
        const labelUList& addr_;
        Op op_;

        matrixCoeffsMultiplyFunctor
        (
            const Field<Type>& psi,
            const scalarField& coeffs,
            const labelUList& addr,
            const Op& op
        )
        :
            psi_(psi),
            coeffs_(coeffs),
            addr_(addr),
            op_(op)
        {}

        Type operator()(const label& cell, const label& face)
        {
            return op_(coeffs_[face]*psi_[addr_[face]]);
        }
    };
    
    template<class Type, class Op>
    struct matrixCoeffsFunctor
    {
        const scalarField& coeffs_;
        Op op_;

        matrixCoeffsFunctor
        (
            const scalarField& coeffs,
            const Op& op
        )
        :
            coeffs_(coeffs),
            op_(op)
        {}

        Type operator()(const label& cell, const label& face)
        {
            return op_(coeffs_[face]);
        }
    };
    
    template<class Type,class Input, class OwnFun,class NeiFun,class OwnOp,class NeiOp>
    inline void matrixOperation
    (
        Input in,
        List<Type>& out,
        const lduAddressing& addr,
        OwnFun o,
        NeiFun n,
        OwnOp oo,
        NeiOp no
    )
    {
        const labelRange indices(out.size());
        std::transform
        (
            std::execution::par_unseq,
            indices.begin(),
            indices.end(),
            in,
            out.begin(),
            lduAddressingFunctor<Type,OwnFun,NeiFun,OwnOp,NeiOp>
            (
                addr.ownerStartAddr(),
                addr.losortStartAddr(),
                addr.losortAddr(),
                o,
                n,
                oo,
                no
            )
        );
    }
    
    template<class Type,class Input, class OwnFun,class NeiFun>
    inline void matrixOperation
    (
        Input in,
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
