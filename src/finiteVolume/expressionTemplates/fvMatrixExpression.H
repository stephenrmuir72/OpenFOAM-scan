/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025-2026 M. Janssens
    Copyright (C) 2025-2026 OpenCFD Ltd.
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

InNamespace
    Foam::Expression

Description
    Expression templates for fvMatrix.

    - diag,upper,lower are scalar expressions
    - boundaryCoeffs,internalCoeffs will be of the underlying field type
    - faceFluxCorrection will be of the face flux field type (underlying field
      type)

    TBD:
    - useImplicit_ etc.
    - lduMatrix::lowerCSR ?

SourceFiles
    fvMatrixExpression.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvMatrixExpression_H
#define Foam_fvMatrixExpression_H

#include "ListExpression.H"
#include "GeometricFieldExpression.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "fvm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type> class fvMatrix;

namespace Expression
{

/*---------------------------------------------------------------------------*\
                     Class fvMatrixExpression Declaration
\*---------------------------------------------------------------------------*/

//- Expression of fvMatrix
//
// Bit different:
// - immediately stores state (hasDiag, hasUpper etc.)
// - delays execution of building coefficients
// - big problem: individual components can be built differently so
//   different template argument ... Needs to be templated on all outputs
//   separately!
template
<
    typename E,
    typename DiagExpr,
    typename UpperExpr,
    typename LowerExpr,
    typename BoundaryExpr,
    typename FaceFluxExpr,
    typename SourceExpr
>
class fvMatrixExpression
{
protected:

    const bool hasDiag_;
    const bool hasUpper_;
    const bool hasLower_;
    const bool hasFaceFluxCorrection_;
    dimensionSet dimensions_;

public:
    static constexpr bool is_leaf = false;

    // Constructors

        fvMatrixExpression
        (
            const bool hasDiag,
            const bool hasUpper,
            const bool hasLower,
            const bool hasFaceFluxCorrection,
            const dimensionSet& dimensions
        )
        :
            hasDiag_(hasDiag),
            hasUpper_(hasUpper),
            hasLower_(hasLower),
            hasFaceFluxCorrection_(hasFaceFluxCorrection),
            dimensions_(dimensions)
        {
            #ifdef DEBUG_ET
            Pout<< "fvMatrixExpression :" << nl
                << "    E              : " << error::demangle<E>() << nl;
            if (hasDiag_)
            {
                Pout<< "    DiagExpr       : "
                    << error::demangle<DiagExpr>() << nl;
            }
            {
                Pout<< "    SourceExpr     : "
                    << error::demangle<SourceExpr>() << nl;
            }
            if (hasUpper_)
            {
                Pout<< "    UpperExpr      : "
                    << error::demangle<UpperExpr>() << nl;
            }
            if (hasLower_)
            {
                Pout<< "    LowerExpr      : "
                    << error::demangle<LowerExpr>() << nl;
            }
            {
                Pout<< "    BoundaryExpr   : "
                    << error::demangle<BoundaryExpr>() << nl;
            }
            if (hasFaceFluxCorrection_)
            {
                Pout<< "    FaceFluxExpr   : "
                    << error::demangle<FaceFluxExpr>() << nl;
            }
            #endif
        }


    // lduMatrix part

        DiagExpr diag() const
        {
            return static_cast<const E&>(*this).diag();
        }

        DiagExpr diag()
        {
            return static_cast<E&>(*this).diag();
        }

        UpperExpr upper() const
        {
            return static_cast<const E&>(*this).upper();
        }

        UpperExpr upper()
        {
            return static_cast<E&>(*this).upper();
        }

        LowerExpr lower() const
        {
            return static_cast<const E&>(*this).lower();
        }

        LowerExpr lower()
        {
            return static_cast<E&>(*this).lower();
        }

        bool hasDiag() const
        {
             return hasDiag_;
        }
        bool hasUpper() const
        {
             return hasUpper_;
        }
        bool hasLower() const
        {
             return hasLower_;
        }
        bool diagonal() const
        {
             return (hasDiag() && !hasLower() && !hasUpper());
        }
        bool symmetric() const
        {
             return (hasDiag() && !hasLower() && hasUpper());
        }
        bool asymmetric() const
        {
             return (hasDiag() && hasLower() && hasUpper());
        }
        bool hasFaceFluxCorrection() const
        {
             return hasFaceFluxCorrection_;
        }


    // fvMatrix

        SourceExpr source()
        {
            return static_cast<E&>(*this).source();
        }

        SourceExpr source() const
        {
            return static_cast<const E&>(*this).source();
        }

        BoundaryExpr internalCoeffs(const label i)
        {
            return static_cast<E&>(*this).internalCoeffs(i);
        }

        BoundaryExpr internalCoeffs(const label i) const
        {
            return static_cast<const E&>(*this).internalCoeffs(i);
        }

        BoundaryExpr boundaryCoeffs(const label i)
        {
            return static_cast<E&>(*this).boundaryCoeffs(i);
        }

        BoundaryExpr boundaryCoeffs(const label i) const
        {
            return static_cast<const E&>(*this).boundaryCoeffs(i);
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return static_cast<const E&>(*this).faceFluxCorrection();
        }

        FaceFluxExpr faceFluxCorrection()
        {
            return static_cast<E&>(*this).faceFluxCorrection();
        }

        const auto& dimensions() const noexcept
        {
            return dimensions_;
        }

        auto& dimensions() noexcept
        {
            return dimensions_;
        }


    // Other

        //- Helper to evaluate (=construct) an fvMatrix
        template<class Matrix>
        Matrix& evaluate(Matrix& m) const
        {
            // Do any direct evaluation into the matrix
            static_cast<const E&>(*this).evaluate(m, scalar(1));
            {
                // Add to source.
                (source() + m.source().expr()).evaluate(m.source());
            }
            if (hasDiag())
            {
                if (m.hasDiag())
                {
                    (diag() + m.diag().expr()).evaluate(m.diag());
                }
                else
                {
                    diag().evaluate(m.diag());
                }
            }
            if (hasUpper())
            {
                if (m.hasUpper())
                {
                    (upper() + m.upper().expr()).evaluate(m.upper());
                }
                else
                {
                    upper().evaluate(m.upper());
                }
            }
            if (hasLower())
            {
                if (m.hasLower())
                {
                    (lower() + m.lower().expr()).evaluate(m.lower());
                }
                else
                {
                    lower().evaluate(m.lower());
                }
            }
            if (hasFaceFluxCorrection())
            {
                if (m.hasFaceFluxCorrection())
                {
                    auto& mffCorr = *m.faceFluxCorrectionPtr();
                    (faceFluxCorrection() + mffCorr.expr()).evaluate(mffCorr);
                }
                else
                {
                    faceFluxCorrection().evaluate(*m.faceFluxCorrectionPtr());
                }
            }

            // Do boundary
            if (hasUpper() || hasLower())
            {
                auto& intCoeffs = m.internalCoeffs();
                auto& bouCoeffs = m.boundaryCoeffs();
                const label n = intCoeffs.size();

                typedef
                    typename std::remove_reference_t<decltype(m.source())>::
                    value_type value_type;
                typedef Field<value_type> FieldType;

#if defined(__STDPAR__)
                constexpr bool foam_offload = true;
#else
                constexpr bool foam_offload = false;
#endif
                if constexpr (foam_offload)
                {
                    //- Parallel via ListsRefWrap. See GeometricFieldExpression.
                    //- Different behaviour for if adding v.s. writing

                    DynamicList<label> intID;
                    intID.reserve(n);
                    bool adding = false;
                    for (label i = 0; i < n; ++i)
                    {
                        intID.emplace_back(i);
                        if (intCoeffs.set(i) || bouCoeffs.set(i))
                        {
                            adding = true;
                        }
                        if (!intCoeffs.set(i))
                        {
                            const auto intExpr = internalCoeffs(i);
                            intCoeffs.set(i, new FieldType(intExpr.size()));
                        }
                        if (!bouCoeffs.set(i))
                        {
                            const auto bouExpr = boundaryCoeffs(i);
                            bouCoeffs.set(i, new FieldType(bouExpr.size()));
                        }
                    }

                    // Internal : register 'chunks' for destination container
                    ListsRefWrap<FieldType> intExpr
                    (
                        intCoeffs,      // destination
                        static_cast<const labelUList&>(intID)
                    );
                    // Boundary : register 'chunks' for destination container
                    ListsRefWrap<FieldType> bouExpr
                    (
                        bouCoeffs,      // destination
                        static_cast<const labelUList&>(intID)
                    );

                    if (adding)
                    {
                        // Create new expression for adding components.

                        typedef List_add
                        <
                            BoundaryExpr,
                            ListConstRefWrap<typename BoundaryExpr::value_type>
                        > NewBoundaryExpr;

                        // Split into internal and boundary expressions
                        std::vector<NewBoundaryExpr> intExprs;
                        std::vector<NewBoundaryExpr> bouExprs;

                        intExprs.reserve(n);
                        bouExprs.reserve(n);
                        for (label i = 0; i < n; ++i)
                        {
                            const auto intExpr = internalCoeffs(i);
                            const auto intExpr2 = intExpr + intCoeffs[i].expr();
                            intExprs.emplace_back(intExpr2);

                            const auto bouExpr = boundaryCoeffs(i);
                            const auto bouExpr2 = bouExpr + bouCoeffs[i].expr();
                            bouExprs.emplace_back(bouExpr2);
                        }

                        // Internal : evaluate all at once
                        ListsRefWrap<NewBoundaryExpr>(intExprs).evaluate(intExpr);

                        // Boundary : evaluate all at once
                        ListsRefWrap<NewBoundaryExpr>(bouExprs).evaluate(bouExpr);
                    }
                    else
                    {
                        // Split into internal and boundary expressions
                        std::vector<BoundaryExpr> intExprs;
                        std::vector<BoundaryExpr> bouExprs;

                        intExprs.reserve(n);
                        bouExprs.reserve(n);
                        for (label i = 0; i < n; ++i)
                        {
                            const auto intExpr = internalCoeffs(i);
                            intExprs.emplace_back(internalCoeffs(i));
                            const auto bouExpr = boundaryCoeffs(i);
                            bouExprs.emplace_back(boundaryCoeffs(i));
                        }

                        // Internal : evaluate all at once
                        ListsRefWrap<BoundaryExpr>(intExprs).evaluate(intExpr);

                        // Boundary : evaluate all at once
                        ListsRefWrap<BoundaryExpr>(bouExprs).evaluate(bouExpr);
                    }
                }
                else
                {
                    //- Sequential
                    for (label i = 0; i < n; ++i)
                    {
                        const auto intExpr = internalCoeffs(i);
                        if (intCoeffs.set(i))
                        {
                            // Add
                            (intExpr + intCoeffs[i].expr()).
                                evaluate(intCoeffs[i]);
                        }
                        else
                        {
                            intCoeffs.set(i, new FieldType(intExpr.size()));
                            intExpr.evaluate(intCoeffs[i]);
                        }
                        const auto bouExpr = boundaryCoeffs(i);
                        if (bouCoeffs.set(i))
                        {
                            // Add
                            (bouExpr + bouCoeffs[i].expr()).
                                evaluate(bouCoeffs[i]);
                        }
                        else
                        {
                            bouCoeffs.set(i, new FieldType(bouExpr.size()));
                            bouExpr.evaluate(bouCoeffs[i]);
                        }
                    }
                }
            }
            const_cast<dimensionSet&>(m.dimensions()) = dimensions();

            return m;
        }
};


/*---------------------------------------------------------------------------*\
                          Class Default Declaration
\*---------------------------------------------------------------------------*/

//- Allow accessing outside of range, returning a default value
template<class E>
class Default
:
    public ListExpression<Default<E>>
{
    typename std::conditional<E::is_leaf, const E&, const E>::type expr_;

    const label size_;

    const typename E::value_type val_;

public:

    // ! Store as copy (since holds reference)
    static constexpr bool is_leaf = false;

    typedef typename E::value_type value_type;

    // construct from component. Pass through.
    Default(const E& expr)
    :
        expr_(expr),
        size_(expr_.size()),
        val_(Zero)
    {}

    // construct from components
    Default(const label size, const value_type val, const E& expr)
    :
        expr_(expr),
        size_(size),
        val_(val)
    {}

    // construct from components
    template<class... Args>
    Default(const label size, const value_type val, const Args&... args)
    :
        expr_(args...),
        size_(size),
        val_(val)
    {}

    auto size() const noexcept
    {
        return size_;
    }

    value_type operator[](const label i) const
    {
        if (i < 0 || i >= expr_.size())
        {
            return val_;
        }
        return expr_[i];
    }

    bool is_alias(const void* lst) const noexcept
    {
        return expr_.is_alias(lst);
    }

    typedef typename E::const_iterator const_iterator;
    // struct const_iterator
    // {
    //     // Minimalistic const_iter implementation
    //     typedef typename E::const_iterator ConstIter;
    //     typedef typename std::iterator_traits<ConstIter> ConstIterTraits;
    //     using difference_type = typename ConstIterTraits::difference_type;
    //     using value_type = typename ConstIterTraits::value_type;
    //     using pointer = typename ConstIterTraits::pointer;
    //     using reference = typename ConstIterTraits::reference;
    //     using iterator_category = std::random_access_iterator_tag;

    //     ConstIter base_;

    //     const_iterator(ConstIter it) : base_(it) {}
    //     const_iterator() : base_() {}
    //     ConstIter& operator=(ConstIter it)
    //     {
    //         base_ = it;
    //         return base_;
    //     }

    // // Indexing operators

    //     value_type operator*() const
    //     {
    //         return *base_;
    //     }
    //     //value_type operator->() const
    //     //{
    //     //    return *base_;
    //     //}
    //     value_type operator[](const difference_type i) const
    //     {
    //        return *(base_+i);
    //     }

    // // Increment/Decrement

    //     const_iterator& operator++()
    //     {
    //         ++base_;
    //         return *this;
    //     }
    //     const_iterator operator++(int)
    //     {
    //         const_iterator old(*this);
    //         ++base_;
    //         return old;
    //     }
    //     const_iterator& operator--()
    //     {
    //         --base_;
    //         return *this;
    //     }
    //     const_iterator operator--(int)
    //     {
    //         const_iterator old(*this);
    //         --base_;
    //         return old;
    //     }

    // // Arithmetic operators

    //     const_iterator operator+(const difference_type i) const
    //     {
    //         return const_iterator(base_+i);
    //     }
    //     const_iterator operator-(const difference_type i) const
    //     {
    //         return const_iterator(base_-i);
    //     }
    //     difference_type operator-(const const_iterator& it) const
    //     {
    //         return (base_ - it.base_);
    //     }

    // // Comparison operators

    //     bool operator==(const const_iterator& it) const
    //     {
    //         return (base_ == it.base_);
    //     }
    //     bool operator!=(const const_iterator& it) const
    //     {
    //         return (base_ != it.base_);
    //     }
    //     bool operator<(const const_iterator& it) const
    //     {
    //         return (base_ < it.base_);
    //     }
    //     bool operator<=(const const_iterator& it) const
    //     {
    //         return (base_ <= it.base_);
    //     }
    //     bool operator>(const const_iterator& it) const
    //     {
    //         return (base_ > it.base_);
    //     }
    //     bool operator>=(const const_iterator& it) const
    //     {
    //         return (base_ >= it.base_);
    //     }

    // // Compound assignment

    //     const_iterator& operator+=(const difference_type n)
    //     {
    //         base_ += n;
    //         return *this;
    //     }
    //     const_iterator& operator-=(const difference_type n)
    //     {
    //         base_ -= n;
    //         return *this;
    //     }
    // };
    auto cbegin() const
    {
        return const_iterator(E::cbegin());
    }
    auto cend() const
    {
        return const_iterator(E::cend());
    }
};


/*---------------------------------------------------------------------------*\
                       Class fvMatrixRefWrap Declaration
\*---------------------------------------------------------------------------*/

//- Expression wrap of non-const reference to fvMatrix
template<class Matrix>
class fvMatrixRefWrap
:
    public fvMatrixExpression
    <
        fvMatrixRefWrap<Matrix>,
        ListRefWrap<scalar>, // diag
        ListRefWrap<scalar>, // upper
        ListRefWrap<scalar>, // lower
        ListRefWrap<typename Matrix::psiFieldType::value_type>, // boundary
        GeometricFieldRefWrap<typename Matrix::faceFluxFieldType>,// fflux
        ListRefWrap<typename Matrix::psiFieldType::value_type>  // source
    >
{
public:

    static constexpr bool is_leaf = false;  //true;

    //- The fvMatrix type
    typedef Matrix this_type;

    //- Type to return for internal field
    typedef typename Matrix::psiFieldType::value_type value_type;

    //- Type to return for containers
    typedef ListRefWrap<scalar> DiagExpr;
    typedef ListRefWrap<scalar> UpperExpr;
    typedef ListRefWrap<scalar> LowerExpr;
    typedef ListRefWrap<value_type> BoundaryExpr;
    typedef GeometricFieldRefWrap
    <
        typename Matrix::faceFluxFieldType
    > FaceFluxExpr;
    typedef ListRefWrap<value_type> SourceExpr;


private:

    this_type& elems_;


public:

    //- Copy construct
    fvMatrixRefWrap(this_type& elems)
    :
        fvMatrixExpression
        <
            fvMatrixRefWrap<Matrix>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            BoundaryExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            elems.hasDiag(),
            elems.hasUpper(),
            elems.hasLower(),
            elems.hasFaceFluxCorrection(),
            elems.dimensions()
        ),
        elems_(elems)
    {}

    //- Move construct
    fvMatrixRefWrap(this_type&& elems)
    :
        fvMatrixExpression
        <
            fvMatrixRefWrap<Matrix>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            BoundaryExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            elems.hasDiag(),
            elems.hasUpper(),
            elems.hasLower(),
            elems.hasFaceFluxCorrection(),
            elems.dimensions()
        ),
        elems_(elems)
    {}

    // Construct from Matrix, forcing its evaluation.
    template<typename E>
    fvMatrixRefWrap
    (
        this_type& elems,
        const fvMatrixExpression
        <
            E,
            typename E::DiagExpr,
            typename E::UpperExpr,
            typename E::LowerExpr,
            typename E::BoundaryExpr,
            typename E::FaceFluxExpr,
            typename E::SourceExpr
        >& expr
    )
    :
        fvMatrixExpression
        <
            fvMatrixRefWrap<Matrix>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            BoundaryExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            expr.hasDiag(),
            expr.hasUpper(),
            expr.hasLower(),
            elems.hasFaceFluxCorrection(),
            elems.dimensions()
        ),
        elems_(elems)
    {
        expr.evaluate(elems_);
    }

    //- Assignment
    template<typename E>
    void operator=
    (
        const fvMatrixExpression
        <
            E,
            typename E::DiagExpr,
            typename E::UpperExpr,
            typename E::LowerExpr,
            typename E::BoundaryExpr,
            typename E::FaceFluxExpr,
            typename E::SourceExpr
        >& expr
    )
    {
        expr.evaluate(elems_);
    }

    //- Evaluate and return as GeoField. Rename to evaluate to make it clear
    //- it takes time? Or leave as indexing for convenience?
    template<typename E>
    this_type& evaluate
    (
        const fvMatrixExpression
        <
            E,
            typename E::DiagExpr,
            typename E::UpperExpr,
            typename E::LowerExpr,
            typename E::BoundaryExpr,
            typename E::FaceFluxExpr,
            typename E::SourceExpr
        >& expr
    )
    {
        return expr.evaluate(elems_);
    }


    // lduMatrix part

        DiagExpr diag() const
        {
            return elems_.diag();
        }

        DiagExpr diag()
        {
            return elems_.diag();
        }

        UpperExpr upper() const
        {
            return elems_.upper();
        }

        UpperExpr upper()
        {
            return elems_.upper();
        }

        LowerExpr lower() const
        {
            return elems_.lower();
        }

        LowerExpr lower()
        {
            return elems_.lower();
        }


    // fvMatrix

        SourceExpr source() const
        {
            return elems_.source();
        }

        SourceExpr source()
        {
            return elems_.source();
        }

        BoundaryExpr internalCoeffs(const label i)
        {
            return
            (
                elems_.internalCoeffs().set(i)
              ? elems_.internalCoeffs()[i]
              : List<value_type>::null()
            );
        }

        BoundaryExpr internalCoeffs(const label i) const
        {
            return
            (
                elems_.internalCoeffs().set(i)
              ? elems_.internalCoeffs()[i]
              : List<value_type>::null()
            );
        }

        BoundaryExpr boundaryCoeffs(const label i)
        {
            return
            (
                elems_.boundaryCoeffs().set(i)
              ? elems_.boundaryCoeffs()[i]
              : List<value_type>::null()
            );
        }

        BoundaryExpr boundaryCoeffs(const label i) const
        {
            return
            (
                elems_.boundaryCoeffs().set(i)
              ? elems_.boundaryCoeffs()[i]
              : List<value_type>::null()
            );
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return *const_cast<this_type&>
            (
                elems_
            ).faceFluxCorrectionPtr().get();
        }

        FaceFluxExpr faceFluxCorrection()
        {
            return *elems_.faceFluxCorrectionPtr().get();
        }

        // template<class Matrix>
        Matrix& evaluate(Matrix& m, const scalar scale) const
        {
            return m;
        }
};


/*---------------------------------------------------------------------------*\
                    Class fvMatrixConstRefWrap Declaration
\*---------------------------------------------------------------------------*/

//- Expression wrap of const reference to fvMatrix
// TBD: have Default for faceFluxFieldType as well.
template<class Matrix>
class fvMatrixConstRefWrap
:
    public fvMatrixExpression
    <
        fvMatrixConstRefWrap<Matrix>,
        Default<ListConstRefWrap<scalar>>, // diag
        Default<ListConstRefWrap<scalar>>, // upper
        Default<ListConstRefWrap<scalar>>, // lower
        Default<ListConstRefWrap<typename Matrix::psiFieldType::value_type>>,
        GeometricFieldConstRefWrap<typename Matrix::faceFluxFieldType>, // fflux
        ListConstRefWrap<typename Matrix::psiFieldType::value_type>  // source
    >
{
public:

    static constexpr bool is_leaf = false;  //true;

    //- The fvMatrix type
    typedef Matrix this_type;

    //- Type to return for internal field
    typedef typename Matrix::psiFieldType::value_type value_type;

    //- Type to return for containers
    typedef Default<ListConstRefWrap<scalar>> DiagExpr;
    typedef Default<ListConstRefWrap<scalar>> UpperExpr;
    typedef Default<ListConstRefWrap<scalar>> LowerExpr;
    typedef Default<ListConstRefWrap<value_type>> BoundaryExpr;
    typedef GeometricFieldConstRefWrap
    <
        typename Matrix::faceFluxFieldType
    > FaceFluxExpr;
    typedef ListConstRefWrap<value_type> SourceExpr;

private:

    const this_type& elems_;


public:

    // Construct from components
    fvMatrixConstRefWrap(const this_type& elems)
    :
        fvMatrixExpression
        <
            fvMatrixConstRefWrap<Matrix>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            BoundaryExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            elems.hasDiag(),
            elems.hasUpper(),
            elems.hasLower(),
            elems.hasFaceFluxCorrection(),
            elems.dimensions()
        ),
        elems_(elems)
    {}


    // lduMatrix part

        DiagExpr diag() const
        {
            // Ensure diag exists. hasDiag() can be true for a whole
            // expression even if the underlying lduMatrix has no diag.
            if (!elems_.hasDiag())
            {
                #ifdef DEBUG_ET
                Pout<< error::demangle<fvMatrixConstRefWrap<Matrix>>()
                    << " : dummy diag" << nl;
                #endif
                const label sz = elems_.lduAddr().size();
                return DiagExpr(sz, Zero, List<scalar>::null());
            }
            return DiagExpr(elems_.diag());
        }

        UpperExpr upper() const
        {
            if (!elems_.hasUpper())
            {
                #ifdef DEBUG_ET
                Pout<< error::demangle<fvMatrixConstRefWrap<Matrix>>()
                    << " : dummy upper" << nl;
                #endif
                const label sz = elems_.lduAddr().lowerAddr().size();
                return UpperExpr(sz, Zero, List<scalar>::null());
            }
            return UpperExpr(elems_.upper());
        }

        LowerExpr lower() const
        {
            if (!elems_.hasLower())
            {
                #ifdef DEBUG_ET
                Pout<< error::demangle<fvMatrixConstRefWrap<Matrix>>()
                    << " : dummy lower" << nl;
                #endif
                const label sz = elems_.lduAddr().lowerAddr().size();
                return LowerExpr(sz, Zero, List<scalar>::null());
            }
            return LowerExpr(elems_.lower());
        }


    // fvMatrix

        SourceExpr source() const
        {
            return elems_.source();
        }

        BoundaryExpr internalCoeffs(const label i) const
        {
            const label sz = elems_.psi().mesh().boundary()[i].size();
            return
            (
                elems_.internalCoeffs().set(i)
              ? BoundaryExpr(elems_.internalCoeffs()[i])
              : BoundaryExpr(sz, Zero, List<value_type>::null())
            );
        }

        BoundaryExpr boundaryCoeffs(const label i) const
        {
            const label sz = elems_.psi().mesh().boundary()[i].size();
            return
            (
                elems_.boundaryCoeffs().set(i)
              ? BoundaryExpr(elems_.boundaryCoeffs()[i])
              : BoundaryExpr(sz, Zero, List<value_type>::null())
            );
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            typename Matrix::faceFluxFieldType& ffc =
                *const_cast<this_type&>(elems_).faceFluxCorrectionPtr();
            return ffc;
        }

        Matrix& evaluate(Matrix& m, const scalar scale) const
        {
            return m;
        }
};


/*---------------------------------------------------------------------------*\
                    Class fvMatrixConstTmpWrap Declaration
\*---------------------------------------------------------------------------*/

//- Expression wrap of const reference to fvMatrix
template<class Matrix>
class fvMatrixConstTmpWrap
:
    public fvMatrixExpression
    <
        fvMatrixConstTmpWrap<Matrix>,
        Default<ListConstRefWrap<scalar>>, // diag
        Default<ListConstRefWrap<scalar>>, // upper
        Default<ListConstRefWrap<scalar>>, // lower
        Default<ListConstRefWrap<typename Matrix::psiFieldType::value_type>>,
        GeometricFieldConstRefWrap<typename Matrix::faceFluxFieldType>, // fflux
        ListConstRefWrap<typename Matrix::psiFieldType::value_type>  // source
    >
{
public:

    static constexpr bool is_leaf = false;  //true;

    //- The fvMatrix type
    typedef tmp<Matrix> this_type;

    //- Type to return for internal field
    typedef typename Matrix::psiFieldType::value_type value_type;

    //- Type to return for containers
    typedef Default<ListConstRefWrap<scalar>> DiagExpr;
    typedef Default<ListConstRefWrap<scalar>> UpperExpr;
    typedef Default<ListConstRefWrap<scalar>> LowerExpr;
    typedef Default<ListConstRefWrap<value_type>> BoundaryExpr;
    typedef GeometricFieldConstRefWrap
    <
        typename Matrix::faceFluxFieldType
    > FaceFluxExpr;
    typedef ListConstRefWrap<value_type> SourceExpr;


private:

    const this_type elems_;


public:

    // Construct from components
    fvMatrixConstTmpWrap(const this_type& elems)
    :
        fvMatrixExpression
        <
            fvMatrixConstTmpWrap<Matrix>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            BoundaryExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            elems().hasDiag(),
            elems().hasUpper(),
            elems().hasLower(),
            elems().hasFaceFluxCorrection(),
            elems().dimensions()
        ),
        elems_(elems)
    {}


    // lduMatrix part

        DiagExpr diag() const
        {
            // Ensure diag exists. hasDiag() can be true for a whole
            // expression even if the underlying lduMatrix has no diag.
            if (!elems_().hasDiag())
            {
                #ifdef DEBUG_ET
                Pout<< error::demangle<fvMatrixConstRefWrap<Matrix>>()
                    << " : dummy diag" << nl;
                #endif
                const label sz = elems_().lduAddr().size();
                return DiagExpr(sz, Zero, List<scalar>::null());
            }
            return DiagExpr(elems_().diag());
        }

        UpperExpr upper() const
        {
            if (!elems_().hasUpper())
            {
                #ifdef DEBUG_ET
                Pout<< error::demangle<fvMatrixConstRefWrap<Matrix>>()
                    << " : dummy upper" << nl;
                #endif
                const label sz = elems_().lduAddr().lowerAddr().size();
                return UpperExpr(sz, Zero, List<scalar>::null());
            }
            return UpperExpr(elems_().upper());
        }

        LowerExpr lower() const
        {
            if (!elems_().hasLower())
            {
                #ifdef DEBUG_ET
                Pout<< error::demangle<fvMatrixConstRefWrap<Matrix>>()
                    << " : dummy lower" << nl;
                #endif
                const label sz = elems_().lduAddr().lowerAddr().size();
                return LowerExpr(sz, Zero, List<scalar>::null());
            }
            return LowerExpr(elems_().lower());
        }


    // fvMatrix

        SourceExpr source() const
        {
            return SourceExpr(elems_().source());
        }

        BoundaryExpr internalCoeffs(const label i) const
        {
            const label sz = elems_().psi().mesh().boundary()[i].size();
            return
            (
                elems_().internalCoeffs().set(i)
              ? BoundaryExpr(elems_().internalCoeffs()[i])
              : BoundaryExpr(sz, Zero, List<value_type>::null())
            );
        }

        BoundaryExpr boundaryCoeffs(const label i) const
        {
            const label sz = elems_().psi().mesh().boundary()[i].size();
            return
            (
                elems_().boundaryCoeffs().set(i)
              ? BoundaryExpr(elems_().boundaryCoeffs()[i])
              : BoundaryExpr(sz, Zero, List<value_type>::null())
            );
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return *const_cast<Matrix&>
            (
                elems_()
            ).faceFluxCorrectionPtr();
        }

        Matrix& evaluate(Matrix& m, const scalar scale) const
        {
            return m;
        }
};


/*---------------------------------------------------------------------------*\
                    Class fvMatrixConstTmpWrap Declaration
\*---------------------------------------------------------------------------*/

//- Expression wrap of const reference to fvMatrix
template<class Matrix>
class fvMatrixConstTmpWrap
:
    public fvMatrixExpression
    <
        fvMatrixConstTmpWrap<Matrix>,
        ListConstRefWrap<typename Matrix::psiFieldType::value_type>, // diag
        ListConstRefWrap<typename Matrix::psiFieldType::value_type>, // upper
        ListConstRefWrap<typename Matrix::psiFieldType::value_type>, // lower
        GeometricFieldConstRefWrap<typename Matrix::faceFluxFieldType>, // fflux
        ListConstRefWrap<typename Matrix::psiFieldType::value_type>  // source
    >
{
public:

    static constexpr bool is_leaf = false;  //true;

    //- The fvMatrix type
    typedef tmp<Matrix> this_type;

    //- Type to return for internal field
    typedef typename Matrix::psiFieldType::value_type value_type;

    //- Type to return for containers
    typedef ListConstRefWrap<value_type> DiagExpr;
    typedef ListConstRefWrap<value_type> UpperExpr;
    typedef ListConstRefWrap<value_type> LowerExpr;
    typedef GeometricFieldConstRefWrap
    <
        typename Matrix::faceFluxFieldType
    > FaceFluxExpr;
    typedef ListConstRefWrap<value_type> SourceExpr;


private:

    const this_type elems_;


public:

    // Construct from components
    fvMatrixConstTmpWrap(const this_type& elems)
    :
        fvMatrixExpression
        <
            fvMatrixConstTmpWrap<Matrix>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            elems().hasDiag(),
            elems().hasUpper(),
            elems().hasLower(),
            elems().hasFaceFluxCorrection(),
            elems().dimensions()
        ),
        elems_(elems)
    {}


    // lduMatrix part

        DiagExpr diag() const
        {
            return elems_().diag();
        }

        UpperExpr upper() const
        {
            return elems_().upper();
        }

        LowerExpr lower() const
        {
            return elems_().lower();
        }


    // fvMatrix

        SourceExpr source() const
        {
            return elems_().source();
        }

        UpperExpr internalCoeffs(const label i) const
        {
            return
            (
                elems_().internalCoeffs().set(i)
              ? elems_().internalCoeffs()[i]
              : List<value_type>::null()
            );
        }

        UpperExpr boundaryCoeffs(const label i) const
        {
            return
            (
                elems_().boundaryCoeffs().set(i)
              ? elems_().boundaryCoeffs()[i]
              : List<value_type>::null()
            );
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return *const_cast<Matrix&>
            (
                elems_()
            ).faceFluxCorrectionPtr();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Expressions on fvMatrices
// ~~~~~~~~~~~~~~~~~~~~~~~~~

//- Helper to filter dimensions for +- operations. Used since Su,Sp (when
//  constructed from field) are dimless.
inline dimensionSet filterDims(const dimensionSet& u, const dimensionSet& v)
{
    if (u == dimless)
    {
        return v;
    }
    else if (v == dimless)
    {
        return u;
    }
    else
    {
        return u;
    }
}


// operator+
template<typename E1, typename E2>
class fvm_add
:
    public fvMatrixExpression
    <
        fvm_add<E1, E2>,
        List_add                    // DiagExpr
        <
            typename E1::DiagExpr,
            typename E2::DiagExpr
        >,
        List_add                    // UpperExpr
        <
            typename E1::UpperExpr,
            typename E2::UpperExpr
        >,
        List_add                    // LowerExpr
        <
            typename E1::LowerExpr,
            typename E2::LowerExpr
        >,
        List_add                    // BoundaryExpr
        <
            typename E1::BoundaryExpr,
            typename E2::BoundaryExpr
        >,
        GF_add                      // FaceFluxExpr
        <
            typename E1::FaceFluxExpr,
            typename E2::FaceFluxExpr
        >,
        List_add                    // SourceExpr
        <
            typename E1::SourceExpr,
            typename E2::SourceExpr
        >
    >
{
    // cref if leaf, copy otherwise
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type v_;

public:
    static constexpr bool is_leaf = false;

    //- Type to return for internal field
    typedef typename E1::value_type value_type;

    //- Type to return for containers
    typedef List_add
    <
        typename E1::DiagExpr,
        typename E2::DiagExpr
    > DiagExpr;
    typedef List_add
    <
        typename E1::UpperExpr,
        typename E2::UpperExpr
    > UpperExpr;
    typedef List_add
    <
        typename E1::LowerExpr,
        typename E2::LowerExpr
    > LowerExpr;
    typedef List_add
    <
        typename E1::BoundaryExpr,
        typename E2::BoundaryExpr
    > BoundaryExpr;
    typedef GF_add
    <
        typename E1::FaceFluxExpr,
        typename E2::FaceFluxExpr
    > FaceFluxExpr;
    typedef List_add
    <
        typename E1::SourceExpr,
        typename E2::SourceExpr
    > SourceExpr;

    fvm_add(const E1& u, const E2& v)
    :
        fvMatrixExpression
        <
            fvm_add<E1, E2>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            BoundaryExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            u.hasDiag() || v.hasDiag(),
            u.hasUpper() || v.hasUpper(),
            u.hasLower() || v.hasLower(),
            u.hasFaceFluxCorrection() || v.hasFaceFluxCorrection(),
            filterDims(u.dimensions(), v.dimensions())
        ),
        u_(u),
        v_(v)
    {}

    // lduMatrix part

        DiagExpr diag() const
        {
            //return u_.diag() + v_.diag();
            return DiagExpr(u_.diag(), v_.diag());
        }

        UpperExpr upper() const
        {
            return UpperExpr(u_.upper(), v_.upper());
        }

        LowerExpr lower() const
        {
            return LowerExpr(u_.lower(), v_.lower());
        }

    // fvMatrix

        SourceExpr source() const
        {
            return SourceExpr(u_.source(), v_.source());
        }

        BoundaryExpr boundaryCoeffs(const label i) const
        {
            return BoundaryExpr(u_.boundaryCoeffs(i), v_.boundaryCoeffs(i));
        }

        BoundaryExpr internalCoeffs(const label i) const
        {
            return BoundaryExpr(u_.internalCoeffs(i), v_.internalCoeffs(i));
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return FaceFluxExpr(u_.faceFluxCorrection(), v_.faceFluxCorrection());
        }

        template<class Matrix>
        Matrix& evaluate(Matrix& m, const scalar scale) const
        {
            u_.evaluate(m, scale);
            v_.evaluate(m, scale);
            return m;
        }
};
template<typename E1, typename E2>
fvm_add<E1, E2>
operator+
(
    fvMatrixExpression
    <
        E1,
        typename E1::DiagExpr,
        typename E1::UpperExpr,
        typename E1::LowerExpr,
        typename E1::BoundaryExpr,
        typename E1::FaceFluxExpr,
        typename E1::SourceExpr
    > const& u,
    fvMatrixExpression
    <
        E2,
        typename E2::DiagExpr,
        typename E2::UpperExpr,
        typename E2::LowerExpr,
        typename E2::BoundaryExpr,
        typename E2::FaceFluxExpr,
        typename E2::SourceExpr
    > const& v
)
{
    return fvm_add<E1, E2>
    (
        static_cast<const E1&>(u),
        static_cast<const E2&>(v)
    );
}
template<typename E1, typename E2>
fvm_add<E1, fvMatrixConstTmpWrap<fvMatrix<typename E1::value_type>>>
operator+
(
    fvMatrixExpression
    <
        E1,
        typename E1::DiagExpr,
        typename E1::UpperExpr,
        typename E1::LowerExpr,
        typename E1::BoundaryExpr,
        typename E1::FaceFluxExpr,
        typename E1::SourceExpr
    > const& u,
    tmp<fvMatrix<typename E1::value_type>> const& v
)
{
    typedef typename E1::value_type value_type;
    return u + fvMatrixConstTmpWrap<fvMatrix<value_type>>(v);
}
template<typename E1, typename E2>
fvm_add<E1, E2>
operator+
(
    tmp<fvMatrix<typename E1::value_type>> const& u,
    fvMatrixExpression
    <
        E1,
        typename E1::DiagExpr,
        typename E1::UpperExpr,
        typename E1::LowerExpr,
        typename E1::BoundaryExpr,
        typename E1::FaceFluxExpr,
        typename E1::SourceExpr
    > const& v
)
{
    typedef typename E1::value_type value_type;
    return fvMatrixConstTmpWrap<fvMatrix<value_type>>(u) + v;
}


// operator-, operator==
template<typename E1, typename E2>
class fvm_subtract
:
    public fvMatrixExpression
    <
        fvm_subtract<E1, E2>,
        List_subtract                   // diag
        <
            typename E1::DiagExpr,
            typename E2::DiagExpr
        >,
        List_subtract                   // upper
        <
            typename E1::UpperExpr,
            typename E2::UpperExpr
        >,
        List_subtract                   // lower
        <
            typename E1::LowerExpr,
            typename E2::LowerExpr
        >,
        List_subtract                   // boundary
        <
            typename E1::BoundaryExpr,
            typename E2::BoundaryExpr
        >,
        GF_subtract                     // faceflux
        <
            typename E1::FaceFluxExpr,
            typename E2::FaceFluxExpr
        >,
        List_subtract                   // source
        <
            typename E1::SourceExpr,
            typename E2::SourceExpr
        >
    >
{
    // cref if leaf, copy otherwise
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type v_;

public:
    static constexpr bool is_leaf = false;

    //- Type to return for internal field
    typedef typename E1::value_type value_type;

    //- Type to return for containers
    typedef List_subtract
    <
        typename E1::DiagExpr,
        typename E2::DiagExpr
    > DiagExpr;
    typedef List_subtract
    <
        typename E1::UpperExpr,
        typename E2::UpperExpr
    > UpperExpr;
    typedef List_subtract
    <
        typename E1::LowerExpr,
        typename E2::LowerExpr
    > LowerExpr;
    typedef List_subtract
    <
        typename E1::BoundaryExpr,
        typename E2::BoundaryExpr
    > BoundaryExpr;
    typedef GF_subtract
    <
        typename E1::FaceFluxExpr,
        typename E2::FaceFluxExpr
    > FaceFluxExpr;
    typedef List_subtract
    <
        typename E1::SourceExpr,
        typename E2::SourceExpr
    > SourceExpr;

    fvm_subtract(const E1& u, const E2& v)
    :
        fvMatrixExpression
        <
            fvm_subtract<E1, E2>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            BoundaryExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            u.hasDiag() || v.hasDiag(),
            u.hasUpper() || v.hasUpper(),
            u.hasLower() || v.hasLower(),
            u.hasFaceFluxCorrection() || v.hasFaceFluxCorrection(),
            filterDims(u.dimensions(), v.dimensions())
        ),
        u_(u),
        v_(v)
    {}


    // lduMatrix part

        DiagExpr diag() const
        {
            return DiagExpr(u_.diag(), v_.diag());
        }

        UpperExpr upper() const
        {
            return UpperExpr(u_.upper(), v_.upper());
        }

        LowerExpr lower() const
        {
            return LowerExpr(u_.lower(), v_.lower());
        }


    // fvMatrix

        SourceExpr source() const
        {
            return SourceExpr(u_.source(), v_.source());
        }

        BoundaryExpr internalCoeffs(const label i) const
        {
            return BoundaryExpr(u_.internalCoeffs(i), v_.internalCoeffs(i));
        }

        BoundaryExpr boundaryCoeffs(const label i) const
        {
            return BoundaryExpr(u_.boundaryCoeffs(i), v_.boundaryCoeffs(i));
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return FaceFluxExpr(u_.faceFluxCorrection(), v_.faceFluxCorrection());
        }

        template<class Matrix>
        Matrix& evaluate(Matrix& m, const scalar scale) const
        {
            u_.evaluate(m, scale);
            v_.evaluate(m, -1.0*scale);
            return m;
        }
};
template<typename E1, typename E2>
fvm_subtract<E1, E2>
operator-
(
    fvMatrixExpression
    <
        E1,
        typename E1::DiagExpr,
        typename E1::UpperExpr,
        typename E1::LowerExpr,
        typename E1::BoundaryExpr,
        typename E1::FaceFluxExpr,
        typename E1::SourceExpr
    > const& u,
    fvMatrixExpression
    <
        E2,
        typename E2::DiagExpr,
        typename E2::UpperExpr,
        typename E2::LowerExpr,
        typename E2::BoundaryExpr,
        typename E2::FaceFluxExpr,
        typename E2::SourceExpr
    > const& v
)
{
    return fvm_subtract<E1, E2>
    (
        static_cast<const E1&>(u),
        static_cast<const E2&>(v)
    );
}
template<typename E1, typename E2>
fvm_subtract<E1, E2>
operator==
(
    fvMatrixExpression
    <
        E1,
        typename E1::DiagExpr,
        typename E1::UpperExpr,
        typename E1::LowerExpr,
        typename E1::BoundaryExpr,
        typename E1::FaceFluxExpr,
        typename E1::SourceExpr
    > const& u,
    fvMatrixExpression
    <
        E2,
        typename E2::DiagExpr,
        typename E2::UpperExpr,
        typename E2::LowerExpr,
        typename E2::BoundaryExpr,
        typename E2::FaceFluxExpr,
        typename E2::SourceExpr
    > const& v
)
{
    // workaround: Su,Sp etc return dimless
    if (v.dimensions() == dimless)
    {
        const_cast<dimensionSet&>(v.dimensions()).reset(u.dimensions());
    }
    else if (u.dimensions() == dimless)
    {
        const_cast<dimensionSet&>(u.dimensions()).reset(v.dimensions());
    }
    return fvm_subtract<E1, E2>
    (
        static_cast<const E1&>(u),
        static_cast<const E2&>(v)
    );
}

// Forward Declarations
template<typename E1>
auto Su
(
    const E1&,
    const GeometricField<typename E1::value_type, fvPatchField, volMesh>&
);

template<typename E1, typename E2>
fvm_subtract<E1, E2>
operator==
(
    fvMatrixExpression
    <
        E1,
        typename E1::DiagExpr,
        typename E1::UpperExpr,
        typename E1::LowerExpr,
        typename E1::BoundaryExpr,
        typename E1::FaceFluxExpr,
        typename E1::SourceExpr
    > const& u,
    tmp
    <
        GeometricField<typename E1::value_type, fvPatchField, volMesh>
    > const& v
)
{
    typedef typename E1::value_type value_type;
    typedef GeometricField<value_type, fvPatchField, volMesh> FldType;
    return fvm_subtract<E1, E2>(u, Su<FldType>(v, v()));
}
template<typename E1, typename E2>
fvm_subtract<E1, E2>
operator-
(
    fvMatrixExpression
    <
        E1,
        typename E1::DiagExpr,
        typename E1::UpperExpr,
        typename E1::LowerExpr,
        typename E1::BoundaryExpr,
        typename E1::FaceFluxExpr,
        typename E1::SourceExpr
    > const& u,
    tmp<fvMatrix<typename E1::value_type>> const& v
)
{
    typedef typename E1::value_type value_type;
    return fvm_subtract<E1, E2>
    (
        u,
        fvMatrixConstTmpWrap<fvMatrix<value_type>>(v)
    );
}
template<typename E1, typename E2>
fvm_subtract<E1, E2>
operator-
(
    tmp<fvMatrix<typename E1::value_type>> const& u,
    fvMatrixExpression
    <
        E1,
        typename E1::DiagExpr,
        typename E1::UpperExpr,
        typename E1::LowerExpr,
        typename E1::BoundaryExpr,
        typename E1::FaceFluxExpr,
        typename E1::SourceExpr
    > const& v
)
{
    typedef typename E1::value_type value_type;
    return fvm_subtract<E1, E2>
    (
        fvMatrixConstTmpWrap<fvMatrix<value_type>>(u),
        v
    );
}


// Negation : operator-
template<typename E1>
class fvm_negate
:
    public fvMatrixExpression
    <
        fvm_negate<E1>,
        List_negate<typename E1::DiagExpr>,
        List_negate<typename E1::UpperExpr>,
        List_negate<typename E1::LowerExpr>,
        List_negate<typename E1::BoundaryExpr>,
        GF_negate<typename E1::FaceFluxExpr>,
        List_negate<typename E1::SourceExpr>
    >
{
    // cref if leaf, copy otherwise
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;

public:
    static constexpr bool is_leaf = false;

    //- Type to return for internal field
    typedef typename E1::value_type value_type;

    //- Type to return for containers
    typedef List_negate<typename E1::DiagExpr> DiagExpr;
    typedef List_negate<typename E1::UpperExpr> UpperExpr;
    typedef List_negate<typename E1::LowerExpr> LowerExpr;
    typedef List_negate<typename E1::BoundaryExpr> BoundaryExpr;
    typedef GF_negate<typename E1::FaceFluxExpr> FaceFluxExpr;
    typedef List_negate<typename E1::SourceExpr> SourceExpr;

    fvm_negate(const E1& u)
    :
        fvMatrixExpression
        <
            fvm_negate<E1>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            BoundaryExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            u.hasDiag(),
            u.hasUpper(),
            u.hasLower(),
            u.hasFaceFluxCorrection(),
            u.dimensions()
        ),
        u_(u)
    {}


    // lduMatrix part

        // Note: conversion constructor from diag to List_negate<> takes
        // care of negation
        DiagExpr diag() const
        {
            return u_.diag();
        }

        UpperExpr upper() const
        {
            return u_.upper();
        }

        LowerExpr lower() const
        {
            return u_.lower();
        }


    // fvMatrix

        SourceExpr source() const
        {
            return u_.source();
        }

        BoundaryExpr internalCoeffs(const label i) const
        {
            return u_.internalCoeffs(i);
        }

        BoundaryExpr boundaryCoeffs(const label i) const
        {
            return u_.boundaryCoeffs(i);
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return u_.faceFluxCorrection();
        }

        template<class Matrix>
        Matrix& evaluate(Matrix& m, const scalar scale) const
        {
            u_.evaluate(m, -1.0*scale);
            return m;
        }
};
template<typename E1>
auto operator-
(
    fvMatrixExpression
    <
        E1,
        typename E1::DiagExpr,
        typename E1::UpperExpr,
        typename E1::LowerExpr,
        typename E1::BoundaryExpr,
        typename E1::FaceFluxExpr,
        typename E1::SourceExpr
    > const& u
)
{
    return fvm_negate<E1>(static_cast<const E1&>(u));
}


// Source manipulation
// ~~~~~~~~~~~~~~~~~~~

//- fvMatrix, internal parts and dimensions only. No upper/lower/boundaryCoeffs
//  etc. Used for source manipulation.
template<class E1, class E2>
class fvMatrixInternal
:
    public fvMatrixExpression
    <
        fvMatrixInternal<E1, E2>,
        E1, // Diag
        UniformListWrap<scalar>, // Upper (= 0)
        UniformListWrap<scalar>, // Lower (= 0)
        UniformListWrap<typename E2::value_type>, // Boundary (= 0)
        GeometricFieldConstRefWrap
        <
            GeometricField<typename E2::value_type, fvsPatchField, surfaceMesh>
        >,  // FaceFlux, Not used
        E2  // Source
    >
{
public:

    static constexpr bool is_leaf = false;  //true;

    //- Type to return for internal field
    typedef typename E2::value_type value_type;

    //- Type to return for containers
    typedef E1 DiagExpr;

    //- Uniform 0 for upper/lower
    typedef UniformListWrap<scalar> UpperExpr;
    typedef UniformListWrap<scalar> LowerExpr;
    typedef UniformListWrap<value_type> BoundaryExpr;

    typedef GeometricField<value_type, fvsPatchField, surfaceMesh>
            faceFluxFieldType;
    typedef GeometricFieldConstRefWrap<faceFluxFieldType> FaceFluxExpr;
    typedef E2 SourceExpr;


private:

    //- Underlying matrix - used for sizes only
    const fvMatrix<value_type>& m_;

    const E1 diag_;
    const SourceExpr source_;


public:

    // Construct from component expressions
    fvMatrixInternal
    (
        const fvMatrix<value_type>& m,
        const DiagExpr& diag,
        const SourceExpr& source
    )
    :
        fvMatrixExpression
        <
            fvMatrixInternal<DiagExpr, SourceExpr>,
            DiagExpr,
            UpperExpr,       // upper not yet used
            LowerExpr,       // lower not yet used
            BoundaryExpr,    // boundary not yet used
            FaceFluxExpr,    // face flux not used
            SourceExpr
        >
        (
            true,       //elems.hasDiag(),
            false,      //elems.hasUpper(),
            false,      //elems.hasLower(),
            false,      //elems.hasFaceFluxCorrection(),
            m.dimensions()  //elems.dimensions()
        ),
        m_(m),
        diag_(diag),
        source_(source)
    {}

    // Construct from component expressions
    fvMatrixInternal
    (
        const fvMatrix<value_type>& m,
        const DiagExpr&& diag,
        const SourceExpr&& source
    )
    :
        fvMatrixExpression
        <
            fvMatrixInternal<DiagExpr, SourceExpr>,
            DiagExpr,
            UpperExpr,       // upper not yet used
            LowerExpr,       // lower not yet used
            BoundaryExpr,    // boundary not yet used
            FaceFluxExpr,       // faceflux not yet used
            SourceExpr
        >
        (
            true,       //elems.hasDiag(),
            false,      //elems.hasUpper(),
            false,      //elems.hasLower(),
            false,      //elems.hasFaceFluxCorrection(),
            m.dimensions()  //elems.dimensions()
        ),
        m_(m),
        diag_(diag),
        source_(source)
    {}


    // lduMatrix part

        DiagExpr diag() const
        {
            return diag_;
        }

        UpperExpr upper() const
        {
            return UpperExpr(m_.upper().size(), scalar(0));
        }

        LowerExpr lower() const
        {
            return LowerExpr(m_.lower().size(), scalar(0));
        }


    // fvMatrix

        SourceExpr source() const
        {
            return source_;
        }

        BoundaryExpr internalCoeffs(const label i) const
        {
            return BoundaryExpr
            (
                m_.internalCoeffs()[i].size(),
                Zero
            );
        }

        BoundaryExpr boundaryCoeffs(const label i) const
        {
            return BoundaryExpr
            (
                m_.boundaryCoeffs()[i].size(),
                Zero
            );
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return
            (
                *const_cast<fvMatrix<value_type>&>(m_).
                    faceFluxCorrectionPtr().get()
            );
        }

        template<class Matrix>
        Matrix& evaluate(Matrix& m, const scalar scale) const
        {
            return m;
        }
};

//- Sources using matrix - provides dimensions
template<typename E1, typename Type>
auto SuSp
(
    const E1& expression,
    const fvMatrix<Type>& m
)
{
    // Equivalent of fvm::SuSp

    //- Wrap of constant as a list expression
    typedef Expression::UniformListWrap<scalar> constant;

    // Wrap constants
    const constant constantZero(m.psi().size(), 0.0);
    const constant constantMinOne(m.psi().size(), -1.0);

    // Wrap mesh volume
    const ListConstRefWrap<scalar> V(m.psi().mesh().V());
    return fvMatrixInternal
    (
        m,  // m.dimensions(),
        V*max(expression, constantZero),        // diag
        constantMinOne                          // source
       *V
       *min(expression, constantZero)
       *m.psi().internalField().expr()
    );
}
// Provide ListExpression for diag; dimensions from matrix
template<typename E1, typename Type>
auto Sp
(
    const E1& expression,
    const fvMatrix<Type>& m
)
{
    // Equivalent of fvm::Sp

    typedef typename fvMatrix<Type>::value_type value_type;
    //- Wrap of constant as a list expression
    typedef Expression::UniformListWrap<value_type> constant;

    // Wrap zero
    const constant constantZero(m.psi().size(), 0.0);

    // Wrap mesh volume
    const ListConstRefWrap<scalar> V(m.psi().mesh().V());

    return fvMatrixInternal
    (
        m,              // m.dimensions(),
        V*expression,   // diag
        constantZero    // source
    );
}
template<typename E1>
auto Su
(
    const E1& expression,
    const fvMatrix<typename E1::value_type>& m
)
{
    // Equivalent of fvm::Su

    //- Wrap of constant as a list expression
    typedef Expression::UniformListWrap<scalar> constant;

    // Wrap zero
    const constant constantZero(m.psi().size(), Zero);

    // Wrap mesh volume
    const ListConstRefWrap<scalar> V(m.psi().mesh().V());

    return fvMatrixInternal
    (
        m,              // m.dimensions(),
        constantZero,   // Still provides an expression for diag. TBD.
        -V*expression
    );
}


//- Variant of fvMatrixInternal:
//      - diag, source only
//      - returns dimless (handled in outer expressions)
//  Used for Su, Sp, SuSp to avoid needing matrix info
template<class E1, class E2>
class fvMatrixInternal2
:
    public fvMatrixExpression
    <
        fvMatrixInternal2<E1, E2>,
        E1, // Diag
        UniformListWrap<scalar>, // Upper, Not used
        UniformListWrap<scalar>, // Lower, Not used
        UniformListWrap<typename E2::value_type>, // Boundary, Not used
        UniformGeometricFieldWrap2<typename E2::value_type>,  // FFlux, Not used
        E2  // Source
    >
{
public:

    static constexpr bool is_leaf = false;  //true;

    //- Type to return for internal field
    typedef typename E2::value_type value_type;

    //- Type to return for containers
    typedef E1 DiagExpr;
    typedef UniformListWrap<scalar> UpperExpr;
    typedef UniformListWrap<scalar> LowerExpr;
    typedef UniformListWrap<value_type> BoundaryExpr;
    typedef UniformGeometricFieldWrap2<value_type> FaceFluxExpr;
    typedef E2 SourceExpr;


private:

    const GeometricField<value_type, fvPatchField, volMesh>& psi_;
    const E1 diag_;
    const SourceExpr source_;


public:

    // Construct from component expressions
    fvMatrixInternal2
    (
        const GeometricField<value_type, fvPatchField, volMesh>& psi,
        const DiagExpr& diag,
        const SourceExpr& source
    )
    :
        fvMatrixExpression
        <
            fvMatrixInternal2<DiagExpr, SourceExpr>,
            DiagExpr,
            UpperExpr,       // upper not yet used
            LowerExpr,       // lower not yet used
            BoundaryExpr,    // boundary not yet used
            FaceFluxExpr,    // face flux not used
            SourceExpr
        >
        (
            true,       //elems.hasDiag(),
            false,      //elems.hasUpper(),
            false,      //elems.hasLower(),
            false,      //elems.hasFaceFluxCorrection(),
            dimless     //elems.dimensions()
        ),
        psi_(psi),
        diag_(diag),
        source_(source)
    {}


    // lduMatrix part

        DiagExpr diag() const
        {
            return diag_;
        }

        UpperExpr upper() const
        {
            const label sz = psi_.mesh().lduAddr().lowerAddr().size();
            return UpperExpr(sz, scalar(0));
        }

        LowerExpr lower() const
        {
            const label sz = psi_.mesh().lduAddr().lowerAddr().size();
            return LowerExpr(sz, scalar(0));
        }


    // fvMatrix

        SourceExpr source() const
        {
            return source_;
        }

        BoundaryExpr internalCoeffs(const label i) const
        {
            const label sz = psi_.mesh().boundary()[i].size();
            return BoundaryExpr(sz, Zero);
        }

        BoundaryExpr boundaryCoeffs(const label i) const
        {
            const label sz = psi_.mesh().boundary()[i].size();
            return BoundaryExpr(sz, Zero);
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return FaceFluxExpr(Zero);
        }

        template<class Matrix>
        Matrix& evaluate(Matrix& m, const scalar scale) const
        {
            return m;
        }
};
template<typename E1, typename Type>
auto SuSp
(
    const E1& expression,
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    //- Wrap of constant as a list expression
    typedef Expression::UniformListWrap<scalar> constant;

    // Wrap constants
    const constant constantZero(fld.size(), 0.0);
    const constant constantMinOne(fld.size(), -1.0);

    // Wrap mesh volume
    const ListConstRefWrap<scalar> V(fld.mesh().V());

    return fvMatrixInternal2
    (
        fld,
        V*max(expression, constantZero),        // diag
        V*min(expression, constantZero)
         *constantMinOne
         *fld.internalField().expr()            // source
    );
}
template<typename DiagField, typename Type>
auto SuSp
(
    const tmp<DiagField> te,
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    return SuSp(te().internalField().expr(), fld);
}
template<typename E1, typename Type>
auto Sp
(
    const E1& diag,
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    // Wrap zero
    const Expression::UniformListWrap<scalar> constantZero(fld.size(), 0.0);

    // Wrap mesh volume
    const ListConstRefWrap<scalar> V(fld.mesh().V());

    return fvMatrixInternal2(fld, V*diag, constantZero);
}
template<typename GeoField>
auto Sp
(
    const tmp<GeoField> tdiag,
    const GeoField& fld
)
{
    return Sp(tdiag().internalField().expr(), fld);
}
template<typename E1>
auto Su
(
    const E1& source,
    const GeometricField<typename E1::value_type, fvPatchField, volMesh>& fld
)
{
    // Wrap zero
    const Expression::UniformListWrap
    <
        typename E1::value_type
    > constantZero(fld.size(), Zero);

    // Wrap mesh volume
    const ListConstRefWrap<scalar> V(fld.mesh().V());

    return fvMatrixInternal2(fld, constantZero, -V*source);
}
template<typename GeoField>
auto Su
(
    const tmp<GeoField> tsource,
    const GeoField& fld
)
{
    return Su(tsource().internalField().expr(), fld);
}


// Some discretisation using expressions
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<typename E1, typename E2>
class interpolate
:
    public GeometricFieldExpression
    <
        interpolate<E1, E2>,

        //- Type to return for internal field
        List_interpolate
        <
            IndirectConstWrap<typename E1::IntExpr>,
            IndirectConstWrap<typename E1::IntExpr>,
            typename E2::IntExpr
        >,

        //- Type to return for uncoupled patch field
        typename E1::UncoupledPatchExpr,

        //- Type to return for coupled patch field
        List_interpolate
        <
            IndirectConstWrap<typename E1::IntExpr>,
            typename E1::CoupledPatchExpr,
            typename E2::CoupledPatchExpr
        >,

        //- Type of actual values
        typename E1::value_type
    >
{
    //- Internal values
    typename std::conditional
    <
        E1::is_leaf,
        const E1&,
        const E1
    >::type cellVals_;

    //- Interpolation weights
    typename std::conditional
    <
        E2::is_leaf,
        const E2&,
        const E2
    >::type faceWeights_;

    //- Addressing container
    const fvMesh& mesh_;

public:

    static constexpr bool is_leaf = false;

    //- Type to return for internal field
    typedef List_interpolate
    <
        IndirectConstWrap<typename E1::IntExpr>,
        IndirectConstWrap<typename E1::IntExpr>,
        typename E2::IntExpr
    > IntExpr;

    //- Type to return for uncoupled patch field
    typedef typename E1::UncoupledPatchExpr UncoupledPatchExpr;

    //- Type to return for coupled patch field
    typedef List_interpolate
    <
        IndirectConstWrap<typename E1::IntExpr>,
        typename E1::CoupledPatchExpr,
        typename E2::CoupledPatchExpr
    > CoupledPatchExpr;

    //- Type of actual values
    typedef typename E1::value_type value_type;


    // Construct from components
    interpolate
    (
        const E1& cellVals,
        const E2& faceWeights,
        const fvMesh& mesh
    )
    :
        GeometricFieldExpression
        <
            interpolate<E1, E2>,
            IntExpr,
            UncoupledPatchExpr,
            CoupledPatchExpr,
            value_type
        >
        (
            cellVals.dimensions()*faceWeights.dimensions(),
            faceWeights.oriented()
        ),
        cellVals_(cellVals),
        faceWeights_(faceWeights),
        mesh_(mesh)
    {}

    value_type operator[](const label i) const
    {
        const auto ownVal = cellVals_[mesh_.owner()[i]];
        const auto neiVal = cellVals_[mesh_.neighbour()[i]];
        return faceWeights_[i]*(ownVal-neiVal) + neiVal;
    }

    label size() const
    {
        return faceWeights_.size(); // number of internal faces
    }

    IntExpr internalField() const
    {
        return IntExpr
        (
            IndirectConstWrap<typename E1::IntExpr>
            (
                cellVals_.internalField(),
                mesh_.owner()
            ),
            IndirectConstWrap<typename E1::IntExpr>
            (
                cellVals_.internalField(),
                mesh_.neighbour()
            ),
            faceWeights_.internalField()
        );
    }

    UncoupledPatchExpr patchField(const label i) const
    {
        return UncoupledPatchExpr(cellVals_.patchField(i));
    }

    CoupledPatchExpr coupledPatchField(const label i) const
    {
        return CoupledPatchExpr
        (
            IndirectConstWrap<typename E1::IntExpr>
            (
                cellVals_.internalField(),
                mesh_.boundaryMesh()[i].faceCells()
            ),
            cellVals_.access
            (
                [&](const auto& fld, const label i)
                {
                    return ListConstTmpWrap<Field<value_type>>
                    (
                        fld.boundaryField()[i].patchNeighbourField()
                    );
                },
                i
            ),
            faceWeights_.coupledPatchField(i)
        );
    }
};


template<typename E1>
class linearInterpolate
:
    public interpolate
    <
        E1,
        GeometricFieldConstRefWrap<surfaceScalarField>
    >
{
public:

    // Construct from components
    linearInterpolate
    (
        const E1& cellVals,
        const fvMesh& mesh
    )
    :
        interpolate
        <
            E1,
            GeometricFieldConstRefWrap<surfaceScalarField>
        >
        (
            cellVals,
            mesh.surfaceInterpolation::weights().expr(),
            mesh
        )
    {}
};


template<class Type, class E1, class E2>
void fvmLaplacianUncorrected
(
    fvMatrix<Type>& fvm,
    const E1& gammaMagSf,
    const E2& deltaCoeffs
)
{
    typedef Expression::ListConstRefWrap<Type> expr;
    typedef Expression::ListRefWrap<Type> evaluator;

    const auto& vf = fvm.psi();

    // Set dimensions. Or should we just check?
    const_cast<dimensionSet&>(fvm.dimensions()).reset
    (
        deltaCoeffs.dimensions()
      * gammaMagSf.dimensions()
      * vf.dimensions()
    );

    Expression::ListRefWrap<scalar>
    (
        fvm.upper(),
        deltaCoeffs.internalField()*gammaMagSf.internalField()
    );
    fvm.negSumDiag();

    forAll(vf.boundaryField(), patchi)
    {
        const auto& pvf = vf.boundaryField()[patchi];
        auto& intCoeffs = fvm.internalCoeffs()[patchi];
        auto& bouCoeffs = fvm.boundaryCoeffs()[patchi];

        if (pvf.coupled())
        {
            // Evaluate deltaCoeffs
            scalarField pDeltaCoeffs;
            ListRefWrap<scalar>(pDeltaCoeffs, deltaCoeffs.patchField(patchi));

            //intCoeffs = pGamma*pvf.gradientInternalCoeffs(pDeltaCoeffs);
            evaluator
            (
                intCoeffs,
                expr(pvf.gradientInternalCoeffs(pDeltaCoeffs)())
               *gammaMagSf.patchField(patchi)
            );
            //bouCoeffs = -pGamma*pvf.gradientBoundaryCoeffs(pDeltaCoeffs);
            evaluator
            (
                bouCoeffs,
               -expr(pvf.gradientBoundaryCoeffs(pDeltaCoeffs)())
               *gammaMagSf.patchField(patchi)
            );
        }
        else
        {
            //intCoeffs = pGamma*pvf.gradientInternalCoeffs();
            evaluator
            (
                intCoeffs,
                expr(pvf.gradientInternalCoeffs()())
               *gammaMagSf.patchField(patchi)
            );
            //bouCoeffs = -pGamma*pvf.gradientBoundaryCoeffs();
            evaluator
            (
                bouCoeffs,
               -expr(pvf.gradientBoundaryCoeffs()())
               *gammaMagSf.patchField(patchi)
            );
        }
    }
}

//XXXXXXXX
template<class Type, class EvalOp>
class DirectOp
:
    public fvMatrixExpression
    <
        DirectOp<Type, EvalOp>,
        UniformListWrap<scalar>, // Diag
        UniformListWrap<scalar>, // Upper
        UniformListWrap<scalar>, // Lower
        UniformListWrap<Type>, // Boundary, Not used
        UniformGeometricFieldWrap2<Type>,  // FFlux, Not used
        UniformListWrap<Type>
    >
{
public:

    static constexpr bool is_leaf = false;  //true;

    //- Type to return for internal field
    typedef Type value_type;

    //- Type to return for containers
    typedef UniformListWrap<scalar> DiagExpr;
    typedef UniformListWrap<scalar> UpperExpr;
    typedef UniformListWrap<scalar> LowerExpr;
    typedef UniformListWrap<value_type> BoundaryExpr;
    typedef UniformGeometricFieldWrap2<value_type> FaceFluxExpr;
    typedef UniformListWrap<value_type> SourceExpr;


private:

    const GeometricField<value_type, fvPatchField, volMesh>& psi_;

    const EvalOp evalOp_;

public:

    // Construct from component expressions
    DirectOp
    (
        const GeometricField<value_type, fvPatchField, volMesh>& psi,
        const EvalOp& evalOp,
        const dimensionSet& dimensions
    )
    :
        fvMatrixExpression
        <
            DirectOp<Type, EvalOp>,
            DiagExpr,
            UpperExpr,       // upper not yet used
            LowerExpr,       // lower not yet used
            BoundaryExpr,    // boundary not yet used
            FaceFluxExpr,    // face flux not used
            SourceExpr
        >
        (
            false,      //elems.hasDiag(),
            false,      //elems.hasUpper(),
            false,      //elems.hasLower(),
            false,      //elems.hasFaceFluxCorrection(),
            dimensions  //elems.dimensions()
        ),
        psi_(psi),
        evalOp_(evalOp)
    {}


    // lduMatrix part

        DiagExpr diag() const
        {
            const label sz = psi_.mesh().lduAddr().size();
            return DiagExpr(sz, Zero);
        }

        UpperExpr upper() const
        {
            const label sz = psi_.mesh().lduAddr().lowerAddr().size();
            return UpperExpr(sz, Zero);
        }

        LowerExpr lower() const
        {
            const label sz = psi_.mesh().lduAddr().lowerAddr().size();
            return LowerExpr(sz, Zero);
        }


    // fvMatrix

        SourceExpr source() const
        {
            const label sz = psi_.mesh().lduAddr().size();
            return SourceExpr(sz, Zero);
        }

        BoundaryExpr internalCoeffs(const label i) const
        {
            const label sz = psi_.mesh().boundary()[i].size();
            return BoundaryExpr(sz, Zero);
        }

        BoundaryExpr boundaryCoeffs(const label i) const
        {
            const label sz = psi_.mesh().boundary()[i].size();
            return BoundaryExpr(sz, Zero);
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return FaceFluxExpr(Zero);
        }

        template<class Matrix>
        Matrix& evaluate(Matrix& m, const scalar scale) const
        {
            // Dimensions supplied in constructor - feed through to matrix
            const_cast<dimensionSet&>(m.dimensions()).reset(this->dimensions());
            // Fill in matrix using supplied operator
            evalOp_(m, scale, psi_);
            return m;
        }
};
struct ddtOp
{
    template<class Type>
    void operator()
    (
        fvMatrix<Type>& fvm,
        const scalar scale,
        const GeometricField<Type, fvPatchField, volMesh>& vf
    ) const
    {
        fvm::ddt<Type>(fvm, scale, vf);
    }
};
template<class Type>
class ddt
:
    public DirectOp<Type, ddtOp>
{
public:

    ddt(const GeometricField<Type, fvPatchField, volMesh>& psi)
    :
        DirectOp<Type, ddtOp>
        (
            psi,
            ddtOp(),
            psi.dimensions()*dimVol/dimTime
        )
    {}
};
struct laplacianOp
{
    template<class Type>
    void operator()
    (
        fvMatrix<Type>& fvm,
        const scalar scale,
        const GeometricField<Type, fvPatchField, volMesh>& vf
    ) const
    {
        fvm = scale*fvm::laplacian(vf);
    }
};
template<class Type>
class laplacian
:
    public DirectOp<Type, laplacianOp>
{
public:

    laplacian(const GeometricField<Type, fvPatchField, volMesh>& psi)
    :
        DirectOp<Type, laplacianOp>
        (
            psi,
            laplacianOp(),
            psi.dimensions()*dimLength*dimLength/dimTime
        )
    {}
};//XXXXXXXX


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Expression
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
