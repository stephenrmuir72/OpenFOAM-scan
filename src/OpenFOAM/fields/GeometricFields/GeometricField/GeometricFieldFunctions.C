/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2025 OpenCFD Ltd.
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

#include "GeometricFieldReuseFunctions.H"
#include "GeometricFieldOps.H"

#define TEMPLATE \
    template<class Type, template<class> class PatchField, class GeoMesh>
#include "GeometricFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template
<
    class ResultType,
    class Type,
    template<class> class PatchField,
    class GeoMesh,
    class UnaryOp,
    class FieldOp
>
void applyUnaryBoundary
(
    GeometricBoundaryField<ResultType, PatchField, GeoMesh>& result,
    const GeometricBoundaryField<Type, PatchField, GeoMesh>& source,
    const UnaryOp& unaryOp,
    const FieldOp& fieldOp
)
{
    if (result.hasDeferred() || source.hasDeferred())
    {
        result.transform(source, unaryOp);
    }
    else
    {
        fieldOp(result, source);
    }
}

template<class Type, template<class> class PatchField, class GeoMesh>
void component
(
    GeometricField
    <
        typename GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
    >& result,
    const GeometricField<Type, PatchField, GeoMesh>& f1,
    const direction d
)
{
    component(result.primitiveFieldRef(), f1.primitiveField(), d);
    applyUnaryBoundary
    (
        result.boundaryFieldRef(),
        f1.boundaryField(),
        [d](const auto& value) { return Foam::component(value, d); },
        [d](auto& resultBoundary, const auto& sourceBoundary)
        {
            Foam::component(resultBoundary, sourceBoundary, d);
        }
    );
    result.oriented() = f1.oriented();
    if (GeometricBoundaryField<Type, PatchField, GeoMesh>::debug)
    {
        result.boundaryField().check();
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void T
(
     GeometricField<Type, PatchField, GeoMesh>& result,
     const GeometricField<Type, PatchField, GeoMesh>& f1
)
{
    T(result.primitiveFieldRef(), f1.primitiveField());
    applyUnaryBoundary
    (
        result.boundaryFieldRef(),
        f1.boundaryField(),
        [](const auto& value) { return Foam::T(value); },
        [](auto& resultBoundary, const auto& sourceBoundary)
        {
            Foam::T(resultBoundary, sourceBoundary);
        }
    );
    result.oriented() = f1.oriented();
    if (GeometricBoundaryField<Type, PatchField, GeoMesh>::debug)
    {
        result.boundaryField().check();
    }
}


template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh,
    direction r
>
void pow
(
    GeometricField
    <typename powProduct<Type, r>::type, PatchField, GeoMesh>& result,
    const GeometricField<Type, PatchField, GeoMesh>& f1
)
{
    pow(result.primitiveFieldRef(), f1.primitiveField(), r);
    applyUnaryBoundary
    (
        result.boundaryFieldRef(),
        f1.boundaryField(),
        [](const auto& value) { return Foam::pow(value, r); },
        [](auto& resultBoundary, const auto& sourceBoundary)
        {
            Foam::pow(resultBoundary, sourceBoundary, r);
        }
    );
    result.oriented() = pow(f1.oriented(), r);
    result.correctLocalBoundaryConditions();
    if (GeometricBoundaryField<Type, PatchField, GeoMesh>::debug)
    {
        result.boundaryField().check();
    }
}


template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh,
    direction r
>
tmp<GeometricField<typename powProduct<Type, r>::type, PatchField, GeoMesh>>
pow
(
    const GeometricField<Type, PatchField, GeoMesh>& f1,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type resultType;

    auto tres =
        reuseTmpGeometricField<resultType, Type, PatchField, GeoMesh>::New
        (
            f1,
            "pow(" + f1.name() + ',' + Foam::name(r) + ')',
            pow(f1.dimensions(), r)
        );

    pow<Type, r, PatchField, GeoMesh>(tres.ref(), f1);

    return tres;
}


template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh,
    direction r
>
tmp<GeometricField<typename powProduct<Type, r>::type, PatchField, GeoMesh>>
pow
(
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf1,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type resultType;

    const auto& f1 = tf1();

    auto tres =
        reuseTmpGeometricField<resultType, Type, PatchField, GeoMesh>::New
        (
            tf1,
            "pow(" + f1.name() + ',' + Foam::name(r) + ')',
            pow(f1.dimensions(), r)
        );

    pow<Type, r, PatchField, GeoMesh>(tres.ref(), f1);

    tf1.clear();
    return tres;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void sqr
(
    GeometricField
    <
        typename outerProduct<Type, Type>::type, PatchField, GeoMesh
    >& result,
    const GeometricField<Type, PatchField, GeoMesh>& f1
)
{
    sqr(result.primitiveFieldRef(), f1.primitiveField());
    applyUnaryBoundary
    (
        result.boundaryFieldRef(),
        f1.boundaryField(),
        [](const auto& value) { return Foam::sqr(value); },
        [](auto& resultBoundary, const auto& sourceBoundary)
        {
            Foam::sqr(resultBoundary, sourceBoundary);
        }
    );
    result.oriented() = sqr(f1.oriented());
    result.correctLocalBoundaryConditions();
    if (GeometricBoundaryField<Type, PatchField, GeoMesh>::debug)
    {
        result.boundaryField().check();
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp
<
    GeometricField
    <
        typename outerProduct<Type, Type>::type,
        PatchField,
        GeoMesh
    >
>
sqr(const GeometricField<Type, PatchField, GeoMesh>& f1)
{
    typedef typename outerProduct<Type, Type>::type resultType;

    auto tres =
        reuseTmpGeometricField<resultType, Type, PatchField, GeoMesh>::New
        (
            f1,
            "sqr(" + f1.name() + ')',
            sqr(f1.dimensions())
        );

    sqr(tres.ref(), f1);

    return tres;
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp
<
    GeometricField
    <
        typename outerProduct<Type, Type>::type,
        PatchField,
        GeoMesh
    >
>
sqr(const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf1)
{
    typedef typename outerProduct<Type, Type>::type resultType;

    const auto& f1 = tf1();

    auto tres =
        reuseTmpGeometricField<resultType, Type, PatchField, GeoMesh>::New
        (
            tf1,
            "sqr(" + f1.name() + ')',
            sqr(f1.dimensions())
        );

    sqr(tres.ref(), f1);

    tf1.clear();
    return tres;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void magSqr
(
    GeometricField<typename typeOfMag<Type>::type, PatchField, GeoMesh>& result,
    const GeometricField<Type, PatchField, GeoMesh>& f1
)
{
    magSqr(result.primitiveFieldRef(), f1.primitiveField());
    applyUnaryBoundary
    (
        result.boundaryFieldRef(),
        f1.boundaryField(),
        [](const auto& value) { return Foam::magSqr(value); },
        [](auto& resultBoundary, const auto& sourceBoundary)
        {
            Foam::magSqr(resultBoundary, sourceBoundary);
        }
    );
    result.oriented() = magSqr(f1.oriented());
    result.correctLocalBoundaryConditions();
    if (GeometricBoundaryField<Type, PatchField, GeoMesh>::debug)
    {
        result.boundaryField().check();
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<typename typeOfMag<Type>::type, PatchField, GeoMesh>>
magSqr
(
    const GeometricField<Type, PatchField, GeoMesh>& f1
)
{
    typedef typename typeOfMag<Type>::type resultType;

    auto tres =
        reuseTmpGeometricField<resultType, Type, PatchField, GeoMesh>::New
        (
            f1,
            "magSqr(" + f1.name() + ')',
            sqr(f1.dimensions())
        );

    magSqr(tres.ref(), f1);

    return tres;
}

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<typename typeOfMag<Type>::type, PatchField, GeoMesh>>
magSqr
(
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf1
)
{
    auto tres = magSqr(tf1.cref());
    tf1.clear();

    return tres;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void mag
(
    GeometricField<typename typeOfMag<Type>::type, PatchField, GeoMesh>& result,
    const GeometricField<Type, PatchField, GeoMesh>& f1
)
{
    mag(result.primitiveFieldRef(), f1.primitiveField());
    applyUnaryBoundary
    (
        result.boundaryFieldRef(),
        f1.boundaryField(),
        [](const auto& value) { return Foam::mag(value); },
        [](auto& resultBoundary, const auto& sourceBoundary)
        {
            Foam::mag(resultBoundary, sourceBoundary);
        }
    );
    result.oriented() = mag(f1.oriented());
    result.correctLocalBoundaryConditions();
    if (GeometricBoundaryField<Type, PatchField, GeoMesh>::debug)
    {
        result.boundaryField().check();
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<typename typeOfMag<Type>::type, PatchField, GeoMesh>>
mag
(
    const GeometricField<Type, PatchField, GeoMesh>& f1
)
{
    typedef typename typeOfMag<Type>::type resultType;

    auto tres =
        reuseTmpGeometricField<resultType, Type, PatchField, GeoMesh>::New
        (
            f1,
            "mag(" + f1.name() + ')',
            f1.dimensions()
        );

    mag(tres.ref(), f1);

    return tres;
}

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<typename typeOfMag<Type>::type, PatchField, GeoMesh>>
mag
(
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf1
)
{
    auto tres = mag(tf1.cref());
    tf1.clear();

    return tres;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void cmptAv
(
    GeometricField
    <
        typename GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
    >& result,
    const GeometricField<Type, PatchField, GeoMesh>& f1
)
{
    cmptAv(result.primitiveFieldRef(), f1.primitiveField());
    applyUnaryBoundary
    (
        result.boundaryFieldRef(),
        f1.boundaryField(),
        [](const auto& value) { return Foam::cmptAv(value); },
        [](auto& resultBoundary, const auto& sourceBoundary)
        {
            Foam::cmptAv(resultBoundary, sourceBoundary);
        }
    );
    result.oriented() = cmptAv(f1.oriented());
    if (GeometricBoundaryField<Type, PatchField, GeoMesh>::debug)
    {
        result.boundaryField().check();
    }
}

template<class Type, template<class> class PatchField, class GeoMesh>
tmp
<
    GeometricField
    <
        typename GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
    >
>
cmptAv(const GeometricField<Type, PatchField, GeoMesh>& f1)
{
    typedef typename
        GeometricField<Type, PatchField, GeoMesh>::cmptType resultType;

    auto tres =
        reuseTmpGeometricField<resultType, Type, PatchField, GeoMesh>::New
        (
            f1,
            "cmptAv(" + f1.name() + ')',
            f1.dimensions()
        );

    cmptAv(tres.ref(), f1);

    return tres;
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp
<
    GeometricField
    <
        typename GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
    >
>
cmptAv(const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf1)
{
    auto tres = cmptAv(tf1.cref());
    tf1.clear();

    return tres;
}


#define UNARY_REDUCTION_FUNCTION_WITH_BOUNDARY(ReturnType, Func, BinaryOp)     \
                                                                               \
template<class Type, template<class> class PatchField, class GeoMesh>          \
dimensioned<ReturnType> Func                                                   \
(                                                                              \
    const GeometricField<Type, PatchField, GeoMesh>& f1,                       \
    const label comm                                                           \
)                                                                              \
{                                                                              \
    return dimensioned<ReturnType>                                             \
    (                                                                          \
        #Func "(" + f1.name() + ')',                                           \
        f1.dimensions(),                                                       \
        returnReduce                                                           \
        (                                                                      \
            Foam::Func                                                         \
            (                                                                  \
                Foam::Func(f1.primitiveField()),                               \
                Foam::Func(f1.boundaryField())                                 \
            ),                                                                 \
            BinaryOp<ReturnType>(),                                            \
            UPstream::msgType(),                                               \
            comm                                                               \
        )                                                                      \
    );                                                                         \
}                                                                              \
                                                                               \
template<class Type, template<class> class PatchField, class GeoMesh>          \
dimensioned<ReturnType> Func                                                   \
(                                                                              \
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf1,                 \
    const label comm                                                           \
)                                                                              \
{                                                                              \
    dimensioned<ReturnType> res = Func(tf1(), comm);                           \
    tf1.clear();                                                               \
    return res;                                                                \
}


UNARY_REDUCTION_FUNCTION_WITH_BOUNDARY(Type, max, maxOp)
UNARY_REDUCTION_FUNCTION_WITH_BOUNDARY(Type, min, minOp)
UNARY_REDUCTION_FUNCTION_WITH_BOUNDARY(MinMax<Type>, minMax, plusOp)
UNARY_REDUCTION_FUNCTION_WITH_BOUNDARY(scalarMinMax, minMaxMag, plusOp)

#undef UNARY_REDUCTION_FUNCTION_WITH_BOUNDARY


// Forward to DimensionedField directly (same name)
#define UNARY_REDUCTION_FUNCTION(ReturnType, Func)                             \
                                                                               \
template<class Type, template<class> class PatchField, class GeoMesh>          \
dimensioned<ReturnType> Func                                                   \
(                                                                              \
    const GeometricField<Type, PatchField, GeoMesh>& f1,                       \
    const label comm                                                           \
)                                                                              \
{                                                                              \
    return Func(f1.internalField(), comm);                                     \
}                                                                              \
                                                                               \
template<class Type, template<class> class PatchField, class GeoMesh>          \
dimensioned<ReturnType> Func                                                   \
(                                                                              \
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf1,                 \
    const label comm                                                           \
)                                                                              \
{                                                                              \
    auto result = Func(tf1(), comm);                                           \
    tf1.clear();                                                               \
    return result;                                                             \
}

UNARY_REDUCTION_FUNCTION(Type, sum)
UNARY_REDUCTION_FUNCTION(Type, average)
UNARY_REDUCTION_FUNCTION(typename typeOfMag<Type>::type, sumMag)

#undef UNARY_REDUCTION_FUNCTION


BINARY_FUNCTION(Type, Type, Type, max)
BINARY_FUNCTION(Type, Type, Type, min)
BINARY_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_FUNCTION(Type, Type, Type, cmptDivide)

BINARY_TYPE_FUNCTION(Type, Type, Type, max)
BINARY_TYPE_FUNCTION(Type, Type, Type, min)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptDivide)


// ------------------------------------------------------------------------- //

// Clamp Methods

template<class Type, template<class> class PatchField, class GeoMesh>
void clamp
(
    GeometricField<Type, PatchField, GeoMesh>& result,
    const GeometricField<Type, PatchField, GeoMesh>& f1,
    Foam::zero_one
)
{
    const MinMax<Type> range(Foam::zero_one{});

    clamp(result.primitiveFieldRef(), f1.primitiveField(), range);
    clamp(result.boundaryFieldRef(), f1.boundaryField(), range);
    result.oriented() = f1.oriented();
    result.correctLocalBoundaryConditions();
    if (GeometricBoundaryField<Type, PatchField, GeoMesh>::debug)
    {
        result.boundaryField().check();
    }
}

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>>
clamp
(
    const GeometricField<Type, PatchField, GeoMesh>& f1,
    Foam::zero_one
)
{
    auto tres =
        reuseTmpGeometricField<Type, Type, PatchField, GeoMesh>::New
        (
            f1,
            "clamp01(" + f1.name() + ')',
            f1.dimensions()
        );

    clamp(tres.ref(), f1, Foam::zero_one{});

    return tres;
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>>
clamp
(
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf1,
    Foam::zero_one
)
{
    const auto& f1 = tf1();

    auto tres =
        reuseTmpGeometricField<Type, Type, PatchField, GeoMesh>::New
        (
            tf1,
            "clamp01(" + f1.name() + ')',
            f1.dimensions()
        );

    clamp(tres.ref(), f1, Foam::zero_one{});

    tf1.clear();
    return tres;
}

BINARY_TYPE_FUNCTION_FS(Type, Type, MinMax<Type>, clamp)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

TERNARY_FUNCTION(Type, Type, Type, scalar, lerp)
TERNARY_TYPE_FUNCTION_FFS(Type, Type, Type, scalar, lerp)


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(Type, Type, -, negate, transform)

BINARY_OPERATOR(Type, Type, scalar, *, '*', multiply)
BINARY_OPERATOR(Type, scalar, Type, *, '*', multiply)
BINARY_OPERATOR(Type, Type, scalar, /, '|', divide)

BINARY_TYPE_OPERATOR_SF(Type, scalar, Type, *, '*', multiply)
BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, *, '*', multiply)

BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, /, '|', divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define PRODUCT_OPERATOR(product, Op, OpFunc)                                  \
                                                                               \
template                                                                       \
<class Type1, class Type2, template<class> class PatchField, class GeoMesh>    \
void OpFunc                                                                    \
(                                                                              \
    GeometricField                                                             \
    <typename product<Type1, Type2>::type, PatchField, GeoMesh>& result,       \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
)                                                                              \
{                                                                              \
    if                                                                         \
    (                                                                          \
        !FieldOps::assignFlat                                                   \
        (                                                                      \
            result,                                                            \
            f1,                                                                \
            f2,                                                                \
            [](const auto& a, const auto& b) { return a Op b; }                \
        )                                                                      \
    )                                                                          \
    {                                                                          \
        Foam::OpFunc                                                          \
        (                                                                      \
            result.primitiveFieldRef(),                                        \
            f1.primitiveField(),                                               \
            f2.primitiveField()                                                \
        );                                                                     \
        if                                                                     \
        (                                                                      \
            result.boundaryField().hasDeferred()                              \
         || f1.boundaryField().hasDeferred()                                   \
         || f2.boundaryField().hasDeferred()                                   \
        )                                                                      \
        {                                                                      \
            result.boundaryFieldRef().apply                                   \
            (                                                                  \
                f1.boundaryField(),                                            \
                f2.boundaryField(),                                            \
                [](const auto& a, const auto& b) { return a Op b; }            \
            );                                                                 \
        }                                                                      \
        else                                                                   \
        {                                                                      \
            Foam::OpFunc                                                      \
            (                                                                  \
                result.boundaryFieldRef(),                                     \
                f1.boundaryField(),                                            \
                f2.boundaryField()                                             \
            );                                                                 \
        }                                                                      \
    }                                                                          \
                                                                               \
    result.oriented() = (f1.oriented() Op f2.oriented());                      \
    result.correctLocalBoundaryConditions();                                   \
    if (GeometricBoundaryField<Type1, PatchField, GeoMesh>::debug)             \
    {                                                                          \
        result.boundaryField().check();                                        \
    }                                                                          \
}                                                                              \
                                                                               \
                                                                               \
template                                                                       \
<class Type1, class Type2, template<class> class PatchField, class GeoMesh>    \
tmp                                                                            \
<                                                                              \
    GeometricField<typename product<Type1, Type2>::type, PatchField, GeoMesh>  \
>                                                                              \
operator Op                                                                    \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type resultType;                   \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<resultType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            f1,                                                                \
           '(' + f1.name() + #Op + f2.name() + ')',                            \
            (f1.dimensions() Op f2.dimensions())                               \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref(), f1, f2);                                          \
                                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template                                                                       \
<class Type1, class Type2, template<class> class PatchField, class GeoMesh>    \
tmp                                                                            \
<                                                                              \
    GeometricField<typename product<Type1, Type2>::type, PatchField, GeoMesh>  \
>                                                                              \
operator Op                                                                    \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2                 \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type resultType;                   \
                                                                               \
    const auto& f2 = tf2();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<resultType, Type2, PatchField, GeoMesh>::New    \
        (                                                                      \
            tf2,                                                               \
            '(' + f1.name() + #Op + f2.name() + ')',                           \
            (f1.dimensions() Op f2.dimensions())                               \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref(), f1, f2);                                          \
                                                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<class Type1, class Type2, template<class> class PatchField, class GeoMesh>    \
tmp                                                                            \
<                                                                              \
    GeometricField<typename product<Type1, Type2>::type, PatchField, GeoMesh>  \
>                                                                              \
operator Op                                                                    \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type resultType;                   \
                                                                               \
    const auto& f1 = tf1();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<resultType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            tf1,                                                               \
            '(' + f1.name() + #Op + f2.name() + ')',                           \
            (f1.dimensions() Op f2.dimensions())                               \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref(), f1, f2);                                          \
                                                                               \
    tf1.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<class Type1, class Type2, template<class> class PatchField, class GeoMesh>    \
tmp                                                                            \
<                                                                              \
    GeometricField<typename product<Type1, Type2>::type, PatchField, GeoMesh>  \
>                                                                              \
operator Op                                                                    \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2                 \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type resultType;                   \
                                                                               \
    const auto& f1 = tf1();                                                    \
    const auto& f2 = tf2();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpTmpGeometricField                                              \
        <resultType, Type1, Type1, Type2, PatchField, GeoMesh>::New            \
        (                                                                      \
            tf1,                                                               \
            tf2,                                                               \
            '(' + f1.name() + #Op + f2.name() + ')',                           \
            (f1.dimensions() Op f2.dimensions())                               \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref(), f1, f2);                                          \
                                                                               \
    tf1.clear();                                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<class Form, class Type, template<class> class PatchField, class GeoMesh>      \
void OpFunc                                                                    \
(                                                                              \
    GeometricField                                                             \
    <typename product<Type, Form>::type, PatchField, GeoMesh>& result,         \
    const GeometricField<Type, PatchField, GeoMesh>& f1,                       \
    const dimensioned<Form>& dvs                                               \
)                                                                              \
{                                                                              \
    Foam::OpFunc(result.primitiveFieldRef(), f1.primitiveField(), dvs.value());\
    if                                                                       \
    (                                                                          \
        result.boundaryField().hasDeferred()                                  \
     || f1.boundaryField().hasDeferred()                                       \
    )                                                                          \
    {                                                                          \
        result.boundaryFieldRef().apply                                       \
        (                                                                      \
            f1.boundaryField(),                                                \
            dvs.value(),                                                       \
            [](const auto& a, const auto& b) { return a Op b; }                \
        );                                                                     \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        Foam::OpFunc                                                          \
        (                                                                      \
            result.boundaryFieldRef(),                                         \
            f1.boundaryField(),                                                \
            dvs.value()                                                        \
        );                                                                     \
    }                                                                          \
    result.oriented() = f1.oriented();                                         \
    if (GeometricBoundaryField<Type, PatchField, GeoMesh>::debug)              \
    {                                                                          \
        result.boundaryField().check();                                        \
    }                                                                          \
}                                                                              \
                                                                               \
template                                                                       \
<class Form, class Type, template<class> class PatchField, class GeoMesh>      \
tmp<GeometricField<typename product<Type, Form>::type, PatchField, GeoMesh>>   \
operator Op                                                                    \
(                                                                              \
    const GeometricField<Type, PatchField, GeoMesh>& f1,                       \
    const dimensioned<Form>& dvs                                               \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type resultType;                     \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<resultType, Type, PatchField, GeoMesh>::New     \
        (                                                                      \
            f1,                                                                \
            '(' + f1.name() + #Op + dvs.name() + ')',                          \
            (f1.dimensions() Op dvs.dimensions())                              \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref(), f1, dvs);                                         \
                                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class Type, template<class> class PatchField,                              \
    class GeoMesh                                                              \
>                                                                              \
tmp<GeometricField<typename product<Form, Type>::type, PatchField, GeoMesh>>   \
operator Op                                                                    \
(                                                                              \
    const GeometricField<Type, PatchField, GeoMesh>& f1,                       \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    return f1 Op dimensioned<Form>(static_cast<const Form&>(vs));              \
}                                                                              \
                                                                               \
                                                                               \
template                                                                       \
<class Form, class Type, template<class> class PatchField, class GeoMesh>      \
tmp<GeometricField<typename product<Type, Form>::type, PatchField, GeoMesh>>   \
operator Op                                                                    \
(                                                                              \
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf1,                 \
    const dimensioned<Form>& dvs                                               \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type resultType;                     \
                                                                               \
    const auto& f1 = tf1();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<resultType, Type, PatchField, GeoMesh>::New     \
        (                                                                      \
            tf1,                                                               \
            '(' + f1.name() + #Op + dvs.name() + ')',                          \
            (f1.dimensions() Op dvs.dimensions())                              \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref(), f1, dvs);                                         \
                                                                               \
    tf1.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class Type, template<class> class PatchField,                              \
    class GeoMesh                                                              \
>                                                                              \
tmp<GeometricField<typename product<Form, Type>::type, PatchField, GeoMesh>>   \
operator Op                                                                    \
(                                                                              \
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf1,                 \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    return tf1 Op dimensioned<Form>(static_cast<const Form&>(vs));             \
}                                                                              \
                                                                               \
                                                                               \
template                                                                       \
<class Form, class Type, template<class> class PatchField, class GeoMesh>      \
void OpFunc                                                                    \
(                                                                              \
    GeometricField                                                             \
    <typename product<Form, Type>::type, PatchField, GeoMesh>& result,         \
    const dimensioned<Form>& dvs,                                              \
    const GeometricField<Type, PatchField, GeoMesh>& f2                        \
)                                                                              \
{                                                                              \
    Foam::OpFunc(result.primitiveFieldRef(), dvs.value(), f2.primitiveField());\
    if                                                                       \
    (                                                                          \
        result.boundaryField().hasDeferred()                                  \
     || f2.boundaryField().hasDeferred()                                       \
    )                                                                          \
    {                                                                          \
        result.boundaryFieldRef().apply                                       \
        (                                                                      \
            dvs.value(),                                                       \
            f2.boundaryField(),                                                \
            [](const auto& a, const auto& b) { return a Op b; }                \
        );                                                                     \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        Foam::OpFunc                                                          \
        (                                                                      \
            result.boundaryFieldRef(),                                         \
            dvs.value(),                                                       \
            f2.boundaryField()                                                 \
        );                                                                     \
    }                                                                          \
    result.oriented() = f2.oriented();                                         \
    if (GeometricBoundaryField<Type, PatchField, GeoMesh>::debug)              \
    {                                                                          \
        result.boundaryField().check();                                        \
    }                                                                          \
}                                                                              \
                                                                               \
template                                                                       \
<class Form, class Type, template<class> class PatchField, class GeoMesh>      \
tmp<GeometricField<typename product<Form, Type>::type, PatchField, GeoMesh>>   \
operator Op                                                                    \
(                                                                              \
    const dimensioned<Form>& dvs,                                              \
    const GeometricField<Type, PatchField, GeoMesh>& f2                        \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type resultType;                     \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<resultType, Type, PatchField, GeoMesh>::New     \
        (                                                                      \
            f2,                                                                \
            '(' + dvs.name() + #Op + f2.name() + ')',                          \
            (dvs.dimensions() Op f2.dimensions())                              \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref(), dvs, f2);                                         \
                                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class Type, template<class> class PatchField,                              \
    class GeoMesh                                                              \
>                                                                              \
tmp<GeometricField<typename product<Form, Type>::type, PatchField, GeoMesh>>   \
operator Op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const GeometricField<Type, PatchField, GeoMesh>& f2                        \
)                                                                              \
{                                                                              \
    return dimensioned<Form>(static_cast<const Form&>(vs)) Op f2;              \
}                                                                              \
                                                                               \
template                                                                       \
<class Form, class Type, template<class> class PatchField, class GeoMesh>      \
tmp<GeometricField<typename product<Form, Type>::type, PatchField, GeoMesh>>   \
operator Op                                                                    \
(                                                                              \
    const dimensioned<Form>& dvs,                                              \
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf2                  \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type resultType;                     \
                                                                               \
    const auto& f2 = tf2();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<resultType, Type, PatchField, GeoMesh>::New     \
        (                                                                      \
            tf2,                                                               \
            '(' + dvs.name() + #Op + f2.name() + ')',                          \
            (dvs.dimensions() Op f2.dimensions())                              \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref(), dvs, f2);                                         \
                                                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
template                                                                       \
<                                                                              \
    class Form,                                                                \
    class Cmpt,                                                                \
    direction nCmpt,                                                           \
    class Type, template<class> class PatchField,                              \
    class GeoMesh                                                              \
>                                                                              \
tmp<GeometricField<typename product<Form, Type>::type, PatchField, GeoMesh>>   \
operator Op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tf2                  \
)                                                                              \
{                                                                              \
    return dimensioned<Form>(static_cast<const Form&>(vs)) Op tf2;             \
}


PRODUCT_OPERATOR(typeOfSum, +, add)
PRODUCT_OPERATOR(typeOfSum, -, subtract)

PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
