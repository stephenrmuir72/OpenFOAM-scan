/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019,2023 OpenCFD Ltd.
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
    Foam::FieldOps

Description
    Various utility functions to work on geometric fields

SourceFiles
    GeometricFieldOps.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_GeometricFieldOps_H
#define Foam_GeometricFieldOps_H

#include "FieldOps.H"
#include "GeometricField.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FieldOps
{

/*---------------------------------------------------------------------------*\
                        Namespace FieldOps Declarations
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class PatchField, class GeoMesh>
label flatSize
(
    const GeometricField<Type, PatchField, GeoMesh>& field
)
{
    if constexpr
    (
        GeometricBoundaryField<Type, PatchField, GeoMesh>::supportsDeferred
    )
    {
        return
            GeoMesh::size(field.mesh())
          + GeoMesh::boundary_size(field.mesh());
    }

    return GeoMesh::size(field.mesh());
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool hasFlatStorage
(
    const GeometricField<Type, PatchField, GeoMesh>& field
)
{
    if constexpr
    (
        GeometricBoundaryField<Type, PatchField, GeoMesh>::supportsDeferred
    )
    {
        return
            FieldBase::unifiedGeometricField
         && field.primitiveField().capacity() >= flatSize(field);
    }

    return false;
}


template<class Type, template<class> class PatchField, class GeoMesh>
UList<Type> flatValues
(
    GeometricField<Type, PatchField, GeoMesh>& field
)
{
    auto& values = field.primitiveFieldRef();
    return UList<Type>(values.begin(), flatSize(field));
}


template<class Type, template<class> class PatchField, class GeoMesh>
UList<Type> flatValues
(
    const GeometricField<Type, PatchField, GeoMesh>& field
)
{
    const auto& values = field.primitiveField();
    return UList<Type>
    (
        const_cast<Type*>(values.begin()),
        flatSize(field)
    );
}


template
<
    class Tout, class T1, class UnaryOp,
    template<class> class PatchField, class GeoMesh
>
bool assignFlat
(
    GeometricField<Tout, PatchField, GeoMesh>& result,
    const GeometricField<T1, PatchField, GeoMesh>& input,
    const UnaryOp& op
)
{
    if
    (
        !hasFlatStorage(result)
     || !hasFlatStorage(input)
     || flatSize(result) != flatSize(input)
    )
    {
        return false;
    }

    auto resultValues = flatValues(result);
    const auto inputValues = flatValues(input);

    forAll(resultValues, i)
    {
        resultValues[i] = op(inputValues[i]);
    }

    return true;
}


template
<
    class Tout, class T1, class T2, class BinaryOp,
    template<class> class PatchField, class GeoMesh
>
bool assignFlat
(
    GeometricField<Tout, PatchField, GeoMesh>& result,
    const GeometricField<T1, PatchField, GeoMesh>& lhs,
    const GeometricField<T2, PatchField, GeoMesh>& rhs,
    const BinaryOp& op
)
{
    if
    (
        !hasFlatStorage(result)
     || !hasFlatStorage(lhs)
     || !hasFlatStorage(rhs)
     || flatSize(result) != flatSize(lhs)
     || flatSize(result) != flatSize(rhs)
    )
    {
        return false;
    }

    auto resultValues = flatValues(result);
    const auto lhsValues = flatValues(lhs);
    const auto rhsValues = flatValues(rhs);

    forAll(resultValues, i)
    {
        resultValues[i] = op(lhsValues[i], rhsValues[i]);
    }

    return true;
}


template
<
    class Tout, class T1, class T2, class BinaryOp,
    template<class> class PatchField, class GeoMesh
>
bool assignFlat
(
    GeometricField<Tout, PatchField, GeoMesh>& result,
    const T1& lhs,
    const GeometricField<T2, PatchField, GeoMesh>& rhs,
    const BinaryOp& op
)
{
    if
    (
        !hasFlatStorage(result)
     || !hasFlatStorage(rhs)
     || flatSize(result) != flatSize(rhs)
    )
    {
        return false;
    }

    auto resultValues = flatValues(result);
    const auto rhsValues = flatValues(rhs);

    forAll(resultValues, i)
    {
        resultValues[i] = op(lhs, rhsValues[i]);
    }

    return true;
}


template
<
    class Tout, class T1, class T2, class BinaryOp,
    template<class> class PatchField, class GeoMesh
>
bool assignFlat
(
    GeometricField<Tout, PatchField, GeoMesh>& result,
    const GeometricField<T1, PatchField, GeoMesh>& lhs,
    const T2& rhs,
    const BinaryOp& op
)
{
    if
    (
        !hasFlatStorage(result)
     || !hasFlatStorage(lhs)
     || flatSize(result) != flatSize(lhs)
    )
    {
        return false;
    }

    auto resultValues = flatValues(result);
    const auto lhsValues = flatValues(lhs);

    forAll(resultValues, i)
    {
        resultValues[i] = op(lhsValues[i], rhs);
    }

    return true;
}


template
<
    class Tout, class T1, class T2, class T3, class TernaryOp,
    template<class> class PatchField, class GeoMesh
>
bool assignFlat
(
    GeometricField<Tout, PatchField, GeoMesh>& result,
    const GeometricField<T1, PatchField, GeoMesh>& first,
    const GeometricField<T2, PatchField, GeoMesh>& second,
    const GeometricField<T3, PatchField, GeoMesh>& third,
    const TernaryOp& op
)
{
    if
    (
        !hasFlatStorage(result)
     || !hasFlatStorage(first)
     || !hasFlatStorage(second)
     || !hasFlatStorage(third)
     || flatSize(result) != flatSize(first)
     || flatSize(result) != flatSize(second)
     || flatSize(result) != flatSize(third)
    )
    {
        return false;
    }

    auto resultValues = flatValues(result);
    const auto firstValues = flatValues(first);
    const auto secondValues = flatValues(second);
    const auto thirdValues = flatValues(third);

    forAll(resultValues, i)
    {
        resultValues[i] = op
        (
            firstValues[i],
            secondValues[i],
            thirdValues[i]
        );
    }

    return true;
}


template
<
    class Tout, class T1, class T2, class T3, class TernaryOp,
    template<class> class PatchField, class GeoMesh
>
bool assignFlat
(
    GeometricField<Tout, PatchField, GeoMesh>& result,
    const GeometricField<T1, PatchField, GeoMesh>& first,
    const GeometricField<T2, PatchField, GeoMesh>& second,
    const T3& third,
    const TernaryOp& op
)
{
    if
    (
        !hasFlatStorage(result)
     || !hasFlatStorage(first)
     || !hasFlatStorage(second)
     || flatSize(result) != flatSize(first)
     || flatSize(result) != flatSize(second)
    )
    {
        return false;
    }

    auto resultValues = flatValues(result);
    const auto firstValues = flatValues(first);
    const auto secondValues = flatValues(second);

    forAll(resultValues, i)
    {
        resultValues[i] = op(firstValues[i], secondValues[i], third);
    }

    return true;
}

//- Populate a geometric field as the result of a unary operation on an input.
//  It is permissible for inputs/outputs to refer to the same field(s),
//  but partially overlapping regions are ill-defined.
template
<
    class Tout, class T1, class UnaryOp,
    template<class> class PatchField, class GeoMesh
>
void assign
(
    GeometricField<Tout, PatchField, GeoMesh>& result,
    const GeometricField<T1, PatchField, GeoMesh>& a,
    const UnaryOp& op
)
{
    if constexpr
    (
        GeometricBoundaryField<Tout, PatchField, GeoMesh>::supportsDeferred
    )
    {
        if (!assignFlat(result, a, op))
        {
            FieldOps::assign
            (
                result.primitiveFieldRef(),
                a.primitiveField(),
                op
            );

            auto& bfld = result.boundaryFieldRef();

            const label len = bfld.size();

            for (label i = 0; i < len; ++i)
            {
                FieldOps::assign
                (
                    bfld[i],
                    a.boundaryField()[i],
                    op
                );
            }
        }
    }
    else
    {
        FieldOps::assign
        (
            result.primitiveFieldRef(),
            a.primitiveField(),
            op
        );
    }

    result.correctLocalBoundaryConditions();
    if (GeometricBoundaryField<Tout, PatchField, GeoMesh>::debug)
    {
        result.boundaryField().check();
    }
}


//- Populate a geometric field from the binary operation on two inputs.
//  It is permissible for inputs/outputs to refer to the same field(s),
//  but partially overlapping regions are ill-defined.
template
<
    class Tout, class T1, class T2, class BinaryOp,
    template<class> class PatchField, class GeoMesh
>
void assign
(
    GeometricField<Tout, PatchField, GeoMesh>& result,
    const GeometricField<T1, PatchField, GeoMesh>& a,
    const GeometricField<T2, PatchField, GeoMesh>& b,
    const BinaryOp& bop
)
{
    if constexpr
    (
        GeometricBoundaryField<Tout, PatchField, GeoMesh>::supportsDeferred
    )
    {
        if (!assignFlat(result, a, b, bop))
        {
            FieldOps::assign
            (
                result.primitiveFieldRef(),
                a.primitiveField(),
                b.primitiveField(),
                bop
            );

            auto& bfld = result.boundaryFieldRef();

            const label len = bfld.size();

            for (label i = 0; i < len; ++i)
            {
                FieldOps::assign
                (
                    bfld[i],
                    a.boundaryField()[i],
                    b.boundaryField()[i],
                    bop
                );
            }
        }
    }
    else
    {
        FieldOps::assign
        (
            result.primitiveFieldRef(),
            a.primitiveField(),
            b.primitiveField(),
            bop
        );
    }

    result.correctLocalBoundaryConditions();
    if (GeometricBoundaryField<Tout, PatchField, GeoMesh>::debug)
    {
        result.boundaryField().check();
    }
}


//- Emulate a ternary operation, selecting values from a or b
//- depending on the binary predicate.
template
<
    class T, class BinaryOp,
    template<class> class PatchField, class GeoMesh
>
void ternary
(
    GeometricField<T, PatchField, GeoMesh>& result,
    const GeometricField<T, PatchField, GeoMesh>& a,
    const GeometricField<T, PatchField, GeoMesh>& b,
    const BinaryOp& bop
)
{
    if constexpr
    (
        GeometricBoundaryField<T, PatchField, GeoMesh>::supportsDeferred
    )
    {
        if
        (
            !assignFlat
            (
                result,
                a,
                b,
                [&bop](const T& lhs, const T& rhs)
                {
                    return bop(lhs, rhs) ? lhs : rhs;
                }
            )
        )
        {
            FieldOps::ternary
            (
                result.primitiveFieldRef(),
                a.primitiveField(),
                b.primitiveField(),
                bop
            );

            auto& bfld = result.boundaryFieldRef();

            const label len = bfld.size();

            for (label i = 0; i < len; ++i)
            {
                FieldOps::ternary
                (
                    bfld[i],
                    a.boundaryField()[i],
                    b.boundaryField()[i],
                    bop
                );
            }
        }
    }
    else
    {
        FieldOps::ternary
        (
            result.primitiveFieldRef(),
            a.primitiveField(),
            b.primitiveField(),
            bop
        );
    }
    result.correctLocalBoundaryConditions();
    if (GeometricBoundaryField<T, PatchField, GeoMesh>::debug)
    {
        result.boundaryField().check();
    }
}


//- Emulate a ternary operation, selecting field values from a or b
//- depending on the conditional.
//
//  Since boolean fields are not normally used, a flip operation is
//  a general requirement.
template
<
    class T, class BoolType, class FlipOp,
    template<class> class PatchField, class GeoMesh
>
void ternarySelect
(
    GeometricField<T, PatchField, GeoMesh>& result,
    const GeometricField<BoolType, PatchField, GeoMesh>& cond,
    const GeometricField<T, PatchField, GeoMesh>& a,
    const GeometricField<T, PatchField, GeoMesh>& b,
    const FlipOp& flip
)
{
    if constexpr
    (
        GeometricBoundaryField<T, PatchField, GeoMesh>::supportsDeferred
    )
    {
        if
        (
            !assignFlat
            (
                result,
                cond,
                a,
                b,
                [&flip](const auto& selector, const T& lhs, const T& rhs)
                {
                    return flip(selector) ? lhs : rhs;
                }
            )
        )
        {
            FieldOps::ternarySelect
            (
                result.primitiveFieldRef(),
                cond.primitiveField(),
                a.primitiveField(),
                b.primitiveField(),
                flip
            );

            auto& bfld = result.boundaryFieldRef();

            const label len = bfld.size();

            for (label i = 0; i < len; ++i)
            {
                FieldOps::ternarySelect
                (
                    bfld[i],
                    cond.boundaryField()[i],
                    a.boundaryField()[i],
                    b.boundaryField()[i],
                    flip
                );
            }
        }
    }
    else
    {
        FieldOps::ternarySelect
        (
            result.primitiveFieldRef(),
            cond.primitiveField(),
            a.primitiveField(),
            b.primitiveField(),
            flip
        );
    }
    result.correctLocalBoundaryConditions();
    if (GeometricBoundaryField<T, PatchField, GeoMesh>::debug)
    {
        result.boundaryField().check();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FieldOps
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
