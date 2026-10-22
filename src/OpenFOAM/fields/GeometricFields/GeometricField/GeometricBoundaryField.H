/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017,2022 OpenFOAM Foundation
    Copyright (C) 2015-2026 OpenCFD Ltd.
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

Class
    Foam::GeometricBoundaryField

Description
    Generic GeometricBoundaryField class.

SourceFiles
    GeometricBoundaryField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_GeometricBoundaryField_H
#define Foam_GeometricBoundaryField_H

#include "dimensionedTypes.H"
#include "DimensionedField.H"
#include "FieldField.H"
#include "lduInterfaceFieldPtrsList.H"
#include "LduInterfaceFieldPtrsList.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointMesh;
template<class Type> class fvPatchField;

/*---------------------------------------------------------------------------*\
                   Class GeometricBoundaryField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricBoundaryField
:
    public FieldField<PatchField, Type>
{
public:

    // Public Typedefs

        //- The template instance used for this GeometricBoundaryField type
        typedef GeometricBoundaryField<Type, PatchField, GeoMesh> this_type;

        //- The boundary mesh type for the boundary fields
        typedef typename GeoMesh::BoundaryMesh BoundaryMesh;

        //- The internal field type associated with the boundary fields
        typedef DimensionedField<Type, GeoMesh> Internal;

        //- This boundary field type
        typedef GeometricBoundaryField<Type, PatchField, GeoMesh> Boundary;

        //- The patch field type for the boundary fields
        typedef PatchField<Type> Patch;

        //- Whether calculated patches may use flat internal storage
        static constexpr bool supportsDeferred =
            !std::is_same_v<GeoMesh, pointMesh>;


private:

    // Private Data

        //- Reference to BoundaryMesh for which this field is defined
        const BoundaryMesh& bmesh_;

        //- Internal field used to materialise deferred calculated patches
        const Internal* deferredInternal_ = nullptr;

        //- Patch-field type to use when materialising deferred patches
        word deferredPatchFieldType_ = PatchField<Type>::calculatedType();

            //- Per-patch types when deferred patches were specified as a list
            wordList deferredPatchFieldTypes_;

            //- Return the patch-field type to use when materialising a patch
            const word& deferredPatchFieldType(const label patchi) const
            {
                return deferredPatchFieldTypes_.empty()
                ? deferredPatchFieldType_
                : deferredPatchFieldTypes_[patchi];
            }

        //- Return the flat-storage view for a deferred patch
        UList<Type> deferredValues(const label patchi) const
        {
            if (!deferredInternal_)
            {
                FatalErrorInFunction
                    << "Cannot access an unset patch field at index "
                    << patchi << abort(FatalError);
            }

            const label offset =
                GeoMesh::size(deferredInternal_->mesh())
              + bmesh_[patchi].offset();

            const auto& values =
                static_cast<const DynamicField<Type>&>(*deferredInternal_);

            return UList<Type>
            (
                const_cast<Type*>(&values[offset]),
                bmesh_[patchi].size()
            );
        }

        //- Refresh deferred extrapolated-calculated values from owner cells
        void updateDeferredExtrapolated() const;


    // Private Member Functions

        //- Check if patch fields have been consistently evaluated.
        //  If not:
        //      exit = true  : FatalError
        //      exit = false : return bool
        template<class CheckPatchField>
        bool checkConsistency(const scalar tol, const bool exitIfBad) const;


public:

    //- Enable debug
    static int debug;


    // Constructors

        //- Construct from a BoundaryMesh, setting patches later
        explicit GeometricBoundaryField(const BoundaryMesh& bmesh);

        //- Construct from a BoundaryMesh, reference to the internal field
        //- and a patch type
        GeometricBoundaryField
        (
            const BoundaryMesh& bmesh,
            const Internal& iField,
            const word& patchFieldType = PatchField<Type>::calculatedType(),
            const bool deferCalculatedPatches = true
        );

        //- Construct from a BoundaryMesh, reference to the internal field
        //- and a wordList of patch types and optional the actual patch
        //- types (to override constraint patches)
        GeometricBoundaryField
        (
            const BoundaryMesh& bmesh,
            const Internal& iField,
            const wordList& wantedPatchTypes,
            const wordList& actualPatchTypes = wordList()
        );

        //- Construct from a BoundaryMesh, reference to the internal field
        //- and list of patch fields (to be cloned)
        GeometricBoundaryField
        (
            const BoundaryMesh& bmesh,
            const Internal& iField,
            //! Patch fields to be cloned.
            //! [For advanced usage: can be empty == set patches later]
            const UPtrList<PatchField<Type>>& pflds
        );

        //- Construct as copy, setting the reference to the internal field
        GeometricBoundaryField
        (
            const Internal& iField,
            const GeometricBoundaryField& btf
        );

        //- Construct as copy, setting the reference to the internal field
        //- and resetting type of field for given patch IDs
        GeometricBoundaryField
        (
            const Internal& iField,
            const GeometricBoundaryField& btf,
            const labelUList& patchIDs,
            const word& patchFieldName
        );

        //- Construct from dictionary
        GeometricBoundaryField
        (
            const BoundaryMesh& bmesh,
            const Internal& iField,
            const dictionary& dict
        );

        //- No copy construct
        //  Dangerous because it may be set to a field which gets deleted.
        //  Need new type of BoundaryField, one which is part of a geometric
        //  field for which snGrad etc. may be called and a free standing
        //  BoundaryField for which such operations are unavailable.
        GeometricBoundaryField(const GeometricBoundaryField& bf) = delete;
        // :
        //     FieldField<PatchField, Type>(bf),
        //     bmesh_(bf.bmesh_)
        // {}

        //- No move construct: likely just as dangerous as copy construct
        GeometricBoundaryField(GeometricBoundaryField&&) = delete;

        // FUTURE?
        // //- Clone, setting reference to the internal field
        // tmp<GeometricBoundaryField> clone(const Internal& iField) const
        // {
        //     return tmp<GeometricBoundaryField>::New(iField, *this);
        // }


    // Member Functions

        //- Materialise and return the patch field at the given index
        Patch& patchField(const label patchi);

        //- Materialise and return the patch field at the given index
        const Patch& patchField(const label patchi) const
        {
            return constCast().patchField(patchi);
        }

        //- True when at least one null slot denotes a deferred calculated patch
        bool hasDeferred() const noexcept
        {
            if (!deferredInternal_)
            {
                return false;
            }

            forAll(*this, patchi)
            {
                if (!this->get(patchi)) return true;
            }

            return false;
        }

        //- Values of a patch, using its flat-storage slice when deferred
        UList<Type> patchValues(const label patchi)
        {
            if constexpr (supportsDeferred)
            {
                if (Patch* pfld = this->get(patchi))
                {
                    return UList<Type>(pfld->begin(), pfld->size());
                }

                return deferredValues(patchi);
            }
            else
            {
                FatalErrorInFunction
                    << "Flat patch values are unavailable for point fields."
                    << " Trying to access values of patch " << patchi
                    << " of field "
                    << this->patchField(patchi).internalField().name()
                    << abort(FatalError);
                return {};
            }
        }

        //- Values of a patch, using its flat-storage slice when deferred
        UList<Type> patchValues(const label patchi) const
        {
            if constexpr (supportsDeferred)
            {
                if (const Patch* pfld = this->get(patchi))
                {
                    return UList<Type>
                    (
                        const_cast<Type*>(pfld->begin()),
                        pfld->size()
                    );
                }

                return deferredValues(patchi);
            }
            else
            {
                FatalErrorInFunction
                    << "Flat patch values are unavailable for point fields."
                    << " Trying to access values of patch " << patchi
                    << " of field "
                    << this->patchField(patchi).internalField().name()
                    << abort(FatalError);
                return {};
            }
        }

        //- Apply an operation to corresponding patch values
        template<class Type1, class BinaryOp>
        void apply
        (
            const GeometricBoundaryField<Type1, PatchField, GeoMesh>& rhs,
            const BinaryOp& op
        )
        {
            forAll(*this, patchi)
            {
                auto resultValues = patchValues(patchi);
                const auto rhsValues = rhs.patchValues(patchi);

                #ifdef FOAM_OFFLOAD
                apply_policy
                (
                    static_cast<std::size_t>(resultValues.size()),
                    [&](auto policy)
                    {
                        std::transform
                        (
                            policy,
                            resultValues.begin(),
                            resultValues.end(),
                            rhsValues.begin(),
                            resultValues.begin(),
                            [&](const auto& resultValue, const auto& rhsValue)
                            {
                                return op(resultValue, rhsValue);
                            }
                        );
                    }
                );
                #else
                forAll(resultValues, i)
                {
                    resultValues[i] = op(resultValues[i], rhsValues[i]);
                }
                #endif
            }
        }

        //- Apply an operation to corresponding values from two fields
        template<class Type1, class Type2, class BinaryOp>
        void apply
        (
            const GeometricBoundaryField<Type1, PatchField, GeoMesh>& lhs,
            const GeometricBoundaryField<Type2, PatchField, GeoMesh>& rhs,
            const BinaryOp& op
        )
        {
            forAll(*this, patchi)
            {
                auto resultValues = patchValues(patchi);
                const auto lhsValues = lhs.patchValues(patchi);
                const auto rhsValues = rhs.patchValues(patchi);

                #ifdef FOAM_OFFLOAD
                apply_policy
                (
                    static_cast<std::size_t>(resultValues.size()),
                    [&](auto policy)
                    {
                        std::transform
                        (
                            policy,
                            lhsValues.begin(),
                            lhsValues.end(),
                            rhsValues.begin(),
                            resultValues.begin(),
                            [&](const auto& lhsValue, const auto& rhsValue)
                            {
                                return op(lhsValue, rhsValue);
                            }
                        );
                    }
                );
                #else
                forAll(resultValues, i)
                {
                    resultValues[i] = op(lhsValues[i], rhsValues[i]);
                }
                #endif
            }
        }

        //- Apply an operation to a field and uniform value
        template<class Type1, class BinaryOp>
        void apply
        (
            const Type1& lhs,
            const GeometricBoundaryField<Type1, PatchField, GeoMesh>& rhs,
            const BinaryOp& op
        )
        {
            forAll(*this, patchi)
            {
                auto resultValues = patchValues(patchi);
                const auto rhsValues = rhs.patchValues(patchi);

                #ifdef FOAM_OFFLOAD
                apply_policy
                (
                    static_cast<std::size_t>(resultValues.size()),
                    [&](auto policy)
                    {
                        std::transform
                        (
                            policy,
                            rhsValues.begin(),
                            rhsValues.end(),
                            resultValues.begin(),
                            [&](const auto& rhsValue)
                            {
                                return op(lhs, rhsValue);
                            }
                        );
                    }
                );
                #else
                forAll(resultValues, i)
                {
                    resultValues[i] = op(lhs, rhsValues[i]);
                }
                #endif
            }
        }

        //- Apply an operation to a field and uniform value
        template<class Type1, class BinaryOp>
        void apply
        (
            const GeometricBoundaryField<Type1, PatchField, GeoMesh>& lhs,
            const Type1& rhs,
            const BinaryOp& op
        )
        {
            forAll(*this, patchi)
            {
                auto resultValues = patchValues(patchi);
                const auto lhsValues = lhs.patchValues(patchi);

                #ifdef FOAM_OFFLOAD
                apply_policy
                (
                    static_cast<std::size_t>(resultValues.size()),
                    [&](auto policy)
                    {
                        std::transform
                        (
                            policy,
                            lhsValues.begin(),
                            lhsValues.end(),
                            resultValues.begin(),
                            [&](const auto& lhsValue)
                            {
                                return op(lhsValue, rhs);
                            }
                        );
                    }
                );
                #else
                forAll(resultValues, i)
                {
                    resultValues[i] = op(lhsValues[i], rhs);
                }
                #endif
            }
        }

        //- Transform corresponding patch values without materialising fields
        template<class Type1, class UnaryOp>
        void transform
        (
            const GeometricBoundaryField<Type1, PatchField, GeoMesh>& rhs,
            const UnaryOp& op
        )
        {
            forAll(*this, patchi)
            {
                auto resultValues = patchValues(patchi);
                const auto rhsValues = rhs.patchValues(patchi);

                #ifdef FOAM_OFFLOAD
                apply_policy
                (
                    static_cast<std::size_t>(resultValues.size()),
                    [&](auto policy)
                    {
                        std::transform
                        (
                            policy,
                            rhsValues.begin(),
                            rhsValues.end(),
                            resultValues.begin(),
                            [&](const auto& rhsValue)
                            {
                                return op(rhsValue);
                            }
                        );
                    }
                );
                #else
                forAll(resultValues, i)
                {
                    resultValues[i] = op(rhsValues[i]);
                }
                #endif
            }
        }

        template<class Type1, class Type2, class BinaryOp>
        void apply
        (
            const GeometricBoundaryField<Type1, PatchField, GeoMesh>& rhs,
            const Type2& value,
            const BinaryOp& op
        )
        {
            forAll(*this, patchi)
            {
                auto resultValues = patchValues(patchi);
                const auto rhsValues = rhs.patchValues(patchi);

                #ifdef FOAM_OFFLOAD
                apply_policy
                (
                    static_cast<std::size_t>(resultValues.size()),
                    [&](auto policy)
                    {
                        std::transform
                        (
                            policy,
                            rhsValues.begin(),
                            rhsValues.end(),
                            resultValues.begin(),
                            [&](const auto& rhsValue)
                            {
                                return op(rhsValue, value);
                            }
                        );
                    }
                );
                #else
                forAll(resultValues, i)
                {
                    resultValues[i] = op(rhsValues[i], value);
                }
                #endif
            }
        }

        template<class Type1, class Type2, class BinaryOp>
        void apply
        (
            const Type1& value,
            const GeometricBoundaryField<Type2, PatchField, GeoMesh>& rhs,
            const BinaryOp& op
        )
        {
            forAll(*this, patchi)
            {
                auto resultValues = patchValues(patchi);
                const auto rhsValues = rhs.patchValues(patchi);

                #ifdef FOAM_OFFLOAD
                apply_policy
                (
                    static_cast<std::size_t>(resultValues.size()),
                    [&](auto policy)
                    {
                        std::transform
                        (
                            policy,
                            rhsValues.begin(),
                            rhsValues.end(),
                            resultValues.begin(),
                            [&](const auto& rhsValue)
                            {
                                return op(value, rhsValue);
                            }
                        );
                    }
                );
                #else
                forAll(resultValues, i)
                {
                    resultValues[i] = op(value, rhsValues[i]);
                }
                #endif
            }
        }

        //- Return non-const reference to this boundary field
        this_type& constCast() const noexcept
        {
            return const_cast<this_type&>(*this);
        }

        //- Read the boundary field
        void readField(const Internal& iField, const dictionary& dict);

        //- Update the boundary condition coefficients
        void updateCoeffs();

        //- Evaluate boundary conditions for each patch field.
        //- Uses specified or default comms.
        void evaluate
        (
            const UPstream::commsTypes commsType = UPstream::defaultCommsType
        );

        //- Evaluate boundary conditions for patch fields matching the
        //- given predicate. Uses specified or default comms.
        template<class UnaryPredicate>
        void evaluate_if
        (
            const UnaryPredicate& pred,
            const UPstream::commsTypes commsType = UPstream::defaultCommsType
        );

        //- Evaluate boundary conditions after change in local values.
        //- Uses specified or default comms.
        void evaluateLocal
        (
            const UPstream::commsTypes commsType = UPstream::defaultCommsType
        );

        //- Evaluate boundary conditions on coupled patches of the given type,
        //- using specified or default comms.
        //
        //  \tparam CoupledPatchType The coupled \b patch type that should
        //     be evaluated. For a \c void type,
        //     all coupled patch fields will be evaluated.
        template<class CoupledPatchType>
        void evaluateCoupled
        (
            const UPstream::commsTypes commsType = UPstream::defaultCommsType
        );

        //- Evaluate boundary conditions for selected patches.
        //- Uses non-blocking comms.
        void evaluateSelected(const labelUList& patchIDs);

        //- Return a list of the patch types
        wordList types() const;

        //- Return boundary field of values neighbouring the boundary
        tmp<GeometricBoundaryField> boundaryInternalField() const;

        //- Return a list of pointers for each patch field with only those
        //- pointing to interfaces being set
        LduInterfaceFieldPtrsList<Type> interfaces() const;

        //- Return a list of pointers for each patch field with only those
        //- pointing to interfaces being set
        lduInterfaceFieldPtrsList scalarInterfaces() const;

        //- Helper: check if field has been evaluated. See instantiations.
        bool check() const;


    // Write

        //- Write boundary field as dictionary entry
        void writeEntry(const word& keyword, Ostream& os) const;

        //- Write dictionary entries of the individual boundary fields.
        void writeEntries(Ostream& os) const;


    // Member Operators

        //- Preserve FieldField access to component indices
        using FieldField<PatchField, Type>::operator[];

        //- Access a patch field, materialising it if deferred
        Patch& operator[](const label patchi)
        {
            return patchField(patchi);
        }

        //- Access a patch field, materialising it if deferred
        const Patch& operator[](const label patchi) const
        {
            return patchField(patchi);
        }

        //- Copy assignment from GeometricBoundaryField
        void operator=(const GeometricBoundaryField& bf);

        //- Copy assignment from FieldField\<PatchField, Type\>
        void operator=(const FieldField<PatchField, Type>& bf);

        //- Assignment to uniform value
        void operator=(const Type& val);

        //- Assignment to zero
        void operator=(Foam::zero);

        //- Forced assignment from GeometricBoundaryField
        void operator==(const GeometricBoundaryField& bf);

        //- Forced assignment from FieldField\<PatchField, Type\>
        void operator==(const FieldField<PatchField, Type>& bf);

        //- Forced assignment to uniform value
        void operator==(const Type& val);

        //- Forced assignment to zero
        void operator==(Foam::zero) { *this == Type(Foam::zero{}); }

        // Prevent automatic comparison rewriting (c++20)
        bool operator!=(const GeometricBoundaryField&) = delete;
        bool operator!=(const FieldField<PatchField, Type>&) = delete;
        bool operator!=(const Type&) = delete;
        bool operator!=(Foam::zero) = delete;
};


template<class Type, template<class> class PatchField, class GeoMesh>
Ostream& operator<<
(
    Ostream&,
    const GeometricBoundaryField<Type, PatchField, GeoMesh>&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "GeometricBoundaryField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
