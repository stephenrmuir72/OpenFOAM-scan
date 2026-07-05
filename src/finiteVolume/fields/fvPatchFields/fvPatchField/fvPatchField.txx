/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2025 OpenCFD Ltd.
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

#include "IOobject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "fvPatchFieldMapper.H"
#include "volMesh.H"


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
void Foam::fvPatchField<Type>::rebindToIF() noexcept
{
    const fvPatch& p = fvPatchFieldBase::patch();

    // Flat-storage (rebind) design: a patch field normally slices the
    // boundary extension of its DimensionedField so internal + boundary data
    // are contiguous.  Two cases cannot rebind and must own p.size() storage
    // so the patch field still behaves like a normal Field:
    //   1. Null DimensionedField (NewCalculatedType for arithmetic
    //      temporaries; slicedFvPatchField which then shallowCopies).
    //   2. A DimensionedField without a boundary extension (e.g. the
    //      SlicedDimensionedField dummy mesh.V() used by nearWallDist).
    if (&internalField_ == &DimensionedField<Type, volMesh>::null())
    {
        ownAllocation_ = true;
        static_cast<List<Type>&>(*this).resize_nocopy(p.size());
        return;
    }

    const label nCells = internalField_.mesh().nCells();
    const label offset = nCells + p.offset();

    const auto& df =
        static_cast<const DynamicField<Type>&>(internalField_);

    if
    (
        FieldBase::unifiedGeometricField
     && df.capacity() >= offset + p.size()
    )
    {
        // Rebind: point into the DimensionedField boundary extension.
        // ownAllocation_ stays false; destructor will null the pointer.
        static_cast<UList<Type>&>(*this) =
            UList<Type>
            (
                const_cast<Type*>(&internalField_[offset]),
                p.size()
            );
    }
    else
    {
        // No boundary extension available: own p.size() storage.
        ownAllocation_ = true;
        static_cast<List<Type>&>(*this).resize_nocopy(p.size());
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fvPatchField<Type>::releaseStorage() noexcept
{
    if (ownAllocation_)
    {
        // Free the own heap allocation
        static_cast<List<Type>&>(*this).clear();
        ownAllocation_ = false;
    }
    else
    {
        // Detach a rebind slice without freeing the DimensionedField extension
        static_cast<UList<Type>&>(*this) = UList<Type>();
    }
}


template<class Type>
bool Foam::fvPatchField<Type>::readValueEntry
(
    const dictionary& dict,
    IOobjectOption::readOption readOpt
)
{
    if (!IOobjectOption::isAnyRead(readOpt)) return false;
    const auto& p = fvPatchFieldBase::patch();


    const auto* eptr = dict.findEntry("value", keyType::LITERAL);

    if (eptr)
    {
        Field<Type>::assign(*eptr, p.size());
        return true;
    }

    if (IOobjectOption::isReadRequired(readOpt))
    {
        FatalIOErrorInFunction(dict)
            << "Required entry 'value' : missing for patch " << p.name()
            << " in dictionary " << dict.relativeName() << nl
            << exit(FatalIOError);
    }

    return false;
}


template<class Type>
void Foam::fvPatchField<Type>::extrapolateInternal()
{
    fvPatchFieldBase::patch().patchInternalField(internalField_, *this);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchFieldBase(p),
    Field<Type>(),
    internalField_(iF)
{
    rebindToIF();
}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const word& patchType
)
:
    fvPatchFieldBase(p, patchType),
    Field<Type>(),
    internalField_(iF)
{
    rebindToIF();
}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Type& value
)
:
    fvPatchFieldBase(p),
    Field<Type>(),
    internalField_(iF)
{
    rebindToIF();
    Field<Type>::operator=(value);
}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& pfld
)
:
    fvPatchFieldBase(p),
    Field<Type>(),
    internalField_(iF)
{
    rebindToIF();
    if (pfld.size()) Field<Type>::operator=(pfld);
}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    Field<Type>&& pfld
)
:
    fvPatchFieldBase(p),
    Field<Type>(),
    internalField_(iF)
{
    rebindToIF();
    if (pfld.size()) Field<Type>::operator=(pfld);
}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    IOobjectOption::readOption requireValue
)
:
    fvPatchFieldBase(p, dict),
    Field<Type>(),
    internalField_(iF)
{
    rebindToIF();
    readValueEntry(dict, requireValue);
}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchFieldBase(ptf, p),
    Field<Type>(),
    internalField_(iF)
{
    rebindToIF();
    // For unmapped faces set to internal field value (zero-gradient)
    if (notNull(iF) && mapper.hasUnmapped())
    {
        this->extrapolateInternal();
    }
    Field<Type>::map(ptf, mapper);
}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatchField<Type>& pfld,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Type& value
)
:
    fvPatchFieldBase(pfld, p),
    Field<Type>(),
    internalField_(iF)
{
    rebindToIF();
    Field<Type>::operator=(value);
}


template<class Type>
Foam::fvPatchField<Type>::fvPatchField
(
    const fvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchFieldBase(ptf),
    Field<Type>(),
    internalField_(iF)
{
    rebindToIF();
    Field<Type>::operator=(ptf);
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class Type>
Foam::fvPatchField<Type>::~fvPatchField() noexcept
{
    if (!ownAllocation_)
    {
        // Null the UList so List<Type> destructor does not attempt to free
        // memory it does not own.  Covers three cases:
        //   1. Rebind into DimensionedField boundary extension.
        //   2. Null internalField (SlicedGeometricField, slicedFvPatchField)
        //      where v_ was set externally via shallowCopy.
        //   3. rebindToIF() exited early (null reference).
        static_cast<UList<Type>&>(*this) = UList<Type>();
    }
    // else: own allocation from rebindToIF() fallback → let List<T> free it.
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvPatchField<Type>::check(const fvPatchField<Type>& rhs) const
{
    fvPatchFieldBase::checkPatch(rhs);
}


template<class Type>
void Foam::fvPatchField<Type>::snGrad(UList<Type>& result) const
{
    // Get patch internal field, store temporarily in result
    this->patchInternalField(result);
    const auto& pif = result;

    const Field<Type>& pfld = *this;
    const auto& dc = patch().deltaCoeffs();

    const label len = result.size();

    labelRange faceRange(len);

    apply_policy(static_cast<std::size_t>(len), [&](auto policy)
    {
        std::for_each
        (
            policy,
            faceRange.begin(),
            faceRange.end(),
            [&](const label i)
            {
                result[i] = dc[i]*(pfld[i] - pif[i]);
            }
        );
    });
}

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvPatchField<Type>::snGrad() const
{
    auto tfld = tmp<Field<Type>>::New(this->size());
    this->snGrad(static_cast<UList<Type>&>(tfld.ref()));
    return tfld;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fvPatchField<Type>::patchInternalField() const
{
    return patch().patchInternalField(internalField_);
}


template<class Type>
void Foam::fvPatchField<Type>::patchInternalField(UList<Type>& pfld) const
{
    patch().patchInternalField(internalField_, pfld);
}


template<class Type>
void Foam::fvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    if (!this->size() && !mapper.distributed())
    {
        if (this->size())
        {
            Field<Type>::operator=(this->patchInternalField()());
        }
    }
    else
    {
        Field<Type>::autoMap(mapper);

        // For unmapped faces set to internal field value (zero-gradient)
        if (mapper.hasUnmapped())
        {
            Field<Type> pif(this->patchInternalField());

            if
            (
                mapper.direct()
             && notNull(mapper.directAddressing())
             && mapper.directAddressing().size()
            )
            {
                const labelUList& mapAddressing = mapper.directAddressing();

                forAll(mapAddressing, i)
                {
                    if (mapAddressing[i] < 0)
                    {
                        Field<Type>::operator[](i) = pif[i];
                    }
                }
            }
            else if (!mapper.direct() && mapper.addressing().size())
            {
                const labelListList& mapAddressing = mapper.addressing();

                forAll(mapAddressing, i)
                {
                    const labelList& localAddrs = mapAddressing[i];

                    if (!localAddrs.size())
                    {
                        Field<Type>::operator[](i) = pif[i];
                    }
                }
            }
        }
    }
}


template<class Type>
void Foam::fvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    Field<Type>::rmap(ptf, addr);
}


template<class Type>
void Foam::fvPatchField<Type>::updateCoeffs()
{
    fvPatchFieldBase::setUpdated(true);
}


template<class Type>
void Foam::fvPatchField<Type>::updateWeightedCoeffs(const scalarField& weights)
{
    // Default behaviour ignores the weights
    if (!updated())
    {
        updateCoeffs();

        fvPatchFieldBase::setUpdated(true);
    }
}


template<class Type>
void Foam::fvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!updated())
    {
        updateCoeffs();
    }

    fvPatchFieldBase::setUpdated(false);
    fvPatchFieldBase::setManipulated(false);
}


template<class Type>
void Foam::fvPatchField<Type>::manipulateMatrix(fvMatrix<Type>& matrix)
{
    fvPatchFieldBase::setManipulated(true);
}


template<class Type>
void Foam::fvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& matrix,
    const scalarField& weights
)
{
    fvPatchFieldBase::setManipulated(true);
}


template<class Type>
void Foam::fvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& matrix,
    const label iMatrix,
    const direction cmp
)
{
    fvPatchFieldBase::setManipulated(true);
}


template<class Type>
void Foam::fvPatchField<Type>::write(Ostream& os) const
{
    os.writeEntry("type", type());

    if (!patchType().empty())
    {
        os.writeEntry("patchType", patchType());
    }
    if (useImplicit())
    {
        os.writeEntry("useImplicit", "true");
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    Field<Type>::operator=(ul);
}


template<class Type>
void Foam::fvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchFieldBase::checkPatch(ptf);
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator+=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchFieldBase::checkPatch(ptf);
    Field<Type>::operator+=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator-=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchFieldBase::checkPatch(ptf);
    Field<Type>::operator-=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator*=
(
    const fvPatchField<scalar>& ptf
)
{
    fvPatchFieldBase::checkPatch(ptf);
    Field<Type>::operator*=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator/=
(
    const fvPatchField<scalar>& ptf
)
{
    fvPatchFieldBase::checkPatch(ptf);
    Field<Type>::operator/=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator+=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator-=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    Field<Type>::operator*=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    Field<Type>::operator/=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


template<class Type>
void Foam::fvPatchField<Type>::operator+=
(
    const Type& t
)
{
    Field<Type>::operator+=(t);
}


template<class Type>
void Foam::fvPatchField<Type>::operator-=
(
    const Type& t
)
{
    Field<Type>::operator-=(t);
}


template<class Type>
void Foam::fvPatchField<Type>::operator*=
(
    const scalar s
)
{
    Field<Type>::operator*=(s);
}


template<class Type>
void Foam::fvPatchField<Type>::operator/=
(
    const scalar s
)
{
    Field<Type>::operator/=(s);
}


template<class Type>
void Foam::fvPatchField<Type>::operator==
(
    const fvPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<class Type>
void Foam::fvPatchField<Type>::operator==
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const fvPatchField<Type>& ptf)
{
    ptf.write(os);

    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
