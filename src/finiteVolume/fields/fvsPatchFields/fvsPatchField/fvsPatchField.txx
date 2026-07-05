/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
#include "surfaceMesh.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fvsPatchField<Type>::rebindToIF() noexcept
{
    const fvPatch& p = fvsPatchFieldBase::patch();

    // Flat-storage (rebind) design: cannot rebind for a null DimensionedField
    // (NewCalculatedType arithmetic temporaries, slicedFvsPatchField) or one
    // without a boundary extension; own p.size() storage in those cases.
    if (&internalField_ == &DimensionedField<Type, surfaceMesh>::null())
    {
        ownAllocation_ = true;
        static_cast<List<Type>&>(*this).resize_nocopy(p.size());
        return;
    }

    const label nInternal = internalField_.mesh().nInternalFaces();
    const label offset = nInternal + p.offset();

    const auto& df =
        static_cast<const DynamicField<Type>&>(internalField_);

    if
    (
        FieldBase::unifiedGeometricField
     && df.capacity() >= offset + p.size()
    )
    {
        static_cast<UList<Type>&>(*this) =
            UList<Type>
            (
                const_cast<Type*>(&internalField_[offset]),
                p.size()
            );
    }
    else
    {
        ownAllocation_ = true;
        static_cast<List<Type>&>(*this).resize_nocopy(p.size());
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fvsPatchField<Type>::releaseStorage() noexcept
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
bool Foam::fvsPatchField<Type>::readValueEntry
(
    const dictionary& dict,
    IOobjectOption::readOption readOpt
)
{
    if (!IOobjectOption::isAnyRead(readOpt)) return false;
    const auto& p = fvsPatchFieldBase::patch();


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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fvsPatchField<Type>::fvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    fvsPatchFieldBase(p),
    Field<Type>(),
    internalField_(iF)
{
    rebindToIF();
}


template<class Type>
Foam::fvsPatchField<Type>::fvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const Type& value
)
:
    fvsPatchFieldBase(p),
    Field<Type>(),
    internalField_(iF)
{
    rebindToIF();
    Field<Type>::operator=(value);
}


template<class Type>
Foam::fvsPatchField<Type>::fvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const Field<Type>& pfld
)
:
    fvsPatchFieldBase(p),
    Field<Type>(),
    internalField_(iF)
{
    rebindToIF();
    if (pfld.size()) Field<Type>::operator=(pfld);
}


template<class Type>
Foam::fvsPatchField<Type>::fvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    Field<Type>&& pfld
)
:
    fvsPatchFieldBase(p),
    Field<Type>(),
    internalField_(iF)
{
    rebindToIF();
    if (pfld.size()) Field<Type>::operator=(pfld);
}


template<class Type>
Foam::fvsPatchField<Type>::fvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict,
    IOobjectOption::readOption requireValue
)
:
    fvsPatchFieldBase(p, dict),
    Field<Type>(),
    internalField_(iF)
{
    rebindToIF();
    readValueEntry(dict, requireValue);
}


template<class Type>
Foam::fvsPatchField<Type>::fvsPatchField
(
    const fvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvsPatchFieldBase(ptf, p),
    Field<Type>(),
    internalField_(iF)
{
    rebindToIF();
    Field<Type>::map(ptf, mapper);
}


template<class Type>
Foam::fvsPatchField<Type>::fvsPatchField
(
    const fvsPatchField<Type>& pfld,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const Type& value
)
:
    fvsPatchFieldBase(pfld, p),
    Field<Type>(),
    internalField_(iF)
{
    rebindToIF();
    Field<Type>::operator=(value);
}


template<class Type>
Foam::fvsPatchField<Type>::fvsPatchField
(
    const fvsPatchField<Type>& pfld,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    fvsPatchFieldBase(pfld),
    Field<Type>(),
    internalField_(iF)
{
    rebindToIF();
    Field<Type>::operator=(pfld);
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class Type>
Foam::fvsPatchField<Type>::~fvsPatchField() noexcept
{
    if (!ownAllocation_)
    {
        static_cast<UList<Type>&>(*this) = UList<Type>();
    }
    // else: own allocation from rebindToIF() fallback → let List<T> free it.
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvsPatchField<Type>::check(const fvsPatchField<Type>& rhs) const
{
    fvsPatchFieldBase::checkPatch(rhs);
}


template<class Type>
void Foam::fvsPatchField<Type>::autoMap(const fvPatchFieldMapper& m)
{
    Field<Type>::autoMap(m, internalField_.is_oriented());
}


template<class Type>
void Foam::fvsPatchField<Type>::rmap
(
    const fvsPatchField<Type>& ptf,
    const labelList& addr
)
{
    Field<Type>::rmap(ptf, addr);
}


template<class Type>
void Foam::fvsPatchField<Type>::write(Ostream& os) const
{
    os.writeEntry("type", type());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvsPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    Field<Type>::operator=(ul);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator=
(
    const fvsPatchField<Type>& ptf
)
{
    fvsPatchFieldBase::checkPatch(ptf);
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator+=
(
    const fvsPatchField<Type>& ptf
)
{
    fvsPatchFieldBase::checkPatch(ptf);
    Field<Type>::operator+=(ptf);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator-=
(
    const fvsPatchField<Type>& ptf
)
{
    fvsPatchFieldBase::checkPatch(ptf);
    Field<Type>::operator-=(ptf);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator*=
(
    const fvsPatchField<scalar>& ptf
)
{
    fvsPatchFieldBase::checkPatch(ptf);
    Field<Type>::operator*=(ptf);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator/=
(
    const fvsPatchField<scalar>& ptf
)
{
    fvsPatchFieldBase::checkPatch(ptf);
    Field<Type>::operator/=(ptf);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator+=(tf);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator-=(tf);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    Field<Type>::operator*=(tf);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    Field<Type>::operator/=(tf);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator+=
(
    const Type& t
)
{
    Field<Type>::operator+=(t);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator-=
(
    const Type& t
)
{
    Field<Type>::operator-=(t);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator*=
(
    const scalar s
)
{
    Field<Type>::operator*=(s);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator/=
(
    const scalar s
)
{
    Field<Type>::operator/=(s);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator==
(
    const fvsPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<class Type>
void Foam::fvsPatchField<Type>::operator==
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const fvsPatchField<Type>& ptf)
{
    ptf.write(os);

    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
