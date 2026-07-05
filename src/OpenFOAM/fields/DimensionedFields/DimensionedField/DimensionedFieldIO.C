/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2025 OpenCFD Ltd.
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

#include "DimensionedField.H"
#include "IOstreams.H"
#include "localIOdictionary.H"
#include "SubList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, class GeoMesh>
void Foam::DimensionedField<Type, GeoMesh>::readField
(
    const dictionary& fieldDict,
    const word& fieldDictEntry
)
{
    dimensions_.readEntry("dimensions", fieldDict);

    // Note: oriented state may have already been set on construction
    // - if so - do not reset by re-reading
    // - required for backwards compatibility in case of restarting from
    //   an old run where the oriented state may not have been set
    if (oriented_.oriented() != orientedType::ORIENTED)
    {
        oriented_.read(fieldDict);  // The "oriented" entry (if present)
    }

    const label meshSize = GeoMesh::size(mesh_);
    const label totalSize = meshSize + GeoMesh::boundary_size(mesh_);

    // The primitive field
    auto& fld = static_cast<DynamicField<Type>&>(*this);

    // Resize primitive field to make space for the internal field.
    // When boundary patch fields already exist (rebind design), they point
    // into the extension [meshSize..totalSize-1].  Use reserve_exact +
    // setAddressableSize instead of clear() so the existing allocation (and
    // any live patch field pointers into it) are preserved.
    // NOTE: Field<Type>::assign() uses non-virtual List::resize_nocopy
    // internally, which would call clear() and free the extension.  Read into
    // a temporary Field and then deepCopy to avoid any List resize/free.
    if (FieldBase::unifiedGeometricField && fld.capacity() < totalSize)
    {
        // No patch fields yet (first/early call): safe to free and realloc.
        fld.clear();
        fld.reserve_exact(totalSize);
    }
    // Set addressable size to meshSize via DynamicField's own resize_nocopy,
    // which preserves capacity (unlike List::resize_nocopy which would free).
    fld.resize_nocopy(meshSize);

    {
        // Read the meshSize internal values into a temporary Field.
        // Field<Type>(entry, len) calls assign() on its OWN heap allocation,
        // so it never touches the DimensionedField's extension.
        const auto* eptr =
            fieldDict.findEntry(fieldDictEntry, keyType::LITERAL);

        if (!eptr)
        {
            FatalIOErrorInFunction(fieldDict)
                << "Required entry '" << fieldDictEntry
                << "' missing in dictionary "
                << fieldDict.relativeName() << nl
                << exit(FatalIOError);
        }

        Field<Type> tmp(*eptr, meshSize);

        // Copy meshSize values in-place; extension is untouched.
        UList<Type>(fld.data(), meshSize).deepCopy(tmp);
    }
}


template<class Type, class GeoMesh>
void Foam::DimensionedField<Type, GeoMesh>::readField
(
    const word& fieldDictEntry
)
{
    dictionary dict
    (
        localIOdictionary::readContents
        (
            IOobject
            (
                this->name(),
                this->instance(),
                this->local(),
                this->db(),
                IOobjectOption::MUST_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            ),
            typeName
        )
    );

    this->close();

    readField(dict, fieldDictEntry);
}


template<class Type, class GeoMesh>
bool Foam::DimensionedField<Type, GeoMesh>::readIfPresent
(
    const word& fieldDictEntry
)
{
    if
    (
        this->isReadRequired()
     || (this->isReadOptional() && this->headerOk())
    )
    {
        readField(fieldDictEntry);
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const word& fieldDictEntry,
    const bool extraCapacity
)
:
    regIOobject(io),
    mesh_(mesh)
{
    if (extraCapacity && FieldBase::unifiedGeometricField)
    {
        DynamicField<Type>::reserve_exact
        (
            GeoMesh::size(mesh) + GeoMesh::boundary_size(mesh)
        );
    }
    readField(fieldDictEntry);
}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const dictionary& fieldDict,
    const word& fieldDictEntry,
    const bool extraCapacity
)
:
    regIOobject(io),
    mesh_(mesh)
{
    if (extraCapacity && FieldBase::unifiedGeometricField)
    {
        DynamicField<Type>::reserve_exact
        (
            GeoMesh::size(mesh) + GeoMesh::boundary_size(mesh)
        );
    }
    readField(fieldDict, fieldDictEntry);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
bool Foam::DimensionedField<Type, GeoMesh>::writeData
(
    Ostream& os,
    const word& fieldDictEntry
) const
{
    os.writeEntry("dimensions", dimensions());
    os << nl;

    if (oriented_.writeEntry(os))
    {
        os << nl;
    }

    DynamicField<Type>::writeEntry(fieldDictEntry, os);

    os.check(FUNCTION_NAME);
    return os.good();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const DimensionedField<Type, GeoMesh>& fld
)
{
    fld.writeData(os);

    return os;
}


template<class Type, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const tmp<DimensionedField<Type, GeoMesh>>& tfld
)
{
    tfld().writeData(os);
    tfld.clear();

    return os;
}


// ************************************************************************* //
