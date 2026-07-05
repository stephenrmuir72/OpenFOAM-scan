/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "SlicedGeometricField.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //

template
<
    class Type,
    template<class> class PatchField,
    template<class> class SlicedPatchField,
    class GeoMesh
>
bool
Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
isBoundaryAddressing
(
    const Mesh& mesh,
    const label fieldSize
)
{
    label maxAddress(0);

    if (!mesh.boundary().empty())
    {
        const auto& p = mesh.boundary().back();
        maxAddress = (p.start() + p.size());
    }

    // If field size appear to not include internal field
    return (fieldSize < maxAddress);
}


template
<
    class Type,
    template<class> class PatchField,
    template<class> class SlicedPatchField,
    class GeoMesh
>
Foam::tmp<Foam::FieldField<PatchField, Type>>
Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
makeBoundary
(
    const Mesh& mesh,
    const Field<Type>& completeOrBoundaryField,
    const bool preserveCouples,
    const bool preserveProcessorOnly,
    const bool isBoundaryOnly
) const
{
    typedef typename
        SlicedPatchField<Type>::processorPatchType
        processorPatchType;

    auto tbf = tmp<FieldField<PatchField, Type>>::New(mesh.boundary().size());
    auto& bf = tbf.ref();

    forAll(mesh.boundary(), patchi)
    {
        const auto& p = mesh.boundary()[patchi];

        if
        (
            preserveCouples && p.coupled()
         && (!preserveProcessorOnly || isA<processorPatchType>(p))
        )
        {
            // For coupled patched construct the correct patch field type
            bf.set
            (
                patchi,
                PatchField<Type>::New(p.type(), p, *this)
            );

            // Initialize the values on the coupled patch to those of the slice
            // of the given field.
            // Note: these will usually be over-ridden by the boundary field
            // evaluation e.g. in the case of processor and cyclic patches.
            bf[patchi] = SlicedPatchField<Type>
            (
                p,
                DimensionedField<Type, GeoMesh>::null(),
                completeOrBoundaryField,
                isBoundaryOnly
            );
        }
        else
        {
            bf.set
            (
                patchi,
                new SlicedPatchField<Type>
                (
                    p,
                    // Bind to this (now sliced) field: parallel-safe because
                    // makeBoundary is called after the internal slice is set.
                    *this,
                    completeOrBoundaryField,
                    isBoundaryOnly
                )
            );
        }
    }

    return tbf;
}


template
<
    class Type,
    template<class> class PatchField,
    template<class> class SlicedPatchField,
    class GeoMesh
>
Foam::tmp<Foam::FieldField<PatchField, Type>>
Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
makeBoundary
(
    const Mesh& mesh,
    const FieldField<PatchField, Type>& bField,
    const bool preserveCouples
) const
{
    auto tbf = tmp<FieldField<PatchField, Type>>::New(mesh.boundary().size());
    auto& bf = tbf.ref();

    forAll(mesh.boundary(), patchi)
    {
        const auto& p = mesh.boundary()[patchi];

        if (preserveCouples && p.coupled())
        {
            // For coupled patched construct the correct patch field type
            bf.set
            (
                patchi,
                PatchField<Type>::New(p.type(), p, *this)
            );

            // Assign field
            bf[patchi] == bField[patchi];
        }
        else
        {
            // Create unallocated copy of patch field
            bf.set
            (
                patchi,
                new SlicedPatchField<Type>
                (
                    p,
                    // Bind to this (now sliced) field: parallel-safe because
                    // makeBoundary is called after the internal slice is set.
                    *this
                )
            );

            bf[patchi].UList<Type>::shallowCopy(bField[patchi]);
        }
    }

    return tbf;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    class Type,
    template<class> class PatchField,
    template<class> class SlicedPatchField,
    class GeoMesh
>
Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
SlicedGeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& dims,
    const Field<Type>& completeField,
    const bool preserveCouples
)
:
    GeometricField<Type, PatchField, GeoMesh>
    (
        io,
        mesh,
        dims,
        Field<Type>(),
        // Empty boundary: populated below once the internal slice is valid
        UPtrList<PatchField<Type>>()
    )
{
    // The internal field was given a boundary-extension capacity by its
    // constructor (for the flat-storage rebind design); a sliced field owns
    // no such storage.  Release it (capacity -> 0) before pointing at the
    // external data, then build the boundary against the valid slice so patch
    // fields (incl. parallel processor patches) bind correctly.
    static_cast<DynamicField<Type>&>(*this).clearStorage();

    // Set internalField to the slice of the complete field
    UList<Type>::shallowCopy
    (
        SubList<Type>(completeField, GeoMesh::size(mesh))
    );

    this->boundaryFieldRef().transfer
    (
        makeBoundary(mesh, completeField, preserveCouples).ref()
    );

    correctBoundaryConditions();
}


template
<
    class Type,
    template<class> class PatchField,
    template<class> class SlicedPatchField,
    class GeoMesh
>
Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
SlicedGeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& dims,
    const Field<Type>& completeIField,
    const Field<Type>& completeBField,
    const bool preserveCouples,
    const bool preserveProcessorOnly
)
:
    GeometricField<Type, PatchField, GeoMesh>
    (
        io,
        mesh,
        dims,
        Field<Type>(),
        // Empty boundary: populated below once the internal slice is valid
        UPtrList<PatchField<Type>>()
    )
{
    // Release the boundary-extension capacity (see other constructor) before
    // pointing at the external data, then build the boundary against the slice.
    static_cast<DynamicField<Type>&>(*this).clearStorage();

    // Set internalField to the slice of the complete field
    UList<Type>::shallowCopy
    (
        SubList<Type>(completeIField, GeoMesh::size(mesh))
    );

    this->boundaryFieldRef().transfer
    (
        makeBoundary
        (
            mesh,
            completeBField,
            preserveCouples,
            preserveProcessorOnly,
            isBoundaryAddressing(mesh, completeBField.size())
        ).ref()
    );

    correctBoundaryConditions();
}


template
<
    class Type,
    template<class> class PatchField,
    template<class> class SlicedPatchField,
    class GeoMesh
>
Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
SlicedGeometricField
(
    const IOobject& io,
    const GeometricField<Type, PatchField, GeoMesh>& gf,
    const bool preserveCouples
)
:
    GeometricField<Type, PatchField, GeoMesh>
    (
        io,
        gf.mesh(),
        gf.dimensions(),
        Field<Type>(),
        // Empty boundary: populated below once the internal slice is valid
        UPtrList<PatchField<Type>>()
    )
{
    // Release the boundary-extension capacity before pointing at external data
    static_cast<DynamicField<Type>&>(*this).clearStorage();

    // Set internalField to the internal field
    UList<Type>::shallowCopy(gf.primitiveField());

    this->boundaryFieldRef().transfer
    (
        makeBoundary(gf.mesh(), gf.boundaryField(), preserveCouples).ref()
    );

    correctBoundaryConditions();
}


template
<
    class Type,
    template<class> class PatchField,
    template<class> class SlicedPatchField,
    class GeoMesh
>
Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
SlicedGeometricField
(
    const SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>& gf
)
:
    GeometricField<Type, PatchField, GeoMesh>
    (
        gf,
        gf.mesh(),
        gf.dimensions(),
        Field<Type>(),
        // Empty boundary: populated below once the internal slice is valid
        UPtrList<PatchField<Type>>()
    )
{
    // Release the boundary-extension capacity before pointing at external data
    static_cast<DynamicField<Type>&>(*this).clearStorage();

    // Set internalField to the internal field
    UList<Type>::shallowCopy(gf.primitiveField());

    // preserveCouples = true
    this->boundaryFieldRef().transfer
    (
        makeBoundary(gf.mesh(), gf.boundaryField(), true).ref()
    );
}


template
<
    class Type,
    template<class> class PatchField,
    template<class> class SlicedPatchField,
    class GeoMesh
>
Foam::tmp
<
    Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>
>
Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
clone() const
{
    return tmp
    <
        SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>
    >::New
    (
        *this
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template
<
    class Type,
    template<class> class PatchField,
    template<class> class SlicedPatchField,
    class GeoMesh
>
Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
~SlicedGeometricField()
{
    // Set internalField to nullptr to avoid deletion of underlying field
    UList<Type>::shallowCopy(nullptr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    class Type,
    template<class> class PatchField,
    template<class> class SlicedPatchField,
    class GeoMesh
>
void Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
correctBoundaryConditions()
{
    GeometricField<Type, PatchField, GeoMesh>::correctBoundaryConditions();
}


// ************************************************************************* //
