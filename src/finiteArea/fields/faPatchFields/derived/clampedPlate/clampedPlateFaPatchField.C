/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>

\*---------------------------------------------------------------------------*/

#include "clampedPlateFaPatchField.H"
#include "areaFields.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
clampedPlateFaPatchField<Type>::clampedPlateFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(p, iF)
{}


template<class Type>
clampedPlateFaPatchField<Type>::clampedPlateFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    faPatchField<Type>(p, iF, dict, IOobjectOption::NO_READ)
{
    this->extrapolateInternal();  // Zero-gradient patch values
}


template<class Type>
clampedPlateFaPatchField<Type>::clampedPlateFaPatchField
(
    const clampedPlateFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    faPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
clampedPlateFaPatchField<Type>::clampedPlateFaPatchField
(
    const clampedPlateFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void clampedPlateFaPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    // Clamp the edge (patch value) and the adjacent internal cells to zero.
    // The operation is type-generic (Zero is valid for every Type), so the
    // single template serves scalar/vector/tensor plate fields alike.
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=(Zero);

    auto& iF = const_cast<Field<Type>&>(this->primitiveField());

    // faPatch::edgeFaces() maps to fvPatch::faceCells() in the
    // reduced-order (2D) fvMesh representation
    for (const label facei : this->patch().faceCells())
    {
        iF[facei] = Zero;
    }

    faPatchField<Type>::evaluate();
}


template<class Type>
tmp<Field<Type>> clampedPlateFaPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type>>::New(this->size(), Zero);
}


template<class Type>
tmp<Field<Type>> clampedPlateFaPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return *this;
}


template<class Type>
tmp<Field<Type>>
clampedPlateFaPatchField<Type>::gradientInternalCoeffs() const
{
    return -Type(pTraits<Type>::one)*this->patch().deltaCoeffs();
}


template<class Type>
tmp<Field<Type>>
clampedPlateFaPatchField<Type>::gradientBoundaryCoeffs() const
{
    return this->patch().deltaCoeffs()*(*this);
}


template<class Type>
void clampedPlateFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    faPatchField<Type>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
