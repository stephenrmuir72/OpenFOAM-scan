/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
    Copyright (C) Keysight Technologies
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

#include "meshPointPatchDisplacementVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "pointConstraint.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshPointPatchDisplacementVectorField::
meshPointPatchDisplacementVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    mpPatch_(refCast<const meshPointPatch>(p))
{}


Foam::meshPointPatchDisplacementVectorField::
meshPointPatchDisplacementVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    mpPatch_(refCast<const meshPointPatch>(p))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


Foam::meshPointPatchDisplacementVectorField::
meshPointPatchDisplacementVectorField
(
    const meshPointPatchDisplacementVectorField& ppf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ppf, p, iF, mapper),
    mpPatch_(refCast<const meshPointPatch>(p))
{}


// Foam::meshPointPatchDisplacementVectorField::
// meshPointPatchDisplacementVectorField
// (
//     const meshPointPatchDisplacementVectorField& ppf
// )
// :
//     fixedValuePointPatchField<vector>(ppf),
//     mpPatch_(refCast<const meshPointPatch>(ppf.patch()))
// {}


Foam::meshPointPatchDisplacementVectorField::
meshPointPatchDisplacementVectorField
(
    const meshPointPatchDisplacementVectorField& ppf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ppf, iF),
    mpPatch_(refCast<const meshPointPatch>(ppf.patch()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshPointPatchDisplacementVectorField::updateCoeffs()
{
    if (debug)
    {
        PoutInFunction<< " applying constraints for patch " << patch().name()
            << endl;
    }
    vectorField displacement(this->patchInternalField());

    // Apply constrained value

    // Constraints from patch
    const auto& pcs = mpPatch_.constraints();
    // Offset (e.g. position on surface) from valuePointPatchField
    const auto& wanted = *this;

    forAll(pcs, i)
    {
        const vector diff(displacement[i] - wanted[i]);
        displacement[i] = wanted[i] + pcs[i].constrainDisplacement(diff);
    }

    // Insert calculated field into internal field
    auto& iF = const_cast<Field<vector>&>(this->primitiveField());
    this->setInInternalField(iF, displacement);

    // Bypass fixedValue (since inserts *this)
    pointPatchField<vector>::updateCoeffs();
}


void Foam::meshPointPatchDisplacementVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!updated())
    {
        updateCoeffs();
    }

    pointPatchField<vector>::evaluate(commsType);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePointPatchTypeField
(
    pointPatchVectorField,
    meshPointPatchDisplacementVectorField
);

} // End namespace Foam

// ************************************************************************* //
