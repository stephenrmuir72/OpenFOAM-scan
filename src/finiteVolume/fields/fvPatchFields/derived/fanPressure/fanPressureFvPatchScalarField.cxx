/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd
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

#include "fanPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "TableFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::fanPressureFvPatchScalarField::fanFlowDirection
>
Foam::fanPressureFvPatchScalarField::fanFlowDirectionNames_
({
    { fanFlowDirection::ffdIn, "in" },
    { fanFlowDirection::ffdOut, "out" },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(p, iF),
    fanCurve_(nullptr),
    nonDimensional_(false),
    fanNonDimensional_
    (
        p,
        iF.dimensions(),
        iF.name(),
        iF.objectPath()
    ),
    relaxFactor_(-1),
    pdFan0_(0),
    timeIndex_(-1),
    direction_(ffdOut)
{}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const this_bctype& rhs,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper),
    fanCurve_(rhs.fanCurve_.clone()),
    nonDimensional_(rhs.nonDimensional_),
    fanNonDimensional_
    (
        rhs.fanNonDimensional_,
        p,
        iF.dimensions(),
        iF.name(),
        iF.objectPath()
    ),
    relaxFactor_(rhs.relaxFactor_),
    pdFan0_(rhs.pdFan0_),
    timeIndex_(rhs.timeIndex_),
    direction_(rhs.direction_)
{}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict),
    fanCurve_(nullptr),
    nonDimensional_(dict.getOrDefault("nonDimensional", false)),
    fanNonDimensional_
    (
        p,
        iF.dimensions(),
        iF.name(),
        iF.objectPath(),
        dict,
        this->db(),
        nonDimensional_
    ),
    relaxFactor_(-1),
    pdFan0_(0),
    timeIndex_(this->db().time().timeIndex()),
    direction_(fanFlowDirectionNames_.get("direction", dict))
{
    // Validate the under-relaxation factor. The value is optional and
    // defaults to -1 (disabled); when specified it must satisfy 0 < relax <= 1.
    dict.readCheckIfPresent
    (
        "relax",
        relaxFactor_,
        [](const scalar v)
        {
            return v > 0 && v <= 1;
        }
    );

    // Backwards compatibility
    if (dict.found("file"))
    {
        fanCurve_.reset
        (
            new Function1Types::TableFile<scalar>("fanCurve", dict, &this->db())
        );
    }
    else
    {
        fanCurve_.reset(Function1<scalar>::New("fanCurve", dict, &this->db()));
    }

    // Read optional saved relaxation state for restart
    if (dict.found("deltaP0"))
    {
        pdFan0_ = dict.get<scalar>("deltaP0");
    }
}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const this_bctype& rhs,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(rhs, iF),
    fanCurve_(rhs.fanCurve_.clone()),
    nonDimensional_(rhs.nonDimensional_),
    fanNonDimensional_
    (
        rhs.fanNonDimensional_,
        this->patch(),
        iF.dimensions(),
        iF.name(),
        iF.objectPath()
    ),
    relaxFactor_(rhs.relaxFactor_),
    pdFan0_(rhs.pdFan0_),
    timeIndex_(rhs.timeIndex_),
    direction_(rhs.direction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fanPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Retrieve flux field
    const auto& phip = patch().lookupPatchField<surfaceScalarField>(phiName());

    const int dir = 2*direction_ - 1;

    // Average volumetric flow rate
    scalar volFlowRate = 0;

    if (phip.internalField().dimensions() == dimVolume/dimTime)
    {
        volFlowRate = dir*gSum(phip);
    }
    else if (phip.internalField().dimensions() == dimMass/dimTime)
    {
        const auto& rhop = patch().lookupPatchField<volScalarField>(rhoName());
        volFlowRate = dir*gSum(phip/rhop);
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct\n"
            << "    on patch " << patch().name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath() << nl
            << exit(FatalError);
    }

    // Preserve the physical volume flow rate for reporting
    const scalar Q = volFlowRate;

    // Pressure drop for this flow rate
    scalar pdFan;

    if (nonDimensional_)
    {
        const scalar timeValue =
            this->db().time().timeOutputValue();

        volFlowRate =
            fanNonDimensional_.nonDimensionalFlowRate
            (
                volFlowRate,
                timeValue
            );

        pdFan = fanCurve_->value(max(volFlowRate, scalar(0)));

        pdFan =
            fanNonDimensional_.dimensionalPressure
            (
                pdFan,
                rhoName(),
                timeValue
            );
    }
    else
    {
        pdFan = fanCurve_->value(max(volFlowRate, scalar(0)));
    }

    if (relaxFactor_ > 0)
    {
        pdFan = lerp(pdFan0_, pdFan, relaxFactor_);

        if (timeIndex_ != this->db().time().timeIndex())
        {
            pdFan0_ = pdFan;
            timeIndex_ = this->db().time().timeIndex();
        }
    }

    if (debug)
    {
        Info<< type() << " on patch " << patch().name()
            << " of field " << internalField().name() << ':' << nl
            << "    volume flow rate = " << Q << " m3/s"
            << ", deltaP = " << pdFan;

        if (internalField().dimensions() == dimPressure)
        {
            Info<< " Pa" << endl;
        }
        else
        {
            Info<< " Pa/(kg/m3)" << endl;
        }
    }

    this->parent_bctype::updateCoeffs
    (
        p0() - dir*pdFan,
        patch().lookupPatchField<volVectorField>(UName())
    );
}


void Foam::fanPressureFvPatchScalarField::write(Ostream& os) const
{
    this->parent_bctype::write(os);
    fanCurve_->writeData(os);
    os.writeEntry("direction", fanFlowDirectionNames_[direction_]);

    if (relaxFactor_ > 0)
    {
        os.writeEntry("relax", relaxFactor_);
        os.writeEntry("deltaP0", pdFan0_);
    }

    if (nonDimensional_)
    {
        os.writeEntry("nonDimensional", "true");
        fanNonDimensional_.write(os);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fanPressureFvPatchScalarField
    );
};


// ************************************************************************* //
