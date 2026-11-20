/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Keysight Technologies
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

#include "turbulentEpsilonOmegaInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * //

namespace Foam
{

//- Return true for epsilon dimensions, false for omega dimensions
static bool isEpsilonField(const DimensionedField<scalar, volMesh>& iF)
{
    if (iF.dimensions() == dimensionSet(0, 2, -3, 0, 0, 0, 0))
    {
        return true;
    }

    if (iF.dimensions() == dimensionSet(0, 0, -1, 0, 0, 0, 0))
    {
        return false;
    }

    FatalErrorInFunction
        << "Field " << iF.name() << " has dimensions " << iF.dimensions()
        << " but turbulentEpsilonOmegaInlet requires epsilon [m2/s3] or "
        << "omega [1/s]" << exit(FatalError);

    return false;
}

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{

turbulentEpsilonOmegaInletFvPatchScalarField::
turbulentEpsilonOmegaInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(p, iF),
    kName_("k"),
    viscosityRatioMethod_(false),
    epsilon_(isEpsilonField(iF)),
    mixingLength_(0.0),
    viscosityRatio_(0.0),
    Cmu_(0.09)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


turbulentEpsilonOmegaInletFvPatchScalarField::
turbulentEpsilonOmegaInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    const word& method
)
:
    parent_bctype(p, iF),
    kName_(dict.getOrDefault<word>("k", "k")),
    viscosityRatioMethod_(method == "viscosityRatio"),
    epsilon_(isEpsilonField(iF)),
    mixingLength_(0.0),
    viscosityRatio_(0.0),
    Cmu_(dict.getCheckOrDefault<scalar>("Cmu", 0.09, scalarMinMax::ge(SMALL)))
{
    if (viscosityRatioMethod_)
    {
        viscosityRatio_ =
            dict.getCheck<scalar>("viscosityRatio", scalarMinMax::ge(SMALL));
    }
    else if (method == "mixingLength")
    {
        mixingLength_ =
            dict.getCheck<scalar>("mixingLength", scalarMinMax::ge(SMALL));
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Unknown method " << method << nl
            << "    Valid methods are mixingLength and viscosityRatio"
            << exit(FatalIOError);
    }

    this->phiName_ = dict.getOrDefault<word>("phi", "phi");

    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


turbulentEpsilonOmegaInletFvPatchScalarField::
turbulentEpsilonOmegaInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    turbulentEpsilonOmegaInletFvPatchScalarField
    (
        p,
        iF,
        dict,
        dict.get<word>("method")
    )
{}


turbulentEpsilonOmegaInletFvPatchScalarField::
turbulentEpsilonOmegaInletFvPatchScalarField
(
    const this_bctype& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(ptf, p, iF, mapper),
    kName_(ptf.kName_),
    viscosityRatioMethod_(ptf.viscosityRatioMethod_),
    epsilon_(ptf.epsilon_),
    mixingLength_(ptf.mixingLength_),
    viscosityRatio_(ptf.viscosityRatio_),
    Cmu_(ptf.Cmu_)
{}


turbulentEpsilonOmegaInletFvPatchScalarField::
turbulentEpsilonOmegaInletFvPatchScalarField
(
    const this_bctype& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(ptf, iF),
    kName_(ptf.kName_),
    viscosityRatioMethod_(ptf.viscosityRatioMethod_),
    epsilon_(ptf.epsilon_),
    mixingLength_(ptf.mixingLength_),
    viscosityRatio_(ptf.viscosityRatio_),
    Cmu_(ptf.Cmu_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentEpsilonOmegaInletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Lookup Cmu corresponding to the turbulence model selected
    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    Cmu_ = turbModel.coeffDict().getOrDefault<scalar>("Cmu", Cmu_);

    const auto& kp =
        patch().lookupPatchField<volScalarField>(kName_);

    const auto& phip =
        patch().lookupPatchField<surfaceScalarField>(this->phiName_);

    if (viscosityRatioMethod_)
    {
        const tmp<scalarField> tnup = turbModel.nu(patch().index());
        const scalarField& nup = tnup();

        if (epsilon_)
        {
            this->refValue() = (Cmu_*sqr(kp))/(viscosityRatio_*nup);
        }
        else
        {
            this->refValue() = kp/(viscosityRatio_*nup);
        }
    }
    else
    {
        if (epsilon_)
        {
            this->refValue() =
                (pow(Cmu_, 0.75)/mixingLength_)*pow(kp, 1.5);
        }
        else
        {
            this->refValue() =
                sqrt(kp)/(pow(Cmu_, 0.25)*mixingLength_);
        }
    }

    this->valueFraction() = neg(phip);

    this->parent_bctype::updateCoeffs();
}


void turbulentEpsilonOmegaInletFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry
    (
        "method",
        word(viscosityRatioMethod_ ? "viscosityRatio" : "mixingLength")
    );
    if (viscosityRatioMethod_)
    {
        os.writeEntry("viscosityRatio", viscosityRatio_);
    }
    else
    {
        os.writeEntry("mixingLength", mixingLength_);
    }
    os.writeEntryIfDifferent<scalar>("Cmu", 0.09, Cmu_);
    os.writeEntry("phi", this->phiName_);
    os.writeEntry("k", kName_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentEpsilonOmegaInletFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
