/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "totalPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::totalPressureFvPatchScalarField::totalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(p, iF),
    UName_("U"),
    phiName_("phi"),
    rhoName_("rho"),
    psiName_("none"),
    gamma_(0.0),
    p0_(p.size(), Zero),
    flowResistanceModel_(),
    flowResistanceRelax_(1),
    flowResistancePressureCorrection_(0),
    flowResistancePressureCorrectionInitialised_(false)
{}


Foam::totalPressureFvPatchScalarField::totalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict, IOobjectOption::NO_READ),
    UName_(dict.getOrDefault<word>("U", "U")),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    psiName_(dict.getOrDefault<word>("psi", "none")),
    gamma_(psiName_ != "none" ? dict.get<scalar>("gamma") : 1),
    p0_("p0", dict, p.size()),
    flowResistanceModel_(),
    flowResistanceRelax_(1),
    flowResistancePressureCorrection_(0),
    flowResistancePressureCorrectionInitialised_(false)
{
    if (!this->readValueEntry(dict))
    {
        fvPatchField<scalar>::operator=(p0_);
    }

    // Read optional flow resistance
    if (dict.found("flowResistanceCoeffs"))
    {
        const dictionary& resistanceDict =
            dict.subDict("flowResistanceCoeffs");

        flowResistanceModel_.reset
        (
            new flowResistanceModel(resistanceDict)
        );

        flowResistanceRelax_ =
            resistanceDict.getOrDefault<scalar>("relax", 1);

        if (flowResistanceRelax_ <= 0 || flowResistanceRelax_ > 1)
        {
            FatalIOErrorInFunction(resistanceDict)
                << "The flow resistance relaxation factor must be in "
                   "the range (0, 1]"
                << exit(FatalIOError);
        }

        if (resistanceDict.found("relaxedPressureCorrection"))
        {
            flowResistancePressureCorrection_ =
                resistanceDict.get<scalar>("relaxedPressureCorrection");
            flowResistancePressureCorrectionInitialised_ = true;
        }
    }
}


Foam::totalPressureFvPatchScalarField::totalPressureFvPatchScalarField
(
    const this_bctype& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_),
    p0_(ptf.p0_, mapper),
    flowResistanceModel_(ptf.flowResistanceModel_.clone()),
    flowResistanceRelax_(ptf.flowResistanceRelax_),
    flowResistancePressureCorrection_
    (
        ptf.flowResistancePressureCorrection_
    ),
    flowResistancePressureCorrectionInitialised_
    (
        ptf.flowResistancePressureCorrectionInitialised_
    )
{}


Foam::totalPressureFvPatchScalarField::totalPressureFvPatchScalarField
(
    const this_bctype& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(tppsf, iF),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    p0_(tppsf.p0_),
    flowResistanceModel_(tppsf.flowResistanceModel_.clone()),
    flowResistanceRelax_(tppsf.flowResistanceRelax_),
    flowResistancePressureCorrection_
    (
        tppsf.flowResistancePressureCorrection_
    ),
    flowResistancePressureCorrectionInitialised_
    (
        tppsf.flowResistancePressureCorrectionInitialised_
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::totalPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    this->parent_bctype::autoMap(m);
    p0_.autoMap(m);
}


void Foam::totalPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    this->parent_bctype::rmap(ptf, addr);

    const auto& tiptf = refCast<const this_bctype>(ptf);

    p0_.rmap(tiptf.p0_, addr);
}


void Foam::totalPressureFvPatchScalarField::updateCoeffs
(
    const scalarField& p0p,
    const vectorField& Up
)
{
    if (updated())
    {
        return;
    }

    const auto& phip = patch().lookupPatchField<surfaceScalarField>(phiName_);

    if (internalField().dimensions() == dimPressure)
    {
        if (psiName_ == "none")
        {
            // Variable density and low-speed compressible flow

            const auto& rho =
                patch().lookupPatchField<volScalarField>(rhoName_);

            operator==(p0p - 0.5*rho*(neg(phip))*magSqr(Up));
        }
        else
        {
            // High-speed compressible flow

            const auto& psip =
                patch().lookupPatchField<volScalarField>(psiName_);

            if (gamma_ > 1)
            {
                scalar gM1ByG = (gamma_ - 1)/gamma_;

                operator==
                (
                    p0p
                   /pow
                    (
                        (1.0 + 0.5*psip*gM1ByG*(neg(phip))*magSqr(Up)),
                        1.0/gM1ByG
                    )
                );
            }
            else
            {
                operator==(p0p/(1.0 + 0.5*psip*(neg(phip))*magSqr(Up)));
            }
        }

    }
    else if (internalField().dimensions() == dimPressure/dimDensity)
    {
        // Incompressible flow
        operator==(p0p - 0.5*(neg(phip))*magSqr(Up));
    }
    else
    {
        FatalErrorInFunction
            << " Incorrect pressure dimensions " << internalField().dimensions()
            << nl
            << "    Should be " << dimPressure
            << " for compressible/variable density flow" << nl
            << "    or " << dimPressure/dimDensity
            << " for incompressible flow," << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    this->parent_bctype::updateCoeffs();
}


void Foam::totalPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (flowResistanceModel_)
    {
        // Flow resistance active — modify p0 by the pressure loss

        const vectorField& Up =
            patch().lookupPatchField<volVectorField>(UName());
        const scalarField& phip =
            patch().lookupPatchField<surfaceScalarField>(phiName_);

        // Loss coefficient K
        const scalar K = flowResistanceModel_->lossCoefficient();

        const scalar A = gSum(patch().magSf());

        // Uniform pressure drop: Δp = K * 0.5 * rho * UnAvg^2
        // (incompressible: no rho factor)
        scalar deltaP;
        scalar UnAvg;

        if (internalField().dimensions() == dimPressure)
        {
            const scalarField& rhop =
                patch().lookupPatchField<volScalarField>(rhoName_);

            UnAvg = gSum(phip/rhop)/A;
            const scalar rhoAvg = gSum(rhop*patch().magSf())/A;

            deltaP = K * 0.5 * rhoAvg * sqr(UnAvg);
        }
        else
        {
            UnAvg = gSum(phip)/A;

            deltaP = K * 0.5 * sqr(UnAvg);
        }

        // Relax the signed correction so that flow reversal is also smooth
        const scalar targetPressureCorrection = sign(UnAvg)*deltaP;

        if (!flowResistancePressureCorrectionInitialised_)
        {
            flowResistancePressureCorrection_ = targetPressureCorrection;
            flowResistancePressureCorrectionInitialised_ = true;
        }
        else
        {
            flowResistancePressureCorrection_ =
                (1 - flowResistanceRelax_)
               *flowResistancePressureCorrection_
              + flowResistanceRelax_*targetPressureCorrection;
        }

        // Apply one uniform correction based on the net patch flow direction
        const scalarField p0Adj
        (
            p0_ + flowResistancePressureCorrection_
        );

        updateCoeffs(p0Adj, Up);

        if (debug)
        {
            Info<< "Patch average velocity = " << UnAvg << nl
                << "Target pressure correction = "
                << targetPressureCorrection << nl
                << "Relaxed pressure correction = "
                << flowResistancePressureCorrection_ << endl;
        }
    }
    else
    {
        // Standard behaviour — no flow resistance
        updateCoeffs(p0(), patch().lookupPatchField<volVectorField>(UName()));
    }
}


void Foam::totalPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntry("rho", rhoName_);
    os.writeEntry("psi", psiName_);
    os.writeEntry("gamma", gamma_);
    p0_.writeEntry("p0", os);

    // Write flowResistanceCoeffs if active
    if (flowResistanceModel_)
    {
        os.beginBlock("flowResistanceCoeffs");
        flowResistanceModel_->write(os);
        os.writeEntryIfDifferent<scalar>
        (
            "relax",
            1,
            flowResistanceRelax_
        );
        if
        (
            flowResistanceRelax_ < 1
         && flowResistancePressureCorrectionInitialised_
        )
        {
            os.writeEntry
            (
                "relaxedPressureCorrection",
                flowResistancePressureCorrection_
            );
        }
        os.endBlock();
    }

    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        totalPressureFvPatchScalarField
    );
}

// ************************************************************************* //
