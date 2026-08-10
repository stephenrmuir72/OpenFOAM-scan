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
    nuName_("nu"),
    muName_("thermo:mu"),
    psiName_("none"),
    gamma_(0.0),
    p0_(p.size(), Zero),
    flowResistanceModel_()
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
    nuName_(dict.getOrDefault<word>("nu", "nu")),
    muName_(dict.getOrDefault<word>("mu", "thermo:mu")),
    psiName_(dict.getOrDefault<word>("psi", "none")),
    gamma_(psiName_ != "none" ? dict.get<scalar>("gamma") : 1),
    p0_("p0", dict, p.size()),
    flowResistanceModel_()
{
    if (!this->readValueEntry(dict))
    {
        fvPatchField<scalar>::operator=(p0_);
    }

    // Read optional flow resistance
    if (dict.found("flowResistanceCoeffs"))
    {
        flowResistanceModel_ = flowResistanceModel::New
        (
            dict.subDict("flowResistanceCoeffs")
        );
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
    nuName_(ptf.nuName_),
    muName_(ptf.muName_),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_),
    p0_(ptf.p0_, mapper),
    flowResistanceModel_(ptf.flowResistanceModel_.clone())
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
    nuName_(tppsf.nuName_),
    muName_(tppsf.muName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    p0_(tppsf.p0_),
    flowResistanceModel_(tppsf.flowResistanceModel_.clone())
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

        const scalar A = gSum(patch().magSf());

        if (A <= 0)
        {
            FatalErrorInFunction
                << "Zero or negative patch area (" << A << ") on patch "
                << this->patch().name() << " of field "
                << this->internalField().name()
                << "; cannot evaluate the flow resistance correction."
                << exit(FatalError);
        }

        // Velocity-flux sum and area-averaged density (1 incompressible).
        // Mirror the same dimension validation used by the
        // updateCoeffs(p0p, Up) overload below, so an invalid pressure
        // dimension is caught here rather than silently treated as
        // incompressible.
        scalar Q;
        scalar rhoAvg;

        if (internalField().dimensions() == dimPressure)
        {
            const scalarField& rhop =
                patch().lookupPatchField<volScalarField>(rhoName_);

            Q = gSum(phip/rhop);
            rhoAvg = gSum(rhop*patch().magSf())/A;
        }
        else if (internalField().dimensions() == dimPressure/dimDensity)
        {
            Q = gSum(phip);
            rhoAvg = 1;
        }
        else
        {
            FatalErrorInFunction
                << " Incorrect pressure dimensions "
                << internalField().dimensions() << nl
                << "    Should be " << dimPressure
                << " for compressible/variable density flow" << nl
                << "    or " << dimPressure/dimDensity
                << " for incompressible flow," << nl
                << "    on patch " << this->patch().name()
                << " of field " << this->internalField().name()
                << " in file " << this->internalField().objectPath()
                << exit(FatalError);
        }

        // Patch-averaged normal velocity
        const scalar UnAvg = Q/A;

        // Kinematic viscosity: user-specified, otherwise looked up
        // directly from the registered nu/mu field (nu = mu/rho for
        // compressible flow), avoiding a hard dependency on a
        // turbulenceModel (which would introduce a circular library
        // dependency, since TurbulenceModels itself links against
        // finiteVolume). Only required for the Reynolds-number-
        // dependent resistance types (perforated/slotted/screen); no
        // viscosity field need be registered for the 'constant' type.
        scalar nuAvg = flowResistanceModel_->userNu();

        if (flowResistanceModel_->needsNu())
        {
            const word nuName(IOobject::groupName(nuName_, internalField().group()));
            const word muName(IOobject::groupName(muName_, internalField().group()));

            if (db().foundObject<volScalarField>(nuName))
            {
                const auto& nu =
                    patch().lookupPatchField<volScalarField>(nuName);

                nuAvg = gSum(nu*patch().magSf())/A;
            }
            else if (db().foundObject<volScalarField>(muName))
            {
                const auto& mu =
                    patch().lookupPatchField<volScalarField>(muName);

                const scalarField& rhop =
                    patch().lookupPatchField<volScalarField>(rhoName_);

                nuAvg = gSum((mu/rhop)*patch().magSf())/A;
            }
            else
            {
                FatalErrorInFunction
                    << "Cannot find a registered '" << nuName
                    << "' or '" << muName << "' field to evaluate the "
                       "kinematic viscosity for the flow resistance "
                       "model. Specify the 'nu' entry in "
                       "flowResistanceCoeffs instead."
                    << exit(FatalError);
            }
        }

        // Patch-averaged normal velocity and relaxed pressure correction
        const scalar correction =
            flowResistanceModel_->updateCorrection(Q, A, rhoAvg, nuAvg);

        // Apply one uniform correction based on the net patch flow
        // direction
        updateCoeffs(p0_ + correction, Up);

        if (debug)
        {
            Info<< "Patch average velocity = " << UnAvg << nl
                << "Target pressure correction = "
                << flowResistanceModel_->targetCorrection() << nl
                << "Relaxed pressure correction = "
                << flowResistanceModel_->relaxedCorrection() << endl;
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
    os.writeEntryIfDifferent<word>("nu", "nu", nuName_);
    os.writeEntryIfDifferent<word>("mu", "thermo:mu", muName_);
    os.writeEntry("psi", psiName_);
    os.writeEntry("gamma", gamma_);
    p0_.writeEntry("p0", os);

    // Write flowResistanceCoeffs if active
    if (flowResistanceModel_)
    {
        os.beginBlock("flowResistanceCoeffs");
        flowResistanceModel_->write(os);
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
