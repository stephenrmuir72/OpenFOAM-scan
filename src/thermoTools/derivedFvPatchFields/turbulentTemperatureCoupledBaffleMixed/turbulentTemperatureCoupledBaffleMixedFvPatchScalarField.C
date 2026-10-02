/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "turbulentTemperatureCoupledBaffleMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(p, iF),
    temperatureCoupledBase(patch()),  // default method (fluidThermo)
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(p, iF),
        *this
    ),
    TnbrName_("undefined-Tnbr")
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 1.0;
}


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const this_bctype& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(p, iF),
        *this,
        ptf
    ),
    TnbrName_(ptf.TnbrName_),
    thicknessLayers_(ptf.thicknessLayers_),
    thicknessLayer_(ptf.thicknessLayer_.clone(p.patch())),
    kappaLayers_(ptf.kappaLayers_),
    kappaLayer_(ptf.kappaLayer_.clone(p.patch()))
{}


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF),
    temperatureCoupledBase(patch(), dict),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(p, iF),
        *this,
        dict
    ),
    TnbrName_(dict.get<word>("Tnbr"))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    WarningInFunction
        << "This BC has been superseded by "
        << "compressible::turbulentTemperatureRadCoupledMixed "
        << "which has more functionalities and it can handle "
        << "the assemble coupled option for energy. "
        << endl;

    // Read list of layers
    if (dict.readIfPresent("thicknessLayers", thicknessLayers_))
    {
        dict.readEntry("kappaLayers", kappaLayers_);
    }
    // Read single additional PatchFunction1
    thicknessLayer_ = PatchFunction1<scalar>::NewIfPresent
    (
        p.patch(),
        "thicknessLayer",
        dict
    );
    kappaLayer_ = PatchFunction1<scalar>::NewIfPresent
    (
        p.patch(),
        "kappaLayer",
        dict
    );


    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    if (this->readMixedEntries(dict))
    {
        // Full restart
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = Zero;
        valueFraction() = 1.0;
    }

// This blocks (crashes) with more than two worlds!
//
///    // Store patch value as initial guess when running in database mode
///    mappedPatchFieldBase<scalar>::initRetrieveField
///    (
///        this->internalField().name(),
///        *this
///    );
///    mappedPatchFieldBase<scalar>::initRetrieveField
///    (
///        this->internalField().name() + "_weights",
///        this->patch().deltaCoeffs()
///    );
}


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const this_bctype& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(wtcsf, iF),
    temperatureCoupledBase(patch(), wtcsf),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(patch(), iF),
        *this,
        wtcsf
    ),
    TnbrName_(wtcsf.TnbrName_),
    thicknessLayers_(wtcsf.thicknessLayers_),
    thicknessLayer_(wtcsf.thicknessLayer_.clone(patch().patch())),
    kappaLayers_(wtcsf.kappaLayers_),
    kappaLayer_(wtcsf.kappaLayer_.clone(patch().patch()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    this->parent_bctype::autoMap(mapper);
    temperatureCoupledBase::autoMap(mapper);
    //mappedPatchFieldBase<scalar>::autoMap(mapper);
    if (thicknessLayer_)
    {
        thicknessLayer_().autoMap(mapper);
        kappaLayer_().autoMap(mapper);
    }
}


void turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    this->parent_bctype::rmap(ptf, addr);

    const auto& tiptf = refCast<const this_bctype>(ptf);

    temperatureCoupledBase::rmap(tiptf, addr);
    //mappedPatchFieldBase<scalar>::rmap(ptf, addr);
    if (thicknessLayer_ && tiptf.thicknessLayer_)
    {
        thicknessLayer_().rmap(tiptf.thicknessLayer_(), addr);
        kappaLayer_().rmap(tiptf.kappaLayer_(), addr);
    }
}


scalarField
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::layerResistance()
const
{
    scalarField resistance(patch().size(), Zero);
    if (thicknessLayer_)
    {
        const scalar t = db().time().timeOutputValue();
        resistance += thicknessLayer_().value(t)/kappaLayer_().value(t);
    }
    forAll(thicknessLayers_, iLayer)
    {
        resistance += thicknessLayers_[iLayer]/kappaLayers_[iLayer];
    }
    return resistance;
}


tmp<Foam::scalarField>
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::kappa
(
    const scalarField& Tp
) const
{
    // Get kappa from relevant thermo
    tmp<scalarField> tk(temperatureCoupledBase::kappa(Tp));

    // Optionally modify with explicit resistance
    if (thicknessLayer_ || thicknessLayers_.size())
    {
        scalarField KDelta(tk*patch().deltaCoeffs());

        // Harmonic averaging of kappa*deltaCoeffs
        {
            KDelta = 1.0/KDelta;
            KDelta += layerResistance();
            KDelta = 1.0/KDelta;
        }

        // Update kappa from KDelta
        tk = KDelta/patch().deltaCoeffs();
    }

    return tk;
}


void turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::incrMsgType();

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        mappedPatchFieldBase<scalar>::mapper
        (
            patch(),
            this->internalField()
        );

    const scalarField& Tp = *this;
    const auto myFlux = wallCoupledFlux(Tp);
    const scalarField myResistance(layerResistance());

    scalarField nbrIntFld;
    scalarField nbrKDelta;
    scalarField nbrNormalKDelta;
    scalarField nbrCorrection;
    scalarField nbrResistance;
    scalarField nbrDirectional(patch().size(), myFlux.directional ? 1 : 0);
    if (mpp.sameWorld())
    {
        // Same world so lookup
        const auto& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
        const label nbrPatchID = mpp.samplePolyPatch().index();
        const auto& nbrPatch = nbrMesh.boundary()[nbrPatchID];

        const auto& nbrField = refCast<const this_bctype>
        (
            nbrPatch.lookupPatchField<volScalarField>(TnbrName_)
        );

        // Swap to obtain full local values of neighbour K*delta
        nbrIntFld = nbrField.patchInternalField();
        nbrKDelta = nbrField.kappa(nbrField)*nbrPatch.deltaCoeffs();
        auto nbrFlux = nbrField.wallCoupledFlux(nbrField);
        nbrNormalKDelta = std::move(nbrFlux.conductance);
        nbrCorrection = std::move(nbrFlux.correction);
        nbrResistance = nbrField.layerResistance();
        nbrDirectional = scalarField
        (
            nbrPatch.size(),
            nbrFlux.directional ? 1 : 0
        );
    }
    else
    {
        // Different world so use my region,patch. Distribution below will
        // do the reordering.
        nbrIntFld = patchInternalField();
        nbrKDelta = kappa(Tp)*patch().deltaCoeffs();
        nbrNormalKDelta = myFlux.conductance;
        nbrCorrection = myFlux.correction;
        nbrResistance = myResistance;
    }
    distribute(this->internalField().name() + "_value", nbrIntFld);
    distribute(this->internalField().name() + "_weights", nbrKDelta);
    distribute
    (
        this->internalField().name() + "_normalWeights",
        nbrNormalKDelta
    );
    distribute
    (
        this->internalField().name() + "_correction",
        nbrCorrection
    );
    distribute
    (
        this->internalField().name() + "_layerResistance",
        nbrResistance
    );
    distribute
    (
        this->internalField().name() + "_directional",
        nbrDirectional
    );

    const bool useTensor =
        myFlux.directional || (max(nbrDirectional) > 0.5);

    scalarField myKDelta;
    scalarField myCorrection;
    if (useTensor)
    {
        myKDelta = myFlux.conductance;
        myCorrection = myFlux.correction;

        const scalarField rawNbrKDelta(nbrNormalKDelta);
        nbrKDelta =
            1/(1/rawNbrKDelta + myResistance + nbrResistance);
        nbrCorrection *= nbrKDelta/rawNbrKDelta;
    }
    else
    {
        myKDelta = kappa(Tp)*patch().deltaCoeffs();
        myCorrection.setSize(patch().size(), Zero);
    }


    // The mapped conductance includes both layers; refGrad balances the
    // explicit face flux on both sides.

    this->refValue() = nbrIntFld;
    this->refGrad() = Zero;
    if (useTensor)
    {
        this->refGrad() =
            -(myCorrection + nbrCorrection)
           *patch().deltaCoeffs()/myKDelta;
    }
    this->valueFraction() = nbrKDelta/(nbrKDelta + myKDelta);

    this->parent_bctype::updateCoeffs();

    if (debug)
    {
        const scalar Q =
            useTensor
          ? gSum
            (
                patch().magSf()
               *(myKDelta*snGrad()/patch().deltaCoeffs() + myCorrection)
            )
          : gSum(kappa(Tp)*patch().magSf()*snGrad());

        auto limits = gMinMax(*this);
        auto avg = gAverage(*this);

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << mpp.sampleRegion() << ':'
            << mpp.samplePatch() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature "
            << " min:" << limits.min()
            << " max:" << limits.max()
            << " avg:" << avg
            << endl;
    }

    UPstream::msgType(oldTag);  // Restore tag
}


void turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& m,
    const label iMatrix,
    const direction cmpt
)
{
    FatalErrorInFunction
        << "This BC does not support energy coupling "
        << "Use compressible::turbulentTemperatureRadCoupledMixed "
        << "which has more functionalities and it can handle "
        << "the assemble coupled option for energy. "
        << abort(FatalError);
}


tmp<Field<scalar>>
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::coeffs
(
    fvMatrix<scalar>& matrix,
    const Field<scalar>& coeffs,
    const label mat
) const
{
    FatalErrorInFunction
        << "This BC does not support energy coupling "
        << "Use compressible::turbulentTemperatureRadCoupledMixed "
        << "which has more functionalities and it can handle "
        << "the assemble coupled option for energy. "
        << abort(FatalError);

    return nullptr;
}


void turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    this->parent_bctype::write(os);
    os.writeEntry("Tnbr", TnbrName_);
    if (thicknessLayer_)
    {
        thicknessLayer_().writeData(os);
        kappaLayer_().writeData(os);
    }
    if (thicknessLayers_.size())
    {
        thicknessLayers_.writeEntry("thicknessLayers", os);
        kappaLayers_.writeEntry("kappaLayers", os);
    }
    temperatureCoupledBase::write(os);
    mappedPatchFieldBase<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
