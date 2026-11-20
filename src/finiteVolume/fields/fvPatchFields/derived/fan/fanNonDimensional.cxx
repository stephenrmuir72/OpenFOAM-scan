/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  | www.openfoam.com
-------------------------------------------------------------------------------
    Copyright (C) 2026 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fanNonDimensional.H"
#include "fvPatch.H"
#include "volFields.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fanNonDimensional::fanNonDimensional
(
    const fvPatch& patch,
    const dimensionSet& dimensions,
    const word& fieldName,
    const fileName& objectPath
)
:
    patch_(patch),
    dimensions_(dimensions),
    fieldName_(fieldName),
    objectPath_(objectPath),
    enabled_(false),
    rpm_(nullptr),
    dm_(nullptr),
    legacy_(false),
    rhoCorrection_(false),
    rhoRefUserProvided_(false),
    rhoRef_(defaultRho())
{}


Foam::fanNonDimensional::fanNonDimensional
(
    const fvPatch& patch,
    const dimensionSet& dimensions,
    const word& fieldName,
    const fileName& objectPath,
    const dictionary& dict,
    const objectRegistry& db,
    const bool enabled
)
:
    patch_(patch),
    dimensions_(dimensions),
    fieldName_(fieldName),
    objectPath_(objectPath),
    enabled_(enabled),
    rpm_(nullptr),
    dm_(nullptr),
    legacy_(dict.getOrDefault("legacy", false)),
    rhoCorrection_(dict.getOrDefault("rhoCorrection", false)),
    rhoRefUserProvided_(false),
    rhoRef_(defaultRho())
{
    rhoRefUserProvided_ = dict.readIfPresent("rhoRef", rhoRef_);

    if (enabled_)
    {
        if (legacy_)
        {
            WarningInFunction
                << "The legacy non-dimensional fan formulation is incorrect "
                   "and should only be used to reproduce existing cases" << nl
                << "    on patch " << patch_.name()
                << " of field " << fieldName_
                << " in file " << objectPath_
                << endl;
        }

        rpm_.reset(Function1<scalar>::New("rpm", dict, &db));
        dm_.reset(Function1<scalar>::New("dm", dict, &db));

        if (dimensions_ == dimPressure && !legacy_)
        {
            if (!rhoCorrection_)
            {
                if (rhoRefUserProvided_)
                {
                    if (rhoRef_ < SMALL)
                    {
                        FatalIOErrorInFunction(dict)
                            << "rhoRef cannot be zero or negative, but is "
                            << rhoRef_ << nl
                            << "    on patch " << patch_.name()
                            << " of field " << fieldName_
                            << " in file " << objectPath_
                            << exit(FatalIOError);
                    }
                }
                else
                {
                    WarningInFunction
                        << "rhoRef not specified; using default rated density "
                        << rhoRef_ << " [kg/m^3]" << nl
                        << "    on patch " << patch_.name()
                        << " of field " << fieldName_
                        << " in file " << objectPath_
                        << endl;
                }
            }
        }
        else if (rhoCorrection_ || rhoRefUserProvided_)
        {
            WarningInFunction
                << "Entries 'rhoCorrection' and/or 'rhoRef' are ignored for "
                   "kinematic pressure fields (dimPressure/dimDensity)" << nl
                << "    on patch " << patch_.name()
                << " of field " << fieldName_
                << " in file " << objectPath_
                << endl;
        }
    }
    else if (rhoCorrection_ || rhoRefUserProvided_)
    {
        WarningInFunction
            << "Entries 'rhoCorrection' and/or 'rhoRef' are ignored unless "
               "non-dimensional operation is enabled" << nl
            << "    on patch " << patch_.name()
            << " of field " << fieldName_
            << " in file " << objectPath_
            << endl;
    }
}


Foam::fanNonDimensional::fanNonDimensional
(
    const fanNonDimensional& rhs,
    const fvPatch& patch,
    const dimensionSet& dimensions,
    const word& fieldName,
    const fileName& objectPath
)
:
    patch_(patch),
    dimensions_(dimensions),
    fieldName_(fieldName),
    objectPath_(objectPath),
    enabled_(rhs.enabled_),
    rpm_(rhs.rpm_.clone()),
    dm_(rhs.dm_.clone()),
    legacy_(rhs.legacy_),
    rhoCorrection_(rhs.rhoCorrection_),
    rhoRefUserProvided_(rhs.rhoRefUserProvided_),
    rhoRef_(rhs.rhoRef_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fanNonDimensional::nonDimensionalFlowRate
(
    const scalar volumeFlowRate,
    const scalar timeValue
) const
{
    const scalar rpm = rpm_->value(timeValue);
    const scalar meanDiam = dm_->value(timeValue);
    const scalar pi = constant::mathematical::pi;

    if (legacy_)
    {
        return
            120.0*volumeFlowRate
           /stabilise
            (
                sqr(pi)*pi*sqr(meanDiam)*meanDiam*rpm,
                VSMALL
            );
    }

    return
        60.0*volumeFlowRate
       /stabilise
        (
            pi*sqr(meanDiam)*meanDiam*rpm,
            VSMALL
        );
}


Foam::scalar Foam::fanNonDimensional::dimensionalPressure
(
    const scalar nonDimensionalPressure,
    const word& rhoName,
    const scalar timeValue
) const
{
    const scalar rpm = rpm_->value(timeValue);
    const scalar meanDiam = dm_->value(timeValue);
    const scalar pi = constant::mathematical::pi;

    if (legacy_)
    {
        return
            nonDimensionalPressure
           *sqr(sqr(pi)*meanDiam*rpm)
           /1800.0;
    }

    if (dimensions_ == dimPressure)
    {
        scalar rhoEff = rhoRef_;

        if (rhoCorrection_)
        {
            const auto& rhop =
                patch_.lookupPatchField<volScalarField>(rhoName);

            rhoEff = gWeightedAverage(patch_.magSf(), rhop);
        }

        return
            nonDimensionalPressure
           *rhoEff*sqr(pi*meanDiam*rpm)
           /3600.0;
    }
    else if (dimensions_ == dimPressure/dimDensity)
    {
        return
            nonDimensionalPressure
           *sqr(pi*meanDiam*rpm)
           /3600.0;
    }

    FatalErrorInFunction
        << "Incorrect pressure dimensions "
        << dimensions_ << nl
        << "    Should be " << dimPressure
        << " for compressible/variable density flow" << nl
        << "    or " << dimPressure/dimDensity
        << " for incompressible flow," << nl
        << "    on patch " << patch_.name()
        << " of field " << fieldName_
        << " in file " << objectPath_
        << exit(FatalError);

    return 0;
}


void Foam::fanNonDimensional::write(Ostream& os) const
{
    if (!enabled_)
    {
        return;
    }

    rpm_->writeData(os);
    dm_->writeData(os);
    os.writeEntryIfDifferent<bool>("legacy", false, legacy_);
    os.writeEntryIfDifferent<bool>("rhoCorrection", false, rhoCorrection_);

    if (rhoRefUserProvided_)
    {
        os.writeEntry("rhoRef", rhoRef_);
    }
}


// ************************************************************************* //
