/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 OpenFOAM Foundation
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

#include "flowResistanceModel.H"

// * * * * * * * * * * * * * * * * Static Data Members  * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(flowResistanceModel, 0);
    defineRunTimeSelectionTable(flowResistanceModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowResistanceModel::flowResistanceModel(const dictionary& dict)
:
    relax_(dict.getOrDefault<scalar>("relax", 1)),
    userNu_(dict.getOrDefault<scalar>("nu", -1)),
    targetCorrection_(0),
    relaxedCorrection_(0),
    relaxedCorrectionInitialised_(false),
    patchAreaValidated_(false)
{
    // Validate the relaxation factor
    if (relax_ <= 0 || relax_ > 1)
    {
        FatalIOErrorInFunction(dict)
            << "The flow resistance relaxation factor must be in "
               "the range (0, 1]"
            << exit(FatalIOError);
    }

    // Resume from a previously written relaxed correction
    if (dict.found("relaxedPressureCorrection"))
    {
        relaxedCorrection_ =
            dict.get<scalar>("relaxedPressureCorrection");
        relaxedCorrectionInitialised_ = true;
    }
}


Foam::flowResistanceModel::flowResistanceModel
(
    const flowResistanceModel& m
)
:
    relax_(m.relax_),
    userNu_(m.userNu_),
    targetCorrection_(m.targetCorrection_),
    relaxedCorrection_(m.relaxedCorrection_),
    relaxedCorrectionInitialised_(m.relaxedCorrectionInitialised_),
    patchAreaValidated_(m.patchAreaValidated_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::flowResistanceModel::validatePatchArea(const scalar A) const
{
    // Default: no geometry to validate
}


Foam::scalar Foam::flowResistanceModel::updateCorrection
(
    const scalar Q,
    const scalar A,
    const scalar rhoAvg,
    const scalar nuAvg
)
{
    // Validate the resistance geometry against the patch area on the
    // first update (no-op for resistance types without hole geometry)
    if (!patchAreaValidated_)
    {
        validatePatchArea(A);
        patchAreaValidated_ = true;
    }

    const scalar UnAvg = Q/A;
    const scalar deltaP =
        lossCoefficient(UnAvg, nuAvg)*0.5*rhoAvg*sqr(UnAvg);

    // Signed target: inflow (UnAvg < 0) lowers p0, outflow raises it.
    //
    // Sign convention: Q (and hence UnAvg = Q/A) is derived from the
    // patch flux phi = U.Sf, and Sf is outward-pointing on a boundary
    // patch, so UnAvg > 0 means net outflow and UnAvg < 0 means net
    // inflow. p0_ represents the total pressure on the far (reservoir)
    // side of the resistive element:
    //  - On inflow the resistance dissipates part of the driving
    //    pressure before it reaches the patch, so the effective total
    //    pressure applied to the patch (p0_ + correction) must be
    //    lower than the nominal reservoir p0_, i.e. correction < 0.
    //  - On outflow the domain interior must overcome both the
    //    reservoir pressure and the resistance's back-pressure, so the
    //    effective total pressure is higher, i.e. correction > 0.
    // deltaP = zeta*0.5*rhoAvg*UnAvg^2 >= 0 always, so multiplying by
    // sign(UnAvg) gives exactly this behaviour and the correction
    // varies continuously (approximately linearly in UnAvg for
    // constant zeta) through UnAvg = 0.
    targetCorrection_ = sign(UnAvg)*deltaP;

    if (!relaxedCorrectionInitialised_)
    {
        // First update takes the full step
        relaxedCorrection_ = targetCorrection_;
        relaxedCorrectionInitialised_ = true;
    }
    else
    {
        relaxedCorrection_ =
            lerp(relaxedCorrection_, targetCorrection_, relax_);
    }

    return relaxedCorrection_;
}


void Foam::flowResistanceModel::write(Ostream& os) const
{
    os.writeEntry("resistanceType", type());

    writeCoeffs(os);

    os.writeEntryIfDifferent<scalar>("relax", 1, relax_);

    os.writeEntryIfDifferent<scalar>("nu", -1, userNu_);

    if (relax_ < 1 && relaxedCorrectionInitialised_)
    {
        os.writeEntry("relaxedPressureCorrection", relaxedCorrection_);
    }
}


// ************************************************************************* //
