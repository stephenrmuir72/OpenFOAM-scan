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

#include "holeFlowResistance.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowResistanceModels::holeFlowResistance::holeFlowResistance
(
    const dictionary& dict
)
:
    flowResistanceModel(dict),
    openAreaFraction_(-1),
    longitudinalPitch_(-1),
    transversePitch_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::flowResistanceModels::holeFlowResistance::resolveOpenAreaFraction
(
    const dictionary& dict,
    const dictionary& typeDict,
    const word& geometryDesc
)
{
    // The open area fraction is either specified directly or derived
    // from the geometry and the pitches; the two options are exclusive
    const bool havePitches =
        typeDict.found("longitudinalPitch")
     || typeDict.found("transversePitch");

    if (typeDict.found("openAreaFraction") && havePitches)
    {
        FatalIOErrorInFunction(dict)
            << "openAreaFraction must not be specified together "
               "with 'longitudinalPitch'/'transversePitch'"
            << exit(FatalIOError);
    }

    scalar f = -1;

    if (typeDict.found("openAreaFraction"))
    {
        f = typeDict.get<scalar>("openAreaFraction");

        if (!(f > 0 && f <= 1))
        {
            FatalIOErrorInFunction(dict)
                << "openAreaFraction must be in (0, 1]"
                << exit(FatalIOError);
        }
    }
    else
    {
        longitudinalPitch_ = typeDict.get<scalar>("longitudinalPitch");
        transversePitch_ = typeDict.get<scalar>("transversePitch");

        if (!(longitudinalPitch_ > 0 && transversePitch_ > 0))
        {
            FatalIOErrorInFunction(dict)
                << "longitudinalPitch and transversePitch must be positive"
                << exit(FatalIOError);
        }

        f = holeArea()/(longitudinalPitch_*transversePitch_);

        if (!(f > 0 && f <= 1))
        {
            FatalIOErrorInFunction(dict)
                << "Computed openAreaFraction " << f
                << " is not in (0, 1]. "
                   "Check " << geometryDesc << " geometry and pitch."
                << exit(FatalIOError);
        }
    }

    openAreaFraction_ = f;

    return f;
}


Foam::scalar
Foam::flowResistanceModels::holeFlowResistance::shapeFactor() const
{
    return 1;
}


Foam::scalar
Foam::flowResistanceModels::holeFlowResistance::zetaTurbulent() const
{
    const scalar f = openAreaFraction_;
    const scalar zetaT = sqr(0.707*pow(1 - f, 0.375) + (1 - f))/sqr(f);

    return zetaT*shapeFactor();
}


Foam::scalar
Foam::flowResistanceModels::holeFlowResistance::lossCoefficient
(
    const scalar,
    const scalar
) const
{
    const scalar f = openAreaFraction_;

    // Defensive invariant: openAreaFraction_ is validated in the constructor
    if (!(f > 0 && f <= 1))
    {
        FatalErrorInFunction
            << "Valid openAreaFraction required for resistance"
            << abort(FatalError);
    }

    // Re-independent hole-flow resistance: only the fully turbulent
    // loss coefficient is used (Eq. 13, including the slot shape
    // factor for the slotted type)
    return zetaTurbulent();
}


void
Foam::flowResistanceModels::holeFlowResistance::validatePatchArea
(
    const scalar A
) const
{
    const scalar Ahole = holeArea();

    if (openAreaFraction_*A < Ahole)
    {
        FatalErrorInFunction
            << "The open area fraction " << openAreaFraction_
            << " and hole area " << Ahole
            << " imply fewer than one hole on a patch of area "
            << A << ". Check openAreaFraction, hole geometry and "
               "pitch."
            << abort(FatalError);
    }
}


// ************************************************************************* //
