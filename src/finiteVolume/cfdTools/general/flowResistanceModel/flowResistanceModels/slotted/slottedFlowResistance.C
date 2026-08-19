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

#include "slottedFlowResistance.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * //

namespace Foam
{
namespace flowResistanceModels
{
    defineTypeNameAndDebug(Slotted, 0);
    addToRunTimeSelectionTable(flowResistanceModel, Slotted, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowResistanceModels::Slotted::Slotted(const dictionary& dict)
:
    holeFlowResistance(dict),
    slotLength_(-1),
    slotWidth_(-1)
{
    const word coeffsName = "slottedCoeffs";

    // The slot length and width are always required
    if (!dict.found(coeffsName))
    {
        FatalIOErrorInFunction(dict)
            << "resistanceType 'slotted' requires the '"
            << coeffsName << "' sub-dict containing 'slotLength' "
               "and 'slotWidth'"
            << exit(FatalIOError);
    }

    if (dict.found("openAreaFraction"))
    {
        FatalIOErrorInFunction(dict)
            << "openAreaFraction must be specified inside '"
            << coeffsName << "'"
            << exit(FatalIOError);
    }

    const dictionary& typeDict = dict.subDict(coeffsName);

    slotLength_ = typeDict.get<scalar>("slotLength");
    slotWidth_ = typeDict.get<scalar>("slotWidth");

    if (!(slotWidth_ > 0))
    {
        FatalIOErrorInFunction(dict)
            << "slotWidth must be positive"
            << exit(FatalIOError);
    }

    if (!(slotLength_ > slotWidth_))
    {
        FatalIOErrorInFunction(dict)
            << "slotLength must be greater than slotWidth (a > b) "
               "to give a positive straight section"
            << exit(FatalIOError);
    }

    if
    (
        !typeDict.found("openAreaFraction")
     &&
        (
            typeDict.found("longitudinalPitch")
         || typeDict.found("transversePitch")
        )
    )
    {
        const scalar longitudinalPitch =
            typeDict.get<scalar>("longitudinalPitch");
        const scalar transversePitch =
            typeDict.get<scalar>("transversePitch");

        if (!(longitudinalPitch >= slotLength_))
        {
            FatalIOErrorInFunction(dict)
                << "longitudinalPitch must be greater than or equal to "
                   "slotLength to avoid overlapping slots"
                << exit(FatalIOError);
        }

        if (!(transversePitch >= slotWidth_))
        {
            FatalIOErrorInFunction(dict)
                << "transversePitch must be greater than or equal to "
                   "slotWidth to avoid overlapping slots"
                << exit(FatalIOError);
        }
    }

    if (dict.found("nu"))
    {
        FatalIOErrorInFunction(dict)
            << "The 'nu' entry is not permitted for resistanceType "
               "'slotted'. The loss coefficient is independent of the "
               "Reynolds number, so a kinematic viscosity is not required. "
               "Remove the 'nu' entry from flowResistanceCoeffs."
            << exit(FatalIOError);
    }

    resolveOpenAreaFraction(dict, typeDict, "slot");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::flowResistanceModels::Slotted::holeArea() const
{
    return
        slotWidth_*(slotLength_ - slotWidth_)
      + constant::mathematical::pi*sqr(slotWidth_)/4.0;
}


Foam::scalar Foam::flowResistanceModels::Slotted::shapeFactor() const
{
    const scalar ab = slotLength_/slotWidth_;

    return
    (
        ab > 10
      ? 0.81
      : 1 - 0.05*pow(ab - 1, 0.6)
    );
}


void Foam::flowResistanceModels::Slotted::writeCoeffs(Ostream& os) const
{
    os.beginBlock("slottedCoeffs");

    os.writeEntry("slotLength", slotLength_);
    os.writeEntry("slotWidth", slotWidth_);

    if (longitudinalPitch_ > 0)
    {
        os.writeEntry("longitudinalPitch", longitudinalPitch_);
        os.writeEntry("transversePitch", transversePitch_);
    }
    else
    {
        os.writeEntry("openAreaFraction", openAreaFraction_);
    }

    os.endBlock();
}


// ************************************************************************* //
