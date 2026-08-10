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

#include "perforatedFlowResistance.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * //

namespace Foam
{
namespace flowResistanceModels
{
    defineTypeNameAndDebug(Perforated, 0);
    addToRunTimeSelectionTable(flowResistanceModel, Perforated, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowResistanceModels::Perforated::Perforated(const dictionary& dict)
:
    holeFlowResistance(dict),
    holeShape_(),
    holeSize_(-1)
{
    const word coeffsName = "perforatedCoeffs";

    // The hole shape and size are always required
    if (!dict.found(coeffsName))
    {
        FatalIOErrorInFunction(dict)
            << "resistanceType 'perforated' requires the '"
            << coeffsName << "' sub-dict containing 'holeShape' "
               "and 'holeSize'"
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

    holeShape_ = typeDict.get<word>("holeShape");
    holeSize_ = typeDict.get<scalar>("holeSize");

    if (holeShape_ != "circular" && holeShape_ != "square")
    {
        FatalIOErrorInFunction(dict)
            << "holeShape must be 'circular' or 'square'"
            << exit(FatalIOError);
    }

    resolveOpenAreaFraction(dict, typeDict, "hole");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::flowResistanceModels::Perforated::holeArea() const
{
    return
    (
        holeShape_ == "circular"
      ? constant::mathematical::pi*sqr(holeSize_)/4.0
      : sqr(holeSize_)
    );
}


Foam::scalar Foam::flowResistanceModels::Perforated::holePerimeter() const
{
    return
    (
        holeShape_ == "circular"
      ? constant::mathematical::pi*holeSize_
      : 4.0*holeSize_
    );
}


void Foam::flowResistanceModels::Perforated::writeCoeffs(Ostream& os) const
{
    os.beginBlock("perforatedCoeffs");

    os.writeEntry("holeShape", holeShape_);
    os.writeEntry("holeSize", holeSize_);

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
