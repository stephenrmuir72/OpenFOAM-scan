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

#include "screenFlowResistance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * //

namespace Foam
{
namespace flowResistanceModels
{
    defineTypeNameAndDebug(Screen, 0);
    addToRunTimeSelectionTable(flowResistanceModel, Screen, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowResistanceModels::Screen::Screen(const dictionary& dict)
:
    flowResistanceModel(dict),
    openAreaFraction_(-1),
    strandShape_(),
    strandSize_(-1),
    widthSeparation_(-1),
    lengthSeparation_(-1)
{
    const word coeffsName = "screenCoeffs";

    // The strand shape and size are always required (the size enters
    // the strand Reynolds number)
    if (!dict.found(coeffsName))
    {
        FatalIOErrorInFunction(dict)
            << "resistanceType 'screen' requires the '"
            << coeffsName << "' sub-dict containing 'strandShape' "
               "and 'strandSize'"
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

    // Strand shape defaults to 'round' when not specified
    strandShape_ = typeDict.getOrDefault<word>("strandShape", "round");
    strandSize_ = typeDict.get<scalar>("strandSize");

    if (strandShape_ != "round" && strandShape_ != "square"
     && strandShape_ != "flat")
    {
        FatalIOErrorInFunction(dict)
            << "strandShape must be 'round', 'square' or 'flat'"
            << exit(FatalIOError);
    }

    if (strandSize_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "strandSize must be positive"
            << exit(FatalIOError);
    }

    // The open area fraction is either specified directly or derived
    // from the strand geometry and the separations; the two options
    // are exclusive
    const bool haveSeparations =
        typeDict.found("widthSeparation")
     || typeDict.found("lengthSeparation");

    if (typeDict.found("openAreaFraction") && haveSeparations)
    {
        FatalIOErrorInFunction(dict)
            << "openAreaFraction must not be specified together "
               "with 'widthSeparation'/'lengthSeparation'"
            << exit(FatalIOError);
    }

    if (typeDict.found("openAreaFraction"))
    {
        openAreaFraction_ = typeDict.get<scalar>("openAreaFraction");

        if (openAreaFraction_ <= 0 || openAreaFraction_ > 1)
        {
            FatalIOErrorInFunction(dict)
                << "openAreaFraction must be in (0, 1]"
                << exit(FatalIOError);
        }
    }
    else
    {
        widthSeparation_ = typeDict.get<scalar>("widthSeparation");
        lengthSeparation_ = typeDict.get<scalar>("lengthSeparation");

        openAreaFraction_ =
            (1 - strandSize_/widthSeparation_)
           *(1 - strandSize_/lengthSeparation_);

        if (openAreaFraction_ <= 0 || openAreaFraction_ > 1)
        {
            FatalIOErrorInFunction(dict)
                << "Computed openAreaFraction " << openAreaFraction_
                << " is not in (0, 1]. Check strand geometry and "
                   "separation."
                << exit(FatalIOError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::flowResistanceModels::Screen::kShape() const
{
    if (strandShape_ == "round")
    {
        return 1.3;
    }
    else if (strandShape_ == "square")
    {
        return 1.8;
    }
    else
    {
        return 2.1;
    }
}


Foam::scalar Foam::flowResistanceModels::Screen::zetaTurbulent() const
{
    const scalar f = openAreaFraction_;

    return kShape()*(1 - f) + sqr((1 - f)/f);
}


Foam::scalar Foam::flowResistanceModels::Screen::lossCoefficient
(
    const scalar UnAvg,
    const scalar nuAvg
) const
{
    const scalar f = openAreaFraction_;

    if (f <= 0 || f > 1)
    {
        FatalErrorInFunction
            << "Valid openAreaFraction required for resistance"
            << abort(FatalError);
    }

    if (nuAvg <= 0)
    {
        FatalErrorInFunction
            << "Kinematic viscosity must be positive to evaluate the "
               "strand Reynolds number, got " << nuAvg << ". Check the "
               "registered nu/mu field or specify a positive 'nu' "
               "entry in flowResistanceCoeffs."
            << abort(FatalError);
    }

    // Strand Reynolds number, clamped to a minimum of 1 to keep the
    // laminar term finite at (near) zero velocity
    const scalar Re =
        max(mag(UnAvg)*strandSize_/(f*nuAvg), scalar(1));

    // Cache for debug reporting (base class)
    Re_ = Re;

    const scalar zetaT = zetaTurbulent();

    if (Re < 50)
    {
        return 22/Re + zetaT;
    }
    else if (Re < 1000)
    {
        return (1 + 20/Re)*zetaT;
    }
    else
    {
        return zetaT;
    }
}


void Foam::flowResistanceModels::Screen::writeCoeffs(Ostream& os) const
{
    os.beginBlock("screenCoeffs");

    os.writeEntry("strandShape", strandShape_);
    os.writeEntry("strandSize", strandSize_);

    if (widthSeparation_ > 0)
    {
        os.writeEntry("widthSeparation", widthSeparation_);
        os.writeEntry("lengthSeparation", lengthSeparation_);
    }
    else
    {
        os.writeEntry("openAreaFraction", openAreaFraction_);
    }

    os.endBlock();
}


// ************************************************************************* //
