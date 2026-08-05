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
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowResistanceModel::flowResistanceModel(const dictionary& dict)
:
    resistanceType_(dict.get<word>("resistanceType")),
    openAreaFraction_(dict.getOrDefault<scalar>("openAreaFraction", -1)),
    lossCoefficient_(dict.getOrDefault<scalar>("lossCoefficient", -1)),
    holeShape_(),
    holeSize_(-1),
    longitudinalPitch_(-1),
    transversePitch_(-1),
    geometryUsed_(false)
{
    // Validate known resistanceType
    const bool knownType =
        resistanceType_ == "constant"
     || resistanceType_ == "perforated"
     || resistanceType_ == "slotted"
     || resistanceType_ == "screen";

    if (!knownType)
    {
        FatalIOErrorInFunction(dict)
            << "Unknown resistanceType \"" << resistanceType_
            << "\". Supported: constant, perforated, slotted, screen"
            << exit(FatalIOError);
    }

    if (resistanceType_ != "constant")
    {
        const word coeffsName = resistanceType_ + "Coeffs";

        if (dict.found(coeffsName))
        {
            // Geometry sub-dict present — compute f from it
            if (dict.found("openAreaFraction"))
            {
                FatalIOErrorInFunction(dict)
                    << "openAreaFraction must not be specified together "
                       "with " << coeffsName
                    << exit(FatalIOError);
            }

            const dictionary& typeDict = dict.subDict(coeffsName);

            if (resistanceType_ == "perforated")
            {
                holeShape_ = typeDict.get<word>("holeShape");
                holeSize_ = typeDict.get<scalar>("holeSize");
                longitudinalPitch_ =
                    typeDict.get<scalar>("longitudinalPitch");
                transversePitch_ =
                    typeDict.get<scalar>("transversePitch");

                if
                (
                    holeShape_ != "circular"
                 && holeShape_ != "square"
                )
                {
                    FatalIOErrorInFunction(dict)
                        << "holeShape must be 'circular' or 'square'"
                        << exit(FatalIOError);
                }

                const scalar Ahole =
                (
                    holeShape_ == "circular"
                  ? constant::mathematical::pi*sqr(holeSize_)/4.0
                  : sqr(holeSize_)
                );

                openAreaFraction_ =
                    Ahole/(longitudinalPitch_*transversePitch_);

                if (openAreaFraction_ <= 0 || openAreaFraction_ > 1)
                {
                    FatalIOErrorInFunction(dict)
                        << "Computed openAreaFraction "
                        << openAreaFraction_
                        << " is not in (0, 1]. "
                           "Check hole geometry and pitch."
                        << exit(FatalIOError);
                }

                geometryUsed_ = true;
            }
            else
            {
                // slotted or screen — placeholder
                WarningInFunction
                    << resistanceType_
                    << " not yet implemented" << endl;
                openAreaFraction_ = 0;
                geometryUsed_ = true;
            }
        }
        else
        {
            // No geometry sub-dict — openAreaFraction required
            if (openAreaFraction_ < 0)
            {
                FatalIOErrorInFunction(dict)
                    << "resistanceType '" << resistanceType_
                    << "' requires either 'openAreaFraction' or '"
                    << coeffsName << "'"
                    << exit(FatalIOError);
            }

            if (openAreaFraction_ <= 0 || openAreaFraction_ > 1)
            {
                FatalIOErrorInFunction(dict)
                    << "openAreaFraction must be in (0, 1]"
                    << exit(FatalIOError);
            }
        }
    }

    if (resistanceType_ == "constant")
    {
        if (lossCoefficient_ < 0)
        {
            FatalIOErrorInFunction(dict)
                << "resistanceType 'constant' requires 'lossCoefficient'"
                << exit(FatalIOError);
        }
    }
    else if (lossCoefficient_ >= 0)
    {
        FatalIOErrorInFunction(dict)
            << "'lossCoefficient' is only valid for "
               "resistanceType 'constant'"
            << exit(FatalIOError);
    }
}


Foam::flowResistanceModel::flowResistanceModel
(
    const flowResistanceModel& m
)
:
    resistanceType_(m.resistanceType_),
    openAreaFraction_(m.openAreaFraction_),
    lossCoefficient_(m.lossCoefficient_),
    holeShape_(m.holeShape_),
    holeSize_(m.holeSize_),
    longitudinalPitch_(m.longitudinalPitch_),
    transversePitch_(m.transversePitch_),
    geometryUsed_(m.geometryUsed_)
{}


Foam::autoPtr<Foam::flowResistanceModel>
Foam::flowResistanceModel::clone() const
{
    return autoPtr<flowResistanceModel>(new flowResistanceModel(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::flowResistanceModel::lossCoefficient() const
{
    scalar K = 0;

    if (lossCoefficient_ >= 0)
    {
        K = lossCoefficient_;
    }
    else if (resistanceType_ == "perforated")
    {
        const scalar f = openAreaFraction_;

        if (f <= 0 || f > 1)
        {
            FatalErrorInFunction
                << "Valid openAreaFraction required for perforated resistance"
                << abort(FatalError);
        }

        K = sqr(0.707*(1 - f))/sqr(f);
    }
    else if (resistanceType_ == "slotted")
    {
        WarningInFunction
            << "Slotted resistance type not yet implemented" << endl;
    }
    else if (resistanceType_ == "screen")
    {
        WarningInFunction
            << "Screen resistance type not yet implemented" << endl;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown resistanceType: " << resistanceType_
            << ". Supported: constant, perforated, slotted, screen"
            << abort(FatalError);
    }

    return K;
}


void Foam::flowResistanceModel::write(Ostream& os) const
{
    os.writeEntry("resistanceType", resistanceType_);

    if (geometryUsed_)
    {
        const word coeffsName = resistanceType_ + "Coeffs";
        os.beginBlock(coeffsName);

        if (resistanceType_ == "perforated")
        {
            os.writeEntry("holeShape", holeShape_);
            os.writeEntry("holeSize", holeSize_);
            os.writeEntry("longitudinalPitch", longitudinalPitch_);
            os.writeEntry("transversePitch", transversePitch_);
        }

        os.endBlock();
    }
    else if (openAreaFraction_ > 0)
    {
        os.writeEntry("openAreaFraction", openAreaFraction_);
    }

    if (resistanceType_ == "constant")
    {
        os.writeEntry("lossCoefficient", lossCoefficient_);
    }
}


// ************************************************************************* //
