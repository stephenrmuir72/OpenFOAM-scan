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
#include "interpolateXY.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scalar Foam::flowResistanceModel::resolveOpenAreaFraction
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

        if (f <= 0 || f > 1)
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

        f = holeArea()/(longitudinalPitch_*transversePitch_);

        if (f <= 0 || f > 1)
        {
            FatalIOErrorInFunction(dict)
                << "Computed openAreaFraction " << f
                << " is not in (0, 1]. "
                   "Check " << geometryDesc << " geometry and pitch."
                << exit(FatalIOError);
        }
    }

    return f;
}


Foam::flowResistanceModel::flowResistanceModel(const dictionary& dict)
:
    resistanceType_(dict.get<word>("resistanceType")),
    openAreaFraction_(dict.getOrDefault<scalar>("openAreaFraction", -1)),
    lossCoefficient_(dict.getOrDefault<scalar>("lossCoefficient", -1)),
    holeShape_(),
    holeSize_(-1),
    longitudinalPitch_(-1),
    transversePitch_(-1),
    slotLength_(-1),
    slotWidth_(-1),
    geometryUsed_(false),
    placeholderCoeffsDict_(),
    relax_(dict.getOrDefault<scalar>("relax", 1)),
    userNu_(dict.getOrDefault<scalar>("nu", -1)),
    targetCorrection_(0),
    relaxedCorrection_(0),
    relaxedCorrectionInitialised_(false),
    patchAreaValidated_(false)
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

    if (resistanceType_ == "perforated")
    {
        const word coeffsName = resistanceType_ + "Coeffs";

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

        openAreaFraction_ =
            resolveOpenAreaFraction(dict, typeDict, "hole");

        geometryUsed_ = true;
    }
    else if (resistanceType_ == "slotted")
    {
        const word coeffsName = resistanceType_ + "Coeffs";

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

        if (slotWidth_ <= 0)
        {
            FatalIOErrorInFunction(dict)
                << "slotWidth must be positive"
                << exit(FatalIOError);
        }

        if (slotLength_ <= slotWidth_)
        {
            FatalIOErrorInFunction(dict)
                << "slotLength must be greater than slotWidth (a > b) "
                   "to give a positive straight section"
                << exit(FatalIOError);
        }

        openAreaFraction_ =
            resolveOpenAreaFraction(dict, typeDict, "slot");

        geometryUsed_ = true;
    }
    else if (resistanceType_ != "constant")
    {
        // screen — placeholder
        const word coeffsName = resistanceType_ + "Coeffs";

        if (dict.found(coeffsName))
        {
            // Geometry sub-dict present
            if (dict.found("openAreaFraction"))
            {
                FatalIOErrorInFunction(dict)
                    << "openAreaFraction must not be specified together "
                       "with " << coeffsName
                    << exit(FatalIOError);
            }

            WarningInFunction
                << resistanceType_
                << " not yet implemented" << endl;
            openAreaFraction_ = 0;
            geometryUsed_ = true;

            // Preserve the user-supplied sub-dict verbatim so that it
            // is not silently dropped when write() is called (e.g. on
            // a restart write of the boundary condition)
            placeholderCoeffsDict_ = dict.subDict(coeffsName);
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
    resistanceType_(m.resistanceType_),
    openAreaFraction_(m.openAreaFraction_),
    lossCoefficient_(m.lossCoefficient_),
    holeShape_(m.holeShape_),
    holeSize_(m.holeSize_),
    longitudinalPitch_(m.longitudinalPitch_),
    transversePitch_(m.transversePitch_),
    slotLength_(m.slotLength_),
    slotWidth_(m.slotWidth_),
    geometryUsed_(m.geometryUsed_),
    placeholderCoeffsDict_(m.placeholderCoeffsDict_),
    relax_(m.relax_),
    userNu_(m.userNu_),
    targetCorrection_(m.targetCorrection_),
    relaxedCorrection_(m.relaxedCorrection_),
    relaxedCorrectionInitialised_(m.relaxedCorrectionInitialised_),
    patchAreaValidated_(m.patchAreaValidated_)
{}


Foam::autoPtr<Foam::flowResistanceModel>
Foam::flowResistanceModel::clone() const
{
    return autoPtr<flowResistanceModel>(new flowResistanceModel(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::flowResistanceModel::holeArea() const
{
    if (resistanceType_ == "perforated")
    {
        return
        (
            holeShape_ == "circular"
          ? constant::mathematical::pi*sqr(holeSize_)/4.0
          : sqr(holeSize_)
        );
    }
    else if (resistanceType_ == "slotted")
    {
        return
            slotWidth_*(slotLength_ - slotWidth_)
          + constant::mathematical::pi*sqr(slotWidth_)/4.0;
    }

    FatalErrorInFunction
        << "No single-hole area for resistanceType '" << resistanceType_
        << "'"
        << abort(FatalError);

    return 0;
}


Foam::scalar Foam::flowResistanceModel::holePerimeter() const
{
    if (resistanceType_ == "perforated")
    {
        return
        (
            holeShape_ == "circular"
          ? constant::mathematical::pi*holeSize_
          : 4.0*holeSize_
        );
    }
    else if (resistanceType_ == "slotted")
    {
        return 2.0*(slotLength_ - slotWidth_)
             + constant::mathematical::pi*slotWidth_;
    }

    FatalErrorInFunction
        << "No single-hole perimeter for resistanceType '"
        << resistanceType_ << "'"
        << abort(FatalError);

    return 0;
}


Foam::scalar Foam::flowResistanceModel::hydraulicDiameter() const
{
    return 4.0*holeArea()/holePerimeter();
}


Foam::scalar Foam::flowResistanceModel::zetaTurbulent() const
{
    const scalar f = openAreaFraction_;
    const scalar zetaT = sqr(0.707*sqrt(1 - f) + (1 - f))/sqr(f);

    if (resistanceType_ == "slotted")
    {
        const scalar ab = slotLength_/slotWidth_;

        const scalar kShape =
        (
            ab > 10
          ? 0.81
          : 1 - 0.05*pow(ab - 1, 0.6)
        );

        return zetaT*kShape;
    }

    return zetaT;
}


Foam::scalar Foam::flowResistanceModel::eps0(const scalar Re) const
{
    static const scalarList ReTable
    {
        10, 15, 20, 25, 40, 60, 100, 200, 400, 1000, 2000, 4000,
        10000, 20000, 100000, 200000, 1000000
    };

    static const scalarList epsTable
    {
        0.32, 0.33, 0.34, 0.34, 0.36, 0.37, 0.4, 0.42, 0.46, 0.53,
        0.59, 0.64, 0.74, 0.81, 0.94, 0.96, 0.98
    };

    // log(Re) table, precomputed once, for the log-linear interpolation
    static const scalarList logReTable
    (
        []()
        {
            scalarList logRe(ReTable.size());
            forAll(ReTable, i)
            {
                logRe[i] = log(ReTable[i]);
            }
            return logRe;
        }()
    );

    // Clamp to the table range, then log-linearly interpolate
    // (interpolateXY holds the boundary value beyond the table range)
    const scalar ReC = max(min(Re, ReTable.last()), ReTable.first());

    return interpolateXY(log(ReC), logReTable, epsTable);
}


Foam::scalar Foam::flowResistanceModel::zetaPhi
(
    const scalar Re,
    const scalar f
) const
{
    static const scalar ReTable[] =
    {
        25, 40, 60, 100, 200, 400, 1000, 2000, 4000, 10000, 20000,
        100000, 200000
    };

    static const scalar fTable[] =
    {
        0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0
    };

    static const scalar zetaTable[13][11] =
    {
        {1.94, 1.78, 1.57, 1.35, 1.1, 0.85, 0.58, 0.4, 0.2, 0.03, 0},
        {1.38, 1.36, 1.16, 0.99, 0.75, 0.56, 0.37, 0.24, 0.13, 0.03, 0},
        {1.14, 1.05, 0.88, 0.79, 0.55, 0.3, 0.23, 0.13, 0.08, 0.02, 0},
        {0.89, 0.85, 0.75, 0.57, 0.34, 0.19, 0.11, 0.06, 0.03, 0, 0},
        {0.69, 0.67, 0.57, 0.4, 0.19, 0.1, 0.06, 0.03, 0.01, 0, 0},
        {0.64, 0.57, 0.43, 0.28, 0.12, 0.06, 0.03, 0.02, 0, 0, 0},
        {0.39, 0.36, 0.3, 0.19, 0.07, 0.03, 0.02, 0.01, 0, 0, 0},
        {0.3, 0.26, 0.22, 0.14, 0.05, 0.02, 0.01, 0, 0, 0, 0},
        {0.22, 0.2, 0.17, 0.1, 0.03, 0.01, 0, 0, 0, 0, 0},
        {0.15, 0.13, 0.1, 0.06, 0.02, 0.01, 0, 0, 0, 0, 0},
        {0.11, 0.09, 0.07, 0.04, 0.01, 0, 0, 0, 0, 0, 0},
        {0.04, 0.03, 0.02, 0.02, 0.01, 0, 0, 0, 0, 0, 0},
        {0.01, 0.01, 0.01, 0.01, 0.01, 0, 0, 0, 0, 0, 0}
    };

    const label nRe = 13;
    const label nf = 11;

    // Clamp to the table range (the caller only uses Re <= 1e5)
    const scalar ReC = max(min(Re, ReTable[nRe - 1]), ReTable[0]);
    const scalar fC = max(min(f, fTable[nf - 1]), fTable[0]);

    // Re bracket (log-linear) and f bracket (linear)
    label i = 0;
    while (i < nRe - 2 && ReC >= ReTable[i + 1])
    {
        ++i;
    }

    label j = 0;
    while (j < nf - 2 && fC >= fTable[j + 1])
    {
        ++j;
    }

    const scalar tRe =
    (
        (log(ReC) - log(ReTable[i]))
      / (log(ReTable[i + 1]) - log(ReTable[i]))
    );

    const scalar tf =
    (
        (fC - fTable[j])/(fTable[j + 1] - fTable[j])
    );

    const scalar zeta00 = zetaTable[i][j];
    const scalar zeta01 = zetaTable[i][j + 1];
    const scalar zeta10 = zetaTable[i + 1][j];
    const scalar zeta11 = zetaTable[i + 1][j + 1];

    // Linear in f at each Re row, then log-linear in Re
    const scalar zeta0 = zeta00 + tf*(zeta01 - zeta00);
    const scalar zeta1 = zeta10 + tf*(zeta11 - zeta10);

    return zeta0 + tRe*(zeta1 - zeta0);
}


Foam::scalar Foam::flowResistanceModel::lossCoefficient
(
    const scalar UnAvg,
    const scalar nuAvg
) const
{
    if (lossCoefficient_ >= 0)
    {
        return lossCoefficient_;
    }

    if (resistanceType_ == "perforated" || resistanceType_ == "slotted")
    {
        const scalar f = openAreaFraction_;

        if (f <= 0 || f > 1)
        {
            FatalErrorInFunction
                << "Valid openAreaFraction required for resistance"
                << abort(FatalError);
        }

        // Hole Reynolds number, clamped to a minimum of 1 to keep the
        // laminar term finite at (near) zero velocity
        const scalar Re =
            max(mag(UnAvg)*hydraulicDiameter()/(f*nuAvg), scalar(1));

        const scalar zetaT = zetaTurbulent();

        if (Re < 10)
        {
            return 33/Re/sqr(f);
        }
        else if (Re < 25)
        {
            return 33/Re/sqr(f) + eps0(Re)*zetaT;
        }
        else if (Re <= 1e5)
        {
            return zetaPhi(Re, f)/sqr(f) + eps0(Re)*zetaT;
        }
        else
        {
            return zetaT;
        }
    }
    else if (resistanceType_ == "screen")
    {
        WarningInFunction
            << "Screen resistance type not yet implemented" << endl;

        return 0;
    }

    FatalErrorInFunction
        << "Unknown resistanceType: " << resistanceType_
        << ". Supported: constant, perforated, slotted, screen"
        << abort(FatalError);

    return 0;
}


Foam::scalar Foam::flowResistanceModel::updateCorrection
(
    const scalar Q,
    const scalar A,
    const scalar rhoAvg,
    const scalar nuAvg
)
{
    // Validate the hole count against the patch area on the first update
    if (!patchAreaValidated_)
    {
        if (resistanceType_ == "perforated" || resistanceType_ == "slotted")
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

        patchAreaValidated_ = true;
    }

    const scalar UnAvg = Q/A;
    const scalar deltaP =
        lossCoefficient(UnAvg, nuAvg)*0.5*rhoAvg*sqr(UnAvg);

    // Signed target: inflow (UnAvg < 0) lowers p0, outflow raises it
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
            (1 - relax_)*relaxedCorrection_
          + relax_*targetCorrection_;
    }

    return relaxedCorrection_;
}


void Foam::flowResistanceModel::write(Ostream& os) const
{
    os.writeEntry("resistanceType", resistanceType_);

    if (geometryUsed_)
    {
        const word coeffsName = resistanceType_ + "Coeffs";

        if (resistanceType_ == "perforated" || resistanceType_ == "slotted")
        {
            os.beginBlock(coeffsName);

            if (resistanceType_ == "perforated")
            {
                os.writeEntry("holeShape", holeShape_);
                os.writeEntry("holeSize", holeSize_);
            }
            else
            {
                os.writeEntry("slotLength", slotLength_);
                os.writeEntry("slotWidth", slotWidth_);
            }

            if (longitudinalPitch_ > 0)
            {
                os.writeEntry("longitudinalPitch", longitudinalPitch_);
                os.writeEntry("transversePitch", transversePitch_);
            }
            else if (openAreaFraction_ > 0)
            {
                os.writeEntry("openAreaFraction", openAreaFraction_);
            }

            os.endBlock();
        }
        else
        {
            // Placeholder resistance type (e.g. 'screen'): re-emit the
            // user-supplied sub-dict verbatim rather than reconstructing
            // it from (currently unused) member data
            os.writeKeyword(coeffsName);
            placeholderCoeffsDict_.write(os, true);
        }
    }
    else if (openAreaFraction_ > 0)
    {
        os.writeEntry("openAreaFraction", openAreaFraction_);
    }

    if (resistanceType_ == "constant")
    {
        os.writeEntry("lossCoefficient", lossCoefficient_);
    }

    os.writeEntryIfDifferent<scalar>("relax", 1, relax_);

    os.writeEntryIfDifferent<scalar>("nu", -1, userNu_);

    if (relax_ < 1 && relaxedCorrectionInitialised_)
    {
        os.writeEntry("relaxedPressureCorrection", relaxedCorrection_);
    }
}


// ************************************************************************* //
