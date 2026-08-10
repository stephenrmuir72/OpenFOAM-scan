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

#include "holeFlowResistance.H"
#include "interpolateXY.H"
#include "scalarList.H"

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

    openAreaFraction_ = f;

    return f;
}


Foam::scalar
Foam::flowResistanceModels::holeFlowResistance::hydraulicDiameter() const
{
    return 4.0*holeArea()/holePerimeter();
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
    const scalar zetaT = sqr(0.707*sqrt(1 - f) + (1 - f))/sqr(f);

    return zetaT*shapeFactor();
}


Foam::scalar
Foam::flowResistanceModels::holeFlowResistance::eps0(const scalar Re)
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


Foam::scalar
Foam::flowResistanceModels::holeFlowResistance::zetaPhi
(
    const scalar Re,
    const scalar f
)
{
    static const scalar ReTable[] =
    {
        25, 40, 60, 100, 200, 400, 1000, 2000, 4000, 10000, 20000,
        100000, 200000
    };

    static const scalarList fTable
    {
        0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0
    };

    static const List<scalarList> zetaTable
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

    // Clamp to the table range (the caller only uses Re <= 1e5)
    const scalar ReC = max(min(Re, ReTable[nRe - 1]), ReTable[0]);
    const scalar fC = max(min(f, fTable.last()), fTable.first());

    // Re bracket (log-linear); the f interpolation of each bracketing
    // row is delegated to interpolateXY (linear in f)
    label i = 0;
    while (i < nRe - 2 && ReC >= ReTable[i + 1])
    {
        ++i;
    }

    const scalar tRe =
    (
        (log(ReC) - log(ReTable[i]))
      / (log(ReTable[i + 1]) - log(ReTable[i]))
    );

    const scalar zeta0 = interpolateXY(fC, fTable, zetaTable[i]);
    const scalar zeta1 = interpolateXY(fC, fTable, zetaTable[i + 1]);

    // Linear in f at each Re row (via interpolateXY above), then
    // log-linear in Re
    return zeta0 + tRe*(zeta1 - zeta0);
}


Foam::scalar
Foam::flowResistanceModels::holeFlowResistance::lossCoefficient
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
