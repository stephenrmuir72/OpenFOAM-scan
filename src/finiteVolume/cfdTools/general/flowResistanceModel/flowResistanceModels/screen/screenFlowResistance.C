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
    openAreaFraction_(dict.getOrDefault<scalar>("openAreaFraction", -1)),
    placeholderCoeffsDict_(),
    geometryUsed_(false)
{
    const word coeffsName = "screenCoeffs";

    if (dict.found(coeffsName))
    {
        if (dict.found("openAreaFraction"))
        {
            FatalIOErrorInFunction(dict)
                << "openAreaFraction must not be specified together "
                   "with " << coeffsName
                << exit(FatalIOError);
        }

        WarningInFunction
            << "screen resistanceType not yet implemented" << endl;

        geometryUsed_ = true;

        // Preserve the user-supplied sub-dict verbatim so that it is
        // not silently dropped when write() is called (e.g. on a
        // restart write of the boundary condition)
        placeholderCoeffsDict_ = dict.subDict(coeffsName);
    }
    else
    {
        // No geometry sub-dict — openAreaFraction required
        if (openAreaFraction_ < 0)
        {
            FatalIOErrorInFunction(dict)
                << "resistanceType 'screen' requires either "
                   "'openAreaFraction' or '" << coeffsName << "'"
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::flowResistanceModels::Screen::lossCoefficient
(
    const scalar UnAvg,
    const scalar nuAvg
) const
{
    return 0;
}


void Foam::flowResistanceModels::Screen::writeCoeffs(Ostream& os) const
{
    if (geometryUsed_)
    {
        // Re-emit the user-supplied sub-dict verbatim rather than
        // reconstructing it from (currently unused) member data
        os.writeKeyword("screenCoeffs");
        placeholderCoeffsDict_.write(os, true);
    }
    else if (openAreaFraction_ > 0)
    {
        os.writeEntry("openAreaFraction", openAreaFraction_);
    }
}


// ************************************************************************* //
