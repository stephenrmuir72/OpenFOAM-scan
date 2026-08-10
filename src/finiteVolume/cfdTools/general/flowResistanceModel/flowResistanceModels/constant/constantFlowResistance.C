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

#include "constantFlowResistance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * //

namespace Foam
{
namespace flowResistanceModels
{
    defineTypeNameAndDebug(Constant, 0);
    addToRunTimeSelectionTable(flowResistanceModel, Constant, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowResistanceModels::Constant::Constant(const dictionary& dict)
:
    flowResistanceModel(dict),
    lossCoefficient_(dict.get<scalar>("lossCoefficient"))
{
    if (lossCoefficient_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "lossCoefficient must not be negative"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::flowResistanceModels::Constant::lossCoefficient
(
    const scalar UnAvg,
    const scalar nuAvg
) const
{
    return lossCoefficient_;
}


void Foam::flowResistanceModels::Constant::writeCoeffs(Ostream& os) const
{
    os.writeEntry("lossCoefficient", lossCoefficient_);
}


// ************************************************************************* //
