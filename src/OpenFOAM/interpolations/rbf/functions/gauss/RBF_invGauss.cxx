/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2009 Frank Bos, TU Delft
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

#include "RBF_invGauss.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBFTypes
{

defineTypeName(invGauss);
addToRunTimeSelectionTable(RBFFunction, invGauss, dictionary);

} // End namespace RBFTypes
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBFTypes::invGauss::invGauss(const dictionary& dict)
:
    RBFFunction(dict.get<scalar>("radius"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RBFTypes::invGauss::weights
(
    UList<scalar>& wts,
    const UList<point>& points,
    const point& refPoint
) const
{
    RBFFunction::checkFieldSizes(wts, points, FUNCTION_NAME);

    const label loop_len = wts.size();

    for (label i = 0; i < loop_len; ++i)
    {
        wts[i] = weight(points[i], refPoint);
    }
}


// ************************************************************************* //
