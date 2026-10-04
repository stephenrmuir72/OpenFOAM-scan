/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Keysight Technologies
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "custom.H"
#include "Function1.H"
#include "addToRunTimeSelectionTable.H"

#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * //

defineTypeNameAndDebug(custom, 0);
addToRunTimeSelectionTable(temporalKernel, custom, dictionary);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

custom::custom(const dictionary& dict)
:
    temporalKernel(typeName, dict),
    releaseRatePtr_(Function1<scalar>::New("releaseRate", dict)),
    totalIntegral_(releaseRatePtr_->integrate(0, 1))
{
    if (!(totalIntegral_ > 0) || !std::isfinite(totalIntegral_))
    {
        FatalIOErrorInFunction(dict_)
            << "The integral of releaseRate over normalized age [0, 1] must "
            << "be positive and finite, received " << totalIntegral_
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

scalar custom::cumulative
(
    const scalar age,
    const scalar duration
) const
{
    if (!(duration > VSMALL) || !std::isfinite(duration))
    {
        FatalErrorInFunction
            << "Duration must be positive and finite, received " << duration
            << abort(FatalError);
    }

    if (!std::isfinite(age))
    {
        FatalErrorInFunction
            << "Age must be finite, received " << age
            << abort(FatalError);
    }

    if (age <= 0)
    {
        return 0;
    }

    if (age >= duration)
    {
        return 1;
    }

    const scalar x = age/duration;
    const scalar fraction = releaseRatePtr_->integrate(0, x)/totalIntegral_;

    if
    (
        !std::isfinite(fraction)
     || fraction < 0
     || fraction > 1
    )
    {
        FatalIOErrorInFunction(dict_)
            << "Invalid cumulative fraction from releaseRate at normalized "
            << "age " << x << ": " << fraction
            << exit(FatalIOError);
    }

    return fraction;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
