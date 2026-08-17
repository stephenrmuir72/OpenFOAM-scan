/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
      \\/     M anipulation  |
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

#include "morphologyBasedAggregationScale.H"
#include "addToRunTimeSelectionTable.H"
#include "morphologyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * //

defineTypeNameAndDebug(morphologyBasedAggregationScale, 0);
addToRunTimeSelectionTable
(
    aggregationScaleModel,
    morphologyBasedAggregationScale,
    dictionary
);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

morphologyBasedAggregationScale::morphologyBasedAggregationScale
(
    const liquidFilmBase& film,
    const dictionary& dict
)
:
    aggregationScaleModel(type(), film, dict),
    Ccapillary_(dict.getOrDefault<scalar>("Ccapillary", 1.0)),
    Cligament_(dict.getOrDefault<scalar>("Cligament", 2.915)),
    Csheet_(dict.getOrDefault<scalar>("Csheet", 1.0)),
    Ccoherence_(dict.getOrDefault<scalar>("Ccoherence", 1.0))
{
    if (Ccapillary_ <= VSMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Ccapillary must be positive" << exit(FatalIOError);
    }

    if (Cligament_ <= VSMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Cligament must be positive" << exit(FatalIOError);
    }

    if (Csheet_ <= VSMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Csheet must be positive" << exit(FatalIOError);
    }

    if (Ccoherence_ <= VSMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Ccoherence must be positive" << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Static Member Functions  * * * * * * * * * //

scalar morphologyBasedAggregationScale::inertioCapillaryTime
(
    const scalar rho,
    const scalar sigma,
    const scalar diameter,
    const scalar coefficient
)
{
    if (rho <= VSMALL || sigma <= VSMALL || diameter <= VSMALL)
    {
        FatalErrorInFunction
            << "Inertio-capillary time requires positive rho, sigma, and "
            << "diameter. Received rho=" << rho
            << ", sigma=" << sigma
            << ", diameter=" << diameter
            << abort(FatalError);
    }

    return coefficient*Foam::sqrt(rho*Foam::pow(diameter, 3)/sigma);
}


scalar morphologyBasedAggregationScale::ligamentBreakupTime
(
    const scalar rho,
    const scalar mu,
    const scalar sigma,
    const scalar diameter,
    const scalar coefficient
)
{
    if (mu < 0)
    {
        FatalErrorInFunction
            << "Ligament breakup time requires non-negative viscosity. "
            << "Received mu=" << mu << abort(FatalError);
    }

    const scalar tauCapillary =
        inertioCapillaryTime(rho, sigma, diameter, coefficient);
    const scalar Oh =
        mu/Foam::sqrt(rho*sigma*diameter);

    return tauCapillary*Foam::pow(1.0 + 3.0*Oh, 1.0/6.0);
}


scalar morphologyBasedAggregationScale::sheetBreakupTime
(
    const scalar diameter,
    const scalar relativeSpeed,
    const scalar coefficient
)
{
    if (diameter <= VSMALL || relativeSpeed <= VSMALL)
    {
        FatalErrorInFunction
            << "Sheet breakup time requires positive diameter and relative "
            << "primary/film speed. Received diameter=" << diameter
            << ", relativeSpeed=" << relativeSpeed
            << abort(FatalError);
    }

    return coefficient*diameter/relativeSpeed;
}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

void morphologyBasedAggregationScale::calculate
(
    const labelField& morphology,
    const scalarField& physicalDiameter,
    const vectorField& Up,
    const boolList& activeFace,
    scalarField& eventTime,
    scalarField& coherenceLength
) const
{
    const label nFaces = mesh_.nFaces();

    if
    (
        morphology.size() != nFaces
     || physicalDiameter.size() != nFaces
     || Up.size() != nFaces
     || activeFace.size() != nFaces
     || eventTime.size() != nFaces
     || coherenceLength.size() != nFaces
    )
    {
        FatalErrorInFunction
            << "All aggregation scale fields must have finite-area face size "
            << nFaces << abort(FatalError);
    }

    const scalarField& rho = film_.rho().primitiveField();
    const scalarField& mu = film_.mu().primitiveField();
    const scalarField& sigma = film_.sigma().primitiveField();
    const vectorField& Uf = film_.Uf().primitiveField();

    eventTime = GREAT;
    coherenceLength = 0;

    forAll(morphology, facei)
    {
        if (!activeFace[facei])
        {
            continue;
        }

        if
        (
            morphology[facei] < 0
         || morphology[facei] >= morphologyModel::nMorphologyTypes
        )
        {
            FatalErrorInFunction
                << "Invalid active morphology ID " << morphology[facei]
                << " on finite-area face " << facei
                << abort(FatalError);
        }

        const scalar diameter = physicalDiameter[facei];

        if (diameter <= VSMALL)
        {
            FatalErrorInFunction
                << "Active aggregation face " << facei
                << " has non-positive physical diameter " << diameter
                << abort(FatalError);
        }

        if (rho[facei] <= VSMALL)
        {
            FatalErrorInFunction
                << "Active aggregation face " << facei
                << " has non-positive liquid density " << rho[facei]
                << abort(FatalError);
        }

        if (mu[facei] < 0)
        {
            FatalErrorInFunction
                << "Active aggregation face " << facei
                << " has negative liquid viscosity " << mu[facei]
                << abort(FatalError);
        }

        if (sigma[facei] <= VSMALL)
        {
            FatalErrorInFunction
                << "Active aggregation face " << facei
                << " has non-positive surface tension " << sigma[facei]
                << abort(FatalError);
        }

        const scalar magUrel = mag(Up[facei] - Uf[facei]);
        scalar tau = 0;

        switch
        (
            static_cast<morphologyModel::morphologyType>(morphology[facei])
        )
        {
            case morphologyModel::DIRECT_DROPLET:
            case morphologyModel::DRIP:
                tau =
                    inertioCapillaryTime
                    (
                        rho[facei],
                        sigma[facei],
                        diameter,
                        Ccapillary_
                    );
                break;

            case morphologyModel::LIGAMENT:
                tau =
                    ligamentBreakupTime
                    (
                        rho[facei],
                        mu[facei],
                        sigma[facei],
                        diameter,
                        Cligament_
                    );
                break;

            case morphologyModel::SHEET:
                tau =
                    sheetBreakupTime
                    (
                        diameter,
                        magUrel,
                        Csheet_
                    );
                break;

            default:
                FatalErrorInFunction
                    << "Unsupported active morphology ID "
                    << morphology[facei] << abort(FatalError);
        }

        const scalar transportSpeed = max(mag(Uf[facei]), magUrel);

        eventTime[facei] = tau;
        coherenceLength[facei] =
            Ccoherence_*max(diameter, transportSpeed*tau);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
