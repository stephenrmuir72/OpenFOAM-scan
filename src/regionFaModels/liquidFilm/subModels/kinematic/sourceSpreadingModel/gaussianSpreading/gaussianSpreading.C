/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 OpenCFD Ltd.
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

#include "gaussianSpreading.H"
#include "addToRunTimeSelectionTable.H"
#include "HashSet.H"
#include "areaFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gaussianSpreading, 0);
addToRunTimeSelectionTable
(
    sourceSpreadingModel,
    gaussianSpreading,
    dictionary
);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
void gaussianSpreading::spreadImpl
(
    Field<Type>& field,
    const scalarField& diameter
) const
{
    const labelListList& faceFaces = mesh_.patch().faceFaces();
    const scalarField& S = mesh_.S().field();
    const vectorField& Cf = mesh_.areaCentres().primitiveField();

    const scalar invTwoSigmaSqr = 1.0/(2.0*sqr(sigma_));

    Field<Type> spreadField(field.size(), Zero);

    forAll(field, facei)
    {
        if (mag(field[facei]) > SMALL)
        {
            const vector& centre = Cf[facei];

            // BFS outward for nLayers rings
            labelHashSet stencil;
            stencil.insert(facei);

            labelHashSet currentRing;
            currentRing.insert(facei);

            for (label layer = 0; layer < nLayers_; ++layer)
            {
                labelHashSet nextRing;

                for (const label fi : currentRing)
                {
                    for (const label nbr : faceFaces[fi])
                    {
                        if (stencil.insert(nbr))
                        {
                            nextRing.insert(nbr);
                        }
                    }
                }

                currentRing = std::move(nextRing);
            }

            // Compute Gaussian weights
            scalar totalWeight = 0;
            for (const label fi : stencil)
            {
                const scalar distSqr = magSqr(Cf[fi] - centre);
                totalWeight += S[fi] * Foam::exp(-distSqr*invTwoSigmaSqr);
            }

            if (totalWeight > SMALL)
            {
                const Type totalSource = field[facei];
                for (const label fi : stencil)
                {
                    const scalar distSqr = magSqr(Cf[fi] - centre);
                    const scalar w =
                        S[fi] * Foam::exp(-distSqr*invTwoSigmaSqr);
                    spreadField[fi] += totalSource * (w / totalWeight);
                }
            }
            else
            {
                spreadField[facei] += field[facei];
            }
        }
        else
        {
            spreadField[facei] += field[facei];
        }
    }

    field = spreadField;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

gaussianSpreading::gaussianSpreading
(
    const faMesh& mesh,
    const dictionary& dict
)
:
    sourceSpreadingModel(type(), mesh, dict),
    sigma_(dict.get<scalar>("sigma")),
    nLayers_(dict.getOrDefault<label>("nLayers", 3))
{
    if (sigma_ <= 0)
    {
        FatalErrorInFunction
            << "sigma must be positive, received " << sigma_
            << abort(FatalError);
    }

    if (nLayers_ < 1)
    {
        FatalErrorInFunction
            << "nLayers must be >= 1, received " << nLayers_
            << abort(FatalError);
    }

    Info<< "    Gaussian source spreading:"
        << " sigma=" << sigma_
        << " nLayers=" << nLayers_
        << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void gaussianSpreading::spread
(
    scalarField& field,
    const scalarField& diameter
) const
{
    spreadImpl(field, diameter);
}


void gaussianSpreading::spread
(
    vectorField& field,
    const scalarField& diameter
) const
{
    spreadImpl(field, diameter);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
