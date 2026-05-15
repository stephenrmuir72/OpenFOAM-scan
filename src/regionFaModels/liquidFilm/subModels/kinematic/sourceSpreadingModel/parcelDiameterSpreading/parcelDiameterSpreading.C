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

#include "parcelDiameterSpreading.H"
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

defineTypeNameAndDebug(parcelDiameterSpreading, 0);
addToRunTimeSelectionTable
(
    sourceSpreadingModel,
    parcelDiameterSpreading,
    dictionary
);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
void parcelDiameterSpreading::spreadImpl
(
    Field<Type>& field,
    const scalarField& diameter
) const
{
    const labelListList& faceFaces = mesh_.patch().faceFaces();
    const scalarField& S = mesh_.S().field();
    const vectorField& Cf = mesh_.areaCentres().primitiveField();

    Field<Type> spreadField(field.size(), Zero);

    forAll(field, facei)
    {
        if (mag(field[facei]) > SMALL && diameter[facei] > SMALL)
        {
            const scalar radius = 0.5*diameter[facei];
            const scalar radiusSqr = sqr(radius);
            const vector& centre = Cf[facei];

            // BFS outward, collecting faces within radius
            labelHashSet stencil;
            stencil.insert(facei);

            labelHashSet currentRing;
            currentRing.insert(facei);

            bool expanding = true;
            while (expanding)
            {
                expanding = false;
                labelHashSet nextRing;

                for (const label fi : currentRing)
                {
                    for (const label nbr : faceFaces[fi])
                    {
                        if (stencil.insert(nbr))
                        {
                            if (magSqr(Cf[nbr] - centre) <= radiusSqr)
                            {
                                nextRing.insert(nbr);
                                expanding = true;
                            }
                        }
                    }
                }

                currentRing = std::move(nextRing);
            }

            // Filter stencil to faces within radius
            labelHashSet filtered;
            for (const label fi : stencil)
            {
                if (magSqr(Cf[fi] - centre) <= radiusSqr)
                {
                    filtered.insert(fi);
                }
            }

            // Area-weighted distribution
            scalar totalArea = 0;
            for (const label fi : filtered)
            {
                totalArea += S[fi];
            }

            if (totalArea > SMALL)
            {
                const Type totalSource = field[facei];
                for (const label fi : filtered)
                {
                    spreadField[fi] += totalSource * (S[fi] / totalArea);
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

parcelDiameterSpreading::parcelDiameterSpreading
(
    const faMesh& mesh,
    const dictionary& dict
)
:
    sourceSpreadingModel(type(), mesh, dict)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void parcelDiameterSpreading::spread
(
    scalarField& field,
    const scalarField& diameter
) const
{
    spreadImpl(field, diameter);
}


void parcelDiameterSpreading::spread
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
