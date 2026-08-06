/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2024 OpenCFD Ltd.
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

#include "filmSeparation.H"
#include "filmSeparationModel.H"
#include "separationAggregationModel.H"
#include "liquidFilmBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Local Helper Functions  * * * * * * * * * * * * * //

namespace
{

//- Build the IOobject for a restartable accumulation buffer.  When the buffer
//  is active it auto-writes and is read on restart; otherwise it is an
//  unregistered, non-writing placeholder (the buffers are unused when
//  accumulation is disabled).
Foam::IOobject accumulatorIO
(
    const Foam::word& name,
    const Foam::regionModels::areaSurfaceFilmModels::liquidFilmBase& film,
    const bool active
)
{
    return Foam::IOobject
    (
        name,
        film.regionMesh().time().timeName(),
        film.regionMesh().thisDb(),
        active ? Foam::IOobject::READ_IF_PRESENT : Foam::IOobject::NO_READ,
        active ? Foam::IOobject::AUTO_WRITE : Foam::IOobject::NO_WRITE,
        active ? Foam::IOobject::REGISTER : Foam::IOobject::NO_REGISTER
    );
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(filmSeparation, 0);
addToRunTimeSelectionTable
(
    injectionModel,
    filmSeparation,
    dictionary
);
}
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::vector
Foam::regionModels::areaSurfaceFilmModels::filmSeparation::randomDirection()
{
    // Three independent normal deviates give a direction uniformly distributed
    // on the unit sphere; the previous uniform-in-cube sample biased the jitter
    // towards the cube diagonals.
    vector dir
    (
        rndGen_.GaussNormal<scalar>(),
        rndGen_.GaussNormal<scalar>(),
        rndGen_.GaussNormal<scalar>()
    );

    const scalar magDir = mag(dir);

    return (magDir > VSMALL ? dir/magDir : vector(1, 0, 0));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::areaSurfaceFilmModels::filmSeparation::filmSeparation
(
    liquidFilmBase& film,
    const dictionary& dict
)
:
    injectionModel(type(), film, dict),
    filmSeparationModelPtr_(filmSeparationModel::New(film, coeffDict_)),
    separationAggregationModelPtr_
    (
        separationAggregationModel::New(film, coeffDict_)
    ),
    massAccumulator_
    (
        accumulatorIO
        (
            "massAccumulator",
            film,
            separationAggregationModelPtr_.valid()
        ),
        scalarField
        (
            separationAggregationModelPtr_.valid()
          ? film.regionMesh().nFaces()
          : 0,
            Zero
        )
    ),
    momentumAccumulator_
    (
        accumulatorIO
        (
            "momentumAccumulator",
            film,
            separationAggregationModelPtr_.valid()
        ),
        vectorField
        (
            separationAggregationModelPtr_.valid()
          ? film.regionMesh().nFaces()
          : 0,
            Zero
        )
    ),
    diameterAccumulator_
    (
        accumulatorIO
        (
            "diameterAccumulator",
            film,
            separationAggregationModelPtr_.valid()
        ),
        scalarField
        (
            separationAggregationModelPtr_.valid()
          ? film.regionMesh().nFaces()
          : 0,
            Zero
        )
    ),
    accumulationTime_
    (
        accumulatorIO
        (
            "accumulationTime",
            film,
            separationAggregationModelPtr_.valid()
        ),
        scalarField
        (
            separationAggregationModelPtr_.valid()
          ? film.regionMesh().nFaces()
          : 0,
            Zero
        )
    ),
    velocityJitter_
    (
        coeffDict_.getOrDefault<scalar>("velocityJitter", 0.0)
    ),
    rndGen_(coeffDict_.getOrDefault<label>("seed", 0) + Pstream::myProcNo()),
    UfInjectionPtr_(nullptr)
{
    if (separationAggregationModelPtr_)
    {
        // Register a velocity field on the FA mesh for Lagrangian reading
        UfInjectionPtr_.reset
        (
            new areaVectorField
            (
                IOobject
                (
                    "UfInjection",
                    film.regionMesh().time().timeName(),
                    film.regionMesh().thisDb(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::REGISTER
                ),
                film.regionMesh(),
                dimensionedVector(dimVelocity, Zero)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionModels::areaSurfaceFilmModels::filmSeparation::~filmSeparation()
{}  // filmSeparationModel was forward declared


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::regionModels::areaSurfaceFilmModels::filmSeparation::correct
(
    scalarField& availableMass,
    scalarField& massToInject,
    scalarField& diameterToInject
)
{
    const faMesh& mesh = film().regionMesh();

    // Calculate the mass ratio of film separation
    tmp<scalarField> tmassRatio = filmSeparationModelPtr_->separatedMassRatio();
    const auto& massRatio = tmassRatio.cref();

    if (!separationAggregationModelPtr_.valid())
    {
        massToInject = massRatio*availableMass;
        diameterToInject = massRatio*film().h();
        availableMass -= massToInject;

        addToInjectedMass(sum(massToInject));
    }
    else
    {
        const scalar deltaT = mesh.time().deltaTValue();
        const areaVectorField& Uf = film().Uf();

        scalarField effectiveAvailMass
        (
            max(availableMass - massAccumulator_, scalar(0))
        );

        scalarField separatedMass(massRatio*effectiveAvailMass);

        forAll(separatedMass, faceI)
        {
            if (separatedMass[faceI] > VSMALL)
            {
                massAccumulator_[faceI] += separatedMass[faceI];
                momentumAccumulator_[faceI] += separatedMass[faceI]*Uf[faceI];
                diameterAccumulator_[faceI] +=
                    separatedMass[faceI]*film().h()[faceI];
            }

            if (massAccumulator_[faceI] > VSMALL)
            {
                accumulationTime_[faceI] += deltaT;
            }
        }

        UfInjectionPtr_().primitiveFieldRef() = Uf.primitiveField();

        separationAggregationModelPtr_->aggregate
        (
            massAccumulator_,
            momentumAccumulator_,
            diameterAccumulator_,
            accumulationTime_,
            massToInject,
            diameterToInject,
            UfInjectionPtr_().primitiveFieldRef()
        );

        if (velocityJitter_ > VSMALL)
        {
            forAll(massToInject, faceI)
            {
                if (massToInject[faceI] > VSMALL)
                {
                    vector& meanVelocity =
                        UfInjectionPtr_().primitiveFieldRef()[faceI];

                    meanVelocity +=
                        velocityJitter_*mag(meanVelocity)*randomDirection();
                }
            }
        }

        availableMass =
            max(availableMass - massAccumulator_ - massToInject, scalar(0));

        addToInjectedMass(sum(massToInject));
    }

    injectionModel::correct();


    if (debug && mesh.time().writeTime())
    {
        {
            areaScalarField areaSeparated
            (
                mesh.newIOobject("separated"),
                mesh,
                dimensionedScalar(dimMass, Zero)
            );
            areaSeparated.primitiveFieldRef() = massRatio;
            areaSeparated.write();
        }

        {
            areaScalarField areaMassToInject
            (
                mesh.newIOobject("massToInject"),
                mesh,
                dimensionedScalar(dimMass, Zero)
            );
            areaMassToInject.primitiveFieldRef() = massToInject;
            areaMassToInject.write();
        }

        if (separationAggregationModelPtr_.valid())
        {
            areaScalarField areaAccumulated
            (
                mesh.newIOobject("massAccumulated"),
                mesh,
                dimensionedScalar(dimMass, Zero)
            );
            areaAccumulated.primitiveFieldRef() = massAccumulator_;
            areaAccumulated.write();
        }
    }
}


// ************************************************************************* //
