/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2025 OpenCFD Ltd.
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
#include "addToRunTimeSelectionTable.H"

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


const Foam::Enum
<
    Foam::regionModels::areaSurfaceFilmModels::filmSeparation::diameterModelType
>
Foam::regionModels::areaSurfaceFilmModels::filmSeparation::
diameterModelTypeNames
({
    { diameterModelType::FILM_THICKNESS, "filmThickness" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar
Foam::regionModels::areaSurfaceFilmModels::filmSeparation::computeDiameter
(
    const label faceI
) const
{
    switch (diameterModel_)
    {
        case diameterModelType::FILM_THICKNESS:
        {
            return film().h()[faceI];
        }
        default:
        {
            return film().h()[faceI];
        }
    }
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
    accumulationEnabled_
    (
        coeffDict_.getOrDefault<bool>("accumulation", false)
    ),
    massAccumulator_(),
    momentumAccumulator_(),
    accumulationTime_(),
    minParcelMass_
    (
        coeffDict_.getOrDefault<scalar>("minParcelMass", 1e-8)
    ),
    maxAccumulationTime_
    (
        coeffDict_.getOrDefault<scalar>("maxAccumulationTime", GREAT)
    ),
    diameterModel_
    (
        diameterModelTypeNames.getOrDefault
        (
            "diameterModel",
            coeffDict_,
            diameterModelType::FILM_THICKNESS
        )
    ),
    velocityJitter_
    (
        coeffDict_.getOrDefault<scalar>("velocityJitter", 0.0)
    ),
    rndGen_(coeffDict_.getOrDefault<label>("seed", 0)),
    UfInjectionPtr_(nullptr)
{
    if (accumulationEnabled_)
    {
        const label nFaces = film.regionMesh().nFaces();

        massAccumulator_.resize(nFaces, Zero);
        momentumAccumulator_.resize(nFaces, Zero);
        accumulationTime_.resize(nFaces, Zero);

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

        Info<< "    filmSeparation: mass accumulation enabled" << nl
            << "        minParcelMass       = " << minParcelMass_ << nl
            << "        maxAccumulationTime = " << maxAccumulationTime_ << nl
            << "        diameterModel       = "
            << diameterModelTypeNames[diameterModel_] << nl
            << "        velocityJitter      = " << velocityJitter_ << endl;
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

    if (!accumulationEnabled_)
    {
        // --- Original (non-accumulated) path ---

        massToInject = massRatio*availableMass;
        diameterToInject = massRatio*film().h();
        availableMass -= massRatio*availableMass;

        addToInjectedMass(sum(massToInject));
    }
    else
    {
        // --- Mass accumulation path ---

        const scalar deltaT = mesh.time().deltaTValue();
        const areaVectorField& Uf = film().Uf();

        // Effective available mass: subtract already-buffered mass to
        // prevent re-separation of the same mass across time steps
        scalarField effectiveAvailMass
        (
            max(availableMass - massAccumulator_, scalar(0))
        );

        // Compute the mass separating this step
        scalarField separatedMass(massRatio*effectiveAvailMass);

        // Accumulate mass, momentum, and time per face
        forAll(separatedMass, faceI)
        {
            if (separatedMass[faceI] > VSMALL)
            {
                massAccumulator_[faceI] += separatedMass[faceI];
                momentumAccumulator_[faceI] +=
                    separatedMass[faceI]*Uf[faceI];
            }

            // Advance time for any face with accumulated mass
            if (massAccumulator_[faceI] > VSMALL)
            {
                accumulationTime_[faceI] += deltaT;
            }
        }

        // Initialise UfInjection to the instantaneous film velocity
        UfInjectionPtr_().primitiveFieldRef() = Uf.primitiveField();

        scalar totalInjectedMass = 0;

        // Check thresholds and inject
        forAll(massAccumulator_, faceI)
        {
            const bool thresholdMet =
            (
                massAccumulator_[faceI] >= minParcelMass_
             || accumulationTime_[faceI] >= maxAccumulationTime_
            );

            if (thresholdMet && massAccumulator_[faceI] > VSMALL)
            {
                // --- Inject accumulated mass ---

                massToInject[faceI] = massAccumulator_[faceI];
                diameterToInject[faceI] = computeDiameter(faceI);
                totalInjectedMass += massAccumulator_[faceI];

                // Mass-weighted average velocity
                vector meanVelocity =
                    momentumAccumulator_[faceI]/massAccumulator_[faceI];

                // Optional stochastic jitter
                if (velocityJitter_ > VSMALL)
                {
                    vector randomDir
                    (
                        rndGen_.position<scalar>(-1, 1),
                        rndGen_.position<scalar>(-1, 1),
                        rndGen_.position<scalar>(-1, 1)
                    );

                    randomDir.normalise();

                    meanVelocity +=
                        velocityJitter_*mag(meanVelocity)*randomDir;
                }

                UfInjectionPtr_().primitiveFieldRef()[faceI] = meanVelocity;

                // Reset buffers
                massAccumulator_[faceI] = 0;
                momentumAccumulator_[faceI] = Zero;
                accumulationTime_[faceI] = 0;
            }
            else
            {
                // Suppress injection for this face
                massToInject[faceI] = 0;
                diameterToInject[faceI] = 0;
            }
        }

        // Reduce available mass for subsequent injection models
        availableMass -= separatedMass;

        addToInjectedMass(totalInjectedMass);
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

        if (accumulationEnabled_)
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
