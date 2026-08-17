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
#include "diameterModel.H"
#include "filmSeparationModel.H"
#include "morphologyModel.H"
#include "PstreamReduceOps.H"
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


Foam::word defaultDiameterModel
(
    const Foam::regionModels::areaSurfaceFilmModels::morphologyModel::
        morphologyType morphology
)
{
    using morphologyModel =
        Foam::regionModels::areaSurfaceFilmModels::morphologyModel;

    switch (morphology)
    {
        case morphologyModel::DIRECT_DROPLET:
        case morphologyModel::DRIP:
            return "filmThickness";

        case morphologyModel::LIGAMENT:
            return "RayleighPlateau";

        case morphologyModel::SHEET:
            return "KelvinHelmholtz";

        default:
            return "filmThickness";
    }
}


Foam::vector sumMomentum(const Foam::vectorField& momentum)
{
    Foam::vector result(Foam::Zero);

    forAll(momentum, facei)
    {
        result += momentum[facei];
    }

    return result;
}


Foam::vector summedMassVelocity
(
    const Foam::scalarField& mass,
    const Foam::vectorField& velocity
)
{
    Foam::vector result(Foam::Zero);

    forAll(mass, facei)
    {
        result += mass[facei]*velocity[facei];
    }

    return result;
}


Foam::word morphologyMassPropertyName(const Foam::label morphology)
{
    using morphologyModel =
        Foam::regionModels::areaSurfaceFilmModels::morphologyModel;

    const auto type =
        static_cast<morphologyModel::morphologyType>(morphology);

    return Foam::word
    (
        "morphologySeparatedMass_"
      + morphologyModel::morphologyTypeNames[type]
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
    morphologyModelPtr_(morphologyModel::NewIfPresent(film, coeffDict_)),
    diameterModels_(morphologyModel::nMorphologyTypes),
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
    morphologyMassTotal_
    (
        accumulatorIO("morphologyMassTotal", film, true),
        scalarField(morphologyModel::nMorphologyTypes, Zero)
    ),
    morphologySeparatedMass_(morphologyModel::nMorphologyTypes, Zero),
    massConservationResidual_(Zero),
    momentumConservationResidual_(Zero),
    velocityJitter_
    (
        coeffDict_.getOrDefault<scalar>("velocityJitter", 0.0)
    ),
    rndGen_(coeffDict_.getOrDefault<label>("seed", 0) + Pstream::myProcNo()),
    UfInjectionPtr_(nullptr)
{
    boolList possibleMorphologies(morphologyModel::nMorphologyTypes, false);

    if (morphologyModelPtr_.valid())
    {
        possibleMorphologies = morphologyModelPtr_->possibleTypes();
    }
    else
    {
        possibleMorphologies.set(morphologyModel::DIRECT_DROPLET);
    }

    for
    (
        label morphology = 0;
        morphology < possibleMorphologies.size();
        ++morphology
    )
    {
        if (!possibleMorphologies[morphology])
        {
            continue;
        }

        dictionary diameterDict(coeffDict_);
        word modelType
        (
            morphologyModelPtr_.valid()
          ? defaultDiameterModel
            (
                static_cast<morphologyModel::morphologyType>(morphology)
            )
          : coeffDict_.getOrDefault<word>
            (
                "diameterModel",
                "filmThickness"
            )
        );

        if (morphologyModelPtr_.valid())
        {
            if (const dictionary* diameterDictPtr =
                coeffDict_.findDict("diameterModels"))
            {
                const word& morphologyName =
                    morphologyModel::morphologyTypeNames
                    [
                        static_cast<morphologyModel::morphologyType>
                        (
                            morphology
                        )
                    ];

                if
                (
                    const dictionary* modelDictPtr =
                        diameterDictPtr->findDict(morphologyName)
                )
                {
                    diameterDict.merge(*modelDictPtr);
                    modelType = modelDictPtr->getOrDefault<word>
                    (
                        "type",
                        modelDictPtr->getOrDefault<word>
                        (
                            "diameterModel",
                            modelType
                        )
                    );
                }
            }
        }

        diameterModels_.set
        (
            morphology,
            diameterModel::New(film, diameterDict, modelType)
        );
    }

    if
    (
        morphologyModelPtr_.valid()
     && separationAggregationModelPtr_.valid()
     && separationAggregationModelPtr_->overridesDiameter()
    )
    {
        FatalIOErrorInFunction(coeffDict_)
            << "The fixedDiameter option of separationAggregation is a "
            << "legacy aggregate-level override and cannot be used with an "
            << "explicit morphology model. Select a physical model in "
            << "diameterModels instead."
            << exit(FatalIOError);
    }

    if (separationAggregationModelPtr_.valid())
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

    if (morphologyModelPtr_.valid())
    {
        Info<< "    filmSeparation: physical morphology model "
            << morphologyModelPtr_->type() << endl;
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
    const areaVectorField& Uf = film().Uf();

    // Calculate the mass ratio of film separation
    tmp<scalarField> tmassRatio = filmSeparationModelPtr_->separatedMassRatio();
    const auto& massRatio = tmassRatio.cref();

    tmp<labelField> tMorphology
    (
        morphologyModelPtr_.valid()
      ? morphologyModelPtr_->morphology()
      : tmp<labelField>::New
        (
            mesh.nFaces(),
            label(morphologyModel::DIRECT_DROPLET)
        )
    );
    const labelField& morphology = tMorphology.cref();

    scalarField physicalDiameter(mesh.nFaces(), Zero);

    for
    (
        label morphologyType = 0;
        morphologyType < morphologyModel::nMorphologyTypes;
        ++morphologyType
    )
    {
        const diameterModel* diameterModelPtr =
            diameterModels_.set(morphologyType);

        if (!diameterModelPtr)
        {
            continue;
        }

        tmp<scalarField> tDiameter = diameterModelPtr->diameter();
        const scalarField& diameter = tDiameter.cref();

        if (diameter.size() != mesh.nFaces())
        {
            FatalErrorInFunction
                << "Physical diameter model "
                << diameterModelPtr->type()
                << " returned " << diameter.size() << " entries for "
                << mesh.nFaces() << " finite-area faces."
                << abort(FatalError);
        }

        forAll(morphology, faceI)
        {
            if (morphology[faceI] == morphologyType)
            {
                physicalDiameter[faceI] = diameter[faceI];
            }
        }
    }

    forAll(morphology, faceI)
    {
        if
        (
            morphology[faceI] < 0
         || morphology[faceI] >= morphologyModel::nMorphologyTypes
         || !diameterModels_.set(morphology[faceI])
        )
        {
            FatalErrorInFunction
                << "Invalid or unconfigured morphology ID "
                << morphology[faceI] << " on finite-area face " << faceI
                << abort(FatalError);
        }
    }

    if (!separationAggregationModelPtr_.valid())
    {
        massToInject = massRatio*availableMass;
        diameterToInject = Zero;

        const bool legacyDirectDiameter =
        (
            !morphologyModelPtr_.valid()
         && diameterModels_[morphologyModel::DIRECT_DROPLET].
                scalesWithSeparationRatio()
        );

        scalar separatedMassTotal = 0;
        vector separatedMomentum(Zero);

        forAll(massToInject, faceI)
        {
            if (massToInject[faceI] > VSMALL)
            {
                const label morphologyType = morphology[faceI];

                diameterToInject[faceI] =
                (
                    legacyDirectDiameter
                  ? massRatio[faceI]*physicalDiameter[faceI]
                  : physicalDiameter[faceI]
                );

                morphologySeparatedMass_[morphologyType] +=
                    massToInject[faceI];
                separatedMassTotal += massToInject[faceI];
                separatedMomentum += massToInject[faceI]*Uf[faceI];
            }
        }

        availableMass -= massToInject;

        massConservationResidual_ +=
            separatedMassTotal - sum(massToInject);
        momentumConservationResidual_ +=
            separatedMomentum - summedMassVelocity(massToInject, Uf);

        addToInjectedMass(sum(massToInject));
    }
    else
    {
        const scalar deltaT = mesh.time().deltaTValue();
        const scalar bufferedMassBefore = sum(massAccumulator_);
        const vector bufferedMomentumBefore = sumMomentum(momentumAccumulator_);

        scalarField effectiveAvailMass
        (
            max(availableMass - massAccumulator_, scalar(0))
        );

        scalarField separatedMass(massRatio*effectiveAvailMass);
        scalar separatedMassTotal = 0;
        vector separatedMomentum(Zero);

        forAll(separatedMass, faceI)
        {
            if (separatedMass[faceI] > VSMALL)
            {
                massAccumulator_[faceI] += separatedMass[faceI];
                momentumAccumulator_[faceI] += separatedMass[faceI]*Uf[faceI];
                diameterAccumulator_[faceI] +=
                    separatedMass[faceI]*physicalDiameter[faceI];
                morphologySeparatedMass_[morphology[faceI]] +=
                    separatedMass[faceI];
                separatedMassTotal += separatedMass[faceI];
                separatedMomentum += separatedMass[faceI]*Uf[faceI];
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

        massConservationResidual_ +=
        (
            bufferedMassBefore + separatedMassTotal
          - sum(massAccumulator_) - sum(massToInject)
        );
        momentumConservationResidual_ +=
        (
            bufferedMomentumBefore + separatedMomentum
          - sumMomentum(momentumAccumulator_)
          - summedMassVelocity
            (
                massToInject,
                UfInjectionPtr_().primitiveField()
            )
        );

        addToInjectedMass(sum(massToInject));
    }

    if (debug)
    {
        const vectorField& injectionVelocity =
        (
            UfInjectionPtr_
          ? UfInjectionPtr_().primitiveField()
          : film().Uf().primitiveField()
        );

        label localFaces = 0;
        scalar localMass = 0;
        vector localMomentum(Zero);

        forAll(massToInject, facei)
        {
            if (massToInject[facei] > VSMALL)
            {
                ++localFaces;
                localMass += massToInject[facei];
                localMomentum +=
                    massToInject[facei]*injectionVelocity[facei];

                if (debug >= 2)
                {
                    Pout<< "[filmSeparation:injection] t="
                        << mesh.time().timeOutputValue()
                        << " face=" << facei
                        << " mass=" << massToInject[facei]
                        << " diameter=" << diameterToInject[facei]
                        << " velocity=" << injectionVelocity[facei]
                        << " momentum="
                        << massToInject[facei]*injectionVelocity[facei]
                        << endl;
                }
            }
        }

        Info<< "[filmSeparation:injection] t="
            << mesh.time().timeOutputValue()
            << " faces=" << returnReduce(localFaces, sumOp<label>())
            << " mass="
            << returnReduce(localMass, sumOp<scalar>())
            << " momentum="
            << returnReduce(localMomentum, sumOp<vector>())
            << " bufferedMass="
            << (
                    separationAggregationModelPtr_
                  ? returnReduce(sum(massAccumulator_), sumOp<scalar>())
                  : scalar(0)
               )
            << endl;
    }

    injectionModel::correct();

    if (writeTime())
    {
        scalarField morphologyMass
        (
            morphologyModel::nMorphologyTypes,
            Zero
        );

        scalar totalMorphologyMass = 0;

        forAll(morphologyMass, morphologyType)
        {
            const scalar stepMass =
                returnReduce
                (
                    morphologySeparatedMass_[morphologyType],
                    sumOp<scalar>()
                );

            morphologyMassTotal_[morphologyType] += stepMass;
            morphologyMass[morphologyType] =
                morphologyMassTotal_[morphologyType];
            totalMorphologyMass += morphologyMass[morphologyType];

            setModelProperty<scalar>
            (
                morphologyMassPropertyName(morphologyType),
                morphologyMass[morphologyType]
            );
            morphologySeparatedMass_[morphologyType] = 0;
        }

        const scalar massResidual =
            returnReduce(massConservationResidual_, sumOp<scalar>());
        const vector momentumResidual =
            returnReduce(momentumConservationResidual_, sumOp<vector>());

        setModelProperty<scalar>("massTransferClosure", massResidual);
        setModelProperty<vector>("momentumTransferClosure", momentumResidual);

        Info<< "    filmSeparation morphology-resolved separated mass:" << nl;

        forAll(morphologyMass, morphologyType)
        {
            const auto type =
                static_cast<morphologyModel::morphologyType>(morphologyType);

            Info<< "        "
                << morphologyModel::morphologyTypeNames[type] << " = "
                << morphologyMass[morphologyType] << " kg"
                << " (fraction "
                << morphologyMass[morphologyType]
                   /max(totalMorphologyMass, VSMALL) << ')' << nl;
        }

        Info<< "        transfer mass closure = " << massResidual << " kg" << nl
            << "        transfer momentum closure = " << momentumResidual
            << " kg.m/s" << endl;

        massConservationResidual_ = 0;
        momentumConservationResidual_ = Zero;
    }


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

        {
            areaScalarField areaDiameterToInject
            (
                mesh.newIOobject("diameterToInject"),
                mesh,
                dimensionedScalar(dimLength, Zero)
            );
            areaDiameterToInject.primitiveFieldRef() = diameterToInject;
            areaDiameterToInject.write();
        }

        {
            areaVectorField areaVelocityToInject
            (
                mesh.newIOobject("velocityToInject"),
                mesh,
                dimensionedVector(dimVelocity, Zero)
            );

            areaVelocityToInject.primitiveFieldRef() =
            (
                separationAggregationModelPtr_.valid()
              ? UfInjectionPtr_().primitiveField()
              : Uf.primitiveField()
            );
            areaVelocityToInject.write();
        }

        {
            areaScalarField areaMorphology
            (
                mesh.newIOobject("morphology"),
                mesh,
                dimensionedScalar(dimless, Zero)
            );

            forAll(morphology, faceI)
            {
                areaMorphology.primitiveFieldRef()[faceI] = morphology[faceI];
            }

            areaMorphology.write();
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
