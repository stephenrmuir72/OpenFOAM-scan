/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "temporalSpreading.H"
#include "runTimeSelectionTables.H"
#include "PstreamReduceOps.H"
#include "temporalKernel.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

Foam::word eventFieldName
(
    const Foam::word& channel,
    const Foam::label id
)
{
    return "temporalSpreading" + channel + '_' + Foam::name(id);
}

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(temporalSpreading, 0);
defineRunTimeSelectionTable(temporalSpreading, dictionary);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

IOobject temporalSpreading::stateIO(const faMesh& mesh)
{
    return IOobject
    (
        "temporalSpreadingState",
        mesh.time().timeName(),
        mesh.thisDb(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE,
        IOobject::REGISTER
    );
}


IOobject temporalSpreading::eventIO
(
    const word& fieldName,
    const faMesh& mesh,
    const label id,
    const bool readExisting
)
{
    return IOobject
    (
        eventFieldName(fieldName, id),
        mesh.time().timeName(),
        mesh.thisDb(),
        readExisting ? IOobject::MUST_READ : IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        IOobject::REGISTER
    );
}


void temporalSpreading::restoreState()
{
    if (!state_.found("eventIds"))
    {
        return;
    }

    const word savedModelType(state_.get<word>("modelType"));
    const word savedKernel(state_.get<word>("temporalKernel"));

    if (savedModelType != modelType_)
    {
        FatalErrorInFunction
            << "Cannot restart temporalSpreading state with model type '"
            << savedModelType << "' using configured type '" << modelType_
            << "'. Remove the pending state only after it has fully released."
            << abort(FatalError);
    }

    if (savedKernel != kernelPtr_->type())
    {
        FatalErrorInFunction
            << "Cannot restart temporalSpreading state with kernel '"
            << savedKernel << "' using configured kernel '"
            << kernelPtr_->type()
            << "'. Remove the pending state only after it has fully released."
            << abort(FatalError);
    }

    const labelList eventIds(state_.get<labelList>("eventIds"));
    const scalarList eventStartTimes(state_.get<scalarList>("eventStartTimes"));

    if (eventIds.size() != eventStartTimes.size())
    {
        FatalErrorInFunction
            << "Inconsistent temporalSpreading restart metadata: "
            << eventIds.size() << " event ids but "
            << eventStartTimes.size() << " start times."
            << abort(FatalError);
    }

    nextEventId_ = state_.get<label>("nextEventId");

    label maxEventId = -1;

    forAll(eventIds, eventi)
    {
        const label id = eventIds[eventi];

        if (id < 0 || !std::isfinite(eventStartTimes[eventi]))
        {
            FatalErrorInFunction
                << "Invalid temporalSpreading restart event id or start time."
                << abort(FatalError);
        }

        for (label previousi = 0; previousi < eventi; ++previousi)
        {
            if (eventIds[previousi] == id)
            {
                FatalErrorInFunction
                    << "Duplicate temporalSpreading restart event id " << id
                    << abort(FatalError);
            }
        }

        events_.emplace_back(mesh_, id, eventStartTimes[eventi], true);
        maxEventId = Foam::max(maxEventId, id);
    }

    if (nextEventId_ <= maxEventId)
    {
        FatalErrorInFunction
            << "Invalid next temporalSpreading event id " << nextEventId_
            << "; it must exceed the largest active event id " << maxEventId
            << abort(FatalError);
    }
}


void temporalSpreading::updateState()
{
    label nEvents = 0;

    forAll(events_, eventi)
    {
        if (events_.set(eventi))
        {
            ++nEvents;
        }
    }

    labelList eventIds(nEvents);
    scalarList eventStartTimes(nEvents);
    label eventi = 0;

    forAll(events_, eventSlot)
    {
        if (events_.set(eventSlot))
        {
            eventIds[eventi] = events_[eventSlot].id();
            eventStartTimes[eventi] = events_[eventSlot].startTime();
            ++eventi;
        }
    }

    state_.set("modelType", modelType_);
    state_.set("temporalKernel", kernelPtr_->type());
    state_.set("nextEventId", nextEventId_);
    state_.set("eventIds", eventIds);
    state_.set("eventStartTimes", eventStartTimes);
}


void temporalSpreading::checkResolution
(
    const areaScalarField& duration,
    const areaScalarField& mass,
    const areaVectorField& momentum,
    const areaScalarField& normalImpulse
) const
{
    if (!(minStepsPerDuration_ > 0))
    {
        return;
    }

    bool localActive = false;
    scalar localMinDuration = GREAT;

    forAll(duration, facei)
    {
        if
        (
            mag(mass[facei]) > VSMALL
         || mag(momentum[facei]) > VSMALL
         || mag(normalImpulse[facei]) > VSMALL
        )
        {
            localActive = true;
            localMinDuration = Foam::min(localMinDuration, duration[facei]);
        }
    }

    if (!returnReduceOr(localActive))
    {
        return;
    }

    const scalar minDuration = returnReduce(localMinDuration, minOp<scalar>());
    const scalar deltaT = mesh_.time().deltaTValue();
    const scalar steps = minDuration/deltaT;

    if (steps < minStepsPerDuration_)
    {
        WarningInFunction
            << "temporalSpreading duration is resolved by " << steps
            << " CFD steps, below impactTimeResolution.minStepsPerDuration "
            << minStepsPerDuration_ << ". The physical duration is unchanged."
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

temporalSpreading::event::event
(
    const faMesh& mesh,
    const label id,
    const scalar startTime,
    const bool readExisting
)
:
    id_(id),
    startTime_(startTime),
    massPtr_
    (
        readExisting
      ? new areaScalarField
        (
            temporalSpreading::eventIO("Mass", mesh, id, true),
            mesh
        )
      : new areaScalarField
        (
            temporalSpreading::eventIO("Mass", mesh, id, false),
            mesh,
            dimensionedScalar(dimMass, Zero)
        )
    ),
    momentumPtr_
    (
        readExisting
      ? new areaVectorField
        (
            temporalSpreading::eventIO("Momentum", mesh, id, true),
            mesh
        )
      : new areaVectorField
        (
            temporalSpreading::eventIO("Momentum", mesh, id, false),
            mesh,
            dimensionedVector(dimMass*dimVelocity, Zero)
        )
    ),
    normalImpulsePtr_
    (
        readExisting
      ? new areaScalarField
        (
            temporalSpreading::eventIO
            (
                "NormalImpulse",
                mesh,
                id,
                true
            ),
            mesh
        )
      : new areaScalarField
        (
            temporalSpreading::eventIO
            (
                "NormalImpulse",
                mesh,
                id,
                false
            ),
            mesh,
            dimensionedScalar(dimMass*dimVelocity, Zero)
        )
    ),
    durationPtr_
    (
        readExisting
      ? new areaScalarField
        (
            temporalSpreading::eventIO("Duration", mesh, id, true),
            mesh
        )
      : new areaScalarField
        (
            temporalSpreading::eventIO("Duration", mesh, id, false),
            mesh,
            dimensionedScalar(dimTime, Zero)
        )
    )
{
    if
    (
        readExisting
     &&
        (
            mass().dimensions() != dimMass
         || momentum().dimensions() != dimMass*dimVelocity
         || normalImpulse().dimensions() != dimMass*dimVelocity
         || duration().dimensions() != dimTime
        )
    )
    {
        FatalErrorInFunction
            << "Invalid dimensions in temporal-spreading restart event "
            << id
            << abort(FatalError);
    }
}


temporalSpreading::temporalSpreading
(
    const word& modelType,
    const faMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    modelType_(modelType),
    kernelPtr_(temporalKernel::New(dict)),
    minStepsPerDuration_(-1),
    state_(stateIO(mesh), word::null),
    nextEventId_(0),
    events_()
{
    if (dict.readIfPresent("minStepsPerDuration", minStepsPerDuration_))
    {
        if (!(minStepsPerDuration_ > 0) || !std::isfinite(minStepsPerDuration_))
        {
            FatalIOErrorInFunction(dict)
                << "minStepsPerDuration must be finite and > 0, received "
                << minStepsPerDuration_
                << exit(FatalIOError);
        }
    }

    restoreState();
    updateState();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

temporalSpreading::~temporalSpreading()
{}  // temporalKernel was forward declared


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void temporalSpreading::duration
(
    areaScalarField& duration,
    const areaScalarField& mass,
    const areaVectorField& momentum,
    const areaScalarField& normalImpulse,
    const areaScalarField& diameter,
    const areaScalarField& rho,
    const areaScalarField& sigma
) const
{
    calculateDuration
    (
        duration,
        mass,
        momentum,
        normalImpulse,
        diameter,
        rho,
        sigma
    );

    forAll(duration, facei)
    {
        if
        (
            mag(mass[facei]) > VSMALL
         || mag(momentum[facei]) > VSMALL
         || mag(normalImpulse[facei]) > VSMALL
        )
        {
            if (!(duration[facei] > VSMALL) || !std::isfinite(duration[facei]))
            {
                FatalErrorInFunction
                    << "Invalid temporal-spreading duration " << duration[facei]
                    << " on active face " << facei
                    << abort(FatalError);
            }
        }
    }

    checkResolution(duration, mass, momentum, normalImpulse);
}


void temporalSpreading::enqueue
(
    const areaScalarField& mass,
    const areaVectorField& momentum,
    const areaScalarField& normalImpulse,
    const areaScalarField& duration,
    const scalar startTime
)
{
    bool localActive = false;

    forAll(mass, facei)
    {
        localActive =
            localActive
         || mag(mass[facei]) > VSMALL
         || mag(momentum[facei]) > VSMALL
         || mag(normalImpulse[facei]) > VSMALL;
    }

    if (!returnReduceOr(localActive))
    {
        return;
    }

    if (!std::isfinite(startTime))
    {
        FatalErrorInFunction
            << "Non-finite temporal-spreading event start time " << startTime
            << abort(FatalError);
    }

    event& pending =
        events_.emplace_back(mesh_, nextEventId_++, startTime, false);

    pending.mass().primitiveFieldRef() = mass.primitiveField();
    pending.momentum().primitiveFieldRef() = momentum.primitiveField();
    pending.normalImpulse().primitiveFieldRef() =
        normalImpulse.primitiveField();
    pending.duration().primitiveFieldRef() = duration.primitiveField();

    pending.mass().correctBoundaryConditions();
    pending.momentum().correctBoundaryConditions();
    pending.normalImpulse().correctBoundaryConditions();
    pending.duration().correctBoundaryConditions();

    updateState();
}


void temporalSpreading::release
(
    areaScalarField& mass,
    areaVectorField& momentum,
    areaScalarField& normalImpulse,
    const scalar beginTime,
    const scalar endTime
)
{
    if
    (
        !std::isfinite(beginTime)
     || !std::isfinite(endTime)
     || endTime < beginTime
    )
    {
        FatalErrorInFunction
            << "Invalid temporal-spreading release interval ["
            << beginTime << ", " << endTime << ']'
            << abort(FatalError);
    }

    mass.primitiveFieldRef() = Zero;
    momentum.primitiveFieldRef() = Zero;
    normalImpulse.primitiveFieldRef() = Zero;

    forAll(events_, eventi)
    {
        if (!events_.set(eventi))
        {
            continue;
        }

        event& pending = events_[eventi];
        areaScalarField& pendingMass = pending.mass();
        areaVectorField& pendingMomentum = pending.momentum();
        areaScalarField& pendingNormalImpulse = pending.normalImpulse();
        const areaScalarField& pendingDuration = pending.duration();

        bool localPending = false;

        forAll(pendingMass, facei)
        {
            const bool active =
                mag(pendingMass[facei]) > 0
             || mag(pendingMomentum[facei]) > 0
             || mag(pendingNormalImpulse[facei]) > 0;

            if (!active)
            {
                continue;
            }

            const scalar eventDuration = pendingDuration[facei];

            if
            (
                !(eventDuration > VSMALL)
             || !std::isfinite(eventDuration)
            )
            {
                FatalErrorInFunction
                    << "Invalid pending temporal-spreading duration "
                    << eventDuration << " on event " << pending.id()
                    << ", face " << facei
                    << abort(FatalError);
            }

            const scalar age0 = beginTime - pending.startTime();
            const scalar age1 = endTime - pending.startTime();
            const bool complete = age1 >= eventDuration;

            scalar fraction = 0;

            if (complete)
            {
                fraction = 1;
            }
            else
            {
                const scalar F0 = kernelPtr_->cumulative(age0, eventDuration);
                const scalar F1 = kernelPtr_->cumulative(age1, eventDuration);

                if (F1 < F0 || F0 < 0 || F1 > 1)
                {
                    FatalErrorInFunction
                        << "Invalid temporal-kernel cumulative interval for "
                        << "event " << pending.id() << ", face " << facei
                        << ": F0=" << F0 << ", F1=" << F1
                        << abort(FatalError);
                }

                if (F1 >= 1)
                {
                    fraction = 1;
                }
                else if (F0 < 1)
                {
                    fraction = (F1 - F0)/(1 - F0);
                    fraction =
                        Foam::max(scalar(0), Foam::min(scalar(1), fraction));
                }
            }

            const scalar releasedMass = fraction*pendingMass[facei];
            const vector releasedMomentum = fraction*pendingMomentum[facei];
            const scalar releasedNormalImpulse =
                fraction*pendingNormalImpulse[facei];

            mass[facei] += releasedMass;
            momentum[facei] += releasedMomentum;
            normalImpulse[facei] += releasedNormalImpulse;

            if (fraction == 1)
            {
                pendingMass[facei] = Zero;
                pendingMomentum[facei] = Zero;
                pendingNormalImpulse[facei] = Zero;
            }
            else
            {
                pendingMass[facei] -= releasedMass;
                pendingMomentum[facei] -= releasedMomentum;
                pendingNormalImpulse[facei] -= releasedNormalImpulse;

                localPending =
                    localPending
                 || mag(pendingMass[facei]) > 0
                 || mag(pendingMomentum[facei]) > 0
                 || mag(pendingNormalImpulse[facei]) > 0;
            }
        }

        if (!returnReduceOr(localPending))
        {
            events_.set(eventi, nullptr);
        }
    }

    mass.correctBoundaryConditions();
    momentum.correctBoundaryConditions();
    normalImpulse.correctBoundaryConditions();

    updateState();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
