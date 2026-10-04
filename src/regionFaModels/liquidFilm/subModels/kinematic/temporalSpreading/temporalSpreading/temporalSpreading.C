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
#include "Function1.H"
#include "OStringStream.H"
#include "TableBase.H"
#include "SHA1Digest.H"

#include <cmath>

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

Foam::string temporalKernelSignature
(
    const Foam::word& kernelType,
    const Foam::dictionary& dict
)
{
    if (kernelType != "custom")
    {
        return kernelType;
    }

    Foam::autoPtr<Foam::Function1<Foam::scalar>> releaseRate
    (
        Foam::Function1<Foam::scalar>::New("releaseRate", dict)
    );
    Foam::OStringStream os;
    releaseRate->writeData(os);

    if (releaseRate->type() == "tableFile")
    {
        const auto* table =
            dynamic_cast
            <
                const Foam::Function1Types::TableBase<Foam::scalar>*
            >
            (&releaseRate());

        if (!table)
        {
            FatalErrorInFunction
                << "Function1 tableFile does not provide table data"
                << abort(Foam::FatalError);
        }

        // TableFile writes its path, not the table data loaded from it.
        const auto x = table->x();
        const auto y = table->y();
        Foam::OStringStream values;
        values << x() << y();
        const Foam::string tableValues(values.str());
        const Foam::SHA1Digest tableDigest
        (
            tableValues.c_str(),
            tableValues.size()
        );
        os << Foam::nl << "loadedTableDigest " << tableDigest;
    }

    const Foam::string config(os.str());
    const Foam::SHA1Digest configDigest(config.c_str(), config.size());
    return Foam::string(configDigest.str().c_str());
}

Foam::scalar releaseFraction
(
    const Foam::regionModels::areaSurfaceFilmModels::temporalKernel& kernel,
    const Foam::scalar age0,
    const Foam::scalar age1,
    const Foam::scalar duration
)
{
    if (age1 >= duration)
    {
        return 1;
    }

    const Foam::scalar F0 = kernel.cumulative(age0, duration);
    const Foam::scalar F1 = kernel.cumulative(age1, duration);
    if
    (
        !std::isfinite(F0)
     || !std::isfinite(F1)
     || F0 < 0
     || F1 < F0
     || F1 > 1
    )
    {
        FatalErrorInFunction
            << "Invalid temporal-kernel cumulative interval: F0="
            << F0 << ", F1=" << F1 << abort(Foam::FatalError);
    }
    if (F1 >= 1)
    {
        return 1;
    }
    return F0 < 1 ? Foam::max(Foam::scalar(0), Foam::min(
        Foam::scalar(1), (F1 - F0)/(1 - F0)
    )) : 0;
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
        if (state_.found("stateVersion") || state_.found("sparseFaces"))
        {
            FatalErrorInFunction
                << "Incomplete temporalSpreadingState: missing eventIds"
                << abort(FatalError);
        }
        return;
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

    const labelList sparseFaces
    (
        state_.found("sparseFaces")
      ? state_.get<labelList>("sparseFaces") : labelList()
    );
    const bool hasPendingState =
        !eventIds.empty() || !sparseFaces.empty();

    const word savedModelType(state_.get<word>("modelType"));
    const word savedKernel(state_.get<word>("temporalKernel"));

    if (hasPendingState && savedModelType != modelType_)
    {
        FatalErrorInFunction
            << "Cannot restart temporalSpreading state with model type '"
            << savedModelType << "' using configured type '" << modelType_
            << "'. Remove the pending state only after it has fully released."
            << abort(FatalError);
    }

    if (hasPendingState && savedKernel != kernelPtr_->type())
    {
        FatalErrorInFunction
            << "Cannot restart temporalSpreading state with kernel '"
            << savedKernel << "' using configured kernel '"
            << kernelPtr_->type()
            << "'. Remove the pending state only after it has fully released."
            << abort(FatalError);
    }

    if (savedKernel == "custom" && hasPendingState)
    {
        if (!state_.found("temporalKernelSignature"))
        {
            FatalErrorInFunction
                << "Cannot resume pending custom temporal-spreading events "
                << "because the saved state has no releaseRate signature. "
                << "Resume with the original executable and configuration "
                << "until all events release, then write a "
                << "new restart."
                << abort(FatalError);
        }

        const string savedKernelSignature
        (
            state_.get<string>("temporalKernelSignature")
        );
        if (savedKernelSignature != temporalKernelSignature_)
        {
            FatalErrorInFunction
                << "Cannot resume pending custom temporal-spreading events "
                << "with a changed releaseRate signature. Keep the "
                << "original configuration until all events release, then "
                << "write a new restart."
                << abort(FatalError);
        }
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

    if (!state_.found("stateVersion"))
    {
        return;
    }

    const label stateVersion = state_.get<label>("stateVersion");
    if (!temporalSpreadingValidation::supportsSparseStateVersion(stateVersion))
    {
        FatalErrorInFunction
            << "Unsupported temporalSpreadingState version"
            << abort(FatalError);
    }

    const wordList required
    ({
        "sparseFaces", "sparseStartTimes", "sparseDurations",
        "sparseMasses", "sparseMomenta", "sparseNormalImpulses",
        "sparseMassDimensions", "sparseMomentumDimensions",
        "sparseNormalImpulseDimensions", "sparseDurationDimensions"
    });
    for (const word& entry : required)
    {
        if (!state_.found(entry))
        {
            FatalErrorInFunction
                << "Incomplete temporalSpreadingState: missing " << entry
                << abort(FatalError);
        }
    }
    if
    (
        state_.get<dimensionSet>("sparseMassDimensions") != dimMass
     || state_.get<dimensionSet>("sparseMomentumDimensions")
        != dimMass*dimVelocity
     || state_.get<dimensionSet>("sparseNormalImpulseDimensions")
        != dimMass*dimVelocity
     || state_.get<dimensionSet>("sparseDurationDimensions") != dimTime
    )
    {
        FatalErrorInFunction
            << "Invalid sparse temporal-spreading state dimensions"
            << abort(FatalError);
    }

    const labelList faces(state_.get<labelList>("sparseFaces"));
    const scalarList starts(state_.get<scalarList>("sparseStartTimes"));
    const scalarList durations(state_.get<scalarList>("sparseDurations"));
    const scalarList masses(state_.get<scalarList>("sparseMasses"));
    const vectorField momenta(state_.get<vectorField>("sparseMomenta"));
    const scalarList impulses
    (
        state_.get<scalarList>("sparseNormalImpulses")
    );

    if
    (
        faces.size() != starts.size()
     || faces.size() != durations.size()
     || faces.size() != masses.size()
     || faces.size() != momenta.size()
     || faces.size() != impulses.size()
    )
    {
        FatalErrorInFunction
            << "Inconsistent sparse temporal-spreading state lengths"
            << abort(FatalError);
    }

    if (!temporalSpreadingValidation::canRestoreSparseState
        (stateVersion, faces.size()))
    {
        FatalErrorInFunction
            << "Cannot restore temporal-spreading state version 2 with "
            << faces.size() << " pending sparse contribution(s): the "
            << "state does not contain finite-area face-addressing "
            << "metadata. Resume with the original executable and "
            << "unchanged mesh until these contributions are released, "
            << "then write a new restart."
            << abort(FatalError);
    }

    if (stateVersion == 2)
    {
        sparseFaceLabels_ = mesh_.faceLabels();
    }
    else
    {
        if (!state_.found("sparseFaceLabels"))
        {
            FatalErrorInFunction
                << "Incomplete temporalSpreadingState: missing "
                << "sparseFaceLabels"
                << abort(FatalError);
        }

        sparseFaceLabels_ = state_.get<labelList>("sparseFaceLabels");
        if
        (
            !faces.empty()
         && !temporalSpreadingValidation::faceAddressingMatches
            (sparseFaceLabels_, mesh_.faceLabels())
        )
        {
            FatalErrorInFunction
                << "Cannot restore " << faces.size()
                << " pending sparse temporal-spreading contribution(s) "
                << "because finite-area face addressing changed."
                << abort(FatalError);
        }
    }

    forAll(faces, i)
    {
        const contribution transfer
        {
            faces[i], masses[i], momenta[i], impulses[i], durations[i]
        };
        if
        (
            transfer.face < 0
         || transfer.face >= mesh_.nFaces()
         || !(transfer.duration > VSMALL)
         || !std::isfinite(transfer.duration)
         || !std::isfinite(starts[i])
         || !std::isfinite(transfer.mass)
         || !std::isfinite(mag(transfer.momentum))
         || !std::isfinite(transfer.normalImpulse)
         || (
                mag(transfer.mass) == 0
             && mag(transfer.momentum) == 0
             && mag(transfer.normalImpulse) == 0
            )
        )
        {
            FatalErrorInFunction
                << "Invalid sparse temporal-spreading record " << i
                << abort(FatalError);
        }
        sparse_.append(pendingContribution{transfer, starts[i]});
    }
}


void temporalSpreading::checkSparseFaceAddressing() const
{
    if
    (
        !sparse_.empty()
     && !temporalSpreadingValidation::faceAddressingMatches
        (sparseFaceLabels_, mesh_.faceLabels())
    )
    {
        FatalErrorInFunction
            << "Cannot use " << sparse_.size()
            << " pending sparse temporal-spreading contribution(s) after "
            << "finite-area face addressing changed. The pending local face "
            << "indices are no longer safe to use."
            << abort(FatalError);
    }
}


void temporalSpreading::updateState()
{
    if (sparse_.empty())
    {
        sparseFaceLabels_ = mesh_.faceLabels();
    }
    else
    {
        checkSparseFaceAddressing();
    }

    events_.squeezeNull();

    label nEvents = 0;

    forAll(events_, eventi)
    {
        if (events_.get(eventi))
        {
            ++nEvents;
        }
    }

    labelList eventIds(nEvents);
    scalarList eventStartTimes(nEvents);
    label eventi = 0;

    forAll(events_, eventSlot)
    {
        if (events_.get(eventSlot))
        {
            eventIds[eventi] = events_[eventSlot].id();
            eventStartTimes[eventi] = events_[eventSlot].startTime();
            ++eventi;
        }
    }

    state_.set("modelType", modelType_);
    state_.set("temporalKernel", kernelPtr_->type());
    state_.set("temporalKernelSignature", temporalKernelSignature_);
    state_.set("nextEventId", nextEventId_);
    state_.set("eventIds", eventIds);
    state_.set("eventStartTimes", eventStartTimes);
    state_.set("stateVersion", 4);
    state_.set("sparseFaceLabels", sparseFaceLabels_);
    state_.set("sparseMassDimensions", dimMass);
    state_.set("sparseMomentumDimensions", dimMass*dimVelocity);
    state_.set("sparseNormalImpulseDimensions", dimMass*dimVelocity);
    state_.set("sparseDurationDimensions", dimTime);

    labelList faces(sparse_.size());
    scalarList starts(sparse_.size());
    scalarList durations(sparse_.size());
    scalarList masses(sparse_.size());
    vectorField momenta(sparse_.size());
    scalarList impulses(sparse_.size());
    forAll(sparse_, i)
    {
        const pendingContribution& pending = sparse_[i];
        faces[i] = pending.transfer.face;
        starts[i] = pending.startTime;
        durations[i] = pending.transfer.duration;
        masses[i] = pending.transfer.mass;
        momenta[i] = pending.transfer.momentum;
        impulses[i] = pending.transfer.normalImpulse;
    }
    state_.set("sparseFaces", faces);
    state_.set("sparseStartTimes", starts);
    state_.set("sparseDurations", durations);
    state_.set("sparseMasses", masses);
    state_.set("sparseMomenta", momenta);
    state_.set("sparseNormalImpulses", impulses);
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
            << " CFD steps, below "
            << "temporalSpreading.minStepsPerDuration "
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
    temporalKernelSignature_
    (
        temporalKernelSignature(kernelPtr_->type(), dict)
    ),
    minStepsPerDuration_(-1),
    state_(stateIO(mesh), word::null),
    nextEventId_(0),
    events_(),
    sparse_(),
    sparseFaceLabels_()
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
}

scalar temporalSpreading::calculateDurationAt
(
    const label face,
    const scalar mass,
    const vector& momentum,
    const scalar normalImpulse,
    const scalar diameter,
    const scalar rho,
    const scalar sigma
) const
{
    FatalErrorInFunction
        << "Selected temporalSpreading model " << modelType_
        << " does not support individual impact durations"
        << abort(FatalError);
    return 0;
}


scalar temporalSpreading::durationAt
(
    const label face,
    const scalar mass,
    const vector& momentum,
    const scalar normalImpulse,
    const scalar diameter,
    const scalar rho,
    const scalar sigma
) const
{
    if (face < 0 || face >= mesh_.nFaces())
    {
        FatalErrorInFunction
            << "Invalid duration source face " << face << abort(FatalError);
    }
    const scalar duration = calculateDurationAt
    (
        face, mass, momentum, normalImpulse, diameter, rho, sigma
    );
    if (!(duration > VSMALL) || !std::isfinite(duration))
    {
        FatalErrorInFunction
            << "Invalid duration " << duration << " on impact face " << face
            << abort(FatalError);
    }
    return duration;
}


void temporalSpreading::enqueueSparse
(
    const UList<contribution>& transfers,
    const scalar startTime
)
{
    if (!std::isfinite(startTime))
    {
        FatalErrorInFunction
            << "Non-finite temporal-spreading start time " << startTime
            << abort(FatalError);
    }

    checkSparseFaceAddressing();

    scalar localMinDuration = GREAT;
    bool localActive = false;
    for (const contribution& transfer : transfers)
    {
        if
        (
            transfer.face < 0
         || transfer.face >= mesh_.nFaces()
         || !(transfer.duration > VSMALL)
         || !std::isfinite(transfer.duration)
         || !std::isfinite(transfer.mass)
         || !std::isfinite(mag(transfer.momentum))
         || !std::isfinite(transfer.normalImpulse)
        )
        {
            FatalErrorInFunction
                << "Invalid sparse temporal-spreading contribution at face "
                << transfer.face << abort(FatalError);
        }
        if
        (
            mag(transfer.mass) <= VSMALL
         && mag(transfer.momentum) <= VSMALL
         && mag(transfer.normalImpulse) <= VSMALL
        )
        {
            continue;
        }
        if (sparse_.empty())
        {
            sparseFaceLabels_ = mesh_.faceLabels();
        }
        sparse_.append(pendingContribution{transfer, startTime});
        localActive = true;
        localMinDuration = Foam::min(localMinDuration, transfer.duration);
    }
    if (minStepsPerDuration_ > 0 && returnReduceOr(localActive))
    {
        const scalar minDuration =
            returnReduce(localMinDuration, minOp<scalar>());
        const scalar steps = minDuration/mesh_.time().deltaTValue();
        if (steps < minStepsPerDuration_)
        {
            WarningInFunction
                << "temporalSpreading duration is resolved by " << steps
                << " CFD steps, below "
                << "temporalSpreading.minStepsPerDuration "
                << minStepsPerDuration_
                << ". The physical duration is unchanged." << endl;
        }
    }
    updateState();
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

    checkResolution(duration, mass, momentum, normalImpulse);

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

    checkSparseFaceAddressing();

    mass.primitiveFieldRef() = Zero;
    momentum.primitiveFieldRef() = Zero;
    normalImpulse.primitiveFieldRef() = Zero;

    forAll(events_, eventi)
    {
        if (!events_.get(eventi))
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

            const scalar fraction = releaseFraction
            (
                *kernelPtr_,
                beginTime - pending.startTime(),
                endTime - pending.startTime(),
                eventDuration
            );

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

    label retained = 0;
    for (pendingContribution& pending : sparse_)
    {
        contribution& transfer = pending.transfer;
        const scalar fraction = releaseFraction
        (
            *kernelPtr_,
            beginTime - pending.startTime,
            endTime - pending.startTime,
            transfer.duration
        );
        const scalar releasedMass = fraction*transfer.mass;
        const vector releasedMomentum = fraction*transfer.momentum;
        const scalar releasedImpulse = fraction*transfer.normalImpulse;
        mass[transfer.face] += releasedMass;
        momentum[transfer.face] += releasedMomentum;
        normalImpulse[transfer.face] += releasedImpulse;

        if (fraction != 1)
        {
            transfer.mass -= releasedMass;
            transfer.momentum -= releasedMomentum;
            transfer.normalImpulse -= releasedImpulse;
            if
            (
                mag(transfer.mass) > 0
             || mag(transfer.momentum) > 0
             || mag(transfer.normalImpulse) > 0
            )
            {
                sparse_[retained++] = pending;
            }
        }
    }
    sparse_.setSize(retained);

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
