/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O  peration     |
    \\  /    A  nd           |
     \\/     M anipulation   |
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

#include "spatialSpreading.H"
#include "faceCoverage.H"
#include "footprintDistance.H"
#include "footprintKernel.H"
#include "footprintShape.H"
#include "Pstream.H"
#include "Map.H"

#include <cmath>

// * * * * * * * * * * * * * * * Local Data * * * * * * * * * * * * * * * * //

namespace
{

struct sourceData
{
    Foam::label sourceProc = -1;
    Foam::label sourceId = -1;
    Foam::label sourceFace = -1;
    Foam::label component = -1;
    Foam::label hasFootprint = 0;
    Foam::point centre = Foam::Zero;
    Foam::vector normal = Foam::Zero;
    Foam::vector majorAxis = Foam::Zero;
    Foam::vector minorAxis = Foam::Zero;
    Foam::scalar majorRadius = 0;
    Foam::scalar minorRadius = 0;
    Foam::scalar searchRadius = 0;
    Foam::scalar mass = 0;
    Foam::vector momentum = Foam::Zero;
    Foam::scalar normalImpulse = 0;
    Foam::scalar temporalWeight = 0;
    Foam::scalar temporalMoment = 0;
    Foam::scalar duration = 0;

    bool operator==(const sourceData& rhs) const
    {
        return
        (
            sourceProc == rhs.sourceProc
         && sourceId == rhs.sourceId
         && sourceFace == rhs.sourceFace
         && component == rhs.component
         && hasFootprint == rhs.hasFootprint
         && centre == rhs.centre
         && normal == rhs.normal
         && majorAxis == rhs.majorAxis
         && minorAxis == rhs.minorAxis
         && majorRadius == rhs.majorRadius
         && minorRadius == rhs.minorRadius
         && searchRadius == rhs.searchRadius
         && mass == rhs.mass
         && momentum == rhs.momentum
         && normalImpulse == rhs.normalImpulse
         && temporalWeight == rhs.temporalWeight
         && temporalMoment == rhs.temporalMoment
         && duration == rhs.duration
        );
    }

    bool operator!=(const sourceData& rhs) const
    {
        return !(*this == rhs);
    }
};


struct preparedSource
{
    sourceData source;
    Foam::DynamicList<Foam::label> targetFaces;
    Foam::DynamicList<Foam::scalar> rawWeights;
    Foam::scalar localWeightSum = 0;
    Foam::scalar globalWeightSum = 0;
};


Foam::Ostream& operator<<(Foam::Ostream& os, const sourceData& source)
{
    os  << source.sourceProc
        << source.sourceId
        << source.sourceFace
        << source.component
        << source.hasFootprint
        << source.centre
        << source.normal
        << source.majorAxis
        << source.minorAxis
        << source.majorRadius
        << source.minorRadius
        << source.searchRadius
        << source.mass
        << source.momentum
        << source.normalImpulse
        << source.temporalWeight
        << source.temporalMoment
        << source.duration;

    return os;
}


Foam::Istream& operator>>(Foam::Istream& is, sourceData& source)
{
    is  >> source.sourceProc
        >> source.sourceId
        >> source.sourceFace
        >> source.component
        >> source.hasFootprint
        >> source.centre
        >> source.normal
        >> source.majorAxis
        >> source.minorAxis
        >> source.majorRadius
        >> source.minorRadius
        >> source.searchRadius
        >> source.mass
        >> source.momentum
        >> source.normalImpulse
        >> source.temporalWeight
        >> source.temporalMoment
        >> source.duration;

    return is;
}


bool sourceIntersectsBounds
(
    const sourceData& source,
    const Foam::point& boundsMin,
    const Foam::point& boundsMax,
    const bool geodesic
)
{
    if (!source.hasFootprint)
    {
        return false;
    }

    // Conservative rank filtering: reject only bounds that cannot intersect
    // the geodesic sphere or the footprint's tangent-coordinate envelope.
    if (geodesic)
    {
        const Foam::point nearest
        (
            Foam::max
            (
                boundsMin.x(),
                Foam::min(source.centre.x(), boundsMax.x())
            ),
            Foam::max
            (
                boundsMin.y(),
                Foam::min(source.centre.y(), boundsMax.y())
            ),
            Foam::max
            (
                boundsMin.z(),
                Foam::min(source.centre.z(), boundsMax.z())
            )
        );
        return
            Foam::magSqr(nearest - source.centre)
         <= Foam::sqr(source.searchRadius);
    }

    Foam::scalar minMajor = Foam::GREAT;
    Foam::scalar maxMajor = -Foam::GREAT;
    Foam::scalar minMinor = Foam::GREAT;
    Foam::scalar maxMinor = -Foam::GREAT;

    for (Foam::label corneri = 0; corneri < 8; ++corneri)
    {
        const Foam::point corner
        (
            (corneri & 1) ? boundsMax.x() : boundsMin.x(),
            (corneri & 2) ? boundsMax.y() : boundsMin.y(),
            (corneri & 4) ? boundsMax.z() : boundsMin.z()
        );
        const Foam::vector delta(corner - source.centre);
        const Foam::scalar major = delta & source.majorAxis;
        const Foam::scalar minor = delta & source.minorAxis;
        minMajor = Foam::min(minMajor, major);
        maxMajor = Foam::max(maxMajor, major);
        minMinor = Foam::min(minMinor, minor);
        maxMinor = Foam::max(maxMinor, minor);
    }

    return
    (
        minMajor <= source.searchRadius
     && maxMajor >= -source.searchRadius
     && minMinor <= source.searchRadius
     && maxMinor >= -source.searchRadius
    );
}


Foam::vector tangentToFace
(
    const Foam::vector& momentum,
    const Foam::vector& normal
)
{
    return momentum - normal*(momentum & normal);
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

defineTypeNameAndDebug(spatialSpreading, 0);
defineRunTimeSelectionTable(spatialSpreading, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

spatialSpreading::spatialSpreading
(
    const word& modelType,
    const faMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    data_(mesh),
    sourceExchangePhases_(0),
    weightSumExchangePhases_(0),
    weightFactorExchangePhases_(0),
    sourceDescriptorsSent_(0),
    sourceDescriptorsReceived_(0),
    maxSourcesRetained_(0),
    candidateFacesVisited_(0)
{
    shapePtr_ = footprintShape::New(dict);
    distancePtr_ = footprintDistance::New(dict);
    kernelPtr_ = footprintKernel::New(dict);
    coveragePtr_ =faceCoverage::New(dict);
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

spatialSpreading::~spatialSpreading() = default;


// * * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

void spatialSpreading::spreadSources
(
    areaScalarField& mass,
    areaVectorField& momentum,
    areaScalarField& normalImpulse,
    const areaScalarField& diameter,
    const scalar beta,
    areaScalarField* temporalWeight,
    areaScalarField* temporalMoment,
    const UList<sourceInput>* impacts,
    const areaScalarField* legacyDuration,
    DynamicList<targetContribution>* transfers
) const
{
    if ((temporalWeight == nullptr) != (temporalMoment == nullptr))
    {
        FatalErrorInFunction
            << "Temporal weight and moment fields must be supplied together"
            << abort(FatalError);
    }

    const scalarField& localMass = mass.primitiveField();
    const vectorField& localMomentum = momentum.primitiveField();
    const scalarField& localNormal = normalImpulse.primitiveField();
    const scalarField& localDiameter = diameter.primitiveField();
    const scalarField* localTemporalWeight =
        temporalWeight ? &temporalWeight->primitiveField() : nullptr;
    const scalarField* localTemporalMoment =
        temporalMoment ? &temporalMoment->primitiveField() : nullptr;
    const scalarField* localDuration =
        legacyDuration ? &legacyDuration->primitiveField() : nullptr;
    if (transfers && !legacyDuration)
    {
        FatalErrorInFunction
            << "Timed transfers require legacy source-face durations"
            << abort(FatalError);
    }

    bool localActive = impacts && !impacts->empty();
    forAll(localMass, facei)
    {
        localActive =
            localActive
         || Foam::mag(localMass[facei]) > VSMALL
         || Foam::mag(localMomentum[facei]) > VSMALL
         || Foam::mag(localNormal[facei]) > VSMALL;
    }

    const label comm = mesh_.comm();
    const label myProc = UPstream::myProcNo(comm);
    const label nProcs = UPstream::nProcs(comm);

    label globalActive = localActive;
    reduce(globalActive, maxOp<label>(), UPstream::msgType(), comm);
    if (!globalActive)
    {
        return;
    }

    const bool needFacePoints = coveragePtr_->needsFacePoints();
    data_.update(needFacePoints);
    const labelList& components = data_.components();
    const vectorField& targetNormals = data_.normals();
    const bool geodesic = distancePtr_->usesGeodesicSearch();
    const scalar extraRadius =
        needFacePoints ? data_.maxFaceRadius() : scalar(0);

    point localBoundsMin(Zero);
    point localBoundsMax(Zero);
    const bool haveBounds =
        data_.centreBounds(localBoundsMin, localBoundsMax);
    List<List<point>> procBounds(nProcs);
    if (haveBounds)
    {
        List<point> bounds(2);
        bounds[0] = localBoundsMin;
        bounds[1] = localBoundsMax;
        procBounds[myProc] = bounds;
    }
    // Exchange only rank bounds here; full source records follow only for
    // ranks whose bounds can intersect a footprint search.
    Pstream::allGatherList(procBounds, UPstream::msgType(), comm);
    ++sourceExchangePhases_;

    DynamicList<sourceData> localSources;

    auto setFootprint = [&](sourceData& source, const scalar sourceDiameter)
    {
        if (!std::isfinite(sourceDiameter))
        {
            FatalErrorInFunction
                << "Non-finite parcel diameter " << sourceDiameter
                << " for source face " << source.sourceFace
                << abort(FatalError);
        }

        const scalar minorRadius = 0.5*beta*sourceDiameter;
        if (!std::isfinite(minorRadius))
        {
            FatalErrorInFunction
                << "Non-finite footprint radius " << minorRadius
                << " for source face " << source.sourceFace
                << " (diameter=" << sourceDiameter << ", beta=" << beta
                << ')'
                << abort(FatalError);
        }

        source.searchRadius = extraRadius;
        if (!std::isfinite(source.searchRadius))
        {
            FatalErrorInFunction
                << "Non-finite footprint search radius "
                << source.searchRadius << " for source face "
                << source.sourceFace
                << abort(FatalError);
        }

        if (minorRadius > SMALL)
        {
            const footprintShape::footprintGeometry fp =
                shapePtr_->footprintAt(data_, source.sourceFace, minorRadius);

            if
            (
                !std::isfinite(fp.majorRadius)
             || !std::isfinite(fp.minorRadius)
            )
            {
                FatalErrorInFunction
                    << "Non-finite footprint radii for source face "
                    << source.sourceFace << " (major=" << fp.majorRadius
                    << ", minor=" << fp.minorRadius << ')'
                    << abort(FatalError);
            }

            source.hasFootprint = 1;
            source.centre = fp.centre;
            source.normal = fp.normal;
            source.majorAxis = fp.majorAxis;
            source.minorAxis = fp.minorAxis;
            source.majorRadius = fp.majorRadius;
            source.minorRadius = fp.minorRadius;
            source.searchRadius += fp.enclosingRadius();
            if (!std::isfinite(source.searchRadius))
            {
                FatalErrorInFunction
                    << "Non-finite footprint search radius "
                    << source.searchRadius << " for source face "
                    << source.sourceFace
                    << abort(FatalError);
            }
        }
    };

    auto appendTransfer =
        [&](const sourceData& source, const label face, const scalar weight)
    {
        if
        (
            !transfers
         || (
                mag(source.mass) <= VSMALL
             && mag(source.momentum) <= VSMALL
             && mag(source.normalImpulse) <= VSMALL
            )
        )
        {
            return;
        }
        if (!(source.duration > VSMALL) || !std::isfinite(source.duration))
        {
            FatalErrorInFunction
                << "Invalid source duration " << source.duration
                << " for source " << source.sourceProc << ':'
                << source.sourceId << abort(FatalError);
        }
        transfers->append
        (
            targetContribution
            {
                source.sourceProc,
                source.sourceId,
                face,
                weight*source.mass,
                weight*tangentToFace(source.momentum, targetNormals[face]),
                weight*source.normalImpulse,
                source.duration
            }
        );
    };

    forAll(localMass, sourceFacei)
    {
        if
        (
            Foam::mag(localMass[sourceFacei]) <= VSMALL
         && Foam::mag(localMomentum[sourceFacei]) <= VSMALL
         && Foam::mag(localNormal[sourceFacei]) <= VSMALL
        )
        {
            continue;
        }

        sourceData source;
        source.sourceProc = myProc;
        source.sourceId = localSources.size();
        source.sourceFace = sourceFacei;
        source.component = components[sourceFacei];
        source.mass = localMass[sourceFacei];
        source.momentum = localMomentum[sourceFacei];
        source.normalImpulse = localNormal[sourceFacei];
        source.temporalWeight =
            localTemporalWeight ? (*localTemporalWeight)[sourceFacei] : 0;
        source.temporalMoment =
            localTemporalMoment ? (*localTemporalMoment)[sourceFacei] : 0;
        source.duration =
            localDuration ? (*localDuration)[sourceFacei] : 0;

        setFootprint(source, localDiameter[sourceFacei]);

        localSources.append(source);
    }

    if (impacts)
    {
        for (const sourceInput& impact : *impacts)
        {
            if
            (
                impact.face < 0
             || impact.face >= mesh_.nFaces()
             || !(impact.mass > 0)
             || !std::isfinite(impact.mass)
             || !std::isfinite(mag(impact.momentum))
             || !std::isfinite(impact.normalImpulse)
             || !(impact.diameter > 0)
             || !std::isfinite(impact.diameter)
            )
            {
                FatalErrorInFunction
                    << "Invalid parcel impact on face " << impact.face
                    << " (mass=" << impact.mass
                    << ", diameter=" << impact.diameter << ')'
                    << exit(FatalError);
            }

            sourceData source;
            source.sourceProc = myProc;
            source.sourceId = localSources.size();
            source.sourceFace = impact.face;
            source.component = components[impact.face];
            source.mass = impact.mass;
            source.momentum = impact.momentum;
            source.normalImpulse = impact.normalImpulse;
            source.duration = impact.duration;
            setFootprint(source, impact.diameter);
            localSources.append(source);
        }
    }

    List<DynamicList<sourceData>> outgoingSources(nProcs);
    for (const sourceData& source : localSources)
    {
        if (!source.hasFootprint)
        {
            continue;
        }

        for (label proci = 0; proci < nProcs; ++proci)
        {
            if
            (
                proci != myProc
             && procBounds[proci].size() == 2
             && sourceIntersectsBounds
                (
                    source,
                    procBounds[proci][0],
                    procBounds[proci][1],
                    geodesic
                )
            )
            {
                outgoingSources[proci].append(source);
                ++sourceDescriptorsSent_;
            }
        }
    }

    PstreamBuffers sourceBufs(comm);
    for (label proci = 0; proci < nProcs; ++proci)
    {
        if (proci == myProc || outgoingSources[proci].empty())
        {
            continue;
        }

        UOPstream toProc(proci, sourceBufs);
        toProc << List<sourceData>(outgoingSources[proci]);
    }
    sourceBufs.finishedSends();
    ++sourceExchangePhases_;

    DynamicList<sourceData> sources;
    for (const sourceData& source : localSources)
    {
        sources.append(source);
    }
    for (label proci = 0; proci < nProcs; ++proci)
    {
        if (proci == myProc || !sourceBufs.recvDataCount(proci))
        {
            continue;
        }

        List<sourceData> remoteSources;
        UIPstream fromProc(proci, sourceBufs);
        fromProc >> remoteSources;
        sourceDescriptorsReceived_ += remoteSources.size();
        for (const sourceData& source : remoteSources)
        {
            sources.append(source);
        }
    }
    maxSourcesRetained_ =
        Foam::max(maxSourcesRetained_, label(sources.size()));

    scalarField nextMass(mesh_.nFaces(), Zero);
    vectorField nextMomentum(mesh_.nFaces(), Zero);
    scalarField nextNormal(mesh_.nFaces(), Zero);
    scalarField nextTemporalWeight(mesh_.nFaces(), Zero);
    scalarField nextTemporalMoment(mesh_.nFaces(), Zero);

    DynamicList<label> geodesicSourceProcs;
    DynamicList<label> geodesicSourceIds;
    DynamicList<label> geodesicSourceFaces;
    DynamicList<scalar> geodesicSearchRadii;
    labelList geodesicSourceIndex(sources.size(), -1);
    List<labelList> geodesicCandidates;
    List<scalarField> geodesicDistances;

    if (geodesic)
    {
        forAll(sources, sourcei)
        {
            if (sources[sourcei].hasFootprint)
            {
                geodesicSourceIndex[sourcei] = geodesicSourceProcs.size();
                geodesicSourceProcs.append(sources[sourcei].sourceProc);
                geodesicSourceIds.append(sources[sourcei].sourceId);
                geodesicSourceFaces.append(sources[sourcei].sourceFace);
                geodesicSearchRadii.append(sources[sourcei].searchRadius);
            }
        }

        label globalHasGeodesicSources = !geodesicSourceProcs.empty();
        reduce
        (
            globalHasGeodesicSources,
            maxOp<label>(),
            UPstream::msgType(),
            comm
        );
        if (globalHasGeodesicSources)
        {
            data_.geodesicCandidatesBatch
            (
                geodesicSourceProcs,
                geodesicSourceIds,
                geodesicSourceFaces,
                scalarField(geodesicSearchRadii),
                geodesicCandidates,
                geodesicDistances
            );
        }
    }

    DynamicList<preparedSource> prepared;
    List<Map<label>> preparedIndices(nProcs);

    forAll(sources, sourcei)
    {
        const sourceData& source = sources[sourcei];
        if (!source.hasFootprint)
        {
            if (source.sourceProc == myProc)
            {
                nextMass[source.sourceFace] += source.mass;
                nextMomentum[source.sourceFace] += tangentToFace
                (
                    source.momentum,
                    targetNormals[source.sourceFace]
                );
                nextNormal[source.sourceFace] += source.normalImpulse;
                nextTemporalWeight[source.sourceFace] +=
                    source.temporalWeight;
                nextTemporalMoment[source.sourceFace] +=
                    source.temporalMoment;
                appendTransfer(source, source.sourceFace, 1);
            }
            continue;
        }

        footprintShape::footprintGeometry fp;
        fp.centre = source.centre;
        fp.normal = source.normal;
        fp.majorAxis = source.majorAxis;
        fp.minorAxis = source.minorAxis;
        fp.majorRadius = source.majorRadius;
        fp.minorRadius = source.minorRadius;

        labelList candidateFaces;
        scalarField selectedDistances;
        if (geodesic)
        {
            const label geoIndex = geodesicSourceIndex[sourcei];
            if (geoIndex < 0)
            {
                FatalErrorInFunction
                    << "Missing batched geodesic source index for "
                    << source.sourceProc << ':' << source.sourceId
                    << exit(FatalError);
            }
            candidateFaces = geodesicCandidates[geoIndex];
            selectedDistances = geodesicDistances[geoIndex];
        }
        else
        {
            DynamicList<label> localCandidates;
            DynamicList<scalar> localDistances;
            distancePtr_->findCandidates
            (
                data_,
                *shapePtr_,
                fp,
                source.sourceProc,
                source.sourceFace,
                source.searchRadius,
                localCandidates,
                localDistances
            );
            candidateFaces = localCandidates;
            selectedDistances = localDistances;
        }
        candidateFacesVisited_ += candidateFaces.size();

        preparedSource localSource;
        localSource.source = source;

        forAll(candidateFaces, candidatei)
        {
            const label targetFace = candidateFaces[candidatei];

            if (components[targetFace] != source.component)
            {
                continue;
            }

            scalar x = 0;
            scalar y = 0;
            shapePtr_->planarOffset
            (
                fp,
                data_.centres()[targetFace],
                x,
                y
            );

            const scalar planarDistance =
                Foam::sqrt(Foam::sqr(x) + Foam::sqr(y));
            const scalar boundaryRadius =
                shapePtr_->boundaryRadius
                (
                    x,
                    y,
                    fp.majorRadius,
                    fp.minorRadius
                );

            const scalar selectedDistance = selectedDistances[candidatei];
            const faceCoverage::coverageResult coverage =
                coveragePtr_->evaluate
                (
                    data_,
                    fp,
                    targetFace,
                    selectedDistance,
                    boundaryRadius
                );

            if (!coverage.included)
            {
                continue;
            }

            scalar kernelDistance = selectedDistance;
            if (coveragePtr_->usesProjectedCentroid())
            {
                const vector centroidDelta(coverage.centre - fp.centre);
                const scalar centroidPlanarDistance = Foam::sqrt
                (
                    Foam::sqr(centroidDelta & fp.majorAxis)
                  + Foam::sqr(centroidDelta & fp.minorAxis)
                );

                kernelDistance =
                    distancePtr_->kernelDistance
                    (
                        selectedDistance,
                        planarDistance,
                        centroidPlanarDistance
                    );
            }

            const scalar rawWeight =
                kernelPtr_->weight(coverage.measure, kernelDistance);

            if (rawWeight > VSMALL)
            {
                localSource.targetFaces.append(targetFace);
                localSource.rawWeights.append(rawWeight);
                localSource.localWeightSum += rawWeight;
            }
        }

        if (localSource.localWeightSum > VSMALL)
        {
            if
            (
                preparedIndices[source.sourceProc].contains
                (
                    source.sourceId
                )
            )
            {
                FatalErrorInFunction
                    << "Duplicate prepared source "
                    << source.sourceProc << ':' << source.sourceId
                    << exit(FatalError);
            }

            preparedIndices[source.sourceProc].insert
            (
                source.sourceId,
                prepared.size()
            );
            prepared.append(std::move(localSource));
        }
    }

    scalarField globalWeightSum(localSources.size(), Zero);
    List<DynamicList<label>> sendWeightFaces(nProcs);
    List<DynamicList<scalar>> sendWeightValues(nProcs);

    for (const preparedSource& source : prepared)
    {
        if (source.source.sourceProc == myProc)
        {
            globalWeightSum[source.source.sourceId] +=
                source.localWeightSum;
        }
        else
        {
            const label owner = source.source.sourceProc;
            sendWeightFaces[owner].append(source.source.sourceId);
            sendWeightValues[owner].append(source.localWeightSum);
        }
    }

    // Batch local normalization sums back to the source owners.
    PstreamBuffers weightSumBufs(comm);
    for (label proci = 0; proci < nProcs; ++proci)
    {
        if (proci == myProc || sendWeightFaces[proci].empty())
        {
            continue;
        }

        UOPstream toProc(proci, weightSumBufs);
        toProc
            << labelList(sendWeightFaces[proci])
            << scalarField(sendWeightValues[proci]);
    }
    weightSumBufs.finishedSends();
    ++weightSumExchangePhases_;

    List<DynamicList<label>> reporterFaces(nProcs);
    for (label proci = 0; proci < nProcs; ++proci)
    {
        if (proci == myProc || !weightSumBufs.recvDataCount(proci))
        {
            continue;
        }

        labelList remoteFaces;
        scalarField remoteWeights;
        UIPstream fromProc(proci, weightSumBufs);
        fromProc >> remoteFaces >> remoteWeights;

        if (remoteFaces.size() != remoteWeights.size())
        {
            FatalErrorInFunction
                << "Inconsistent source weight sums from processor "
                << proci << exit(FatalError);
        }

        forAll(remoteFaces, i)
        {
            const label sourceId = remoteFaces[i];
            if (sourceId < 0 || sourceId >= localSources.size())
            {
                FatalErrorInFunction
                    << "Invalid source ID " << sourceId
                    << " reported to owner processor " << myProc
                    << exit(FatalError);
            }
            globalWeightSum[sourceId] += remoteWeights[i];
            reporterFaces[proci].append(sourceId);
        }
    }

    // Return one global normalizer to each rank with a positive local weight.
    PstreamBuffers weightFactorBufs(comm);
    for (label proci = 0; proci < nProcs; ++proci)
    {
        if (proci == myProc || reporterFaces[proci].empty())
        {
            continue;
        }

        scalarField globalWeights(reporterFaces[proci].size(), Zero);
        forAll(reporterFaces[proci], i)
        {
            globalWeights[i] = globalWeightSum[reporterFaces[proci][i]];
        }

        UOPstream toProc(proci, weightFactorBufs);
        toProc << labelList(reporterFaces[proci]) << globalWeights;
    }
    weightFactorBufs.finishedSends();
    ++weightFactorExchangePhases_;

    for (preparedSource& source : prepared)
    {
        if (source.source.sourceProc == myProc)
        {
            source.globalWeightSum =
                globalWeightSum[source.source.sourceId];
        }
    }

    for (label proci = 0; proci < nProcs; ++proci)
    {
        if (proci == myProc || !weightFactorBufs.recvDataCount(proci))
        {
            continue;
        }

        labelList remoteFaces;
        scalarField remoteWeights;
        UIPstream fromProc(proci, weightFactorBufs);
        fromProc >> remoteFaces >> remoteWeights;

        if (remoteFaces.size() != remoteWeights.size())
        {
            FatalErrorInFunction
                << "Inconsistent source weight factors from processor "
                << proci << exit(FatalError);
        }

        forAll(remoteFaces, i)
        {
            if (!preparedIndices[proci].contains(remoteFaces[i]))
            {
                FatalErrorInFunction
                    << "Missing prepared source "
                    << proci << ':' << remoteFaces[i]
                    << " for returned weight factor"
                    << exit(FatalError);
            }

            preparedSource& source =
                prepared[preparedIndices[proci].at(remoteFaces[i])];
            source.globalWeightSum = remoteWeights[i];
        }
    }

    for (const sourceData& source : localSources)
    {
        if
        (
            source.hasFootprint
         && globalWeightSum[source.sourceId] <= VSMALL
        )
        {
            nextMass[source.sourceFace] += source.mass;
            nextMomentum[source.sourceFace] += tangentToFace
            (
                source.momentum,
                targetNormals[source.sourceFace]
            );
            nextNormal[source.sourceFace] += source.normalImpulse;
            nextTemporalWeight[source.sourceFace] += source.temporalWeight;
            nextTemporalMoment[source.sourceFace] += source.temporalMoment;
            appendTransfer(source, source.sourceFace, 1);
        }
    }

    for (const preparedSource& source : prepared)
    {
        if (source.globalWeightSum <= VSMALL)
        {
            FatalErrorInFunction
                << "Missing global weight for source "
                << source.source.sourceProc << ':'
                << source.source.sourceId
                << exit(FatalError);
        }

        forAll(source.targetFaces, targeti)
        {
            const scalar fraction =
                source.rawWeights[targeti]/source.globalWeightSum;
            const label targetFace = source.targetFaces[targeti];

            nextMass[targetFace] += fraction*source.source.mass;
            nextMomentum[targetFace] += fraction*tangentToFace
            (
                source.source.momentum,
                targetNormals[targetFace]
            );
            nextNormal[targetFace] +=
                fraction*source.source.normalImpulse;
            nextTemporalWeight[targetFace] +=
                fraction*source.source.temporalWeight;
            nextTemporalMoment[targetFace] +=
                fraction*source.source.temporalMoment;
            appendTransfer(source.source, targetFace, fraction);
        }
    }

    mass.primitiveFieldRef() = nextMass;
    momentum.primitiveFieldRef() = nextMomentum;
    normalImpulse.primitiveFieldRef() = nextNormal;

    mass.correctBoundaryConditions();
    momentum.correctBoundaryConditions();
    normalImpulse.correctBoundaryConditions();

    if (temporalWeight)
    {
        temporalWeight->primitiveFieldRef() = nextTemporalWeight;
        temporalMoment->primitiveFieldRef() = nextTemporalMoment;
        temporalWeight->correctBoundaryConditions();
        temporalMoment->correctBoundaryConditions();
    }
}


void spatialSpreading::spread
(
    areaScalarField& mass,
    areaVectorField& momentum,
    areaScalarField& normalImpulse,
    const areaScalarField& diameter,
    const UList<sourceInput>& impacts
) const
{
    FatalErrorInFunction
        << "The selected spatialSpreading model does not support "
        << "individual parcel impacts."
        << abort(FatalError);
}


void spatialSpreading::spread
(
    areaScalarField& mass,
    areaVectorField& momentum,
    areaScalarField& normalImpulse,
    areaScalarField& temporalWeight,
    areaScalarField& temporalMoment,
    const areaScalarField& diameter
) const
{
    FatalErrorInFunction
        << "The selected spatialSpreading model does not support "
        << "temporal duration moments."
        << abort(FatalError);
}

void spatialSpreading::spread
(
    areaScalarField& mass,
    areaVectorField& momentum,
    areaScalarField& normalImpulse,
    const areaScalarField& diameter,
    const areaScalarField& legacyDuration,
    const UList<sourceInput>& impacts,
    DynamicList<targetContribution>& transfers
) const
{
    FatalErrorInFunction
        << "The selected spatialSpreading model does not support "
        << "individually timed target contributions."
        << abort(FatalError);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
