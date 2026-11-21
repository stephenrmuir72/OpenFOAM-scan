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

#include "spatialSpreading.H"
#include "faceCoverage.H"
#include "footprintDistance.H"
#include "footprintKernel.H"
#include "footprintShape.H"

// * * * * * * * * * * * * * * * Local Data * * * * * * * * * * * * * * * * //

namespace
{

// Data unit containing the data of sources waiting to be spread
struct sourceData
{
    bool hasFootprint = false;
    Foam::label sourceProc = -1;
    Foam::label sourceId = -1;
    Foam::label sourceFacei = -1;
    Foam::label component = -1;
    Foam::regionModels::areaSurfaceFilmModels::footprintShape::
        footprintGeometry footprint;
    Foam::scalar searchRadius = 0;
    Foam::scalar mass = 0;
    Foam::vector momentum = Foam::zero{};
    Foam::scalar normalImpulse = 0;
};

// Data unit containing 'sourceData' with target faces and weights
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
    os  << source.hasFootprint
        << source.sourceProc
        << source.sourceId
        << source.sourceFacei
        << source.component
        << source.footprint.centre
        << source.footprint.normal
        << source.footprint.majorAxis
        << source.footprint.minorAxis
        << source.footprint.majorRadius
        << source.footprint.minorRadius
        << source.searchRadius
        << source.mass
        << source.momentum
        << source.normalImpulse;

    return os;
}


Foam::Istream& operator>>(Foam::Istream& is, sourceData& source)
{
    is  >> source.hasFootprint
        >> source.sourceProc
        >> source.sourceId
        >> source.sourceFacei
        >> source.component
        >> source.footprint.centre
        >> source.footprint.normal
        >> source.footprint.majorAxis
        >> source.footprint.minorAxis
        >> source.footprint.majorRadius
        >> source.footprint.minorRadius
        >> source.searchRadius
        >> source.mass
        >> source.momentum
        >> source.normalImpulse;

    return is;
}


bool sourceIntersectsBounds
(
    const sourceData& source,
    const Foam::boundBox& bounds
)
{
    if (!source.hasFootprint)
    {
        return false;
    }

    // Conservative rank filtering: reject only bounds that cannot intersect
    // the footprint's tangent-coordinate envelope.
    Foam::scalar minMajor = Foam::GREAT;
    Foam::scalar maxMajor = -Foam::GREAT;
    Foam::scalar minMinor = Foam::GREAT;
    Foam::scalar maxMinor = -Foam::GREAT;

    for (Foam::label corneri = 0; corneri < 8; ++corneri)
    {
        const Foam::point corner
        (
            bounds.octCorner(corneri)
        );
        const Foam::vector delta(corner - source.footprint.centre);
        const Foam::scalar major = delta.inner(source.footprint.majorAxis);
        const Foam::scalar minor = delta.inner(source.footprint.minorAxis);
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
    data_(mesh)
{
    shapePtr_ = footprintShape::New(dict_);
    distancePtr_ = footprintDistance::New(dict_);
    kernelPtr_ = footprintKernel::New(dict_);
    coveragePtr_ = faceCoverage::New(dict_);
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

spatialSpreading::~spatialSpreading() = default;


// * * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

void spatialSpreading::spreadImpl
(
    areaScalarField& mass,
    areaVectorField& momentum,
    areaScalarField& normalImpulse,
    const areaScalarField& diameter,
    const scalar beta,
    const UList<sourceInput>& impacts
) const
{
    const scalarField& localMass = mass.primitiveField();
    const vectorField& localMomentum = momentum.primitiveField();
    const scalarField& localNormal = normalImpulse.primitiveField();
    const scalarField& localDiameter = diameter.primitiveField();


    // Check whether any processor has anything to spread
    bool localActive = !impacts.empty();
    forAll(localMass, facei)
    {
        localActive =
            localActive
         || Foam::mag(localMass[facei]) > VSMALL
         || Foam::mag(localMomentum[facei]) > VSMALL
         || Foam::mag(localNormal[facei]) > VSMALL;
    }

    const auto comm = mesh_.comm();
    const auto myProci = UPstream::myProcNo(comm);
    const auto numProc = UPstream::nProcs(comm);

    if
    (
        const bool globalActive =
            returnReduce(localActive, orOp<bool>(), UPstream::msgType(), comm);
        !globalActive
    )
    {
        // If nobody has a source to spread, return now.
        return;
    }


    // Prepare surface information and exchange small maps
    const bool needFacePoints = coveragePtr_->needsFacePoints();
    data_.prepare(needFacePoints);
    const scalar extraRadius = needFacePoints ? data_.maxFaceRadius() : 0;
    const labelList& components = data_.components();
    const vectorField& targetNormals = data_.normals();


    List<boundBox> procBounds(numProc);
    procBounds[myProci] = data_.centreBounds();
    // Exchange only rank bounds here; full source records follow only for
    // ranks whose bounds can intersect a footprint search later
    Pstream::allGatherList(procBounds, UPstream::msgType(), comm);


    // Set footprint data per source
    DynamicList<sourceData> localSources;

    auto setFootprint = [&](sourceData& source, const scalar sourceDiameter)
    {
        if (debug && !std::isfinite(sourceDiameter))
        {
            FatalErrorInFunction
                << "Non-finite parcel diameter " << sourceDiameter
                << " for source face " << source.sourceFacei
                << abort(FatalError);
        }

        const scalar minorRadius = 0.5*beta*sourceDiameter;
        if (debug && !std::isfinite(minorRadius))
        {
            FatalErrorInFunction
                << "Non-finite footprint radius " << minorRadius
                << " for source face " << source.sourceFacei
                << " (diameter=" << sourceDiameter << ", beta=" << beta
                << ')'
                << abort(FatalError);
        }

        source.searchRadius = extraRadius;
        if (debug && !std::isfinite(source.searchRadius))
        {
            FatalErrorInFunction
                << "Non-finite footprint search radius "
                << source.searchRadius << " for source face "
                << source.sourceFacei
                << abort(FatalError);
        }

        // If radius is too small, the source is left without a footprint and
        // will stay on its original face without any spread
        if (minorRadius > SMALL)
        {
            // Instruct footprint-shape model to create footprint at source face
            const footprintShape::footprintGeometry fp =
                shapePtr_->footprintAt(data_, source.sourceFacei, minorRadius);

            if
            (
                debug
             && (
                    !std::isfinite(fp.majorRadius)
                 || !std::isfinite(fp.minorRadius)
                )
            )
            {
                FatalErrorInFunction
                    << "Non-finite footprint radii for source face "
                    << source.sourceFacei << " (major=" << fp.majorRadius
                    << ", minor=" << fp.minorRadius << ')'
                    << abort(FatalError);
            }

            source.hasFootprint = true;
            source.footprint = fp;
            source.searchRadius += fp.enclosingRadius();

            if (debug && !std::isfinite(source.searchRadius))
            {
                FatalErrorInFunction
                    << "Non-finite footprint search radius "
                    << source.searchRadius << " for source face "
                    << source.sourceFacei
                    << abort(FatalError);
            }
        }
    };

    auto appendSource =
    [&]
    (
        const label sourceFacei,
        const scalar sourceMass,
        const vector& sourceMomentum,
        const scalar sourceNormalImpulse,
        const scalar sourceDiameter
    )
    {
        sourceData source;
        source.sourceProc = myProci;
        source.sourceId = localSources.size();
        source.sourceFacei = sourceFacei;
        source.component = components[sourceFacei];
        source.mass = sourceMass;
        source.momentum = sourceMomentum;
        source.normalImpulse = sourceNormalImpulse;

        setFootprint(source, sourceDiameter);
        localSources.push_back(source);
    };


    // Turn input values into source records
    forAll(localMass, sourceFacei)
    {
        if
        (
            Foam::mag(localMass[sourceFacei]) <= VSMALL
         && Foam::mag(localMomentum[sourceFacei]) <= VSMALL
         && Foam::mag(localNormal[sourceFacei]) <= VSMALL
        )
        {
            // There is nothing meaningful there to spread; skip
            continue;
        }

        appendSource
        (
            sourceFacei,
            localMass[sourceFacei],
            localMomentum[sourceFacei],
            localNormal[sourceFacei],
            localDiameter[sourceFacei]
        );
    }


    // Optionally append individual parcel impacts to localSources
    if (!impacts.empty())
    {
        for (const sourceInput& impact : impacts)
        {
            if
            (
                debug
             && (
                    impact.face < 0
                 || impact.face >= mesh_.nFaces()
                 || !(impact.mass > 0)
                 || !std::isfinite(impact.mass)
                 || !std::isfinite(mag(impact.momentum))
                 || !std::isfinite(impact.normalImpulse)
                 || !(impact.diameter > 0)
                 || !std::isfinite(impact.diameter)
                )
            )
            {
                FatalErrorInFunction
                    << "Invalid parcel impact on face " << impact.face
                    << " (mass=" << impact.mass
                    << ", diameter=" << impact.diameter << ')'
                    << exit(FatalError);
            }

            appendSource
            (
                impact.face,
                impact.mass,
                impact.momentum,
                impact.normalImpulse,
                impact.diameter
            );
        }
    }


    // Send each footprint source only to procs that may need it
    List<DynamicList<sourceData>> outgoingSources(numProc);
    for (const sourceData& source : localSources)
    {
        if (!source.hasFootprint)
        {
            // Owner will handle this source; skip
            continue;
        }

        for (label proci = 0; proci < numProc; ++proci)
        {
            // Check if footprint's search area reachs bounding box per proc
            if
            (
                proci != myProci
             && sourceIntersectsBounds
                (
                    source,
                    procBounds[proci]
                )
            )
            {
                outgoingSources[proci].push_back(source);
            }
        }
    }

    PstreamBuffers sourceBufs(comm);
    for (label proci = 0; proci < numProc; ++proci)
    {
        if (proci == myProci || outgoingSources[proci].empty())
        {
            continue;
        }

        UOPstream toProc(proci, sourceBufs);
        toProc << List<sourceData>(outgoingSources[proci]);
    }
    sourceBufs.finishedSends();


    DynamicList<sourceData> sources;
    for (const sourceData& source : localSources)
    {
        sources.push_back(source);
    }
    for (label proci = 0; proci < numProc; ++proci)
    {
        if (proci == myProci || !sourceBufs.recvDataCount(proci))
        {
            continue;
        }

        List<sourceData> remoteSources;
        UIPstream fromProc(proci, sourceBufs);
        fromProc >> remoteSources;
        for (const sourceData& source : remoteSources)
        {
            sources.push_back(source);
        }
    }


    scalarField nextMass(mesh_.nFaces(), zero{});
    vectorField nextMomentum(mesh_.nFaces(), zero{});
    scalarField nextNormal(mesh_.nFaces(), zero{});

    auto addSourceToFace = [&]
    (
        const sourceData& source,
        const label targetFacei,
        const scalar fraction
    )
    {
        nextMass[targetFacei] += fraction*source.mass;

        vector tangentToFace(source.momentum);
        tangentToFace.removeCollinear(targetNormals[targetFacei]);
        nextMomentum[targetFacei] += fraction*tangentToFace;

        nextNormal[targetFacei] += fraction*source.normalImpulse;
    };


    // Prepare and spread sources
    DynamicList<preparedSource> preparedSources;
    List<Map<label>> preparedIndices(numProc);

    forAll(sources, sourcei)
    {
        const sourceData& source = sources[sourcei];
        if (!source.hasFootprint)
        {
            if (source.sourceProc == myProci)
            {
                addSourceToFace(source, source.sourceFacei, scalar(1));
            }
            continue;
        }

        const footprintShape::footprintGeometry& fp = source.footprint;

        // Find local candidate faces and dists with footprint-distance model
        labelList candidateFaces;
        scalarField candidateDistances;
        distancePtr_->findCandidates
        (
            data_,
            fp,
            source.sourceProc,
            source.searchRadius,
            candidateFaces,
            candidateDistances
        );

        preparedSource sourceInPrep;
        sourceInPrep.source = source;

        forAll(candidateFaces, candidatei)
        {
            const label targetFacei = candidateFaces[candidatei];

            if (components[targetFacei] != source.component)
            {
                // Skip candidate if its connected-piece label does not match
                // the source's to stop spreading across disconnected parts
                continue;
            }

            // Find face's place on the footprint
            scalar x = 0;
            scalar y = 0;
            footprintShape::planarOffset
            (
                fp,
                data_.centres()[targetFacei],
                x,
                y
            );

            // How far the footprint reaches in x, y dirs - target face centre
            const scalar planarDistance = Foam::hypot(x, y);
            const scalar footprintRadius =
                shapePtr_->footprintRadius
                (
                    x,
                    y,
                    fp.majorRadius,
                    fp.minorRadius
                );

            // Calc how much of the footprint covers this face
            const scalar candidateDistance = candidateDistances[candidatei];
            const faceCoverage::coverageResult coverage =
                coveragePtr_->evaluate
                (
                    data_,
                    fp,
                    targetFacei,
                    candidateDistance,
                    footprintRadius
                );

            if (!coverage.included)
            {
                // If this face outside the coverage; skip
                continue;
            }

            scalar kernelDistance = candidateDistance;
            // If the coverage has a centre point which is not this face centre,
            // use the projected centre for footprint-kernel calculation
            if (coveragePtr_->usesProjectedCentroid())
            {
                scalar centroidX = 0;
                scalar centroidY = 0;
                footprintShape::planarOffset
                (
                    fp,
                    coverage.centre,
                    centroidX,
                    centroidY
                );
                const scalar centroidPlanarDistance =
                    Foam::hypot(centroidX, centroidY);

                kernelDistance =
                    distancePtr_->kernelDistance
                    (
                        candidateDistance,
                        planarDistance,
                        centroidPlanarDistance
                    );
            }

            const scalar rawWeight =
                kernelPtr_->weight(coverage.measure, kernelDistance);

            if (rawWeight > VSMALL)
            {
                sourceInPrep.targetFaces.push_back(targetFacei);
                sourceInPrep.rawWeights.push_back(rawWeight);
                sourceInPrep.localWeightSum += rawWeight;
            }
        }

        // Save the source if it found a target
        if (sourceInPrep.localWeightSum > VSMALL)
        {
            if
            (
                debug
             && preparedIndices[source.sourceProc].contains(source.sourceId)
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
                preparedSources.size()
            );
            preparedSources.push_back(std::move(sourceInPrep));
        }
    }


    scalarField globalWeightSum(localSources.size(), zero{});
    List<DynamicList<label>> sendWeightFaces(numProc);
    List<DynamicList<scalar>> sendWeightValues(numProc);

    for (const preparedSource& psrc : preparedSources)
    {
        if (psrc.source.sourceProc == myProci)
        {
            globalWeightSum[psrc.source.sourceId] += psrc.localWeightSum;
        }
        else
        {
            const label owner = psrc.source.sourceProc;
            sendWeightFaces[owner].push_back(psrc.source.sourceId);
            sendWeightValues[owner].push_back(psrc.localWeightSum);
        }
    }


    // Batch local normalization sums back to the source owners
    PstreamBuffers weightSumBufs(comm);
    for (label proci = 0; proci < numProc; ++proci)
    {
        if (proci == myProci || sendWeightFaces[proci].empty())
        {
            continue;
        }

        UOPstream toProc(proci, weightSumBufs);
        toProc
            << labelList(sendWeightFaces[proci])
            << scalarField(sendWeightValues[proci]);
    }
    weightSumBufs.finishedSends();

    // Record which source IDs each proc reported; the owner will need this
    // list to send the final totals back
    List<DynamicList<label>> reporterFaces(numProc);
    for (label proci = 0; proci < numProc; ++proci)
    {
        if (proci == myProci || !weightSumBufs.recvDataCount(proci))
        {
            continue;
        }

        labelList remoteFaces;
        scalarField remoteWeights;
        UIPstream fromProc(proci, weightSumBufs);
        fromProc >> remoteFaces >> remoteWeights;

        if (debug && (remoteFaces.size() != remoteWeights.size()))
        {
            FatalErrorInFunction
                << "Inconsistent source weight sums from processor "
                << proci << exit(FatalError);
        }

        forAll(remoteFaces, i)
        {
            const label sourceId = remoteFaces[i];
            if (debug && (sourceId < 0 || sourceId >= localSources.size()))
            {
                FatalErrorInFunction
                    << "Invalid source ID " << sourceId
                    << " reported to owner processor " << myProci
                    << exit(FatalError);
            }
            globalWeightSum[sourceId] += remoteWeights[i];
            reporterFaces[proci].push_back(sourceId);
        }
    }


    // Return one global normalizer to each rank with a positive local weight
    PstreamBuffers weightFactorBufs(comm);
    for (label proci = 0; proci < numProc; ++proci)
    {
        if (proci == myProci || reporterFaces[proci].empty())
        {
            continue;
        }

        scalarField globalWeights(reporterFaces[proci].size(), zero{});
        // Put the totals in the same order
        forAll(reporterFaces[proci], i)
        {
            globalWeights[i] = globalWeightSum[reporterFaces[proci][i]];
        }

        UOPstream toProc(proci, weightFactorBufs);
        toProc << labelList(reporterFaces[proci]) << globalWeights;
    }
    weightFactorBufs.finishedSends();

    for (preparedSource& psrc : preparedSources)
    {
        if (psrc.source.sourceProc == myProci)
        {
            psrc.globalWeightSum = globalWeightSum[psrc.source.sourceId];
        }
    }

    for (label proci = 0; proci < numProc; ++proci)
    {
        if (proci == myProci || !weightFactorBufs.recvDataCount(proci))
        {
            continue;
        }

        labelList remoteFaces;
        scalarField remoteWeights;
        UIPstream fromProc(proci, weightFactorBufs);
        fromProc >> remoteFaces >> remoteWeights;

        if (debug && (remoteFaces.size() != remoteWeights.size()))
        {
            FatalErrorInFunction
                << "Inconsistent source weight factors from processor "
                << proci << exit(FatalError);
        }

        forAll(remoteFaces, i)
        {
            if (debug && !preparedIndices[proci].contains(remoteFaces[i]))
            {
                FatalErrorInFunction
                    << "Missing prepared source "
                    << proci << ':' << remoteFaces[i]
                    << " for returned weight factor"
                    << exit(FatalError);
            }

            preparedSource& psrc =
                preparedSources[preparedIndices[proci].at(remoteFaces[i])];
            psrc.globalWeightSum = remoteWeights[i];
        }
    }


    // Apply weights to sources waiting to be spread
    for (const sourceData& source : localSources)
    {
        // Use the original face if no target gets a useful weight
        if
        (
            source.hasFootprint
         && globalWeightSum[source.sourceId] <= VSMALL
        )
        {
            addSourceToFace(source, source.sourceFacei, scalar(1));
        }
    }

    for (const preparedSource& psrc : preparedSources)
    {
        if (debug && (psrc.globalWeightSum <= VSMALL))
        {
            FatalErrorInFunction
                << "Missing global weight for source "
                << psrc.source.sourceProc << ':'
                << psrc.source.sourceId
                << exit(FatalError);
        }

        forAll(psrc.targetFaces, i)
        {
            const scalar fraction = psrc.rawWeights[i]/psrc.globalWeightSum;
            const label targetFacei = psrc.targetFaces[i];

            addSourceToFace(psrc.source, targetFacei, fraction);
        }
    }

    mass.primitiveFieldRef() = nextMass;
    momentum.primitiveFieldRef() = nextMomentum;
    normalImpulse.primitiveFieldRef() = nextNormal;

    mass.correctBoundaryConditions();
    momentum.correctBoundaryConditions();
    normalImpulse.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
