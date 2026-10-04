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

#include "spatialSpreadingData.H"
#include "areaFields.H"
#include "Map.H"
#include "pointIOField.H"
#include <queue>
#include <tuple>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

using batchQueueEntry =
    std::tuple<Foam::scalar, Foam::label, Foam::label>;
using batchQueue =
    std::priority_queue
    <
        batchQueueEntry,
        std::vector<batchQueueEntry>,
        std::greater<batchQueueEntry>
    >;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

spatialSpreadingData::spatialSpreadingData(const faMesh& mesh)
:
    mesh_(mesh),
    adjacency_(),
    boundaryEdges_(),
    boundaryEdgesValid_(false),
    localComponents_(),
    centreTreePtr_(nullptr),
    centreBoundsMin_(Zero),
    centreBoundsMax_(Zero),
    centreBoundsValid_(false),
    facePoints_(),
    faceRadii_(),
    maxFaceRadius_(0),
    topologyValid_(false),
    componentsValid_(false),
    facePointsValid_(false),
    nFaces_(-1),
    nPoints_(-1),
    pointEventNo_(-1),
    pointsInstance_(),
    facesInstance_(),
    tangentBoundCornerChecks_(0),
    centreBoundsBuilds_(0),
    centreBoundsFaceVisits_(0),
    geodesicExchangePhases_(0)
{}


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

bool spatialSpreadingData::topologyCurrent() const
{
    return
    (
        topologyValid_
     && nFaces_ == mesh_.nFaces()
     && nPoints_ == mesh_.nPoints()
     && facesInstance_ == mesh_.facesInstance()
    );
}


bool spatialSpreadingData::current() const
{
    return
    (
        topologyCurrent()
     && pointsInstance_ == mesh_.pointsInstance()
     && pointEventNo_ == currentPointEventNo()
    );
}


label spatialSpreadingData::currentPointEventNo() const
{
    const pointIOField& points =
        static_cast<const pointIOField&>(mesh_.mesh().points());
    return points.eventNo();
}


void spatialSpreadingData::updateTopology() const
{
    adjacency_.setSize(mesh_.nFaces());

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    forAll(owner, edgei)
    {
        adjacency_[owner[edgei]].append(neighbour[edgei]);
        adjacency_[neighbour[edgei]].append(owner[edgei]);
    }

    forAll(adjacency_, facei)
    {
        Foam::sort(adjacency_[facei]);

        label nUnique = 0;
        for (const label nbr : adjacency_[facei])
        {
            if (!nUnique || adjacency_[facei][nUnique - 1] != nbr)
            {
                adjacency_[facei][nUnique++] = nbr;
            }
        }
        adjacency_[facei].setSize(nUnique);
    }

    centreTreePtr_.reset(nullptr);
    boundaryEdges_.clear();
    boundaryEdgesValid_ = false;
    centreBoundsValid_ = false;
    localComponents_.clear();
    componentsValid_ = false;
    facePoints_.clear();
    faceRadii_.clear();
    maxFaceRadius_ = 0;
    facePointsValid_ = false;
    topologyValid_ = true;

    nFaces_ = mesh_.nFaces();
    nPoints_ = mesh_.nPoints();
    pointEventNo_ = currentPointEventNo();
    pointsInstance_ = mesh_.pointsInstance();
    facesInstance_ = mesh_.facesInstance();
}


void spatialSpreadingData::updateBoundaryEdges() const
{
    if (boundaryEdgesValid_)
    {
        return;
    }

    boundaryEdges_.setSize(mesh_.nFaces());
    const List<labelPair>& boundaryConnections =
        mesh_.boundaryConnections();
    const labelUList& boundaryFaces = mesh_.patch().boundaryFaces();
    forAll(boundaryConnections, edgei)
    {
        const label facei = boundaryFaces[edgei];
        if (facei >= 0 && facei < boundaryEdges_.size())
        {
            boundaryEdges_[facei].append(edgei);
        }
    }

    boundaryEdgesValid_ = true;
}


void spatialSpreadingData::updateCentreBounds() const
{
    const vectorField& centres = this->centres();
    centreBoundsValid_ = !centres.empty();
    if (!centreBoundsValid_)
    {
        return;
    }

    centreBoundsMin_ = centres[0];
    centreBoundsMax_ = centres[0];

    for (const point& centre : centres)
    {
        centreBoundsMin_.x() =
            Foam::min(centreBoundsMin_.x(), centre.x());
        centreBoundsMin_.y() =
            Foam::min(centreBoundsMin_.y(), centre.y());
        centreBoundsMin_.z() =
            Foam::min(centreBoundsMin_.z(), centre.z());
        centreBoundsMax_.x() =
            Foam::max(centreBoundsMax_.x(), centre.x());
        centreBoundsMax_.y() =
            Foam::max(centreBoundsMax_.y(), centre.y());
        centreBoundsMax_.z() =
            Foam::max(centreBoundsMax_.z(), centre.z());
    }

    ++centreBoundsBuilds_;
    centreBoundsFaceVisits_ += centres.size();
}


void spatialSpreadingData::updateGeometry() const
{
    centreTreePtr_.reset(nullptr);
    centreBoundsValid_ = false;
    facePoints_.clear();
    faceRadii_.clear();
    maxFaceRadius_ = 0;
    facePointsValid_ = false;
    pointEventNo_ = currentPointEventNo();
    pointsInstance_ = mesh_.pointsInstance();
}


void spatialSpreadingData::updateFacePoints() const
{
    facePoints_.resize(mesh_.nFaces());
    faceRadii_.setSize(mesh_.nFaces(), Zero);
    scalar localMaxFaceRadius = 0;

    const faceList& faces = mesh_.faces();
    const pointField& points = mesh_.points();
    const vectorField& centres = this->centres();

    forAll(faces, facei)
    {
        const face& f = faces[facei];
        pointField& pointsForFace = facePoints_[facei];
        pointsForFace.setSize(f.size());

        scalar radius = 0;
        forAll(f, pointi)
        {
            pointsForFace[pointi] = points[f[pointi]];
            radius =
                Foam::max
                (
                    radius,
                    mag(pointsForFace[pointi] - centres[facei])
                );
        }

        faceRadii_[facei] = radius;
        localMaxFaceRadius = Foam::max(localMaxFaceRadius, radius);
    }

    reduce
    (
        localMaxFaceRadius,
        maxOp<scalar>(),
        UPstream::msgType(),
        mesh_.comm()
    );
    maxFaceRadius_ = localMaxFaceRadius;

    facePointsValid_ = true;
}


void spatialSpreadingData::updateComponents() const
{
    const label nFaces = mesh_.nFaces();
    localComponents_.setSize(nFaces, -1);

    const globalIndex& globalFaces =
        mesh_.mesh().globalData().globalMeshFaceAddr();

    labelList queue(nFaces);
    forAll(localComponents_, seedFace)
    {
        if (localComponents_[seedFace] >= 0)
        {
            continue;
        }

        label head = 0;
        label tail = 0;
        queue[tail++] = seedFace;

        const label rootToken =
            globalFaces.toGlobal(mesh_.faceLabels()[seedFace]);

        localComponents_[seedFace] = rootToken;

        while (head < tail)
        {
            const label facei = queue[head++];

            for (const label nbr : adjacency_[facei])
            {
                if (localComponents_[nbr] < 0)
                {
                    localComponents_[nbr] = rootToken;
                    queue[tail++] = nbr;
                }
            }
        }
    }

    const List<labelPair>& boundaryConnections =
        mesh_.boundaryConnections();
    const labelUList& boundaryFaces = mesh_.patch().boundaryFaces();
    const label comm = mesh_.comm();
    const label nProcs = UPstream::nProcs(comm);
    const label myProc = UPstream::myProcNo(comm);

    bool converged = false;
    while (!converged)
    {
        List<DynamicList<label>> sendFaces(nProcs);
        List<DynamicList<label>> sendTokens(nProcs);
        bool localChanged = false;

        forAll(boundaryConnections, edgei)
        {
            const label localFacei = boundaryFaces[edgei];
            const label nbrProc = boundaryConnections[edgei].first();
            const label nbrMeshFace = boundaryConnections[edgei].second();

            if (nbrProc == myProc)
            {
                const label nbrFacei = mesh_.whichFace(nbrMeshFace);
                if
                (
                    nbrFacei >= 0
                 && localComponents_[localFacei]
                    > localComponents_[nbrFacei]
                )
                {
                    localComponents_[localFacei] =
                        localComponents_[nbrFacei];
                    localChanged = true;
                }
            }
            else if (nbrProc >= 0 && nbrProc < nProcs)
            {
                sendFaces[nbrProc].append(nbrMeshFace);
                sendTokens[nbrProc].append(localComponents_[localFacei]);
            }
        }

        PstreamBuffers pBufs(comm);

        for (label proci = 0; proci < nProcs; ++proci)
        {
            if (proci == myProc)
            {
                continue;
            }

            UOPstream toProc(proci, pBufs);
            toProc
                << labelList(sendFaces[proci])
                << labelList(sendTokens[proci]);
        }

        pBufs.finishedSends();

        for (label proci = 0; proci < nProcs; ++proci)
        {
            if (proci == myProc)
            {
                continue;
            }

            labelList remoteFaces;
            labelList remoteTokens;
            UIPstream fromProc(proci, pBufs);
            fromProc >> remoteFaces >> remoteTokens;

            if (remoteFaces.size() != remoteTokens.size())
            {
                FatalErrorInFunction
                    << "Inconsistent distributed component data from "
                    << proci << exit(FatalError);
            }

            forAll(remoteFaces, i)
            {
                const label localFacei = mesh_.whichFace(remoteFaces[i]);
                if
                (
                    localFacei >= 0
                 && localComponents_[localFacei] > remoteTokens[i]
                )
                {
                    localComponents_[localFacei] = remoteTokens[i];
                    localChanged = true;
                }
            }
        }

        forAll(adjacency_, facei)
        {
            for (const label nbr : adjacency_[facei])
            {
                const label token =
                    Foam::min(localComponents_[facei],
                              localComponents_[nbr]);

                if (localComponents_[facei] != token)
                {
                    localComponents_[facei] = token;
                    localChanged = true;
                }
                if (localComponents_[nbr] != token)
                {
                    localComponents_[nbr] = token;
                    localChanged = true;
                }
            }
        }

        label changed = localChanged;
        reduce(changed, maxOp<label>(), UPstream::msgType(), comm);
        converged = !changed;
    }

    componentsValid_ = true;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void spatialSpreadingData::update(const bool needFacePoints) const
{
    if (!topologyCurrent())
    {
        updateTopology();
    }
    else if (!current())
    {
        updateGeometry();
    }

    if (!centreBoundsValid_)
    {
        updateCentreBounds();
    }

    if (needFacePoints && !facePointsValid_)
    {
        updateFacePoints();
    }

    if (!centreTreePtr_)
    {
        const vectorField& centres = this->centres();

        if (!centres.empty())
        {
            treeBoundBox bounds(centres);
            bounds.inflate(ROOTSMALL);

            centreTreePtr_.reset
            (
                new indexedOctree<treeDataPoint>
                (
                    treeDataPoint(centres),
                    bounds,
                    20,
                    100,
                    100
                )
            );
        }
    }
}


const vectorField& spatialSpreadingData::centres() const
{
    return mesh_.areaCentres().primitiveField();
}


bool spatialSpreadingData::centreBounds(point& min, point& max) const
{
    update(false);
    if (!centreBoundsValid_)
    {
        return false;
    }

    min = centreBoundsMin_;
    max = centreBoundsMax_;
    return true;
}


const vectorField& spatialSpreadingData::normals() const
{
    return mesh_.faceAreaNormals().primitiveField();
}


const scalarField& spatialSpreadingData::areas() const
{
    return mesh_.S().field();
}


const List<pointField>& spatialSpreadingData::facePoints() const
{
    update(true);
    return facePoints_;
}


const scalarField& spatialSpreadingData::faceRadii() const
{
    update(true);
    return faceRadii_;
}


scalar spatialSpreadingData::maxFaceRadius() const
{
    update(true);
    return maxFaceRadius_;
}


const labelListList& spatialSpreadingData::adjacency() const
{
    update(false);
    return adjacency_;
}


const labelList& spatialSpreadingData::components() const
{
    update(false);
    if (!componentsValid_)
    {
        updateComponents();
    }
    return localComponents_;
}


labelList spatialSpreadingData::candidatesInSphere
(
    const point& centre,
    const scalar radius
) const
{
    update(false);

    if (centreTreePtr_)
    {
        labelList candidates
        (
            centreTreePtr_->findSphere(centre, Foam::sqr(radius))
        );
        Foam::sort(candidates);
        return candidates;
    }

    labelList candidates(mesh_.nFaces());
    forAll(candidates, facei)
    {
        candidates[facei] = facei;
    }
    return candidates;
}


labelList spatialSpreadingData::candidatesInTangentPlane
(
    const point& centre,
    const vector& normal,
    const vector& majorAxis,
    const vector& minorAxis,
    const scalar radius
) const
{
    update(false);

    if (!centreTreePtr_ || mesh_.nFaces() == 0)
    {
        return labelList();
    }

    scalar minNormalDistance = GREAT;
    scalar maxNormalDistance = -GREAT;
    // Projection extrema of an AABB occur at its eight corners.
    for (label corneri = 0; corneri < 8; ++corneri)
    {
        const point corner
        (
            (corneri & 1) ? centreBoundsMax_.x() : centreBoundsMin_.x(),
            (corneri & 2) ? centreBoundsMax_.y() : centreBoundsMin_.y(),
            (corneri & 4) ? centreBoundsMax_.z() : centreBoundsMin_.z()
        );
        const scalar normalDistance = (corner - centre) & normal;
        minNormalDistance = Foam::min(minNormalDistance, normalDistance);
        maxNormalDistance = Foam::max(maxNormalDistance, normalDistance);
    }
    tangentBoundCornerChecks_ += 8;

    const scalar queryRadius = Foam::max(radius, ROOTSMALL);
    pointField boundsPoints(8);
    label pointi = 0;

    for (const scalar normalDistance :
        {minNormalDistance, maxNormalDistance})
    {
        for (const scalar majorSign : {-1.0, 1.0})
        {
            for (const scalar minorSign : {-1.0, 1.0})
            {
                boundsPoints[pointi++] =
                    centre
                  + queryRadius*
                    (
                        majorSign*majorAxis
                      + minorSign*minorAxis
                    )
                  + normalDistance*normal;
            }
        }
    }

    labelList candidates
    (
        centreTreePtr_->findBox(treeBoundBox(boundsPoints))
    );
    Foam::sort(candidates);
    return candidates;
}


void spatialSpreadingData::geodesicCandidates
(
    const label sourceProc,
    const label sourceFace,
    const scalar searchRadius,
    DynamicList<label>& candidateFaces,
    DynamicList<scalar>& selectedDistances
) const
{
    labelList sourceProcs(1, sourceProc);
    labelList sourceFaces(1, sourceFace);
    scalarField searchRadii(1, searchRadius);
    List<labelList> candidates(1);
    List<scalarField> distances(1);

    geodesicCandidatesBatch
    (
        sourceProcs,
        sourceFaces,
        searchRadii,
        candidates,
        distances
    );

    candidateFaces.clear();
    selectedDistances.clear();
    candidateFaces.append(candidates[0]);
    selectedDistances.append(distances[0]);
}


void spatialSpreadingData::geodesicCandidatesBatch
(
    const labelUList& sourceProcs,
    const labelUList& sourceFaces,
    const scalarField& searchRadii,
    List<labelList>& candidateFaces,
    List<scalarField>& selectedDistances
) const
{
    geodesicCandidatesBatch
    (
        sourceProcs, sourceFaces, sourceFaces, searchRadii,
        candidateFaces, selectedDistances
    );
}


void spatialSpreadingData::geodesicCandidatesBatch
(
    const labelUList& sourceProcs,
    const labelUList& sourceIds,
    const labelUList& sourceFaces,
    const scalarField& searchRadii,
    List<labelList>& candidateFaces,
    List<scalarField>& selectedDistances
) const
{
    if
    (
        sourceProcs.size() != sourceIds.size()
     || sourceProcs.size() != sourceFaces.size()
     || sourceProcs.size() != searchRadii.size()
    )
    {
        FatalErrorInFunction
            << "Inconsistent batched geodesic source data"
            << exit(FatalError);
    }

    update(false);
    updateBoundaryEdges();

    const label nSources = sourceProcs.size();
    const label comm = mesh_.comm();
    const label nProcs = UPstream::nProcs(comm);
    const label myProc = UPstream::myProcNo(comm);
    const List<labelPair>& boundaryConnections =
        mesh_.boundaryConnections();
    const pointField& haloCentres = mesh_.haloFaceCentres();
    const vectorField& centres = this->centres();

    List<Map<scalar>> distances(nSources);
    List<Map<label>> sourceIndices(nProcs);
    batchQueue queue;

    forAll(sourceProcs, sourcei)
    {
        const label sourceProc = sourceProcs[sourcei];
        const label sourceId = sourceIds[sourcei];
        const label sourceFace = sourceFaces[sourcei];

        if (sourceProc < 0 || sourceProc >= nProcs)
        {
            FatalErrorInFunction
                << "Invalid geodesic source processor " << sourceProc
                << exit(FatalError);
        }

        if (sourceId < 0 || sourceIndices[sourceProc].contains(sourceId))
        {
            FatalErrorInFunction
                << "Duplicate geodesic source "
                << sourceProc << ':' << sourceId
                << exit(FatalError);
        }
        sourceIndices[sourceProc].insert(sourceId, sourcei);

        if (sourceProc == myProc)
        {
            if (sourceFace < 0 || sourceFace >= mesh_.nFaces())
            {
                FatalErrorInFunction
                    << "Source face " << sourceFace
                    << " is invalid on processor " << myProc
                    << exit(FatalError);
            }

            distances[sourcei].insert(sourceFace, 0);
            queue.emplace(0, sourcei, sourceFace);
        }
    }

    bool converged = false;
    while (!converged)
    {
        List<DynamicList<label>> sendSourceProcs(nProcs);
        List<DynamicList<label>> sendSourceFaces(nProcs);
        List<DynamicList<label>> sendFaces(nProcs);
        List<DynamicList<scalar>> sendDistances(nProcs);
        bool localChanged = false;

        while (!queue.empty())
        {
            const batchQueueEntry entry = queue.top();
            queue.pop();

            const scalar distance = std::get<0>(entry);
            const label sourcei = std::get<1>(entry);
            const label facei = std::get<2>(entry);
            Map<scalar>& sourceDistances = distances[sourcei];

            if
            (
                distance > sourceDistances.at(facei)
             || distance > searchRadii[sourcei]
            )
            {
                continue;
            }

            for (const label nbr : adjacency_[facei])
            {
                const scalar nextDistance =
                    distance + mag(centres[nbr] - centres[facei]);

                if (nextDistance > searchRadii[sourcei])
                {
                    continue;
                }

                if
                (
                    !sourceDistances.contains(nbr)
                 || nextDistance < sourceDistances.at(nbr)
                )
                {
                    if (sourceDistances.contains(nbr))
                    {
                        sourceDistances.at(nbr) = nextDistance;
                    }
                    else
                    {
                        sourceDistances.insert(nbr, nextDistance);
                    }
                    queue.emplace(nextDistance, sourcei, nbr);
                    localChanged = true;
                }
            }

            for (const label edgei : boundaryEdges_[facei])
            {
                const label nbrProc = boundaryConnections[edgei].first();
                const label nbrMeshFace =
                    boundaryConnections[edgei].second();

                if (nbrProc == myProc)
                {
                    const label nbrFacei = mesh_.whichFace(nbrMeshFace);
                    if (nbrFacei < 0)
                    {
                        continue;
                    }

                    const scalar nextDistance =
                        distance
                      + mag(centres[nbrFacei] - centres[facei]);

                    if
                    (
                        nextDistance <= searchRadii[sourcei]
                     &&
                        (
                            !sourceDistances.contains(nbrFacei)
                         || nextDistance < sourceDistances.at(nbrFacei)
                        )
                    )
                    {
                        if (sourceDistances.contains(nbrFacei))
                        {
                            sourceDistances.at(nbrFacei) = nextDistance;
                        }
                        else
                        {
                            sourceDistances.insert(nbrFacei, nextDistance);
                        }
                        queue.emplace(nextDistance, sourcei, nbrFacei);
                        localChanged = true;
                    }
                }
                else if (nbrProc >= 0 && nbrProc < nProcs)
                {
                    const scalar nextDistance =
                        distance
                      + mag(haloCentres[edgei] - centres[facei]);

                    if (nextDistance <= searchRadii[sourcei])
                    {
                        sendSourceProcs[nbrProc].append(sourceProcs[sourcei]);
                        sendSourceFaces[nbrProc].append(sourceIds[sourcei]);
                        sendFaces[nbrProc].append(nbrMeshFace);
                        sendDistances[nbrProc].append(nextDistance);
                        localChanged = true;
                    }
                }
            }
        }

        // Exchange all sources' frontier updates in one message phase.
        PstreamBuffers pBufs(comm);
        for (label proci = 0; proci < nProcs; ++proci)
        {
            if (proci == myProc || sendFaces[proci].empty())
            {
                continue;
            }

            UOPstream toProc(proci, pBufs);
            toProc
                << labelList(sendSourceProcs[proci])
                << labelList(sendSourceFaces[proci])
                << labelList(sendFaces[proci])
                << scalarField(sendDistances[proci]);
        }
        pBufs.finishedSends();
        ++geodesicExchangePhases_;

        for (label proci = 0; proci < nProcs; ++proci)
        {
            if (proci == myProc || !pBufs.recvDataCount(proci))
            {
                continue;
            }

            labelList remoteSourceProcs;
            labelList remoteSourceFaces;
            labelList remoteFaces;
            scalarField remoteDistances;
            UIPstream fromProc(proci, pBufs);
            fromProc
                >> remoteSourceProcs
                >> remoteSourceFaces
                >> remoteFaces
                >> remoteDistances;

            if
            (
                remoteSourceProcs.size() != remoteSourceFaces.size()
             || remoteSourceProcs.size() != remoteFaces.size()
             || remoteSourceProcs.size() != remoteDistances.size()
            )
            {
                FatalErrorInFunction
                    << "Inconsistent batched geodesic data from "
                    << proci << exit(FatalError);
            }

            forAll(remoteFaces, i)
            {
                const label sourceProc = remoteSourceProcs[i];
                const label sourceId = remoteSourceFaces[i];
                if
                (
                    sourceProc < 0 || sourceProc >= nProcs
                 || !sourceIndices[sourceProc].contains(sourceId)
                )
                {
                    FatalErrorInFunction
                        << "Missing geodesic source descriptor "
                        << sourceProc << ':' << sourceId
                        << " on processor " << myProc
                        << exit(FatalError);
                }

                const label sourcei = sourceIndices[sourceProc].at(sourceId);
                const scalar nextDistance = remoteDistances[i];
                const label localFacei = mesh_.whichFace(remoteFaces[i]);

                if
                (
                    localFacei >= 0
                 && nextDistance <= searchRadii[sourcei]
                 &&
                    (
                        !distances[sourcei].contains(localFacei)
                     || nextDistance
                        < distances[sourcei].at(localFacei)
                    )
                )
                {
                    if (distances[sourcei].contains(localFacei))
                    {
                        distances[sourcei].at(localFacei) = nextDistance;
                    }
                    else
                    {
                        distances[sourcei].insert(localFacei, nextDistance);
                    }
                    queue.emplace(nextDistance, sourcei, localFacei);
                    localChanged = true;
                }
            }
        }

        label changed = localChanged || !queue.empty();
        reduce(changed, maxOp<label>(), UPstream::msgType(), comm);
        converged = !changed;
    }

    candidateFaces.setSize(nSources);
    selectedDistances.setSize(nSources);
    forAll(distances, sourcei)
    {
        const labelList faceLabels = distances[sourcei].sortedToc();
        DynamicList<label> sourceCandidates(faceLabels.size());
        DynamicList<scalar> sourceDistances(faceLabels.size());

        for (const label facei : faceLabels)
        {
            const scalar distance = distances[sourcei].at(facei);
            if (distance <= searchRadii[sourcei])
            {
                sourceCandidates.append(facei);
                sourceDistances.append(distance);
            }
        }

        candidateFaces[sourcei] = sourceCandidates;
        selectedDistances[sourcei] = sourceDistances;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
