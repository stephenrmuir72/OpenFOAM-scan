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

#include "physicalSpreadingGeometry.H"
#include "HashSet.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

physicalSpreadingGeometry::physicalSpreadingGeometry(const faMesh& mesh)
:
    mesh_(mesh),
    nFaces_(-1),
    nPoints_(-1),
    pointsInstance_(),
    facesInstance_(),
    hasFacePoints_(false),
    globalFacesPtr_(nullptr),
    centreTreePtr_(nullptr),
    centres_(),
    normals_(),
    areas_(),
    facePoints_(),
    adjacency_(),
    components_(),
    faceRadii_(),
    maxFaceRadius_(0)
{}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

void physicalSpreadingGeometry::update(const bool needFacePoints) const
{
    if
    (
        globalFacesPtr_
     && !mesh_.moving()
     && nFaces_ == mesh_.nFaces()
     && nPoints_ == mesh_.nPoints()
     && pointsInstance_ == mesh_.pointsInstance()
     && facesInstance_ == mesh_.facesInstance()
     && (hasFacePoints_ || !needFacePoints)
    )
    {
        return;
    }

    centreTreePtr_.reset(nullptr);

    const label comm = mesh_.comm();
    const label nProcs = UPstream::nProcs(comm);
    const label myProc = UPstream::myProcNo(comm);

    globalFacesPtr_.reset(new globalIndex(mesh_.nFaces(), comm));
    const globalIndex& globalFaces = *globalFacesPtr_;
    const label nGlobalFaces = globalFaces.totalSize();

    const vectorField localCentres(mesh_.areaCentres().primitiveField());
    const vectorField localNormals(mesh_.faceAreaNormals().primitiveField());
    const scalarField localAreas(mesh_.S().field());
    const labelList localPolyFaces(mesh_.faceLabels());

    List<pointField> localFacePoints;
    if (needFacePoints)
    {
        localFacePoints.resize(mesh_.nFaces());

        const faceList& localFaces = mesh_.faces();
        const pointField& localPoints = mesh_.points();

        forAll(localFaces, facei)
        {
            const face& f = localFaces[facei];
            pointField& facePoints = localFacePoints[facei];
            facePoints.setSize(f.size());

            forAll(f, pointi)
            {
                facePoints[pointi] = localPoints[f[pointi]];
            }
        }
    }

    List<vectorField> procCentres(nProcs);
    List<vectorField> procNormals(nProcs);
    List<scalarField> procAreas(nProcs);
    List<List<pointField>> procFacePoints
    (
        needFacePoints ? nProcs : 0
    );
    List<labelList> procPolyFaces(nProcs);

    procCentres[myProc] = localCentres;
    procNormals[myProc] = localNormals;
    procAreas[myProc] = localAreas;
    if (needFacePoints)
    {
        procFacePoints[myProc] = localFacePoints;
    }
    procPolyFaces[myProc] = localPolyFaces;

    Pstream::allGatherList(procCentres, UPstream::msgType(), comm);
    Pstream::allGatherList(procNormals, UPstream::msgType(), comm);
    Pstream::allGatherList(procAreas, UPstream::msgType(), comm);
    if (needFacePoints)
    {
        Pstream::allGatherList
        (
            procFacePoints,
            UPstream::msgType(),
            comm
        );
    }
    Pstream::allGatherList(procPolyFaces, UPstream::msgType(), comm);

    centres_.setSize(nGlobalFaces);
    normals_.setSize(nGlobalFaces);
    areas_.setSize(nGlobalFaces);
    facePoints_.setSize(needFacePoints ? nGlobalFaces : 0);
    faceRadii_.setSize(nGlobalFaces, Zero);
    maxFaceRadius_ = 0;

    forAll(procCentres, proci)
    {
        const label start = globalFaces.localStart(proci);

        if
        (
            procCentres[proci].size() != procNormals[proci].size()
         || procCentres[proci].size() != procAreas[proci].size()
         || procCentres[proci].size() != procPolyFaces[proci].size()
         || (
                needFacePoints
             && procCentres[proci].size() != procFacePoints[proci].size()
            )
        )
        {
            FatalErrorInFunction
                << "Inconsistent gathered finite-area geometry on processor "
                << proci << exit(FatalError);
        }

        forAll(procCentres[proci], facei)
        {
            const label globalFacei = start + facei;
            centres_[globalFacei] = procCentres[proci][facei];

            const vector normal(procNormals[proci][facei]);
            const scalar magNormal = mag(normal);
            if (magNormal <= ROOTVSMALL)
            {
                FatalErrorInFunction
                    << "Zero normal for global finite-area face "
                    << globalFacei << exit(FatalError);
            }
            normals_[globalFacei] = normal/magNormal;

            areas_[globalFacei] = procAreas[proci][facei];

            if (needFacePoints)
            {
                facePoints_[globalFacei] = procFacePoints[proci][facei];

                scalar faceRadius = 0;
                for (const point& p : facePoints_[globalFacei])
                {
                    faceRadius =
                        Foam::max(faceRadius, mag(p - centres_[globalFacei]));
                }

                faceRadii_[globalFacei] = faceRadius;
                maxFaceRadius_ = Foam::max(maxFaceRadius_, faceRadius);
            }
        }
    }

    List<labelHashSet> localAdjacency(mesh_.nFaces());
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    forAll(owner, edgei)
    {
        const label own = owner[edgei];
        const label nei = neighbour[edgei];
        localAdjacency[own].insert(globalFaces.toGlobal(myProc, nei));
        localAdjacency[nei].insert(globalFaces.toGlobal(myProc, own));
    }

    const List<labelPair>& boundaryConnections = mesh_.boundaryConnections();
    const labelUList& boundaryFaces = mesh_.patch().boundaryFaces();

    forAll(boundaryConnections, edgei)
    {
        const label localFacei = boundaryFaces[edgei];
        const label nbrProc = boundaryConnections[edgei].first();
        const label nbrMeshFace = boundaryConnections[edgei].second();

        if (nbrProc < 0 || nbrProc >= nProcs)
        {
            continue;
        }

        const label nbrLocalFace = procPolyFaces[nbrProc].find(nbrMeshFace);
        if (nbrLocalFace < 0)
        {
            continue;
        }

        const label nbrGlobalFace = globalFaces.toGlobal(nbrProc, nbrLocalFace);
        const label ownGlobalFace = globalFaces.toGlobal(myProc, localFacei);

        if (nbrGlobalFace != ownGlobalFace)
        {
            localAdjacency[localFacei].insert(nbrGlobalFace);
        }
    }

    labelListList localAdjacencyLists(mesh_.nFaces());
    forAll(localAdjacency, facei)
    {
        localAdjacencyLists[facei] = localAdjacency[facei].sortedToc();
    }

    List<labelListList> procAdjacency(nProcs);
    procAdjacency[myProc] = localAdjacencyLists;
    Pstream::allGatherList(procAdjacency, UPstream::msgType(), comm);

    List<labelHashSet> adjacencySets(nGlobalFaces);
    forAll(procAdjacency, proci)
    {
        const label start = globalFaces.localStart(proci);

        forAll(procAdjacency[proci], facei)
        {
            const label globalFacei = start + facei;

            for (const label nbr : procAdjacency[proci][facei])
            {
                if (nbr >= 0 && nbr < nGlobalFaces && nbr != globalFacei)
                {
                    adjacencySets[globalFacei].insert(nbr);
                    adjacencySets[nbr].insert(globalFacei);
                }
            }
        }
    }

    adjacency_.setSize(nGlobalFaces);
    forAll(adjacency_, facei)
    {
        adjacency_[facei] = adjacencySets[facei].sortedToc();
    }

    components_.setSize(nGlobalFaces, -1);
    labelList queue(nGlobalFaces);
    label componenti = 0;

    forAll(components_, seedFace)
    {
        if (components_[seedFace] >= 0)
        {
            continue;
        }

        label head = 0;
        label tail = 0;
        queue[tail++] = seedFace;
        components_[seedFace] = componenti;

        while (head < tail)
        {
            const label facei = queue[head++];

            for (const label nbr : adjacency_[facei])
            {
                if (components_[nbr] < 0)
                {
                    components_[nbr] = componenti;
                    queue[tail++] = nbr;
                }
            }
        }

        ++componenti;
    }

    if (!centres_.empty())
    {
        treeBoundBox overallBb(centres_);
        overallBb.inflate(Foam::max(ROOTSMALL, maxFaceRadius_));

        centreTreePtr_.reset
        (
            new indexedOctree<treeDataPoint>
            (
                treeDataPoint(centres_),
                overallBb,
                20,
                100,
                100
            )
        );
    }

    nFaces_ = mesh_.nFaces();
    nPoints_ = mesh_.nPoints();
    pointsInstance_ = mesh_.pointsInstance();
    facesInstance_ = mesh_.facesInstance();
    hasFacePoints_ = needFacePoints;
}


labelList physicalSpreadingGeometry::candidatesInSphere
(
    const point& centre,
    const scalar radius
) const
{
    if (centreTreePtr_)
    {
        labelList candidates
        (
            centreTreePtr_->findSphere(centre, Foam::sqr(radius))
        );
        Foam::sort(candidates);
        return candidates;
    }

    labelList candidates(centres_.size());
    forAll(candidates, facei)
    {
        candidates[facei] = facei;
    }

    return candidates;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
