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
    localComponents_(),
    centreTreePtr_(nullptr),
    centreBounds_(),
    facePoints_(),
    faceRadii_(),
    maxFaceRadius_(0),
    centreSearchBuilt_(false),
    componentsValid_(false),
    facePointsValid_(false)
{}


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

void spatialSpreadingData::buildCentreSearch() const
{
    const vectorField& centres = this->centres();
    centreBounds_ = treeBoundBox(centres);

    if (centreBounds_.empty())
    {
        return;
    }

    treeBoundBox bounds(centreBounds_);
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


void spatialSpreadingData::buildFacePoints() const
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


void spatialSpreadingData::buildComponents() const
{
    const label nFaces = mesh_.nFaces();
    localComponents_.setSize(nFaces, -1);

    // Preserve the local component graph defined by the finite-area LDU pairs.
    const labelListList& faceEdges = mesh_.patch().faceEdges();
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

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

            for (const label edgei : faceEdges[facei])
            {
                // Only internal edges have local owner/neighbour entries.
                if (edgei >= owner.size())
                {
                    continue;
                }

                label nbr;
                if (owner[edgei] == facei)
                {
                    nbr = neighbour[edgei];
                }
                else if (neighbour[edgei] == facei)
                {
                    nbr = owner[edgei];
                }
                else
                {
                    // Other faces on a multiply-connected edge are not
                    // connected by the finite-area LDU addressing.
                    continue;
                }

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

        forAll(owner, edgei)
        {
            const label ownerFacei = owner[edgei];
            const label neighbourFacei = neighbour[edgei];
            const label token =
                Foam::min
                (
                    localComponents_[ownerFacei],
                    localComponents_[neighbourFacei]
                );

            if (localComponents_[ownerFacei] != token)
            {
                localComponents_[ownerFacei] = token;
                localChanged = true;
            }
            if (localComponents_[neighbourFacei] != token)
            {
                localComponents_[neighbourFacei] = token;
                localChanged = true;
            }
        }

        label changed = localChanged;
        reduce(changed, maxOp<label>(), UPstream::msgType(), comm);
        converged = !changed;
    }

    componentsValid_ = true;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void spatialSpreadingData::prepare(const bool needFacePoints) const
{
    if (!centreSearchBuilt_)
    {
        buildCentreSearch();
        centreSearchBuilt_ = true;
    }

    if (needFacePoints && !facePointsValid_)
    {
        buildFacePoints();
    }

    if (!componentsValid_)
    {
        buildComponents();
    }
}


const vectorField& spatialSpreadingData::centres() const
{
    return mesh_.areaCentres().primitiveField();
}


const treeBoundBox& spatialSpreadingData::centreBounds() const
{
    if (!centreSearchBuilt_)
    {
        FatalErrorInFunction
            << "Centre search data has not been prepared"
            << exit(FatalError);
    }

    return centreBounds_;
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
    if (!facePointsValid_)
    {
        FatalErrorInFunction
            << "Face-point data was not requested during preparation"
            << exit(FatalError);
    }

    return facePoints_;
}


const scalarField& spatialSpreadingData::faceRadii() const
{
    if (!facePointsValid_)
    {
        FatalErrorInFunction
            << "Face-radius data was not requested during preparation"
            << exit(FatalError);
    }

    return faceRadii_;
}


scalar spatialSpreadingData::maxFaceRadius() const
{
    if (!facePointsValid_)
    {
        FatalErrorInFunction
            << "Face-radius data was not requested during preparation"
            << exit(FatalError);
    }

    return maxFaceRadius_;
}


const labelList& spatialSpreadingData::components() const
{
    if (!componentsValid_)
    {
        FatalErrorInFunction
            << "Component data has not been prepared"
            << exit(FatalError);
    }

    return localComponents_;
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
    if (!centreSearchBuilt_)
    {
        FatalErrorInFunction
            << "Centre search data has not been prepared"
            << exit(FatalError);
    }

    if (!centreTreePtr_ || mesh_.nFaces() == 0)
    {
        return labelList();
    }

    scalar minNormalDistance = GREAT;
    scalar maxNormalDistance = -GREAT;
    // Projection extrema of an AABB occur at its eight corners.
    for (label corneri = 0; corneri < 8; ++corneri)
    {
        const point corner(centreBounds_.octCorner(corneri));
        const scalar normalDistance = (corner - centre).inner(normal);
        minNormalDistance = Foam::min(minNormalDistance, normalDistance);
        maxNormalDistance = Foam::max(maxNormalDistance, normalDistance);
    }

    const scalar queryRadius = Foam::max(radius, ROOTSMALL);
    treeBoundBox bounds;

    for (const scalar normalDistance :
        {minNormalDistance, maxNormalDistance})
    {
        for (const scalar majorSign : {-1.0, 1.0})
        {
            for (const scalar minorSign : {-1.0, 1.0})
            {
                bounds.add
                (
                    centre
                  + queryRadius*
                    (
                        majorSign*majorAxis
                      + minorSign*minorAxis
                    )
                  + normalDistance*normal
                );
            }
        }
    }

    labelList candidates
    (
        centreTreePtr_->findBox(bounds)
    );
    Foam::sort(candidates);
    return candidates;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
