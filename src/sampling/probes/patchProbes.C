/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "patchProbes.H"
#include "volFields.H"
#include "IOmanip.H"
#include "treeBoundBox.H"
#include "treeDataFace.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchProbes, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        patchProbes,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchProbes::findElements(const fvMesh& mesh)
{
    const int communicator = UPstream::worldComm;
    const auto myProci = UPstream::myProcNo(communicator);
    const auto numProc = UPstream::nProcs(communicator);

    // Test for valid proc number
    const auto valid_proc = [numProc](int proci) noexcept
    {
        return (proci >= 0 && proci < numProc);
    };

    (void)mesh.tetBasePtIs();

    const polyBoundaryMesh& bm = mesh.boundaryMesh();

    // The extracted sampling information (output):
    // - cellIds_   : cells, not used
    // - faceIds_   : mesh face id
    // - patchIds_  : patch corresponding to faceIds_
    // - procIds_   : processor
    //
    // NOTE: faceIds_ and patchIds_ are filled with local
    // information and re-filtered later

    elementList_.resize_nocopy(this->size());
    elementList_ = -1;
    faceList_.resize_nocopy(this->size());
    faceList_ = -1;
    patchIDList_.resize_nocopy(this->size());
    patchIDList_ = -1;
    processor_.resize_nocopy(this->size());
    processor_ = -1;

    // Nearest selection key: (sqr(distance), rank)
    // Construct as missed sample
    // == (large distance, proc = max_ranks)
    List<Tuple2<scalar, int>> nearest
    (
        this->size(),
        { VGREAT, UPstream::max_ranks() }
    );

    // Initialize all hit points as point::max
    pointField hitPoints(this->size(), point::max);

    const labelList patchIDs(bm.patchSet(patchNames_).sortedToc());

    label nFaces = 0;
    forAll(patchIDs, i)
    {
        nFaces += bm[patchIDs[i]].size();
    }

    if (nFaces > 0)
    {
        // Collect mesh faces and bounding box
        labelList bndFaces(nFaces);
        treeBoundBox overallBb;

        nFaces = 0;
        forAll(patchIDs, i)
        {
            const polyPatch& pp = bm[patchIDs[i]];
            forAll(pp, i)
            {
                bndFaces[nFaces++] = pp.start()+i;
                const face& f = pp[i];

                // Without reduction.
                overallBb.add(pp.points(), f);
            }
        }

        Random rndGen(123456);
        overallBb.inflate(rndGen, 1e-4, ROOTVSMALL);


        const indexedOctree<treeDataFace> boundaryTree
        (
            // patch faces only
            treeDataFace(mesh, std::move(bndFaces)),

            overallBb,      // overall search domain
            8,              // maxLevel
            10,             // leafsize
            3.0             // duplicity
        );

        const auto& treeData = boundaryTree.shapes();

        const auto& samples = probeLocations();

        forAll(samples, probei)
        {
            const point& sample = samples[probei];

            pointIndexHit pHit = boundaryTree.findNearest
            (
                sample,
                Foam::sqr(boundaryTree.bb().mag())
            );

            if (!pHit.hit())
            {
                pHit = boundaryTree.findNearest(sample, Foam::sqr(GREAT));
            }

            const point& p = pHit.point();

            if (pHit.hit())
            {
                label facei = treeData.objectIndex(pHit.index());
                label patchi = bm.whichPatch(facei);

                if (isA<emptyPolyPatch>(bm[patchi]))
                {
                    WarningInFunction
                        << " The sample point: " << sample
                        << " belongs to " << patchi
                        << " which is an empty patch. This is not permitted. "
                        << " This sample will not be included "
                        << endl;
                }
                else
                {
                    faceList_[probei] = facei;
                    patchIDList_[probei] = patchi;
                    hitPoints[probei] = p;
                    nearest[probei] = { sample.distSqr(p), myProci };
                }
            }
        }
    }

    // Find nearest - globally consistent
    {
        UPstream::minLocation(nearest, communicator);
    }

    // Update hitPoints, with extra sanity
    forAll(nearest, i)
    {
        const auto proci = nearest[i].second();
        const point& sample = probeLocations()[i];

        // No hit on my proc
        // - set as point::max (for later reduction)
        //   or preserve original sample point if no hit anywhere

        if (proci != myProci)
        {
            hitPoints[i] = valid_proc(proci) ? point::max : sample;
        }
    }

    // Selects the best point, the original sample or point::max
    #if (OPENFOAM >= 2506)
    Foam::reduce
    (
        hitPoints.data(),
        hitPoints.size(),
        minOp<point>(),
        UPstream::msgType(),
        UPstream::worldComm
    );
    #else
    Foam::reduce
    (
        reinterpret_cast<point::cmptType*>(hitPoints.data()),
        (point::nComponents * hitPoints.size()),
        minOp<point::cmptType>(),
        UPstream::msgType(),
        UPstream::worldComm
    );
    #endif

    if (debug)
    {
        InfoInFunction << nl;

        const auto& samples = probeLocations();

        forAll(nearest, samplei)
        {
            const auto proci = nearest[samplei].second();
            const auto facei = faceList_[samplei];
            const point& sample = samples[samplei];

            Info<< "    " << samplei << " coord:" << sample
                << " found on processor:"
                << (proci < UPstream::nProcs() ? proci : -1)
                << " in local face:" << facei
                << " with location:" << hitPoints[samplei]
                << endl;
        }
    }

    // Update actual probe locations and store old ones
    {
        pointField& locations = this->probeLocations();

        oldPoints_ = locations;
        locations = std::move(hitPoints);
    }


    // Re-filter the local faceIds_/patchIds_
    // to only retain the globally nearest sample point.

    forAll(nearest, samplei)
    {
        const auto proci = nearest[samplei].second();

        // Valid (local) face to sample?
        if (proci == myProci && faceList_[samplei] >= 0)
        {
            processor_[samplei] = proci;
        }
        else
        {
            faceList_[samplei] = -1;
            patchIDList_[samplei] = -1;
            processor_[samplei] = -1;
        }
    }

    // NB: never reduce (cellIds_, faceIds_) : they are processor-local

    Foam::reduce
    (
        processor_.data(),
        processor_.size(),
        maxOp<label>(),
        UPstream::msgType(),
        communicator
    );

    Foam::reduce
    (
        patchIDList_.data(),
        patchIDList_.size(),
        maxOp<label>(),
        UPstream::msgType(),
        communicator
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchProbes::patchProbes
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const bool loadFromFiles,
    const bool readFields
)
:
    probes(name, runTime, dict, loadFromFiles, false)
{
    if (readFields)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchProbes::performAction(unsigned request)
{
    if (!pointField::empty() && request && prepare(request))
    {
        performAction(scalarFields_, request);
        performAction(vectorFields_, request);
        performAction(sphericalTensorFields_, request);
        performAction(symmTensorFields_, request);
        performAction(tensorFields_, request);

        performAction(surfaceScalarFields_, request);
        performAction(surfaceVectorFields_, request);
        performAction(surfaceSphericalTensorFields_, request);
        performAction(surfaceSymmTensorFields_, request);
        performAction(surfaceTensorFields_, request);
    }
    return true;
}


bool Foam::patchProbes::execute()
{
    if (onExecute_)
    {
        return performAction(ACTION_ALL & ~ACTION_WRITE);
    }

    return true;
}


bool Foam::patchProbes::write()
{
    return performAction(ACTION_ALL);
}


bool Foam::patchProbes::read(const dictionary& dict)
{
    if (!dict.readIfPresent("patches", patchNames_))
    {
        patchNames_.resize(1);
        patchNames_.first() = dict.get<word>("patch");
    }

    return probes::read(dict);
}


// ************************************************************************* //
