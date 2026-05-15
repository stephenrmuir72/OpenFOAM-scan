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

#include "physicalSpreading.H"
#include "faceCoverage.H"
#include "footprintDistance.H"
#include "footprintKernel.H"
#include "footprintShape.H"
#include "DynamicList.H"
#include "IPstream.H"
#include "OPstream.H"
#include "Pstream.H"
#include "PstreamBuffers.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

struct sourceContribution
{
    Foam::label sourceGlobalFace;
    Foam::label targetGlobalFace;
    Foam::scalar mass;
    Foam::vector momentum;
    Foam::scalar normalImpulse;
};


bool sourceContributionLess
(
    const sourceContribution& a,
    const sourceContribution& b,
    const Foam::vectorField& centres
)
{
    if (a.targetGlobalFace != b.targetGlobalFace)
    {
        return a.targetGlobalFace < b.targetGlobalFace;
    }

    const Foam::vector& aCentre = centres[a.sourceGlobalFace];
    const Foam::vector& bCentre = centres[b.sourceGlobalFace];

    if (aCentre.x() != bCentre.x())
    {
        return aCentre.x() < bCentre.x();
    }
    if (aCentre.y() != bCentre.y())
    {
        return aCentre.y() < bCentre.y();
    }
    if (aCentre.z() != bCentre.z())
    {
        return aCentre.z() < bCentre.z();
    }

    return a.sourceGlobalFace < b.sourceGlobalFace;
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

defineTypeNameAndDebug(physicalSpreading, 0);
defineRunTimeSelectionTable(physicalSpreading, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

physicalSpreading::physicalSpreading
(
    const word& modelType,
    const faMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    geometry_(mesh),
    shapePtr_(footprintShape::New(dict)),
    distancePtr_(footprintDistance::New(dict)),
    kernelPtr_(footprintKernel::New(dict)),
    coveragePtr_(faceCoverage::New(dict))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

physicalSpreading::~physicalSpreading()
{}  // Various classes were forward declared


// * * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

void physicalSpreading::spreadSources
(
    areaScalarField& mass,
    areaVectorField& momentum,
    areaScalarField& normalImpulse,
    const areaScalarField& diameter,
    const scalar beta
) const
{
    const scalarField& localMass = mass.primitiveField();
    const vectorField& localMomentum = momentum.primitiveField();
    const scalarField& localNormal = normalImpulse.primitiveField();
    const scalarField& localDiameter = diameter.primitiveField();

    bool localActive = false;
    forAll(localMass, facei)
    {
        localActive =
            localActive
         || mag(localMass[facei]) > VSMALL
         || mag(localMomentum[facei]) > VSMALL
         || mag(localNormal[facei]) > VSMALL;
    }

    if (!returnReduceOr(localActive))
    {
        return;
    }

    const bool needFacePoints = coveragePtr_->needsFacePoints();
    geometry_.update(needFacePoints);

    const globalIndex& globalFaces = geometry_.globalFaces();
    const label comm = mesh_.comm();
    const label nProcs = UPstream::nProcs(comm);
    const label myProc = UPstream::myProcNo(comm);
    const label localStart = globalFaces.localStart(myProc);

    const vectorField& globalCentres = geometry_.centres();
    const labelList& globalComponents = geometry_.components();
    const scalar maxFaceRadius = geometry_.maxFaceRadius();

    DynamicList<label> localSourceIds;
    DynamicList<label> localTargetIds;
    DynamicList<scalar> localMassContributions;
    DynamicList<vector> localMomentumContributions;
    DynamicList<scalar> localNormalContributions;

    forAll(localMass, sourceFacei)
    {
        if
        (
            mag(localMass[sourceFacei]) <= VSMALL
         && mag(localMomentum[sourceFacei]) <= VSMALL
         && mag(localNormal[sourceFacei]) <= VSMALL
        )
        {
            continue;
        }

        const label sourceGlobalFace = localStart + sourceFacei;
        const scalar minorRadius = 0.5*beta*localDiameter[sourceFacei];

        auto appendContribution =
        [
            &localSourceIds,
            &localTargetIds,
            &localMassContributions,
            &localMomentumContributions,
            &localNormalContributions,
            &localMass,
            &localMomentum,
            &localNormal,
            sourceFacei,
            sourceGlobalFace
        ]
        (
            const label targetGlobalFace,
            const scalar fraction
        )
        {
            localSourceIds.append(sourceGlobalFace);
            localTargetIds.append(targetGlobalFace);
            localMassContributions.append
            (
                fraction*localMass[sourceFacei]
            );
            localMomentumContributions.append
            (
                fraction*localMomentum[sourceFacei]
            );
            localNormalContributions.append
            (
                fraction*localNormal[sourceFacei]
            );
        };

        if (!(minorRadius > SMALL))
        {
            appendContribution(sourceGlobalFace, 1);
            continue;
        }

        const footprintShape::footprintGeometry fp =
            shapePtr_->footprintAt(geometry_, sourceGlobalFace, minorRadius);

        const scalar enclosingRadius = fp.enclosingRadius();
        const scalar searchRadius =
            enclosingRadius + (needFacePoints ? maxFaceRadius : scalar(0));

        DynamicList<label> candidateFaces;
        DynamicList<scalar> selectedDistances;
        distancePtr_->findCandidates
        (
            geometry_,
            *shapePtr_,
            fp,
            sourceGlobalFace,
            searchRadius,
            candidateFaces,
            selectedDistances
        );

        DynamicList<label> targetFaces;
        DynamicList<scalar> rawWeights;
        scalar weightSum = 0;

        const label sourceComponent = globalComponents[sourceGlobalFace];

        forAll(candidateFaces, candidatei)
        {
            const label targetGlobalFace = candidateFaces[candidatei];

            if (globalComponents[targetGlobalFace] != sourceComponent)
            {
                continue;
            }

            scalar x = 0;
            scalar y = 0;
            shapePtr_->planarOffset
            (
                fp,
                globalCentres[targetGlobalFace],
                x,
                y
            );
            const scalar planarDistance =
                Foam::sqrt(Foam::sqr(x) + Foam::sqr(y));
            const scalar boundaryRadius =
                shapePtr_->boundaryRadius(x, y, fp.majorRadius, fp.minorRadius);

            const scalar selectedDistance = selectedDistances[candidatei];

            const faceCoverage::coverageResult coverage =
                coveragePtr_->evaluate
                (
                    geometry_,
                    fp,
                    targetGlobalFace,
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
                const vector centroidDelta = coverage.centre - fp.centre;
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
                targetFaces.append(targetGlobalFace);
                rawWeights.append(rawWeight);
                weightSum += rawWeight;
            }
        }

        if (weightSum <= VSMALL)
        {
            appendContribution(sourceGlobalFace, 1);
            continue;
        }

        forAll(targetFaces, targeti)
        {
            appendContribution
            (
                targetFaces[targeti],
                rawWeights[targeti]/weightSum
            );
        }
    }

    List<DynamicList<label>> sendSourceIds(nProcs);
    List<DynamicList<label>> sendTargetIds(nProcs);
    List<DynamicList<scalar>> sendMassContributions(nProcs);
    List<DynamicList<vector>> sendMomentumContributions(nProcs);
    List<DynamicList<scalar>> sendNormalContributions(nProcs);

    forAll(localTargetIds, contributioni)
    {
        const label targetProc =
            globalFaces.whichProcID
            (
                myProc,
                localTargetIds[contributioni]
            );

        sendSourceIds[targetProc].append(localSourceIds[contributioni]);
        sendTargetIds[targetProc].append(localTargetIds[contributioni]);
        sendMassContributions[targetProc].append
        (
            localMassContributions[contributioni]
        );
        sendMomentumContributions[targetProc].append
        (
            localMomentumContributions[contributioni]
        );
        sendNormalContributions[targetProc].append
        (
            localNormalContributions[contributioni]
        );
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
            << labelList(sendSourceIds[proci])
            << labelList(sendTargetIds[proci])
            << scalarField(sendMassContributions[proci])
            << vectorField(sendMomentumContributions[proci])
            << scalarField(sendNormalContributions[proci]);
    }

    pBufs.finishedSends();

    DynamicList<sourceContribution> contributions;

    auto appendContributions = [&contributions]
    (
        const labelList& sourceIds,
        const labelList& targetIds,
        const scalarField& masses,
        const vectorField& momenta,
        const scalarField& normals
    )
    {
        if
        (
            sourceIds.size() != targetIds.size()
         || targetIds.size() != masses.size()
         || targetIds.size() != momenta.size()
         || targetIds.size() != normals.size()
        )
        {
            FatalErrorInFunction
                << "Inconsistent physical-spreading contribution lists"
                << abort(FatalError);
        }

        forAll(targetIds, contributioni)
        {
            sourceContribution contribution;
            contribution.sourceGlobalFace = sourceIds[contributioni];
            contribution.targetGlobalFace = targetIds[contributioni];
            contribution.mass = masses[contributioni];
            contribution.momentum = momenta[contributioni];
            contribution.normalImpulse = normals[contributioni];
            contributions.append(contribution);
        }
    };

    appendContributions
    (
        labelList(sendSourceIds[myProc]),
        labelList(sendTargetIds[myProc]),
        scalarField(sendMassContributions[myProc]),
        vectorField(sendMomentumContributions[myProc]),
        scalarField(sendNormalContributions[myProc])
    );

    for (label proci = 0; proci < nProcs; ++proci)
    {
        if (proci == myProc)
        {
            continue;
        }

        labelList sourceIds;
        labelList targetIds;
        scalarField masses;
        vectorField momenta;
        scalarField normals;

        UIPstream fromProc(proci, pBufs);
        fromProc >> sourceIds >> targetIds >> masses >> momenta >> normals;

        appendContributions
        (
            sourceIds,
            targetIds,
            masses,
            momenta,
            normals
        );
    }

    std::sort
    (
        contributions.begin(),
        contributions.end(),
        [&globalCentres]
        (
            const sourceContribution& a,
            const sourceContribution& b
        )
        {
            return sourceContributionLess(a, b, globalCentres);
        }
    );

    scalarField nextMass(mesh_.nFaces(), Zero);
    vectorField nextMomentum(mesh_.nFaces(), Zero);
    scalarField nextNormal(mesh_.nFaces(), Zero);
    const label localEnd = localStart + mesh_.nFaces();

    for (const sourceContribution& contribution : contributions)
    {
        if
        (
            contribution.targetGlobalFace < localStart
         || contribution.targetGlobalFace >= localEnd
        )
        {
            FatalErrorInFunction
                << "Received a contribution for non-local target face "
                << contribution.targetGlobalFace
                << ", local range is " << localStart << ".."
                << (localEnd - 1)
                << abort(FatalError);
        }

        const label targetFacei =
            contribution.targetGlobalFace - localStart;
        nextMass[targetFacei] += contribution.mass;
        nextMomentum[targetFacei] += contribution.momentum;
        nextNormal[targetFacei] += contribution.normalImpulse;
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
