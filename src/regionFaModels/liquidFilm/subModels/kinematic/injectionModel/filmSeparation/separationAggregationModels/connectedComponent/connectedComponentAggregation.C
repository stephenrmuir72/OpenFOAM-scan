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

#include "connectedComponentAggregation.H"
#include "addToRunTimeSelectionTable.H"
#include "faMesh.H"
#include "globalIndex.H"
#include "Map.H"
#include "DynamicList.H"
#include "Pstream.H"
#include "PstreamBuffers.H"
#include "PstreamReduceOps.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * Static Data Members and Selectors * * * * * * * * * * //

defineTypeNameAndDebug(connectedComponentAggregation, 0);
addToRunTimeSelectionTable
(
    separationAggregationModel,
    connectedComponentAggregation,
    dictionary
);


const Enum<connectedComponentAggregation::diameterModelType>
connectedComponentAggregation::diameterModelTypeNames
({
    {
        diameterModelType::MASS_WEIGHTED_FILM_THICKNESS,
        "massWeightedFilmThickness"
    },
    { diameterModelType::FIXED_DIAMETER, "fixedDiameter" },
});


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

static void aggregateConnectedFaces
(
    const labelList& globalFaces,
    const scalarField& masses,
    const vectorField& momenta,
    const scalarField& diameterMoments,
    const scalarField& times,
    const scalarField& eventTimes,
    const scalarField& coherenceLengths,
    const labelField& morphologies,
    const pointField& faceCentres,
    const labelList& edgeStarts,
    const labelList& edgeEnds,
    const connectedComponentAggregation::diameterModelType diameterModel,
    const scalar fixedDiameter,
    DynamicList<label>& injectGlobalFaces,
    DynamicList<scalar>& injectMasses,
    DynamicList<scalar>& injectDiameters,
    DynamicList<vector>& injectVelocities,
    DynamicList<label>& resetGlobalFaces
)
{
    injectGlobalFaces.clear();
    injectMasses.clear();
    injectDiameters.clear();
    injectVelocities.clear();
    resetGlobalFaces.clear();

    const label nActive = globalFaces.size();

    if
    (
        masses.size() != nActive
     || momenta.size() != nActive
     || diameterMoments.size() != nActive
     || times.size() != nActive
     || eventTimes.size() != nActive
     || coherenceLengths.size() != nActive
     || morphologies.size() != nActive
     || faceCentres.size() != nActive
     || edgeStarts.size() != edgeEnds.size()
    )
    {
        FatalErrorInFunction
            << "Inconsistent connected-face aggregation input sizes"
            << abort(FatalError);
    }

    if (!nActive)
    {
        return;
    }

    Map<label> globalToNode(2*nActive);

    forAll(globalFaces, nodei)
    {
        globalToNode.insert(globalFaces[nodei], nodei);
    }

    List<DynamicList<label>> neighbours(nActive);

    forAll(edgeStarts, edgei)
    {
        const label nodeA = globalToNode.lookup(edgeStarts[edgei], -1);
        const label nodeB = globalToNode.lookup(edgeEnds[edgei], -1);

        if
        (
            nodeA != -1
         && nodeB != -1
         && nodeA != nodeB
         && morphologies[nodeA] == morphologies[nodeB]
         && mag(faceCentres[nodeA] - faceCentres[nodeB])
            <= min(coherenceLengths[nodeA], coherenceLengths[nodeB])
        )
        {
            if (!neighbours[nodeA].found(nodeB))
            {
                neighbours[nodeA].append(nodeB);
            }

            if (!neighbours[nodeB].found(nodeA))
            {
                neighbours[nodeB].append(nodeA);
            }
        }
    }

    auto nodeCentreLess = [&](const label nodeA, const label nodeB)
    {
        const point& centreA = faceCentres[nodeA];
        const point& centreB = faceCentres[nodeB];

        if (centreA.x() < centreB.x()) return true;
        if (centreB.x() < centreA.x()) return false;
        if (centreA.y() < centreB.y()) return true;
        if (centreB.y() < centreA.y()) return false;
        if (centreA.z() < centreB.z()) return true;
        if (centreB.z() < centreA.z()) return false;

        return globalFaces[nodeA] < globalFaces[nodeB];
    };

    for (auto& nei : neighbours)
    {
        Foam::sort(nei, nodeCentreLess);
    }

    labelList component(nActive, -1);
    label nComponents = 0;
    DynamicList<label> front;

    forAll(globalFaces, startNodei)
    {
        if (component[startNodei] == -1)
        {
            front.clear();
            front.append(startNodei);
            component[startNodei] = nComponents;

            for (label fronti = 0; fronti < front.size(); ++fronti)
            {
                const label nodei = front[fronti];
                const DynamicList<label>& nodeNbrs = neighbours[nodei];

                forAll(nodeNbrs, nbri)
                {
                    const label nbrNodei = nodeNbrs[nbri];

                    if (component[nbrNodei] == -1)
                    {
                        component[nbrNodei] = nComponents;
                        front.append(nbrNodei);
                    }
                }
            }

            ++nComponents;
        }
    }

    List<DynamicList<label>> components(nComponents);

    forAll(component, nodei)
    {
        components[component[nodei]].append(nodei);
    }

    forAll(components, compi)
    {
        DynamicList<label>& nodes = components[compi];
        Foam::sort(nodes, nodeCentreLess);

        scalar mass = 0;
        scalar ageMoment = 0;
        scalar eventRateMoment = 0;

        forAll(nodes, i)
        {
            const label nodei = nodes[i];

            if (eventTimes[nodei] <= VSMALL)
            {
                FatalErrorInFunction
                    << "Connected component face " << globalFaces[nodei]
                    << " has a non-positive physical event time "
                    << eventTimes[nodei] << abort(FatalError);
            }

            mass += masses[nodei];
            ageMoment += masses[nodei]*times[nodei];
            eventRateMoment += masses[nodei]/eventTimes[nodei];
        }

        if (mass <= VSMALL || eventRateMoment <= VSMALL)
        {
            continue;
        }

        const scalar meanAge = ageMoment/mass;
        const scalar meanEventTime = mass/eventRateMoment;

        // Compare times, not an unnormalised mass-time moment.
        if (meanAge < meanEventTime)
        {
            continue;
        }

        scalar diameterMoment = 0;
        vector momentum(Zero);
        point weightedCentre(Zero);

        forAll(nodes, i)
        {
            const label nodei = nodes[i];
            diameterMoment += diameterMoments[nodei];
            momentum += momenta[nodei];
            weightedCentre += masses[nodei]*faceCentres[nodei];
        }

        weightedCentre /= mass;

        label representativeNodei = nodes[0];
        scalar bestDistance =
            magSqr(faceCentres[representativeNodei] - weightedCentre);

        for (const label nodei : nodes)
        {
            const scalar d2 = magSqr(faceCentres[nodei] - weightedCentre);

            if
            (
                d2 < bestDistance
             || (
                    d2 == bestDistance
                 && nodeCentreLess(nodei, representativeNodei)
                )
            )
            {
                representativeNodei = nodei;
                bestDistance = d2;
            }
        }

        scalar diameter = fixedDiameter;

        if
        (
            diameterModel
         == connectedComponentAggregation::MASS_WEIGHTED_FILM_THICKNESS
        )
        {
            diameter = diameterMoment/mass;
        }

        for (const label groupi : nodes)
        {
            resetGlobalFaces.append(globalFaces[groupi]);
        }

        injectGlobalFaces.append(globalFaces[representativeNodei]);
        injectMasses.append(mass);
        injectDiameters.append(diameter);
        injectVelocities.append(momentum/mass);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

connectedComponentAggregation::connectedComponentAggregation
(
    const liquidFilmBase& film,
    const dictionary& dict
)
:
    separationAggregationModel(type(), film, dict),
    diameterModel_
    (
        diameterModelTypeNames.getOrDefault
        (
            "diameterModel",
            dict,
            MASS_WEIGHTED_FILM_THICKNESS
        )
    ),
    fixedDiameter_(dict.getOrDefault<scalar>("fixedDiameter", 0))
{
    // boundaryConnections() processor ids and the global face numbering used by
    // aggregate() are defined on the world communicator, so the finite-area
    // mesh must use it too (otherwise toGlobal(proc, face) is mis-offset).
    if (mesh_.comm() != UPstream::worldComm)
    {
        FatalErrorInFunction
            << "connectedComponentAggregation requires the finite-area mesh to "
            << "use the world communicator." << nl
            << exit(FatalError);
    }

    if (diameterModel_ == FIXED_DIAMETER && fixedDiameter_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "The fixedDiameter entry must be positive when "
            << "diameterModel is fixedDiameter" << nl
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void connectedComponentAggregation::aggregate
(
    scalarField& massAccumulator,
    vectorField& momentumAccumulator,
    scalarField& diameterAccumulator,
    scalarField& accumulationTime,
    const scalarField& eventTime,
    const scalarField& coherenceLength,
    const labelField& morphology,
    scalarField& massToInject,
    scalarField& diameterToInject,
    vectorField& velocityToInject
) const
{
    const faMesh& mesh = mesh_;
    const label nFaces = mesh.nFaces();
    const label comm = mesh.comm();
    const diameterModelType diameterModel = diameterModel_;
    const scalar fixedDiameter = fixedDiameter_;

    if
    (
        massAccumulator.size() != nFaces
     || momentumAccumulator.size() != nFaces
     || diameterAccumulator.size() != nFaces
     || accumulationTime.size() != nFaces
     || eventTime.size() != nFaces
     || coherenceLength.size() != nFaces
     || morphology.size() != nFaces
     || massToInject.size() != nFaces
     || diameterToInject.size() != nFaces
     || velocityToInject.size() != nFaces
    )
    {
        FatalErrorInFunction
            << "All aggregation fields must have finite-area face size "
            << nFaces << abort(FatalError);
    }

    massToInject = 0;
    diameterToInject = 0;

    boolList activeFace(nFaces, false);
    label nLocalActive = 0;

    forAll(massAccumulator, facei)
    {
        if (massAccumulator[facei] > VSMALL)
        {
            if (coherenceLength[facei] <= VSMALL)
            {
                FatalErrorInFunction
                    << "Active connected-component face " << facei
                    << " has a non-positive physical coherence reach "
                    << coherenceLength[facei] << abort(FatalError);
            }

            activeFace[facei] = true;
            ++nLocalActive;
        }
    }

    if (!returnReduceOr(nLocalActive > 0, comm))
    {
        return;
    }

    // Global face numbering on worldComm so that boundaryConnections() proc
    // ids (also worldComm-based) index the same offsets; mesh.comm() is
    // asserted equal to worldComm at construction.
    const globalIndex globalFaceNumber
    (
        mesh.mesh().nFaces(),
        UPstream::worldComm
    );
    const labelList& meshFaceLabels = mesh.faceLabels();
    const pointField& centres = mesh.areaCentres().primitiveField();

    Map<label> globalToLocal(2*nFaces);
    forAll(meshFaceLabels, facei)
    {
        globalToLocal.insert
        (
            globalFaceNumber.toGlobal(meshFaceLabels[facei]),
            facei
        );
    }

    DynamicList<label> localGlobalFaces(nLocalActive);
    DynamicList<scalar> localMasses(nLocalActive);
    DynamicList<vector> localMomenta(nLocalActive);
    DynamicList<scalar> localDiameters(nLocalActive);
    DynamicList<scalar> localTimes(nLocalActive);
    DynamicList<scalar> localEventTimes(nLocalActive);
    DynamicList<scalar> localCoherenceLengths(nLocalActive);
    DynamicList<label> localMorphologies(nLocalActive);
    DynamicList<point> localCentres(nLocalActive);
    DynamicList<label> localEdgeStarts;
    DynamicList<label> localEdgeEnds;

    forAll(activeFace, facei)
    {
        if (activeFace[facei])
        {
            localGlobalFaces.append
            (
                globalFaceNumber.toGlobal(meshFaceLabels[facei])
            );
            localMasses.append(massAccumulator[facei]);
            localMomenta.append(momentumAccumulator[facei]);
            localDiameters.append(diameterAccumulator[facei]);
            localTimes.append(accumulationTime[facei]);
            localEventTimes.append(eventTime[facei]);
            localCoherenceLengths.append(coherenceLength[facei]);
            localMorphologies.append(morphology[facei]);
            localCentres.append(centres[facei]);
        }
    }

    const labelListList& faceFaces = mesh.patch().faceFaces();

    forAll(activeFace, facei)
    {
        if (activeFace[facei])
        {
            const label globalFacei =
                globalFaceNumber.toGlobal(meshFaceLabels[facei]);

            const labelList& nbrFaces = faceFaces[facei];
            for (const label nbrFacei : nbrFaces)
            {
                if (activeFace[nbrFacei])
                {
                    const label globalNbrFacei =
                        globalFaceNumber.toGlobal(meshFaceLabels[nbrFacei]);

                    if (globalFacei < globalNbrFacei)
                    {
                        localEdgeStarts.append(globalFacei);
                        localEdgeEnds.append(globalNbrFacei);
                    }
                }
            }
        }
    }

    const labelList& boundaryFaces = mesh.patch().boundaryFaces();
    const List<labelPair>& boundaryConnections = mesh.boundaryConnections();

    forAll(boundaryConnections, edgei)
    {
        const label facei = boundaryFaces[edgei];

        if (activeFace[facei])
        {
            const label globalFacei =
                globalFaceNumber.toGlobal(meshFaceLabels[facei]);

            const labelPair& connection = boundaryConnections[edgei];

            // Neighbour face indices are polyMesh-local, not faMesh-local.
            const label globalNbrFacei =
                globalFaceNumber.toGlobal
                (
                    connection.first(),
                    connection.second()
                );

            localEdgeStarts.append(globalFacei);
            localEdgeEnds.append(globalNbrFacei);
        }
    }

    labelList localGlobalFaceList(std::move(localGlobalFaces));
    scalarField localMassList(std::move(localMasses));
    vectorField localMomentumList(std::move(localMomenta));
    scalarField localDiameterList(std::move(localDiameters));
    scalarField localTimeList(std::move(localTimes));
    scalarField localEventTimeList(std::move(localEventTimes));
    scalarField localCoherenceLengthList(std::move(localCoherenceLengths));
    labelList localMorphologyList(std::move(localMorphologies));
    pointField localCentreList(std::move(localCentres));
    labelList localEdgeStartList(std::move(localEdgeStarts));
    labelList localEdgeEndList(std::move(localEdgeEnds));

    const label nProcs = UPstream::nProcs(comm);
    const label myRank = UPstream::myProcNo(comm);

    List<labelList> allGlobalFaces(nProcs);
    List<scalarField> allMasses(nProcs);
    List<vectorField> allMomenta(nProcs);
    List<scalarField> allDiameters(nProcs);
    List<scalarField> allTimes(nProcs);
    List<scalarField> allEventTimes(nProcs);
    List<scalarField> allCoherenceLengths(nProcs);
    List<labelList> allMorphologies(nProcs);
    List<pointField> allCentres(nProcs);
    List<labelList> allEdgeStarts(nProcs);
    List<labelList> allEdgeEnds(nProcs);

    allGlobalFaces[myRank] = std::move(localGlobalFaceList);
    allMasses[myRank] = std::move(localMassList);
    allMomenta[myRank] = std::move(localMomentumList);
    allDiameters[myRank] = std::move(localDiameterList);
    allTimes[myRank] = std::move(localTimeList);
    allEventTimes[myRank] = std::move(localEventTimeList);
    allCoherenceLengths[myRank] = std::move(localCoherenceLengthList);
    allMorphologies[myRank] = std::move(localMorphologyList);
    allCentres[myRank] = std::move(localCentreList);
    allEdgeStarts[myRank] = std::move(localEdgeStartList);
    allEdgeEnds[myRank] = std::move(localEdgeEndList);

    if (Pstream::parRun())
    {
        Pstream::gatherList(allGlobalFaces, UPstream::msgType(), comm);
        Pstream::gatherList(allMasses, UPstream::msgType(), comm);
        Pstream::gatherList(allMomenta, UPstream::msgType(), comm);
        Pstream::gatherList(allDiameters, UPstream::msgType(), comm);
        Pstream::gatherList(allTimes, UPstream::msgType(), comm);
        Pstream::gatherList(allEventTimes, UPstream::msgType(), comm);
        Pstream::gatherList(allCoherenceLengths, UPstream::msgType(), comm);
        Pstream::gatherList(allMorphologies, UPstream::msgType(), comm);
        Pstream::gatherList(allCentres, UPstream::msgType(), comm);
        Pstream::gatherList(allEdgeStarts, UPstream::msgType(), comm);
        Pstream::gatherList(allEdgeEnds, UPstream::msgType(), comm);
    }

    DynamicList<label> injectGlobalFaces;
    DynamicList<scalar> injectMasses;
    DynamicList<scalar> injectDiameters;
    DynamicList<vector> injectVelocities;
    DynamicList<label> resetGlobalFaces;

    if (Pstream::master(comm))
    {
        label nActive = 0;
        label nEdges = 0;

        forAll(allGlobalFaces, proci)
        {
            nActive += allGlobalFaces[proci].size();
            nEdges += allEdgeStarts[proci].size();
        }

        labelList globalFaces(nActive);
        scalarField masses(nActive);
        vectorField momenta(nActive);
        scalarField diameters(nActive);
        scalarField times(nActive);
        scalarField eventTimes(nActive);
        scalarField coherenceLengths(nActive);
        labelField morphologies(nActive);
        pointField faceCentres(nActive);
        labelList edgeStarts(nEdges);
        labelList edgeEnds(nEdges);

        label nodei = 0;
        forAll(allGlobalFaces, proci)
        {
            const labelList& procFaces = allGlobalFaces[proci];
            forAll(procFaces, prociFacei)
            {
                globalFaces[nodei] = procFaces[prociFacei];
                masses[nodei] = allMasses[proci][prociFacei];
                momenta[nodei] = allMomenta[proci][prociFacei];
                diameters[nodei] = allDiameters[proci][prociFacei];
                times[nodei] = allTimes[proci][prociFacei];
                eventTimes[nodei] = allEventTimes[proci][prociFacei];
                coherenceLengths[nodei] =
                    allCoherenceLengths[proci][prociFacei];
                morphologies[nodei] =
                    allMorphologies[proci][prociFacei];
                faceCentres[nodei] = allCentres[proci][prociFacei];

                ++nodei;
            }
        }

        label edgei = 0;
        forAll(allEdgeStarts, proci)
        {
            const labelList& starts = allEdgeStarts[proci];
            const labelList& ends = allEdgeEnds[proci];

            forAll(starts, prociEdgei)
            {
                edgeStarts[edgei] = starts[prociEdgei];
                edgeEnds[edgei] = ends[prociEdgei];
                ++edgei;
            }
        }

        aggregateConnectedFaces
        (
            globalFaces,
            masses,
            momenta,
            diameters,
            times,
            eventTimes,
            coherenceLengths,
            morphologies,
            faceCentres,
            edgeStarts,
            edgeEnds,
            diameterModel,
            fixedDiameter,
            injectGlobalFaces,
            injectMasses,
            injectDiameters,
            injectVelocities,
            resetGlobalFaces
        );
    }

    labelList injectGlobalFaceList(std::move(injectGlobalFaces));
    scalarField injectMassList(std::move(injectMasses));
    scalarField injectDiameterList(std::move(injectDiameters));
    vectorField injectVelocityList(std::move(injectVelocities));
    labelList resetGlobalFaceList(std::move(resetGlobalFaces));

    if (Pstream::parRun())
    {
        Pstream::broadcast(injectGlobalFaceList, comm);
        Pstream::broadcast(injectMassList, comm);
        Pstream::broadcast(injectDiameterList, comm);
        Pstream::broadcast(injectVelocityList, comm);
        Pstream::broadcast(resetGlobalFaceList, comm);
    }

    forAll(resetGlobalFaceList, i)
    {
        const label facei = globalToLocal.lookup(resetGlobalFaceList[i], -1);

        if (facei != -1)
        {
            massAccumulator[facei] = 0;
            momentumAccumulator[facei] = Zero;
            diameterAccumulator[facei] = 0;
            accumulationTime[facei] = 0;
        }
    }

    forAll(injectGlobalFaceList, i)
    {
        const label facei = globalToLocal.lookup(injectGlobalFaceList[i], -1);

        if (facei != -1)
        {
            massToInject[facei] = injectMassList[i];
            diameterToInject[facei] = injectDiameterList[i];
            velocityToInject[facei] = injectVelocityList[i];
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
