/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024-2025 OpenCFD Ltd.
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

Description
    All to do with snapping to the surface with buffer layer

\*----------------------------------------------------------------------------*/

#include "snappySnapDriver.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "fvMesh.H"
#include "Time.H"
#include "OBJstream.H"
#include "mapPolyMesh.H"
#include "snapParameters.H"
#include "unitConversion.H"
#include "PatchTools.H"
#include "profiling.H"
#include "addPatchCellLayer.H"
#include "snappyLayerDriver.H"
#include "weightedPosition.H"

#include "localPointRegion.H"
#include "pointConstraints.H"
#include "displacementMotionSolver.H"
#include "meshPointPatch.H"
#include "processorPointPatch.H"
#include "dummyTransform.H"
#include "faceSet.H"
#include "motionSmoother.H"
// #include "tetDecomposer.H"
// #include "tetMatcher.H"

//- Feature splitting
#include "refinementFeatures.H"
#include "refinementSurfaces.H"
#include "searchableSurfaces.H"
#include "processorMeshes.H"
#include "surfaceFeatures.H"
#include "treeDataEdge.H"
#include "treeDataPoint.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::snappySnapDriver::avg
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& pp,
    const bitSet& isPatchMasterEdge,
    const pointField& localPoints
)
{
    const labelListList& pointEdges = pp.pointEdges();
    const labelList& meshPoints = pp.meshPoints();
    const edgeList& edges = pp.edges();

    Field<weightedPosition> wps
    (
        pointEdges.size(),
        pTraits<weightedPosition>::zero
    );

    // weightedPosition::syncPoints(mesh, meshPoints, wps);
    forAll(pointEdges, verti)
    {
        weightedPosition& wp = wps[verti];

        for (const label edgei : pointEdges[verti])
        {
            if (isPatchMasterEdge[edgei])
            {
                const label otherVerti = edges[edgei].otherVertex(verti);

                wp.first() += 1.0;
                wp.second() += localPoints[otherVerti];
            }
        }
    }

    weightedPosition::syncPoints(mesh, meshPoints, wps);

    tmp<pointField> tavg(new pointField(wps.size()));
    pointField& avg = tavg.ref();

    forAll(wps, verti)
    {
        const weightedPosition& wp = wps[verti];

        if (mag(wp.first()) < VSMALL)
        {
            // Set to zero?
            avg[verti] = Zero;
        }
        else
        {
            avg[verti] = wp.second()/wp.first();
        }
    }

    return tavg;
}


Foam::tmp<Foam::pointField>
Foam::snappySnapDriver::smoothLambdaMuPatchDisplacement
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& pp,
    const pointField& ppLocalPoints,
    const List<labelPair>& baffles
)
{
    const labelList meshEdges
    (
        pp.meshEdges(mesh.edges(), mesh.pointEdges())
    );
    const bitSet isPatchMasterEdge
    (
        meshRefinement::getMasterEdges
        (
            mesh,
            meshEdges
        )
    );

    pointField newLocalPoints(ppLocalPoints);

    const label iters = 90;
    const scalar lambda = 0.33;
    const scalar mu = 0.34;

    for (label iter = 0; iter < iters; iter++)
    {
        // Lambda
        newLocalPoints =
            (1 - lambda)*newLocalPoints
          + lambda*avg(mesh, pp, isPatchMasterEdge, newLocalPoints);

        // Mu
        newLocalPoints =
            (1 + mu)*newLocalPoints
          - mu*avg(mesh, pp, isPatchMasterEdge, newLocalPoints);
    }
    return newLocalPoints-ppLocalPoints;
}


// Surface feature walking
// ~~~~~~~~~~~~~~~~~~~~~~~

void Foam::snappySnapDriver::pointPoints
(
    const indirectPrimitivePatch& pp,
    const bitSet& isCandidateFace,
    const List<pointConstraint>& patchConstraints,
    const label pointi,
    DynamicList<label>& pointPoints,
    DynamicList<label>& pointEdges,      // pp edge or -1 for diagonal
    DynamicList<label>& pointFaces
) const
{
    pointPoints.clear();
    pointEdges.clear();
    pointFaces.clear();

    // Note: tried adding only points that are not yet feat-edge attracted
    // (patchConstraints[f[fp]].first() <= 1)
    // but that gives worse snapping.

    // Points on existing edges
    for (const label edgei : pp.pointEdges()[pointi])
    {
        // Do not walk off candidate faces
        bool isCandidate = false;
        for (const label facei : pp.edgeFaces()[edgei])
        {
            if (isCandidateFace[facei])
            {
                isCandidate = true;
                break;
            }
        }

        if (!isCandidate)
        {
            continue;
        }

        const label otherPointi = pp.edges()[edgei].otherVertex(pointi);
        if (!pointPoints.found(otherPointi))
        {
            pointPoints.append(otherPointi);
            pointEdges.append(edgei);
            pointFaces.append(-1);
        }
    }

    // Diagonals
    const auto& pFaces = pp.pointFaces()[pointi];
    for (const label facei : pFaces)
    {
        // Do not walk off candidate faces
        if (!isCandidateFace[facei])
        {
            continue;
        }

        const auto& f = pp.localFaces()[facei];
        const label fp = f.find(pointi);
        const label nextPointi = f[f.fcIndex(fp)];
        const label prevPointi = f[f.rcIndex(fp)];

        forAll(f, fp)
        {
            const label otherPointi = f[fp];
            if
            (
                otherPointi != pointi
             && otherPointi != prevPointi
             && otherPointi != nextPointi
            )
            {
                if (!pointPoints.found(otherPointi))
                {
                    pointPoints.append(otherPointi);
                    pointEdges.append(-1);
                    pointFaces.append(facei);
                }
            }
        }
    }
}


Foam::point Foam::snappySnapDriver::avgPoint
(
    const indirectPrimitivePatch& pp,
    const pointField& points,
    const label pointi
) const
{
    const auto& pointEdges = pp.pointEdges();
    const auto& edges = pp.edges();

    point avgPoint(Zero);
    label nPoints = 0;
    for (const label edgei : pointEdges[pointi])
    {
        avgPoint += points[edges[edgei].otherVertex(pointi)];
        nPoints++;
    }
    return avgPoint/nPoints;
}


void Foam::snappySnapDriver::smoothAttraction
(
    const indirectPrimitivePatch& pp,
    const List<pointConstraint>& patchConstraints,
    const label ppPointi,
    const label maxConstraint,
    const point& attractPt,
    pointField& points
) const
{
    const auto& pointEdges = pp.pointEdges();
    const auto& edges = pp.edges();


    // Attract point
    points[ppPointi] = attractPt;

    for (label i = 0; i < 3; i++)
    {
        bitSet isMovedPoint(pp.nPoints());
        isMovedPoint.set(ppPointi);


        // Walk out and average other points
        DynamicList<label> front;
        front.append(ppPointi);

        for (label iter = 0; iter < 4; iter++)
        {
            DynamicList<label> newFront;
            for (const label frontPointi : front)
            {
                for (const label edgei : pointEdges[frontPointi])
                {
                    const edge& e = edges[edgei];
                    const label nextPointi = e.otherVertex(frontPointi);

                    if (isMovedPoint.set(nextPointi))
                    {
                        if
                        (
                            patchConstraints[nextPointi].first()
                         <= maxConstraint
                        )
                        {
                            const point avgPt(avgPoint(pp, points, nextPointi));
                            points[nextPointi] =
                                0.5*points[nextPointi]+0.5*avgPt;
                        }
                        newFront.append(nextPointi);
                    }
                }
            }

            front = newFront;
        }
    }
}


void Foam::snappySnapDriver::mapPatchAddressing
(
    const mapPolyMesh& morphMap,
    // Old patch to mesh addresssing
    const labelUList& oldMeshPoints,
    const labelUList& oldMeshFaces,
    // New patch to mesh addressing
    const labelUList& newMeshPoints,
    const labelUList& newMeshFaces,
    // new-to-old point map
    labelList& ppMap,
    // old-to-new point map
    labelList& ppReverseMap,
    // new-to-old face map
    labelList& ppFaceMap
)
{
    // Construct old mesh-to-patch addressing
    const Map<label> oldMeshPointMap(invertToMap(oldMeshPoints));

    ppMap.resize_nocopy(newMeshPoints.size());
    ppMap = -1;
    ppReverseMap.resize_nocopy(oldMeshPoints.size());
    ppReverseMap = -1;

    forAll(newMeshPoints, i)
    {
        const label oldMeshPointi = morphMap.pointMap()[newMeshPoints[i]];
        const auto mpFnd = oldMeshPointMap.find(oldMeshPointi);
        if (mpFnd)
        {
            ppMap[i] = mpFnd();
            ppReverseMap[mpFnd()] = i;
        }
    }
    // faceMap (new-to-old) for pp faces
    ppFaceMap.resize_nocopy(newMeshFaces.size());
    ppFaceMap = -1;

    const Map<label> oldMeshFaceMap(invertToMap(oldMeshFaces));
    forAll(newMeshFaces, patchFacei)
    {
        const label newFacei = newMeshFaces[patchFacei];
        const label oldMeshFacei = morphMap.faceMap()[newFacei];
        const auto mpFnd = oldMeshFaceMap.find(oldMeshFacei);
        if (mpFnd)
        {
            ppFaceMap[patchFacei] = mpFnd();
        }
    }
}


/*
// Surface routines
// ~~~~~~~~~~~~~~~~

void Foam::snappySnapDriver::determineFeatureEdgeRegion
(
    const surfaceFeatures& surfFeatures,
    labelList& edgeToRegion
) const
{
    const triSurface& surf = surfFeatures.surface();
    const edgeList& edges = surf.edges();
    const labelListList& pointEdges = surf.pointEdges();

    bitSet isFeaturePoint(surf.nPoints());
    isFeaturePoint.set(surfFeatures.featurePoints());

    bitSet isFeatureEdge(surf.nEdges());
    isFeatureEdge.set(surfFeatures.featureEdges());

    // Walk from feature points
    label regioni = 0;

    for (const label pointi : surfFeatures.featurePoints())
    {
        const labelList& pEdges = surf.pointEdges()[pointi];

        for (const label edgei : pEdges)
        {
            if (isFeatureEdge[edgei] && edgeToRegion[edgei] == -1)
            {
                // Start walking from point, edge

                label stringPointi = pointi;
                label stringEdgei = edgei;

                // Walk until reaching already set feature edge (?loop?)
                // or end of string
                while (true)
                {
                    const edge& e = edges[stringEdgei];
                    // Mark edge
                    edgeToRegion[stringEdgei] = regioni;
                    // Find next point
                    label nextPointi = e.otherVertex(stringPointi);

                    if (isFeaturePoint[nextPointi])
                    {
                        // Found end of string
                        //Pout<< "Found (feature) end of string " << regioni
                        //    << endl;
                        break;
                    }

                    // Find next edge that is feature and not stringEdgei
                    stringEdgei = -1;
                    for (const label nextEdgei : pointEdges[nextPointi])
                    {
                        if
                        (
                            isFeatureEdge[nextEdgei]
                         && edgeToRegion[nextEdgei] == -1
                        )
                        {
                            // Found continuation
                            stringPointi = nextPointi;
                            stringEdgei = nextEdgei;
                            break;
                        }
                    }

                    if (stringEdgei == -1)
                    {
                        //Pout<< "Found (non-feature) end of string " << regioni
                        //    << endl;
                        break;
                    }
                }

                regioni++;
            }
        }
    }
}
*/


void Foam::snappySnapDriver::updateMap
(
    const labelUList& oldToNew,
    Map<label>& map
)
{
    Map<label> newMap(map.capacity());
    forAllIters(map, iter)
    {
        if (iter.key() == -1)
        {
            FatalErrorInFunction<< "problem" << exit(FatalError);
        }
        const label newPointi = oldToNew[iter()];
        if (newPointi != -1)
        {
            newMap.insert(iter.key(), newPointi);
        }
    }
    map = std::move(newMap);
}


void Foam::snappySnapDriver::writeConstraints
(
    const fileName& prefix,
    const pointField& patchPoints,
    const pointField& patchToNearest,
    const List<pointConstraint>& patchConstraints
) const
{
    List<autoPtr<OBJstream>> osPtrs(4);
    const fileName dir(meshRefiner_.mesh().time().timePath());
    osPtrs[0] = autoPtr<OBJstream>::New(dir/"preSnapOutside_no-attr.obj");
    osPtrs[1] = autoPtr<OBJstream>::New(dir/"preSnapOutside_surf-attr.obj");
    osPtrs[2] = autoPtr<OBJstream>::New(dir/"preSnapOutside_line-attr.obj");
    osPtrs[3] = autoPtr<OBJstream>::New(dir/"preSnapOutside_point-attr.obj");
    forAll(patchConstraints, pointi)
    {
        auto& os = osPtrs[patchConstraints[pointi].first()]();
        os.write(linePointRef(patchPoints[pointi], patchToNearest[pointi]));
    }

    Pout<< "Writing constraints to OBJ files to time dir:"
        << meshRefiner_.mesh().time().timeName() << endl;
    forAll(osPtrs, i)
    {
        auto& os = osPtrs[i]();
        Pout<< "nConstraints:" << i
            << " nPoints:" << os.nVertices()/2
            << " out of:" << patchConstraints.size()
            << " to file: " << os.name().name()
            << endl;
    }
}


void Foam::snappySnapDriver::writeConstraints
(
    const fileName& prefix,
    const indirectPrimitivePatch& pp,
    const pointField& patchPoints,
    const pointField& patchToNearest,
    const List<pointConstraint>& patchConstraints
) const
{
    const fileName dir(meshRefiner_.mesh().time().timePath());
    OBJstream os(dir/prefix + ".obj");
    const auto& localFaces = pp.localFaces();

    forAll(localFaces, facei)
    {
        const auto& f = localFaces[facei];
        forAll(f, fp)
        {
            const label pointi = f[fp];
            const label constrainti = patchConstraints[pointi].first();
            if (constrainti == 3)
            {
                // Draw line to attraction
                os.write
                (
                    linePointRef
                    (
                        patchPoints[pointi],
                        patchToNearest[pointi]
                    )
                );
            }
            const label nextPointi = f[f.fcIndex(fp)];
            const label nextConstrainti =
                patchConstraints[nextPointi].first();
            if (constrainti >= 2 && nextConstrainti >= 2)
            {
                // Draw along face edge
                os.write
                (
                    linePointRef
                    (
                        patchPoints[pointi],
                        patchPoints[nextPointi]
                    )
                );
            }
        }
    }
    Pout<< "Feature walk with points:" << os.nVertices()/2
        << " written to file: " << os.name().name()
        << endl;
}


/*
void Foam::snappySnapDriver::determineSurfaceFeatures
(
    const triSurfaceMesh& surfMesh,
    const surfaceFeatures& surfFeatures,

    vectorField& surfEdgeNormal,
    labelList& surfFeatPointOrder,
    labelList& surfFeatEdgeRegion,
    autoPtr<indexedOctree<treeDataEdge>>& surfEdgeTreePtr
) const
{
    // Average normal
    surfEdgeNormal.resize_nocopy(surfMesh.nEdges());
    {
        const auto& edgeFaces = surfMesh.edgeFaces();
        const auto& faceNormals = surfMesh.faceNormals();
        forAll(edgeFaces, edgei)
        {
            const auto& eFaces = edgeFaces[edgei];
            if (eFaces.size() == 2)
            {
                surfEdgeNormal[edgei] = faceNormals[eFaces[0]];
                surfEdgeNormal[edgei] += faceNormals[eFaces[1]];
                surfEdgeNormal[edgei].normalise();
            }
            else
            {
                // Use face normal for open edges or non-manifold ones
                surfEdgeNormal[edgei] = faceNormals[eFaces[0]];
            }
        }
    }

    // Get the topological order (= number of outgoing feature edges) of the
    // feature points. E.g. 3 for corner point of cube.
    // TBD: distributedTriSurfaceMesh? Unless we have buffer layers? Or
    //      store nEdges on every local point (like the normal) ...
    surfFeatPointOrder.resize_nocopy(surfMesh.nPoints());
    surfFeatPointOrder = 0;
    {
        const auto& edges = surfMesh.edges();
        for (const label edgei : surfFeatures.featureEdges())
        {
            const edge& e = edges[edgei];
            surfFeatPointOrder[e[0]]++;
            surfFeatPointOrder[e[1]]++;
        }
    }

    // Assign zoning to all the feature edges. Maybe this should move to
    // surfaceFeatures?
    surfFeatEdgeRegion.resize_nocopy(surfMesh.nEdges());
    surfFeatEdgeRegion = -1;
    determineFeatureEdgeRegion(surfFeatures, surfFeatEdgeRegion);

    surfEdgeTreePtr.reset
    (
        new indexedOctree<treeDataEdge>
        (
            // Feature edges only
            treeDataEdge
            (
                surfMesh.edges(),
                surfMesh.localPoints(),
                surfFeatures.featureEdges()
            ),
            treeBoundBox(surfMesh.bounds()),  // bb
            surfMesh.maxTreeDepth(),     // maxLevel
            10,                 // leafsize
            3.0                 // duplicity
        )
    );
}
*/

/*
// Check that patch point is surrounded by faces attracted to all the
// surface regions.
bool Foam::snappySnapDriver::hasAllRegions
(
    const triSurface& s,
    const label sPointi,

    const indirectPrimitivePatch& pp,
    const labelList& faceToRegion,
    const label ppPointi
) const
{
    // Work out regions on patch
    labelHashSet patchRegions;
    {
        for (const label facei : pp.pointFaces()[ppPointi])
        {
            patchRegions.insert(faceToRegion[facei]);
        }
    }

    // Check that all the surface regions are present
    for (const label facei : s.pointFaces()[sPointi])
    {
        if (!patchRegions.found(s[facei].region()))
        {
            return false;
        }
    }
    return true;
}
void Foam::snappySnapDriver::matchFeaturePoints
(
    const indirectPrimitivePatch& pp,
    const pointField& ppLocalPoints,
    const indexedOctree<treeDataPoint>& ppPointTree,
    const labelList& patchFaceToSurfaceRegion,

    const triSurface& surfMesh,
    const surfaceFeatures& surfFeatures,
    const labelList& surfFeatPointOrder,


    // Index of surface feature point (or -1)
    labelList& patchToPoint,
    // Point on nearest surface
    vectorField& patchToNearest,
    // Normal on nearest surface
    vectorField& patchToNormal,
    // Region on nearest surface
    labelList& patchToRegion,
    // Constraints at surface (1 = surface, 2 = featedge, 3 = featpoint)
    List<pointConstraint>& patchConstraints
) const
{
    //const fvMesh& mesh = meshRefiner_.mesh();

    labelHashSet ppPoints;
    DynamicList<label> candidates;
    for (const label surfPti : surfFeatures.featurePoints())
    {
        const point& surfPt = surfMesh.localPoints()[surfPti];

        // Find nearest mesh point
        const pointIndexHit nearInfo = ppPointTree.findNearest(surfPt, 1.0);
        if (nearInfo.hit())
        {
            ppPoints.clear();
            const label ppi = nearInfo.index();
            ppPoints.insert(ppi);
            // Add surround: points on connected edges
            const auto& pEdges = pp.pointEdges()[ppi];
            for (const label edgei : pEdges)
            {
                const edge& e = pp.edges()[edgei];
                ppPoints.insert(e.otherVertex(ppi));
            }
            // Optional: points on connected faces
            const auto& pFaces = pp.pointFaces()[ppi];
            for (const label facei : pFaces)
            {
                const auto& f = pp.localFaces()[facei];
                ppPoints.insert(f);
            }
        }


        //// Find initial nearish mesh points
        //ppPointTree.findSphere(surfPt, 1.0, ppPoints);
        // Pout<< "For surface feat point:" << surfPt
        //     << " degree:" << surfFeatPointOrder[surfPti]
        //     << " found patch points:"
        //     << UIndirectList<point>(ppLocalPoints, ppPoints.sortedToc())
        //     << endl;

        const auto& pFaces = surfMesh.pointFaces()[surfPti];
        // Pout<< "    pFaces:" << pFaces << endl;
        // Pout<< "    pFC:"
        //     << UIndirectList<point>(surfMesh.faceCentres(), pFaces) << endl;

        DynamicList<label> surfRegions;
        for (const label trii : pFaces)
        {
            surfRegions.push_uniq(surfMesh[trii].region());
        }

        // Find points that have correct order and correct surrounding face
        // attraction
        candidates.clear();
        for (const label ppi : ppPoints)
        {
            if (pp.pointEdges()[ppi].size() >= surfFeatPointOrder[surfPti])
            {
                // Check that surrounding face centres snap to the correct
                // surface region. This might be problematic ...
                const bool hasSurfaceRegions = hasAllRegions
                (
                    surfMesh,
                    surfPti,      // vertex on surface

                    pp,
                    patchFaceToSurfaceRegion,   // from pp to surface
                    ppi                         // vertex on pp
                );

                if (hasSurfaceRegions)
                {
                    candidates.append(ppi);
                }
            }
        }

        // Find best fit
        label minPointi = -1;
        scalar minDistSqr = Foam::sqr(GREAT);
        for (const label ppi : candidates)
        {
            const scalar d2 = magSqr(ppLocalPoints[ppi]-surfPt);
            if (d2 < minDistSqr)
            {
                // Pout<< "    Found nearest patch point:"
                //     << ppLocalPoints[ppi]
                //     << " distance:" << Foam::sqrt(d2) << endl;

                minDistSqr = d2;
                minPointi = ppi;
            }
        }

        if (minPointi != -1)
        {
            patchToPoint[minPointi] = surfPti;
            patchToNearest[minPointi] = surfPt;
            patchToNormal[minPointi] = surfMesh.pointNormals()[surfPti];
            patchConstraints[minPointi] = pointConstraint(3, surfPt);
            patchToRegion[minPointi] = -1;
        }
        else
        {
            // Fall-back - use nearest
            // Pout<< "    Fallback : search for nearest any patch point" << endl;

            // Find nearish mesh point
            //nearInfo = ppPointTree.findNearest(surfPt, 1.0);
            if (nearInfo.hit())
            {
                patchToPoint[nearInfo.index()] = surfPti;
                patchToNearest[nearInfo.index()] = surfPt;
                patchConstraints[nearInfo.index()] =
                    pointConstraint(3, nearInfo.hitPoint());
                patchToRegion[nearInfo.index()] = -1;
            }
        }
    }
}
*/


void Foam::snappySnapDriver::matchFeaturePoints
(
    const indirectPrimitivePatch& pp,
    const pointField& ppLocalPoints,
    const indexedOctree<treeDataPoint>& ppPointTree,
    const labelList& patchFaceToSurfaceRegion,

    // Index of surface/feature-edge-mesh
    labelList& patchToFeature,
    // Index of feature point (or -1)
    labelList& patchToPoint,
    // Point on nearest surface
    vectorField& patchToNearest,
    // Normal on nearest surface
    vectorField& patchToNormal,
    // Region on nearest surface
    labelList& patchToRegion,
    // Constraints at surface (1 = surface, 2 = featedge, 3 = featpoint)
    List<pointConstraint>& patchConstraints
) const
{
    const auto& features = meshRefiner_.features();

    labelHashSet ppPoints;
    labelHashSet ppRegions;
    DynamicList<label> candidates;

    forAll(features, feati)
    {
        const extendedEdgeMesh& surfFeatures = features[feati];

        const auto& surfPointEdges = surfFeatures.pointEdges();
        forAll(surfPointEdges, surfPti)
        {
            const auto& surfPEdges = surfPointEdges[surfPti];
            const point& surfPt = surfFeatures.points()[surfPti];

            // Find nearest mesh point. TBD: limit search radius.
            const pointIndexHit nearInfo =
                ppPointTree.findNearest(surfPt, GREAT);

            if (!nearInfo.hit() || patchToPoint[nearInfo.index()] != -1)
            {
                continue;
            }

            // Get candidate patch points
            ppPoints.clear();
            const label ppi = nearInfo.index();
            ppPoints.insert(ppi);
            // Add surround: points on connected edges
            const auto& pEdges = pp.pointEdges()[ppi];
            for (const label edgei : pEdges)
            {
                const edge& e = pp.edges()[edgei];
                ppPoints.insert(e.otherVertex(ppi));
            }
            // Optional: points on connected faces (so also points across diag)
            const auto& pFaces = pp.pointFaces()[ppi];
            for (const label facei : pFaces)
            {
                const auto& f = pp.localFaces()[facei];
                ppPoints.insert(f);
            }

            // Find points that have correct order and correct surrounding face
            // attraction
            candidates.clear();
            for (const label ppi : ppPoints)
            {
                if (pp.pointEdges()[ppi].size() >= surfPEdges.size())
                {
                    // Pout<< "    Candidate patch point:" << ppLocalPoints[ppi]
                    //     << " with degree:" << pp.pointEdges()[ppi].size()
                    //     << endl;
                    // Pout<< "    Connected feature edges:" << surfPEdges << endl;
                    // Check that surrounding face centres snap to the correct
                    // surface region. This might be problematic on
                    // extendedEdgeMesh since info not stored there. TBD.
                    // ppRegions.clear();
                    // for (const label facei : pp.pointFaces()[ppi])
                    // {
                    //     ppRegions.insert(patchFaceToSurfaceRegion[facei]);
                    // }
                    // if (ppRegions.size() >= surfPEdges.size())
                    {
                        candidates.append(ppi);
                    }
                }
            }

            // Find best fit
            label minPointi = -1;
            scalar minDistSqr = Foam::sqr(GREAT);
            for (const label ppi : candidates)
            {
                const scalar d2 = magSqr(ppLocalPoints[ppi]-surfPt);
                if (d2 < minDistSqr)
                {
                    // Pout<< "    Found nearest patch point:"
                    //     << ppLocalPoints[ppi]
                    //     << " distance:" << Foam::sqrt(d2) << endl;

                    minDistSqr = d2;
                    minPointi = ppi;
                }
            }

            if (minPointi != -1)
            {
                //- Hard: what to do if pp point already assigned?
                //  - keep closest?
                //  - if not what is constraint? Slide-on-surface?
                // if (patchToPoint[minPointi] != -1)
                // {
                //     const scalar d2 = magSqr
                //     (
                //         ppLocalPoints[minPointi]-patchToNearest[minPointi]
                //     );
                //     if (d2 < minDistSqr)
                //     {
                //         // Existing is better
                //         continue;
                //     }
                //     FatalErrorInFunction
                //         << "point " << minPointi
                //         << " already assigned to feature point "
                //         << patchToPoint[minPointi]
                //         << exit(FatalError);
                // }
                patchToFeature[minPointi] = feati;
                patchToPoint[minPointi] = surfPti;
                patchToNearest[minPointi] = surfPt;
                //patchToNormal[minPointi] = surfMesh.pointNormals()[surfPti];
                patchConstraints[minPointi] = pointConstraint(3, surfPt);
                patchToRegion[minPointi] = -1;

                // // Attract connected points a bit as well.
                // smoothAttraction
                // (
                //     pp,
                //     patchConstraints,
                //     minPointi,      // attract point
                //     1,              // max constraint to attract
                //     surfPt,         // attract to
                //     ppLocalPoints   // points to update
                // );
            }
            else
            {
                // Fall-back - use previously found nearest mesh point (discard
                // surrounding points). Note: is exactly the same as above
                // - or should we avoid attraction and leave it to feature
                // edge matching later on?
                if (nearInfo.hit())
                {
                    patchToFeature[nearInfo.index()] = feati;
                    patchToPoint[nearInfo.index()] = surfPti;
                    patchToNearest[nearInfo.index()] = surfPt;
                    patchConstraints[nearInfo.index()] =
                        pointConstraint(3, surfPt);
                    patchToRegion[nearInfo.index()] = -1;
                }
            }
        }
    }
}


void Foam::snappySnapDriver::determineCandidateFaces
(
    const scalar featCosAngle,

    const indirectPrimitivePatch& pp,

    // Point on nearest surface
    const pointField& patchToNearest,
    // Normal on nearest surface
    const vectorField& patchToNormal,
    // Region on nearest surface
    const labelList& patchToRegion,
    // Constraints at surface (1 = surface, 2 = featedge, 3 = featpoint)
    const List<pointConstraint>& patchConstraints,

    bitSet& isCandidateFace
) const
{
    const auto& localFaces = pp.localFaces();
    forAll(localFaces, facei)
    {
        const auto& localF = localFaces[facei];

        // Topological : check if surface region differs across face
        forAll(localF, fp)
        {
            if (patchToRegion[localF[fp]] != patchToRegion[localF[0]])
            {
                // Different vertices attracted to different regions
                isCandidateFace.set(facei);
                break;
            }
            if (patchConstraints[localF[fp]].first() >= 2)
            {
                // Vertex already marked as feature edge/point
                isCandidateFace.set(facei);
                break;
            }
        }

        // Geometric : check if surface normal differs across face
        if (!isCandidateFace[facei] && featCosAngle > -1)
        {
            forAll(localF, i)
            {
                const vector& ni = patchToNormal[localF[i]];

                for (label j = i+1; j < localF.size(); j++)
                {
                    const vector& nj = patchToNormal[localF[j]];

                    // Differing surface normal
                    if ((ni & nj) < featCosAngle)
                    {
                        isCandidateFace.set(facei);
                        break;
                    }
                }
                if (isCandidateFace[facei])
                {
                    break;
                }
            }
        }
    }
}


Foam::pointIndexHit Foam::snappySnapDriver::findNearestPointOnEdge
(
    const point& pt,
    vector& edgeNormal
) const
{
    labelList nearFeature;
    List<pointIndexHit> nearInfo;
    vectorField nearNormal;
    meshRefiner_.features().findNearestEdge
    (
        pointField(1, pt),
        scalarField(1, Foam::sqr(GREAT)),
        nearFeature,        // which edgeMesh
        nearInfo,
        nearNormal
    );

    edgeNormal = nearNormal[0];
    return pointIndexHit(nearInfo[0]);
}


void Foam::snappySnapDriver::matchFeatureEdges
(
    // Mesh faces to snap
    const indirectPrimitivePatch& pp,
    // Coordinates of pp points (updated during snapping)
    pointField& localPoints,

    // Whether pp face has different attraction
    const bitSet& isCandidateFace,

    // Starting patch feature point
    const label startPPointi,
    // Starting surface feature
    const label surfFeati,
    // Starting surface feature point
    const label surfFeatPointi,


    // Per patch face: diagonal cut
    Map<edge>& patchFaceToDiagonal,

    // Per patch point:

    // Index of surface edge
    labelList& patchToEdge,
    // Point on nearest surface
    vectorField& patchToNearest,
    // Normal on nearest surface
    vectorField& patchToNormal,
    // Constraints at surface (1 = surface, 2 = featedge, 3 = featpoint)
    List<pointConstraint>& patchConstraints
) const
{
    // Start walking from feature point (constraints == 3) along the pp
    // edges (+ diagonal edges) in the direction of the feature edges
    // emanating from the feature point.

    const auto& features = meshRefiner_.features();

    patchToEdge.setSize(localPoints.size(), -1);
    patchToNearest.setSize(localPoints.size(), vector::uniform(GREAT));
    patchToNormal.setSize(localPoints.size(), vector::uniform(GREAT));
    patchConstraints.resize_nocopy(localPoints.size());

    if (meshRefiner_.features().size() == 0)
    {
        // No feature edges
        WarningInFunction
            << "No feature edges defined. Cannot match feature edges."
            << endl;
        return;
    }


    static label iter = 0;
    autoPtr<OBJstream> osPtr
    (
        (debug & meshRefinement::MESH)
      ? new OBJstream
        (
            meshRefiner_.mesh().time().timePath()
          / "walk_"+Foam::name(iter)+".obj"
        )
      : nullptr
    );


    // Mark points that are on candidateFaces
    bitSet isCandidatePoint(pp.nPoints());
    for (const label facei : isCandidateFace)
    {
        isCandidatePoint.set(pp.localFaces()[facei]);
    }

    const extendedEdgeMesh& surfFeatures = features[surfFeati];
    const auto& surfPoints = surfFeatures.points();
    const auto& surfEdges = surfFeatures.edges();

    for (const label surfEdgei : surfFeatures.pointEdges()[surfFeatPointi])
    {
        // Mark feature edge as having been done. Not used yet?
        // surfFeatEdgeRegionToStart.insert(regioni, startPPointi);

        // Get the direction of the edge
        point surfPt(surfPoints[surfFeatPointi]);
        vector surfEdgeDir
        (
            linePointRef
            (
                surfPt,
                surfPoints[surfEdges[surfEdgei].otherVertex(surfFeatPointi)]
            ).unitVec()
        );

        // Attract & smooth a bit. Snaps startPPointi to surfPt (in localPoints)
        // and smooths a few surrounding (featureEdge- or surface-) points
        smoothAttraction
        (
            pp,
            patchConstraints,
            startPPointi,
            2,  // Do not affect feature points
            surfPt,
            localPoints
        );

        // Current
        label ppPointi = startPPointi;
        // label ppFacei = -1;


        // From
        //  - patch point startPPointi
        //  - surf point surfFeatPointi
        // start walking along surface faces updating
        //  - surface : surfPt,surfEdgeDir
        //  - patch   : ppPointi, ppFacei

        DynamicList<label> pPoints;
        DynamicList<label> pEdges;
        DynamicList<label> pFaces;

        while (true)
        {
            // Walk from ppPointi to nextPpPointi
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            // Best next point on patch
            label nextPpPointi = -1;
            label nextPpEdgei = -1;
            label nextPpFacei = -1;

            // Nearest point on surface edge and its edge direction
            pointIndexHit nearestSurf(vector::uniform(GREAT));
            vector nearestNormal(vector::uniform(GREAT));

            // Get surrounding points (also diagonal across faces)
            // - point
            // - edge  : -1 for diagonal
            // - face  : >=0 for diagonal
            pointPoints
            (
                pp,
                isCandidateFace,
                patchConstraints,
                ppPointi,
                pPoints,
                pEdges,
                pFaces
            );

            // Algo 1 : try any cut that makes least angle with featedge
            if (nextPpPointi == -1)
            {
                scalar maxCos = Foam::cos(degToRad(60.0));

                forAll(pPoints, i)
                {
                    const label otherPointi = pPoints[i];
                    const label otherEdgei = pEdges[i];
                    const label otherFacei = pFaces[i];

                    // Get other point (in new coordinates)
                    const point& otherPt = localPoints[otherPointi];
                    const vector eVec
                    (
                        (otherPt-localPoints[ppPointi]).normalise()
                    );
                    // Check if in better direction
                    if ((eVec & surfEdgeDir) > maxCos)
                    {
                        // TBD: find on this feature only? Or do we like it
                        // being able to switch?
                        vector edgeNormal;
                        const pointIndexHit info = findNearestPointOnEdge
                        (
                            otherPt,
                            edgeNormal
                        );
                        // if (info.hit())
                        // {
                        //     Pout<< "otherPt:" << otherPt << " info2:" << info
                        //         << " edgeNormal:" << edgeNormal << endl;
                        // }

                        if (info.hit())
                        {
                            maxCos = (eVec & surfEdgeDir);
                            nearestSurf = info;
                            nearestNormal = edgeNormal;
                            nextPpPointi = otherPointi;
                            nextPpEdgei = otherEdgei;
                            nextPpFacei = otherFacei;
                        }
                    }
                }
            }

            // Algo 3: given best guess (nextPpPointi) check if desired cut
            //         has attraction already which would cause bad triangle
            {
                if (nextPpFacei != -1)
                {
                    // Diagonal cut
                    for (const label pointi : pp.localFaces()[nextPpFacei])
                    {
                        if (pointi != ppPointi && pointi != nextPpPointi)
                        {
                            if (patchConstraints[pointi].first() >= 2)
                            {
                                // Check if snapping to both would cause
                                // triangle problem.
                                const triPointRef tri
                                (
                                    patchToNearest[ppPointi],
                                    nearestSurf.point(),
                                    patchToNearest[pointi]
                                );
                                if (tri.quality() < 0.01)
                                {
                                    // Break diag attraction
                                    nextPpPointi = -1;
                                    break;
                                }
                            }
                        }
                    }
                }
                else if (nextPpEdgei != -1)
                {
                    // Cut along edge
                    for (const label facei : pp.edgeFaces()[nextPpEdgei])
                    {
                        for (const label pointi : pp.localFaces()[facei])
                        {
                            if (pointi != ppPointi && pointi != nextPpPointi)
                            {
                                if (patchConstraints[pointi].first() >= 2)
                                {
                                    // Check if snapping to both would cause
                                    // triangle problem.
                                    const triPointRef tri
                                    (
                                        patchToNearest[ppPointi],
                                        nearestSurf.point(),
                                        patchToNearest[pointi]
                                    );
                                    if (tri.quality() < 0.01)
                                    {
                                        // Break edge attraction
                                        nextPpPointi = -1;
                                        break;
                                    }
                                }
                            }
                        }
                        if (nextPpPointi == -1)
                        {
                            break;
                        }
                    }
                }
            }

            if (nextPpPointi == -1)
            {
                break;
            }

            // Pout<< "From ppPoint:" << ppPointi << nl
            //     << "    at:" << localPoints[ppPointi] << nl
            //     << "walked face to:" << nextPpPointi << nl
            //     << "    at:" << localPoints[nextPpPointi]
            //     << endl;


            // Mark any diagonal cuts
            if (nextPpFacei != -1)
            {
                const face& f = pp.localFaces()[nextPpFacei];
                patchFaceToDiagonal.insert
                (
                    nextPpFacei,
                    edge(f.find(ppPointi), f.find(nextPpPointi))
                );
            }


            if (patchConstraints[nextPpPointi].first() >= 2)
            {
                // Reached end of walk
                if (osPtr)
                {
                    osPtr->write
                    (
                        linePointRef
                        (
                            pp.localPoints()[ppPointi],
                            pp.localPoints()[nextPpPointi]
                        )
                    );
                }

                break;
            }


            // Index of feature, not edge!
            patchToEdge[nextPpPointi] = nearestSurf.index();
            // Point on nearest surface
            patchToNearest[nextPpPointi] = nearestSurf.point();
            // Normal on nearest surface
            patchToNormal[nextPpPointi] = nearestNormal;


            // Constraints at surface (1 = surface, 2 = featedge, 3 =
            // featpoint). Approximate feature edge direction.
            vector eVec((nearestSurf.point()-surfPt).normalise());
            patchConstraints[nextPpPointi] = pointConstraint(2, eVec);

            // Attract & smooth a bit. Modifies localPoints:
            // - nextPpPointi gets set to exact surface point
            // - a few surrounding points get dragged along
            smoothAttraction
            (
                pp,
                patchConstraints,
                nextPpPointi,                       // pp point to attract
                1,  // Only affect non-feature points
                patchToNearest[nextPpPointi],       // surf point
                localPoints
            );


            if (osPtr)
            {
                osPtr->write
                (
                    linePointRef
                    (
                        pp.localPoints()[ppPointi],
                        pp.localPoints()[nextPpPointi]
                    )
                );
            }
            // Step to next point
            ppPointi = nextPpPointi;
            // ppFacei = nextPpFacei;
            surfEdgeDir = patchConstraints[nextPpPointi].second();
            surfPt = nearestSurf.point();
        }
    }

    if (osPtr)
    {
        Pout<< osPtr->name() << " written " << osPtr->nVertices() << " points"
            << endl;
    }

    if (debug & meshRefinement::MESH)
    {
        OBJstream os
        (
            meshRefiner_.mesh().time().timePath()
          / "featEdges_matching_" + Foam::name(iter) + ".obj"
        );
        os.write(pp.localFaces(), localPoints, false);
        Pout<< "Write pp after featEdge matching iter " << iter
            << " to " << os.name()
            << " nVertices:" << os.nVertices() << endl;
    }

    iter++;
}


void Foam::snappySnapDriver::matchSurfaces
(
    // Mesh faces to snap
    const indirectPrimitivePatch& pp,
    // Coordinates of pp points
    pointField& localPoints,

    //- Per patch point

    // Index of surface 'face' (or -1)
    labelList& patchToFace,
    // Region of surface 'face' (or -1)
    labelList& patchToRegion,
    // Point on nearest surface
    vectorField& patchToNearest,
    // Normal on nearest surface
    vectorField& patchToNormal,
    // Constraints at surface (1 = surface, 2 = featedge, 3 = featpoint)
    List<pointConstraint>& patchConstraints
) const
{
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    List<pointIndexHit> info;
    vectorField normal;
    labelList region;
    labelList surface;
    surfaces.findNearestRegion
    (
        surfaces.surfaces(),
        localPoints,
        scalarField(localPoints.size(), sqr(GREAT)),

        surface,
        info,
        region,
        normal
    );

    patchToFace.setSize(localPoints.size(), -1);
    patchToRegion.setSize(localPoints.size(), -1);
    patchToNearest.setSize(localPoints.size(), vector::uniform(GREAT));
    patchToNormal.setSize(localPoints.size(), vector::uniform(GREAT));
    patchConstraints.resize_nocopy(localPoints.size());

    forAll(patchConstraints, pointi)
    {
        if (patchConstraints[pointi].first() <= 1)  // redo nearest surface
        {
            patchToFace[pointi] = info[pointi].index();
            patchToRegion[pointi] = region[pointi];
            patchToNearest[pointi] = info[pointi].point();
            patchToNormal[pointi] = normal[pointi];
            patchConstraints[pointi] = pointConstraint(1, normal[pointi]);
            localPoints[pointi] = info[pointi].point();
        }
    }
}


void Foam::snappySnapDriver::matchToSurface
(
    const bool doFeatures,
    const scalar featCosAngle,

    const indirectPrimitivePatch& pp,   // addressing only
    pointField& ppLocalPoints,          // points (modified during feat-edges)
    const vectorField& ppPointNormals,  // normal
    // const indexedOctree<treeDataPoint>& ppPointTree,

    // Per string index a starting pp point
    Map<label>& surfFeatEdgeRegionToStart,


    //- Per patch face

    // Region on nearest surface
    labelList& patchFaceToSurfaceRegion,

    // Face to diagonal cut
    Map<edge>& patchFaceToDiagonal,

    //- Per patch point

    // Point on nearest surface
    pointField& patchToNearest,
    // Normal on nearest surface
    vectorField& patchToNormal,
    // Index of nearest surface
    labelList& patchToSurface,
    // Region on nearest surface
    labelList& patchToRegion,
    // Surface 'face' index
    labelList& patchToFace,
    // Surface 'edge' index
    labelList& patchToEdge,
    // Surface 'feature' index
    labelList& patchToFeature,
    // Surface feature point (if patch point is attracted)
    labelList& patchToPoint,
    // Constraints at surface (1 = surface, 2 = featedge, 3 = featpoint)
    List<pointConstraint>& patchConstraints
) const
{
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    if (debug & meshRefinement::MESH)
    {
        OBJstream os(meshRefiner_.mesh().time().timePath()/"pp_smoothed.obj");
        os.write(pp.localFaces(), ppLocalPoints, false);
        Pout<< "Write pp after smoothing to " << os.name()
            << " nVertices:" << os.nVertices() << endl;
    }

    patchToNearest.setSize(pp.nPoints(), vector::uniform(GREAT));
    patchToNormal.setSize(pp.nPoints(), vector::uniform(GREAT));
    patchToSurface.setSize(pp.nPoints(), -1);
    patchToRegion.setSize(pp.nPoints(), -1);
    patchToFeature.setSize(pp.nPoints(), -1);
    patchToPoint.setSize(pp.nPoints(), -1);
    patchConstraints.setSize(pp.nPoints());

    // Overall bb
    treeBoundBox overallBb(ppLocalPoints);
    // Grow a bit since points move
    overallBb.inflate(1.5);

    // Additional pp search structures
    const indexedOctree<treeDataPoint> ppPointTree
    (
        treeDataPoint(ppLocalPoints),
        overallBb,
        10, // maxLevel
        10, // leafSize
        3.0 // duplicity
    );


    // Patch face centre to nearest surface location
    {
        const auto& localFaces = pp.localFaces();

        // Calculate local geometry with provided (smoothed) points
        pointField ppFaceCentres(pp.size());
        vectorField ppFaceNormals(pp.size());
        forAll(pp, i)
        {
            ppFaceCentres[i] = localFaces[i].centre(ppLocalPoints);
            ppFaceNormals[i] = localFaces[i].areaNormal(ppLocalPoints);
        }

        labelList surface;
        List<pointIndexHit> info;
        vectorField normal;
        labelList region;
        surfaces.findNearestRegion
        (
            surfaces.surfaces(),
            ppFaceCentres,
            scalarField(pp.size(), sqr(GREAT)),

            surface,
            info,
            region,
            normal
        );
        patchFaceToSurfaceRegion.setSize(pp.size(), -1);
        forAll(info, i)
        {
            if (info[i].hit())
            {
                //const scalar cosAngle = (normal[i] & ppFaceNormals[i]);
                //if (cosAngle > 0)
                {
                    // Store nearest surface point
                    patchFaceToSurfaceRegion[i] = region[i];
                }
            }
        }
    }


    // Patch vertex to nearest surface location
    {
        labelList surface;
        List<pointIndexHit> info;
        vectorField normal;
        labelList region;
        surfaces.findNearestRegion
        (
            surfaces.surfaces(),
            ppLocalPoints,
            scalarField(pp.nPoints(), sqr(GREAT)),

            surface,
            info,
            region,
            normal
        );
        forAll(info, i)
        {
            if (info[i].hit())
            {
                //const scalar cosAngle = (normal[i] & ppPointNormals[i]);
                //// Initial filtering: Ignore hits on 'wrong' side
                //if (cosAngle > 0)
                {
                    // Store nearest surface point & normal
                    patchToNearest[i] = info[i].hitPoint();
                    patchConstraints[i].applyConstraint(normal[i]);
                    patchToNormal[i] = normal[i];
                    //patchToPoint[i] = ? Only if nearest is a feature point
                    patchToRegion[i] = region[i];
                }
            }
        }
    }



    // Match feature points
    // ~~~~~~~~~~~~~~~~~~~~

    if (doFeatures)
    {
        // Feature point:
        // - find nearest mesh point to each surface feature point
        // - attract a bit and attract edge connected points to the emanating
        //   feature edge directions
        matchFeaturePoints
        (
            pp,
            ppLocalPoints,
            ppPointTree,
            patchFaceToSurfaceRegion,

            // Index of surface
            patchToFeature,
            // Index of surface feature point (or -1)
            patchToPoint,
            // Point on nearest surface
            patchToNearest,
            // Normal on nearest surface
            patchToNormal,
            // Region on nearest surface
            patchToRegion,
            // Constraints at surface (1 = surface, 2 = featedge, 3 = featpoint)
            patchConstraints
        );

        // Wip. See if points neighbouring feature points can be dragged a bit
        // to improve initial guess for feature edge matching.
        forAll(patchConstraints, ppPointi)
        {
            if (patchConstraints[ppPointi].first() == 3)
            {
                // Attract point
                const vector d =
                    patchToNearest[ppPointi] - ppLocalPoints[ppPointi];
                ppLocalPoints[ppPointi] += d;

                // Attract neighbouring points so follow on feature
                // edge matching goes better.
                for (const label edgei : pp.pointEdges()[ppPointi])
                {
                    const edge& e = pp.edges()[edgei];
                    const label otherPpPointi = e.otherVertex(ppPointi);
                    if (patchConstraints[otherPpPointi].first() < 3)
                    {
                        ppLocalPoints[otherPpPointi] += 0.5*d;
                    }
                }
            }
        }
        if (debug & meshRefinement::MESH)
        {
            OBJstream os
            (
                meshRefiner_.mesh().time().timePath()
              / "pp_smoothed_featPoints.obj"
            );
            os.write(pp.localFaces(), ppLocalPoints, false);
            Pout<< "Write pp after feat-point attraction to " << os.name()
                << " nVertices:" << os.nVertices() << endl;
        }
    }



    // Match feature edges
    // ~~~~~~~~~~~~~~~~~~~

    if (doFeatures)
    {
        // Detect faces with points (or neighbouring faces?) on differing
        // regions.
        // These are the candidate faces for feature-edge matching
        bitSet isCandidateFace(pp.size());
        determineCandidateFaces
        (
            featCosAngle,

            pp,
            // Point on nearest surface
            patchToNearest,
            // Normal on nearest surface
            patchToNormal,
            // Region on nearest surface
            patchToRegion,
            // Constraints at surface (1 = surface, 2 = featedge, 3 = featpoint)
            patchConstraints,

            isCandidateFace
        );


        // Match feature edges
        // ~~~~~~~~~~~~~~~~~~~

        bitSet isDiagonalCutFace(pp.size());

        forAll(patchConstraints, ppPointi)
        {
            // Search for featpoint (nConstraints = 3) to start walk from
            if (patchConstraints[ppPointi].first() < 3)
            {
                continue;
            }

            matchFeatureEdges
            (
                // Mesh faces to snap
                pp,
                // Per patch point the current (attracted) coordinate
                ppLocalPoints,
                // Whether pp face has different attraction
                isCandidateFace,

                // Starting patch feature point
                ppPointi,
                // Starting surface feature
                patchToFeature[ppPointi],
                // Starting surface feature point
                patchToPoint[ppPointi],


                // Per patch face:

                // diagonal cut
                patchFaceToDiagonal,

                // Per patch point:

                // Index of surface edge
                patchToEdge,
                // Location on nearest edge
                patchToNearest,
                // Average face normal on nearest edge
                patchToNormal,
                // Constraints at surface (1 = surface, 2 = featedge, 3 = featpoint)
                patchConstraints
            );
        }
    }


    // Match remaining to nearest surface
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    matchSurfaces
    (
        pp,
        ppLocalPoints,

        // Per patch point:

        // Index of surface 'face' (or -1)
        patchToFace,
        // Region on nearest surface
        patchToRegion,
        // Point on nearest surface
        patchToNearest,
        // Normal on nearest surface
        patchToNormal,
        // Constraints at surface (1 = surface, 2 = featedge, 3 = featpoint)
        patchConstraints
    );
}


void Foam::snappySnapDriver::matchToSurfaceMerged
(
    const bool doFeatures,
    const scalar featCosAngle,

    const indirectPrimitivePatch& pp,   // addressing only
    pointField& ppLocalPoints,          // points (modified during feat-edges)
    const vectorField& ppPointNormals,  // normal

    // Per string index a starting pp point
    Map<label>& surfFeatEdgeRegionToStart,


    //- Per patch face

    // Region on nearest surface
    labelList& patchFaceToSurfaceRegion,

    // Face to diagonal cut
    Map<edge>& patchFaceToDiagonal,

    //- Per patch point

    // Point on nearest surface
    pointField& patchToNearest,
    // Normal on nearest surface
    vectorField& patchToNormal,
    // Index of nearest surface
    labelList& patchToSurface,
    // Region on nearest surface
    labelList& patchToRegion,
    // Surface 'face' index
    labelList& patchToFace,
    // Surface 'edge' index
    labelList& patchToEdge,
    // Surface 'feature' index
    labelList& patchToFeature,
    // Surface feature point (if patch point is attracted)
    labelList& patchToPoint,
    // Constraints at surface (1 = surface, 2 = featedge, 3 = featpoint)
    List<pointConstraint>& patchConstraints
) const
{
    if (Pstream::parRun())
    {
        // The matching uses patch-based data only and not mesh-based data
        // so we can just send over the patch to the master and do all
        // matching there. The only exception is in the point-averaging
        // during walking the edges but since we're only doing it on
        // the master processor it does not matter.

        Info<< "Collecting all patch faces onto master to do analysis ..."
            << nl << endl;

        // Send patch to master
        autoPtr<globalIndex> globalFacesPtr;
        autoPtr<mapDistribute> pointMapPtr;
        labelList ppToMerged;
        autoPtr<PrimitivePatch<faceList, pointField>> mergedPatchPtr =
            gatherAndMerge
            (
                meshRefiner_.mesh(),
                pp,

                globalFacesPtr,
                pointMapPtr,
                ppToMerged
            );

        // Convert to indirectPrimitivePatch just for convenience. Keep same
        // face ordering so point ordering should be same as well ...
        const indirectPrimitivePatch mergedPatch
        (
            IndirectList<face>
            (
                mergedPatchPtr(),
                identity(mergedPatchPtr().size())
            ),
            mergedPatchPtr().points()
        );
        const label nPts = mergedPatch.nPoints();

        // Map pp data onto merged
        pointField mergedPoints(ppLocalPoints);
        pointMapPtr().reverseDistribute(nPts, mergedPoints);
        vectorField mergedNormals(ppPointNormals);
        pointMapPtr().reverseDistribute(nPts, mergedNormals);

        // Do matching

        // patch face to nearest surface
        labelList mergedFaceToSurfaceRegion;
        // patch face to diagonal split info
        Map<edge> mergedFaceToDiagonal;
        // Point on nearest surface
        pointField mergedToNearest(nPts, vector::uniform(GREAT));
        // Normal on nearest surface
        vectorField mergedToNormal(nPts, vector::uniform(GREAT));
        // Index of nearest surface
        labelList mergedToSurface(nPts, -1);
        // Region on nearest surface
        labelList mergedToRegion(nPts, -1);
        // Patch point to surface 'face'
        labelList mergedToFace(nPts, -1);
        // Surface edge 'region' (if patch point is attracted to it)
        labelList mergedToEdge(nPts, -1);
        // Surface 'feature' index
        labelList mergedToFeature(nPts, -1);
        // Surface feature point (if patch point is attracted)
        labelList mergedToPoint(nPts, -1);
        // Constraints at surface (1 = surface, 2 = featedge, 3 = featpoint)
        List<pointConstraint> mergedConstraints(nPts);

        if (Pstream::master())
        {
            matchToSurface
            (
                doFeatures,
                featCosAngle,

                mergedPatch,                // addressing only
                mergedPoints,               // points (modified during matching!)
                mergedNormals,              // normal

                // Per string index a starting pp point
                surfFeatEdgeRegionToStart,


                //- Per patch face

                // Region on nearest surface
                mergedFaceToSurfaceRegion,
                // Diagonal split info
                mergedFaceToDiagonal,

                //- Per merged point

                // Point on nearest surface
                mergedToNearest,
                // Normal on nearest surface
                mergedToNormal,
                // Index of nearest surface
                mergedToSurface,
                // Region on nearest surface
                mergedToRegion,
                // Surface 'face' index
                mergedToFace,
                // Surface 'edge' index
                mergedToEdge,
                // Surface 'feature' index
                mergedToFeature,
                // Surface feature point (if patch point is attracted)
                mergedToPoint,
                // Constraints at surface (1 = surface, 2 = featedge, 3 = featpoint)
                mergedConstraints
            );
        }


        // Map back to pp
        // ~~~~~~~~~~~~~~

        //- Input data

            // Send master data to processors
            pointMapPtr().distribute(mergedPoints);
            // Shuffle according to ppToMerged
            ppLocalPoints = UIndirectList<point>(mergedPoints, ppToMerged);

        //- Surface-based data

            // TBD. Currently not used anyway
            // surfFeatEdgeRegionToStart

        //- Face-based data

            patchFaceToSurfaceRegion.setSize(pp.size(), -1);
            globalFacesPtr().scatter
            (
                mergedFaceToSurfaceRegion,
                patchFaceToSurfaceRegion
            );

            // patchFaceToDiagonal is the local diagonal cut per face. Since
            // we don't rotate the faces, the diagonal cut should be the same.
            {
                // Expand Map into List for scatter and then back into Map
                const edge null(-1, -1);
                edgeList mergedFaceToDiagonalList
                (
                    mergedPatch.size(),
                    null
                );
                forAllConstIters(mergedFaceToDiagonal, iter)
                {
                    mergedFaceToDiagonalList[iter.key()] = iter();
                }
                edgeList ppFaceToDiagonalList(pp.size());
                globalFacesPtr().scatter
                (
                    mergedFaceToDiagonalList,
                    ppFaceToDiagonalList
                );
                forAll(ppFaceToDiagonalList, facei)
                {
                    if (ppFaceToDiagonalList[facei] != null)
                    {
                        patchFaceToDiagonal.insert
                        (
                            facei,
                            ppFaceToDiagonalList[facei]
                        );
                    }
                }
            }

        //- Point-based data

            // Point on nearest surface
            pointMapPtr().distribute(mergedToNearest);
            patchToNearest = UIndirectList<point>(mergedToNearest, ppToMerged);
            // Normal on nearest surface
            pointMapPtr().distribute(mergedToNormal);
            patchToNormal = UIndirectList<vector>(mergedToNormal, ppToMerged);
            // Index of nearest surface
            pointMapPtr().distribute(mergedToSurface);
            patchToSurface = UIndirectList<label>(mergedToSurface, ppToMerged);
            // Region on nearest surface
            pointMapPtr().distribute(mergedToRegion);
            patchToRegion = UIndirectList<label>(mergedToRegion, ppToMerged);
            // Surface 'face' index
            pointMapPtr().distribute(mergedToFace);
            patchToFace = UIndirectList<label>(mergedToFace, ppToMerged);
            // Surface 'edge' index
            pointMapPtr().distribute(mergedToEdge);
            patchToEdge = UIndirectList<label>(mergedToEdge, ppToMerged);
            // Surface 'feature' index
            pointMapPtr().distribute(mergedToFeature);
            patchToFeature = UIndirectList<label>(mergedToFeature, ppToMerged);
            // Surface feature point (if patch point is attracted)
            pointMapPtr().distribute(mergedToPoint);
            patchToPoint = UIndirectList<label>(mergedToPoint, ppToMerged);
            // Constraints at surface (1 = surface, 2 = featedge, 3 = featpoint)
            pointMapPtr().distribute(mergedConstraints);
            patchConstraints =
                UIndirectList<pointConstraint>(mergedConstraints, ppToMerged);
    }
    else
    {
        // Just do matching directly
        matchToSurface
        (
            doFeatures,
            featCosAngle,

            pp,                     // addressing only
            ppLocalPoints,          // points (modified during matching!)
            ppPointNormals,         // normal

            // Per string index a starting pp point
            surfFeatEdgeRegionToStart,

            //- Per patch face

            // Region on nearest surface
            patchFaceToSurfaceRegion,
            // Diagonal split info
            patchFaceToDiagonal,

            //- Per patch point

            // Point on nearest surface
            patchToNearest,
            // Normal on nearest surface
            patchToNormal,
            // Index of nearest surface
            patchToSurface,
            // Region on nearest surface
            patchToRegion,
            // Surface 'face' index
            patchToFace,
            // Surface 'edge' index
            patchToEdge,
            // Surface 'feature' index
            patchToFeature,
            // Surface feature point (if patch point is attracted)
            patchToPoint,
            // Constraints at surface (1 = surface, 2 = featedge, 3 = featpoint)
            patchConstraints
        );
    }
}


// Topological helper routines
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Foam::snappySnapDriver::splitFace
(
    const face& f,
    const UList<label>& cuts,
    DynamicList<face>& faces
) const
{
    faces.clear();

    label fp = cuts[0];

    DynamicList<label> verts(3);

    forAll(cuts, cuti)
    {
        verts.append(f[fp]);
        while (true)
        {
            fp = f.fcIndex(fp);
            verts.append(f[fp]);

            if (fp == cuts[cuts.fcIndex(cuti)])
            {
                faces.append(face(verts));
                verts.clear();
                break;
            }
        }
    }
}


Foam::label Foam::snappySnapDriver::findSplitFaceRegion
(
    const bool useFaceRegion,
    const indirectPrimitivePatch& pp,
    const labelList& patchFaceToSurfaceRegion,
    const List<pointConstraint>& patchConstraints,
    const labelList& patchToRegion,
    const label facei,
    const label startFp,
    const label endFp
) const
{
    // Walks from startFp to endFp on facei
    // and finds the 'best' point to take surface
    // region from.
    const auto& faces = pp.localFaces();
    const face& f = faces[facei];
    const label n =
    (
        endFp >= startFp
      ? endFp-startFp+1
      : endFp+f.size()-startFp+1
    );
    // Find any face-constrained points
    label fp = startFp;
    for (label i = 0; i < n; i++)
    {
        const label ppPointi = f[fp];
        if
        (
            patchConstraints[ppPointi].first() <= 1
        &&  patchToRegion[ppPointi] != -1
        )
        {
            // Found face-constrained point with
            // attraction
            return patchToRegion[ppPointi];
        }
        fp = f.fcIndex(fp);
    }


    if (useFaceRegion)
    {
        const auto& fEdges = pp.faceEdges()[facei];
        fp = startFp;
        for (label i = 0; i < n-1; i++)
        {
            const label edgei = fEdges[i];
            for (const label nbrFacei : pp.edgeFaces()[edgei])
            {
                if (nbrFacei != facei)
                {
                    // //- Recurse to find connected face that has valid
                    // //  region.
                    // const label nbrRegioni = findSplitFaceRegion
                    // (
                    //     false,  // recurse,
                    //     pp,
                    //     patchConstraints,
                    //     patchToRegion,
                    //     nbrFacei,
                    //     0,  //startFp,
                    //     faces[nbrFacei].size()-1
                    // );
                    const label nbrRegioni = patchFaceToSurfaceRegion[nbrFacei];
                    if (nbrRegioni != -1)
                    {
                        return nbrRegioni;
                    }
                }
            }
            fp = f.fcIndex(fp);
        }
    }
    return -1;
}


Foam::tmp<Foam::scalarField>
Foam::snappySnapDriver::wantedThickness
(
    const indirectPrimitivePatch& pp,
    const scalar cellSizeFraction
) const
{
    fvMesh& mesh = meshRefiner_.mesh();

    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    const labelList& owner = mesh.faceOwner();

    // Undistorted edge length
    const scalar edge0Len =
        meshRefiner_.meshCutter().level0EdgeLength();


    tmp<scalarField> tthickness(tmp<scalarField>::New(pp.nPoints()));
    scalarField& thickness = tthickness.ref();

    labelList maxPointLevel(pp.nPoints(), labelMin);

    forAll(pp, i)
    {
        label ownLevel = cellLevel[owner[pp.addressing()[i]]];

        const face& f = pp.localFaces()[i];

        forAll(f, fp)
        {
            maxPointLevel[f[fp]] = max(maxPointLevel[f[fp]], ownLevel);
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        maxPointLevel,
        maxEqOp<label>(),
        labelMin            // null value
    );

    forAll(thickness, pointi)
    {
        const scalar edgeLen = edge0Len/(1<<maxPointLevel[pointi]);
        thickness[pointi] = cellSizeFraction*edgeLen;
    }
    return tthickness;
}


const Foam::pointMesh& Foam::snappySnapDriver::makePointMesh
(
    const indirectPrimitivePatch& pp,
    const pointConstraintList& pointConstraints,
    const labelList& patchIDs,
    const word& allSurfPatchName,
    const word& allEdgePatchName,
    const word& allPointPatchName
) const
{
    fvMesh& mesh = meshRefiner_.mesh();

    if (pointConstraints.size() != pp.nPoints())
    {
        FatalErrorInFunction<< "pointConstraints:" << pointConstraints.size()
            << " pp:" << pp.nPoints() << exit(FatalError);
    }


    // Expand pointConstraints to all meshPoints
    pointConstraintList meshPointConstraints(mesh.nPoints());
    UIndirectList<pointConstraint>(meshPointConstraints, pp.meshPoints()) =
        pointConstraints;

    const auto combineConstraints = [&]
    (
        pointConstraint& x,
        const pointConstraint& y
    )
    {
        x.combine(y);
    };


    syncTools::syncPointList
    (
        mesh,
        meshPointConstraints,
        combineConstraints,
        pointConstraint(),
        dummyTransform()
    );


    // Sort point constraints:
    // - 0 : does not happen?
    // - 1 : attract to surface
    // - 2 : attract to feature edge
    // - 3 : attract to feature point

    DynamicList<label> surfaceMeshPoints(pointConstraints.size());
    DynamicList<label> featEdgeMeshPoints(pointConstraints.size());
    DynamicList<label> featPointMeshPoints(pointConstraints.size());
    forAll(meshPointConstraints, pointi)
    {
        if (meshPointConstraints[pointi].first() == 1)
        {
            surfaceMeshPoints.append(pointi);
        }
        else if (meshPointConstraints[pointi].first() == 2)
        {
            featEdgeMeshPoints.append(pointi);
        }
        else if (meshPointConstraints[pointi].first() == 3)
        {
            featPointMeshPoints.append(pointi);
        }
    }

    // Lookup / construct (from polyPatches) the pointMesh
    const auto& pMesh = pointMesh::New(meshRefiner_.mesh());

    bool changed = false;

    pointBoundaryMesh& pointBm =
        const_cast<pointBoundaryMesh&>(pMesh.boundary());

    // Note: there are two choices on how to handle the single-patch points:
    //       - introduce additional meshPoint patch ('allSurfPatchName')
    //         and have patch ordering take care of overriding normal,
    //         fixedValue behaviour.
    //       - override original facePointPatches to be meshPointPatch
    //         with optionally only the internal points (so each point is
    //         handled only by a single pointPatch)
    for (const label patchi : patchIDs)
    {
        const pointPatch& ppp = pointBm[patchi];
        if (!isA<meshPointPatch>(ppp))
        {
            const auto* fpPatchPtr = isA<facePointPatch>(ppp);
            if (fpPatchPtr)
            {
                // Get primitive patch
                const auto& pp = fpPatchPtr->patch();

                // Keep only 'internal' points. Could leave feature points
                // in but already handled in additional pointPatches added
                // below so filter out any non-planar points. Problem
                // of this filtering though is that the pointDisplacement field
                // cannot be visualised in paraview anymore ...
                const auto& mp = pp.meshPoints();
                DynamicList<label> patchMeshPoints(mp.size());
                for (const label meshPointi : mp)
                {
                    if (meshPointConstraints[meshPointi].first() == 1)
                    {
                        patchMeshPoints.append(meshPointi);
                    }
                }

                // Override patch
                const_cast<pointBoundaryMesh&>(pointBm).set
                (
                    patchi,
                    new meshPointPatch
                    (
                        pp.name(),
                        patchMeshPoints,        // mp
                        List<pointConstraint>
                        (
                            meshPointConstraints,
                            patchMeshPoints     // mp
                        ),
                        patchi,
                        pointBm,
                        meshPointPatch::typeName
                    )
                );
                changed = true;
            }
        }
    }

    // Check if already has constraint patches
    label surfPatchi = -1;
    if (allSurfPatchName.size())
    {
        surfPatchi = pointBm.findPatchID(allSurfPatchName);
        if (surfPatchi != -1)
        {
            // Delete patch. TBD: clear patchGroup
            pointBm.set(surfPatchi, nullptr);
            changed = true;
        }
    }
    const label edgePatchi = pointBm.findPatchID(allEdgePatchName);
    if (edgePatchi != -1)
    {
        // Delete patch. TBD: clear patchGroup
        pointBm.set(edgePatchi, nullptr);
        changed = true;
    }
    const label pointPatchi = pointBm.findPatchID(allPointPatchName);
    if (pointPatchi != -1)
    {
        // Delete patch. TBD: clear patchGroup
        pointBm.set(pointPatchi, nullptr);
        changed = true;
    }

    // Add additional point patches in order:
    // - original polyPatch based
    // - (optional) surface-plane based constraint
    // - featEdge-based constraints
    // - featPoint-based constraints (so can override edge-based constraints)
    // - processor boundaries (or should this be all constraint patches, e.g.
    //   symmetry plane. Note: hopefully this is already handled in the
    //   feature extraction ...)

    const label nSurfMeshPoints =
        returnReduce(surfaceMeshPoints.size(), sumOp<label>(), mesh.comm());
    if (allSurfPatchName.size() && nSurfMeshPoints)
    {
        Info<< "Adding surface plane patch " << allSurfPatchName
            << " with " << nSurfMeshPoints << " points" << endl;
        if (surfPatchi != -1)
        {
            // Override patch
            const_cast<pointBoundaryMesh&>(pointBm).set
            (
                surfPatchi,
                new meshPointPatch
                (
                    allSurfPatchName,
                    surfaceMeshPoints,
                    List<pointConstraint>
                    (
                        meshPointConstraints,
                        surfaceMeshPoints
                    ),
                    surfPatchi,
                    pointBm,
                    meshPointPatch::typeName
                )
            );
        }
        else
        {
            // Append
            const_cast<pointBoundaryMesh&>(pointBm).push_back
            (
                new meshPointPatch
                (
                    allSurfPatchName,
                    surfaceMeshPoints,
                    List<pointConstraint>
                    (
                        meshPointConstraints,
                        surfaceMeshPoints
                    ),
                    pointBm.size(),
                    pointBm,
                    meshPointPatch::typeName
                )
            );
        }
        changed = true;
    }
    const label nEdgeMeshPoints =
        returnReduce(featEdgeMeshPoints.size(), sumOp<label>(), mesh.comm());
    if (nEdgeMeshPoints)
    {
        Info<< "Adding surface edge patch " << allEdgePatchName
            << " with " << nEdgeMeshPoints << " points" << endl;
        if (edgePatchi != -1)
        {
            // Override patch
            const_cast<pointBoundaryMesh&>(pointBm).set
            (
                edgePatchi,
                new meshPointPatch
                (
                    allEdgePatchName,
                    featEdgeMeshPoints,
                    List<pointConstraint>
                    (
                        meshPointConstraints,
                        featEdgeMeshPoints
                    ),
                    edgePatchi,
                    pointBm,
                    meshPointPatch::typeName
                )
            );
        }
        else
        {
            // Append
            const_cast<pointBoundaryMesh&>(pointBm).push_back
            (
                new meshPointPatch
                (
                    allEdgePatchName,
                    featEdgeMeshPoints,
                    List<pointConstraint>
                    (
                        meshPointConstraints,
                        featEdgeMeshPoints
                    ),
                    pointBm.size(),
                    pointBm,
                    meshPointPatch::typeName
                )
            );
        }
        changed = true;
    }
    const label nFeatMeshPoints =
        returnReduce(featPointMeshPoints.size(), sumOp<label>(), mesh.comm());
    if (nFeatMeshPoints)
    {
        Info<< "Adding surface point patch " << allPointPatchName
            << " with " << nFeatMeshPoints << " points" << endl;
        if (pointPatchi != -1)
        {
            // Override patch
            const_cast<pointBoundaryMesh&>(pointBm).set
            (
                pointPatchi,
                new meshPointPatch
                (
                    allPointPatchName,
                    featPointMeshPoints,
                    List<pointConstraint>
                    (
                        meshPointConstraints,
                        featPointMeshPoints
                    ),
                    pointPatchi,
                    pointBm,
                    meshPointPatch::typeName
                )
            );
        }
        else
        {
            // Append
            const_cast<pointBoundaryMesh&>(pointBm).push_back
            (
                new meshPointPatch
                (
                    allPointPatchName,
                    featPointMeshPoints,
                    List<pointConstraint>
                    (
                        meshPointConstraints,
                        featPointMeshPoints
                    ),
                    pointBm.size(),
                    pointBm,
                    meshPointPatch::typeName
                )
            );
        }
        changed = true;
    }

    // Shuffle into order
    labelList oldToNew(pointBm.size());
    label newPatchi = 0;
    forAll(pointBm, patchi)
    {
        if (!isA<processorPointPatch>(pointBm[patchi]))
        {
            oldToNew[patchi] = newPatchi++;
        }
    }
    forAll(pointBm, patchi)
    {
        if (isA<processorPointPatch>(pointBm[patchi]))
        {
            oldToNew[patchi] = newPatchi++;
        }
    }
    pointBm.reorder(oldToNew, true);

    if (changed)
    {
        const_cast<pointMesh&>(pMesh).setInstance(meshRefiner_.timeName());
        const_cast<pointMesh&>(pMesh).writeOpt(IOobject::AUTO_WRITE);
        pointBm.writeOpt(IOobject::AUTO_WRITE);
    }
    return pMesh;
}


Foam::autoPtr<Foam::displacementMotionSolver>
Foam::snappySnapDriver::makeMotionSolver
(
    const pointMesh& pMesh,
    const autoPtr<indirectPrimitivePatch>& pp,
    const dictionary& snapDict,
    const labelList& adaptPatchIDs
) const
{
    // Note: pointConstraints already in pointBoundaryMesh as
    //       meshPointPatch 'constraint' patches
    fvMesh& mesh = meshRefiner_.mesh();

    tmp<pointVectorField> tallDisp
    (
        meshRefinement::makeDisplacementField
        (
            pMesh,
            adaptPatchIDs,
            fixedValuePointPatchVectorField::typeName
        )
    );

    // Make sure the pointDisplacement is not registered (since
    // displacementMotionSolver itself holds it)
    tallDisp.ref().checkOut();


    // Add faceZone for additional face constraints. To be handled by the
    // smoothers. Avoids having to create baffles.
    if (pp)
    {
        surfaceZonesInfo::addFaceZone
        (
            "boundaryFaces",
            pp().addressing(),
            boolList(pp().addressing().size(), false),
            mesh
        );
    }

    autoPtr<displacementMotionSolver> motionPtr
    (
        displacementMotionSolver::New
        (
            snapDict.get<word>("solver"),
            mesh,
            IOdictionary
            (
                IOobject
                (
                    "dynamicMeshDict",  //"motionSolverDict",
                    pMesh.thisDb().time().constant(),
                    pMesh.thisDb(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                snapDict
            ),
            tallDisp(),
            pointIOField
            (
                IOobject
                (
                    "points0",
                    pMesh.thisDb().time().constant(),
                    polyMesh::meshSubDir,
                    pMesh.thisDb(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh.points()
            )
        )
    );
    return motionPtr;
}


void Foam::snappySnapDriver::setDisplacement
(
    const indirectPrimitivePatch& pp,
    const pointField& patchDisp,    // displacement w.r.t. current mesh
    const labelList& adaptPatchIDs,
    const pointField& points0,
    pointVectorField& fld           // displacement w.r.t. points0
)
{
    const pointMesh& pMesh = fld.mesh();
    const pointField& points = pMesh().points();
    const labelList& meshPoints = pp.meshPoints();

    if
    (
        (points0.size() != points.size())
     || (points0.size() != fld.size())
     || (points0.size() != pMesh.size())
     || (meshPoints.size() != patchDisp.size())
    )
    {
        FatalErrorInFunction
            << "Sizing :"
            << " points0.size():" << points0.size()
            << " points.size():" << points.size()
            << " fld.size():" << fld.size()
            << " patchDisp.size():" << patchDisp.size()
            << " meshPoints.size():" << meshPoints.size()
            << " mesh.nPoints():" << pMesh.size()
            << exit(FatalError);
    }


    // Problem is that the patchDisp might not be consistent
    // across multi-patch points so expand to mesh points

    auto& iF = fld.internalFieldRef();
    forAll(meshPoints, patchPointi)
    {
        const label meshPointi = meshPoints[patchPointi];
        iF[meshPointi] =
            patchDisp[patchPointi]+points[meshPointi]-points0[meshPointi];
    }

    // Assign to bc
    pointVectorField::Boundary& bfld = fld.boundaryFieldRef();
    for (const label patchi : adaptPatchIDs)
    {
        bfld[patchi] == pointField(iF, bfld[patchi].patch().meshPoints());
    }

    // Reset internalField since has to be up to date with mesh points
    iF.vectorField::operator=(vector::zero);

    // Apply multi-patch constraints. Problem: most patches originate from
    // meshing so are type 'wall'. The pointConstraints are only on any points
    // remaining from the starting mesh.
    const pointConstraints& pcs = pointConstraints::New(pMesh);
    pcs.constrainDisplacement(fld, true);
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::snappySnapDriver::addBufferLayers
(
    const indirectPrimitivePatch& pp,
    const pointField& thickness,
    // Layer mesh modifier
    addPatchCellLayer& addLayer
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    // Introduce single layer of cells. Straight from snappyLayerDriver

    // Global face indices engine
    const globalIndex globalFaces(mesh.nFaces());

    // Determine extrudePatch.edgeFaces in global numbering (so across
    // coupled patches). This is used only to string up edges
    // between coupled
    // faces (all edges between same (global)face indices get extruded).
    labelListList edgeGlobalFaces
    (
        addPatchCellLayer::globalEdgeFaces
        (
            mesh,
            globalFaces,
            pp
        )
    );



    // Use global edge - face connectivity to
    // - disable any non-manifold extrusions
    // - boundary edges can be extruded - is handled by addPatchCellLayer
    // - what about setting numLayers to 0 to disable layers? Should that
    //   affect buffer layer addition? Guess not since buffer layer is
    //   to aid snapping ...
    labelList nFaceLayers(pp.size(), 1);
    labelList nPointLayers(pp.nPoints(), 1);
    forAll(edgeGlobalFaces, edgei)
    {
        if (edgeGlobalFaces[edgei].size() > 2)
        {
            const edge& e = pp.edges()[edgei];

            nPointLayers[e[0]] = 0;
            nPointLayers[e[1]] = 0;
        }
        //else if (edgeGlobalFaces[edgei].size() == 1)
        //{
        //    const edge& e = pp.edges()[edgei];
        //    const edge meshE(pp.meshPoints()[e[0]], pp.meshPoints()[e[1]]);
        //
        //    nPointLayers[e[0]] = 0;
        //    nPointLayers[e[1]] = 0;
        //}
        else if (edgeGlobalFaces[edgei].size() == 0)
        {
            const edge& e = pp.edges()[edgei];
            const edge meshE(pp.meshPoints()[e[0]], pp.meshPoints()[e[1]]);
            FatalErrorInFunction << "Edge:" << meshE
                << " At:" << meshE.line(mesh.points())
                << " has no faces!" << exit(FatalError);
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        nPointLayers,
        minEqOp<label>(),
        labelMax            // null value
    );

    // Check for extrusion of baffles i.e. all edges of a face having the
    // same two neighbouring faces, one of which is the current face.
    // This is detected locally already before - this test is for the extremely
    // rare occurrence where the baffle faces are on different processors.
    {
        label nBaffleFaces = 0;

        const labelListList& faceEdges = pp.faceEdges();
        forAll(pp, facei)
        {
            const labelList& fEdges = faceEdges[facei];

            const labelList& globFaces0 = edgeGlobalFaces[fEdges[0]];
            if (globFaces0.size() == 2)
            {
                const edge e0(globFaces0[0], globFaces0[1]);
                bool isBaffle = true;
                for (label fp = 1; fp < fEdges.size(); ++fp)
                {
                    const labelList& globFaces = edgeGlobalFaces[fEdges[fp]];
                    if
                    (
                        (globFaces.size() != 2)
                     || (edge(globFaces[0], globFaces[1]) != e0)
                    )
                    {
                        isBaffle = false;
                        break;
                    }
                }

                if (isBaffle)
                {
                    const face& f = pp.localFaces()[facei];
                    for (const auto& fp : f)
                    {
                        nPointLayers[fp] = 0;
                    }
                    ++nBaffleFaces;
                }
            }
        }

        if (debug)
        {
            reduce(nBaffleFaces, sumOp<label>());

            if (nBaffleFaces)
            {
                Info<< "Set number of buffer layers to zero for all points on "
                    << nBaffleFaces << " baffle faces" << endl;
            }
        }
    }

    forAll(pp.localFaces(), facei)
    {
        const face& f = pp.localFaces()[facei];
        const UIndirectList<label> pLayers(nPointLayers, f);
        if (!pLayers.found(label(1)))
        {
            nFaceLayers[facei] = 0;
        }
    }



    // Determine patches for extruded boundary edges. Adds any
    // additional processor patches (since extruding coupled edge can cause
    // additional connectivity)

    labelList edgePatchID;
    labelList edgeZoneID;
    boolList edgeFlip;
    labelList inflateFaceID;
    snappyLayerDriver::determineSidePatches
    (
        meshRefiner_,
        globalFaces,
        edgeGlobalFaces,
        pp,

        edgePatchID,
        edgeZoneID,
        edgeFlip,
        inflateFaceID
    );


    // Mesh topo change engine. Insert current mesh.
    polyTopoChange meshMod(mesh);


    // Add topo regardless of whether extrudeStatus is extruderemove.
    // Not add layer if patchDisp is zero.

    addLayer.setRefinement
    (
        globalFaces,
        edgeGlobalFaces,

        scalarField(pp.nPoints(), 1),  // expansion ratio
        pp,
        bitSet(pp.size()),          // no flip

        edgePatchID,    // boundary patch for extruded boundary edges
        edgeZoneID,     // zone for extruded edges
        edgeFlip,
        inflateFaceID,


        labelList(0),   // exposed patchIDs, not used for adding layers
        nFaceLayers,
        nPointLayers,
        thickness,  //patchAttraction,

        meshMod
    );

    // Hack to remove meshPhi - mapped incorrectly. TBD.
    mesh.moving(false);
    mesh.clearOut();

    // Apply the stored topo changes to the current mesh.
    autoPtr<mapPolyMesh> mapPtr = meshMod.changeMesh(mesh, false);
    mapPolyMesh& map = *mapPtr;

    // Update fields
    mesh.updateMesh(map);

    // Move mesh (since morphing does not do this)
    if (map.hasMotionPoints())
    {
        mesh.movePoints(map.preMotionPoints());
    }
    else
    {
        // Hack to remove meshPhi - mapped incorrectly. TBD.
        mesh.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh.setInstance(meshRefiner_.timeName());

    // Update numbering on layer
    addLayer.updateMesh
    (
        map,
        identity(pp.size()),
        identity(pp.nPoints())
    );

    // Update intersections
    const labelListList addedCells(addLayer.addedCells());
    bitSet isChangedFace(mesh.nFaces());
    for (const labelList& faceToCells : addedCells)
    {
        for (const label celli : faceToCells)
        {
            isChangedFace.set(mesh.cells()[celli]);
        }
    }
    meshRefiner_.updateMesh(map, isChangedFace.toc());

    if (debug & meshRefinement::MESH)
    {
        const_cast<Time&>(mesh.time())++;
        Info<< "Writing mesh-with-layer to time "
            << meshRefiner_.timeName() << endl;
        meshRefiner_.write
        (
            meshRefinement::debugType(debug),
            meshRefinement::writeType
            (
                meshRefinement::writeLevel()
              | meshRefinement::WRITEMESH
            ),
            meshRefiner_.timeName()
        );
    }
    return mapPtr;
}


// void Foam::snappySnapDriver::doSnapBufferLayers
// (
//     const dictionary& snapDict,
//     const dictionary& motionDict,
//     const meshRefinement::FaceMergeType mergeType,
//     const scalar featureCos,
//     const scalar planarAngle,
//     const snapParameters& snapParams
// )
// {
//     addProfiling(snap, "snappyHexMesh::snap");
//     fvMesh& mesh = meshRefiner_.mesh();

//     Info<< nl
//         << "Morphing phase" << nl
//         << "--------------" << nl
//         << endl;


//     // Name of pointPatch for any points attracted to feature edges
//     //const word allEdgePatchName("boundaryEdges");
//     // Name of pointPatch for any points attracted to feature points
//     //const word allPointPatchName("boundaryPoints");


//     // faceZone handling
//     // ~~~~~~~~~~~~~~~~~
//     //
//     // We convert all faceZones into baffles during snapping so we can use
//     // a standard mesh motion (except for the mesh checking which for baffles
//     // created from internal faces should check across the baffles). The state
//     // is stored in two variables:
//     //      baffles : pairs of boundary faces
//     //      duplicateFace : from mesh face to its baffle colleague (or -1 for
//     //                      normal faces)
//     // There are three types of faceZones according to the faceType property:
//     //
//     // internal
//     // --------
//     // - baffles: need to be checked across
//     // - duplicateFace: from face to duplicate face. Contains
//     //   all faces on faceZone to prevents merging patch faces.
//     //
//     // baffle
//     // ------
//     // - baffles: no need to be checked across
//     // - duplicateFace: contains all faces on faceZone to prevent
//     //   merging patch faces.
//     //
//     // boundary
//     // --------
//     // - baffles: no need to be checked across. Also points get duplicated
//     //            so will no longer be baffles
//     // - duplicateFace: contains no faces on faceZone since both sides can
//     //   merge faces independently.



//     // faceZones of type internal
//     const labelList internalFaceZones
//     (
//         meshRefiner_.getZones
//         (
//             List<surfaceZonesInfo::faceZoneType>
//             (
//                 1,
//                 surfaceZonesInfo::INTERNAL
//             )
//         )
//     );


//     // Create baffles (pairs of faces that share the same points)
//     // Baffles stored as owner and neighbour face that have been created.
//     {
//         List<labelPair> baffles;
//         labelList originatingFaceZone;
//         meshRefiner_.createZoneBaffles
//         (
//             identity(mesh.faceZones().size()),
//             baffles,
//             originatingFaceZone
//         );
//     }


//     // Duplicate points on faceZones of type boundary
//     meshRefiner_.dupNonManifoldBoundaryPoints();

//     // Extract baffles across internal faceZones (for checking mesh quality
//     // across)
//     labelPairList internalBaffles
//     (
//         meshRefiner_.subsetBaffles
//         (
//             mesh,
//             internalFaceZones,
//             localPointRegion::findDuplicateFacePairs(mesh)
//         )
//     );


//     bool doFeatures = true; //false;
//     label nFeatIter = 1;
//     if (snapParams.nFeatureSnap() > 0)
//     {
//         doFeatures = true;

//         if (!dryRun_)
//         {
//             nFeatIter = snapParams.nFeatureSnap();
//         }

//         Info<< "Snapping to features in " << nFeatIter
//             << " iterations ..." << endl;
//     }

//     // Get the labels of added patches.
//     const labelList adaptPatchIDs(meshRefiner_.meshedPatches());

//     autoPtr<indirectPrimitivePatch> ppPtr
//     (
//         meshRefinement::makePatch
//         (
//             mesh,
//             adaptPatchIDs
//         )
//     );


//     if (debug)
//     {
//         const fileName dir(mesh.time().path());
//         mkDir(dir);
//         OBJstream str
//         (
//             dir
//           / "pp_initial_" + meshRefiner_.timeName() + ".obj"
//         );
//         str.write
//         (
//             ppPtr().localFaces(),
//             ppPtr().localPoints(),
//             false
//         );
//     }


//     // Maximum distance to attract to nearest feature on surface
//     scalarField snapDist(calcSnapDistance(mesh, snapParams, ppPtr()));


//     // Smooth patch points. Equivalent of preSmoothPatch but using mesh
//     // motion solver.
//     vectorField patchDisp
//     (
//         smoothLambdaMuPatchDisplacement
//         (
//             mesh,
//             ppPtr(),
//             ppPtr().localPoints(),
//             internalBaffles
//         )
//     );
//     pointField ppLocalPoints(ppPtr().localPoints()+patchDisp);

//     const bool smoothInternal = true;
//     if (smoothInternal)
//     {
//         // Create pointMesh with the correct patches
//         //const pointMesh& pMesh = makePointMesh
//         //(
//         //    ppPtr(),
//         //    patchConstraints,
//         //    allEdgePatchName,
//         //    allPointPatchName
//         //);
//         const pointMesh& pMesh = pointMesh::New(mesh);

//         autoPtr<displacementMotionSolver> motionPtr
//         (
//             makeMotionSolver
//             (
//                 pMesh,
//                 ppPtr,
//                 snapDict,
//                 adaptPatchIDs
//             )
//         );

//         // Insert as bc to motionSmoother. Note that this takes displacement
//         // relative to points0
//         setDisplacement
//         (
//             ppPtr(),
//             patchDisp,
//             adaptPatchIDs,
//             motionPtr().points0(),
//             motionPtr().pointDisplacement()
//         );
//         // Solve internal displacement
//         tmp<pointField> tnewPoints(motionPtr->newPoints());

//         // Move points
//         mesh.movePoints(tnewPoints);

//         // Update pp for new mesh points. Ok as long as we also update geometry
//         // and wanted displacement (usually zero if mesh motion has succeeded)
//         ppLocalPoints = pointField(mesh.points(), ppPtr().meshPoints());
//         patchDisp -= (ppLocalPoints-ppPtr().localPoints());

//         if (debug & meshRefinement::MESH)
//         {
//             const_cast<Time&>(mesh.time())++;
//             Info<< "Writing smoothed mesh to time "
//                 << meshRefiner_.timeName() << endl;
//             meshRefiner_.write
//             (
//                 meshRefinement::debugType(debug),
//                 meshRefinement::writeType
//                 (
//                     meshRefinement::writeLevel()
//                   | meshRefinement::WRITEMESH
//                 ),
//                 meshRefiner_.timeName()
//             );
//         }
//     }


//     // Construct iterative mesh mover.
//     Info<< "Constructing mesh displacer ..." << endl;
//     Info<< "Using mesh parameters " << motionDict << nl << endl;

//     autoPtr<displacementMotionSolver> motionPtr
//     (
//         makeMotionSolver
//         (
//             pointMesh::New(mesh),
//             ppPtr,
//             snapDict,
//             adaptPatchIDs
//         )
//     );
//     autoPtr<motionSmoother> meshMoverPtr
//     (
//         new motionSmoother
//         (
//             mesh,
//             ppPtr(),
//             adaptPatchIDs,
//             motionPtr().pointDisplacement(),
//             motionDict,
//             dryRun_
//         )
//     );


//     // Check initial mesh
//     Info<< "Checking initial mesh ..." << endl;
//     labelHashSet wrongFaces(mesh.nFaces()/100);
//     motionSmoother::checkMesh(false, mesh, motionDict, wrongFaces, dryRun_);
//     const label nInitErrors = returnReduce
//     (
//         wrongFaces.size(),
//         sumOp<label>()
//     );

//     Info<< "Detected " << nInitErrors << " illegal faces"
//         << " (concave, zero area or negative cell pyramid volume)"
//         << endl;


//     Info<< "Checked initial mesh in = "
//         << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;



//     //- Only if in feature attraction mode:
//     // Point on nearest feature
//     vectorField patchFeaturePoint(ppPtr().nPoints(), Zero);
//     // Constraints at feature
//     List<pointConstraint> patchConstraints(ppPtr().nPoints());

//     for (label iter = 0; iter < nFeatIter; iter++)
//     {
//         Info<< nl
//             << "Morph iteration " << iter << nl
//             << "-----------------" << endl;


//         // Calculate displacement at every patch point if we need it:
//         // - if automatic near-surface detection
//         // - if face splitting active
//         pointField nearestPoint(ppPtr().nPoints(), vector::max);
//         vectorField nearestNormal(ppPtr().nPoints(), Zero);


//         const bool strictRegionSnap
//         (
//             iter < nFeatIter/2
//           ? snapParams.strictRegionSnap()
//           : Switch(true)
//         );

//         vectorField disp = calcNearestSurface
//         (
//             strictRegionSnap,               // attract points to region only
//             meshRefiner_,
//             globalToMasterPatch_,           // for if strictRegionSnap
//             globalToSlavePatch_,            // for if strictRegionSnap
//             ppPtr(),
//             ppLocalPoints,
//             snapDist,                       // max snap distance

//             nearestPoint,
//             nearestNormal
//         );

//         // Override displacement at thin gaps
//         if (snapParams.detectNearSurfacesSnap())
//         {
//             detectNearSurfaces
//             (
//                 Foam::cos(degToRad(planarAngle)),// planar cos for gaps
//                 ppPtr(),
//                 ppLocalPoints,
//                 nearestPoint,   // surfacepoint from nearest test
//                 nearestNormal,  // surfacenormal from nearest test

//                 disp
//             );
//         }

//         // Override displacement with feature edge attempt
//         if (doFeatures)
//         {
//             //- Any faces to split
//             DynamicList<label> splitFaces;
//             //- Indices in face to split across
//             DynamicList<labelPair> splits;
//             //- Patches for split face
//             DynamicList<labelPair> splitPatches;

//             // Offset to project to nearest feature. Use in combination with
//             // patchConstraints.
//             vectorField patchAttraction;

//             disp = calcNearestSurfaceFeature
//             (
//                 snapParams,
//                 false,              // no alignMeshEdges
//                 true,               // check >=3 patch points

//                 iter,
//                 featureCos,
//                 scalar(iter+1)/nFeatIter,

//                 snapDist,
//                 disp,               // nearest surface
//                 nearestNormal,
//                 ppPtr(),
//                 ppLocalPoints,

//                 patchAttraction,    // offset wrt ppLocalPoints to nearest
//                                     // feature edge/point
//                 patchConstraints,   // feature type + constraint

//                 splitFaces,
//                 splits,
//                 splitPatches
//             );


//             // Freeze points on exposed points/faces
//             freezeExposedPoints
//             (
//                 meshRefiner_,
//                 "frozenFaces",      // faceZone name
//                 "frozenPoints",     // pointZone name
//                 ppPtr(),
//                 disp
//             );

//             patchFeaturePoint = ppLocalPoints+patchAttraction;

//             if (debug)
//             {
//                 OBJstream str
//                 (
//                     mesh.time().path()
//                   / "calcNearestSurfaceFeature"
//                   + meshRefiner_.timeName()
//                   + ".obj"
//                 );
//                 forAll(ppLocalPoints, pointi)
//                 {
//                     const point& pt = ppLocalPoints[pointi];
//                     str.write(linePointRef(pt, pt+disp[pointi]));
//                 }
//             }



//             // Split any faces:
//             // - does not move/add any points
//             // - tries to move new faces to correct patch

//             if (returnReduce(splitFaces.size(), sumOp<label>()))
//             {
//                 polyTopoChange meshMod(mesh);

//                 // Insert the mesh changes
//                 meshRefiner_.doSplitFaces
//                 (
//                     splitFaces,
//                     splits,
//                     splitPatches,
//                     meshMod
//                 );

//                 // Save old meshPoints before changing mesh
//                 const Map<label> oldMeshPointMap(ppPtr->meshPointMap());
//                 Pout<< "old pp points:" << ppPtr->nPoints()
//                     << " oldMeshPointMap:" << oldMeshPointMap.size()
//                     << endl;

//                 // Remove any unnecessary fields
//                 meshMoverPtr.clear();
//                 motionPtr.clear();
//                 ppPtr.clear();
//                 mesh.clearOut();
//                 mesh.moving(false);

//                 // Change the mesh (no inflation)
//                 autoPtr<mapPolyMesh> mapPtr = meshMod.changeMesh(mesh, false);
//                 mapPolyMesh& map = *mapPtr;

//                 // Update fields
//                 mesh.updateMesh(map);

//                 // Move mesh (since morphing might not do this)
//                 if (map.hasMotionPoints())
//                 {
//                     mesh.movePoints(map.preMotionPoints());
//                 }
//                 else
//                 {
//                     mesh.clearOut();
//                 }

//                 // Reset the instance for if in overwrite mode
//                 mesh.setInstance(meshRefiner_.timeName());
//                 meshRefiner_.setInstance(mesh.facesInstance());

//                 // Update intersections on split faces
//                 {
//                     DynamicList<label> changedFaces(splitFaces.size());
//                     Map<label> splitFacesMap(splitFaces.size());
//                     forAll(splitFaces, i)
//                     {
//                         splitFacesMap.insert(splitFaces[i], i);
//                     }
//                     forAll(map.faceMap(), facei)
//                     {
//                         if (splitFacesMap.find(map.faceMap()[facei]))
//                         {
//                             changedFaces.append(facei);
//                         }
//                     }
//                     // Update intersections on changed faces
//                     meshRefiner_.updateMesh
//                     (
//                         map,
//                         meshRefiner_.growFaceCellFace(changedFaces)
//                     );
//                 }

//                 if (debug&meshRefinement::MESH)
//                 {
//                     const_cast<Time&>(mesh.time())++;
//                     Info<< "Writing split-faces mesh to time "
//                         << meshRefiner_.timeName() << endl;
//                     meshRefiner_.write
//                     (
//                         meshRefinement::debugType(debug),
//                         meshRefinement::writeType
//                         (
//                             meshRefinement::writeLevel()
//                           | meshRefinement::WRITEMESH
//                         ),
//                         mesh.time().path()/meshRefiner_.timeName()
//                     );
//                 }



//                 // Update local mesh data
//                 // ~~~~~~~~~~~~~~~~~~~~~~

//                 Info<< "Updating for face-splitting" << endl;

//                 // baffles
//                 forAll(internalBaffles, i)
//                 {
//                     labelPair& baffle = internalBaffles[i];
//                     baffle.first() = map.reverseFaceMap()[baffle.first()];
//                     baffle.second() = map.reverseFaceMap()[baffle.second()];
//                 }

//                 // re-do patch (since faces might have been split)
//                 ppPtr = meshRefinement::makePatch(mesh, adaptPatchIDs);
//                 motionPtr = makeMotionSolver
//                 (
//                     pointMesh::New(mesh),
//                     ppPtr,
//                     snapDict,
//                     adaptPatchIDs
//                 );
//                 meshMoverPtr.reset
//                 (
//                     new motionSmoother
//                     (
//                         mesh,
//                         ppPtr(),
//                         adaptPatchIDs,
//                         motionPtr().pointDisplacement(),
//                         motionDict,
//                         dryRun_
//                     )
//                 );

//                 const auto& mp = ppPtr->meshPoints();
//                 // pointMap (new-to-old) for pp points. Note: no points changed
//                 // but local point ordering might have changed since faces
//                 // split.
//                 labelList ppMap(mp.size());
//                 forAll(mp, i)
//                 {
//                     ppMap[i] = oldMeshPointMap[mp[i]];
//                 }
//                 // patchDisp
//                 meshRefinement::updateList(ppMap, vector::zero, patchDisp);
//                 // snapDist
//                 meshRefinement::updateList(ppMap, scalar(0), snapDist);
//                 // patchFeaturePoint
//                 meshRefinement::updateList
//                 (
//                     ppMap,
//                     vector::zero,
//                     patchFeaturePoint
//                 );
//                 // patchConstraints
//                 meshRefinement::updateList
//                 (
//                     ppMap,
//                     pointConstraint(),
//                     patchConstraints
//                 );
// //                // nearestPoint
// //                meshRefinement::updateList(ppMap, vector::zero, nearestPoint);
// //                // nearestNormal
// //                meshRefinement::updateList(ppMap, vector::zero, nearestNormal);
//                 // disp
//                 meshRefinement::updateList(ppMap, vector::zero, disp);
//                 // ppLocalPoints
//                 meshRefinement::updateList(ppMap, vector::zero, ppLocalPoints);

//                 Info<< "DONE Updating for face-splitting" << endl;
//             }
//         }


//         // Attract/slide the mesh a bit
//         {
//             //// Create pointMesh with the correct patches
//             //// - points on single polyPatches stay as is
//             //// - points on two polyPatches go to allEdgePatchName
//             //// - points on >two polyPatches go to allPointPatchName
//             //const pointMesh& pMesh = makePointMesh
//             //(
//             //    ppPtr(),
//             //    patchConstraints,
//             //    allEdgePatchName,
//             //    allPointPatchName
//             //);
//             //const pointMesh& pMesh = pointMesh::New(mesh);
//             //
//             //autoPtr<displacementMotionSolver> motionPtr
//             //(
//             //    makeMotionSolver
//             //    (
//             //        pMesh,
//             //        snapDict,
//             //        adaptPatchIDs
//             //        //patchConstraints
//             //    )
//             //);

//             // Insert as bc to motionSmoother. Note that this takes displacement
//             // relative to points0
//             setDisplacement
//             (
//                 ppPtr(),
//                 disp,
//                 adaptPatchIDs,
//                 motionPtr().points0(),
//                 motionPtr().pointDisplacement()
//             );

//             // Solve internal displacement
//             tmp<pointField> tnewPoints(motionPtr->newPoints());

//             // Move points
//             if (false)
//             {
//                 // 1. Directly move points
//                 mesh.movePoints(tnewPoints);
//                 // Optional? Only geometry used is ppLocalPoints which we keep
//                 // and update 'by hand'.
//                 //ppPtr().movePoints(tnewPoints);
//             }
//             else
//             {
//                 // 2. Use motionSmoother
//                 // Set initial distribution of displacement field (on patches)
//                 // from patchDisp and make displacement consistent with b.c.
//                 // on displacement pointVectorField.
//                 const vectorField newDisp(tnewPoints()-meshMoverPtr().oldPoints());

//                 meshMoverPtr().displacement().vectorField::operator=(newDisp);
//                 meshMoverPtr().setDisplacement(disp);

//                 // Apply internal displacement to mesh.
//                 const bool meshOk = scaleMesh
//                 (
//                     snapParams,
//                     nInitErrors,
//                     internalBaffles,
//                     meshMoverPtr()
//                 );

//                 if (!meshOk)
//                 {
//                     WarningInFunction
//                         << "Did not successfully snap mesh."
//                         << " Continuing to snap to resolve easy" << nl
//                         << "    surfaces but the"
//                         << " resulting mesh will not satisfy your quality"
//                         << " constraints" << nl << endl;
//                 }

//                 // Use current mesh as base mesh
//                 meshMoverPtr().correct();
//             }

//             // Update pp for new mesh points. Ok as long as we also update
//             // geometry and wanted displacement (usually zero if mesh motion
//             // has succeeded)
//             ppLocalPoints = pointField(mesh.points(), ppPtr().meshPoints());
//             patchDisp -= (ppLocalPoints-ppPtr().localPoints());

//             if (debug & meshRefinement::MESH)
//             {
//                 const_cast<Time&>(mesh.time())++;
//                 Info<< "Writing partially moved mesh to time "
//                     << meshRefiner_.timeName() << endl;
//                 meshRefiner_.write
//                 (
//                     meshRefinement::debugType(debug),
//                     meshRefinement::writeType
//                     (
//                         meshRefinement::writeLevel()
//                       | meshRefinement::WRITEMESH
//                     ),
//                     meshRefiner_.timeName()
//                 );
//             }
//         }

//         // Split problematic cells
//         {
//             Info<< nl << "Checking moved mesh ..." << endl;
//             faceSet wrongFaces(mesh, "wrongFaces", mesh.nFaces()/1000);
//             motionSmoother::checkMesh
//             (
//                 false,
//                 mesh,
//                 motionDict,
//                 identity(mesh.nFaces()),
//                 internalBaffles,
//                 wrongFaces,
//                 false           // dryRun_
//             );
//             const label nWrong(returnReduce(wrongFaces.size(), sumOp<label>()));
//             Info<< "Detected " << nWrong
//                 << " illegal faces"
//                 << " (concave, zero area or negative cell pyramid volume)"
//                 << endl;

//             //if (nWrong)
//             if (false)
//             {
//                 bitSet decomposeCell(mesh.nCells());
//                 for (const label facei : wrongFaces)
//                 {
//                     const label own = mesh.faceOwner()[facei];
//                     if (!tetMatcher::test(mesh, own))
//                     {
//                         decomposeCell.set(own);
//                     }
//                     if (mesh.isInternalFace(facei))
//                     {
//                         const label nei = mesh.faceNeighbour()[facei];
//                         if (!tetMatcher::test(mesh, nei))
//                         {
//                             decomposeCell.set(nei);
//                         }
//                     }
//                 }
//                 Pout<< "spliyyinG :" << decomposeCell.count() << " cells"
//                     << endl;

//                 tetDecomposer tetDecomp(mesh);
//                 polyTopoChange meshMod(mesh);
//                 tetDecomp.setRefinement
//                 (
//                     tetDecomposer::FACE_CENTRE_TRIS,
//                     decomposeCell,
//                     meshMod
//                 );

//                 // Save old meshPoints before changing mesh
//                 const Map<label> oldMeshPointMap(ppPtr->meshPointMap());

//                 Pout<< "old pp points:" << ppPtr->nPoints()
//                     << " oldMeshPointMap:" << oldMeshPointMap.size()
//                     << endl;

//                 // Remove any unnecessary fields
//                 meshMoverPtr.clear();
//                 motionPtr.clear();
//                 ppPtr.clear();
//                 mesh.clearOut();
//                 mesh.moving(false);

//                 // Change the mesh (no inflation)
//                 autoPtr<mapPolyMesh> mapPtr = meshMod.changeMesh(mesh, false);
//                 mapPolyMesh& map = *mapPtr;

//                 // Update fields
//                 mesh.updateMesh(map);

//                 // Move mesh (since morphing does not do this)
//                 if (map.hasMotionPoints())
//                 {
//                     mesh.movePoints(map.preMotionPoints());
//                 }
//                 else
//                 {
//                     mesh.clearOut();
//                 }

//                 // Reset the instance for if in overwrite mode
//                 mesh.setInstance(meshRefiner_.timeName());
//                 meshRefiner_.setInstance(mesh.facesInstance());

//                 //- Update numbering on tet-decomposition engine
//                 tetDecomp.updateMesh(map);

//                 bitSet isChangedFace(mesh.nFaces());
//                 forAll(map.cellMap(), celli)
//                 {
//                     if (decomposeCell[map.cellMap()[celli]])
//                     {
//                         isChangedFace.set(mesh.cells()[celli]);
//                     }
//                 }
//                 syncTools::syncFaceList
//                 (
//                     mesh,
//                     isChangedFace,
//                     orEqOp<unsigned int>()
//                 );

//                 Pout<< "isChangedFace :" << decomposeCell.count() << " faces"
//                     << endl;


//                 // Update intersection info
//                 meshRefiner_.updateMesh(map, isChangedFace.toc());


//                 if (debug&meshRefinement::MESH)
//                 {
//                     const_cast<Time&>(mesh.time())++;
//                     Info<< "Writing tet-decomp mesh to time "
//                         << meshRefiner_.timeName() << endl;
//                     meshRefiner_.write
//                     (
//                         meshRefinement::debugType(debug),
//                         meshRefinement::writeType
//                         (
//                             meshRefinement::writeLevel()
//                           | meshRefinement::WRITEMESH
//                         ),
//                         mesh.time().path()/meshRefiner_.timeName()
//                     );
//                 }


//                 // Update local mesh data
//                 // ~~~~~~~~~~~~~~~~~~~~~~

//                 Info<< "Updating for tet-decomp" << endl;

//                 // baffles
//                 forAll(internalBaffles, i)
//                 {
//                     labelPair& baffle = internalBaffles[i];
//                     baffle.first() = map.reverseFaceMap()[baffle.first()];
//                     baffle.second() = map.reverseFaceMap()[baffle.second()];
//                 }

//                 // re-do patch (since faces might have been split)
//                 ppPtr = meshRefinement::makePatch(mesh, adaptPatchIDs);
//     Pout<< "new pp points:" << ppPtr->nPoints()
//         << " new meshPointMap:" << ppPtr->meshPointMap().size()
//         << endl;
//                 const auto& mp = ppPtr->meshPoints();
//                 // pointMap (new-to-old) for pp points. Might have new
//                 // face-centre points - these get mapped from any point on
//                 // the originating face.
//                 labelList ppMap(mp.size(), -1);
//                 forAll(mp, i)
//                 {
//                     const label oldMeshPointi = map.pointMap()[mp[i]];
//                     const auto mpFnd = oldMeshPointMap.find(oldMeshPointi);
//                     if (mpFnd)
//                     {
//                         ppMap[i] = mpFnd();
//                     }
//                 }

//                 // patchDisp
//                 meshRefinement::updateList(ppMap, vector::zero, patchDisp);
//                 // snapDist
//                 meshRefinement::updateList(ppMap, scalar(0), snapDist);
//                 // patchFeaturePoint
//                 meshRefinement::updateList
//                 (
//                     ppMap,
//                     vector::zero,
//                     patchFeaturePoint
//                 );
//                 // patchConstraints
//                 meshRefinement::updateList
//                 (
//                     ppMap,
//                     pointConstraint(),
//                     patchConstraints
//                 );
// //                // nearestPoint
// //                meshRefinement::updateList(ppMap, vector::zero, nearestPoint);
// //                // nearestNormal
// //                meshRefinement::updateList(ppMap, vector::zero, nearestNormal);
// //                // disp
// //                meshRefinement::updateList(ppMap, vector::zero, disp);

//                 // Maximum distance to attract to nearest feature on surface
//                 snapDist = calcSnapDistance(mesh, snapParams, ppPtr());
//                 ppLocalPoints = pointField(mesh.points(), ppPtr().meshPoints());

//                 Info<< "DONE Updating for tet-decomp" << endl;
//             }
//         }
//     }

// //XXXXXX
// /*
//     {
//         // Introduce single layer of cells. Straight from snappyLayerDriver

//         // Global face indices engine
//         const globalIndex globalFaces(mesh.nFaces());

//         // Determine extrudePatch.edgeFaces in global numbering (so across
//         // coupled patches). This is used only to string up edges
//         // between coupled
//         // faces (all edges between same (global)face indices get extruded).
//         labelListList edgeGlobalFaces
//         (
//             addPatchCellLayer::globalEdgeFaces
//             (
//                 mesh,
//                 globalFaces,
//                 ppPtr()
//             )
//         );

//         // Determine patches for extruded boundary edges. Calculates if any
//         // additional processor patches need to be constructed.

//         labelList edgePatchID;
//         labelList edgeZoneID;
//         boolList edgeFlip;
//         labelList inflateFaceID;
//         snappyLayerDriver::determineSidePatches
//         (
//             meshRefiner_,
//             globalFaces,
//             edgeGlobalFaces,
//             ppPtr(),

//             edgePatchID,
//             edgeZoneID,
//             edgeFlip,
//             inflateFaceID
//         );


//         // Mesh topo change engine. Insert current mesh.
//         polyTopoChange meshMod(mesh);

//         // Layer mesh modifier
//         addPatchCellLayer addLayer(mesh);

//         // Extrude very thin layer of cells
//         pointField extrusion(PatchTools::pointNormals(mesh, ppPtr()));
//         const tmp<scalarField> thickness
//         (
//             wantedThickness
//             (
//                 ppPtr(),
//                 1e-3        // cellSizeFraction
//             )
//         );
//         extrusion *= thickness;

//         // Add topo regardless of whether extrudeStatus is extruderemove.
//         // Not add layer if patchDisp is zero.

//         const label ppNFaces = ppPtr().size();
//         const label ppNPoints = ppPtr().nPoints();

//         addLayer.setRefinement
//         (
//             globalFaces,
//             edgeGlobalFaces,

//             scalarField(ppNPoints, 1),  // expansion ratio
//             ppPtr(),
//             bitSet(ppPtr().size()),     // no flip

//             edgePatchID,    // boundary patch for extruded boundary edges
//             edgeZoneID,     // zone for extruded edges
//             edgeFlip,
//             inflateFaceID,


//             labelList(0),   // exposed patchIDs, not used for adding layers
//             labelList(ppNFaces, 1),
//             labelList(ppNPoints, 1),
//             extrusion,  //patchAttraction,

//             meshMod
//         );

//         // Save old meshPoints before changing mesh
//         const Map<label> oldMeshPointMap(ppPtr->meshPointMap());
//         Pout<< "old pp points:" << ppPtr->nPoints()
//             << " oldMeshPointMap:" << oldMeshPointMap.size()
//             << endl;

//         // Remove any unnecessary fields
//         meshMoverPtr.clear();
//         motionPtr.clear();
//         ppPtr.clear();
//         mesh.clearOut();
//         mesh.moving(false);

//         // Apply the stored topo changes to the current mesh.
//         autoPtr<mapPolyMesh> mapPtr = meshMod.changeMesh(mesh, false);
//         mapPolyMesh& map = *mapPtr;

//         // Update fields
//         mesh.updateMesh(map);

//         // Move mesh (since morphing does not do this)
//         if (map.hasMotionPoints())
//         {
//             mesh.movePoints(map.preMotionPoints());
//         }
//         else
//         {
//             // Hack to remove meshPhi - mapped incorrectly. TBD.
//             mesh.clearOut();
//         }

//         // Reset the instance for if in overwrite mode
//         mesh.setInstance(meshRefiner_.timeName());

//         // Re-do the patch
//         ppPtr = meshRefinement::makePatch(mesh, adaptPatchIDs);
//     Pout<< "new pp points:" << ppPtr->nPoints()
//         << " new meshPointMap:" << ppPtr->meshPointMap().size()
//         << endl;
//         // Map the old-patch-point data to the new patch points
//         // pointMap (new-to-old) for new pp points
//         labelList ppMap(ppPtr->nPoints(), -1);
//         {
//             const labelListList& oldAddedPoints = addLayer.addedPoints();
//             forAll(oldAddedPoints, oldPatchPointi)
//             {
//                 const label oldPointi = oldAddedPoints[oldPatchPointi].last();
//                 const label newPointi = map.reversePointMap()[oldPointi];
//                 const label newPatchPointi = ppPtr().meshPointMap()[newPointi];

//                 ppMap[newPatchPointi] = oldPatchPointi;
//             }
//         }

//         // Update attraction
//         meshRefinement::updateList(ppMap, vector::zero, patchDisp);
//         // snapDist
//         meshRefinement::updateList(ppMap, scalar(0), snapDist);
//         // patchFeaturePoint
//         meshRefinement::updateList
//         (
//             ppMap,
//             vector::zero,
//             patchFeaturePoint
//         );
//         // patchConstraints
//         meshRefinement::updateList
//         (
//             ppMap,
//             pointConstraint(),
//             patchConstraints
//         );
//         // ppLocalPoints
//         //meshRefinement::updateList(ppMap, vector::zero, ppLocalPoints);
//         ppLocalPoints = pointField(mesh.points(), ppPtr().meshPoints());


//         // Update numbering on layer
//         addLayer.updateMesh
//         (
//             map,
//             identity(ppNFaces),
//             identity(ppNPoints)
//         );

//         // Update intersections
//         const labelListList addedCells(addLayer.addedCells());
//         bitSet isChangedFace(mesh.nFaces());
//         for (const labelList& faceToCells : addedCells)
//         {
//             for (const label celli : faceToCells)
//             {
//                 isChangedFace.set(mesh.cells()[celli]);
//             }
//         }
//         meshRefiner_.updateMesh(map, isChangedFace.toc());


//         if (debug & meshRefinement::MESH)
//         {
//             const_cast<Time&>(mesh.time())++;
//             Info<< "Writing mesh-with-layer to time "
//                 << meshRefiner_.timeName() << endl;
//             meshRefiner_.write
//             (
//                 meshRefinement::debugType(debug),
//                 meshRefinement::writeType
//                 (
//                     meshRefinement::writeLevel()
//                   | meshRefinement::WRITEMESH
//                 ),
//                 meshRefiner_.timeName()
//             );
//         }
//         if (debug)
//         {
//             OBJstream str
//             (
//                 mesh.time().path()
//               / "new_projection"
//               + meshRefiner_.timeName()
//               + ".obj"
//             );
//             forAll(ppLocalPoints, pointi)
//             {
//                 const point& pt = ppLocalPoints[pointi];
//                 str.write(linePointRef(pt, patchFeaturePoint[pointi]));
//             }
//             Pout<< "** writing mapped attraction to " << str.name() << endl;
//         }
//     }
// }


void Foam::snappySnapDriver::doSnapBufferLayers
(
    const dictionary& snapDict,
    const dictionary& motionDict,
    const meshRefinement::FaceMergeType mergeType,
    const scalar featureCos,
    const scalar planarAngle,
    const snapParameters& snapParams
)
{
    addProfiling(snap, "snappyHexMesh::snap");
    fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl
        << "Morphing phase" << nl
        << "--------------" << nl
        << endl;


    // Name of pointPatch for any points attracted to surface
    const word allSurfPatchName("surfacePoints");
    // Name of pointPatch for any points attracted to feature edges
    const word allEdgePatchName("boundaryEdges");
    // Name of pointPatch for any points attracted to feature points
    const word allPointPatchName("boundaryPoints");


    // Dummy value
    const vector GGG(vector::uniform(GREAT));


    // faceZone handling
    // ~~~~~~~~~~~~~~~~~
    //
    // We convert all faceZones into baffles during snapping so we can use
    // a standard mesh motion (except for the mesh checking which for baffles
    // created from internal faces should check across the baffles). The state
    // is stored in two variables:
    //      baffles : pairs of boundary faces
    //      duplicateFace : from mesh face to its baffle colleague (or -1 for
    //                      normal faces)
    // There are three types of faceZones according to the faceType property:
    //
    // internal
    // --------
    // - baffles: need to be checked across
    // - duplicateFace: from face to duplicate face. Contains
    //   all faces on faceZone to prevents merging patch faces.
    //
    // baffle
    // ------
    // - baffles: no need to be checked across
    // - duplicateFace: contains all faces on faceZone to prevent
    //   merging patch faces.
    //
    // boundary
    // --------
    // - baffles: no need to be checked across. Also points get duplicated
    //            so will no longer be baffles
    // - duplicateFace: contains no faces on faceZone since both sides can
    //   merge faces independently.

    // New faceZone algorithm:
    // - keep faceZones as internal faces
    // - extract as primitivePatch
    // - smooth primitivePatch
    // - determine correspondence to surface
    // - optionally split (possibly still internal) mesh faces
    // - add buffer layers:
    //      - note: since cannot add to multiple sides we have to split mesh:
    //      - convert internal faces into baffles
    //      - duplicate points
    //      - determine extrusion direction
    //      - extrude
    //      - merge points
    //      - merge faces
    // - have mesh motion smoother do all (using pointMesh patches, not
    //   polyMesh patches)

    // Get labels of patches where optional buffer layers are added
    DynamicList<label> bufPatchIDs;
    DynamicList<label> bufZoneIDs;
    if
    (
        meshRefiner_.meshType() == meshRefinement::CASTELLATEDBUFFERLAYER
     || meshRefiner_.meshType() == meshRefinement::CASTELLATEDBUFFERLAYER2
    )
    {
        bufPatchIDs.setCapacity(globalToMasterPatch_.size());

        const auto& addLayers = meshRefiner_.surfaces().addBufferLayers();
        // Normal patches
        forAll(globalToMasterPatch_, globalRegioni)
        {
            if (addLayers[globalRegioni])
            {
                bufPatchIDs.push_uniq(globalToMasterPatch_[globalRegioni]);
                bufPatchIDs.push_uniq(globalToSlavePatch_[globalRegioni]);
            }
        }
        // Temporary patches from faceZones
        for (const auto& fz : mesh.faceZones())
        {
            label mpI, spI;
            surfaceZonesInfo::faceZoneType type;
            if (meshRefiner_.getFaceZoneInfo(fz.name(), mpI, spI, type))
            {
                bufPatchIDs.push_uniq(mpI);
                bufPatchIDs.push_uniq(spI);
                bufZoneIDs.push_uniq(fz.index());
            }
        }
    }


    // faceZones of type internal
    const labelList internalFaceZones
    (
        meshRefiner_.getZones
        (
            List<surfaceZonesInfo::faceZoneType>
            (
                1,
                surfaceZonesInfo::INTERNAL
            )
        )
    );
    const labelList baffleFaceZones
    (
        meshRefiner_.getZones
        (
            List<surfaceZonesInfo::faceZoneType>
            (
                1,
                surfaceZonesInfo::BAFFLE
            )
        )
    );
    const labelList boundaryFaceZones
    (
        meshRefiner_.getZones
        (
            List<surfaceZonesInfo::faceZoneType>
            (
                1,
                surfaceZonesInfo::BOUNDARY
            )
        )
    );
    if (boundaryFaceZones.size())
    {
        FatalErrorInFunction
            << "In snap mode "
            << meshRefinement::MeshTypeNames[meshRefiner_.meshType()]
            << " faceZone cannnot be of type \"boundary\"."
            << " Use \"baffle\" instead"
            << exit(FatalError);
    }


    // Get the labels of added patches.
    labelList adaptPatchIDs(meshRefiner_.meshedPatches());

    // Create addressing engine. Note: refers to mesh.points() for e.g.
    // faceCentres, faceNormals so don't use these but calculate patch
    // geometry on-the-fly. Alternatively have whole
    // copy of mesh.points() and modify those.
    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            mesh.points(),
            adaptPatchIDs,
            bufZoneIDs
        )
    );
    pointField ppLocalPoints(ppPtr().localPoints());

    // Maximum distance to attract to nearest feature on surface
    scalarField snapDist(calcSnapDistance(mesh, snapParams, ppPtr()));


    // Smooth patch points. Equivalent of preSmoothPatch. Note that ppPtr is
    // still non-baffled.
    vectorField patchDisp
    (
        smoothLambdaMuPatchDisplacement
        (
            mesh,
            ppPtr(),
            ppLocalPoints,
            List<labelPair>::null() //internalBaffles // Not used!
        )
    );
    ppLocalPoints += patchDisp;

    // Calculate current pointNormals
    const vectorField ppPointNormals
    (
        PatchTools::pointNormals
        (
            mesh,
            ppPtr(),
            ppLocalPoints
        )
    );

    // Save smooth pp points for debugging
    pointField oldPPLocalPoints = ppLocalPoints;




    // Determine correspondence of mesh outside to surface
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Info<< "Determining mesh boundary to surface correspondence ..."
        << endl;

    //- Mesh-to-surface mapping

    // - for every mesh vertex determine nearest surface location. Only accept
    //   if normal is similar
    // - for every surface feature point:
    //      - find nearby mesh vertices
    //      - ignore if its point order is < surface point order
    //      - ignore if surrounding mesh faces are not attracted to surrounding
    //        surface planes
    //      - if all fails just use nearest
    // - determine faces containing points that attract to different surface
    //   topology.
    // - try to walk from feature point to feature point across these faces
    // - if walking across faces split the face.


    // patch face to nearest surface
    labelList patchFaceToSurfaceRegion;
    // patch face to diagonal split info
    Map<edge> patchFaceToDiagonal;
    // Patch point to surface 'face'
    labelList patchToFace(ppPtr().nPoints(), -1);
    // Point on nearest surface
    pointField patchToNearest(ppPtr().nPoints(), vector::uniform(GREAT));
    // Normal on nearest surface
    vectorField patchToNormal(ppPtr().nPoints(), vector::uniform(GREAT));
    // Index of nearest surface
    labelList patchToSurface(ppPtr().nPoints(), -1);
    // Region on nearest surface
    labelList patchToRegion(ppPtr().nPoints(), -1);
    // Surface 'feature' index
    labelList patchToFeature(ppPtr().nPoints(), -1);
    // Surface feature point (if patch point is attracted)
    labelList patchToPoint(ppPtr().nPoints(), -1);
    // Surface edge 'region' (if patch point is attracted to it)
    labelList patchToEdge(ppPtr().nPoints(), -1);
    // Inverse of patchToEdge : from surface feat edge back to pp
    Map<label> surfFeatEdgeRegionToStart;
    // Constraints at surface (1 = surface, 2 = featedge, 3 = featpoint)
    List<pointConstraint> patchConstraints(ppPtr().nPoints());

    matchToSurfaceMerged
    (
        (snapParams.nFeatureSnap() > 0),    // do geom/region features
        featureCos,

        ppPtr(),                // addressing only
        ppLocalPoints,          // points (modified during matching!)
        ppPointNormals,         // normal

       // Per string index a starting pp point
        surfFeatEdgeRegionToStart,


        //- Per patch face

        // Region on nearest surface
        patchFaceToSurfaceRegion,
        // Diagonal split info
        patchFaceToDiagonal,

        //- Per patch point

        // Point on nearest surface
        patchToNearest,
        // Normal on nearest surface
        patchToNormal,
        // Index of nearest surface
        patchToSurface,
        // Region on nearest surface
        patchToRegion,
        // Surface 'face' index
        patchToFace,
        // Surface 'edge' index
        patchToEdge,
        // Surface 'feature' index
        patchToFeature,
        // Surface feature point (if patch point is attracted)
        patchToPoint,
        // Constraints at surface (1 = surface, 2 = featedge, 3 = featpoint)
        patchConstraints
    );


    // Dump the attraction
    if (debug & meshRefinement::MESH)
    {
        writeConstraints
        (
            "pp_analysis",
            oldPPLocalPoints,  //ppPtr().localPoints()
            patchToNearest,
            patchConstraints
        );
    }



    // Split mesh faces diagonally to have continuous feat-edge walk
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Info<< "Splitting boundary feature faces ..." << endl;
    {
        // Topo changes container
        polyTopoChange meshMod(mesh);

        // Work arrays
        DynamicList<label> zones;
        DynamicList<bool> flips;

        forAll(ppPtr().localFaces(), facei)
        {
            const face& localF = ppPtr().localFaces()[facei];
            const label meshFacei = ppPtr().addressing()[facei];
            const label nei =
            (
                mesh.isInternalFace(meshFacei)
              ? mesh.faceNeighbour()[meshFacei]
              : -1
            );

            const auto& iter = patchFaceToDiagonal.find(facei);
            if (iter.good())
            {
                const labelList cuts(iter());
                // cuts is a diagonal/fan split
                // - say split from 1 to 3
                // - face0 : f[1], f[2], f[3]
                // - face1 : f[3], f[0], f[1]
                DynamicList<face> faces;
                splitFace(mesh.faces()[meshFacei], cuts, faces);

                meshMod.faceZones(meshFacei, zones, flips);
                {
                    // Use cuts[0]. Find the non-split
                    label patchi = -1;
                    if (!mesh.isInternalFace(meshFacei))
                    {
                        const label regioni = findSplitFaceRegion
                        (
                            true,
                            ppPtr(),
                            patchFaceToSurfaceRegion,
                            patchConstraints,
                            patchToRegion,
                            facei,
                            cuts[0],
                            cuts[1]
                        );
                        patchi = globalToMasterPatch_[regioni];
                    }
                    meshMod.modifyFace
                    (
                        faces[0],
                        meshFacei,
                        mesh.faceOwner()[meshFacei],
                        nei,
                        false,          //flipFaceFlux
                        patchi,         //patchID
                        zones,          //zoneID
                        flips           //zoneFlip
                    );
                }
                for (label i = 1; i < faces.size(); i++)
                {
                    label patchi = -1;
                    if (!mesh.isInternalFace(meshFacei))
                    {
                        const label regioni = findSplitFaceRegion
                        (
                            true,
                            ppPtr(),
                            patchFaceToSurfaceRegion,
                            patchConstraints,
                            patchToRegion,
                            facei,
                            cuts[i],
                            cuts[cuts.fcIndex(i)]
                        );
                        patchi = globalToMasterPatch_[regioni];
                    }
                    meshMod.addFace
                    (
                        faces[i],
                        mesh.faceOwner()[meshFacei],
                        nei,
                        -1,             //masterPointID
                        -1,             //masterEdgeID
                        meshFacei,      //masterFaceID
                        false,          //flipFaceFlux
                        patchi,         //patchID
                        zones,          //zoneID
                        flips           //zoneFlip
                    );
                }
            }
            else if (!mesh.isInternalFace(meshFacei))
            {
                // Repatch boundary face?
                const label regioni = findSplitFaceRegion
                (
                    true,
                    ppPtr(),
                    patchFaceToSurfaceRegion,
                    patchConstraints,
                    patchToRegion,
                    facei,
                    0,
                    localF.size()-1
                );
                const label patchi = globalToMasterPatch_[regioni];

                meshMod.faceZones(meshFacei, zones, flips);
                meshMod.modifyFace
                (
                    mesh.faces()[meshFacei],
                    meshFacei,
                    mesh.faceOwner()[meshFacei],
                    nei,
                    false,          //flipFaceFlux
                    patchi,         //patchID
                    zones,          //zoneID
                    flips           //zoneFlip
                );
            }
        }

        autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

        mesh.updateMesh(morphMap);

        // Move mesh (since morphing might not do this)
        if (morphMap().hasMotionPoints())
        {
            mesh.movePoints(morphMap().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes. No other way to do this?
            mesh.clearOut();
        }

        // Update intersections. Recalculate intersections for exposed faces.
        const labelList newSplitFaces = renumber
        (
            morphMap().reverseFaceMap(),
            ppPtr().addressing()
        );
        meshRefiner_.updateMesh(morphMap(), newSplitFaces);

        // Update local data (though only split faces; no new points/edges)
        {
            // baffles are in mesh face indices
            // const auto& reverseFaceMap = morphMap().reverseFaceMap();
            // for (auto& baffle : baffles)
            // {
            //     baffle.first() = reverseFaceMap[baffle.first()];
            //     baffle.second() = reverseFaceMap[baffle.second()];
            // }

            const auto oldMeshPoints(ppPtr().meshPoints());
            const auto oldMeshFaces(ppPtr().addressing());
            ppPtr = meshRefinement::makePatch
            (
                mesh,
                mesh.points(),
                adaptPatchIDs,
                bufZoneIDs
            );
            labelList ppMap;        // new-to-old point map
            labelList ppReverseMap; // old-to-new point map
            labelList ppFaceMap;    // new-to-old face map
            mapPatchAddressing
            (
                morphMap(),
                // Old patch to mesh addresssing
                oldMeshPoints,
                oldMeshFaces,
                ppPtr->meshPoints(),    //newMeshPoints,
                ppPtr().addressing(),   //newMeshFaces,

                ppMap,
                ppReverseMap,
                ppFaceMap
            );

            // ppLocalPoints
            meshRefinement::updateList(ppMap, GGG, ppLocalPoints);
            meshRefinement::updateList(ppMap, GGG, oldPPLocalPoints);
            // snapDist
            meshRefinement::updateList(ppMap, scalar(0), snapDist);
            // patchDisp
            meshRefinement::updateList(ppMap, vector::zero, patchDisp);
            // patchFaceToSurfaceRegion
            meshRefinement::updateList(ppFaceMap, -1, patchFaceToSurfaceRegion);
            // patchToFace
            meshRefinement::updateList(ppMap, -1, patchToFace);
            // patchToNearest
            meshRefinement::updateList(ppMap, GGG, patchToNearest);
            // patchToNormal
            meshRefinement::updateList(ppMap, GGG, patchToNormal);
            // patchToRegion
            meshRefinement::updateList(ppMap, -1, patchToRegion);
            // patchToPoint
            meshRefinement::updateList(ppMap, -1, patchToPoint);
            // patchToEdge
            meshRefinement::updateList(ppMap, -1, patchToEdge);
            // surfFeatEdgeRegionToStart
            updateMap(ppReverseMap, surfFeatEdgeRegionToStart);
            // patchConstraints
            meshRefinement::updateList(ppMap, pointConstraint(), patchConstraints);
        }

        if (debug & meshRefinement::MESH)
        {
            const_cast<Time&>(mesh.time())++;
            Info<< "Writing split-face mesh to time "
                << meshRefiner_.timeName() << endl;
            meshRefiner_.write
            (
                meshRefinement::debugType(debug),
                meshRefinement::writeType
                (
                    meshRefinement::writeLevel()
                  | meshRefinement::WRITEMESH
                ),
                meshRefiner_.timeName()
            );
            writeConstraints
            (
                "pp_splitFaces",
                ppPtr(),                //addressing only
                ppPtr().localPoints(),  //ppLocalPoints
                patchToNearest,
                patchConstraints
            );
        }

        topoSet::removeFiles(mesh);
        processorMeshes::removeFiles(mesh);
    }


    // Add buffer layers where needed
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  - thickness just enough not to get geom matching problems
    //    (but will produce illegal cells on warped geom)
    //  - inflation is done later on
    //  - maybe should be done in one go? Ie introduce new points at
    //    desired location directly.

    Info<< "Adding buffer layers ..." << endl;
    {
        // Create baffles (pairs of faces that share the same points)
        // Baffles stored as owner and neighbour face that have been created.
        // Note that 'internal' faceZones are converted to baffles only
        // so we can add layers to both sides. This should not be necessary
        // if addPatchCellLayer can do this.
        List<labelPair> baffles;
        {
            labelList originatingFaceZone;
            autoPtr<mapPolyMesh> mapPtr = meshRefiner_.createZoneBaffles
            (
                identity(mesh.faceZones().size()),
                baffles,
                originatingFaceZone
            );

            // Map all attraction info
            // ~~~~~~~~~~~~~~~~~~~~~~~
            if (mapPtr)
            {
                const auto oldMeshPoints(ppPtr().meshPoints());
                const auto oldMeshFaces(ppPtr().addressing());
                ppPtr = meshRefinement::makePatch
                (
                    mesh,
                    mesh.points(),
                    adaptPatchIDs,
                    bufZoneIDs
                );

                labelList ppMap;        // new-to-old point map
                labelList ppReverseMap; // old-to-new point map
                labelList ppFaceMap;    // new-to-old face map
                mapPatchAddressing
                (
                    mapPtr(),
                    // Old patch to mesh addresssing
                    oldMeshPoints,
                    oldMeshFaces,
                    ppPtr().meshPoints(),   //newMeshPoints,
                    ppPtr().addressing(),   //newMeshFaces,

                    ppMap,                  // new-to-old patch points
                    ppReverseMap,           // old-to-new patch points
                    ppFaceMap               // new-to-old patch faces
                );

                // ppLocalPoints
                meshRefinement::updateList(ppMap, GGG, ppLocalPoints);
                meshRefinement::updateList(ppMap, GGG, oldPPLocalPoints);
                // snapDist
                meshRefinement::updateList(ppMap, scalar(0), snapDist);
                // patchDisp
                meshRefinement::updateList(ppMap, vector::zero, patchDisp);
                // patchFaceToSurfaceRegion
                meshRefinement::updateList(ppFaceMap, -1, patchFaceToSurfaceRegion);
                // patchToFace
                meshRefinement::updateList(ppMap, -1, patchToFace);
                // patchToNearest
                meshRefinement::updateList(ppMap, GGG, patchToNearest);
                // patchToNormal
                meshRefinement::updateList(ppMap, GGG, patchToNormal);
                // patchToRegion
                meshRefinement::updateList(ppMap, -1, patchToRegion);
                // patchToPoint
                meshRefinement::updateList(ppMap, -1, patchToPoint);
                // patchToEdge
                meshRefinement::updateList(ppMap, -1, patchToEdge);
                // surfFeatEdgeRegionToStart
                updateMap(ppReverseMap, surfFeatEdgeRegionToStart);
                // patchConstraints
                meshRefinement::updateList(ppMap, pointConstraint(), patchConstraints);
            }
        }

        // Duplicate points on faceZones of type boundary. Only needed so
        // addPatchCellLayer can walk on the correct side and pointNormals
        // can be calculated. Should be not needed
        // once addPatchCellLayer can add to multiple sides ...
        labelList pointToMaster;
        {
            labelList numLayers(mesh.boundaryMesh().size(), 0);
            UIndirectList<label>(numLayers, bufPatchIDs) = 1;

            autoPtr<mapPolyMesh> mapPtr =
                snappyLayerDriver::dupFaceZonePoints
                (
                    meshRefiner_,
                    bufPatchIDs,    // patch indices (now including faceZones)
                    numLayers,      // num layers per patch
                    baffles,        // updated baffles
                    pointToMaster // per mesh point -1 or the master mesh point
                );

            if (mapPtr)
            {
                const auto oldMeshPoints(ppPtr().meshPoints());
                const auto oldMeshFaces(ppPtr().addressing());
                ppPtr = meshRefinement::makePatch
                (
                    mesh,
                    mesh.points(),
                    adaptPatchIDs,
                    bufZoneIDs
                );

                labelList ppMap;        // new-to-old point map
                labelList ppReverseMap; // old-to-new point map
                labelList ppFaceMap;    // new-to-old face map
                mapPatchAddressing
                (
                    mapPtr(),
                    // Old patch to mesh addresssing
                    oldMeshPoints,
                    oldMeshFaces,
                    ppPtr().meshPoints(),   //newMeshPoints,
                    ppPtr().addressing(),   //newMeshFaces,

                    ppMap,                  // new-to-old patch points
                    ppReverseMap,           // old-to-new patch points
                    ppFaceMap               // new-to-old patch faces
                );

                // ppLocalPoints
                meshRefinement::updateList(ppMap, GGG, ppLocalPoints);
                meshRefinement::updateList(ppMap, GGG, oldPPLocalPoints);
                // snapDist
                meshRefinement::updateList(ppMap, scalar(0), snapDist);
                // patchDisp
                meshRefinement::updateList(ppMap, vector::zero, patchDisp);
                // patchFaceToSurfaceRegion
                meshRefinement::updateList(ppFaceMap, -1, patchFaceToSurfaceRegion);
                // patchToFace
                meshRefinement::updateList(ppMap, -1, patchToFace);
                // patchToNearest
                meshRefinement::updateList(ppMap, GGG, patchToNearest);
                // patchToNormal
                meshRefinement::updateList(ppMap, GGG, patchToNormal);
                // patchToRegion
                meshRefinement::updateList(ppMap, -1, patchToRegion);
                // patchToPoint
                meshRefinement::updateList(ppMap, -1, patchToPoint);
                // patchToEdge
                meshRefinement::updateList(ppMap, -1, patchToEdge);
                // surfFeatEdgeRegionToStart
                updateMap(ppReverseMap, surfFeatEdgeRegionToStart);
                // patchConstraints
                meshRefinement::updateList(ppMap, pointConstraint(), patchConstraints);
            }
        }

        // Layer mesh modifier
        // - use intrusion, not extrusion (similar to splitting
        //   cells)
        // - so points outside stay the same / at the same location
        // - and same for faces
        addPatchCellLayer addLayer(mesh, true, false);

        // Get faces to add buffer layers to
        autoPtr<indirectPrimitivePatch> bufPatchPtr
        (
            meshRefinement::makePatch
            (
                mesh,
                mesh.points(),
                bufPatchIDs,
                bufZoneIDs
            )
        );

        pointField thickness
        (
            wantedThickness(bufPatchPtr(), 1e-1)  //cellSizeFraction
          * PatchTools::pointNormals(mesh, bufPatchPtr())
        );

        // Make sure to adhere to constraints
        {
            // Bit convoluted: have constraint patches on mesh but want
            // to see their effect on overall outside-patch ('bufPatch').
            // - construct pointMesh from polyMesh and fixedValue bc
            //   displacement field
            // - expand outside-patch to mesh
            // - override fixedValue bcs and use multi-patch constraints
            // - map mesh point values back to outside-patch
            const pointMesh& pMesh = pointMesh::New(mesh);
            const labelList& mp = bufPatchPtr().meshPoints();

            tmp<pointVectorField> tdisp
            (
                meshRefinement::makeDisplacementField
                (
                    pMesh,
                    adaptPatchIDs,
                    zeroGradientPointPatchVectorField::typeName
                )
            );
            // Set internal field
            UIndirectList<point>(tdisp.ref(), mp) = thickness;
            // Take over onto boundary field. Needed since constraint
            // patch might be before the fixedValue patches so its
            // value gets overwritten
            for (auto& ppf : tdisp.ref().boundaryFieldRef())
            {
                ppf == ppf.patchInternalField();
            }
            // Adhere to multi-point constraints
            const pointConstraints& pcs = pointConstraints::New(pMesh);

            pcs.constrainDisplacement(tdisp.ref(), false);
            thickness = UIndirectList<point>(tdisp(), mp);
        }



        // Print a bit
        {
            // See snappyLayerDriver::calculateLayerThickness.
            const auto& pbm = mesh.boundaryMesh();
            label maxLen = 0;
            for (const label patchi : bufPatchIDs)
            {
                maxLen = max(maxLen, label(pbm[patchi].name().size()));
            }

            const int oldPrecision = Info.stream().precision();

            Info<< nl
                << setf(ios_base::left) << setw(maxLen) << "patch"
                << setw(0) << " faces    layers thickness[m]" << nl
                << setf(ios_base::left) << setw(maxLen) << "-----"
                << setw(0) << " -----    ------ ------------" << endl;

            for (const label patchi : bufPatchIDs)
            {
                Info<< setf(ios_base::left) << setw(maxLen)
                    << pbm[patchi].name() << setprecision(3)
                    << " " << setw(8)
                    << returnReduce(pbm[patchi].size(), sumOp<scalar>())
                    << " " << setw(6) << 1
                    << " " << setw(8) << gAverage(mag(thickness))
                    << endl;
            }
            Info<< setprecision(oldPrecision) << endl;
        }


        // Do mesh changes : introduce point, faces, cells
        autoPtr<mapPolyMesh> mapPtr = addBufferLayers
        (
            bufPatchPtr(),
            -thickness,
            addLayer
        );

        // Update numbering on baffles. Note: uses
        // geometric tolerance (still?) - could be avoided since now we have
        // addPatchCellLayer still intact. However would like to avoid
        // use of addPatchCellLayer altogether.
        if (mapPtr)
        {
            // baffles are in mesh face indices
            const auto& reverseFaceMap = mapPtr().reverseFaceMap();
            for (auto& baffle : baffles)
            {
                baffle.first() = reverseFaceMap[baffle.first()];
                baffle.second() = reverseFaceMap[baffle.second()];
            }
            // Note: points have been added as duplicate of boundary
            //       points so pointToMaster needs to be mapped using a
            //       old-to-new map.
            {
                const auto& oldToNew = mapPtr().reversePointMap();
                labelList newPointToMaster(mesh.nPoints(), -1);
                forAll(pointToMaster, oldPointi)
                {
                    const label dupPointi = pointToMaster[oldPointi];
                    const label newPointi = oldToNew[oldPointi];
                    if
                    (
                        dupPointi != -1
                     && newPointi != -1
                     && oldToNew[dupPointi] != -1
                    )
                    {
                        newPointToMaster[newPointi] = oldToNew[dupPointi];
                    }
                }
                pointToMaster = std::move(newPointToMaster);
            }

            const auto oldMeshPoints(ppPtr().meshPoints());
            const auto oldMeshFaces(ppPtr().addressing());
            ppPtr = meshRefinement::makePatch
            (
                mesh,
                mesh.points(),
                adaptPatchIDs,
                bufZoneIDs
            );
            labelList ppMap;        // new-to-old point map
            labelList ppReverseMap; // old-to-new point map
            labelList ppFaceMap;    // new-to-old face map
            mapPatchAddressing
            (
                mapPtr(),
                // Old patch to mesh addresssing
                oldMeshPoints,
                oldMeshFaces,
                ppPtr().meshPoints(),   //newMeshPoints,
                ppPtr().addressing(),   //newMeshFaces,

                ppMap,
                ppReverseMap,
                ppFaceMap
            );

            // ppLocalPoints
            meshRefinement::updateList(ppMap, GGG, ppLocalPoints);
            meshRefinement::updateList(ppMap, GGG, oldPPLocalPoints);
            // snapDist
            meshRefinement::updateList(ppMap, scalar(0), snapDist);
            // patchDisp
            meshRefinement::updateList(ppMap, vector::zero, patchDisp);
            // patchFaceToSurfaceRegion
            meshRefinement::updateList(ppFaceMap, -1, patchFaceToSurfaceRegion);
            // patchToFace
            meshRefinement::updateList(ppMap, -1, patchToFace);
            // patchToNearest
            meshRefinement::updateList(ppMap, GGG, patchToNearest);
            // patchToNormal
            meshRefinement::updateList(ppMap, GGG, patchToNormal);
            // patchToRegion
            meshRefinement::updateList(ppMap, -1, patchToRegion);
            // patchToPoint
            meshRefinement::updateList(ppMap, -1, patchToPoint);
            // patchToEdge
            meshRefinement::updateList(ppMap, -1, patchToEdge);
            // surfFeatEdgeRegionToStart
            updateMap(ppReverseMap, surfFeatEdgeRegionToStart);
            // patchConstraints
            meshRefinement::updateList(ppMap, pointConstraint(), patchConstraints);
        }
        // Dump the attraction
        if (debug & meshRefinement::MESH)
        {
            writeConstraints
            (
                "pp_after_buffer_layers",
                ppPtr().localPoints(),  //ppLocalPoints
                patchToNearest,
                patchConstraints
            );
        }
        // Merge any duplicated points
        mapPtr = meshRefiner_.mergePoints(pointToMaster);
        if (mapPtr)
        {
            // baffles are in mesh face indices
            const auto& reverseFaceMap = mapPtr().reverseFaceMap();
            for (auto& baffle : baffles)
            {
                baffle.first() = reverseFaceMap[baffle.first()];
                baffle.second() = reverseFaceMap[baffle.second()];
            }
            // duplicated points
            pointToMaster.clear();
            //meshRefinement::updateList(pointToMaster, -1, mapPtr().pointMap());

            const auto oldMeshPoints(ppPtr().meshPoints());
            const auto oldMeshFaces(ppPtr().addressing());
            ppPtr = meshRefinement::makePatch
            (
                mesh,
                mesh.points(),
                adaptPatchIDs,
                bufZoneIDs
            );
            labelList ppMap;        // new-to-old point map
            labelList ppReverseMap; // old-to-new point map
            labelList ppFaceMap;    // new-to-old face map
            mapPatchAddressing
            (
                mapPtr(),
                // Old patch to mesh addresssing
                oldMeshPoints,
                oldMeshFaces,
                ppPtr().meshPoints(),   //newMeshPoints,
                ppPtr().addressing(),   //newMeshFaces,

                ppMap,
                ppReverseMap,
                ppFaceMap
            );

            // ppLocalPoints
            meshRefinement::updateList(ppMap, GGG, ppLocalPoints);
            meshRefinement::updateList(ppMap, GGG, oldPPLocalPoints);
            // snapDist
            meshRefinement::updateList(ppMap, scalar(0), snapDist);
            // patchDisp
            meshRefinement::updateList(ppMap, vector::zero, patchDisp);
            // patchFaceToSurfaceRegion
            meshRefinement::updateList(ppFaceMap, -1, patchFaceToSurfaceRegion);
            // patchToFace
            meshRefinement::updateList(ppMap, -1, patchToFace);
            // patchToNearest
            meshRefinement::updateList(ppMap, GGG, patchToNearest);
            // patchToNormal
            meshRefinement::updateList(ppMap, GGG, patchToNormal);
            // patchToRegion
            meshRefinement::updateList(ppMap, -1, patchToRegion);
            // patchToPoint
            meshRefinement::updateList(ppMap, -1, patchToPoint);
            // patchToEdge
            meshRefinement::updateList(ppMap, -1, patchToEdge);
            // surfFeatEdgeRegionToStart
            updateMap(ppReverseMap, surfFeatEdgeRegionToStart);
            // patchConstraints
            meshRefinement::updateList(ppMap, pointConstraint(), patchConstraints);
        }
        // Merge any introduced baffles (from faceZones of faceType 'internal')
        mapPtr = meshRefiner_.mergeZoneBaffles
        (
            true,   // internal zones
            false   // baffle zones
        );

        if (mapPtr)
        {
            // baffles are in mesh face indices
            const auto& reverseFaceMap = mapPtr().reverseFaceMap();
            DynamicList<labelPair> newBaffles(baffles.size());
            for (auto& baffle : baffles)
            {
                labelPair newBaffle
                (
                    reverseFaceMap[baffle.first()],
                    reverseFaceMap[baffle.second()]
                );
                if (newBaffle.first() != -1 && newBaffle.second() != -1)
                {
                    newBaffles.append(newBaffle);
                }
            }
            baffles = std::move(newBaffles);
            const auto oldMeshPoints(ppPtr().meshPoints());
            const auto oldMeshFaces(ppPtr().addressing());
            ppPtr = meshRefinement::makePatch
            (
                mesh,
                mesh.points(),
                adaptPatchIDs,
                bufZoneIDs
            );
            labelList ppMap;        // new-to-old point map
            labelList ppReverseMap; // old-to-new point map
            labelList ppFaceMap;    // new-to-old face map
            mapPatchAddressing
            (
                mapPtr(),
                // Old patch to mesh addresssing
                oldMeshPoints,
                oldMeshFaces,
                ppPtr().meshPoints(),   //newMeshPoints,
                ppPtr().addressing(),   //newMeshFaces,

                ppMap,
                ppReverseMap,
                ppFaceMap
            );

            // ppLocalPoints
            meshRefinement::updateList(ppMap, GGG, ppLocalPoints);
            meshRefinement::updateList(ppMap, GGG, oldPPLocalPoints);
            // snapDist
            meshRefinement::updateList(ppMap, scalar(0), snapDist);
            // patchDisp
            meshRefinement::updateList(ppMap, vector::zero, patchDisp);
            // patchFaceToSurfaceRegion
            meshRefinement::updateList(ppFaceMap, -1, patchFaceToSurfaceRegion);
            // patchToFace
            meshRefinement::updateList(ppMap, -1, patchToFace);
            // patchToNearest
            meshRefinement::updateList(ppMap, GGG, patchToNearest);
            // patchToNormal
            meshRefinement::updateList(ppMap, GGG, patchToNormal);
            // patchToRegion
            meshRefinement::updateList(ppMap, -1, patchToRegion);
            // patchToPoint
            meshRefinement::updateList(ppMap, -1, patchToPoint);
            // patchToEdge
            meshRefinement::updateList(ppMap, -1, patchToEdge);
            // surfFeatEdgeRegionToStart
            updateMap(ppReverseMap, surfFeatEdgeRegionToStart);
            // patchConstraints
            meshRefinement::updateList(ppMap, pointConstraint(), patchConstraints);

            if (debug & meshRefinement::MESH)
            {
                const_cast<Time&>(mesh.time())++;
                Info<< "Writing baffle-merged mesh to time "
                    << meshRefiner_.timeName() << endl;
                meshRefiner_.write
                (
                    meshRefinement::debugType(debug),
                    meshRefinement::writeType
                    (
                        meshRefinement::writeLevel()
                      | meshRefinement::WRITEMESH
                    ),
                    meshRefiner_.timeName()
                );
                writeConstraints
                (
                    "pp_after_baffle_merging",
                    ppPtr().localPoints(),  //ppLocalPoints
                    patchToNearest,
                    patchConstraints
                );
            }
        }
    }

    {
        // Create pointMesh with the correct patches
        const pointMesh& pMesh = makePointMesh
        (
            ppPtr(),
            patchConstraints,
            labelList::null(),  //adaptPatchIDs,  // convert into meshPointPatch
            allSurfPatchName,   //word::null,
            allEdgePatchName,
            allPointPatchName
        );
        const auto& pBnd = pMesh.boundary();
        const label allSurfPatchi = pBnd.findPatchID(allSurfPatchName);
        if (allSurfPatchi != -1)
        {
            adaptPatchIDs.push_uniq(allSurfPatchi);
        }
        const label allEdgePatchi = pBnd.findPatchID(allEdgePatchName);
        if (allEdgePatchi != -1)
        {
            adaptPatchIDs.push_uniq(allEdgePatchi);
        }
        const label allPointPatchi = pBnd.findPatchID(allPointPatchName);
        if (allPointPatchi != -1)
        {
            adaptPatchIDs.push_uniq(allPointPatchi);
        }

        autoPtr<displacementMotionSolver> motionPtr
        (
            makeMotionSolver
            (
                pMesh,
                ppPtr,
                snapDict,
                adaptPatchIDs
            )
        );

        if (true)
        {
            // Move points to smoothed patch points first to improve snapping
            // - this is a bit of a hack to make life easier for the
            // motionSolver.
            Info<< "Pre-snapping to smoothed patch points" << endl;
            if (debug)
            {
                OBJstream os(mesh.time().timePath()/"pp_OLD_smoothed.obj");
                os.write(ppPtr().localFaces(), oldPPLocalPoints, false);
                Pout<< "Write pp after smoothing to " << os.name()
                    << " nVertices:" << os.nVertices() << endl;
            }

            // Could introduce baffles for internal face-zone faces here.
            // However this then requires coding to check mesh quality through
            // these baffles. Much easier to put the face-zones faces into
            // a faceZone instead. Does not seem to work ok though...

            // Insert as bc to motionSmoother. Note:
            // - wanted displacement is wrt current mesh (mesh not moved yet
            //   so can use localPoints)
            // - calculates bc on pointDisplacement relative to points0
            setDisplacement
            (
                ppPtr(),
                oldPPLocalPoints-ppPtr().localPoints(),
                adaptPatchIDs,
                motionPtr().points0(),
                motionPtr().pointDisplacement()
            );
            // Solve internal displacement
            tmp<pointField> tnewPoints(motionPtr->newPoints());

            // Move points
            mesh.movePoints(tnewPoints);

            if (debug & meshRefinement::MESH)
            {
                const_cast<Time&>(mesh.time())++;
                Info<< "Writing smoothed mesh to time "
                    << meshRefiner_.timeName() << endl;
                meshRefiner_.write
                (
                    meshRefinement::debugType(debug),
                    meshRefinement::writeType
                    (
                        meshRefinement::writeLevel()
                    | meshRefinement::WRITEMESH
                    ),
                    meshRefiner_.timeName()
                );
            }
        }


        // Insert attraction as bc on pointDisplacement.
        // (Note that this takes displacement relative to points0)
        // pointDisplacement :
        //  - non-meshed patches : slip
        //  - meshed patches : slip using constraints from meshPointPatches
        // Note: previous motionSmoother step has moved mesh points so use
        // up to date pp points
        setDisplacement
        (
            ppPtr(),
            patchToNearest-pointField(mesh.points(), ppPtr().meshPoints()),
            adaptPatchIDs,
            motionPtr().points0(),
            motionPtr().pointDisplacement()
        );

        if (debug & meshRefinement::MESH)
        {
            const_cast<Time&>(mesh.time())++;
            Info<< "Writing pre-snapping mesh to time "
                << meshRefiner_.timeName() << endl;
            meshRefiner_.write
            (
                meshRefinement::debugType(debug),
                meshRefinement::writeType
                (
                    meshRefinement::writeLevel()
                  | meshRefinement::WRITEMESH
                ),
                meshRefiner_.timeName()
            );
            OBJstream os(mesh.time().timePath()/"preSnapOutside.obj");
            Info<< "Writing wanted outside to " << os.name() << endl;
            os.write(ppPtr().localFaces(), patchToNearest, false);
        }

        // Solve internal displacement
        tmp<pointField> tnewPoints(motionPtr->newPoints());

        // Move points
        mesh.movePoints(tnewPoints);

        // Reset moving flag to avoid problems with topo changes
        mesh.moving(false);

        // Update pp for new mesh points. Ok as long as we also update geometry
        // and wanted displacement (usually zero if mesh motion has succeeded)
        ppLocalPoints = pointField(mesh.points(), ppPtr().meshPoints());
        patchDisp -= (ppLocalPoints-ppPtr().localPoints());

        if (debug & meshRefinement::MESH)
        {
            const_cast<Time&>(mesh.time())++;
            Info<< "Writing smoothed LAYER mesh to time "
                << meshRefiner_.timeName() << endl;
            meshRefiner_.write
            (
                meshRefinement::debugType(debug),
                meshRefinement::writeType
                (
                    meshRefinement::writeLevel()
                  | meshRefinement::WRITEMESH
                ),
                meshRefiner_.timeName()
            );
        }
    }



//         // Split problematic cells
//         {
//             Info<< nl << "Checking moved mesh ..." << endl;
//             faceSet wrongFaces(mesh, "wrongFaces", mesh.nFaces()/1000);
//             motionSmoother::checkMesh
//             (
//                 false,
//                 mesh,
//                 motionDict,
//                 identity(mesh.nFaces()),
//                 internalBaffles,
//                 wrongFaces,
//                 false           // dryRun_
//             );
//             const label nWrong(returnReduce(wrongFaces.size(), sumOp<label>()));
//             Info<< "Detected " << nWrong
//                 << " illegal faces"
//                 << " (concave, zero area or negative cell pyramid volume)"
//                 << endl;

//             //if (nWrong)
//             if (false)
//             {
//                 bitSet decomposeCell(mesh.nCells());
//                 for (const label facei : wrongFaces)
//                 {
//                     const label own = mesh.faceOwner()[facei];
//                     if (!tetMatcher::test(mesh, own))
//                     {
//                         decomposeCell.set(own);
//                     }
//                     if (mesh.isInternalFace(facei))
//                     {
//                         const label nei = mesh.faceNeighbour()[facei];
//                         if (!tetMatcher::test(mesh, nei))
//                         {
//                             decomposeCell.set(nei);
//                         }
//                     }
//                 }
//                 Pout<< "spliyyinG :" << decomposeCell.count() << " cells"
//                     << endl;

//                 tetDecomposer tetDecomp(mesh);
//                 polyTopoChange meshMod(mesh);
//                 tetDecomp.setRefinement
//                 (
//                     tetDecomposer::FACE_CENTRE_TRIS,
//                     decomposeCell,
//                     meshMod
//                 );

//                 // Save old meshPoints before changing mesh
//                 const Map<label> oldMeshPointMap(ppPtr->meshPointMap());

//                 Pout<< "old pp points:" << ppPtr->nPoints()
//                     << " oldMeshPointMap:" << oldMeshPointMap.size()
//                     << endl;

//                 // Remove any unnecessary fields
//                 meshMoverPtr.clear();
//                 motionPtr.clear();
//                 ppPtr.clear();
//                 mesh.clearOut();
//                 mesh.moving(false);

//                 // Change the mesh (no inflation)
//                 autoPtr<mapPolyMesh> mapPtr = meshMod.changeMesh(mesh, false);
//                 mapPolyMesh& map = *mapPtr;

//                 // Update fields
//                 mesh.updateMesh(map);

//                 // Move mesh (since morphing does not do this)
//                 if (map.hasMotionPoints())
//                 {
//                     mesh.movePoints(map.preMotionPoints());
//                 }
//                 else
//                 {
//                     mesh.clearOut();
//                 }

//                 // Reset the instance for if in overwrite mode
//                 mesh.setInstance(meshRefiner_.timeName());
//                 meshRefiner_.setInstance(mesh.facesInstance());

//                 //- Update numbering on tet-decomposition engine
//                 tetDecomp.updateMesh(map);

//                 bitSet isChangedFace(mesh.nFaces());
//                 forAll(map.cellMap(), celli)
//                 {
//                     if (decomposeCell[map.cellMap()[celli]])
//                     {
//                         isChangedFace.set(mesh.cells()[celli]);
//                     }
//                 }
//                 syncTools::syncFaceList
//                 (
//                     mesh,
//                     isChangedFace,
//                     orEqOp<unsigned int>()
//                 );

//                 Pout<< "isChangedFace :" << decomposeCell.count() << " faces"
//                     << endl;


//                 // Update intersection info
//                 meshRefiner_.updateMesh(map, isChangedFace.toc());


//                 if (debug&meshRefinement::MESH)
//                 {
//                     const_cast<Time&>(mesh.time())++;
//                     Info<< "Writing tet-decomp mesh to time "
//                         << meshRefiner_.timeName() << endl;
//                     meshRefiner_.write
//                     (
//                         meshRefinement::debugType(debug),
//                         meshRefinement::writeType
//                         (
//                             meshRefinement::writeLevel()
//                           | meshRefinement::WRITEMESH
//                         ),
//                         mesh.time().path()/meshRefiner_.timeName()
//                     );
//                 }


//                 // Update local mesh data
//                 // ~~~~~~~~~~~~~~~~~~~~~~

//                 Info<< "Updating for tet-decomp" << endl;

//                 // baffles
//                 forAll(internalBaffles, i)
//                 {
//                     labelPair& baffle = internalBaffles[i];
//                     baffle.first() = map.reverseFaceMap()[baffle.first()];
//                     baffle.second() = map.reverseFaceMap()[baffle.second()];
//                 }

//                 // re-do patch (since faces might have been split)
//                 ppPtr = meshRefinement::makePatch(mesh, adaptPatchIDs);
//     Pout<< "new pp points:" << ppPtr->nPoints()
//         << " new meshPointMap:" << ppPtr->meshPointMap().size()
//         << endl;
//                 const auto& mp = ppPtr->meshPoints();
//                 // pointMap (new-to-old) for pp points. Might have new
//                 // face-centre points - these get mapped from any point on
//                 // the originating face.
//                 labelList ppMap(mp.size(), -1);
//                 forAll(mp, i)
//                 {
//                     const label oldMeshPointi = map.pointMap()[mp[i]];
//                     const auto mpFnd = oldMeshPointMap.find(oldMeshPointi);
//                     if (mpFnd)
//                     {
//                         ppMap[i] = mpFnd();
//                     }
//                 }

//                 // patchDisp
//                 meshRefinement::updateList(ppMap, vector::zero, patchDisp);
//                 // snapDist
//                 meshRefinement::updateList(ppMap, scalar(0), snapDist);
//                 // patchFeaturePoint
//                 meshRefinement::updateList
//                 (
//                     ppMap,
//                     vector::zero,
//                     patchFeaturePoint
//                 );
//                 // patchConstraints
//                 meshRefinement::updateList
//                 (
//                     ppMap,
//                     pointConstraint(),
//                     patchConstraints
//                 );
// //                // nearestPoint
// //                meshRefinement::updateList(ppMap, vector::zero, nearestPoint);
// //                // nearestNormal
// //                meshRefinement::updateList(ppMap, vector::zero, nearestNormal);
// //                // disp
// //                meshRefinement::updateList(ppMap, vector::zero, disp);

//                 // Maximum distance to attract to nearest feature on surface
//                 snapDist = calcSnapDistance(mesh, snapParams, ppPtr());
//                 ppLocalPoints = pointField(mesh.points(), ppPtr().meshPoints());

//                 Info<< "DONE Updating for tet-decomp" << endl;
//             }
//         }
//     }

    if (debug & meshRefinement::MESH)
    {
        const_cast<Time&>(mesh.time())++;
    }
}


// ************************************************************************* //
