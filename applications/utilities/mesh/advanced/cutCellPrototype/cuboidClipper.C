/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F i e l d       | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O p e r a t i o n |
    \\  /    A n d           |
     \\/     M a n i p u l a t i o n  |
-------------------------------------------------------------------------------
    Description
        Utility-owned clipping of a cut-cell interface for an axis-aligned
        six-region cuboid.
\*---------------------------------------------------------------------------*/

#include "cuboidClipper.H"

#include "boundBox.H"
#include "face.H"
#include "mapPolyMesh.H"
#include "polyBoundaryMesh.H"
#include "meshCutter.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "polyTopoChange.H"
#include "triSurface.H"
#include "wallPolyPatch.H"

namespace Foam
{

namespace
{

struct cuboidPlane
{
    label region;
    label axis;
    scalar coordinate;
    scalar sign;
    vector normal;
};


label dominantAxis(const vector& v)
{
    const vector av(cmptMag(v));

    if (av.x() >= av.y() && av.x() >= av.z())
    {
        return 0;
    }
    if (av.y() >= av.z())
    {
        return 1;
    }
    return 2;
}


scalar component(const point& p, const label axis)
{
    return p[axis];
}


point makePoint
(
    const label fixedAxis,
    const scalar fixedValue,
    const label firstAxis,
    const scalar firstValue,
    const label secondAxis,
    const scalar secondValue
)
{
    point p(Zero);
    p[fixedAxis] = fixedValue;
    p[firstAxis] = firstValue;
    p[secondAxis] = secondValue;
    return p;
}


vector polygonArea(const pointField& points)
{
    vector area(Zero);
    const point& p0 = points[0];

    for (label i = 1; i + 1 < points.size(); ++i)
    {
        area += (points[i] - p0) ^ (points[i + 1] - p0);
    }

    return 0.5*area;
}


label findPoint
(
    const pointField& points,
    const point& target,
    const scalar tolerance
)
{
    const scalar tolSqr = sqr(tolerance);

    forAll(points, pointi)
    {
        if (magSqr(points[pointi] - target) <= tolSqr)
        {
            return pointi;
        }
    }

    return -1;
}


List<cuboidPlane> readCuboidPlanes
(
    const triSurface& surface,
    const scalar tolerance
)
{
    if (surface.patches().size() != 6)
    {
        FatalErrorInFunction
            << "Multi-region clipping requires exactly six STL regions; "
            << "found " << surface.patches().size()
            << ". Arbitrary multi-region surfaces remain unsupported."
            << exit(FatalError);
    }

    const pointField& surfacePoints = surface.localPoints();
    const List<labelledTri>& triangles = surface.localFaces();
    List<cuboidPlane> planes(6);
    boolList found(6, false);
    const scalar normalTolerance = 1e-6;

    forAll(surface.patches(), regioni)
    {
        label firstFace = -1;
        label nRegionFaces = 0;

        forAll(triangles, facei)
        {
            if (triangles[facei].region() == regioni)
            {
                if (firstFace < 0)
                {
                    firstFace = facei;
                }
                ++nRegionFaces;
            }
        }

        if (nRegionFaces == 0)
        {
            FatalErrorInFunction
                << "Cuboid region " << regioni << " is empty."
                << exit(FatalError);
        }

        const labelledTri& first = triangles[firstFace];
        const point& p0 = surfacePoints[first[0]];
        const point& p1 = surfacePoints[first[1]];
        const point& p2 = surfacePoints[first[2]];
        vector normal = (p1 - p0) ^ (p2 - p0);

        if (mag(normal) <= VSMALL)
        {
            FatalErrorInFunction
                << "Cuboid region " << regioni
                << " contains a degenerate triangle."
                << exit(FatalError);
        }

        normal /= mag(normal);
        const label axis = dominantAxis(normal);

        if (mag(normal[axis]) < 1 - normalTolerance)
        {
            FatalErrorInFunction
                << "Cuboid region " << regioni
                << " is not axis aligned; generalized clipping is "
                << "limited to an axis-aligned six-region cuboid."
                << exit(FatalError);
        }

        const scalar coordinate = component(p0, axis);
        const label side = normal[axis] > 0 ? 1 : -1;
        const label slot = 2*axis + (side > 0);

        forAll(triangles, facei)
        {
            if (triangles[facei].region() != regioni)
            {
                continue;
            }

            const labelledTri& tri = triangles[facei];
            const vector triNormal =
                ((surfacePoints[tri[1]] - surfacePoints[tri[0]])
               ^ (surfacePoints[tri[2]] - surfacePoints[tri[0]]));

            if
            (
                mag(triNormal) <= VSMALL
             || mag((triNormal/mag(triNormal)) & normal) < 1 - normalTolerance
             || mag(component(surfacePoints[tri[0]], axis) - coordinate)
                    > tolerance
            )
            {
                FatalErrorInFunction
                    << "STL region " << regioni
                    << " is not one planar axis-aligned cuboid face."
                    << exit(FatalError);
            }
        }

        if (found[slot])
        {
            FatalErrorInFunction
                << "Duplicate axis-aligned cuboid face for region "
                << regioni << exit(FatalError);
        }

        found[slot] = true;
        planes[slot] = {regioni, axis, coordinate, scalar(side), normal};
    }

    forAll(found, slot)
    {
        if (!found[slot])
        {
            FatalErrorInFunction
                << "The six STL regions do not form one face for every "
                << "axis and side." << exit(FatalError);
        }
    }

    return planes;
}


boundBox cellBounds
(
    const polyMesh& mesh,
    const label celli,
    const label addedCelli
)
{
    DynamicList<point> points;

    const labelListList& cellPoints = mesh.cellPoints();

    for (const label pointi : cellPoints[celli])
    {
        points.append(mesh.points()[pointi]);
    }
    for (const label pointi : cellPoints[addedCelli])
    {
        points.append(mesh.points()[pointi]);
    }

    points.shrink();
    return boundBox(pointField(points), false);
}


bool clippedFacePoints
(
    const face& originalFace,
    const pointField& points,
    const boundBox& surfaceBounds,
    const scalar tolerance,
    pointField& clipped
)
{
    if (originalFace.size() != 4)
    {
        return false;
    }

    label fixedAxis = -1;
    for (label axis = 0; axis < 3; ++axis)
    {
        scalar minValue = GREAT;
        scalar maxValue = -GREAT;

        forAll(originalFace, fp)
        {
            minValue = min(minValue, points[originalFace[fp]][axis]);
            maxValue = max(maxValue, points[originalFace[fp]][axis]);
        }

        if (maxValue - minValue <= tolerance)
        {
            fixedAxis = axis;
            break;
        }
    }

    if (fixedAxis < 0)
    {
        return false;
    }

    const scalar fixedValue = points[originalFace[0]][fixedAxis];

    if
    (
        fixedValue < surfaceBounds.min()[fixedAxis] - tolerance
     || fixedValue > surfaceBounds.max()[fixedAxis] + tolerance
    )
    {
        return false;
    }

    label firstAxis = -1;
    label secondAxis = -1;
    for (label axis = 0; axis < 3; ++axis)
    {
        if (axis != fixedAxis)
        {
            if (firstAxis < 0)
            {
                firstAxis = axis;
            }
            else
            {
                secondAxis = axis;
            }
        }
    }

    scalar firstMin = GREAT;
    scalar firstMax = -GREAT;
    scalar secondMin = GREAT;
    scalar secondMax = -GREAT;

    forAll(originalFace, fp)
    {
        const point& p = points[originalFace[fp]];
        firstMin = min(firstMin, p[firstAxis]);
        firstMax = max(firstMax, p[firstAxis]);
        secondMin = min(secondMin, p[secondAxis]);
        secondMax = max(secondMax, p[secondAxis]);
    }

    firstMin = max(firstMin, surfaceBounds.min()[firstAxis]);
    firstMax = min(firstMax, surfaceBounds.max()[firstAxis]);
    secondMin = max(secondMin, surfaceBounds.min()[secondAxis]);
    secondMax = min(secondMax, surfaceBounds.max()[secondAxis]);

    if
    (
        firstMax - firstMin <= tolerance
     || secondMax - secondMin <= tolerance
    )
    {
        return false;
    }

    clipped.setSize(4);
    clipped[0] =
        makePoint
        (
            fixedAxis, fixedValue,
            firstAxis, firstMin,
            secondAxis, secondMin
        );
    clipped[1] =
        makePoint
        (
            fixedAxis, fixedValue,
            firstAxis, firstMax,
            secondAxis, secondMin
        );
    clipped[2] =
        makePoint
        (
            fixedAxis, fixedValue,
            firstAxis, firstMax,
            secondAxis, secondMax
        );
    clipped[3] =
        makePoint
        (
            fixedAxis, fixedValue,
            firstAxis, firstMin,
            secondAxis, secondMax
        );

    return true;
}



label pointFor
(
    polyTopoChange& meshMod,
    const polyMesh& mesh,
    const point& p,
    const scalar tolerance,
    DynamicList<point>& addedPoints,
    DynamicList<label>& addedPointLabels
)
{
    const label existing = findPoint(mesh.points(), p, tolerance);
    if (existing >= 0)
    {
        return existing;
    }

    forAll(addedPoints, i)
    {
        if (magSqr(addedPoints[i] - p) <= sqr(tolerance))
        {
            return addedPointLabels[i];
        }
    }

    const label pointi = meshMod.addPoint(p, -1, -1, true);
    addedPoints.append(p);
    addedPointLabels.append(pointi);
    return pointi;
}


void updateFaceRegions
(
    const mapPolyMesh& morphMap,
    const Map<label>& oldFaceRegions,
    const DynamicList<label>& newActionFaces,
    const DynamicList<label>& newActionRegions,
    Map<label>& faceRegions
)
{
    faceRegions.clear();

    forAllConstIters(oldFaceRegions, iter)
    {
        const label oldFacei = iter.key();
        const label newFacei = morphMap.reverseFaceMap()[oldFacei];

        if (newFacei >= 0)
        {
            faceRegions.insert(newFacei, iter.val());
        }
    }

    forAll(newActionFaces, i)
    {
        const label newFacei =
            morphMap.reverseFaceMap()[newActionFaces[i]];

        if (newFacei < 0)
        {
            FatalErrorInFunction
                << "Clipped cuboid face action " << newActionFaces[i]
                << " disappeared during topology change."
                << exit(FatalError);
        }

        faceRegions.insert(newFacei, newActionRegions[i]);
    }
}


struct rebuiltFace
{
    face geometry;
    label owner;
    label neighbour;
    label patch;
    label zone;
    bool zoneFlip;
};


bool sameFaceGeometry
(
    const face& a,
    const face& b
)
{
    if (a.size() != b.size())
    {
        return false;
    }

    forAll(a, i)
    {
        if (!b.found(a[i]))
        {
            return false;
        }
    }

    return true;
}


label findRebuiltFace
(
    const List<rebuiltFace>& faces,
    const face& candidate
)
{
    forAll(faces, facei)
    {
        if (sameFaceGeometry(faces[facei].geometry, candidate))
        {
            return facei;
        }
    }

    return -1;
}


bool sameFacePoints
(
    const pointField& a,
    const pointField& b,
    const scalar tolerance
)
{
    if (a.size() != b.size())
    {
        return false;
    }

    forAll(a, i)
    {
        bool found = false;

        forAll(b, j)
        {
            if (magSqr(a[i] - b[j]) <= sqr(tolerance))
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            return false;
        }
    }

    return true;
}


label sourceFaceZone
(
    const polyMesh& mesh,
    const pointField& candidatePoints,
    const scalar tolerance,
    bool& zoneFlip
)
{
    zoneFlip = false;

    forAll(mesh.faceZones(), zonei)
    {
        const faceZone& zone = mesh.faceZones()[zonei];

        forAll(zone, i)
        {
            const label facei = zone[i];
            const face& sourceFace = mesh.faces()[facei];
            pointField sourcePoints(sourceFace.size());

            forAll(sourceFace, fp)
            {
                sourcePoints[fp] = mesh.points()[sourceFace[fp]];
            }

            if (!sameFacePoints(candidatePoints, sourcePoints, tolerance))
            {
                continue;
            }

            const vector sourceArea = polygonArea(sourcePoints);
            const vector candidateArea = polygonArea(candidatePoints);
            const bool sameOrientation =
                (sourceArea & candidateArea) >= 0;

            zoneFlip =
                zone.flipMap()[i] ^ !sameOrientation;
            return zonei;
        }
    }

    return -1;
}


label pointForDirect
(
    polyTopoChange& meshMod,
    const point& p,
    const scalar tolerance
)
{
    forAll(meshMod.points(), pointi)
    {
        if (magSqr(meshMod.points()[pointi] - p) <= sqr(tolerance))
        {
            return pointi;
        }
    }

    return meshMod.addPoint(p, -1, -1, true);
}


face boxSubFace
(
    polyTopoChange& meshMod,
    const label fixedAxis,
    const scalar fixedValue,
    const label side,
    const label firstAxis,
    const scalar firstMin,
    const scalar firstMax,
    const label secondAxis,
    const scalar secondMin,
    const scalar secondMax,
    const scalar tolerance
)
{
    pointField points(4);
    points[0] = makePoint
    (
        fixedAxis, fixedValue,
        firstAxis, firstMin,
        secondAxis, secondMin
    );
    points[1] = makePoint
    (
        fixedAxis, fixedValue,
        firstAxis, firstMax,
        secondAxis, secondMin
    );
    points[2] = makePoint
    (
        fixedAxis, fixedValue,
        firstAxis, firstMax,
        secondAxis, secondMax
    );
    points[3] = makePoint
    (
        fixedAxis, fixedValue,
        firstAxis, firstMin,
        secondAxis, secondMax
    );

    face result(4);
    forAll(result, i)
    {
        result[i] = pointForDirect(meshMod, points[i], tolerance);
    }

    const vector area = polygonArea(points);
    vector outward(Zero);
    outward[fixedAxis] = side;

    if ((area & outward) < 0)
    {
        result = result.reverseFace();
    }

    return result;
}


void appendUnique
(
    DynamicList<scalar>& values,
    const scalar value,
    const scalar tolerance
)
{
    forAll(values, i)
    {
        if (mag(values[i] - value) <= tolerance)
        {
            return;
        }
    }

    values.append(value);
}


void sortUnique
(
    DynamicList<scalar>& values,
    const scalar tolerance
)
{
    for (label i = 1; i < values.size(); ++i)
    {
        const scalar value = values[i];
        label j = i;

        while (j > 0 && values[j - 1] > value)
        {
            values[j] = values[j - 1];
            --j;
        }

        values[j] = value;
    }

    DynamicList<scalar> unique(values.size());
    forAll(values, i)
    {
        appendUnique(unique, values[i], tolerance);
    }
    values.transfer(unique);
}


label patchForBoxFace
(
    const boundBox& cellBounds,
    const boundBox& cuboidBounds,
    const label axis,
    const label side,
    const scalar tolerance,
    const labelList& regionBySlot
)
{
    const scalar coordinate =
        side < 0 ? cellBounds.min()[axis] : cellBounds.max()[axis];
    const scalar cuboidCoordinate =
        side < 0 ? cuboidBounds.min()[axis] : cuboidBounds.max()[axis];

    if (mag(coordinate - cuboidCoordinate) <= tolerance)
    {
        return regionBySlot[2*axis + (side > 0)];
    }

    return -1;
}


void replaceBoundaryPatches
(
    polyMesh& mesh,
    const wordList& patchNames
)
{
    mesh.clearOut();

    polyBoundaryMesh& boundary =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    boundary.clear();

    forAll(patchNames, patchi)
    {
        dictionary patchDict;
        patchDict.set("type", wallPolyPatch::typeName);
        patchDict.set("nFaces", 0);
        patchDict.set("startFace", mesh.nFaces());

        boundary.push_back
        (
            polyPatch::New
            (
                patchNames[patchi],
                patchDict,
                patchi,
                boundary
            )
        );
    }
}

} // End anonymous namespace


void cuboidClipper::applyDirect
(
    polyMesh& mesh,
    const triSurface& surface,
    const scalar tolerance,
    const label keepZone
)
{
    const List<cuboidPlane> planes = readCuboidPlanes(surface, tolerance);
    const boundBox cuboidBounds(surface.localPoints(), false);

    wordList patchNames(surface.patches().size());
    labelList regionBySlot(6, -1);
    forAll(surface.patches(), regioni)
    {
        patchNames[regioni] = surface.patches()[regioni].name();
    }
    forAll(planes, slot)
    {
        regionBySlot[slot] = planes[slot].region;
    }

    if (keepZone < 0 || keepZone >= mesh.cellZones().size())
    {
        FatalErrorInFunction
            << "Invalid retained cellZone " << keepZone
            << exit(FatalError);
    }

    polyTopoChange meshMod(patchNames.size(), false);

    labelList oldToNewCell(mesh.nCells(), -1);
    List<boundBox> retainedBounds(mesh.nCells());
    DynamicList<scalar> globalCoordinates[3];
    label nRetained = 0;

    for (label celli = 0; celli < mesh.nCells(); ++celli)
    {
        const boundBox originalBounds = cellBounds(mesh, celli, celli);
        boundBox clippedBounds
        (
            max(originalBounds.min(), cuboidBounds.min()),
            min(originalBounds.max(), cuboidBounds.max())
        );

        if
        (
            clippedBounds.max().x() - clippedBounds.min().x() <= tolerance
         || clippedBounds.max().y() - clippedBounds.min().y() <= tolerance
         || clippedBounds.max().z() - clippedBounds.min().z() <= tolerance
        )
        {
            continue;
        }

        retainedBounds[celli] = clippedBounds;
        for (label axis = 0; axis < 3; ++axis)
        {
            appendUnique
            (
                globalCoordinates[axis],
                clippedBounds.min()[axis],
                tolerance
            );
            appendUnique
            (
                globalCoordinates[axis],
                clippedBounds.max()[axis],
                tolerance
            );
        }

        label sourceZone = mesh.cellZones().whichZone(celli);
        if (sourceZone < 0)
        {
            sourceZone = keepZone;
        }

        oldToNewCell[celli] =
            meshMod.addCell(-1, -1, -1, celli, sourceZone);
        ++nRetained;
    }

    for (label axis = 0; axis < 3; ++axis)
    {
        sortUnique(globalCoordinates[axis], tolerance);
    }

    DynamicList<rebuiltFace> rebuiltFaces(6*nRetained);

    for (label celli = 0; celli < mesh.nCells(); ++celli)
    {
        if (oldToNewCell[celli] < 0)
        {
            continue;
        }

        const boundBox clippedBounds = retainedBounds[celli];
        const point cellCentre =
            0.5*(clippedBounds.min() + clippedBounds.max());

        for (label axis = 0; axis < 3; ++axis)
        {
            for (label side = -1; side <= 1; side += 2)
            {
                label firstAxis = -1;
                label secondAxis = -1;
                for (label otherAxis = 0; otherAxis < 3; ++otherAxis)
                {
                    if (otherAxis != axis)
                    {
                        if (firstAxis < 0)
                        {
                            firstAxis = otherAxis;
                        }
                        else
                        {
                            secondAxis = otherAxis;
                        }
                    }
                }

                const scalar fixedValue =
                    side < 0
                  ? clippedBounds.min()[axis]
                  : clippedBounds.max()[axis];

                DynamicList<scalar> firstCoordinates(globalCoordinates[firstAxis]);
                DynamicList<scalar> secondCoordinates(globalCoordinates[secondAxis]);
                for (label i = firstCoordinates.size() - 1; i >= 0; --i)
                {
                    if
                    (
                        firstCoordinates[i]
                            < clippedBounds.min()[firstAxis] - tolerance
                     || firstCoordinates[i]
                            > clippedBounds.max()[firstAxis] + tolerance
                    )
                    {
                        firstCoordinates.remove(i);
                    }
                }

                for (label i = secondCoordinates.size() - 1; i >= 0; --i)
                {
                    if
                    (
                        secondCoordinates[i]
                            < clippedBounds.min()[secondAxis] - tolerance
                     || secondCoordinates[i]
                            > clippedBounds.max()[secondAxis] + tolerance
                    )
                    {
                        secondCoordinates.remove(i);
                    }
                }

                for (label firstI = 0; firstI + 1 < firstCoordinates.size(); ++firstI)
                {
                    for
                    (
                        label secondI = 0;
                        secondI + 1 < secondCoordinates.size();
                        ++secondI
                    )
                    {
                        face current =
                            boxSubFace
                            (
                                meshMod,
                                axis,
                                fixedValue,
                                side,
                                firstAxis,
                                firstCoordinates[firstI],
                                firstCoordinates[firstI + 1],
                                secondAxis,
                                secondCoordinates[secondI],
                                secondCoordinates[secondI + 1],
                                tolerance
                            );

                        pointField currentPoints(current.size());
                        forAll(current, fp)
                        {
                            currentPoints[fp] = meshMod.points()[current[fp]];
                        }

                        vector area = polygonArea(currentPoints);
                        const point faceCentre =
                            0.25
                           *(
                                currentPoints[0]
                              + currentPoints[1]
                              + currentPoints[2]
                              + currentPoints[3]
                            );

                        if ((area & (faceCentre - cellCentre)) < 0)
                        {
                            current = current.reverseFace();
                            forAll(current, fp)
                            {
                                currentPoints[fp] =
                                    meshMod.points()[current[fp]];
                            }
                        }

                        const label patchi =
                            patchForBoxFace
                            (
                                clippedBounds,
                                cuboidBounds,
                                axis,
                                side,
                                tolerance,
                                regionBySlot
                            );

                        bool zoneFlip = false;
                        const label zonei =
                            sourceFaceZone
                            (
                                mesh,
                                currentPoints,
                                tolerance,
                                zoneFlip
                            );

                        const label existing =
                            findRebuiltFace(rebuiltFaces, current);

                        if (existing >= 0)
                        {
                            rebuiltFace& shared = rebuiltFaces[existing];

                            if (shared.patch >= 0 || patchi >= 0)
                            {
                                FatalErrorInFunction
                                    << "Cuboid intersection generated a face "
                                    << "shared by two cells and a boundary "
                                    << "patch."
                                    << exit(FatalError);
                            }

                            if (shared.neighbour >= 0)
                            {
                                FatalErrorInFunction
                                    << "More than two retained cells share "
                                    << "face " << existing << exit(FatalError);
                            }

                            if (shared.zone != zonei)
                            {
                                FatalErrorInFunction
                                    << "Retained cuboid face " << existing
                                    << "has inconsistent faceZone mapping."
                                    << exit(FatalError);
                            }

                            shared.neighbour = oldToNewCell[celli];
                        }
                        else
                        {
                            rebuiltFace newFace;
                            newFace.geometry = current;
                            newFace.owner = oldToNewCell[celli];
                            newFace.neighbour = -1;
                            newFace.patch = patchi;
                            newFace.zone = zonei;
                            newFace.zoneFlip = zoneFlip;
                            rebuiltFaces.append(newFace);
                        }
                    }
                }
            }
        }
    }

    forAll(rebuiltFaces, facei)
    {
        if
        (
            rebuiltFaces[facei].patch < 0
         && rebuiltFaces[facei].neighbour < 0
        )
        {
            FatalErrorInFunction
                << "Retained cuboid cell face " << facei
                << " is neither shared nor on a cuboid patch."
                << exit(FatalError);
        }
    }

    for (label patchi = -1; patchi < patchNames.size(); ++patchi)
    {
        forAll(rebuiltFaces, facei)
        {
            const rebuiltFace& current = rebuiltFaces[facei];
            if (current.patch != patchi)
            {
                continue;
            }

            meshMod.addFace
            (
                current.geometry,
                current.owner,
                current.neighbour,
                -1,
                -1,
                -1,
                false,
                current.patch,
                current.zone,
                current.zoneFlip
            );
        }
    }

    labelList patchMap(patchNames.size(), -1);
    replaceBoundaryPatches(mesh, patchNames);
    autoPtr<mapPolyMesh> morphMap =
        meshMod.changeMesh(mesh, patchMap, false, false, false, false);

    if (morphMap().hasMotionPoints())
    {
        mesh.movePoints(morphMap().preMotionPoints());
    }

    Info<< "Direct cuboid reconstruction retained " << nRetained
        << " cells and generated " << rebuiltFaces.size()
        << " coherent faces across " << patchNames.size()
        << " patches." << endl;
}


void cuboidClipper::apply
(
    polyMesh& mesh,
    const triSurface& surface,
    const meshCutter& cutter,
    const labelList& originalCutCells,
    const labelList& mappedCutCells,
    const List<labelList>& originalCellFaces,
    const List<faceList>& originalCellFaceGeometry,
    const mapPolyMesh& cutMap,
    const List<labelList>& cutCellRegions,
    const scalar tolerance,
    const bool retainedOnly,
    Map<label>& faceRegions,
    Map<label>& retainedCells
)
{
    bool needsClipping = false;

    forAll(cutCellRegions, i)
    {
        if (cutCellRegions[i].size() > 1)
        {
            needsClipping = true;
            break;
        }
    }

    if (!needsClipping)
    {
        return;
    }

    if (!retainedOnly)
    {
        FatalErrorInFunction
            << "Multi-region cuboid clipping currently constructs the "
            << "retained polyhedron only. Re-run with -extractInside; "
            << "keeping the complementary non-planar parent cell is "
            << "explicitly unsupported in this milestone."
            << exit(FatalError);
    }

    const List<cuboidPlane> planes = readCuboidPlanes(surface, tolerance);
    const boundBox surfaceBounds(surface.localPoints(), false);
    const Map<label> oldFaceRegions(faceRegions);
    polyTopoChange meshMod(mesh);
    DynamicList<label> actionFaces;
    DynamicList<label> actionRegions;
    DynamicList<point> addedPoints;
    DynamicList<label> addedPointLabels;
    label nClipped = 0;

    forAll(originalCutCells, i)
    {
        if (cutCellRegions[i].size() <= 1)
        {
            continue;
        }

        const label celli = mappedCutCells[i];
        const label addedCelli = cutter.addedCells()[celli];
        const label interfaceFacei = cutter.addedFaces()[celli];
        const boundBox bounds = cellBounds(mesh, celli, addedCelli);
        DynamicList<label> replacementFaces;
        scalar newArea = 0;
        point retainedSample(Zero);
        for (label axis = 0; axis < 3; ++axis)
        {
            const scalar lower =
                max(bounds.min()[axis], surfaceBounds.min()[axis]);
            const scalar upper =
                min(bounds.max()[axis], surfaceBounds.max()[axis]);

            if (upper - lower <= tolerance)
            {
                FatalErrorInFunction
                    << "No positive-volume retained cuboid sample exists "
                    << "for multi-region cut cell "
                    << originalCutCells[i] << exit(FatalError);
            }

            retainedSample[axis] = 0.5*(lower + upper);
        }

        const scalar cellDistance =
            magSqr(mesh.cellCentres()[celli] - retainedSample);
        const scalar addedDistance =
            magSqr(mesh.cellCentres()[addedCelli] - retainedSample);

        if (mag(cellDistance - addedDistance) <= sqr(tolerance))
        {
            FatalErrorInFunction
                << "Could not identify one retained cuboid side for "
                << "multi-region cut cell " << originalCutCells[i]
                << ". Descendant centres are "
                << mesh.cellCentres()[celli] << " and "
                << mesh.cellCentres()[addedCelli]
                << "; unsupported ambiguous classification."
                << exit(FatalError);
        }

        const label retainedCell =
            cellDistance < addedDistance ? celli : addedCelli;
        retainedCells.insert(retainedCell, 1);
        label nModifiedOuterFaces = 0;

        forAll(originalCellFaces[i], parentFacei)
        {
            const label oldFacei = originalCellFaces[i][parentFacei];
            label currentFacei = -1;

            for (const label facei : mesh.cells()[retainedCell])
            {
                if (cutMap.faceMap()[facei] == oldFacei)
                {
                    if (currentFacei >= 0)
                    {
                        FatalErrorInFunction
                            << "Parent face " << oldFacei
                            << " produced multiple retained descendant "
                            << "faces in cell " << originalCutCells[i]
                            << "; general face subdivision is unsupported."
                            << exit(FatalError);
                    }
                    currentFacei = facei;
                }
            }

            if (currentFacei < 0)
            {
                continue;
            }

            pointField clippedPoints;
            if
            (
                !clippedFacePoints
                (
                    originalCellFaceGeometry[i][parentFacei],
                    mesh.points(),
                    surfaceBounds,
                    tolerance,
                    clippedPoints
                )
            )
            {
                if (mesh.isInternalFace(currentFacei))
                {
                    FatalErrorInFunction
                        << "Retained cuboid clipping would remove internal "
                        << "face " << currentFacei << " from cell "
                        << originalCutCells[i]
                        << "; neighboring-cell face reconstruction is "
                        << "unsupported in this milestone."
                        << exit(FatalError);
                }

                meshMod.removeFace(currentFacei, -1);
                ++nModifiedOuterFaces;
                continue;
            }

            face replacement(4);
            forAll(replacement, fp)
            {
                replacement[fp] = pointFor
                (
                    meshMod,
                    mesh,
                    clippedPoints[fp],
                    tolerance,
                    addedPoints,
                    addedPointLabels
                );
            }

            vector replacementArea = polygonArea(clippedPoints);
            if ((replacementArea & mesh.faceAreas()[currentFacei]) < 0)
            {
                replacement = replacement.reverseFace();
                replacementArea = -replacementArea;
            }

            const label patchi =
                mesh.isInternalFace(currentFacei)
              ? -1
              : mesh.boundaryMesh().whichPatch(currentFacei);
            const label zonei = mesh.faceZones().whichZone(currentFacei);
            bool zoneFlip = false;

            if (zonei >= 0)
            {
                const faceZone& zone = mesh.faceZones()[zonei];
                zoneFlip = zone.flipMap()[zone.whichFace(currentFacei)];
            }

            meshMod.modifyFace
            (
                replacement,
                currentFacei,
                mesh.faceOwner()[currentFacei],
                mesh.isInternalFace(currentFacei)
              ? mesh.faceNeighbour()[currentFacei] : -1,
                false,
                patchi,
                zonei,
                zoneFlip
            );
            ++nModifiedOuterFaces;
        }

        forAll(planes, planei)
        {
            const cuboidPlane& plane = planes[planei];

            if (!cutCellRegions[i].found(plane.region))
            {
                continue;
            }

            if
            (
                plane.coordinate <= bounds.min()[plane.axis] + tolerance
             || plane.coordinate >= bounds.max()[plane.axis] - tolerance
            )
            {
                FatalErrorInFunction
                    << "Cuboid region " << plane.region
                    << " is coplanar with a parent-cell boundary in cell "
                    << originalCutCells[i]
                    << "; vertex/coplanar clipping is unsupported."
                    << exit(FatalError);
            }

            label firstAxis = -1;
            label secondAxis = -1;
            for (label axis = 0; axis < 3; ++axis)
            {
                if (axis != plane.axis)
                {
                    if (firstAxis < 0)
                    {
                        firstAxis = axis;
                    }
                    else
                    {
                        secondAxis = axis;
                    }
                }
            }

            const scalar firstMin =
                max(bounds.min()[firstAxis], surfaceBounds.min()[firstAxis]);
            const scalar firstMax =
                min(bounds.max()[firstAxis], surfaceBounds.max()[firstAxis]);
            const scalar secondMin =
                max(bounds.min()[secondAxis], surfaceBounds.min()[secondAxis]);
            const scalar secondMax =
                min(bounds.max()[secondAxis], surfaceBounds.max()[secondAxis]);

            if
            (
                firstMax - firstMin <= tolerance
             || secondMax - secondMin <= tolerance
            )
            {
                FatalErrorInFunction
                    << "Cuboid region " << plane.region
                    << " has a degenerate intersection with cell "
                    << originalCutCells[i] << exit(FatalError);
            }

            const point p0 =
                makePoint
                (
                    plane.axis, plane.coordinate,
                    firstAxis, firstMin,
                    secondAxis, secondMin
                );
            const point p1 =
                makePoint
                (
                    plane.axis, plane.coordinate,
                    firstAxis, firstMax,
                    secondAxis, secondMin
                );
            const point p2 =
                makePoint
                (
                    plane.axis, plane.coordinate,
                    firstAxis, firstMax,
                    secondAxis, secondMax
                );
            const point p3 =
                makePoint
                (
                    plane.axis, plane.coordinate,
                    firstAxis, firstMin,
                    secondAxis, secondMax
                );

            face replacement(4);
            replacement[0] = pointFor
            (
                meshMod, mesh, p0, tolerance,
                addedPoints, addedPointLabels
            );
            replacement[1] = pointFor
            (
                meshMod, mesh, p1, tolerance,
                addedPoints, addedPointLabels
            );
            replacement[2] = pointFor
            (
                meshMod, mesh, p2, tolerance,
                addedPoints, addedPointLabels
            );
            replacement[3] = pointFor
            (
                meshMod, mesh, p3, tolerance,
                addedPoints, addedPointLabels
            );

            pointField replacementPoints(4);
            replacementPoints[0] = p0;
            replacementPoints[1] = p1;
            replacementPoints[2] = p2;
            replacementPoints[3] = p3;
            vector replacementArea = polygonArea(replacementPoints);

            const point faceCentre =
                0.25
               *(
                    p0
                  + p1
                  + p2
                  + p3
                );

            if
            (
                (replacementArea & (retainedSample - faceCentre)) > 0
            )
            {
                replacement = replacement.reverseFace();
                replacementArea = -replacementArea;
            }

            newArea += mag(replacementArea);
            replacementFaces.append(plane.region);

            const label actionFacei =
                meshMod.addFace
                (
                    replacement,
                    mesh.faceOwner()[interfaceFacei],
                    mesh.faceNeighbour()[interfaceFacei],
                    -1,
                    -1,
                    interfaceFacei,
                    false,
                    -1,
                    -1,
                    false
                );

            actionFaces.append(actionFacei);
            actionRegions.append(plane.region);
        }

        if (replacementFaces.size() < 2)
        {
            FatalErrorInFunction
                << "Multi-region cut cell " << originalCutCells[i]
                << " did not produce at least two region-specific faces."
                << exit(FatalError);
        }

        if (nModifiedOuterFaces == 0 || newArea <= VSMALL)
        {
        FatalErrorInFunction
            << "Cuboid clipping found no valid retained polyhedron "
            << "boundary for multi-region cut cell "
            << originalCutCells[i] << exit(FatalError);
        }

        meshMod.removeFace(interfaceFacei, -1);
        ++nClipped;
    }

    if (nClipped == 0)
    {
        FatalErrorInFunction
            << "Multi-region cut provenance was detected, but no cuboid "
            << "interface could be constructed." << exit(FatalError);
    }

    autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

    if (morphMap().hasMotionPoints())
    {
        mesh.movePoints(morphMap().preMotionPoints());
    }

    updateFaceRegions
    (
        morphMap(),
        oldFaceRegions,
        actionFaces,
        actionRegions,
        faceRegions
    );

    const Map<label> oldRetainedCells(retainedCells);
    retainedCells.clear();

    forAllConstIters(oldRetainedCells, iter)
    {
        const label newCelli =
            morphMap().reverseCellMap()[iter.key()];

        if (newCelli < 0)
        {
            FatalErrorInFunction
                << "Retained cuboid cell " << iter.key()
                << " disappeared during clipping topology change."
                << exit(FatalError);
        }

        retainedCells.insert(newCelli, 1);
    }

    Info<< "Cuboid clipping replaced " << nClipped
        << " multi-region cut interface"
        << (nClipped == 1 ? "" : "s")
        << " with region-specific planar faces." << endl;
}

} // End namespace Foam

// ************************************************************************* //
