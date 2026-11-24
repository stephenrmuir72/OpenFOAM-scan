/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2009 Frank Bos, TU Delft
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
    The motion is defined by moving patches and static patches.
    The static patches have zero displacement.

    movingIDs_ : point ids from the moving patches (polyMesh meshPoints)
    staticIDs_ : point ids from the static patches (polyMesh meshPoints)
    interpIDs_ : all points ids after removing movingIDs_, staticIDs_
    controlIDs_ : a subset (or all) of the movingIDs_, staticIDs_
        to be used for the RBF interpolations. Can be coarsened
        during construction (uses a simple stride, not geometric coarsening).
    uniqControlIDs_ : (parallel) same as controlIDs_ after topological
        point merging

    For the displacement solver:
    - the RBF interpolator is attached to the globalControlPoints_
      (which means it is parallel consistent) and interpolates
      values with respect to the interpPoints_
    - the solve() stage advances the solid body function, but only
      for the displacements corresponding to the movingIDs_.
    - the curPoints() does the heavy lifting of appplying the displacements
      (from the solve method) and using the RBF interpolator to get all of
      the points.

\*---------------------------------------------------------------------------*/

#include "RBFMotionSolver.H"
#include "emptyPolyPatch.H"
#include "globalIndex.H"
#include "globalMeshData.H"
#include "syncTools.H"
#include "transformField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RBFMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        points0MotionSolver,
        RBFMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::RBFMotionSolver::setControlPoints()
{
    const auto& pts0 = this->points0();

    // Expire old values
    ctrlPointsNumberingPtr_.reset(nullptr);

    if (!UPstream::parRun())
    {
        // serial : trivial
        globalControlPoints_ = UIndirectList<point>(pts0, controlIDs_);

        Info<< "Selected " << controlIDs_.size()
            << " control points" << endl;
        return;
    }

    // parallel:

    // Use a topological merge of control points to eliminate any
    // duplicates across processor boundaries. This is needed regardless
    // of control point coarsening or bounding box filtering.

    // The rank-local id to global index numbering (not needed here)
    labelList pointToGlobal;

    ctrlPointsNumberingPtr_ = mesh().globalData().mergePoints
    (
        controlIDs_,        // [in] rank-local control points
        Map<label>(),       // [in] meshPointMap (unused)
        pointToGlobal,      // [out] ignore here
        uniqControlIDs_     // [out] my rank-local subset
    );
    pointToGlobal.clear();  // discard (unused)

    const label nTotalPoints = ctrlPointsNumberingPtr_->totalSize();

    // Use globalIndex to gather + broadcast.
    {
        const auto& gi = ctrlPointsNumberingPtr_();

        // Concatenate all unique control points (on master)
        gi.gather
        (
            UIndirectList<point>(pts0, uniqControlIDs_),
            globalControlPoints_,
            UPstream::msgType(),
            UPstream::commsTypes::nonBlocking,
            UPstream::worldComm
        );

        // Broadcast the list.
        // - not necessary to broadcast the size separately,
        //   it should be identical to globalIndex totalSize()

        if
        (
            UPstream::master(UPstream::worldComm)
         && (globalControlPoints_.size() != nTotalPoints)
        )
        {
            FatalErrorInFunction
                << "logic error. have " << globalControlPoints_.size()
                << " control points, expecting " << nTotalPoints
                << exit(FatalError);
        }

        if (UPstream::is_subrank(UPstream::worldComm))
        {
            globalControlPoints_.resize_nocopy(nTotalPoints);
        }

        UPstream::broadcast
        (
            globalControlPoints_.data(),
            globalControlPoints_.size(),
            UPstream::worldComm
        );

        Info<< "Selected " << nTotalPoints
            << " unique control points" << endl;
    }
}


void Foam::RBFMotionSolver::makeControlIDs()
{
    // The polyBoundaryMesh
    const auto& pbm = mesh().boundaryMesh();

    // Get the patch id selection
    labelHashSet movingPatchIds(pbm.patchSet(movingPatches_, false));
    labelHashSet staticPatchIds(pbm.patchSet(staticPatches_, false));

    // Only retain non-processor and non-empty types
    {
        const label nNonProcessor = pbm.nNonProcessor();
        const auto real_patch = [&pbm,nNonProcessor](label patchi) -> bool
        {
            if (isType<emptyPolyPatch>(pbm[patchi]))
            {
                return false;
            }
            else
            {
                return (patchi < nNonProcessor);
            }
        };

        movingPatchIds.filterKeys(real_patch);
        staticPatchIds.filterKeys(real_patch);
    }

    // Remove conflicting selections
    if
    (
        labelHashSet duplicates = (movingPatchIds & staticPatchIds);
        (duplicates.size() > 0)
    )
    {
        // Keep as static, remove from moving
        movingPatchIds -= duplicates;

        InfoInFunction
            << "Removed " << duplicates.size()
            << " conflicting moving/static selections from moving patches" << nl
            << "(" << nl;

        for (label patchi : duplicates.sortedToc())
        {
            Info<< "    " << pbm[patchi].name() << nl;
        }
        Info<< ")" << endl;
    }

    // Require something to move
    if (movingPatchIds.empty())
    {
        FatalErrorInFunction
            << "No moving patches selected!" << nl
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Moving patch selections" << nl << "(" << nl;

        for (label patchi : movingPatchIds.sortedToc())
        {
            Info<< "    " << pbm[patchi].name() << nl;
        }
        Info<< ")" << endl;

        Info<< "Static patch selections" << nl << "(" << nl;

        for (label patchi : staticPatchIds.sortedToc())
        {
            Info<< "    " << pbm[patchi].name() << nl;
        }
        Info<< ")" << endl;
    }


    // The point selection states in order of priority
    // (static has higher priority than moving)
    enum pointState : int8_t
    {
        UNVISITED = 0, ZONE_POINT,
        MOVING_POINT, COARSE_MOVING_POINT,
        STATIC_POINT, COARSE_STATIC_POINT
    };

    // The point status
    List<int8_t> markedPoints(mesh().nPoints(), int8_t(UNVISITED));

    // Pass 1: static points (from UNVISITED to STATIC_POINT)
    // - highest priority so done first
    label nStaticPoints = 0;
    for (const auto patchi : staticPatchIds)
    {
        for (const label pointi : pbm[patchi].meshPoints())
        {
            if (markedPoints[pointi] == UNVISITED)
            {
                markedPoints[pointi] = STATIC_POINT;
                ++nStaticPoints;
            }
        }
    }

    // Pass 2: moving points (from UNVISITED to STATIC_POINT)
    // - lower priority so done second
    label nMovingPoints = 0;
    for (const auto patchi : movingPatchIds)
    {
        for (const label pointi : pbm[patchi].meshPoints())
        {
            if (markedPoints[pointi] == UNVISITED)
            {
                markedPoints[pointi] = MOVING_POINT;
                ++nMovingPoints;
            }
        }
    }

    if (UPstream::parRun())
    {
        // Ensure moving/static are consistently marked since
        // they may be touching each other on a processor boundary.
        syncTools::syncPointList
        (
            mesh(),
            markedPoints,
            maxEqOp<int8_t>(),
            int8_t(UNVISITED)
        );

        nMovingPoints = 0;
        nStaticPoints = 0;

        for (const auto status : markedPoints)
        {
            /**/ if (status == STATIC_POINT) { ++nStaticPoints; }
            else if (status == MOVING_POINT) { ++nMovingPoints; }
        }
    }

    // Sizing is now known for moving and static.
    // Globally consistent, no duplicate counting.

    movingIDs_.resize_nocopy(nMovingPoints);
    staticIDs_.resize_nocopy(nStaticPoints);

    // The interpolated points are normallywhatever is left over that is
    // not moving or static. However, with zone motion this will a subset
    // of that.
    // The zoneMotion::pointIDs() are already parallel consistent and can
    // be used without any further synchronization.

    if (zoneMotion::moveAllCells())
    {
        label nInterp
        (
            markedPoints.size() - movingIDs_.size() - staticIDs_.size()
        );

        interpIDs_.resize_nocopy(nInterp);
    }
    else
    {
        // Only select a subset
        label nInterp = 0;

        for (const label pointi : zoneMotion::pointIDs())
        {
            if (markedPoints[pointi] == UNVISITED)
            {
                markedPoints[pointi] = ZONE_POINT;
                ++nInterp;
            }
        }
        interpIDs_.resize_nocopy(nInterp);
    }

    // The point ids by category
    // (note: may contain duplicates across processors)

    {
        label nInterp = 0;
        label nMoving = 0;
        label nStatic = 0;

        // The interpolated points are UNVISITED (all points) or ZONE_POINT
        const int8_t interpPointType =
        (
            zoneMotion::moveAllCells() ? UNVISITED : ZONE_POINT
        );

        forAll(markedPoints, pointi)
        {
            if (markedPoints[pointi] == interpPointType)
            {
                // (UNVISITED | ZONE_POINT)
                interpIDs_[nInterp++] = pointi;
            }
            else if (markedPoints[pointi] == MOVING_POINT)
            {
                movingIDs_[nMoving++] = pointi;
            }
            else if (markedPoints[pointi] == STATIC_POINT)
            {
                staticIDs_[nStatic++] = pointi;
            }
        }

        // Some feedback
        {
            label sizes[3] = { nInterp, nMoving, nStatic };

            Foam::reduce
            (
                sizes, 3,
                sumOp<label>(),
                UPstream::msgType(),
                UPstream::worldComm
            );

            Info<< "Number of interpolated points   : " << sizes[0] << nl
                << "Number of moving boundary points: " << sizes[1] << nl
                << "Number of static boundary points: " << sizes[2] << endl;
        }
    }


    // Pass 3: control points, optionally coarsened.
    // Moving points are first in the list, static (if any) come after.
    // The moving/static ID lists are unique, so any point will only
    // be selected at most once.

    if
    (
        const bool uses_static = includeStaticPatches_ && !staticIDs_.empty();
        (coarseningRatio_ > 1 || uses_static)
    )
    {
        const label ratio = std::max(label(1), coarseningRatio_);

        // The number of control points used is known from the points and
        // the coarsening ratio. Selecting every ratio'th point (starting
        // with the first) gives ceil(n/ratio) points.

        label nControlPoints = ((movingIDs_.size() + ratio - 1) / ratio);
        if (uses_static)
        {
            nControlPoints  += ((staticIDs_.size() + ratio - 1) / ratio);
        }

        controlIDs_.resize_nocopy(nControlPoints);

        // Reset, used as the running insertion index
        nControlPoints = 0;

        // Moving points first
        {
            for (label i = 0; i < movingIDs_.size(); i += ratio)
            {
                const label pointi = movingIDs_[i];
                markedPoints[pointi] = COARSE_MOVING_POINT;
                controlIDs_[nControlPoints++] = pointi;
            }
        }

        // Then static points
        if (uses_static)
        {
            for (label i = 0; i < staticIDs_.size(); i += ratio)
            {
                const label pointi = staticIDs_[i];
                markedPoints[pointi] = COARSE_STATIC_POINT;
                controlIDs_[nControlPoints++] = pointi;
            }
        }
    }
    else
    {
        // Moving points only, no coarsening: just copy them
        controlIDs_ = movingIDs_;
    }

    const pointField& pts0 = this->points0();

    // Reset with zero displacement, updated later...
    motion_.resize_fill(movingIDs_.size(), Foam::zero{});

    // Set the interpolated points
    interpPoints_ = UIndirectList<point>(pts0, interpIDs_);

    // Assign the points0 locations of the unique control points
    // (serial or parallel).
    setControlPoints();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBFMotionSolver::RBFMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    points0MotionSolver(mesh, dict, typeName),
    movingPatches_(coeffDict().get<wordRes>("movingPatches")),
    staticPatches_(coeffDict().get<wordRes>("staticPatches")),
    coarseningRatio_
    (
        coeffDict().get<label>("coarseningRatio")
    ),
    includeStaticPatches_
    (
        coeffDict().get<bool>("includeStaticPatches")
    ),
    frozenInterpolation_
    (
        coeffDict().get<bool>("frozenInterpolation")
    ),
    useFunctionObject_
    (
        coeffDict().getOrDefault("functionObject", false)
    ),
    interpolation_
    (
        coeffDict().subDict("interpolation"),
        globalControlPoints_,
        interpPoints_
    )
{
    if (useFunctionObject_)
    {
        Info<< "Motion driven by functionObject" << endl;
    }
    else
    {
        Info<< "Motion driven by solidBodyMotionFunction" << endl;
        SBMFPtr_ = solidBodyMotionFunction::New(coeffDict(), mesh.time());
    }

    makeControlIDs();

    if (mesh.nSolutionD() == 2)
    {
    Info<<"\nRBF Note: 2D Mesh\n"
        <<" When one of the empty patches is at X,Y, or Z = 0 it's possible\n"
        <<" all of the selected control points X,Y,Z coordinate will be 0.\n"
        <<" This causes all of the entries in a row of the RBF interpolation\n"
        <<" matrix to be zero and the matrix becomes non-invertable. If this\n"
        <<" occurs, translate the mesh slightly in the empty direction."
        << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// Non-default in header (incomplete types)
Foam::RBFMotionSolver::~RBFMotionSolver() = default;


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::RBFMotionSolver::curPoints() const
{
    const auto& pts0 = this->points0();

    // NOTE: unlike other motion solvers (eg, sixDoFRigidBodyMotionSolver)
    // the zonal selection is actually done within makeControlIDs()
    // and is already reflected in interpIDs_.

    // Prepare new points: initially as a zero-displacement field

    auto tcurPoints = tmp<pointField>::New(mesh().nPoints(), Foam::zero{});
    auto& curPoints = tcurPoints.ref();

    // 1. Extract the displacement of the moving points into the
    //    corresponding locations.
    //    Note: motion_ has the same length and order as movingIDs_

    forAll(movingIDs_, i)
    {
        curPoints[movingIDs_[i]] = motion_[i];
    }

    // 2. Displacement for staticIDs_ == zero, but nothing to do here
    //    since displacement == zero when constructed above.

    // The interpolated motion
    tmp<vectorField> tinterpolatedMotion;

    if (UPstream::parRun())
    {
        // Displacements of control points (multi-rank)
        vectorField displacements;

        {
            const auto& gi = ctrlPointsNumberingPtr_();

            gi.gather
            (
                // Gather the unique (rank-local) displacements
                UIndirectList<point>(curPoints, uniqControlIDs_),
                displacements,
                UPstream::msgType(),
                UPstream::commsTypes::nonBlocking,
                UPstream::worldComm
            );

            if
            (
                UPstream::master(UPstream::worldComm)
             && (displacements.size() != gi.totalSize())
            )
            {
                FatalErrorInFunction
                    << "logic error. have " << displacements.size()
                    << " displacements, expecting " << gi.totalSize()
                    << exit(FatalError);
            }

            // Broadcast
            if (UPstream::is_subrank(UPstream::worldComm))
            {
                displacements.resize_nocopy(gi.totalSize());
            }

            UPstream::broadcast
            (
                displacements.data(),
                displacements.size(),
                UPstream::worldComm
            );
        }

        // 3. call interpolation
        tinterpolatedMotion = interpolation_.interpolate(displacements);
    }
    else
    {
        // 3. call interpolation
        tinterpolatedMotion = interpolation_.interpolate
        (
            UIndirectList<point>(curPoints, controlIDs_)
        );
    }
    const auto& interpolatedMotion = tinterpolatedMotion();

    // 4. Insert RBF interpolated motion
    forAll(interpIDs_, i)
    {
        curPoints[interpIDs_[i]] = interpolatedMotion[i];
    }

    // 5. Add old point positions
    curPoints += pts0;

    twoDCorrectPoints(tcurPoints.ref());

    return tcurPoints;
}


void Foam::RBFMotionSolver::solve()
{
    if (useFunctionObject_)
    {
        return;
    }

    const pointField& pts0 = this->points0();

    // This should be a no-op most of the time
    motion_.resize_nocopy(movingIDs_.size());

    // Step 1:
    // calculate the new positions for the moving points

    Foam::transformPoints
    (
        motion_,
        SBMFPtr_().transformation(),
        UIndirectList<point>(pts0, movingIDs_)
    );

    // Step 2:
    // make into relative displacements

    forAll(movingIDs_, i)
    {
        const label pointi = movingIDs_[i];
        motion_[i] -= pts0[pointi];
    }

    // Keep interpolation factors (frozen) or force them to be
    // recalculated on subsequent calls
    if (!frozenInterpolation_)
    {
        interpolation_.clearOut();
    }
}


void Foam::RBFMotionSolver::movePoints(const pointField&)
{
    // No local data to update
}


void Foam::RBFMotionSolver::updateMesh(const mapPolyMesh& mpm)
{
    // Map points0 etc.
    points0MotionSolver::updateMesh(mpm);

    // Zone points may depend on the (changed) mesh
    zoneMotion::rebuildPointIDs(coeffDict(), mesh());

    // Recalculate control point IDs
    makeControlIDs();
}


// ************************************************************************* //
