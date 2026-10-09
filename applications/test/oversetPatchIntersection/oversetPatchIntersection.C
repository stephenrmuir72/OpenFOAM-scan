/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | www.openfoam.com
    \\  /    A nd           |
     \\/     M anipulation  |
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

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "columnFvMesh.H"
#include "patchFaceIntersection.H"

#include <iostream>

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const vector zeroTolerance(vector::zero);
    const boundBox searchBb(point(0, 0, 0), point(1, 1, 1));

    const triPoints patchTriangle
    (
        point(0, 0, 0.5),
        point(1, 0, 0.5),
        point(0, 1, 0.5)
    );

    const boundBox falsePositiveVoxel
    (
        point(0.75, 0.75, 0.4),
        point(1, 1, 0.6)
    );
    if
    (
        patchFaceIntersection::triangleIntersectsBox
        (
            patchTriangle,
            falsePositiveVoxel,
            zeroTolerance
        )
    )
    {
        std::cerr << "Triangle/voxel test accepted an AABB-only overlap\n";
        return 1;
    }

    const boundBox intersectingVoxel
    (
        point(0.1, 0.1, 0.4),
        point(0.2, 0.2, 0.6)
    );
    if
    (
        !patchFaceIntersection::triangleIntersectsBox
        (
            patchTriangle,
            intersectingVoxel,
            zeroTolerance
        )
    )
    {
        std::cerr << "Triangle/voxel test missed a true intersection\n";
        return 1;
    }

    const triPoints separatedTriangle
    (
        point(0.6, 0.6, 0.5),
        point(1, 0.6, 0.5),
        point(0.6, 1, 0.5)
    );
    if
    (
        patchFaceIntersection::trianglesIntersect
        (
            patchTriangle,
            separatedTriangle,
            0
        )
    )
    {
        std::cerr << "Triangle/triangle test accepted an AABB-only overlap\n";
        return 1;
    }

    const triPoints intersectingTriangle
    (
        point(0.1, 0.1, 0.5),
        point(0.2, 0.1, 0.5),
        point(0.1, 0.2, 0.5)
    );
    if
    (
        !patchFaceIntersection::trianglesIntersect
        (
            patchTriangle,
            intersectingTriangle,
            0
        )
    )
    {
        std::cerr << "Triangle/triangle test missed a true intersection\n";
        return 1;
    }

    const pointField facePoints
    ({
        point(0, 0, 0.5),
        point(1, 0, 0.5),
        point(0, 1, 0.5)
    });
    const face patchFace({0, 1, 2});
    const labelVector nDivs(4, 4, 1);
    const labelVector trackingStrides(1, 4, 16);
    const labelVector inverseStrides(4, 1, 16);

    PackedList<2> trackingVoxels(16);
    trackingVoxels = 0;
    patchFaceIntersection::fillFaceVoxels
    (
        trackingVoxels,
        searchBb,
        nDivs,
        trackingStrides,
        patchFace,
        facePoints,
        zeroTolerance,
        1
    );

    if (trackingVoxels[15] || !trackingVoxels[3])
    {
        std::cerr << "Tracking voxel rasterization did not match face geometry\n";
        return 1;
    }

    PackedList<2> inverseVoxels(16);
    inverseVoxels = 0;
    patchFaceIntersection::fillFaceVoxels
    (
        inverseVoxels,
        searchBb,
        nDivs,
        inverseStrides,
        patchFace,
        facePoints,
        zeroTolerance,
        1
    );

    if (inverseVoxels[15] || !inverseVoxels[12])
    {
        std::cerr << "Inverse-distance voxel rasterization used the wrong layout\n";
        return 1;
    }

    label voxelCount = 0;
    if
    (
        !patchFaceIntersection::checkedVoxelCount
        (
            labelVector(4, 4, 1),
            16,
            voxelCount
        )
     || voxelCount != 16
    )
    {
        std::cerr << "Voxel budget rejected a valid grid\n";
        return 1;
    }

    if
    (
        patchFaceIntersection::checkedVoxelCount
        (
            labelVector(4, 4, 2),
            16,
            voxelCount
        )
     || patchFaceIntersection::checkedVoxelCount
        (
            labelVector(100000, 100000, 100000),
            1024,
            voxelCount
        )
    )
    {
        std::cerr << "Voxel budget accepted an over-budget grid\n";
        return 1;
    }

    labelVector currentDivs(1, 1, 1);
    const labelVector geometricD(1, 1, -1);
    labelVector nextDivs;
    label lastVoxelCount = 0;
    while
    (
        patchFaceIntersection::refineVoxelDivisions
        (
            currentDivs,
            geometricD,
            16,
            nextDivs,
            voxelCount
        )
    )
    {
        currentDivs = nextDivs;
        lastVoxelCount = voxelCount;
    }

    if (currentDivs != labelVector(4, 4, 1) || lastVoxelCount != 16)
    {
        std::cerr << "Voxel refinement did not stop at the configured budget\n";
        return 1;
    }

    if
    (
        patchFaceIntersection::refineVoxelDivisions
        (
            labelVector(1, 1, 1),
            labelVector(-1, -1, -1),
            16,
            nextDivs,
            voxelCount
        )
     || nextDivs != labelVector(1, 1, 1)
     || voxelCount != 1
    )
    {
        std::cerr << "Voxel refinement continued without an active dimension\n";
        return 1;
    }

    const pointField aabbOnlyPatch
    ({
        point(1.4, 0.9, 0.9),
        point(0.9, 1.4, 0.9),
        point(0.9, 0.9, 1.4)
    });
    const boundBox cellBb(mesh.cellBb(0));
    const boundBox patchBb(aabbOnlyPatch, false);

    if (!cellBb.overlaps(patchBb))
    {
        std::cerr << "Cell/patch AABB regression is not exercising overlap\n";
        return 1;
    }

    if
    (
        patchFaceIntersection::intersectsCell
        (
            aabbOnlyPatch,
            mesh,
            0,
            0
        )
    )
    {
        std::cerr << "Cell intersection test accepted an AABB-only overlap\n";
        return 1;
    }

    const pointField patchInsideCell
    ({
        point(0.2, 0.2, 0.2),
        point(0.3, 0.2, 0.2),
        point(0.2, 0.3, 0.2)
    });

    if
    (
        !patchFaceIntersection::intersectsCell
        (
            patchInsideCell,
            mesh,
            0,
            0
        )
    )
    {
        std::cerr << "Cell intersection test missed a patch inside the cell\n";
        return 1;
    }

    const pointField patchCrossingCell
    ({
        point(0.5, -1, 0.5),
        point(0.5, 2, 0.5),
        point(0.5, 0.5, 2)
    });

    if
    (
        !patchFaceIntersection::intersectsCell
        (
            patchCrossingCell,
            mesh,
            0,
            0
        )
    )
    {
        std::cerr << "Cell intersection test missed a crossing patch\n";
        return 1;
    }

    std::cout << "Passed exact overset patch intersection checks\n";
    return 0;
}

// ************************************************************************* //
