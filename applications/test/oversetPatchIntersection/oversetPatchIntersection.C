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
#include "inverseDistanceCellCellStencil.H"
#include "trackingInverseDistanceCellCellStencil.H"
#include "voxelMeshSearch.H"
#include "tetIndices.H"

#include <iostream>
#include <chrono>

using namespace Foam;

class inverseVoxelHarness : public cellCellStencils::inverseDistance
{
public:
    using inverseDistance::overlaps;
    using inverseDistance::markBoundaries;
};

class trackingVoxelHarness : public cellCellStencils::trackingInverseDistance
{
public:
    using trackingInverseDistance::markBoundaries;
};

bool exhaustiveIntersectsCell
(
    const pointField& trianglePoints,
    const fvMesh& mesh,
    const label celli,
    const scalar tolerance
)
{
    DynamicList<triPoints> boundaryTriangles;
    for (const label facei : mesh.cells()[celli])
    {
        for
        (
            label tetPointi = 1;
            tetPointi < mesh.faces()[facei].size() - 1;
            ++tetPointi
        )
        {
            const triPointRef tri
            (
                tetIndices(celli, facei, tetPointi).faceTri(mesh)
            );
            boundaryTriangles.append(triPoints(tri.a(), tri.b(), tri.c()));
        }
    }
    for (label pointi = 0; pointi < trianglePoints.size(); pointi += 3)
    {
        const triPoints tri
        (
            trianglePoints[pointi],
            trianglePoints[pointi + 1],
            trianglePoints[pointi + 2]
        );
        for (const point& p : tri)
        {
            if (mesh.pointInCell(p, celli, polyMesh::CELL_TETS))
            {
                return true;
            }
        }
        for (const triPoints& boundaryTri : boundaryTriangles)
        {
            if
            (
                patchFaceIntersection::trianglesIntersect
                (
                    tri, boundaryTri, tolerance
                )
            )
            {
                return true;
            }
        }
    }
    return false;
}

int main(int argc, char *argv[])
{
    argList::addBoolOption("negativeBase", "Check core negative tet-base fallback");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (args.found("negativeBase"))
    {
        // Simulate the core's failed base-point search on an actual mesh.
        const_cast<labelIOList&>(mesh.tetBasePtIs()) = -1;
        const pointField crossing
        ({
            point(0.5, -2, 0.5), point(0.5, 2, 0.5), point(0.5, 2, 3)
        });
        const pointField separated
        ({
            point(3, -2, 0.5), point(3, 2, 0.5), point(3, 2, 3)
        });
        if (min(mesh.tetBasePtIs()) >= 0)
        {
            std::cerr << "Negative tet-base fixture was not loaded\n";
            return 1;
        }
        for (const point& p : crossing)
        {
            if (mesh.pointInCell(p, 0, polyMesh::CELL_TETS))
            {
                std::cerr << "Negative-base crossing has an inside vertex\n";
                return 1;
            }
        }
        for (const scalar tolerance : {scalar(0), scalar(1e-6)})
        {
            if
            (
                !exhaustiveIntersectsCell(crossing, mesh, 0, tolerance)
             || !patchFaceIntersection::intersectsCell
                    (crossing, mesh, 0, tolerance)
             || patchFaceIntersection::intersectsCell
                    (separated, mesh, 0, tolerance)
            )
            {
                std::cerr << "Negative tet-base fallback disagrees with core\n";
                return 1;
            }
        }
        Info<< "Passed negative tet-base fallback checks" << endl;
        return 0;
    }

    const vector zeroTolerance(vector::zero);
    const boundBox searchBb(point(0, 0, 0), point(1, 1, 1));

    const patchFaceIntersection::triangleSearch emptyIndex(pointField(), 0);
    patchFaceIntersection::triangleSearch::queryStats emptyWork;
    if
    (
        !emptyIndex.query(searchBb, &emptyWork).empty()
     || emptyIndex.intersectsCell(mesh, 0)
     || emptyWork.nodes
     || emptyWork.triangleTests
    )
    {
        std::cerr << "Empty triangle index returned candidates\n";
        return 1;
    }

    for (const scalar boundary : {scalar(0), scalar(0.5), scalar(1)})
    {
        const labelVector divisions(2, 1, 1);
        const treeBoundBox query
        (
            point(boundary, 0, 0),
            point(boundary, 1, 1)
        );
        for (label voxel = 0; voxel < 2; ++voxel)
        {
            PackedList<2> voxels(2);
            voxels = 0;
            voxels[voxel] = 1;
            const bool expected =
                boundary == 0.5 || (boundary == 0 ? voxel == 0 : voxel == 1);
            if
            (
                inverseVoxelHarness::overlaps
                (
                    searchBb, divisions, voxels, query, 1
                ) != expected
             || voxelMeshSearch::overlaps
                (
                    searchBb, divisions, query, voxels, 1u
                ) != expected
            )
            {
                std::cerr << "Voxel query dropped a closed boundary contact\n";
                return 1;
            }
        }
    }

    const boundBox nearWallBox
    (
        point(0.005, 0.1, 0.1), point(0.006, 0.9, 0.9)
    );
    const vector nearWallTolerance(0.001, 0.01, 0.001);
    PackedList<2> inverseNearWall(1), trackingNearWall(1);
    labelList patchTypes(mesh.nCells(), cellCellStencil::OTHER);
    const labelList cellMap(identity(mesh.nCells()));
    labelVector nearWallDivs(1, 1, 1);
    inverseVoxelHarness::markBoundaries
    (
        mesh, nearWallTolerance, nearWallBox, nearWallDivs,
        inverseNearWall, cellMap, patchTypes
    );
    if
    (
        !trackingVoxelHarness::markBoundaries
        (
            0, 16*1024*1024, mesh, nearWallTolerance, nearWallBox, nearWallDivs,
            trackingNearWall, cellMap, patchTypes
        )
     || !inverseNearWall[0]
     || !trackingNearWall[0]
    )
    {
        std::cerr << "Face precheck dropped an anisotropic near-contact\n";
        return 1;
    }

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

    const pointField externalUnequalEdge
    ({
        point(1.00015, 0.2, 0.2),
        point(1.00015, 0.20001, 0.2),
        point(1.00015, 0.2, 0.20001)
    });
    if
    (
        patchFaceIntersection::intersectsCell
        (
            externalUnequalEdge,
            mesh,
            0,
            1e-4
        )
    )
    {
        std::cerr << "External unequal-edge patch pierced the cell\n";
        return 1;
    }

    const pointField externalNearHit
    ({
        point(1.00005, 0.2, 0.2),
        point(1.00005, 0.4, 0.2),
        point(1.00005, 0.2, 0.4)
    });
    if
    (
        !patchFaceIntersection::intersectsCell
        (
            externalNearHit,
            mesh,
            0,
            1e-4
        )
    )
    {
        std::cerr << "Near external patch contact was missed\n";
        return 1;
    }

    const pointField warpedCellPatch
    ({
        point(-0.5, 1.5, 1.25),
        point(1.5, -0.5, 1.25),
        point(1.5, -0.4, 1.25)
    });
    for (const point& p : warpedCellPatch)
    {
        if (mesh.pointInCell(p, 0, polyMesh::CELL_TETS))
        {
            std::cerr << "Warped crossing fixture has an inside vertex\n";
            return 1;
        }
    }
    if
    (
        !patchFaceIntersection::intersectsCell
        (
            warpedCellPatch,
            mesh,
            0,
            0
        )
    )
    {
        std::cerr << "CELL_TETS missed the warped-hexahedron patch\n";
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

    const triPoints unequalEdgeTriangle
    (
        point(0.5, -0.00025, 0.5),
        point(0.5, -0.00015, 0.5),
        point(0.5, -0.00015, 0.50001)
    );
    const triPoints cellBoundaryTriangle
    (
        point(0.5, 0, 0.5),
        point(0.5, 1, 0.5),
        point(0.5, 0, 1)
    );
    if
    (
        patchFaceIntersection::trianglesIntersect
        (
            unequalEdgeTriangle,
            cellBoundaryTriangle,
            1e-4
        )
    )
    {
        std::cerr << "Unequal-edge triangle test accepted a false piercing\n";
        return 1;
    }

    const pointField indexedTriangles
    ({
        point(0.05, 0.05, 0.05),
        point(0.15, 0.05, 0.05),
        point(0.05, 0.15, 0.05),
        point(0.8, 0.8, 0.8),
        point(0.9, 0.8, 0.8),
        point(0.8, 0.9, 0.8)
    });
    const patchFaceIntersection::triangleSearch triangleIndex
    (
        indexedTriangles,
        0
    );
    const labelList localCandidates
    (
        triangleIndex.query
        (
            boundBox(point(0, 0, 0), point(0.2, 0.2, 0.2))
        )
    );
    if
    (
        localCandidates.size() != 1
     || localCandidates[0] != 0
    )
    {
        std::cerr << "Triangle octree returned non-local candidates\n";
        return 1;
    }

    const vector anisotropicTolerance(0.01, 0.001, 0.001);
    const pointField anisotropicContact
    ({
        point(0.2, -0.005, 0.2),
        point(0.4, -0.005, 0.2),
        point(0.2, -0.005, 0.4)
    });
    const patchFaceIntersection::triangleSearch anisotropicSearch
    (
        anisotropicContact,
        anisotropicTolerance
    );
    if
    (
        anisotropicSearch.trianglesOverlapping(mesh.cellBb(0)).size() != 3
     || !anisotropicSearch.intersectsCell(mesh, 0)
    )
    {
        std::cerr << "Anisotropic tolerance omitted a scalar-tolerance hit\n";
        return 1;
    }

    const face anisotropicFace({0, 1, 2});
    PackedList<2> anisotropicVoxels(1);
    anisotropicVoxels = 0;
    patchFaceIntersection::fillFaceVoxels
    (
        anisotropicVoxels,
        searchBb,
        labelVector(1, 1, 1),
        labelVector(1, 1, 1),
        anisotropicFace,
        anisotropicContact,
        anisotropicTolerance,
        1
    );
    if (!anisotropicVoxels[0])
    {
        std::cerr << "Anisotropic tolerance omitted a voxel near-contact\n";
        return 1;
    }

    const label triangleSetSizes[] = {1000, 10000, 100000, 1000000};
    for (const label nTriangles : triangleSetSizes)
    {
        pointField separatedTriangles(3*nTriangles);
        for (label trianglei = 0; trianglei < nTriangles; ++trianglei)
        {
            const scalar x = 2.0*trianglei;
            separatedTriangles[3*trianglei] = point(x, 0, 0);
            separatedTriangles[3*trianglei + 1] = point(x + 0.25, 0, 0);
            separatedTriangles[3*trianglei + 2] = point(x, 0.25, 0);
        }

        const label targetTriangle = nTriangles/2;
        const scalar targetX = 2.0*targetTriangle;
        const auto start = std::chrono::steady_clock::now();
        const patchFaceIntersection::triangleSearch scaledIndex
        (
            separatedTriangles,
            0
        );
        const auto built = std::chrono::steady_clock::now();
        patchFaceIntersection::triangleSearch::queryStats work;
        const boundBox queryBox
        (
            point(targetX - 0.1, -0.1, -0.1),
            point(targetX + 0.4, 0.4, 0.1)
        );
        const labelList candidates
        (
            scaledIndex.query(queryBox, &work)
        );
        if
        (
            candidates.size() != 1
         || candidates[0] != targetTriangle
         || work.nodes > 80
         || work.triangleTests > 20
         || scaledIndex.referenceCount() != nTriangles
         || scaledIndex.nodeCount() <= 1
        )
        {
            std::cerr << "Triangle index query scaled with unrelated "
                << nTriangles << " triangles\n";
            return 1;
        }
        label totalCandidates = 0;
        const auto queried = std::chrono::steady_clock::now();
        for (label queryi = 0; queryi < 5000; ++queryi)
        {
            totalCandidates += scaledIndex.query(queryBox).size();
        }
        const auto done = std::chrono::steady_clock::now();
        if (totalCandidates != 5000)
        {
            std::cerr << "Repeated planar queries changed candidates\n";
            return 1;
        }
        Info<< "Planar triangle index: triangles=" << nTriangles
            << " nodes=" << scaledIndex.nodeCount()
            << " references=" << scaledIndex.referenceCount()
            << " queryNodes=" << work.nodes
            << " triangleTests=" << work.triangleTests
            << " storageBytes=" << scaledIndex.storageBytes()
            << " buildMs="
            << std::chrono::duration<double, std::milli>(built - start).count()
            << " queryUs="
            << std::chrono::duration<double, std::micro>(done - queried).count()/5000
            << endl;
    }

    for (const bool planar : {false, true})
    {
        pointField triangles(3*256);
        for (label i = 0; i < 256; ++i)
        {
            const point base
            (
                0.2*(i % 16), 0.2*(i/16), planar ? 0 : 0.1*(i % 5)
            );
            triangles[3*i] = base;
            triangles[3*i + 1] = base + vector(0.3, 0, 0);
            triangles[3*i + 2] = base + vector(0, 0.3, 0);
        }
        const scalar tolerance = 0.01;
        const patchFaceIntersection::triangleSearch index(triangles, tolerance);
        for (label queryi = 0; queryi < 40; ++queryi)
        {
            const point low(0.08*queryi, 0.07*queryi, -0.01);
            boundBox bb(low, low + vector(0.1, 0.1, 0.2));
            const labelList candidates(index.query(bb));
            bb.grow(tolerance);
            label expected = 0;
            for (label i = 0; i < 256; ++i)
            {
                const boundBox triangleBb
                (
                    triPoints
                    (
                        triangles[3*i], triangles[3*i + 1], triangles[3*i + 2]
                    ).box()
                );
                const bool overlaps = triangleBb.overlaps(bb);
                expected += overlaps;
                if (candidates.found(i) != overlaps)
                {
                    std::cerr << "Indexed query differs from exhaustive bounds\n";
                    return 1;
                }
            }
            if (expected != candidates.size())
            {
                std::cerr << "Triangle query returned duplicate candidates\n";
                return 1;
            }
        }
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

    if
    (
        trackingVoxels[15]
     || !trackingVoxels[0]
     || !trackingVoxels[3]
    )
    {
        std::cerr << "Tracking voxel rasterization did not match face geometry\n";
        return 1;
    }

    const pointField boundaryPlanePoints
    ({
        point(0.5, 0, 0),
        point(0.5, 1, 0),
        point(0.5, 0, 1)
    });
    const face boundaryPlaneFace({0, 1, 2});
    PackedList<2> boundaryPlaneVoxels(2);
    boundaryPlaneVoxels = 0;
    patchFaceIntersection::fillFaceVoxels
    (
        boundaryPlaneVoxels,
        searchBb,
        labelVector(2, 1, 1),
        labelVector(1, 2, 2),
        boundaryPlaneFace,
        boundaryPlanePoints,
        zeroTolerance,
        1
    );
    if (!boundaryPlaneVoxels[0] || !boundaryPlaneVoxels[1])
    {
        std::cerr << "Internal boundary touch was not rasterized closed\n";
        return 1;
    }

    const pointField upperPlanePoints
    ({
        point(1, 0, 0),
        point(1, 1, 0),
        point(1, 0, 1)
    });
    const face upperPlaneFace({0, 1, 2});
    if
    (
        !patchFaceIntersection::triangleIntersectsBox
        (
            triPoints
            (
                upperPlanePoints[0],
                upperPlanePoints[1],
                upperPlanePoints[2]
            ),
            boundBox(point(0.5, 0, 0), point(1, 1, 1)),
            zeroTolerance
        )
    )
    {
        std::cerr << "Upper boundary SAT rejected direct contact\n";
        return 1;
    }
    PackedList<2> upperPlaneVoxels(2);
    upperPlaneVoxels = 0;
    patchFaceIntersection::fillFaceVoxels
    (
        upperPlaneVoxels,
        searchBb,
        labelVector(2, 1, 1),
        labelVector(1, 2, 2),
        upperPlaneFace,
        upperPlanePoints,
        zeroTolerance,
        1
    );
    if (!upperPlaneVoxels[1] || upperPlaneVoxels[0])
    {
        std::cerr << "Upper search-box boundary was not rasterized correctly\n";
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

    if
    (
        inverseVoxels[15]
     || !inverseVoxels[0]
     || !inverseVoxels[12]
    )
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

    const pointField triangleSets[]
    {
        pointField(), externalUnequalEdge, externalNearHit, warpedCellPatch,
        aabbOnlyPatch, patchInsideCell, patchCrossingCell, indexedTriangles
    };
    for (const scalar tolerance : {scalar(0), scalar(1e-4)})
    {
        for (const pointField& triangles : triangleSets)
        {
            const patchFaceIntersection::triangleSearch index(triangles, tolerance);
            if
            (
                index.intersectsCell(mesh, 0)
             != exhaustiveIntersectsCell(triangles, mesh, 0, tolerance)
            )
            {
                std::cerr << "Indexed cell intersection differs from exhaustive predicates\n";
                return 1;
            }
        }
    }

    std::cout << "Passed exact overset patch intersection checks\n";
    return 0;
}

// ************************************************************************* //
