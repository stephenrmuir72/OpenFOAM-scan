/*---------------------------------------------------------------------------*\
Application
    Test-AMIInterpolation

Description
    Regression checks for faceAreaWeightAMI and cyclicAMI.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "faceAreaWeightAMI.H"
#include "AMICache.H"
#include "PstreamReduceOps.H"
#include "cyclicAMIFvPatchField.H"
#include "fvMesh.H"
#include "volFields.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"

#include <chrono>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>

using namespace Foam;

namespace
{

const scalar transferTol = 1e-12;
const scalar identityTol = 1e-11;


void report
(
    const word& test,
    const word& id,
    const scalar observed,
    const scalar expected,
    const scalar threshold,
    label& nFailures
)
{
    const scalar error = mag(observed - expected);
    const bool finite = std::isfinite(observed) && std::isfinite(error);
    const bool pass = finite && error <= threshold;

    if (Pstream::master())
    {
        Info<< "AMI_ASSERT"
            << " test=" << test
            << " id=" << id
            << " status=" << (pass ? "pass" : "fail")
            << " observed=" << observed
            << " expected=" << expected
            << " error=" << error
            << " threshold=" << threshold
            << endl;
    }

    if (!pass)
    {
        ++nFailures;
    }
}


void reportBool
(
    const word& test,
    const word& id,
    const bool observed,
    const bool expected,
    label& nFailures
)
{
    report
    (
        test,
        id,
        observed ? 1 : 0,
        expected ? 1 : 0,
        0,
        nFailures
    );
}


void reportScalar
(
    const word& test,
    const word& id,
    const scalar observed,
    const scalar expected,
    const scalar relTol,
    const scalar scale,
    label& nFailures
)
{
    report
    (
        test,
        id,
        observed,
        expected,
        relTol*max(mag(expected), scale),
        nFailures
    );
}


void reportVector
(
    const word& test,
    const word& id,
    const vector& observed,
    const vector& expected,
    const scalar relTol,
    const scalar scale,
    label& nFailures
)
{
    reportScalar
    (
        test,
        id,
        mag(observed - expected),
        0,
        relTol,
        scale,
        nFailures
    );
}


void reportUnitTransfer
(
    const word& test,
    const word& prefix,
    faceAreaWeightAMI& ami,
    const primitivePatch& src,
    const primitivePatch& tgt,
    const scalar tolerance,
    label& nFailures
)
{
    scalar sourceCoverageError = 0;
    scalar targetCoverageError = 0;
    for (const scalar coverage : ami.srcWeightsSum())
    {
        sourceCoverageError =
            max(sourceCoverageError, mag(coverage - 1));
    }
    for (const scalar coverage : ami.tgtWeightsSum())
    {
        targetCoverageError =
            max(targetCoverageError, mag(coverage - 1));
    }
    sourceCoverageError =
        returnReduce(sourceCoverageError, maxOp<scalar>());
    targetCoverageError =
        returnReduce(targetCoverageError, maxOp<scalar>());

    const scalarField sourceOne(src.size(), 1);
    const scalarField targetOne(tgt.size(), 1);
    const scalarField mappedToSource(ami.interpolateToSource(targetOne));
    const scalarField mappedToTarget(ami.interpolateToTarget(sourceOne));

    scalar sourceConstantError = 0;
    scalar targetConstantError = 0;
    for (const scalar value : mappedToSource)
    {
        sourceConstantError = max(sourceConstantError, mag(value - 1));
    }
    for (const scalar value : mappedToTarget)
    {
        targetConstantError = max(targetConstantError, mag(value - 1));
    }
    sourceConstantError =
        returnReduce(sourceConstantError, maxOp<scalar>());
    targetConstantError =
        returnReduce(targetConstantError, maxOp<scalar>());

    reportScalar
    (
        test,
        prefix + "_SOURCE_COVERAGE",
        sourceCoverageError,
        0,
        tolerance,
        1,
        nFailures
    );
    reportScalar
    (
        test,
        prefix + "_TARGET_COVERAGE",
        targetCoverageError,
        0,
        tolerance,
        1,
        nFailures
    );
    reportScalar
    (
        test,
        prefix + "_CONSTANT_SOURCE",
        sourceConstantError,
        0,
        tolerance,
        1,
        nFailures
    );
    reportScalar
    (
        test,
        prefix + "_CONSTANT_TARGET",
        targetConstantError,
        0,
        tolerance,
        1,
        nFailures
    );
}


void reportTiming
(
    const word& test,
    const word& phase,
    const label iterations,
    const double totalSeconds
)
{
    if (Pstream::master())
    {
        Info<< "AMI_TIMING"
            << " test=" << test
            << " phase=" << phase
            << " iterations=" << iterations
            << " total_seconds=" << totalSeconds
            << " mean_seconds=" << totalSeconds/iterations
            << endl;
    }
}


bool benchmarkGeometry
(
    const dictionary& amiDict,
    const primitivePatch& src,
    const primitivePatch& tgt,
    const label iterations
)
{
    double constructorSeconds = 0;
    double calculateSeconds = 0;
    double sourceInterpolationSeconds = 0;
    double targetInterpolationSeconds = 0;
    scalar checksum = 0;

    const scalarField sourceValues(src.size(), 1);
    const scalarField targetValues(tgt.size(), 1);

    for (label iteration = 0; iteration < iterations; ++iteration)
    {
        auto start = std::chrono::steady_clock::now();
        faceAreaWeightAMI ami(amiDict);
        auto constructed = std::chrono::steady_clock::now();

        if (!ami.calculate(src, tgt))
        {
            return false;
        }
        auto calculated = std::chrono::steady_clock::now();

        const scalarField mappedToSource
        (
            ami.interpolateToSource(targetValues)
        );
        auto sourceInterpolated = std::chrono::steady_clock::now();

        const scalarField mappedToTarget
        (
            ami.interpolateToTarget(sourceValues)
        );
        auto targetInterpolated = std::chrono::steady_clock::now();

        constructorSeconds +=
            std::chrono::duration<double>(constructed - start).count();
        calculateSeconds +=
            std::chrono::duration<double>(calculated - constructed).count();
        sourceInterpolationSeconds +=
            std::chrono::duration<double>
            (
                sourceInterpolated - calculated
            ).count();
        targetInterpolationSeconds +=
            std::chrono::duration<double>
            (
                targetInterpolated - sourceInterpolated
            ).count();

        if (!mappedToSource.empty())
        {
            checksum += mappedToSource[0];
        }
        if (!mappedToTarget.empty())
        {
            checksum += mappedToTarget[0];
        }
    }

    reportTiming("T01", "constructor", iterations, constructorSeconds);
    reportTiming("T01", "mapping_rebuild", iterations, calculateSeconds);
    reportTiming
    (
        "T01",
        "interpolate_to_source",
        iterations,
        sourceInterpolationSeconds
    );
    reportTiming
    (
        "T01",
        "interpolate_to_target",
        iterations,
        targetInterpolationSeconds
    );

    if (Pstream::master())
    {
        Info<< "AMI_TIMING"
            << " test=T01 phase=checksum"
            << " iterations=" << iterations
            << " value=" << checksum
            << endl;
    }

    return true;
}


void appendQuad
(
    DynamicList<point>& points,
    DynamicList<face>& faces,
    const scalar x0,
    const scalar x1,
    const scalar y0,
    const scalar y1,
    const scalar z,
    const bool reverse
)
{
    const label start = points.size();
    points.append(point(x0, y0, z));
    points.append(point(x1, y0, z));
    points.append(point(x1, y1, z));
    points.append(point(x0, y1, z));

    face f(4);
    forAll(f, i)
    {
        f[i] = start + (reverse ? 3 - i : i);
    }
    faces.append(f);
}


void reverseFaces(DynamicList<face>& faces)
{
    DynamicList<face> reversed(faces.size());
    for (label i = faces.size() - 1; i >= 0; --i)
    {
        reversed.append(faces[i]);
    }
    faces.swap(reversed);
}


void makeGeometry
(
    const word& test,
    const label srcRank,
    const label tgtRank,
    const scalar gap,
    const scalar scale,
    const scalar offset,
    const bool swap,
    const bool unmatchedSource,
    const bool multiTarget,
    DynamicList<point>& srcPts,
    DynamicList<face>& srcFaces,
    DynamicList<point>& tgtPts,
    DynamicList<face>& tgtFaces
)
{
    const bool ownSrc =
        !UPstream::parRun() || UPstream::myProcNo() == srcRank;
    const bool ownTgt =
        !UPstream::parRun() || UPstream::myProcNo() == tgtRank;

    if (ownSrc)
    {
        if (test == "T05" && multiTarget)
        {
            // Source A spans x in [0,2] and overlaps two target faces.
            // Target T2 spans two source faces (A and B).
            appendQuad(srcPts, srcFaces, 0, 2, 0, 1, 0, false);
            appendQuad(srcPts, srcFaces, 2, 3, 0, 1, 0, false);
        }
        else if (test == "T05")
        {
            appendQuad(srcPts, srcFaces, 0, 1, 0, 1, 0, false);
            appendQuad(srcPts, srcFaces, 2, 3, 0, 1, 0, false);

            // Source face with no target overlap (requireMatch negative)
            if (unmatchedSource)
            {
                appendQuad(srcPts, srcFaces, 10, 11, 0, 1, 0, false);
            }
        }
        else
        {
            appendQuad(srcPts, srcFaces, 0, 1, 0, 1, 0, false);
        }
    }

    if (ownTgt)
    {
        if (test == "T05" && multiTarget)
        {
            appendQuad(tgtPts, tgtFaces, 0, 1.5, 0, 1, gap, true);
            appendQuad(tgtPts, tgtFaces, 1.5, 3, 0, 1, gap, true);
        }
        else if (test == "T05")
        {
            appendQuad(tgtPts, tgtFaces, 2, 3, 0, 1, gap, true);
            appendQuad(tgtPts, tgtFaces, 0, 1, 0, 1, gap, true);
        }
        else if (test == "T06")
        {
            tgtPts =
            {
                point(0, 0, gap),
                point(0.0005, 0, gap),
                point(1, 0, gap),
                point(1, 1, gap),
                point(0, 1, gap),
                point(0, 0.0005, gap),
                point(0.0005, 0.0005, gap)
            };
            tgtFaces.append(face{0, 5, 6, 1});
            tgtFaces.append(face{1, 6, 5, 4, 3, 2});
        }
        else if (test == "T12")
        {
            // The empty target patch is the invalid/empty-input fixture.
        }
        else
        {
            tgtPts =
            {
                point(0, 0, gap),
                point(0.5, 0, gap),
                point(1, 0, gap),
                point(0, 1, gap),
                point(0.5, 1, gap),
                point(1, 1, gap)
            };
            tgtFaces.append(face{0, 3, 4, 1});
            tgtFaces.append(face{1, 4, 5, 2});
        }
    }

    if (test == "T10")
    {
        for (point& p : tgtPts)
        {
            p.x() += 0.1;
        }
    }

    if (swap)
    {
        srcPts.swap(tgtPts);
        srcFaces.swap(tgtFaces);
    }

    for (point& p : srcPts)
    {
        p = scale*p + point::uniform(offset);
    }
    for (point& p : tgtPts)
    {
        p = scale*p + point::uniform(offset);
    }
}


void reportScaleMetrics
(
    const word& test,
    const scalar scale,
    const scalar offset,
    const faceAreaWeightAMI& ami,
    const primitivePatch& src,
    const primitivePatch& tgt,
    label& nFailures
)
{
    const scalar x0 = offset;
    const scalar xHalf = offset + 0.5*scale;
    const scalar x1 = offset + scale;
    const scalar y0 = offset;
    const scalar y1 = offset + scale;
    const scalar z = offset;
    const scalar sourceArea = (x1 - x0)*(y1 - y0);

    scalar localSourceArea = 0;
    const scalarField sourceAreas(src.magFaceAreas());
    for (const scalar area : sourceAreas)
    {
        localSourceArea += area;
    }

    scalar localTargetArea = 0;
    const scalarField targetAreas(tgt.magFaceAreas());
    for (const scalar area : targetAreas)
    {
        localTargetArea += area;
    }

    const scalar globalSourceArea =
        returnReduce(localSourceArea, sumOp<scalar>());
    const scalar globalTargetArea =
        returnReduce(localTargetArea, sumOp<scalar>());

    scalar localSourceOverlapArea = 0;
    scalar localTargetOverlapArea = 0;
    forAll(ami.srcWeights(), facei)
    {
        for (const scalar weight : ami.srcWeights()[facei])
        {
            localSourceOverlapArea += sourceAreas[facei]*weight;
        }
    }
    forAll(ami.tgtWeights(), facei)
    {
        for (const scalar weight : ami.tgtWeights()[facei])
        {
            localTargetOverlapArea += targetAreas[facei]*weight;
        }
    }

    const scalar globalSourceOverlapArea =
        returnReduce(localSourceOverlapArea, sumOp<scalar>());
    const scalar globalTargetOverlapArea =
        returnReduce(localTargetOverlapArea, sumOp<scalar>());

    const point expectedCentres[2] =
    {
        point(x0 + 0.5*(xHalf - x0), y0 + 0.5*(y1 - y0), z),
        point(xHalf + 0.5*(x1 - xHalf), y0 + 0.5*(y1 - y0), z)
    };

    scalar localCentroidError = 0;
    const pointListList& sourceCentroids = ami.srcCentroids();
    forAll(sourceCentroids, facei)
    {
        if (sourceCentroids[facei].size() != 2)
        {
            localCentroidError = GREAT;
            continue;
        }

        for (const point& expectedCentre : expectedCentres)
        {
            scalar nearestCentroidError = GREAT;
            for (const point& centroid : sourceCentroids[facei])
            {
                nearestCentroidError = min
                (
                    nearestCentroidError,
                    mag(centroid - expectedCentre)/scale
                );
            }
            localCentroidError = max
            (
                localCentroidError,
                nearestCentroidError
            );
        }
    }
    const scalar globalCentroidError =
        returnReduce(localCentroidError, maxOp<scalar>());

    const scalar areaTolerance = max
    (
        1e-8,
        64*std::numeric_limits<scalar>::epsilon()*mag(offset)/scale
    );
    const scalar centroidTolerance = max
    (
        1e-9,
        32*std::numeric_limits<scalar>::epsilon()*mag(offset)/scale
    );

    report
    (
        test,
        "T11_SOURCE_AREA_SCALE_RATIO",
        globalSourceArea/sourceArea,
        1,
        areaTolerance,
        nFailures
    );
    report
    (
        test,
        "T11_TARGET_AREA_SCALE_RATIO",
        globalTargetArea/sourceArea,
        1,
        areaTolerance,
        nFailures
    );
    report
    (
        test,
        "T11_SOURCE_OVERLAP_AREA_RATIO",
        globalSourceOverlapArea/sourceArea,
        1,
        areaTolerance,
        nFailures
    );
    report
    (
        test,
        "T11_TARGET_OVERLAP_AREA_RATIO",
        globalTargetOverlapArea/sourceArea,
        1,
        areaTolerance,
        nFailures
    );
    report
    (
        test,
        "T11_CENTROID_NORMALIZED_ERROR",
        globalCentroidError,
        0,
        centroidTolerance,
        nFailures
    );

}


void runGeometryTest
(
    const word& test,
    const argList& args,
    const Time& runTime,
    label& nFailures
)
{
    const scalar gap = args.getOrDefault<scalar>("gap", 0);
    const scalar scale = args.getOrDefault<scalar>("scale", 1);
    const scalar offset = args.getOrDefault<scalar>("offset", 0);
    const label srcRank = args.getOrDefault<label>("src-rank", 0);
    const label tgtRank = args.getOrDefault<label>("tgt-rank", 1);

    DynamicList<point> srcPts, tgtPts;
    DynamicList<face> srcFaces, tgtFaces;
    makeGeometry
    (
        test,
        srcRank,
        tgtRank,
        gap,
        scale,
        offset,
        args.found("swap"),
        args.found("unmatched-source"),
        args.found("t05-multi"),
        srcPts,
        srcFaces,
        tgtPts,
        tgtFaces
    );

    if (test == "T05")
    {
        const word order = args.getOrDefault<word>("t05-order", "as-built");
        if (order == "src-reversed" || order == "both-reversed")
        {
            reverseFaces(srcFaces);
        }
        if (order == "tgt-reversed" || order == "both-reversed")
        {
            reverseFaces(tgtFaces);
        }
        if (!(order == "as-built" || order == "src-reversed"
           || order == "tgt-reversed" || order == "both-reversed"))
        {
            FatalErrorInFunction
                << "Unknown -t05-order " << order << nl
                << exit(FatalError);
        }
    }

    const pointField srcPoints(srcPts);
    const pointField tgtPoints(tgtPts);
    primitivePatch src(SubList<face>(srcFaces), srcPoints);
    primitivePatch tgt(SubList<face>(tgtFaces), tgtPoints);

    dictionary amiDict;
    amiDict.add("requireMatch", test != "T12");
    amiDict.add
    (
        "restartUncoveredSourceFace",
        !args.found("no-restart")
    );

    if (test == "T02")
    {
        amiDict.add("cacheSize", 4);
        amiDict.add("rotationAxis", vector(0, 0, 1));
        amiDict.add("rotationCentre", point(Zero));
    }

    if (test == "T07")
    {
        amiDict.add("maxDistance2", args.getOrDefault<scalar>
        (
            "maxDistance2",
            0.01
        ));
    }
    else if (args.found("maxDistance2"))
    {
        amiDict.add("maxDistance2", args.get<scalar>("maxDistance2"));
    }

    faceAreaWeightAMI ami(amiDict);
    const bool calculated = ami.calculate(src, tgt);

    if (test == "T01" && args.found("timing"))
    {
        const label iterations =
            args.getOrDefault<label>("benchmark-iterations", 1);
        if (!benchmarkGeometry(amiDict, src, tgt, iterations))
        {
            ++nFailures;
        }
    }

    reportBool(test, test + "_CALCULATED", calculated, true, nFailures);

    const scalarField srcCoverage(ami.srcWeightsSum());
    const scalarField tgtCoverage(ami.tgtWeightsSum());

    if (test == "T08")
    {
        reportUnitTransfer
        (
            test,
            "T08",
            ami,
            src,
            tgt,
            transferTol,
            nFailures
        );

        const bool expectedDistributed =
            UPstream::parRun() && srcRank != tgtRank;
        reportBool
        (
            test,
            "T08_DISTRIBUTED",
            ami.distributed(),
            expectedDistributed,
            nFailures
        );
        report
        (
            test,
            "T08_SINGLE_PATCH_PROC",
            ami.singlePatchProc(),
            expectedDistributed ? -1 : 0,
            0,
            nFailures
        );
        return;
    }

    if (test == "T11")
    {
        const scalar tolerance = max
        (
            transferTol,
            8*std::numeric_limits<scalar>::epsilon()
           *mag(offset)/scale
        );
        reportUnitTransfer
        (
            test,
            "T11",
            ami,
            src,
            tgt,
            tolerance,
            nFailures
        );
        reportScaleMetrics
        (
            test,
            scale,
            offset,
            ami,
            src,
            tgt,
            nFailures
        );
        return;
    }

    if (test == "T07")
    {
        scalar localSrcError = 0;
        scalar localTgtError = 0;

        if (!srcCoverage.empty())
        {
            for (const scalar coverage : srcCoverage)
            {
                localSrcError = max(localSrcError, mag(coverage - 1));
            }
        }
        if (!tgtCoverage.empty())
        {
            for (const scalar coverage : tgtCoverage)
            {
                localTgtError = max(localTgtError, mag(coverage - 1));
            }
        }

        const scalar srcError =
            returnReduce(localSrcError, maxOp<scalar>());
        const scalar tgtError =
            returnReduce(localTgtError, maxOp<scalar>());

        reportScalar
        (
            test,
            "T07_SOURCE_COVERAGE",
            1 - srcError,
            1,
            transferTol,
            1,
            nFailures
        );
        reportScalar
        (
            test,
            "T07_TARGET_COVERAGE",
            1 - tgtError,
            1,
            transferTol,
            1,
            nFailures
        );
        return;
    }

    if (test == "T12")
    {
        reportScalar
        (
            test,
            "T12_EMPTY_TARGET_COVERAGE",
            srcCoverage.empty() ? 0 : srcCoverage[0],
            0,
            transferTol,
            1,
            nFailures
        );
        return;
    }

    if (test == "T06")
    {
        reportScalar
        (
            test,
            "T06_SOURCE_COVERAGE",
            srcCoverage[0],
            1,
            transferTol,
            1,
            nFailures
        );
        reportScalar
        (
            test,
            "T06_FINE_FACE_COVERAGE",
            tgtCoverage[0],
            1,
            transferTol,
            1,
            nFailures
        );
        return;
    }

    if (test == "T05")
    {
        reportScalar
        (
            test,
            "T05_SOURCE_COVERAGE",
            returnReduce(min(srcCoverage), minOp<scalar>()),
            1,
            transferTol,
            1,
            nFailures
        );
        reportScalar
        (
            test,
            "T05_TARGET_COVERAGE",
            returnReduce(min(tgtCoverage), minOp<scalar>()),
            1,
            transferTol,
            1,
            nFailures
        );

        // Source-face physical mapping. Local addressing and weights are read
        // directly. Target centroids on other ranks are reached through
        // interpolateToSource, so the check is valid in parallel.
        // Multi-target fixture (-t05-multi): source A spans x in [0,2] and
        // overlaps T1 (x in [0,1.5], weight 0.75) and T2 (x in [1.5,3],
        // weight 0.25). Source B spans x in [2,3] and maps only to T2.
        const bool t05Multi = args.found("t05-multi");
        const pointField srcCtr(src.faceCentres());
        const pointField tgtCtr(tgt.faceCentres());
        scalarField tgtX(tgtCtr.size());
        scalarField tgtY(tgtCtr.size());
        forAll(tgtCtr, tgtFacei)
        {
            tgtX[tgtFacei] = tgtCtr[tgtFacei].x();
            tgtY[tgtFacei] = tgtCtr[tgtFacei].y();
        }
        const scalarField mappedX(ami.interpolateToSource(tgtX));
        const scalarField mappedY(ami.interpolateToSource(tgtY));
        const scalarListList& wght = ami.srcWeights();

        bool physical = true;
        forAll(wght, facei)
        {
            const point& s = srcCtr[facei];
            if (t05Multi)
            {
                // Target centroid x is 0.75 for T1 and 2.25 for T2, y is 0.5
                const scalar expectedX =
                    s.x() < 2 ? 0.75*0.75 + 0.25*2.25 : 2.25;
                if (mag(mappedX[facei] - expectedX) > transferTol
                 || mag(mappedY[facei] - 0.5) > transferTol)
                {
                    physical = false;
                }
            }
            else
            {
                // Exactly one target with the same x-y footprint
                label nMatched = 0;
                forAll(wght[facei], j)
                {
                    if (wght[facei][j] > transferTol)
                    {
                        ++nMatched;
                    }
                }
                if (nMatched != 1
                 || mag(mappedX[facei] - s.x()) > transferTol
                 || mag(mappedY[facei] - s.y()) > transferTol)
                {
                    physical = false;
                }
            }
        }
        physical = returnReduceAnd(physical);
        reportBool(test, "T05_PHYSICAL_MAPPING", physical, true, nFailures);
        return;
    }

    if (test == "T10")
    {
        reportScalar
        (
            test,
            "T10_SOURCE_COVERAGE",
            srcCoverage[0],
            0.9,
            transferTol,
            1,
            nFailures
        );
        reportScalar
        (
            test,
            "T10_TARGET_COVERAGE_0",
            tgtCoverage[0],
            1,
            transferTol,
            1,
            nFailures
        );
        reportScalar
        (
            test,
            "T10_TARGET_COVERAGE_1",
            tgtCoverage[1],
            0.8,
            transferTol,
            1,
            nFailures
        );
        return;
    }

    if (test == "T13")
    {
        scalarField targetValues(tgt.size(), 0);
        targetValues[0] = 0;
        targetValues[1] = 2;
        const scalarField mapped(ami.interpolateToSource(targetValues));

        reportScalar
        (
            test,
            "T13_AREA_AVERAGE",
            mapped[0],
            1,
            transferTol,
            2,
            nFailures
        );
        return;
    }

    scalarField sourceOne(src.size(), 1);
    scalarField targetOne(tgt.size(), 1);
    const scalarField toSource(ami.interpolateToSource(targetOne));
    const scalarField toTarget(ami.interpolateToTarget(sourceOne));

    scalar sourceCoverageError = 0;
    scalar sourceConstantError = 0;
    forAll(srcCoverage, i)
    {
        sourceCoverageError =
            max(sourceCoverageError, mag(srcCoverage[i] - 1));
        sourceConstantError =
            max(sourceConstantError, mag(toSource[i] - 1));
    }

    reportScalar
    (
        test,
        "T01_SOURCE_COVERAGE",
        sourceCoverageError,
        0,
        transferTol,
        1,
        nFailures
    );
    reportScalar
    (
        test,
        "T01_CONSTANT_SOURCE",
        sourceConstantError,
        0,
        transferTol,
        1,
        nFailures
    );

    scalar targetCoverageError = 0;
    scalar targetConstantError = 0;
    forAll(tgtCoverage, i)
    {
        targetCoverageError =
            max(targetCoverageError, mag(tgtCoverage[i] - 1));
        targetConstantError =
            max(targetConstantError, mag(toTarget[i] - 1));
    }

    reportScalar
    (
        test,
        "T01_TARGET_COVERAGE",
        targetCoverageError,
        0,
        transferTol,
        1,
        nFailures
    );
    reportScalar
    (
        test,
        "T01_CONSTANT_TARGET",
        targetConstantError,
        0,
        transferTol,
        1,
        nFailures
    );

    bool nonNegative = true;
    for (const auto& faceWeights : ami.srcWeights())
    {
        for (const scalar weight : faceWeights)
        {
            nonNegative = nonNegative && std::isfinite(weight) && weight >= 0;
        }
    }
    reportBool(test, "T01_NONNEGATIVE_WEIGHTS", nonNegative, true, nFailures);

    scalarField sourceValues(src.size(), 3);
    scalarField targetValues(tgt.size(), 0);
    if (targetValues.size() == 2)
    {
        targetValues[0] = 1;
        targetValues[1] = 4;
    }

    const scalarField mappedToSource
    (
        ami.interpolateToSource(targetValues)
    );
    const scalarField mappedToTarget
    (
        ami.interpolateToTarget(sourceValues)
    );

    scalar sourceIntegral = 0;
    scalar targetIntegral = 0;
    const scalarField srcAreas(src.magFaceAreas());
    const scalarField tgtAreas(tgt.magFaceAreas());

    forAll(srcAreas, i)
    {
        sourceIntegral +=
            srcAreas[i]*sourceValues[i]*mappedToSource[i];
    }
    forAll(tgtAreas, i)
    {
        targetIntegral +=
            tgtAreas[i]*mappedToTarget[i]*targetValues[i];
    }

    reportScalar
    (
        test,
        "T01_WEIGHTED_ADJOINT",
        sourceIntegral,
        targetIntegral,
        identityTol,
        max(mag(sourceIntegral), mag(targetIntegral)),
        nFailures
    );
}


void runStateTransitionTest(label& nFailures)
{
    const label nProcs = UPstream::nProcs();
    const label sourceRanks[3] = {0, 0, nProcs >= 4 ? 1 : 0};
    const label targetRanks[3] =
    {
        nProcs >= 2 ? 1 : 0,
        0,
        nProcs >= 4 ? 3 : (nProcs >= 2 ? 1 : 0)
    };
    const scalar offsets[3] = {0, 0.125, 0.25};
    const word prefixes[3] =
    {
        "T09_INITIAL_DISTRIBUTED",
        "T09_LOCAL_AFTER_DISTRIBUTED",
        "T09_DISTRIBUTED_AFTER_LOCAL"
    };

    dictionary amiDict;
    amiDict.add("requireMatch", true);
    amiDict.add("restartUncoveredSourceFace", true);
    faceAreaWeightAMI ami(amiDict);

    for (label state = 0; state < 3; ++state)
    {
        DynamicList<point> srcPts, tgtPts;
        DynamicList<face> srcFaces, tgtFaces;
        makeGeometry
        (
            "T09",
            sourceRanks[state],
            targetRanks[state],
            0,
            1,
            offsets[state],
            false,
            false,
            false,
            srcPts,
            srcFaces,
            tgtPts,
            tgtFaces
        );

        const pointField srcPoints(srcPts);
        const pointField tgtPoints(tgtPts);
        primitivePatch src(SubList<face>(srcFaces), srcPoints);
        primitivePatch tgt(SubList<face>(tgtFaces), tgtPoints);

        ami.upToDate(false);
        const bool calculated = ami.calculate(src, tgt);
        const bool expectedDistributed =
            UPstream::parRun()
         && sourceRanks[state] != targetRanks[state];
        const label expectedSinglePatchProc =
            expectedDistributed ? -1 : sourceRanks[state];
        const bool localMapsPresent =
            ami.hasSrcMap() || ami.hasTgtMap();
        const bool mapsPresent = returnReduceOr(localMapsPresent);

        reportBool
        (
            "T09",
            prefixes[state] + "_CALCULATED",
            calculated,
            true,
            nFailures
        );
        reportBool
        (
            "T09",
            prefixes[state] + "_DISTRIBUTED",
            ami.distributed(),
            expectedDistributed,
            nFailures
        );
        report
        (
            "T09",
            prefixes[state] + "_SINGLE_PATCH_PROC",
            ami.singlePatchProc(),
            expectedSinglePatchProc,
            0,
            nFailures
        );
        reportBool
        (
            "T09",
            prefixes[state] + "_MAPS_PRESENT",
            mapsPresent,
            expectedDistributed,
            nFailures
        );
        reportUnitTransfer
        (
            "T09",
            prefixes[state],
            ami,
            src,
            tgt,
            transferTol,
            nFailures
        );
    }
}


void runCacheTest
(
    const Time& runTime,
    label& nFailures,
    const bool timing
)
{
    DynamicList<point> srcPts, tgtPts;
    DynamicList<face> srcFaces, tgtFaces;
    makeGeometry
    (
        "T02",
        0,
        1,
        0,
        1,
        0,
        false,
        false,
        false,
        srcPts,
        srcFaces,
        tgtPts,
        tgtFaces
    );

    const pointField srcPoints(srcPts);
    const pointField tgtPoints(tgtPts);
    primitivePatch src(SubList<face>(srcFaces), srcPoints);
    primitivePatch tgt(SubList<face>(tgtFaces), tgtPoints);

    dictionary amiDict;
    amiDict.add("requireMatch", true);
    amiDict.add("cacheSize", 4);
    amiDict.add("rotationAxis", vector(0, 0, 1));
    amiDict.add("rotationCentre", point(Zero));

    faceAreaWeightAMI ami(amiDict);
    ami.calculate(src, tgt);

    AMICache& cache = ami.cache();
    cache.addToCache(ami, point(1, 0, 0));
    ami.srcWeights()[0][0] = 1;
    ami.srcWeights()[0][1] = 0;
    cache.addToCache(ami, point(0, 1, 0));

    auto restoreStarted = std::chrono::steady_clock::now();
    const bool restored = ami.restoreCache(point(1, 1, 0), runTime);
    auto restoreFinished = std::chrono::steady_clock::now();
    if (timing)
    {
        reportTiming
        (
            "T02",
            "cache_restore",
            1,
            std::chrono::duration<double>
            (
                restoreFinished - restoreStarted
            ).count()
        );
    }
    reportBool("T02", "T02_CACHE_RESTORED", restored, true, nFailures);

    constexpr bool weightIsScalar =
        std::is_same_v<decltype(cache.weight()), scalar>;
    reportBool
    (
        "T02",
        "T02_CACHE_WEIGHT_TYPE",
        weightIsScalar,
        true,
        nFailures
    );

    const scalar weight = static_cast<scalar>(cache.weight());
    reportScalar
    (
        "T02",
        "T02_CACHE_FRACTION",
        weight,
        0.5,
        transferTol,
        1,
        nFailures
    );

    scalarField marker(2);
    marker[0] = 0;
    marker[1] = 2;
    auto interpolationStarted = std::chrono::steady_clock::now();
    const scalarField midpoint(ami.interpolateToSource(marker));
    auto interpolationFinished = std::chrono::steady_clock::now();
    if (timing)
    {
        reportTiming
        (
            "T02",
            "cached_interpolation",
            1,
            std::chrono::duration<double>
            (
                interpolationFinished - interpolationStarted
            ).count()
        );
    }
    reportScalar
    (
        "T02",
        "T02_CACHE_MIDPOINT",
        midpoint[0],
        1.5,
        transferTol,
        2,
        nFailures
    );
}


point cachePoint(const scalar theta)
{
    return point(std::cos(theta), std::sin(theta), 0);
}


bool setCacheOperator
(
    faceAreaWeightAMI& ami,
    const bool toSource,
    const bool selectSecond
)
{
    scalarListList& weights =
        toSource ? ami.srcWeights() : ami.tgtWeights();
    scalarField& weightsSum =
        toSource ? ami.srcWeightsSum() : ami.tgtWeightsSum();

    bool valid = weights.size() == weightsSum.size();
    label localRows = 0;

    forAll(weights, facei)
    {
        if (weights[facei].size() != 2)
        {
            valid = false;
            continue;
        }

        weights[facei][0] = selectSecond ? 0 : 1;
        weights[facei][1] = selectSecond ? 1 : 0;
        weightsSum[facei] = 1;
        ++localRows;
    }

    const label globalRows =
        returnReduce(localRows, sumOp<label>());
    const label anyInvalid =
        returnReduce(valid ? 0 : 1, maxOp<label>());

    return globalRows > 0 && anyInvalid == 0;
}


scalar cacheResultSum(const scalarField& result)
{
    scalar localSum = 0;
    for (const scalar value : result)
    {
        localSum += value;
    }
    return returnReduce(localSum, sumOp<scalar>());
}


bool cacheResultFinite(const scalarField& result)
{
    bool localFinite = true;
    for (const scalar value : result)
    {
        localFinite = localFinite && std::isfinite(value);
    }
    return returnReduce(localFinite ? 1 : 0, minOp<label>()) == 1;
}


scalar cacheResultError
(
    const scalarField& result,
    const scalar expected
)
{
    scalar localError = 0;
    for (const scalar value : result)
    {
        localError = max(localError, mag(value - expected));
    }
    return returnReduce(localError, maxOp<scalar>());
}


scalarField cacheInterpolate
(
    faceAreaWeightAMI& ami,
    const bool toSource,
    const scalarField& values
)
{
    if (toSource)
    {
        return scalarField(ami.interpolateToSource(values));
    }
    return scalarField(ami.interpolateToTarget(values));
}


scalar checkCacheQuery
(
    const word& direction,
    const word& query,
    faceAreaWeightAMI& ami,
    const Time& runTime,
    const point& queryPoint,
    const scalar expectedWeight,
    const bool exact,
    const scalarField& variableValues,
    const scalarField& constantValues,
    const scalar expectedResult,
    label& nFailures
)
{
    const word prefix("T02_" + direction + "_" + query);
    const AMICache& cache = ami.cache();
    const bool restored = ami.restoreCache(queryPoint, runTime);

    const label restoredRanks =
        returnReduce(restored ? 1 : 0, minOp<label>());
    reportBool
    (
        "T02",
        prefix + "_RESTORED",
        restoredRanks == 1,
        true,
        nFailures
    );

    const bool localState =
        cache.index0() >= 0
     && (
            exact
          ? cache.index1() == -1 && !cache.applyInterpolate()
          : cache.index1() >= 0
         && cache.index1() != cache.index0()
         && cache.applyInterpolate()
        );
    const label stateRanks =
        returnReduce(localState ? 1 : 0, minOp<label>());
    const label minIndex0 =
        returnReduce(cache.index0(), minOp<label>());
    const label maxIndex0 =
        returnReduce(cache.index0(), maxOp<label>());
    const label minIndex1 =
        returnReduce(cache.index1(), minOp<label>());
    const label maxIndex1 =
        returnReduce(cache.index1(), maxOp<label>());
    reportBool
    (
        "T02",
        prefix + "_INDEX_STATE",
        stateRanks == 1
     && minIndex0 == maxIndex0
     && minIndex1 == maxIndex1,
        true,
        nFailures
    );

    const scalar minWeight =
        returnReduce(cache.weight(), minOp<scalar>());
    const scalar maxWeight =
        returnReduce(cache.weight(), maxOp<scalar>());
    reportScalar
    (
        "T02",
        prefix + "_WEIGHT",
        max(mag(minWeight - expectedWeight), mag(maxWeight - expectedWeight)),
        0,
        transferTol,
        1,
        nFailures
    );
    const bool localWeightValid =
        std::isfinite(cache.weight())
     && cache.weight() >= 0
     && cache.weight() <= 1;
    reportBool
    (
        "T02",
        prefix + "_WEIGHT_BOUNDS",
        returnReduce(localWeightValid ? 1 : 0, minOp<label>()) == 1,
        true,
        nFailures
    );

    const scalarField variableResult
    (
        cacheInterpolate(ami, direction == "SOURCE", variableValues)
    );
    const scalarField constantResult
    (
        cacheInterpolate(ami, direction == "SOURCE", constantValues)
    );

    const label globalResultCount =
        returnReduce(label(variableResult.size()), sumOp<label>());
    reportBool
    (
        "T02",
        prefix + "_RESULT_PRESENT",
        globalResultCount > 0,
        true,
        nFailures
    );
    reportScalar
    (
        "T02",
        prefix + "_RESULT",
        cacheResultError(variableResult, expectedResult),
        0,
        transferTol,
        1,
        nFailures
    );
    reportScalar
    (
        "T02",
        prefix + "_CONSTANT",
        cacheResultError(constantResult, 7),
        0,
        transferTol,
        7,
        nFailures
    );
    reportBool
    (
        "T02",
        prefix + "_FINITE",
        cacheResultFinite(variableResult)
     && cacheResultFinite(constantResult),
        true,
        nFailures
    );

    return cacheResultSum(variableResult);
}


void runCacheFractionDirection
(
    const Time& runTime,
    const label srcRank,
    const label tgtRank,
    const bool toSource,
    label& nFailures
)
{
    const word direction(toSource ? "SOURCE" : "TARGET");
    DynamicList<point> srcPts, tgtPts;
    DynamicList<face> srcFaces, tgtFaces;
    makeGeometry
    (
        "T02",
        srcRank,
        tgtRank,
        0,
        1,
        0,
        !toSource,
        false,
        false,
        srcPts,
        srcFaces,
        tgtPts,
        tgtFaces
    );

    const pointField srcPoints(srcPts);
    const pointField tgtPoints(tgtPts);
    primitivePatch src(SubList<face>(srcFaces), srcPoints);
    primitivePatch tgt(SubList<face>(tgtFaces), tgtPoints);

    dictionary amiDict;
    amiDict.add("requireMatch", true);
    amiDict.add("cacheSize", 4);
    amiDict.add("nThetaStencilMax", 4);
    amiDict.add("rotationAxis", vector(0, 0, 1));
    amiDict.add("rotationCentre", point(Zero));

    faceAreaWeightAMI ami(amiDict);
    const bool calculated = ami.calculate(src, tgt);
    const bool allCalculated =
        returnReduce(calculated ? 1 : 0, minOp<label>()) == 1;
    reportBool
    (
        "T02",
        "T02_" + direction + "_CALCULATED",
        allCalculated,
        true,
        nFailures
    );
    if (!allCalculated)
    {
        return;
    }

    const scalar pi = std::acos(-1.0);
    const scalar lowerAngle = 0;
    const scalar upperAngle = 0.5*pi;
    scalarField variableValues(toSource ? tgt.size() : src.size(), 0);
    if (variableValues.size() == 2)
    {
        variableValues[1] = 2;
    }
    scalarField constantValues(variableValues.size(), 7);

    const bool lowerSet = setCacheOperator(ami, toSource, false);
    ami.addToCache(cachePoint(lowerAngle));
    const bool upperSet = setCacheOperator(ami, toSource, true);
    ami.addToCache(cachePoint(upperAngle));
    const bool operatorsSet =
        returnReduce(lowerSet && upperSet ? 1 : 0, minOp<label>()) == 1;
    reportBool
    (
        "T02",
        "T02_" + direction + "_OPERATORS_SET",
        operatorsSet,
        true,
        nFailures
    );
    if (!operatorsSet)
    {
        return;
    }

    const scalar lowerResult =
        checkCacheQuery
        (
            direction,
            "Q000",
            ami,
            runTime,
            cachePoint(lowerAngle),
            0,
            true,
            variableValues,
            constantValues,
            toSource ? 2 : 0,
            nFailures
        );
    const scalar upperResult =
        checkCacheQuery
        (
            direction,
            "Q1000",
            ami,
            runTime,
            cachePoint(upperAngle),
            0,
            true,
            variableValues,
            constantValues,
            toSource ? 0 : 2,
            nFailures
        );
    reportScalar
    (
        "T02",
        "T02_" + direction + "_ENDPOINT_SEPARATION",
        mag(lowerResult - upperResult),
        2,
        transferTol,
        2,
        nFailures
    );

    for (const scalar fraction : {scalar(0.25), scalar(0.5), scalar(0.75)})
    {
        const word query
        (
            fraction == 0.25 ? "Q0250"
          : fraction == 0.5  ? "Q0500"
                             : "Q0750"
        );
        const scalar expectedResult =
            (1 - fraction)*lowerResult + fraction*upperResult;
        checkCacheQuery
        (
            direction,
            query,
            ami,
            runTime,
            cachePoint(fraction*upperAngle),
            fraction,
            false,
            variableValues,
            constantValues,
            expectedResult,
            nFailures
        );
    }

    const scalar wrapLowerAngle = 1.5*pi;
    const scalar wrapUpperAngle = 0;
    ami.cache().clearCache();
    const bool wrapLowerSet =
        setCacheOperator(ami, toSource, toSource);
    ami.addToCache(cachePoint(wrapLowerAngle));
    const bool wrapUpperSet = setCacheOperator(ami, toSource, !toSource);
    ami.addToCache(cachePoint(wrapUpperAngle));
    const bool wrapOperatorsSet =
        returnReduce(wrapLowerSet && wrapUpperSet ? 1 : 0, minOp<label>()) == 1;
    reportBool
    (
        "T02",
        "T02_" + direction + "_WRAP_OPERATORS_SET",
        wrapOperatorsSet,
        true,
        nFailures
    );
    if (!wrapOperatorsSet)
    {
        return;
    }

    const scalar wrapLowerResult =
        checkCacheQuery
        (
            direction,
            "W000",
            ami,
            runTime,
            cachePoint(wrapLowerAngle),
            0,
            true,
            variableValues,
            constantValues,
            0,
            nFailures
        );
    const scalar wrapUpperResult =
        checkCacheQuery
        (
            direction,
            "W1000",
            ami,
            runTime,
            cachePoint(wrapUpperAngle),
            0,
            true,
            variableValues,
            constantValues,
            2,
            nFailures
        );
    const scalar wrapFraction = 0.5;
    const scalar wrapExpected =
        (1 - wrapFraction)*wrapLowerResult
      + wrapFraction*wrapUpperResult;
    checkCacheQuery
    (
        direction,
        "W0500",
        ami,
        runTime,
        cachePoint(1.75*pi),
        wrapFraction,
        false,
        variableValues,
        constantValues,
        wrapExpected,
        nFailures
    );
}


void runCacheFractionTest
(
    const Time& runTime,
    const argList& args,
    label& nFailures
)
{
    constexpr bool weightIsScalar =
        std::is_same_v<decltype(std::declval<AMICache>().weight()), scalar>;
    reportBool
    (
        "T02",
        "T02_CACHE_WEIGHT_TYPE",
        weightIsScalar,
        true,
        nFailures
    );

    const label srcRank = args.getOrDefault<label>("src-rank", 0);
    const label tgtRank = args.getOrDefault<label>("tgt-rank", 1);
    runCacheFractionDirection
    (
        runTime,
        srcRank,
        tgtRank,
        true,
        nFailures
    );
    runCacheFractionDirection
    (
        runTime,
        srcRank,
        tgtRank,
        false,
        nFailures
    );
}


// Independent Rodrigues rotation used as the T03 oracle.  Deliberately
// not the production cyclicAMI transform code.
vector rotateOracle
(
    const vector& v,
    const vector& axis,
    const scalar angleDeg
)
{
    const scalar theta = angleDeg*constant::mathematical::pi/180.0;
    const vector k = axis/mag(axis);

    return
        v*std::cos(theta)
      + (k ^ v)*std::sin(theta)
      + k*(k & v)*(1 - std::cos(theta));
}


//- Check each receiver face of a cyclicAMI neighbour field against the
//  low-weight oracle. Faces below lowWeightCorrection (strict '<', as in
//  AMIInterpolation::weightedSum) must hold the receiver default; all other
//  faces must hold the rotated donor value.
template<class Type>
void checkLowWeightFaces
(
    const word& test,
    const word& patchTag,
    const cyclicAMIFvPatch& amiPatch,
    const tmp<Field<Type>>& tNeighbour,
    const Type& expectedDefault,
    const Type& expectedDonor,
    label& nFailures
)
{
    const cyclicAMIPolyPatch& cpp = amiPatch.cyclicAMIPatch();
    const scalarField& weightsSum = cpp.weightsSum();
    const scalar lowWeight =
    (
        cpp.owner()
      ? cpp.AMI().lowWeightCorrection()
      : cpp.neighbPatch().AMI().lowWeightCorrection()
    );
    const Field<Type>& neighbour = tNeighbour();

    scalar maxDeviation = 0;
    label nLow = 0;
    label nDonor = 0;

    if (neighbour.size() == weightsSum.size())
    {
        forAll(neighbour, facei)
        {
            const bool low = weightsSum[facei] < lowWeight;
            const Type& expected = low ? expectedDefault : expectedDonor;
            maxDeviation =
                max(maxDeviation, mag(neighbour[facei] - expected));
            if (low)
            {
                ++nLow;
            }
            else
            {
                ++nDonor;
            }
        }
    }
    else
    {
        maxDeviation = std::numeric_limits<scalar>::max();
    }

    reduce(maxDeviation, maxOp<scalar>());
    reduce(nLow, sumOp<label>());
    reduce(nDonor, sumOp<label>());
    label nFaces = neighbour.size();
    reduce(nFaces, sumOp<label>());

    reportScalar
    (
        test,
        "T04_FACE_FRAME_" + patchTag,
        maxDeviation,
        0,
        identityTol,
        1,
        nFailures
    );
    reportBool
    (
        test,
        "T04_FACES_CHECKED_" + patchTag,
        nFaces > 0 && nLow + nDonor == nFaces,
        true,
        nFailures
    );
}


//- Rotated donor oracle, written out explicitly (R & T & R.T()) so the test
//  does not reuse the transform code under test.
inline vector rotatedDonor(const tensor& R, const vector& v)
{
    return R & v;
}

inline scalar rotatedDonor(const tensor&, const scalar s)
{
    return s;
}

inline tensor rotatedDonor(const tensor& R, const tensor& T)
{
    return R & T & R.T();
}

inline symmTensor rotatedDonor(const tensor& R, const symmTensor& T)
{
    return symm(R & tensor(T) & R.T());
}

template<class Type>
void checkLowWeightField
(
    const word& test,
    const word& fieldName,
    GeometricField<Type, fvPatchField, volMesh>& field,
    const label patchA,
    const label patchB,
    const cyclicAMIFvPatch& amiPatchA,
    const cyclicAMIFvPatch& amiPatchB,
    const Type& value,
    label& nFailures
)
{
    // Receiver-side low-weight defaults stay at the unrotated value. Each
    // side's donor expectation uses that patch's own forwardT, because the
    // two sides' calcTransforms may choose opposite rotation signs.
    checkLowWeightFaces
    (
        test,
        fieldName + "_AMI1",
        amiPatchA,
        field.boundaryField()[patchA].patchNeighbourField(),
        value,
        rotatedDonor(amiPatchA.cyclicAMIPatch().forwardT()[0], value),
        nFailures
    );

    checkLowWeightFaces
    (
        test,
        fieldName + "_AMI2",
        amiPatchB,
        field.boundaryField()[patchB].patchNeighbourField(),
        value,
        rotatedDonor(amiPatchB.cyclicAMIPatch().forwardT()[0], value),
        nFailures
    );
}


void runMeshTest
(
    const word& test,
    const Time& runTime,
    label& nFailures,
    const bool timing
)
{
    fvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    if (test == "T03")
    {
        auto& boundary =
            const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
        for (auto& patch : boundary)
        {
            patch.changeTopology();
        }

        polyTopoChange topo(mesh);
        topo.movePoints(mesh.points());
        for (auto& patch : boundary)
        {
            patch.setTopology(topo);
        }

        const auto rebuildStarted = std::chrono::steady_clock::now();
        autoPtr<mapPolyMesh> map = topo.changeMesh(mesh, true);
        mesh.updateMesh(map());

        if (map().hasMotionPoints())
        {
            mesh.movePoints(map().preMotionPoints());
        }
        const auto rebuildFinished = std::chrono::steady_clock::now();
        if (timing)
        {
            reportTiming
            (
                "T03",
                "topology_rebuild",
                1,
                std::chrono::duration<double>
                (
                    rebuildFinished - rebuildStarted
                ).count()
            );
        }
    }

    volVectorField U
    (
        IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );
    U.correctBoundaryConditions();

    label patchA = -1;
    label patchB = -1;
    forAll(mesh.boundary(), patchi)
    {
        if (mesh.boundary()[patchi].name() == "AMI1")
        {
            patchA = patchi;
        }
        else if (mesh.boundary()[patchi].name() == "AMI2")
        {
            patchB = patchi;
        }
    }

    reportBool(test, test + "_AMI_PATCHES_PRESENT", patchA >= 0 && patchB >= 0, true, nFailures);
    if (patchA < 0 || patchB < 0)
    {
        return;
    }

    const auto& amiPatchA =
        refCast<const cyclicAMIFvPatch>(mesh.boundary()[patchA]);
    const auto& amiPatchB =
        refCast<const cyclicAMIFvPatch>(mesh.boundary()[patchB]);

    if (test == "T03")
    {
        // Optional variant description, generated per run by the runner.
        // Absent file reproduces the 90 degree fixture.
        const IOdictionary t03Variant
        (
            IOobject
            (
                "t03Variant",
                runTime.system(),
                runTime,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        );

        const bool rotational =
            t03Variant.lookupOrDefault<word>("transform", word("rotational"))
         == "rotational";
        const vector rotationAxis
        (
            t03Variant.lookupOrDefault<vector>
            (
                "rotationAxis",
                vector(0, 0, 1)
            )
        );
        const vector rotationCentre
        (
            t03Variant.lookupOrDefault<vector>
            (
                "rotationCentre",
                vector(0, 0, 0)
            )
        );
        const scalar rotationAngle =
            rotational
          ? t03Variant.lookupOrDefault<scalar>("rotationAngle", 90.0)
          : 0.0;
        const vector separationVector
        (
            t03Variant.lookupOrDefault<vector>
            (
                "separationVector",
                vector(0, 0, 0)
            )
        );

        // Owner face normal in the fixture; neighbour is the rotated image
        // with its outward normal reversed.
        const vector ownerNormal(0, -1, 0);
        const vector expectedNeighbourArea =
            -rotateOracle(ownerNormal, rotationAxis, rotationAngle);

        const polyPatch& ownerPatch = amiPatchA.patch();
        const polyPatch& nbrPatch = amiPatchB.patch();

        // Face counts are summed over ranks, so every rank takes the same
        // branch.  Each raw/stored quantity is a sum over local faces, which
        // equals the single global face because the count is checked above.
        label nOwnerFaces = ownerPatch.size();
        label nNbrFaces = nbrPatch.size();
        reduce(nOwnerFaces, sumOp<label>());
        reduce(nNbrFaces, sumOp<label>());
        report(test, "T03_OWNER_FACE_EXISTS", nOwnerFaces, 1, 0, nFailures);
        report(test, "T03_NEIGHBOUR_FACE_EXISTS", nNbrFaces, 1, 0, nFailures);

        vector ownerRaw(0, 0, 0);
        vector ownerCentre(0, 0, 0);
        forAll(ownerPatch.localFaces(), fi)
        {
            ownerRaw +=
                ownerPatch.localFaces()[fi].areaNormal
                (
                    ownerPatch.localPoints()
                );
            ownerCentre +=
                ownerPatch.localFaces()[fi].centre(ownerPatch.localPoints());
        }
        reduce(ownerRaw, sumOp<vector>());
        reduce(ownerCentre, sumOp<vector>());

        vector nbrRaw(0, 0, 0);
        vector nbrCentre(0, 0, 0);
        forAll(nbrPatch.localFaces(), fi)
        {
            nbrRaw +=
                nbrPatch.localFaces()[fi].areaNormal
                (
                    nbrPatch.localPoints()
                );
            nbrCentre +=
                nbrPatch.localFaces()[fi].centre(nbrPatch.localPoints());
        }
        reduce(nbrRaw, sumOp<vector>());
        reduce(nbrCentre, sumOp<vector>());

        vector ownerStored(0, 0, 0);
        forAll(amiPatchA.Sf(), fi)
        {
            ownerStored += amiPatchA.Sf()[fi];
        }
        reduce(ownerStored, sumOp<vector>());

        vector nbrStored(0, 0, 0);
        forAll(amiPatchB.Sf(), fi)
        {
            nbrStored += amiPatchB.Sf()[fi];
        }
        reduce(nbrStored, sumOp<vector>());

        const vector expectedNbrCentre =
            rotateOracle
            (
                ownerCentre - rotationCentre,
                rotationAxis,
                rotationAngle
            )
          + rotationCentre
          + separationVector;

        reportVector
        (
            test,
            "T03_OWNER_AREA_GEOMETRY",
            ownerRaw,
            ownerNormal,
            identityTol,
            1,
            nFailures
        );
        reportVector
        (
            test,
            "T03_OWNER_AREA_FRAME",
            ownerStored,
            ownerNormal,
            identityTol,
            1,
            nFailures
        );
        reportVector
        (
            test,
            "T03_NEIGHBOUR_AREA_GEOMETRY",
            nbrRaw,
            expectedNeighbourArea,
            identityTol,
            1,
            nFailures
        );
        reportVector
        (
            test,
            "T03_NEIGHBOUR_AREA_FRAME",
            nbrStored,
            expectedNeighbourArea,
            identityTol,
            1,
            nFailures
        );
        report
        (
            test,
            "T03_NEIGHBOUR_MAGNITUDE",
            mag(nbrStored),
            1,
            identityTol,
            nFailures
        );
        reportVector
        (
            test,
            "T03_NEIGHBOUR_CENTRE_GEOMETRY",
            nbrCentre,
            expectedNbrCentre,
            identityTol,
            1,
            nFailures
        );

        // Cell volumes and per-cell closure over the whole local mesh.
        scalar minVolume = std::numeric_limits<scalar>::max();
        forAll(mesh.V(), celli)
        {
            minVolume = min(minVolume, mesh.V()[celli]);
        }
        reduce(minVolume, minOp<scalar>());
        reportBool
        (
            test,
            "T03_CELL_VOLUMES_POSITIVE",
            minVolume > 0,
            true,
            nFailures
        );

        const surfaceVectorField& Sf = mesh.Sf();
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();
        vectorField closure(mesh.nCells(), Zero);
        forAll(nei, facei)
        {
            closure[own[facei]] += Sf[facei];
            closure[nei[facei]] -= Sf[facei];
        }
        forAll(mesh.boundary(), patchi)
        {
            const vectorField& Sfp = Sf.boundaryField()[patchi];
            const labelUList& faceCells = mesh.boundary()[patchi].faceCells();
            forAll(faceCells, i)
            {
                closure[faceCells[i]] += Sfp[i];
            }
        }
        scalar maxClosure = 0;
        forAll(closure, celli)
        {
            maxClosure = max(maxClosure, mag(closure[celli]));
        }
        reduce(maxClosure, maxOp<scalar>());
        report(test, "T03_CELL_CLOSURE", maxClosure, 0, 1e-10, nFailures);

        // Distributed branch is only taken when the AMI is distributed.
        bool distributed = false;
        if (amiPatchA.size())
        {
            distributed = amiPatchA.AMI().distributed();
        }
        reduce(distributed, orOp<bool>());
        reportBool
        (
            test,
            "T03_AMI_DISTRIBUTED",
            distributed,
            UPstream::parRun(),
            nFailures
        );

        return;
    }

    // Internal and donor values are uniform, so faces with a low weight
    // (defaults) and faces with a donor value can be told apart. The scalar
    // field is a control: the rotational transform is the identity, so it
    // cannot expose the frame defect.
    checkLowWeightField<vector>
    (
        test, "U", U, patchA, patchB, amiPatchA, amiPatchB,
        vector(1, 2, 3), nFailures
    );

    volScalarField S
    (
        IOobject("S", runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );
    S.correctBoundaryConditions();
    checkLowWeightField<scalar>
    (
        test, "S", S, patchA, patchB, amiPatchA, amiPatchB,
        scalar(1), nFailures
    );

    volTensorField tField
    (
        IOobject("T", runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );
    tField.correctBoundaryConditions();
    checkLowWeightField<tensor>
    (
        test, "T", tField, patchA, patchB, amiPatchA, amiPatchB,
        tensor(1, 2, 3, 4, 5, 6, 7, 8, 9), nFailures
    );

    volSymmTensorField stField
    (
        IOobject("ST", runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );
    stField.correctBoundaryConditions();
    checkLowWeightField<symmTensor>
    (
        test, "ST", stField, patchA, patchB, amiPatchA, amiPatchB,
        symmTensor(1, 2, 3, 4, 5, 6), nFailures
    );
}

} // End namespace


int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    argList::addOption("test", "Tnn", "Test identifier");
    argList::addOption("gap", "distance");
    argList::addOption("scale", "length");
    argList::addOption("offset", "coordinate");
    argList::addOption("src-rank", "rank");
    argList::addOption("tgt-rank", "rank");
    argList::addOption("maxDistance2", "squared-distance");
    argList::addOption("benchmark-iterations", "count");
    argList::addOption
    (
        "t05-order",
        "order",
        "T05 face order: as-built | src-reversed | tgt-reversed | both-reversed"
    );
    argList::addBoolOption("timing");
    argList::addBoolOption("no-restart");
    argList::addBoolOption("swap");
    argList::addBoolOption("unmatched-source");
    argList::addBoolOption("t05-multi");
    argList::addBoolOption("topology-test");
    argList::addBoolOption("state-cycle");
    argList::addBoolOption("s04-fractions");

    argList args(argc, argv);
    Time runTime(Time::controlDictName, args);
    Info().precision(16);

    const word test = args.getOrDefault<word>("test", "T01");
    const label benchmarkIterations =
        args.getOrDefault<label>("benchmark-iterations", 1);
    if
    (
        args.found("benchmark-iterations")
     && (!args.found("timing") || benchmarkIterations < 1 || test != "T01")
    )
    {
        FatalErrorInFunction
            << "-benchmark-iterations requires -timing and -test T01, "
            << "and must be positive"
            << exit(FatalError);
    }
    label nFailures = 0;

    if (test == "T01" || test == "T05" || test == "T06" || test == "T08"
     || test == "T07" || test == "T10" || test == "T11"
     || test == "T12" || test == "T13")
    {
        runGeometryTest(test, args, runTime, nFailures);
    }
    else if (test == "T02")
    {
        if (args.found("s04-fractions"))
        {
            runCacheFractionTest(runTime, args, nFailures);
        }
        else
        {
            runCacheTest(runTime, nFailures, args.found("timing"));
        }
    }
    else if (test == "T09")
    {
        if (!args.found("state-cycle"))
        {
            FatalErrorInFunction
                << "T09 requires -state-cycle"
                << exit(FatalError);
        }
        runStateTransitionTest(nFailures);
    }
    else if (test == "T03" || test == "T04")
    {
        runMeshTest(test, runTime, nFailures, args.found("timing"));
    }
    else
    {
        FatalErrorInFunction
            << "Unknown test identifier " << test
            << exit(FatalError);
    }

    const label totalFailures =
        returnReduce(nFailures, maxOp<label>());

    if (Pstream::master())
    {
        Info<< "AMI_SUMMARY test=" << test
            << " status=" << (totalFailures ? "fail" : "pass")
            << " failures=" << totalFailures
            << endl;
    }

    return totalFailures ? 1 : 0;
}


// ************************************************************************* //
