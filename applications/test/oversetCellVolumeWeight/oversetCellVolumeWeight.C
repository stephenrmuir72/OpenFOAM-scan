/*---------------------------------------------------------------------------*\
    OpenFOAM: The Open Source CFD Toolbox
    This file is distributed under the GNU General Public License, version 3
    or later, WITHOUT ANY WARRANTY. See <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "columnFvMesh.H"
#include "cellVolumeWeightCellCellStencil.H"
#include "cellVolumeWeightMethod.H"
#include "meshToMesh.H"
#include "tetOverlapVolume.H"
#include "clockTime.H"

#include <cmath>
#include <limits>
#include <algorithm>

using namespace Foam;

void require(const bool condition, const char* description)
{
    if (!returnReduceAnd(condition))
    {
        FatalErrorInFunction << description << exit(FatalError);
    }
}

class volumeHarness : public cellCellStencils::cellVolumeWeight
{
public:
    volumeHarness(const fvMesh& mesh, const dictionary& dict)
    :
        cellVolumeWeight(mesh, dict, false)
    {}

    using cellVolumeWeight::normaliseVolumeWeights;
    using cellVolumeWeight::affineWeights;
    using cellVolumeWeight::decodeRestartType;
    using cellVolumeWeight::weightedDonorVolume;
    using cellVolumeWeight::distanceWeights;

    void seedHistoricalDonors(const bool preferred = false)
    {
        forAll(cellTypes_, celli)
        {
            const point& p = mesh_.cellCentres()[celli];
            const bool background = zoneID()[celli] == 0
                && p.x() > 1 && p.x() < 3
                && p.y() > 1 && p.y() < 3
                && p.z() > 1 && p.z() < 3;
            cellTypes_[celli] = (preferred ? zoneID()[celli] == 2 : background)
                ? HOLE : CALCULATED;
        }
    }

    void checkFallbackExclusions()
    {
        const globalIndex gi(mesh_.nCells());
        labelList ids(4, labelMax);
        forAll(zoneID(), celli)
        {
            if (zoneID()[celli] == 0)
            {
                ids[label(mesh_.cellCentres()[celli].x())] = gi.toGlobal(celli);
            }
        }
        Pstream::listReduce(ids, minOp<label>());
        for (const bool cycle : {false, true})
        {
            labelList types(mesh_.nCells(), CALCULATED);
            labelListList stencil(mesh_.nCells());
            scalarListList weights(mesh_.nCells());
            labelList zones(mesh_.nCells(), -1);
            boolList fallback(mesh_.nCells(), false);
            boolList excluded(mesh_.nCells(), false);
            const boolList historical(mesh_.nCells(), false);
            forAll(types, celli)
            {
                if (zoneID()[celli] != 0) continue;
                const label i = label(mesh_.cellCentres()[celli].x());
                if (i >= 2) continue;
                types[celli] = INTERPOLATED;
                excluded[celli] = true;
                fallback[celli] = true;
                stencil[celli] = labelList(1, cycle ? ids[1-i] : ids[2]);
                weights[celli] = scalarList(1, 1);
            }
            finaliseDonors
            (
                gi, historical, fallback, List<DynamicList<donorCandidate>>(),
                types, stencil, weights, zones, &excluded
            );
            bool valid = true;
            forAll(types, celli)
            {
                if (zoneID()[celli] == 0 && mesh_.cellCentres()[celli].x() < 2)
                {
                    valid &= types[celli] == (cycle ? POROUS : INTERPOLATED);
                    valid &= cycle ? stencil[celli].empty()
                        : stencil[celli] == labelList(1, ids[2]);
                }
            }
            require(valid, "Ordinary donor exclusions changed history or admitted cycles");
        }
    }

    void checkGeometryCoverage()
    {
        const globalIndex gi(mesh_.nCells());
        for (const scalar fraction :
             {scalar(.5), scalar(.75), scalar(.75-1e-8),
              scalar(.75+1e-8), scalar(1), scalar(2)})
        {
            labelList types(mesh_.nCells(), INTERPOLATED);
            labelListList stencil(mesh_.nCells());
            scalarListList weights(mesh_.nCells());
            labelList donorZones(mesh_.nCells(), -1);
            List<DynamicList<donorCandidate>> candidates(mesh_.nCells());
            forAll(types, celli)
            {
                const scalar self = tetOverlapVolume().cellCellOverlapVolumeMinDecomp
                (
                    mesh_, celli, mesh_, celli, treeBoundBox(mesh_.cellBb(celli))
                );
                donorCandidate candidate;
                candidate.zone = 1;
                candidate.donors = labelList(1, gi.toGlobal(celli));
                candidate.volumes = scalarList(1, fraction*self);
                candidates[celli].append(candidate);
            }
            const boolList flags(mesh_.nCells(), false);
            finaliseDonors
            (
                gi, flags, flags, candidates, types, stencil, weights, donorZones
            );
            bool valid = true;
            const bool accepted = fraction >= .75 && fraction <= 1;
            forAll(types, celli)
            {
                valid &= types[celli] == (accepted ? INTERPOLATED : POROUS);
                valid &= accepted ? weights[celli] == scalarList(1, 1)
                    : weights[celli].empty();
            }
            require(valid, "Geometric coverage confused with FV volume/weight sum");
        }
    }

    void checkPriority()
    {
        const labelList types({OTHER, PATCH, OVERSET});
        for (const labelList& order :
             {labelList({0, 1, 2}), labelList({2, 1, 0}),
              labelList({1, 2, 0}), labelList({0, 2, 1}),
              labelList({1, 0, 2}), labelList({2, 0, 1}),
              labelList({1, 1, 2}), labelList({1, 0}),
              labelList({0, 1}), labelList({0}), labelList()})
        {
            labelList result(1, -1);
            interpolatePatchTypes(labelListList(1, order), types, result);
            const label expected = order.contains(2) ? OVERSET
                : order.contains(1) ? PATCH : order.empty() ? -1 : OTHER;
            require(result[0] == expected, "Patch priority depends on donor order");
        }
    }

    void checkReadiness()
    {
        const globalIndex gi(mesh_.nCells());
        labelList ids(4, labelMax);
        forAll(zoneID(), celli)
        {
            if (zoneID()[celli] == 0)
            {
                ids[label(mesh_.cellCentres()[celli].x())] = gi.toGlobal(celli);
            }
        }
        Pstream::listReduce(ids, minOp<label>());
        for (const bool failure : {false, true})
        {
            labelList types(mesh_.nCells(), CALCULATED);
            labelListList stencil(mesh_.nCells());
            scalarListList weights(mesh_.nCells());
            labelList donorZones(mesh_.nCells(), -1);
            boolList oldHoles(mesh_.nCells(), false);
            const boolList fallback(mesh_.nCells(), false);
            forAll(types, celli)
            {
                if (zoneID()[celli] != 0) continue;
                const label i = label(mesh_.cellCentres()[celli].x());
                if (i == 3) oldHoles[celli] = true;
                if (failure && i == 0) oldHoles[celli] = true;
                if (i < 2)
                {
                    types[celli] = INTERPOLATED;
                    stencil[celli] = i == 0 ? labelList({ids[1], ids[2], ids[3]})
                        : labelList(1, ids[2]);
                    weights[celli] = scalarList(stencil[celli].size(), 1);
                }
                if (failure && i >= 2) types[celli] = HOLE;
                if (failure && i == 3)
                {
                    stencil[celli] = labelList(1, ids[2]);
                    weights[celli] = scalarList(1, 1);
                }
            }
            finaliseDonors
            (
                gi, oldHoles, fallback, List<DynamicList<donorCandidate>>(),
                types, stencil, weights, donorZones
            );
            bool valid = true;
            forAll(types, celli)
            {
                if (zoneID()[celli] != 0) continue;
                const label i = label(mesh_.cellCentres()[celli].x());
                if (failure)
                {
                    valid &= types[celli] == (i == 1 ? POROUS : HOLE)
                        && stencil[celli].empty() && weights[celli].empty();
                }
                else if (i == 0)
                {
                    valid &= types[celli] == INTERPOLATED
                        && stencil[celli] == (allowInterpolatedDonors_
                            ? labelList({ids[1], ids[2]}) : labelList(1, ids[2]))
                        && mag(sum(weights[celli])-1) < 1e-14;
                }
            }
            require(valid, "Readiness or protected-hole policy changed");
        }
    }

    void checkDonorChains()
    {
        const globalIndex gi(mesh_.nCells());
        labelList ids(4, labelMax);
        forAll(zoneID(), celli)
        {
            if (zoneID()[celli] == 0)
            {
                ids[label(mesh_.cellCentres()[celli].x())] = gi.toGlobal(celli);
            }
        }
        Pstream::listReduce(ids, minOp<label>());
        labelList types(mesh_.nCells(), CALCULATED);
        labelListList stencil(mesh_.nCells());
        scalarListList weights(mesh_.nCells());
        labelList donorZones(mesh_.nCells(), -1);
        forAll(types, celli)
        {
            if (zoneID()[celli] == 0)
            {
                const label i = label(mesh_.cellCentres()[celli].x());
                types[celli] = i == 3 ? HOLE : INTERPOLATED;
                if (i < 3)
                {
                    stencil[celli] = labelList(1, ids[i+1]);
                    weights[celli] = scalarList(1, 1e-30);
                }
            }
        }
        const boolList flags(mesh_.nCells(), false);
        finaliseDonors
        (
            gi, flags, flags, List<DynamicList<donorCandidate>>(),
            types, stencil, weights, donorZones
        );
        bool valid = true;
        forAll(types, celli)
        {
            if (zoneID()[celli] == 0)
            {
                const label i = label(mesh_.cellCentres()[celli].x());
                valid &= types[celli] == (i == 3 ? HOLE : POROUS)
                    && stencil[celli].empty() && weights[celli].empty();
            }
        }
        require(valid, "Final stencils retained donors demoted later");
    }

    void checkCoverage()
    {
        const globalIndex gi(mesh_.nCells());
        labelList ids(4, labelMax);
        forAll(zoneID(), celli)
        {
            if (zoneID()[celli] == 0)
            {
                ids[label(mesh_.cellCentres()[celli].x())] = gi.toGlobal(celli);
            }
        }
        Pstream::listReduce(ids, minOp<label>());
        for (label scenario = 0; scenario < 10; ++scenario)
        {
            labelList types(mesh_.nCells(), CALCULATED);
            labelListList stencil(mesh_.nCells());
            scalarListList weights(mesh_.nCells());
            labelList donorZones(mesh_.nCells(), -1);
            List<DynamicList<donorCandidate>> candidates(mesh_.nCells());
            forAll(types, celli)
            {
                if (zoneID()[celli] != 0) continue;
                const label i = label(mesh_.cellCentres()[celli].x());
                if (scenario == 1 && i == 1) types[celli] = HOLE;
                if (i != 0) continue;
                types[celli] = INTERPOLATED;
                donorCandidate preferred;
                preferred.zone = 1;
                preferred.donors = labelList({ids[1], ids[2]});
                preferred.volumes = scenario
                    ? scalarList({0.95*mesh_.V()[celli], 0.05*mesh_.V()[celli]})
                    : scalarList({0.005*mesh_.V()[celli], 0.005*mesh_.V()[celli]});
                donorCandidate alternative;
                alternative.zone = 2;
                alternative.donors = labelList(1, ids[3]);
                alternative.volumes = scalarList(1, mesh_.V()[celli]);
                stencil[celli] = preferred.donors;
                weights[celli] = preferred.volumes;
                candidates[celli].append(preferred);
                if (scenario < 2)
                {
                    candidates[celli].append(alternative);
                }
                else
                {
                    const scalar volume = scenario == 2 ? .9
                        : scenario == 3 ? .9-1e-8 : scenario == 4 ? .9+1e-8
                        : scenario == 5 ? -1 : scenario == 6 ? 0
                        : scenario == 7 ? std::numeric_limits<scalar>::infinity()
                        : scenario == 8 ? std::numeric_limits<scalar>::quiet_NaN()
                        : 2;
                    candidates[celli][0].donors = labelList(1, ids[2]);
                    candidates[celli][0].volumes = scalarList(1, volume*mesh_.V()[celli]);
                }
            }
            const boolList flags(mesh_.nCells(), false);
            finaliseDonors
            (
                gi, flags, flags, candidates, types, stencil, weights, donorZones
            );
            bool valid = true;
            forAll(types, celli)
            {
                if (zoneID()[celli] == 0 && mesh_.cellCentres()[celli].x() < 1)
                {
                    if (scenario < 2)
                    {
                        valid &= types[celli] == INTERPOLATED
                            && donorZones[celli] == 2
                            && stencil[celli] == labelList(1, ids[3]);
                    }
                    else if (scenario == 2 || scenario == 4)
                    {
                        valid &= types[celli] == INTERPOLATED
                            && stencil[celli] == labelList(1, ids[2]);
                    }
                    else
                    {
                        valid &= types[celli] == POROUS
                            && stencil[celli].empty() && weights[celli].empty();
                    }
                }
            }
            require(valid, "Coverage/retained-volume policy did not retry zones");
        }
    }
};

autoPtr<polyMesh> boxesMesh
(
    const word& name,
    const Time& time,
    const List<boundBox>& boxes,
    const bool partition,
    const bool skewed = false,
    const scalar warp = 0,
    const label faceStart = 0,
    const bool collinearTop = false
)
{
    DynamicList<point> points;
    DynamicList<face> faces;
    DynamicList<label> owners;
    const labelListList boxFaces
    ({
        labelList({0,4,7,3}), labelList({1,2,6,5}),
        labelList({0,1,5,4}), labelList({3,7,6,2}),
        labelList({0,3,2,1}), labelList({4,5,6,7})
    });
    label celli = 0;
    forAll(boxes, boxi)
    {
        if (partition && boxi % Pstream::nProcs() != Pstream::myProcNo())
        {
            continue;
        }
        const point& a = boxes[boxi].min();
        const point& b = boxes[boxi].max();
        const label base = points.size();
        for (const point& p : pointList
             ({point(a.x(),a.y(),a.z()), point(b.x(),a.y(),a.z()),
               point(b.x(),b.y(),a.z()), point(a.x(),b.y(),a.z()),
               point(a.x(),a.y(),b.z()), point(b.x(),a.y(),b.z()),
               point(b.x(),b.y(),b.z()), point(a.x(),b.y(),b.z())}))
        {
            points.append(p);
        }
        points[base+6].z() += warp*(b.z()-a.z());
        if (collinearTop)
        {
            points.append((points[base+4]+points[base+5])/2);
        }
        for (const labelList& facePoints : boxFaces)
        {
            const labelList vertices = collinearTop && facePoints[0] == 4
                ? labelList({4,8,5,6,7}) : facePoints;
            face f(vertices.size());
            forAll(f, i)
            {
                f[i] = base + vertices[(i+faceStart)%vertices.size()];
            }
            faces.append(f);
            owners.append(celli);
        }
        ++celli;
    }
    if (skewed)
    {
        for (point& p : points)
        {
            p.x() += .3*p.y();
            p.y() += .2*p.z();
        }
    }
    autoPtr<polyMesh> mesh
    (
        new polyMesh
        (
            IOobject
            (
                name, time.timeName(), time,
                IOobject::NO_READ, IOobject::NO_WRITE, false
            ),
            pointField(points), faceList(faces), labelList(owners),
            labelList(), false
        )
    );
    polyPatchList patches(1);
    patches.set
    (
        0, new polyPatch
        (
            "boxBoundary", mesh->nFaces(), 0, 0,
            mesh->boundaryMesh(), word::null
        )
    );
    mesh->addPatches(patches);
    return mesh;
}

scalar exactVolume(const boundBox& a, const boundBox& b)
{
    scalar volume = 1;
    for (direction c = 0; c < 3; ++c)
    {
        volume *= max(scalar(0), min(a.max()[c], b.max()[c])
                                 - max(a.min()[c], b.min()[c]));
    }
    return volume;
}

label boxIndex(const List<boundBox>& boxes, const point& centre)
{
    forAll(boxes, i)
    {
        if (mag(boxes[i].centre()-centre) < 1e-8*boxes[i].mag()) return i;
    }
    FatalErrorInFunction << "Unmatched reference box centre " << centre
        << exit(FatalError);
    return -1;
}

void checkMapper(const Time& time, const scalar scale, const bool corrected,
                 const meshToMesh::procMapMethod distribution, const bool normalise,
                 const bool skewed = false, const bool reversed = false)
{
    List<boundBox> source
    ({
        boundBox(point(0,0,0), point(1,1,1)),
        boundBox(point(2,0,0), point(3,1,1))
    });
    List<boundBox> target
    ({
        boundBox(point(.2,.2,.2), point(.8,.8,.8)),
        boundBox(point(.1,.1,.1), point(.1001,.1001,.1001)),
        boundBox(point(.9,.3,.2), point(2.1,.4,.3)),
        boundBox(point(.4,.2,.2), point(.6,.2000002,.7)),
        boundBox(point(0,1,0), point(1,2,1))
    });
    for (auto& box : source) { box.min() *= scale; box.max() *= scale; }
    for (auto& box : target) { box.min() *= scale; box.max() *= scale; }
    if (reversed) source.swap(target);
    autoPtr<polyMesh> src = boxesMesh("sourceBoxes", time, source, true, skewed);
    autoPtr<polyMesh> tgt = boxesMesh("targetBoxes", time, target, true, skewed);
    auto referencePoint = [&](point p)
    {
        if (skewed)
        {
            p.y() -= .2*p.z();
            p.x() -= .3*p.y();
        }
        return p;
    };
    auto aspectRatio = [](const boundBox& box)
    {
        const vector span = box.span();
        return max(span.x(), max(span.y(), span.z()))
            /min(span.x(), min(span.y(), span.z()));
    };
    meshToMesh mapper
    (
        src(), tgt(), corrected
            ? meshToMesh::interpolationMethod::imCorrectedCellVolumeWeight
            : meshToMesh::interpolationMethod::imCellVolumeWeight,
        HashTable<word>(0), wordList(0), distribution, normalise
    );
    for (label direction = 0; direction < 2; ++direction)
    {
        const polyMesh& receiver = direction ? tgt() : src();
        const polyMesh& donor = direction ? src() : tgt();
        const auto& receiverBoxes = direction ? target : source;
        const auto& donorBoxes = direction ? source : target;
        const auto& addr = direction
            ? mapper.tgtToSrcCellAddr() : mapper.srcToTgtCellAddr();
        const auto& weights = direction
            ? mapper.tgtToSrcCellWght() : mapper.srcToTgtCellWght();
        const auto& offsets = direction
            ? mapper.tgtToSrcCellVec() : mapper.srcToTgtCellVec();
        pointField donorCentres(donor.cellCentres());
        const auto& map = direction ? mapper.srcMap() : mapper.tgtMap();
        if (map) map->distribute(donorCentres);
        const vector gradient = vector(.7,-.2,.9)/scale;
        scalarField linearValues(donor.nCells());
        vectorField vectorValues(donor.nCells());
        forAll(linearValues, celli)
        {
            linearValues[celli] = 1.25 + (gradient & donor.cellCentres()[celli]);
            vectorValues[celli] = donor.cellCentres()[celli]/scale + vector(.1,.2,.3);
        }
        const vectorField scalarGradients(donor.nCells(), gradient);
        const tensorField vectorGradients(donor.nCells(), tensor::I/scale);
        scalarField mappedScalar(receiver.nCells(), Zero);
        vectorField mappedVector(receiver.nCells(), Zero);
        if (corrected)
        {
            if (direction)
            {
                mapper.mapSrcToTgt<scalar>
                (
                    linearValues, scalarGradients, plusEqOp<scalar>(), mappedScalar
                );
                mapper.mapSrcToTgt<vector>
                (
                    vectorValues, vectorGradients, plusEqOp<vector>(), mappedVector
                );
            }
            else
            {
                mapper.mapTgtToSrc<scalar>
                (
                    linearValues, scalarGradients, plusEqOp<scalar>(), mappedScalar
                );
                mapper.mapTgtToSrc<vector>
                (
                    vectorValues, vectorGradients, plusEqOp<vector>(), mappedVector
                );
            }
        }
        bool valid = addr.size() == receiver.nCells();
        scalar maxVolumeError = 0;
        forAll(addr, celli)
        {
            const label a = boxIndex
            (
                receiverBoxes, referencePoint(receiver.cellCentres()[celli])
            );
            const boundBox& acceptor = receiverBoxes[a];
            label nExpected = 0;
            scalar total = 0;
            for (const auto& box : donorBoxes)
            {
                const scalar volume = exactVolume(acceptor, box);
                if (volume > 1e-6*min(acceptor.volume(), box.volume()))
                {
                    ++nExpected;
                    total += volume;
                }
            }
            valid &= addr[celli].size() == nExpected
                && weights[celli].size() == nExpected;
            if (corrected) valid &= offsets[celli].size() == nExpected;
            scalar expectedScalar = 0;
            vector expectedVector(Zero);
            scalar cellTolerance = 1e-10;
            labelHashSet seen;
            forAll(addr[celli], j)
            {
                const label slot = addr[celli][j];
                const label b = boxIndex(donorBoxes, referencePoint(donorCentres[slot]));
                valid &= seen.insert(b);
                const scalar raw = exactVolume(acceptor, donorBoxes[b]);
                const scalar expected = normalise ? raw/total : raw;
                // Thin-cell clipping loses digits in coordinate subtraction.
                const scalar tolerance = max
                (
                    scalar(1e-10),
                    8*std::numeric_limits<scalar>::epsilon()
                        *max(aspectRatio(acceptor), aspectRatio(donorBoxes[b]))
                );
                cellTolerance = max(cellTolerance, tolerance);
                const scalar volumeError = mag(weights[celli][j]-expected)/expected;
                maxVolumeError = max(maxVolumeError, volumeError);
                if (volumeError > tolerance)
                {
                    Pout<< "CVW overlap error: receiver " << a << " donor " << b
                        << " value " << weights[celli][j] << " reference " << expected
                        << " relative error " << mag(weights[celli][j]-expected)/expected
                        << endl;
                }
                valid &= mag(weights[celli][j]-expected)
                    <= tolerance*max(expected, VSMALL);
                if (corrected)
                {
                    point centroid(Zero);
                    for (Foam::direction c = 0; c < 3; ++c)
                    {
                        centroid[c] =
                            (max(acceptor.min()[c], donorBoxes[b].min()[c])
                             + min(acceptor.max()[c], donorBoxes[b].max()[c]))/2;
                    }
                    if (skewed)
                    {
                        centroid.x() += .3*centroid.y();
                        centroid.y() += .2*centroid.z();
                    }
                    valid &= mag(donorCentres[slot] + offsets[celli][j]-centroid)
                        < tolerance*scale;
                    expectedScalar += expected*(1.25 + (gradient & centroid));
                    expectedVector += expected*(centroid/scale + vector(.1,.2,.3));
                }
            }
            if (corrected)
            {
                valid &= mag(mappedScalar[celli]-expectedScalar)
                        <= cellTolerance*max(mag(expectedScalar), VSMALL)
                    && mag(mappedVector[celli]-expectedVector)
                        <= cellTolerance*max(mag(expectedVector), VSMALL);
            }
        }
        reduce(maxVolumeError, maxOp<scalar>());
        if (maxVolumeError > 1e-10 && Pstream::master())
        {
            Info<< "CVW thin-cell relative volume error: " << maxVolumeError << endl;
        }
        require(valid, "Mapper differs from analytic box intersections/centroids");
    }
}

void checkRefinement(const Time& time, const fvMesh& mesh)
{
    dictionary dict;
    dict.add("weightingScheme", word("linearExact"));
    dict.add("allowNegativeWeights", true);
    volumeHarness affine(mesh, dict);
    scalar lastVolumeError = 0, lastAffineError = 0;
    for (const label n : {8, 16, 32})
    {
        const scalar h = 1.0/n;
        const label donor = n/2;
        List<boundBox> target(n);
        forAll(target, i)
        {
            target[i] = boundBox(point(i*h,0,0), point((i+1)*h,1,1));
        }
        const List<boundBox> source
        ({
            boundBox(point((donor+.35)*h,0,0), point((donor+.85)*h,1,1))
        });
        autoPtr<polyMesh> src = boxesMesh("refineSource", time, source, false);
        autoPtr<polyMesh> tgt = boxesMesh("refineTarget", time, target, false);
        cellVolumeWeightMethod method(src(), tgt());
        labelListList a, b;
        scalarListList wa, wb;
        pointListList va, vb;
        method.calculate(a, wa, va, b, wb, vb);
        require(a[0].size() == 1 && a[0][0] == donor,
            "Refinement reference overlaps incorrect");
        const point sample = src->cellCentres()[0];
        const scalar x = sample.x();
        const scalar xc = tgt->cellCentres()[donor].x();
        const scalar volumeError = mag(xc*xc-x*x);
        pointList centres
        ({
            tgt->cellCentres()[donor-1], tgt->cellCentres()[donor],
            tgt->cellCentres()[donor+1]
        });
        scalarList weights({0,1,0});
        require(affine.affineWeights(sample, centres, weights),
            "Refinement affine correction rejected");
        scalar reconstructed = 0;
        forAll(weights, i) reconstructed += weights[i]*sqr(centres[i].x());
        const scalar affineError = mag(reconstructed-x*x);
        if (lastVolumeError)
        {
            require(lastVolumeError/volumeError > 1.8
                && lastAffineError/affineError > 3.8,
                "Manufactured interpolation errors did not converge");
        }
        lastVolumeError = volumeError;
        lastAffineError = affineError;
        if (Pstream::master())
        {
            Info<< "CVW refinement: " << n << " donor cells, volume error "
                << volumeError << ", affine error " << affineError << endl;
        }
    }
}

void checkSliverThreshold(const Time& time)
{
    const List<boundBox> source
    ({
        boundBox(point(0,0,0), point(1,1,1))
    });
    const List<boundBox> target
    ({
        boundBox(point(1-.5e-6,.2,0), point(2-.5e-6,.3,1)),
        boundBox(point(1-2e-6,.6,0), point(2-2e-6,.7,1))
    });
    autoPtr<polyMesh> src = boxesMesh("sliverSource", time, source, false);
    autoPtr<polyMesh> tgt = boxesMesh("sliverTarget", time, target, false);
    cellVolumeWeightMethod method(src(), tgt());
    labelListList a, b;
    scalarListList wa, wb;
    pointListList va, vb;
    method.calculate(a, wa, va, b, wb, vb);
    const scalar accepted = exactVolume(source[0], target[1]);
    const scalar discarded = exactVolume(source[0], target[0]);
    const scalar roundoff =
        64*std::numeric_limits<scalar>::epsilon()*target[0].volume();
    if (Pstream::master())
    {
        Info<< "CVW slivers: addressing " << a[0] << ", accepted "
            << method.V() << " (reference " << accepted << "), discarded "
            << method.discardedVolume() << " (reference " << discarded
            << "), exact tests " << method.candidateTests() << endl;
    }
    require(a[0] == labelList(1, 1) && b[0].empty()
        && b[1] == labelList(1, 0) && method.candidateTests() == 2
        && mag(method.V()-accepted) < 1e-10*accepted + roundoff
        && mag(method.discardedVolume()-discarded) < 1e-10*discarded + roundoff,
        "Sliver threshold or discarded-volume accounting is incorrect");
}

void checkSearchScaling(const Time& time)
{
    for (const label size : {16, 64, 256})
    {
        List<boundBox> source(size), target(size);
        forAll(source, i)
        {
            source[i] = boundBox(point(2*i,0,0), point(2*i+1,1,1));
            target[i] = boundBox(point(2*i+.2,.2,.2), point(2*i+.8,.8,.8));
        }
        autoPtr<polyMesh> src = boxesMesh("scalingSource", time, source, false);
        autoPtr<polyMesh> tgt = boxesMesh("scalingTarget", time, target, false);
        cellVolumeWeightMethod method(src(), tgt());
        labelListList a, b;
        scalarListList wa, wb;
        pointListList va, vb;
        const clockTime timer;
        method.calculate(a, wa, va, b, wb, vb);
        require(method.candidateTests() == size, "Duplicate/nonlinear candidate work");
        Info<< "CVW work scaling: " << size << " cells, "
            << method.candidateTests() << " exact tests, "
            << timer.elapsedTime() << " s" << endl;
    }
}

void checkWarpedGeometry(const Time& time, const bool mapping)
{
    for (const scalar warp : {scalar(1e-9), scalar(1e-6), scalar(.04), scalar(-.04)})
    for (const label start : {label(0), label(1)})
    for (const scalar scale : {scalar(1e-6), scalar(1), scalar(1e3)})
    {
        const point origin = scale*point(10, -20, 30);
        const List<boundBox> boxes
        (
            1, boundBox(origin, origin+scale*point(1,1,1))
        );
        autoPtr<polyMesh> source
        (
            boxesMesh("warpedSource", time, boxes, true, false, warp, start)
        );
        autoPtr<polyMesh> target
        (
            boxesMesh("warpedTarget", time, boxes, true, false, warp, start)
        );
        scalar geometricVolume = 0;
        point geometricCentre(Zero);
        bool valid = true;
        if (source->nCells())
        {
            scalar fanVolume = 0;
            vector moment(Zero);
            const point& reference = source->points()[0];
            for (const label facei : source->cells()[0])
            {
                const face& f = source->faces()[facei];
                const point& a = source->points()[f[0]];
                for (label i = 1; i+1 < f.size(); ++i)
                {
                    const point& b = source->points()[f[i]];
                    const point& c = source->points()[f[i+1]];
                    const scalar volume =
                        ((a-reference) & ((b-reference)^(c-reference)))/6;
                    geometricVolume += volume;
                    moment += volume*((reference+a+b+c)/4);
                    const point& centre = source->cellCentres()[0];
                    const scalar fan =
                        ((a-centre) & ((b-centre)^(c-centre)))/6;
                    valid &= fan > 0 && std::isfinite(fan);
                    fanVolume += fan;
                }
            }
            geometricCentre = moment/geometricVolume;
            const auto self = tetOverlapVolume().cellCellOverlapMomentMinDecomp
            (
                source(), 0, source(), 0, treeBoundBox(source->cellBb(0))
            );
            valid &= mag(fanVolume-geometricVolume) < 1e-10*geometricVolume
                && mag(self.first()-geometricVolume) < 1e-10*geometricVolume
                && mag(self.second()/self.first()-geometricCentre) < 1e-10*scale;
            if (warp == .04 && start == 0 && scale == 1)
            {
                Info<< "CVW warped reference: FV=" << source->cellVolumes()[0]
                    << " fan=" << fanVolume << " boundary=" << geometricVolume
                    << " self=" << self.first() << " excess="
                    << self.first()/source->cellVolumes()[0]-1 << endl;
            }
        }
        require(valid, "Warped fixture has invalid tessellation/volume/moment");
        if (mapping)
        {
            for (const auto method :
                 {meshToMesh::interpolationMethod::imCellVolumeWeight,
                  meshToMesh::interpolationMethod::imCorrectedCellVolumeWeight})
            {
                meshToMesh mapper
                (
                    source(), target(), method, HashTable<word>(), wordList(),
                    meshToMesh::procMapMethod::pmAABB, false
                );
                if (source->nCells())
                {
                    valid &= mapper.srcToTgtCellAddr()[0].size() == 1
                        && mapper.tgtToSrcCellAddr()[0].size() == 1
                        && mag(mapper.srcToTgtCellWght()[0][0]-geometricVolume)
                            < 1e-10*geometricVolume
                        && mapper.srcToTgtCellWght()[0][0]
                            == mapper.tgtToSrcCellWght()[0][0];
                    if (method == meshToMesh::interpolationMethod::imCorrectedCellVolumeWeight)
                    {
                        valid &= mag(mapper.srcToTgtCellVec()[0][0]
                            -(geometricCentre-target->cellCentres()[0])) < 1e-10*scale;
                    }
                }
                require(valid, "Warped mapper changed volume/moment symmetry");
            }
            for (const bool partial : {false, true})
            for (const bool reversed : {false, true})
            {
                const point lo = partial ? point(.8,.2,.2) : point(.2,.2,.2);
                const point hi = partial ? point(1.2,.8,.8) : point(.8,.8,.8);
                autoPtr<polyMesh> piece
                (
                    boxesMesh
                    (
                        "warpedPiece", time,
                        List<boundBox>(1, boundBox(origin+scale*lo, origin+scale*hi)),
                        true
                    )
                );
                const polyMesh& from = reversed ? piece() : source();
                const polyMesh& to = reversed ? source() : piece();
                meshToMesh mapper
                (
                    from, to,
                    meshToMesh::interpolationMethod::imCorrectedCellVolumeWeight,
                    HashTable<word>(), wordList(),
                    meshToMesh::procMapMethod::pmAABB, false
                );
                const scalar expectedVolume = (partial ? .072 : .216)*pow3(scale);
                const point centre = origin+scale
                    *(partial ? point(.9,.5,.5) : point(.5,.5,.5));
                const vector gradient = vector(.7,-.2,.9)/scale;
                scalarField values(from.nCells());
                vectorField vectors(from.nCells());
                forAll(values, celli)
                {
                    values[celli] = 1.25+(gradient & from.cellCentres()[celli]);
                    vectors[celli] = from.cellCentres()[celli]/scale;
                }
                scalarField mapped(to.nCells(), Zero);
                vectorField mappedVectors(to.nCells(), Zero);
                mapper.mapSrcToTgt<scalar>
                (
                    values, vectorField(from.nCells(), gradient),
                    plusEqOp<scalar>(), mapped
                );
                mapper.mapSrcToTgt<vector>
                (
                    vectors, tensorField(from.nCells(), tensor::I/scale),
                    plusEqOp<vector>(), mappedVectors
                );
                if (to.nCells())
                {
                    const scalar expected = expectedVolume*(1.25+(gradient & centre));
                    valid &= mag(mapper.tgtToSrcCellWght()[0][0]-expectedVolume)
                        < 1e-10*expectedVolume
                        && mag(mapped[0]-expected) < 1e-10*mag(expected)
                        && mag(mappedVectors[0]-expectedVolume*centre/scale)
                            < 1e-10*expectedVolume*mag(centre/scale);
                }
                require(valid, "Warped partial/reversed scalar/vector transfer failed");
            }
        }
    }
    Info<< "Passed CVW warped geometry" << endl;
}

int main(int argc, char *argv[])
{
    argList::addOption("invalidSetting", "name", "Check invalid configuration rejection");
    argList::addOption("restartValue", "value", "Verify decoded restart classification");
    argList::addBoolOption("linearExact", "Exercise full affine update");
    argList::addBoolOption("coverage", "Exercise full coverage update");
    argList::addBoolOption("fixedLayer", "Exercise fixed layer selection");
    argList::addBoolOption("warpedProbe", "Compare warped geometry representations");
    argList::addBoolOption("warpedMapper", "Exercise warped plain/corrected mapping");
    argList::addBoolOption("warpedUpdate", "Exercise warped quality and mesh motion");
    argList::addBoolOption("invalidGeometry", "Reject an inverted first-vertex fan");
    argList::addBoolOption("nonFiniteRestart", "Reject non-finite restart decoding");
    argList::addOption("donorFallback", "mode", "Exercise optional donor fallback");
    argList::addBoolOption("fallbackExercise", "Exercise ordinary recovery in 3D");
    argList::addBoolOption("boundedAffine", "Exercise bounded affine recovery");
    argList::addBoolOption("alternateExercise", "Exercise preferred-zone failure in 3D");
    argList::addBoolOption("holeLeakProbe", "Diagnose the analytic closed-cube solid fixture");
    argList::addOption("updates", "count", "Number of 3D exercise updates");
    argList::addOption("processorMapMethod", "method", "AABB or LOD distribution");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    Info().precision(17);
    if (args.found("holeLeakProbe"))
    {
        volumeHarness update(mesh, dictionary());
        update.update();
        label inside = 0, active = 0;
        bool planar = true;
        const boundBox solid(point(1.8,1.8,1.8), point(3.2,3.2,3.2));
        forAll(update.cellTypes(), celli)
        {
            const scalar geometric =
                tetOverlapVolume::cellVolumeMinDecomp(mesh, celli).first();
            planar &= mag(geometric/mesh.V()[celli]-1) < 1e-10;
            if (update.zoneID()[celli] == 0 && solid.contains(mesh.cellCentres()[celli]))
            {
                ++inside;
                active += update.cellTypes()[celli] != volumeHarness::HOLE;
            }
        }
        reduce(inside, sumOp<label>());
        reduce(active, sumOp<label>());
        require(planar && inside == 1, "Invalid analytic closed-solid fixture");
        Info<< "CVW closed-solid probe: " << inside << " interior receivers, "
            << active << " active." << endl;
        if (active)
        {
            WarningInFunction
                << "Active receiver inside the explicitly defined closed solid. "
                << "Patch/fringe reachability is not an inside/outside classifier; "
                << "this diagnostic does not change hole classification."
                << nl << endl;
        }
        return 0;
    }
    if (args.found("fallbackExercise") || args.found("alternateExercise"))
    {
        dictionary config;
        const word mode(args.get<word>("donorFallback"));
        config.add("donorFallback", mode);
        config.add("layerDonorMetric", word("overlapWeighted"));
        const bool alternate = args.found("alternateExercise");
        config.add("layerRelax", alternate ? scalar(1) : scalar(.5));
        if (args.found("processorMapMethod"))
        {
            config.add("processorMapMethod", args.get<word>("processorMapMethod"));
        }
        if (args.found("linearExact") || args.found("boundedAffine"))
        {
            config.add("weightingScheme", word("linearExact"));
            config.add("allowNegativeWeights", args.found("linearExact"));
        }
        if (args.found("coverage"))
        {
            config.add("donorSelection", word("coverage"));
            config.add("minCoverage", scalar(.5));
        }
        volumeHarness update(mesh, config);
        const label iterations = args.found("updates") ? args.get<label>("updates") : 2;
        require(iterations > 0, "Positive update count required");
        for (label iteration = 0; iteration < iterations; ++iteration)
        {
            update.seedHistoricalDonors(alternate);
            update.update();
            label recovered = 0, partial = 0;
            bool valid = true;
            scalarList donorZones(mesh.nCells());
            labelList historical(mesh.nCells(), Zero);
            forAll(historical, celli)
            {
                const point& p = mesh.cellCentres()[celli];
                const bool background = update.zoneID()[celli] == 0
                    && p.x() > 1 && p.x() < 3
                    && p.y() > 1 && p.y() < 3
                    && p.z() > 1 && p.z() < 3;
                historical[celli] = alternate ? update.zoneID()[celli] == 2 : background;
                donorZones[celli] = update.zoneID()[celli];
            }
            update.cellInterpolationMap().distribute(historical);
            update.cellInterpolationMap().distribute(donorZones);
            for (const label celli : update.interpolationCells())
            {
                const scalarList& weights = update.cellInterpolationWeights()[celli];
                valid &= !weights.empty() && mag(sum(weights)-1) < 1e-12;
                for (const scalar weight : weights)
                {
                    valid &= std::isfinite(weight)
                        && (args.found("linearExact") || weight >= 0);
                }
                for (const label slot : update.cellStencil()[celli])
                {
                    valid &= !historical[slot];
                    if (update.zoneID()[celli] == 1)
                    {
                        valid &= donorZones[slot] == (alternate ? 0 : 1);
                    }
                }
                if (update.zoneID()[celli] == 1)
                {
                    ++recovered;
                    const scalar blend = update.cellInterpolationWeight()[celli];
                    partial += blend > 0 && blend < 1;
                }
            }
            reduce(recovered, sumOp<label>());
            reduce(partial, sumOp<label>());
            require(valid, "Ordinary recovery retained unsafe donors/invalid weights");
            require(alternate
                ? (mode != "legacy" || args.found("coverage")
                    ? recovered > 0 : recovered == 0)
                : (mode == "sameZone" ? recovered > 0 && partial == recovered
                    : recovered == 0), "Donor fallback policy/blend was not preserved");
            Info<< (alternate ? "CVW alternate fallback: " : "CVW ordinary fallback: ")
                << mode << " recovered "
                << recovered << ", partial " << partial << endl;
        }
        Info<< (alternate ? "Passed CVW alternate fallback"
            : "Passed CVW ordinary fallback") << endl;
        return 0;
    }
    if (args.found("warpedProbe") || args.found("warpedMapper"))
    {
        checkWarpedGeometry(runTime, args.found("warpedMapper"));
        return 0;
    }
    if (args.found("invalidGeometry"))
    {
        autoPtr<polyMesh> invalid
        (
            boxesMesh
            (
                "invalidFan", runTime,
                List<boundBox>(1, boundBox(point(0,0,0), point(1,1,1))),
                true, false, -1.2
            )
        );
        if (invalid->nCells())
        {
            tetOverlapVolume::cellVolumeMinDecomp(invalid(), 0);
        }
        require(false, "Invalid geometry was accepted");
    }
    if (args.found("nonFiniteRestart"))
    {
        bool rounded = false, clamped = false;
        volumeHarness::decodeRestartType
        (
            std::numeric_limits<scalar>::quiet_NaN(), rounded, clamped
        );
        require(false, "Non-finite restart was accepted");
    }
    if (args.found("warpedUpdate"))
    {
        for (const scalar warp : {scalar(.04), scalar(-.08)})
        {
            pointField points(mesh.points());
            for (point& p : points)
            {
                p.z() += warp*p.x()*p.y()*p.z()/4;
            }
            mesh.movePoints(points);
            dictionary config;
            config.add("donorSelection", word("coverage"));
            config.add("minCoverage", scalar(.75));
            volumeHarness warped(mesh, config);
            warped.checkGeometryCoverage();
            volumeHarness update(mesh, dictionary());
            update.update();
            bool valid = true;
            require
            (
                returnReduce(update.interpolationCells().size(), sumOp<label>()) > 0,
                "Warped update exercised no interpolation"
            );
            for (const label celli : update.interpolationCells())
            {
                const scalarList& weights = update.cellInterpolationWeights()[celli];
                valid &= !weights.empty() && mag(sum(weights)-1) < 1e-12;
            }
            require(valid, "Warped update lost partition of unity");
        }
        for (const scalar scale : {scalar(1e-6), scalar(1e9)})
        {
            mesh.movePoints(scale*mesh.points());
            dictionary config;
            config.add("minCoverage", scalar(.75));
            volumeHarness scaled(mesh, config);
            scaled.checkGeometryCoverage();
        }
        Info<< "Passed CVW warped update" << endl;
        return 0;
    }

    if (args.found("invalidSetting"))
    {
        dictionary invalid;
        const word key(args.get<word>("invalidSetting"));
        if (key == "minCoverage" || key == "minRetainedFraction" ||
            key == "overlapTolerance") invalid.add(key, scalar(-1));
        else if (key == "stencilExpansionRings") invalid.add(key, label(2));
        else invalid.add(key, word("invalid"));
        volumeHarness rejected(mesh, invalid);
        FatalErrorInFunction << "Invalid setting was accepted" << exit(FatalError);
    }

    dictionary defaults;
    volumeHarness harness(mesh, defaults);
    if (args.found("restartValue"))
    {
        const scalar value = args.get<scalar>("restartValue");
        const label expected = label(std::round(max
        (
            scalar(volumeHarness::CALCULATED),
            min(scalar(volumeHarness::POROUS), value)
        )));
        bool valid = true;
        for (const label type : harness.cellTypes()) valid &= type == expected;
        require(valid, "Incorrect restart classification");
        Info<< "Passed CVW restart decoding" << endl;
        return 0;
    }
    harness.checkPriority();
    harness.checkDonorChains();
    harness.checkReadiness();
    harness.checkFallbackExclusions();
    labelList order({0,1,2,3});
    const pointList centres
    ({point(1,0,0), point(0,0,0), point(0,1,0), point(0,0,1)});
    const scalarList donorVolumes({2,8,1,100});
    do
    {
        require(mag(harness.weightedDonorVolume
        (
            pointList(UIndirectList<point>(centres, order)),
            scalarList(UIndirectList<scalar>(donorVolumes, order)),
            scalarList(4, .25)
        )-27.75) < 1e-14, "Weighted metric depends on donor permutation");
        const pointList donors(UIndirectList<point>(centres, order));
        scalarList weights;
        harness.distanceWeights(point(.2,.3,.4), donors, weights);
        const scalar score = harness.weightedDonorVolume
        (
            donors, scalarList(UIndirectList<scalar>(donorVolumes, order)), weights
        );
        scalarList referenceWeights;
        harness.distanceWeights(point(.2,.3,.4), centres, referenceWeights);
        require(mag(score-harness.weightedDonorVolume
            (centres, donorVolumes, referenceWeights)) < 1e-13,
            "Weighted IDW metric depends on neighbour permutation");
    } while (std::next_permutation(order.begin(), order.end()));
    dictionary strict;
    strict.add("allowInterpolatedDonors", false);
    volumeHarness strictHarness(mesh, strict);
    strictHarness.checkReadiness();
    dictionary quality;
    quality.add("donorSelection", word("coverage"));
    quality.add("minCoverage", scalar(.9));
    quality.add("minRetainedFraction", scalar(.9));
    volumeHarness qualityHarness(mesh, quality);
    qualityHarness.checkCoverage();
    for (const scalar factor : {scalar(1e-300), scalar(1), scalar(1e300)})
    {
        scalarList weights({factor, 3*factor});
        require(harness.normaliseVolumeWeights(weights)
            && mag(weights[0]-.25) < 1e-14, "Volume normalization is scale dependent");
    }
    for (const scalarList& invalid :
         {scalarList(), scalarList({0,0}), scalarList({-1,2}),
          scalarList({std::numeric_limits<scalar>::infinity()})})
    {
        scalarList weights(invalid);
        require(!harness.normaliseVolumeWeights(weights), "Invalid volume weights accepted");
    }

    for (const bool signedWeights : {false, true})
    {
        dictionary config;
        config.add("weightingScheme", word("linearExact"));
        config.add("allowNegativeWeights", signedWeights);
        volumeHarness affine(mesh, config);
        for (const scalar anisotropy : {scalar(1), scalar(1e-6), scalar(1e6)})
        {
            for (const point& translation :
                 {point(0,0,0), point(10,-20,30), point(1e6,-2e6,3e6)})
            {
                const pointList centres
                ({
                    translation, translation + point(anisotropy,0,0),
                    translation + point(0,1,0), translation + point(0,0,1)
                });
                const point sample = translation + point(.2*anisotropy,.3,.4);
                scalarList weights(4, .25);
                require(affine.affineWeights(sample, centres, weights),
                    "Feasible affine correction rejected");
                vector value(Zero);
                forAll(weights, i) value += weights[i]*centres[i];
                const scalar tolerance = 1e-10*max(anisotropy, scalar(1))
                    + 128*std::numeric_limits<scalar>::epsilon()*mag(sample);
                require(mag(value-sample) < tolerance, "Affine moment is incorrect");
                const scalarList baseline(weights);
                const bool accepted = affine.affineWeights
                (
                    translation + point(2*anisotropy,2,2), centres, weights
                );
                require(accepted == signedWeights
                    && (accepted || weights == baseline), "Bounded/signed policy is incorrect");
            }
        }
        scalarList hitWeights(1, 1);
        require(affine.affineWeights(point(1,2,3), pointList(1, point(1,2,3)), hitWeights),
            "Exact-hit one-donor stencil rejected");
        require(!affine.affineWeights(point(1.1,2,3), pointList(1, point(1,2,3)), hitWeights)
            && hitWeights[0] == 1, "Infeasible one-donor stencil accepted");
        scalarList planarWeights(3, 1.0/3);
        const pointList planarCentres({point(0,0,0), point(1,0,0), point(0,1,0)});
        require(affine.affineWeights(point(.2,.3,0), planarCentres, planarWeights),
            "Feasible planar stencil rejected");
        require(!affine.affineWeights(point(.2,.3,.1), planarCentres, planarWeights),
            "Infeasible planar stencil accepted");
        scalarList weights({.5,.5});
        require(affine.affineWeights(point(0,0,0),
            pointList({point(-1,0,0),point(1,0,0)}), weights),
            "Feasible collinear stencil rejected");
        require(!affine.affineWeights(point(0,.1,0),
            pointList({point(-1,0,0),point(1,0,0)}), weights),
            "Infeasible collinear stencil accepted");
    }
    for (const scalar scale : {scalar(1), scalar(1e-6), scalar(1e3)})
    {
        for (const bool corrected : {false, true})
        {
            for (const bool normalise : {false, true})
            {
                for (const auto distribution :
                     {meshToMesh::procMapMethod::pmAABB, meshToMesh::procMapMethod::pmLOD})
                {
                    checkMapper(runTime, scale, corrected, distribution, normalise);
                }
            }
        }
    }
    checkMapper(runTime, 1, true, meshToMesh::procMapMethod::pmAABB, false, true);
    checkMapper(runTime, 1, false, meshToMesh::procMapMethod::pmLOD, true, true);
    checkMapper(runTime, 1, true, meshToMesh::procMapMethod::pmAABB, false, true, true);
    checkMapper(runTime, 1e-6, false, meshToMesh::procMapMethod::pmLOD, true, false, true);
    checkRefinement(runTime, mesh);
    checkSliverThreshold(runTime);
    checkSearchScaling(runTime);
    {
        const List<boundBox> boxes(1, boundBox(point(0,0,0), point(1,1,1)));
        autoPtr<polyMesh> src
        (
            boxesMesh("collinearFace", runTime, boxes, true, false, 0, 0, true)
        );
        autoPtr<polyMesh> tgt(boxesMesh("plainFace", runTime, boxes, true));
        cellVolumeWeightMethod method(src(), tgt());
        labelListList a, b;
        scalarListList wa, wb;
        pointListList va, vb;
        method.calculate(a, wa, va, b, wb, vb);
        require(!src->nCells() ||
            (a[0].size() == 1 && mag(wa[0][0]-1) < 1e-12),
            "Harmless collinear face triangle was rejected");
    }

    dictionary updateConfig;
    if (args.found("donorFallback"))
    {
        updateConfig.add("donorFallback", args.get<word>("donorFallback"));
    }
    if (args.found("linearExact"))
    {
        updateConfig.add("weightingScheme", word("linearExact"));
        updateConfig.add("allowNegativeWeights", true);
    }
    if (args.found("coverage"))
    {
        updateConfig.add("donorSelection", word("coverage"));
        updateConfig.add("minCoverage", scalar(.01));
        updateConfig.add("minRetainedFraction", scalar(.5));
        updateConfig.add("layerDonorMetric", word("overlapWeighted"));
        updateConfig.add("processorMapMethod", word("LOD"));
    }
    if (args.found("fixedLayer"))
    {
        updateConfig.add("holeLayers", label(2));
        updateConfig.add("useLayer", label(1));
    }
    volumeHarness update(mesh, updateConfig);
    for (label iteration = 0; iteration < 3; ++iteration)
    {
        update.update();
        scalarList field(mesh.nCells(), scalar(2.75));
        update.cellInterpolationMap().distribute(field);
        const vector gradient(.7,-.2,.9);
        scalarList linearField(mesh.nCells());
        forAll(linearField, celli)
        {
            linearField[celli] = 1.25 + (gradient & mesh.cellCentres()[celli]);
        }
        update.cellInterpolationMap().distribute(linearField);
        pointField donorCentres(mesh.cellCentres());
        update.cellInterpolationMap().distribute(donorCentres);
        bool valid = true;
        scalarList signature(12*5, Zero);
        label nInterpolated = update.interpolationCells().size();
        for (const label celli : update.interpolationCells())
        {
            const auto& stencil = update.cellStencil()[celli];
            const auto& weights = update.cellInterpolationWeights()[celli];
            scalar value = 0;
            scalar linearValue = 0;
            vector moment(Zero);
            forAll(weights, j) value += weights[j]*field[stencil[j]];
            forAll(weights, j) linearValue += weights[j]*linearField[stencil[j]];
            forAll(weights, j) moment += weights[j]*donorCentres[stencil[j]];
            valid &= std::isfinite(value) && mag(value-2.75) < 1e-10
                && stencil.size() == weights.size();
            if (args.found("linearExact"))
            {
                valid &= mag(moment-mesh.cellCentres()[celli]) < 1e-10;
                valid &= mag(linearValue
                    - (1.25 + (gradient & mesh.cellCentres()[celli]))) < 1e-10;
            }
            else
            {
                for (const scalar weight : weights) valid &= weight >= 0;
            }
        }
        forAll(update.cellTypes(), celli)
        {
            const label index = harness.zoneID()[celli] == 0
                ? label(mesh.cellCentres()[celli].x())
                : 4 + label((mesh.cellCentres()[celli].x()-.1)/.475);
            signature[5*index] = update.cellTypes()[celli]+1;
            signature[5*index+1] = sum(update.cellInterpolationWeights()[celli]);
            forAll(update.cellStencil()[celli], j)
            {
                const scalar weight = update.cellInterpolationWeights()[celli][j];
                const point& centre = donorCentres[update.cellStencil()[celli][j]];
                for (Foam::direction c = 0; c < 3; ++c)
                {
                    signature[5*index+2+c] += weight*centre[c];
                }
            }
        }
        Pstream::listReduce(signature, sumOp<scalar>());
        if (Pstream::master())
        {
            Info<< "CVW signature:";
            for (const scalar value : signature) Info<< ' ' << value;
            Info<< endl;
        }
        reduce(nInterpolated, sumOp<label>());
        require(valid && nInterpolated > 0,
            "Full CVW update violates constant/affine preservation");
    }
    if (Pstream::master()) Info<< "Passed cellVolumeWeight regressions" << endl;
    return 0;
}
