/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "localPointRegion.H"
#include "syncTools.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "globalIndex.H"
#include "indirectPrimitivePatch.H"
#include "dummyTransform.H"
#include "mergePoints.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(localPointRegion, 0);

// Reduction class to get minimum value over face.
class minEqOpFace
{
public:

    void operator()(face& x, const face& y) const
    {
        if (x.size())
        {
            label j = 0;
            forAll(x, i)
            {
                x[i] = min(x[i], y[j]);

                j = y.rcIndex(j);
            }
        }
    }
};

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Are two lists identical either in forward or in reverse order.
bool Foam::localPointRegion::isDuplicate
(
    const face& f0,
    const face& f1,
    const bool forward
)
{
    if (f0.size() != f1.size())
    {
        return false;
    }

    label fp1 = f1.find(f0[0]);

    if (fp1 == -1)
    {
        return false;
    }

    forAll(f0, fp0)
    {
        if (f0[fp0] != f1[fp1])
        {
            return false;
        }

        if (forward)
        {
            fp1 = f1.fcIndex(fp1);
        }
        else
        {
            fp1 = f1.rcIndex(fp1);
        }
    }
    return true;
}


// Count regions per point
void Foam::localPointRegion::countPointRegions
(
    const polyMesh& mesh,
    const boolList& candidatePoint,
    const Map<label>& candidateFace,
    faceList& minRegion
)
{
    // Almost all will have only one so only
    // populate Map if more than one.
    labelList minPointRegion(mesh.nPoints(), -1);
    // From global point to local (multi-region) point numbering
    meshPointMap_.resize(candidateFace.size()/100);
    // From local (multi-region) point to regions
    DynamicList<labelList> pointRegions(meshPointMap_.size());

    // From faces with any duplicated point on it to local face
    meshFaceMap_.resize(meshPointMap_.size());

    forAllConstIters(candidateFace, iter)
    {
        const label facei = iter.key();

        if (!mesh.isInternalFace(facei))
        {
            const face& f = mesh.faces()[facei];

            if (minRegion[facei].empty())
            {
                FatalErrorInFunction
                    << "Face from candidateFace without minRegion set." << endl
                    << "Face:" << facei << " fc:" << mesh.faceCentres()[facei]
                    << " verts:" << f << abort(FatalError);
            }

            forAll(f, fp)
            {
                const label pointi = f[fp];

                // Even points which were not candidates for splitting might
                // be on multiple baffles that are being split so check.

                if (candidatePoint[pointi])
                {
                    const label region = minRegion[facei][fp];

                    if (minPointRegion[pointi] == -1)
                    {
                        minPointRegion[pointi] = region;
                    }
                    else if (minPointRegion[pointi] != region)
                    {
                        // Multiple regions for this point. Add.
                        const auto iter = meshPointMap_.cfind(pointi);
                        if (iter.good())
                        {
                            labelList& regions = pointRegions[iter.val()];
                            if (!regions.found(region))
                            {
                                label sz = regions.size();
                                regions.setSize(sz+1);
                                regions[sz] = region;
                            }
                        }
                        else
                        {
                            const label localPointi = meshPointMap_.size();
                            meshPointMap_.insert(pointi, localPointi);
                            labelList regions(2);
                            regions[0] = minPointRegion[pointi];
                            regions[1] = region;
                            pointRegions.append(regions);
                        }

                        const label meshFaceMapI = meshFaceMap_.size();
                        meshFaceMap_.insert(facei, meshFaceMapI);
                    }
                }
            }
        }
    }
    minPointRegion.clear();

    // Add internal faces that use any duplicated point. Can only have one
    // region!
    forAllConstIters(candidateFace, iter)
    {
        const label facei = iter.key();

        if (mesh.isInternalFace(facei))
        {
            const face& f = mesh.faces()[facei];

            forAll(f, fp)
            {
                // Note: candidatePoint test not really necessary but
                // speeds up rejection.
                if (candidatePoint[f[fp]] && meshPointMap_.found(f[fp]))
                {
                    const label meshFaceMapI = meshFaceMap_.size();
                    meshFaceMap_.insert(facei, meshFaceMapI);
                }
            }
        }
    }


    // Transfer to member data
    pointRegions.shrink();
    pointRegions_.setSize(pointRegions.size());
    forAll(pointRegions, i)
    {
        pointRegions_[i].transfer(pointRegions[i]);
    }

    // Compact minRegion
    faceRegions_.setSize(meshFaceMap_.size());
    forAllConstIters(meshFaceMap_, iter)
    {
        faceRegions_[iter()].labelList::transfer(minRegion[iter.key()]);

        //// Print a bit
        //{
        //    label facei = iter.key();
        //    const face& f = mesh.faces()[facei];
        //    Pout<< "Face:" << facei << " fc:" << mesh.faceCentres()[facei]
        //        << " verts:" << f << endl;
        //    forAll(f, fp)
        //    {
        //        Pout<< "    " << f[fp] << " min:" << faceRegions_[iter()][fp]
        //            << endl;
        //    }
        //    Pout<< endl;
        //}
    }

    // Compact region numbering
    // ? TBD.
}


void Foam::localPointRegion::calcPointRegions
(
    const polyMesh& mesh,
    const labelPairList& baffles,
    boolList& candidatePoint
)
{
    const label nBnd = mesh.nBoundaryFaces();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();


    syncTools::syncPointList
    (
        mesh,
        candidatePoint,
        orEqOp<bool>(),
        false               // nullValue
    );


    // Mark any face/boundaryFace/cell with a point on a candidate point.
    // - candidateFace does not necessary have to be a baffle!
    // - candidateFace is synchronised (since candidatePoint is)
    Map<label> candidateFace(2*nBnd);
    label candidateFacei = 0;

    Map<label> candidateCell(nBnd);
    label candidateCelli = 0;

    forAll(mesh.faces(), facei)
    {
        const face& f = mesh.faces()[facei];

        forAll(f, fp)
        {
            if (candidatePoint[f[fp]])
            {
                // Mark face
                if (candidateFace.insert(facei, candidateFacei))
                {
                    candidateFacei++;
                }

                // Mark cells
                if (candidateCell.insert(faceOwner[facei], candidateCelli))
                {
                    candidateCelli++;
                }

                if (mesh.isInternalFace(facei))
                {
                    label nei = faceNeighbour[facei];
                    if (candidateCell.insert(nei, candidateCelli))
                    {
                        candidateCelli++;
                    }
                }

                break;
            }
        }
    }



    // Get global indices for cells
    globalIndex globalCells(mesh.nCells());


    // Determine for every candidate face per point the minimum region
    // (global cell) it is connected to. (candidateFaces are the
    // only ones using a
    // candidate point so the only ones that can be affected)
    faceList minRegion(mesh.nFaces());
    forAllConstIters(candidateFace, iter)
    {
        label facei = iter.key();
        const face& f = mesh.faces()[facei];

        if (mesh.isInternalFace(facei))
        {
            label globOwn = globalCells.toGlobal(faceOwner[facei]);
            label globNei = globalCells.toGlobal(faceNeighbour[facei]);
            minRegion[facei].setSize(f.size(), min(globOwn, globNei));
        }
        else
        {
            label globOwn = globalCells.toGlobal(faceOwner[facei]);
            minRegion[facei].setSize(f.size(), globOwn);
        }
    }

    // Now minimize over all faces that are connected through internal
    // faces or through cells. This loop iterates over the max number of
    // cells connected to a point (=8 for hex mesh)

    while (true)
    {
        // Transport minimum from face across cell
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        Map<label> minPointValue(128);
        label nChanged = 0;
        forAllConstIters(candidateCell, iter)
        {
            minPointValue.clear();

            label celli = iter.key();
            const cell& cFaces = mesh.cells()[celli];

            // Determine minimum per point
            forAll(cFaces, cFacei)
            {
                label facei = cFaces[cFacei];

                if (minRegion[facei].size())
                {
                    const face& f = mesh.faces()[facei];

                    forAll(f, fp)
                    {
                        const label pointi = f[fp];
                        auto iter = minPointValue.find(pointi);

                        if (!iter.good())
                        {
                            minPointValue.insert(pointi, minRegion[facei][fp]);
                        }
                        else
                        {
                            label currentMin = iter();
                            iter() = min(currentMin, minRegion[facei][fp]);
                        }
                    }
                }
            }

            // Set face minimum from point minimum
            forAll(cFaces, cFacei)
            {
                label facei = cFaces[cFacei];

                if (minRegion[facei].size())
                {
                    const face& f = mesh.faces()[facei];

                    forAll(f, fp)
                    {
                        label minVal = minPointValue[f[fp]];

                        if (minVal != minRegion[facei][fp])
                        {
                            minRegion[facei][fp] = minVal;
                            nChanged++;
                        }
                    }
                }
            }
        }

        //Pout<< "nChanged:" << nChanged << endl;

        if (!returnReduceOr(nChanged))
        {
            break;
        }


        // Transport minimum across coupled faces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        SubList<face> l
        (
            minRegion,
            mesh.nBoundaryFaces(),
            mesh.nInternalFaces()
        );
        syncTools::syncBoundaryFaceList
        (
            mesh,
            l,
            minEqOpFace(),
            Foam::dummyTransform()  // dummy transformation
        );
        forAll(baffles, i)
        {
            label f0 = baffles[i].first();
            label f1 = baffles[i].second();
            minEqOpFace()(minRegion[f0], minRegion[f1]);
            minRegion[f1] = minRegion[f0];
        }
    }


    // Count regions per point
    countPointRegions(mesh, candidatePoint, candidateFace, minRegion);
    minRegion.clear();


    //// Print points with multiple regions. These points need to be duplicated.
    //forAllConstIters(meshPointMap_, iter)
    //{
    //    Pout<< "point:" << iter.key()
    //        << " coord:" << mesh.points()[iter.key()]
    //        << " regions:" << pointRegions_[iter.val()] << endl;
    //}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::localPointRegion::localPointRegion(const polyMesh& mesh)
:
    //nRegions_(0),
    meshPointMap_(0),
    pointRegions_(0),
    meshFaceMap_(0),
    faceRegions_(0)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Get any point on the outside which is on a non-coupled boundary
    boolList candidatePoint(mesh.nPoints(), false);

    forAll(patches, patchi)
    {
        if (!patches[patchi].coupled())
        {
            const polyPatch& pp = patches[patchi];

            forAll(pp.meshPoints(), i)
            {
                candidatePoint[pp.meshPoints()[i]] = true;
            }
        }
    }

    calcPointRegions(mesh, labelPairList(0), candidatePoint);
}


Foam::localPointRegion::localPointRegion
(
    const polyMesh& mesh,
    const labelList& candidatePoints
)
:
    //nRegions_(0),
    meshPointMap_(0),
    pointRegions_(0),
    meshFaceMap_(0),
    faceRegions_(0)
{
    boolList candidatePoint(mesh.nPoints(), false);

    forAll(candidatePoints, i)
    {
        candidatePoint[candidatePoints[i]] = true;
    }

    calcPointRegions(mesh, labelPairList(0), candidatePoint);
}


Foam::localPointRegion::localPointRegion
(
    const polyMesh& mesh,
    const labelPairList& baffles,
    const labelList& candidatePoints
)
:
    //nRegions_(0),
    meshPointMap_(0),
    pointRegions_(0),
    meshFaceMap_(0),
    faceRegions_(0)
{
    boolList candidatePoint(mesh.nPoints(), false);

    forAll(candidatePoints, i)
    {
        candidatePoint[candidatePoints[i]] = true;
    }

    calcPointRegions(mesh, baffles, candidatePoint);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return a list (in allPatch indices) with either -1 or the face label
// of the face that uses the same vertices.
Foam::labelList Foam::localPointRegion::findDuplicateFaces
(
    const primitiveMesh& mesh,
    const labelList& boundaryFaces
)
{
    // Addressing engine for all boundary faces.
    indirectPrimitivePatch allPatch
    (
        IndirectList<face>(mesh.faces(), boundaryFaces),
        mesh.points()
    );

    labelList duplicateFace(allPatch.size(), -1);
    // label nDuplicateFaces = 0;

    // Find all duplicate faces.
    forAll(allPatch, bFacei)
    {
        const face& f = allPatch.localFaces()[bFacei];

        // Get faces connected to f[0].
        // Check whether share all points with f.
        const labelList& pFaces = allPatch.pointFaces()[f[0]];

        forAll(pFaces, i)
        {
            label otherFacei = pFaces[i];

            if (otherFacei > bFacei)
            {
                const face& otherF = allPatch.localFaces()[otherFacei];

                if (isDuplicate(f, otherF, true))
                {
                    FatalErrorInFunction
                        << "Face:" << bFacei + mesh.nInternalFaces()
                        << " has local points:" << f << " at:"
                        << UIndirectList<point>(allPatch.localPoints(), f)
                        << " which are in same order as face:"
                        << otherFacei + mesh.nInternalFaces()
                        << " with local points:" << otherF
                        << abort(FatalError);
                }
                else if (isDuplicate(f, otherF, false))
                {
                    const label meshFace0 = boundaryFaces[bFacei];
                    const label meshFace1 = boundaryFaces[otherFacei];

                    if
                    (
                        duplicateFace[bFacei] != -1
                     || duplicateFace[otherFacei] != -1
                    )
                    {
                        FatalErrorInFunction
                            << "One of two duplicate faces already marked"
                            << " as duplicate." << nl
                            << "This means that three or more faces share"
                            << " the same points and this is illegal." << nl
                            << "Face:" << meshFace0
                            << " with local points:" << f << " at:"
                            << UIndirectList<point>(allPatch.localPoints(), f)
                            << " which are in same order as face:"
                            << meshFace1
                            << " with local points:" << otherF
                            << abort(FatalError);
                    }

                    duplicateFace[bFacei] = otherFacei;
                    duplicateFace[otherFacei] = bFacei;
                    // ++nDuplicateFaces;
                }
            }
        }
    }

    return duplicateFace;
}


Foam::List<Foam::labelPair> Foam::localPointRegion::findDuplicateFacePairs
(
    const polyMesh& mesh
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Faces to test: all boundary faces
    labelList testFaces
    (
        identity(mesh.nBoundaryFaces(), mesh.nInternalFaces())
    );

    // Find corresponding baffle face (or -1)
    const labelList duplicateFace(findDuplicateFaces(mesh, testFaces));

    // Convert into list of coupled face pairs (mesh face labels).
    DynamicList<labelPair> baffles(testFaces.size());

    forAll(duplicateFace, i)
    {
        label otherFacei = duplicateFace[i];

        if (otherFacei != -1 && i < otherFacei)
        {
            label meshFace0 = testFaces[i];
            label patch0 = patches.whichPatch(meshFace0);
            label meshFace1 = testFaces[otherFacei];
            label patch1 = patches.whichPatch(meshFace1);

            // Check for illegal topology. Should normally not happen!
            if
            (
                (patch0 != -1 && isA<processorPolyPatch>(patches[patch0]))
             || (patch1 != -1 && isA<processorPolyPatch>(patches[patch1]))
            )
            {
                FatalErrorInFunction
                    << "One of two duplicate faces is on"
                    << " processorPolyPatch."
                    << "This is not allowed." << nl
                    << "Face:" << meshFace0
                    << " fc:" << mesh.faceCentres()[meshFace0]
                    << " is on patch:" << patches[patch0].name()
                    << nl
                    << "Face:" << meshFace1
                    << " fc:" << mesh.faceCentres()[meshFace1]
                    << " is on patch:" << patches[patch1].name()
                    << abort(FatalError);
            }
            else
            {
                baffles.append(labelPair(meshFace0, meshFace1));
            }
        }
    }

    return List<labelPair>(std::move(baffles));
}


// Are two lists identical either in forward or in reverse order.
bool Foam::localPointRegion::isDuplicate
(
    const scalar mergeDist,
    const pointList& f0,
    const pointList& f1,
    const bool forward
)
{
    if (f0.size() != f1.size())
    {
        return false;
    }

    label fp1 = ListOps::find_if
    (
        f1,
        [f0, mergeDist](const point& p)
        {
            return mag(p - f0[0]) < mergeDist;
        }
    );

    if (fp1 == -1)
    {
        return false;
    }

    forAll(f0, fp0)
    {
        if (mag(f0[fp0] - f1[fp1]) >= mergeDist)
        {
            return false;
        }

        if (forward)
        {
            fp1 = f1.fcIndex(fp1);
        }
        else
        {
            fp1 = f1.rcIndex(fp1);
        }
    }
    return true;
}


void Foam::localPointRegion::findDuplicateFacePairs
(
    const scalar mergeDist,
    const polyMesh& mesh,
    const globalIndex& globalFaces,
    const indirectPrimitivePatch& allPatch,
    DynamicList<labelPair>& localBaffles,   // incidices into allPatch
    DynamicList<labelPair>& globalBaffles   // indices in globalFaces
)
{
    // Return lists of face pairs sharing the same points:
    // - localBaffles : indices into boundaryFaces
    // - globalBaffles : indices in global face numbering

    const labelList& boundaryFaces = allPatch.addressing();
    if (globalFaces.localSize() != boundaryFaces.size())
    {
        FatalErrorInFunction
            << "Size of globalFaces " << globalFaces.localSize()
            << " does not match size of boundaryFaces " << boundaryFaces.size()
            << abort(FatalError);
    }
    const auto& localFaces = allPatch.localFaces();

    localBaffles.setCapacity(allPatch.size());
    globalBaffles.setCapacity(allPatch.size());


    // Find local duplicates
    // ~~~~~~~~~~~~~~~~~~~~~

    bitSet isLocalBaffle(allPatch.size());

    // Find all duplicate faces.
    forAll(localFaces, bFacei)
    {
        const face& f = localFaces[bFacei];

        // Get faces connected to f[0].
        // Check whether share all points with f.
        const labelList& pFaces = allPatch.pointFaces()[f[0]];

        forAll(pFaces, i)
        {
            label otherFacei = pFaces[i];

            if (otherFacei > bFacei)
            {
                const face& otherF = allPatch.localFaces()[otherFacei];
                if (isDuplicate(f, otherF, false))
                {
                    localBaffles.append(labelPair(bFacei, otherFacei));
                    isLocalBaffle.set(bFacei);
                    isLocalBaffle.set(otherFacei);
                }
            }
        }
    }


    // Find remote duplicates
    // ~~~~~~~~~~~~~~~~~~~~~~
    // - compare all face centres
    // - for faces within mergeDist: compare the points


    // Geometric merging for face centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // (opposite boundary faces have opposite numbering but hopefully the same
    //  face centre. Note: this is the same test as done by processor
    //  face matching so should be robust enough)


    // Gather all face centres in processor order
    pointField allFaceCentres;
    globalFaces.gather(allPatch.faceCentres(), allFaceCentres);

    labelList allToMerged;
    mergePoints
    (
        allFaceCentres,
        mergeDist,
        false,      //const bool verbose,
        allToMerged
    );

    // Extract duplicate faces by inverting allToMerged to get number of
    // duplicates per master face. This is similar to ListOps'
    // invertOneToMany except it leaves out non-duplicates. Could probably
    // be done with two lists instead since we're only interested in
    // two faces at the same location.
    labelList nDup(allFaceCentres.size(), 0);
    forAll(allToMerged, globali)
    {
        const label master = allToMerged[globali];
        if (nDup[master] >= 2)
        {
            FatalErrorInFunction
                << "More than two faces at same location. This is not allowed."
                << abort(FatalError);
        }
        nDup[master]++;
    }

    // Collect duplicate faces (same master face)
    List<labelList> dupFaces(allFaceCentres.size());
    forAll(nDup, globali)
    {
        if (nDup[globali] == 2)
        {
            dupFaces[globali].resize_nocopy(nDup[globali]);
        }
        nDup[globali] = 0;
    }
    forAll(dupFaces, globali)
    {
        const label master = allToMerged[globali];
        if (dupFaces[master].size())
        {
            dupFaces[master][nDup[master]++] = globali;

            // if (nDup[master] == 2)
            // {
            //     Pout<< "Found duplicate faces at global face:" << globali
            //         << " globals:" << dupFaces[master]
            //         << " p0:" << globalFaces.whichProcID(dupFaces[master][0])
            //         << " p1:" << globalFaces.whichProcID(dupFaces[master][1])
            //         << " with face centres:"
            //         << UIndirectList<point>(allFaceCentres, dupFaces[master])
            //         << endl;
            // }
        }
    }

    // At this point:
    // - per master face (not same as global face index)
    // - we have the list of global faces


    // Get data from master in global face numbering i.e. expand using
    // allToMerged. Note: this is a bit inefficient since we send the same
    // face centre multiple times but should be ok since only for boundary
    // faces and only for those that are duplicates.
    List<labelList> allDupFaces
    (
        UIndirectList<labelList>(dupFaces, allToMerged)
    );

    // Send back to originating processor.
    List<labelList> dupFacesScatter(allPatch.size());
    globalFaces.scatter(allDupFaces, dupFacesScatter);

    // Send back the vertices of the duplicate faces to master processor
    List<pointList> faceCoords(dupFacesScatter.size());
    forAll(dupFacesScatter, facei)
    {
        if (dupFacesScatter[facei].size())
        {
            const auto& faceIDs = dupFacesScatter[facei];
            forAll(faceIDs, i)
            {
                if (globalFaces.isLocal(faceIDs[i]))
                {
                    const label localFacei = globalFaces.toLocal(faceIDs[i]);
                    if (!isLocalBaffle[localFacei])
                    {
                        // Face is local so can get the points
                        faceCoords[facei] =
                            allPatch[localFacei].points(allPatch.points());
                        break;
                    }
                }
            }
        }
    }
    globalFaces.gatherInplace(faceCoords);

    // Do comparison on master processor and mark dups
    forAll(allToMerged, globali)
    {
        const auto& faceIDs = dupFaces[allToMerged[globali]];
        for (label i = 0; i < faceIDs.size(); i++)
        {
            const auto& f0Verts = faceCoords[faceIDs[i]];
            if (f0Verts.empty())
            {
                break;
            }
            for (label j = i+1; j < faceIDs.size(); j++)
            {
                const auto& f1Verts = faceCoords[faceIDs[j]];
                if (f1Verts.empty())
                {
                    break;
                }

                if (isDuplicate(mergeDist, f0Verts, f1Verts, false))
                {
                    // Mark as duplicate of master
                    // Pout<< "Global face0:" << faceIDs[i]
                    //     << " has local points:" << f0Verts << nl
                    //     << "Global face1:" << faceIDs[j]
                    //     << " has local points:" << f1Verts
                    //     << endl;
                }
                else
                {
                    // Unmark face as duplicate
                    dupFaces[allToMerged[globali]].clear();
                }
            }
        }
    }

    // Re-scatter the global face IDs to originating processors
    // (since we've unmarked the non-duplicate ones)
    allDupFaces = UIndirectList<labelList>(dupFaces, allToMerged);
    globalFaces.scatter(allDupFaces, dupFacesScatter);

    // Extract as face pairs
    forAll(dupFacesScatter, facei)
    {
        const labelList& faceIDs = dupFacesScatter[facei];
        if (faceIDs.size() == 2)
        {
            const label globFace0 = faceIDs[0];
            const label globFace1 = faceIDs[1];
            if (globalFaces.isLocal(globFace0))
            {
                globalBaffles.append(labelPair(globFace0, globFace1));
            }
            else if (globalFaces.isLocal(globFace1))
            {
                globalBaffles.append(labelPair(globFace1, globFace0));
            }
            else
            {
                // Should not happen since at least one should be local
                FatalErrorInFunction << "Problem" << abort(FatalError);
            }
        }
        else if (faceIDs.size())
        {
            // Should not happen since we only allow pairs
            FatalErrorInFunction
                << "More than two faces at same location. This is not allowed."
                << abort(FatalError);
        }
    }
    // Now both processors know the same pair of faces. Look up the
    // one that is local and work out where the other one is.
}


void Foam::localPointRegion::updateMesh(const mapPolyMesh& map)
{
    {
        Map<label> newMap(2*meshFaceMap_.size());

        forAllConstIters(meshFaceMap_, iter)
        {
            const label newFacei = map.reverseFaceMap()[iter.key()];

            if (newFacei >= 0)
            {
                newMap.insert(newFacei, iter.val());
            }
        }
        meshFaceMap_.transfer(newMap);
    }
    {
        Map<label> newMap(2*meshPointMap_.size());

        forAllConstIters(meshPointMap_, iter)
        {
            const label newPointi = map.reversePointMap()[iter.key()];

            if (newPointi >= 0)
            {
                newMap.insert(newPointi, iter.val());
            }
        }

        meshPointMap_.transfer(newMap);
    }
}


// ************************************************************************* //
