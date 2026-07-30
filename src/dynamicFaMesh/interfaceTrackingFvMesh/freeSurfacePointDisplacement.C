/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Zeljko Tukovic, FSB Zagreb.
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

#include "interfaceTrackingFvMesh.H"
#include "primitivePatchInterpolation.H"
#include "processorPolyPatch.H"
#include "wallFvPatch.H"
#include "coordinateSystem.H"
#include "unitConversion.H"
#include "scalarMatrices.H"
#include "tensor2D.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::vectorField>
Foam::interfaceTrackingFvMesh::pointDisplacement()
{
    const pointField& points = aMesh().patch().localPoints();
    const labelListList& pointFaces = aMesh().patch().pointFaces();

    auto tdisplacement = tmp<vectorField>::New(points.size(), Zero);
    auto& displacement = tdisplacement.ref();

    // Calculate displacement of internal points.
    // pointAreaNormals() are reordered into patch()-local numbering so they
    // align with the patch-local point loops below.
    const tmp<vectorField> tPointNormals(patchPointAreaNormals());
    const vectorField& pointNormals = tPointNormals();
    labelList internalPoints = aMesh().internalPoints();

    for (const label curPoint : internalPoints)
    {
        const labelList& curPointFaces = pointFaces[curPoint];

        vectorField lsPoints(curPointFaces.size(), Zero);

        for (label i=0; i<curPointFaces.size(); i++)
        {
            label curFace = curPointFaces[i];

            lsPoints[i] = controlPoints()[curFace];
        }

        vectorField pointAndNormal
        (
            lsPlanePointAndNormal
            (
                lsPoints,
                points[curPoint],
                pointNormals[curPoint]
            )
        );

        vector& P = pointAndNormal[0];
        vector& N = pointAndNormal[1];

        displacement[curPoint] =
            pointsDisplacementDir()[curPoint]
           *((P - points[curPoint])&N)
           /(pointsDisplacementDir()[curPoint]&N);
    }

    // Mirror control points
    FieldField<Field, vector> patchMirrorPoints(aMesh().boundary().size());

    // Old faMesh points
    vectorField oldPoints(aMesh().nPoints(), Zero);
    const labelList& meshPoints = aMesh().patch().meshPoints();
    forAll(oldPoints, pI)
    {
        oldPoints[pI] =
            mesh().oldPoints()[meshPoints[pI]];
    }

    // In finiteArea the reduced-order boundary faces ARE the finite-area
    // boundary edges (2-point faces in reduced-order/patch-local point
    // numbering). A boundary patch's edges are therefore its faces, ordered by
    // the finite-volume patch (patchSlice) ordering, its owning finite-area
    // face is the patch faceCells(), and the neighbour-patch coupling data is
    // keyed the same way.
    const faceList& faEdges = aMesh().faces();

    // Reduced-order face point labels -> patch().localPoints() numbering
    const labelList& faPointMap = aMesh().patchPointMap();

    forAll(patchMirrorPoints, patchI)
    {
        // Use the polyPatch view: a constraint fvPatch (e.g. the 2D empty
        // patch) reports size()/faceCells() as empty, but its reduced-order
        // boundary faces (rim edges) are real and must be mirrored.
        const polyPatch& fap = aMesh().boundary()[patchI].patch();
        const label patchStart = fap.start();
        const label patchSize = fap.size();

        patchMirrorPoints.set(patchI, new vectorField(patchSize, Zero));

        vectorField N
        (
            aMesh().ngbPolyPatchFaceNormals(patchI)
        );

        // Finite-area face (reduced-order cell) owning each boundary edge
        const labelUList& eFaces = fap.faceCells();

        // Correct N according to specified contact angle
        if (contactAnglePtr_)
        {
            label ngbPolyPatchID = aMesh().ngbPolyPatchIndex(patchI);

            if (ngbPolyPatchID != -1)
            {
                if
                (
                    mesh().boundary()[ngbPolyPatchID].type()
                 == wallFvPatch::typeName
                )
                {
                    scalar rotAngle = degToRad
                    (
                        gAverage
                        (
                            90
                          - contactAnglePtr_->boundaryField()[patchI]
                        )
                    );

                    const vectorField& pEdgN =
                        aMesh().edgeAreaNormals().boundaryField()[patchI];

                    vectorField rotationAxis( N^pEdgN );

                    forAll(rotationAxis, edgeI)
                    {
                        const face& faEdge = faEdges[patchStart + edgeI];
                        vector e =
                            oldPoints[faPointMap[faEdge[1]]]
                          - oldPoints[faPointMap[faEdge[0]]];

                        // Adjust direction
                        rotationAxis[edgeI] =
                            e*(e&rotationAxis[edgeI])
                           /mag((e&rotationAxis[edgeI]));
                    }
                    rotationAxis /= mag(rotationAxis) + SMALL;

                    vectorField rotationAxis2 = rotationAxis;
                    forAll(rotationAxis2, edgeI)
                    {
                        rotationAxis2[edgeI] =
                            (N[edgeI]^facesDisplacementDir()[eFaces[edgeI]]);

                        // Adjust direction
                        rotationAxis2[edgeI] =
                            rotationAxis2[edgeI]
                           *(rotationAxis2[edgeI]&rotationAxis[edgeI])
                           /(
                                mag((rotationAxis2[edgeI]&rotationAxis[edgeI]))
                              + SMALL
                            );
                    }
                    rotationAxis2 /= mag(rotationAxis2) + SMALL;

                    // Rodrigues' rotation formula
                    N = N*cos(rotAngle)
                      + rotationAxis*(rotationAxis & N)*(1 - cos(rotAngle))
                      + (rotationAxis^N)*sin(rotAngle);

                    N /= mag(N) + SMALL;

                    N = (rotationAxis^N);

                    N = (N^rotationAxis2);

                    N /= mag(N) + SMALL;
                }
            }
        }

        vectorField peCentres(patchSize, Zero);
        forAll(peCentres, edgeI)
        {
            const face& faEdge = faEdges[patchStart + edgeI];
            peCentres[edgeI] =
                0.5*(points[faPointMap[faEdge[0]]]
              + points[faPointMap[faEdge[1]]]);
        }

        vectorField delta
        (
            vectorField(controlPoints(), eFaces)
          - peCentres
        );

        patchMirrorPoints[patchI] =
            peCentres + ((I - 2*N*N)&delta);
    }

    // Calculate displacement of boundary points
    labelList boundaryPoints = aMesh().boundaryPoints();

    const labelListList& pointEdges = aMesh().patch().pointEdges();

    // Parent-patch edge -> reduced-order boundary face (global mesh-face id),
    // -1 for an internal edge. Bridges the parent-patch point/edge topology to
    // the reduced-order boundary patches that carry the mirror points.
    const labelList& patchEdgeBndFace = aMesh().patchEdgeToBoundaryFace();
    const polyBoundaryMesh& faBoundary = aMesh().boundaryMesh();

    for (const label curPoint : boundaryPoints)
    {
        if (motionPointsMask()[curPoint] == 1)
        {
            // Calculating mirror points
            const labelList& curPointEdges = pointEdges[curPoint];

            vectorField mirrorPoints(2, Zero);

            label counter = -1;

            forAll(curPointEdges, edgeI)
            {
                const label curEdge = curPointEdges[edgeI];
                const label bndFacei = patchEdgeBndFace[curEdge];

                // A finite-area boundary edge (coincides with a reduced-order
                // boundary face); internal edges map to -1.
                if (bndFacei != -1)
                {
                    const label patchID = faBoundary.whichPatch(bndFacei);
                    const label edgeID =
                        bndFacei - faBoundary[patchID].start();

                    mirrorPoints[++counter] =
                        patchMirrorPoints[patchID][edgeID];
                }
            }

            // Calculating LS plane fit
            const labelList& curPointFaces = pointFaces[curPoint];

            vectorField lsPoints
            (
                curPointFaces.size() + mirrorPoints.size(),
                Zero
            );

            counter = -1;

            for (label i=0; i<curPointFaces.size(); i++)
            {
                label curFace = curPointFaces[i];

                lsPoints[++counter] = controlPoints()[curFace];
            }

            for (label i=0; i<mirrorPoints.size(); i++)
            {
                lsPoints[++counter] = mirrorPoints[i];
            }

            vectorField pointAndNormal
            (
                lsPlanePointAndNormal
                (
                    lsPoints,
                    points[curPoint],
                    pointNormals[curPoint]
                )
            );

            vector& P = pointAndNormal[0];
            vector& N = pointAndNormal[1];

            displacement[curPoint] =
                pointsDisplacementDir()[curPoint]
               *((P - points[curPoint])&N)
               /(pointsDisplacementDir()[curPoint]&N);
        }
    }

    // Calculate displacement of processor patch points.
    //
    // The legacy path exchanged the least-squares stencil control points
    // across processorFaPatches (via nonGlobalPatchPoints()/neighbPoints() and
    // the per-point boundaryProcessorFaPatchPoints mirror construction). That
    // finite-area processor point coupling is not yet ported to finiteArea -
    // guard it loudly rather than silently drop the neighbour contribution.
    // Serial runs have no processor patches and never enter this path.
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "Parallel free-surface point displacement (processorFaPatch"
            << " least-squares stencil exchange) is not yet ported to"
            << " finiteArea. Run serially." << nl << abort(FatalError);
    }


    // Calculate displacement of global processor patch points
    if (aMesh().globalData().nGlobalPoints() > 0)
    {
        const labelList& spLabels =
            aMesh().globalData().sharedPointLabels();

        const labelList& addr = aMesh().globalData().sharedPointAddr();

        for (label k=0; k<aMesh().globalData().nGlobalPoints(); k++)
        {
            List<List<vector>> procLsPoints(Pstream::nProcs());

            label curSharedPointIndex = addr.find(k);

            if (curSharedPointIndex != -1)
            {
                label curPoint = spLabels[curSharedPointIndex];

                const labelList& curPointFaces = pointFaces[curPoint];

                procLsPoints[Pstream::myProcNo()] =
                    List<vector>(curPointFaces.size());

                forAll(curPointFaces, faceI)
                {
                    label curFace = curPointFaces[faceI];

                    procLsPoints[Pstream::myProcNo()][faceI] =
                        controlPoints()[curFace];
                }
            }

            Pstream::allGatherList(procLsPoints);

            if (curSharedPointIndex != -1)
            {
                label curPoint = spLabels[curSharedPointIndex];

                label nAllPoints = 0;
                forAll(procLsPoints, procI)
                {
                    nAllPoints += procLsPoints[procI].size();
                }

                vectorField allPoints(nAllPoints, Zero);

                label counter = 0;
                forAll(procLsPoints, procI)
                {
                    forAll(procLsPoints[procI], pointI)
                    {
                        allPoints[counter++] =
                            procLsPoints[procI][pointI];
                    }
                }

                vectorField pointAndNormal
                (
                    lsPlanePointAndNormal
                    (
                        allPoints,
                        points[curPoint],
                        pointNormals[curPoint]
                    )
                );

                const vector& P = pointAndNormal[0];
                const vector& N = pointAndNormal[1];

                displacement[curPoint] =
                    pointsDisplacementDir()[curPoint]
                   *((P - points[curPoint])&N)
                   /(pointsDisplacementDir()[curPoint]&N);
            }
        }
    }

    return tdisplacement;
}


Foam::tmp<Foam::vectorField>
Foam::interfaceTrackingFvMesh::lsPlanePointAndNormal
(
    const vectorField& points,
    const vector& origin,
    const vector& axis
) const
{
    // LS in local CS
    vector dir = (points[0] - origin);
    dir -= axis*(axis&dir);
    dir /= mag(dir);
    coordinateSystem cs("cs", origin, axis, dir);

    vectorField localPoints(cs.localPosition(points));
    vector avgLocalPoint = average(localPoints);

    // scalarField W = 1.0/(mag(points - origin) + SMALL);
    scalarField W(points.size(), scalar(1));

    const label nCoeffs = 2;
    scalarRectangularMatrix M(points.size(), nCoeffs, Zero);

    scalar L = 2*max(mag(localPoints-avgLocalPoint));
    for (label i=0; i<localPoints.size(); i++)
    {
        M[i][0] = (localPoints[i].x() - avgLocalPoint.x())/L;
        M[i][1] = (localPoints[i].y() - avgLocalPoint.y())/L;
    }

    // Applying weights
    for (label i=0; i<M.n(); i++)
    {
        for (label j=0; j<M.m(); j++)
        {
            M[i][j] *= W[i];
        }
    }

    tensor2D lsM(Zero);
    // scalarSquareMatrix lsM(nCoeffs, Zero);

    for (label i=0; i<nCoeffs; i++)
    {
        for (label j=0; j<nCoeffs; j++)
        {
            for (label k=0; k<M.n(); k++)
            {
                lsM[i*nCoeffs+j] += M[k][i]*M[k][j];
                // lsM(i,j) += M[k][i]*M[k][j];
            }
        }
    }

    // Calculate inverse
    tensor2D invLsM = inv(lsM);

    scalarRectangularMatrix curInvMatrix(nCoeffs, points.size(), Zero);

    for (label i=0; i<nCoeffs; i++)
    {
        for (label j=0; j<M.n(); j++)
        {
            for (label k=0; k<nCoeffs; k++)
            {
                curInvMatrix[i][j] += invLsM[i*nCoeffs+k]*M[j][k]*W[j];
                // curInvMatrix[i][j] += invLsM[i][k]*M[j][k]*W[j];
            }
        }
    }

    scalarField coeffs(nCoeffs, Zero);
    scalarField source(points.size(), Zero);

    for (label i=0; i<points.size(); i++)
    {
        source[i] = (localPoints[i].z() - avgLocalPoint.z())/L;
    }

    for (label i=0; i<nCoeffs; i++)
    {
        for (label j=0; j<source.size(); j++)
        {
            coeffs[i] += curInvMatrix[i][j]*source[j];
        }
    }

    vector n0(-coeffs[0], -coeffs[1], 1.0);
    n0 = cs.globalVector(n0);
    n0 /= mag(n0);

    vector p0 = avgLocalPoint;
    p0 = cs.globalPosition(p0);

    auto tpointAndNormal = tmp<vectorField>::New(2);
    auto& pointAndNormal = tpointAndNormal.ref();

    pointAndNormal[0] = p0;
    pointAndNormal[1] = n0;

    return tpointAndNormal;
}


// ************************************************************************* //
