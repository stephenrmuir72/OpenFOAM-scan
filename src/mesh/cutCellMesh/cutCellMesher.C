/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
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

#include "cutCellMesher.H"
#include "cutCellCutter.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "labelIOList.H"
#include "cutCellMerge.H"
#include "removeFaces.H"
#include "labelListIOList.H"
#include "cellZone.H"
#include "faceZone.H"
#include "treeDataTriSurface.H"
#include "treeBoundBox.H"
#include "DynamicList.H"
#include "boundBox.H"
#include "pointIOField.H"
#include "boolList.H"
#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Local Functions  * * * * * * * * * * * * * * //

namespace
{

//- Signed distance of a point from the plane through 'refPoint' with normal n
inline scalar planeDist(const vector& n, const point& refPoint, const point& x)
{
    return n & (x - refPoint);
}


//- True when two planes are the same plane, allowing for the last bits
inline bool samePlane
(
    const vector& n1,
    const point& p1,
    const vector& n2,
    const point& p2,
    const scalar tol
)
{
    return
        (n1 & n2) > 1 - tol
     && mag(n1 & (p2 - p1)) <= tol;
}


//- A fixed order for planes, so that a result does not depend on the order the
//- cells that produced them happen to be in
inline bool planeLess(const plane& a, const plane& b)
{
    const vector na(a.normal());
    const vector nb(b.normal());

    for (direction i = 0; i < 3; i++)
    {
        if (mag(na[i] - nb[i]) > SMALL)
        {
            return na[i] < nb[i];
        }
    }

    return (a.normal() & a.refPoint()) < (b.normal() & b.refPoint());
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::cutCellMesher::cutCellMesher
(
    polyMesh& mesh,
    const PtrList<cutCellSurface>& surfaces
)
:
    mesh_(mesh),
    surfaces_(surfaces),
    tol_(1e-9*mag(mesh.bounds().max() - mesh.bounds().min())),
    info_(),
    backgroundCell_()
{
    if (tol_ <= 0)
    {
        tol_ = SMALL;
    }
}


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

Foam::label Foam::cutCellMesher::cellZoneID(const word& name)
{
    cellZoneMesh& zones = mesh_.cellZones();

    //- Finds an existing zone by name, or creates an empty one
    zones(name, true);

    return zones.findZoneID(name);
}


Foam::label Foam::cutCellMesher::faceZoneID(const word& name)
{
    faceZoneMesh& zones = mesh_.faceZones();

    zones(name, true);

    return zones.findZoneID(name);
}


Foam::cutCellPolyhedron Foam::cutCellMesher::cellPolyhedron
(
    const polyMesh& mesh,
    const label celli
)
{
    const labelList& cPoints = mesh.cellPoints()[celli];

    //- Mesh point label to local index
    Map<label> toLocal(2*cPoints.size());

    pointField pts(cPoints.size());

    forAll(cPoints, i)
    {
        toLocal.insert(cPoints[i], i);
        pts[i] = mesh.points()[cPoints[i]];
    }

    const faceList& fcs = mesh.faces();

    //- The cell's faces, outward oriented. A cell's face list is not ordered
    //- relative to the cell, so the owner decides whether to reverse it.
    faceList cellFcs(mesh.cells()[celli].size());

    forAll(cellFcs, i)
    {
        const label facei = mesh.cells()[celli][i];
        const face& f = fcs[facei];

        face& lf = cellFcs[i];
        lf.setSize(f.size());

        if (mesh.faceOwner()[facei] == celli)
        {
            forAll(f, k)
            {
                lf[k] = toLocal[f[k]];
            }
        }
        else
        {
            forAll(f, k)
            {
                lf[k] = toLocal[f[f.size() - 1 - k]];
            }
        }
    }

    return cutCellPolyhedron(pts, cellFcs);
}


void Foam::cutCellMesher::snapPoints()
{
    snapDelta_.setSize(mesh_.points().size(), Zero);

    if (vertexSnapFraction_ <= 0)
    {
        return;
    }

    if (surfaces_.size() != 1)
    {
        FatalErrorInFunction
            << "Cutting with " << surfaces_.size() << " surfaces is not"
            << " supported in this slice. Cells the surface does not cross still"
            << " need an inside/outside label, so give exactly one geometry"
            << " entry." << exit(FatalError);
    }

    const cutCellSurface& surf = surfaces_[0];

    const pointField& pts = mesh_.points();

    const scalar f = vertexSnapFraction_;

    //- The domain is the user's, so a point on a boundary patch never moves.
    //- That is also what keeps the total volume: a closed mesh's volume depends
    //- on its boundary alone, so interior moves leave it unchanged.
    boolList onBoundary(pts.size(), false);

    forAll(mesh_.boundaryMesh(), patchi)
    {
        const faceList& patchFaces = mesh_.boundaryMesh()[patchi];

        forAll(patchFaces, pf)
        {
            const face& pF = patchFaces[pf];

            forAll(pF, k)
            {
                onBoundary[pF[k]] = true;
            }
        }
    }

    //- The planes pulling each point, as a chain per point. A pull is one
    //- (point, plane) pair: a cut crossing an edge near its endpoint.
    Map<label> firstPull(pts.size()/64 + 16);
    DynamicList<label> nextPull(pts.size()/64 + 16);
    DynamicList<label> pullPoint(pts.size()/64 + 16);
    DynamicList<vector> pullNormal(pts.size()/64 + 16);
    DynamicList<point> pullRef(pts.size()/64 + 16);

    auto addPull = [&](const label pointi, const vector& n, const point& ref)
    {
        //- One facet is several triangles, so the same plane can pull the same
        //- point more than once. It is one pull, or the resolution would weigh
        //- that plane twice.
        for
        (
            label pulli = (firstPull.found(pointi) ? firstPull[pointi] : -1);
            pulli != -1;
            pulli = nextPull[pulli]
        )
        {
            if
            (
                (pullNormal[pulli] & n) > 1 - 1e-6
             && mag(pullNormal[pulli] & (pullRef[pulli] - ref)) <= tol_
            )
            {
                return;
            }
        }

        const label pulli = pullPoint.size();

        pullPoint.append(pointi);
        pullNormal.append(n);
        pullRef.append(ref);
        nextPull.append(firstPull.found(pointi) ? firstPull[pointi] : -1);
        firstPull.set(pointi, pulli);
    };

    //- A cut crossing within a fraction of the edge's length of one of its
    //- endpoints is collapsed onto that endpoint. The measure is the fraction
    //- along the edge rather than a distance, which is what 05 section 7.1.2
    //- describes and what makes the test independent of the refinement level.
    forAll(mesh_.cells(), celli)
    {
        const cutCellPolyhedron cell(cellPolyhedron(mesh_, celli));

        const boundBox bb(cell.points(), false);

        const labelList cand
        (
            surf.search().tree().findBox(treeBoundBox(bb.min(), bb.max()))
        );

        if (cand.empty())
        {
            continue;
        }

        const pointField& cPts = cell.points();

        forAll(cand, i)
        {
            const pointField tri(surf.triangle(cand[i]));

            const vector nRaw = (tri[1] - tri[0]) ^ (tri[2] - tri[0]);

            if (mag(nRaw) <= VSMALL)
            {
                continue;
            }

            const vector n = nRaw/mag(nRaw);
            const point ref = tri[0];

            const pointField poly(triangleBoxPolygon(tri, bb, tol_));

            if (poly.size() < 3)
            {
                continue;
            }

            //- Twice the signed area of a triangle about the facet's normal
            auto triArea2 = [&](const point& a, const point& b, const point& c)
            {
                return ((b - a) ^ (c - a)) & n;
            };

            scalar polyArea2 = 0;

            forAll(poly, k)
            {
                polyArea2 += triArea2(poly[0], poly[k], poly[k + 1 < poly.size() ? k + 1 : 0]);
            }

            //- Every edge of the cell, each of them seen twice, once from each
            //- of the faces that share it
            forAll(cell.faces(), cfi)
            {
                const face& cF = cell.faces()[cfi];

                forAll(cF, k)
                {
                    const label a = cF[k];
                    const label b = cF[cF.fcIndex(k)];

                    const scalar da = n & (cPts[a] - ref);
                    const scalar db = n & (cPts[b] - ref);

                    //- The plane has to cross the edge, not touch it
                    if (da*db >= 0)
                    {
                        continue;
                    }

                    const scalar t = da/(da - db);

                    label endLocal = -1;

                    if (t <= f)
                    {
                        endLocal = a;
                    }
                    else if (t >= 1 - f)
                    {
                        endLocal = b;
                    }
                    else
                    {
                        continue;
                    }

                    //- Only a crossing the facet really contains counts: a plane
                    //- can pass near a vertex with its facet somewhere else, and
                    //- then there is nothing there to collapse.
                    const point x(cPts[a] + t*(cPts[b] - cPts[a]));

                    scalar area2 = 0;

                    forAll(poly, pk)
                    {
                        area2 += triArea2
                        (
                            x, poly[pk], poly[pk + 1 < poly.size() ? pk + 1 : 0]
                        );
                    }

                    if (mag(area2 - polyArea2) > 1e-6*max(mag(polyArea2), VSMALL))
                    {
                        continue;
                    }

                    addPull(mesh_.cellPoints()[celli][endLocal], n, ref);
                }
            }
        }
    }

    //- Resolve each pulled point once, so that every cell it belongs to sees
    //- the same move, and then move it
    pointField newPts(pts);

    forAllConstIters(firstPull, iter)
    {
        const label pointi = iter.key();

        if (onBoundary[pointi])
        {
            continue;
        }

        DynamicList<vector> normals;
        DynamicList<point> refs;

        for (label pulli = iter(); pulli != -1; pulli = nextPull[pulli])
        {
            normals.append(pullNormal[pulli]);
            refs.append(pullRef[pulli]);
        }

        const point& P = pts[pointi];

        point target = P;

        if (normals.size() == 1)
        {
            //- Exactly onto the plane, which is the common case
            target = P - (normals[0] & (P - refs[0]))*normals[0];
        }
        else
        {
            //- The smallest move that best satisfies all of them: for two
            //- planes the closest point on their line, for three their corner.
            //- A small multiple of the identity keeps the solve well posed when
            //- they have no common point, and leaves the answer within
            //- round-off of the exact one when they have.
            symmTensor A = Zero;
            vector b = Zero;

            forAll(normals, i)
            {
                A += sqr(normals[i]);
                b += normals[i]*(normals[i] & refs[i]);
            }

            const scalar eps = 1e-12;

            const tensor M
            (
                A.xx() + eps, A.xy(), A.xz(),
                A.xy(), A.yy() + eps, A.yz(),
                A.xz(), A.yz(), A.zz() + eps
            );

            //- The solve is a division by that determinant, and it is the
            //- division that can trap: two pulls a fraction of a degree
            //- apart give a determinant at the scale of the formula's own
            //- rounding, and dividing by it raised SIGFPE and aborted a run
            //- whose geometry is sound. Measured: `07_sphere_cut`'s 63
            //- solves bottom out at 1.4e-14, a pair about seven degrees
            //- apart, while `10_box_regions` gave 1.0e-17 and then exactly 0
            //- for a pair a fifth of a degree apart. The normals are unit
            //- vectors, so those determinants are scale free and one
            //- threshold fits every mesh; it sits at about 1.8 degrees. An
            //- ill-conditioned pair leaves the point where it is, which is
            //- what the move bound does with the runaway answer it would
            //- otherwise produce, and what the collapse's policy already
            //- says: the sliver stays and is reported.
            if (mag(det(M)) <= 1e-15)
            {
                info_.nSnapRejected++;
                continue;
            }

            target = inv(M) & (b + eps*P);
        }

        const vector delta(target - P);

        //- A resolved point can be asked to move far. A single pull cannot: the
        //- crossing it collapses onto lies within vertexSnapFraction of the
        //- edge's endpoint, so the move is at most that fraction of that edge.
        //- The multi-plane solve has no such bound - two nearly parallel planes
        //- meet along a line that runs away as their dihedral angle closes, and
        //- a coarsely tessellated curved surface puts such pairs in one cell.
        //- Measured on a 320-triangle sphere in a 16-cubed background: a move of
        //- 0.0209 where the cell is 0.0625, a third of a cell, which distorted
        //- the cells until no piece contained a point of its own cell and the
        //- assembly aborted with no mesh written at all.
        //-
        //- The collapse exists to remove a sliver of the order of the trigger,
        //- so a move of the order of the trigger is the largest that can be
        //- doing that job. Beyond it the point is left where it is: the sliver
        //- stays and is reported, for the agglomeration to absorb, which is a
        //- worse mesh than a collapse would have given but a mesh nonetheless.
        if (maxSnapMoveFraction_ > 0)
        {
            scalar shortestEdge = GREAT;

            forAll(mesh_.pointEdges()[pointi], pei)
            {
                shortestEdge = min
                (
                    shortestEdge,
                    mesh_.edges()[mesh_.pointEdges()[pointi][pei]].mag(pts)
                );
            }

            if (mag(delta) > maxSnapMoveFraction_*shortestEdge)
            {
                info_.nSnapRejected++;
                continue;
            }
        }

        if (mag(delta) > SMALL)
        {
            newPts[pointi] = target;
            snapDelta_[pointi] = delta;
            info_.nSnapped++;
            info_.maxSnapDistance = max(info_.maxSnapDistance, mag(delta));
        }
    }

    if (info_.nSnapped)
    {
        mesh_.movePoints(newPts);
    }
}


void Foam::cutCellMesher::findCutCells(List<cellCut>& cuts)
{
    //- Every cell needs an inside/outside label, including the ones the
    //- surface does not cross, so a single closed surface is required here
    if (surfaces_.size() != 1)
    {
        FatalErrorInFunction
            << "Cutting with " << surfaces_.size() << " surfaces is not"
            << " supported in this slice. Cells the surface does not cross still"
            << " need an inside/outside label, so give exactly one geometry"
            << " entry." << exit(FatalError);
    }

    const cutCellSurface& surf = surfaces_[0];

    const label insideZone = cellZoneID(surf.name() + "_inside");
    const label outsideZone = cellZoneID(surf.name() + "_outside");

    //- Inside/outside of the cell centres, for the cells the surface does not
    //- cross. Ray casting against the same triangulation the cutting uses, so
    //- the two cannot disagree.
    pointField centres(mesh_.nCells());

    forAll(centres, celli)
    {
        centres[celli] = mesh_.cellCentres()[celli];
    }

    const boolList cellInside(surf.search().calcInside(centres));

    //- Region zone of every cell, including the ones the surface does not
    //- cross, so that the whole mesh ends up labelled
    cellRegionZone_.setSize(mesh_.nCells());

    forAll(cellInside, celli)
    {
        cellRegionZone_[celli] = (cellInside[celli] ? insideZone : outsideZone);
    }

    DynamicList<cellCut> found(cuts.size());

    label maxPlanes = 0;

    forAll(mesh_.cells(), celli)
    {
        const cutCellPolyhedron cell(cellPolyhedron(mesh_, celli));

        const boundBox bb(cell.points(), false);

        const labelList cand(surf.search().tree().findBox(treeBoundBox(bb.min(), bb.max())));

        if (cand.empty())
        {
            continue;
        }

        List<pointField> tris(cand.size());

        forAll(cand, i)
        {
            tris[i] = surf.triangle(cand[i]);
        }

        cutCellResult result;

        if (celli == dumpCell_)
        {
            //- The one cell, with the triangles crossing it, written so that a
            //- failing cell can be replayed without running the case again.
            dictionary d;
            d.add("tol", tol_);
            d.add("celli", celli);
            d.add("cellPoints", cell.points());
            d.add("cellFaces", cell.faces());
            d.add("tris", tris);

            OFstream os(mesh_.time().path()/"cutCellDump");
            d.write(os, false);
        }

        //- Boundary association by inheritance: each triangle carries the
        //- region of the surface it belongs to, and a face a facet produces
        //- inherits it. This replaces matching geometry to patches afterwards.
        labelList patchOf(cand.size());

        forAll(cand, i)
        {
            //- The region comes from the surface's own record, not from
            //- surface()[i].region(): that is 0 for every face, because the
            //- triSurface built from our triFaces puts every face in region 0.
            //- A format with no regions leaves the record empty, and then every
            //- face is region 0 and the caller makes a single zone.
            patchOf[i] = surf.regionOf().empty()
              ? 0
              : surf.regionOf()[cand[i]];
        }

        const cutState state
        (
            cutMethod_ == "loop"
          ? cutCellByLoop(cell, tris, patchOf, tol_, result, celli == dumpCell_)
          : cutCellByTriangles(cell, tris, patchOf, tol_, result)
        );

        if (state == cutState::untouched)
        {
            if (loopCut_)
            {
                static label nUntouched = 0;

                if (nUntouched < 6)
                {
                    nUntouched++;

                    Info<< "LOOPCELL " << celli << " untouched, polygons "
                        << result.nPolygons << endl;
                }
            }

            if (result.nPolygons > 0)
            {
                //- Triangles cross the cell, but no plane of theirs divides it.
                //- When one of the cell's points was collapsed onto a plane,
                //- that is exactly what this looks like: the plane passes
                //- through the point, so there is nothing left to divide. It is
                //- an expected outcome of the collapse, counted and reported
                //- rather than failing the run. Otherwise it is a case we
                //- cannot handle and is reported as one.
                bool collapsed = false;

                for (const label p : mesh_.cellPoints()[celli])
                {
                    if (mag(snapDelta_[p]) > 0)
                    {
                        collapsed = true;
                        break;
                    }
                }

                if (collapsed)
                {
                    info_.nDegenerate++;
                    info_.degenerate.append(celli);
                    info_.degenerateReason.append
                    (
                        "one of its points was collapsed onto the surface, so"
                        " the plane passes through it and there is nothing left"
                        " to divide"
                    );
                }
                else
                {
                    info_.nUnsupported++;
                    info_.unsupported.append(celli);
                    info_.unsupportedReason.append
                    (
                        "the surface touches it without dividing it"
                    );
                }
            }
            continue;
        }

        if (state == cutState::refused)
        {
            //- Two pieces would have met along a line rather than over an area,
            //- so there is no face to put between them. The cell is left whole.
            //- A grazing surface is a geometric outcome rather than a case we
            //- cannot handle, so it is counted and listed, and the run does not
            //- fail because of it.
            info_.nDegenerate++;
            info_.degenerate.append(celli);

            if (result.refusal == cutRefusal::malformedPiece)
            {
                info_.nRefusedMalformed++;
                info_.degenerateReason.append
                (
                    "its pieces are malformed, so no face of one is a face of"
                    " the other"
                );
            }
            else if (result.refusal == cutRefusal::degenerateSharedFace)
            {
                info_.nRefusedDegenerateFace++;
                info_.degenerateReason.append
                (
                    "the face between its pieces has no area"
                );
            }
            else if (result.refusal == cutRefusal::faceWithHole)
            {
                info_.nRefusedFaceWithHole++;
                info_.degenerateReason.append
                (
                    "a chain closes on itself inside one of its faces, so that"
                    " face's far side is a face with a hole"
                );
            }

            continue;
        }

        if (result.pieces.size() < 2)
        {
            if (loopCut_)
            {
                static label nLoopWhole = 0;

                if (nLoopWhole < 6)
                {
                    nLoopWhole++;

                    Info<< "LOOPCELL " << celli << " left whole ("
                        << result.pieces.size() << " pieces)" << endl;
                }
            }

            info_.nUnsupported++;
            info_.unsupported.append(celli);
            info_.unsupportedReason.append
            (
                "the surface crosses it but produces no pieces"
            );
            continue;
        }

        //- Classify every piece against the surface. The pieces of one cell are
        //- not in general one region each, because a facet can cover part of
        //- the cross-section only.
        pointField pieceCentres(result.pieces.size());

        forAll(result.pieces, pi)
        {
            pieceCentres[pi] = result.pieces[pi].poly.centre();
        }

        const boolList pieceInside(surf.search().calcInside(pieceCentres));

        cellCut c;
        c.celli = celli;
        c.planes = result.planes;
        c.faceRegions = result.faceRegions;
        c.faceRegionPiece = result.faceRegionPiece;
        c.pieces = result.pieces;
        c.sharedFaces = result.sharedFaces;
        c.pieceCell.setSize(c.pieces.size());
        c.pieceZone.setSize(c.pieces.size());

        label nIn = 0;
        label nOut = 0;

        forAll(c.pieces, pi)
        {
            c.pieceZone[pi] = (pieceInside[pi] ? insideZone : outsideZone);

            if (pieceInside[pi])
            {
                nIn++;
            }
            else
            {
                nOut++;
            }
        }

        //- Whether the surface's own planes should decide this cell, and there
        //- are two reasons to ask, both of them evidence rather than
        //- preference.
        //-
        //- The first is that the centre's verdict left every piece in one
        //- region, so the cell is not separated at all. Measured on
        //- 10_box_regions: 17 of its 698 crossed cells were left whole this way
        //- and carried 5.3% of the interface area of an exact planar solid away
        //- with them, silently.
        //-
        //- The second is that it left a **shared face with both of its pieces in
        //- one region**. A shared face lies on the surface, so its two pieces are
        //- on opposite sides of it and cannot be one region: where they are, the
        //- verdict is demonstrably wrong for this cell. Measured on the same
        //- case, cell 2945: pieces 0 and 1 were both put in the inside zone
        //- while the surface separates them, piece 1 standing 0.027 - half a
        //- cell - in front of the first plane's facet, and the two faces between
        //- them were then dropped as internal to a merged cell, taking 0.0031 of
        //- interface area with them.
        //-
        //- The cause in both cases is the same: a piece's mean centre is not a
        //- sound witness, because the loop builds a piece by joining the
        //- polygons of the chains, so a piece can be non-convex and its mean can
        //- lie on the other side of the surface from the piece itself. On the 17
        //- cells above, the solid's own geometry agrees with the query: both
        //- pieces' means really are inside the box.
        //-
        //- So the pieces are asked the surface's own question instead. The
        //- facets crossing the cell give outward-oriented planes, and a piece is
        //- inside when every one of its own points is behind every one of those
        //- planes. That is the rule the arrangement has always used to order its
        //- pieces - the sign of the largest distance over a piece's points
        //- rather than the sign at its centre, which a sliver makes meaningless -
        //- and it is asked *only* where the evidence above says the centre has
        //- failed, so a cell that separates cleanly cannot be changed by it.
        bool decideByPlanes = (nIn == 0 || nOut == 0);
        bool byDisagreement = false;

        if (!decideByPlanes)
        {
            forAll(c.sharedFaces, si)
            {
                const cutCellPieceFace& shared = c.sharedFaces[si];

                if
                (
                    shared.pieceA < 0 || shared.pieceA >= c.pieceZone.size()
                 || shared.pieceB < 0 || shared.pieceB >= c.pieceZone.size()
                )
                {
                    continue;
                }

                if
                (
                    shared.patch >= 0
                 && c.pieceZone[shared.pieceA] == c.pieceZone[shared.pieceB]
                )
                {
                    decideByPlanes = true;
                    byDisagreement = true;
                    break;
                }
            }
        }

        if (decideByPlanes)
        {
            //- The margin the pieces were built with: cutCellLoop merges a
            //- piece's points within ten times the caller's tolerance
            const scalar surfaceTol = 10*tol_;

            labelList fromSurface(c.pieces.size(), -1);
            label nInSurface = 0;
            label nOutSurface = 0;

            forAll(c.pieces, pi)
            {
                bool behindAll = !result.planes.empty();

                forAll(result.planes, pli)
                {
                    const vector& n = result.planes[pli].normal();
                    const point& p0 = result.planes[pli].refPoint();

                    scalar deepest = -GREAT;

                    forAll(c.pieces[pi].poly.points(), qi)
                    {
                        deepest = max
                        (
                            deepest,
                            (n & (c.pieces[pi].poly.points()[qi] - p0))
                        );
                    }

                    if (deepest > surfaceTol)
                    {
                        behindAll = false;
                        break;
                    }
                }

                fromSurface[pi] = (behindAll ? 1 : 0);

                if (behindAll)
                {
                    nInSurface++;
                }
                else
                {
                    nOutSurface++;
                }
            }

            if (nInSurface > 0 && nOutSurface > 0)
            {
                nIn = nInSurface;
                nOut = nOutSurface;

                forAll(c.pieces, pi)
                {
                    c.pieceZone[pi] = (fromSurface[pi] == 1)
                      ? insideZone
                      : outsideZone;
                }

                if (byDisagreement)
                {
                    info_.nZoneRepaired++;
                }
            }
            else if (byDisagreement)
            {
                //- The planes could not separate it either, so the centre's
                //- verdict stands and the faces this cell shares stay internal.
                //- Counted and reported rather than dropped in silence.
                info_.nZoneUnresolved++;
            }
            else
            {
                //- The surface passes through the cell without separating it, so
                //- every piece is in one region. Counted, since the pieces are
                //- still valid cells.
                info_.nSingleRegion++;
            }
        }

        //- The interface is the shared faces whose pieces are in different
        //- regions, so its total area is known before the mesh exists
        forAll(c.sharedFaces, si)
        {
            const cutCellPieceFace& shared = c.sharedFaces[si];

            if (c.pieceZone[shared.pieceA] != c.pieceZone[shared.pieceB])
            {
                const cutCellPolyhedron& piece = c.pieces[shared.pieceA].poly;
                const face& f = piece.faces()[shared.faceInA];

                pointField poly(f.size());

                forAll(f, i)
                {
                    poly[i] = piece.points()[f[i]];
                }

                info_.interfaceArea += polygonArea(poly);
            }
        }

        //- The pieces of one region that touch each other are one cell, so a
        //- cell the surface crosses twice is two cells, not four
        c.pieceGroup = groupPieces(c.pieces.size(), c.sharedFaces, c.pieceZone);

        const label nCellGroups = max(c.pieceGroup) + 1;

        info_.nCut++;

        info_.maxPieces = max(info_.maxPieces, result.pieces.size());

        if (loopCut_)
        {
            static label nLoopReported = 0;

            if (nLoopReported < 6)
            {
                nLoopReported++;

                Info<< "LOOPCELL " << celli << " pieces "
                    << result.pieces.size() << endl;
            }
        }

        if (result.pieces.size() > 2)
        {
            info_.nManyPieces++;

            if (info_.manyPieces.size() < 10)
            {
                info_.manyPieces.append(celli);
            }
        }


        info_.nPieces += c.pieces.size();
        info_.nGroups += nCellGroups;
        info_.nMerged += c.pieces.size() - nCellGroups;

        maxPlanes = max(maxPlanes, c.planes.size());

        found.append(c);
    }

    //- Cut cells by the number of planes crossing them
    info_.byPlanes.setSize(maxPlanes + 1, Zero);

    forAll(found, ci)
    {
        info_.byPlanes[found[ci].planes.size()]++;
    }

    info_.nCandidate = info_.nCut + info_.nUnsupported + info_.nDegenerate;

    cuts.clear();

    forAll(found, ci)
    {
        cuts.append(found[ci]);
    }
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::cutCellMesher::assemble
(
    const List<cellCut>& cuts
)
{
    const pointField& pts = mesh_.points();
    const faceList& fcs = mesh_.faces();
    const cellList& cells = mesh_.cells();
    const label nPoints = mesh_.nPoints();

    //- The zones must exist before the change: polyTopoChange validates the
    //- zoneIDs it is given against the mesh
    const label insideZone = cellZoneID(surfaces_[0].name() + "_inside");

    //- One faceZone per surface region when the surface has more than one, so a
    //- face is filed under the region it came from and a case can name its
    //- patches from that. The zone takes the region's own name, which is what
    //- snappy's faceZoneNaming REGION does: faceZoneNames_ = regionNames. A
    //- single-region surface keeps the original name, so the cases that predate
    //- this cannot move.
    const word& surfName = surfaces_[0].name();
    const wordList& regionNames = surfaces_[0].regionNames();

    //- The region count comes from the faces, not from surface().patches():
    //- that returns only the patches the file declares, and a file carrying its
    //- regions as tokens declares none, so reading it there found one patch for
    //- a three-region surface and the per-region path never ran.
    label maxRegion = -1;

    forAll(cuts, ci)
    {
        forAll(cuts[ci].sharedFaces, si)
        {
            maxRegion = max(maxRegion, cuts[ci].sharedFaces[si].patch);
        }
    }

    const bool multiRegion = (maxRegion > 0) && (maxRegion < regionNames.size());


    Info<< "  interface zones: maxRegion " << maxRegion
        << ", names " << regionNames.size()
        << ", multiRegion " << multiRegion << endl;

    labelList interfaceZoneOf(maxRegion + 1, -1);

    forAll(interfaceZoneOf, ri)
    {
        interfaceZoneOf[ri] = faceZoneID
        (
            multiRegion
          ? regionNames[ri]
          : word(surfName + "_interface")
        );
    }

    //- The fallback for a face whose region is unknown.
    const label interfaceZone = (multiRegion ? -1 : interfaceZoneOf[0]);

    //- All points, faces and cells of the mesh are taken as they are, and the
    //- cut is expressed as modifications on top of them. That is what keeps the
    //- labels of untouched cells, and so the parent map, simple.
    polyTopoChange topoChange(mesh_);

    //- Where a cut cell is, by background cell label
    Map<label> cutOf(2*cuts.size());

    forAll(cuts, ci)
    {
        cutOf.insert(cuts[ci].celli, ci);
    }


    //- New points, keyed by the segment they lie on, so that every cell meeting
    //- a position refers to one point. One chain per segment, since two planes
    //- crossing one segment are two points.
    //- One point per location. Keyed on the coordinate rather than on the segment
    //- between two points, because the same location is reached under different
    //- point labels depending on which route asks for it: an edge crossing
    //- through newPoint, the same crossing from a face split, or an interior
    //- point of a piece that touches its neighbour along a line. Only
    //- coincidences matter, so the key is the coordinate rounded to tol_, each
    //- bucket carries a chain, and the look up scans the 26 neighbouring buckets
    //- as well so that a pair whose coordinates straddle a bucket boundary is
    //- still found.
    //- How close two points must be to count as the same point. Derived from
    //- tol_ rather than from the mesh bounds, because tol_ is decided once in
    //- the constructor and guarded against being zero, while reading the bounds
    //- here gave zero and trapped. tol_ is the geometric scale, 1e-9 of the mesh
    //- diagonal, and identity is a thousandth of that, 1e-12 of the diagonal.
    //- It has to be far tighter than tol_: the two routes that reach one
    //- crossing differ by round off, around 1e-16 here, and merging at tol_
    //- would collapse real features a nanometre across, which showed up as a
    //- change in the total volume.
    const scalar pointTol(1e-3*tol_);

    Map<label> firstPoint(2*cuts.size());
    DynamicList<point> newCoord(2*cuts.size());
    DynamicList<label> newLabel(2*cuts.size());
    DynamicList<label> newNext(2*cuts.size());

    auto coordOf = [&](const label p) -> point
    {
        return (p < nPoints ? pts[p] : newCoord[p - nPoints]);
    };

    //- A hash of the coordinate quantised at pointTol. Quantising at that scale
    //- on a mesh a unit across gives indices around 1e11, which do not fit in a
    //- 32 bit label, so the hash is folded on the scalar side and only the
    //- result is a label. A collision costs a comparison, never a wrong answer,
    //- because the look up tests the coordinate itself.
    auto bucket = [&](const point& x) -> label
    {
        const scalar h
        (
            std::fmod
            (
                std::floor(x.x()/pointTol)*73856093.0
              + std::floor(x.y()/pointTol)*19349663.0
              + std::floor(x.z()/pointTol)*83492791.0,
                2147483647.0
            )
        );

        return label(h);
    };

    //- The mesh point in each bucket. The cache above holds only the points this
    //- assembly made, so without this a crossing that lands exactly on a mesh
    //- point would be added as a second point at the same coordinates. Measured
    //- on a 320-triangle sphere in a 16-cubed background: 18 duplicate points,
    //- every one of them at the sphere's pole, where five facets meet and their
    //- crossings land on existing grid points. It is one defect and not two,
    //- because a duplicated point leaves the faces either side of it not sharing
    //- it, which is what the 860 edges without two faces in that mesh are.
    //-
    //- One point per bucket is enough: two distinct mesh points cannot be within
    //- pointTol of each other in a mesh that is not already degenerate.
    Map<label> meshPointInBucket(pts.size());

    forAll(pts, pointi)
    {
        meshPointInBucket.set(bucket(pts[pointi]), pointi);
    }

    //- The label of the point already at x, or -1
    auto findPoint = [&](const point& x) -> label
    {
        for (label di = -1; di <= 1; di++)
        {
            for (label dj = -1; dj <= 1; dj++)
            {
                for (label dk = -1; dk <= 1; dk++)
                {
                    const point probe
                    (
                        x + point
                        (
                            pointTol*scalar(di),
                            pointTol*scalar(dj),
                            pointTol*scalar(dk)
                        )
                    );

                    Map<label>::const_iterator iter =
                        firstPoint.find(bucket(probe));

                    if (iter == firstPoint.end())
                    {
                        continue;
                    }

                    for (label n = iter(); n >= 0; n = newNext[n])
                    {
                        if (mag(newCoord[n] - x) <= pointTol)
                        {
                            return newLabel[n];
                        }
                    }
                }
            }
        }

        //- And the mesh's own points, which the cache above does not hold
        for (label di = -1; di <= 1; di++)
        {
            for (label dj = -1; dj <= 1; dj++)
            {
                for (label dk = -1; dk <= 1; dk++)
                {
                    const point probe
                    (
                        x + point
                        (
                            pointTol*scalar(di),
                            pointTol*scalar(dj),
                            pointTol*scalar(dk)
                        )
                    );

                    Map<label>::const_iterator iter =
                        meshPointInBucket.find(bucket(probe));

                    if (iter == meshPointInBucket.end())
                    {
                        continue;
                    }

                    if (mag(pts[iter()] - x) <= pointTol)
                    {
                        return iter();
                    }
                }
            }
        }

        return -1;
    };

    //- The point at a location, made once however many routes ask for it
    auto pointAt = [&](const point& x) -> label
    {
        const label found = findPoint(x);

        if (found >= 0)
        {
            return found;
        }

        const label ptLabel = topoChange.addPoint(x, -1, -1, true);

        const label n = newCoord.size();

        newCoord.append(x);
        newLabel.append(ptLabel);
        newNext.append(-1);

        const label key = bucket(x);

        Map<label>::iterator iter = firstPoint.find(key);

        if (iter == firstPoint.end())
        {
            firstPoint.insert(key, n);
        }
        else
        {
            newNext[n] = iter();
            iter() = n;
        }

        info_.nNewPoints++;

        return ptLabel;
    };

    auto newPoint = [&]
    (
        const label a,
        const label b,
        const vector& n,
        const point& refPoint
    ) -> label
    {
        const point pa(coordOf(a));
        const point pb(coordOf(b));

        const scalar da = planeDist(n, refPoint, pa);
        const scalar db = planeDist(n, refPoint, pb);

        const scalar w = da/(da - db);

        return pointAt(pa + w*(pb - pa));
    };


    //- Split a face polygon by a plane, in mesh point labels. False when the
    //- face is not split, that is when all of its points lie on one side of the
    //- plane or on the plane itself.
    auto splitFace = [&]
    (
        const face& f,
        const vector& n,
        const point& refPoint,
        face& minusPoly,
        face& plusPoly
    ) -> bool
    {
        DynamicList<label> m(f.size() + 2);
        DynamicList<label> p(f.size() + 2);

        bool hasMinus = false;
        bool hasPlus = false;

        forAll(f, i)
        {
            const label a = f[i];
            const label b = f[f.fcIndex(i)];

            const scalar da = planeDist(n, refPoint, coordOf(a));
            const scalar db = planeDist(n, refPoint, coordOf(b));

            if (da <= tol_)
            {
                m.append(a);
            }

            if (da >= -tol_)
            {
                p.append(a);
            }

            if (da < -tol_)
            {
                hasMinus = true;
            }

            if (da > tol_)
            {
                hasPlus = true;
            }

            if
            (
                (da > tol_ && db < -tol_)
             || (da < -tol_ && db > tol_)
            )
            {
                const label x = newPoint(a, b, n, refPoint);

                m.append(x);
                p.append(x);
            }
        }

        if (!hasMinus || !hasPlus)
        {
            return false;
        }

        //- A plane that grazes the face within the tolerance keeps a point on
        //- both sides, so one of the two polygons is a strip whose vertices
        //- coincide and which has no area. A face with no area is not a face of
        //- a mesh and leaves the two cells it borders open, so the split is
        //- refused and the face stays whole with respect to this plane. The
        //- subdivision is computed once per face from the union of the planes
        //- of every cell that meets it, so both cells still see the same
        //- sub-faces and the mesh stays conformal.
        pointField mCoords(m.size());
        pointField pCoords(p.size());

        forAll(m, i)
        {
            mCoords[i] = coordOf(m[i]);
        }

        forAll(p, i)
        {
            pCoords[i] = coordOf(p[i]);
        }

        if
        (
            polygonDegenerate(mCoords, tol_)
         || polygonDegenerate(pCoords, tol_)
        )
        {
            return false;
        }

        minusPoly.setSize(m.size());

        forAll(m, i)
        {
            minusPoly[i] = m[i];
        }

        plusPoly.setSize(p.size());

        forAll(p, i)
        {
            plusPoly[i] = p[i];
        }

        return true;
    };


    //- Subdivide the faces the surface crosses, once each, by the planes of
    //- every cut cell that meets them, so that all the cells on either side see
    //- the same sub-faces. This is what keeps the mesh conformal when two cells
    //- sharing a face are cut by different planes.
    struct faceSplit
    {
        //- The sub-faces, in mesh point labels
        List<face> subFaces;

        //- The points on the face and where they are
        List<label> points;
        List<point> coords;
    };

    Map<label> splitOf(2*cuts.size());
    DynamicList<faceSplit> splits(2*cuts.size());

    {
        //- Which cuts meet each face
        Map<DynamicList<label>> faceCuts(2*cuts.size());

        forAll(cuts, ci)
        {
            for (const label facei : cells[cuts[ci].celli])
            {
                if (!faceCuts.found(facei))
                {
                    faceCuts.insert(facei, DynamicList<label>(2));
                }

                faceCuts[facei].append(ci);
            }
        }

        forAllConstIters(faceCuts, iter)
        {
            const label facei = iter.key();

            //- The union of their planes, since every cell on the face has to
            //- see a face split by all of them
            DynamicList<plane> facePlanes(4);

            for (const label ci : iter())
            {
                for (const plane& pln : cuts[ci].planes)
                {
                    bool known = false;

                    forAll(facePlanes, pi)
                    {
                        if
                        (
                            samePlane
                            (
                                pln.normal(), pln.refPoint(),
                                facePlanes[pi].normal(),
                                facePlanes[pi].refPoint(),
                                tol_
                            )
                        )
                        {
                            known = true;
                            break;
                        }
                    }

                    if (!known)
                    {
                        facePlanes.append(pln);
                    }
                }
            }

            std::sort(facePlanes.begin(), facePlanes.end(), planeLess);

            //- Clip the face by each plane in turn
            DynamicList<face> frontier(2);
            frontier.append(fcs[facei]);

            forAll(facePlanes, pi)
            {
                DynamicList<face> next(frontier.size() + 1);

                forAll(frontier, fi)
                {
                    face m, p;

                    if
                    (
                        splitFace
                        (
                            frontier[fi],
                            facePlanes[pi].normal(),
                            facePlanes[pi].refPoint(),
                            m, p
                        )
                    )
                    {
                        next.append(m);
                        next.append(p);
                    }
                    else
                    {
                        next.append(frontier[fi]);
                    }
                }

                frontier.transfer(next);
            }

            if (frontier.size() < 2)
            {
                continue;
            }

            faceSplit split;
            split.subFaces.setSize(frontier.size());

            forAll(frontier, fi)
            {
                split.subFaces[fi] = frontier[fi];
            }

            //- Every face meeting an edge has to carry the same points on it. A
            //- face is split only by the planes of the cut cells on its own side,
            //- so a crossing made by a neighbouring cell's plane on a shared edge
            //- is missing here, and the two faces' sub-polygons then leave a
            //- sliver of area attributed to the wrong pair of cells. The point is
            //- inserted rather than the split repeated: it lies on the edge, so
            //- no shape changes, and its label comes from newPoint, the same
            //- cache the split uses, so the face on the other side of the edge
            //- receives exactly the same one.
            for (const label edgei : mesh_.faceEdges()[facei])
            {
                const edge& e = mesh_.edges()[edgei];
                const label p0 = e[0];
                const label p1 = e[1];
                const point q0(coordOf(p0));
                const point q1(coordOf(p1));

                const vector d(q1 - q0);
                const scalar d2(d & d);

                if (d2 <= VSMALL)
                {
                    continue;
                }

                //- Where a point sits along this edge, and whether it is on it
                auto param = [&](const point& x) -> scalar
                {
                    return ((x - q0) & d)/d2;
                };

                auto onEdge = [&](const point& x) -> bool
                {
                    const scalar t = param(x);

                    if (t < -SMALL || t > 1 + SMALL)
                    {
                        return false;
                    }

                    return mag(x - (q0 + t*d)) <= pointTol;
                };

                //- The planes of every cut cell containing this edge
                DynamicList<plane> edgePlanes(4);

                for (const label celli : mesh_.edgeCells()[edgei])
                {
                    if (!cutOf.found(celli))
                    {
                        continue;
                    }

                    for (const plane& pln : cuts[cutOf[celli]].planes)
                    {
                        bool known = false;

                        forAll(edgePlanes, pi)
                        {
                            if
                            (
                                samePlane
                                (
                                    pln.normal(), pln.refPoint(),
                                    edgePlanes[pi].normal(),
                                    edgePlanes[pi].refPoint(),
                                    tol_
                                )
                            )
                            {
                                known = true;
                                break;
                            }
                        }

                        if (!known)
                        {
                            edgePlanes.append(pln);
                        }
                    }
                }

                if (edgePlanes.empty())
                {
                    continue;
                }

                //- Their crossings, with the parameter along the edge
                DynamicList<label> crossLabel(edgePlanes.size());
                DynamicList<scalar> crossT(edgePlanes.size());

                forAll(edgePlanes, pi)
                {
                    const plane& pln = edgePlanes[pi];

                    const scalar d0 = planeDist(pln.normal(), pln.refPoint(), q0);
                    const scalar d1 = planeDist(pln.normal(), pln.refPoint(), q1);

                    if (mag(d0 - d1) <= VSMALL)
                    {
                        continue;
                    }

                    const scalar t = d0/(d0 - d1);

                    if (t <= 0 || t >= 1)
                    {
                        continue;
                    }

                    crossLabel.append
                    (
                        newPoint(p0, p1, pln.normal(), pln.refPoint())
                    );

                    crossT.append(t);
                }

                if (crossLabel.empty())
                {
                    continue;
                }

                //- Few entries, so a straight insertion sort
                for (label i = 1; i < crossLabel.size(); i++)
                {
                    for (label j = i; j > 0 && crossT[j] < crossT[j - 1]; j--)
                    {
                        std::swap(crossT[j], crossT[j - 1]);
                        std::swap(crossLabel[j], crossLabel[j - 1]);
                    }
                }

                //- Insert them into every sub-face edge that lies along this one
                forAll(split.subFaces, sfi)
                {
                    const face& sub = split.subFaces[sfi];

                    DynamicList<label> grown(sub.size() + crossLabel.size());

                    forAll(sub, k)
                    {
                        const label a = sub[k];
                        const label b = sub[sub.fcIndex(k)];

                        grown.append(a);

                        const point pa(coordOf(a));
                        const point pb(coordOf(b));

                        if (!onEdge(pa) || !onEdge(pb))
                        {
                            continue;
                        }

                        const scalar ta = param(pa);
                        const scalar tb = param(pb);

                        const scalar lo = min(ta, tb);
                        const scalar hi = max(ta, tb);

                        forAll(crossT, ci)
                        {
                            if (crossT[ci] <= lo || crossT[ci] >= hi)
                            {
                                continue;
                            }

                            if (crossLabel[ci] != a && crossLabel[ci] != b)
                            {
                                grown.append(crossLabel[ci]);
                            }
                        }
                    }

                    if (grown.size() != sub.size())
                    {
                        face newSub(grown.size());

                        forAll(grown, i)
                        {
                            newSub[i] = grown[i];
                        }

                        split.subFaces[sfi] = newSub;
                    }
                }
            }

            //- The points on the face, so that a point of a piece that lies on
            //- this face resolves to the same label
            forAll(split.subFaces, fi)
            {
                forAll(split.subFaces[fi], i)
                {
                    const label p = split.subFaces[fi][i];

                    bool known = false;

                    forAll(split.points, j)
                    {
                        if (split.points[j] == p)
                        {
                            known = true;
                            break;
                        }
                    }

                    if (!known)
                    {
                        split.points.append(p);
                        split.coords.append(coordOf(p));
                    }
                }
            }

            splitOf.insert(facei, splits.size());
            splits.append(split);
        }

        //- A face the surface does not cross can still share a line with one that
        //- is split. Step 3 visits only the faces in faceCuts - those bordering a
        //- cut cell - so a face between two uncut cells is never visited, keeps an
        //- edge whole, and disagrees with the neighbour that was split. Grow the
        //- set until it is closed, adding each newly reached face the points
        //- already computed for the line it shares.
        //
        //- The match cannot be by edge identity: a split neighbour's edges are the
        //- HALVES of the shared line, so a whole edge spanning the same line shares
        //- no edge with it. Match by whether the edge CONTAINS a split point.
        label nGrown = 0;
        label nSweep = 0;

        for (label sweep = 0; sweep < 8; sweep++)
        {
            //- Every point a split face puts on one of its own edges, keyed by
            //- that mesh edge. A split point lies on a mesh edge, so the only two
            //- faces that can take it are the two faces of that edge; keeping the
            //- points by edge lets a candidate face look at the few points on its
            //- own edge instead of at every split point the cut has made.
            //-
            //- Visited in the same order as the flat list it replaces, and a point
            //- is appended to its edge only the first time it is met, so each
            //- edge's list is that flat list restricted to the edge: the points a
            //- face takes, and the order it takes them, are unchanged.
            Map<DynamicList<Tuple2<point, label>>> edgeSplitPts;

            forAllConstIters(splitOf, iter)
            {
                const faceSplit& sp = splits[iter()];
                const face& sf = mesh_.faces()[iter.key()];
                const labelList& sfEdges = mesh_.faceEdges()[iter.key()];

                //- Every point the split carries. A split point is an ENDPOINT of
                //- the sub-face edges, not interior to one, so it has to be taken
                //- from the point list rather than found by an edge test.
                forAll(sp.points, k)
                {
                    const point& x = sp.coords[k];
                    const label p = sp.points[k];

                    forAll(sf, i)
                    {
                        const point& pa = pts[sf[i]];
                        const point& pb = pts[sf[sf.fcIndex(i)]];
                        const vector ab(pb - pa);
                        const scalar ab2 = ab & ab;

                        if (ab2 <= SMALL)
                        {
                            continue;
                        }

                        if (mag(x - pa) <= tol_ || mag(x - pb) <= tol_)
                        {
                            continue;
                        }

                        const scalar t = ((x - pa) & ab)/ab2;

                        if (t <= SMALL || t >= 1 - SMALL)
                        {
                            continue;
                        }

                        if (mag(x - (pa + t*ab)) > tol_)
                        {
                            continue;
                        }

                        DynamicList<Tuple2<point, label>>& onEdge =
                            edgeSplitPts(sfEdges[i]);

                        bool known = false;

                        forAll(onEdge, q)
                        {
                            if (onEdge[q].second() == p)
                            {
                                known = true;
                                break;
                            }
                        }

                        if (!known)
                        {
                            onEdge.append(Tuple2<point, label>(x, p));
                        }

                        break;
                    }
                }
            }

            if (!edgeSplitPts.size())
            {
                break;
            }

            label addedThisSweep = 0;

            //- Only a face that shares a mesh edge with an already-split face can
            //- carry one of these points: a split point lies on a mesh edge, so
            //- the face either side of that edge is the only place it can be
            //- taken. Scanning every face of the mesh tested all 13056 of them
            //- against every split point and was 244 s of a 258 s run at 25.5
            //- facets per cell. The grown set is unchanged - a face with no split
            //- point strictly inside one of its edges is never taken - and the
            //- candidates are still visited in label order, so the mesh the
            //- assembly builds is the same one.
            boolList candidate(mesh_.nFaces(), false);

            forAllConstIters(splitOf, splIter)
            {
                for (const label edgei : mesh_.faceEdges()[splIter.key()])
                {
                    for (const label nei : mesh_.edgeFaces()[edgei])
                    {
                        candidate[nei] = true;
                    }
                }
            }

            forAll(candidate, facei)
            {
                if (!candidate[facei] || splitOf.found(facei))
                {
                    continue;
                }

                const face& f = mesh_.faces()[facei];
                const labelList& fEdges = mesh_.faceEdges()[facei];

                DynamicList<label> take(4);

                forAll(f, i)
                {
                    //- Only the split points on this face's own edge, which is the
                    //- few the edge index holds
                    if (!edgeSplitPts.found(fEdges[i]))
                    {
                        continue;
                    }

                    const point pa(pts[f[i]]);
                    const point pb(pts[f[f.fcIndex(i)]]);
                    const vector ab(pb - pa);
                    const scalar ab2 = ab & ab;

                    if (ab2 <= SMALL)
                    {
                        continue;
                    }

                    const DynamicList<Tuple2<point, label>>& onEdge =
                        edgeSplitPts[fEdges[i]];

                    forAll(onEdge, q)
                    {
                        const point& x = onEdge[q].first();

                        if (mag(x - pa) <= tol_ || mag(x - pb) <= tol_)
                        {
                            continue;
                        }

                        const scalar t = ((x - pa) & ab)/ab2;

                        if (t <= SMALL || t >= 1 - SMALL)
                        {
                            continue;
                        }

                        if (mag(x - (pa + t*ab)) > tol_)
                        {
                            continue;
                        }

                        if (findIndex(take, onEdge[q].second()) < 0)
                        {
                            take.append(onEdge[q].second());
                        }
                    }
                }

                if (!take.size())
                {
                    continue;
                }

                //- The face with the taken points inserted along its edges
                DynamicList<label> expanded(2*f.size());

                forAll(f, i)
                {
                    const label first = f[i];
                    const label second = f[f.fcIndex(i)];

                    const point pa(pts[first]);
                    const point pb(pts[second]);
                    const vector ab(pb - pa);
                    const scalar ab2 = ab & ab;

                    DynamicList<Tuple2<scalar, label>> onEdge(4);

                    if (ab2 > SMALL)
                    {
                        forAll(take, k)
                        {
                            const label id = take[k];

                            if (id == first || id == second)
                            {
                                continue;
                            }

                            const point x(coordOf(id));
                            const scalar t = ((x - pa) & ab)/ab2;

                            if (t <= SMALL || t >= 1 - SMALL)
                            {
                                continue;
                            }

                            if (mag(x - (pa + t*ab)) <= tol_)
                            {
                                onEdge.append(Tuple2<scalar, label>(t, id));
                            }
                        }
                    }

                    std::sort
                    (
                        onEdge.begin(),
                        onEdge.end(),
                        [](const Tuple2<scalar, label>& x,
                           const Tuple2<scalar, label>& y)
                        {
                            return x.first() < y.first();
                        }
                    );

                    expanded.append(first);

                    forAll(onEdge, k)
                    {
                        expanded.append(onEdge[k].second());
                    }
                }

                faceSplit grown;
                grown.subFaces.setSize(1);
                grown.subFaces[0].transfer(expanded);

                grown.points.setSize(take.size());
                grown.coords.setSize(take.size());

                forAll(take, k)
                {
                    grown.points[k] = take[k];
                    grown.coords[k] = coordOf(take[k]);
                }

                splitOf.insert(facei, splits.size());
                splits.append(grown);

                addedThisSweep++;
                nGrown++;
            }

            nSweep = sweep + 1;

            if (!addedThisSweep)
            {
                break;
            }
        }

        Info<< "  faces grown to stay conformal: " << nGrown
            << " over " << nSweep << " sweep(s)" << endl;
    }


    //- A point of a piece, as a mesh point label. A point of a cut lies either
    //- on an edge of the cell, which the cells sharing that edge all use, or on
    //- a face of it, which the face's subdivision has already made, or inside
    //- it, where it belongs to this cell alone.

    auto pointOf = [&](const cellCut& c, const point& x) -> label
    {
        const labelList& cPts = mesh_.cellPoints()[c.celli];

        for (const label p : cPts)
        {
            if (mag(pts[p] - x) <= tol_)
            {
                return p;
            }
        }

        for (const label facei : cells[c.celli])
        {
            const face& f = fcs[facei];

            forAll(f, i)
            {
                const label a = f[i];
                const label b = f[f.fcIndex(i)];

                const point pa(pts[a]);
                const point pb(pts[b]);

                const vector ab(pb - pa);
                const scalar t = ((x - pa) & ab)/(ab & ab);

                if (t < -SMALL || t > 1 + SMALL)
                {
                    continue;
                }

                if (mag(x - (pa + t*ab)) <= tol_)
                {
                    //- Which of the cell's planes crosses the edge where this
                    //- point is, if any
                    forAll(c.planes, pi)
                    {
                        const vector& n = c.planes[pi].normal();
                        const point& refPoint = c.planes[pi].refPoint();

                        const scalar da = planeDist(n, refPoint, pa);
                        const scalar db = planeDist(n, refPoint, pb);

                        if
                        (
                            (da > tol_ && db < -tol_)
                         || (da < -tol_ && db > tol_)
                        )
                        {
                            const scalar w = da/(da - db);

                            if (mag((pa + w*(pb - pa)) - x) <= tol_)
                            {
                                return newPoint(a, b, n, refPoint);
                            }
                        }
                    }
                }
            }
        }

        //- Whether any of the cell's faces claims the point by its plane. A point
        //- that a face claims is that face's business, even when the face's
        //- subdivision has nothing at the same place: the next route must not
        //- take it over. Measured, taking it over gave the loop's `09` mesh 96
        //- fewer points and a new `invalidVertexLabels`, on faces no collapse
        //- had bent.
        bool faceClaimed = false;

        for (const label facei : cells[c.celli])
        {
            const scalar magArea = mag(mesh_.faceAreas()[facei]);

            if (magArea <= VSMALL)
            {
                continue;
            }

            const vector n(mesh_.faceAreas()[facei]/magArea);
            const point refPoint(mesh_.faceCentres()[facei]);

            const scalar d = mag(planeDist(n, refPoint, x));

            if (d <= tol_)
            {
                faceClaimed = true;

                if (splitOf.found(facei))
                {
                    const faceSplit& split = splits[splitOf[facei]];

                    forAll(split.points, j)
                    {
                        if (mag(split.coords[j] - x) <= tol_)
                        {
                            return split.points[j];
                        }
                    }
                }
            }
        }

        //- A face the collapse has bent is not a plane, so a piece's vertex on it
        //- can miss the plane test above by far more than any identity tolerance
        //- - 4.55e-4 against 1.7e-9 on 06_oblique_box. The vertex is still that
        //- face's: the cutter's successive cuts and the face's subdivision draw
        //- different chords across the same bent face, so the piece carries one
        //- crossing where the subdivision has already made the other, and they
        //- are 1.34e-5 apart (07 section 36.10). Both are found by a nearest
        //- choice, and each choice is guarded by the margin to its runner-up
        //- rather than by a distance, so no constant enters: measured there the
        //- nearest face wins by 6.7 over the next and its nearest subdivision
        //- point by about 300 over the next, against the 2 asked of both.
        //- Reached only when **no** face claims the point, which is where a bent
        //- face puts a piece's vertex: its plane test then fails for every face of
        //- the cell. The guard is what keeps a face that does claim its point from
        //- being taken over - without it the loop's `09` mesh lost 96 points and
        //- gained `invalidVertexLabels`, on faces no collapse had bent.
        if (!faceClaimed)
        {
            scalar bestFaceD = GREAT;
            scalar nextFaceD = GREAT;
            label bestFace = -1;

            for (const label facei : cells[c.celli])
            {
                const scalar magArea = mag(mesh_.faceAreas()[facei]);

                if (magArea <= VSMALL || !splitOf.found(facei))
                {
                    continue;
                }

                const scalar d =
                    mag
                    (
                        planeDist
                        (
                            mesh_.faceAreas()[facei]/magArea,
                            mesh_.faceCentres()[facei],
                            x
                        )
                    );

                if (d < bestFaceD)
                {
                    nextFaceD = bestFaceD;
                    bestFaceD = d;
                    bestFace = facei;
                }
                else if (d < nextFaceD)
                {
                    nextFaceD = d;
                }
            }

            if (bestFace != -1 && 2*bestFaceD < nextFaceD)
            {
                const faceSplit& split = splits[splitOf[bestFace]];

                scalar bestPtD = GREAT;
                scalar nextPtD = GREAT;
                label bestPt = -1;

                forAll(split.points, j)
                {
                    const scalar d = mag(split.coords[j] - x);

                    if (d < bestPtD)
                    {
                        nextPtD = bestPtD;
                        bestPtD = d;
                        bestPt = split.points[j];
                    }
                    else if (d < nextPtD)
                    {
                        nextPtD = d;
                    }
                }

                if (bestPt != -1 && 2*bestPtD < nextPtD)
                {
                    return bestPt;
                }
            }
        }

        //- Inside the cell, where the point belongs to this cell alone. A cell
        //- has few of these, so a linear scan is enough.
        //- Inside the cell, and shared with any other piece that touches it
        return pointAt(x);
    };


    //- The mesh cell of a cut cell on the side of a point, which is the piece
    //- whose side against every plane of the cell the point agrees with
    auto cellAt = [&](const label celli, const point& x) -> label
    {
        if (!cutOf.found(celli))
        {
            return celli;
        }

        const cellCut& c = cuts[cutOf[celli]];

        //- The loop's pieces are the regions the surface separates, not
        //- combinations of half-spaces, so a sign against each plane does not
        //- identify one: a piece can straddle a plane. Containment does, and
        //- needs no convexity of the piece to be asked.
        if (loopCut_)
        {
            //- What the lookup saw, to be printed if it finds nothing
            scalar bestPlane = GREAT;
            label bestFace = -1;
            label bestRegion = -1;
            scalar bestArea = 0;

            //- The piece of the region containing the query, taken from the
            //- regions of the cell's faces rather than by asking a piece
            //- whether it contains the point. The regions tile each face by
            //- construction, so one of them contains it, and asking a
            //- polyhedron instead failed as soon as the mesher's own plane
            //- splitting and the cutter's regions disagreed about a face:
            //- "No piece of cell 1111 contains (...) although the point is
            //- inside that cell".
            forAll(c.faceRegions, ffi)
            {
                forAll(c.faceRegions[ffi], ri)
                {
                    const pointField& r = c.faceRegions[ffi][ri];

                    if (r.size() < 3)
                    {
                        continue;
                    }

                    const vector rn =
                        normalised((r[1] - r[0]) ^ (r[2] - r[0]));

                    if (mag(rn & (x - r[0])) > 100*tol_)
                    {
                        continue;
                    }

                    point rc = Zero;

                    forAll(r, k)
                    {
                        rc += r[k];
                    }

                    rc /= scalar(r.size());

                    bool in = true;

                    forAll(r, k)
                    {
                        const point& a = r[k];
                        const point& b = r[r.fcIndex(k)];

                        const vector e = (b - a) ^ rn;
                        const vector inward = (((rc - a) & e) >= 0 ? e : -e);

                        if ((inward & (x - a)) < -100*tol_)
                        {
                            in = false;
                            break;
                        }
                    }

                    if (in)
                    {
                        return c.pieceCell[c.faceRegionPiece[ffi][ri]];
                    }

                    const scalar pd = mag(rn & (x - r[0]));

                    if (pd < bestPlane)
                    {
                        bestPlane = pd;
                        bestFace = ffi;
                        bestRegion = ri;
                        bestArea = polygonArea(r);
                    }
                }
            }

            static label nNotFound = 0;

            if (nNotFound < 5)
            {
                nNotFound++;

                Info<< "CELLAT cell " << celli << " at " << x
                    << " no region: faces " << c.faceRegions.size()
                    << " nearest face " << bestFace << " region " << bestRegion
                    << " plane distance " << bestPlane
                    << " area " << bestArea << endl;
            }
        }

        forAll(c.pieces, pi)
        {
            bool match = true;

            forAll(c.planes, k)
            {
                const scalar d = planeDist
                (
                    c.planes[k].normal(), c.planes[k].refPoint(), x
                );

                //- A query within tolerance of a plane lies on it, so it does
                //- not decide which piece it belongs to. A collapsed cut passes
                //- through a mesh point, and a query at that point has to match
                //- the piece on either side rather than a side chosen by the
                //- sign of a distance that is zero.
                if (mag(d) <= tol_)
                {
                    continue;
                }

                if ((d > 0 ? 1 : -1) != c.pieces[pi].side[k])
                {
                    match = false;
                    break;
                }
            }

            if (match)
            {
                return c.pieceCell[pi];
            }
        }

        FatalErrorInFunction
            << "No piece of cell " << c.celli << " contains " << x
            << ", although the point is inside that cell." << exit(FatalError);

        return c.celli;
    };

    //- Set a face of the new mesh, either an existing face or a new one.
    //-
    //- polyMesh keeps internal faces in upper-triangular order, meaning the
    //- lower cell label is the owner. Cutting can reverse the order of a
    //- background face, because a sub-cell added for the positive side is
    //- labelled above every background cell, so reorder and reverse the face
    //- when that happens. The topology change does this itself only when it
    //- renumbers cells, and this cut removes no cell, so it is done here.
    //- Which faces the emission reached. A face of a cut cell never reached is
    //- left with no new position and mapPolyMesh aborts.
    labelHashSet emittedFaces(4*cuts.size());
    labelHashSet removedFaces(4*cuts.size());
    label nDropSize = 0;
    label nDropSame = 0;
    label nNoRegions = 0;

    auto setFace = [&]
    (
        const face& f,
        const label targetFace,     //- existing face to modify, or -1 to add
        label own,
        label nei,
        const label patchi,
        const label zonei,
        bool zoneFlip,
        const char* tag = ""            //- probe: which emission path this is
    )
    {
        emittedFaces.insert(targetFace);

        //- Probe: for one cell, name every new face and the path that made it.
        if (loopCut_ && targetFace < 0 && dumpCell_ >= 0
            && (own == dumpCell_ || nei == dumpCell_))
        {
            Info<< "PROBEFACE " << tag << " own " << own << " nei " << nei
                << " zone " << zonei << " size " << f.size() << " labels";
            forAll(f, i)
            {
                Info<< ' ' << f[i];
            }
            Info<< endl;
        }

        if (targetFace == 2463 || (nei < 0 && patchi < 0))
        {
            static label nOdd = 0;

            if (nOdd < 8)
            {
                nOdd++;

                Info<< "SETFACE target " << targetFace << " size " << f.size()
                    << " own " << own << " nei " << nei
                    << " patch " << patchi << " zone " << zonei
                    << " | isInternal " << (targetFace >= 0
                        ? mesh_.isInternalFace(targetFace) : -1)
                    << " faceNei " << (targetFace >= 0
                        ? mesh_.faceNeighbour()[targetFace] : -1)
                    << " faceOwner " << (targetFace >= 0
                        ? mesh_.faceOwner()[targetFace] : -1)
                    << " nFaces " << mesh_.nFaces()
                    << " nInternal " << mesh_.nInternalFaces() << endl;
            }
        }

        face nf(f);

        bool flipFlux = false;

        if (nei >= 0 && nei < own)
        {
            nf.flip();
            std::swap(own, nei);
            zoneFlip = !zoneFlip;
            flipFlux = true;
        }

        if (targetFace >= 0)
        {
            topoChange.modifyFace
            (
                nf, targetFace, own, nei, flipFlux, patchi, zonei, zoneFlip
            );
        }
        else
        {
            topoChange.addFace
            (
                nf, own, nei, -1, -1, -1, flipFlux, patchi, zonei, zoneFlip
            );
        }
    };


    //- 1. One cell per group of pieces. The pieces of one region that share a
    //- face are one cell, so a cell crossed by two planes is two cells, not
    //- four. The group holding the first piece, which for a single plane is the
    //- one inside the surface, keeps the label of the background cell.
    forAll(cuts, ci)
    {
        cellCut& c = const_cast<cellCut&>(cuts[ci]);

        c.pieceCell.setSize(c.pieces.size());

        //- The mesh cell of each group, and its zone, taken from the first
        //- piece of the group
        labelList groupCell(max(c.pieceGroup) + 1, -1);

        forAll(c.pieces, pi)
        {
            const label g = c.pieceGroup[pi];

            if (groupCell[g] < 0)
            {
                if (g == 0)
                {
                    groupCell[g] = c.celli;

                    topoChange.modifyCell
                    (
                        c.celli, labelList(1, c.pieceZone[pi])
                    );
                }
                else
                {
                    groupCell[g] = topoChange.addCell
                    (
                        -1, -1, -1, c.celli, c.pieceZone[pi]
                    );
                }
            }

            c.pieceCell[pi] = groupCell[g];
        }
    }


    //- 2. Cells the surface does not cross keep their label, and take the
    //- region their centre falls in
    {
        labelHashSet isCut(2*cuts.size());

        forAll(cuts, ci)
        {
            isCut.insert(cuts[ci].celli);
        }

        forAll(cellRegionZone_, celli)
        {
            if (!isCut.found(celli))
            {
                topoChange.modifyCell
                (
                    celli, labelList(1, cellRegionZone_[celli])
                );
            }
        }
    }


    //- 3. The faces of the cells that were cut: every sub-face goes to the
    //- piece of the cell on its own side, and to the piece of the neighbouring
    //- cell on the other side, or to that cell when it was not cut. A face
    //- between two cut cells is reached from both, so it is done once.
    labelHashSet doneFaces(2*cuts.size());

    forAll(cuts, ci)
    {
        const cellCut& c = cuts[ci];

        for (const label facei : cells[c.celli])
        {
            if (doneFaces.found(facei))
            {
                continue;
            }

            doneFaces.insert(facei);

            const label own = mesh_.faceOwner()[facei];

            //- faceNeighbour is only defined for internal faces, so it must not
            //- be read for a boundary face
            const label nei = mesh_.isInternalFace(facei)
              ? mesh_.faceNeighbour()[facei]
              : -1;

            const label patchi = mesh_.boundaryMesh().whichPatch(facei);

            //- Keep whatever faceZone the face was in
            const label zonei = mesh_.faceZones().whichZone(facei);

            bool zoneFlip = false;

            if (zonei >= 0)
            {
                const faceZone& fz = mesh_.faceZones()[zonei];
                zoneFlip = fz.flipMap()[fz.whichFace(facei)];
            }

            //- On the loop path the face's sub-faces ARE the regions the cutter
            //- built, carrying their pieces, so nothing is looked up by a query.
            //- Asking a subdivision the cutter does not know about cannot work:
            //- the mesher's sub-faces come from the facet planes and the
            //- cutter's regions from the loop's chains, so a sub-face crosses a
            //- region boundary and its centre lands outside every region - "no
            //- region: plane distance 0", which is a point in the plane and
            //- outside the polygon.
            if (loopCut_)
            {
                const label cutA = (cutOf.found(own) ? own : -1);
                const label cutB = ((nei >= 0 && cutOf.found(nei)) ? nei : -1);

                if ((cutA >= 0 || cutB >= 0) && nei >= 0)
                {
                    static label nB2 = 0;
                    nB2++;
                    if (dumpCell_ >= 0 && nB2 % 500 == 0)
                    {
                        Info<< "REGB2 loopCut faces reaching block 2: " << nB2 << endl;
                    }

                    const label anchor = (cutA >= 0 ? cutA : cutB);
                    const cellCut& cc = cuts[cutOf[anchor]];

                    label ffi = -1;

                    forAll(cells[anchor], k)
                    {
                        if (cells[anchor][k] == facei)
                        {
                            ffi = k;
                            break;
                        }
                    }

                    const bool hasRegions =
                    (
                        ffi >= 0
                     && ffi < cc.faceRegions.size()
                     && !cc.faceRegions[ffi].empty()
                    );

                    if (!hasRegions)
                    {
                        nNoRegions++;
                    }

                    if (hasRegions)
                    {
                        static label nB3 = 0;
                        static label nB3Other = 0;
                        nB3++;
                        if ((cutA >= 0 ? cutB : cutA) >= 0) nB3Other++;
                        if (dumpCell_ >= 0 && nB3 % 500 == 0)
                        {
                            Info<< "REGB3 region block entered " << nB3
                                << " times, other cut " << nB3Other << endl;
                        }

                        //- The neighbour's regions on this face, to give each
                        //- sub-face its piece there as well
                        const label other = (cutA >= 0 ? cutB : cutA);
                        const cellCut* oc = nullptr;
                        label offi = -1;

                        if (other >= 0)
                        {
                            oc = &cuts[cutOf[other]];

                            forAll(cells[other], k)
                            {
                                if (cells[other][k] == facei)
                                {
                                    offi = k;
                                    break;
                                }
                            }

                            static label nOffi = 0;

                            if (dumpCell_ >= 0 && nOffi < 5)
                            {
                                nOffi++;

                                Info<< "OFFI face " << facei
                                    << " anchor " << anchor
                                    << " other " << other
                                    << " nFacesOfOther " << cells[other].size()
                                    << " offi " << offi
                                    << " faceRegionsOfOther "
                                    << oc->faceRegions.size()
                                    << " regionsOfAnchor " << cc.faceRegions.size() << endl;
                            }
                        }

                        //- Do the two cells agree on this face's regions? They
                        //- must, or the sub-faces of one do not meet those of
                        //- the other and the mesh gets edges without two faces.
                        if (oc && offi >= 0)
                        {
                            static label nCompared = 0;
                            static label nMismatch = 0;

                            scalar areaA = 0;
                            scalar areaB = 0;

                            forAll(cc.faceRegions[ffi], rk)
                            {
                                areaA += polygonArea(cc.faceRegions[ffi][rk]);
                            }

                            forAll(oc->faceRegions[offi], rk)
                            {
                                areaB += polygonArea(oc->faceRegions[offi][rk]);
                            }

                            nCompared++;

                            const bool same =
                            (
                                cc.faceRegions[ffi].size()
                             == oc->faceRegions[offi].size()
                             && mag(areaA - areaB) <= 1e-6*max(areaA, SMALL)
                            );

                            if (!same)
                            {
                                nMismatch++;

                                if (dumpCell_ >= 0 && nMismatch <= 6)
                                {
                                    Info<< "REGDIFF face " << facei
                                        << " anchor regions " << cc.faceRegions[ffi].size()
                                        << " area " << areaA
                                        << " other regions " << oc->faceRegions[offi].size()
                                        << " area " << areaB << endl;
                                }
                            }

                            if (dumpCell_ >= 0 && nCompared % 2000 == 0)
                            {
                                Info<< "REGDIFF compared " << nCompared
                                    << " mismatched " << nMismatch << endl;
                            }
                        }

                        bool keptAny = false;

                        forAll(cc.faceRegions[ffi], ri)
                        {
                            const pointField& r = cc.faceRegions[ffi][ri];

                            if (r.size() < 3)
                            {
                                nDropSize++;
                                continue;
                            }

                            point rc = Zero;

                            forAll(r, k)
                            {
                                rc += r[k];
                            }

                            rc /= scalar(r.size());

                            face sub(r.size());

                            forAll(r, k)
                            {
                                sub[k] = pointOf(cc, r[k]);
                            }

                            const label pieceA =
                                cc.pieceCell[cc.faceRegionPiece[ffi][ri]];

                            //- Never -1 on an internal face: setFace reads a
                            //- negative neighbour as a boundary face and
                            //- polyTopoChange then cannot place it, aborting with
                            //- "Did not determine new position for face ...".
                            //- A cut neighbour starts at its own first piece and
                            //- the region match below overwrites it.
                            label pieceB = -1;

                            if (nei >= 0)
                            {
                                //- The anchor is the owner when the owner is
                                //- cut. When only the NEIGHBOUR is cut the
                                //- anchor is the neighbour, and the other side
                                //- of the face is the owner - which is uncut
                                //- here. Taking the neighbour's own first piece
                                //- instead made pieceB == pieceA, dropped the
                                //- sub-face, and removed the face: 120 faces,
                                //- and with them the volume they enclose.
                                if (cutA < 0)
                                {
                                    pieceB = own;
                                }
                                else
                                {
                                    pieceB =
                                        (cutB < 0
                                      ? nei
                                      : cuts[cutOf[cutB]].pieceCell[0]);
                                }
                            }

                            if (oc && offi >= 0)
                            {
                                //- The same polygon, so its centre names the
                                //- region the neighbour grouped it into
                                scalar best = GREAT;

                                forAll(oc->faceRegions[offi], rj)
                                {
                                    const pointField& ro = oc->faceRegions[offi][rj];

                                    if (ro.size() < 3)
                                    {
                                        continue;
                                    }

                                    point roc = Zero;

                                    forAll(ro, k)
                                    {
                                        roc += ro[k];
                                    }

                                    roc /= scalar(ro.size());

                                    const scalar d = mag(roc - rc);

                                    if (d < best)
                                    {
                                        best = d;
                                        pieceB =
                                            oc->pieceCell
                                            [
                                                oc->faceRegionPiece[offi][rj]
                                            ];
                                    }
                                }
                            }

                            //- A sub-face with the same cell on both sides is
                            //- inside a piece, not a face of the mesh. Dropping
                            //- it is right; emitting it aborts in polyTopoChange
                            //- with "owner 1113 neighbour 1113".
                            if (pieceA >= 0 && pieceA == pieceB)
                            {
                                nDropSame++;

                                if (dumpCell_ >= 0 && nDropSame <= 6)
                                {
                                    Info<< "DROPSAME face " << facei
                                        << " own " << own << " nei " << nei
                                        << " cutA " << cutA << " cutB " << cutB
                                        << " anchor " << anchor << " other " << other
                                        << " ffi " << ffi << " ri " << ri
                                        << " nReg " << cc.faceRegions[ffi].size()
                                        << " pieceA " << pieceA
                                        << " pieceB " << pieceB
                                        << " anchorCell " << cc.celli
                                        << " anchorPieces " << cc.pieceCell
                                        << " otherCell " << (oc ? oc->celli : -1)
                                        << " otherPieces "
                                        << (oc ? oc->pieceCell : labelList())
                                        << endl;
                                }

                                continue;
                            }

                            setFace
                            (
                                sub,
                                (keptAny ? -1 : facei),
                                pieceA,
                                pieceB,
                                patchi,
                                zonei,
                                zoneFlip,
                                "loop-sub"
                            );

                            keptAny = true;
                        }

                        //- Every sub-face was inside one piece, so the old face
                        //- has no successor and must go, or mapPolyMesh has no
                        //- position for it
                        if (!keptAny)
                        {
                            removedFaces.insert(facei);
                            topoChange.removeFace(facei, -1);
                        }

                        continue;
                    }
                }
            }

            if (!splitOf.found(facei))
            {
                //- A face the surface does not cross lies wholly on one side, so
                //- it moves from the background cell to the piece on that side
                const point x(mesh_.faceCentres()[facei]);

                setFace
                (
                    fcs[facei],
                    facei,
                    cellAt(own, x),
                    (nei >= 0 ? cellAt(nei, x) : -1),
                    patchi,
                    zonei,
                    zoneFlip,
                    "whole"
                );

                continue;
            }

            const faceSplit& split = splits[splitOf[facei]];

            forAll(split.subFaces, sfi)
            {
                const face& sub = split.subFaces[sfi];

                //- The centre of the sub-face, for the side tests. Sub-faces
                //- are convex, so the average is inside them.
                point x(Zero);

                forAll(sub, i)
                {
                    x += coordOf(sub[i]);
                }

                x /= scalar(sub.size());

                setFace
                (
                    sub,
                    (sfi == 0 ? facei : -1),
                    cellAt(own, x),
                    (nei >= 0 ? cellAt(nei, x) : -1),
                    patchi,
                    zonei,
                    zoneFlip,
                    "split"
                );
            }
        }
    }


    //- 4. The faces the pieces of one cell share: the interface where the two
    //- pieces are in different regions, an ordinary internal face where they
    //- are not
    forAll(cuts, ci)
    {
        const cellCut& c = cuts[ci];

        //- Probe: how many times does section 4 reach the cell, and with how
        //- many shared faces each time?
        if (loopCut_ && dumpCell_ >= 0 && c.celli == dumpCell_)
        {
            Info<< "PROBECUT ci " << ci << " cell " << c.celli
                << " sharedFaces " << c.sharedFaces.size()
                << " pieces " << c.pieces.size()
                << " nCells " << mesh_.nCells() << endl;

        }

        forAll(c.sharedFaces, si)
        {
            const cutCellPieceFace& shared = c.sharedFaces[si];

            if (loopCut_ && dumpCell_ >= 0 && c.celli == dumpCell_)
            {
                const face& pf = c.pieces[shared.pieceA].poly.faces()
                    [shared.faceInA];
                const point p0 = c.pieces[shared.pieceA].poly.points()[pf[0]];

                Info<< "   PROBESHARED si " << si
                    << " pieceA " << shared.pieceA
                    << " pieceB " << shared.pieceB
                    << " faceInA " << shared.faceInA
                    << " size " << pf.size()
                    << " p0 " << p0 << endl;
            }

            if (c.pieceGroup[shared.pieceA] == c.pieceGroup[shared.pieceB])
            {
                //- The two pieces are one cell now, so this face is inside it
                //- and is not a face of the mesh at all
                const_cast<report&>(info_).nMergedFaces++;
                continue;
            }

            const bool separates =
            (
                c.pieceZone[shared.pieceA] != c.pieceZone[shared.pieceB]
            );

            //- Pieces that share a face and are in one region are one cell, so
            //- a shared face that survives has to be an interface face
            if (!separates)
            {
                FatalErrorInFunction
                    << "Cells " << c.pieceCell[shared.pieceA] << " and "
                    << c.pieceCell[shared.pieceB] << " of one region in cell "
                    << c.celli << " share a face that the merge did not take."
                    << exit(FatalError);
            }

            //- Take the polygon from the piece inside the surface where the
            //- face separates the regions, so its normal points out of the
            //- solid, which is the direction the interface zone records
            label from = shared.pieceA;
            label to = shared.pieceB;
            label fromFace = shared.faceInA;

            if (separates && c.pieceZone[shared.pieceB] == insideZone)
            {
                from = shared.pieceB;
                to = shared.pieceA;
                fromFace = shared.faceInB;
            }

            const cutCellPolyhedron& poly = c.pieces[from].poly;
            const face& f = poly.faces()[fromFace];

            face meshFace(f.size());

            forAll(f, i)
            {
                meshFace[i] = pointOf(c, poly.points()[f[i]]);
            }

            //- The face's own region decides the zone. A face whose region is
            //- unknown falls back to the single zone when there is one, so a
            //- single-region surface behaves exactly as it did before.
            label zone = interfaceZone;

            if
            (
                shared.patch >= 0
             && shared.patch < interfaceZoneOf.size()
            )
            {
                zone = interfaceZoneOf[shared.patch];
            }

            setFace
            (
                meshFace,
                -1,
                c.pieceCell[from],
                c.pieceCell[to],
                -1,
                (separates ? zone : -1),
                false,
                "iface"
            );
        }
    }


    //- Faces grown for conformality belong to no cut cell, so the emit loop above
    //- never reaches them - it is driven by cuts. Their cells are uncut, so the
    //- sub-face keeps the face's own owner and neighbour.
    {
        label nGrownEmitted = 0;
        label nGrownInternal = 0;

        forAllConstIters(splitOf, iter)
        {
            const label facei = iter.key();

            if (doneFaces.found(facei))
            {
                continue;
            }

            const label own = mesh_.faceOwner()[facei];
            const label nei = mesh_.isInternalFace(facei)
              ? mesh_.faceNeighbour()[facei]
              : -1;

            const label patchi = mesh_.boundaryMesh().whichPatch(facei);
            const label zonei = mesh_.faceZones().whichZone(facei);

            bool zoneFlip = false;

            if (zonei >= 0)
            {
                const faceZone& fz = mesh_.faceZones()[zonei];
                zoneFlip = fz.flipMap()[fz.whichFace(facei)];
            }

            const faceSplit& split = splits[iter()];

            forAll(split.subFaces, sfi)
            {
                setFace
                (
                    split.subFaces[sfi],
                    (sfi == 0 ? facei : -1),
                    own,
                    nei,
                    patchi,
                    zonei,
                    zoneFlip,
                    "grown"
                );
            }

            nGrownEmitted++;

            if (nei >= 0)
            {
                nGrownInternal++;
            }
        }

        Info<< "  grown faces emitted: " << nGrownEmitted
            << " (internal " << nGrownInternal << ")" << endl;
    }


    {
        label nMissed = 0;
        label nRemoved = 0;

        forAll(cuts, ci)
        {
            for (const label facei : cells[cuts[ci].celli])
            {
                if (!emittedFaces.found(facei))
                {
                    if (removedFaces.found(facei))
                    {
                        nRemoved++;
                    }

                    if (nMissed < 6)
                    {
                        Info<< "MISSED face " << facei
                            << " owner " << mesh_.faceOwner()[facei]
                            << " patch " << mesh_.boundaryMesh().whichPatch(facei)
                            << " inSplit " << splitOf.found(facei)
                            << " removed " << removedFaces.found(facei) << endl;
                    }

                    nMissed++;
                }
            }
        }

        Info<< "faces of cut cells not emitted: " << nMissed
            << " (removed: " << nRemoved << "), sub-faces dropped by size "
            << nDropSize << ", by same piece " << nDropSame
            << ", cut-cell faces with no regions " << nNoRegions << endl;
    }

    //- 5. Apply, and take the mesh with it
    autoPtr<mapPolyMesh> map = topoChange.changeMesh(mesh_, false);

    mesh_.updateMesh(map());

    //- 6. Drop the points nothing uses. The subdivision makes a crossing point
    //- for every sub-face it builds, and on the loop path the polygon emitted for
    //- a face the cutter has regions for is the cutter's region, not the
    //- subdivision's sub-face, so most of those points are never referenced by a
    //- face. They are points and nothing else - no face, no cell, no zone - so
    //- removing them takes nothing out of the mesh but the points themselves.
    boolList used(mesh_.points().size(), false);
    label nUsed = 0;

    forAll(mesh_.faces(), facei)
    {
        for (const label p : mesh_.faces()[facei])
        {
            if (!used[p])
            {
                used[p] = true;
                nUsed++;
            }
        }
    }

    const label nUnused = mesh_.points().size() - nUsed;

    if (nUnused == 0)
    {
        return map;
    }

    polyTopoChange dropUnused(mesh_);

    forAll(used, p)
    {
        if (!used[p])
        {
            dropUnused.removePoint(p, -1);
        }
    }

    autoPtr<mapPolyMesh> dropMap = dropUnused.changeMesh(mesh_, false);

    mesh_.updateMesh(dropMap());

    Info<< "points no face used, dropped from the mesh: " << nUnused << nl;

    return dropMap;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mapPolyMesh> Foam::cutCellMesher::cut()
{
    Info<< nl << "Cutting with surface" << nl << "--------------------" << nl
        << endl;

    info_ = report();
    info_.nCells = mesh_.nCells();

    //- Collapse the cut crossings that fall too close to a mesh point, which
    //- moves those points. This has to happen before anything reads the mesh
    //- geometry, so that the cut and every query about it see one set of
    //- positions.
    snapPoints();

    //- The parent volumes have to be taken before the cut, since afterwards the
    //- cells are the pieces and their parents are gone
    backgroundVolumes_ = mesh_.cellVolumes();

    //- Every cell starts as its own background cell
    backgroundCell_.setSize(mesh_.nCells());

    forAll(backgroundCell_, celli)
    {
        backgroundCell_[celli] = celli;
    }

    List<cellCut> cuts;
    findCutCells(cuts);

    //- The collapse is what makes a cell graze. It puts a mesh point exactly on
    //- the plane of a facet, so the surface passes through a corner of every cell
    //- that point is a corner of, and such a cell's two pieces would meet along a
    //- line, with no face to put between them. A mesh point is shared, so the
    //- choice is per point and not per cell: the moved corners of the refused
    //- cells are put back and the cut is taken again.
    //-
    //- The collapse is kept only where that clears the refusals. Where it does
    //- not, the refusals were not the collapse's doing - a piece that is malformed
    //- or a face with a hole - and putting the points back would buy nothing while
    //- costing the slivers the collapse removed, so the collapsed positions are
    //- restored exactly and the original cut stands. Measured: on the 320-triangle
    //- sphere 12 points and 72 refused cells, and the reversion cuts 272 of 272
    //- where the collapse cut 200 and failed the case's gate; on `06_oblique_box`
    //- 4 refused cells remain either way, so its cut is left as it was.
    if (info_.nSnapped && info_.nDegenerate)
    {
        const pointField snappedPts(mesh_.points());
        const pointField snappedDeltas(snapDelta_);
        const label snappedCount = info_.nSnapped;
        const scalar snapDistance = info_.maxSnapDistance;
        const label nSnapRejected = info_.nSnapRejected;

        boolList restore(mesh_.points().size(), false);
        label nRestore = 0;

        forAll(info_.degenerate, i)
        {
            for (const label p : mesh_.cellPoints()[info_.degenerate[i]])
            {
                if (mag(snapDelta_[p]) > 0 && !restore[p])
                {
                    restore[p] = true;
                    nRestore++;
                }
            }
        }

        if (nRestore)
        {
            pointField pts(mesh_.points());

            forAll(restore, p)
            {
                if (restore[p])
                {
                    pts[p] -= snapDelta_[p];
                    snapDelta_[p] = Zero;
                }
            }

            mesh_.movePoints(pts);

            label nStillSnapped = 0;
            scalar maxStillSnapped = 0;

            forAll(snapDelta_, p)
            {
                const scalar d = mag(snapDelta_[p]);

                if (d > 0)
                {
                    nStillSnapped++;
                    maxStillSnapped = max(maxStillSnapped, d);
                }
            }

            info_ = report();
            info_.nSnapped = nStillSnapped;
            info_.maxSnapDistance = maxStillSnapped;
            info_.nSnapRejected = nSnapRejected;

            findCutCells(cuts);

            //- A refusal that is a face with a hole is not the collapse's
            //- doing, so it must not veto the reversion of the refusals that
            //- are: with the veto left whole, one geometric refusal kept the
            //- collapse's own damage in place everywhere. Measured on
            //- 06_oblique_box at the tool's own collapse: 50 refusals, of which
            //- 4 were holes and 46 were the collapse's, all 46 reverting once
            //- the 4 stopped vetoing - cells cut 68 to 112 and the interface
            //- 0.648324617528 to 1.1315630264, which is what the same case cuts
            //- with the collapse off.
            if (info_.nDegenerate - info_.nRefusedFaceWithHole > 0)
            {
                //- Not the collapse's doing. Put the collapsed positions back and
                //- take the original cut again, so that the result is the one the
                //- collapse was applied for rather than one that lost it and
                //- gained nothing.
                mesh_.movePoints(snappedPts);
                snapDelta_ = snappedDeltas;

                info_ = report();
                info_.nSnapped = snappedCount;
                info_.maxSnapDistance = snapDistance;
                info_.nSnapRejected = nSnapRejected;

                findCutCells(cuts);
            }
        }
    }

    if (info_.nSnapped)
    {
        Info<< "Intersections collapsed onto a point: " << info_.nSnapped
            << ", largest move " << info_.maxSnapDistance << nl;
    }

    if (info_.nSnapRejected)
    {
        Info<< "  collapses refused as too large      : " << info_.nSnapRejected
            << nl;
    }

    Info<< "Cells cut into >2 pieces     : " << info_.nManyPieces
        << " (most " << info_.maxPieces << ") " << info_.manyPieces << nl;

    Info<< "Cells crossed by the surface : " << info_.nCandidate << nl
        << "Cells cut                    : " << info_.nCut << nl
        << "Cells refused as degenerate  : " << info_.nDegenerate << nl
        << "  of those, piece malformed   : " << info_.nRefusedMalformed << nl
        << "  of those, face has no area  : " << info_.nRefusedDegenerateFace << nl
        << "  of those, face has a hole  : " << info_.nRefusedFaceWithHole << nl
        << "Cells reported, left as is   : " << info_.nUnsupported << endl;

    forAll(info_.degenerate, i)
    {
        Info<< "    cell " << info_.degenerate[i] << ": "
            << info_.degenerateReason[i] << ", and it was left whole" << endl;
    }

    forAll(info_.unsupported, i)
    {
        Info<< "    cell " << info_.unsupported[i] << ": "
            << info_.unsupportedReason[i] << endl;
    }

    if (cuts.empty())
    {
        Info<< "Nothing to cut" << endl;
        return autoPtr<mapPolyMesh>();
    }

    //- Every group of pieces is a cell of the cut mesh. One group per cut cell
    //- keeps the label of that cell, so the rest are the cells that are added.
    backgroundCell_.setSize(mesh_.nCells() + info_.nGroups - info_.nCut);

    autoPtr<mapPolyMesh> map = assemble(cuts);

    //- A sub-cell's parent is the background cell it came from. One piece keeps
    //- the label of that cell, so the rest are the ones to record.
    forAll(cuts, ci)
    {
        const cellCut& c = cuts[ci];

        forAll(c.pieceCell, pi)
        {
            backgroundCell_[c.pieceCell[pi]] = c.celli;
        }
    }

    //- backgroundCell_ was grown with setSize, which leaves the new entries
    //- uninitialised, and only the cells appearing in pieceCell have been
    //- written. A cell the fill missed would carry garbage into its parent list
    //- and from there into backgroundParents, so report it here rather than let
    //- a -1 reach the output file.
    label nNoParent = 0;

    forAll(backgroundCell_, celli)
    {
        const label b = backgroundCell_[celli];

        if (b < 0 || b >= backgroundVolumes_.size())
        {
            if (nNoParent < 5)
            {
                Info<< "  no background cell for cell " << celli
                    << ": parent reads " << b << endl;
            }

            ++nNoParent;
            backgroundCell_[celli] = 0;
        }
    }

    if (nNoParent)
    {
        Info<< "  cells with no background parent: " << nNoParent
            << " of " << backgroundCell_.size()
            << " - the piece fill missed them" << endl;
    }

    //- The parent list starts as the parent map, one background cell each
    backgroundParents_.setSize(mesh_.nCells());

    forAll(backgroundParents_, celli)
    {
        backgroundParents_[celli] = labelList(1, backgroundCell_[celli]);
    }

    if (sliverVolumeRatio_ > 0 && minVolumeRatio_ > 0)
    {
        mergeSlivers();
    }

    Info<< "Cells after cutting          : " << mesh_.nCells() << nl
        << "Points added                 : " << info_.nNewPoints << nl
        << "Interface area               : " << info_.interfaceArea << endl;

    if (info_.nMerged)
    {
        Info<< "Sub-cells merged by region   : " << info_.nMerged
            << " of " << info_.nPieces << ", in " << info_.nGroups
            << " cell(s), dropping " << info_.nMergedFaces << " face(s)" << endl;
    }

    //- Cut cells by the number of planes of the surface inside them, which says
    //- how much of the mesh needed the arrangement
    for (label n = 1; n < info_.byPlanes.size(); n++)
    {
        if (info_.byPlanes[n])
        {
            Info<< "Cut by " << n << " plane(s)              : "
                << info_.byPlanes[n] << " cell(s)" << endl;
        }
    }

    if (info_.nSingleRegion)
    {
        Info<< "Cut cells whose pieces are all one region: "
            << info_.nSingleRegion << endl;
    }

    if (info_.nZoneRepaired)
    {
        Info<< "Cut cells whose zones the surface's planes re-decided, a"
            << " shared face having both its pieces in one region: "
            << info_.nZoneRepaired << endl;
    }

    if (info_.nZoneUnresolved)
    {
        Info<< "Cut cells whose zones could not be separated by the surface's"
            << " planes either, so the centre's verdict stands: "
            << info_.nZoneUnresolved << endl;
    }

    return map;
}


void Foam::cutCellMesher::writeSnapDisplacement() const
{
    //- Nothing moved, so there is nothing to record. An all-zero field would
    //- also be written in a compact form that the readers have to expand.
    if (snapDelta_.empty() || info_.nSnapped == 0)
    {
        return;
    }

    pointIOField snapDisplacement
    (
        IOobject
        (
            "snapDisplacement",
            mesh_.time().timeName(),
            "polyMesh",             // local, so it sits with the mesh
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        snapDelta_
    );

    snapDisplacement.write();
}


void Foam::cutCellMesher::mergeSlivers()
{
    const scalarField& V = mesh_.cellVolumes();
    const pointField& C = mesh_.cellCentres();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    const auto& zoneMesh = mesh_.cellZones();
    const auto& faceZoneMesh = mesh_.faceZones();

    //- Parents' volume of each cell, so eligibility means the same thing for a
    //- cell that has already been merged
    scalarField parentVolume(V.size(), 0);

    forAll(parentVolume, celli)
    {
        forAll(backgroundParents_[celli], pi)
        {
            parentVolume[celli] +=
                backgroundVolumes_[backgroundParents_[celli][pi]];
        }
    }

    DynamicList<label> pierced;

    //- Why a sliver ends without a candidate, counted per reason so the run
    //- says what it could not do instead of dropping it silently. The last two
    //- exist only because of removeFaces and go with it.
    label cntNotInternal = 0;
    label cntZonedFace = 0;
    label cntOtherRegion = 0;
    label cntMultiShared = 0;
    label cntEdgeRefused = 0;
    label cntInvalid = 0;
    label cntNoCandidate = 0;
    label cntRankingDiffers = 0;
    label cntOnlyRef2 = 0;

    forAll(V, celli)
    {
        const label zonei = zoneMesh.whichZone(celli);

        if (zonei < 0 || parentVolume[celli] <= 0)
        {
            continue;
        }

        if (V[celli] > sliverVolumeRatio_*parentVolume[celli])
        {
            continue;
        }

        label bestFace = -1;
        label bestConcave = labelMax;
        scalar bestRatio = -GREAT;
        scalar bestArea = -1;
        scalar bestHealthyArea = -1;
        label bestHealthyFace = -1;

        forAll(mesh_.cells()[celli], k)
        {
            const label facei = mesh_.cells()[celli][k];

            if (!mesh_.isInternalFace(facei))
            {
                ++cntNotInternal;
                continue;
            }

            //- Never remove a face of a zone: the interface lives in one, and
            //- removing it would put the two regions in one cell
            if (faceZoneMesh.whichZone(facei) >= 0)
            {
                ++cntZonedFace;
                continue;
            }

            const label nb = (own[facei] == celli ? nei[facei] : own[facei]);

            //- Only within one region
            if (zoneMesh.whichZone(nb) != zonei)
            {
                ++cntOtherRegion;
                continue;
            }

            //- The two filters that were here - a single shared face, and three
            //- faces at every edge of it - existed only because removeFaces
            //- refused those cases. The assembly has no such requirement, and the
            //- interface is protected by the same-region test above.

            //- The boundary the merged cell would have: the faces of both,
            //- minus the one between them
            List<pointField> boundary(0);

            forAll(mesh_.cells()[celli], i)
            {
                const label fi = mesh_.cells()[celli][i];

                if (fi != facei)
                {
                    boundary.append
                    (
                        pointField(mesh_.faces()[fi].points(mesh_.points()))
                    );
                }
            }

            forAll(mesh_.cells()[nb], i)
            {
                const label fi = mesh_.cells()[nb][i];

                if (fi != facei)
                {
                    boundary.append
                    (
                        pointField(mesh_.faces()[fi].points(mesh_.points()))
                    );
                }
            }

            const point ref =
            (
                C[celli]*V[celli] + C[nb]*V[nb]
            )/max(V[celli] + V[nb], VSMALL);

            const mergeQuality q(evaluateMerge(boundary, ref));

            if (!q.valid)
            {
                ++cntInvalid;
                continue;
            }

            const scalar area = mag(mesh_.faceAreas()[facei]);

            //- REF-2 ranks by the largest shared face area with a healthy
            //- neighbour, healthy meaning at or above the sliver threshold
            //- of its own parent
            if
            (
                V[nb] > sliverVolumeRatio_*parentVolume[nb]
             && area > bestHealthyArea
            )
            {
                bestHealthyArea = area;
                bestHealthyFace = facei;
            }

            const bool better =
            (
                q.concaveFaces < bestConcave
             || (q.concaveFaces == bestConcave && q.minPyramidRatio > bestRatio)
             || (q.concaveFaces == bestConcave && q.minPyramidRatio == bestRatio
                 && area > bestArea)
            );

            if (better)
            {
                bestConcave = q.concaveFaces;
                bestRatio = q.minPyramidRatio;
                bestArea = area;
                bestFace = facei;
            }
        }

        if (bestFace >= 0)
        {
            pierced.append(bestFace);
        }
        else
        {
            ++cntNoCandidate;
        }

        if (bestFace >= 0 && bestHealthyFace != bestFace)
        {
            ++cntRankingDiffers;

            if (cntRankingDiffers <= 5)
            {
                Info<< "  ranking differs on cell " << celli
                    << ": quality picks face " << bestFace
                    << " area " << bestArea
                    << ", REF-2 picks face " << bestHealthyFace
                    << " area " << bestHealthyArea << endl;
            }
        }
        else if (bestFace < 0 && bestHealthyFace >= 0)
        {
            ++cntOnlyRef2;
        }
    }

    Info<< nl << "Sliver merge candidates" << nl
        << "    selected                     : " << pierced.size() << nl
        << "    no candidate                 : " << cntNoCandidate << nl
        << "    ranking differs              : " << cntRankingDiffers
        << " (REF-2 would pick where quality did not: " << cntOnlyRef2 << ")" << nl
        << "    rejected, not an internal face    : " << cntNotInternal << nl
        << "    rejected, face is in a zone       : " << cntZonedFace << nl
        << "    rejected, neighbour in another zone: " << cntOtherRegion << nl
        << "    rejected, more than one shared face: " << cntMultiShared << nl
        << "    rejected, edge has under three faces: " << cntEdgeRefused
        << "   <- a removeFaces requirement" << nl
        << "    rejected, merged cell invalid      : " << cntInvalid << nl
        << "    (the last two reasons disappear with polyTopoChange)" << endl;

    if (pierced.empty())
    {
        return;
    }

    const label nBefore = mesh_.nCells();

    //- Resolve the pierced faces into groups. The cells either side of a pierced
    //- face become one cell, and groups chain when a sliver's partner has a pick
    //- of its own. The master is the largest-volume cell, so the label and the
    //- cellZone that survive belong to a substantial cell.
    labelList cellRegionMaster(mesh_.nCells());
    labelList cellRegion(mesh_.nCells(), -1);
    labelList group(mesh_.nCells());

    forAll(group, c)
    {
        group[c] = c;
    }

    scalarField groupVolume(V);

    auto groupRoot = [&](label c) -> label
    {
        while (group[c] != c)
        {
            group[c] = group[group[c]];
            c = group[c];
        }

        return c;
    };

    forAll(pierced, i)
    {
        const label f = pierced[i];

        label a = groupRoot(own[f]);
        label b = groupRoot(nei[f]);

        if (a == b)
        {
            continue;
        }

        if (groupVolume[a] < groupVolume[b])
        {
            std::swap(a, b);
        }

        group[b] = a;
        groupVolume[a] += groupVolume[b];
    }

    //- How big each group is. A cell that ends up merged records its group,
    //- the master included: the parent routing below uses cellRegion >= 0 to
    //- recognise a merged cell, and leaving the master at -1 left its parent
    //- list empty, which wrote a -1 into backgroundParents.
    labelList groupSize(mesh_.nCells(), 0);

    forAll(group, c)
    {
        ++groupSize[groupRoot(c)];
    }

    label nRegions = 0;

    forAll(group, c)
    {
        const label r = groupRoot(c);
        cellRegionMaster[c] = r;

        if (groupSize[r] > 1)
        {
            cellRegion[c] = r;
        }

        if (r == c && groupSize[r] > 1)
        {
            ++nRegions;
        }
    }

    if (!nRegions)
    {
        return;
    }

    //- Grow the groups that are still short of their target. REF-2 step 6: a
    //- group absorbs the neighbour it shares the largest boundary face with,
    //- within its region, and keeps going until it reaches minVolumeRatio times
    //- the volume of the background cells it came from. The target accumulates
    //- as the group grows, and the master volume accumulates with it.
    scalarField groupTarget(mesh_.nCells(), 0);

    //- Summed over the group's members: the grouping ran before this existed, so
    //- a target set only for the root would be minVolumeRatio times one member's
    //- parents, which a group of six cells reaches immediately.
    forAll(group, c)
    {
        groupTarget[groupRoot(c)] += minVolumeRatio_*parentVolume[c];
    }

    label nShort = 0;
    label nGrown = 0;

    const label maxGrow = 100;

    for (label iter = 0; iter < maxGrow; iter++)
    {
        bool grew = false;

        forAll(group, c)
        {
            if (groupRoot(c) != c)
            {
                continue;
            }

            if (groupVolume[c] >= groupTarget[c])
            {
                continue;
            }

            //- The best boundary face of the group that leads into the same
            //- region, by area, which is REF-2's rule for growing
            label bestFace = -1;
            label bestOther = -1;
            scalar bestArea = -1;

            forAll(group, m)
            {
                if (groupRoot(m) != c)
                {
                    continue;
                }

                const label zonei = zoneMesh.whichZone(m);

                if (zonei < 0)
                {
                    continue;
                }

                forAll(mesh_.cells()[m], k)
                {
                    const label f = mesh_.cells()[m][k];

                    if (!mesh_.isInternalFace(f))
                    {
                        continue;
                    }

                    const label other = (own[f] == m ? nei[f] : own[f]);
                    const label ro = groupRoot(other);

                    if (ro == c || zoneMesh.whichZone(other) != zonei)
                    {
                        continue;
                    }

                    const scalar area = mag(mesh_.faceAreas()[f]);

                    if (area > bestArea)
                    {
                        bestArea = area;
                        bestFace = f;
                        bestOther = other;
                    }
                }
            }

            if (bestFace < 0)
            {
                continue;
            }

            //- Absorb it, keeping this group's root so the surviving label is
            //- the one the sliver already had
            const label otherRoot = groupRoot(bestOther);
            group[otherRoot] = c;
            groupVolume[c] += groupVolume[otherRoot];
            groupTarget[c] += groupTarget[otherRoot];
            ++nGrown;
            grew = true;
        }

        if (!grew)
        {
            break;
        }
    }

    forAll(group, c)
    {
        if (groupRoot(c) == c && groupVolume[c] < groupTarget[c])
        {
            ++nShort;
        }
    }

    Info<< "Sliver merge growth          : " << nGrown << " cell(s) absorbed"
        << ", " << nShort << " group(s) still short of the target" << endl;

    //- Parents of every region, its master's first so that entry zero stays
    //- the primary parent
    //- The growth re-roots groups, so the master of every cell has to be
    //- recomputed here: taken during the grouping it is stale, and it splits a
    //- grown group's parents across several old roots. Each new cell then listed
    //- only its share, and backgroundParents under-reported - a cell that
    //- absorbed twelve background cells listed two.
    forAll(cellRegionMaster, celli)
    {
        cellRegionMaster[celli] = groupRoot(celli);
    }

    List<labelList> masterParents(mesh_.nCells(), labelList(0));

    forAll(cellRegionMaster, celli)
    {
        const label m = cellRegionMaster[celli];

        if (cellRegion[celli] < 0 || m < 0 || m >= mesh_.nCells())
        {
            continue;
        }

        if (celli == m)
        {
            masterParents[m] = backgroundParents_[celli];
        }
    }

    forAll(cellRegionMaster, celli)
    {
        const label m = cellRegionMaster[celli];

        if (cellRegion[celli] < 0 || m < 0 || m >= mesh_.nCells())
        {
            continue;
        }

        forAll(backgroundParents_[celli], pi)
        {
            const label p = backgroundParents_[celli][pi];

            if (!masterParents[m].found(p))
            {
                masterParents[m].append(p);
            }
        }
    }

    polyTopoChange topoChange(mesh_);

    //- The topology, expressed directly. A master keeps the faces whose other
    //- side is outside its group; the internal ones go. removeCell only marks the
    //- cell, so a face of a removed cell has to be re-owned explicitly or the map
    //- leaves it with a cell that no longer exists and it vanishes - which for a
    //- boundary face would put a hole in the domain, and for an interface face
    //- would break the CHT contract.
    //- The cells. A cell's face list is not declared anywhere: it is composed
    //- from the faces, so the master needs no modifyCell - there is no overload
    //- that sets faces, only zones. The non-masters are removed into it, and the
    //- faces below decide what the merged cell is bounded by.
    forAll(group, c)
    {
        const label r = groupRoot(c);

        if (r != c)
        {
            topoChange.removeCell(c, r);
        }
    }

    //- Note the local faceZoneMesh above shadows the type of the same name,
    //- so the type cannot be written here
    const auto& fzMesh = mesh_.faceZones();

    forAll(mesh_.faces(), f)
    {
        const bool isInternal = mesh_.isInternalFace(f);
        const label o = own[f];
        const label n = (isInternal ? nei[f] : -1);
        const label ro = groupRoot(o);
        const label rn = (isInternal ? groupRoot(n) : -1);

        //- Both sides are one cell now, so the face is inside it and goes. This
        //- covers more than the pierced faces: two cells of a group may share
        //- another face as well
        if (isInternal && ro == rn)
        {
            topoChange.removeFace(f, -1);
            continue;
        }

        //- Nothing moved
        if (o == ro && (!isInternal || n == rn))
        {
            continue;
        }

        const label zonei = fzMesh.whichZone(f);

        //- The flip belongs to the zone and is indexed by the face's position
        //- within it, not by the face label
        bool zoneFlip = false;

        if (zonei >= 0)
        {
            const faceZone& fz = fzMesh[zonei];
            zoneFlip = fz.flipMap()[fz.whichFace(f)];
        }

        if (isInternal)
        {
            //- An internal face has to stay ordered owner < neighbour. The
            //- winding is a geometric ordering as seen from the owner, so when
            //- the two masters swap places the face has to be reversed too.
            if (ro < rn)
            {
                topoChange.modifyFace
                (
                    mesh_.faces()[f], f, ro, rn, false, -1, zonei, zoneFlip
                );
            }
            else
            {
                face reversed(mesh_.faces()[f]);
                reversed = reversed.reverseFace();

                topoChange.modifyFace
                (
                    reversed, f, rn, ro, false, -1, zonei, zoneFlip
                );
            }
        }
        else
        {
            topoChange.modifyFace
            (
                mesh_.faces()[f], f, ro, -1, false,
                mesh_.boundaryMesh().whichPatch(f), zonei, zoneFlip
            );
        }
    }

    autoPtr<mapPolyMesh> mergeMap = topoChange.changeMesh(mesh_, false);

    mesh_.updateMesh(mergeMap());

    //- Route the parent lists through the change: a merged cell takes its
    //- region's parents, every other cell keeps its own
    //- Route the parents driven by the new cells. Iterating the forward map
    //- missed 61 of 556 cells on the oblique case and left their parent at -1:
    //- a new index is not guaranteed to be visited. The reverse map names an old
    //- cell behind each new one, and that cell's group holds the parents, so
    //- every new cell is covered.
    //- cellMap is indexed by the NEW cell and holds an old one: removeCell sets
    //- cellMap_[removed] = -2 and encodes the merge target in reverseCellMap_.
    //- Iterating its index therefore covers every new cell, one old cell each.
    const labelList& cellMap = mergeMap().cellMap();

    labelList newBackground(mesh_.nCells(), -1);
    List<labelList> newParents(mesh_.nCells(), labelList(0));

    label nNoOld = 0;

    forAll(newBackground, newCelli)
    {
        if
        (
            newCelli >= cellMap.size()
         || cellMap[newCelli] < 0
         || cellMap[newCelli] >= cellRegionMaster.size()
        )
        {
            ++nNoOld;
            continue;
        }

        const label oldCelli = cellMap[newCelli];
        const label m = cellRegionMaster[oldCelli];

        if
        (
            cellRegion[oldCelli] >= 0
         && m >= 0
         && m < masterParents.size()
         && !masterParents[m].empty()
        )
        {
            newParents[newCelli] = masterParents[m];
            newBackground[newCelli] = masterParents[m][0];
        }
        else
        {
            newParents[newCelli] = backgroundParents_[oldCelli];
            newBackground[newCelli] = backgroundCell_[oldCelli];
        }
    }

    if (nNoOld)
    {
        Info<< "  new cells with no old cell behind them: " << nNoOld << endl;
    }

    label nStillNeg = 0;

    forAll(newBackground, i)
    {
        if (newBackground[i] < 0)
        {
            if (nStillNeg < 5)
            {
                Info<< "  new cell " << i << " still has no background parent" << endl;
            }

            ++nStillNeg;
        }
    }

    if (nStillNeg)
    {
        Info<< "  cells without a background parent after the merge: "
            << nStillNeg << endl;
    }

    //- How many cells the merge left below the trigger, and the worst of them.
    //- The counters above are per neighbour and per reason; this is the count
    //- that matters, and it is what makes the next change measurable.
    {
        const scalarField& vol = mesh_.cellVolumes();

        label nStill = 0;

        forAll(newBackground, celli)
        {
            scalar parentV = 0;

            forAll(newParents[celli], pi)
            {
                const label p = newParents[celli][pi];

                if (p >= 0 && p < backgroundVolumes_.size())
                {
                    parentV += backgroundVolumes_[p];
                }
            }

            if (parentV <= 0)
            {
                continue;
            }

            const scalar ratio = vol[celli]/parentV;

            if (ratio <= sliverVolumeRatio_)
            {
                if (nStill < 5)
                {
                    const label oldCelli =
                    (
                        celli < cellMap.size() && cellMap[celli] >= 0
                      ? cellMap[celli]
                      : -1
                    );

                    Info<< "  below the trigger: cell " << celli
                        << " ratio " << ratio
                        << " parents " << newParents[celli].size()
                        << " vol " << vol[celli]
                        << " parentSum " << parentV;

                    if (oldCelli >= 0)
                    {
                        //- What the growth would have compared against, for the
                        //- one old cell the map names behind this new one
                        Info<< " oldCell " << oldCelli
                            << " oldParentVolume " << parentVolume[oldCelli];
                    }

                    Info<< endl;
                }

                ++nStill;
            }
        }

        Info<< "Cells still below the trigger : " << nStill
            << " of " << mesh_.nCells()
            << " (trigger " << sliverVolumeRatio_ << ")" << endl;
    }

    backgroundCell_ = newBackground;
    backgroundParents_ = newParents;
    info_.nAbsorbed = nBefore - mesh_.nCells();

    //- Anything still at -1 was never written by the routing above
    label nNeverWritten = 0;

    forAll(newBackground, i)
    {
        if (newBackground[i] < 0)
        {
            if (nNeverWritten < 5)
            {
                Info<< "  new cell " << i << " has no background parent" << endl;
            }

            ++nNeverWritten;
        }
    }

    Info<< "Cells absorbed as slivers    : " << info_.nAbsorbed
        << " in " << nRegions << " region(s)" << endl;
}


void Foam::cutCellMesher::writeBackgroundParents() const
{
    labelListIOList backgroundParents
    (
        IOobject
        (
            "backgroundParents",
            mesh_.time().timeName(),
            "polyMesh",
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        backgroundParents_
    );

    backgroundParents.write();
}


void Foam::cutCellMesher::writeBackgroundCell() const
{
    labelIOList backgroundCell
    (
        IOobject
        (
            "backgroundCell",
            mesh_.time().timeName(),
            "polyMesh",             // local, so it sits with the mesh
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        backgroundCell_
    );

    backgroundCell.write();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
