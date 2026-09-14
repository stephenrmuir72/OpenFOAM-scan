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
#include "cellZone.H"
#include "faceZone.H"
#include "treeDataTriSurface.H"
#include "treeBoundBox.H"
#include "DynamicList.H"
#include "boundBox.H"
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

        const cutState state = cutCellByTriangles(cell, tris, tol_, result);

        if (state == cutState::untouched)
        {
            if (result.nPolygons > 0)
            {
                //- Triangles cross the cell, but no plane of theirs divides it.
                //- Reported rather than silently left alone.
                info_.nUnsupported++;
                info_.unsupported.append(celli);
                info_.unsupportedReason.append
                (
                    "the surface touches it without dividing it"
                );
            }
            continue;
        }

        if (state == cutState::refused || result.pieces.size() < 2)
        {
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

        if (nIn == 0 || nOut == 0)
        {
            //- The surface passes through the cell without separating it, so
            //- every piece is in one region. Counted, since the pieces are
            //- still valid cells.
            info_.nSingleRegion++;
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

        info_.nCut++;
        info_.nPieces += c.pieces.size();

        maxPlanes = max(maxPlanes, c.planes.size());

        found.append(c);
    }

    //- Cut cells by the number of planes crossing them
    info_.byPlanes.setSize(maxPlanes + 1, Zero);

    forAll(found, ci)
    {
        info_.byPlanes[found[ci].planes.size()]++;
    }

    info_.nCandidate = info_.nCut + info_.nUnsupported;

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
    const label interfaceZone = faceZoneID(surfaces_[0].name() + "_interface");

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
    Map<label> firstPoint(2*cuts.size());
    DynamicList<point> newCoord(2*cuts.size());
    DynamicList<label> newLabel(2*cuts.size());
    DynamicList<label> newNext(2*cuts.size());
    DynamicList<point> newA(2*cuts.size());
    DynamicList<point> newB(2*cuts.size());

    auto coordOf = [&](const label p) -> point
    {
        return (p < nPoints ? pts[p] : newCoord[p - nPoints]);
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

        const point x(pa + w*(pb - pa));

        const label key = min(a, b);
        const bool had = firstPoint.found(key);
        const label prev = had ? firstPoint[key] : -1;

        for (label i = prev; i >= 0; i = newNext[i])
        {
            //- The same segment, whichever way round it was met. Two faces
            //- sharing an edge traverse it in opposite orders.
            const bool sameSegment =
            (
                (mag(newA[i] - pa) <= tol_ && mag(newB[i] - pb) <= tol_)
             || (mag(newA[i] - pb) <= tol_ && mag(newB[i] - pa) <= tol_)
            );

            if (sameSegment && mag(newCoord[i] - x) <= tol_)
            {
                return newLabel[i];
            }
        }

        const label ptLabel = topoChange.addPoint(x, -1, -1, true);

        const label i = newCoord.size();

        newCoord.append(x);
        newLabel.append(ptLabel);
        newNext.append(prev);
        newA.append(pa);
        newB.append(pb);

        if (had)
        {
            firstPoint.set(key, i);
        }
        else
        {
            firstPoint.insert(key, i);
        }

        info_.nNewPoints++;

        return ptLabel;
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
    }


    //- A point of a piece, as a mesh point label. A point of a cut lies either
    //- on an edge of the cell, which the cells sharing that edge all use, or on
    //- a face of it, which the face's subdivision has already made, or inside
    //- it, where it belongs to this cell alone.
    DynamicList<point> insideCoord(4);
    DynamicList<label> insideLabel(4);

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

        for (const label facei : cells[c.celli])
        {
            const scalar magArea = mag(mesh_.faceAreas()[facei]);

            if (magArea <= VSMALL)
            {
                continue;
            }

            const vector n(mesh_.faceAreas()[facei]/magArea);
            const point refPoint(mesh_.faceCentres()[facei]);

            if (mag(planeDist(n, refPoint, x)) <= tol_ && splitOf.found(facei))
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

        //- Inside the cell, where the point belongs to this cell alone. A cell
        //- has few of these, so a linear scan is enough.
        forAll(insideCoord, i)
        {
            if (mag(insideCoord[i] - x) <= tol_)
            {
                return insideLabel[i];
            }
        }

        const label p = topoChange.addPoint(x, -1, -1, true);

        insideCoord.append(x);
        insideLabel.append(p);

        info_.nNewPoints++;

        return p;
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

        forAll(c.pieces, pi)
        {
            bool match = true;

            forAll(c.planes, k)
            {
                const scalar d = planeDist
                (
                    c.planes[k].normal(), c.planes[k].refPoint(), x
                );

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
    auto setFace = [&]
    (
        const face& f,
        const label targetFace,     //- existing face to modify, or -1 to add
        label own,
        label nei,
        const label patchi,
        const label zonei,
        bool zoneFlip
    )
    {
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


    //- 1. One piece keeps the label of the background cell, the others are
    //- added. The pieces are ordered negative first against the planes, so for
    //- a cell cut by one plane the first piece is the one inside the surface,
    //- which is the label the mesh has always kept.
    forAll(cuts, ci)
    {
        cellCut& c = const_cast<cellCut&>(cuts[ci]);

        forAll(c.pieces, pi)
        {
            if (pi == 0)
            {
                c.pieceCell[pi] = c.celli;

                topoChange.modifyCell
                (
                    c.celli, labelList(1, c.pieceZone[pi])
                );
            }
            else
            {
                c.pieceCell[pi] = topoChange.addCell
                (
                    -1, -1, -1, c.celli, c.pieceZone[pi]
                );
            }
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
                    zoneFlip
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
                    zoneFlip
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

        forAll(c.sharedFaces, si)
        {
            const cutCellPieceFace& shared = c.sharedFaces[si];

            const bool separates =
            (
                c.pieceZone[shared.pieceA] != c.pieceZone[shared.pieceB]
            );

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

            setFace
            (
                meshFace,
                -1,
                c.pieceCell[from],
                c.pieceCell[to],
                -1,
                (separates ? interfaceZone : -1),
                false
            );
        }
    }


    //- 5. Apply, and take the mesh with it
    autoPtr<mapPolyMesh> map = topoChange.changeMesh(mesh_, false);

    mesh_.updateMesh(map());

    return map;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mapPolyMesh> Foam::cutCellMesher::cut()
{
    Info<< nl << "Cutting with surface" << nl << "--------------------" << nl
        << endl;

    info_ = report();
    info_.nCells = mesh_.nCells();

    //- Every cell starts as its own background cell
    backgroundCell_.setSize(mesh_.nCells());

    forAll(backgroundCell_, celli)
    {
        backgroundCell_[celli] = celli;
    }

    List<cellCut> cuts;
    findCutCells(cuts);

    Info<< "Cells crossed by the surface : " << info_.nCandidate << nl
        << "Cells cut                    : " << info_.nCut << nl
        << "Cells reported, left as is   : " << info_.nUnsupported << endl;

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

    //- Every piece is a cell of the cut mesh. One piece per cut cell keeps the
    //- label of that cell, so the rest are the cells that will be added.
    backgroundCell_.setSize(mesh_.nCells() + info_.nPieces - info_.nCut);

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

    Info<< "Cells after cutting          : " << mesh_.nCells() << nl
        << "Points added                 : " << info_.nNewPoints << nl
        << "Interface area               : " << info_.interfaceArea << endl;

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

    return map;
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
