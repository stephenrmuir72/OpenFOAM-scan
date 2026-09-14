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

//- Key of the undirected mesh edge between two points
inline label edgeKey(const label a, const label b, const label nPoints)
{
    return min(a, b)*max(label(1), nPoints) + max(a, b);
}


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


//- Classification of a polygon's points against a plane
//  side = -1 if any point is strictly below, +1 if any is strictly above
void polygonSides
(
    const pointField& pts,
    const face& f,
    const vector& n,
    const point& refPoint,
    const scalar tol,
    bool& hasMinus,
    bool& hasPlus
)
{
    forAll(f, i)
    {
        const scalar d = planeDist(n, refPoint, pts[f[i]]);

        if (d < -tol)
        {
            hasMinus = true;
        }

        if (d > tol)
        {
            hasPlus = true;
        }
    }
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

        cutCellPolyhedron minus, plus;
        pointField interface;
        cutCellDetail detail;

        const cutState state = cutCellByTriangles
        (
            cell, tris, tol_, minus, plus, interface, detail
        );

        if (state == cutState::refused)
        {
            info_.nUnsupported++;
            info_.unsupported.append(celli);
            info_.unsupportedReason.append
            (
                "more than one plane crosses it, or the plane does not span it"
            );
            continue;
        }

        if (state == cutState::untouched)
        {
            if (detail.nPolygons > 0)
            {
                //- Triangles cross the cell, but the plane they lie on does not
                //- divide it. Reported rather than silently left alone.
                info_.nUnsupported++;
                info_.unsupported.append(celli);
                info_.unsupportedReason.append
                (
                    "a plane crosses it but does not divide it"
                );
            }
            continue;
        }

        //- Crossed by one plane. The plane normal is the outward normal of the
        //- surface, so the positive piece is the outside.
        pointField pieces(2);
        pieces[0] = minus.centre();
        pieces[1] = plus.centre();

        const boolList pieceInside(surf.search().calcInside(pieces));

        if (pieceInside[0] == pieceInside[1])
        {
            FatalErrorInFunction
                << "Both pieces of cell " << celli << " classify as "
                << (pieceInside[0] ? "inside" : "outside") << " the surface."
                << " The cut is not separating them, which means the geometry"
                << " and the labelling disagree." << exit(FatalError);
        }

        if (!pieceInside[0])
        {
            FatalErrorInFunction
                << "The negative piece of cell " << celli << " classifies as"
                << " outside the surface, so the surface is not outward"
                << " oriented with respect to the cut. Refusing to label the"
                << " regions the wrong way round." << exit(FatalError);
        }

        cellCut c;
        c.celli = celli;
        c.normal = detail.normal;
        c.refPoint = detail.refPoint;
        c.minusZone = insideZone;
        c.plusZone = outsideZone;
        c.plusCelli = -1;

        found.append(c);

        info_.nCut++;
        info_.interfaceArea += polygonArea(interface);
    }

    //- A face can only be crossed by one plane if it is to be split into one
    //- pair of sub-faces, which is what keeps the mesh conformal. Where two
    //- cells cut by different planes share a face, both planes cross it, and
    //- neither cell can be cut in this slice.
    labelHashSet conflicted;

    {
        const pointField& pts = mesh_.points();
        const faceList& fcs = mesh_.faces();

        //- Which cuts cross each face, in the order they were found
        Map<DynamicList<label>> faceCuts;

        forAll(found, ci)
        {
            const cellCut& c = found[ci];

            for (const label facei : mesh_.cells()[c.celli])
            {
                bool hasMinus = false;
                bool hasPlus = false;

                polygonSides
                (
                    pts, fcs[facei], c.normal, c.refPoint, tol_,
                    hasMinus, hasPlus
                );

                if (hasMinus && hasPlus)
                {
                    if (!faceCuts.found(facei))
                    {
                        faceCuts.insert(facei, DynamicList<label>(2));
                    }

                    faceCuts[facei].append(ci);
                }
            }
        }

        forAllConstIters(faceCuts, iter)
        {
            const DynamicList<label>& onFace = iter();

            forAll(onFace, i)
            {
                for (label j = i + 1; j < onFace.size(); j++)
                {
                    const cellCut& a = found[onFace[i]];
                    const cellCut& b = found[onFace[j]];

                    if
                    (
                        !samePlane
                        (
                            a.normal, a.refPoint, b.normal, b.refPoint, tol_
                        )
                    )
                    {
                        conflicted.insert(onFace[i]);
                        conflicted.insert(onFace[j]);
                    }
                }
            }
        }
    }

    cuts.clear();

    forAll(found, ci)
    {
        if (conflicted.found(ci))
        {
            info_.nCut--;
            info_.nUnsupported++;
            info_.unsupported.append(found[ci].celli);
            info_.unsupportedReason.append
            (
                "a plane through a neighbouring cell crosses the same face, so"
                " the face cannot be split once"
            );
        }
        else
        {
            cuts.append(found[ci]);
        }
    }

    info_.nCandidate = info_.nCut + info_.nUnsupported;
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

    //- The interface faceZone must exist before the change: polyTopoChange
    //- validates the zoneIDs it is given against the mesh
    const label interfaceZone = faceZoneID(surfaces_[0].name() + "_interface");

    //- All points, faces and cells of the mesh are taken as they are, and the
    //- cut is expressed as modifications on top of them. That is what keeps the
    //- labels of untouched cells, and so the parent map, simple.
    polyTopoChange topoChange(mesh_);


    //- New points, keyed by the mesh edge they lie on. One chain per edge, since
    //- two different planes crossing one edge are two different points.
    Map<label> edgeFirstPoint(2*cuts.size());
    List<label> nextPoint(2*cuts.size());
    List<point> allPointCoord(2*cuts.size());
    List<label> allPointLabel(2*cuts.size());

    auto newPoint = [&]
    (
        const label a,
        const label b,
        const vector& n,
        const point& refPoint
    ) -> label
    {
        //- Where this plane crosses this mesh edge. Every cell that asks for the
        //- same edge and plane gets the same place on it, so the cells sharing
        //- the edge share the point, which is what makes the mesh conformal.
        const scalar da = planeDist(n, refPoint, pts[a]);
        const scalar db = planeDist(n, refPoint, pts[b]);

        const scalar w = da/(da - db);

        const point x(pts[a] + w*(pts[b] - pts[a]));

        const label key = edgeKey(a, b, nPoints);

        const bool had = edgeFirstPoint.found(key);
        const label prev = had ? edgeFirstPoint[key] : -1;

        for (label i = prev; i >= 0; i = nextPoint[i])
        {
            if (mag(allPointCoord[i] - x) <= tol_)
            {
                return allPointLabel[i];
            }
        }

        const label ptLabel = topoChange.addPoint(x, -1, -1, true);

        const label i = allPointCoord.size();

        allPointCoord.append(x);
        allPointLabel.append(ptLabel);
        nextPoint.append(prev);

        if (had)
        {
            edgeFirstPoint[key] = i;
        }
        else
        {
            edgeFirstPoint.insert(key, i);
        }

        info_.nNewPoints++;

        return ptLabel;
    };


    //- Split a mesh face by a cut plane, in mesh point labels. False when the
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

            const scalar da = planeDist(n, refPoint, pts[a]);
            const scalar db = planeDist(n, refPoint, pts[b]);

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


    //- 2. Split the faces, once each, from the face's own points so that both
    //- of the cells it borders see the same two sub-polygons
    Map<label> faceSplit(2*cuts.size());
    List<face> splitMinus(2*cuts.size());
    List<face> splitPlus(2*cuts.size());

    forAll(cuts, ci)
    {
        const cellCut& c = cuts[ci];

        for (const label facei : cells[c.celli])
        {
            if (faceSplit.found(facei))
            {
                continue;
            }

            face m, p;

            if (splitFace(fcs[facei], c.normal, c.refPoint, m, p))
            {
                faceSplit.insert(facei, splitMinus.size());
                splitMinus.append(m);
                splitPlus.append(p);
            }
        }
    }


    //- 3. One piece keeps the label of the background cell, the other is added.
    //- The negative piece is the inside, the positive the outside.
    forAll(cuts, ci)
    {
        cellCut& c = const_cast<cellCut&>(cuts[ci]);

        c.plusCelli = topoChange.addCell
        (
            -1, -1, -1, c.celli, c.plusZone
        );

        topoChange.modifyCell(c.celli, labelList(1, c.minusZone));
    }


    //- 3. Cells the surface does not cross keep their label, and take the region
    //- their centre falls in. The cells that were cut take their zones from the
    //- two pieces instead.
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
                topoChange.modifyCell(celli, labelList(1, cellRegionZone_[celli]));
            }
        }
    }


    //- 4. The faces of every cell that was cut. A face belongs to a piece rather
    //- than to the background cell, so its owner and neighbour are set from
    //- which side of the plane the piece is on. Untouched cells keep their
    //- labels, so this is the only place the pieces are wired up.
    labelHashSet done;

    auto sideCell = [&](const label celli, const int side) -> label
    {
        forAll(cuts, ci)
        {
            if (cuts[ci].celli == celli)
            {
                return (side < 0 ? cuts[ci].celli : cuts[ci].plusCelli);
            }
        }

        return celli;
    };

    //- Set a face of the new mesh, either an existing face or a new one.
    //-
    //- polyMesh keeps internal faces in upper-triangular order, meaning the
    //- lower cell label is the owner. Cutting can reverse the order of a
    //- background face, because the sub-cell added on the positive side is
    //- labelled above every background cell, so reorder and reverse the face
    //- when that happens. The topology change does this itself only when it
    //- renumbers cells, and this cut removes no cell, so it has to be done here.
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

    forAll(cuts, ci)
    {
        const cellCut& c = cuts[ci];

        for (const label facei : cells[c.celli])
        {
            if (done.found(facei))
            {
                continue;
            }

            done.insert(facei);

            const label own = mesh_.faceOwner()[facei];

            //- faceNeighbour is only defined for the internal faces, so it must
            //- not be read for a boundary face. Reading it there returns
            //- whatever lies past the end of the list, which shows up much
            //- later as a corrupt neighbour in the assembled mesh.
            const label nei =
            (
                mesh_.isInternalFace(facei)
              ? mesh_.faceNeighbour()[facei]
              : -1
            );

            const label patchi = mesh_.boundaryMesh().whichPatch(facei);

            //- Keep whatever faceZone the face was in
            const label zonei = mesh_.faceZones().whichZone(facei);

            bool zoneFlip = false;

            if (zonei >= 0)
            {
                const faceZone& fz = mesh_.faceZones()[zonei];
                zoneFlip = fz.flipMap()[fz.whichFace(facei)];
            }

            if (faceSplit.found(facei))
            {
                const label sI = faceSplit[facei];

                //- Each sub-face keeps the orientation of the face it came from,
                //- so its own and neighbour are the pieces on its own side
                setFace
                (
                    splitMinus[sI],
                    facei,
                    sideCell(own, -1),
                    (nei >= 0 ? sideCell(nei, -1) : -1),
                    patchi,
                    zonei,
                    zoneFlip
                );

                setFace
                (
                    splitPlus[sI],
                    -1,
                    sideCell(own, +1),
                    (nei >= 0 ? sideCell(nei, +1) : -1),
                    patchi,
                    zonei,
                    zoneFlip
                );
            }
            else
            {
                //- A face the plane does not cross lies wholly on one side, so
                //- it moves from the background cell to the piece on that side
                const scalar d = planeDist(c.normal, c.refPoint, mesh_.faceCentres()[facei]);

                const int side = (d > 0 ? +1 : -1);

                setFace
                (
                    fcs[facei],
                    facei,
                    sideCell(own, side),
                    (nei >= 0 ? sideCell(nei, side) : -1),
                    patchi,
                    zonei,
                    zoneFlip
                );
            }
        }
    }


    //- 5. The interface, one internal face per cut cell, oriented from the
    //- inside piece to the outside piece and recorded in the faceZone
    forAll(cuts, ci)
    {
        const cellCut& c = cuts[ci];

        //- Recompute the interface polygon of this cell, in mesh point labels.
        //- The cut is redone here rather than remembered, because the points it
        //- needs are the shared ones, which only exist now.
        const cutCellPolyhedron cell(cellPolyhedron(mesh_, c.celli));

        const boundBox bb(cell.points(), false);
        const labelList cand(surfaces_[0].search().tree().findBox(treeBoundBox(bb.min(), bb.max())));

        List<pointField> tris(cand.size());
        forAll(cand, i)
        {
            tris[i] = surfaces_[0].triangle(cand[i]);
        }

        cutCellPolyhedron minus, plus;
        pointField interface;
        cutCellDetail detail;

        const cutState state = cutCellByTriangles
        (
            cell, tris, tol_, minus, plus, interface, detail
        );

        if (state != cutState::cut)
        {
            FatalErrorInFunction
                << "Cell " << c.celli << " was accepted and is no longer"
                << " accepted on the second pass. The cut is not reproducible,"
                << " so the mesh would be inconsistent." << exit(FatalError);
        }

        const labelList& cPoints = mesh_.cellPoints()[c.celli];

        face f(detail.interfaceSource.size());

        forAll(detail.interfaceSource, i)
        {
            const cutPointSource& src = detail.interfaceSource[i];

            if (src.second() < 0)
            {
                f[i] = cPoints[src.first()];
            }
            else
            {
                f[i] = newPoint
                (
                    cPoints[src.first()],
                    cPoints[src.second()],
                    c.normal,
                    c.refPoint
                );
            }
        }

        //- Oriented from the inside piece to the outside piece, which is the
        //- direction the cut polygon comes out in, and recorded in the zone
        setFace(f, -1, c.celli, c.plusCelli, -1, interfaceZone, false);
    }


    //- 6. Apply, and take the mesh with it
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

    //- Room for the cells that will be added
    backgroundCell_.setSize(mesh_.nCells() + cuts.size());

    autoPtr<mapPolyMesh> map = assemble(cuts);

    //- A sub-cell's parent is the background cell it came from. The negative
    //- piece keeps the label of that cell, the positive piece was added, so
    //- only the added ones still have to be recorded.
    forAll(cuts, ci)
    {
        backgroundCell_[cuts[ci].plusCelli] = cuts[ci].celli;
    }

    Info<< "Cells after cutting          : " << mesh_.nCells() << nl
        << "Points added                 : " << info_.nNewPoints << nl
        << "Interface area               : " << info_.interfaceArea << endl;

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
