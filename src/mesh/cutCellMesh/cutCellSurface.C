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

#include "cutCellSurface.H"
#include "IFstream.H"
#include "dictionary.H"
#include "Time.H"
#include "boundBox.H"
#include "error.H"
#include "geometricSurfacePatch.H"
#include "surfZoneIdentifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

void Foam::cutCellSurface::boxTriangles
(
    const word& name,
    const dictionary& dict,
    pointField& points,
    List<triFace>& tris
)
{
    const point mn(dict.lookup("min"));
    const point mx(dict.lookup("max"));

    if (mag(mx - mn) < SMALL)
    {
        FatalErrorInFunction
            << "Geometry entry '" << name << "' has no extent: min " << mn
            << " max " << mx << exit(FatalError);
    }

    //- The eight corners, indexed by the bits of their coordinates
    points.setSize(8);

    for (label i = 0; i < 8; i++)
    {
        points[i] = point
        (
            (i & 1) ? mx.x() : mn.x(),
            (i & 2) ? mx.y() : mn.y(),
            (i & 4) ? mx.z() : mn.z()
        );
    }

    const point c(0.5*(mn + mx));

    tris.clear();

    //- Six faces, each the four corners on one side of one axis, split into two
    //- triangles. The winding is settled by comparing the triangle normal with
    //- the direction away from the box centre, so the orientation is checked
    //- rather than reasoned out by hand.
    static const label uu[4] = {0, 1, 1, 0};
    static const label vv[4] = {0, 0, 1, 1};

    for (direction axis = 0; axis < 3; axis++)
    {
        const direction a1 = (axis + 1) % 3;
        const direction a2 = (axis + 2) % 3;

        for (label side = 0; side < 2; side++)
        {
            label quad[4];

            for (label k = 0; k < 4; k++)
            {
                label bits = 0;

                if (side)   bits |= (1 << axis);
                if (uu[k])  bits |= (1 << a1);
                if (vv[k])  bits |= (1 << a2);

                quad[k] = bits;
            }

            auto addTri = [&](const label i0, const label i1, const label i2)
            {
                const vector n
                (
                    (points[i1] - points[i0]) ^ (points[i2] - points[i0])
                );

                if ((n & (c - points[i0])) > 0)
                {
                    //- Winding points inward, so reverse it
                    tris.append(triFace(i0, i2, i1));
                }
                else
                {
                    tris.append(triFace(i0, i1, i2));
                }
            };

            addTri(quad[0], quad[1], quad[2]);
            addTri(quad[0], quad[2], quad[3]);
        }
    }
}


void Foam::cutCellSurface::fileTriangles
(
    const word& name,
    const dictionary& dict,
    const Time& runTime,
    pointField& points,
    List<triFace>& tris
)
{
    const fileName path(surfacePath(name, dict, runTime));

    if (!isFile(path))
    {
        FatalErrorInFunction
            << "Cannot find surface file " << path << " for geometry entry '"
            << name << "'" << exit(FatalError);
    }

    const triSurface s(path);

    points = s.points();
    tris.setSize(s.size());

    forAll(s, i)
    {
        tris[i] = triFace(s[i][0], s[i][1], s[i][2]);
    }
}


void Foam::cutCellSurface::checkClosedOriented
(
    const triSurface& s,
    const word& name
)
{
    //- Two triangles share an edge when the undirected pair of their point
    //- labels agrees. For a closed, consistently oriented surface each such
    //- edge is then met exactly twice, and the two triangles traverse it in
    //- opposite directions.
    Map<label> nEdges(2*s.size());
    Map<label> firstSign(2*s.size());

    forAll(s, i)
    {
        const labelledTri& t = s[i];

        for (direction e = 0; e < 3; e++)
        {
            const label a = t[e];
            const label b = t[(e + 1) % 3];
            const label sgn = (a < b) ? 1 : -1;

            const label key =
                min(a, b)*max(label(1), label(s.points().size())) + max(a, b);

            Map<label>::iterator iter = nEdges.find(key);

            if (iter == nEdges.end())
            {
                nEdges.insert(key, 1);
                firstSign.insert(key, sgn);
            }
            else if (iter() == 1)
            {
                if (firstSign[key] == sgn)
                {
                    FatalErrorInFunction
                        << "Surface for geometry entry '" << name
                        << "' is not consistently oriented: the edge from point "
                        << min(a, b) << " to point " << max(a, b)
                        << " is traversed the same way by two triangles"
                        << exit(FatalError);
                }

                iter()++;
            }
            else
            {
                FatalErrorInFunction
                    << "Surface for geometry entry '" << name
                    << "' is not closed: the edge from point " << min(a, b)
                    << " to point " << max(a, b)
                    << " is shared by more than two triangles"
                    << exit(FatalError);
            }
        }
    }

    const label nPoints = s.points().size();

    forAllConstIters(nEdges, iter)
    {
        if (iter() != 2)
        {
            FatalErrorInFunction
                << "Surface for geometry entry '" << name << "' is not closed:"
                << " the edge from point " << iter.key()/nPoints
                << " to point " << iter.key() % nPoints << " is met by "
                << iter() << " triangle(s), not 2" << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fileName Foam::cutCellSurface::surfacePath
(
    const word& name,
    const dictionary& dict,
    const Time& runTime
)
{
    const fileName file
    (
        dict.lookupOrDefault<fileName>("file", fileName(name + ".stl"))
    );

    return runTime.constant()/"triSurface"/file;
}


void Foam::cutCellSurface::readRegions
(
    const word& name,
    const dictionary& dict,
    const Time& runTime,
    wordList& names,
    labelList& regionOf
)
{
    names.clear();
    regionOf.clear();

    const word type(dict.lookup("type"));

    //- A dictionary shape is generated here and carries no regions, so the
    //- caller falls back to a single zone.
    if (type != "triSurfaceMesh" && type != "closedTriSurfaceMesh")
    {
        return;
    }

    const fileName path(surfacePath(name, dict, runTime));

    //- The region of every face is the file's own record and the reader keeps
    //- it whatever the format is: the solids of an ASCII STL, and the per-face
    //- 0x<token> of the native triSurface format.
    const triSurface s(path);

    regionOf.setSize(s.size());

    forAll(s, i)
    {
        regionOf[i] = s[i].region();
    }

    //- The reader generates a name for a region the file did not name itself:
    //- patch0, patch1, ... from geometricSurfacePatch, and zone0, ... from
    //- surfZoneIdentifier. So a name that is one of those means the file said
    //- nothing, and the native triSurface format is the one whose names sit on
    //- its faces as tokens instead.
    bool declared = false;

    forAll(s.patches(), patchi)
    {
        const word& patchn = s.patches()[patchi].name();

        if
        (
            patchn != geometricSurfacePatch::defaultName(patchi)
         && patchn != surfZoneIdentifier::defaultName(patchi)
        )
        {
            declared = true;
            break;
        }
    }

    if (declared)
    {
        //- An STL's solids, or any other format that names its own regions.
        forAll(s.patches(), patchi)
        {
            names.append(s.patches()[patchi].name());
        }
    }
    else if (path.ext() == "tri")
    {
        //- The names are the tokens in the file, in order of first appearance,
        //- prefixed with "zone" as the reader does: a line ending 0xwall names
        //- its region zonewall. That assignment of ids is the one the region
        //- indices above use, so the two are consistent, and the gate checks it.
        IFstream is(path);

        string line;

        while (is.good() && is.getLine(line))
        {
            const string::size_type at = line.rfind("0x");

            if (at == string::npos)
            {
                continue;
            }

            const word token(line.substr(at + 2));

            if (token.empty())
            {
                continue;
            }

            const word zoneName("zone" + token);

            if (findIndex(names, zoneName) < 0)
            {
                names.append(zoneName);
            }
        }
    }

    Info<< "  surface regions : " << names.size();

    forAll(names, i)
    {
        Info<< " \"" << names[i] << "\"";
    }

    Info<< endl;
}


Foam::triSurface Foam::cutCellSurface::makeSurface
(
    const word& name,
    const dictionary& dict,
    const Time& runTime
)
{
    const word type(dict.lookup("type"));

    pointField points;
    List<triFace> tris;

    if (type == "searchableBox")
    {
        boxTriangles(name, dict, points, tris);
    }
    else if (type == "triSurfaceMesh" || type == "closedTriSurfaceMesh")
    {
        fileTriangles(name, dict, runTime, points, tris);
    }
    else
    {
        FatalErrorInFunction
            << "Cannot cut with geometry entry '" << name << "' of type '"
            << type << "'." << nl
            << "  Supported in this slice: searchableBox, triSurfaceMesh,"
            << " closedTriSurfaceMesh." << nl
            << "  Curved dictionary shapes are generated offline as STL, see"
            << " 07 section 6.11." << exit(FatalError);
    }

    const triSurface s(tris, points);

    checkClosedOriented(s, name);

    //- The sign of the volume is what distinguishes an outward oriented closed
    //- surface from an inward oriented one
    const scalar vol = signedVolume(s);

    if (vol <= 0)
    {
        FatalErrorInFunction
            << "Surface for geometry entry '" << name << "' has signed volume "
            << vol << ", so it is inward oriented. Cutting needs the triangles"
            << " to face outwards, and it needs to know which side of the"
            << " surface is inside." << exit(FatalError);
    }

    Info<< "Cut surface " << name << ": " << s.size() << " triangles, closed,"
        << " outward oriented, volume " << vol << endl;

    return s;
}


Foam::cutCellSurface::cutCellSurface
(
    const word& name,
    const dictionary& dict,
    const Time& runTime
)
:
    name_(name),
    surface_(makeSurface(name, dict, runTime)),
    search_(surface_)
{
    readRegions(name, dict, runTime, regionNames_, regionOf_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pointField Foam::cutCellSurface::triangle(const label i) const
{
    const labelledTri& t = surface_[i];

    pointField pts(3);

    forAll(pts, k)
    {
        pts[k] = surface_.points()[t[k]];
    }

    return pts;
}


Foam::scalar Foam::cutCellSurface::signedVolume(const triSurface& s)
{
    //- Referenced to the centre of the bounding box rather than the origin, so
    //- that a surface far from the origin does not lose the answer to
    //- cancellation between terms of order the coordinate cubed.
    const pointField& pts = s.points();

    if (pts.empty())
    {
        return 0;
    }

    const point ref(boundBox(pts, false).centre());

    scalar sumT = 0;

    forAll(s, i)
    {
        const labelledTri& t = s[i];

        const vector a(pts[t[0]] - ref);
        const vector b(pts[t[1]] - ref);
        const vector c(pts[t[2]] - ref);

        sumT += a & (b ^ c);
    }

    return sumT/6.0;
}


Foam::scalar Foam::cutCellSurface::signedVolume() const
{
    return signedVolume(surface_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
