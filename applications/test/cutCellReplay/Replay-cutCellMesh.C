/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.com
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

Application
    Replay-cutCellMesh

Description
    Replay one cell through the loop cut, from a dump the mesher wrote with the
    dumpCell keyword. Temporary: it exists to fix a failing cell without running
    a case, and once that cell is a fixture in Test-cutCellMesh it goes.

\*---------------------------------------------------------------------------*/

#include "cutCellCutter.H"
#include "cutCellPolyhedron.H"
#include "dictionary.H"
#include "IFstream.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        Info<< "usage: Replay-cutCellMesh <dump>" << endl;
        return 1;
    }

    IFstream is(argv[1]);
    dictionary d(is);

    const scalar tol(d.get<scalar>("tol"));
    const label celli(d.get<label>("celli"));
    const pointField pts(d.get<pointField>("cellPoints"));
    const faceList fcs(d.get<faceList>("cellFaces"));
    const List<pointField> tris(d.get<List<pointField>>("tris"));

    Info<< "cell " << celli << ": " << pts.size() << " points, "
        << fcs.size() << " faces, " << tris.size() << " triangles, tol "
        << tol << nl << endl;

    const cutCellPolyhedron cell(pts, fcs);

    labelList patchOf(tris.size(), 0);

    //- The loop's own report, per cell edge
    cutCellLoopDebug(true);

    cutCellResult result;

    const cutState state = cutCellByLoop(cell, tris, patchOf, tol, result);

    Info<< nl << "state ";

    switch (state)
    {
        case cutState::cut: Info<< "cut"; break;
        case cutState::untouched: Info<< "untouched"; break;
        case cutState::refused: Info<< "refused"; break;
    }

    Info<< ", pieces " << result.pieces.size()
        << ", shared faces " << result.sharedFaces.size() << endl;

    scalar total = 0;

    forAll(result.pieces, pi)
    {
        const cutCellPolyhedron& p = result.pieces[pi].poly;

        Info<< "  piece " << pi << ": faces " << p.nFaces()
            << " badEdges " << p.nBadEdges()
            << " volume " << p.volume() << endl;

        total += p.volume();
    }

    Info<< "  volumes sum to " << total
        << ", cell volume " << cell.volume()
        << ", difference " << mag(total - cell.volume()) << endl;

    return 0;
}


// ************************************************************************* //
