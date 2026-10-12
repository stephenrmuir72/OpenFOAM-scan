/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
    Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "tetOverlapVolume.H"
#include "tetrahedron.H"
#include "polyMesh.H"
#include "OFstream.H"
#include "treeBoundBox.H"
#include "indexedOctree.H"
#include "treeDataCell.H"

#include <cmath>
#include <limits>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tetOverlapVolume, 0);
}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

Foam::treeBoundBox Foam::tetOverlapVolume::pyrBb
(
    const pointField& points,
    const face& f,
    const point& fc
)
{
    treeBoundBox bb(fc);
    bb.add(points, f);

    return bb;
}


// * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

Foam::Tuple2<Foam::scalar, Foam::scalar>
Foam::tetOverlapVolume::cellVolumeMinDecomp
(
    const primitiveMesh& mesh,
    const label celli
)
{
    const scalar fvVolume = mesh.cellVolumes()[celli];
    if (!std::isfinite(fvVolume) || fvVolume <= 0)
    {
        FatalErrorInFunction
            << "Invalid finite-volume cell volume " << fvVolume
            << " at cell " << celli << exit(FatalError);
    }
    const point& centre = mesh.cellCentres()[celli];
    const pointField& points = mesh.points();
    treeBoundBox bounds(centre);
    long double volume = 0;
    label nTets = 0;
    scalar minAltitude = GREAT;
    for (const label facei : mesh.cells()[celli])
    {
        const face& f = mesh.faces()[facei];
        bounds.add(points, f);
        for (label i = 1; i+1 < f.size(); ++i)
        {
            const label a = mesh.faceOwner()[facei] == celli ? f[i] : f[i+1];
            const label b = mesh.faceOwner()[facei] == celli ? f[i+1] : f[i];
            const tetPointRef tet(centre, points[f[0]], points[a], points[b]);
            const scalar tetVolume = tet.mag();
            ++nTets;
            // The clipping engine also skips exactly degenerate face triangles.
            if
            (
                tetVolume == 0
             && magSqr((points[a]-points[f[0]])^(points[b]-points[f[0]])) == 0
            )
            {
                continue;
            }
            if (!std::isfinite(tetVolume) || tetVolume <= 0)
            {
                FatalErrorInFunction
                    << "Invalid CVW tessellation at cell " << celli
                    << ", face " << facei << ", fan triangle " << i
                    << ": oriented volume " << tetVolume
                    << ". A positive first-vertex fan about the FV cell "
                    << "centre is required." << exit(FatalError);
            }
            volume += tetVolume;
            for (label side = 0; side < 4; ++side)
            {
                const scalar areaSqr = magSqr(tet.tri(side).areaNormal());
                if (!std::isfinite(areaSqr) || areaSqr < ROOTVSMALL)
                {
                    FatalErrorInFunction
                        << "Unresolved CVW tetrahedron area at cell " << celli
                        << ", face " << facei << exit(FatalError);
                }
                const scalar area = Foam::sqrt(areaSqr);
                minAltitude = min(minAltitude, 3*tetVolume/area);
            }
        }
    }
    const scalar geometricVolume = scalar(volume);
    const vector span(bounds.span());
    const scalar minSpan = min
    (
        minAltitude, min(span.x(), min(span.y(), span.z()))
    );
    if
    (
        !std::isfinite(geometricVolume) || geometricVolume <= 0
     || !std::isfinite(minSpan) || minSpan <= 0
    )
    {
        FatalErrorInFunction
            << "Unresolved CVW geometric volume at cell " << celli
            << ": volume " << geometricVolume << ", minimum span "
            << minSpan << exit(FatalError);
    }
    const scalar coordinateScale = max(mag(bounds.min()), mag(bounds.max()));
    const scalar relativeRoundoff =
        64*std::numeric_limits<scalar>::epsilon()*(nTets+1)
       *max(scalar(1), coordinateScale/minSpan);
    if
    (
        !std::isfinite(relativeRoundoff) || relativeRoundoff >= 1
    )
    {
        FatalErrorInFunction
            << "Unresolved CVW geometric volume at cell " << celli
            << ": volume " << geometricVolume << ", relative allowance "
            << relativeRoundoff << exit(FatalError);
    }
    return Tuple2<scalar, scalar>(geometricVolume, relativeRoundoff);
}


bool Foam::tetOverlapVolume::cellCellOverlapMinDecomp
(
    const primitiveMesh& meshA,
    const label cellAI,
    const primitiveMesh& meshB,
    const label cellBI,
    const treeBoundBox& cellBbB,
    const scalar threshold
) const
{
    hasOverlapOp overlapCheckOp(threshold);
    cellCellOverlapMinDecomp<hasOverlapOp>
    (
        meshA,
        cellAI,
        meshB,
        cellBI,
        cellBbB,
        overlapCheckOp
    );

    return overlapCheckOp.ok_;
}


Foam::scalar Foam::tetOverlapVolume::cellCellOverlapVolumeMinDecomp
(
    const primitiveMesh& meshA,
    const label cellAI,

    const primitiveMesh& meshB,
    const label cellBI,
    const treeBoundBox& cellBbB
) const
{
    sumOverlapOp overlapSumOp;
    cellCellOverlapMinDecomp<sumOverlapOp>
    (
        meshA,
        cellAI,
        meshB,
        cellBI,
        cellBbB,
        overlapSumOp
    );

    return overlapSumOp.iop_.vol_;
}


Foam::Tuple2<Foam::scalar, Foam::point>
Foam::tetOverlapVolume::cellCellOverlapMomentMinDecomp
(
    const primitiveMesh& meshA,
    const label cellAI,

    const primitiveMesh& meshB,
    const label cellBI,
    const treeBoundBox& cellBbB
) const
{
    sumOverlapMomentOp overlapSumOp;
    cellCellOverlapMinDecomp<sumOverlapMomentOp>
    (
        meshA,
        cellAI,
        meshB,
        cellBI,
        cellBbB,
        overlapSumOp
    );

    return overlapSumOp.iop_.vol_;
}


Foam::labelList Foam::tetOverlapVolume::overlappingCells
(
    const polyMesh& fromMesh,
    const polyMesh& toMesh,
    const label iTo
) const
{
    const indexedOctree<treeDataCell>& treeA = fromMesh.cellTree();

    treeBoundBox bbB(toMesh.points(), toMesh.cellPoints()[iTo]);

    return treeA.findBox(bbB);
}


// ************************************************************************* //
