/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 OpenFOAM Foundation
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

Class
    Foam::mergeQuality

Description
    The shape of a cell formed by merging two cells, used to rank candidate
    merges rather than to gate them.

    The quality measures come from the sliver logic of the prototype branch
    `cutCellPrototype` (see 07 section 4.1), which evaluated a merged cell on
    the mesh. Here the evaluation takes geometry only - the boundary polygons
    of the merged cell and a reference point inside it - so that it can be
    tested without building a mesh.

    Convexity is deliberately not a criterion. Merging a sliver into a
    neighbour always produces a non-convex cell, so a strict convexity gate
    would forbid every merge; `concaveFaces` is a score instead.

SourceFiles
    cutCellMerge.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cutCellMerge_H
#define Foam_cutCellMerge_H

#include "pointField.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class plane;

/*---------------------------------------------------------------------------*\
                        Class mergeQuality Declaration
\*---------------------------------------------------------------------------*/

//- The shape of a merged cell
struct mergeQuality
{
    //- The boundary closes and encloses a positive volume
    bool valid;

    //- Signed volume of the merged cell
    scalar volume;

    //- Smallest pyramid volume over the merged volume. Small means one face is
    //- much closer to the reference point than the others, which is what a
    //- sliver next to the reference point looks like.
    scalar minPyramidRatio;

    //- Faces whose outward normal points at another boundary face's centre.
    //- A non-zero count is expected for merged cut cells and only ranks them.
    label concaveFaces;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- The shape of the cell whose boundary is `boundary`, taking `reference` as
//- the point the pyramids are built from. Returns valid = false when the
//- boundary encloses no volume.
mergeQuality evaluateMerge
(
    const List<pointField>& boundary,
    const point& reference
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
