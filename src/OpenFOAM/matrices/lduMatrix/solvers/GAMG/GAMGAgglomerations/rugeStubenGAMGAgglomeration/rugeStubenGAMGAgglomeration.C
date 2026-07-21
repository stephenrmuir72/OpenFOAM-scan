/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 OpenCFD Ltd.
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

#include "rugeStubenGAMGAgglomeration.H"
#include "RugeStubenCoarsening.H"
#include "lduMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rugeStubenGAMGAgglomeration, 0);

    addToRunTimeSelectionTable
    (
        GAMGAgglomeration,
        rugeStubenGAMGAgglomeration,
        lduMatrix
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rugeStubenGAMGAgglomeration::rugeStubenGAMGAgglomeration
(
    const lduMatrix& matrix,
    const dictionary& controlDict
)
:
    pairGAMGAgglomeration(matrix.mesh(), controlDict),
    strengthThreshold_
    (
        controlDict.getOrDefault<scalar>("strengthThreshold", 0.25)
    )
{
    agglomerate(matrix);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rugeStubenGAMGAgglomeration::agglomerate(const lduMatrix& matrix)
{
    if (matrix.symmetric())
    {
        // Classical Ruge-Stuben strength-of-connection + two-pass C/F
        // splitting on the real finest-level matrix. Purely rank-local
        // (mirrors pairGAMGAgglomeration's own local-only pairing
        // decision - processor coupling is handled afterwards by the
        // generic GAMGInterface/GAMGProcAgglomeration machinery, not by
        // coordinating the partition decision itself across ranks).
        RugeStubenCoarsening coarsening
        (
            matrix,
            strengthThreshold_,
            matrix.mesh().comm()
        );

        const List<RugeStubenCoarsening::pointType>& pointTypes =
            coarsening.pointTypes();
        const labelListList& strongNbrs = coarsening.strongNeighbours();

        const label nFineCells = pointTypes.size();

        auto tfineToCoarse = tmp<labelField>::New(nFineCells, label(-1));
        labelField& fineToCoarse = tfineToCoarse.ref();

        label nCoarseCells = 0;

        // C-points become coarse cells in their own right
        forAll(pointTypes, celli)
        {
            if (pointTypes[celli] == RugeStubenCoarsening::pointType::C)
            {
                fineToCoarse[celli] = nCoarseCells++;
            }
        }

        // F-points are clustered into the coarse cell of their first
        // strong C-neighbour (the classical two-pass RS splitting
        // guarantees every F-point has at least one). This gives a
        // "nearest strong C-neighbour" partition - not the weighted
        // direct-interpolation of RugeStubenCoarsening::
        // interpolationAddr()/interpolationWeights() - but it is
        // sufficient to drive the standard (unweighted) GAMGSolver
        // matrix-agglomeration/interpolation pipeline unchanged.
        forAll(pointTypes, celli)
        {
            if (pointTypes[celli] == RugeStubenCoarsening::pointType::F)
            {
                label parent = -1;

                for (const label nbr : strongNbrs[celli])
                {
                    if (pointTypes[nbr] == RugeStubenCoarsening::pointType::C)
                    {
                        parent = fineToCoarse[nbr];
                        break;
                    }
                }

                if (parent == -1)
                {
                    // Should not normally occur for a valid two-pass RS
                    // splitting (every F-point is guaranteed a strong
                    // C-neighbour); fall back to a singleton coarse
                    // cell rather than failing.
                    parent = nCoarseCells++;
                }

                fineToCoarse[celli] = parent;
            }
        }

        if
        (
            continueAgglomerating
            (
                nCellsInCoarsestLevel_,
                nFineCells,
                nCoarseCells,
                matrix.mesh().comm()
            )
        )
        {
            nCells_[0] = nCoarseCells;
            restrictAddressing_.set(0, tfineToCoarse);

            // Create the coarse mesh addressing for level 0
            agglomerateLduAddressing(0);

            // Agglomerate the face-magnitude field (a generic strength
            // proxy, not real matrix coefficients) so the inherited
            // pair-agglomeration algorithm can continue from level 1
            // onwards, exactly as algebraicPairGAMGAgglomeration does.
            scalarField faceWeights
            (
                meshLevels_[0].lduAddr().upperAddr().size(),
                Zero
            );
            restrictFaceField(faceWeights, mag(matrix.upper())(), 0);

            // Continue with the standard pair-agglomeration algorithm
            // for all subsequent levels (no real coarse matrix
            // coefficients are available to GAMGAgglomeration beyond
            // level 0).
            pairGAMGAgglomeration::agglomerate
            (
                nCellsInCoarsestLevel_,
                1,
                faceWeights,
                true
            );

            return;
        }
    }

    // Either the matrix is not a suitable (symmetric) candidate for RS
    // splitting, or RS splitting did not reduce the cell count enough
    // to continue past level 0 - fall back to standard, plain
    // pair-agglomeration starting at level 0 (same behaviour as
    // algebraicPairGAMGAgglomeration).
    if (matrix.hasLower())
    {
        pairGAMGAgglomeration::agglomerate
        (
            nCellsInCoarsestLevel_,
            0,
            max(mag(matrix.upper()), mag(matrix.lower())),
            true
        );
    }
    else
    {
        pairGAMGAgglomeration::agglomerate
        (
            nCellsInCoarsestLevel_,
            0,
            mag(matrix.upper()),
            true
        );
    }
}


// ************************************************************************* //
