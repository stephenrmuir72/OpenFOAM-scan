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

Class
    Foam::rugeStubenGAMGAgglomeration

Description
    Prototype GAMGAgglomeration that derives its finest-level (level 0)
    cell partition from a classical Ruge-Stuben strength-of-connection
    analysis and two-pass C/F splitting (Foam::RugeStubenCoarsening),
    instead of pairGAMGAgglomeration's own pairwise measure. Every
    F-point is clustered into the coarse cell of its first strong
    C-neighbour (guaranteed to exist by the classical two-pass
    splitting), giving a disjoint single-parent partition compatible
    with the standard (unweighted, injection-based) GAMGAgglomeration /
    GAMGSolver pipeline - including GAMGProcAgglomeration processor
    coarsening.

    Scope (prototype):
      - Only the finest level (level 0), where the real lduMatrix
        coefficients are available, uses RS strength-of-connection/
        C-F splitting. GAMGAgglomeration itself never has access to
        real coarse-level matrix coefficients (only agglomerated
        face-weight proxies), so all subsequent levels fall back to
        the inherited pairGAMGAgglomeration::agglomerate() pairwise
        algorithm on those proxies - exactly as
        algebraicPairGAMGAgglomeration already does end-to-end.
      - Only symmetric (SPD) matrices trigger the RS-based level-0
        split (matching RugeStubenCoarsening's own scope); asymmetric
        or diagonal-only matrices silently fall back to plain
        pair-agglomeration from level 0.
      - The C/F splitting itself is purely rank-local (mirrors
        pairGAMGAgglomeration's own local-only pairing decision).
        Parallel correctness/convergence comes from the existing
        GAMGInterface/GAMGProcAgglomeration machinery, not from
        coordinating the level-0 partition decision across ranks.
      - This does not (yet) use RugeStubenCoarsening's weighted,
        multi-parent direct-interpolation data
        (interpolationAddr()/interpolationWeights()); it only reuses
        the C/F splitting to build a classical single-parent
        restrictAddressing_ partition. A weighted Galerkin (RAP)
        coarse-matrix path now exists in GAMGAgglomeration
        (coarseFaceLookup()/fill-in-aware topology) and GAMGSolver
        (agglomerateMatrixWeighted()), gated behind
        GAMGAgglomeration::hasInterpolation() - but this class
        deliberately does not populate interpolationAddr_/
        interpolationWeights_ yet: empirical testing found the
        resulting coarse matrix's aggregate row-sum/conditioning
        behaviour was not good enough to use in practice when
        composed with GAMGSolver's deeper V-cycle machinery, and the
        root cause was not fully isolated. Revisit only once that is
        resolved.

SourceFiles
    rugeStubenGAMGAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef rugeStubenGAMGAgglomeration_H
#define rugeStubenGAMGAgglomeration_H

#include "pairGAMGAgglomeration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class rugeStubenGAMGAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class rugeStubenGAMGAgglomeration
:
    public pairGAMGAgglomeration
{
    // Private Data

        //- Strength-of-connection threshold theta, in [0,1], used for
        //- the finest-level Ruge-Stuben C/F splitting
        scalar strengthThreshold_;


    // Private Member Functions

        //- Build the finest-level (level 0) partition from a classical
        //- Ruge-Stuben strength-of-connection + C/F splitting (when
        //- applicable), then hand off to the inherited
        //- pairGAMGAgglomeration::agglomerate() for all subsequent
        //- (proxy face-weight only) levels.
        void agglomerate(const lduMatrix& matrix);

        //- No copy construct
        rugeStubenGAMGAgglomeration
        (
            const rugeStubenGAMGAgglomeration&
        ) = delete;

        //- No copy assignment
        void operator=(const rugeStubenGAMGAgglomeration&) = delete;


public:

    //- Runtime type information
    TypeName("rugeStuben");


    // Constructors

        //- Construct given matrix and controls
        rugeStubenGAMGAgglomeration
        (
            const lduMatrix& matrix,
            const dictionary& controlDict
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
