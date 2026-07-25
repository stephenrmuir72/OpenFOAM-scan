/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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

Class
    Foam::RugeStubenAMG

Group
    grpLduMatrixSolvers

Description
    Classical Ruge-Stuben algebraic multigrid solver (prototype).

    Characteristics:
      - Requires a symmetric (SPD) matrix; asymmetric matrices are
        rejected with a FatalError.
      - Coarsening: classical two-pass Ruge-Stuben C/F splitting with
        classical direct interpolation (see RugeStubenCoarsening).
      - Coarse matrix creation: Galerkin (RAP) triple product.
      - Smoother: Gauss-Seidel.
      - Type of cycle: V-cycle.
      - Coarsest level approximated with additional Gauss-Seidel
        sweeps (see Scope note below).

    Parallel running: every coarsening level (including the finest)
    is built with the interface-aware RugeStubenCoarsening
    constructor, which performs a distributed CLJP-style C/F
    splitting and a distributed Galerkin (RAP) coarse-matrix
    assembly, producing new coarse-level processor interfaces
    (reusing processorGAMGInterface/processorGAMGInterfaceField) as
    needed. The whole V-cycle (smoothing, residual, restriction,
    prolongation) therefore runs genuinely distributed at every
    level - no gather-to-master step remains. The coarsest level is
    solved with additional Gauss-Seidel sweeps (not a distributed
    direct/Krylov solve) - see the Scope note below.

    Scope (prototype): no fvSolution-tuned defaults; intended to
    validate hierarchy/unit correctness rather than production
    performance. Rare multi-rank "corner" (2-hop) Galerkin
    contributions are dropped (see RugeStubenCoarsening); this
    affects convergence RATE only, not correctness, since smoothing
    and residual evaluation always use the exact level matrix. The
    coarsest level is solved with a fixed number of additional
    Gauss-Seidel sweeps rather than a dedicated distributed direct or
    Krylov solve. Coarsening-level rank imbalance (a rank reaching
    zero local coarse cells while others still have many) is not
    specifically handled. See RugeStubenCoarsening for the
    coarsening/interpolation/Galerkin-assembly algorithm.

SourceFiles
    RugeStubenAMG.C

\*---------------------------------------------------------------------------*/

#ifndef RugeStubenAMG_H
#define RugeStubenAMG_H

#include "lduMatrix.H"
#include "lduPrimitiveMesh.H"
#include "RugeStubenCoarsening.H"
#include "PtrDynList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class RugeStubenAMG Declaration
\*---------------------------------------------------------------------------*/

class RugeStubenAMG
:
    public lduMatrix::solver
{
    // Private Data

        //- Strength-of-connection threshold theta, in [0,1]
        scalar strengthThreshold_;

        //- Stop coarsening once a level has this many cells or fewer
        //- (global total across all ranks)
        label nCellsInCoarsestLevel_;

        //- Number of pre-smoothing sweeps
        label nPreSweeps_;

        //- Number of post-smoothing sweeps
        label nPostSweeps_;

        //- Number of Gauss-Seidel sweeps used to (approximately) solve
        //- the coarsest level
        label nCoarsestSweeps_;

        //- Hierarchy of coarsening levels. levels_[i] coarsens the
        //- level-i matrix (matrix_ if i==0, else
        //- levels_[i-1].coarseMatrix()) to level i+1. Built with the
        //- interface-aware constructor at every level (including
        //- level 0), so each level's coarse interfaces/coefficients
        //- are available via levels_[i].coarseInterfaces()/
        //- coarseInterfaceBouCoeffs() for use by level i+1.
        PtrDynList<RugeStubenCoarsening> levels_;


    // Private Member Functions

        //- Return the matrix at the given level (0 = finest)
        const lduMatrix& matrixLevel(const label leveli) const;

        //- Return the processor-interface fields at the given level
        //- (0 = finest, using the real mesh interfaces_; else the
        //- coarse interfaces built by levels_[leveli-1])
        const lduInterfaceFieldPtrsList& interfacesLevel
        (
            const label leveli
        ) const;

        //- Return the interface boundary coefficients at the given
        //- level, parallel to interfacesLevel(leveli)
        const FieldField<Field, scalar>& interfaceBouCoeffsLevel
        (
            const label leveli
        ) const;

        //- Smooth psi at the given level for nSweeps Gauss-Seidel sweeps
        void smoothLevel
        (
            const label leveli,
            scalarField& psi,
            const scalarField& source,
            const label nSweeps
        ) const;

        //- Calculate the residual at the given level
        void residualLevel
        (
            const label leveli,
            const scalarField& psi,
            const scalarField& source,
            scalarField& res
        ) const;

        //- Recursive V-cycle starting at the given level
        void Vcycle
        (
            const label leveli,
            scalarField& psi,
            const scalarField& source
        ) const;

        //- No copy construct
        RugeStubenAMG(const RugeStubenAMG&) = delete;

        //- No copy assignment
        void operator=(const RugeStubenAMG&) = delete;


public:

    //- Runtime type information
    TypeName("RugeStubenAMG");


    // Constructors

        //- Construct from matrix components and solver controls
        RugeStubenAMG
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    //- Destructor
    virtual ~RugeStubenAMG() = default;


    // Member Functions

        //- Read the control parameters from the controlDict_
        virtual void readControls();

        //- Number of coarsening levels actually built
        //- (0 = matrix was already at or below nCellsInCoarsestLevel)
        label nLevels() const noexcept
        {
            return levels_.size();
        }

        //- Access to a coarsening level (for inspection/testing)
        const RugeStubenCoarsening& level(const label i) const
        {
            return levels_[i];
        }

        //- Solve the matrix with this solver
        virtual solverPerformance solve
        (
            scalarField& psi,
            const scalarField& source,
            const direction cmpt = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
