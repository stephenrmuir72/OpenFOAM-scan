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
    Foam::RugeStubenCoarsening

Description
    Builds a single classical Ruge-Stuben AMG coarsening level from a
    symmetric (SPD) lduMatrix:

      -# Strength-of-connection graph (theta-scaled, per-row max
         opposite-sign-to-diagonal off-diagonal test - handles both the
         positive-diagonal/negative-off-diagonal convention and its
         mirror image, e.g. fvm::laplacian's negative-diagonal/
         positive-off-diagonal assembly).
      -# Classical two-pass C/F splitting:
         - Pass 1: greedy, measure-based (lambda = |S^T|) selection,
           following the standard sequential RS coarsening algorithm
           (see Briggs, Henson & McCormick, "A Multigrid Tutorial",
           2nd ed., Section 8.2, Algorithm "Coarsening").
         - Pass 2: enforces that every strongly-connected F-F pair
           shares a common strong C-point, promoting points to C as
           required.

         An alternative, genuinely distributed CLJP-style
         (Cleary-Luby-Jones-Plassmann) C/F splitting is available via
         a second constructor that additionally takes the matrix's
         processor-interface data. Rather than the sequential
         "pick-the-single-global-max" greedy rule above, CLJP assigns
         every point a fixed weight w_i = lambda_i + tiebreak_i and
         repeatedly selects *all* currently-undecided local maxima
         (w.r.t. still-undecided strong neighbours, including ghost
         neighbours across processor interfaces) as C in parallel -
         a Luby-style maximal-independent-set algorithm requiring only
         nearest-neighbour communication per round. See
         calcPointTypeParallel() for details, and Cleary, Falgout,
         Henson & Jones, "Coarse-grid selection for parallel algebraic
         multigrid", 1998.
      -# Classical direct interpolation weights (Ruge-Stuben 1987),
         with separate scaling of the negative and positive
         off-diagonal parts, restricted to strong coarse neighbours.
      -# Galerkin (RAP) coarse-matrix assembly, assembled from scratch
         via a sparse row-triple-product (no existing OpenFOAM utility
         performs this for weighted, non-injection interpolation).

    Scope (prototype):
      - The coarse-level lduPrimitiveMesh is constructed without
        interfaces in both constructors.
      - The distributed CLJP-style constructor only distributes the
        C/F *splitting* (calcPointType) across processor interfaces.
        The subsequent interpolation-weight and Galerkin coarse-matrix
        assembly stages still only consider purely local row data;
        cross-processor interpolation/Galerkin coupling is not yet
        implemented, so the resulting coarse level does not (yet)
        capture cross-rank coarse coupling. This is a deliberate,
        documented scope limit rather than an oversight.
      - Symmetric (SPD) matrices only; asymmetric matrices are
        rejected with a FatalError.
      - Not optimised: the greedy C/F selection is a simple O(n^2)
        scan rather than a bucket-sorted priority queue. Adequate for
        validating hierarchy/unit correctness on small matrices, not
        intended for production-scale performance.

SourceFiles
    RugeStubenCoarsening.C

\*---------------------------------------------------------------------------*/

#ifndef RugeStubenCoarsening_H
#define RugeStubenCoarsening_H

#include "lduMatrix.H"
#include "lduPrimitiveMesh.H"
#include "Tuple2.H"
#include "autoPtr.H"
#include "FieldField.H"
#include "lduInterfaceFieldPtrsList.H"
#include "boolList.H"
#include "globalIndex.H"
#include "PtrList.H"
#include "processorGAMGInterface.H"
#include "processorGAMGInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class RugeStubenCoarsening Declaration
\*---------------------------------------------------------------------------*/

class RugeStubenCoarsening
{
public:

    //- Classification of a point as coarse (C) or fine (F)
    enum class pointType : uint8_t
    {
        F = 0,
        C = 1
    };


private:

    // Private Data

        //- Strength threshold theta, in [0,1]
        scalar strengthThreshold_;

        //- Per-cell classification (size = nFineCells)
        List<pointType> pointType_;

        //- Strong-dependency graph: strongNeighbours_[i] is the sorted
        //- set S_i of points that point i strongly depends on
        //- (sgn(a_ii)*(-a_ij) >= theta * max_k(sgn(a_ii)*(-a_ik)),
        //- tested against row i)
        labelListList strongNeighbours_;

        //- Fine-cell to coarse-cell map. -1 for F-points.
        labelList fineToCoarse_;

        //- Number of coarse cells (= number of C-points)
        label nCoarseCells_;

        //- Sparse interpolation (prolongation) addressing per fine
        //- cell, in coarse-cell indices. Single-entry {coarse index}
        //- for C-points (identity); multiple entries (over C_i) for
        //- F-points.
        labelListList interpolationAddr_;

        //- Sparse interpolation weights per fine cell, parallel to
        //- interpolationAddr_
        List<scalarField> interpolationWeights_;

        //- The coarse-level addressing (no interfaces)
        autoPtr<lduPrimitiveMesh> coarseMeshPtr_;

        //- The coarse-level Galerkin (RAP) matrix
        autoPtr<lduMatrix> coarseMatrixPtr_;


        // Parallel-only boundary bookkeeping (populated by the
        // interface-aware constructor only; left empty by the plain
        // serial constructor).

            //- Per interface, per face: the equivalent-internal
            //- boundary coefficient (-bouCoeffs[facei])
            List<scalarField> boundaryCoeff_;

            //- Per interface, per face: whether the local faceCell
            //- strongly depends on its ghost neighbour
            List<boolList> ownerStrong_;

            //- Per interface, per face: the ghost neighbour's coarse
            //- GLOBAL index if it is a C-point, -1 if it is F
            List<labelList> boundaryGhostCoarseIndex_;

            //- Per fine cell: every (patchi, facei) processor-boundary
            //- connection touching this cell (strong or weak) - unlike
            //- calcPointTypeParallel()'s local-only cellBndFaces, which
            //- records strong connections only
            List<DynamicList<Tuple2<label, label>>> cellBoundaryFaces_;

            //- This level's coarse-cell global numbering
            autoPtr<globalIndex> coarseGlobalIndexPtr_;

            //- Sparse GHOST interpolation addressing (coarse GLOBAL
            //- indices), parallel to interpolationAddr_/
            //- interpolationWeights_, for fine cells that strongly
            //- depend on a ghost C-neighbour. Empty for purely-local
            //- fine cells (the overwhelming majority).
            labelListList interpolationGhostAddr_;
            List<scalarField> interpolationGhostWeights_;

            //- Coarse-level processor interfaces built by
            //- calcCoarseMatrixParallel(), reusing the existing GAMG
            //- processor-interface/interfaceField pair rather than a
            //- bespoke class - empty unless the interface-aware
            //- constructor was used and cross-processor coarse
            //- coupling was found. The underlying processorGAMGInterface
            //- objects themselves are owned solely by coarseMeshPtr_
            //- (via lduPrimitiveMesh::addInterfaces()/
            //- primitiveInterfaces_) - NOT duplicated here - to avoid
            //- a double-free at destruction time. coarseInterfaceFields_
            //- holds a non-owning reference to each one (see
            //- calcCoarseMatrixParallel()).
            PtrList<processorGAMGInterfaceField> coarseInterfaceFields_;
            lduInterfaceFieldPtrsList coarseInterfaces_;
            FieldField<Field, scalar> coarseInterfaceBouCoeffs_;


    // Private Member Functions

        //- Build, for every cell, the list of (neighbour cell,
        //- off-diagonal coefficient) pairs from both ends of every
        //- face (valid for symmetric matrices only)
        static List<DynamicList<Tuple2<label, scalar>>> buildRowAdjacency
        (
            const lduMatrix& matrix
        );

        //- Build the strong-dependency graph S_i
        void calcStrongNeighbours
        (
            const lduMatrix& matrix,
            const List<DynamicList<Tuple2<label, scalar>>>& rowAdjacency
        );

        //- Classical two-pass Ruge-Stuben C/F splitting
        void calcPointType();

        //- Pass 2 of calcPointType()/calcPointTypeParallel(): enforce
        //- that every strongly-connected (purely local) F-F pair
        //- shares a common strong C-point, promoting points to C as
        //- required. Iterated to a fixed point (bounded by nCells
        //- passes). Operates on pointType_/strongNeighbours_ only, so
        //- (for calcPointTypeParallel) it does not detect/fix F-F
        //- violations whose only common strong neighbour would be a
        //- ghost (cross-processor) point - see
        //- applyFFCommonCNeighbourFixupParallel() for the extension
        //- that covers that case.
        void applyFFCommonCNeighbourFixup();

        //- Extends applyFFCommonCNeighbourFixup() to also detect and
        //- fix F-F pairs that span a processor interface. Checks two
        //- cross-rank sources of a common strong C-neighbour, each
        //- requiring only nearest-neighbour communication:
        //-   (1) another LOCAL cell reaching the SAME ghost through a
        //-       (possibly different) strong boundary face - detected
        //-       purely from locally-exchanged ghost GLOBAL indices;
        //-   (2) a cell that is a strong LOCAL C-neighbour of the
        //-       ghost (from the ghost-owning rank's perspective) and
        //-       is ALSO one of this rank's own ghost connections
        //-       (a "shared corner" cell) - detected via one ragged
        //-       exchange, per round, of GLOBAL indices.
        //- Alternates with applyFFCommonCNeighbourFixup() (local pass)
        //- until no rank changes anything (global reduce), bounded by
        //- nCells+2 rounds. Any F-F pair whose only common
        //- C-neighbour is reachable through a longer (2+ hop) chain
        //- is not detected - a documented, deliberately-scoped
        //- approximation.
        void applyFFCommonCNeighbourFixupParallel
        (
            const List<DynamicList<Tuple2<label, label>>>& cellBndFaces,
            const lduInterfaceFieldPtrsList& interfaces,
            const label comm
        );

        //- For every processor interface, determine whether each
        //- boundary face is a "strong" connection from the owning
        //- (local) row's perspective, i.e. whether the local
        //- faceCell strongly depends on its ghost (across-interface)
        //- neighbour. Folds interfaceBouCoeffs into an extended,
        //- boundary-aware row measure (max magnitude opposite-sign-
        //- to-diagonal off-diagonal, following the same convention-
        //- agnostic sign test as calcStrongNeighbours) used only for
        //- this boundary-face test; does not alter strongNeighbours_.
        //- Persists boundaryCoeff_, ownerStrong_ and
        //- cellBoundaryFaces_ members for later (parallel-only) use.
        void calcBoundaryStrength
        (
            const lduMatrix& matrix,
            const List<DynamicList<Tuple2<label, scalar>>>& rowAdjacency,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        );

        //- Distributed CLJP-style (Cleary-Luby-Jones-Plassmann)
        //- parallel C/F splitting. Requires strongNeighbours_ to
        //- already be computed (local part of S). Uses ghost
        //- exchange across processor interfaces (nearest-neighbour
        //- communication only, no global gather) to classify points
        //- whose strong-connection stencil crosses a processor
        //- boundary consistently with points classified purely
        //- locally.
        void calcPointTypeParallel
        (
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const label comm
        );

        //- Build classical direct interpolation weights
        void calcInterpolation
        (
            const lduMatrix& matrix,
            const List<DynamicList<Tuple2<label, scalar>>>& rowAdjacency
        );

        //- Build classical direct interpolation weights, extended to
        //- also account for strong ghost (across-processor-interface)
        //- C-neighbours - populates interpolationGhostAddr_/
        //- interpolationGhostWeights_ in addition to the usual
        //- (local-only) interpolationAddr_/interpolationWeights_.
        //- Requires boundaryCoeff_/ownerStrong_/
        //- cellBoundaryFaces_/boundaryGhostCoarseIndex_ to already be
        //- populated (calcBoundaryStrength() + calcCoarseGlobalIndex()).
        void calcInterpolationParallel
        (
            const lduMatrix& matrix,
            const List<DynamicList<Tuple2<label, scalar>>>& rowAdjacency
        );

        //- Build this level's coarse-cell global numbering
        //- (coarseGlobalIndexPtr_) and, for every boundary face,
        //- exchange the ghost neighbour's coarse global index (or -1
        //- if it is F) into boundaryGhostCoarseIndex_. Requires
        //- fineToCoarse_/nCoarseCells_ to already be final.
        void calcCoarseGlobalIndex
        (
            const lduInterfaceFieldPtrsList& interfaces,
            const label comm
        );

        //- Assemble the Galerkin coarse matrix and coarse addressing
        void calcCoarseMatrix
        (
            const lduMatrix& matrix,
            const List<DynamicList<Tuple2<label, scalar>>>& rowAdjacency,
            const label comm
        );

        //- Assemble the Galerkin coarse matrix, extended to also
        //- account for cross-processor RAP contributions arising from
        //- boundary (processor-interface) fine-level edges, and to
        //- build the resulting coarse-level processor interfaces
        //- (reusing processorGAMGInterface/processorGAMGInterfaceField).
        //- See the .C file for the documented scope limitation on
        //- rare multi-rank "corner" (2-hop) contributions.
        void calcCoarseMatrixParallel
        (
            const lduMatrix& matrix,
            const List<DynamicList<Tuple2<label, scalar>>>& rowAdjacency,
            const lduInterfaceFieldPtrsList& interfaces,
            const label comm
        );

        //- No copy construct
        RugeStubenCoarsening(const RugeStubenCoarsening&) = delete;

        //- No copy assignment
        void operator=(const RugeStubenCoarsening&) = delete;


public:

    // Constructors

        //- Construct from a symmetric lduMatrix and strength threshold,
        //- performing strength analysis, C/F splitting, interpolation
        //- weight calculation and Galerkin coarse-matrix assembly.
        RugeStubenCoarsening
        (
            const lduMatrix& matrix,
            const scalar strengthThreshold,
            const label comm
        );

        //- Construct from a symmetric lduMatrix, strength threshold
        //- and processor-interface coupling data, performing a
        //- distributed CLJP-style parallel C/F splitting instead of
        //- the classical serial greedy algorithm above (see
        //- calcPointTypeParallel()). Interpolation weights and the
        //- Galerkin coarse matrix are still assembled from purely
        //- local row data (see the class-level Scope note).
        RugeStubenCoarsening
        (
            const lduMatrix& matrix,
            const scalar strengthThreshold,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const label comm
        );


    // Member Functions

        //- Strength threshold used
        scalar strengthThreshold() const noexcept
        {
            return strengthThreshold_;
        }

        //- Per-cell classification (C or F)
        const List<pointType>& pointTypes() const noexcept
        {
            return pointType_;
        }

        //- Strong-dependency sets S_i
        const labelListList& strongNeighbours() const noexcept
        {
            return strongNeighbours_;
        }

        //- Fine-to-coarse map (-1 for F-points)
        const labelList& fineToCoarse() const noexcept
        {
            return fineToCoarse_;
        }

        //- Number of coarse cells
        label nCoarseCells() const noexcept
        {
            return nCoarseCells_;
        }

        //- Number of fine cells
        label nFineCells() const noexcept
        {
            return pointType_.size();
        }

        //- Sparse interpolation addressing (coarse indices) per fine
        //- cell.
        const labelListList& interpolationAddr() const noexcept
        {
            return interpolationAddr_;
        }

        //- Sparse interpolation weights per fine cell, parallel to
        //- interpolationAddr()
        const List<scalarField>& interpolationWeights() const noexcept
        {
            return interpolationWeights_;
        }

        //- The coarse-level lduMesh (addressing only, no interfaces)
        const lduPrimitiveMesh& coarseMesh() const
        {
            return *coarseMeshPtr_;
        }

        //- The coarse-level Galerkin matrix
        const lduMatrix& coarseMatrix() const
        {
            return *coarseMatrixPtr_;
        }

        //- The coarse-level processor interfaces (empty unless the
        //- interface-aware constructor was used and cross-processor
        //- coarse coupling was found)
        const lduInterfaceFieldPtrsList& coarseInterfaces() const noexcept
        {
            return coarseInterfaces_;
        }

        //- Coupling coefficients for coarseInterfaces(), parallel to it
        const FieldField<Field, scalar>& coarseInterfaceBouCoeffs()
        const noexcept
        {
            return coarseInterfaceBouCoeffs_;
        }

        //- Restrict a fine-level field to the coarse level using R = P^T
        void restrictField
        (
            scalarField& coarseFld,
            const scalarField& fineFld,
            const lduInterfaceFieldPtrsList& fineInterfaces
        ) const;

        //- Prolong (interpolate) a coarse-level field to the fine
        //- level using the direct-interpolation weights P
        void prolongField
        (
            scalarField& fineFld,
            const scalarField& coarseFld,
            const lduInterfaceFieldPtrsList& fineInterfaces
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
