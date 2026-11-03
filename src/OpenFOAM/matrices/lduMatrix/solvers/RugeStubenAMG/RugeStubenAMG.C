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

\*---------------------------------------------------------------------------*/

#include "RugeStubenAMG.H"
#include "GaussSeidelSmoother.H"
#include "processorLduInterface.H"
#include "UPstream.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RugeStubenAMG, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<RugeStubenAMG>
        addRugeStubenAMGSymMatrixConstructorToTable_;

    lduMatrix::solver::addasymMatrixConstructorToTable<RugeStubenAMG>
        addRugeStubenAMGAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RugeStubenAMG::RugeStubenAMG
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    ),
    strengthThreshold_(0.25),
    nCellsInCoarsestLevel_(1),
    nPreSweeps_(2),
    nPostSweeps_(2),
    levels_()
{
    readControls();

    const label comm = matrix_.mesh().comm();

    // Eagerly build the full coarsening hierarchy. Every level
    // (including level 0) uses the interface-aware constructor, so
    // processor-boundary coupling is carried through to every coarser
    // level via the previous level's coarseInterfaces()/
    // coarseInterfaceBouCoeffs(). The loop termination is based on the
    // GLOBAL (all-rank) cell count so every rank builds the same
    // number of levels (required since coarse-level communication
    // needs a consistent hierarchy depth across ranks).
    while (true)
    {
        const lduMatrix& curMatrix = matrixLevel(levels_.size());
        const label nCells = curMatrix.diag().size();
        const label nCellsGlobal = returnReduce(nCells, sumOp<label>(), comm);

        if (nCellsGlobal <= nCellsInCoarsestLevel_)
        {
            break;
        }

        autoPtr<RugeStubenCoarsening> levelPtr
        (
            new RugeStubenCoarsening
            (
                curMatrix,
                strengthThreshold_,
                interfaceBouCoeffsLevel(levels_.size()),
                interfacesLevel(levels_.size()),
                comm
            )
        );

        const label nCoarseLocal = levelPtr().nCoarseCells();
        const label nCoarseGlobal =
            returnReduce(nCoarseLocal, sumOp<label>(), comm);

        if (nCoarseGlobal >= nCellsGlobal)
        {
            // No coarsening progress possible: stop here
            break;
        }

        // Stop BEFORE appending a level that leaves any rank with zero
        // coarse cells. A rank fully coarsened away can no longer carry
        // its processor-interface couplings, which risks silently
        // fragmenting the coarsest level into disconnected per-rank
        // islands (e.g. two diagonally-adjacent ranks in a structured
        // decomposition losing their only bridging rank at the coarsest
        // level, leaving them with no path between them in the ldu
        // graph). Since the coarsest level is exactly where the global/
        // near-null-space error mode needs to be resolved (critical for
        // singular/near-singular systems such as an all-Neumann pressure
        // equation), silently losing global connectivity there causes
        // severe convergence degradation. Keep the previous, still
        // fully-populated-per-rank and globally-connected level as the
        // coarsest instead.
        const label nRanksWithZeroCoarseCells =
            returnReduce(label(nCoarseLocal == 0 ? 1 : 0), sumOp<label>(), comm);

        if (nRanksWithZeroCoarseCells > 0)
        {
            break;
        }

        // RSDEBUG: temporary per-level row-sum-zero (constant-
        // preservation) diagnostic. Checks diag[i] + all off-diagonal
        // (internal AND interface, using the internal-equivalent sign
        // convention coeff = -bouCoeffs) sums to (near) zero for the
        // newly-built coarse level, to see whether violations grow
        // with coarsening depth.
        if (debug)
        {
            const lduMatrix& cm = levelPtr().coarseMatrix();
            const lduAddressing& cmAddr = levelPtr().coarseMesh().lduAddr();
            const labelUList& cOwner = cmAddr.lowerAddr();
            const labelUList& cNeighbour = cmAddr.upperAddr();
            const scalarField& cDiag = cm.diag();
            const scalarField& cUpper = cm.upper();
            const scalarField& cLower =
            (
                cm.hasLower() ? cm.lower() : cm.upper()
            );

            scalarField rowSum(cDiag);

            forAll(cOwner, facei)
            {
                rowSum[cOwner[facei]] += cUpper[facei];
                rowSum[cNeighbour[facei]] += cLower[facei];
            }

            const lduInterfaceFieldPtrsList& cInterfaces =
                levelPtr().coarseInterfaces();
            const FieldField<Field, scalar>& cBouCoeffs =
                levelPtr().coarseInterfaceBouCoeffs();

            forAll(cInterfaces, patchi)
            {
                if (cInterfaces.set(patchi))
                {
                    const labelUList& faceCells =
                        cInterfaces[patchi].interface().faceCells();
                    const scalarField& bc = cBouCoeffs[patchi];

                    forAll(faceCells, facei)
                    {
                        rowSum[faceCells[facei]] += -bc[facei];
                    }
                }
            }

            scalar maxRel = 0;
            scalar sumRel = 0;
            label nViol = 0;
            const label n = rowSum.size();

            forAll(rowSum, i)
            {
                const scalar scale = max(mag(cDiag[i]), SMALL);
                const scalar rel = mag(rowSum[i])/scale;
                if (rel > 1e-6) ++nViol;
                sumRel += rel;
                if (rel > maxRel) maxRel = rel;
            }

            Pout<< "RSDEBUG-LEVELROWSUM level=" << levels_.size()
                << " nCells=" << n
                << " nViol=" << nViol
                << " meanRel=" << (n > 0 ? sumRel/n : 0)
                << " maxRel=" << maxRel
                << endl;
        }

        levels_.append(levelPtr.ptr());
    }

    if (levels_.empty())
    {
        FatalErrorInFunction
            << "No coarse levels created: matrix already at or below"
               " nCellsInCoarsestLevel (" << nCellsInCoarsestLevel_ << ")."
               " Choose another solver or reduce nCellsInCoarsestLevel."
            << exit(FatalError);
    }

    entry* coarseEntry = controlDict_.findEntry
    (
        "coarsestSolver",
        keyType::LITERAL_RECURSIVE
    );

    if (coarseEntry && coarseEntry->isDict())
    {
        coarsestSolverPtr_ = lduMatrix::solver::New
        (
            "coarsestLevelCorr",
            matrixLevel(levels_.size()),
            interfaceBouCoeffsLevel(levels_.size()),
            interfaceIntCoeffs_,
            interfacesLevel(levels_.size()),
            coarseEntry->dict()
        );
    }
    else
    {
        dictionary coarseControls;
        const bool asymmetric = matrixLevel(levels_.size()).asymmetric();

        coarseControls.add
        (
            "solver",
            word(asymmetric ? "PBiCGStab" : "PCG")
        );
        coarseControls.add
        (
            "preconditioner",
            word(asymmetric ? "DILU" : "DIC")
        );
        coarseControls.add("tolerance", tolerance_);
        coarseControls.add("relTol", relTol_);

        coarsestSolverPtr_ = lduMatrix::solver::New
        (
            "coarsestLevelCorr",
            matrixLevel(levels_.size()),
            interfaceBouCoeffsLevel(levels_.size()),
            interfaceIntCoeffs_,
            interfacesLevel(levels_.size()),
            coarseControls
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * //

const Foam::lduMatrix& Foam::RugeStubenAMG::matrixLevel(const label leveli) const
{
    return (leveli == 0 ? matrix_ : levels_[leveli - 1].coarseMatrix());
}


const Foam::lduInterfaceFieldPtrsList& Foam::RugeStubenAMG::interfacesLevel
(
    const label leveli
) const
{
    return (leveli == 0 ? interfaces_ : levels_[leveli - 1].coarseInterfaces());
}


const Foam::FieldField<Foam::Field, Foam::scalar>&
Foam::RugeStubenAMG::interfaceBouCoeffsLevel
(
    const label leveli
) const
{
    return
    (
        leveli == 0
      ? interfaceBouCoeffs_
      : levels_[leveli - 1].coarseInterfaceBouCoeffs()
    );
}


void Foam::RugeStubenAMG::smoothLevel
(
    const label leveli,
    scalarField& psi,
    const scalarField& source,
    const label nSweeps
) const
{
    const lduMatrix& A = matrixLevel(leveli);

    GaussSeidelSmoother smoother
    (
        fieldName_,
        A,
        interfaceBouCoeffsLevel(leveli),
        interfaceIntCoeffs_,
        interfacesLevel(leveli),
        dictionary::null
    );

    smoother.smooth(psi, source, 0, nSweeps);
}


void Foam::RugeStubenAMG::residualLevel
(
    const label leveli,
    const scalarField& psi,
    const scalarField& source,
    scalarField& res
) const
{
    const lduMatrix& A = matrixLevel(leveli);
    res.setSize(psi.size());

    A.residual
    (
        res, psi, source, interfaceBouCoeffsLevel(leveli), interfacesLevel(leveli), 0
    );
}


void Foam::RugeStubenAMG::Vcycle
(
    const label leveli,
    scalarField& psi,
    const scalarField& source
) const
{
    if (leveli == levels_.size())
    {
        coarsestSolverPtr_->scalarSolve(psi, source);
        return;
    }

    // Pre-smoothing
    smoothLevel(leveli, psi, source, nPreSweeps_);

    // Fine-level residual
    scalarField res;
    residualLevel(leveli, psi, source, res);

    // Restrict residual to the coarse level
    const RugeStubenCoarsening& lvl = levels_[leveli];

    scalarField coarseSource;
    lvl.restrictField(coarseSource, res, interfacesLevel(leveli));

    scalarField coarseCorr(lvl.nCoarseCells(), Zero);

    // Recurse to the next coarser level
    Vcycle(leveli + 1, coarseCorr, coarseSource);

    // Prolong and apply the coarse-grid correction
    scalarField fineCorr;
    lvl.prolongField(fineCorr, coarseCorr, interfacesLevel(leveli));

    psi += fineCorr;

    // Post-smoothing
    smoothLevel(leveli, psi, source, nPostSweeps_);
}




// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RugeStubenAMG::readControls()
{
    lduMatrix::solver::readControls();

    strengthThreshold_ =
        controlDict_.getOrDefault<scalar>("strengthThreshold", 0.25);

    nCellsInCoarsestLevel_ =
        controlDict_.getOrDefault<label>("nCellsInCoarsestLevel", 1);

    nPreSweeps_ = controlDict_.getOrDefault<label>("nPreSweeps", 2);
    nPostSweeps_ = controlDict_.getOrDefault<label>("nPostSweeps", 2);
}


Foam::solverPerformance Foam::RugeStubenAMG::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{
    solverPerformance solverPerf(typeName, fieldName_);

    scalarField Apsi(psi.size());
    scalarField tmpField(psi.size());

    matrix_.Amul(Apsi, psi, interfaceBouCoeffs_, interfaces_, cmpt);

    scalar normFactor = this->normFactor(psi, source, Apsi, tmpField);

    scalarField rA(psi.size());
    matrix_.residual(rA, psi, source, interfaceBouCoeffs_, interfaces_, cmpt);

    solverPerf.initialResidual() = gSumMag(rA, matrix_.mesh().comm())/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    if
    (
        minIter_ > 0
     || !solverPerf.checkConvergence(tolerance_, relTol_, log_)
    )
    {
        do
        {
            Vcycle(0, psi, source);

            matrix_.residual
            (
                rA, psi, source, interfaceBouCoeffs_, interfaces_, cmpt
            );

            solverPerf.finalResidual() =
                gSumMag(rA, matrix_.mesh().comm())/normFactor;

            ++solverPerf.nIterations();
        }
        while
        (
            solverPerf.nIterations() < maxIter_
         && !(
                solverPerf.checkConvergence(tolerance_, relTol_, log_)
             && solverPerf.nIterations() >= minIter_
            )
        );
    }

    return solverPerf;
}


// ************************************************************************* //
