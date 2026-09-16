/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2021,2023 OpenCFD Ltd.
    Copyright (C) 2023 Huawei (Yu Ankun)
    Copyright (C) 2023 OpenCFD Ltd.
    Copyright (C) 2026 Keysight Technologies
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

#include "GAMGSolver.H"
#include "SubField.H"
#include "PrecisionAdaptor.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::GAMGSolver::solve
(
    scalarField& psi_s,
    const scalarField& source,
    const direction cmpt
) const
{
    // Agglomerate current matrix
    forAll(agglomeration_, fineLevelIndex)
    {
        // Agglomerate on to coarse level mesh
        const_cast<GAMGSolver*>(this)->agglomerateMatrix
        (
            fineLevelIndex,
            agglomeration_.meshLevel(fineLevelIndex + 1),
            agglomeration_.interfaceLevel(fineLevelIndex + 1)
        );
    }

    PrecisionAdaptor<solveScalar, scalar> tpsi(psi_s);
    solveScalarField& psi = tpsi.ref();

    ConstPrecisionAdaptor<solveScalar, scalar> tsource(source);

    // Setup class containing solver performance data
    solverPerformance solverPerf(typeName, fieldName_);

    // Calculate A.psi used to calculate the initial residual
    solveScalarField Apsi(psi.size());
    matrix_.Amul(Apsi, psi, interfaceBouCoeffs_, interfaces_, cmpt);

    // Create the storage for the finestCorrection which may be used as a
    // temporary in normFactor
    solveScalarField finestCorrection(psi.size());

    // Calculate normalisation factor
    solveScalar normFactor =
        this->normFactor(psi, tsource(), Apsi, finestCorrection);

    if ((log_ >= 2) || (debug >= 2))
    {
        Pout<< "   Normalisation factor = " << normFactor << endl;
    }

    // Calculate initial finest-grid residual field
    solveScalarField finestResidual(psi.size());
    matrix_.residual
    (
        finestResidual,
        psi,
        source,
        interfaceBouCoeffs_,
        interfaces_,
        cmpt
    );

    matrix().setResidualField
    (
        ConstPrecisionAdaptor<scalar, solveScalar>(finestResidual)(),
        fieldName_,
        true
    );

    // Calculate normalised residual for convergence test
    solverPerf.initialResidual() = gSumMag
    (
        finestResidual,
        matrix().mesh().comm()
    )/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();


    // Check convergence, solve if not converged
    if
    (
        minIter_ > 0
     || !solverPerf.checkConvergence(tolerance_, relTol_, log_)
    )
    {
        // Initialise coarse data structures
        if (coarseCorrFields_.size() != matrixLevels_.size())
        {
            initVcycle
            (
                coarseCorrFields_,
                coarseSources_,
                smoothers_,
                scratch1_,
                scratch2_
            );
        }

        do
        {
            Vcycle
            (
                smoothers_,
                psi,
                source,
                Apsi,
                finestCorrection,
                finestResidual,

                (scratch1_.size() ? scratch1_ : Apsi),
                (scratch2_.size() ? scratch2_ : finestCorrection),

                coarseCorrFields_,
                coarseSources_,
                cmpt
            );

            // Calculate finest level residual field
            matrix_.residual
            (
                finestResidual,
                psi,
                source,
                interfaceBouCoeffs_,
                interfaces_,
                cmpt
            );

            solverPerf.finalResidual() = gSumMag
            (
                finestResidual,
                matrix().mesh().comm()
            )/normFactor;

            if ((log_ >= 2) || (debug >= 2))
            {
                solverPerf.print(Info.masterStream(matrix().mesh().comm()));
            }
        } while
        (
            (
              ++solverPerf.nIterations() < maxIter_
            && !solverPerf.checkConvergence(tolerance_, relTol_, log_)
            )
         || solverPerf.nIterations() < minIter_
        );
    }

    matrix().setResidualField
    (
        ConstPrecisionAdaptor<scalar, solveScalar>(finestResidual)(),
        fieldName_,
        false
    );

    return solverPerf;
}


void Foam::GAMGSolver::Vcycle
(
    const PtrList<lduMatrix::smoother>& smoothers,
    solveScalarField& psi,
    const scalarField& source,
    solveScalarField& Apsi,
    solveScalarField& finestCorrection,
    solveScalarField& finestResidual,

    solveScalarField& scratch1,
    solveScalarField& scratch2,

    PtrList<solveScalarField>& coarseCorrFields,
    PtrList<solveScalarField>& coarseSources,
    const direction cmpt
) const
{
    //debug = 2;

    const label coarsestLevel = matrixLevels_.size() - 1;

    // Restrict finest grid residual for the next level up.
    agglomeration_.restrictField(coarseSources[0], finestResidual, 0, true);

    if (nPreSweeps_ && ((log_ >= 2) || (debug >= 2)))
    {
        Pout<< "Pre-smoothing scaling factors: ";
    }


    // Residual restriction (going to coarser levels)
    for (label leveli = 0; leveli < coarsestLevel; leveli++)
    {
        if (coarseSources_.set(leveli + 1))
        {
            // If the optional pre-smoothing sweeps are selected
            // smooth the coarse-grid field for the restricted source
            if (nPreSweeps_)
            {
                coarseCorrFields_[leveli] = 0.0;

                smoothers_[leveli + 1].scalarSmooth
                (
                    coarseCorrFields_[leveli],
                    coarseSources_[leveli],  //coarseSource,
                    cmpt,
                    min
                    (
                        nPreSweeps_ +  preSweepsLevelMultiplier_*leveli,
                        maxPreSweeps_
                    )
                );

                // Scale coarse-grid correction field
                // but not on the coarsest level because it evaluates to 1
                if (scaleCorrection_ && leveli < coarsestLevel - 1)
                {
                    solveScalarField::subField ACf
                    (
                        scratch1_,
                        coarseCorrFields_[leveli].size()
                    );

                    scale
                    (
                        coarseCorrFields_[leveli],
                        const_cast<solveScalarField&>
                        (
                            static_cast<const solveScalarField&>(ACf)
                        ),
                        matrixLevels_[leveli],
                        interfaceLevelsBouCoeffs_[leveli],
                        interfaceLevels_[leveli],
                        coarseSources_[leveli],
                        cmpt
                    );
                }

                // Correct the residual with the new solution
                // residual can be used by fusing Amul with b-Amul
                matrixLevels_[leveli].residual
                (
                    coarseSources_[leveli],
                    coarseCorrFields_[leveli],
                    ConstPrecisionAdaptor<scalar, solveScalar>
                    (
                        coarseSources_[leveli]
                    )(),
                    interfaceLevelsBouCoeffs_[leveli],
                    interfaceLevels_[leveli],
                    cmpt
                );
            }

            // Residual is equal to source
            agglomeration_.restrictField
            (
                coarseSources_[leveli + 1],
                coarseSources_[leveli],
                leveli + 1,
                true
            );
        }
    }

    if (nPreSweeps_ && ((log_ >= 2) || (debug >= 2)))
    {
        Pout<< endl;
    }


    // Solve Coarsest level with either an iterative or direct solver
    if (coarseCorrFields.set(coarsestLevel))
    {
        solveCoarsestLevel
        (
            coarseCorrFields_[coarsestLevel],
            coarseSources_[coarsestLevel]
        );
    }

    if ((log_ >= 2) || (debug >= 2))
    {
        Pout<< "Post-smoothing scaling factors: ";
    }

    // Smoothing and prolongation of the coarse correction fields
    // (going to finer levels)

    solveScalarField dummyField(0);

    // Work storage for prolongation
    solveScalarField work;

    for (label leveli = coarsestLevel - 1; leveli >= 0; leveli--)
    {
        if (coarseCorrFields.set(leveli))
        {
            // Create a field for the pre-smoothed correction field
            // as a sub-field of the finestCorrection which is not
            // currently being used
            solveScalarField::subField preSmoothedCoarseCorrField
            (
                scratch2_,
                coarseCorrFields_[leveli].size()
            );

            // Only store the preSmoothedCoarseCorrField if pre-smoothing is
            // used
            if (nPreSweeps_)
            {
                preSmoothedCoarseCorrField = coarseCorrFields_[leveli];
            }


            // Prolong correction to leveli
            const auto& cf = agglomeration_.prolongField
            (
                coarseCorrFields_[leveli],   // current level
                work,
                (
                    coarseCorrFields_.set(leveli + 1)
                  ? coarseCorrFields_[leveli + 1]
                  : dummyField              // dummy value
                ),
                leveli + 1
            );


            // Create A.psi for this coarse level as a sub-field of Apsi
            solveScalarField::subField ACf
            (
                scratch1_,
                coarseCorrFields_[leveli].size()
            );
            auto& ACfRef = const_cast<solveScalarField&>
            (
                static_cast<const solveScalarField&>(ACf)
            );

            if (interpolateCorrection_)
            {
                // Normal operation : have both coarse level and fine
                // level. No processor agglomeration
                interpolate
                (
                    coarseCorrFields_[leveli],
                    ACfRef,
                    matrixLevels_[leveli],
                    interfaceLevelsBouCoeffs_[leveli],
                    interfaceLevels_[leveli],
                    agglomeration_.restrictAddressing(leveli + 1),
                    cf,
                    cmpt
                );
            }

            // Scale coarse-grid correction field
            // but not on the coarsest level because it evaluates to 1
            if
            (
                scaleCorrection_
             && (interpolateCorrection_ || leveli < coarsestLevel - 1)
            )
            {
                scale
                (
                    coarseCorrFields_[leveli],
                    ACfRef,
                    matrixLevels_[leveli],
                    interfaceLevelsBouCoeffs_[leveli],
                    interfaceLevels_[leveli],
                    coarseSources_[leveli],
                    cmpt
                );
            }

            // Only add the preSmoothedCoarseCorrField if pre-smoothing is
            // used
            if (nPreSweeps_)
            {
                coarseCorrFields_[leveli] += preSmoothedCoarseCorrField;
            }

            smoothers[leveli + 1].scalarSmooth
            (
                coarseCorrFields_[leveli],
                coarseSources_[leveli],  //coarseSource,
                cmpt,
                min
                (
                    nPostSweeps_ + postSweepsLevelMultiplier_*leveli,
                    maxPostSweeps_
                )
            );
        }
    }

    // Prolong the finest level correction
    agglomeration_.prolongField
    (
        finestCorrection,
        coarseCorrFields_[0],
        0,
        true
    );

    if (interpolateCorrection_)
    {
        interpolate
        (
            finestCorrection,
            Apsi,
            matrix_,
            interfaceBouCoeffs_,
            interfaces_,
            agglomeration_.restrictAddressing(0),
            coarseCorrFields_[0],
            cmpt
        );
    }

    if (scaleCorrection_)
    {
        // Scale the finest level correction
        scale
        (
            finestCorrection,
            Apsi,
            matrix_,
            interfaceBouCoeffs_,
            interfaces_,
            finestResidual,
            cmpt
        );
    }

    forAll(psi, i)
    {
        psi[i] += finestCorrection[i];
    }

    smoothers[0].smooth
    (
        psi,
        source,
        cmpt,
        nFinestSweeps_
    );
}


void Foam::GAMGSolver::initVcycle
(
    PtrList<solveScalarField>& coarseCorrFields,
    PtrList<solveScalarField>& coarseSources,
    PtrList<lduMatrix::smoother>& smoothers,
    solveScalarField& scratch1,
    solveScalarField& scratch2
) const
{
    label maxSize = matrix_.diag().size();

    coarseCorrFields.setSize(matrixLevels_.size());
    coarseSources.setSize(matrixLevels_.size());
    smoothers.setSize(matrixLevels_.size() + 1);

    // Create the smoother for the finest level
    smoothers.set
    (
        0,
        lduMatrix::smoother::New
        (
            fieldName_,
            matrix_,
            interfaceBouCoeffs_,
            interfaceIntCoeffs_,
            interfaces_,
            controlDict_
        )
    );

    forAll(matrixLevels_, leveli)
    {
        if (agglomeration_.nCells(leveli) >= 0)
        {
            label nCoarseCells = agglomeration_.nCells(leveli);

            coarseSources.set(leveli, new solveScalarField(nCoarseCells));
        }

        if (matrixLevels_.set(leveli))
        {
            const lduMatrix& mat = matrixLevels_[leveli];

            label nCoarseCells = mat.diag().size();

            maxSize = max(maxSize, nCoarseCells);

            coarseCorrFields.set(leveli, new solveScalarField(nCoarseCells));

            smoothers.set
            (
                leveli + 1,
                lduMatrix::smoother::New
                (
                    fieldName_,
                    matrixLevels_[leveli],
                    interfaceLevelsBouCoeffs_[leveli],
                    interfaceLevelsIntCoeffs_[leveli],
                    interfaceLevels_[leveli],
                    controlDict_
                )
            );
        }
    }

    if (maxSize > matrix_.diag().size())
    {
        // Allocate some scratch storage
        scratch1.setSize(maxSize);
        scratch2.setSize(maxSize);
    }
    else
    {
        // Always allocate so SubField references inside Vcycle are valid
        scratch1.setSize(matrix_.diag().size());
        scratch2.setSize(matrix_.diag().size());
    }
}


Foam::dictionary Foam::GAMGSolver::PCGsolverDict
(
    const scalar tol,
    const scalar relTol
) const
{
    dictionary dict(IStringStream("solver PCG; preconditioner DIC;")());
    dict.add("tolerance", tol);
    dict.add("relTol", relTol);

    return dict;
}


Foam::dictionary Foam::GAMGSolver::PBiCGStabSolverDict
(
    const scalar tol,
    const scalar relTol
) const
{
    dictionary dict(IStringStream("solver PBiCGStab; preconditioner DILU;")());
    dict.add("tolerance", tol);
    dict.add("relTol", relTol);

    return dict;
}


void Foam::GAMGSolver::solveCoarsestLevel
(
    solveScalarField& coarsestCorrField,
    const solveScalarField& coarsestSource
) const
{
    const label coarsestLevel = matrixLevels_.size() - 1;

    const label coarseComm = matrixLevels_[coarsestLevel].mesh().comm();

    if (directSolveCoarsest_)
    {
        PrecisionAdaptor<scalar, solveScalar> tcorrField(coarsestCorrField);

        coarsestLUMatrixPtr_->solve
        (
            tcorrField.ref(),
            ConstPrecisionAdaptor<scalar, solveScalar>(coarsestSource)()
        );
    }
    //else if
    //(
    //    agglomeration_.processorAgglomerate()
    // && procMatrixLevels_.set(coarsestLevel)
    //)
    //{
    //    //const labelList& agglomProcIDs = agglomeration_.agglomProcIDs
    //    //(
    //    //    coarsestLevel
    //    //);
    //    //
    //    //scalarField allSource;
    //    //
    //    //globalIndex cellOffsets;
    //    //if (Pstream::myProcNo(coarseComm) == agglomProcIDs[0])
    //    //{
    //    //    cellOffsets.offsets() =
    //    //        agglomeration_.cellOffsets(coarsestLevel);
    //    //}
    //    //
    //    //cellOffsets.gather
    //    //(
    //    //    coarseComm,
    //    //    agglomProcIDs,
    //    //    coarsestSource,
    //    //    allSource
    //    //);
    //    //
    //    //scalarField allCorrField;
    //    //solverPerformance coarseSolverPerf;
    //
    //    label solveComm = agglomeration_.procCommunicator(coarsestLevel);
    //
    //    coarsestCorrField = 0;
    //    solverPerformance coarseSolverPerf;
    //
    //    if (Pstream::myProcNo(solveComm) != -1)
    //    {
    //        const lduMatrix& allMatrix = procMatrixLevels_[coarsestLevel];
    //
    //        {
    //            Pout<< "** Master:Solving on comm:" << solveComm
    //                << " with procs:" << UPstream::procID(solveComm) << endl;
    //
    //            if (allMatrix.asymmetric())
    //            {
    //                coarseSolverPerf = PBiCGStab
    //                (
    //                    "coarsestLevelCorr",
    //                    allMatrix,
    //                    procInterfaceLevelsBouCoeffs_[coarsestLevel],
    //                    procInterfaceLevelsIntCoeffs_[coarsestLevel],
    //                    procInterfaceLevels_[coarsestLevel],
    //                    PBiCGStabSolverDict(tolerance_, relTol_)
    //                ).solve
    //                (
    //                    coarsestCorrField,
    //                    coarsestSource
    //                );
    //            }
    //            else
    //            {
    //                coarseSolverPerf = PCG
    //                (
    //                    "coarsestLevelCorr",
    //                    allMatrix,
    //                    procInterfaceLevelsBouCoeffs_[coarsestLevel],
    //                    procInterfaceLevelsIntCoeffs_[coarsestLevel],
    //                    procInterfaceLevels_[coarsestLevel],
    //                    PCGsolverDict(tolerance_, relTol_)
    //                ).solve
    //                (
    //                    coarsestCorrField,
    //                    coarsestSource
    //                );
    //            }
    //        }
    //    }
    //
    //    Pout<< "done master solve." << endl;
    //
    //    //// Scatter to all processors
    //    //coarsestCorrField.setSize(coarsestSource.size());
    //    //cellOffsets.scatter
    //    //(
    //    //    coarseComm,
    //    //    agglomProcIDs,
    //    //    allCorrField,
    //    //    coarsestCorrField
    //    //);
    //
    //    if (debug >= 2)
    //    {
    //        coarseSolverPerf.print(Info.masterStream(coarseComm));
    //    }
    //
    //    Pout<< "procAgglom: coarsestSource   :" << coarsestSource << endl;
    //    Pout<< "procAgglom: coarsestCorrField:" << coarsestCorrField << endl;
    //}
    else
    {
        coarsestCorrField = 0;
        const solverPerformance coarseSolverPerf
        (
            coarsestSolverPtr_->scalarSolve
            (
                coarsestCorrField,
                coarsestSource
            )
        );

        if ((log_ >= 2) || debug)
        {
            coarseSolverPerf.print(Info.masterStream(coarseComm));
        }
    }
}


// ************************************************************************* //
