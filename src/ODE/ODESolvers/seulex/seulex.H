/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2026 OpenFOAM Contributors
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
    Foam::seulex

Group
    grpODESolvers

Description
    An extrapolation-algorithm, based on the linearly implicit Euler method
    with step size control and order selection.

    With the optional switch \c sparse (on by default), the linear systems
    \c (1/dx)I-J are factorised by a sparse LU on the sparsity pattern of
    the Jacobian (seulexSparseLU) whenever the symbolic factorisation needs
    fewer than half the multiply-adds of the dense LU, falling back to the
    dense LU with partial pivoting if a pivot is too small.

    With the optional switch \c autonomous (off by default), for ODE systems
    whose derivatives do not depend on the time argument (such as chemistry),
    the initial right-hand side \c f(y0) of a step is evaluated once and
    shared across extrapolation columns, and the right-hand side evaluated by
    the divergence test is reused for the next substep.

    Reference:
    \verbatim
        Hairer, E., Nørsett, S. P., & Wanner, G. (1996).
        Solving Ordinary Differential Equations II:
        Stiff and Differential-Algebraic Problems, second edition",
        Springer-Verlag, Berlin.
    \endverbatim

SourceFiles
    seulex.C

\*---------------------------------------------------------------------------*/

#ifndef seulex_H
#define seulex_H

#include "ODESolver.H"
#include "scalarMatrices.H"
#include "labelField.H"
#include "boolList.H"
#include "seulexSparseLU.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class seulex Declaration
\*---------------------------------------------------------------------------*/

class seulex
:
    public ODESolver
{
    // Private data

        // Static constants

            static const label kMaxx_ = 12;
            static const label iMaxx_ = kMaxx_ + 1;

            static const scalar
                stepFactor1_, stepFactor2_, stepFactor3_,
                stepFactor4_, stepFactor5_,
                kFactor1_, kFactor2_;

        //- Allow sparse LU factorisation when the pattern is sparse (default: true)
        const bool useSparse_;

        //- The derivatives do not depend on the time argument: right-hand
        //  sides at the same state are shared (switch autonomous, default: false)
        const bool autonomous_;

        // Evaluated constants

            scalar jacRedo_;
            labelField nSeq_;
            scalarField cpu_;
            scalarSquareMatrix coeff_;

        // Temporary storage
        // held to avoid dynamic memory allocation between calls
        // and to transfer internal values between functions

            mutable scalar theta_;
            mutable label kTarg_;
            mutable scalarRectangularMatrix table_;

            mutable scalarField dfdx_;
            mutable scalarSquareMatrix dfdy_;
            mutable scalarSquareMatrix a_;
            mutable labelList pivotIndices_;

            // Fields space for "solve" function
            mutable scalarField dxOpt_, temp_;
            mutable scalarField y0_, ySequence_, scale_;

            // Fields used in "seul" function
            mutable scalarField dy_, yTemp_, dydx_;

            //- f(y0) of the current step attempt (autonomous systems)
            mutable scalarField f0_;

        // Sparse linear algebra

            //- Possible non-zeros of the Jacobian, n x n row-major
            mutable boolList pattern_;

            //- Ordering and symbolic factorisation for pattern_
            mutable seulexSparseLU slu_;

            //- The sparse factorisation is cheaper for this pattern
            mutable bool sparse_;

            //- The current factorisation is the dense one
            mutable bool denseNow_;

            //- Pattern rebuilds; past a limit (the system keeps changing
            //  size, as with TDAC) the solver stays dense
            mutable label nBuild_;


    // Private Member Functions

        //- Add the Jacobian's non-zeros to the pattern; rebuild the
        //  symbolic factorisation when it grew
        void updatePattern(const scalar x, const scalarField& y) const;

        //- Factorise (1/dx) I - J
        void factorise(const scalar dx) const;

        //- Solve with the current factorisation
        void linSolve(scalarField& b) const;

        //- Computes the j-th line of the extrapolation table
        bool seul
        (
            const scalar x0,
            const scalarField& y0,
            const scalar dxTot,
            const label k,
            scalarField& y,
            const scalarField& scale
        ) const;

        //- Polynomial extrapolation
        void extrapolate
        (
            const label k,
            scalarRectangularMatrix& table,
            scalarField& y
        ) const;


public:

    //- Runtime type information
    TypeName("seulex");


    // Constructors

        //- Construct from ODESystem
        seulex(const ODESystem& ode, const dictionary& dict);


    //- Destructor
    virtual ~seulex() = default;


    // Member Functions

        //- Resize the ODE solver
        virtual bool resize();

        //- Solve the ODE system and the update the state
        virtual void solve
        (
            scalar& x,
            scalarField& y,
            stepState& step
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
