/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2009 Frank Bos, TU Delft.
    Copyright (C) 2009 Dubravko Matijasevic, FSB Zagreb
    Copyright (C) 2013-2020 Darrin Stephens, Applied CCM Pty Ltd
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

Class
    Foam::RBFInterpolation

Description
    Interpolation class which uses Radial Basis Functions (RBF) to interpolate
    the fluid displacements for given boundary displacements.

    The coefficient vectors, alpha and beta, are determined by solving
    the system:

    \verbatim
        | db | = | Mbb Pb | | alpha |
        | 0  |   | Pb^T 0 | | beta  |
    \endverbatim

    where
    \vartable
        db    | Given boundary displacements
        Mbb   | Boundary RBF correlation matrix (Nb x Nb), containing RBF
              | evaluations at the boundary nodes
        Pb    | Linear polynomial matrix (Nb x 4)
    \endvartable

    Those coefficients are calculated every timestep, with the current
    boundary displacements db, with the inverse of Mbb. Using those
    coefficients, the RBF is evaluated at all fluid points every
    timestep.

    The efficiency of this method is increased by:
    -# using control points, which are a subset of the moving boundary
       points. Those control points are selected by a coarsening function.
    -# neglecting the outer boundary points, since a cut-off function
       is used toward the outer boundaries.

Usage
    Example of the dictionary specification:
    \verbatim
    {
        RBF             <word>;
        radius          <scalar>;
        focalPoint      <point>;
        innerDistance   <scalar>;
        outerDistance   <scalar>;
        polynomials     <bool>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property    | Description                         | Type   | Reqd | Deflt
      RBF         | RBF type name (eg, WendlandC2..)    | word   | yes  | -
      radius      | Radius (support/shape) of the RBF   | scalar | yes  | -
      focalPoint  | Focal point of the cut-off radii    | point  | yes  | -
      innerDistance | Inner cut-off distance            | scalar | yes  | -
      outerDistance | Outer cut-off distance            | scalar | yes  | -
      polynomials | Add linear polynomial terms         | bool   | no   | false
    \endtable

SourceFiles
    RBFInterpolation.cxx
    RBFInterpolation.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_RBFInterpolation_H
#define Foam_RBFInterpolation_H

#include "RBFFunction.H"
#include "pointField.H"
#include "IndirectList.H"
#include "SquareMatrix.H"
#include "SubList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;

/*---------------------------------------------------------------------------*\
                      Class RBFInterpolation Declaration
\*---------------------------------------------------------------------------*/

class RBFInterpolation
{
public:

    // Forward Declarations
    class InterpolationMatrix;

private:

    // Private Data

        //- Reference to control points (all ranks)
        const pointField& controlPoints_;

        //- Reference to the interpolation points (local rank only)
        const pointField& interpPoints_;

        //- RBF function
        autoPtr<RBFFunction> RBFptr_;

        //- Interpolation matrix
        mutable std::unique_ptr<InterpolationMatrix> matrixPtr_;

        //- Focal point for morphing region
        point focalPoint_;

        //- Inner morphing distance
        scalar innerDist_;

        //- Outer morphing distance
        scalar outerDist_;

        //- Add polynomials to RBF matrix
        bool withPolynomials_;


    // Private Member Functions

        //- Calculate the interpolation matrix (matrixPtr_)
        void calcInterpolationMatrix() const;

        //- Return the interpolation matrix (demand-driven)
        const InterpolationMatrix& interpMatrix() const
        {
            if (!matrixPtr_) calcInterpolationMatrix();
            return *matrixPtr_;
        }

        //- Interpolate control field
        template<class Type, class CtrlField>
        tmp<Field<Type>> interpolate_impl(const CtrlField& ctrlField) const;


public:

    // Generated Methods

        //- No copy assignment
        void operator=(const RBFInterpolation&) = delete;


    // Constructors

        //- Copy construct
        RBFInterpolation(const RBFInterpolation&);

        //- Construct from components
        RBFInterpolation
        (
            const dictionary& dict,
            //! Reference to control points (multi-rank)
            const pointField& controlPoints,
            //! Reference to the interpolation points (local rank only)
            const pointField& interpPoints
        );


    //- Destructor
    ~RBFInterpolation();

        //- Clear out interpolation matrix
        void clearOut();


    // Member Functions

        //- The RBF function
        const RBFFunction& RBF() const { return RBFptr_(); }

        //- Interpolate the control field,
        //- which has a size correspondence to the controlPoints
        template<class Type>
        tmp<Field<Type>> interpolate(const Field<Type>& ctrlField) const;

        //- Interpolate the control field,
        //- which has a size correspondence to the controlPoints
        template<class Type>
        tmp<Field<Type>>
        interpolate(const UIndirectList<Type>& ctrlField) const;
};


/*---------------------------------------------------------------------------*\
            Class RBFInterpolation::InterpolationMatrix Declaration
\*---------------------------------------------------------------------------*/

//- A simple symmetric square matrix (scalar) with packed storage
//- of the upper triangle only (row-wise) and minimal functionality.
class RBFInterpolation::InterpolationMatrix
{
    // Private Data

        //- Number of rows/columns
        label nCols_;

        //- Packed upper-triangle values, size n*(n+1)/2
        List<scalar> v_;


    // Private Member Functions

        //- Storage index for (i,j), using symmetry
        label index(label i, label j) const noexcept
        {
            if (i > j) std::swap(i, j);
            return (i*nCols_ - (i*(i - 1))/2 + (j - i));
        }

public:

    // Constructors

        //- Default construct (empty matrix)
        InterpolationMatrix() noexcept : nCols_(0) {}

        //- Construct empty (uninitialised) storage for an n x n matrix
        explicit InterpolationMatrix(label n)
        :
            nCols_(n),
            v_((n*(n + 1))/2)
        {}

        //- Copy construct from the upper-triangle of a square matrix
        explicit InterpolationMatrix(const SquareMatrix<scalar>& mat)
        :
            nCols_(mat.n()),
            v_((nCols_*(nCols_ + 1))/2)
        {
            label k = 0;
            for (label i = 0; i < nCols_; ++i)
            {
                for (label j = i; j < nCols_; ++j)
                {
                    v_[k++] = mat(i, j);
                }
            }
        }


    // Member Functions

        //- Number of rows/columns
        label n() const noexcept { return nCols_; }

        //- Number of stored (upper-triangle) values
        label size() const noexcept { return v_.size(); }

        //- Empty matrix?
        bool empty() const noexcept { return v_.empty(); }

        //- Const access to the packed storage
        const scalar* cdata() const noexcept { return v_.cdata(); }

        //- Access to the packed storage
        scalar* data() noexcept { return v_.data(); }


    // Row-based access

        //- Const pointer to start of the row data in the packed storage
        const scalar* rowData(label irow) const noexcept
        {
            return (v_.cdata() + index(irow, irow));
        }

        //- Pointer to start of the row data in the packed storage
        scalar* rowData(label irow) noexcept
        {
            return (v_.data() + index(irow, irow));
        }


        //- Offset in the packed storage of the start of the given row.
        //  Row i holds only the trailing portion (columns i..n-1).
        label rowStart(label irow) const { return index(irow, irow); }

        //- Number of stored values in the trailing portion of the row,
        //- i.e. columns irow..n-1
        label rowSize(label irow) const noexcept { return (nCols_ - irow); }

        //- Trailing portion of the row (columns irow..n-1),
        //- no subscript checking
        const SubList<scalar> rowTail(label irow) const
        {
            return v_.slice(rowStart(irow), rowSize(irow));
        }

        //- Trailing portion of the row (columns irow..n-1),
        //- no subscript checking
        SubList<scalar> rowTail(label irow)
        {
            return v_.slice(rowStart(irow), rowSize(irow));
        }


    // Member Operators

        //- The (i,j) value; symmetric, so (j,i) is the same storage
        scalar& operator()(label i, label j)
        {
            return v_[index(i, j)];
        }

        //- The (i,j) value; symmetric, so (j,i) is the same storage
        const scalar& operator()(label i, label j) const
        {
            return v_[index(i, j)];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "RBFInterpolation.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
