/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2009 Frank Bos, TU Delft
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
    Foam::RBFTypes::invGauss

Description
    Gauss radial basis function with the non-compact Gaussian kernel

    \verbatim
    exp(-(d r)^2)
    \endverbatim

    where \c d is the radial distance and \c radius is the width parameter
    (larger \c radius == \b narrower filter).

    Available dictionary entries:
    \table
      Property  | Description              | Type     | Reqd | Deflt
      radius    | The inverse filter width | scalar   | yes  | -
    \endtable

Note
    This formulation corresponds to the \em "Gauss" method from the
    extend-project. However here it it uses the "invGauss" name
    to reflect that the filter width is a multiplier and not normalization
    etc.

SourceFiles
    RBF_invGauss.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_RBFTypes_invGauss_H
#define Foam_RBFTypes_invGauss_H

#include "RBFFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBFTypes
{

/*---------------------------------------------------------------------------*\
                          Class invGauss Declaration
\*---------------------------------------------------------------------------*/

class invGauss
:
    public RBFFunction
{
    // Private Member Functions

        //- Single point weight [non-virtual]
        inline scalar weight(const point& p, const point& p0) const
        {
            return std::exp(-p.distSqr(p0)*Foam::sqr(radius_));
        }

        //- Threshold radius beyond which the weight is below the fraction
        //- \p tol of its peak value (1 at r=0).
        //  Inverse of weight(): solves exp(-r^2 R^2) = tol
        inline scalar threshold_radius(scalar tol) const
        {
            if (tol <= VSMALL) return GREAT;
            if (tol >= 1) return 0;
            return std::sqrt(-std::log(tol))/radius_;
        }

public:

    //- Runtime type information (no debug)
    TypeNameNoDebug("invGauss");


    // Constructors

        //- Default construct (unit radius) or with specified radius
        explicit invGauss(scalar radius = 1) noexcept
        :
            RBFFunction(radius)
        {}

        //- Construct from dictionary
        explicit invGauss(const dictionary& dict);

        //- Create and return a clone
        virtual autoPtr<RBFFunction> clone() const
        {
            return autoPtr<RBFFunction>(new invGauss(this->radius()));
        }


    // Member Functions

        //- The filter is non-compact
        virtual bool compact() const noexcept { return false; }

        //- Threshold radius beyond which the weight is below the fraction
        //- \p tol of its peak value (1 at r=0).
        virtual scalar cutoff(scalar tol) const
        {
            return threshold_radius(tol);
        }

        //- Calculate the weights given the points.
        //- The input \p points and output \p wts sizes must match.
        virtual void weights
        (
            UList<scalar>& wts,
            const UList<point>& points,
            const point& refPoint
        ) const;

        //- Inherit the two-parameter weights() from base
        using RBFFunction::weights;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBFTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
