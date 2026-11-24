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
    Foam::RBFTypes::ThinPlateSpline

Description
    Thin plate spline radial basis function.

    Thin plate spline kernel, d^2 log(d), where d is the distance.
    Used for RBF interpolation weights; no radius parameter.

    Available dictionary entries:
    \table
      Property  | Description          | Type     | Reqd | Deflt
      radius    | The filter radius    | scalar   | yes  | -
    \endtable

SourceFiles
    RBF_ThinPlateSpline.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_RBFTypes_ThinPlateSpline_H
#define Foam_RBFTypes_ThinPlateSpline_H

#include "RBFFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBFTypes
{

/*---------------------------------------------------------------------------*\
                  Class RBFTypes::ThinPlateSpline Declaration
\*---------------------------------------------------------------------------*/

class ThinPlateSpline
:
    public RBFFunction
{
    // Private Member Functions

        //- Single point weight.
        //  Note: it does not actually use radius at all
        inline scalar weight(const point& p, const point& p0) const
        {
            scalar w(0);
            if (scalar dist = p.dist(p0); dist > SMALL)
            {
                w = Foam::sqr(dist)*std::log(dist);
            }
            return w;
        }

public:

    //- Runtime type information (no debug)
    TypeNameNoDebug("TPS");


    // Constructors

        //- Default construct (unit radius) or with specified radius
        explicit ThinPlateSpline(scalar radius = 1) noexcept
        :
            RBFFunction(radius)
        {}

        //- Construct from dictionary
        ThinPlateSpline(const dictionary& dict);

        //- Return a clone
        virtual autoPtr<RBFFunction> clone() const
        {
            return autoPtr<RBFFunction>(new ThinPlateSpline(this->radius()));
        }


    // Member Functions

        //- The filter is non-compact
        virtual bool compact() const noexcept { return false; }

        //- The filter weight grows with distance (instead of decaying)
        //- so a cutoff distance is not meaningful. Returns GREAT.
        virtual scalar cutoff(scalar tol) const { return GREAT; }


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
