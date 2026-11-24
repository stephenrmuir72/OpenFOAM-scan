/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::RBFTypes::WendlandC6

Description
    Wendland C6 compactly supported radial basis function:

    \verbatim
    (1 - q)^8 (32q^3 + 25q^2 + 8q +1)

    for q = (d/radius) when (d < radius), zero otherwise
    \endverbatim

    where \c d is the radial distance and \c radius is the support radius.

    Available dictionary entries:
    \table
      Property  | Description          | Type     | Reqd | Deflt
      radius    | The filter radius    | scalar   | yes  | -
    \endtable

SourceFiles
    RBF_WendlandC6.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_RBFTypes_WendlandC6_H
#define Foam_RBFTypes_WendlandC6_H

#include "RBFFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBFTypes
{

/*---------------------------------------------------------------------------*\
                    Class RBFTypes::WendlandC6 Declaration
\*---------------------------------------------------------------------------*/

class WendlandC6
:
    public RBFFunction
{
    // Private Member Functions

        //- Single point weight [non-virtual]
        inline scalar weight(const point& p, const point& p0) const
        {
            if (const scalar q = p.dist(p0)/radius_; q < 1)
            {
                return
                (
                    Foam::pow4(Foam::sqr(1 - q))   // ie, (1-q)^8
                  * (32*Foam::pow3(q) + 25*Foam::sqr(q) + 8*q + 1)
                );
            }
            else
            {
                return 0;
            }
        }

public:

    //- Runtime type information (no debug)
    TypeNameNoDebug("W6");


    // Constructors

        //- Default construct (unit radius) or with specified radius
        explicit WendlandC6(scalar radius = 1) noexcept
        :
            RBFFunction(radius)
        {}

        //- Construct from dictionary
        explicit WendlandC6(const dictionary& dict);

        //- Return a clone
        virtual autoPtr<RBFFunction> clone() const
        {
            return autoPtr<RBFFunction>(new WendlandC6(this->radius()));
        }


    // Member Functions

        //- The filter is compact (identically zero outside of the radius)
        virtual bool compact() const noexcept { return true; }


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
