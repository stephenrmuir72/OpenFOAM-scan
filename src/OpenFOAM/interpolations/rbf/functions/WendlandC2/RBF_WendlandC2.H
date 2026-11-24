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
    Foam::RBFTypes::WendlandC2

Description
    Wendland C2 compactly supported radial basis function:

    \verbatim
    (1 - q)^4 (1 + 4 q)

    for q = (d/radius) when (d < radius), zero otherwise
    \endverbatim

    where \c d is the radial distance and \c radius is the support radius.

    Available dictionary entries:
    \table
      Property  | Description          | Type     | Reqd | Deflt
      radius    | The filter radius    | scalar   | yes  | -
    \endtable

SourceFiles
    RBF_WendlandC2.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_RBFTypes_WendlandC2_H
#define Foam_RBFTypes_WendlandC2_H

#include "RBFFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBFTypes
{

/*---------------------------------------------------------------------------*\
                    Class RBFTypes::WendlandC2 Declaration
\*---------------------------------------------------------------------------*/

class WendlandC2
:
    public RBFFunction
{
    // Private Member Functions

        //- Single point weight [non-virtual]
        inline scalar weight(const point& p, const point& p0) const
        {
            if (const scalar q = p.dist(p0)/radius_; q < 1)
            {
                return Foam::pow4(1 - q)*(4*q + 1);
            }
            else
            {
                return 0;
            }
        }

public:

    //- Runtime type information (no debug)
    TypeNameNoDebug("WendlandC2");


    // Constructors

        //- Default construct (unit radius) or with specified radius
        explicit WendlandC2(scalar radius = 1) noexcept
        :
            RBFFunction(radius)
        {}

        //- Construct from dictionary
        explicit WendlandC2(const dictionary& dict);

        //- Return a clone
        virtual autoPtr<RBFFunction> clone() const
        {
            return autoPtr<RBFFunction>(new WendlandC2(this->radius()));
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
