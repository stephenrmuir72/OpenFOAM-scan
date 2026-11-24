/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::RBFTypes::Cosine

Description
    Cosine radial basis function with the compact kernel

    \verbatim
    1/2 (1 + cos(pi q))

    for q = (d/radius) when (d < radius), zero otherwise
    \endverbatim

    where \c d is the radial distance and \c radius is the support radius.

    Available dictionary entries:
    \table
      Property  | Description          | Type     | Reqd | Deflt
      radius    | The filter radius    | scalar   | yes  | -
    \endtable

SourceFiles
    RBF_Cosine.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_RBFTypes_Cosine_H
#define Foam_RBFTypes_Cosine_H

#include "RBFFunction.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBFTypes
{

/*---------------------------------------------------------------------------*\
                            Class Cosine Declaration
\*---------------------------------------------------------------------------*/

class Cosine
:
    public RBFFunction
{
    // Private Member Functions

        //- Single point weight [non-virtual]
        inline scalar weight(const point& p, const point& p0) const
        {
            if (const scalar q = p.dist(p0)/radius_; q < 1)
            {
                return 0.5*(1 + std::cos(q * Foam::constant::mathematical::pi));
            }
            else
            {
                return 0;
            }
        }

public:

    //- Runtime type information (no debug)
    TypeNameNoDebug("cosine");


    // Constructors

        //- Default construct (unit radius) or with specified radius
        explicit Cosine(scalar radius = 1)
        :
            RBFFunction(radius)
        {}

        //- Construct from dictionary
        explicit Cosine(const dictionary& dict);

        //- Create and return a clone
        virtual autoPtr<RBFFunction> clone() const
        {
            return autoPtr<RBFFunction>(new Cosine(this->radius()));
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
