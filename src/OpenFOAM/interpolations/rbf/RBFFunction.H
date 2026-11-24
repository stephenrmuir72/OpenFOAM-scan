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
    Foam::RBFFunction

Description
    Radial basis function virtual base class

    Run-time selectable interface for the RBF kernels.
    Used to compute interpolation weights from point distances.

SourceFiles
    RBFFunction.cxx

Namespace
    Foam::RBFTypes

Description
    Namespace for RBF functions etc.

\*---------------------------------------------------------------------------*/

#ifndef Foam_RBFFunction_H
#define Foam_RBFFunction_H

#include "error.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "tmp.H"
#include "point.H"
#include "Field.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;

/*---------------------------------------------------------------------------*\
                         Class RBFFunction Declaration
\*---------------------------------------------------------------------------*/

class RBFFunction
{
protected:

    // Protected Data

        //- Radius (unused by the weighting functions)
        scalar radius_;


    // Protected Member Functions

        //- Check field sizes or report FatalError
        template<class Type1, class Type2>
        static void checkFieldSizes
        (
            const UList<Type1>& wts,
            const UList<Type2>& pts,
            const char* funcName
        )
        {
            if (wts.size() != pts.size())
            {
                FatalErrorIn(funcName)
                    << "Size mismatch: weights " << wts.size()
                    << " points " << pts.size()
                    << exit(FatalError);
            }
        }

public:

    //- Runtime type information (no debug)
    TypeNameNoDebug("RBFFunction");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            RBFFunction,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Selectors

        //- Return a pointer to the selected RBF function
        static autoPtr<RBFFunction> New
        (
            const word& RBFname,
            const dictionary& dict
        );


    // Constructors

        //- Default construct (unit radius) or with specified radius
        explicit constexpr RBFFunction(scalar radius = 1) noexcept
        :
            radius_(radius)
        {}

        //- Create and return a clone
        virtual autoPtr<RBFFunction> clone() const = 0;


    //- Destructor
    virtual ~RBFFunction() = default;


    // Properties

        //- Return the filter radius [non-virtual]
        scalar radius() const noexcept { return radius_; }

        //- Set the filter radius, return previous value [non-virtual]
        scalar radius(scalar r) noexcept
        {
            scalar old(radius_);
            radius_ = r;
            return old;
        }

        //- Is the filter compact? (identically zero outside of the radius).
        //  Many basis functions are compact, so that is also the default.
        //  Non-compact functions should override this.
        virtual bool compact() const noexcept { return false; }

        //- For \em non-compact filters, the cutoff distance beyond which
        //- the filter weight is below the fraction \p tol of its peak value
        //
        //  Compact functions simply return the radius() and ignore the
        //  specified tolerance.
        //  Non-compact functions should override with their characteristics.
        virtual scalar cutoff(scalar tol) const { return radius_; }


    // Member Functions

        //- Calculate the weights given the points.
        //- The input \p points and output \p wts sizes must match.
        virtual void weights
        (
            UList<scalar>& wts,
            const UList<point>& points,
            const point& refPoint
        ) const = 0;

        //- Return RBF weights
        virtual tmp<Field<scalar>> weights
        (
            const UList<point>& points,
            const point& refPoint
        ) const
        {
            auto tresult = tmp<Field<scalar>>::New(points.size());
            weights(tresult.ref(), points, refPoint);
            return tresult;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
