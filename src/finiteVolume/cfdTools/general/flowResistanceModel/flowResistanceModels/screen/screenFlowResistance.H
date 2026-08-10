/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::flowResistanceModels::Screen

Description
    Placeholder for a screen (mesh/fabric) flow resistance. Not yet
    implemented: always returns a zero loss coefficient. Any
    user-supplied \c screenCoeffs sub-dict is preserved verbatim so
    that it round-trips through write() without being silently
    dropped, e.g. on a restart write of the boundary condition.

    \verbatim
    flowResistanceCoeffs
    {
        resistanceType      screen;

        screenCoeffs
        {
            // Not yet implemented; entries preserved verbatim
        }
    }
    \endverbatim

SourceFiles
    screenFlowResistance.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_flowResistanceModels_Screen_H
#define Foam_flowResistanceModels_Screen_H

#include "flowResistanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace flowResistanceModels
{

/*---------------------------------------------------------------------------*\
                          Class Screen Declaration
\*---------------------------------------------------------------------------*/

class Screen
:
    public flowResistanceModel
{
    // Private Data

        //- Open area fraction f (0..1), if specified directly.
        //  Negative = not set (only the screenCoeffs sub-dict form
        //  was used, or neither was provided)
        scalar openAreaFraction_;

        //- Verbatim copy of the "screenCoeffs" sub-dict, so that any
        //  user-supplied content (currently unused, since screen is
        //  not yet implemented) is preserved on write() rather than
        //  silently dropped
        dictionary placeholderCoeffsDict_;

        //- True if the screenCoeffs sub-dict was used (rather than a
        //  direct openAreaFraction entry)
        bool geometryUsed_;


protected:

    // Protected Member Functions

        //- Not yet implemented: always returns 0
        virtual scalar lossCoefficient
        (
            const scalar UnAvg,
            const scalar nuAvg
        ) const;

        //- Write the "screenCoeffs" sub-dict verbatim, or the
        //  "openAreaFraction" entry
        virtual void writeCoeffs(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("screen");


    // Constructors

        //- Construct from dictionary (the flowResistanceCoeffs sub-dict)
        explicit Screen(const dictionary& dict);

        //- Copy constructor
        Screen(const Screen&) = default;

        //- Construct and return a clone
        virtual autoPtr<flowResistanceModel> clone() const
        {
            return autoPtr<flowResistanceModel>(new Screen(*this));
        }


    //- Destructor
    virtual ~Screen() = default;


    // Member Functions

        //- The screen resistance is not yet implemented and so never
        //  requires the kinematic viscosity
        virtual bool needsNu() const
        {
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace flowResistanceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
