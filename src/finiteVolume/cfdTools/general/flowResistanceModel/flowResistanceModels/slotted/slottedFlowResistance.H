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
    Foam::flowResistanceModels::Slotted

Description
    Flow resistance for a slotted opening (regular array of slots)
    using the fully turbulent loss coefficient with the slot shape
    factor k_shape. The loss coefficient is independent of the
    Reynolds number.

    \verbatim
    flowResistanceCoeffs
    {
        resistanceType      slotted;

        slottedCoeffs
        {
            slotLength          0.01;
            slotWidth           0.005;
            longitudinalPitch   0.02;
            transversePitch     0.015;
        }
    }
    \endverbatim

SourceFiles
    slottedFlowResistance.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_flowResistanceModels_Slotted_H
#define Foam_flowResistanceModels_Slotted_H

#include "holeFlowResistance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace flowResistanceModels
{

/*---------------------------------------------------------------------------*\
                          Class Slotted Declaration
\*---------------------------------------------------------------------------*/

class Slotted
:
    public holeFlowResistance
{
    // Private Data

        //- Slot length [m] (a: straight length + cap diameter)
        scalar slotLength_;

        //- Slot width [m] (b: cap diameter = slot width)
        scalar slotWidth_;


protected:

    // Protected Member Functions

        //- Area of a single slot [m^2]
        virtual scalar holeArea() const;

        //- Slot shape factor k_shape (Eq. 16)
        virtual scalar shapeFactor() const;

        //- Write the "slottedCoeffs" sub-dict
        virtual void writeCoeffs(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("slotted");


    // Constructors

        //- Construct from dictionary (the flowResistanceCoeffs sub-dict)
        explicit Slotted(const dictionary& dict);

        //- Copy constructor
        Slotted(const Slotted&) = default;

        //- Construct and return a clone
        virtual autoPtr<flowResistanceModel> clone() const
        {
            return autoPtr<flowResistanceModel>(new Slotted(*this));
        }


    //- Destructor
    virtual ~Slotted() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace flowResistanceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
