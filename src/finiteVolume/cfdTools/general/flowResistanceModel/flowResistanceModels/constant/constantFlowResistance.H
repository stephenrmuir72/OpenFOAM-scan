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
    Foam::flowResistanceModels::Constant

Description
    User-specified constant loss coefficient K, independent of the
    Reynolds number.

    \verbatim
    flowResistanceCoeffs
    {
        resistanceType      constant;
        lossCoefficient     2.5;
        relax               0.3;
    }
    \endverbatim

SourceFiles
    constantFlowResistance.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_flowResistanceModels_Constant_H
#define Foam_flowResistanceModels_Constant_H

#include "flowResistanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace flowResistanceModels
{

/*---------------------------------------------------------------------------*\
                          Class Constant Declaration
\*---------------------------------------------------------------------------*/

class Constant
:
    public flowResistanceModel
{
    // Private Data

        //- User-specified loss coefficient K
        scalar lossCoefficient_;


protected:

    // Protected Member Functions

        //- Return the user-specified loss coefficient
        virtual scalar lossCoefficient
        (
            const scalar UnAvg,
            const scalar nuAvg
        ) const;

        //- Write the "lossCoefficient" entry
        virtual void writeCoeffs(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from dictionary (the flowResistanceCoeffs sub-dict)
        explicit Constant(const dictionary& dict);

        //- Copy constructor
        Constant(const Constant&) = default;

        //- Construct and return a clone
        virtual autoPtr<flowResistanceModel> clone() const
        {
            return autoPtr<flowResistanceModel>(new Constant(*this));
        }


    //- Destructor
    virtual ~Constant() = default;


    // Member Functions

        //- A constant loss coefficient never requires the kinematic
        //  viscosity
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
