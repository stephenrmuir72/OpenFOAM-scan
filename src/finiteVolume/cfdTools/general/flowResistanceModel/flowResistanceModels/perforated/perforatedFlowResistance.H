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
    Foam::flowResistanceModels::Perforated

Description
    Reynolds-number-dependent flow resistance for a perforated plate
    (regular array of circular or square holes).

    \verbatim
    flowResistanceCoeffs
    {
        resistanceType      perforated;

        perforatedCoeffs
        {
            holeShape           circular;
            holeSize            0.005;
            openAreaFraction    0.5;
        }
    }
    \endverbatim

    or, with the open area fraction derived from the hole geometry
    and pitches:
    \verbatim
    flowResistanceCoeffs
    {
        resistanceType      perforated;

        perforatedCoeffs
        {
            holeShape           circular;
            holeSize            0.005;
            longitudinalPitch   0.01;
            transversePitch     0.01;
        }
    }
    \endverbatim

SourceFiles
    perforatedFlowResistance.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_flowResistanceModels_Perforated_H
#define Foam_flowResistanceModels_Perforated_H

#include "holeFlowResistance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace flowResistanceModels
{

/*---------------------------------------------------------------------------*\
                         Class Perforated Declaration
\*---------------------------------------------------------------------------*/

class Perforated
:
    public holeFlowResistance
{
    // Private Data

        //- Hole shape: circular or square
        word holeShape_;

        //- Hole size [m] (diameter or side)
        scalar holeSize_;


protected:

    // Protected Member Functions

        //- Area of a single hole [m^2]
        virtual scalar holeArea() const;

        //- Perimeter of a single hole [m]
        virtual scalar holePerimeter() const;

        //- Write the "perforatedCoeffs" sub-dict
        virtual void writeCoeffs(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("perforated");


    // Constructors

        //- Construct from dictionary (the flowResistanceCoeffs sub-dict)
        explicit Perforated(const dictionary& dict);

        //- Copy constructor
        Perforated(const Perforated&) = default;

        //- Construct and return a clone
        virtual autoPtr<flowResistanceModel> clone() const
        {
            return autoPtr<flowResistanceModel>(new Perforated(*this));
        }


    //- Destructor
    virtual ~Perforated() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace flowResistanceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
