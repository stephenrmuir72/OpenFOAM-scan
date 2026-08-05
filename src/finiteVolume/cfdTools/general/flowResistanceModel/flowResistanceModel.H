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
    Foam::flowResistanceModel

Description
    Computes the loss coefficient K for planar flow resistance elements
    (e.g., perforated plates, slotted openings, wire screens) using
    geometry-based correlations or a user-specified constant.

    The pressure drop across the resistance is:

    \f[
        \Delta p = K \frac{1}{2} \rho |U_n|^2
    \f]

    where \c K is the loss coefficient and \c U_n is the face-normal
    velocity.

    Supported resistance types:
    - \c perforated: K computed from open area fraction using
      \f$ K = \frac{[0.707(1-f)]^2}{f^2} \f$
    - \c slotted: placeholder (returns 0)
    - \c screen: placeholder (returns 0)
    - \c constant: user-specified K via \c lossCoefficient entry

    The open area fraction f can be specified directly or derived from
    hole geometry (shape, size, and pitch).

Usage
    Example with direct open area fraction:
    \verbatim
    flowResistanceCoeffs
    {
        resistanceType      perforated;
        openAreaFraction    0.5;
    }
    \endverbatim

    Example with geometry-derived open area fraction:
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

    Example with constant K:
    \verbatim
    flowResistanceCoeffs
    {
        resistanceType      constant;
        lossCoefficient     2.5;
    }
    \endverbatim

SourceFiles
    flowResistanceModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_flowResistanceModel_H
#define Foam_flowResistanceModel_H

#include "dictionary.H"
#include "autoPtr.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class flowResistanceModel Declaration
\*---------------------------------------------------------------------------*/

class flowResistanceModel
{
    // Private Data

        //- Type of flow resistance geometry
        word resistanceType_;

        //- Open area fraction f (0..1). Negative = not set.
        scalar openAreaFraction_;

        //- User-specified loss coefficient K. Negative = not set.
        scalar lossCoefficient_;

        //- Hole shape for perforated geometry: circular or square
        word holeShape_;

        //- Hole size for perforated geometry [m] (diameter or side)
        scalar holeSize_;

        //- Longitudinal pitch for perforated geometry [m] (S_l)
        scalar longitudinalPitch_;

        //- Transverse pitch for perforated geometry [m] (S_t)
        scalar transversePitch_;

        //- True if openAreaFraction was derived from geometry inputs
        bool geometryUsed_;


public:

    // Constructors

        //- Construct from dictionary (the flowResistanceCoeffs sub-dict)
        flowResistanceModel(const dictionary& dict);

        //- Copy constructor
        flowResistanceModel(const flowResistanceModel&);

        //- Clone
        autoPtr<flowResistanceModel> clone() const;


    // Member Functions

        //- Return the loss coefficient K
        scalar lossCoefficient() const;

        //- Write the dictionary entries to stream
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
