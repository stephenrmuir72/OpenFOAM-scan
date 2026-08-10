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
    Flow resistance of a screen (mesh of strands), Section D.1.3 of the
    vent-opening formulation. The strands may be round wires (diameter
    \c d), square wires (side \c w) or thin flat ribbons (width \c w);
    the strand shape selects the shape factor

    \f[
        k_{shape} = 1.3,\ 1.8,\ 2.1  \quad \text{(round, square, flat)}
    \f]

    The open area fraction \c f is either specified directly via the
    \c openAreaFraction entry or derived from the strand size and the
    width/length separations:

    \f[
        f = \left(1 - \frac{d}{S_w}\right)\left(1 - \frac{d}{S_l}\right)
    \f]

    The two options are mutually exclusive. The velocity through the
    screen is \f$ U_{hole} = U_n/f \f$ and the strand Reynolds number is

    \f[
        Re = \frac{U_{hole}\,d}{\nu}
    \f]

    with the loss coefficient evaluated in three Reynolds regimes:

    \f[
        \zeta = \begin{cases}
            \frac{22}{Re} + \zeta_{turb} & Re < 50 \\
            \left(1 + \frac{20}{Re}\right)\zeta_{turb}
                & 50 \le Re < 1000 \\
            \zeta_{turb} & Re \ge 1000
        \end{cases}
    \f]

    with

    \f[
        \zeta_{turb} = k_{shape}(1-f) + \left(\frac{1-f}{f}\right)^2
    \f]

    The strand Reynolds number is clamped to a minimum of 1 to keep the
    laminar term finite at (near) zero velocity.

    \verbatim
    flowResistanceCoeffs
    {
        resistanceType      screen;

        screenCoeffs
        {
            strandShape         round;
            strandSize          0.001;
            widthSeparation     0.005;
            lengthSeparation    0.005;
            // or
            // openAreaFraction 0.6;
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

        //- Open area fraction f (0..1]
        scalar openAreaFraction_;

        //- Strand shape: 'round', 'square' or 'flat'
        word strandShape_;

        //- Strand size [m]: wire diameter d or wire/ribbon width w
        scalar strandSize_;

        //- Width separation [m] (S_w). Negative if not used, i.e.
        //  openAreaFraction was specified directly.
        scalar widthSeparation_;

        //- Length separation [m] (S_l). Negative if not used.
        scalar lengthSeparation_;


    // Private Member Functions

        //- Shape factor k_shape for the strand shape (1.3 round,
        //  1.8 square, 2.1 flat)
        scalar kShape() const;

        //- Turbulent loss coefficient zeta_turb
        scalar zetaTurbulent() const;


protected:

    // Protected Member Functions

        //- Loss coefficient zeta (Section D.1.3)
        virtual scalar lossCoefficient
        (
            const scalar UnAvg,
            const scalar nuAvg
        ) const;

        //- Write the "screenCoeffs" sub-dict
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

        //- A screen resistance is Reynolds-number dependent and so
        //  requires the kinematic viscosity unless the user has
        //  supplied an override
        virtual bool needsNu() const
        {
            return userNu_ < 0;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace flowResistanceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
