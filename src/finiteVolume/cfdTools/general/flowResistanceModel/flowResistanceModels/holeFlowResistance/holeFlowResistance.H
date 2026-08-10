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
    Foam::flowResistanceModels::holeFlowResistance

Description
    Abstract intermediate base class for Reynolds-number-dependent
    flow resistances through a regular array of holes (\c perforated,
    \c slotted). Provides the shared open-area-fraction resolution and
    the Re-dependent loss coefficient correlation:

    \f[
        \zeta_{turb} = \frac{[0.707\sqrt{1-f} + 1-f]^2}{f^2}
    \f]

    evaluated in four Reynolds regimes (Eqs. 9-12), where \c f is the
    open area fraction and the hole Reynolds number is

    \f[
        Re = \frac{U_n d_h}{f \nu}
    \f]

    with \c d_h = 4A_h/P_h the hydraulic diameter of a single hole.
    Derived classes provide the single-hole area/perimeter (and any
    shape factor applied to \c zeta_turb, e.g. the slotted \c k_shape).

    The open area fraction f is either specified directly via the
    \c openAreaFraction entry or derived from the hole geometry and
    the pitches:

    \f[
        f = \frac{A_h}{S_l S_t}
    \f]

    The two options are mutually exclusive. At the first update the
    number of holes implied by f and the hole geometry is validated
    against the patch area: at least one hole must fit on the patch.

SourceFiles
    holeFlowResistance.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_flowResistanceModels_holeFlowResistance_H
#define Foam_flowResistanceModels_holeFlowResistance_H

#include "flowResistanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace flowResistanceModels
{

/*---------------------------------------------------------------------------*\
                    Class holeFlowResistance Declaration
\*---------------------------------------------------------------------------*/

class holeFlowResistance
:
    public flowResistanceModel
{
protected:

    // Protected Data

        //- Open area fraction f (0..1)
        scalar openAreaFraction_;

        //- Longitudinal pitch [m] (S_l). Negative if not used, i.e.
        //  openAreaFraction was specified directly.
        scalar longitudinalPitch_;

        //- Transverse pitch [m] (S_t). Negative if not used.
        scalar transversePitch_;


private:

    // Private Member Functions

        //- Empirical coefficient eps0(Re) (Eqs. 10-11), log-linear
        //  interpolation of the Section D table
        static scalar eps0(const scalar Re);

        //- Empirical coefficient zeta_phi(Re, f) (Eq. 11), log-linear
        //  in Re and linear in f interpolation of the Section D table
        static scalar zetaPhi(const scalar Re, const scalar f);


protected:

    // Protected Member Functions

        //- Area of a single hole [m^2]
        virtual scalar holeArea() const = 0;

        //- Perimeter of a single hole [m]
        virtual scalar holePerimeter() const = 0;

        //- Hydraulic diameter of a single hole [m] = 4A/P
        scalar hydraulicDiameter() const;

        //- Shape factor applied to the fully turbulent loss
        //  coefficient (Eq. 16 for the slotted type; 1 for perforated)
        virtual scalar shapeFactor() const;

        //- Turbulent loss coefficient zeta_turb (Eq. 13), including
        //  the shape factor
        scalar zetaTurbulent() const;

        //- Resolve the open area fraction f from the given
        //  \c typeDict (the "<resistanceType>Coeffs" sub-dict):
        //  either taken directly from the 'openAreaFraction' entry,
        //  or derived from the hole/slot geometry (via holeArea(),
        //  which must already be usable at this point) and the
        //  'longitudinalPitch'/'transversePitch' entries. The two
        //  options are validated as mutually exclusive and the
        //  result is validated to lie in (0, 1]. \c geometryDesc is
        //  used to word the "check ... geometry and pitch" error
        //  message (e.g. "hole" or "slot").
        scalar resolveOpenAreaFraction
        (
            const dictionary& dict,
            const dictionary& typeDict,
            const word& geometryDesc
        );

        //- Loss coefficient zeta for the given area-averaged normal
        //  velocity and kinematic viscosity (Eqs. 9-12)
        virtual scalar lossCoefficient
        (
            const scalar UnAvg,
            const scalar nuAvg
        ) const;

        //- Validate that at least one hole fits on a patch of the
        //  given area
        virtual void validatePatchArea(const scalar A) const;


public:

    // Constructors

        //- Construct from dictionary (the flowResistanceCoeffs
        //  sub-dict). The open area fraction is not resolved by this
        //  constructor; derived classes must call
        //  resolveOpenAreaFraction() once their geometry is read.
        explicit holeFlowResistance(const dictionary& dict);

        //- Copy constructor
        holeFlowResistance(const holeFlowResistance&) = default;


    //- Destructor
    virtual ~holeFlowResistance() = default;


    // Member Functions

        //- The open area fraction f
        scalar openAreaFraction() const
        {
            return openAreaFraction_;
        }

        //- A hole-based resistance requires the kinematic viscosity
        //  unless the user has supplied an override
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
