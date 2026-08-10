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
    Computes the Reynolds-number-dependent loss coefficient zeta for
    planar flow resistance elements (perforated plates, slotted
    openings) using the correlations of the vent-openings formulation,
    or a user-specified constant.

    The pressure drop across the resistance is:

    \f[
        \Delta p = \zeta \frac{1}{2} \rho U_n^2
    \f]

    where \c zeta is the loss coefficient, \c U_n is the face-normal
    velocity and \c rho is the density (1 for incompressible flow).
    The patch-averaged normal velocity is obtained from the velocity
    flux sum \c Q = sum(U.Sf) over the patch faces and the patch
    area \c A:

    \f[
        U_n = \frac{Q}{A}
    \f]

    The signed correction \c sign(U_n)*deltaP is under-relaxed with
    the optional \c relax entry (0 < relax <= 1, default 1) so that
    flow reversal is also smooth:

    \f[
        corr_n = (1 - relax) corr_{n-1} + relax * target_n
    \f]

    The relaxed correction is computed by \c updateCorrection() and can
    be written as \c relaxedPressureCorrection for restart.

    For \c perforated and \c slotted the loss coefficient is
    Reynolds-number-dependent. With \c f the open area fraction and
    \c d_h = 4A_h/P_h the hydraulic diameter of a single hole, the hole
    Reynolds number is

    \f[
        Re = \frac{U_n d_h}{f \nu}
    \f]

    where \c nu is the kinematic viscosity (from the turbulence model,
    or from the optional \c nu entry). The loss coefficient is then
    evaluated in four Reynolds regimes. The fully turbulent coefficient
    is

    \f[
        \zeta_{turb} = \frac{[0.707\sqrt{1-f} + 1-f]^2}{f^2}
    \f]

    multiplied by a slot shape factor for the \c slotted type.

    Supported resistance types:
    - \c perforated: Re-dependent zeta (four Reynolds regimes)
    - \c slotted: Re-dependent zeta with the slot shape factor k_shape
    - \c screen: placeholder (returns 0)
    - \c constant: user-specified zeta via \c lossCoefficient entry

    For \c perforated and \c slotted the hole geometry must always be
    provided in the respective \c perforatedCoeffs or \c slottedCoeffs
    sub-dict. The open area fraction f is either specified directly via
    the \c openAreaFraction entry or derived from the geometry and the
    pitches:

    \f[
        f = \frac{A_h}{S_l S_t}
    \f]

    where \c A_h is the area of a single hole:
    - perforated: \f$ A_h = \frac{\pi d^2}{4} \f$ for circular holes or
      \f$ A_h = d^2 \f$ for square holes
    - slotted: \f$ A_h = b(a-b) + \frac{\pi b^2}{4} \f$ for slot length
      \f$ a \f$ and slot width \f$ b \f$ (a > b)

    The two options are mutually exclusive.

    At the first update the number of holes implied by f and the hole
    geometry is validated against the patch area: at least one hole must
    fit on the patch.

Usage
    Example with direct open area fraction:
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

    Example with slotted geometry:
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

    Example with constant K:
    \verbatim
    flowResistanceCoeffs
    {
        resistanceType      constant;
        lossCoefficient     2.5;
        relax               0.3;
    }
    \endverbatim

    For \c perforated and \c slotted the kinematic viscosity is obtained
    from the turbulence model (\c nu = \c mu/rho). The optional \c nu
    entry overrides this, e.g. when no turbulence model is present:

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

        nu                  1.5e-5;
    }
    \endverbatim

    A restarting run may supply the optional \c relaxedPressureCorrection
    entry to resume from a previously written relaxed correction.

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

        //- Longitudinal pitch [m] (S_l)
        scalar longitudinalPitch_;

        //- Transverse pitch [m] (S_t)
        scalar transversePitch_;

        //- Slot length for slotted geometry [m] (a: straight length +
        //  cap diameter)
        scalar slotLength_;

        //- Slot width for slotted geometry [m] (b: cap diameter =
        //  slot width)
        scalar slotWidth_;

        //- True if hole geometry inputs were used
        bool geometryUsed_;

        //- Under-relaxation factor for the signed pressure correction
        scalar relax_;

        //- User-specified kinematic viscosity [m^2/s]. Negative = not set.
        scalar userNu_;

        //- Target (unrelaxed) signed pressure correction from last update
        scalar targetCorrection_;

        //- Relaxed signed pressure correction
        scalar relaxedCorrection_;

        //- Whether the relaxed correction has been initialised
        bool relaxedCorrectionInitialised_;

        //- Whether the patch area has been validated against the hole
        //  geometry (i.e. that at least one hole fits on the patch)
        bool patchAreaValidated_;


    // Private Member Functions

        //- Area of a single hole [m^2] for the current resistance
        //  geometry (perforated or slotted)
        scalar holeArea() const;

        //- Perimeter of a single hole [m] for the current resistance
        //  geometry (perforated or slotted)
        scalar holePerimeter() const;

        //- Hydraulic diameter of a single hole [m] = 4A/P
        scalar hydraulicDiameter() const;

        //- Turbulent loss coefficient zeta_turb (Eq. 13), including the
        //  slot shape factor k_shape (Eq. 16) for the slotted type
        scalar zetaTurbulent() const;

        //- Loss coefficient zeta for the current resistance geometry:
        //  user-specified for 'constant', otherwise Re-dependent
        //  (Eqs. 9-12) for perforated/slotted
        scalar lossCoefficient
        (
            const scalar UnAvg,
            const scalar nuAvg
        ) const;


public:

    // Constructors

        //- Construct from dictionary (the flowResistanceCoeffs sub-dict)
        flowResistanceModel(const dictionary& dict);

        //- Copy constructor
        flowResistanceModel(const flowResistanceModel&);

        //- Clone
        autoPtr<flowResistanceModel> clone() const;


    // Member Functions

        //- The user-specified kinematic viscosity [m^2/s]
        //  (negative if not set)
        scalar userNu() const
        {
            return userNu_;
        }

        //- True if this resistance requires the kinematic viscosity
        //  (i.e. it is Reynolds-number dependent) and no user-specified
        //  \c nu overrides it. Used by the caller to decide whether a
        //  turbulence model lookup is required.
        bool needsNu() const
        {
            return
            (
                userNu_ < 0
             && (resistanceType_ == "perforated" || resistanceType_ == "slotted")
            );
        }

        //- Empirical coefficient eps0(Re) (Eqs. 10-11), log-linear
        //  interpolation of the Section D table
        scalar eps0(const scalar Re) const;

        //- Empirical coefficient zeta_phi(Re, f) (Eq. 11), log-linear
        //  in Re and linear in f interpolation of the Section D table
        scalar zetaPhi(const scalar Re, const scalar f) const;

        //- Compute and return the relaxed signed pressure correction
        //  for the given velocity flux sum Q, patch area A,
        //  area-averaged density rhoAvg (1 for incompressible flow)
        //  and area-averaged kinematic viscosity nuAvg:
        //  UnAvg = Q/A, deltaP = zeta(Re)*0.5*rhoAvg*UnAvg^2,
        //  target = sign(UnAvg)*deltaP, then apply the relaxation.
        //  The hole Reynolds number is Re = |UnAvg|*d_h/(f*nuAvg).
        scalar updateCorrection
        (
            const scalar Q,
            const scalar A,
            const scalar rhoAvg,
            const scalar nuAvg
        );

        //- The target (unrelaxed) signed correction from the last update
        scalar targetCorrection() const
        {
            return targetCorrection_;
        }

        //- The current relaxed signed correction
        scalar relaxedCorrection() const
        {
            return relaxedCorrection_;
        }

        //- Write the dictionary entries to stream
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
