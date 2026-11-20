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

Namespace
    Foam::flowResistanceModels

Description
    A namespace for various flow resistance model implementations.

Class
    Foam::flowResistanceModel

Description
    Runtime-selectable model computing the loss coefficient zeta for
    planar flow resistance elements attached to a \c totalPressure
    boundary condition.

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

    Available resistance types:
    - \c perforated: Re-independent fully turbulent zeta
    - \c slotted: Re-independent fully turbulent zeta with the slot
      shape factor k_shape
    - \c screen: Re-dependent zeta (three Reynolds regimes, strand
      shape factor k_shape)
    - \c constant: user-specified zeta via \c lossCoefficient entry

    For \c screen the kinematic viscosity is looked up directly from
    the registered 'nu' field (incompressible flow) or derived from the
    registered 'mu' and 'rho' fields (compressible flow:
    \c nu = \c mu/rho). The optional \c nu entry overrides this, e.g.
    when no such field is registered. The \c nu entry is only permitted
    for \c screen; it is rejected for \c perforated, \c slotted and
    \c constant.

    \verbatim
    flowResistanceCoeffs
    {
        resistanceType      screen;

        screenCoeffs
        {
            strandShape         round;
            strandSize          0.001;
            openAreaFraction    0.6;
        }

        nu                  1.5e-5;
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

    A restarting run may supply the optional \c relaxedPressureCorrection
    entry to resume from a previously written relaxed correction.

    Reference:
    \verbatim
        Idelchik, I. E. (2008).
        Handbook of Hydraulic Resistance (3rd ed.).
        Begell House.
        (Re-dependent loss coefficient correlations and eps0(Re)/
        zeta_phi(Re, f) tables for perforated plates, slotted
        openings and screens/grids.)
    \endverbatim

SourceFiles
    flowResistanceModel.C
    flowResistanceModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_flowResistanceModel_H
#define Foam_flowResistanceModel_H

#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class flowResistanceModel Declaration
\*---------------------------------------------------------------------------*/

class flowResistanceModel
{
protected:

    // Protected Data

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

        //- Whether the patch area has been validated against the
        //  resistance geometry (i.e. that at least one hole fits on
        //  the patch), independent of relaxedCorrectionInitialised_
        //  so that it is still performed on a restarted run
        bool patchAreaValidated_;

        //- Reynolds number from the last lossCoefficient() evaluation,
        //  for debug reporting. mutable: assigned from the const
        //  lossCoefficient() implementation. Zero when not applicable
        //  (e.g. the constant resistance type).
        mutable scalar Re_;


    // Protected Member Functions

        //- Loss coefficient zeta for the given area-averaged normal
        //  velocity and kinematic viscosity
        virtual scalar lossCoefficient
        (
            const scalar UnAvg,
            const scalar nuAvg
        ) const = 0;

        //- Validate the resistance geometry/open area fraction against
        //  the patch area (i.e. that at least one hole fits on the
        //  patch), called once on the first updateCorrection(). The
        //  default implementation does nothing; hole-based resistance
        //  types override this.
        virtual void validatePatchArea(const scalar A) const;

        //- Write the resistanceType-specific dictionary entries
        //  (i.e. the "<resistanceType>Coeffs" sub-dict and/or
        //  "openAreaFraction"/"lossCoefficient" as applicable)
        virtual void writeCoeffs(Ostream& os) const = 0;


public:

    //- Runtime type information
    TypeName("flowResistanceModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        flowResistanceModel,
        dictionary,
        (const dictionary& dict),
        (dict)
    );


    // Constructors

        //- Construct from dictionary (the flowResistanceCoeffs sub-dict)
        explicit flowResistanceModel(const dictionary& dict);

        //- Copy constructor
        flowResistanceModel(const flowResistanceModel&);

        //- Construct and return a clone
        virtual autoPtr<flowResistanceModel> clone() const = 0;

        //- No copy assignment
        void operator=(const flowResistanceModel&) = delete;


    //- Selector
    static autoPtr<flowResistanceModel> New(const dictionary& dict);


    //- Destructor
    virtual ~flowResistanceModel() = default;


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
        virtual bool needsNu() const = 0;

        //- Compute and return the relaxed signed pressure correction
        //  for the given velocity flux sum Q, patch area A,
        //  area-averaged density rhoAvg (1 for incompressible flow)
        //  and area-averaged kinematic viscosity nuAvg:
        //  UnAvg = Q/A, deltaP = zeta(Re)*0.5*rhoAvg*UnAvg^2,
        //  target = sign(UnAvg)*deltaP, then apply the relaxation.
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
