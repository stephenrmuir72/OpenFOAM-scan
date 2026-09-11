/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  | www.openfoam.com
-------------------------------------------------------------------------------
    Copyright (C) 2026 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::fanNonDimensional

Description
    Common non-dimensional fan-curve formulation used by the \c fan and
    \c fanPressure boundary conditions.

    The formulation follows ISO 13348:2007.

    The fan flow coefficient (non-dimensional flow rate) is calculated as:

    \verbatim
        Phi = 60*Q/(PI*dm^3*rpm)
    \endverbatim

    where \c Q is the volume flow rate, \c dm is the impeller diameter and 
    \c rpm is the fan speed.

    The pressure coefficient (non-dimensional \c deltaP)

    For fields with dimensions \c dimPressure, where pressure is in Pa:

    \verbatim
        Psi = 3600*deltaP/(rho*(sqr(PI*rpm*dm)))
    \endverbatim

    The density \c rho used to convert the pressure coefficient back to Pa is
    \c rhoRef by default. If \c rhoCorrection is enabled, the area-weighted
    patch density is used instead. This correction accounts for density
    effects when applying the fan curve.

    For fields with dimensions \c dimPressure/dimDensity, the field is a
    kinematic pressure and density does not appear in the formulation:

    \verbatim
        Psi = 3600*deltaP/(sqr(PI*rpm*dm))
    \endverbatim

    In this case, \c deltaP is a pressure drop per unit density and has
    dimensions \c dimPressure/dimDensity (m^2/s^2).

    If \c legacy is enabled, the previous formulation is used:

    \verbatim
        Phi = 120*Q/(PI^3*dm^3*rpm)
        Psi = 1800*deltaP/(PI^4*(dm*rpm)^2)
    \endverbatim

    This formulation is incorrect and should only be used to reproduce
    existing cases.

    This is an implementation helper and is not a runtime-selectable boundary
    condition.

Usage
    The helper reads the following entries when non-dimensional operation is
    enabled:
    \table
        Property       | Description                | Required  | Default
        rpm            | Fan rpm                    | yes       |
        dm             | Impeller diameter          | yes       |
        legacy         | Use the previous formulation | no | false
        rhoCorrection  | Account for density effects | no | false
        rhoRef         | Rated density (when rhoCorrection = false) | no | 1.2
    \endtable

\*---------------------------------------------------------------------------*/

#ifndef Foam_fanNonDimensional_H
#define Foam_fanNonDimensional_H

#include "Function1.H"
#include "dimensionSet.H"
#include "fileName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;
class fvPatch;
class objectRegistry;
class Ostream;

/*---------------------------------------------------------------------------*\
                    Class fanNonDimensional Declaration
\*---------------------------------------------------------------------------*/

class fanNonDimensional
{
    // Private Static Member Functions

        //- The default reference density (1.2)
        static constexpr scalar defaultRho() noexcept { return 1.2; }


    // Private Data

        //- Patch on which the fan is applied
        const fvPatch& patch_;

        //- Dimensions of the field to which the fan is applied
        const dimensionSet dimensions_;

        //- Field name used in diagnostics
        const word fieldName_;

        //- Field path used in diagnostics
        const fileName objectPath_;

        //- Whether non-dimensional operation is enabled
        const bool enabled_;

        //- Fan rpm
        autoPtr<Function1<scalar>> rpm_;

        //- Fan mean diameter
        autoPtr<Function1<scalar>> dm_;

        //- Use the previous non-dimensional formulation
        bool legacy_;

        //- Use patch-averaged density instead of rated density
        bool rhoCorrection_;

        //- Whether rhoRef was explicitly specified
        bool rhoRefUserProvided_;

        //- Rated/reference density
        scalar rhoRef_;


public:

    // Constructors

        //- Construct with non-dimensional operation disabled
        fanNonDimensional
        (
            const fvPatch&,
            const dimensionSet&,
            const word& fieldName,
            const fileName& objectPath
        );

        //- Construct from dictionary
        fanNonDimensional
        (
            const fvPatch&,
            const dimensionSet&,
            const word& fieldName,
            const fileName& objectPath,
            const dictionary&,
            const objectRegistry&,
            const bool enabled
        );

        //- Construct by mapping onto a new patch
        fanNonDimensional
        (
            const fanNonDimensional&,
            const fvPatch&,
            const dimensionSet&,
            const word& fieldName,
            const fileName& objectPath
        );

        //- No copy without a patch
        fanNonDimensional(const fanNonDimensional&) = delete;


    // Member Functions

        //- Convert a physical volume flow rate into a table input
        scalar nonDimensionalFlowRate
        (
            const scalar volumeFlowRate,
            const scalar timeValue
        ) const;

        //- Convert a table pressure value into a dimensional pressure drop
        scalar dimensionalPressure
        (
            const scalar nonDimensionalPressure,
            const word& rhoName,
            const scalar timeValue
        ) const;

        //- Write the non-dimensional entries
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
