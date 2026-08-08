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
    Foam::regionModels::areaSurfaceFilmModels::aerodynamicDiameter

Description
    Critical-Weber physical droplet diameter:

    \f[
        d = \frac{We_{crit}\,\sigma}{\rho_g \, |U_p - U_f|^2}.
    \f]

    where:
    \vartable
      d            | Diameter [m]
      We_{crit}    | Critical Weber number [-]
      \sigma       | Surface tension [kg/s^2]
      \rho_g       | Primary-phase density [kg/m^3]
      U_p          | Primary-phase velocity [m/s]
      U_f          | Film velocity [m/s]
    \endvartable

    The model falls back to the local film thickness where the relative
    velocity vanishes.

Usage
    \verbatim
    diameterModels
    {
        <morphology>
        {
            // Mandatory entries
            type        aerodynamic;
            rhop        <scalar>;

            // Optional entries
            WeCrit      <scalar>;
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      type         | Type name: aerodynamic             | word | yes  | -
      rhop         | Primary-phase density            | scalar | yes  | -
      WeCrit       | Critical Weber number            | scalar | no   | 12.0
    \endtable

SourceFiles
    aerodynamicDiameter.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_aerodynamicDiameter_H
#define regionFaModels_aerodynamicDiameter_H

#include "diameterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class aerodynamicDiameter Declaration
\*---------------------------------------------------------------------------*/

class aerodynamicDiameter
:
    public diameterModel
{
    // Private Data

        //- Primary-phase density
        scalar rhop_;

        //- Critical Weber number
        scalar WeCrit_;


public:

    //- Runtime type information
    TypeName("aerodynamic");


    // Static Member Functions

        //- Return the critical-Weber physical diameter for one face
        static scalar dParcel
        (
            const scalar sigma,
            const scalar magSqrUrel,
            const scalar rhop,
            const scalar WeCrit,
            const scalar hFallback
        );


    // Constructors

        //- Construct from the film and model dictionary
        aerodynamicDiameter
        (
            const liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~aerodynamicDiameter() = default;


    // Member Functions

        //- Return the clipped physical diameter for every film face
        virtual tmp<scalarField> diameter() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
