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
    Foam::regionModels::areaSurfaceFilmModels::RayleighPlateau

Description
    Capillary ligament breakup model:

    \f[
        d = C_{lig} \, h(1 + 3\,Oh)^{1/6}, \qquad
        Oh = \frac{\mu}{\sqrt{\rho_l\, \sigma \, h}}.
    \f]

    where:
    \vartable
      d            | Diameter [m]
      C_{lig}      | Ligament-to-drop diameter ratio
      Oh           | Ohnesorge correction
      h            | Film thickness
      \mu          | Film dynamic viscosity
      \rho_l       | Film density
      \sigma       | Surface tension
    \endvartable


Usage
    \verbatim
    diameterModels
    {
        <morphology>
        {
            // Mandatory entries
            type                    RayleighPlateau;

            // Optional entries
            Cligament               <scalar>;
            viscousCorrection       <bool>;
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      type         | Type name: RayleighPlateau         | word | yes  | -
      Cligament    | Ligament-to-drop diameter ratio  | scalar | no   | 1.89
      viscousCorrection | Ohnesorge correction          | bool | no   | false
    \endtable

SourceFiles
    RayleighPlateau.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_RayleighPlateau_H
#define regionFaModels_RayleighPlateau_H

#include "diameterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                     Class RayleighPlateau Declaration
\*---------------------------------------------------------------------------*/

class RayleighPlateau
:
    public diameterModel
{
    // Private Data

        //- Ligament-to-drop diameter ratio
        scalar Cligament_;

        //- Apply the Ohnesorge correction
        bool viscousCorrection_;


public:

    //- Runtime type information
    TypeName("RayleighPlateau");


    // Static Member Functions

        //- Return the Rayleigh-Plateau physical diameter for one face
        static scalar dParcel
        (
            const scalar h,
            const scalar rho,
            const scalar mu,
            const scalar sigma,
            const scalar Cligament,
            const bool viscousCorrection
        );


    // Constructors

        //- Construct from the film and model dictionary
        RayleighPlateau
        (
            const liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~RayleighPlateau() = default;


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
