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
    Foam::regionModels::areaSurfaceFilmModels::KelvinHelmholtz

Description
    Kelvin-Helmholtz sheet-breakup model.  Gas shear selects the most unstable
    sheet wavelength, which forms a ligament that then undergoes capillary
    breakup into a physical droplet size.

Usage
    \verbatim
    diameterModels
    {
        <morphology>
        {
            // Mandatory entries
            type        KelvinHelmholtz;
            rhop        <scalar>;

            // Optional entries
            Cdrop       <scalar>;
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      type         | Type name: KelvinHelmholtz         | word | yes  | -
      rhop         | Primary-phase density            | scalar | yes  | -
      Cdrop        | Ligament-to-drop diameter ratio  | scalar | no   | 1.88
    \endtable

SourceFiles
    KelvinHelmholtz.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_KelvinHelmholtz_H
#define regionFaModels_KelvinHelmholtz_H

#include "diameterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                     Class KelvinHelmholtz Declaration
\*---------------------------------------------------------------------------*/

class KelvinHelmholtz
:
    public diameterModel
{
    // Private Data

        //- Primary-phase density
        scalar rhop_;

        //- Ligament-to-drop diameter ratio
        scalar Cdrop_;


public:

    //- Runtime type information
    TypeName("KelvinHelmholtz");


    // Static Member Functions

        //- Return the Kelvin-Helmholtz physical diameter for one face
        static scalar dParcel
        (
            const scalar h,
            const scalar rho,
            const scalar mu,
            const scalar sigma,
            const scalar magSqrUrel,
            const scalar rhop,
            const scalar Cdrop,
            const scalar hFallback
        );


    // Constructors

        //- Construct from the film and model dictionary
        KelvinHelmholtz
        (
            const liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~KelvinHelmholtz() = default;


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
