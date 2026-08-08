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
    Foam::regionModels::areaSurfaceFilmModels::regimeMapMorphology

Description
    Local deterministic morphology classifier based on the gas Weber number:

    \f[
        We_g = \frac{\rho_g |U_p-U_f|^2 h}{\sigma}.
    \f]

    Low, intermediate, and high aerodynamic forcing select the configured
    low-forcing state, ligament, and sheet, respectively.  Thresholds are
    calibration coefficients and are not universal physical constants.

Usage
    \verbatim
    morphology
    {
        // Mandatory entries
        type regimeMap;

        regimeMapCoeffs
        {
            // Mandatory entries
            rhop                    1.2;
            WeDripToLigament        1;
            WeLigamentToSheet       10;

            // Optional entries
            lowWeMorphology         drip;
        }
    }
    \endverbatim

    The optional \c lowWeMorphology entry defaults to \c drip and must be
    either \c directDroplet or \c drip.  The threshold entries must satisfy
    \f$0 \mathrel{\leq} WeDripToLigament < WeLigamentToSheet\f$.

SourceFiles
    regimeMapMorphology.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_regimeMapMorphology_H
#define regionFaModels_regimeMapMorphology_H

#include "morphologyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                   Class regimeMapMorphology Declaration
\*---------------------------------------------------------------------------*/

class regimeMapMorphology
:
    public morphologyModel
{
    // Private Data

        //- Primary-phase density [kg/m3]
        scalar rhop_;

        //- Transition from the low-forcing state to ligament
        scalar WeDripToLigament_;

        //- Transition from ligament to sheet
        scalar WeLigamentToSheet_;

        //- Low-forcing morphology, directDroplet or drip
        morphologyType lowWeMorphology_;


public:

    //- Runtime type information
    TypeName("regimeMap");


    // Constructors

        //- Construct from the film and model dictionary
        regimeMapMorphology
        (
            const liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~regimeMapMorphology() = default;


    // Member Functions

        //- Classify each film face by its local gas Weber number
        virtual tmp<labelField> morphology() const;

        //- Return all morphology types reachable by the configured map
        virtual boolList possibleTypes() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
