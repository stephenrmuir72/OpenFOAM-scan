/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 OpenCFD Ltd.
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
    Foam::pureZoneMixtureCorrection

Description
    Apply zone-mixture diffusion correction to a solid enthalpy equation.

    For zone-based mixtures where kappa != alpha * Cp, the standard enthalpy
    equation solves:
        ddt(rho*h) - laplacian(alpha, h) = 0

    This function adds explicit correction terms so that at convergence the
    equation correctly represents:
        ddt(rho*h) - laplacian(kappa, T) = 0

    The correction relies on cancellation of fvm:: and fvc:: laplacian terms
    at convergence:
        fvm::laplacian(betav*coeff, h) - fvc::laplacian(betav*coeff, h) -> 0

    leaving:
        - fvc::laplacian(kappaf, T)

Usage
    \code
    pureZoneMixtureCorrection(hEqn, thermo, betav, h);

    // For anisotropic:
    pureZoneMixtureCorrection
    (
        hEqn, thermo, betav, h,
        thermo.isotropic() ? nullptr : &taniAlpha()
    );
    \endcode

SourceFiles
    pureZoneMixtureCorrection.H (this file)

\*---------------------------------------------------------------------------*/

#ifndef pureZoneMixtureCorrection_H
#define pureZoneMixtureCorrection_H

#include "fvScalarMatrix.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "solidThermo.H"
#include "calcHarmonicKappaf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Function pureZoneMixtureCorrection
\*---------------------------------------------------------------------------*/

//- Apply zone-mixture correction to a solid enthalpy equation.
//
//  \param hEqn        The enthalpy equation matrix to modify
//  \param thermo      The solid thermo (provides alpha, kappa, T, Cp, etc.)
//  \param betav       The volume porosity field
//  \param h           The enthalpy field
//  \param aniAlphaPtr Pointer to the anisotropic tensor field
//                     (nullptr for isotropic materials)
inline void pureZoneMixtureCorrection
(
    fvScalarMatrix& hEqn,
    const solidThermo& thermo,
    const volScalarField& betav,
    const volScalarField& h,
    const volSymmTensorField* aniAlphaPtr = nullptr
)
{
    if (!thermo.isZoneMixture()) return;

    if (thermo.isotropic())
    {
        // +laplacian(betav*alpha, h) - laplacian(betav*kappa, T)
        //
        // At convergence: fvm::laplacian(betav*alpha, h) from the main
        // equation and this fvc::laplacian(betav*alpha, h) cancel,
        // leaving -fvc::laplacian(betav*kappa, T).
        hEqn += fvc::laplacian(betav*thermo.alpha(), h, "laplacian(alpha,h)");
        hEqn -= fvc::laplacian(betav*thermo.kappa(), thermo.T(), "laplacian(kappa,T)");
    }
    else
    {
        // +laplacian(betav*taniAlpha, h) - laplacian(tkappaf, T)
        //
        // At convergence: fvm::laplacian(betav*taniAlpha, h) from the main
        // equation and this fvc::laplacian(betav*taniAlpha, h) cancel,
        // leaving -fvc::laplacian(tkappaf, T).

        const auto& aniAlpha = *aniAlphaPtr;

        const auto tkappaf =
            fvc::interpolate(betav)
          * calcHarmonicKappaf(h.mesh(), aniAlpha, thermo.Cp());

        hEqn += fvc::laplacian(betav*aniAlpha, h, "laplacian(alpha,h)");
        hEqn -= fvc::laplacian(tkappaf, thermo.T(), "laplacian(kappa,T)");
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
