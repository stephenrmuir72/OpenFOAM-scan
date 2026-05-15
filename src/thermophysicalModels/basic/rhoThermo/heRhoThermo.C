/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "heRhoThermo.H"

#include <cmath>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

inline bool validTemperature(const Foam::scalar T)
{
    return std::isfinite(T) && T > 0;
}


inline void checkTemperature
(
    const Foam::scalar T,
    const Foam::scalar p,
    const Foam::scalar he,
    const Foam::label celli
)
{
    if (FOAM_UNLIKELY(!validTemperature(T)))
    {
        using namespace Foam;

        FatalErrorInFunction
            << "Invalid temperature in cell " << celli
            << ": T=" << T
            << ", p=" << p
            << ", he=" << he << nl
            << Foam::abort(Foam::FatalError);
    }
}


inline void checkTemperature
(
    const Foam::scalar T,
    const Foam::scalar p,
    const Foam::word& patchName,
    const Foam::label patchi,
    const Foam::label facei
)
{
    if (FOAM_UNLIKELY(!validTemperature(T)))
    {
        using namespace Foam;

        FatalErrorInFunction
            << "Invalid temperature on patch " << patchi
            << " (" << patchName << "), face " << facei
            << ": T=" << T
            << ", p=" << p << nl
            << Foam::abort(Foam::FatalError);
    }
}


inline void checkTemperature
(
    const Foam::scalar T,
    const Foam::scalar p,
    const Foam::scalar he,
    const Foam::word& patchName,
    const Foam::label patchi,
    const Foam::label facei
)
{
    if (FOAM_UNLIKELY(!validTemperature(T)))
    {
        using namespace Foam;

        FatalErrorInFunction
            << "Invalid temperature on patch " << patchi
            << " (" << patchName << "), face " << facei
            << ": T=" << T
            << ", p=" << p
            << ", he=" << he << nl
            << Foam::abort(Foam::FatalError);
    }
}

}  // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::calculate
(
    const volScalarField& p,
    volScalarField& T,
    volScalarField& he,
    volScalarField& psi,
    volScalarField& rho,
    volScalarField& mu,
    volScalarField& alpha,
    const bool doOldTimes
)
{
    // Note: update oldTimes before current time so that if T.oldTime() is
    // created from T, it starts from the unconverted T
    if (doOldTimes && (p.nOldTimes() || T.nOldTimes()))
    {
        calculate
        (
            p.oldTime(),
            T.oldTime(),
            he.oldTime(),
            psi.oldTime(),
            rho.oldTime(),
            mu.oldTime(),
            alpha.oldTime(),
            true
        );
    }

    const scalarField& hCells = he.primitiveField();
    const scalarField& pCells = p.primitiveField();

    scalarField& TCells = T.primitiveFieldRef();
    scalarField& psiCells = psi.primitiveFieldRef();
    scalarField& rhoCells = rho.primitiveFieldRef();
    scalarField& muCells = mu.primitiveFieldRef();
    scalarField& alphaCells = alpha.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        if (this->updateT())
        {
            TCells[celli] = mixture_.THE
            (
                hCells[celli],
                pCells[celli],
                TCells[celli]
            );
        }

        checkTemperature(TCells[celli], pCells[celli], hCells[celli], celli);

        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
        rhoCells[celli] = mixture_.rho(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
    }

    const volScalarField::Boundary& pBf = p.boundaryField();
    volScalarField::Boundary& TBf = T.boundaryFieldRef();
    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();
    volScalarField::Boundary& rhoBf = rho.boundaryFieldRef();
    volScalarField::Boundary& heBf = he.boundaryFieldRef();
    volScalarField::Boundary& muBf = mu.boundaryFieldRef();
    volScalarField::Boundary& alphaBf = alpha.boundaryFieldRef();

    forAll(pBf, patchi)
    {
        const fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                checkTemperature
                (
                    pT[facei],
                    pp[facei],
                    pp.patch().name(),
                    patchi,
                    facei
                );

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                if (this->updateT())
                {
                    pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);
                }

                checkTemperature
                (
                    pT[facei],
                    pp[facei],
                    phe[facei],
                    pp.patch().name(),
                    patchi,
                    facei
                );

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heRhoThermo<BasicPsiThermo, MixtureType>::heRhoThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName)
{
    calculate
    (
        this->p_,
        this->T_,
        this->he_,
        this->psi_,
        this->rho_,
        this->mu_,
        this->alpha_,
        true                    // Create old time fields
    );
}


template<class BasicPsiThermo, class MixtureType>
Foam::heRhoThermo<BasicPsiThermo, MixtureType>::heRhoThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName
)
:
    heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName, dictName)
{
    calculate
    (
        this->p_,
        this->T_,
        this->he_,
        this->psi_,
        this->rho_,
        this->mu_,
        this->alpha_,
        true                    // Create old time fields
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heRhoThermo<BasicPsiThermo, MixtureType>::~heRhoThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::correct()
{
    DebugInFunction << endl;

    calculate
    (
        this->p_,
        this->T_,
        this->he_,
        this->psi_,
        this->rho_,
        this->mu_,
        this->alpha_,
        false           // No need to update old times
    );

    DebugInFunction << "Finished" << endl;
}

// ************************************************************************* //
