/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

    Some parts here added have been inspired by the work done at the following
    repository https://github.com/SimFlowCFD/RapidCFD-dev
\*---------------------------------------------------------------------------*/

#include "adjustPhi.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

#ifdef FOAM_OFFLOAD
namespace Foam
{
    template<bool positive>
    struct adjustPhiFilter : public std::unary_function<scalar,scalar>
    {
        scalar operator()(const scalar& s)
        {
            if (positive)
                return s>0?s:0;
            else
                return s<0?s:0;
        } 
};

    struct adjustPhiCorrectFunctor : public std::unary_function<scalar,scalar>
    {
        const scalar massCorr;
        adjustPhiCorrectFunctor(const scalar _massCorr):massCorr(_massCorr){}
        scalar operator()(const scalar& phi)
        {
            if (phi > 0.0)
                return phi * massCorr;
            else
                return phi;
        }
    };
}
#endif

bool Foam::adjustPhi
(
    surfaceScalarField& phi,
    const volVectorField& U,
    volScalarField& p
)
{
    if (p.needReference())
    {
        scalar massIn = 0.0;
        scalar fixedMassOut = 0.0;
        scalar adjustableMassOut = 0.0;

        surfaceScalarField::Boundary& bphi =
            phi.boundaryFieldRef();

        forAll(bphi, patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            const fvsPatchScalarField& phip = bphi[patchi];

            if (!phip.coupled())
            {
                #ifdef FOAM_OFFLOAD
                apply_policy(phip.size(), [&](auto policy)                                    
                { 
                    massIn =  std::transform_reduce
                        (
                            policy,
                            phip.begin(),
                            phip.end(),
                            0.0, // initial value
                            std::minus<scalar>(),
                            [](const scalar& val) 
                            { 
                                return adjustPhiFilter<false>()(val); 
                            }
                         );
                });
                #endif
                if (Up.fixesValue() && !isA<inletOutletFvPatchVectorField>(Up))
                {
                    #ifdef FOAM_OFFLOAD
                    apply_policy(phip.size(), [&](auto policy)
                    {  
                        fixedMassOut = std::transform_reduce
                            (
                                policy,
                                phip.begin(),
                                phip.end(),
                                0.0, // initial value
                                std::plus<scalar>(),
                                [](const scalar& val) 
                                {
                                     return adjustPhiFilter<true>()(val); 
                                } // mapping function
                            );
                    });
                    #else
                    forAll(phip, i)
                    {
                        if (phip[i] < 0.0)
                        {
                            massIn -= phip[i];
                        }
                        else
                        {
                            fixedMassOut += phip[i];
                        }
                    }
                    #endif

                }
                else
                {
                    #ifdef FOAM_OFFLOAD
                    apply_policy(phip.size(), [&](auto policy)
                    {
                        adjustableMassOut = std::transform_reduce
                        (
                            policy,
                            phip.begin(),
                            phip.end(),
                            0.0, // initial value
                            std::plus<scalar>(), // sum operation
                            [](const scalar& val)
                            { 
                                return adjustPhiFilter<true>()(val); 
                            }
                        ); 
                    });
                    #else
                    forAll(phip, i)
                    {
                        if (phip[i] < 0.0)
                        {
                            massIn -= phip[i];
                        }
                        else
                        {
                            adjustableMassOut += phip[i];
                        }
                    }
                    #endif
                }
            }
        }

        // Calculate the total flux in the domain, used for normalisation
        scalar totalFlux = VSMALL + sum(mag(phi)).value();

        reduce(massIn, sumOp<scalar>());
        reduce(fixedMassOut, sumOp<scalar>());
        reduce(adjustableMassOut, sumOp<scalar>());

        scalar massCorr = 1.0;
        scalar magAdjustableMassOut = mag(adjustableMassOut);

        if
        (
            magAdjustableMassOut > VSMALL
         && magAdjustableMassOut/totalFlux > SMALL
        )
        {
            massCorr = (massIn - fixedMassOut)/adjustableMassOut;
        }
        else if (mag(fixedMassOut - massIn)/totalFlux > 1e-8)
        {
            FatalErrorInFunction
                << "Continuity error cannot be removed by adjusting the"
                   " outflow.\nPlease check the velocity boundary conditions"
                   " and/or run potentialFoam to initialise the outflow." << nl
                << "Total flux              : " << totalFlux << nl
                << "Specified mass inflow   : " << massIn << nl
                << "Specified mass outflow  : " << fixedMassOut << nl
                << "Adjustable mass outflow : " << adjustableMassOut << nl
                << exit(FatalError);
        }

        forAll(bphi, patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            fvsPatchScalarField& phip = bphi[patchi];

            if (!phip.coupled())
            {
                if
                (
                    !Up.fixesValue()
                 || isA<inletOutletFvPatchVectorField>(Up)
                )
                {
                    #ifdef FOAM_OFFLOAD
                    apply_policy(phip.size(), [&](auto policy)
                    {
                        std::transform
                        (   
                            policy,
                            phip.begin(),
                            phip.end(),
                            phip.begin(),
                            adjustPhiCorrectFunctor(massCorr)
                        );
                    });
                    #else
                    forAll(phip, i)
                    {
                        if (phip[i] > 0.0)
                        {
                            phip[i] *= massCorr;
                        }
                    }
                    #endif
                }
            }
        }

        return mag(massIn)/totalFlux < SMALL
            && mag(fixedMassOut)/totalFlux < SMALL
            && mag(adjustableMassOut)/totalFlux < SMALL;
    }

    return false;
}


// ************************************************************************* //
