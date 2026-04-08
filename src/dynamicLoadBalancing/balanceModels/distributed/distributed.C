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

\*---------------------------------------------------------------------------*/

#include "distributed.H"
#include "fvMeshDistribute.H"
#include "mapDistributePolyMesh.H"
#include "dynamicRefineFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"
#include "decompositionMethod.H"
#include "decompositionModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace balanceModels
{
    defineTypeNameAndDebug(distributed, 0);
    addToRunTimeSelectionTable(balanceModel, distributed, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::balanceModels::distributed::transformToIntegers
(
    const scalarField& w,
    const label targetSum
)
{
    const label n = w.size();
    labelList out(n, 0);

    if (n == 0 || targetSum <= 0)
    {
        return out;
    }

    // Sum of positive weights
    scalar sumW = 0.0;
    for (label i = 0; i < n; ++i)
    {
        if (w[i] > SMALL)
        {
            sumW += w[i];
        }
    }

    if (sumW <= SMALL)
    {
        return out;
    }


    // Normalise the weights
    scalarField wn(n, Foam::zero{});
    for (label i = 0; i < n; ++i)
    {
        if (w[i] > SMALL)
        {
            wn[i] = w[i]/sumW;
        }
        else
        {
            wn[i] = 0.0;
        }
    }


    // Base allocation by floor
    scalarList raw(n, Foam::zero{});
    for (label i = 0; i < n; ++i)
    {
        raw[i] = wn[i]*targetSum;
        out[i] = static_cast<label>(std::floor(raw[i]));
    }


    // Compute remainders (fractional parts) with indices
    std::vector<std::pair<scalar,label>> rem;
    rem.reserve(n);
    for (label i = 0; i < n; ++i)
    {
        const scalar ri = raw[i] - scalar(out[i]);
        rem.emplace_back(ri, i);
    }


    // Seats left to distribute
    label assigned = 0;
    for (label i = 0; i < n; ++i)
    {
        assigned += out[i];
    }
    label left = targetSum - assigned;


    // Distribute by largest remainders (stable tie-break by smaller index)
    std::stable_sort
    (
        rem.begin(), rem.end(),
        [](const std::pair<scalar,label>& a, const std::pair<scalar,label>& b)
        {
            if (a.first == b.first) return a.second < b.second;
            return a.first > b.first; // descending
        }
    );

    for (label k = 0; k < left; ++k)
    {
        ++out[rem[k].second];
    }

    return out;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::balanceModels::distributed::distributed
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    balanceModel(mesh, dict),
    distributionDictPtr_
    (
        new IOdictionary
        (
            IOobject
            (
                "decomposeParDict",
                mesh.time().system(),
                mesh.time(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        )
    )
{
    const labelList defaultWeights(UPstream::nProcs(), 1);

    dictionary coeffs;
    coeffs.add("processorWeights", defaultWeights);
    distributionDictPtr_->set("coeffs", coeffs);

    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::balanceModels::distributed::read(const dictionary& dict)
{
    if (!balanceModel::read(dict))
    {
        return false;
    }

    return true;
}


bool Foam::balanceModels::distributed::balance(const scalarField& procWeights)
{
    DebugInFunction
        << "Number of cells (before balance): " << mesh().nCells() << endl;

    // Check if the size of procWeights matches nProcs
    if (procWeights.size() != UPstream::nProcs())
    {
        FatalErrorInFunction
            << "processorWeights size: " << procWeights.size()
            << ", does not match nProcs: " << UPstream::nProcs() << "."
            << exit(FatalError);
    }


    // Check if weights are valid
    scalar sumW = gSum(procWeights);
    scalar minW = gMin(procWeights);
    if (!std::isfinite(sumW) || sumW <= SMALL || minW <= SMALL)
    {
        FatalErrorInFunction
            << "Invalid processorWeights "
            << "(sum:" << sumW << ", min:" << minW << ")."
            << exit(FatalError);
    }


    // Update the internal decomposition dictionary with processor weights
    dictionary coeffs;
    coeffs.set
    (
        "processorWeights",
        transformToIntegers(procWeights)
    );
    distributionDictPtr_->set("coeffs", coeffs);
    DebugInFunction<< "distributionPtr: " << *distributionDictPtr_ << endl;


    // Create decomposition model
    auto decomposerPtr_ = decompositionMethod::New
    (
        *distributionDictPtr_,
        mesh().name()
    );


    // Create new decomposition distribution, new processor assignment per cell
    const labelList distribution
    (
        decomposerPtr_->decompose
        (
            mesh(),
            tmp<scalarField>::New()  // cell weights - redundant
        )
    );


    // Create mesh distribution engine
    fvMeshDistribute distributor(mesh());


    // Do actual sending/receiving of mesh
    autoPtr<mapDistributePolyMesh> map = distributor.distribute(distribution);

    if (isA<dynamicRefineFvMesh>(mesh()))
    {
        auto& refineMesh = dynamicCast<dynamicRefineFvMesh>(mesh());

        hexRef8& meshCutter = refineMesh.meshCutter();

        meshCutter.distribute(map());
    }

    // Update motion solver state (points0_ and zone-based pointIDs_)
    // after redistribution, if the mesh uses a motion solver.
    if (isA<dynamicMotionSolverFvMesh>(mesh()))
    {
        auto& dmsMesh = dynamicCast<dynamicMotionSolverFvMesh>(mesh());
        dmsMesh.distribute(map());
    }

    DebugInFunction
        << "Number of cells (after balance): " << mesh().nCells() << endl;

    // What am I going to do with the 'map'?
    return true;
}


// ************************************************************************* //
