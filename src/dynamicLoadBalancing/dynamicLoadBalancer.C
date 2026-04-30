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

#include "dynamicLoadBalancer.H"
#include "balanceModel.H"
#include "performanceMetric.H"
#include "profilingPstream.H"
#include "Function1.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicLoadBalancer, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::dynamicLoadBalancer::isBalancingDue() const
{
    // Requires that  isImbalanced() has been called before

    const scalar t = mesh_.time().value();
    const label timeIndex = mesh_.time().timeIndex();
    const label freq = max(label(1), frequencyPtr_().value(t));

    DebugInFunction
        << "t: " << t << nl
        << "timeIndex: " << timeIndex << nl
        << "freq: " << freq << nl
        << endl;

    bool skipBalance = false;
    if (skipBalancePtr_)
    {
        skipBalance = bool(skipBalancePtr_().value(t));
    }

    bool isBalancingDue =
        (
            UPstream::nProcs() > 1
         && timeIndex > 1
         && curTimeIndex_ != timeIndex
         && (freq <= 1 || (timeIndex % freq == 0))
         && !skipBalance
        );

    if (isBalancingDue)
    {
        if (checkChronicImbalance_)
        {
            const scalar tol = recomputeImbalanceTolerance_;
            const scalar backoffPeriod = backoffPeriod_;

            // Already reduced across processors
            const scalar currImbl = performanceMetricPtr_->imbalanceValue();
            const scalar prevImbl = performanceMetricPtr_->prevImbalanceValue();

            if (mag(currImbl - prevImbl) < tol)  // this clashes with imbalanceThreshold
            {
                if ((t - lastBalanceTime_) < backoffPeriod)
                {
                    Info<< "    Chronic imbalance: skipping balancing" << endl;

                    isBalancingDue = false;
                }
            }
        }
    }

    DebugInFunction
        << "isBalancingDue: " << isBalancingDue << endl;

    return isBalancingDue;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicLoadBalancer::dynamicLoadBalancer
(
    fvMesh& mesh
)
:
    mesh_(mesh),
    dict_
    (
        IOdictionary::readContents
        (
            IOobject
            (
                "dynamicLoadBalanceDict",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ
            )
        )
    ),
    balanceModelPtr_(balanceModel::New(mesh, dict_)),
    performanceMetricPtr_(performanceMetric::New(mesh, dict_)),
    frequencyPtr_(nullptr),
    skipBalancePtr_(nullptr),
    recomputeImbalanceTolerance_(0),
    backoffPeriod_(0),
    checkChronicImbalance_(false),
    lastBalanceTime_(0.0),
    curTimeIndex_(-1)
{
    profilingPstream::enable();
    read(dict_);

    if (!UPstream::parRun())
    {
        WarningInFunction
            << "Dynamic load balancing requires running in parallel."
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicLoadBalancer::~dynamicLoadBalancer()
{
    profilingPstream::disable();
}  // Various models were forward declared


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicLoadBalancer::read(const dictionary& dict)
{
    if
    (
        !balanceModelPtr_->read(dict)
     || !performanceMetricPtr_->read(dict)
    )
    {
        return false;
    }

    skipBalancePtr_ =
        Function1<label>::NewIfPresent("skipBalance", dict, &mesh_);


    if
    (
        const auto& chronic = dict.subOrEmptyDict("chronicImbalanceSettings");
       !chronic.empty()
    )
    {
        recomputeImbalanceTolerance_ =
            chronic.getScalar("recomputeImbalanceTolerance");
        backoffPeriod_ = chronic.getScalar("backoffPeriod");
        checkChronicImbalance_ = true;
    }

    frequencyPtr_ = Function1<label>::New("frequency", dict, &mesh_);

    bool ok = true;
    ok = ok && (frequencyPtr_);

    return ok;
}


bool Foam::dynamicLoadBalancer::balance()
{
    const bool imbalanced = performanceMetricPtr_->isImbalanced();
    const bool due = isBalancingDue();
    const bool anyDueImbalanced = returnReduce(due && imbalanced, orOp<bool>());

    if (!anyDueImbalanced)
    {
        DebugInFunction << "Dynamic load balancing is not active." << endl;

        // Snapshot current imbalance for the next-step chronic check.
        performanceMetricPtr_->updatePrevImbalanceValue();
        return false;
    }

    curTimeIndex_ = mesh_.time().timeIndex();
    lastBalanceTime_ = mesh_.time().value();

    DebugInFunction << "Dynamic load balancing will be performed." << endl;

    bool balanced =
        balanceModelPtr_->balance(performanceMetricPtr_->processorWeights());

    DebugInFunction << "Dynamic load balancing is performed." << endl;

    performanceMetricPtr_->resetTimeCounters();

    // Snapshot current imbalance for the next-step chronic check, after
    // the chronic comparison in isBalancingDue() has already consumed the
    // previous value.
    performanceMetricPtr_->updatePrevImbalanceValue();

    return balanced;
}


// ************************************************************************* //
