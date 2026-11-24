/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 IH-Cantabria
    Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "streamFunctionFentonWaveModel.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(streamFunctionFenton, 0);
    addToRunTimeSelectionTable
    (
        waveModel,
        streamFunctionFenton,
        patch
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::waveModels::streamFunctionFenton::eta
(
    const scalar h,
    const scalar kx,
    const scalar ky,
    const scalar T,
    const scalar x,
    const scalar y,
    const scalar omega,
    const scalar t,
    const scalar phase
) const
{
   const scalar k = sqrt(kx*kx + ky*ky);
   scalar strfnAux = 0.0;
   
   // Ejs are dimensionless, physical eta = (1/k) * Sum(E_j * cos(...))
   forAll(Ejs_, iterSF)
   {
        strfnAux +=
            Ejs_[iterSF]*cos((iterSF + 1) *(kx*x + ky*y - omega*t + phase));
   }

   return (1.0/k)*strfnAux;
}


Foam::vector Foam::waveModels::streamFunctionFenton::Uf
(
    const scalar h,
    const scalar kx,
    const scalar ky,
    const scalar T,
    const scalar x,
    const scalar y,
    const scalar omega,
    const scalar t,
    const scalar phase,
    const scalar z
) const
{
    const scalar k = sqrt(kx*kx + ky*ky);
    const scalar phaseTot = kx*x + ky*y - omega*t + phase;

    scalar u = 0.0;
    scalar w = 0.0;

    // Bjs are dimensionless.
    forAll(Bjs_, iterSF2)
    {
        scalar j = iterSF2 + 1.0;
        u +=
            j*Bjs_[iterSF2]*cosh(j*k*z)
           /cosh(j*k*h)*cos(j*phaseTot);

        w +=
            j*Bjs_[iterSF2]*sinh(j*k*z)
           /cosh(j*k*h)*sin(j*phaseTot);
    }

    // Convert dimensionless velocity sum to physical
    // u_phys = C - uMean + sqrt(g/k) * u_dim
    // Note: waveLength_/T is phase speed C
    
    u = waveLength_/T - uMean_ + sqrt(mag(g_)/k)*u;
    w = sqrt(mag(g_)/k)*w;

    scalar v = u*sin(waveAngle_);
    u *= cos(waveAngle_);

    return vector(u, v, w);
}

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::waveModels::streamFunctionFenton::evaluateFunction
(
    const scalarList& x, 
    scalarList& f, 
    const scalar kd, 
    const scalar kH
) const
{
    // Mapping variables from flat vector x
    // 0: tau, 1: celerity, 2: ce, 3: uBar, 4: q, 5: rkg
    scalar tau = x[0];
    scalar celerity = x[1];
    scalar ce = x[2];
    scalar uBar = x[3];
    scalar q = x[4];
    scalar rkg = x[5];

    // Offsets
    label etaOffset = 6;
    label bOffset = 6 + (nModes_ + 1);

    // Helper accessors
    auto getEta = [&](label m) { return x[etaOffset + m]; };
    auto getB = [&](label j) { return x[bOffset + (j-1)]; }; // j 1-based

    // 1. Scalar Constraints
    f[0] = celerity*tau - mathematical::twoPi;
    f[1] = ce + uBar - celerity;
    f[2] = uBar - celerity - q/kd;
    
    scalar sumEtaInner = 0.0;
    for(label j=1; j<nModes_; ++j) sumEtaInner += getEta(j);
    f[3] = getEta(0) + getEta(nModes_) + 2.0*sumEtaInner;
    
    f[4] = getEta(0) - getEta(nModes_) - kH;

    // 2. Kinematic Boundary Condition (BBC1)
    // EqKinematic for m = 0 to N
    for (label m = 0; m <= nModes_; ++m)
    {
        scalar termSum = 0.0;
        for (label j = 1; j <= nModes_; ++j)
        {
            scalar arg = j*(getEta(m) + kd);
            scalar cosTerm = cos((j*m*mathematical::pi)/scalar(nModes_));
            termSum += getB(j) * (sinh(arg)/cosh(j*kd)) * cosTerm;
        }
        f[5 + m] = -getEta(m)*uBar + termSum - q;
    }

    // 3. Dynamic Boundary Condition (BBC2)
    // EqDynamic for m = 0 to N
    label dynamicOffset = 5 + (nModes_ + 1);
    for (label m = 0; m <= nModes_; ++m)
    {
        scalar sumU = 0.0;
        scalar sumV = 0.0;
        
        for (label j = 1; j <= nModes_; ++j)
        {
            scalar arg = j*(getEta(m) + kd);
            scalar angle = (j*m*mathematical::pi)/scalar(nModes_);
            scalar coshTerm = cosh(arg)/cosh(j*kd);
            scalar sinhTerm = sinh(arg)/cosh(j*kd);
            
            sumU += j * getB(j) * coshTerm * cos(angle);
            sumV += j * getB(j) * sinhTerm * sin(angle);
        }

        scalar uTerm = (-uBar + sumU);
        scalar termSq = uTerm*uTerm + sumV*sumV;
        
        f[dynamicOffset + m] = getEta(m) - rkg + 0.5 * termSq;
    }
}

void Foam::waveModels::streamFunctionFenton::calculateCoeffs()
{
    Info<< "Calculating Stream Function coefficients (Fenton 1988)..." << endl;

    Info << " waterDepthRef_ = " << waterDepthRef_ 
         << ", waveHeight_ = " << waveHeight_ 
         << ", waveLength_ = " << waveLength_ << endl;

    const scalar waveK = mathematical::twoPi/waveLength_;
    const scalar kd = waveK * waterDepthRef_;
    const scalar kH = waveK * waveHeight_;

    // Total variables: 6 scalars + (N+1) etas + N Bs
    const label numVars = 6 + (nModes_ + 1) + nModes_;
    
    scalarList x(numVars, 0.0);
    scalarList f(numVars, 0.0);

    // --- 1. Initial Guess (Linear Wave Theory) ---
    Info << " Step 1. Initializing Stream Function solver..." << endl;
    scalar tanhKd = tanh(kd);
    scalar sqrtTanhKd = sqrt(tanhKd);

    // Scalars
    x[0] = mathematical::twoPi / sqrtTanhKd; // tau
    x[1] = sqrtTanhKd;                       // celerity
    x[2] = 0.0;                              // ce
    x[3] = sqrtTanhKd;                       // uBar
    x[4] = 0.0;                              // q
    x[5] = tanhKd / 2.0;                     // rkg

    // Etas (Cosine profile)
    label etaOffset = 6;
    for(label i=0; i<=nModes_; ++i)
    {
        x[etaOffset + i] = (kH/2.0) * cos(i*mathematical::pi/scalar(nModes_));
    }

    // Bs (Only B1 is non-zero initially)
    label bOffset = 6 + (nModes_ + 1);
    x[bOffset] = (kH/2.0) / sqrtTanhKd;

    // --- 2. Finite Difference Jacobian Initialization ---
    Info << " Step 2. Initializing Jacobian matrix..." << endl;
    RectangularMatrix<scalar> J(numVars, numVars, 0.0);
    scalar eps = 1e-7;

    evaluateFunction(x, f, kd, kH);
    scalarList f0 = f; 

    for(label j=0; j<numVars; ++j)
    {
        scalar oldVal = x[j];
        x[j] += eps;
        
        scalarList fPert(numVars, 0.0); // [Fix] Initialize to 0
        evaluateFunction(x, fPert, kd, kH);
        
        for(label i=0; i<numVars; ++i)
        {
            J(i, j) = (fPert[i] - f0[i]) / eps;
        }
        
        x[j] = oldVal; 
    }

    // --- 3. Broyden's Method Loop ---
    Info << " Step 3. Solving nonlinear system..." << endl;
    label maxIter = 100; // Increase max iter
    scalar tol = 1e-8;
    scalar relaxation = 0.5; // [Fix] Under-relaxation factor
    
    for(label iter=0; iter<maxIter; ++iter)
    {
        evaluateFunction(x, f, kd, kH);
        
        scalar resNorm = 0.0;
        forAll(f, i) resNorm += f[i]*f[i];
        resNorm = sqrt(resNorm);

        // Print residual every step to monitor convergence
        if (iter % 5 == 0 || iter < 5)
        {
             Info << "  Iter " << iter << ", Residual = " << resNorm << endl;
        }

        if (resNorm < tol)
        {
            Info<< "StreamFunction converged in " << iter << " iterations. Residual: " << resNorm << endl;
            break;
        }
        
        // [Fix] Safety Check for NaN/Inf
        if (!std::isfinite(resNorm))
        {
            FatalErrorInFunction 
                << "Solver diverged (NaN/Inf detected) at iter " << iter 
                << ". Try reducing relaxation or wave height." 
                << exit(FatalError);
        }

        // Solve J * s = -f
        scalarList deltaX(numVars, 0.0);
        scalarList minusF = -1.0 * f;

        {
            QRMatrix<RectangularMatrix<scalar>> QRM
            (
                J,
                QRMatrix<RectangularMatrix<scalar>>::modes::FULL,
                QRMatrix<RectangularMatrix<scalar>>::outputs::BOTH_QR
            );
            
            const RectangularMatrix<scalar>& Q = QRM.Q();
            const RectangularMatrix<scalar>& R = QRM.R();

            // Calculate y = Q^T * (-f)
            scalarList y(numVars, 0.0);
            for(label i=0; i<numVars; ++i)
            {
                for(label k=0; k<numVars; ++k)
                {
                    y[i] += Q(k, i) * minusF[k];
                }
            }

            // Back substitution R * deltaX = y
            for (label i = numVars - 1; i >= 0; --i)
            {
                scalar sum = 0.0;
                for (label j = i + 1; j < numVars; ++j)
                {
                    sum += R(i, j) * deltaX[j];
                }
                
                // [Fix] More robust check for singularity
                if (mag(R(i, i)) > 1e-12) 
                {
                    deltaX[i] = (y[i] - sum) / R(i, i);
                }
                else
                {
                    // Handle singular modes (e.g. unconstrained current)
                    deltaX[i] = 0.0; 
                }
            }
        }

        // Broyden Update
        scalar sDotS = 0.0;
        forAll(deltaX, i) 
        {
            sDotS += deltaX[i]*deltaX[i];
        }
        
        // [Fix] Limit step size if deltaX is huge
        if (sDotS > 100.0) 
        {
            WarningInFunction << "Large step detected, clipping." << endl;
            scalar scale = 10.0 / sqrt(sDotS);
            forAll(deltaX, i) deltaX[i] *= scale;
            sDotS = 0.0; 
            forAll(deltaX, i) sDotS += deltaX[i]*deltaX[i];
        }

        // Update x with Relaxation
        scalarList xNew = x;
        forAll(xNew, i) 
        {
            xNew[i] += deltaX[i] * relaxation; // [Fix] Apply relaxation
        }

        // Compute new F
        scalarList fNew(numVars);
        evaluateFunction(xNew, fNew, kd, kH);

        // Update Jacobian
        if (sDotS > 1e-20) // [Fix] Protect against div by zero
        {
             for(label i=0; i<numVars; ++i)
             {
                 for(label j=0; j<numVars; ++j)
                 {
                     // Broyden update: J += (df * dx^T) / (dx^T * dx)
                     // Here we use simplified update suitable when linear solve was accurate
                     // J_new = J_old + (fNew * deltaX^T) / (deltaX^T * deltaX) * relax?
                     // Standard Broyden uses full step s. 
                     // We use the 'deltaX' (Newton step) as the direction 's'.
                     
                     // Warning: Broyden update is sensitive. 
                     // Using fNew directly approximates (y - J*s) if J*s = -fOld.
                     J(i, j) += (fNew[i] * deltaX[j]) / sDotS;
                 }
             }
        }

        x = xNew;
        
        if (iter == maxIter - 1)
        {
            WarningInFunction << "Maximum iterations reached. Residual: " << resNorm << endl;
        }
    }

    // --- 4. Store Results ---
    Info << " Step 4. Storing calculated coefficients..." << endl;
    
    // scalar ce = x[2];
    uMean_ = (x[1]+x[2]) * sqrt(mag(g_)/waveK); 

    Info << "  T old = " << wavePeriod_ << endl;
    wavePeriod_ = waveLength_ / (x[1] * sqrt(mag(g_)/waveK)); 

    Info << "  T new = " << wavePeriod_ << endl;

    Bjs_.setSize(nModes_);
    for(label j=0; j<nModes_; ++j)
    {
        Bjs_[j] = x[bOffset + j];
    }

    Ejs_.setSize(nModes_);
    scalarList calculatedEtas(nModes_ + 1);
    for(label i=0; i<=nModes_; ++i) calculatedEtas[i] = x[etaOffset + i];

    for(label j=1; j<=nModes_; ++j)
    {
        scalar sum = 0.0;
        sum += 0.5 * (calculatedEtas[0] + calculatedEtas[nModes_] * cos(j * mathematical::pi));
        
        for(label m=1; m<nModes_; ++m)
        {
             sum += calculatedEtas[m] * cos(m * j * mathematical::pi / scalar(nModes_));
        }
        
        Ejs_[j-1] = (2.0/scalar(nModes_)) * sum;
    }
}

void Foam::waveModels::streamFunctionFenton::setLevel
(
    const scalar t,
    const scalar tCoeff,
    scalarField& level
) const
{
    const scalar waveOmega = mathematical::twoPi/wavePeriod_;
    const scalar waveK = mathematical::twoPi/waveLength_;
    const scalar waveKx = waveK*cos(waveAngle_);
    const scalar waveKy = waveK*sin(waveAngle_);

    forAll(level, paddlei)
    {
        const scalar eta =
            this->eta
            (
                waterDepthRef_,
                waveKx,
                waveKy,
                wavePeriod_,
                xPaddle_[paddlei],
                yPaddle_[paddlei],
                waveOmega,
                t,
                wavePhase_
            );

        level[paddlei] = waterDepthRef_ + tCoeff*eta;
    }
}


void Foam::waveModels::streamFunctionFenton::setVelocity
(
    const scalar t,
    const scalar tCoeff,
    const scalarField& level
)
{
    const scalar waveOmega = mathematical::twoPi/wavePeriod_;
    const scalar waveK = mathematical::twoPi/waveLength_;
    const scalar waveKx = waveK*cos(waveAngle_);
    const scalar waveKy = waveK*sin(waveAngle_);

    forAll(U_, facei)
    {
        // Fraction of geometry represented by paddle - to be set
        scalar fraction = 1;

        // Height - to be set
        scalar z = 0;

        setPaddlePropeties(level, facei, fraction, z);

        if (fraction > 0)
        {
            const label paddlei = faceToPaddle_[facei];

            const vector Uf = this->Uf
            (
                waterDepthRef_,
                waveKx,
                waveKy,
                wavePeriod_,
                xPaddle_[paddlei],
                yPaddle_[paddlei],
                waveOmega,
                t,
                wavePhase_,
                z
            );

            U_[facei] = fraction*Uf*tCoeff;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::streamFunctionFenton::streamFunctionFenton
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    regularWaveModel(dict, mesh, patch, false),
    nModes_(10), // Default order
    uMean_(0),
    Bjs_(),
    Ejs_()
{
    if (readFields)
    {
        readDict(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::streamFunctionFenton::readDict(const dictionary& overrideDict)
{
    if (regularWaveModel::readDict(overrideDict))
    {
        overrideDict.readEntry("waveLength", waveLength_);

        if (overrideDict.found("Bjs") && overrideDict.found("Ejs"))
        {
            Info << "Reading Stream Function coefficients from dictionary." << endl;
            overrideDict.readEntry("uMean", uMean_);
            overrideDict.readEntry("Bjs", Bjs_);
            overrideDict.readEntry("Ejs", Ejs_);
        }
        else
        {
            Info << "Calculating Stream Function coefficients..." << endl;
            overrideDict.readIfPresent("nModes", nModes_);
            // Bjs_, Ejs_, uMean, and wavePeriod_ will be calculated.
            calculateCoeffs();
        }

        return true;
    }

    return false;
}


void Foam::waveModels::streamFunctionFenton::info(Ostream& os) const
{
    regularWaveModel::info(os);

    os  << "    Stream function order (N) : " << nModes_ << nl
        << "    uMean : " << uMean_ << nl
        << "    Stream function wavelength : " << waveLength_ << nl
        << "    Bj coefficients : " << Bjs_ << nl
        << "    Ej coefficients : " << Ejs_ << nl;
}


// ************************************************************************* //


