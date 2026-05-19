/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    volPointInterpolationTest

\*---------------------------------------------------------------------------*/

#include "ListExpression.H"
#include "fvMatrixExpression.H"
#include "Time.H"
#include "argList.H"
#include "fvMesh.H"
#include "ListExpression.H"
#include "GeometricFieldExpression.H"
#include "fvCFD.H"
#include "EulerDdtScheme.H"
#include "fvOptions.H"
#include <ratio>
#include <chrono>
#include "gaussLaplacianScheme.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*
- fvMatrixExpression expr = Expr::ddt(p) - Expr::div(phi, U) == Expr::Su;
    -> ddt : fvMatrixInternal
    -> templated as min<<ddtU,divU>,Su>
- fvMatrix<scalar>(p, expr);
    - calls expr.evaluate(m)
        - calls InitExpr(m)?
            - calls virtual ddtScheme to in-place populate and add to m
            - so how do we include the negation of the divU contribution?
                - cannot do it as a template instantiation
                - could do it at constructor time - start off with 1.0 and then
                  multiply by -1.0 for divU contribution
                - could do it as a lambda or scalar passed to the InitExpr, but then we
                  lose the ability to use the expression templates for the divU contribution
*/

template<class T>
void print(const std::string& msg, const fvMatrix<T>& m)
{
    Pout<< msg << " psi: " << m.psi().name() << nl
        << "    size:" << m.lduAddr().size()
        << " dimensions:" << m.dimensions()
        << " diag:" << m.hasDiag()
        << " uppper:" << m.hasUpper()
        << " lower:" << m.hasLower()
        << nl;
    if (m.hasDiag())
    {
        Pout<< "    diag:" << flatOutput(m.diag()) << nl;
    }
    Pout<< "    src :" << flatOutput(m.source()) << nl;
    if (m.hasUpper())
    {
        Pout<< "    upper:" << flatOutput(m.upper()) << nl;
    }
    if (m.hasLower())
    {
        Pout<< "    lower:" << flatOutput(m.lower()) << nl;
    }
    forAll(m.internalCoeffs(), i)
    {
        if (m.internalCoeffs().set(i))
        {
            Pout<< "    i:" << i
                << " internal:" << m.internalCoeffs()[i]
                << endl;
        }
        if (m.boundaryCoeffs().set(i))
        {
            Pout<< "    i:" << i
                << " boundary:" << m.boundaryCoeffs()[i]
                << endl;
        }
    }
}

tmp<volScalarField> someFunction(const volScalarField& fld)
{
    return fld*1.0001;
}


using namespace std::chrono;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Fill with random value
    for (label i=0; i<p.internalField().size(); ++i)
    {
        const_cast<DynamicField<scalar>&>(p.primitiveField())[i] = scalar(i);
    }

    // fvm expression
    {
        const fvMatrix<scalar> m3(-fvm::ddt(p));
        // fvMatrix<scalar> m3(p, p.dimensions()*dimVol/dimTime);
        // fvm::ddt<scalar>(m3, 1.0, p);
        print("-fvm::ddt", m3);
        auto expr(-Expression::ddt(p));
        DebugVar(expr.dimensions());
        DebugVar(expr.hasDiag());
        DebugVar(expr.source().size());
        const fvMatrix<scalar> m2(p, expr);
        print("-Expression::ddt", m2);
    }


    if (false)
    {
        // Test ListsRefWrap

        volScalarField p2("p2", p);
        DebugVar(p2.boundaryFieldRef());

        for (auto& pf : p.boundaryFieldRef())
        {
            scalarField newVals(pf.size());
            forAll(pf, i)
            {
                newVals[i] = scalar(i);
            }
            pf == newVals;
        }

        Expression::ListsConstRefWrap<scalarField> expr(p.boundaryField());

        const auto twoA = expr + expr;
        Expression::ListsRefWrap<scalarField> expr2(p2.boundaryFieldRef());
        Pout<< "**before assignment twoA:" << twoA.size()
            << " expr2:" << expr2.size() << endl;

        expr2 = twoA;

        Pout<< "**after assignment twoA:" << twoA.size()
            << " expr2:" << expr2.size() << endl;
        DebugVar(p2.boundaryField());
    }

    if (false)
    {
        surfaceScalarField result
        (
            IOobject
            (
                "result",
                runTime.timeName(),
                mesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(p.dimensions(), 0)
        );

        //- Supplied weights
        //result = Expression::interpolate
        //(
        //    p.expr(),
        //    tweights().expr(),
        //    mesh
        //);

        //- Mesh-based linear weights
        result = Expression::linearInterpolate(p.expr(), mesh);

        DebugVar(result);
    }



    volScalarField p2
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh
    );

    if (false)
    {
        // Expressions of volFields
        volScalarField result
        (
            IOobject
            (
                "result",
                runTime.timeName(),
                mesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(p.dimensions(), 0)
        );

        {
            Pout<< "No expression templates:" << endl;
            const high_resolution_clock::time_point t1 =
                high_resolution_clock::now();
            result = p + p2;
            const high_resolution_clock::time_point t2 =
                high_resolution_clock::now();
            const duration<double> time_span = t2 - t1;
            Pout<< "Operation time:" << time_span.count() << endl;
        }
        {
            Pout<< "With expression templates:" << endl;
            const high_resolution_clock::time_point t1 =
                high_resolution_clock::now();
            result = p.expr() + p2.expr();
            const high_resolution_clock::time_point t2 =
                high_resolution_clock::now();
            const duration<double> time_span = t2 - t1;
            Pout<< "Operation time:" << time_span.count() << endl;
        }
    }
    if (false)
    {
        // Expressions of (tmp)volFields
        volScalarField result
        (
            IOobject
            (
                "result",
                runTime.timeName(),
                mesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(p.dimensions(), 0)
        );
        auto expression = someFunction(p).expr() + someFunction(p).expr();
        result = expression;
        DebugVar(result);
    }

    if (false)
    {
        // Expressions of volFields
        volScalarField result
        (
            "result",
            mesh,
            sqr(p.expr() + p.expr())
        );
        DebugVar(result);
    }

    /*
    if (false)
    {
        // Fill p with some values
        forAll(p, celli)
        {
            p[celli] = celli;
        }
        p.correctBoundaryConditions();

        // Interpolate to surface field
        surfaceScalarField result
        (
            "result",
            mesh,
            Expression::interpolate    //<volExpr, surfaceExpr>
            (
                p.expr(),
                mesh.surfaceInterpolation::weights().expr(),
                mesh
            )
        );
        DebugVar(result);
    }
    */

    const auto tsomeField = fvc::laplacian(p);
    const auto& someField = tsomeField();
    if (false)
    {
        // Expressions of fvMatrix

        // tmp<fvMatrix<scalar>> tm0(fvm::laplacian(p));

        // // Print a bit
        // {
        //     const fvMatrix<scalar>& m0 = tm0();
        //     print("fvm::laplacian", m0);
        // }

        // if (false)
        // {
        //     // Test templated tmp<fvMatrix> wrapper
        //     Expression::fvMatrixConstTmpWrap2<fvScalarMatrix> mWrap(tm0);
        //     DebugVar(mWrap.dimensions());
        //     DebugVar(mWrap.hasDiag());
        //     DebugVar(mWrap.hasUpper());
        //     DebugVar(mWrap.hasLower());
        // }


        // fv::options& fvOptions(fv::options::New(mesh));

        // const auto opts = fvOptions(p);

        // Do some expression in fvMatrix syntax
        const tmp<fvMatrix<scalar>> tm1
        (
            // fvm::ddt(p)
            fvm::laplacian(p)
         == fvc::Su(someField, p)
          // + opts.expr()
        );
DebugVar("here");

        // Evaluate expression
        const fvMatrix<scalar> m2
        (
            p2,
            (
                fvm::laplacian(p).expr()
             == Expression::Su(someField, p)
            )
        );

        // Print a bit
        print("fvm::laplacian,fvc::Su", tm1());
        print("fvm::laplacian,Expression::Su", m2);
    }


DebugVar("----------------------------");
    {
        const fvMatrix<scalar> m3(fvm::ddt(p) - fvm::laplacian(p) == fvc::Su(someField, p));
        // fv::gaussLaplacianScheme<scalar, scalar> scheme(mesh);
        // const fvMatrix<scalar> m3
        // (
        //     scheme.fvmLaplacian(fvc::interpolate(one), p)
        // );
        print("** fv::ddt,fvm::laplacian,fvc::Su", m3);
    }
    {
        fvMatrix<scalar> m2(p2, Expression::ddt(p) - Expression::laplacian(p) == Expression::Su(someField, p));
        //fvm::laplacian<scalar, scalar>(m2, -1.0, p);
        print("**Expression::ddt,laplacian,Su", m2);
    }

Pout<< "*** fvm::div ***" << nl;
    {
        // volVectorField U2
        // (
        //     IOobject
        //     (
        //         "U",
        //         runTime.timeName(),
        //         mesh.thisDb(),
        //         IOobject::MUST_READ,
        //         IOobject::AUTO_WRITE,
        //         IOobject::NO_REGISTER
        //     ),
        //     mesh
        // );
        surfaceScalarField phi
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh.thisDb(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE,
                IOobject::NO_REGISTER
            ),
            fvc::flux(U)
        );

        const fvMatrix<vector> m3(p*fvm::ddt(U));   // - fvm::div(phi, U));
        print("** fvm::ddt,fvm::div", m3);

        // Have operator* for scalarField*Expression
        fvMatrix<vector> m2(U, p.expr()*Expression::ddt(U));   // - Expression::div(phi, U));

        // p*Expression::ddt(U);
        // becomes an object
        // scale<volScalarField, Expression::ddt>
        print("**Expression::ddt,div", m2);
    }

    return 0;
}


// ************************************************************************* //
