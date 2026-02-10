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
#include "fvOptions.H"
#include <ratio>
#include <chrono>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
    DebugVar(U);
    #include "createPhi.H"
    DebugVar(phi);

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

    {
        // Expressions of fvMatrix

        // tmp<fvMatrix<scalar>> tm0(fvm::laplacian(p));

        // // Print a bit
        // {
        //     const fvMatrix<scalar>& m0 = tm0();
        //     DebugVar(m0.dimensions());
        //     DebugVar(m0.hasDiag());
        //     DebugVar(m0.hasUpper());
        //     DebugVar(m0.hasLower());
        //     forAll(m0.internalCoeffs(), i)
        //     {
        //         DebugVar(i);
        //         if (m0.internalCoeffs().set(i))
        //         {
        //             DebugVar(m0.internalCoeffs()[i]);
        //         }
        //         if (m0.boundaryCoeffs().set(i))
        //         {
        //             DebugVar(m0.boundaryCoeffs()[i]);
        //         }
        //     }
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
         == fvc::Su(fvc::laplacian(p), p)
          // + opts.expr()
        );

        // Do some expression in fvMatrixExpression syntax
        // const auto tlaplacian = fvc::laplacian(p);
        // const UList<scalar>& tl = tlaplacian().internalField();
        // const auto expr =
        // (
        //     // fvm::ddt(p).expr()
        //     fvm::laplacian(p).expr()
        //  == Expression::Su<(fvc::laplacian(p), p)
        // //  + opts.expr()
        // );
// {
//     DebugVar(expr.hasDiag());
//     DebugVar(typeid(expr.diag()).name());
//     DebugVar(expr.diag().size());
//     scalarField diag;
//     expr.diag().evaluate(diag);
//     Pout<< "Expression diag: " << flatOutput(diag) << nl;
// }
// {
//     DebugVar(typeid(expr.source()).name());
//     DebugVar(expr.source().size());
//     scalarField source;
//     expr.source().evaluate(source);
//     Pout<< "Expression source: " << flatOutput(source) << nl;
// }

        // Evaluate expression
        const fvMatrix<scalar> m2
        (
            p2,
            (
                fvm::laplacian(p).expr()
             == Expression::Su(fvc::laplacian(p), p)
            )
        );

DebugVar("here");

// Print a bit
        {
            const fvMatrix<scalar>& m1 = tm1();
            DebugVar(m1.dimensions());
            DebugVar(m1.hasDiag());
            DebugVar(m1.hasUpper());
            DebugVar(m1.hasLower());
            forAll(m1.internalCoeffs(), i)
            {
                DebugVar(i);
                if (m1.internalCoeffs().set(i))
                {
                    DebugVar(m1.internalCoeffs()[i]);
                }
                if (m1.boundaryCoeffs().set(i))
                {
                    DebugVar(m1.boundaryCoeffs()[i]);
                }
            }
        }
        // Print a bit
        {
            DebugVar(m2.dimensions());
            DebugVar(m2.hasDiag());
            DebugVar(m2.hasUpper());
            DebugVar(m2.hasLower());
            forAll(m2.internalCoeffs(), i)
            {
                DebugVar(i);
                if (m2.internalCoeffs().set(i))
                {
                    DebugVar(m2.internalCoeffs()[i]);
                }
                if (m2.boundaryCoeffs().set(i))
                {
                    DebugVar(m2.boundaryCoeffs()[i]);
                }
            }
        }
    }

    return 0;
}


// ************************************************************************* //
