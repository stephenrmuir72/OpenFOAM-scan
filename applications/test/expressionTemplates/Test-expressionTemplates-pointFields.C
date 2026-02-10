/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 M. Janssens
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
    Test-expressionTemplates

\*---------------------------------------------------------------------------*/

#include "fvMatrixExpression.H"
#include "Time.H"
#include "argList.H"
#include "polyMesh.H"
#include "pointMesh.H"
#include "pointFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createPolyMesh.H"


    {
        // Multiplication
        Expression::ListWrap<scalar> expr1({1.0, 2.0, 3.0});
        Expression::ListWrap<vector> expr2
        ({
            vector(1,2,3),
            vector(4,5,6),
            vector(7,8,9)
        });
        auto expr3 = expr1*expr2;
        typedef decltype(expr1)::value_type expr1Type;
        typedef decltype(expr2)::value_type expr2Type;
        //typedef decltype(expr3)::value_type resultType;
        typedef decltype(std::declval<expr1Type>() * std::declval<expr2Type>()) resultType;
        DebugVar(typeid(resultType).name());

        Field<resultType> result(expr3);
        DebugVar(result);
    }

    {
        const scalarField vals1({1.0, 2.0, 3.0});
        const vectorField vals2
        ({
            vector(1,2,3),
            vector(4,5,6),
            vector(7,8,9)
        });
        typedef typename outerProduct
        <
            decltype(vals1)::value_type,
            decltype(vals2)::value_type
        >::type resultType;

    }

    {
        // Use wrappers

        typedef Expression::ListConstRefWrap<scalar> expr;
        typedef Expression::Default<expr> defaultExpr;
        scalarField vals0({1.0, 2.0, 3.0});
        // const expr wvals0(vals0);
        const defaultExpr dvals0(5, -1, vals0);
        DebugVar(dvals0.size());
        DebugVar(dvals0[3]);

        scalarField vals1;
        Expression::ListRefWrap<scalar> wvals1(dvals0.size(), vals1);
        wvals1 = dvals0 + dvals0 + sin(dvals0);
        DebugVar(vals1);
        DebugVar(dvals0.size());
        DebugVar(dvals0[3]);
        // return 0;
    }




    const pointMesh& pMesh = pointMesh::New(mesh);

    // Field, dimensionedScalar as List expression
    {
        scalarField vals({1.0, 2.0, 3.0});
        dimensionedScalar d(dimless, 4.0);

        scalarField result(d.expr(vals.size())*vals.expr());

        Pout<< "result:" << result << endl;
    }

    Info<< "Reading field p\n" << endl;
    pointScalarField p
    (
        IOobject
        (
            "pointDisplacement",
            runTime.timeName(),
            pMesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh
    );

    pointScalarField result
    (
        IOobject
        (
            "result",
            runTime.timeName(),
            pMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        pMesh,
        dimensionedScalar("zero", Foam::sqr(p().dimensions()), 0)
    );

    {
        // Construct value + dimensions with correct sizing
        const auto oneField(dimensionedScalar(dimless, 1.0).expr(p));

        // Make expression
        auto expression = oneField * p.expr() + p.expr();

        // Combine expressions
        auto newExpression = sqr(expression);

        // Assign values
        result = newExpression;

        DebugVar(result);
    }

    {
        // Construct from expression
        pointScalarField result2("result2", pMesh, sqrt(p.expr()));

        DebugVar(result2);
    }



    return 0;
}


// ************************************************************************* //
