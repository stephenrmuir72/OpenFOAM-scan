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

Description
    Testing std::execution::par_unseq policy

\*---------------------------------------------------------------------------*/

#include "UList.H"
#include "argList.H"
#include "Time.H"
#include "OSspecific.H"
#include "FieldOps.H"
#include "clockTime.H"
#include "Random.H"
#include "vectorField.H"
#include "ListPolicy.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:
int main(int argc, char *argv[])
{
    constexpr label defaultCount = 10000000;

    argList::noCheckProcessorDirectories();

    argList::addOption("count", "number of 'cells'");

    #include "setRootCase.H"

    const label nCells = args.getOrDefault<label>("count", defaultCount);
    const label nIntFaces = (nCells-1);

    List<label> cellVals(nCells);

    Info<< "offload_enabled: "
        << Switch(ListPolicy::offload_enabled) << nl << nl;

    Info<< "Testing fill of large list (" << nCells << ")" << nl;
    {
        clockTime timing;
        for (label i = 0; i < 1000; ++i)
        {
            // Do some work to ensure we have a reasonable amount of time in the
            // test to allow for attaching a profiler, etc.
            std::fill_n(cellVals.begin(), cellVals.size(), 123);
        }
        Info<< "Sequential fill took " << timing.timeIncrement() << "s" << nl;
        for (label i = 0; i < 1000; ++i)
        {
            // Do some work to ensure we have a reasonable amount of time in the
            // test to allow for attaching a profiler, etc.
            std::fill_n
            (
                #ifdef __cpp_lib_execution
                std::execution::par_unseq,
                #else
                std::execution::seq,
                #endif
                cellVals.begin(),
                cellVals.size(),
                123
            );
        }
        Info<< "Parallel fill took " << timing.timeIncrement() << "s" << nl;
        for (label i = 0; i < 1000; ++i)
        {
            // Do some work to ensure we have a reasonable amount of time in the
            // test to allow for attaching a profiler, etc.
            std::fill_n(cellVals.begin(), cellVals.size(), 123);
        }
        Info<< "Sequential fill took " << timing.timeIncrement() << "s" << nl;
    }
    Info<< nl << endl;


    // Construct lduAddressing for all internal faces
    // to allow testing of interpolation with non-trivial addressing

    labelList P(nIntFaces);
    labelList N(P.size());
    forAll(P, facei)
    {
        P[facei] = facei;
        N[facei] = P[facei]+1;
    }
    // Interpolation weights for linear interpolation
    const scalarField lambda(P.size(), 0.5);
    // Faked face area vectors for testing
    const vectorField Sfi(P.size(), vector(1, 0, 0));



    Info<< "Testing dotInterpolate (cell to face) using lambda" << nl;
    {
        Field<scalar> sfi(P.size());
        Field<vector> vfi(nCells);
        Random::uniformGeneratorOp<vector> gen(123);
        FieldOps::assign(vfi, vfi, gen);

        const labelRange faceRange(P.size());

        auto dotInterpolate = [&](const label fi)
        {
            sfi[fi] =
                Sfi[fi]
              & (lambda[fi]*(vfi[P[fi]] - vfi[N[fi]]) + vfi[N[fi]]);
        };

        clockTime timing;
        for (label i = 0; i < 1000; ++i)
        {
            std::for_each
            (
                faceRange.begin(),
                faceRange.end(),
                dotInterpolate
            );
        }
        Info<< "Sequential interpolation took " << timing.timeIncrement()
            << "s" << nl;
        for (label i = 0; i < 1000; ++i)
        {
            // Do some work to ensure we have a reasonable amount of time in the
            // test to allow for attaching a profiler, etc.
            std::for_each
            (
                #ifdef __cpp_lib_execution
                std::execution::par_unseq,
                #else
                std::execution::seq,
                #endif
                faceRange.begin(),
                faceRange.end(),
                dotInterpolate
            );
        }
        Info<< "Parallel interpolation took " << timing.timeIncrement()
            << "s" << nl;
    }

    // Sanity check: ListPolicy::for_each_n, transform_n
    {
        const label n = 1000;
        labelList a(n, 0), b(n, 0), c(n, 0);

        ListPolicy::for_each_n(n, [&](label i) { a[i] = i; });
        ListPolicy::transform_n
        (
            n, b.begin(), [](label i) { return 2*i; }
        );
        ListPolicy::transform_n
        (
            n, a.cbegin(), c.begin(),
            [](label i, label x) { return i + x; }
        );

        bool ok = true;
        forAll(a, i)
        {
            ok = ok && (a[i] == i) && (b[i] == 2*i) && (c[i] == 2*i);
        }
        Info<< nl
            << "ListPolicy for_each_n/transform_n: "
            << (ok ? "ok" : "FAILED") << nl;
    }

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
