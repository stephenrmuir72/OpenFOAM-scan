if (finalIter)
{
    mesh.data().setFinalIteration(true);
}

{
    const auto* aniAlphaPtr = thermo.Anialpha();

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; ++nonOrth)
    {
        fvScalarMatrix hEqn
        (
            fvm::ddt(betav*rho, h)
          - (
               aniAlphaPtr
             ? fvm::laplacian(betav*(*aniAlphaPtr), h, "laplacian(alpha,h)")
             : fvm::laplacian(betav*thermo.alpha(), h, "laplacian(alpha,h)")
            )
          ==
            fvOptions(rho, h)
        );

        hEqn.relax();

        fvOptions.constrain(hEqn);

        hEqn.solve(h.select(finalIter));

        fvOptions.correct(h);
    }

    thermo.correct();

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}

if (finalIter)
{
    mesh.data().setFinalIteration(false);
}
