{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; ++nonOrth)
    {
        fvScalarMatrix hEqn
        (
            thermo.isotropic()
          ? -fvm::laplacian(betav*thermo.alpha(), h, "laplacian(alpha,h)")
          : -fvm::laplacian(betav*taniAlpha(), h, "laplacian(alpha,h)")
          ==
            fvOptions(rho, h)
        );

        if (thermo.isZoneMixture())
        {
            if (thermo.isotropic())
            {
                hEqn += fvc::laplacian(betav*thermo.alpha(), h, "laplacian(alpha,h)");
                hEqn -= fvc::laplacian(betav*thermo.kappa(), thermo.T(), "laplacian(kappa,T)");
            }
            else
            {
                hEqn += fvc::laplacian(betav*taniAlpha(), h, "laplacian(alpha,h)");
                hEqn -= fvc::laplacian(betav*taniAlpha()*cp, thermo.T(), "laplacian(kappa,T)");
            }
        }

        hEqn.relax();

        fvOptions.constrain(hEqn);

        if (coupled)
        {
            fvMatrixAssemblyPtr->addFvMatrix(hEqn);
        }
        else
        {
            Info<< "\nSolving for solid region "
                << solidRegions[i].name() << endl;

            hEqn.solve();
            fvOptions.correct(h);
            thermo.correct();

            Info<< "Min/max T:" << min(thermo.T()).value() << ' '
                << max(thermo.T()).value() << endl;
        }
    }
}
