{
    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
     ==
        fvOptions(rho, U)
    );

    #include "UEqnAddPorosity.H"

    UEqn.relax();
    fvOptions.constrain(UEqn);

    while (pimple.correct())
    {
        volScalarField rAU("rAU", 1.0/UEqn.A());
        surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
        volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            fvc::flux(HbyA)
          + MRF.zeroFilter(fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi, Uf))
        );
        MRF.makeRelative(phiHbyA);

        if (p_rgh.needReference())
        {
            fvc::makeRelative(phiHbyA, U);
            adjustPhi(phiHbyA, U, p_rgh);
            fvc::makeAbsolute(phiHbyA, U);
        }

        surfaceScalarField phig
        (
            (
                mixture.surfaceTensionForce()
              - ghf*fvc::snGrad(rho)
            )*rAUf*mesh.magSf()
        );

        phiHbyA += phig;

        constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix p_rghEqn
            (
                fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
            );

            p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

            p_rghEqn.solve(p_rgh.select(pimple.finalInnerIter()));

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA - p_rghEqn.flux();

                p_rgh.relax();

                U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
                U.correctBoundaryConditions();
                fvOptions.correct(U);
            }
        }

        #include "continuityErrs.H"

        fvc::correctUf(Uf, U, phi);
        fvc::makeRelative(phi, U);

        p == p_rgh + rho*gh;

        if (p_rgh.needReference())
        {
            p += dimensionedScalar
            (
                "p",
                p.dimensions(),
                pRefValue - getRefCellValue(p, pRefCell)
            );
            p_rgh = p - rho*gh;
        }
    }
}

