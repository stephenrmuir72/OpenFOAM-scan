solidThermo& thermo = thermos[i];

tmp<volScalarField> trho = thermo.rho();
const volScalarField& rho = trho();

tmp<volScalarField> tcp = thermo.Cp();
const volScalarField& cp = tcp();

tmp<volSymmTensorField> taniAlpha;
if (!thermo.isotropic())
{
    volSymmTensorField& aniAlpha = aniAlphas[i];
    tmp<volVectorField> tkappaByCp = thermo.Kappa()/cp;
    const coordinateSystem& coodSys = coordinates[i];

    aniAlpha.primitiveFieldRef() =
        coodSys.transformPrincipal
        (
            mesh.cellCentres(),
            tkappaByCp()
        );

    aniAlpha.correctBoundaryConditions();

    taniAlpha = tmp<volSymmTensorField>
    (
        new volSymmTensorField(aniAlpha)
    );
}


volScalarField& h = thermo.he();

const volScalarField& betav = betavSolid[i];

tmp<surfaceScalarField> talphaf;
tmp<surfaceScalarField> tkappaf;
if (thermo.isZoneMixture() && !thermo.isotropic())
{
    auto coeffs = calcHarmonicZoneCoeffs(mesh, taniAlpha(), cp);

    talphaf = fvc::interpolate(betav) * coeffs.alphaf;
    tkappaf = fvc::interpolate(betav) * coeffs.kappaf;
}

fv::options& fvOptions = solidHeatSources[i];
