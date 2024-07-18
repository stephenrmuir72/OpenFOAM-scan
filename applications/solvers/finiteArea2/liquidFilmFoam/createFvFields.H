// Volume-to surface mapping object
const word patchName("film");
const label patchi = mesh.boundaryMesh().findPatchID(patchName);
const auto& pp = mesh.boundaryMesh()[patchi];
labelList meshFaces(pp.size());
forAll(pp, i)
{
    meshFaces[i] = pp.start()+i;
}
const volSurfaceMapping vsm(aMesh, meshFaces);

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh.thisDb(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimVelocity, Zero)
);


volScalarField H
(
    IOobject
    (
        "H",
        runTime.timeName(),
        mesh.thisDb(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength, Zero)
);
