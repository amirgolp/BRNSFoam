// Intrinsic permeability scaling field K (mesh-resolution-independent)
Info<< "Reading intrinsic permeability scale K\n" << endl;
volScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Initial Solid-indicator from porosity (eps>0.99 => fluid)
forAll(eps, celli)
{
    Solid[celli] = (eps[celli] > 0.99) ? 0.0 : 1.0;
}
forAll(mesh.boundary(), patchi)
{
    forAll(Solid.boundaryField()[patchi], facei)
    {
        Solid.boundaryFieldRef()[patchi][facei] =
            (eps.boundaryField()[patchi][facei] > 0.99) ? 0.0 : 1.0;
    }
}
