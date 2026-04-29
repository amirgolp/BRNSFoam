// Recompute Kozeny-Carman permeability + Darcy drag from current eps.
// Drag is forced to zero in fluid cells (Solid==0) so the DBS reduces to
// Stokes/NS in the resolved-fluid regions.

volScalarField k0
(
    "k0",
    K*pow3(eps)/max(sqr(1.0 - eps), dimensionedScalar("small", dimless, SMALL))
);

volScalarField M    ("M",    k0/muf);
volScalarField Drag ("Drag", 1.0/M);

forAll(Drag, celli)
{
    if (Solid[celli] == 0)
    {
        Drag[celli] = 0.0;
    }
}

// Refresh porosity-weighted mass flux used by UEqn
rhofPhi      = rhof*phi;
rhoPhiByEps  = rhofPhi*linearInterpolate(1.0/(eps + dimensionedScalar("small", dimless, SMALL)));
