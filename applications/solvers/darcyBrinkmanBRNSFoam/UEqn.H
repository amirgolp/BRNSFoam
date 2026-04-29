fvVectorMatrix UEqn
(
    // Inertial terms (porosity-scaled)
    (1.0/eps)*(fvm::ddt(rhof, U) + fvm::div(rhoPhiByEps, U))

    // Brinkman viscous term (porosity-scaled)
  - (1.0/eps)*fvm::laplacian(muf, U)

    // Implicit Darcy drag (zero in resolved-fluid cells)
  + fvm::Sp(Drag, U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));
    fvOptions.correct(U);
}
