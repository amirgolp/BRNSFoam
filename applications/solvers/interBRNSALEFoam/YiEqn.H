{
    // Henry-mapped formulation with consistent advection.
    //
    // With Dbar = alpha + H*(1-alpha), the species concentration is
    // Yi = Dbar*Y2 where Y2 is continuous across the gas-liquid interface
    // at equilibrium. The transport equation in Y2 is
    //
    //   d(Dbar*Y2)/dt + div(phiDbar*Y2) - laplacian(gamma, Y2) = 0
    //
    // with gamma = D1*alpha + H*D2*(1-alpha). The diffusion + Henry-transfer
    // face flux equals -gamma*snGrad(Y2); solving for Y2 makes the
    // interfacial transfer implicit and unconditionally stable.
    //
    // Advective flux phiDbar is built from the mixture flux phi (already
    // relative to mesh motion in the ALE solver) and the MULES phase-1
    // volumetric flux alphaPhi10 so that the discrete identity
    //
    //   d(Dbar)/dt + div(phiDbar) = 0
    //
    // holds (Y2 = const must remain a solution). Since Dbar = H - (H-1)*alpha
    // and div(phi) = 0 for the incompressible mixture, the choice
    //
    //   phiDbar = H*phi + (1 - H)*alphaPhi10
    //
    // gives the required divergence and reduces to phi*Yi at single-phase
    // boundary faces (water-inlet alphaPhi10 = phi -> phiDbar = phi;
    // gas-inlet alphaPhi10 = 0 -> phiDbar = H*phi, and Y2 = Yi/H there, so
    // the advected mass is phi*Yi in either case).
    //
    // Without this term, alpha is advected by the flow but the species
    // equation is purely diffusive, so a fluid parcel leaves its dissolved
    // Yi behind as the interface sweeps through it; with H >> 1 (or << 1)
    // the resulting Dbar*Y2 recovery rescales the residual concentration
    // and produces unbounded drift in Yi where there is no reactive sink.
    //
    // Note: advection makes the matrix asymmetric, so the fvSolution Yi
    // block must use a symmetry-agnostic solver (smoothSolver +
    // symGaussSeidel, or PBiCGStab + diagonal). DIC/DILU will fail.

    const speciesTable& solutionSpecies = speciesMixture.species();

    forAll(solutionSpecies, i)
    {
        volScalarField& Yi = speciesMixture.Y(i);
        const dimensionedScalar H  = speciesMixture.HY(i);
        const dimensionedScalar D1 = speciesMixture.D1Y(i);
        const dimensionedScalar D2 = speciesMixture.D2Y(i);

        volScalarField Dbar
        (
            IOobject
            (
                "Dbar_" + Yi.name(),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            alpha1 + H*(scalar(1) - alpha1)
        );
        Dbar.oldTime() == alpha1.oldTime() + H*(scalar(1) - alpha1.oldTime());

        volScalarField Y2(Yi);
        Y2.rename("Y2_" + Yi.name());
        Y2.primitiveFieldRef() /= Dbar.primitiveField();
        forAll(Y2.boundaryField(), patchi)
        {
            Y2.boundaryFieldRef()[patchi] =
                Yi.boundaryField()[patchi]/Dbar.boundaryField()[patchi];
        }
        Y2.oldTime();
        Y2.oldTime().primitiveFieldRef() =
            Yi.oldTime().primitiveField()/Dbar.oldTime().primitiveField();

        surfaceScalarField gamma
        (
            "Henry_gamma_" + Yi.name(),
            fvc::interpolate(D1*alpha1 + H*D2*(scalar(1) - alpha1))
        );

        surfaceScalarField phiDbar
        (
            "phiDbar_" + Yi.name(),
            H*phi + (scalar(1) - H)*alphaPhi10
        );

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix Y2Eqn
            (
                fvm::ddt(Dbar, Y2)
              + fvm::div(phiDbar, Y2, "div(phi,Yi)")
              - fvm::laplacian(gamma, Y2)
            );

            if (simple.finalNonOrthogonalIter())
            {
                Y2Eqn.solve(mesh.solver("YiFinal"));
            }
            else
            {
                Y2Eqn.solve(mesh.solver("Yi"));
            }
        }

        Yi.primitiveFieldRef() = Dbar.primitiveField()*Y2.primitiveField();
        Yi.correctBoundaryConditions();

        Info<< "Species concentration = "
            << Yi.weightedAverage(mesh.V()).value()
            << "  Min(Yi) = " << gMin(Yi.internalField())
            << "  Max(Yi) = " << gMax(Yi.internalField())
            << endl;
    }

    // Conditional BRNS reaction or transport-only
    if (brnsEnabled)
    {
        #include "BRNSReaction.H"

        forAll(brnsSpecies, k)
        {
            if (brnsIsSolidPhase[k])
            {
                volScalarField& X = *brnsFields[k];
                Info<< "Solid phase " << brnsSpecies[k] << " avg = "
                    << X.weightedAverage(mesh.V()).value()
                    << "  Min = " << gMin(X.internalField())
                    << "  Max = " << gMax(X.internalField()) << nl;
            }
        }
    }
    else
    {
        // Transport-only mode: just correct the mixture
        speciesMixture.correct();
    }
}
