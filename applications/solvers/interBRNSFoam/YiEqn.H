{
    // Henry-mapped formulation. With Dbar = alpha + H*(1-alpha), the
    // combined diffusion + Henry-transfer face flux equals -gamma*snGrad(Y2)
    // for Y2 = Yi/Dbar and gamma = D1*alpha + H*D2*(1-alpha). Solving for Y2
    // makes the transfer term implicit and unconditionally stable; Yi is
    // recovered as Dbar*Y2 after the solve.

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

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix Y2Eqn
            (
                fvm::ddt(Dbar, Y2)
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
