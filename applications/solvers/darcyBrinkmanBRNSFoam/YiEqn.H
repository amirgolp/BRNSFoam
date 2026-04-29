{
    // Single-phase species transport with porosity weighting:
    //   d(eps*Y)/dt + div(phi*Y) - laplacian(eps*D, Y) = 0
    //
    // No Henry mapping (no second phase). For each species, solve and clip
    // tiny negatives. After the species loop, invoke BRNS at reactingWall
    // patches if enabled.

    forAll(solutionSpecies, i)
    {
        volScalarField& Yi = Y[i];
        const dimensionedScalar& Di = DY[i];

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix YiEqn
            (
                fvm::ddt(eps, Yi)
              + fvm::div(phi, Yi)
              - fvm::laplacian(eps*Di, Yi)
              ==
                fvOptions(Yi)
            );

            YiEqn.relax();
            fvOptions.constrain(YiEqn);

            if (pimple.finalNonOrthogonalIter())
            {
                YiEqn.solve(mesh.solver("YiFinal"));
            }
            else
            {
                YiEqn.solve(mesh.solver("Yi"));
            }
        }

        Yi.max(0.0);
        Yi.correctBoundaryConditions();

        Info<< "Species " << solutionSpecies[i]
            << " avg = " << Yi.weightedAverage(mesh.V()).value()
            << "  Min = " << gMin(Yi.internalField())
            << "  Max = " << gMax(Yi.internalField()) << endl;
    }

    if (brnsEnabled)
    {
        #include "BRNSReaction.H"

        forAll(brnsSpecies, k)
        {
            if (brnsIsSolidPhase[k])
            {
                volScalarField& X = *brnsFields[k];
                Info<< "Solid phase " << brnsSpecies[k]
                    << " avg = " << X.weightedAverage(mesh.V()).value()
                    << "  Min = " << gMin(X.internalField())
                    << "  Max = " << gMax(X.internalField()) << nl;
            }
        }
    }
}
