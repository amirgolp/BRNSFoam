/*---------------------------------------------------------------------------*\

License
    This file is part of BRNSFoam, derived from GeoChemFoam.

    BRNSFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. See <http://www.gnu.org/licenses/>.

Application
    interBRNSFoam

Description
    Solver for two incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes.

    Features:
    - Multi-species transport using CST method
    - Optional BRNS biomass reactions (auto-detected from brnsSpecies list)
    - Transport-only mode when brnsSpecies is empty or missing

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "simpleControl.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "inertMultiComponentMixture.H"
#include "basicTwoPhaseMultiComponentTransportMixture.H"
#include "twoPhaseMultiComponentTransportMixture.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Define monitor time indexes
#define ifMonitor if (runTime.timeIndex() % 10 == 0)

// BRNS Fortran interface
extern "C" {
void invokebrns_(
    double *theCurArray, double *thePreArray, double *outputArray, int *numComp,
    double *time_step, int *boundary_flag, int *return_value,
    double *x_pos, double *y_pos, double *z_pos,
    double *porosity, double *saturation, double *parameterVector
);
}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids\n"
        "using VOF phase-fraction based interface capturing.\n"
        "Supports optional BRNS biomass reactions (auto-detected).\n"
        "Transport-only mode when brnsSpecies is empty or missing."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    basicTwoPhaseMultiComponentTransportMixture& speciesMixture = pSpeciesMixture();

    // BRNS >>> initialize biomass at reactingWall (only if BRNS enabled)
    if (brnsEnabled)
    {
        Info<< "Initializing solid-phase fields on reactingWall patches..." << nl << endl;
        forAll(Surf.boundaryField(), patchi)
        {
            const fvPatchScalarField& sf = Surf.boundaryField()[patchi];
            if (sf.type() == "reactingWall")
            {
                const labelList& faceCells = sf.patch().faceCells();
                forAll(brnsSpecies, k)
                {
                    if (!brnsIsSolidPhase[k]) continue;
                    volScalarField& X = *brnsFields[k];
                    forAll(faceCells, facei)
                    {
                        const label c = faceCells[facei];
                        // Initialize adjacent cell concentration
                        X[c] = solidPhaseInit;
                        // Keep patch value equal to owner cell (equilibration)
                        X.boundaryFieldRef()[patchi][facei] = X[c];
                    }
                    X.correctBoundaryConditions();
                    X.write();
                }
            }
        }
        Info<< "Solid-phase initialization done." << nl << endl;
    }
    // BRNS <<<

    #include "createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "alphaControls.H"

            #include "YiMulesEqn.H"
            #include "alphaEqnSubCycle.H"

            gradalpha1 = mag(fvc::grad(alpha1));

            mixture.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        #include "YiEqn.H"

        runTime.write();

        // Monitor average and max velocity
        ifMonitor
        {
            Info << "\n         Umax = " << max(mag(U)).value() << " m/s  "
                 << "Uavg = " << mag(average(U)).value() << " m/s";
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
