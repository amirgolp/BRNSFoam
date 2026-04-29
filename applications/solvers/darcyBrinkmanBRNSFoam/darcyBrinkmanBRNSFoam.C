/*---------------------------------------------------------------------------*\

License
    This file is part of BRNSFoam, derived from GeoChemFoam.

    BRNSFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. See <http://www.gnu.org/licenses/>.

Application
    darcyBrinkmanBRNSFoam

Description
    Single-phase Darcy-Brinkman-Stokes solver coupled to BRNS biomass /
    biogeochemical kinetics.

    Flow is governed by the volume-averaged DBS momentum equation with a
    Kozeny-Carman drag closure on a porosity field eps. Multi-species
    transport uses a porosity-weighted advection-diffusion equation; BRNS
    is invoked at reactingWall-tagged patches each step (auto-detected
    from the brnsSpecies list, transport-only otherwise).

    Intended for continuum / hybrid pore-Darcy cases where the geometry
    is not fully resolved. For geometry-resolved pore-scale work use
    interBRNSFoam (two-phase) or interBRNSALEFoam (moving mesh).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define ifMonitor if (runTime.timeIndex() % 10 == 0)

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
        "Single-phase Darcy-Brinkman-Stokes solver with optional BRNS\n"
        "biogeochemical reactions on reactingWall patches."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"

    #include "createFields.H"
    #include "createPorousMediaFields.H"

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
                        X[c] = solidPhaseInit;
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

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Refresh porosity-derived drag (may have changed last step via reactions)
        #include "updateVariables.H"

        // --- Pressure-velocity PIMPLE corrector loop (use 1 outer for PISO mode)
        while (pimple.loop())
        {
            #include "UEqn.H"

            while (pimple.correct())
            {
                #include "pEqn.H"
            }
        }

        #include "YiEqn.H"

        runTime.write();

        ifMonitor
        {
            Info<< "\n         Umax = " << max(mag(U)).value() << " m/s  "
                << "Uavg = " << mag(average(U)).value() << " m/s";
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
