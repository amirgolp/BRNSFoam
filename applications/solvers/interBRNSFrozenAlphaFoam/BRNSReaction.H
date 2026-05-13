{
    double porosity = 1.0;

    // temporal vectors reused per cell
    static std::vector<double> conc;       conc.resize(nBrns, 0.0);
    static std::vector<int>    bflag;      bflag.resize(nBrns, 0);
    int numComp = static_cast<int>(nBrns);
    int retVal  = -1;

    label nReacted = 0;

    forAll(Surf.boundaryField(), patchi)
    {
        const fvPatchScalarField& sf = Surf.boundaryField()[patchi];
        if (sf.type() != "reactingWall") continue;

        const labelList& faceCells = sf.patch().faceCells();

        forAll(faceCells, facei)
        {
            const label c = faceCells[facei];

            // Fill concentrations in the exact brnsSpecies order
            forAll(brnsSpecies, k)
            {
                volScalarField& fld = *brnsFields[k];
                const scalar val = fld[c];
                conc[k] = (val > 1e-20) ? double(val) : 0.0;
                bflag[k] = 0; // default: no boundary constraint
            }

            // Optional: only react if aqueous phase dominates this cell
            double sat = std::max(0.0, std::min(1.0, double(alpha1[c])));
            if (sat > 0.5)
            {
                // Cell centre in µm (BRNS coordinate system)
                const vector& cc = mesh.C()[c];
                double pos_x = double(cc.x()) * 1e6;
                double pos_y = double(cc.y()) * 1e6;
                double pos_z = double(cc.z()) * 1e6;

                double dt = double(runTime.deltaT().value());
                invokebrns_(
                    conc.data(),     // current
                    conc.data(),     // previous (same array is ok)
                    conc.data(),     // output  (in-place)
                    &numComp,
                    &dt,
                    bflag.data(),
                    &retVal,
                    &pos_x, &pos_y, &pos_z,
                    &porosity,
                    &sat,
                    /* parameterVector */ nullptr
                );

                if (retVal != 0)
                {
                    WarningInFunction << "invokebrns_ returned code " << retVal
                                      << " at cell " << c << " (patch " << patchi << ")\n";
                }

                // Write back to the correct fields
                forAll(brnsSpecies, k)
                {
                    volScalarField& fld = *brnsFields[k];
                    const scalar newVal = scalar( conc[k] > 1e-20 ? conc[k] : 0.0 );
                    fld[c] = newVal;

                    // keep patch face equilibrated with owner cell for solid phase
                    if (brnsIsSolidPhase[k])
                    {
                        fld.boundaryFieldRef()[patchi][facei] = newVal;
                    }
                }

                ++nReacted;
            }
        }
    }

    if (nReacted)
    {
        Info<< "BRNS: reacted " << nReacted << " wall-adjacent cells." << nl << endl;

        // Update mixture only if your reactive mixture needs it; safe no-op otherwise
        speciesMixture.correct();

        // Make sure boundaries are consistent
        forAll(brnsSpecies, k)
        {
            brnsFields[k]->correctBoundaryConditions();
        }
    }
}
