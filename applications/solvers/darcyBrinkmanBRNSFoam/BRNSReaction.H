{
    // Single-phase BRNS coupling: invoke BRNS once per wall-adjacent cell
    // along reactingWall patches. No saturation gating (single phase, sat=1).
    // Uses the local cell porosity (eps), unlike the two-phase variant
    // which fixes porosity=1.

    static std::vector<double> conc;  conc.resize(nBrns, 0.0);
    static std::vector<int>    bflag; bflag.resize(nBrns, 0);
    int numComp = static_cast<int>(nBrns);
    int retVal  = -1;
    double sat  = 1.0;

    label nReacted = 0;

    forAll(Surf.boundaryField(), patchi)
    {
        const fvPatchScalarField& sf = Surf.boundaryField()[patchi];
        if (sf.type() != "reactingWall") continue;

        const labelList& faceCells = sf.patch().faceCells();

        forAll(faceCells, facei)
        {
            const label c = faceCells[facei];

            forAll(brnsSpecies, k)
            {
                volScalarField& fld = *brnsFields[k];
                const scalar val = fld[c];
                conc[k]  = (val > 1e-20) ? double(val) : 0.0;
                bflag[k] = 0;
            }

            const vector& cc = mesh.C()[c];
            double pos_x = double(cc.x()) * 1e6;
            double pos_y = double(cc.y()) * 1e6;
            double pos_z = double(cc.z()) * 1e6;

            double dt = double(runTime.deltaT().value());
            double porosity = double(eps[c]);

            invokebrns_(
                conc.data(),
                conc.data(),
                conc.data(),
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
                WarningInFunction
                    << "invokebrns_ returned code " << retVal
                    << " at cell " << c << " (patch " << patchi << ")\n";
            }

            forAll(brnsSpecies, k)
            {
                volScalarField& fld = *brnsFields[k];
                const scalar newVal = scalar(conc[k] > 1e-20 ? conc[k] : 0.0);
                fld[c] = newVal;

                if (brnsIsSolidPhase[k])
                {
                    fld.boundaryFieldRef()[patchi][facei] = newVal;
                }
            }

            ++nReacted;
        }
    }

    if (nReacted)
    {
        Info<< "BRNS: reacted " << nReacted << " wall-adjacent cells." << nl << endl;

        forAll(brnsSpecies, k)
        {
            brnsFields[k]->correctBoundaryConditions();
        }
    }
}
