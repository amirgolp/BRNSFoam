Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar muf  (transportProperties.lookup("muf"));
dimensionedScalar rhof (transportProperties.lookup("rhof"));

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Reading porosity field eps\n" << endl;
volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Solid
(
    IOobject
    (
        "Solid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Solid", dimless, 0)
);

surfaceScalarField rhofPhi
(
    IOobject
    (
        "rhofPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rhof*phi
);

surfaceScalarField rhoPhiByEps
(
    IOobject
    (
        "rhoPhiByEps",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rhofPhi*linearInterpolate(1.0/(eps + dimensionedScalar("small", dimless, SMALL)))
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

// =========================================================================
// thermoPhysicalProperties + species + BRNS auto-detection
// =========================================================================
Info<< "Reading thermoPhysicalProperties\n" << endl;

IOdictionary thermoPhysicalProperties
(
    IOobject
    (
        "thermoPhysicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Species list and per-species diffusion coefficients live in
// thermoPhysicalProperties. Single-phase form: solutionSpecies + DY entries.
//   solutionSpecies ( O2 NO3 NH4 ... );
//   D_O2  D [ 0 2 -1 0 0 0 0 ] 1.5e-9;
//   ...
wordList solutionSpecies(thermoPhysicalProperties.lookup("solutionSpecies"));
const label nSpecies = solutionSpecies.size();

PtrList<volScalarField> Y(nSpecies);
PtrList<dimensionedScalar> DY(nSpecies);

forAll(solutionSpecies, i)
{
    const word& nm = solutionSpecies[i];
    Info<< "Reading species field " << nm << nl;

    Y.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                nm,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    DY.set
    (
        i,
        new dimensionedScalar
        (
            thermoPhysicalProperties.lookup("D_" + nm)
        )
    );
}

HashTable<label> soluteIdx;
forAll(solutionSpecies, i) { soluteIdx.insert(solutionSpecies[i], i); }

// BRNS detection
bool brnsEnabled = false;
wordList brnsSpecies;
label nBrns = 0;

if (thermoPhysicalProperties.found("brnsSpecies"))
{
    brnsSpecies = wordList(thermoPhysicalProperties.lookup("brnsSpecies"));
    nBrns = brnsSpecies.size();
    if (nBrns > 0)
    {
        brnsEnabled = true;
        Info<< "BRNS mode ENABLED: " << nBrns << " species" << nl << endl;
    }
}
if (!brnsEnabled)
{
    Info<< "BRNS mode DISABLED (transport-only)" << nl << endl;
}

Info<< "Reading field Surf\n" << endl;
volScalarField Surf
(
    IOobject
    (
        "Surf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

scalar solidPhaseInit =
    thermoPhysicalProperties.lookupOrDefault<scalar>("solidPhaseInit", 0.25);
(void)solidPhaseInit;

wordHashSet solidPhaseSet;
if (thermoPhysicalProperties.found("solidPhaseFields"))
{
    solidPhaseSet =
        wordHashSet(wordList(thermoPhysicalProperties.lookup("solidPhaseFields")));
    Info<< "Solid-phase fields: " << solidPhaseSet << nl << endl;
}

List<volScalarField*> brnsFields;
List<bool> brnsIsSolidPhase;

if (brnsEnabled)
{
    brnsFields.resize(nBrns, nullptr);
    brnsIsSolidPhase.resize(nBrns, false);

    forAll(brnsSpecies, i)
    {
        const word& nm = brnsSpecies[i];
        const bool isSolid = solidPhaseSet.found(nm);
        brnsIsSolidPhase[i] = isSolid;

        if (isSolid)
        {
            Info<< "Reading solid-phase field '" << nm << "'\n";
            brnsFields[i] = new volScalarField
            (
                IOobject
                (
                    nm,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            );
        }
        else
        {
            if (!soluteIdx.found(nm))
            {
                FatalErrorInFunction
                    << "BRNS solute '" << nm << "' not found in solutionSpecies.\n"
                    << "Available: " << solutionSpecies << nl
                    << exit(FatalError);
            }
            brnsFields[i] = &Y[soluteIdx[nm]];
        }
    }
}

#include "createFvOptions.H"
