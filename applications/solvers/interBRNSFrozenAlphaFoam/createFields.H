#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
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


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

volScalarField gradalpha1
(
    IOobject
    (
        "gradalpha1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(fvc::grad(alpha1))
);

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

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

bool phaseTransfer = thermoPhysicalProperties.lookupOrDefault("PhaseTransfer",false);
if (phaseTransfer==true)
{
    FatalErrorInFunction
        << "Phase transfer not supported in interBRNSFoam.\n"
        << abort(FatalError);
}


Info << "Create species mixture\n" << endl;
autoPtr<basicTwoPhaseMultiComponentTransportMixture> pSpeciesMixture
(
    basicTwoPhaseMultiComponentTransportMixture::New(mesh, alpha1)
);

// ============================================================================
// BRNS Auto-Detection
// ============================================================================
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
        Info<< "BRNS mode ENABLED: found " << nBrns << " species in brnsSpecies list" << nl << endl;
    }
    else
    {
        Info<< "BRNS mode DISABLED: brnsSpecies list is empty" << nl << endl;
    }
}
else
{
    Info<< "BRNS mode DISABLED: no brnsSpecies entry found (transport-only mode)" << nl << endl;
}

// BRNS >>> read Surf (to detect reactingWall BC)
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

// Optional: initial value for solid-phase wall-adjacent cells
scalar solidPhaseInit = thermoPhysicalProperties.lookupOrDefault<scalar>("solidPhaseInit", 0.25);
(void)solidPhaseInit; // Used later in brnsEnabled block

// Configurable solid-phase field detection (replaces hardcoded X_ prefix)
wordHashSet solidPhaseSet;
if (thermoPhysicalProperties.found("solidPhaseFields"))
{
    solidPhaseSet = wordHashSet(wordList(thermoPhysicalProperties.lookup("solidPhaseFields")));
    Info<< "Solid-phase fields: " << solidPhaseSet << nl << endl;
}

// Conditional BRNS field setup
List<volScalarField*> brnsFields;
List<bool> brnsIsSolidPhase;

if (brnsEnabled)
{
    // Map brnsSpecies -> volScalarField*
    const speciesTable& solutionSpecies = pSpeciesMixture().species();
    HashTable<label> soluteIdx;
    forAll(solutionSpecies, si) { soluteIdx.insert(solutionSpecies[si], si); }

    brnsFields.resize(nBrns, nullptr);
    brnsIsSolidPhase.resize(nBrns, false);

    forAll(brnsSpecies, i)
    {
        const word& nm = brnsSpecies[i];
        const bool isSolid = solidPhaseSet.found(nm);
        brnsIsSolidPhase[i] = isSolid;

        if (isSolid)
        {
            // Solid-phase field exists as a regular volume field, not in speciesMixture
            Info<< "Reading solid-phase field '" << nm << "'\n" << endl;
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
                    << "BRNS solute '" << nm << "' not found in mixture species.\n"
                    << "Available: " << solutionSpecies << nl
                    << exit(FatalError);
            }
            brnsFields[i] = &pSpeciesMixture().Y(soluteIdx[nm]);
        }
    }
}
// ============================================================================

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

const volScalarField& pc = mixture.pc();

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + pc + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - pc - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());
mesh.setFluxRequired(pc.name());

// MULES compressed flux is registered in case scalarTransport FO needs it.
surfaceScalarField alphaPhiUn
(
    IOobject
    (
        "alphaPhiUn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(phi.dimensions(), Zero)
);

/******************************************************************/
// scheme for phiH
//******************************************************************//

OStringStream buf;

const tokenList& tokens = mesh.divScheme("div(phiH,Yi)");

label len = tokens.size();

forAll (tokens, tok)
{
    buf << mesh.divScheme("div(phiH,Yi)")[tok];
    if (--len)
    {
        buf << ' ';
    }
}

string phiHScheme = buf.str();

#include "createMRF.H"
#include "createFvOptions.H"
