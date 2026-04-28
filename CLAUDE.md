# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GeoChemFoam (v5.1) is an OpenFOAM-based package for pore-scale reactive multiphase transport simulation. It extends OpenFOAM v2212 with specialized solvers for chemical reactions, species transport, and multiphase flow at the micro-scale.

## Repository Topology

This local working directory is **not** the source of truth. Three locations exist:

1. **HPC upstream (authoritative for code + builds):**
   `nhdfamgo@login.cluster.uni-hannover.de:~/BRNSFoam-new`
   This is where the code is built and run. OpenFOAM v2212 and the toolchain only live on the cluster — `wmake` cannot be invoked locally.
   When the user has it SFTP-mounted on the laptop, it appears at:
   `/run/user/1000/gvfs/sftp:host=login.cluster.uni-hannover.de,user=nhdfamgo/home/nhdfamgo/BRNSFoam-new`
   The mount is only present when the user has opened the SFTP connection in Files; if it's missing, ask before assuming it's available.

2. **This local working copy** (`/home/amir/workspace/BRNSFoam-new`):
   Used for code editing only. Builds happen on the HPC. Sync direction is upstream → local for refresh, and local → upstream for changes the user wants deployed. Use `rsync -c --exclude='.git' --exclude='.claude' --exclude='*.o' --exclude='*.dep' --exclude='lnInclude' --exclude='linux64*' --exclude='processor*' --exclude='*.so'` and confirm direction with the user before any sync that overwrites.

3. **GitHub mirror:** `git@github.com:amirgolp/BRNSFoam.git` (origin/main).
   Tracks the local working copy. Not authoritative — the HPC is. Push here only when the user asks.

## Build Commands

```bash
# Set up environment (required before any build/run commands)
source etc/bashrc

# Build everything (ThirdParty libs, src libs, solvers, utilities)
./Allwmake > logMake.out

# Clean all build artifacts
./Allwclean

# Verify installation
./checkInstall.sh
```

### Building Individual Components

```bash
# Build specific library
cd src/transportModels/interfaceProperties && wmake libso

# Build specific solver
cd applications/solvers/multiphase/interFoam && wmake

# Build specific utility
cd applications/utilities/postProcessing/processConcentration && wmake
```

## Build Hierarchy

Build order matters due to library dependencies:
1. `lib/BRNS` - BRNS biomass reaction network library (Fortran)
2. `src/thermophysicalModels/` - Reaction thermo models
3. `src/` - Core libraries (finiteVolume, fvMotionSolver, meshTools)
4. `src/transportModels/` - Two-phase mixture, interface properties (6 libraries)
5. `applications/solvers/` - 18 solver applications
6. `applications/utilities/` - Pre/post-processing tools

### Building BRNS Library

The BRNS (Biogeochemical Reaction Network Simulator) library provides biomass reaction kinetics for biofilm growth simulations. It must be built before compiling BRNS-based solvers (e.g., `interBRNSFoam`, `interBRNSALEFoam`).

**Location:** `lib/BRNS/`

**Prerequisites:**
- gfortran compiler
- BRNSFoam environment loaded (`source etc/bashrc`)

**Build methods:**

```bash
# Method 1: Using Allwmake (builds BRNS + all components)
./Allwmake

# Method 2: Build BRNS library only
cd lib/BRNS && ./Allwmake

# Method 3: Build directly with make
cd lib/BRNS/src && make
```

**Clean BRNS build:**

```bash
cd lib/BRNS && ./Allwclean
# or
cd lib/BRNS/src && make clean
```

**How it works:**
- Compiles 16 Fortran source files (`.f` and `.F`) with `-fPIC` flag
- Links object files into shared library `brns.so`
- Places library in `lib/BRNS/lib/brns.so`
- Automatically added to `LD_LIBRARY_PATH` by environment setup

**Required source files:**
- Core: `basic.f`, `biogeo.f`, `boundaries.f`, `drivervalues.f`, `gaussj.f`, `jacobian.f`, `limits.f`, `rates.f`, `residual.f`, `switches.f`, `newtonsub.f`
- LU decomposition: `LUBKSB.F`, `LUDCMP.F`, `MPROVE.F`
- Interface: `BrnsDll/invokebrns.f`, `BrnsDll/parameters.f`
- Includes: 6 `.inc` files for common blocks and definitions

**Troubleshooting:**
- Undefined BRNS symbols during solver build: Rebuild `brns.so` with `make clean && make`
- Missing library: Verify `lib/BRNS/lib/brns.so` exists before building solvers
- Fortran/C interface issues: Ensure gfortran is compatible with your OpenFOAM build

## Code Architecture

### Solver Hierarchy

Solvers build upon each other with increasing complexity:

```
Basic Solvers
├── laplacianFoam (scalar diffusion)
├── dispersionFoam (dispersion closure)
└── scalarTransportDBSFoam (DBS scalar transport)

Transport Layer
└── multiSpeciesTransportFoam (multi-species diffusion)

Multiphase Layer (VOF-based)
├── interFoam → interOSFoam, diffInterFoam
└── interTransportFoam, interTransferFoam

Reactive Layer
├── reactiveTransportFoam (single-phase + Phreeqc)
├── reactiveTransportDBSFoam, reactiveTransportALEFoam
├── interReactiveTransportFoam (two-phase + Phreeqc)
└── interBRNSFoam (two-phase + BRNS biomass reactions)
```

### Key Libraries

| Library | Purpose |
|---------|---------|
| `reactionThermophysicalModelsGCFOAM` | Phreeqc-based chemical equilibrium |
| `interfacePropertiesGCFOAM` | Surface tension, contact angles |
| `immiscibleIncompressibleTwoPhaseMixtureGCFOAM` | Two-phase flow properties |
| `twoPhaseProperties` | Wetting and contact angle models |

### Solver File Structure

Each solver follows this pattern:
```
solverName/
├── solverName.C          # Main solver loop
├── createFields.H        # Field initialization
├── UEqn.H, pEqn.H        # Momentum, pressure equations
├── alphaEqn.H            # VOF phase fraction (multiphase)
├── YiEqn.H               # Species transport
└── Make/
    ├── files             # Source files list
    └── options           # Compiler flags, library links
```

### Custom Boundary Conditions

Located in `src/thermophysicalModels/reactionThermo/`:
- `reactingWall` - Surface reactions at walls
- `globalConcentrationMixed` - Concentration-based BC
- `reactiveSurfaceConcentrationMixed` - Reactive surface BC

## Running Tutorial Cases

Tutorial cases in `tutorials/` follow a standard workflow:

```bash
cd tutorials/basic/laplacianFoam/Fontainebleau

# 1. Create mesh from geometry
./createMesh.sh

# 2. Initialize fields
./initCase.sh

# 3. Run simulation
./runCase.sh

# 4. Post-process results
./processCase.sh

# Cleanup
./deleteAll.sh  # Remove all generated files
./deleteT.sh    # Remove time directories only
```

Parallel runs are automatic when `processor*` directories exist from decomposition.

## Environment Variables

After `source etc/bashrc`:
- `GCFOAM_DIR` - Installation root
- `GCFOAM_SRC` - Source libraries
- `GCFOAM_APPBIN` - Compiled executables (added to PATH)
- `GCFOAM_LIBBIN` - Compiled libraries (added to LD_LIBRARY_PATH)
- `GCFOAM_TUTORIALS` - Tutorial cases

Aliases: `dir`, `src`, `app`, `tut`, `run` for quick navigation.

## Third-Party Dependencies

- **OpenFOAM v2212** - Base CFD framework
- **BRNS** - Biomass reaction network (Fortran library in `lib/BRNS/`, called via extern C interface)
- **gfortran** - GNU Fortran compiler (required for BRNS library)

## Key File Locations

- Solver sources: `applications/solvers/`
- Post-processing utilities: `applications/utilities/postProcessing/`
- Transport models: `src/transportModels/`
- Reaction chemistry: `src/thermophysicalModels/reactionThermo/`
- Example cases: `tutorials/`
