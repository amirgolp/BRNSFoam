# BRNSFoam

**Pore-scale Biogeochemical Reactive Transport Simulator**

BRNSFoam is an OpenFOAM-based package for simulating biofilm growth and reactive transport at the pore scale. It couples multiphase flow (Volume of Fluid method), multi-component solute transport, and biogeochemical reaction networks (BRNS) for high-resolution simulations in complex porous microstructures.

**Note:** This is a **BRNS-only version**. PhreeqcRM chemical equilibrium solvers have been removed. Supported solvers include `interBRNSFoam`, `interBRNSALEFoam`, and transport-only solvers.

## Prerequisites

### Required Software

1. **OpenFOAM v2212** (or compatible version)
   - Must be sourced before building BRNSFoam
   - Check: `echo $WM_PROJECT_VERSION` should output `2212` or similar

2. **GNU Fortran Compiler (gfortran)**
   - Required for building BRNS library
   - Check: `gfortran --version`
   - Minimum version: 4.8+

3. **Standard build tools**
   - gcc/g++, make, bash
   - Usually available on HPC systems

### Optional Dependencies

- **OpenBLAS** - For BRNS linear algebra operations
- **ParaView** - For visualization (post-processing)

## Installation

### 1. Clone or Download

```bash
cd ~/workspace
# Assuming you already have the BRNSFoam directory
cd BRNSFoam
```

### 2. Source OpenFOAM Environment

```bash
# Load your OpenFOAM installation first
source /path/to/OpenFOAM/etc/bashrc
# Example: source /opt/openfoam2212/etc/bashrc
```

### 3. Source BRNSFoam Environment

```bash
source etc/bashrc
```

This sets up:
- `BRNSFOAM_DIR` - Installation root
- `BRNSFOAM_APPBIN` - Executable binaries (`bin/`)
- `BRNSFOAM_LIBBIN` - Compiled libraries (`lib/`)
- `BRNS_DIR` - BRNS library location (`lib/BRNS`)
- `BRNS_LIBBIN` - BRNS shared library (`lib/BRNS/lib/`)

### 4. Build Everything

```bash
# Build BRNS library + all source libraries + all solvers
./Allwmake 2>&1 | tee build.log
```

Or build components separately:

```bash
# Build BRNS library only
cd lib/BRNS && ./Allwmake && cd ../..

# Build source libraries
cd src && ./Allwmake && cd ..

# Build applications (solvers + utilities)
cd applications && ./Allwmake && cd ..
```

### 5. Verify Installation

Check that key files were created:

```bash
# BRNS library
ls -lh lib/BRNS/lib/brns.so

# Example solver
ls -lh bin/interBRNSFoam

# Source libraries
ls -lh lib/*.so
```

## Build Components

### BRNS Library (`lib/BRNS`)

The BRNS (Biogeochemical Reaction Network Simulator) library handles biomass reaction kinetics:

```bash
cd lib/BRNS/src
make          # Build brns.so
make clean    # Remove object files and library
```

Output: `lib/BRNS/lib/brns.so`

### Source Libraries (`src/`)

Core libraries for reactive transport:
- `reactionThermophysicalModelsBRNSFOAM` - Reaction-coupled thermo models
- `interfacePropertiesBRNSFOAM` - Surface tension, contact angles
- `immiscibleIncompressibleTwoPhaseMixtureBRNSFOAM` - VOF two-phase properties
- `twoPhaseProperties` - Wetting models
- Plus OpenFOAM libraries: `finiteVolume`, `fvMotionSolver`, `meshTools`

### Solvers (`applications/solvers/`)

Key solvers:
- **`interBRNSFoam`** - Two-phase flow with BRNS biomass reactions
- **`interBRNSALEFoam`** - Biofilm growth with mesh motion
- **`laplacianFoam`** - Basic scalar diffusion
- **`scalarTransportDBSFoam`** - DBS transport solver
- Many more (18 total solvers)

### Utilities (`applications/utilities/`)

Pre-processing and post-processing tools for concentration, porosity, permeability, etc.

## Post-Processing: Reaction Rate Fields

After a simulation, BRNS writes reaction rates to `ratesAtFinish.dat` in the case directory. Use `mapRatesToFields.py` to convert these into OpenFOAM `volScalarField` files viewable in ParaView.

### File format

Each line in `ratesAtFinish.dat` corresponds to one reacted cell per BRNS call:

```
x_pos(µm)  y_pos(µm)  z_pos(µm)  r1  r2  ...  rN
```

Coordinates are written in **micrometres** (cell centres from `mesh.C()`, scaled ×1e6 by `BRNSReaction.H`) using full double-precision scientific notation.

### Usage

```bash
# Auto-detect rate names (rate_R1, rate_R2, ...):
python3 mapRatesToFields.py /path/to/case

# Name the rates explicitly:
python3 mapRatesToFields.py . --rate-names "aerobic,anoxic,decay"

# Custom rates file location:
python3 mapRatesToFields.py . --rates-file /path/to/ratesAtFinish.dat

# Override coord scale (default 1e-6 converts µm → m):
python3 mapRatesToFields.py . --coord-scale 1e-6

# Generate an editable config file:
python3 mapRatesToFields.py . --write-config
```

### Output

For each time directory (e.g. `0.05/`):

```
0.05/rate_aerobic    ← volScalarField [mol m⁻³ s⁻¹]
0.05/rate_anoxic
0.05/rate_decay
```

Open the `.foam` file in ParaView — the rate fields appear automatically alongside all other solution fields.

### Config file (`rateConfig.json`)

```json
{
  "rate_names": ["aerobic", "anoxic", "decay"],
  "rates_file": "",
  "tol_coord": 1e-6,
  "coord_scale": 1e-6
}
```

Pass with `--config rateConfig.json`.

---

## Running Simulations

### Tutorial Cases

Example cases are in `tutorials/` directory:

```bash
# Navigate to a tutorial case
cd tutorials/basic/laplacianFoam/Fontainebleau

# Run the case workflow
./createMesh.sh      # Generate mesh
./initCase.sh        # Initialize fields
./runCase.sh         # Run solver
./processCase.sh     # Post-process results

# Clean up
./deleteAll.sh       # Remove all generated files
```

For parallel runs, decompose the domain first, then run with `mpirun`.

## Cleaning Build Artifacts

```bash
# Clean everything (preserves source code)
./Allwclean

# Clean individual components
cd lib/BRNS && ./Allwclean
cd src && ./Allwclean
cd applications && ./Allwclean
```

**Note:** `Allwclean` preserves BRNS source files in `lib/BRNS/src/`

## Troubleshooting

### "BRNS library not found" during solver build

```bash
# Check if brns.so exists
ls -lh lib/BRNS/lib/brns.so

# Rebuild BRNS library
cd lib/BRNS/src
make clean
make
```

### "OpenFOAM not found" errors

```bash
# Ensure OpenFOAM is sourced first
echo $WM_PROJECT_VERSION  # Should show 2212

# Source OpenFOAM, then BRNSFoam
source /path/to/OpenFOAM/etc/bashrc
source etc/bashrc
```

### Fortran compiler errors

```bash
# Check gfortran is available
which gfortran
gfortran --version

# If missing, install or load module on HPC:
module load gcc
module load gfortran
```

### Undefined symbols / linking errors

Make sure build order is correct:
1. BRNS library (`lib/BRNS`)
2. Source libraries (`src/`)
3. Applications (`applications/`)

```bash
# Rebuild in correct order
./Allwclean
./Allwmake 2>&1 | tee build.log
```

## Directory Structure

```
BRNSFoam/
├── README.md                 # This file
├── CLAUDE.md                 # Development guide for Claude Code
├── Allwmake                  # Build script (all components)
├── Allwclean                 # Clean script
├── etc/
│   └── bashrc                # Environment setup
├── lib/
│   └── BRNS/                 # BRNS reaction library
│       ├── src/              # Fortran source files
│       │   ├── Makefile      # Build configuration
│       │   ├── *.f, *.F      # Fortran sources
│       │   └── BrnsDll/      # Interface files
│       └── lib/              # Output: brns.so
├── src/                      # Source libraries
│   ├── finiteVolume/
│   ├── transportModels/
│   ├── thermophysicalModels/
│   └── ...
├── applications/
│   ├── solvers/              # Solver applications
│   └── utilities/            # Pre/post-processing tools
├── tutorials/                # Example cases
├── bin/                      # Compiled executables (generated)
└── lib/                      # Compiled libraries (generated)
```

## Environment Variables

After sourcing `etc/bashrc`:

| Variable | Description | Example |
|----------|-------------|---------|
| `BRNSFOAM_DIR` | Installation root | `/home/user/BRNSFoam` |
| `BRNSFOAM_APPBIN` | Executables directory | `$BRNSFOAM_DIR/bin` |
| `BRNSFOAM_LIBBIN` | Libraries directory | `$BRNSFOAM_DIR/lib` |
| `BRNS_DIR` | BRNS module location | `$BRNSFOAM_DIR/lib/BRNS` |
| `BRNS_LIBBIN` | BRNS library path | `$BRNS_DIR/lib` |

Aliases: `brnsdir`, `brnssrc`, `brnsapp`, `brnstut`

## Citation

If you use BRNSFoam in your research, please cite:

> Golparvar, A., and Thullner, M. (2024). P3D-BRNS v1.0.0: a three-dimensional, multiphase, multicomponent, pore-scale reactive transport modelling package for simulating biogeochemical processes in subsurface environments. *Geoscientific Model Development*, **17**, 881-903.

## Support

For issues and questions:
- Check `CLAUDE.md` for detailed development documentation
- Review tutorial cases in `tutorials/`
- Check build logs in `build.log` if build fails

## License

[Add your license information here]
