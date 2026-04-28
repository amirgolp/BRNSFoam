# Build Test Sequence for HPC

Run these commands on your HPC to test the build process.

## Prerequisites Check

```bash
# 1. Verify OpenFOAM is loaded
echo $WM_PROJECT_VERSION
# Should output: 2212 (or your version)

# 2. Check gfortran
which gfortran
gfortran --version
# Should show GNU Fortran version

# 3. Check wmake
which wmake
# Should show OpenFOAM wmake path
```

If any check fails:
```bash
# Load OpenFOAM
source /path/to/OpenFOAM/etc/bashrc

# Load gfortran (if needed on HPC)
module load gcc
# or
module load gfortran
```

---

## Automated Test (Recommended)

```bash
cd /path/to/BRNSFoam-new

# Run automated test script
./test_build.sh
```

This will:
- ✓ Check all prerequisites
- ✓ Clean previous build
- ✓ Build BRNS library
- ✓ Build source libraries
- ✓ Build all solvers
- ✓ Verify critical components
- ✓ Generate build logs

**Expected output:** "BUILD SUCCESSFUL!" with green checkmarks

**If build fails:** Check the generated log files:
- `brns_build.log` - BRNS library errors
- `src_build.log` - Source library errors
- `apps_build.log` - Solver/application errors

---

## Manual Test (Step-by-Step)

If you prefer to run each step manually:

### Step 1: Load Environment

```bash
cd /path/to/BRNSFoam-new

# Load OpenFOAM first
source /path/to/OpenFOAM/etc/bashrc

# Load BRNSFoam
source etc/bashrc
```

**Verify:**
```bash
echo $BRNSFOAM_DIR
echo $BRNS_LIBBIN
# Should show correct paths
```

### Step 2: Clean Previous Build

```bash
./Allwclean
```

### Step 3: Build BRNS Library

```bash
cd lib/BRNS/src

make clean
make 2>&1 | tee brns_build.log

# Check output
ls -lh ../lib/brns.so
# Should show: brns.so with size ~50-100KB

cd ../../..
```

**Expected output:**
```
gfortran -O4 -x f77-cpp-input -I. -fPIC -c basic.f -o basic.o
gfortran -O4 -x f77-cpp-input -I. -fPIC -c biogeo.f -o biogeo.o
...
gfortran -L. -O4 -shared -s basic.o biogeo.o ... -o brns.so
```

**If this fails:** Check compiler errors in output, ensure gfortran works

### Step 4: Build Source Libraries

```bash
cd src
./Allwmake 2>&1 | tee src_build.log
cd ..
```

**Check libraries were created:**
```bash
ls -lh lib/*.so
```

**Should see:**
- `libreactionThermophysicalModelsBRNSFOAM.so`
- `libinterfacePropertiesBRNSFOAM.so`
- `libimmiscibleIncompressibleTwoPhaseMixtureBRNSFOAM.so`
- `libtwoPhaseProperties.so`
- Plus other OpenFOAM libraries

### Step 5: Build Applications

```bash
cd applications
./Allwmake 2>&1 | tee apps_build.log
cd ..
```

**Check solvers were created:**
```bash
ls -lh bin/
```

**Should see executables like:**
- `interBRNSFoam`
- `interBRNSALEFoam`
- `laplacianFoam`
- `scalarTransportDBSFoam`
- Many others (18+ solvers)

### Step 6: Verify Installation

```bash
# Check BRNS library
ls -lh lib/BRNS/lib/brns.so

# Check critical solver
which interBRNSFoam
interBRNSFoam -help

# Count built components
echo "Libraries: $(find lib -maxdepth 1 -name '*.so' | wc -l)"
echo "Executables: $(find bin -type f | wc -l)"
```

**Expected:**
```
Libraries: 10-15
Executables: 18-25
```

---

## Quick Rebuild Commands

After successful initial build, for incremental changes:

```bash
# Rebuild just BRNS
cd lib/BRNS/src && make clean && make && cd ../../..

# Rebuild just one solver
cd applications/solvers/interBRNSFoam && wmake && cd ../../..

# Rebuild all source libraries
cd src && ./Allwmake && cd ..

# Full rebuild
./Allwclean && ./Allwmake 2>&1 | tee build.log
```

---

## Common Errors and Solutions

### Error: "BRNS library not found"

```bash
# Check if brns.so exists
ls -lh lib/BRNS/lib/brns.so

# If missing, rebuild
cd lib/BRNS/src
make clean
make
cd ../../..
```

### Error: "undefined reference to BRNS functions"

The BRNS library wasn't built before the solver. Rebuild in order:

```bash
./Allwclean
cd lib/BRNS/src && make && cd ../../..
cd src && ./Allwmake && cd ..
cd applications && ./Allwmake && cd ..
```

### Error: "OpenFOAM headers not found"

OpenFOAM environment not loaded:

```bash
source /path/to/OpenFOAM/etc/bashrc
source etc/bashrc
```

### Error: "gfortran: command not found"

Load fortran compiler:

```bash
# On HPC with modules
module load gcc

# Or check available compilers
module avail gcc
module avail gfortran
```

---

## Success Criteria

Build is successful when:

- ✓ No compilation errors in any log
- ✓ `lib/BRNS/lib/brns.so` exists (~50-100KB)
- ✓ At least 10 libraries in `lib/*.so`
- ✓ At least 15 executables in `bin/`
- ✓ `which interBRNSFoam` shows correct path
- ✓ `interBRNSFoam -help` runs without errors

---

## Report Back

After running the test, report:

1. **Test script output:**
   - Did `./test_build.sh` complete successfully?
   - Copy the final summary section

2. **Component counts:**
   ```bash
   find lib -maxdepth 1 -name "*.so" | wc -l
   find bin -type f | wc -l
   ```

3. **Critical files:**
   ```bash
   ls -lh lib/BRNS/lib/brns.so
   ls -lh bin/interBRNSFoam
   ```

4. **Any errors:**
   - Check `brns_build.log`, `src_build.log`, `apps_build.log`
   - Copy relevant error messages
