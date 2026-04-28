#!/bin/bash
#------------------------------------------------------------------------------
# BRNSFoam Build Test Script
#
# Run this script on your HPC to test the complete build process
# Usage: bash test_build.sh
#------------------------------------------------------------------------------

set -e  # Exit on any error

echo "============================================"
echo "BRNSFoam Build Test"
echo "============================================"
echo ""

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print status
print_status() {
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ $1${NC}"
    else
        echo -e "${RED}✗ $1${NC}"
        exit 1
    fi
}

# Function to print info
print_info() {
    echo -e "${YELLOW}→ $1${NC}"
}

#------------------------------------------------------------------------------
# Step 1: Check Prerequisites
#------------------------------------------------------------------------------
echo "Step 1: Checking Prerequisites"
echo "----------------------------------------------"

print_info "Checking OpenFOAM environment..."
if [ -z "$WM_PROJECT_VERSION" ]; then
    echo -e "${RED}ERROR: OpenFOAM environment not loaded!${NC}"
    echo "Please run: source /path/to/OpenFOAM/etc/bashrc"
    exit 1
fi
echo "  OpenFOAM version: $WM_PROJECT_VERSION"
print_status "OpenFOAM environment OK"

print_info "Checking gfortran compiler..."
if ! command -v gfortran &> /dev/null; then
    echo -e "${RED}ERROR: gfortran not found!${NC}"
    echo "Please install gfortran or load module: module load gcc"
    exit 1
fi
gfortran --version | head -n 1
print_status "gfortran found"

print_info "Checking wmake..."
if ! command -v wmake &> /dev/null; then
    echo -e "${RED}ERROR: wmake not found!${NC}"
    echo "Make sure OpenFOAM is properly sourced"
    exit 1
fi
print_status "wmake found"

echo ""

#------------------------------------------------------------------------------
# Step 2: Source BRNSFoam Environment
#------------------------------------------------------------------------------
echo "Step 2: Loading BRNSFoam Environment"
echo "----------------------------------------------"

print_info "Sourcing etc/bashrc..."
source etc/bashrc
print_status "Environment loaded"

echo "  BRNSFOAM_DIR:    $BRNSFOAM_DIR"
echo "  BRNSFOAM_APPBIN: $BRNSFOAM_APPBIN"
echo "  BRNSFOAM_LIBBIN: $BRNSFOAM_LIBBIN"
echo "  BRNS_DIR:        $BRNS_DIR"
echo "  BRNS_LIBBIN:     $BRNS_LIBBIN"

echo ""

#------------------------------------------------------------------------------
# Step 3: Clean Previous Build
#------------------------------------------------------------------------------
echo "Step 3: Cleaning Previous Build"
echo "----------------------------------------------"

print_info "Running Allwclean..."
./Allwclean > /dev/null 2>&1 || true
print_status "Cleaned"

echo ""

#------------------------------------------------------------------------------
# Step 4: Build BRNS Library
#------------------------------------------------------------------------------
echo "Step 4: Building BRNS Library"
echo "----------------------------------------------"

print_info "Building lib/BRNS/lib/brns.so..."
cd lib/BRNS/src

make clean > /dev/null 2>&1 || true
make 2>&1 | tee ../../../brns_build.log | grep -E "gfortran|brns.so|error|Error|fatal"

if [ ! -f "../lib/brns.so" ]; then
    echo -e "${RED}ERROR: brns.so not created!${NC}"
    echo "Check brns_build.log for details"
    exit 1
fi

cd ../../..
print_status "BRNS library built successfully"

ls -lh lib/BRNS/lib/brns.so
echo ""

#------------------------------------------------------------------------------
# Step 5: Build Source Libraries
#------------------------------------------------------------------------------
echo "Step 5: Building Source Libraries"
echo "----------------------------------------------"

print_info "Building src/ libraries..."
cd src
./Allwmake 2>&1 | tee ../src_build.log | grep -E "wmake|error|Error|fatal|up to date"
cd ..
print_status "Source libraries built"

# Check key libraries
print_info "Verifying libraries..."
LIBS=(
    "lib/libreactionThermophysicalModelsBRNSFOAM.so"
    "lib/libinterfacePropertiesBRNSFOAM.so"
    "lib/libimmiscibleIncompressibleTwoPhaseMixtureBRNSFOAM.so"
)

for lib in "${LIBS[@]}"; do
    if [ -f "$lib" ]; then
        echo -e "  ${GREEN}✓${NC} $(basename $lib)"
    else
        echo -e "  ${RED}✗${NC} $(basename $lib) - MISSING"
    fi
done

echo ""

#------------------------------------------------------------------------------
# Step 6: Build Applications
#------------------------------------------------------------------------------
echo "Step 6: Building Applications"
echo "----------------------------------------------"

print_info "Building applications/solvers..."
cd applications
./Allwmake 2>&1 | tee ../apps_build.log | grep -E "wmake|error|Error|fatal|up to date"
cd ..
print_status "Applications built"

# Check key solvers
print_info "Verifying solvers..."
SOLVERS=(
    "bin/interBRNSFoam"
    "bin/interBRNSALEFoam"
    "bin/laplacianFoam"
    "bin/scalarTransportDBSFoam"
)

for solver in "${SOLVERS[@]}"; do
    if [ -f "$solver" ]; then
        echo -e "  ${GREEN}✓${NC} $(basename $solver)"
    else
        echo -e "  ${RED}✗${NC} $(basename $solver) - MISSING"
    fi
done

echo ""

#------------------------------------------------------------------------------
# Step 7: Summary
#------------------------------------------------------------------------------
echo "============================================"
echo "Build Test Summary"
echo "============================================"

# Count files
num_libs=$(find lib -maxdepth 1 -name "*.so" -type f | wc -l)
num_solvers=$(find bin -type f -executable | wc -l)

echo "Libraries built:     $num_libs"
echo "Executables built:   $num_solvers"
echo ""

# Check critical components
all_ok=true

if [ ! -f "lib/BRNS/lib/brns.so" ]; then
    echo -e "${RED}✗ BRNS library missing${NC}"
    all_ok=false
else
    echo -e "${GREEN}✓ BRNS library${NC}"
fi

if [ ! -f "bin/interBRNSFoam" ]; then
    echo -e "${RED}✗ interBRNSFoam solver missing${NC}"
    all_ok=false
else
    echo -e "${GREEN}✓ interBRNSFoam solver${NC}"
fi

if [ ! -f "lib/libreactionThermophysicalModelsBRNSFOAM.so" ]; then
    echo -e "${RED}✗ Reaction thermo library missing${NC}"
    all_ok=false
else
    echo -e "${GREEN}✓ Reaction thermo library${NC}"
fi

echo ""

if [ "$all_ok" = true ]; then
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}   BUILD SUCCESSFUL!${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo ""
    echo "Build logs saved:"
    echo "  - brns_build.log (BRNS library)"
    echo "  - src_build.log (source libraries)"
    echo "  - apps_build.log (applications)"
    echo ""
    echo "Next steps:"
    echo "  1. Try running a tutorial case in tutorials/"
    echo "  2. Check 'which interBRNSFoam' to verify PATH"
    echo ""
else
    echo -e "${RED}========================================${NC}"
    echo -e "${RED}   BUILD FAILED!${NC}"
    echo -e "${RED}========================================${NC}"
    echo ""
    echo "Check build logs for errors:"
    echo "  - brns_build.log"
    echo "  - src_build.log"
    echo "  - apps_build.log"
    echo ""
    exit 1
fi

#------------------------------------------------------------------------------
