#!/bin/bash
#------------------------------------------------------------------------------
# Sync BRNSFoam-new to HPC and Rebuild
#
# INSTRUCTIONS:
# 1. Edit the variables below to match your HPC setup
# 2. Make executable: chmod +x sync_and_rebuild_hpc.sh
# 3. Run: ./sync_and_rebuild_hpc.sh
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# EDIT THESE VARIABLES
#------------------------------------------------------------------------------

# HPC connection details
HPC_USER="nhdfamgo"              # Your HPC username
HPC_HOST="cirrus.ac.uk"           # HPC hostname (e.g., cirrus.ac.uk, archer2.ac.uk)
HPC_PATH="BRNSFoam-new"          # Target directory on HPC (relative to home)

# OpenFOAM path on HPC (needed for rebuild)
OPENFOAM_PATH="/sw/apps/software/arch/MPI/GCC/11.2.0/OpenMPI/4.1.1/OpenFOAM/v2212/OpenFOAM-v2212"

# Optional: Set to "yes" to do dry-run first (recommended)
DRY_RUN="no"

#------------------------------------------------------------------------------
# DO NOT EDIT BELOW THIS LINE
#------------------------------------------------------------------------------

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo "============================================"
echo "BRNSFoam HPC Sync & Rebuild"
echo "============================================"
echo ""
echo "Source:      $(pwd)"
echo "Destination: ${HPC_USER}@${HPC_HOST}:~/${HPC_PATH}/"
echo ""

#------------------------------------------------------------------------------
# Step 1: Sync to HPC
#------------------------------------------------------------------------------

echo -e "${YELLOW}Step 1: Syncing files to HPC...${NC}"
echo "----------------------------------------------"

# Rsync options:
# -a  : archive mode (recursive, preserve permissions)
# -v  : verbose
# -z  : compress during transfer
# -P  : show progress
# --delete : delete files on destination that don't exist in source
# --exclude : exclude build artifacts

RSYNC_OPTS="-avzP --delete"
RSYNC_EXCLUDE=(
    --exclude='*.o'
    --exclude='*.dep'
    --exclude='*.so'
    --exclude='Make/linux64*'
    --exclude='lnInclude/'
    --exclude='bin/'
    --exclude='lib/*.so'
    --exclude='lib/BRNS/lib/brns.so'
    --exclude='*.log'
    --exclude='processor*/'
    --exclude='[0-9]*/'
    --exclude='[0-9]*.*/'
)

if [ "$DRY_RUN" = "yes" ]; then
    echo -e "${YELLOW}DRY RUN MODE - No files will be transferred${NC}"
    RSYNC_OPTS="${RSYNC_OPTS} --dry-run"
fi

# Run rsync
rsync ${RSYNC_OPTS} \
    "${RSYNC_EXCLUDE[@]}" \
    ./ ${HPC_USER}@${HPC_HOST}:~/${HPC_PATH}/

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ Sync completed successfully${NC}"
else
    echo -e "${RED}✗ Sync failed${NC}"
    exit 1
fi

echo ""

if [ "$DRY_RUN" = "yes" ]; then
    echo -e "${YELLOW}DRY RUN COMPLETE - No actual changes made${NC}"
    echo "Set DRY_RUN=\"no\" to perform actual sync"
    exit 0
fi

#------------------------------------------------------------------------------
# Step 2: Rebuild on HPC
#------------------------------------------------------------------------------

echo -e "${YELLOW}Step 2: Rebuilding on HPC...${NC}"
echo "----------------------------------------------"

# Create rebuild script to run on HPC
REBUILD_SCRIPT=$(cat <<'EOF'
#!/bin/bash
set -e

echo "============================================"
echo "Building BRNSFoam on HPC"
echo "============================================"
echo ""

# Source OpenFOAM
echo "Loading OpenFOAM environment..."
source OPENFOAM_PATH_PLACEHOLDER/etc/bashrc
echo "  OpenFOAM version: $WM_PROJECT_VERSION"

# Source BRNSFoam
echo "Loading BRNSFoam environment..."
cd ~/HPC_PATH_PLACEHOLDER
source etc/bashrc

# Clean previous build
echo ""
echo "Cleaning previous build..."
./Allwclean > /dev/null 2>&1 || true

# Build everything
echo ""
echo "Building BRNSFoam..."
./Allwmake 2>&1 | tee build.log

# Check build status
if [ ${PIPESTATUS[0]} -eq 0 ]; then
    echo ""
    echo "============================================"
    echo "Build Verification"
    echo "============================================"

    # Count components
    echo "Libraries built:    $(find lib -maxdepth 1 -name '*.so' 2>/dev/null | wc -l)"
    echo "Executables built:  $(find bin -type f 2>/dev/null | wc -l)"

    # Check critical files
    echo ""
    echo "Critical components:"
    if [ -f "lib/BRNS/lib/brns.so" ]; then
        echo "  ✓ BRNS library"
    else
        echo "  ✗ BRNS library MISSING"
    fi

    if [ -f "bin/interBRNSFoam" ]; then
        echo "  ✓ interBRNSFoam solver"
    else
        echo "  ✗ interBRNSFoam solver MISSING"
    fi

    if [ -f "bin/interBRNSALEFoam" ]; then
        echo "  ✓ interBRNSALEFoam solver"
    else
        echo "  ✗ interBRNSALEFoam solver MISSING"
    fi

    echo ""
    echo "============================================"
    echo "BUILD SUCCESSFUL!"
    echo "============================================"
    exit 0
else
    echo ""
    echo "============================================"
    echo "BUILD FAILED!"
    echo "============================================"
    echo "Check build.log for errors"
    exit 1
fi
EOF
)

# Replace placeholders
REBUILD_SCRIPT="${REBUILD_SCRIPT//OPENFOAM_PATH_PLACEHOLDER/$OPENFOAM_PATH}"
REBUILD_SCRIPT="${REBUILD_SCRIPT//HPC_PATH_PLACEHOLDER/$HPC_PATH}"

# Execute rebuild on HPC
echo "Connecting to HPC and building..."
echo ""

ssh ${HPC_USER}@${HPC_HOST} "bash -s" <<< "$REBUILD_SCRIPT"

if [ $? -eq 0 ]; then
    echo ""
    echo -e "${GREEN}============================================${NC}"
    echo -e "${GREEN}Sync and Rebuild Complete!${NC}"
    echo -e "${GREEN}============================================${NC}"
else
    echo ""
    echo -e "${RED}============================================${NC}"
    echo -e "${RED}Rebuild Failed!${NC}"
    echo -e "${RED}============================================${NC}"
    echo "Check build.log on HPC: ${HPC_USER}@${HPC_HOST}:~/${HPC_PATH}/build.log"
    exit 1
fi
