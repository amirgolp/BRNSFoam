#!/bin/bash
#
# Allrun.sh
#
# End-to-end pipeline for the channelStationary cyclic-buffer case:
#   1. createMesh_cyclicBuffer.sh   blockMesh background + per-component
#                                   STL clip (clip_buffer_stl.py)
#   2. runSnappyHexMesh.sh          snappy + createPatch + subsetMesh +
#                                   setFields + 0.orig snapshot
#   3. <solver>                     interFoam by default; override with
#                                   $SOLVER (e.g. SOLVER=interBRNSFoam)
#
# After step 2 the script reports the cyclicAMI weight quality.  Anything
# outside [0.999, 1.001] means the snap step is reaching across the buffer
# slab and pulling inlet/outlet vertices asymmetrically -- widen buffer_x
# in createMesh_cyclicBuffer.sh and re-Allclean+Allrun.
#
# Environment overrides:
#   SOLVER    name of OpenFOAM solver binary  (default: interFoam)
#   NP        number of MPI ranks             (default: 1, serial)
#
# Examples:
#   ./Allrun.sh                          # serial interFoam
#   SOLVER=interBRNSFoam ./Allrun.sh     # serial interBRNSFoam
#   NP=8 ./Allrun.sh                     # parallel interFoam on 8 ranks
#
# Logs each step to Allrun_*.log in the case root.

set -e
set -u
cd "$(dirname "$0")"

SOLVER="${SOLVER:-interFoam}"
NP="${NP:-1}"

# --- Preconditions ---------------------------------------------------------
# 0/ must contain BC templates for every field.  setFields rewrites their
# internalField -- it does not create the files.  Without templates the
# solver loads uniform 0 alpha.water (no water).  Restore from 0.orig if
# 0/ is missing.
if [ ! -d 0 ] || [ -z "$(ls -A 0 2>/dev/null)" ]; then
    if [ -d 0.orig ]; then
        echo "[Allrun] 0/ missing or empty -- restoring from 0.orig/"
        rm -rf 0
        cp -r 0.orig 0
    else
        echo "[Allrun] ERROR: neither 0/ nor 0.orig/ exists -- cannot run."
        echo "                Ensure field templates are checked into the case."
        exit 1
    fi
fi

for f in createMesh_cyclicBuffer.sh runSnappyHexMesh.sh; do
    if [ ! -x "./$f" ]; then
        echo "[Allrun] ERROR: ./$f missing or not executable."
        exit 1
    fi
done

# --- Step 1: blockMesh background + per-component STL clip -----------------
echo "[Allrun] (1/3) createMesh_cyclicBuffer.sh"
./createMesh_cyclicBuffer.sh 2>&1 | tee Allrun_createMesh.log

# --- Step 2: snappy + createPatch + subsetMesh + setFields ----------------
echo "[Allrun] (2/3) runSnappyHexMesh.sh"
./runSnappyHexMesh.sh 2>&1 | tee Allrun_runSnappy.log

# Note on AMI weight quality:
#   cyclicAMI weights are computed lazily by the FIRST utility that
#   constructs cyclicAMI patches from the mesh.  Empirically that's the
#   solver itself (or decomposePar for parallel) -- createPatch only
#   rewrites the boundary-file patch types, and checkMesh / setFields /
#   subsetMesh don't trigger AMI construction.  So we cannot pre-probe
#   weights here; they will appear in the solver's startup output
#   (Allrun_solver.log) within the first second as:
#       AMI: Patch source sum(weights) min:1 max:1 average:1
#       AMI: Patch target sum(weights) min:...     max:...    average:...
#   Perfect cyclic match gives min=max=1.0.  Anything > 0.001 off means
#   snap is reaching across the buffer slab -- widen buffer_x in
#   createMesh_cyclicBuffer.sh, or drop snapControls.tolerance toward 1.0.

# --- Step 3: solver --------------------------------------------------------
if [ "$NP" -gt 1 ]; then
    echo "[Allrun] (3/3) decomposePar -> $SOLVER -parallel on $NP ranks"
    rm -rf processor*
    decomposePar -force > decomposePar.log 2>&1
    mpirun -np "$NP" "$SOLVER" -parallel 2>&1 | tee Allrun_solver.log
else
    echo "[Allrun] (3/3) $SOLVER (serial)"
    "$SOLVER" 2>&1 | tee Allrun_solver.log
fi

# --- Post-run summary: AMI weights from solver startup ---------------------
# After the solver runs (or is interrupted), pull AMI weights from its log.
# This is the most reliable place to find them -- see comment at step 2a.
if [ -f Allrun_solver.log ]; then
    src_line=$(grep -m1 "AMI: Patch source sum(weights)" Allrun_solver.log || true)
    tgt_line=$(grep -m1 "AMI: Patch target sum(weights)" Allrun_solver.log || true)
    if [ -n "$src_line" ]; then
        echo ""
        echo "[Allrun] cyclicAMI weight quality (from solver startup):"
        echo "         source  ${src_line#*AMI: Patch source }"
        echo "         target  ${tgt_line#*AMI: Patch target }"
    fi
fi

echo "[Allrun] Done."
