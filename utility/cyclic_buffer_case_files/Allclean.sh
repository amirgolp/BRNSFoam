#!/bin/bash
#
# Allclean.sh
#
# Reset the case to a clean pre-mesh state.  Preserves every checked-in
# template -- system/*Dict, constant/transportProperties, the 0/ field
# templates with their BC definitions, constant/triSurface/*.py and the
# source STL -- and removes everything the createMesh + runSnappyHexMesh
# pipeline generates.
#
# After this, `./Allrun.sh` rebuilds the mesh and runs the solver from
# scratch.
#
# Usage:
#   ./Allclean.sh           # standard clean (keep 0/ as-is)
#   ./Allclean.sh --deep    # also overwrite 0/ from 0.orig/
#
# --deep is useful when 0/ has been corrupted (e.g. partial subsetMesh run
# left 0/<field> with the wrong nonuniform size).  0.orig/ is the snapshot
# runSnappyHexMesh.sh saves after a successful setFields.

set -u
cd "$(dirname "$0")"

deep=0
if [ "${1:-}" = "--deep" ] || [ "${1:-}" = "-d" ]; then
    deep=1
fi

echo "[Allclean] Removing generated mesh and runtime artefacts..."

# --- Mesh -------------------------------------------------------------------
rm -rf constant/polyMesh
rm -rf constant/extendedFeatureEdgeMesh

# Generated STL (createMesh re-translates/clips fresh each run).
rm -f  constant/triSurface/Image_meshed.stl
rm -f  constant/triSurface/Image_meshed_buffered.stl

# --- Time directories (numeric, > 0) ---------------------------------------
# Keep 0/ (BC template, mandatory) and 0.orig/ (post-setFields snapshot).
for d in [0-9]*; do
    [ -d "$d" ] || continue
    case "$d" in
        0|0.orig) continue ;;
        *)        rm -rf "$d" ;;
    esac
done

# --- Parallel decomposition -----------------------------------------------
rm -rf processor*

# --- Generated controlDict / blockMeshDict ---------------------------------
# Only remove if a template exists, so we don't kill a hand-edited dict
# that the user has not yet baked into controlDictRun / blockMeshDict2D*.
[ -f system/controlDictRun ]   && rm -f system/controlDict
[ -f system/blockMeshDict2D0 ] && rm -f system/blockMeshDict

# --- Logs and crud ---------------------------------------------------------
rm -f  snappyHexMesh.out createPatch.out
rm -f  checkMesh.out checkMesh_post_subset.out
rm -f  subsetMesh.out setFields.out
rm -f  clip_buffer_stl.out
rm -f  surfaceTransformPoints1.out surfaceTransformPoints2.out
rm -f  decomposePar.log
rm -f  Allrun_createMesh.log Allrun_runSnappy.log Allrun_solver.log
rm -f  Allrun_ami_probe.log
rm -f  log.* *.foam
rm -f  core core.*
rm -rf postProcessing VTK

# --- Optional deep clean of 0/ ---------------------------------------------
if [ "$deep" -eq 1 ]; then
    if [ -d 0.orig ]; then
        echo "[Allclean] --deep: restoring 0/ from 0.orig/"
        rm -rf 0
        cp -r 0.orig 0
    else
        echo "[Allclean] --deep requested but no 0.orig/ to restore from."
        echo "           0/ left untouched."
    fi
fi

echo "[Allclean] Done."
