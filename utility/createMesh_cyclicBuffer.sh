#!/bin/bash
#
# createMesh_cyclicBuffer.sh
# ==========================
#
# Build a mesh suitable for cyclic-x BCs by combining:
#   - the source grain STL (preserved exactly inside the keep-window), with
#   - a pure-pore slab of width $buffer_x at each x-end.
#
# This is a drop-in replacement for `channelStationaryS75/createMesh.sh` and
# uses `clip_buffer_stl.py` instead of `make_periodic_stl.py`.  The voxel
# round-trip in make_periodic_stl.py was destroying the discrete-pebble
# structure (merging grains, pinching off pore throats, leaving disconnected
# pore islands -- see checkMesh "Number of regions: 2" on that case).
#
# Run from the case root.  After it finishes, run runSnappyHexMesh.sh and
# then `createPatch -overwrite` (with the existing createPatchDict) to
# convert inlet/outlet from `patch` -> `cyclic`.
#
# Required co-located files inside the case:
#   constant/triSurface/clip_buffer_stl.py     (copy from utility/)
#   constant/triSurface/raw2stl.py             (only if format=='raw')
#   system/blockMeshDict2D0                    (template with dx/dy/dz/nx/ny/nz/res placeholders)
#   system/snappyHexMeshDict2D                 (template with poreIndex0/1/2 placeholders)
#   system/controlDictInit                     (mesh-build controlDict)
#
# Required Python venv with: numpy, scipy, trimesh, numpy-stl
# (use the same setup_venv.sh as channelStationaryS75 -- no new deps).

###### USERS INPUT ############################################################

Image_name="Grainstones"
dir="$BRNSFOAM_DIR/images/stl"
format='stl'
compressed='yes'

pore_index_X=0.0001522
pore_index_Y=0.0001522
pore_index_Z=0

res=1

# Buffer-zone setup.  Extend blockMesh by buffer_x on each x-end so the inlet
# and outlet patches sit in pure-pore slabs.
# Domain 400 um interior + 2 um at each end = 404 um, n_x = 202 (pitch 2 um).
x_min=-0.000002
x_max=0.000402
y_min=0
y_max=0.00025
z_min=-0.000005
z_max=0.000005

n_x=202
n_y=125
n_z=1

# Physical width of the pure-pore buffer at each x-end (m).  Must be a
# nonnegative multiple of the blockMesh cell pitch -- e.g. 2e-6 for one cell
# at 2 um pitch.  Set to 0 to skip the buffer step entirely (no cyclic).
buffer_x=2e-6

direction=0
NP=6

#### END OF USER INPUT #######################################################

set -e
source ~/.bashrc

if [ -z "${BRNSFOAM_DIR-}" ] || [ ! -d "$BRNSFOAM_DIR/images/stl" ]; then
    echo "ERROR: image library $BRNSFOAM_DIR/images/stl not reachable."
    echo "       BRNSFOAM_DIR=$BRNSFOAM_DIR.  Did etc/bashrc source succeed?"
    exit 1
fi

VENV_DIR="${VENV_DIR:-$HOME/venvs/openfoam_py}"
if [ -f "$VENV_DIR/bin/activate" ]; then
    # shellcheck disable=SC1091
    source "$VENV_DIR/bin/activate"
    if [ -n "${PYTHONPATH-}" ]; then
        unset PYTHONPATH
    fi
    echo "Activated venv: $VENV_DIR"
else
    echo "ERROR: $VENV_DIR not found.  Run setup_venv.sh first." >&2
    exit 1
fi

if [ $format != 'raw' ] && [ $format != 'stl' ]; then
    echo "ERROR: only raw and stl format are implemented"
    exit 1
fi

filename=$Image_name\.$format
if [ $compressed == 'yes' ]; then
    filename=$Image_name\.$format\.tar.gz
fi

cp $dir/$filename constant/triSurface/.
cd constant/triSurface
if [ $compressed == 'yes' ]; then
    tar -xf $filename
fi

if [ $format == 'raw' ]; then
    echo "make stl"
    python raw2stl.py --x_min=$x_min --x_max=$x_max --y_min=$y_min --y_max=$y_max \
                      --z_min=$z_min --z_max=$z_max \
                      --pores_value=$pores_value --solid_value=$solid_value \
                      --image_name=$Image_name \
                      --x_dim=$x_dim --y_dim=$y_dim --z_dim=$z_dim
    rm $Image_name\.*
    surfaceTransformPoints -translate '(-0.5 -0.5 -0.5)' Image_meshed.stl Image_meshed.stl > ../../surfaceTransformPoints1.out
    vector="($res $res $res)"
    surfaceTransformPoints -scale "$vector" Image_meshed.stl Image_meshed.stl > ../../surfaceTransformPoints2.out
    pore_index_X="$(cat pore_indx)"
    pore_index_Y="$(cat pore_indy)"
    pore_index_Z="$(cat pore_indz)"
    pore_index_0=$(expr $pore_index_X*$res | bc)
    pore_index_1=$(expr $pore_index_Y*$res | bc)
    pore_index_2=$(expr $pore_index_Z*$res | bc)
elif [ $format == 'stl' ]; then
    echo "reading stl"
    mv $Image_name\.stl Image_meshed.stl
    rm -f $Image_name\.*
    tranX=$(expr -1*$x_min | bc)
    tranY=$(expr -1*$y_min | bc)
    tranZ=$(expr -1*$z_min | bc)
    vector="($tranX $tranY $tranZ)"
    surfaceTransformPoints -translate "$vector" Image_meshed.stl Image_meshed.stl > ../../surfaceTransformPoints1.out
    vector="($res $res $res)"
    surfaceTransformPoints -scale "$vector" Image_meshed.stl Image_meshed.stl > ../../surfaceTransformPoints2.out
    pore_index_0=$(expr $res*$pore_index_X-$res*$x_min | bc)
    pore_index_1=$(expr $res*$pore_index_Y-$res*$y_min | bc)
    pore_index_2=$(expr $res*$pore_index_Z-$res*$z_min | bc)
fi

# === BUFFER-SLAB STEP ========================================================
# After translate+scale the STL is in the blockMesh frame:
#   x in [0, dx_phys],  y in [0, dy_phys],  z in [0, dz_phys].
# Drop grain components that overlap the buffer slabs so the inlet (x=0) and
# outlet (x=Lx) patches end up as identical full-cross-section pore faces.
# clip_buffer_stl.py handles buffer_x=0 gracefully (only drops grains entirely
# outside blockMesh), so the call is unconditional.
dx_phys=$(echo "($x_max - $x_min) * $res" | bc -l)

echo "Building buffered STL (Lx=$dx_phys m, buffer_x=$buffer_x m)"
python3 clip_buffer_stl.py Image_meshed.stl Image_meshed_buffered.stl \
        $dx_phys $buffer_x \
        > ../../clip_buffer_stl.out
if [ ! -s Image_meshed_buffered.stl ]; then
    echo "ERROR: clip_buffer_stl.py produced no output. See clip_buffer_stl.out."
    exit 1
fi
mv -f Image_meshed_buffered.stl Image_meshed.stl
echo "Buffered STL ready."
# =============================================================================

echo "Coordinates at center of a pore = ($pore_index_0,$pore_index_1,$pore_index_2)"

cd ../..

echo "Create background mesh"
cp system/blockMeshDict2D$direction system/blockMeshDict
dx=$(expr $x_max-1*$x_min | bc)
dy=$(expr $y_max-1*$y_min | bc)
dz=$(expr $z_max-1*$z_min | bc)

sed -i "s/dx/$dx/g"   system/blockMeshDict
sed -i "s/dy/$dy/g"   system/blockMeshDict
sed -i "s/dz/$dz/g"   system/blockMeshDict
sed -i "s/nx/$n_x/g"  system/blockMeshDict
sed -i "s/ny/$n_y/g"  system/blockMeshDict
sed -i "s/nz/$n_z/g"  system/blockMeshDict
sed -i "s/res/$res/g" system/blockMeshDict

cp system/controlDictInit system/controlDict
blockMesh > blockMesh.out

# Restore the flow-time controlDict if a template is present so the user
# isn't left with controlDictInit's reactiveTransportFoam settings after
# meshing (interBRNSFoam will fail with "Entry 'maxAlphaCo' not found"
# otherwise).
if [ -f system/controlDictRun ]; then
    cp system/controlDictRun system/controlDict
    echo "Restored flow-time system/controlDict from system/controlDictRun."
else
    echo "WARNING: system/controlDictRun not found.  system/controlDict still"
    echo "         holds the meshing dict (controlDictInit form).  Edit it"
    echo "         before running the solver, or drop in a controlDictRun."
fi

cp system/snappyHexMeshDict2D system/snappyHexMeshDict
sed -i "s/poreIndex0/$pore_index_0/g" system/snappyHexMeshDict
sed -i "s/poreIndex1/$pore_index_1/g" system/snappyHexMeshDict
sed -i "s/poreIndex2/$pore_index_2/g" system/snappyHexMeshDict

# Wipe any stale processor* directories from a previous mesh build.  Leaving
# them in place caused snappyHexMesh -parallel to read a previous snapped
# mesh (with createPatch-assigned cyclicAMI patches), then segfault inside
# faceAreaWeightAMI::calculate when its AMI weights no longer matched the
# new geometry.  After this wipe, runSnappyHexMesh.sh sees no processor*
# dirs and runs snappy serially on the fresh constant/polyMesh.
if compgen -G "processor*" > /dev/null; then
    echo "Removing stale processor* directories (snappy will run serially)."
    rm -rf processor*
fi

# Intentionally NOT calling decomposePar here.  The case's 0/ field BCs are
# already typed `cyclicAMI` (post-createPatch state), and decomposePar would
# fail with "Attempt to cast type patch to type lduInterface" against the
# fresh blockMesh's plain `patch` inlet/outlet.  Snappy runs in serial on
# ~25k cells in well under a minute, so the parallel path isn't worth the
# field-type juggling at this stage.

echo "BlockMesh and clipped STL ready. Next:"
echo "    ./runSnappyHexMesh.sh       # serial snappy + serial createPatch"
echo "    checkMesh -allTopology      # expect 'Number of regions: 1' this time"
echo "    # then decomposePar to spread fields onto the cyclic mesh,"
echo "    # and run the flow solver in parallel."
