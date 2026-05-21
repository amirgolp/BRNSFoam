#!/bin/bash
#
# runSnappyHexMesh.sh
#
# Runs the full mesh-finishing pipeline after createMesh_cyclicBuffer.sh:
#   1. snappyHexMesh        -- carve grain cells from the buffered STL
#   2. createPatch          -- inlet/outlet patch -> cyclicAMI
#   3. checkMesh            -- topology check, writes per-region cellSets
#   4. subsetMesh           -- drop disconnected pore islands (if any)
#   5. setFields            -- (re-)lay the alpha.water stripe IC
#   6. cp -r 0 0.orig       -- baseline snapshot
#
# Steps 3-6 used to be done by hand and were repeatedly forgotten, with
# the following failure modes:
#   - skip (3) + (4)  ==>  multi-region pressure system is singular,
#                          GAMG p_rgh coarsest-level DIC FPEs at step 1
#   - skip (5)        ==>  alpha.water defaults to uniform 0.75 after a
#                          mesh rebuild, GAMG explodes on the degenerate
#                          interface-free initial pressure system
#   - skip (6)        ==>  no clean baseline to restore for re-spinups
#
# All four are now unconditional after snappy+createPatch.

set -e

# Pick parallel vs serial automatically.  createMesh_cyclicBuffer.sh
# deliberately leaves no processor* dirs, so we'll be in the serial branch
# unless the user has decomposed manually.
if [ -d "processor0" ]; then
    NP="$(find processor* -maxdepth 0 -type d -print | wc -l)"
    echo "Run snappyHexMesh in parallel on $NP processors"
    mpirun -np $NP snappyHexMesh -overwrite -parallel > snappyHexMesh.out
    echo "Run createPatch (parallel) to convert inlet/outlet to cyclicAMI"
    mpirun -np $NP createPatch -overwrite -parallel > createPatch.out
else
    echo "Run snappyHexMesh"
    snappyHexMesh -overwrite > snappyHexMesh.out
    echo "Run createPatch to convert inlet/outlet to cyclicAMI"
    createPatch -overwrite > createPatch.out
fi

# Step 3: topology check.  -allTopology writes region cellSets to
# constant/polyMesh/sets/region{0,1,...} for any disconnected pore islands.
# `|| true` because checkMesh exits non-zero on benign mesh-quality warnings
# (failed faces, AR, skewness) which would otherwise trip `set -e` and skip
# subsetMesh / setFields below.  We only need the region-count string here.
echo "Run checkMesh -allTopology"
checkMesh -allTopology > checkMesh.out 2>&1 || true
n_regions=$(grep -E "Number of regions:" checkMesh.out | awk '{print $4}')
echo "  checkMesh: ${n_regions:-?} region(s)"

# Step 4: drop disconnected pore islands (keep region0 = the main pore
# network that contains snappy's locationInMesh).  Skipped if there's
# only one region.
if [ -f constant/polyMesh/sets/region1 ]; then
    # subsetMesh has to READ every 0/<field> against the current
    # (pre-subset) mesh before it can subset them.  If 0/<field>'s
    # nonuniform list has a different cell count -- e.g. left over from a
    # previous subset run -- subsetMesh aborts with "size N is not equal
    # to the expected length M".
    #
    # Fix: collapse every 0/<field>'s internalField AND every
    # boundaryField/<patch>/value entry to a uniform default before subsetMesh
    # runs.  Uniform values have no explicit size and subset trivially.
    # setFields (step 5) then re-lays alpha.water with the stripe IC at the
    # post-subset cell count, and the solver's first PIMPLE iteration refills
    # any boundary value lists that the BC type cares about.
    #
    # Resetting only internalField is not enough: BCs that carry a nonuniform
    # `value` list (e.g. constantAlphaContactAngle on solidwalls) still match
    # the OLD patch face count.  setFields fatals at "size N is not equal to
    # the expected length M" when it tries to read the patch values, before
    # painting any cell -- so alpha.water silently stays at uniform 0 and the
    # solver runs with no water in the domain.
    echo "Reset 0/ internalField and boundaryField values to uniform (subsetMesh-safe)"
    for f in 0/*; do
        [ -f "$f" ] || continue
        cls=$(awk '/^FoamFile$/{flag=1} flag && /^}/{exit} flag' "$f" \
              | grep -oP 'class\s+\K\w+' | head -1)
        case "$cls" in
            volScalarField|surfaceScalarField)
                uval="uniform 0"
                ;;
            volVectorField|surfaceVectorField)
                uval="uniform (0 0 0)"
                ;;
            labelList)
                # cellToRegion etc. -- stale per-cell label list, just delete
                rm -f "$f"
                continue
                ;;
            *)
                continue
                ;;
        esac
        foamDictionary -entry internalField -set "$uval" "$f" > /dev/null
        # Walk every patch in boundaryField; if it has a `value` entry, force
        # it uniform.  Patches without a `value` entry (zeroGradient, empty,
        # ...) are skipped silently.
        patches=$(foamDictionary -entry boundaryField -keywords "$f" 2>/dev/null) || continue
        for p in $patches; do
            if foamDictionary -entry "boundaryField/${p}/value" "$f" >/dev/null 2>&1; then
                foamDictionary -entry "boundaryField/${p}/value" -set "$uval" "$f" > /dev/null
            fi
        done
    done

    echo "Run subsetMesh -overwrite region0 (drop disconnected pore islands)"
    subsetMesh -overwrite region0 > subsetMesh.out 2>&1
    # Re-checkMesh post-subset to confirm one region and to refresh cellSets.
    # `|| true` for the same reason as the first checkMesh above -- a non-zero
    # exit on quality warnings would skip setFields and leave alpha.water at
    # uniform 0, which silently launches the solver with no water in the
    # domain (interFoam reports Phase-1 volume fraction = 0).
    checkMesh -allTopology > checkMesh_post_subset.out 2>&1 || true
    n_regions=$(grep -E "Number of regions:" checkMesh_post_subset.out | awk '{print $4}')
    echo "  after subsetMesh: ${n_regions:-?} region(s)"
else
    echo "  only one region; subsetMesh not needed"
fi

# Step 5: lay down the alpha.water stripe IC for the (possibly subsetted)
# cell layout.  setFields rewrites 0/alpha.water in-place per setFieldsDict.
echo "Run setFields (alpha.water stripe IC per setFieldsDict)"
setFields > setFields.out 2>&1
sw=$(grep -oP "Phase[- ]?1 volume fraction = \K[0-9.e+-]+" setFields.out | tail -1 || true)
echo "  setFields done${sw:+ (Sw ~ $sw)}"

# Step 6: baseline snapshot.  Overwrites any prior 0.orig.
echo "Save 0/ baseline -> 0.orig"
rm -rf 0.orig
cp -r 0 0.orig

echo
echo "Mesh + IC ready.  Inspect in ParaView before running flow."
echo "  Sw (volume fraction)            : see setFields.out"
echo "  cyclic match (AMI weights)      : see createPatch.out"
echo "  topology (1 region after subset): see checkMesh*.out"
