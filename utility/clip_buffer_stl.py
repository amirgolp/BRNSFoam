#!/usr/bin/env python3
"""
clip_buffer_stl.py
==================

Drop-in replacement for `make_periodic_stl.py` (the voxelize-and-marching-cubes
buffer builder).  Carves a pure-pore slab of width `buffer_x` at each x-end of
the simulation domain WITHOUT round-tripping the geometry through a voxel grid.

What it does
------------
1. Load the source grain-surface STL (already translated/scaled into the
   blockMesh frame by `surfaceTransformPoints` upstream).
2. Split into connected components via face-edge adjacency
   (`scipy.sparse.csgraph.connected_components`).  Each component is one
   grain/pebble.
3. For every component, look at its x-bounding-box and DROP it if it:
      - lies entirely outside the blockMesh x-range  [0, Lx], OR
      - overlaps the left buffer slab               [0, buffer_x], OR
      - overlaps the right buffer slab              [Lx - buffer_x, Lx].
   Surviving components lie entirely inside [buffer_x, Lx - buffer_x].
4. Write the surviving triangles back out as a binary STL.

Why this and not make_periodic_stl.py
-------------------------------------
`make_periodic_stl.py` voxelizes the STL onto a finite lattice and re-extracts
a surface with marching cubes.  At any voxel pitch that isn't very fine
(< ~0.2 µm for the Grainstones STL), this:
  - merges adjacent pebbles whose grain-grain contact is sub-voxel,
  - pinches off narrow pore throats < 2 voxels wide,
  - creates disconnected pore islands (checkMesh: "Number of regions: 2").
On the Grainstones STL at 1 µm voxel pitch we measured pore fraction jumping
from 42 % (raw STL) to 64 % (voxelized) — large enough to destroy the discrete-
pebble structure visible in the source image.

This script preserves the source geometry exactly inside the keep-window.  The
only distortion is at the buffer-slab boundaries, where grains straddling
those boundaries are dropped wholesale (no clipping, no capping).  Dropping
whole components is the only way to guarantee the buffer slab is PURE pore,
which in turn is what makes the inlet/outlet patches identical full-cross-
section faces — the prerequisite for plain `cyclic` BCs without AMI.

Limits
------
- Components that span the keep-window discontinuously in x (e.g. two pebbles
  connected by a thin bridge that crosses the buffer) will be classified as
  one component and dropped if any part overlaps the buffer.  For the
  Grainstones STL each pebble is one component, so this is a non-issue.
- The script does not modify normal orientation.  snappyHexMesh uses
  locationInMesh flood-fill so normal orientation does not matter.

Dependencies
------------
numpy, scipy, trimesh, numpy-stl.  Exactly what `setup_venv.sh` already
installs — no new deps.

Usage
-----
::

    python3 clip_buffer_stl.py INPUT.stl OUTPUT.stl Lx buffer_x

Lx        : simulation domain length in x (metres).  Must match blockMesh.
buffer_x  : pore-only slab width at each x-end (metres).  Typically one or
            two blockMesh cell pitches.

Example for the channelStationaryS75 case (Lx = 404 µm, buffer = 2 µm)::

    python3 clip_buffer_stl.py Image_meshed.stl Image_meshed_buffered.stl \\
        4.04e-4 2e-6
"""
import os
import sys

try:
    import numpy as np
    import trimesh
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components
    from stl import mesh as stlmesh
except ImportError as exc:
    sys.stderr.write(
        "Missing dependency: " + str(exc) + "\n"
        "Run setup_venv.sh first, or install: numpy scipy trimesh numpy-stl\n"
    )
    sys.exit(2)


def clip_buffer(input_stl, output_stl, Lx, buffer_x):
    if buffer_x < 0:
        sys.stderr.write("buffer_x must be >= 0.\n")
        sys.exit(1)
    if 2 * buffer_x >= Lx:
        sys.stderr.write(f"buffer_x ({buffer_x}) too large for Lx ({Lx}).\n")
        sys.exit(1)

    inp_abs = os.path.abspath(input_stl)
    print(f"[1/4] Loading STL: {inp_abs}")
    m = trimesh.load(input_stl, force='mesh')
    if not hasattr(m, 'faces') or len(m.faces) == 0:
        sys.stderr.write("Input is not a single non-empty triangulated mesh.\n")
        sys.exit(1)
    n_faces = len(m.faces)
    print(f"      triangles: {n_faces}")
    print(f"      bounds:    {m.bounds[0]}  ->  {m.bounds[1]}")

    print(f"[2/4] Building face-adjacency graph and labelling components ...")
    fa = m.face_adjacency  # (E, 2) pairs of face indices sharing an edge
    if len(fa) == 0:
        sys.stderr.write("STL has no shared edges (all triangles disjoint).  "
                         "Cannot identify components.\n")
        sys.exit(2)
    data = np.ones(len(fa) * 2, dtype=np.int8)
    rows = np.concatenate([fa[:, 0], fa[:, 1]])
    cols = np.concatenate([fa[:, 1], fa[:, 0]])
    g = csr_matrix((data, (rows, cols)), shape=(n_faces, n_faces))
    n_comp, labels = connected_components(g, directed=False)
    sizes = np.bincount(labels)
    print(f"      components: {n_comp}  (tri counts: min {sizes.min()}, "
          f"max {sizes.max()}, median {int(np.median(sizes))})")

    # STL stores vertices as float32, so a value translated to exactly
    # `buffer_x` round-trips as `buffer_x - O(1e-13)`.  A 1 nm tolerance is
    # well below blockMesh cell pitch and lets grains that sit *exactly* on
    # the buffer boundary count as "inside the keep window" instead of
    # "overlapping the buffer".
    eps = 1e-9
    print(f"[3/4] Filtering components against buffer slabs (eps={eps:g} m) ...")
    print(f"      keep window:        x in [{buffer_x:g}, {Lx - buffer_x:g}]")
    print(f"      left buffer  drop:  x in [0, {buffer_x:g}]")
    print(f"      right buffer drop:  x in [{Lx - buffer_x:g}, {Lx:g}]")
    keep_mask = np.zeros(n_faces, dtype=bool)
    n_keep = n_drop_left = n_drop_right = n_drop_outside = 0
    tri_keep = tri_drop_left = tri_drop_right = tri_drop_outside = 0
    for label in range(n_comp):
        face_idx = np.where(labels == label)[0]
        vert_idx = np.unique(m.faces[face_idx].flatten())
        xmin = float(m.vertices[vert_idx, 0].min())
        xmax = float(m.vertices[vert_idx, 0].max())
        # Entirely outside blockMesh -> ignored anyway, but drop to keep
        # output STL small.
        if xmax < -eps or xmin > Lx + eps:
            n_drop_outside += 1
            tri_drop_outside += len(face_idx)
            continue
        # Overlaps left buffer [0, buffer_x] -- strictly intrudes, not just
        # touches the boundary.
        overlaps_left = (xmin < buffer_x - eps) and (xmax > 0.0 + eps)
        # Overlaps right buffer [Lx - buffer_x, Lx]
        overlaps_right = (xmin < Lx - eps) and (xmax > Lx - buffer_x + eps)
        if overlaps_left:
            n_drop_left += 1
            tri_drop_left += len(face_idx)
            continue
        if overlaps_right:
            n_drop_right += 1
            tri_drop_right += len(face_idx)
            continue
        n_keep += 1
        tri_keep += len(face_idx)
        keep_mask[face_idx] = True

    print(f"      components: keep={n_keep}  drop_left={n_drop_left}  "
          f"drop_right={n_drop_right}  drop_outside={n_drop_outside}")
    print(f"      triangles : keep={tri_keep}  drop_left={tri_drop_left}  "
          f"drop_right={tri_drop_right}  drop_outside={tri_drop_outside}")
    if n_keep == 0:
        sys.stderr.write("All components dropped.  Check Lx / buffer_x.\n")
        sys.exit(2)

    print(f"[4/4] Writing {output_stl}")
    kept_faces = m.faces[keep_mask]
    out = stlmesh.Mesh(np.zeros(len(kept_faces), dtype=stlmesh.Mesh.dtype))
    for i, f in enumerate(kept_faces):
        for j in range(3):
            out.vectors[i][j] = m.vertices[f[j]]
    out.save(output_stl)
    out_abs = os.path.abspath(output_stl)
    print(f"      wrote {out_abs} ({len(kept_faces)} triangles)")
    # Report final x-bounds so the user can verify the buffer is clean.
    kx = m.vertices[np.unique(kept_faces.flatten()), 0]
    print(f"      kept-mesh x-bounds: [{kx.min():g}, {kx.max():g}]")


def main(argv):
    if len(argv) != 5:
        sys.stderr.write(__doc__)
        sys.exit(1)
    clip_buffer(argv[1], argv[2], float(argv[3]), float(argv[4]))


if __name__ == '__main__':
    main(sys.argv)
