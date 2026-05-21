#!/usr/bin/env python3
"""
clip_buffer_stl.py
==================

Drop-in replacement for `make_periodic_stl.py`.  Carves a pure-pore slab of
width `buffer_x` at each x-end of the simulation domain by:

  1. Loading the source grain-surface STL (already translated/scaled into
     the blockMesh frame by upstream `surfaceTransformPoints` calls).
  2. Splitting the STL into connected components (one per grain).
  3. For each component:
       - Drop if entirely outside the blockMesh x-range [0, Lx].
       - Keep intact if entirely inside the keep window [buffer_x,
         Lx - buffer_x].
       - Otherwise: slice with the plane(s) it crosses (keep the in-keep-
         window side), CAPPING each cut.  Capping is done per-component so
         the cap polygon for each grain is just that grain's cross-section
         -- no long earcut "bridge" triangles connecting different grains'
         cross-sections.
  4. Concatenate the surviving pieces and write the result.

Why per-component capping
-------------------------
The earlier monolithic `slice_mesh_plane(whole_mesh, cap=True)` triangulates
the union of ALL cut cross-sections at the buffer plane as a single complex
polygon.  Earcut bridges across the gaps between separate grains, producing
long cap triangles that span from one grain's cut profile to another's --
visible in ParaView as straight wall segments connecting otherwise-distinct
grain outlines.  Snappy treats these as solid walls, putting phantom wall
connectors between physically separate grains.

Per-component slicing closes each grain's cut with a polygon whose only
boundary is that one grain's cross-section -- no inter-grain bridges.

Dependencies
------------
trimesh, numpy, scipy (for connected-components labelling via face-edge
adjacency), and trimesh's cap=True triangulation chain (networkx, shapely,
mapbox-earcut).  setup_venv.sh installs all of them.

Usage
-----
::

    python3 clip_buffer_stl.py INPUT.stl OUTPUT.stl Lx buffer_x

Lx        : simulation domain length in x (metres).  Must match blockMesh.
buffer_x  : pore-only slab width at each x-end (metres).  Typically one
            blockMesh cell pitch.  buffer_x = 0 disables clipping
            entirely.
"""
import os
import sys

try:
    import trimesh
    from trimesh.intersections import slice_mesh_plane
except ImportError as exc:
    sys.stderr.write(
        "Missing core dependency: " + str(exc) + "\n"
        "Run setup_venv.sh first.\n"
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
    n_in = len(m.faces)
    print(f"      triangles: {n_in}")
    print(f"      bounds:    {m.bounds[0]}  ->  {m.bounds[1]}")

    if buffer_x == 0:
        print("[2/4] buffer_x = 0; skipping clip step.")
        m.export(output_stl)
        print(f"[3/4] Wrote {output_stl} unchanged.")
        return

    print(f"[2/4] Splitting into connected components (one per grain) ...")
    components = m.split(only_watertight=False)
    print(f"      components: {len(components)}")

    print(f"[3/4] Per-component slice+cap (right at x={Lx-buffer_x:g}, "
          f"left at x={buffer_x:g}) ...")
    # Slightly looser than float32 STL precision so a vertex landing at
    # exactly buffer_x rounds INSIDE the keep window rather than triggering
    # an unnecessary slice that produces a degenerate cap.
    eps = 1e-9

    kept = []
    n_intact = n_cut = n_drop_outside = n_drop_empty = 0
    tri_intact = tri_cut = 0
    for c in components:
        cmin, cmax = float(c.bounds[0][0]), float(c.bounds[1][0])
        # Entirely outside blockMesh x-range: drop to keep output STL small
        # (snappy ignores them anyway).
        if cmax < -eps or cmin > Lx + eps:
            n_drop_outside += 1
            continue
        # Entirely inside the keep window: pass through unchanged.
        if cmin >= buffer_x - eps and cmax <= Lx - buffer_x + eps:
            kept.append(c)
            n_intact += 1
            tri_intact += len(c.faces)
            continue
        # Otherwise: at least one buffer plane crosses this component.
        # Slice and cap, per-component, so the cap polygon is just this
        # grain's cross-section -- no cross-grain bridges.
        piece = c
        if cmax > Lx - buffer_x + eps:
            try:
                piece = slice_mesh_plane(
                    piece,
                    plane_normal=[-1.0, 0.0, 0.0],
                    plane_origin=[float(Lx - buffer_x), 0.0, 0.0],
                    cap=True,
                )
            except Exception as exc:
                sys.stderr.write(
                    f"  right slice failed on a component (xmin={cmin*1e6:.2f} um, "
                    f"xmax={cmax*1e6:.2f} um): {exc}\n"
                )
                continue
        if piece is None or len(piece.faces) == 0:
            n_drop_empty += 1
            continue
        if cmin < buffer_x - eps:
            try:
                piece = slice_mesh_plane(
                    piece,
                    plane_normal=[1.0, 0.0, 0.0],
                    plane_origin=[float(buffer_x), 0.0, 0.0],
                    cap=True,
                )
            except Exception as exc:
                sys.stderr.write(
                    f"  left slice failed on a component (xmin={cmin*1e6:.2f} um, "
                    f"xmax={cmax*1e6:.2f} um): {exc}\n"
                )
                continue
        if piece is None or len(piece.faces) == 0:
            n_drop_empty += 1
            continue
        kept.append(piece)
        n_cut += 1
        tri_cut += len(piece.faces)

    print(f"      components: intact={n_intact}  sliced+capped={n_cut}  "
          f"dropped_outside={n_drop_outside}  dropped_empty={n_drop_empty}")
    print(f"      triangles : intact={tri_intact}  sliced+capped={tri_cut}")
    if not kept:
        sys.stderr.write("All components dropped.  Check Lx / buffer_x.\n")
        sys.exit(2)

    print(f"[4/4] Concatenating {len(kept)} surviving pieces and writing "
          f"{output_stl}")
    merged = trimesh.util.concatenate(kept)
    merged.export(output_stl)
    print(f"      wrote {os.path.abspath(output_stl)} ({len(merged.faces)} "
          f"triangles)")
    print(f"      final bounds: {merged.bounds[0]}  ->  {merged.bounds[1]}")
    print(f"      preserved {len(merged.faces)} of {n_in} input triangles "
          f"({100.0 * len(merged.faces) / n_in:.1f}%)")


def main(argv):
    if len(argv) != 5:
        sys.stderr.write(__doc__)
        sys.exit(1)
    clip_buffer(argv[1], argv[2], float(argv[3]), float(argv[4]))


if __name__ == '__main__':
    main(sys.argv)
