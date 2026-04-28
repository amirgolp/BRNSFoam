#!/usr/bin/env python3
"""
mapRatesToFields.py
-------------------
Maps BRNS reaction rates from ratesAtFinish.dat into OpenFOAM volScalarField
files for every write timestep, ready for visualization in ParaView.

HOW IT WORKS
  ratesAtFinish.dat format (Fortran format '3f12.4, 512e14.7'):
      x_pos   y_pos   z_pos   r1   r2   ...   rN
  where (x_pos, y_pos, z_pos) are the real cell-centre coordinates written by
  BRNSReaction.H via mesh.C()[c]. Each line = one reacted cell per BRNS call.

  Because we have real coordinates, each line is self-describing:
    - No dependency on patch face ordering
    - Works regardless of how many reactions are defined
    - Reusable across cases with different meshes

USAGE
  # Auto-detect everything (n_rates from file, generic names rate_R1 ...):
  python3 mapRatesToFields.py /path/to/case

  # Name the rates explicitly (must match number of rate columns):
  python3 mapRatesToFields.py . --rate-names "aerobic,anoxic,decay"

  # Override rates file location:
  python3 mapRatesToFields.py . --rates-file /path/to/ratesAtFinish.dat

  # Use a config file (see --write-config):
  python3 mapRatesToFields.py . --config rateConfig.json

  # Write an example config file to edit:
  python3 mapRatesToFields.py . --write-config

OUTPUT
  For each time directory (e.g. 0.531/):
    0.531/rate_aerobic   <- volScalarField, mol/m³/s
    0.531/rate_anoxic
    0.531/rate_decay
  Open the .foam file in ParaView — these fields appear automatically.
"""

import os
import re
import sys
import json
import argparse
import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# OpenFOAM field header template
# ─────────────────────────────────────────────────────────────────────────────

OF_HEADER = (
    "/*--------------------------------*- C++ -*----------------------------------*\\\n"
    "| =========                 |                                                 |\n"
    "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
    "|  \\\\    /   O peration     | Version:  v2106                                 |\n"
    "|   \\\\  /    A nd           | Website:  www.openfoam.com                      |\n"
    "|    \\/     M anipulation  |                                                  |\n"
    "\\*---------------------------------------------------------------------------*/\n"
    "FoamFile\n"
    "{{\n"
    "    version     2.0;\n"
    "    format      ascii;\n"
    "    class       volScalarField;\n"
    "    location    \"{time}\";\n"
    "    object      {name};\n"
    "}}\n"
    "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n"
    "dimensions      [0 -3 -1 0 1 0 0];\n\n"
)

# ─────────────────────────────────────────────────────────────────────────────
# Step 1: Read ratesAtFinish.dat
# ─────────────────────────────────────────────────────────────────────────────

def read_rates_file(rates_file):
    """
    Reads ratesAtFinish.dat. Returns:
      coords : np.array shape (N, 3)  — cell centre (x, y, z) per line
      rates  : np.array shape (N, R)  — R rate values per line
      n_rates: int                    — auto-detected number of rate columns
    """
    print(f"Reading {rates_file} ...")
    rows = []
    with open(rates_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            vals = [float(v) for v in line.split()]
            rows.append(vals)

    if not rows:
        sys.exit("ERROR: ratesAtFinish.dat is empty.")

    n_cols  = len(rows[0])
    n_rates = n_cols - 3          # first 3 = x, y, z
    if n_rates < 1:
        sys.exit(f"ERROR: Expected at least 4 columns, got {n_cols}.")

    arr    = np.array(rows)
    coords = arr[:, :3]
    rates  = arr[:, 3:]
    print(f"  Lines     : {len(rows)}")
    print(f"  Reactions : {n_rates}  (auto-detected)")
    return coords, rates, n_rates


# ─────────────────────────────────────────────────────────────────────────────
# Step 2: Build coordinate → cell-id lookup from the mesh
# ─────────────────────────────────────────────────────────────────────────────

def read_cell_centres(case_dir):
    """
    Reads cell centre coordinates from constant/polyMesh/cellCentres if it
    exists (written by writeMeshObj / postProcess -func writeCellCentres),
    otherwise computes them from points + owner/neighbour connectivity.
    Returns np.array shape (nCells, 3).
    """
    cc_file = os.path.join(case_dir, "constant", "polyMesh", "cellCentres")
    if os.path.isfile(cc_file):
        return _parse_of_vector_field(cc_file)

    # Fallback: compute face centres, then average per cell
    print("  cellCentres not found — computing from faces (slower).")
    return _compute_cell_centres(case_dir)


def _parse_of_vector_field(path):
    with open(path) as f:
        raw = f.read()
    m = re.search(r'\(\s*([\s\S]+?)\s*\)', raw)
    if not m:
        sys.exit(f"ERROR: cannot parse vector field {path}")
    tuples = re.findall(r'\(\s*([^\)]+?)\s*\)', m.group(1))
    return np.array([[float(v) for v in t.split()] for t in tuples])


def _compute_cell_centres(case_dir):
    pm = os.path.join(case_dir, "constant", "polyMesh")

    def read_list(fname):
        with open(os.path.join(pm, fname)) as f:
            raw = f.read()
        m = re.search(r'\(\s*([\s\S]+?)\s*\)\s*;', raw)
        return m.group(1).split()

    # Points
    pts_raw = read_list("points")
    pts = []
    i = 0
    while i < len(pts_raw):
        token = pts_raw[i].strip("()")
        if token == "":
            i += 1; continue
        if "(" in pts_raw[i]:
            pts.append([float(pts_raw[i].strip("()")),
                        float(pts_raw[i+1].strip("()")),
                        float(pts_raw[i+2].strip("()"))])
            i += 3
        else:
            i += 1
    points = np.array(pts)

    # Faces
    with open(os.path.join(pm, "faces")) as f:
        raw = f.read()
    face_blocks = re.findall(r'(\d+)\(([^)]+)\)', raw)
    faces = [[int(v) for v in b[1].split()] for b in face_blocks]

    # Owner
    owner_raw = read_list("owner")
    owners = [int(v) for v in owner_raw]
    n_cells = max(owners) + 1

    # Average face centres per cell
    cell_sum   = np.zeros((n_cells, 3))
    cell_count = np.zeros(n_cells, dtype=int)
    for fi, face in enumerate(faces):
        fc = points[face].mean(axis=0)
        c  = owners[fi]
        cell_sum[c]   += fc
        cell_count[c] += 1

    return cell_sum / np.maximum(cell_count[:, None], 1)


def coordinates_are_degenerate(coords, tol=1e-6):
    """
    Returns True if all coordinate rows are identical within tol.
    This happens when the Fortran format '3f12.4' truncates pore-scale
    coordinates (µm range) to just 4 decimal places, making every cell
    appear at the same rounded location.
    """
    spread = np.max(np.abs(coords - coords[0]), axis=0)
    return bool(np.all(spread < tol))


def build_coord_index(cell_centres):
    """
    Builds a KD-tree for fast nearest-neighbour lookup.
    Falls back to brute-force if scipy is unavailable.
    """
    try:
        from scipy.spatial import cKDTree
        tree = cKDTree(cell_centres)
        def lookup(pts):
            _, ids = tree.query(pts, workers=-1)
            return ids
        print(f"  KD-tree built over {len(cell_centres)} cells (scipy).")
    except ImportError:
        print("  scipy not available — using brute-force lookup (slower).")
        def lookup(pts):
            ids = []
            for p in pts:
                diffs = cell_centres - p
                ids.append(int(np.argmin((diffs**2).sum(axis=1))))
            return np.array(ids)
    return lookup


# ─────────────────────────────────────────────────────────────────────────────
# Step 3a: Patch-ordering fallback (when coordinates are degenerate)
# ─────────────────────────────────────────────────────────────────────────────

def get_reacting_cells(case_dir, patch_name):
    """
    Returns (face_cells, n_faces) where face_cells[i] is the owner cell index
    of face i on the given patch — in the fixed order OpenFOAM iterates them.
    This is the same order BRNSReaction.H uses in its forAll(faceCells) loop.
    """
    pm = os.path.join(case_dir, "constant", "polyMesh")

    with open(os.path.join(pm, "boundary")) as f:
        content = f.read()

    pattern = rf'{re.escape(patch_name)}\s*\{{([^}}]+)\}}'
    m = re.search(pattern, content)
    if not m:
        sys.exit(f"ERROR: patch '{patch_name}' not found in boundary file.")

    block     = m.group(1)
    nFaces    = int(re.search(r'nFaces\s+(\d+)',    block).group(1))
    startFace = int(re.search(r'startFace\s+(\d+)', block).group(1))
    print(f"  Patch '{patch_name}': nFaces={nFaces}, startFace={startFace}")

    with open(os.path.join(pm, "owner")) as f:
        raw = f.read()
    m2 = re.search(r'\(\s*([\s\S]+?)\s*\)', raw)
    owners = [int(x) for x in m2.group(1).split()]
    face_cells = owners[startFace: startFace + nFaces]
    print(f"  → {len(face_cells)} owner cell indices loaded")
    return np.array(face_cells), nFaces


# ─────────────────────────────────────────────────────────────────────────────
# Step 3b: Coordinate-based block size detection (when coords are valid)
# ─────────────────────────────────────────────────────────────────────────────

def detect_block_size_from_coords(coords, tol=1e-6):
    """
    Auto-detect block size by finding the first repeated coordinate
    (i.e. where the same cell appears again = start of next block).
    """
    n_total = len(coords)
    ref = coords[0]
    for n in range(1, n_total):
        if np.linalg.norm(coords[n] - ref) < tol:
            if n_total % n == 0:
                print(f"  Block size auto-detected: {n} lines/call "
                      f"-> {n_total // n} BRNS calls total")
                return n
    print(f"  WARNING: Could not detect block size — treating all {n_total} "
          "lines as one block.")
    return n_total


# ─────────────────────────────────────────────────────────────────────────────
# Step 3c: Read alpha field to reconstruct the reaction conditions
# ─────────────────────────────────────────────────────────────────────────────

def read_alpha_field(case_dir, time_str, field_name, n_cells):
    """
    Parses a volScalarField (like alpha.water) to find the saturation state.
    Handles 'uniform' and 'nonuniform List<scalar>'.
    """
    path = os.path.join(case_dir, time_str, field_name)
    if not os.path.isfile(path):
        # Fallback if alpha is missing in a specific dir
        return np.ones(n_cells)

    with open(path) as f:
        content = f.read()

    # Uniform check
    m_uni = re.search(r'internalField\s+uniform\s+([0-9\.eE\+\-]+)\s*;', content)
    if m_uni:
        return np.full(n_cells, float(m_uni.group(1)))

    # Non-uniform check
    m_nonuni = re.search(r'internalField\s+nonuniform\s+List<scalar>\s*\d+\s*\(\s*([\s\S]+?)\s*\)\s*;', content)
    if m_nonuni:
        return np.array([float(x) for x in m_nonuni.group(1).split()])

    print(f"  WARNING: Could not parse internalField in {path}")
    return np.ones(n_cells)


# ─────────────────────────────────────────────────────────────────────────────
# Step 4: Write OpenFOAM volScalarField
# ─────────────────────────────────────────────────────────────────────────────

def write_of_field(path, time_str, name, n_cells, cell_rates, patch_names):
    """
    Writes a volScalarField with zero everywhere except reacted cells.
    cell_rates: dict {cell_id: rate_value}
    """
    field = np.zeros(n_cells)
    for cid, val in cell_rates.items():
        field[cid] = val

    with open(path, 'w') as f:
        f.write(OF_HEADER.format(time=time_str, name=name))

        f.write(f"internalField   nonuniform List<scalar>\n{n_cells}\n(\n")
        for v in field:
            f.write(f"{v:.8e}\n")
        f.write(");\n\n")

        f.write("boundaryField\n{\n")
        for patch in patch_names:
            if patch in ("frontandback", "front", "back"):
                f.write(f"    {patch}\n    {{\n        type            empty;\n    }}\n")
            else:
                f.write(f"    {patch}\n    {{\n        type            zeroGradient;\n    }}\n")
        f.write("}\n\n")
        f.write("// ************************************************************************* //\n")


# ─────────────────────────────────────────────────────────────────────────────
# Step 5: Read patch names from boundary file
# ─────────────────────────────────────────────────────────────────────────────

def read_patch_names(case_dir):
    boundary_file = os.path.join(case_dir, "constant", "polyMesh", "boundary")
    with open(boundary_file) as f:
        content = f.read()
    return re.findall(r'^\s{4}(\w+)\s*$', content, re.MULTILINE)


# ─────────────────────────────────────────────────────────────────────────────
# Config file helpers
# ─────────────────────────────────────────────────────────────────────────────

DEFAULT_CONFIG = {
    "_comment": "Edit rate_names to match your reactions. Leave empty [] to auto-name them.",
    "rate_names": [],
    "rates_file": "",
    "tol_coord": 1e-6,
}

def write_config(path):
    with open(path, 'w') as f:
        json.dump(DEFAULT_CONFIG, f, indent=2)
    print(f"Example config written to {path}. Edit and re-run with --config {path}")
    sys.exit(0)


def load_config(path):
    with open(path) as f:
        return json.load(f)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Map BRNS rates (ratesAtFinish.dat) to OpenFOAM volScalarFields.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument("case_dir", nargs="?", default=".",
                        help="OpenFOAM case directory (default: .)")
    parser.add_argument("--rates-file", default=None,
                        help="Path to ratesAtFinish.dat")
    parser.add_argument("--rate-names", default=None,
                        help='Comma-separated rate names, e.g. "aerobic,anoxic,decay"')
    parser.add_argument("--tol-coord", type=float, default=1e-6,
                        help="Coordinate tolerance for block-size detection (default: 1e-6)")
    parser.add_argument("--coord-scale", type=float, default=1e-6,
                        help="Scale factor applied to file coords before matching mesh "
                             "(default: 1e-6, converts µm → m as written by BRNSReaction.H)")
    parser.add_argument("--config", default=None,
                        help="Path to JSON config file (overrides other flags)")
    parser.add_argument("--write-config", action="store_true",
                        help="Write an example config file and exit")
    args = parser.parse_args()

    case_dir = os.path.abspath(args.case_dir)

    if args.write_config:
        write_config(os.path.join(case_dir, "rateConfig.json"))

    # Load config file if given (it overrides CLI for rate_names / rates_file)
    cfg = {}
    if args.config:
        cfg = load_config(args.config)

    rates_file   = (cfg.get("rates_file") or args.rates_file
                    or os.path.join(case_dir, "ratesAtFinish.dat"))
    tol          = cfg.get("tol_coord",    args.tol_coord)
    coord_scale  = cfg.get("coord_scale",  args.coord_scale)

    print("=" * 60)
    print(f"Case directory : {case_dir}")
    print(f"Rates file     : {rates_file}")
    print("=" * 60)

    # ── 1. Read dat file ───────────────────────────────────────────────────
    coords, rates_data, n_rates = read_rates_file(rates_file)
    if coord_scale != 1.0:
        coords = coords * coord_scale
        print(f"  Coords scaled by {coord_scale} (file units → mesh units)")

    # ── 2. Rate names ──────────────────────────────────────────────────────
    names_str = cfg.get("rate_names") or []
    if args.rate_names:
        names_str = [n.strip() for n in args.rate_names.split(",")]
    if not names_str:
        names_str = [f"rate_R{i+1}" for i in range(n_rates)]
    if len(names_str) != n_rates:
        sys.exit(f"ERROR: {len(names_str)} rate names given but {n_rates} "
                 f"rate columns detected. Adjust --rate-names.")

    print(f"Rate fields    : {names_str}")

    # ── 3. Read mesh cell centres ──────────────────────────────────────────
    print("\nReading mesh cell centres ...")
    cell_centres = read_cell_centres(case_dir)
    n_cells = len(cell_centres)
    print(f"  Total cells : {n_cells}")

    # ── 4. Align blocks to time directories ───────────────────────────────
    time_dirs = sorted(
        [d for d in os.listdir(case_dir)
         if re.match(r'^\d+\.?\d*$', d) and d != "0"],
        key=float
    )
    print(f"\nTime directories (excl. t=0) : {len(time_dirs)}")

    # ── 5. Choose mapping mode and write fields ───────────────────────────
    n_total = len(coords)
    patch_name = cfg.get("react_patch", "solidwalls")
    patch_names = read_patch_names(case_dir)
    
    is_degenerate = coordinates_are_degenerate(coords, tol=tol)

    print(f"\nWriting timestep(s) × {n_rates} rate field(s) ...")
    if is_degenerate:
        print("  [MODE] Coordinates degenerate -> using DYNAMIC PATCH-ORDERING.")
        print(f"  Reading faceCells of patch '{patch_name}' ...")
        face_cells, _ = get_reacting_cells(case_dir, patch_name)
        
        # We must align rows dynamically using alpha.water
        row_offset = 0
        written_count = 0
        
        for t_idx, time_str in enumerate(time_dirs):
            time_dir = os.path.join(case_dir, time_str)
            
            # Find which cells satisfy the brns trigger (sat > 0.5)
            alpha_field = read_alpha_field(case_dir, time_str, "alpha.water", n_cells)
            active_cells = [c for c in face_cells if min(max(alpha_field[c], 0.0), 1.0) > 0.5]
            block_size = len(active_cells)
            
            if row_offset + block_size > n_total:
                print(f"  Ran out of dat file rows at t={time_str}. Stop.")
                break
                
            block_rates = rates_data[row_offset : row_offset + block_size]
            row_offset += block_size
            written_count += 1
            
            for r_idx, rname in enumerate(names_str):
                cell_vals = {}
                for face_i, cid in enumerate(active_cells):
                    val = float(block_rates[face_i, r_idx])
                    if cid not in cell_vals or val > cell_vals[cid]:
                        cell_vals[cid] = val
                
                out_path = os.path.join(time_dir, rname)
                write_of_field(out_path, time_str, rname, n_cells, cell_vals, patch_names)
                
            if (written_count) % 50 == 0 or t_idx == len(time_dirs) - 1:
                print(f"  [{written_count:>4}/{len(time_dirs)}]  t = {time_str}  ({block_size} cells reacted)")

    else:
        print("  [MODE] Coordinates valid -> using COORDINATE-BASED mapping.")
        lookup       = build_coord_index(cell_centres)
        block_size   = detect_block_size_from_coords(coords, tol=tol)
        all_cell_ids = lookup(coords)
        n_blocks     = n_total // block_size
        n_write      = min(n_blocks, len(time_dirs))
        
        for t_idx in range(n_write):
            time_str = time_dirs[t_idx]
            time_dir = os.path.join(case_dir, time_str)

            start = t_idx * block_size
            end   = start + block_size
            block_cells = all_cell_ids[start:end]
            block_rates = rates_data[start:end]

            for r_idx, rname in enumerate(names_str):
                cell_vals = {}
                for face_i, cid in enumerate(block_cells):
                    val = float(block_rates[face_i, r_idx])
                    if cid not in cell_vals or val > cell_vals[cid]:
                        cell_vals[cid] = val

                out_path = os.path.join(time_dir, rname)
                write_of_field(out_path, time_str, rname, n_cells, cell_vals, patch_names)

            if (t_idx + 1) % 50 == 0 or t_idx == n_write - 1:
                print(f"  [{t_idx+1:>4}/{n_write}]  t = {time_str}")

    print("\nDone. Open your .foam file in ParaView — rate fields appear automatically.")


if __name__ == "__main__":
    main()
