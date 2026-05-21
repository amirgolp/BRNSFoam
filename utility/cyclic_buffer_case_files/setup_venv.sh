#!/bin/bash
#
# setup_venv.sh
#
# One-time setup for the Python utilities used in this case:
#   - clip_buffer_stl.py       (slice the source STL into a buffered form)
#
# Creates a venv at $VENV_DIR and installs:
#   numpy<2.0     pinned for trimesh+scikit-image compatibility
#   scipy         general scientific
#   scikit-image  (legacy from make_periodic_stl.py; harmless)
#   numpy-stl     STL reader/writer
#   trimesh       mesh ops
#   rtree         spatial-index acceleration for trimesh
#   networkx      enclosure-tree traversal during slice-mesh-plane cap
#   shapely       polygon validity / containment during cap
#   mapbox-earcut polygon triangulation for cap
#
# Reuse this venv across all your channelStationary* cases; there's no
# per-case state in it.  Re-run this script if the venv is missing or
# if you want to refresh the dependencies.
#
# Usage:
#   ./setup_venv.sh                  # create at the default location
#   VENV_DIR=/some/path ./setup_venv.sh
#
set -e

VENV_DIR="${VENV_DIR:-$HOME/venvs/openfoam_py}"

echo "Target venv: $VENV_DIR"

PYTHON_BIN="${PYTHON_BIN:-python3}"
echo "Using interpreter: $($PYTHON_BIN --version 2>&1)  ($(command -v $PYTHON_BIN))"

if [ -d "$VENV_DIR" ]; then
    echo "Venv already exists; reusing."
else
    echo "Creating venv ..."
    mkdir -p "$(dirname "$VENV_DIR")"
    "$PYTHON_BIN" -m venv "$VENV_DIR"
fi

# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"

# Cluster module systems (SciPy-bundle, etc.) prepend their site-packages
# to PYTHONPATH and shadow the venv at import time.  Drop it for the rest
# of this script.
if [ -n "${PYTHONPATH-}" ]; then
    echo "Clearing PYTHONPATH (was: $PYTHONPATH) so the venv isn't shadowed."
    unset PYTHONPATH
fi

python -m pip install --upgrade pip wheel

# trimesh's heavy "all" extras pull in pyembree etc. -- we only need the
# core slicing / capping code which lives in the base install plus three
# small triangulation deps (networkx + shapely + mapbox-earcut).
python -m pip install --upgrade \
    "numpy<2.0" \
    "scipy" \
    "scikit-image" \
    "numpy-stl" \
    "trimesh" \
    "rtree" \
    "networkx" \
    "shapely" \
    "mapbox-earcut"

echo
echo "Installed packages relevant to clip_buffer_stl.py:"
python -m pip list 2>/dev/null | grep -Ei '^(numpy|scipy|scikit-image|numpy-stl|trimesh|rtree|networkx|shapely|mapbox-earcut)\b' || true

echo
echo "Done.  To use this venv from createMesh*.sh (or any shell):"
echo "    source $VENV_DIR/bin/activate"
