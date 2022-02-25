#!/bin/sh
set -e -u

blockMesh
touch fluid-openfoam.foam

./run-openfoam.sh "$@"
. ./openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs
