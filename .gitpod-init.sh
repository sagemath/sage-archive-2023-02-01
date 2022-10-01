#!/usr/bin/env bash

# Exit on error
set -e

# Create conda environment
./bootstrap-conda
mamba env create --file src/environment-dev.yml --prefix venv
conda config --append envs_dirs $(pwd)
conda activate $(pwd)/venv

# Build sage
./bootstrap
./configure --enable-build-as-root --with-python=$CONDA_PREFIX/bin/python --prefix=$CONDA_PREFIX
pip install --no-build-isolation -v -v -e ./pkgs/sage-conf ./pkgs/sage-setup
pip install --no-build-isolation -v -v -e ./src
