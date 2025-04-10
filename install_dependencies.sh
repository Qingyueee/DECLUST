#!/bin/bash

echo "===== Installing DECLUST dependencies ====="

# Step 0: Check for conda
if ! command -v conda &> /dev/null
then
    echo "conda not found. Please install Miniconda or Anaconda first."
    exit 1
fi

conda config --add channels defaults
conda config --add channels conda-forge

# Step 1: Install R
echo "Installing R 4.3 via conda-forge..."
conda install -y -c conda-forge r-base=4.3

# Step 2: Pre-install system-level R packages
echo "Installing R system packages via conda..."
conda install -y -c conda-forge \
    r-matrix \
    r-rcpp \
    r-lattice \
    r-reticulate \
    r-bh \
    r-glue \
    r-rlang \
    r-ellipsis

# Step 3: Install scanpy using pip
echo "Installing scanpy via pip..."
pip install scanpy

# Step 4: Install rpy2
echo "Installing rpy2 via conda-forge..."
conda install -y -c conda-forge rpy2

# Step 5: Install R packages via script
echo "Installing R packages from install_R_dependencies.R..."
Rscript install_R_dependencies.R

echo "===== All dependencies installed successfully ====="
