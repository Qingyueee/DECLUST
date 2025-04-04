#!/bin/bash

echo "===== Installing DECLUST dependencies ====="

# Step 0: Check for conda
if ! command -v conda &> /dev/null
then
    echo "conda not found. Please install Miniconda or Anaconda first."
    exit 1
fi

# Step 1: Install R via conda-forge
echo "Installing R 4.3 via conda-forge..."
conda install -y -c conda-forge r-base=4.3

# Step 2: Install scanpy using pip
echo "Installing scanpy via pip..."
pip install scanpy

# Step 3: Install rpy2 using conda-forge
echo "Installing rpy2 via conda-forge..."
conda install -y -c conda-forge rpy2

# Step 4: Install R package 'dplyr' via rpy2
echo "Installing R package 'dplyr' via rpy2..."
python - <<EOF
from rpy2.robjects.packages import importr
utils = importr('utils')
utils.chooseCRANmirror(ind=1)  # choose the first CRAN mirror
utils.install_packages('dplyr')
EOF

echo "===== All dependencies installed successfully ====="
