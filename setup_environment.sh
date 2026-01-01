#!/bin/bash
# QTY μOR - Fentanyl Docking Pipeline
# Environment Setup Script

set -e

echo "=================================================="
echo "Setting up QTY Docking Pipeline Environment"
echo "=================================================="

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: conda not found. Please install Miniconda or Anaconda first."
    echo "Download from: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Create conda environment
ENV_NAME="qty_docking"

echo ""
echo "Creating conda environment: $ENV_NAME"
echo ""

conda create -n $ENV_NAME python=3.10 -y

# Activate environment
echo "Activating environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

# Install core packages from conda-forge
echo ""
echo "Installing core packages..."
echo ""

conda install -c conda-forge \
    numpy \
    scipy \
    pandas \
    matplotlib \
    biopython \
    openmm \
    pdbfixer \
    mdtraj \
    prody \
    pyyaml \
    -y

# Install docking software
echo ""
echo "Installing docking software..."
echo ""

conda install -c conda-forge smina -y || {
    echo "Warning: Could not install smina from conda-forge"
    echo "You may need to install it manually from:"
    echo "https://sourceforge.net/projects/smina/"
}

# Install Open Babel for file conversion
conda install -c conda-forge openbabel -y || {
    echo "Warning: Could not install openbabel"
}

# Optional: Install PyMOL for visualization
echo ""
echo "Installing PyMOL (optional, for visualization)..."
echo ""

conda install -c conda-forge pymol-open-source -y || {
    echo "Warning: Could not install PyMOL"
    echo "You can install it later with: conda install -c conda-forge pymol-open-source"
}

# Create directory structure
echo ""
echo "Creating directory structure..."
echo ""

mkdir -p 01_prepare_structure/output
mkdir -p 02_introduce_mutations/output
mkdir -p 03_relaxation/output
mkdir -p 04_docking/output
mkdir -p 05_md_simulation/output
mkdir -p 06_analysis/output
mkdir -p results

echo ""
echo "=================================================="
echo "Environment setup complete!"
echo "=================================================="
echo ""
echo "To activate the environment, run:"
echo "    conda activate $ENV_NAME"
echo ""
echo "To run the pipeline, use:"
echo "    python run_pipeline.py --mutations L65Q I67T V70T ..."
echo ""
echo "Or with a config file:"
echo "    python run_pipeline.py --config config.yaml"
echo ""
echo "IMPORTANT NOTES:"
echo ""
echo "1. Fentanyl parameterization:"
echo "   The ligand (fentanyl) needs to be parameterized for MD."
echo "   If you need to do this, install AmberTools:"
echo "       conda install -c conda-forge ambertools"
echo "   Then use antechamber to generate parameters."
echo ""
echo "2. For best mutation modeling, consider using PyMOL's"
echo "   mutagenesis wizard or Rosetta for proper side-chain repacking."
echo ""
echo "3. Check config.yaml and update the QTY mutations section"
echo "   with your actual mutations before running."
echo ""
