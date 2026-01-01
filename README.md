# QTY μOR - Fentanyl Docking Pipeline

## Project Overview
This pipeline tests whether a QTY-modified (water-soluble) μ-opioid receptor (OPRM1) variant can still bind fentanyl in aqueous solution.

## Pipeline Steps

### 1. Structure Preparation (`01_prepare_structure/`)
- Download PDB 8EF5 (fentanyl-bound human μOR–Gi complex)
- Extract μOR + fentanyl, remove Gi/scFv/cholesterols
- Keep pocket waters (optional)

### 2. QTY Mutation Introduction (`02_introduce_mutations/`)
- Introduce QTY mutations onto μOR coordinates
- Side-chain repacking

### 3. Structure Relaxation (`03_relaxation/`)
- Local relax/minimization using OpenMM
- Generate ensemble of ~10-30 conformations with different seeds

### 4. Docking (`04_docking/`)
- Dock fentanyl into ensemble using smina/AutoDock Vina
- Define docking box centered on 8EF5 fentanyl pocket
- Analyze pose clustering and key interactions (Asp3.32 salt bridge)

### 5. MD Simulation (`05_md_simulation/`)
- Short aqueous MD (20-100 ns) on top 3-5 docked complexes
- Track: ligand RMSD, Asp3.32 salt bridge persistence, pocket collapse

### 6. Analysis (`06_analysis/`)
- Compare QTY variant to WT control
- Rank by pose stability and MM/GBSA binding estimates

## Requirements
```bash
# Create conda environment
conda create -n qty_docking python=3.10
conda activate qty_docking

# Core packages
conda install -c conda-forge openmm mdtraj biopython numpy scipy matplotlib pandas

# For docking
conda install -c conda-forge vina
# Or install smina separately

# For structure manipulation
pip install prody pdbfixer

# Optional: PyMOL for visualization
conda install -c conda-forge pymol-open-source
```

## Quick Start
```bash
# 1. Set up environment
conda activate qty_docking

# 2. Edit config.yaml with your QTY mutations

# 3. Run the pipeline
python run_pipeline.py
```

## Key Validation Checkpoints

1. **Pose Clustering**: Ligand should dock into same pocket with similar pose across ensemble
2. **Asp3.32 Salt Bridge**: Critical interaction for opioid binding - must be preserved
3. **Ligand RMSD in MD**: Should stay < 3Å from starting pose
4. **Pocket Integrity**: Helix bundle should not collapse in water

## Controls
- **WT μOR**: Run same pipeline starting from 8EF5 as control
- **Negative Control**: Asp3.32→Ala mutation should abolish binding
