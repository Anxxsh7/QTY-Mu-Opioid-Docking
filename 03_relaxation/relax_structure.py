#!/usr/bin/env python3
"""
Step 3: Relaxation and Ensemble Generation
==========================================
Performs local relaxation/minimization on the mutated structure
and generates an ensemble of conformations for docking.

Uses OpenMM for energy minimization and short MD relaxation.

Author: QTY Docking Pipeline
"""

import os
import sys
import argparse
import json
from pathlib import Path
from typing import List, Tuple, Optional
import warnings

import numpy as np

try:
    import openmm as mm
    from openmm import app, unit
    from openmm.app import PDBFile, ForceField, Modeller, Simulation
    from openmm.app import PME, HBonds, NoCutoff
except ImportError:
    print("Error: OpenMM not installed. Run: conda install -c conda-forge openmm")
    sys.exit(1)

try:
    from pdbfixer import PDBFixer
    HAS_PDBFIXER = True
except ImportError:
    HAS_PDBFIXER = False


def prepare_structure_for_openmm(input_pdb: Path, output_pdb: Path,
                                  add_hydrogens: bool = True,
                                  ph: float = 7.4) -> Path:
    """Prepare structure with PDBFixer."""
    
    if not HAS_PDBFIXER:
        print("Warning: PDBFixer not available")
        return input_pdb
    
    print(f"Preparing structure with PDBFixer...")
    
    fixer = PDBFixer(filename=str(input_pdb))
    
    fixer.findMissingResidues()
    
    # Remove terminal missing residues
    keys_to_remove = [k for k in fixer.missingResidues 
                      if k[1] == 0 or k[1] > 500]
    for key in keys_to_remove:
        del fixer.missingResidues[key]
    
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    
    if add_hydrogens:
        fixer.addMissingHydrogens(ph)
    
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    print(f"Prepared structure: {output_pdb}")
    return output_pdb


def get_backbone_atoms(topology) -> List[int]:
    """Get indices of backbone atoms for restraints."""
    backbone_names = {'CA', 'C', 'N', 'O'}
    backbone_indices = []
    
    for atom in topology.atoms():
        if atom.name in backbone_names:
            backbone_indices.append(atom.index)
    
    return backbone_indices


def add_backbone_restraints(system, positions, backbone_indices, 
                            force_constant: float = 100.0):
    """Add harmonic restraints to backbone atoms."""
    
    # force_constant in kJ/mol/nm^2
    k = force_constant * unit.kilojoules_per_mole / unit.nanometer**2
    
    restraint_force = mm.CustomExternalForce(
        "0.5*k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
    )
    restraint_force.addGlobalParameter("k", k)
    restraint_force.addPerParticleParameter("x0")
    restraint_force.addPerParticleParameter("y0")
    restraint_force.addPerParticleParameter("z0")
    
    for idx in backbone_indices:
        pos = positions[idx]
        restraint_force.addParticle(idx, [pos[0], pos[1], pos[2]])
    
    system.addForce(restraint_force)
    print(f"Added restraints to {len(backbone_indices)} backbone atoms")


def minimize_structure(pdb_file: Path, output_pdb: Path,
                       max_iterations: int = 1000,
                       tolerance: float = 10.0,
                       restrain_backbone: bool = True,
                       restraint_force: float = 100.0) -> Path:
    """Perform energy minimization on structure."""
    
    print(f"\nPerforming energy minimization...")
    print(f"  Max iterations: {max_iterations}")
    print(f"  Tolerance: {tolerance} kJ/mol/nm")
    
    pdb = PDBFile(str(pdb_file))
    
    # Use implicit solvent for speed
    forcefield = ForceField('amber14-all.xml', 'implicit/gbn2.xml')
    
    # Create modeller and add hydrogens if needed
    modeller = Modeller(pdb.topology, pdb.positions)
    
    # Create system
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=HBonds
    )
    
    # Add backbone restraints if requested
    if restrain_backbone:
        backbone_indices = get_backbone_atoms(modeller.topology)
        add_backbone_restraints(system, modeller.positions, 
                               backbone_indices, restraint_force)
    
    # Create simulation
    integrator = mm.LangevinMiddleIntegrator(
        300*unit.kelvin, 
        1/unit.picosecond, 
        0.002*unit.picoseconds
    )
    
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    
    # Get initial energy
    state = simulation.context.getState(getEnergy=True)
    initial_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"  Initial energy: {initial_energy:.2f} kJ/mol")
    
    # Minimize
    simulation.minimizeEnergy(
        maxIterations=max_iterations,
        tolerance=tolerance*unit.kilojoules_per_mole/unit.nanometer
    )
    
    # Get final energy
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    final_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"  Final energy: {final_energy:.2f} kJ/mol")
    print(f"  Energy change: {final_energy - initial_energy:.2f} kJ/mol")
    
    # Save minimized structure
    positions = state.getPositions()
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, positions, f)
    
    print(f"  Minimized structure: {output_pdb}")
    return output_pdb


def run_relaxation_md(pdb_file: Path, output_pdb: Path,
                      n_steps: int = 25000,
                      temperature: float = 300.0,
                      restrain_backbone: bool = True,
                      restraint_force: float = 50.0,
                      seed: int = None) -> Path:
    """Run short MD for relaxation/equilibration."""
    
    print(f"\nRunning relaxation MD...")
    print(f"  Steps: {n_steps} ({n_steps * 0.002:.1f} ps)")
    print(f"  Temperature: {temperature} K")
    if seed:
        print(f"  Random seed: {seed}")
    
    pdb = PDBFile(str(pdb_file))
    
    forcefield = ForceField('amber14-all.xml', 'implicit/gbn2.xml')
    modeller = Modeller(pdb.topology, pdb.positions)
    
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=HBonds
    )
    
    if restrain_backbone:
        backbone_indices = get_backbone_atoms(modeller.topology)
        add_backbone_restraints(system, modeller.positions, 
                               backbone_indices, restraint_force)
    
    integrator = mm.LangevinMiddleIntegrator(
        temperature*unit.kelvin, 
        1/unit.picosecond, 
        0.002*unit.picoseconds
    )
    
    if seed is not None:
        integrator.setRandomNumberSeed(seed)
    
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    simulation.context.setVelocitiesToTemperature(temperature*unit.kelvin)
    
    # Run MD
    simulation.step(n_steps)
    
    # Save final structure
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions()
    
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, positions, f)
    
    print(f"  Relaxed structure: {output_pdb}")
    return output_pdb


def generate_ensemble(pdb_file: Path, output_dir: Path,
                      n_structures: int = 20,
                      relaxation_steps: int = 25000,
                      temperature: float = 300.0,
                      seeds: List[int] = None) -> List[Path]:
    """Generate ensemble of structures by running relaxation with different seeds."""
    
    print(f"\n{'='*60}")
    print(f"Generating ensemble of {n_structures} structures")
    print(f"{'='*60}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if seeds is None:
        seeds = list(range(42, 42 + n_structures))
    
    ensemble_files = []
    
    for i, seed in enumerate(seeds[:n_structures]):
        print(f"\n--- Structure {i+1}/{n_structures} (seed={seed}) ---")
        
        output_pdb = output_dir / f"ensemble_{i+1:03d}.pdb"
        
        run_relaxation_md(
            pdb_file, output_pdb,
            n_steps=relaxation_steps,
            temperature=temperature,
            restrain_backbone=True,
            restraint_force=50.0,
            seed=seed
        )
        
        ensemble_files.append(output_pdb)
    
    # Save ensemble info
    ensemble_info = {
        'n_structures': n_structures,
        'relaxation_steps': relaxation_steps,
        'temperature': temperature,
        'seeds': seeds[:n_structures],
        'files': [str(f) for f in ensemble_files]
    }
    
    with open(output_dir / 'ensemble_info.json', 'w') as f:
        json.dump(ensemble_info, f, indent=2)
    
    print(f"\n{'='*60}")
    print(f"Ensemble generation complete!")
    print(f"Structures saved in: {output_dir}")
    print(f"{'='*60}")
    
    return ensemble_files


def main():
    parser = argparse.ArgumentParser(
        description="Relax mutated structure and generate ensemble"
    )
    parser.add_argument('--input-pdb', required=True,
                        help='Input PDB file (mutated structure)')
    parser.add_argument('--output-dir', default='03_relaxation/output',
                        help='Output directory')
    parser.add_argument('--minimize-only', action='store_true',
                        help='Only perform minimization, no MD')
    parser.add_argument('--n-structures', type=int, default=20,
                        help='Number of ensemble structures to generate')
    parser.add_argument('--relaxation-steps', type=int, default=25000,
                        help='MD steps for each relaxation (default: 25000 = 50ps)')
    parser.add_argument('--temperature', type=float, default=300.0,
                        help='Temperature in Kelvin')
    parser.add_argument('--no-restraints', action='store_true',
                        help='Do not restrain backbone during relaxation')
    
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Step 3: Structure Relaxation and Ensemble Generation")
    print("=" * 60)
    
    input_pdb = Path(args.input_pdb)
    
    # Prepare structure
    if HAS_PDBFIXER:
        prepared_pdb = output_dir / "prepared.pdb"
        prepare_structure_for_openmm(input_pdb, prepared_pdb)
        working_pdb = prepared_pdb
    else:
        working_pdb = input_pdb
    
    # Minimize
    minimized_pdb = output_dir / "minimized.pdb"
    minimize_structure(
        working_pdb, minimized_pdb,
        max_iterations=1000,
        restrain_backbone=not args.no_restraints
    )
    
    if args.minimize_only:
        print("\nMinimization complete (--minimize-only flag set)")
        return
    
    # Generate ensemble
    ensemble_dir = output_dir / "ensemble"
    ensemble_files = generate_ensemble(
        minimized_pdb, ensemble_dir,
        n_structures=args.n_structures,
        relaxation_steps=args.relaxation_steps,
        temperature=args.temperature
    )
    
    print(f"\nReady for docking with {len(ensemble_files)} structures")


if __name__ == "__main__":
    main()
    