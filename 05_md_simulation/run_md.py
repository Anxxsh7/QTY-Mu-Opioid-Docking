#!/usr/bin/env python3
"""
Step 5: MD Simulation Setup and Analysis
========================================
Sets up and runs explicit-solvent MD simulations to test
fentanyl binding stability in water (the "soluble" test).

Author: QTY Docking Pipeline
"""

import os
import sys
import argparse
import json
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import warnings

import numpy as np

try:
    import openmm as mm
    from openmm import app, unit
    from openmm.app import (PDBFile, ForceField, Modeller, Simulation,
                            PME, HBonds, StateDataReporter, DCDReporter,
                            PDBReporter)
except ImportError:
    print("Error: OpenMM not installed")
    sys.exit(1)

try:
    from pdbfixer import PDBFixer
    HAS_PDBFIXER = True
except ImportError:
    HAS_PDBFIXER = False

try:
    import mdtraj as md
    HAS_MDTRAJ = True
except ImportError:
    HAS_MDTRAJ = False
    print("Warning: MDTraj not installed - analysis will be limited")


class LigandParameterizer:
    """Handle ligand parameterization for MD."""
    
    @staticmethod
    def check_openff():
        """Check if OpenFF toolkit is available."""
        try:
            from openff.toolkit import Molecule
            return True
        except ImportError:
            return False
    
    @staticmethod
    def parameterize_ligand_gaff(ligand_pdb: Path, output_dir: Path) -> Dict:
        """Parameterize ligand using GAFF (via antechamber)."""
        
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # This would use antechamber/acpype
        # For now, return placeholder
        print("NOTE: Ligand parameterization requires antechamber or OpenFF")
        print("Please parameterize fentanyl separately using:")
        print("  antechamber -i ligand.pdb -fi pdb -o ligand.mol2 -fo mol2 -c bcc -s 2")
        print("  parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod")
        
        return {
            'mol2': output_dir / 'ligand.mol2',
            'frcmod': output_dir / 'ligand.frcmod'
        }


def setup_md_system(complex_pdb: Path, output_dir: Path,
                    box_padding: float = 12.0,
                    ionic_strength: float = 0.15) -> Tuple[any, any, Path]:
    """Set up solvated MD system."""
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\nSetting up MD system...")
    print(f"  Box padding: {box_padding} Å")
    print(f"  Ionic strength: {ionic_strength} M")
    
    # Load and fix structure
    if HAS_PDBFIXER:
        fixer = PDBFixer(filename=str(complex_pdb))
        fixer.findMissingResidues()
        
        # Remove terminal missing residues
        keys_to_remove = [k for k in fixer.missingResidues 
                          if k[1] == 0 or k[1] > 500]
        for key in keys_to_remove:
            del fixer.missingResidues[key]
        
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.4)
        
        topology = fixer.topology
        positions = fixer.positions
    else:
        pdb = PDBFile(str(complex_pdb))
        topology = pdb.topology
        positions = pdb.positions
    
    # Load forcefield
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    
    # Create modeller
    modeller = Modeller(topology, positions)
    
    # Add solvent
    modeller.addSolvent(
        forcefield,
        padding=box_padding * unit.angstrom,
        ionicStrength=ionic_strength * unit.molar,
        neutralize=True
    )
    
    print(f"  Total atoms after solvation: {modeller.topology.getNumAtoms()}")
    
    # Save solvated system
    solvated_pdb = output_dir / "solvated.pdb"
    with open(solvated_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    print(f"  Solvated system: {solvated_pdb}")
    
    # Create system
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=10 * unit.angstrom,
        constraints=HBonds
    )
    
    return system, modeller, solvated_pdb


def get_atom_indices(topology, selection: str) -> List[int]:
    """Get atom indices for a selection."""
    
    indices = []
    
    if selection == 'backbone':
        backbone_names = {'CA', 'C', 'N', 'O'}
        for atom in topology.atoms():
            if atom.name in backbone_names:
                indices.append(atom.index)
    
    elif selection == 'protein':
        for atom in topology.atoms():
            if atom.residue.name in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                                      'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                                      'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                                      'SER', 'THR', 'TRP', 'TYR', 'VAL']:
                indices.append(atom.index)
    
    elif selection.startswith('resnum:'):
        resnum = int(selection.split(':')[1])
        for atom in topology.atoms():
            if atom.residue.index == resnum:
                indices.append(atom.index)
    
    return indices


def add_position_restraints(system, positions, atom_indices: List[int],
                            force_constant: float = 100.0):
    """Add position restraints to atoms."""
    
    k = force_constant * unit.kilojoules_per_mole / unit.nanometer**2
    
    force = mm.CustomExternalForce(
        "0.5*k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
    )
    force.addGlobalParameter("k", k)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    
    for idx in atom_indices:
        pos = positions[idx]
        force.addParticle(idx, [pos[0], pos[1], pos[2]])
    
    system.addForce(force)
    return force


def run_equilibration(system, modeller, output_dir: Path,
                      nvt_steps: int = 50000,
                      npt_steps: int = 250000,
                      temperature: float = 300.0,
                      restraint_force: float = 100.0) -> Path:
    """Run equilibration (NVT then NPT)."""
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n--- Equilibration ---")
    
    # Add restraints to backbone
    backbone_indices = get_atom_indices(modeller.topology, 'backbone')
    restraint = add_position_restraints(
        system, modeller.positions, backbone_indices, restraint_force
    )
    print(f"  Restraining {len(backbone_indices)} backbone atoms")
    
    # NVT equilibration
    print(f"\n  NVT equilibration ({nvt_steps} steps)...")
    
    integrator = mm.LangevinMiddleIntegrator(
        temperature * unit.kelvin,
        1 / unit.picosecond,
        2 * unit.femtoseconds
    )
    
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    
    # Minimize first
    print("    Minimizing...")
    simulation.minimizeEnergy(maxIterations=1000)
    
    # Set velocities
    simulation.context.setVelocitiesToTemperature(temperature * unit.kelvin)
    
    # NVT
    simulation.reporters.append(
        StateDataReporter(
            str(output_dir / 'nvt.log'), 5000,
            step=True, temperature=True, potentialEnergy=True
        )
    )
    simulation.step(nvt_steps)
    
    # Save NVT checkpoint
    nvt_state = simulation.context.getState(getPositions=True, getVelocities=True)
    
    # NPT equilibration
    print(f"\n  NPT equilibration ({npt_steps} steps)...")
    
    # Add barostat
    barostat = mm.MonteCarloBarostat(1 * unit.bar, temperature * unit.kelvin)
    system.addForce(barostat)
    
    # Need new simulation with barostat
    integrator2 = mm.LangevinMiddleIntegrator(
        temperature * unit.kelvin,
        1 / unit.picosecond,
        2 * unit.femtoseconds
    )
    
    simulation2 = Simulation(modeller.topology, system, integrator2)
    simulation2.context.setPositions(nvt_state.getPositions())
    simulation2.context.setVelocities(nvt_state.getVelocities())
    
    simulation2.reporters.append(
        StateDataReporter(
            str(output_dir / 'npt.log'), 5000,
            step=True, temperature=True, potentialEnergy=True,
            density=True, volume=True
        )
    )
    simulation2.step(npt_steps)
    
    # Save equilibrated structure
    eq_state = simulation2.context.getState(getPositions=True)
    eq_pdb = output_dir / "equilibrated.pdb"
    with open(eq_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, eq_state.getPositions(), f)
    
    print(f"  Equilibrated structure: {eq_pdb}")
    
    # Save checkpoint
    simulation2.saveCheckpoint(str(output_dir / 'equilibrated.chk'))
    
    return eq_pdb


def run_production_md(system, topology, positions, output_dir: Path,
                      duration_ns: float = 50.0,
                      temperature: float = 300.0,
                      save_interval_ps: float = 10.0,
                      checkpoint_file: Path = None,
                      run_id: int = 1) -> Dict:
    """Run production MD simulation."""
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    timestep = 2.0  # femtoseconds
    n_steps = int(duration_ns * 1e6 / timestep)
    save_freq = int(save_interval_ps * 1000 / timestep)
    
    print(f"\n--- Production MD (Run {run_id}) ---")
    print(f"  Duration: {duration_ns} ns")
    print(f"  Steps: {n_steps}")
    print(f"  Save interval: {save_interval_ps} ps")
    
    # Set up integrator
    integrator = mm.LangevinMiddleIntegrator(
        temperature * unit.kelvin,
        1 / unit.picosecond,
        timestep * unit.femtoseconds
    )
    
    simulation = Simulation(topology, system, integrator)
    
    if checkpoint_file and checkpoint_file.exists():
        simulation.loadCheckpoint(str(checkpoint_file))
    else:
        simulation.context.setPositions(positions)
        simulation.context.setVelocitiesToTemperature(temperature * unit.kelvin)
    
    # Add reporters
    dcd_file = output_dir / f"trajectory_run{run_id}.dcd"
    log_file = output_dir / f"production_run{run_id}.log"
    
    simulation.reporters.append(DCDReporter(str(dcd_file), save_freq))
    simulation.reporters.append(
        StateDataReporter(
            str(log_file), save_freq,
            step=True, time=True, temperature=True,
            potentialEnergy=True, kineticEnergy=True,
            totalEnergy=True, speed=True
        )
    )
    
    # Run production
    print(f"  Running simulation...")
    simulation.step(n_steps)
    
    # Save final state
    final_state = simulation.context.getState(getPositions=True)
    final_pdb = output_dir / f"final_run{run_id}.pdb"
    with open(final_pdb, 'w') as f:
        PDBFile.writeFile(topology, final_state.getPositions(), f)
    
    print(f"  Trajectory: {dcd_file}")
    print(f"  Final structure: {final_pdb}")
    
    return {
        'trajectory': str(dcd_file),
        'log': str(log_file),
        'final_structure': str(final_pdb),
        'duration_ns': duration_ns
    }


def analyze_trajectory(topology_file: Path, trajectory_file: Path,
                       output_dir: Path,
                       ligand_resname: str = 'ZPE',
                       asp_resnum: int = 147) -> Dict:
    """Analyze MD trajectory for binding stability."""
    
    if not HAS_MDTRAJ:
        print("MDTraj not available - skipping trajectory analysis")
        return {}
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n--- Trajectory Analysis ---")
    
    # Load trajectory
    traj = md.load(str(trajectory_file), top=str(topology_file))
    print(f"  Frames: {traj.n_frames}")
    print(f"  Duration: {traj.time[-1]/1000:.2f} ns")
    
    # Find ligand atoms
    ligand_atoms = traj.topology.select(f"resname {ligand_resname}")
    if len(ligand_atoms) == 0:
        print(f"  Warning: Ligand {ligand_resname} not found!")
        return {}
    
    print(f"  Ligand atoms: {len(ligand_atoms)}")
    
    # Calculate ligand RMSD
    # First, align on protein backbone
    protein_ca = traj.topology.select("name CA")
    traj_aligned = traj.superpose(traj, atom_indices=protein_ca)
    
    # Calculate ligand RMSD relative to first frame
    ligand_xyz = traj_aligned.xyz[:, ligand_atoms, :]
    ligand_ref = ligand_xyz[0]
    
    rmsd = np.sqrt(np.mean(np.sum((ligand_xyz - ligand_ref)**2, axis=2), axis=1))
    rmsd_nm = rmsd  # Already in nm
    rmsd_angstrom = rmsd * 10
    
    print(f"  Ligand RMSD:")
    print(f"    Mean: {np.mean(rmsd_angstrom):.2f} Å")
    print(f"    Max: {np.max(rmsd_angstrom):.2f} Å")
    print(f"    Final: {rmsd_angstrom[-1]:.2f} Å")
    
    # Check Asp salt bridge
    # Find Asp OD1/OD2 and ligand nitrogen
    asp_oxygens = traj.topology.select(f"resid {asp_resnum} and (name OD1 or name OD2)")
    ligand_nitrogens = traj.topology.select(f"resname {ligand_resname} and element N")
    
    salt_bridge_distances = []
    if len(asp_oxygens) > 0 and len(ligand_nitrogens) > 0:
        for frame_idx in range(traj.n_frames):
            min_dist = float('inf')
            for o_idx in asp_oxygens:
                for n_idx in ligand_nitrogens:
                    dist = np.linalg.norm(
                        traj.xyz[frame_idx, o_idx] - traj.xyz[frame_idx, n_idx]
                    )
                    min_dist = min(min_dist, dist)
            salt_bridge_distances.append(min_dist * 10)  # nm to Å
        
        salt_bridge_distances = np.array(salt_bridge_distances)
        salt_bridge_intact = salt_bridge_distances < 4.0
        
        print(f"\n  Asp{asp_resnum} Salt Bridge:")
        print(f"    Mean distance: {np.mean(salt_bridge_distances):.2f} Å")
        print(f"    Intact fraction: {np.mean(salt_bridge_intact)*100:.1f}%")
    else:
        salt_bridge_distances = None
        print(f"\n  Could not analyze Asp{asp_resnum} salt bridge")
    
    # Save analysis results
    analysis = {
        'n_frames': traj.n_frames,
        'duration_ns': float(traj.time[-1] / 1000),
        'ligand_rmsd': {
            'mean_angstrom': float(np.mean(rmsd_angstrom)),
            'max_angstrom': float(np.max(rmsd_angstrom)),
            'final_angstrom': float(rmsd_angstrom[-1]),
            'values': rmsd_angstrom.tolist()
        }
    }
    
    if salt_bridge_distances is not None:
        analysis['asp_salt_bridge'] = {
            'mean_distance': float(np.mean(salt_bridge_distances)),
            'intact_fraction': float(np.mean(salt_bridge_intact)),
            'values': salt_bridge_distances.tolist()
        }
    
    # Save to file
    with open(output_dir / 'md_analysis.json', 'w') as f:
        json.dump(analysis, f, indent=2)
    
    # Save plots
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        fig, axes = plt.subplots(2, 1, figsize=(10, 8))
        
        time_ns = traj.time / 1000
        
        # RMSD plot
        axes[0].plot(time_ns, rmsd_angstrom)
        axes[0].axhline(y=3.0, color='r', linestyle='--', label='3 Å threshold')
        axes[0].set_xlabel('Time (ns)')
        axes[0].set_ylabel('Ligand RMSD (Å)')
        axes[0].set_title('Fentanyl RMSD')
        axes[0].legend()
        
        # Salt bridge plot
        if salt_bridge_distances is not None:
            axes[1].plot(time_ns, salt_bridge_distances)
            axes[1].axhline(y=4.0, color='r', linestyle='--', label='4 Å cutoff')
            axes[1].set_xlabel('Time (ns)')
            axes[1].set_ylabel('Distance (Å)')
            axes[1].set_title(f'Asp{asp_resnum} - Ligand N Distance')
            axes[1].legend()
        
        plt.tight_layout()
        plt.savefig(output_dir / 'md_analysis.png', dpi=150)
        plt.close()
        
        print(f"\n  Analysis plot saved: {output_dir / 'md_analysis.png'}")
        
    except ImportError:
        print("  Matplotlib not available - skipping plots")
    
    return analysis


def main():
    parser = argparse.ArgumentParser(
        description="Run MD simulation for fentanyl binding stability"
    )
    parser.add_argument('--complex-pdb', required=True,
                        help='Docked complex PDB file')
    parser.add_argument('--output-dir', default='05_md_simulation/output',
                        help='Output directory')
    parser.add_argument('--duration-ns', type=float, default=50.0,
                        help='Production MD duration in ns')
    parser.add_argument('--n-replicates', type=int, default=3,
                        help='Number of replicate simulations')
    parser.add_argument('--ligand-resname', default='ZPE',
                        help='Ligand residue name')
    parser.add_argument('--asp-resnum', type=int, default=147,
                        help='Asp3.32 residue number')
    parser.add_argument('--equilibrate-only', action='store_true',
                        help='Only run equilibration')
    parser.add_argument('--analyze-only',
                        help='Only analyze existing trajectory')
    
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Step 5: MD Simulation for Binding Stability")
    print("=" * 60)
    
    # Analysis only mode
    if args.analyze_only:
        traj_dir = Path(args.analyze_only)
        topology = traj_dir / "solvated.pdb"
        for traj_file in traj_dir.glob("trajectory_*.dcd"):
            analysis = analyze_trajectory(
                topology, traj_file, traj_dir / "analysis",
                args.ligand_resname, args.asp_resnum
            )
        return
    
    # Set up system
    complex_pdb = Path(args.complex_pdb)
    system, modeller, solvated_pdb = setup_md_system(
        complex_pdb, output_dir
    )
    
    # Run equilibration
    eq_pdb = run_equilibration(
        system, modeller, output_dir / "equilibration"
    )
    
    if args.equilibrate_only:
        print("\nEquilibration complete (--equilibrate-only flag set)")
        return
    
    # Run production MD replicates
    all_results = []
    
    for rep in range(1, args.n_replicates + 1):
        rep_dir = output_dir / f"replicate_{rep}"
        
        # Load equilibrated state
        eq_pdb_loaded = PDBFile(str(eq_pdb))
        
        result = run_production_md(
            system, modeller.topology, eq_pdb_loaded.positions,
            rep_dir, duration_ns=args.duration_ns,
            run_id=rep
        )
        
        # Analyze trajectory
        analysis = analyze_trajectory(
            solvated_pdb, Path(result['trajectory']),
            rep_dir / "analysis",
            args.ligand_resname, args.asp_resnum
        )
        
        result['analysis'] = analysis
        all_results.append(result)
    
    # Save summary
    summary = {
        'n_replicates': args.n_replicates,
        'duration_ns': args.duration_ns,
        'results': all_results
    }
    
    with open(output_dir / 'md_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Print summary
    print("\n" + "=" * 60)
    print("MD SIMULATION SUMMARY")
    print("=" * 60)
    
    rmsd_finals = []
    salt_bridge_fracs = []
    
    for i, r in enumerate(all_results):
        if 'analysis' in r and r['analysis']:
            rmsd = r['analysis'].get('ligand_rmsd', {}).get('final_angstrom')
            sb = r['analysis'].get('asp_salt_bridge', {}).get('intact_fraction')
            
            if rmsd:
                rmsd_finals.append(rmsd)
            if sb is not None:
                salt_bridge_fracs.append(sb)
            
            print(f"\nReplicate {i+1}:")
            print(f"  Final ligand RMSD: {rmsd:.2f} Å" if rmsd else "  RMSD: N/A")
            print(f"  Salt bridge intact: {sb*100:.1f}%" if sb else "  Salt bridge: N/A")
    
    if rmsd_finals:
        print(f"\n--- Overall ---")
        print(f"  Mean final RMSD: {np.mean(rmsd_finals):.2f} ± {np.std(rmsd_finals):.2f} Å")
        
        if np.mean(rmsd_finals) < 3.0:
            print("  ✓ Ligand appears STABLE in binding pocket")
        else:
            print("  ⚠ Ligand may be DRIFTING from binding pocket")
    
    if salt_bridge_fracs:
        print(f"  Mean salt bridge intact: {np.mean(salt_bridge_fracs)*100:.1f}%")
    
    print("=" * 60)


if __name__ == "__main__":
    main()
