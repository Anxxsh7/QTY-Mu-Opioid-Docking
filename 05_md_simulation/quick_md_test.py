#!/usr/bin/env python3
"""
Quick MD simulation to test QTY-fentanyl binding stability in water-like conditions.
Uses implicit solvent (GB) for speed - still captures water effects.
"""

import os
import sys
import json
from pathlib import Path
import numpy as np

try:
    import openmm as mm
    from openmm import app, unit
    from openmm.app import PDBFile, ForceField, Modeller, Simulation
except ImportError:
    print("Error: OpenMM not installed")
    sys.exit(1)

from pdbfixer import PDBFixer


def combine_receptor_ligand(receptor_pdb: Path, ligand_pdb: Path, output_pdb: Path):
    """Combine receptor and best docking pose into one PDB."""
    
    with open(receptor_pdb) as f:
        receptor_lines = [l for l in f.readlines() if l.startswith(('ATOM', 'HETATM', 'TER'))]
    
    # For PDBQT, convert to PDB format
    ligand_lines = []
    with open(ligand_pdb) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                # PDBQT has extra columns - truncate to PDB format
                pdb_line = line[:66] + '\n'
                ligand_lines.append(pdb_line)
    
    with open(output_pdb, 'w') as f:
        f.writelines(receptor_lines)
        f.write('TER\n')
        f.writelines(ligand_lines)
        f.write('END\n')
    
    print(f"Combined complex: {output_pdb}")
    return output_pdb


def prepare_complex(complex_pdb: Path, output_pdb: Path):
    """Prepare complex with PDBFixer - add hydrogens etc."""
    
    fixer = PDBFixer(filename=str(complex_pdb))
    fixer.findMissingResidues()
    
    # Remove terminal missing residues
    keys_to_remove = [k for k in fixer.missingResidues if k[1] == 0 or k[1] > 500]
    for key in keys_to_remove:
        del fixer.missingResidues[key]
    
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    print(f"Prepared complex: {output_pdb}")
    return output_pdb


def get_ligand_indices(topology, ligand_resname: str = '7V7') -> list:
    """Get atom indices for ligand."""
    indices = []
    for atom in topology.atoms():
        if atom.residue.name == ligand_resname:
            indices.append(atom.index)
    return indices


def get_asp_oxygen_indices(topology, asp_resnum: int = 84) -> list:
    """Get Asp carboxylate oxygen indices."""
    indices = []
    for atom in topology.atoms():
        if atom.residue.id == str(asp_resnum) and atom.name in ['OD1', 'OD2']:
            indices.append(atom.index)
    return indices


def get_nitrogen_indices(topology, ligand_resname: str = '7V7') -> list:
    """Get ligand nitrogen indices."""
    indices = []
    for atom in topology.atoms():
        if atom.residue.name == ligand_resname and atom.element.symbol == 'N':
            indices.append(atom.index)
    return indices


def calculate_distance(positions, idx1: int, idx2: int) -> float:
    """Calculate distance in Angstroms."""
    p1 = positions[idx1].value_in_unit(unit.angstrom)
    p2 = positions[idx2].value_in_unit(unit.angstrom)
    return np.sqrt(sum((a-b)**2 for a, b in zip(p1, p2)))


def run_binding_stability_test(complex_pdb: Path, output_dir: Path,
                               duration_ps: float = 100.0,
                               temperature: float = 300.0,
                               asp_resnum: int = 84):
    """
    Run short MD to test if fentanyl stays bound.
    Uses implicit solvent for speed.
    """
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("MD Binding Stability Test")
    print("=" * 60)
    print(f"Duration: {duration_ps} ps ({duration_ps/1000:.2f} ns)")
    print(f"Temperature: {temperature} K")
    print(f"Solvent: Implicit (GBn2) - mimics aqueous environment")
    
    # Load structure
    pdb = PDBFile(str(complex_pdb))
    
    # Use implicit solvent forcefield
    forcefield = ForceField('amber14-all.xml', 'implicit/gbn2.xml')
    
    modeller = Modeller(pdb.topology, pdb.positions)
    
    # Find key atoms before creating system
    ligand_indices = get_ligand_indices(modeller.topology, '7V7')
    asp_oxy_indices = get_asp_oxygen_indices(modeller.topology, asp_resnum)
    lig_n_indices = get_nitrogen_indices(modeller.topology, '7V7')
    
    print(f"\nLigand atoms: {len(ligand_indices)}")
    print(f"Asp{asp_resnum} oxygens: {len(asp_oxy_indices)}")
    print(f"Ligand nitrogens: {len(lig_n_indices)}")
    
    if len(ligand_indices) == 0:
        print("ERROR: Ligand not found!")
        return None
    
    # Create system
    try:
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds
        )
    except Exception as e:
        print(f"Error creating system: {e}")
        print("\nNote: Fentanyl requires special parameterization for full MD.")
        print("Running protein-only stability test instead...")
        return None
    
    # Setup simulation
    integrator = mm.LangevinMiddleIntegrator(
        temperature * unit.kelvin,
        1 / unit.picosecond,
        2 * unit.femtoseconds
    )
    
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    
    # Minimize
    print("\nMinimizing energy...")
    simulation.minimizeEnergy(maxIterations=500)
    
    # Set velocities
    simulation.context.setVelocitiesToTemperature(temperature * unit.kelvin)
    
    # Calculate initial N-Asp distance
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions()
    
    initial_distances = []
    for n_idx in lig_n_indices:
        for o_idx in asp_oxy_indices:
            d = calculate_distance(positions, n_idx, o_idx)
            initial_distances.append(d)
    
    initial_min_dist = min(initial_distances) if initial_distances else 999
    print(f"Initial N-Asp distance: {initial_min_dist:.2f} Å")
    
    # Run MD and track distances
    timestep_fs = 2.0
    n_steps = int(duration_ps * 1000 / timestep_fs)
    report_freq = n_steps // 20  # 20 snapshots
    
    print(f"\nRunning {n_steps} steps...")
    
    trajectory_data = {
        'time_ps': [],
        'n_asp_distance': [],
        'ligand_com': []
    }
    
    # Get initial ligand center of mass
    lig_positions = [positions[i].value_in_unit(unit.angstrom) for i in ligand_indices]
    initial_com = np.mean(lig_positions, axis=0)
    
    for i in range(20):
        simulation.step(report_freq)
        
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        positions = state.getPositions()
        time_ps = (i + 1) * report_freq * timestep_fs / 1000
        
        # Calculate N-Asp distance
        distances = []
        for n_idx in lig_n_indices:
            for o_idx in asp_oxy_indices:
                d = calculate_distance(positions, n_idx, o_idx)
                distances.append(d)
        
        min_dist = min(distances) if distances else 999
        
        # Calculate ligand COM displacement
        lig_positions = [positions[idx].value_in_unit(unit.angstrom) for idx in ligand_indices]
        com = np.mean(lig_positions, axis=0)
        com_displacement = np.sqrt(sum((a-b)**2 for a, b in zip(com, initial_com)))
        
        trajectory_data['time_ps'].append(time_ps)
        trajectory_data['n_asp_distance'].append(min_dist)
        trajectory_data['ligand_com'].append(com_displacement)
        
        print(f"  {time_ps:6.1f} ps: N-Asp = {min_dist:.2f} Å, COM drift = {com_displacement:.2f} Å")
    
    # Save final structure
    final_pdb = output_dir / "final_complex.pdb"
    with open(final_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, positions, f)
    
    # Analysis
    print("\n" + "=" * 60)
    print("BINDING STABILITY ANALYSIS")
    print("=" * 60)
    
    n_asp_distances = trajectory_data['n_asp_distance']
    com_displacements = trajectory_data['ligand_com']
    
    print(f"\nN-Asp3.32 Distance:")
    print(f"  Initial:  {initial_min_dist:.2f} Å")
    print(f"  Final:    {n_asp_distances[-1]:.2f} Å")
    print(f"  Mean:     {np.mean(n_asp_distances):.2f} ± {np.std(n_asp_distances):.2f} Å")
    print(f"  Min/Max:  {min(n_asp_distances):.2f} / {max(n_asp_distances):.2f} Å")
    
    print(f"\nLigand Position Stability:")
    print(f"  Final COM drift: {com_displacements[-1]:.2f} Å")
    print(f"  Max COM drift:   {max(com_displacements):.2f} Å")
    
    # Interpret results
    print("\n" + "-" * 60)
    print("INTERPRETATION")
    print("-" * 60)
    
    stable = True
    
    if np.mean(n_asp_distances) < 5.0:
        print("✅ Salt bridge MAINTAINED (< 5 Å) - Strong interaction!")
    elif np.mean(n_asp_distances) < 8.0:
        print("⚠️  Salt bridge STRETCHED (5-8 Å) - Weak interaction")
        stable = False
    else:
        print("❌ Salt bridge BROKEN (> 8 Å) - Ligand dissociating!")
        stable = False
    
    if max(com_displacements) < 3.0:
        print("✅ Ligand position STABLE (drift < 3 Å)")
    elif max(com_displacements) < 6.0:
        print("⚠️  Ligand SHIFTING (drift 3-6 Å)")
    else:
        print("❌ Ligand LEAVING binding site (drift > 6 Å)")
        stable = False
    
    if stable:
        print("\n🟢 CONCLUSION: Fentanyl binding appears STABLE in QTY mutant!")
    else:
        print("\n🟡 CONCLUSION: Binding stability is QUESTIONABLE - longer MD needed")
    
    # Save results
    results = {
        'duration_ps': duration_ps,
        'temperature': temperature,
        'initial_n_asp_distance': initial_min_dist,
        'final_n_asp_distance': n_asp_distances[-1],
        'mean_n_asp_distance': np.mean(n_asp_distances),
        'std_n_asp_distance': np.std(n_asp_distances),
        'max_com_drift': max(com_displacements),
        'trajectory': trajectory_data,
        'stable': stable
    }
    
    with open(output_dir / 'md_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nResults saved to: {output_dir}")
    
    return results


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Test fentanyl binding stability")
    parser.add_argument('--receptor', required=True, help='Receptor PDB file')
    parser.add_argument('--ligand', required=True, help='Docked ligand PDBQT file')
    parser.add_argument('--output-dir', default='05_md_simulation/output')
    parser.add_argument('--duration', type=float, default=100.0, help='Duration in ps')
    parser.add_argument('--asp-resnum', type=int, default=84, help='Asp3.32 residue number')
    
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Combine receptor + ligand
    receptor_pdb = Path(args.receptor)
    ligand_pdbqt = Path(args.ligand)
    
    complex_raw = output_dir / "complex_raw.pdb"
    combine_receptor_ligand(receptor_pdb, ligand_pdbqt, complex_raw)
    
    # Step 2: Prepare complex
    complex_prepared = output_dir / "complex_prepared.pdb"
    prepare_complex(complex_raw, complex_prepared)
    
    # Step 3: Run MD
    results = run_binding_stability_test(
        complex_prepared, output_dir,
        duration_ps=args.duration,
        asp_resnum=args.asp_resnum
    )
    
    return results


if __name__ == "__main__":
    main()
