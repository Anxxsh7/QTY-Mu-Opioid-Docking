#!/usr/bin/env python3
"""
Full MD Simulation of QTY μOR + Fentanyl Complex
================================================
Uses AmberTools (antechamber) for fentanyl parameterization
and OpenMM for simulation in explicit water.

Author: QTY Docking Pipeline
"""

import os
import sys
import subprocess
import json
import shutil
from pathlib import Path
from typing import Dict, Tuple, Optional
import numpy as np

try:
    import openmm as mm
    from openmm import app, unit
    from openmm.app import (PDBFile, ForceField, Modeller, Simulation,
                            PME, HBonds, StateDataReporter, DCDReporter)
except ImportError:
    print("Error: OpenMM not installed. Run: pip install openmm")
    sys.exit(1)

try:
    from pdbfixer import PDBFixer
except ImportError:
    print("Error: PDBFixer not installed. Run: pip install pdbfixer")
    sys.exit(1)

try:
    import parmed as pmd
    HAS_PARMED = True
except ImportError:
    HAS_PARMED = False
    print("Warning: parmed not installed - some features limited")


def check_ambertools():
    """Check if AmberTools is available."""
    try:
        result = subprocess.run(['antechamber', '-h'], capture_output=True, text=True)
        return True
    except FileNotFoundError:
        return False


def convert_pdbqt_to_pdb(pdbqt_file: Path, pdb_file: Path) -> Path:
    """Convert PDBQT to PDB format."""
    
    with open(pdbqt_file) as f_in, open(pdb_file, 'w') as f_out:
        atom_num = 1
        for line in f_in:
            if line.startswith(('ATOM', 'HETATM')):
                # Extract atom info from PDBQT
                atom_name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                
                # Determine element from atom name
                element = atom_name[0] if atom_name[0].isalpha() else atom_name[1] if len(atom_name) > 1 else 'C'
                
                # Write PDB format
                pdb_line = f"HETATM{atom_num:5d} {atom_name:4s} LIG A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n"
                f_out.write(pdb_line)
                atom_num += 1
        f_out.write("END\n")
    
    return pdb_file


def prepare_receptor_for_amber(receptor_pdb: Path, output_pdb: Path) -> Path:
    """
    Prepare receptor PDB for AMBER/tleap.
    
    Strategy: Remove ALL hydrogens from the protein.
    tleap will add them with correct names and types.
    This avoids naming conflicts between OpenMM and AMBER conventions.
    """
    print("\n--- Preparing Receptor for AMBER ---")
    
    receptor_pdb = Path(receptor_pdb).resolve()
    output_pdb = Path(output_pdb).resolve()
    
    # Read original PDB
    with open(receptor_pdb) as f:
        lines = f.readlines()
    
    removed_hydrogens = 0
    output_lines = []
    
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            atom_name = line[12:16].strip()
            element = line[76:78].strip() if len(line) > 76 else ''
            
            # Skip all hydrogens - tleap will add them correctly
            # Hydrogen detection: element is H, or atom name starts with H (but not HIS atoms like HA)
            is_hydrogen = (element == 'H' or 
                          (atom_name.startswith('H') and atom_name not in ['HG', 'HG1']))
            
            if is_hydrogen:
                removed_hydrogens += 1
                continue
            
            output_lines.append(line)
        elif line.startswith(('TER', 'END', 'REMARK')):
            output_lines.append(line)
    
    # Write cleaned PDB (heavy atoms only)
    with open(output_pdb, 'w') as f:
        f.writelines(output_lines)
    
    print(f"  Prepared receptor: {output_pdb}")
    print(f"  - Removed {removed_hydrogens} hydrogen atoms (tleap will add them)")
    return output_pdb


def parameterize_fentanyl(ligand_pdb: Path, output_dir: Path, net_charge: int = 1) -> Dict:
    """
    Parameterize fentanyl using AmberTools antechamber.
    
    Fentanyl is protonated at physiological pH, so net_charge = +1
    """
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("\n--- Parameterizing Fentanyl with AmberTools ---")
    
    # Use absolute paths
    ligand_pdb = Path(ligand_pdb).resolve()
    output_dir = Path(output_dir).resolve()
    
    mol2_file = output_dir / "fentanyl.mol2"
    frcmod_file = output_dir / "fentanyl.frcmod"
    
    # Run antechamber to generate mol2 with AM1-BCC charges
    print("  Running antechamber (AM1-BCC charges)...")
    print(f"  Input: {ligand_pdb}")
    antechamber_cmd = [
        'antechamber',
        '-i', str(ligand_pdb),
        '-fi', 'pdb',
        '-o', str(mol2_file),
        '-fo', 'mol2',
        '-c', 'bcc',           # AM1-BCC charges
        '-s', '2',             # Status verbosity
        '-nc', str(net_charge), # Net charge (+1 for protonated fentanyl)
        '-at', 'gaff2'         # GAFF2 atom types
    ]
    
    result = subprocess.run(antechamber_cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"  Error in antechamber: {result.stderr}")
        # Try without BCC (faster but less accurate)
        print("  Retrying with Gasteiger charges...")
        antechamber_cmd[antechamber_cmd.index('bcc')] = 'gas'
        result = subprocess.run(antechamber_cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            raise RuntimeError(f"Antechamber failed: {result.stderr}")
    
    print(f"  Generated: {mol2_file}")
    
    # Run parmchk2 to generate missing parameters
    print("  Running parmchk2...")
    parmchk_cmd = [
        'parmchk2',
        '-i', str(mol2_file),
        '-f', 'mol2',
        '-o', str(frcmod_file),
        '-s', 'gaff2'
    ]
    
    result = subprocess.run(parmchk_cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"  Warning in parmchk2: {result.stderr}")
    
    print(f"  Generated: {frcmod_file}")
    
    # Clean up antechamber temp files
    for temp_file in output_dir.glob("ANTECHAMBER*"):
        temp_file.unlink()
    for temp_file in output_dir.glob("ATOMTYPE*"):
        temp_file.unlink()
    for temp_file in output_dir.glob("sqm*"):
        temp_file.unlink()
    
    return {
        'mol2': mol2_file,
        'frcmod': frcmod_file,
        'success': mol2_file.exists() and frcmod_file.exists()
    }


def create_complex_with_tleap(receptor_pdb: Path, ligand_mol2: Path, 
                               ligand_frcmod: Path, output_dir: Path) -> Dict:
    """
    Use tleap to create parameterized complex with AMBER parameters.
    """
    
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Convert all paths to absolute
    receptor_pdb = Path(receptor_pdb).resolve()
    ligand_mol2 = Path(ligand_mol2).resolve()
    ligand_frcmod = Path(ligand_frcmod).resolve()
    
    print("\n--- Creating Complex with tleap ---")
    
    # Create tleap input file
    tleap_input = output_dir / "tleap.in"
    prmtop_file = output_dir / "complex.prmtop"
    inpcrd_file = output_dir / "complex.inpcrd"
    complex_pdb = output_dir / "complex.pdb"
    
    tleap_script = f"""
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p

# Load ligand parameters
loadamberparams {ligand_frcmod}
LIG = loadmol2 {ligand_mol2}

# Load receptor
receptor = loadpdb {receptor_pdb}

# Create complex
complex = combine {{receptor LIG}}

# Add water box (12 Å padding)
solvatebox complex TIP3PBOX 12.0

# Add ions (0.15 M NaCl)
addionsrand complex Na+ 0
addionsrand complex Cl- 0

# Save files
saveamberparm complex {prmtop_file} {inpcrd_file}
savepdb complex {complex_pdb}

quit
"""
    
    with open(tleap_input, 'w') as f:
        f.write(tleap_script)
    
    print("  Running tleap...")
    result = subprocess.run(
        ['tleap', '-f', str(tleap_input)],
        capture_output=True, text=True
    )
    
    if not prmtop_file.exists():
        print(f"  tleap output: {result.stdout}")
        print(f"  tleap errors: {result.stderr}")
        raise RuntimeError("tleap failed to create topology files")
    
    print(f"  Generated: {prmtop_file}")
    print(f"  Generated: {inpcrd_file}")
    print(f"  Generated: {complex_pdb}")
    
    return {
        'prmtop': prmtop_file,
        'inpcrd': inpcrd_file,
        'pdb': complex_pdb,
        'success': True
    }


def run_md_with_amber_params(prmtop: Path, inpcrd: Path, output_dir: Path,
                              duration_ps: float = 1000.0,
                              temperature: float = 300.0) -> Dict:
    """
    Run MD simulation using AMBER parameter files with OpenMM.
    """
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("\n" + "=" * 70)
    print("MD SIMULATION: QTY μOR + Fentanyl in Water")
    print("=" * 70)
    print(f"Duration: {duration_ps} ps ({duration_ps/1000:.2f} ns)")
    print(f"Temperature: {temperature} K")
    print(f"Topology: {prmtop}")
    
    # Load AMBER files with parmed
    print("\n--- Loading AMBER Parameters ---")
    
    if HAS_PARMED:
        parm = pmd.load_file(str(prmtop), str(inpcrd))
        system = parm.createSystem(
            nonbondedMethod=PME,
            nonbondedCutoff=10 * unit.angstrom,
            constraints=HBonds
        )
        topology = parm.topology
        positions = parm.positions
    else:
        # Use OpenMM's AMBER file readers
        prmtop_obj = app.AmberPrmtopFile(str(prmtop))
        inpcrd_obj = app.AmberInpcrdFile(str(inpcrd))
        
        system = prmtop_obj.createSystem(
            nonbondedMethod=PME,
            nonbondedCutoff=10 * unit.angstrom,
            constraints=HBonds
        )
        topology = prmtop_obj.topology
        positions = inpcrd_obj.positions
    
    n_atoms = topology.getNumAtoms() if hasattr(topology, 'getNumAtoms') else len(list(topology.atoms()))
    print(f"  Loaded system with {n_atoms} atoms")
    
    # Add barostat for NPT
    barostat = mm.MonteCarloBarostat(1 * unit.bar, temperature * unit.kelvin)
    system.addForce(barostat)
    
    # Create integrator
    integrator = mm.LangevinMiddleIntegrator(
        temperature * unit.kelvin,
        1 / unit.picosecond,
        2 * unit.femtoseconds
    )
    
    # Create simulation
    simulation = Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)
    
    # Minimize
    print("\n--- Energy Minimization ---")
    print("  Minimizing...")
    simulation.minimizeEnergy(maxIterations=1000)
    
    min_state = simulation.context.getState(getEnergy=True)
    print(f"  Minimized energy: {min_state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole):.0f} kJ/mol")
    
    # Get initial state for RMSD
    state = simulation.context.getState(getPositions=True)
    initial_positions = state.getPositions()
    
    # Get CA indices for protein RMSD
    ca_indices = []
    for i, atom in enumerate(topology.atoms()):
        if atom.name == 'CA':
            ca_indices.append(i)
    
    # Get ligand indices
    lig_indices = []
    for i, atom in enumerate(topology.atoms()):
        if hasattr(atom.residue, 'name') and atom.residue.name in ['LIG', 'UNL', '7V7']:
            lig_indices.append(i)
    
    initial_ca = np.array([initial_positions[i].value_in_unit(unit.angstrom) for i in ca_indices])
    initial_lig = np.array([initial_positions[i].value_in_unit(unit.angstrom) for i in lig_indices]) if lig_indices else None
    
    print(f"  Tracking {len(ca_indices)} CA atoms, {len(lig_indices)} ligand atoms")
    
    # Set velocities
    simulation.context.setVelocitiesToTemperature(temperature * unit.kelvin)
    
    # Setup reporters
    dcd_file = output_dir / "trajectory.dcd"
    simulation.reporters.append(DCDReporter(str(dcd_file), 5000))  # Every 10 ps
    simulation.reporters.append(StateDataReporter(
        str(output_dir / 'simulation.log'), 2500,
        step=True, time=True, potentialEnergy=True, temperature=True,
        speed=True
    ))
    
    # Run production MD
    print(f"\n--- Production MD ({duration_ps} ps) ---")
    
    timestep_fs = 2.0
    n_steps = int(duration_ps * 1000 / timestep_fs)
    report_interval = max(1, n_steps // 20)
    
    trajectory_data = {
        'time_ps': [],
        'protein_rmsd': [],
        'ligand_rmsd': [],
        'potential_energy': []
    }
    
    for i in range(20):
        simulation.step(report_interval)
        
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        positions = state.getPositions()
        time_ps = (i + 1) * report_interval * timestep_fs / 1000
        
        # Protein RMSD
        current_ca = np.array([positions[idx].value_in_unit(unit.angstrom) for idx in ca_indices])
        protein_rmsd = np.sqrt(np.mean(np.sum((current_ca - initial_ca)**2, axis=1)))
        
        # Ligand RMSD
        if initial_lig is not None and len(lig_indices) > 0:
            current_lig = np.array([positions[idx].value_in_unit(unit.angstrom) for idx in lig_indices])
            ligand_rmsd = np.sqrt(np.mean(np.sum((current_lig - initial_lig)**2, axis=1)))
        else:
            ligand_rmsd = 0.0
        
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        
        trajectory_data['time_ps'].append(time_ps)
        trajectory_data['protein_rmsd'].append(protein_rmsd)
        trajectory_data['ligand_rmsd'].append(ligand_rmsd)
        trajectory_data['potential_energy'].append(pe)
        
        print(f"  {time_ps:7.1f} ps | Protein RMSD: {protein_rmsd:5.2f} Å | Ligand RMSD: {ligand_rmsd:5.2f} Å")
    
    # Save final structure
    final_positions = simulation.context.getState(getPositions=True).getPositions()
    final_pdb = output_dir / "final_complex.pdb"
    with open(final_pdb, 'w') as f:
        PDBFile.writeFile(topology, final_positions, f)
    
    # Analysis
    print("\n" + "=" * 70)
    print("SIMULATION RESULTS")
    print("=" * 70)
    
    protein_rmsds = trajectory_data['protein_rmsd']
    ligand_rmsds = trajectory_data['ligand_rmsd']
    
    print(f"\nProtein Stability (CA RMSD):")
    print(f"  Final:   {protein_rmsds[-1]:.2f} Å")
    print(f"  Mean:    {np.mean(protein_rmsds):.2f} ± {np.std(protein_rmsds):.2f} Å")
    
    if lig_indices:
        print(f"\nLigand Stability (RMSD):")
        print(f"  Final:   {ligand_rmsds[-1]:.2f} Å")
        print(f"  Mean:    {np.mean(ligand_rmsds):.2f} ± {np.std(ligand_rmsds):.2f} Å")
        
        if np.mean(ligand_rmsds) < 3.0:
            print("  ✅ Fentanyl STABLE in binding pocket")
            ligand_stable = True
        elif np.mean(ligand_rmsds) < 5.0:
            print("  ⚠️  Fentanyl showing FLEXIBILITY")
            ligand_stable = True
        else:
            print("  ❌ Fentanyl may be LEAVING binding pocket")
            ligand_stable = False
    else:
        ligand_stable = None
    
    # Save results
    results = {
        'simulation': {
            'duration_ps': duration_ps,
            'temperature_K': temperature,
            'n_atoms': n_atoms,
            'n_ca_atoms': len(ca_indices),
            'n_ligand_atoms': len(lig_indices)
        },
        'protein_rmsd': {
            'final': protein_rmsds[-1],
            'mean': np.mean(protein_rmsds),
            'std': np.std(protein_rmsds)
        },
        'ligand_rmsd': {
            'final': ligand_rmsds[-1] if lig_indices else None,
            'mean': np.mean(ligand_rmsds) if lig_indices else None,
            'std': np.std(ligand_rmsds) if lig_indices else None
        },
        'trajectory': trajectory_data,
        'ligand_stable': ligand_stable,
        'files': {
            'trajectory': str(dcd_file),
            'final_structure': str(final_pdb)
        }
    }
    
    with open(output_dir / 'md_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n📁 Output files saved to: {output_dir}")
    
    return results


def run_full_complex_md(receptor_pdb: Path, ligand_pdbqt: Path, 
                        output_dir: Path, duration_ps: float = 1000.0):
    """
    Full pipeline: parameterize ligand, build complex, run MD.
    """
    
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("QTY μOR + FENTANYL MD SIMULATION PIPELINE")
    print("=" * 70)
    
    # Check AmberTools
    if not check_ambertools():
        print("\n❌ ERROR: AmberTools not found!")
        print("   Install with: conda install -c conda-forge ambertools")
        sys.exit(1)
    
    print("✅ AmberTools found")
    
    # Step 1: Convert PDBQT to PDB
    print("\n[Step 1/5] Converting ligand format...")
    ligand_pdb = output_dir / "fentanyl.pdb"
    convert_pdbqt_to_pdb(ligand_pdbqt, ligand_pdb)
    print(f"  Converted: {ligand_pdb}")
    
    # Step 2: Prepare receptor for AMBER
    print("\n[Step 2/5] Preparing receptor...")
    prepared_receptor = output_dir / "receptor_prepared.pdb"
    prepare_receptor_for_amber(receptor_pdb, prepared_receptor)
    
    # Step 3: Parameterize fentanyl
    print("\n[Step 3/5] Parameterizing fentanyl...")
    param_dir = output_dir / "parameters"
    params = parameterize_fentanyl(ligand_pdb, param_dir)
    
    if not params['success']:
        print("❌ Parameterization failed!")
        sys.exit(1)
    
    # Step 4: Create solvated complex with tleap
    print("\n[Step 4/5] Building solvated complex...")
    complex_files = create_complex_with_tleap(
        prepared_receptor, params['mol2'], params['frcmod'], 
        output_dir / "system"
    )
    
    # Step 5: Run MD
    print("\n[Step 5/5] Running MD simulation...")
    results = run_md_with_amber_params(
        complex_files['prmtop'],
        complex_files['inpcrd'],
        output_dir / "simulation",
        duration_ps=duration_ps
    )
    
    return results


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Full MD simulation of QTY μOR + Fentanyl complex"
    )
    parser.add_argument('--receptor', required=True, help='Receptor PDB file')
    parser.add_argument('--ligand', required=True, help='Ligand PDBQT file')
    parser.add_argument('--output-dir', default='05_md_simulation/output/complex_md')
    parser.add_argument('--duration', type=float, default=1000.0, 
                        help='Duration in ps (default: 1000 = 1 ns)')
    
    args = parser.parse_args()
    
    results = run_full_complex_md(
        Path(args.receptor),
        Path(args.ligand),
        Path(args.output_dir),
        duration_ps=args.duration
    )
    
    return results


if __name__ == "__main__":
    main()
