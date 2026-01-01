#!/usr/bin/env python3
"""
Step 4: Docking Pipeline with smina/AutoDock Vina
=================================================
Docks fentanyl into the ensemble of relaxed QTY structures.
Analyzes pose clustering and key interactions.

Author: QTY Docking Pipeline
"""

import os
import sys
import subprocess
import argparse
import json
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import warnings

import numpy as np

try:
    from Bio.PDB import PDBParser
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
except ImportError:
    print("Error: BioPython not installed")
    sys.exit(1)

warnings.filterwarnings('ignore', category=PDBConstructionWarning)


def check_docking_software() -> str:
    """Check which docking software is available."""
    
    # Check for smina
    try:
        result = subprocess.run(['smina', '--version'], 
                               capture_output=True, text=True)
        if result.returncode == 0:
            print("Found: smina")
            return 'smina'
    except FileNotFoundError:
        pass
    
    # Check for vina
    try:
        result = subprocess.run(['vina', '--version'], 
                               capture_output=True, text=True)
        if result.returncode == 0:
            print("Found: AutoDock Vina")
            return 'vina'
    except FileNotFoundError:
        pass
    
    print("WARNING: Neither smina nor vina found in PATH")
    print("Install smina: conda install -c conda-forge smina")
    print("Or download from: https://sourceforge.net/projects/smina/")
    return None


def pdb_to_pdbqt(input_pdb: Path, output_pdbqt: Path, 
                 is_receptor: bool = True) -> Path:
    """Convert PDB to PDBQT format for docking."""
    
    # Try using Open Babel
    try:
        cmd = ['obabel', str(input_pdb), '-O', str(output_pdbqt)]
        if is_receptor:
            cmd.extend(['-xr'])  # receptor mode
        else:
            cmd.extend(['-p', '7.4'])  # add hydrogens at pH 7.4 for ligand
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0 and output_pdbqt.exists():
            print(f"Converted: {output_pdbqt}")
            return output_pdbqt
    except FileNotFoundError:
        pass
    
    # Try using MGLTools prepare_receptor/prepare_ligand
    try:
        if is_receptor:
            cmd = ['prepare_receptor4.py', '-r', str(input_pdb), 
                   '-o', str(output_pdbqt)]
        else:
            cmd = ['prepare_ligand4.py', '-l', str(input_pdb), 
                   '-o', str(output_pdbqt)]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0 and output_pdbqt.exists():
            print(f"Converted: {output_pdbqt}")
            return output_pdbqt
    except FileNotFoundError:
        pass
    
    print(f"WARNING: Could not convert {input_pdb} to PDBQT")
    print("Install Open Babel: conda install -c conda-forge openbabel")
    return None


def get_box_from_ligand(pdb_file: Path, ligand_resname: str,
                        padding: float = 5.0) -> Dict:
    """Calculate docking box from ligand position."""
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('complex', str(pdb_file))
    
    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname.strip() == ligand_resname:
                    for atom in residue:
                        coords.append(atom.coord)
    
    if not coords:
        raise ValueError(f"Ligand {ligand_resname} not found in {pdb_file}")
    
    coords = np.array(coords)
    center = np.mean(coords, axis=0)
    
    # Box size based on ligand extent + padding
    extent = np.max(coords, axis=0) - np.min(coords, axis=0)
    size = extent + 2 * padding
    
    # Minimum box size
    size = np.maximum(size, [20, 20, 20])
    
    return {
        'center_x': float(center[0]),
        'center_y': float(center[1]),
        'center_z': float(center[2]),
        'size_x': float(size[0]),
        'size_y': float(size[1]),
        'size_z': float(size[2])
    }


def run_docking(receptor_pdbqt: Path, ligand_pdbqt: Path, 
                output_pdbqt: Path, box: Dict,
                exhaustiveness: int = 32,
                num_modes: int = 20,
                energy_range: float = 5.0,
                docking_program: str = 'smina') -> Dict:
    """Run docking with smina or vina."""
    
    cmd = [
        docking_program,
        '--receptor', str(receptor_pdbqt),
        '--ligand', str(ligand_pdbqt),
        '--out', str(output_pdbqt),
        '--center_x', str(box['center_x']),
        '--center_y', str(box['center_y']),
        '--center_z', str(box['center_z']),
        '--size_x', str(box['size_x']),
        '--size_y', str(box['size_y']),
        '--size_z', str(box['size_z']),
        '--exhaustiveness', str(exhaustiveness),
        '--num_modes', str(num_modes),
        '--energy_range', str(energy_range)
    ]
    
    # smina-specific options
    if docking_program == 'smina':
        cmd.extend(['--scoring', 'vinardo'])
    
    print(f"Running {docking_program}...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Docking failed: {result.stderr}")
        return None
    
    # Parse output for scores
    scores = parse_docking_output(result.stdout)
    
    return {
        'output_file': str(output_pdbqt),
        'scores': scores,
        'stdout': result.stdout
    }


def parse_docking_output(output: str) -> List[Dict]:
    """Parse docking output to extract scores."""
    
    scores = []
    in_results = False
    
    for line in output.split('\n'):
        if 'mode' in line.lower() and 'affinity' in line.lower():
            in_results = True
            continue
        
        if in_results:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    mode = int(parts[0])
                    affinity = float(parts[1])
                    scores.append({
                        'mode': mode,
                        'affinity': affinity,
                        'rmsd_lb': float(parts[2]) if len(parts) > 2 else None,
                        'rmsd_ub': float(parts[3]) if len(parts) > 3 else None
                    })
                except (ValueError, IndexError):
                    continue
    
    return scores


def split_docking_output(pdbqt_file: Path, output_dir: Path) -> List[Path]:
    """Split multi-model PDBQT into individual poses."""
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    poses = []
    current_pose = []
    pose_num = 0
    
    with open(pdbqt_file) as f:
        for line in f:
            if line.startswith('MODEL'):
                pose_num += 1
                current_pose = []
            elif line.startswith('ENDMDL'):
                if current_pose:
                    pose_file = output_dir / f"pose_{pose_num:02d}.pdbqt"
                    with open(pose_file, 'w') as pf:
                        pf.writelines(current_pose)
                    poses.append(pose_file)
            else:
                current_pose.append(line)
    
    return poses


def calculate_pose_rmsd(pose1_file: Path, pose2_file: Path) -> float:
    """Calculate RMSD between two poses."""
    
    def get_coords(pdbqt_file):
        coords = []
        with open(pdbqt_file) as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
        return np.array(coords)
    
    coords1 = get_coords(pose1_file)
    coords2 = get_coords(pose2_file)
    
    if len(coords1) != len(coords2):
        return float('inf')
    
    diff = coords1 - coords2
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    return rmsd


def cluster_poses(poses: List[Path], scores: List[Dict], 
                  rmsd_cutoff: float = 2.0) -> List[Dict]:
    """Cluster docking poses by RMSD."""
    
    n_poses = len(poses)
    
    # Calculate RMSD matrix
    rmsd_matrix = np.zeros((n_poses, n_poses))
    for i in range(n_poses):
        for j in range(i+1, n_poses):
            rmsd = calculate_pose_rmsd(poses[i], poses[j])
            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd
    
    # Simple clustering: group poses within rmsd_cutoff of best pose
    clusters = []
    assigned = set()
    
    # Sort by score (best first)
    sorted_indices = sorted(range(n_poses), 
                           key=lambda i: scores[i]['affinity'])
    
    for i in sorted_indices:
        if i in assigned:
            continue
        
        cluster = {
            'representative': i,
            'representative_file': str(poses[i]),
            'score': scores[i]['affinity'],
            'members': [i]
        }
        assigned.add(i)
        
        for j in sorted_indices:
            if j not in assigned and rmsd_matrix[i, j] < rmsd_cutoff:
                cluster['members'].append(j)
                assigned.add(j)
        
        cluster['size'] = len(cluster['members'])
        clusters.append(cluster)
    
    return clusters


def check_asp_interaction(receptor_pdb: Path, ligand_pdbqt: Path,
                          asp_resnum: int = 147,
                          distance_cutoff: float = 4.0) -> Dict:
    """Check if ligand maintains salt bridge with Asp3.32."""
    
    parser = PDBParser(QUIET=True)
    
    # Get Asp carboxylate oxygens
    structure = parser.get_structure('receptor', str(receptor_pdb))
    asp_oxygens = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[1] == asp_resnum:
                    for atom in residue:
                        if atom.name in ['OD1', 'OD2']:
                            asp_oxygens.append(atom.coord)
    
    if not asp_oxygens:
        return {'found': False, 'reason': f'Asp{asp_resnum} not found'}
    
    # Get ligand nitrogen (protonated amine)
    ligand_nitrogens = []
    with open(ligand_pdbqt) as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_type = line[77:79].strip()
                if atom_type.startswith('N'):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    ligand_nitrogens.append(np.array([x, y, z]))
    
    if not ligand_nitrogens:
        return {'found': False, 'reason': 'No nitrogen found in ligand'}
    
    # Check distances
    min_distance = float('inf')
    for asp_o in asp_oxygens:
        for lig_n in ligand_nitrogens:
            dist = np.linalg.norm(np.array(asp_o) - lig_n)
            min_distance = min(min_distance, dist)
    
    has_interaction = min_distance <= distance_cutoff
    
    return {
        'found': True,
        'has_salt_bridge': has_interaction,
        'min_distance': float(min_distance),
        'cutoff': distance_cutoff
    }


def dock_ensemble(ensemble_dir: Path, ligand_pdbqt: Path,
                  output_dir: Path, box: Dict,
                  docking_program: str = 'smina',
                  asp_resnum: int = 147) -> Dict:
    """Dock ligand into all ensemble structures."""
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find ensemble structures
    ensemble_files = sorted(ensemble_dir.glob("*.pdb"))
    if not ensemble_files:
        ensemble_files = sorted(ensemble_dir.glob("*.pdbqt"))
    
    print(f"\nDocking into {len(ensemble_files)} ensemble structures...")
    
    all_results = []
    
    for i, receptor_pdb in enumerate(ensemble_files):
        print(f"\n--- Receptor {i+1}/{len(ensemble_files)}: {receptor_pdb.name} ---")
        
        # Convert to PDBQT if needed
        if receptor_pdb.suffix == '.pdb':
            receptor_pdbqt = output_dir / f"{receptor_pdb.stem}.pdbqt"
            pdb_to_pdbqt(receptor_pdb, receptor_pdbqt, is_receptor=True)
        else:
            receptor_pdbqt = receptor_pdb
        
        if not receptor_pdbqt or not receptor_pdbqt.exists():
            print(f"  Skipping - could not prepare receptor")
            continue
        
        # Run docking
        docking_output = output_dir / f"{receptor_pdb.stem}_docked.pdbqt"
        result = run_docking(
            receptor_pdbqt, ligand_pdbqt, docking_output,
            box, docking_program=docking_program
        )
        
        if result is None:
            print(f"  Docking failed")
            continue
        
        # Split poses
        poses_dir = output_dir / f"{receptor_pdb.stem}_poses"
        poses = split_docking_output(docking_output, poses_dir)
        
        # Check Asp interaction for top pose
        if poses and result['scores']:
            asp_check = check_asp_interaction(
                receptor_pdb, poses[0], asp_resnum
            )
            result['asp_interaction'] = asp_check
            
            if asp_check.get('has_salt_bridge'):
                print(f"  ✓ Asp{asp_resnum} salt bridge: {asp_check['min_distance']:.2f} Å")
            else:
                print(f"  ✗ No Asp{asp_resnum} salt bridge (min dist: {asp_check.get('min_distance', 'N/A')} Å)")
        
        result['receptor'] = str(receptor_pdb)
        result['poses'] = [str(p) for p in poses]
        all_results.append(result)
        
        if result['scores']:
            print(f"  Best score: {result['scores'][0]['affinity']:.2f} kcal/mol")
    
    return {
        'n_receptors': len(ensemble_files),
        'n_successful': len(all_results),
        'results': all_results
    }


def analyze_ensemble_docking(results: Dict, output_dir: Path) -> Dict:
    """Analyze docking results across ensemble."""
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Collect all top scores
    top_scores = []
    asp_interactions = []
    
    for r in results['results']:
        if r['scores']:
            top_scores.append(r['scores'][0]['affinity'])
        
        if 'asp_interaction' in r:
            asp_interactions.append(r['asp_interaction'].get('has_salt_bridge', False))
    
    analysis = {
        'n_structures': results['n_receptors'],
        'n_successful_docking': results['n_successful'],
        'scores': {
            'mean': float(np.mean(top_scores)) if top_scores else None,
            'std': float(np.std(top_scores)) if top_scores else None,
            'min': float(np.min(top_scores)) if top_scores else None,
            'max': float(np.max(top_scores)) if top_scores else None,
            'all': top_scores
        },
        'asp_salt_bridge': {
            'n_with_interaction': sum(asp_interactions),
            'fraction': sum(asp_interactions) / len(asp_interactions) if asp_interactions else 0
        }
    }
    
    # Save analysis
    with open(output_dir / 'docking_analysis.json', 'w') as f:
        json.dump(analysis, f, indent=2)
    
    print("\n" + "=" * 60)
    print("DOCKING ANALYSIS SUMMARY")
    print("=" * 60)
    print(f"Structures docked: {analysis['n_successful_docking']}/{analysis['n_structures']}")
    
    if analysis['scores']['mean']:
        print(f"\nBinding scores (kcal/mol):")
        print(f"  Mean:  {analysis['scores']['mean']:.2f} ± {analysis['scores']['std']:.2f}")
        print(f"  Range: {analysis['scores']['min']:.2f} to {analysis['scores']['max']:.2f}")
    
    asp_frac = analysis['asp_salt_bridge']['fraction']
    print(f"\nAsp3.32 salt bridge preservation: {asp_frac*100:.1f}%")
    
    if asp_frac < 0.5:
        print("  ⚠️  WARNING: Salt bridge lost in majority of poses!")
        print("      This suggests binding may be compromised.")
    else:
        print("  ✓ Salt bridge maintained in most poses - good sign!")
    
    print("=" * 60)
    
    return analysis


def main():
    parser = argparse.ArgumentParser(
        description="Dock fentanyl into QTY μOR ensemble"
    )
    parser.add_argument('--ensemble-dir', required=True,
                        help='Directory with ensemble PDB files')
    parser.add_argument('--ligand', required=True,
                        help='Ligand PDB or PDBQT file')
    parser.add_argument('--output-dir', default='04_docking/output',
                        help='Output directory')
    parser.add_argument('--box-file',
                        help='File with docking box parameters')
    parser.add_argument('--reference-pdb',
                        help='Reference PDB for box calculation (with ligand)')
    parser.add_argument('--ligand-resname', default='ZPE',
                        help='Ligand residue name for box calculation')
    parser.add_argument('--asp-resnum', type=int, default=147,
                        help='Residue number of Asp3.32 (default: 147)')
    parser.add_argument('--exhaustiveness', type=int, default=32,
                        help='Docking exhaustiveness')
    
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Step 4: Ensemble Docking with smina/Vina")
    print("=" * 60)
    
    # Check for docking software
    docking_program = check_docking_software()
    if not docking_program:
        print("\nCannot proceed without docking software.")
        sys.exit(1)
    
    # Get docking box
    if args.box_file:
        box = {}
        with open(args.box_file) as f:
            for line in f:
                if '=' in line and not line.startswith('#'):
                    key, val = line.strip().split('=')
                    box[key.strip()] = float(val.strip())
    elif args.reference_pdb:
        box = get_box_from_ligand(Path(args.reference_pdb), args.ligand_resname)
    else:
        print("ERROR: Must provide --box-file or --reference-pdb")
        sys.exit(1)
    
    print(f"\nDocking box:")
    print(f"  Center: ({box['center_x']:.2f}, {box['center_y']:.2f}, {box['center_z']:.2f})")
    print(f"  Size: ({box['size_x']:.1f}, {box['size_y']:.1f}, {box['size_z']:.1f})")
    
    # Prepare ligand
    ligand_path = Path(args.ligand)
    if ligand_path.suffix == '.pdb':
        ligand_pdbqt = output_dir / f"{ligand_path.stem}.pdbqt"
        pdb_to_pdbqt(ligand_path, ligand_pdbqt, is_receptor=False)
    else:
        ligand_pdbqt = ligand_path
    
    # Run ensemble docking
    results = dock_ensemble(
        Path(args.ensemble_dir), ligand_pdbqt, output_dir,
        box, docking_program=docking_program,
        asp_resnum=args.asp_resnum
    )
    
    # Save raw results
    with open(output_dir / 'docking_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    # Analyze
    analysis = analyze_ensemble_docking(results, output_dir)
    
    print(f"\nResults saved in: {output_dir}")


if __name__ == "__main__":
    main()
