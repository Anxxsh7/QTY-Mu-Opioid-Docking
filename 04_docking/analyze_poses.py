#!/usr/bin/env python3
"""
Analyze fentanyl docking poses - check binding mode and key interactions.
"""

import os
import json
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser
import warnings
warnings.filterwarnings('ignore')

def read_pdbqt_coords(pdbqt_file: Path) -> dict:
    """Read coordinates from PDBQT file."""
    atoms = {}
    with open(pdbqt_file) as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms[name] = np.array([x, y, z])
    return atoms

def find_nitrogen_ligand(pdbqt_file: Path) -> np.ndarray:
    """Find protonated nitrogen in fentanyl."""
    atoms = read_pdbqt_coords(pdbqt_file)
    
    # Fentanyl has a piperidine nitrogen that should be protonated
    # Look for nitrogen atoms
    nitrogens = {k: v for k, v in atoms.items() if k.startswith('N')}
    
    if nitrogens:
        # Return the first nitrogen (should be the piperidine N)
        for name, coord in nitrogens.items():
            print(f"  Found nitrogen: {name} at {coord}")
            return coord
    return None

def get_asp_oxygens(pdb_file: Path, asp_resnum: int) -> list:
    """Get Asp carboxylate oxygen coordinates."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('receptor', str(pdb_file))
    
    oxygens = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[1] == asp_resnum and residue.resname == 'ASP':
                    for atom in residue:
                        if atom.name in ['OD1', 'OD2']:
                            oxygens.append(atom.coord)
                            print(f"  Found Asp{asp_resnum} {atom.name} at {atom.coord}")
    return oxygens

def analyze_docking_pose(receptor_pdb: Path, ligand_pdbqt: Path, asp_resnum: int = 84):
    """Analyze a single docking pose."""
    
    print(f"\nAnalyzing: {ligand_pdbqt.name}")
    
    # Get ligand nitrogen
    lig_n = find_nitrogen_ligand(ligand_pdbqt)
    if lig_n is None:
        print("  No nitrogen found in ligand!")
        return None
    
    # Get Asp oxygens
    asp_oxygens = get_asp_oxygens(receptor_pdb, asp_resnum)
    if not asp_oxygens:
        print(f"  Asp{asp_resnum} not found in receptor!")
        return None
    
    # Calculate distances
    min_dist = float('inf')
    for oxy in asp_oxygens:
        dist = np.linalg.norm(lig_n - oxy)
        if dist < min_dist:
            min_dist = dist
    
    print(f"  N-Asp distance: {min_dist:.2f} Å")
    
    # Salt bridge cutoff is typically 4.0-4.5 Å
    has_salt_bridge = min_dist < 4.5
    print(f"  Salt bridge: {'YES ✓' if has_salt_bridge else 'NO'}")
    
    return {
        'min_distance': min_dist,
        'has_salt_bridge': has_salt_bridge
    }

def main():
    output_dir = Path("04_docking/output")
    ensemble_dir = Path("03_relaxation/output/ensemble")
    
    print("=" * 60)
    print("Detailed Fentanyl-μOR Binding Analysis")
    print("=" * 60)
    
    # Find all docked poses
    results = []
    
    for i in range(1, 6):
        receptor = ensemble_dir / f"ensemble_{i:03d}.pdb"
        docked = output_dir / f"ensemble_{i:03d}_docked.pdbqt"
        
        print(f"\nLooking for: {docked}")
        if docked.exists() and receptor.exists():
            result = analyze_docking_pose(receptor, docked, asp_resnum=84)
            if result:
                results.append(result)
        else:
            print(f"  Files not found!")
    
    if results:
        print("\n" + "=" * 60)
        print("SUMMARY")
        print("=" * 60)
        
        distances = [r['min_distance'] for r in results]
        salt_bridges = sum(1 for r in results if r['has_salt_bridge'])
        
        print(f"Poses analyzed: {len(results)}")
        print(f"N-Asp distances: {min(distances):.2f} - {max(distances):.2f} Å")
        print(f"Mean distance: {np.mean(distances):.2f} Å")
        print(f"Salt bridges: {salt_bridges}/{len(results)} ({100*salt_bridges/len(results):.0f}%)")
        
        if np.mean(distances) < 5.0:
            print("\n✓ Fentanyl is positioned near the critical Asp3.32 residue!")
            print("  This suggests the QTY mutant maintains a functional binding site.")
        elif np.mean(distances) < 8.0:
            print("\n⚠ Fentanyl is somewhat displaced from Asp3.32")
            print("  Binding may be weaker but still possible.")
        else:
            print("\n✗ Fentanyl is far from Asp3.32 - binding mode may be non-native.")

if __name__ == "__main__":
    main()
