#!/usr/bin/env python3
"""
Step 1: Download and Process PDB 8EF5
=====================================
Downloads the fentanyl-bound μOR structure and extracts:
- μOR receptor chain
- Fentanyl ligand
- Optional: pocket waters

Author: QTY Docking Pipeline
"""

import os
import sys
import argparse
import urllib.request
import warnings
from pathlib import Path

import numpy as np

try:
    from Bio.PDB import PDBParser, PDBIO, Select, NeighborSearch
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
except ImportError:
    print("Error: BioPython not installed. Run: pip install biopython")
    sys.exit(1)

try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    HAS_PDBFIXER = True
except ImportError:
    HAS_PDBFIXER = False
    print("Warning: PDBFixer not installed. Some features may be limited.")


# Suppress BioPython warnings
warnings.filterwarnings('ignore', category=PDBConstructionWarning)


class ReceptorLigandSelect(Select):
    """Select only receptor chain and ligand for output."""
    
    def __init__(self, receptor_chain, ligand_resname, keep_waters=False, 
                 water_atoms=None):
        self.receptor_chain = receptor_chain
        self.ligand_resname = ligand_resname
        self.keep_waters = keep_waters
        self.water_atoms = water_atoms or set()
        
    def accept_chain(self, chain):
        return chain.id == self.receptor_chain
    
    def accept_residue(self, residue):
        resname = residue.resname.strip()
        chain = residue.parent.id
        
        # Always keep the ligand
        if resname == self.ligand_resname:
            return True
            
        # Keep receptor residues (standard amino acids)
        if chain == self.receptor_chain:
            if resname in ['HOH', 'WAT', 'TIP3', 'SOL']:
                if self.keep_waters:
                    # Check if this water is in our list of pocket waters
                    for atom in residue:
                        if atom.serial_number in self.water_atoms:
                            return True
                return False
            # Skip other heteroatoms (cholesterol, etc.)
            hetflag = residue.id[0]
            if hetflag.startswith('H_') and resname != self.ligand_resname:
                return False
            return True
        return False


def download_pdb(pdb_id: str, output_dir: Path) -> Path:
    """Download PDB file from RCSB."""
    output_file = output_dir / f"{pdb_id}.pdb"
    
    if output_file.exists():
        print(f"PDB file already exists: {output_file}")
        return output_file
    
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"Downloading {pdb_id} from RCSB...")
    
    try:
        urllib.request.urlretrieve(url, output_file)
        print(f"Downloaded: {output_file}")
    except Exception as e:
        print(f"Error downloading PDB: {e}")
        # Try CIF format
        url_cif = f"https://files.rcsb.org/download/{pdb_id}.cif"
        output_cif = output_dir / f"{pdb_id}.cif"
        urllib.request.urlretrieve(url_cif, output_cif)
        print(f"Downloaded CIF format: {output_cif}")
        return output_cif
    
    return output_file


def find_pocket_waters(structure, ligand_resname: str, receptor_chain: str,
                       distance_cutoff: float = 5.0) -> set:
    """Find water molecules within distance_cutoff of the ligand."""
    
    # Get ligand atoms
    ligand_atoms = []
    water_atoms = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.resname.strip()
                if resname == ligand_resname:
                    ligand_atoms.extend(residue.get_atoms())
                elif resname in ['HOH', 'WAT', 'TIP3', 'SOL']:
                    water_atoms.extend(residue.get_atoms())
    
    if not ligand_atoms:
        print(f"Warning: No ligand found with resname '{ligand_resname}'")
        return set()
    
    if not water_atoms:
        print("No water molecules found in structure")
        return set()
    
    # Find waters near ligand
    ns = NeighborSearch(water_atoms)
    pocket_water_serials = set()
    
    for atom in ligand_atoms:
        nearby = ns.search(atom.coord, distance_cutoff)
        for water_atom in nearby:
            pocket_water_serials.add(water_atom.serial_number)
    
    print(f"Found {len(pocket_water_serials)} water atoms within {distance_cutoff}Å of ligand")
    return pocket_water_serials


def get_ligand_info(structure, ligand_resname: str = None):
    """Get information about ligands in the structure."""
    hetero_residues = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag = residue.id[0]
                if hetflag.startswith('H_'):
                    resname = residue.resname.strip()
                    if resname not in ['HOH', 'WAT', 'TIP3', 'SOL']:
                        hetero_residues.append({
                            'name': resname,
                            'chain': chain.id,
                            'id': residue.id,
                            'n_atoms': len(list(residue.get_atoms()))
                        })
    
    print("\nHetero residues found in structure:")
    for het in hetero_residues:
        print(f"  {het['name']} (chain {het['chain']}, {het['n_atoms']} atoms)")
    
    return hetero_residues


def get_chain_info(structure):
    """Get information about chains in the structure."""
    chain_info = []
    
    for model in structure:
        for chain in model:
            residues = [r for r in chain if r.id[0] == ' ']  # Standard residues only
            if residues:
                chain_info.append({
                    'id': chain.id,
                    'n_residues': len(residues),
                    'start': residues[0].id[1],
                    'end': residues[-1].id[1]
                })
    
    print("\nProtein chains found:")
    for c in chain_info:
        print(f"  Chain {c['id']}: {c['n_residues']} residues ({c['start']}-{c['end']})")
    
    return chain_info


def extract_receptor_ligand(input_pdb: Path, output_pdb: Path,
                           receptor_chain: str, ligand_resname: str,
                           keep_waters: bool = True,
                           water_distance: float = 5.0):
    """Extract receptor and ligand from the full structure."""
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('complex', str(input_pdb))
    
    # Get structure info
    get_chain_info(structure)
    get_ligand_info(structure)
    
    # Find pocket waters if requested
    pocket_waters = set()
    if keep_waters:
        pocket_waters = find_pocket_waters(structure, ligand_resname, 
                                           receptor_chain, water_distance)
    
    # Select and save
    selector = ReceptorLigandSelect(receptor_chain, ligand_resname,
                                    keep_waters, pocket_waters)
    
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_pdb), selector)
    
    print(f"\nExtracted receptor-ligand complex: {output_pdb}")
    
    return output_pdb


def fix_structure(input_pdb: Path, output_pdb: Path):
    """Use PDBFixer to add missing atoms and fix structure issues."""
    
    if not HAS_PDBFIXER:
        print("PDBFixer not available, skipping structure fixing")
        return input_pdb
    
    print("\nFixing structure with PDBFixer...")
    
    fixer = PDBFixer(filename=str(input_pdb))
    
    # Find and add missing residues (careful with termini)
    fixer.findMissingResidues()
    # Don't add missing terminal residues - can cause issues
    keys_to_remove = []
    for key in fixer.missingResidues:
        chain_id, res_id = key
        if res_id == 0 or res_id > 1000:  # Terminal regions
            keys_to_remove.append(key)
    for key in keys_to_remove:
        del fixer.missingResidues[key]
    
    # Find and add missing atoms
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    
    # Add missing hydrogens
    fixer.addMissingHydrogens(7.4)  # pH 7.4
    
    # Save fixed structure
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    print(f"Fixed structure saved: {output_pdb}")
    return output_pdb


def extract_ligand_separately(input_pdb: Path, output_dir: Path, 
                              ligand_resname: str) -> Path:
    """Extract just the ligand to a separate file for docking."""
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('complex', str(input_pdb))
    
    output_pdb = output_dir / f"{ligand_resname}_ligand.pdb"
    
    class LigandSelect(Select):
        def accept_residue(self, residue):
            return residue.resname.strip() == ligand_resname
    
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_pdb), LigandSelect())
    
    print(f"Ligand extracted: {output_pdb}")
    return output_pdb


def get_ligand_center(pdb_file: Path, ligand_resname: str) -> np.ndarray:
    """Get the center of mass of the ligand."""
    
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
        raise ValueError(f"Ligand {ligand_resname} not found")
    
    center = np.mean(coords, axis=0)
    print(f"\nLigand center of mass: {center}")
    return center


def main():
    parser = argparse.ArgumentParser(
        description="Download and process PDB 8EF5 for QTY docking"
    )
    parser.add_argument('--pdb-id', default='8EF5',
                        help='PDB ID to download (default: 8EF5)')
    parser.add_argument('--receptor-chain', default='R',
                        help='Chain ID for μOR receptor (default: R)')
    parser.add_argument('--ligand-name', default='ZPE',
                        help='Ligand residue name (default: ZPE for fentanyl in 8EF5)')
    parser.add_argument('--output-dir', default='01_prepare_structure/output',
                        help='Output directory')
    parser.add_argument('--keep-waters', action='store_true', default=True,
                        help='Keep pocket waters')
    parser.add_argument('--water-distance', type=float, default=5.0,
                        help='Distance cutoff for pocket waters (Å)')
    parser.add_argument('--no-fix', action='store_true',
                        help='Skip PDBFixer step')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Step 1: Download and Process PDB Structure")
    print("=" * 60)
    
    # Download PDB
    pdb_file = download_pdb(args.pdb_id, output_dir)
    
    # Extract receptor + ligand
    extracted_pdb = output_dir / f"{args.pdb_id}_receptor_ligand.pdb"
    extract_receptor_ligand(
        pdb_file, extracted_pdb,
        args.receptor_chain, args.ligand_name,
        args.keep_waters, args.water_distance
    )
    
    # Fix structure if PDBFixer available
    if not args.no_fix:
        fixed_pdb = output_dir / f"{args.pdb_id}_receptor_ligand_fixed.pdb"
        fix_structure(extracted_pdb, fixed_pdb)
        final_pdb = fixed_pdb
    else:
        final_pdb = extracted_pdb
    
    # Extract ligand separately
    extract_ligand_separately(final_pdb, output_dir, args.ligand_name)
    
    # Get ligand center for docking box
    center = get_ligand_center(final_pdb, args.ligand_name)
    
    # Save docking box info
    box_info = output_dir / "docking_box.txt"
    with open(box_info, 'w') as f:
        f.write(f"# Docking box centered on {args.ligand_name} in {args.pdb_id}\n")
        f.write(f"center_x = {center[0]:.3f}\n")
        f.write(f"center_y = {center[1]:.3f}\n")
        f.write(f"center_z = {center[2]:.3f}\n")
        f.write(f"size_x = 25\n")
        f.write(f"size_y = 25\n")
        f.write(f"size_z = 25\n")
    print(f"\nDocking box info saved: {box_info}")
    
    print("\n" + "=" * 60)
    print("Structure preparation complete!")
    print(f"Receptor-ligand complex: {final_pdb}")
    print("=" * 60)
    
    return final_pdb


if __name__ == "__main__":
    main()
