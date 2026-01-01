#!/usr/bin/env python3
"""
Extract fentanyl ligand from original PDB for docking.
"""

import sys
from pathlib import Path

def extract_ligand(input_pdb: Path, output_pdb: Path, ligand_name: str = "7V7"):
    """Extract ligand residue from PDB file."""
    
    with open(input_pdb, 'r') as f:
        lines = f.readlines()
    
    ligand_lines = []
    
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            res_name = line[17:20].strip()
            if res_name == ligand_name:
                ligand_lines.append(line)
    
    if not ligand_lines:
        print(f"ERROR: Ligand {ligand_name} not found!")
        return None
    
    ligand_lines.append("END\n")
    
    with open(output_pdb, 'w') as f:
        f.writelines(ligand_lines)
    
    print(f"Extracted {len(ligand_lines)-1} atoms of {ligand_name}")
    print(f"Ligand saved to: {output_pdb}")
    return output_pdb


if __name__ == "__main__":
    input_pdb = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("01_prepare_pdb/output/8EF5_processed.pdb")
    output_pdb = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("04_docking/output/fentanyl.pdb")
    ligand = sys.argv[3] if len(sys.argv) > 3 else "7V7"
    
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    extract_ligand(input_pdb, output_pdb, ligand)
