#!/usr/bin/env python3
"""
Remove ligand from PDB for relaxation.
We'll dock the ligand back in step 4.
"""

import sys
from pathlib import Path

def remove_ligand(input_pdb: Path, output_pdb: Path, ligand_name: str = "7V7"):
    """Remove ligand residue from PDB file."""
    
    with open(input_pdb, 'r') as f:
        lines = f.readlines()
    
    protein_lines = []
    ligand_count = 0
    
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            res_name = line[17:20].strip()
            if res_name == ligand_name:
                ligand_count += 1
                continue
        elif line.startswith('CONECT'):
            # Skip CONECT records for ligand
            continue
        protein_lines.append(line)
    
    with open(output_pdb, 'w') as f:
        f.writelines(protein_lines)
    
    print(f"Removed {ligand_count} atoms of {ligand_name}")
    print(f"Apo structure saved to: {output_pdb}")
    return output_pdb


if __name__ == "__main__":
    input_pdb = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("02_introduce_mutations/output/qty_mutant.pdb")
    output_pdb = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("03_relaxation/output/qty_mutant_apo.pdb")
    ligand = sys.argv[3] if len(sys.argv) > 3 else "7V7"
    
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    remove_ligand(input_pdb, output_pdb, ligand)
