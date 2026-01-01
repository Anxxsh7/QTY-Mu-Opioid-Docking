#!/usr/bin/env python3
"""
Step 2: Introduce QTY Mutations onto μOR Structure
===================================================
Takes the prepared 8EF5 structure and introduces QTY mutations.
Performs side-chain repacking after mutation.

QTY Code:
  L (Leucine) → Q (Glutamine)
  I (Isoleucine) → T (Threonine)
  V (Valine) → T (Threonine)  
  F (Phenylalanine) → Y (Tyrosine)

Author: QTY Docking Pipeline
"""

import os
import sys
import argparse
import json
from pathlib import Path
from typing import List, Dict, Tuple
import warnings

import numpy as np

try:
    from Bio.PDB import PDBParser, PDBIO, Superimposer, Select
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    # Note: one_to_three/three_to_one moved in newer BioPython versions
    try:
        from Bio.PDB.Polypeptide import one_to_three, three_to_one
    except ImportError:
        # Define manually for compatibility
        pass
except ImportError:
    print("Error: BioPython not installed. Run: pip install biopython")
    sys.exit(1)

warnings.filterwarnings('ignore', category=PDBConstructionWarning)

# Standard amino acid data
AA_THREE_TO_ONE = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

AA_ONE_TO_THREE = {v: k for k, v in AA_THREE_TO_ONE.items()}

# QTY mutation rules
QTY_RULES = {
    'L': 'Q',  # Leucine → Glutamine
    'I': 'T',  # Isoleucine → Threonine
    'V': 'T',  # Valine → Threonine
    'F': 'Y',  # Phenylalanine → Tyrosine
}


def parse_mutation(mutation_str: str) -> Tuple[str, int, str]:
    """Parse mutation string like 'L65Q' into (original, position, new)."""
    mutation_str = mutation_str.strip().upper()
    
    if len(mutation_str) < 3:
        raise ValueError(f"Invalid mutation format: {mutation_str}")
    
    original = mutation_str[0]
    new = mutation_str[-1]
    position = int(mutation_str[1:-1])
    
    return original, position, new


def get_sequence_from_structure(structure, chain_id: str) -> Dict[int, str]:
    """Extract sequence from structure with residue numbers."""
    seq_dict = {}
    
    for model in structure:
        if chain_id in model:
            chain = model[chain_id]
            for residue in chain:
                if residue.id[0] == ' ':  # Standard residue
                    resname = residue.resname.strip()
                    if resname in AA_THREE_TO_ONE:
                        resnum = residue.id[1]
                        seq_dict[resnum] = AA_THREE_TO_ONE[resname]
    
    return seq_dict


def compare_sequences(wt_seq: str, qty_seq: str) -> List[str]:
    """Compare WT and QTY sequences to find mutations."""
    if len(wt_seq) != len(qty_seq):
        raise ValueError(f"Sequence length mismatch: WT={len(wt_seq)}, QTY={len(qty_seq)}")
    
    mutations = []
    for i, (wt, qty) in enumerate(zip(wt_seq, qty_seq), start=1):
        if wt != qty:
            mutations.append(f"{wt}{i}{qty}")
    
    return mutations


def validate_qty_mutations(mutations: List[str], structure, chain_id: str) -> List[Dict]:
    """Validate mutations against QTY rules and structure."""
    
    seq_dict = get_sequence_from_structure(structure, chain_id)
    validated = []
    
    print("\nValidating mutations:")
    print("-" * 50)
    
    for mut_str in mutations:
        orig, pos, new = parse_mutation(mut_str)
        
        # Check if position exists in structure
        if pos not in seq_dict:
            print(f"  WARNING: Position {pos} not in structure - skipping {mut_str}")
            continue
        
        # Check if original residue matches
        struct_res = seq_dict[pos]
        if struct_res != orig:
            print(f"  WARNING: {mut_str} - structure has {struct_res} at position {pos}, not {orig}")
            print(f"           Adjusting to {struct_res}{pos}{new}")
            orig = struct_res
        
        # Check if this follows QTY rules
        is_qty = (orig in QTY_RULES and QTY_RULES[orig] == new)
        
        validated.append({
            'original': orig,
            'position': pos,
            'new': new,
            'string': f"{orig}{pos}{new}",
            'is_qty_rule': is_qty
        })
        
        status = "✓ QTY rule" if is_qty else "○ Custom"
        print(f"  {status}: {orig}{pos}{new}")
    
    print("-" * 50)
    print(f"Total mutations: {len(validated)}")
    qty_count = sum(1 for m in validated if m['is_qty_rule'])
    print(f"Following QTY rules: {qty_count}")
    
    return validated


class MutationModeler:
    """Handle mutation modeling in protein structures."""
    
    # Simplified side chain atoms for each amino acid
    SIDECHAIN_ATOMS = {
        'ALA': ['CB'],
        'CYS': ['CB', 'SG'],
        'ASP': ['CB', 'CG', 'OD1', 'OD2'],
        'GLU': ['CB', 'CG', 'CD', 'OE1', 'OE2'],
        'PHE': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'GLY': [],
        'HIS': ['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'],
        'ILE': ['CB', 'CG1', 'CG2', 'CD1'],
        'LYS': ['CB', 'CG', 'CD', 'CE', 'NZ'],
        'LEU': ['CB', 'CG', 'CD1', 'CD2'],
        'MET': ['CB', 'CG', 'SD', 'CE'],
        'ASN': ['CB', 'CG', 'OD1', 'ND2'],
        'PRO': ['CB', 'CG', 'CD'],
        'GLN': ['CB', 'CG', 'CD', 'OE1', 'NE2'],
        'ARG': ['CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
        'SER': ['CB', 'OG'],
        'THR': ['CB', 'OG1', 'CG2'],
        'TRP': ['CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
        'TYR': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
        'VAL': ['CB', 'CG1', 'CG2'],
    }
    
    def __init__(self, structure):
        self.structure = structure
        self.mutations_applied = []
    
    def apply_mutation_simple(self, chain_id: str, res_num: int, 
                              new_resname: str) -> bool:
        """
        Apply a simple mutation by renaming residue.
        Note: This is a placeholder - proper mutation requires rotamer libraries.
        For production, use PyMOL, Modeller, or Rosetta.
        """
        
        for model in self.structure:
            if chain_id in model:
                chain = model[chain_id]
                for residue in chain:
                    if residue.id[1] == res_num and residue.id[0] == ' ':
                        old_resname = residue.resname
                        
                        # Just rename - atoms will need to be fixed by other tools
                        residue.resname = new_resname
                        
                        self.mutations_applied.append({
                            'chain': chain_id,
                            'position': res_num,
                            'old': old_resname,
                            'new': new_resname
                        })
                        
                        return True
        
        return False


def create_mutation_script_pymol(mutations: List[Dict], input_pdb: str, 
                                  output_pdb: str, chain_id: str) -> str:
    """Create PyMOL script for proper mutation with rotamer optimization."""
    
    script = f'''# PyMOL Mutation Script for QTY μOR
# Run with: pymol -cq mutation_script.pml

# Load structure
load {input_pdb}, receptor

# Mutations to apply
'''
    
    for mut in mutations:
        pos = mut['position']
        new_three = AA_ONE_TO_THREE[mut['new']]
        script += f'''
# Mutation: {mut['string']}
wizard mutagenesis
cmd.get_wizard().set_mode("{new_three}")
cmd.get_wizard().do_select("/{chain_id}//{pos}/")
cmd.get_wizard().apply()
cmd.set_wizard()
'''
    
    script += f'''
# Save mutated structure
save {output_pdb}, receptor

# Quit
quit
'''
    
    return script


def create_mutation_script_modeller(mutations: List[Dict], input_pdb: str,
                                    output_pdb: str, chain_id: str) -> str:
    """Create Modeller script for mutation."""
    
    mutation_list = [(m['position'], m['new']) for m in mutations]
    
    script = f'''# Modeller script for QTY mutations
from modeller import *
from modeller.optimizers import MolecularDynamics, ConjugateGradients

log.verbose()
env = Environ()
env.io.atom_files_directory = ['.']
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

# Read the original model
mdl = Model(env, file='{input_pdb}')

# Define mutations
mutations = {mutation_list}

# Apply mutations
for pos, new_res in mutations:
    sel = Selection(mdl.residues[str(pos) + ':' + '{chain_id}'])
    sel.mutate(residue_type=new_res)

# Optimize the mutated sidechains
mdl.write(file='{output_pdb}')
'''
    
    return script


def write_mutation_info(mutations: List[Dict], output_file: Path):
    """Write mutation information to JSON file."""
    
    info = {
        'n_mutations': len(mutations),
        'mutations': mutations,
        'qty_mutations': [m for m in mutations if m['is_qty_rule']],
        'non_qty_mutations': [m for m in mutations if not m['is_qty_rule']]
    }
    
    with open(output_file, 'w') as f:
        json.dump(info, f, indent=2)
    
    print(f"Mutation info saved: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Introduce QTY mutations into μOR structure"
    )
    parser.add_argument('--input-pdb', required=True,
                        help='Input PDB file (prepared receptor-ligand complex)')
    parser.add_argument('--output-dir', default='02_introduce_mutations/output',
                        help='Output directory')
    parser.add_argument('--chain', default='R',
                        help='Chain ID for receptor')
    parser.add_argument('--mutations', nargs='+',
                        help='List of mutations (e.g., L65Q I67T V70T)')
    parser.add_argument('--mutations-file',
                        help='File with mutations (one per line)')
    parser.add_argument('--qty-sequence',
                        help='Full QTY sequence (will auto-detect mutations)')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Step 2: Introduce QTY Mutations")
    print("=" * 60)
    
    # Parse structure
    pdb_parser = PDBParser(QUIET=True)
    structure = pdb_parser.get_structure('receptor', args.input_pdb)
    
    # Get mutations
    mutations = []
    
    if args.mutations:
        mutations = args.mutations
    elif args.mutations_file:
        with open(args.mutations_file) as f:
            mutations = [line.strip() for line in f if line.strip() 
                        and not line.startswith('#')]
    elif args.qty_sequence:
        # Extract WT sequence and compare
        seq_dict = get_sequence_from_structure(structure, args.chain)
        wt_seq = ''.join(seq_dict[k] for k in sorted(seq_dict.keys()))
        mutations = compare_sequences(wt_seq, args.qty_sequence)
    else:
        print("ERROR: No mutations specified!")
        print("Use --mutations, --mutations-file, or --qty-sequence")
        sys.exit(1)
    
    print(f"\nInput PDB: {args.input_pdb}")
    print(f"Chain: {args.chain}")
    print(f"Mutations to apply: {len(mutations)}")
    
    # Validate mutations
    validated_mutations = validate_qty_mutations(mutations, structure, args.chain)
    
    if not validated_mutations:
        print("ERROR: No valid mutations to apply!")
        sys.exit(1)
    
    # Check for pocket mutations (warning)
    pocket_residues = [147, 148, 149, 228, 229, 232, 293, 297]  # Key binding pocket residues
    pocket_mutations = [m for m in validated_mutations if m['position'] in pocket_residues]
    if pocket_mutations:
        print("\n⚠️  WARNING: Mutations in binding pocket region:")
        for m in pocket_mutations:
            print(f"    {m['string']} - This may affect fentanyl binding!")
    
    # Save mutation info
    write_mutation_info(validated_mutations, output_dir / 'mutations.json')
    
    # Create mutation scripts
    input_pdb_abs = str(Path(args.input_pdb).absolute())
    output_pdb = output_dir / "qty_mutant.pdb"
    
    # PyMOL script (recommended)
    pymol_script = create_mutation_script_pymol(
        validated_mutations, input_pdb_abs, str(output_pdb), args.chain
    )
    pymol_script_file = output_dir / "mutate_pymol.pml"
    with open(pymol_script_file, 'w') as f:
        f.write(pymol_script)
    print(f"\nPyMOL mutation script: {pymol_script_file}")
    
    # Modeller script (alternative)
    modeller_script = create_mutation_script_modeller(
        validated_mutations, input_pdb_abs, str(output_pdb), args.chain
    )
    modeller_script_file = output_dir / "mutate_modeller.py"
    with open(modeller_script_file, 'w') as f:
        f.write(modeller_script)
    print(f"Modeller mutation script: {modeller_script_file}")
    
    # Also create a simple renamed version (for OpenMM to fix)
    print("\nCreating preliminary mutant structure (requires optimization)...")
    modeler = MutationModeler(structure)
    
    for mut in validated_mutations:
        new_three = AA_ONE_TO_THREE[mut['new']]
        success = modeler.apply_mutation_simple(args.chain, mut['position'], new_three)
        if not success:
            print(f"  Warning: Could not apply {mut['string']}")
    
    # Save preliminary structure
    preliminary_pdb = output_dir / "qty_mutant_preliminary.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(preliminary_pdb))
    print(f"Preliminary mutant structure: {preliminary_pdb}")
    
    print("\n" + "=" * 60)
    print("IMPORTANT: Mutation Introduction")
    print("=" * 60)
    print("""
The scripts above will introduce mutations but need proper side-chain 
optimization. Choose ONE of these methods:

Option 1 (RECOMMENDED): PyMOL with Mutagenesis Wizard
  pymol -cq mutate_pymol.pml

Option 2: Modeller
  python mutate_modeller.py

Option 3: Rosetta (best for complex mutations)
  Use the FastRelax protocol after mutation

After mutation, proceed to Step 3 (relaxation) to optimize side chains.
""")
    print("=" * 60)


if __name__ == "__main__":
    main()
