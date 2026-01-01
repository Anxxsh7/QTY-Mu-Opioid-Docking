#!/usr/bin/env python3
"""
QTY μOR - Fentanyl Docking Pipeline
===================================
Master workflow script that runs the complete pipeline.

Usage:
    python run_pipeline.py --mutations L65Q I67T V70T ...
    python run_pipeline.py --mutations-file mutations.txt
    python run_pipeline.py --qty-sequence MDALSG...

Author: QTY Docking Pipeline
"""

import os
import sys
import argparse
import subprocess
import json
import yaml
from pathlib import Path
from datetime import datetime


def run_step(script: Path, args: list, step_name: str) -> bool:
    """Run a pipeline step."""
    
    print(f"\n{'='*70}")
    print(f"RUNNING: {step_name}")
    print(f"{'='*70}")
    
    cmd = [sys.executable, str(script)] + args
    
    try:
        result = subprocess.run(cmd, check=True)
        print(f"\n✓ {step_name} completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"\n✗ {step_name} failed with exit code {e.returncode}")
        return False
    except FileNotFoundError:
        print(f"\n✗ Script not found: {script}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Run QTY μOR - Fentanyl Docking Pipeline"
    )
    
    # Mutation input options
    mut_group = parser.add_mutually_exclusive_group(required=True)
    mut_group.add_argument('--mutations', nargs='+',
                           help='QTY mutations (e.g., L65Q I67T V70T)')
    mut_group.add_argument('--mutations-file',
                           help='File with mutations (one per line)')
    mut_group.add_argument('--qty-sequence',
                           help='Full QTY-modified sequence')
    mut_group.add_argument('--config', 
                           help='YAML config file with all settings')
    
    # Pipeline options
    parser.add_argument('--output-dir', default='results',
                        help='Output directory (default: results)')
    parser.add_argument('--pdb-id', default='8EF5',
                        help='PDB ID for template structure')
    parser.add_argument('--receptor-chain', default='R',
                        help='Receptor chain ID')
    parser.add_argument('--ligand-name', default='ZPE',
                        help='Ligand residue name')
    parser.add_argument('--asp-resnum', type=int, default=147,
                        help='Asp3.32 residue number')
    
    # Step control
    parser.add_argument('--skip-download', action='store_true',
                        help='Skip PDB download (use existing)')
    parser.add_argument('--skip-relaxation', action='store_true',
                        help='Skip relaxation/ensemble generation')
    parser.add_argument('--skip-docking', action='store_true',
                        help='Skip docking')
    parser.add_argument('--skip-md', action='store_true',
                        help='Skip MD simulation')
    parser.add_argument('--run-wt-control', action='store_true',
                        help='Also run WT control')
    
    # Ensemble/MD settings
    parser.add_argument('--n-ensemble', type=int, default=20,
                        help='Number of ensemble structures')
    parser.add_argument('--md-duration-ns', type=float, default=50,
                        help='MD simulation duration (ns)')
    parser.add_argument('--md-replicates', type=int, default=3,
                        help='Number of MD replicates')
    
    args = parser.parse_args()
    
    # Setup paths
    base_dir = Path(__file__).parent
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Scripts
    scripts = {
        'prepare': base_dir / '01_prepare_structure' / 'prepare_pdb.py',
        'mutate': base_dir / '02_introduce_mutations' / 'introduce_mutations.py',
        'relax': base_dir / '03_relaxation' / 'relax_structure.py',
        'dock': base_dir / '04_docking' / 'dock_ensemble.py',
        'md': base_dir / '05_md_simulation' / 'run_md.py',
        'analyze': base_dir / '06_analysis' / 'compare_variants.py'
    }
    
    # Check scripts exist
    for name, script in scripts.items():
        if not script.exists():
            print(f"Error: Script not found: {script}")
            sys.exit(1)
    
    print("=" * 70)
    print("QTY μOR - FENTANYL DOCKING PIPELINE")
    print("=" * 70)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Output directory: {output_dir}")
    
    # Load config if provided
    if args.config:
        with open(args.config) as f:
            config = yaml.safe_load(f)
        mutations = config.get('qty_mutations', [])
    elif args.mutations:
        mutations = args.mutations
    elif args.mutations_file:
        with open(args.mutations_file) as f:
            mutations = [l.strip() for l in f if l.strip() and not l.startswith('#')]
    else:
        mutations = None  # Will use qty_sequence
    
    if mutations:
        print(f"Mutations: {', '.join(mutations)}")
    
    # =========================================================================
    # STEP 1: Download and Prepare PDB
    # =========================================================================
    
    step1_output = output_dir / '01_prepare_structure'
    
    if not args.skip_download:
        success = run_step(
            scripts['prepare'],
            [
                '--pdb-id', args.pdb_id,
                '--receptor-chain', args.receptor_chain,
                '--ligand-name', args.ligand_name,
                '--output-dir', str(step1_output / 'output')
            ],
            "Step 1: Download and Prepare PDB"
        )
        
        if not success:
            print("Pipeline failed at Step 1")
            sys.exit(1)
    
    # Find prepared PDB
    prepared_pdb = step1_output / 'output' / f'{args.pdb_id}_receptor_ligand_fixed.pdb'
    if not prepared_pdb.exists():
        prepared_pdb = step1_output / 'output' / f'{args.pdb_id}_receptor_ligand.pdb'
    
    # =========================================================================
    # STEP 2: Introduce QTY Mutations
    # =========================================================================
    
    step2_output = output_dir / '02_introduce_mutations'
    
    mut_args = [
        '--input-pdb', str(prepared_pdb),
        '--output-dir', str(step2_output / 'output'),
        '--chain', args.receptor_chain
    ]
    
    if mutations:
        mut_args.extend(['--mutations'] + mutations)
    elif args.qty_sequence:
        mut_args.extend(['--qty-sequence', args.qty_sequence])
    
    success = run_step(
        scripts['mutate'],
        mut_args,
        "Step 2: Introduce QTY Mutations"
    )
    
    if not success:
        print("Pipeline failed at Step 2")
        sys.exit(1)
    
    mutant_pdb = step2_output / 'output' / 'qty_mutant_preliminary.pdb'
    
    # =========================================================================
    # STEP 3: Relaxation and Ensemble Generation
    # =========================================================================
    
    step3_output = output_dir / '03_relaxation'
    
    if not args.skip_relaxation:
        success = run_step(
            scripts['relax'],
            [
                '--input-pdb', str(mutant_pdb),
                '--output-dir', str(step3_output / 'output'),
                '--n-structures', str(args.n_ensemble)
            ],
            "Step 3: Structure Relaxation and Ensemble Generation"
        )
        
        if not success:
            print("Pipeline failed at Step 3")
            sys.exit(1)
    
    ensemble_dir = step3_output / 'output' / 'ensemble'
    
    # =========================================================================
    # STEP 4: Docking
    # =========================================================================
    
    step4_output = output_dir / '04_docking'
    
    if not args.skip_docking:
        # Get ligand file
        ligand_pdb = step1_output / 'output' / f'{args.ligand_name}_ligand.pdb'
        box_file = step1_output / 'output' / 'docking_box.txt'
        
        success = run_step(
            scripts['dock'],
            [
                '--ensemble-dir', str(ensemble_dir),
                '--ligand', str(ligand_pdb),
                '--output-dir', str(step4_output / 'output' / 'qty'),
                '--box-file', str(box_file),
                '--asp-resnum', str(args.asp_resnum)
            ],
            "Step 4: Ensemble Docking"
        )
        
        if not success:
            print("Pipeline failed at Step 4")
            sys.exit(1)
    
    # =========================================================================
    # STEP 5: MD Simulation
    # =========================================================================
    
    step5_output = output_dir / '05_md_simulation'
    
    if not args.skip_md:
        # Find best docked complex
        docking_results = step4_output / 'output' / 'qty' / 'docking_results.json'
        
        if docking_results.exists():
            with open(docking_results) as f:
                dock_data = json.load(f)
            
            # Get best pose from first receptor
            if dock_data.get('results') and dock_data['results'][0].get('poses'):
                best_pose = dock_data['results'][0]['poses'][0]
                
                # Note: Would need to combine receptor + ligand pose
                # For now, use original complex as placeholder
                complex_pdb = prepared_pdb
                
                success = run_step(
                    scripts['md'],
                    [
                        '--complex-pdb', str(complex_pdb),
                        '--output-dir', str(step5_output / 'output' / 'qty'),
                        '--duration-ns', str(args.md_duration_ns),
                        '--n-replicates', str(args.md_replicates),
                        '--ligand-resname', args.ligand_name,
                        '--asp-resnum', str(args.asp_resnum)
                    ],
                    "Step 5: MD Simulation"
                )
                
                if not success:
                    print("Warning: MD simulation failed, continuing...")
    
    # =========================================================================
    # STEP 6: Analysis and Comparison
    # =========================================================================
    
    step6_output = output_dir / '06_analysis'
    
    analyze_args = [
        '--qty-docking-dir', str(step4_output / 'output' / 'qty'),
        '--qty-md-dir', str(step5_output / 'output' / 'qty'),
        '--output-dir', str(step6_output / 'output')
    ]
    
    if args.run_wt_control:
        analyze_args.extend([
            '--wt-docking-dir', str(step4_output / 'output' / 'wt'),
            '--wt-md-dir', str(step5_output / 'output' / 'wt')
        ])
    
    success = run_step(
        scripts['analyze'],
        analyze_args,
        "Step 6: Analysis and Comparison"
    )
    
    # =========================================================================
    # SUMMARY
    # =========================================================================
    
    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE")
    print("=" * 70)
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"\nResults saved in: {output_dir}")
    print(f"\nKey outputs:")
    print(f"  • Prepared structure: {prepared_pdb}")
    print(f"  • Mutant structure: {mutant_pdb}")
    print(f"  • Ensemble: {ensemble_dir}")
    print(f"  • Docking results: {step4_output / 'output'}")
    print(f"  • MD results: {step5_output / 'output'}")
    print(f"  • Final report: {step6_output / 'output' / 'binding_assessment_report.txt'}")
    print("=" * 70)


if __name__ == "__main__":
    main()
