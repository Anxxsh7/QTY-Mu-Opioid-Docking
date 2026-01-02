#!/usr/bin/env python3
"""
Visualize MD Simulation Results
================================
Creates PyMOL visualization script and RMSD plots.
"""

import json
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


def create_pymol_script(simulation_dir: Path, output_script: Path):
    """Create PyMOL script to visualize the trajectory."""
    
    pdb_file = simulation_dir / "final_complex.pdb"
    traj_file = simulation_dir / "trajectory.dcd"
    
    script = f'''# PyMOL Visualization Script for QTY μOR + Fentanyl MD
# Run with: pymol {output_script.name}

# Load the final structure
load {pdb_file}, complex

# Color by secondary structure
hide everything
show cartoon, polymer
color marine, polymer

# Show ligand (fentanyl) - residue name LIG
select ligand, resn LIG
show sticks, ligand
color orange, ligand

# Show binding site residues
select binding_site, (polymer within 5 of ligand)
show sticks, binding_site
color cyan, binding_site

# Highlight Asp84 (key salt bridge residue)
select asp84, resi 84 and resn ASP
show sticks, asp84
color red, asp84
label asp84 and name CA, "D84"

# Show water as small dots (optional - comment out if too slow)
# select waters, resn WAT or resn HOH
# show dots, waters
# color lightblue, waters
# set dot_radius, 0.3

# Hide waters for cleaner view
hide everything, resn WAT or resn HOH or resn Na+ or resn Cl-

# Load trajectory (if DCD reading is available)
# Note: Requires PyMOL with trajectory support
# load_traj {traj_file}, complex

# Center and zoom
center ligand
zoom ligand, 15

# Set nice rendering
set cartoon_fancy_helices, 1
set cartoon_side_chain_helper, 1
set ray_shadow, 0
bg_color white

# Save a nice image
# ray 1920, 1080
# png {simulation_dir}/md_snapshot.png

print("Visualization loaded!")
print("- Blue cartoon: Protein backbone")
print("- Orange sticks: Fentanyl")  
print("- Cyan sticks: Binding site residues")
print("- Red sticks: Asp84 (salt bridge)")
'''
    
    with open(output_script, 'w') as f:
        f.write(script)
    
    return output_script


def plot_rmsd(results_json: Path, output_plot: Path):
    """Create RMSD plot from MD results."""
    
    with open(results_json) as f:
        results = json.load(f)
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Get data from trajectory
    time_ps = results['trajectory']['time_ps']
    protein_rmsd = results['trajectory']['protein_rmsd']
    ligand_rmsd = results['trajectory']['ligand_rmsd']
    
    # Protein RMSD
    ax1.plot(time_ps, protein_rmsd, 'b-', linewidth=1.5, label='Protein (Cα)')
    ax1.axhline(np.mean(protein_rmsd), color='b', linestyle='--', alpha=0.5,
                label=f'Mean: {np.mean(protein_rmsd):.2f} Å')
    ax1.fill_between(time_ps, 
                     np.mean(protein_rmsd) - np.std(protein_rmsd),
                     np.mean(protein_rmsd) + np.std(protein_rmsd),
                     alpha=0.2, color='blue')
    ax1.set_ylabel('RMSD (Å)', fontsize=12)
    ax1.set_title('Protein Backbone Stability', fontsize=14)
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, max(protein_rmsd) * 1.2)
    
    # Ligand RMSD (already defined above)
    ax2.plot(time_ps, ligand_rmsd, 'orange', linewidth=1.5, label='Fentanyl')
    ax2.axhline(np.mean(ligand_rmsd), color='orange', linestyle='--', alpha=0.5,
                label=f'Mean: {np.mean(ligand_rmsd):.2f} Å')
    ax2.fill_between(time_ps,
                     np.mean(ligand_rmsd) - np.std(ligand_rmsd),
                     np.mean(ligand_rmsd) + np.std(ligand_rmsd),
                     alpha=0.2, color='orange')
    ax2.set_xlabel('Time (ps)', fontsize=12)
    ax2.set_ylabel('RMSD (Å)', fontsize=12)
    ax2.set_title('Ligand Stability in Binding Pocket', fontsize=14)
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, max(ligand_rmsd) * 1.2)
    
    # Add stability interpretation
    mean_lig_rmsd = np.mean(ligand_rmsd)
    if mean_lig_rmsd < 2:
        stability = "Very Stable Binding"
        color = 'green'
    elif mean_lig_rmsd < 4:
        stability = "Stable Binding"
        color = 'blue'
    elif mean_lig_rmsd < 6:
        stability = "Moderate Flexibility"
        color = 'orange'
    else:
        stability = "Unstable / Dissociating"
        color = 'red'
    
    fig.text(0.5, 0.02, f'Binding Assessment: {stability}', 
             ha='center', fontsize=14, color=color, fontweight='bold')
    
    plt.suptitle('QTY μOR + Fentanyl MD Simulation', fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    
    plt.savefig(output_plot, dpi=150, bbox_inches='tight')
    plt.savefig(output_plot.with_suffix('.pdf'), bbox_inches='tight')
    print(f"  Saved: {output_plot}")
    
    return fig


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Visualize MD simulation results")
    parser.add_argument('--sim-dir', type=str, 
                        default='05_md_simulation/output/complex_md/simulation',
                        help='Simulation output directory')
    parser.add_argument('--show', action='store_true', help='Display plots interactively')
    
    args = parser.parse_args()
    
    sim_dir = Path(args.sim_dir)
    
    print("=" * 60)
    print("MD SIMULATION VISUALIZATION")
    print("=" * 60)
    
    # Check files exist
    if not sim_dir.exists():
        print(f"Error: Directory not found: {sim_dir}")
        return
    
    results_json = sim_dir / "md_results.json"
    
    if not results_json.exists():
        print(f"Error: Results file not found: {results_json}")
        return
    
    # Create RMSD plot
    print("\n[1] Creating RMSD plot...")
    plot_file = sim_dir / "rmsd_plot.png"
    fig = plot_rmsd(results_json, plot_file)
    
    # Create PyMOL script
    print("\n[2] Creating PyMOL visualization script...")
    pymol_script = sim_dir / "visualize.pml"
    create_pymol_script(sim_dir, pymol_script)
    print(f"  Saved: {pymol_script}")
    
    # Print summary
    with open(results_json) as f:
        results = json.load(f)
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"\nProtein RMSD: {results['protein_rmsd']['final']:.2f} Å (final)")
    print(f"Ligand RMSD:  {results['ligand_rmsd']['final']:.2f} Å (final)")
    print(f"\nFiles created:")
    print(f"  📊 {plot_file}")
    print(f"  🔬 {pymol_script}")
    print(f"\nTo view in PyMOL:")
    print(f"  pymol {pymol_script}")
    
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
