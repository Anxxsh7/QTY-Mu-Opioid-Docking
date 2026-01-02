# PyMOL Visualization Script - Fixed Ligand Display
# Uses mol2 file for proper ligand bonds

# Load the complex
load /Users/aneeshbondugula/Downloads/QTY_docking/05_md_simulation/output/complex_md/system/complex.pdb, complex

# Load trajectory
load_traj /Users/aneeshbondugula/Downloads/QTY_docking/05_md_simulation/output/complex_md/simulation/trajectory.dcd, complex

# Hide everything first
hide everything

# Show protein as cartoon
select protein, polymer
show cartoon, protein
color marine, protein

# Select ligand and fix its display
select ligand, resn LIG

# Rebuild bonds for ligand based on distance
unbond ligand, ligand
h_add ligand
bond ligand, ligand, 2  # Connect atoms within 2 Angstroms

# Alternative: use valence to guess bonds
set valence, 1
rebuild ligand

# Show ligand nicely
show sticks, ligand
color orange, ligand
set stick_radius, 0.2, ligand

# Show binding site
select binding_site, (byres polymer within 5 of ligand)
show sticks, binding_site
color cyan, binding_site
set stick_radius, 0.15, binding_site

# Highlight Asp84
select asp84, resi 84 and resn ASP
color red, asp84

# Hide waters
hide everything, resn WAT or resn HOH or resn Na+ or resn Cl-

# Center view
center ligand
zoom ligand, 20

# Settings
set cartoon_fancy_helices, 1
bg_color white
set movie_fps, 10

print "=== VISUALIZATION LOADED ==="
print "Frames: " + str(cmd.count_frames())
