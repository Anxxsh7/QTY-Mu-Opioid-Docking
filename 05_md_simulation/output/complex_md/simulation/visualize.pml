# PyMOL Visualization Script for QTY μOR + Fentanyl MD
# Run with: pymol visualize.pml

# Load the final structure
load 05_md_simulation/output/complex_md/simulation/final_complex.pdb, complex

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
# load_traj 05_md_simulation/output/complex_md/simulation/trajectory.dcd, complex

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
# png 05_md_simulation/output/complex_md/simulation/md_snapshot.png

print("Visualization loaded!")
print("- Blue cartoon: Protein backbone")
print("- Orange sticks: Fentanyl")  
print("- Cyan sticks: Binding site residues")
print("- Red sticks: Asp84 (salt bridge)")
