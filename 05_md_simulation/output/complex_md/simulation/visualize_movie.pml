# PyMOL Visualization Script for QTY μOR + Fentanyl MD Trajectory
# Load structure first, then trajectory

# Load the complex from the system PDB (starting structure)
load /Users/aneeshbondugula/Downloads/QTY_docking/05_md_simulation/output/complex_md/system/complex.pdb, complex

# Load the DCD trajectory
load_traj /Users/aneeshbondugula/Downloads/QTY_docking/05_md_simulation/output/complex_md/simulation/trajectory.dcd, complex

# Hide everything first
hide everything

# Show protein as cartoon
select protein, polymer
show cartoon, protein
color marine, protein

# Show ligand (fentanyl) - residue name LIG
select ligand, resn LIG
show sticks, ligand
color orange, ligand
show spheres, ligand
set sphere_scale, 0.3, ligand

# Show binding site residues
select binding_site, (polymer within 5 of ligand) and not name H*
show sticks, binding_site
color cyan, binding_site

# Highlight Asp84 (key salt bridge residue)
select asp84, resi 84 and resn ASP
show sticks, asp84
color red, asp84

# Hide waters for cleaner view
hide everything, resn WAT or resn HOH or resn Na+ or resn Cl-

# Center on ligand
center ligand
zoom ligand, 20

# Nice rendering settings
set cartoon_fancy_helices, 1
bg_color white
set ray_shadow, 0

# Movie settings
set movie_fps, 10
set cache_frames, 1

print "=== MD TRAJECTORY LOADED ==="
print "Frames loaded: " + str(cmd.count_frames())
print ""
print "Controls:"
print "  - Press PLAY button or use Movie menu"
print "  - Drag slider to scrub through frames"
print "  - Press Home to go to first frame"
print "  - Press End to go to last frame"
