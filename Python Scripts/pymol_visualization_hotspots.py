import os

# Output filename for the PyMOL script
pml_filename = "visualize_hotspots.pml"

# PyMOL commands content
pml_content = """
# STEP 4: HOTSPOT INTERACTION VISUALIZATION

# 1. LOAD PRECESSED STRUCTURE
# Use the local fixed file to ensure consistency with energy calculations
reinitialize
load 6m0j_fixed.pdb

# 2. GENERAL SCENE SETUP
bg_color white
hide everything
# Define chains for easier selection
select RBD, chain E
select ACE2, chain A

# 3. SELECT KEY RESIDUES (HOTSPOTS)
# Based on Alanine Scanning results (Step 3)
# Phe486: Hydrophobic anchor
# Tyr505: H-bond stabilizer
select hotspots, chain E and resi 486+505

# 4. SELECT INTERACTION PARTNERS
# Find residues on ACE2 within 4.0 Angstroms of the hotspots
select neighbors, chain A within 4.0 of hotspots

# 5. VISUAL REPRESENTATION
# Context: Transparent cartoon
show cartoon
set cartoon_transparency, 0.7
color grey90, ACE2
color pale_cyan, RBD

# Details: Show side chains as sticks
show sticks, hotspots
show sticks, neighbors

# 6. COLORING SCHEME
# Hotspots in Red (to highlight importance/destabilization)
color red, hotspots
# Neighbors in Marine Blue (for contrast)
color marine, neighbors

# Apply elemental coloring (N=Blue, O=Red) to sticks
# This is crucial to visualize chemical interactions
util.cnc

# 7. VISUALIZE INTERACTIONS (H-BONDS)
# Calculate polar contacts automatically
dist interactions, hotspots, neighbors, mode=2
# Style the lines
set dash_gap, 0.2
set dash_width, 2.5
set dash_color, black
hide labels, interactions

# 8. LABELS AND ANNOTATIONS
# Label the hotspots so the reader knows which residues they are
label n. CA and hotspots, '%s-%s' % (resn, resi)
set label_color, black
set label_size, 20
set label_position, [2, 2, 2]

# 9. RENDERING SETTINGS (PUBLICATION QUALITY)
zoom hotspots, 8
orient hotspots

# Remove shadows for clearer structural details
set ray_shadows, 0
# Orthoscopic view removes perspective distortion (better for analysis)
set orthoscopic, on
# High quality smoothing
set antialias, 2

"""

# Write the file
with open(pml_filename, "w") as f:
    f.write(pml_content)

print(f"PyMOL script '{pml_filename}' created successfully.")
