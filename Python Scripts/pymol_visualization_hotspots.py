import os

# Content of the PyMOL script (PML)
pml_content = """
# 1. LOAD STRUCTURE
reinitialize
fetch 6m0j, type=pdb
# If you prefer to use your local fixed file, uncomment the line below:
# load 6m0j_fixed.pdbqt

# 2. GENERAL SETUP
bg_color white
hide everything
select RBD, chain E
select ACE2, chain A

# 3. SELECT HOTSPOTS (Phe486 and Tyr505)
# These were identified as the top contributors in the plot
select hotspots, chain E and resi 486+505

# 4. SELECT NEIGHBORS ON ACE2 (Interaction partners)
# Select residues on chain A within 5.0 Angstroms of the hotspots
select neighbors, chain A within 5.0 of hotspots

# 5. REPRESENTATION (Visual Style)
# Show transparent cartoon for context
show cartoon
set cartoon_transparency, 0.6
color grey80, ACE2
color pale_cyan, RBD

# Show Side Chains (Sticks)
show sticks, hotspots
show sticks, neighbors

# 6. COLORING AND HIGHLIGHTS
# Color Hotspots Red (to match the bar chart)
color red, hotspots
# Color ACE2 Neighbors Marine Blue
color marine, neighbors

# Color by element (crucial to see Oxygens/Nitrogens)
# Keeps the base carbon color but paints O red and N blue
util.cnc

# 7. VISUALIZE INTERACTIONS
# Find polar contacts (Hydrogen Bonds)
dist interactions, hotspots, neighbors, mode=2
set dash_gap, 0.2
set dash_width, 2.0
color black, interactions

# 8. CAMERA AND REFINEMENT
# Smooth zoom onto the residues of interest
zoom hotspots, 8
orient hotspots

# Improve visual quality (Shadows and lighting for the report)
set ray_shadows, 0
set specular, 0.5
set antialias, 2
set orthoscopic, on  # Removes perspective distortion (better for scientific figures)

# Message for the user
print("Visualization ready! Press 'Ray' (top right) to render a high-quality image.")
"""

# Output filename
filename = "visualize_hotspots.pml"

# Write the file
with open(filename, "w") as f:
    f.write(pml_content)

print(f"The file '{filename}' has been created.")
