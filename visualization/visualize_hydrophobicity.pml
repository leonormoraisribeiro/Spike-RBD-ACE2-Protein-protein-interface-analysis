
# HYDROPHOBIC COMPLEMENTARITY

reinitialize
load 6m0j_fixed.pdb

bg_color white
hide everything

# 1. DEFINE COLORS
# Color Hydrophobic residues (Oily) Orange
# Color Polar/Charged residues (Water-loving) White
# This helps visualize the 'Hydrophobic Effect' (-4.1 kcal/mol)
select hydrophobic, resn ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO+TYR
select polar, not hydrophobic

# 2. SURFACE REPRESENTATION
# Show the surface of ACE2 to see the "pockets"
show surface, chain A
color white, chain A
set transparency, 0.2, chain A

# Color the hydrophobic pockets on ACE2 surface
color pale_yellow, chain A and hydrophobic

# 3. SHOW THE RBD ANCHOR (Phe486)
# We want to see the Phe486 sticking into the surface
show sticks, chain E and resi 486
color orange, chain E and resi 486

# Show the rest of the RBD interface as a ribbon
show ribbon, chain E
color grey, chain E

# 4. HIGHLIGHT THE INTERACTION
# Select the pocket on ACE2 where Phe486 lands
select pocket, chain A within 5.0 of (chain E and resi 486)
color brightorange, pocket

# 5. LABELS
label n. CA and (chain E and resi 486), "Phe-486 (Anchor)"
set label_color, black
set label_position, [3,3,3]

# 6. RENDER SETTINGS
orient chain E and resi 486
zoom chain E and resi 486, 12
set ray_shadows, 0
set antialias, 2
set orthoscopic, on

