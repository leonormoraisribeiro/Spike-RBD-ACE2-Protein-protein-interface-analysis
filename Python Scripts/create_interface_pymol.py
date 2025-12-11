import sys
import os

# import interface residue list
try:
    from interface_data import INTERFACE_LIST
    print(f" List loaded from'interface_data.py'.")
    print(f" Interface residues: {len(INTERFACE_LIST)}")
except ImportError:
    print(" File 'interface_data.py' not found.")
    print(" Please run the ASA analysis script first.")
    sys.exit(1)

pml_filename = "visualize_interface.pml"
pdb_filename = "6m0j_fixed.pdb" 


# chain A residues
resis_A = [str(r[1]) for r in INTERFACE_LIST if r[0] == 'A']

# chain E residues
resis_E = [str(r[1]) for r in INTERFACE_LIST if r[0] in ['E', 'B']]

# Create PyMol selection strings (e.g., "19+24+27")
sel_string_A = "+".join(resis_A)
sel_string_E = "+".join(resis_E)

print(f"   -> Chain A residues: {len(resis_A)}")
print(f"   -> Chain E residues: {len(resis_E)}")

# write the PyMol script

with open(pml_filename, "w") as f:

    f.write(f"load {pdb_filename}\n")
    f.write("bg_color white\n")
    f.write("hide all\n")
    
    # Base Visualization (Transparent Gray Cartoon)
    f.write("show cartoon\n")
    f.write("color gray80\n")
    f.write("set transparency, 0.4\n") 
    
    # Paint Chain A Interface (Cyan)
    if sel_string_A:
        f.write(f"select interface_A, chain A and resi {sel_string_A}\n")
        f.write("show sticks, interface_A\n")
        f.write("color cyan, interface_A\n")
        # Carbons with different color for contrast
        f.write("color deepteal, interface_A and name C*\n") 

    # Paint Chain E Interface (Orange)
    if sel_string_E:
        chain_id = INTERFACE_LIST[0][0] if resis_E else 'E' # Get the first entry's chain letter if exists
        # If the list has a mix, this filter resolves:
        target_chain = 'E' if 'E' in [r[0] for r in INTERFACE_LIST] else 'B'
        
        f.write(f"select interface_E, chain {target_chain} and resi {sel_string_E}\n")
        f.write("show sticks, interface_E\n")
        f.write("color orange, interface_E\n")
        f.write("color brightorange, interface_E and name C*\n")

    f.write("deselect\n")
    f.write("zoom interface_A or interface_E, 8\n") # Zoom with margin of 8 Angstroms
    f.write("set ray_shadows, 0\n") # Remove shadows for cleaner image

print(f"PyMol script '{pml_filename}' created")