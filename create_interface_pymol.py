import os

# interface_asa_variation.py
def read_rsa(filename):
    asa = {}
    try:
        with open(filename) as f:
            for line in f:
                if line.startswith("RES"):
                    parts = line.split()
                    resname = parts[1]
                    chain = parts[2]
                    resnum = int(parts[3])
                    asa_value = float(parts[4])
                    asa[(chain, resnum, resname)] = asa_value
    except FileNotFoundError:
        print(f"Erro: O ficheiro {filename} não foi encontrado.")
        return {}
    return asa

# load data
asa_bound = read_rsa("6m0j_fixed.rsa")
asa_A = read_rsa("A.rsa") 
asa_B = read_rsa("B.rsa") 

interface_residues = []

# compare ASA values to find interface residues
for (chain, resnum, resname), bound in asa_bound.items():
    if chain == 'A':
        free = asa_A.get((chain, resnum, resname), 0)
    elif chain == 'E':
        free = asa_B.get((chain, resnum, resname), 0)
    else:
        continue
    
    delta = free - bound
    if delta > 1.0: # Interface threshold
        interface_residues.append((chain, resnum))

# create PyMol script to visualize interface residues

# Separate residue numbers by chain
resis_A = [str(r[1]) for r in interface_residues if r[0] == 'A']
resis_E = [str(r[1]) for r in interface_residues if r[0] == 'E']

# Create selection string (e.g., 19+24+27...)
sel_string_A = "+".join(resis_A)
sel_string_E = "+".join(resis_E)

pml_filename = "visualize_interface.pml"
pdb_filename = "6m0j_fixed.pdb"

print(f"Creating PyMol file: {pml_filename}")
print(f"Interface Residues A: {len(resis_A)}")
print(f"Interface Residues E: {len(resis_E)}")

with open(pml_filename, "w") as f:
    # 1. Load structure
    f.write(f"load {pdb_filename}\n")
    f.write("bg_color white\n")
    f.write("hide all\n")
    
    # 2. Base visualization (Light gray Cartoon)
    f.write("show cartoon\n")
    f.write("color gray80\n")
    
    # 3. Select and color Interface of Chain A (ACE2 - e.g., Cyan)
    if sel_string_A:
        f.write(f"select interface_A, chain A and resi {sel_string_A}\n")
        f.write("show sticks, interface_A\n")
        f.write("color cyan, interface_A\n")
        # Color carbon atoms differently for contrast
        f.write("color deepteal, interface_A and name C*\n") 

    # 4. Select and color Interface of Chain E (RBD - e.g., Orange)
    if sel_string_E:
        f.write(f"select interface_E, chain E and resi {sel_string_E}\n")
        f.write("show sticks, interface_E\n")
        f.write("color orange, interface_E\n")
        f.write("color brightorange, interface_E and name C*\n")

    # 5. Zoom and final details
    f.write("deselect\n")
    f.write("zoom interface_A or interface_E, 5\n")
    f.write("set transparency, 0.4\n")