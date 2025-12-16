import sys
from Bio.PDB import PDBParser, NeighborSearch

# --- CONFIGURATION ---
PDB_FILENAME = "6m0j_fixed.pdb"
RSA_COMPLEX = "6m0j_fixed.rsa"
RSA_A = "A.rsa"
RSA_B = "B.rsa"
OUTPUT_PML = "visualize_comparison.pml"

DIST_CUTOFF = 8.0
ASA_THRESHOLD = 0.01

print("--- Generating Comparison Visualization ---")

# --- 1. GET DISTANCE INTERFACE ---
def get_distance_interface():
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', PDB_FILENAME)
    
    atoms_A = list(structure[0]['A'].get_atoms())
    # Handle Chain E or B
    try:
        atoms_E = list(structure[0]['E'].get_atoms())
        chain_E_id = 'E'
    except:
        atoms_E = list(structure[0]['B'].get_atoms())
        chain_E_id = 'B'

    ns = NeighborSearch(atoms_E)
    dist_set = set()

    for atom_a in atoms_A:
        if ns.search(atom_a.get_coord(), DIST_CUTOFF):
            dist_set.add(('A', atom_a.get_parent().id[1]))
            for atom_e in ns.search(atom_a.get_coord(), DIST_CUTOFF):
                dist_set.add((chain_E_id, atom_e.get_parent().id[1]))
    return dist_set, chain_E_id

# --- 2. GET ASA INTERFACE ---
def read_rsa(filename):
    data = {}
    try:
        with open(filename) as f:
            for line in f:
                if line.startswith("RES"):
                    p = line.split()
                    data[(p[2], int(p[3]))] = float(p[4])
    except: pass
    return data

def get_asa_interface(chain_e_id):
    bound = read_rsa(RSA_COMPLEX)
    f_a = read_rsa(RSA_A)
    f_b = read_rsa(RSA_B)
    asa_set = set()
    
    for (chain, res), val_bound in bound.items():
        val_free = 0.0
        if chain == 'A': val_free = f_a.get((chain, res), 0.0)
        else: val_free = f_b.get((chain, res), f_b.get(('E', res), f_b.get(('B', res), 0.0)))
        
        if (val_free - val_bound) > ASA_THRESHOLD:
            asa_set.add((chain, res))
    return asa_set

# --- 3. CALCULATE SETS ---
set_dist, chain_e_id = get_distance_interface()
set_asa = get_asa_interface(chain_e_id)

# The Intersection (Core Interface)
common = set_dist.intersection(set_asa)
# Only in Distance (Periphery/Bystanders)
only_dist = set_dist - set_asa
# Only in ASA (Buried but distant - usually empty)
only_asa = set_asa - set_dist

print(f"Core residues (Both): {len(common)}")
print(f"Peripheral residues (Dist only): {len(only_dist)}")

# --- 4. GENERATE PYMOL SCRIPT ---
def make_sel_string(residue_set, target_chain):
    resis = sorted(list(set([str(r[1]) for r in residue_set if r[0] == target_chain])))
    return "+".join(resis)

with open(OUTPUT_PML, "w") as f:
    f.write(f"load {PDB_FILENAME}\n")
    f.write("bg_color white\n")
    f.write("hide all\n")
    
    # Context (Ghostly Cartoon)
    f.write("show cartoon\n")
    f.write("color gray90\n")
    f.write("set transparency, 0.6\n") # Very transparent

    # GROUP 1: CORE INTERFACE (MAGENTA)
    # These are the residues that are close AND bury surface
    sA = make_sel_string(common, 'A')
    sE = make_sel_string(common, chain_e_id)
    if sA or sE:
        sel = []
        if sA: sel.append(f"(chain A and resi {sA})")
        if sE: sel.append(f"(chain {chain_e_id} and resi {sE})")
        f.write(f"select core_interface, {' or '.join(sel)}\n")
        f.write("show sticks, core_interface\n")
        f.write("color magenta, core_interface\n") # Strong color
        f.write("color hotpink, core_interface and name C*\n")
    
    # GROUP 2: PERIPHERAL/DISTANCE ONLY (CYAN)
    # These are close but don't interact strongly (water mediated, etc)
    sA = make_sel_string(only_dist, 'A')
    sE = make_sel_string(only_dist, chain_e_id)
    if sA or sE:
        sel = []
        if sA: sel.append(f"(chain A and resi {sA})")
        if sE: sel.append(f"(chain {chain_e_id} and resi {sE})")
        f.write(f"select periphery_dist, {' or '.join(sel)}\n")
        f.write("show lines, periphery_dist\n") # Lines are thinner than sticks
        f.write("color lightblue, periphery_dist\n")

    # GROUP 3: ENERGY ONLY (ORANGE) - Usually empty
    sA = make_sel_string(only_asa, 'A')
    sE = make_sel_string(only_asa, chain_e_id)
    if sA or sE:
        sel = []
        if sA: sel.append(f"(chain A and resi {sA})")
        if sE: sel.append(f"(chain {chain_e_id} and resi {sE})")
        f.write(f"select hidden_energy, {' or '.join(sel)}\n")
        f.write("show sticks, hidden_energy\n")
        f.write("color orange, hidden_energy\n")
        f.write("set sphere_scale, 0.3\n")
        f.write("show spheres, hidden_energy\n") # Highlight these if they exist

    # Final Polish
    f.write("deselect\n")
    f.write("zoom core_interface, 5\n")
    f.write("set ray_shadows, 0\n")
    f.write("set stick_radius, 0.25\n")
    
    # Print legend to PyMol console
    f.write("echo [LEGEND] Magenta: Core Interface (Energy + Dist)\n")
    f.write("echo [LEGEND] LightBlue: Periphery (Dist Only)\n")

print(f"Created '{OUTPUT_PML}'. Open in PyMol to see the layers.")