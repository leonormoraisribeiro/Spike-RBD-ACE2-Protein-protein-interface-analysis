import sys
from Bio.PDB import PDBParser, NeighborSearch

# CONFIGURATION 
PDB_FILENAME = "6m0j_fixed.pdb"
OUTPUT_PML = "visualize_distance_interface.pml"
CUTOFF = 8.0  # Geometric cutoff in Angstroms

print(f"--- Calculating Distance Interface ({CUTOFF} A) ---")

# 1. Load Structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure('struct', PDB_FILENAME)
model = structure[0]

# 2. Select Chains (Handle E vs B for RBD)
atoms_A = list(model['A'].get_atoms())

chain_E_id = 'E'
if 'E' in model:
    atoms_E = list(model['E'].get_atoms())
elif 'B' in model:
    chain_E_id = 'B'
    atoms_E = list(model['B'].get_atoms())
else:
    print("Error: RBD chain (E or B) not found.")
    sys.exit(1)

# 3. Neighbor Search
# Create a search tree for chain E atoms
ns = NeighborSearch(atoms_E)
interface_residues = set()

# Check every atom in A against atoms in E
for atom_a in atoms_A:
    neighbors = ns.search(atom_a.get_coord(), CUTOFF)
    if neighbors:
        # Add residue from Chain A
        interface_residues.add(('A', atom_a.get_parent().id[1]))
        
        # Add neighbor residues from Chain E
        for atom_e in neighbors:
            interface_residues.add((chain_E_id, atom_e.get_parent().id[1]))

print(f"Total residues found (Distance): {len(interface_residues)}")

# 4. Generate PyMOL Script
# Organize and sort residues for the selection string
resis_A = sorted(list(set([str(r[1]) for r in interface_residues if r[0] == 'A'])))
resis_E = sorted(list(set([str(r[1]) for r in interface_residues if r[0] == chain_E_id])))

# Create selection strings (e.g., "10+12+15")
sel_string_A = "+".join(resis_A)
sel_string_E = "+".join(resis_E)

print(f"Generating '{OUTPUT_PML}'...")

with open(OUTPUT_PML, "w") as f:
    # Initial setup
    f.write(f"load {PDB_FILENAME}\n")
    f.write("bg_color white\n")
    f.write("hide all\n")
    
    # Base Visualization (Transparent Cartoon)
    f.write("show cartoon\n")
    f.write("color gray90\n")
    f.write("set transparency, 0.5\n")
    
    # Highlight Chain A Interface (Blue-ish)
    if sel_string_A:
        f.write(f"select dist_A, chain A and resi {sel_string_A}\n")
        f.write("show sticks, dist_A\n")
        f.write("color slate, dist_A\n")
    
    # Highlight Chain E Interface (Red-ish)
    if sel_string_E:
        f.write(f"select dist_E, chain {chain_E_id} and resi {sel_string_E}\n")
        f.write("show sticks, dist_E\n")
        f.write("color raspberry, dist_E\n")

    # Camera settings
    f.write("deselect\n")
    f.write("zoom dist_A or dist_E, 5\n")
    f.write(f"echo Geometric Interface (Dist < {CUTOFF} A)\n")
