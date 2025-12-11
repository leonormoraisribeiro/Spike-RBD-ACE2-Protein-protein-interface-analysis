from Bio.PDB import PDBParser, NeighborSearch

# load structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure('6m0j', '6m0j_fixed.pdb')

# separate chains A and E
atoms_A = [atom for atom in structure[0]['A'].get_atoms()]
atoms_E = [atom for atom in structure[0]['E'].get_atoms()] 

# create neighbor search system for chain E
ns = NeighborSearch(atoms_E)

interface_residues = set()
cutoff_distance = 8.0 # Distance we found in pymol
# For each atom in chain A, find neighbors in chain E within cutoff
for atom_a in atoms_A:
    neighbors = ns.search(atom_a.get_coord(), cutoff_distance)
    if neighbors:
        # If there are neighbors, the residue of A belongs to the interface
        interface_residues.add(atom_a.get_parent())
        # The neighbors found in E also belong
        for atom_e in neighbors:
            interface_residues.add(atom_e.get_parent())

print(f"Interface residues (<{cutoff_distance}A): {len(interface_residues)}")