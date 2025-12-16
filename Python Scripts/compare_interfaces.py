import sys
from Bio.PDB import PDBParser, NeighborSearch

# Distance-based interface detection
def get_distance_interface(pdb_file, cutoff=8.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', pdb_file)
    
    # Get atoms for chain A and try to find E or B
    atoms_A = list(structure[0]['A'].get_atoms())
    try:
        atoms_E = list(structure[0]['E'].get_atoms())
        chain_id_E = 'E'
    except KeyError:
        atoms_E = list(structure[0]['B'].get_atoms())
        chain_id_E = 'B'

    # Search neighbors
    ns = NeighborSearch(atoms_E)
    interface_set = set()

    for atom_a in atoms_A:
        if ns.search(atom_a.get_coord(), cutoff):
            # Add residue from A
            interface_set.add(('A', atom_a.get_parent().id[1]))
            # Add neighbors from E
            for atom_e in ns.search(atom_a.get_coord(), cutoff):
                interface_set.add((chain_id_E, atom_e.get_parent().id[1]))
                
    return interface_set

# Helper to read NACCESS RSA files
def read_rsa(filename):
    data = {}
    try:
        with open(filename) as f:
            for line in f:
                if line.startswith("RES"):
                    parts = line.split()
                    # Store as (Chain, ResNum) -> ASA
                    data[(parts[2], int(parts[3]))] = float(parts[4])
    except FileNotFoundError:
        pass
    return data

# Energy/ASA-based interface detection
def get_asa_interface(complex_file, free_a, free_b, threshold=0.01):
    bound = read_rsa(complex_file)
    f_a = read_rsa(free_a)
    f_b = read_rsa(free_b)
    
    interface_set = set()
    
    for (chain, res), val_bound in bound.items():
        val_free = 0.0
        
        if chain == 'A':
            val_free = f_a.get((chain, res), 0.0)
        else:
            # Check B or E chains
            val_free = f_b.get((chain, res), f_b.get(('B', res), f_b.get(('E', res), 0.0)))

        if (val_free - val_bound) > threshold:
            interface_set.add((chain, res))
            
    return interface_set

# Main execution
print("Interface Correlation Analysis")

# Calculate sets
set_dist = get_distance_interface("6m0j_fixed.pdb", cutoff=8.0)
set_asa = get_asa_interface("6m0j_fixed.rsa", "A.rsa", "B.rsa", threshold=0.01)

# Set operations
intersection = set_dist.intersection(set_asa)
only_dist = set_dist - set_asa
only_asa = set_asa - set_dist
union = set_dist.union(set_asa)

# Print stats
print(f"Distance-based residues: {len(set_dist)}")
print(f"Energy-based residues:   {len(set_asa)}")
print(f"Intersection:            {len(intersection)}")

# Calculate Jaccard Index
if union:
    jaccard = len(intersection) / len(union)
    print(f"Jaccard Index:           {jaccard:.4f}")

# List discrepancies
if only_dist:
    print("\nUnique to Distance (Geometric only):")
    print(sorted(list(only_dist)))
    
if only_asa:
    print("\nUnique to Energy (Buried but distant):")
    print(sorted(list(only_asa)))