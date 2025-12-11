import math
import sys

# load forcefield parameters
VDW_PRM_FILE = "vdwprm.txt" 
PDB_FILE = "6m0j_fixed.pdbqt" 
ASA_COMPLEX = "6m0j_fixed.asa"
ASA_CHAIN_A = "A.asa" 
ASA_CHAIN_E = "B.asa" 
ASA_THRESHOLD = 0.01  # Threshold 


class VdwParamset():
    def __init__(self, filename):
        self.at_types = {}
        try:
            with open(filename, 'r') as f:
                for line in f:
                    if line.startswith('#') or len(line) < 10: continue
                    parts = line.split()
                    # Formato: ID Eps Sig ... fsrf
                    self.at_types[parts[0]] = {
                        'eps': float(parts[1]), 
                        'sig': float(parts[2]), 
                        'fsrf': float(parts[4])
                    }
        except:
            print(f"ERROR: {filename} not found."); sys.exit(1)

def get_residue_asa(asa_file):
    """Reads .asa and sums by residue to define the interface"""
    res_asa = {}
    try:
        with open(asa_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    chain = line[21]; res = int(line[22:26])
                    res_asa[(chain, res)] = res_asa.get((chain, res), 0.0) + float(line[54:62])
    except: pass
    return res_asa

def read_atomic_asa(filename):
    """Reads .asa atom by atom for energy calculation"""
    data = {}
    try:
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    name = line[12:16].strip()
                    res = int(line[22:26])
                    chain = line[21]
                    asa = float(line[54:62])
                    data[(chain, res, name)] = asa
    except: pass
    return data

def guess_atom_type(atom_name, pdbqt_type):
    """Uses the PDBQT type if it exists in the forcefield, otherwise guesses by name"""
    if pdbqt_type and pdbqt_type in ff.at_types:
        return pdbqt_type
        
    n = atom_name.strip().upper()
    if n.startswith('C'): return 'C'
    if n.startswith('N'): return 'N'
    if n.startswith('O'): return 'OA'
    if n.startswith('S'): return 'SA'
    return 'C'


print(f"Final energy calculation for interface residues...")

# A. Define Interface
rsa_c = get_residue_asa(ASA_COMPLEX)
rsa_a = get_residue_asa(ASA_CHAIN_A)
rsa_e = get_residue_asa(ASA_CHAIN_E)
interface_list = []

for (chain, res), bound in rsa_c.items():
    free = 0.0
    if chain == 'A': free = rsa_a.get((chain, res), 0.0)
    else: free = rsa_e.get((chain, res), rsa_e.get(('E', res), 0.0))
    
    if (free - bound) > ASA_THRESHOLD:
        interface_list.append((chain, res))

INTERFACE_SET = set(interface_list)
print(f"1. Interface defined: {len(interface_list)} residues.")

# B. Load Data
ff = VdwParamset(VDW_PRM_FILE)
asa_atom_c = read_atomic_asa(ASA_COMPLEX)
asa_atom_a = read_atomic_asa(ASA_CHAIN_A)
asa_atom_e = read_atomic_asa(ASA_CHAIN_E)

# C. Read PDBQT and Calculate (Single efficient loop)
atoms_A = []
atoms_E = []
E_solv = 0.0

with open(PDB_FILE, 'r') as f:
    for line in f:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            chain = line[21]; res = int(line[22:26])
            
            # Filter: Only process interface atoms
            if (chain, res) not in INTERFACE_SET: continue
            
            name = line[12:16].strip()
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            
            # READ CHARGE FROM PDBQT 
            q = 0.0
            if len(line) > 70:
                try: q = float(line[70:76])
                except: pass
            
            # Read Type
            pdbqt_type = line[77:].strip() if len(line) > 77 else ""
            atype = guess_atom_type(name, pdbqt_type)
            
            # Parameters
            p = ff.at_types.get(atype, ff.at_types['C'])
            
            # 1. Solvation Energy
            bound = asa_atom_c.get((chain, res, name), 0.0)
            free = asa_atom_a.get((chain, res, name), 0.0) if chain=='A' else asa_atom_e.get((chain, res, name), 0.0)
            E_solv += p['fsrf'] * (bound - free)
            
            # Save atom for pair interactions
            atom_data = {'x':x, 'y':y, 'z':z, 'q':q, 'eps':p['eps'], 'sig':p['sig']}
            if chain == 'A': atoms_A.append(atom_data)
            else: atoms_E.append(atom_data)

print(f"2. Atoms loaded successfully (A: {len(atoms_A)}, E: {len(atoms_E)})")
print("3. Calculating interactions ...")

# D. Pair Interactions (VdW + Elec)
E_vdw = 0.0
E_elec = 0.0

for a1 in atoms_A:
    for a2 in atoms_E:
        dx, dy, dz = a1['x']-a2['x'], a1['y']-a2['y'], a1['z']-a2['z']
        r2 = dx*dx + dy*dy + dz*dz
        
        if r2 > 144.0: continue # Cutoff 12A
        r = math.sqrt(r2)
        if r < 0.1: continue

        # VdW (Lennard-Jones)
        sig = (a1['sig'] + a2['sig']) * 0.5
        eps = math.sqrt(a1['eps'] * a2['eps'])
        E_vdw += 4.0 * eps * ((sig/r)**12 - (sig/r)**6)
        
        # Electrostatics (Mehler-Solmajer)
        if abs(a1['q']) > 0.0001 and abs(a2['q']) > 0.0001:
            denom = 1 - 7.7839 * math.exp(-0.3153 * r)
            if denom < 0.01: denom = 0.01
            eps_r = (86.9525 / denom) - 8.5525
            if eps_r < 1.0: eps_r = 1.0
            
            E_elec += 332.16 * a1['q'] * a2['q'] / (eps_r * r)


print(f"FINAL RESULTS (Interface with Real Charges)")
print(f"Van der Waals:      {E_vdw:10.4f} kcal/mol")
print(f"Electrostatics:      {E_elec:10.4f} kcal/mol")
print(f"Solvation:         {E_solv:10.4f} kcal/mol")
print(f"DELTA G TOTAL:      {E_vdw + E_elec + E_solv:10.4f} kcal/mol")
