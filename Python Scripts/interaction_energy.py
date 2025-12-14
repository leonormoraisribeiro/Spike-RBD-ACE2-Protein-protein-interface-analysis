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
                    if line.startswith('#') or len(line) < 10:
                        continue
                    parts = line.split()
                    self.at_types[parts[0]] = {
                        'eps': float(parts[1]),
                        'sig': float(parts[2]),
                        'fsrf': float(parts[4])
                    }
        except:
            print(f"ERROR: {filename} not found.")
            sys.exit(1)


def get_residue_asa(asa_file):
    """Reads .asa and sums by residue to define the interface"""
    res_asa = {}
    try:
        with open(asa_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    chain = line[21]
                    res = int(line[22:26])
                    res_asa[(chain, res)] = res_asa.get((chain, res), 0.0) + float(line[54:62])
    except:
        pass
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
    except:
        pass
    return data


def guess_atom_type(atom_name, pdbqt_type, ff):
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

#   Original Function — COMPUTES WT ΔG
def compute_interaction_energy(pdbqt_file,
                               asa_complex="6m0j_fixed.asa",
                               asa_A="A.asa",
                               asa_E="B.asa",
                               return_components=True,
                               verbose=False):

    rsa_c = get_residue_asa(asa_complex)
    rsa_a = get_residue_asa(asa_A)
    rsa_e = get_residue_asa(asa_E)

    INTERFACE = set()

    for (chain, res), bound in rsa_c.items():
        free = rsa_a.get((chain, res), 0.0) if chain == "A" else rsa_e.get((chain, res), 0.0)
        if (free - bound) > ASA_THRESHOLD:
            INTERFACE.add((chain, res))

    if verbose:
        print(f"Interface residues: {len(INTERFACE)}")

    ff = VdwParamset(VDW_PRM_FILE)
    asa_atom_c = read_atomic_asa(asa_complex)
    asa_atom_a = read_atomic_asa(asa_A)
    asa_atom_e = read_atomic_asa(asa_E)

    atoms_A = []
    atoms_E = []
    E_solv = 0.0

    with open(pdbqt_file, "r") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            chain = line[21]
            res = int(line[22:26])

            if (chain, res) not in INTERFACE:
                continue

            name = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            q = float(line[70:76]) if len(line) > 76 else 0.0

            pdbqt_type = line[77:].strip() if len(line) > 77 else ""
            atype = guess_atom_type(name, pdbqt_type, ff)
            p = ff.at_types.get(atype, ff.at_types["C"])

            bound = asa_atom_c.get((chain, res, name), 0.0)
            free = asa_atom_a.get((chain, res, name), 0.0) if chain == 'A' else asa_atom_e.get((chain, res, name), 0.0)
            E_solv += p["fsrf"] * (bound - free)

            atom = {"x": x, "y": y, "z": z, "q": q, "eps": p["eps"], "sig": p["sig"]}

            if chain == "A":
                atoms_A.append(atom)
            else:
                atoms_E.append(atom)

    E_vdw = 0.0
    E_elec = 0.0

    for a1 in atoms_A:
        for a2 in atoms_E:
            dx = a1["x"] - a2["x"]
            dy = a1["y"] - a2["y"]
            dz = a1["z"] - a2["z"]
            r2 = dx*dx + dy*dy + dz*dz

            if r2 < 0.1 or r2 > 144:
                continue

            r = math.sqrt(r2)

            sig = 0.5 * (a1["sig"] + a2["sig"])
            eps = math.sqrt(a1["eps"] * a2["eps"])
            E_vdw += 4 * eps * ((sig/r)**12 - (sig/r)**6)

            if abs(a1["q"]) > 1e-4 and abs(a2["q"]) > 1e-4:
                denom = 1 - 7.7839 * math.exp(-0.3153*r)
                denom = max(denom, 0.01)
                eps_r = max((86.9525 / denom) - 8.5525, 1.0)
                E_elec += 332.16 * a1["q"] * a2["q"] / (eps_r * r)

    total = E_vdw + E_elec + E_solv

    if return_components:
        return total, E_vdw, E_elec, E_solv
    else:
        return total


#   function for alanine scanning 
def compute_interaction_energy_with_ala(pdbqt_file,
                                        chain_mut,
                                        res_mut,
                                        asa_complex="6m0j_fixed.asa",
                                        asa_A="A.asa",
                                        asa_E="B.asa",
                                        return_components=True,
                                        verbose=False):

    allowed = {"N", "CA", "C", "O", "CB"}  # alanine atoms

    rsa_c = get_residue_asa(asa_complex)
    rsa_a = get_residue_asa(asa_A)
    rsa_e = get_residue_asa(asa_E)

    INTERFACE = set()

    for (ch, rn), bound in rsa_c.items():
        free = rsa_a.get((ch, rn), 0.0) if ch == "A" else rsa_e.get((ch, rn), 0.0)
        if (free - bound) > ASA_THRESHOLD:
            INTERFACE.add((ch, rn))

    ff = VdwParamset(VDW_PRM_FILE)
    asa_atom_c = read_atomic_asa(asa_complex)
    asa_atom_a = read_atomic_asa(asa_A)
    asa_atom_e = read_atomic_asa(asa_E)

    atoms_A = []
    atoms_E = []
    E_solv = 0.0

    with open(pdbqt_file, "r") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            ch = line[21]
            rn = int(line[22:26])
            atom = line[12:16].strip()

            if (ch == chain_mut) and (rn == res_mut):
                if atom not in allowed:
                    continue

            if (ch, rn) not in INTERFACE:
                continue

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            q = float(line[70:76]) if len(line) > 76 else 0.0

            pdbqt_type = line[77:].strip()
            atype = guess_atom_type(atom, pdbqt_type, ff)
            p = ff.at_types.get(atype, ff.at_types["C"])

            bound = asa_atom_c.get((ch, rn, atom), 0.0)
            free = asa_atom_a.get((ch, rn, atom), 0.0) if ch == "A" else asa_atom_e.get((ch, rn, atom), 0.0)
            E_solv += p["fsrf"] * (bound - free)

            data = {"x": x, "y": y, "z": z, "q": q, "eps": p["eps"], "sig": p["sig"]}

            if ch == "A":
                atoms_A.append(data)
            else:
                atoms_E.append(data)

    E_vdw = 0.0
    E_elec = 0.0

    for a1 in atoms_A:
        for a2 in atoms_E:
            dx = a1["x"] - a2["x"]
            dy = a1["y"] - a2["y"]
            dz = a1["z"] - a2["z"]
            r2 = dx*dx + dy*dy + dz*dz

            if r2 < 0.1 or r2 > 144:
                continue

            r = math.sqrt(r2)

            sig = 0.5*(a1["sig"]+a2["sig"])
            eps = math.sqrt(a1["eps"]*a2["eps"])
            E_vdw += 4*eps*((sig/r)**12 - (sig/r)**6)

            if abs(a1["q"]) > 1e-4 and abs(a2["q"]) > 1e-4:
                denom = 1 - 7.7839*math.exp(-0.3153*r)
                denom = max(denom, 0.01)
                eps_r = max((86.9525/denom) - 8.5525, 1)
                E_elec += 332.16*a1["q"]*a2["q"]/(eps_r*r)

    total = E_vdw + E_elec + E_solv
    if return_components:
        return total, E_vdw, E_elec, E_solv
    else:
        return total


# run wt as a test
if __name__ == "__main__":
    total, lj, elec, solv = compute_interaction_energy(
        PDB_FILE,
        asa_complex=ASA_COMPLEX,
        asa_A=ASA_CHAIN_A,
        asa_E=ASA_CHAIN_E,
        return_components=True,
        verbose=True
    )

    print("Interaction energy (WT)")
    print(f"Van der Waals:   {lj:10.4f} kcal/mol")
    print(f"Electrostatics:  {elec:10.4f} kcal/mol")
    print(f"Solvation:       {solv:10.4f} kcal/mol")
    print(f"TOTAL ΔG:        {total:10.4f} kcal/mol")
