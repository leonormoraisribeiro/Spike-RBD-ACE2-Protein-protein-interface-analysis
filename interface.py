def read_rsa(filename):
    asa = {}
    with open(filename) as f:
        for line in f:
            if line.startswith("RES"):
                parts = line.split()
                resname = parts[1]
                chain = parts[2]
                resnum = int(parts[3])
                asa_value = float(parts[4])
                asa[(chain, resnum, resname)] = asa_value
    return asa


asa_bound = read_rsa("6m0j_fixed.rsa")
asa_A = read_rsa("A.rsa")
asa_B = read_rsa("B.rsa")

interface = []

for (chain, resnum, resname), bound in asa_bound.items():
    if chain == 'A':
        free = asa_A.get((chain, resnum, resname), 0)
    elif chain == 'B':
        free = asa_B.get((chain, resnum, resname), 0)
    else:
        continue
    
    delta = free - bound
    if delta > 1.0:
        interface.append((chain, resname, resnum, delta))

print("Interface residues:")
for c, rname, rnum, d in interface:
    print(f"{c} {rname}{rnum} ΔASA={d:.2f}")
