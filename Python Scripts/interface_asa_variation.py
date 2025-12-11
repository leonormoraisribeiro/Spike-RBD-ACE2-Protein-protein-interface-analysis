# read RSA files that we created with NACCESS
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

# compare ASA values to find interface residues
for (chain, resnum, resname), bound in asa_bound.items():
    if chain == 'A':
        free = asa_A.get((chain, resnum, resname), 0)
    elif chain == 'E':
        free = asa_B.get((chain, resnum, resname), 0)
    else:
        continue
    # compute delta
    delta = free - bound
    if delta > 0.01:
        interface.append((chain, resname, resnum, delta))

print("Interface residues:")
i = 0
for c, rname, rnum, d in interface:
    print(f"{c} {rname}{rnum} ΔASA={d:.2f}")
    i += 1
print(f"Total interface residues: {i}")

with open("interface_data.py", "w") as f:
    f.write("#This file was automatically generated\n")
    f.write("# Format: ('Chain', ResNum)\n")
    f.write("INTERFACE_LIST = [\n")
    
    interface.sort(key=lambda x: (x[0], x[2]))
    
    for c, rname, rnum, d in interface:
        # Add the rest as a comment to help reading
        f.write(f"    ('{c}', {rnum}),  # {rname} dASA={d:.2f}\n")
    
    f.write("]\n")

print("File 'interface_data.py' created!")