#!/usr/bin/env python3
import argparse
from Bio.PDB import PDBParser, NeighborSearch, Selection
import numpy as np
import csv
import os

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--pdb", required=True, help="PDB file path")
    p.add_argument("--chainA", default="A", help="Chain A id (e.g. ACE2)")
    p.add_argument("--chainB", default="E", help="Chain B id (e.g. RBD)")
    p.add_argument("--cutoff", type=float, default=5.5, help="cut-off distance (Å)")
    p.add_argument("--out", default="interface.csv", help="output CSV")
    p.add_argument("--pymol", default="pymol_selections.pml", help="pymol pml output")
    return p.parse_args()

def collect_residue_atoms(chain):
    residues = []
    for res in chain:
        # skip hetatoms and waters
        hetflag = res.id[0]
        if hetflag != " ":
            continue
        atoms = [a for a in res.get_atoms()]
        if len(atoms) == 0:
            continue
        residues.append((res, atoms))
    return residues

def main():
    args = parse_args()
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", args.pdb)

    model = structure[0]
    chainA = model[args.chainA]
    chainB = model[args.chainB]

    res_atoms_A = collect_residue_atoms(chainA)
    res_atoms_B = collect_residue_atoms(chainB)

    # Build global atom list for neighbor search
    all_atoms_B = [atom for _, atoms in res_atoms_B for atom in atoms]
    ns_B = NeighborSearch(all_atoms_B)

    results = []
    for resA, atomsA in res_atoms_A:
        min_dist = np.inf
        contact_count = 0
        for a in atomsA:
            close = ns_B.search(a.coord, args.cutoff)  # returns atoms within cutoff
            if close:
                # compute true distances to get min
                dists = [np.linalg.norm(a.coord - b.coord) for b in close]
                local_min = min(dists)
                if local_min < min_dist:
                    min_dist = local_min
                contact_count += len(close)
        if min_dist == np.inf:
            min_dist = None
        results.append({
            "chain": args.chainA,
            "resnum": resA.id[1],
            "insertion": resA.id[2].strip() if resA.id[2] != ' ' else '',
            "resname": resA.get_resname(),
            "min_dist": min_dist,
            "contacts": contact_count
        })

    # Repeat for chain B (against A)
    all_atoms_A = [atom for _, atoms in res_atoms_A for atom in atoms]
    ns_A = NeighborSearch(all_atoms_A)
    for resB, atomsB in res_atoms_B:
        min_dist = np.inf
        contact_count = 0
        for b in atomsB:
            close = ns_A.search(b.coord, args.cutoff)
            if close:
                dists = [np.linalg.norm(b.coord - a.coord) for a in close]
                local_min = min(dists)
                if local_min < min_dist:
                    min_dist = local_min
                contact_count += len(close)
        if min_dist == np.inf:
            min_dist = None
        results.append({
            "chain": args.chainB,
            "resnum": resB.id[1],
            "insertion": resB.id[2].strip() if resB.id[2] != ' ' else '',
            "resname": resB.get_resname(),
            "min_dist": min_dist,
            "contacts": contact_count
        })

    # Write CSV
    with open(args.out, "w", newline='') as csvf:
        writer = csv.DictWriter(csvf, fieldnames=["chain","resnum","insertion","resname","min_dist","contacts"])
        writer.writeheader()
        for r in sorted(results, key=lambda x: (x["chain"], int(x["resnum"]))):
            writer.writerow(r)

    print(f"Wrote {args.out} with {len(results)} residues.")

    # Create PyMOL selection script to highlight interface residues (those with min_dist not None)
    with open(args.pymol, "w") as pm:
        pm.write("# PyMOL script to highlight interface residues\n")
        pm.write("hide everything\nshow cartoon\n")
        pm.write(f"load {os.path.abspath(args.pdb)}\n")
        pm.write("bg_color white\n")
        pm.write("color grey70, all\n")
        pm.write("\n")
        # colour interface residues
        for r in results:
            if r["min_dist"] is not None:
                sel = f"/{structure.id}//{r['chain']}/{r['resnum']}"
                label = f"{r['chain']}{r['resnum']}{r['insertion']}"
                pm.write(f"select interface_{r['chain']}_{r['resnum']}, chain {r['chain']} and resi {r['resnum']}\n")
        pm.write("\n")
        pm.write("color red, interface_*\n")
        pm.write("show sticks, interface_*\n")
        pm.write("zoom interface_*\n")
    print(f"Wrote PyMOL script {args.pymol}")

if __name__ == "__main__":
    main()
