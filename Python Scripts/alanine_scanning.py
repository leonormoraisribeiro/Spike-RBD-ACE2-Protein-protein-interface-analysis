from interaction_energy_new import (
    compute_interaction_energy,
    compute_interaction_energy_with_ala
)

# WT PDBQT
PDBQT_ORIGINAL = "6m0j_fixed.pdbqt"

interface_residues = [
    ("A", 19), ("A", 24), ("A", 27), ("A", 28), ("A", 30), ("A", 31),
    ("A", 34), ("A", 35), ("A", 37), ("A", 38), ("A", 41), ("A", 42),
    ("A", 45), ("A", 79), ("A", 82), ("A", 83), ("A", 324), ("A", 325),
    ("A", 326), ("A", 330), ("A", 353), ("A", 354), ("A", 355), ("A", 357),
    ("A", 386), ("A", 393),
    ("E", 403), ("E", 417), ("E", 445), ("E", 446), ("E", 449), ("E", 453),
    ("E", 455), ("E", 456), ("E", 473), ("E", 475), ("E", 476), ("E", 477),
    ("E", 484), ("E", 485), ("E", 486), ("E", 487), ("E", 489), ("E", 490),
    ("E", 493), ("E", 496), ("E", 498), ("E", 500), ("E", 501), ("E", 502),
    ("E", 503), ("E", 505)
]

# ===============================
# WT ENERGY
# ===============================
print("Computing WT interaction energy...\n")

WT_total, WT_lj, WT_elec, WT_solv = compute_interaction_energy(
    PDBQT_ORIGINAL,
    asa_complex="6m0j_fixed.asa",
    asa_A="A.asa",
    asa_E="B.asa",
    return_components=True,
    verbose=True
)

print("\n=== WT ENERGY ===")
print(f"LJ:     {WT_lj:.3f}")
print(f"Elec:   {WT_elec:.3f}")
print(f"Solv:   {WT_solv:.3f}")
print(f"TOTAL:  {WT_total:.3f}\n")

# ===============================
# ALANINE SCANNING (NO FILES)
# ===============================
results = []

for chain, res in interface_residues:
    print(f"\nMutating {chain}{res} → ALA...")

    mut_total, lj, elec, solv = compute_interaction_energy_with_ala(
        PDBQT_ORIGINAL,
        chain,
        res,
        asa_complex="6m0j_fixed.asa",
        asa_A="A.asa",
        asa_E="B.asa",
        return_components=True,
        verbose=False
    )

    ddG = mut_total - WT_total
    results.append((f"{chain}{res}", mut_total, ddG))

    print(f"  ΔG(mut) = {mut_total:.3f} kcal/mol")
    print(f"  ΔΔG     = {ddG:.3f} kcal/mol")


# ===============================
# Save to CSV
# ===============================
import csv

with open("alanine_scan_results.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Residue", "dG_mut", "dDG"])

    for resname, mut_E, ddG in results:
        writer.writerow([resname, f"{mut_E:.3f}", f"{ddG:.3f}"])

print("\nCSV written to alanine_scan_results.csv\n")
