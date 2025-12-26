from interaction_energy import compute_interaction_energy

# ----------------------------
# WT files (in this folder)
# ----------------------------
WT_PDBQT = "6m0j_fixed.pdbqt"
WT_ASA_COMPLEX = "6m0j_fixed.asa"
WT_ASA_A = "A.asa"
WT_ASA_B = "B.asa"

# ----------------------------
# Mutant files (L452R)
# ----------------------------
MUT_PDBQT = "mut_L452R_complex.pdbqt"
MUT_ASA_COMPLEX = "mut_L452R_complex.asa"
MUT_ASA_A = "mut_L452R_A.asa"
MUT_ASA_B = "mut_L452R_B.asa"

print("Final energy calculation for interface residues...")
print("Calculating WT interaction energy...")

G_wt = compute_interaction_energy(
    WT_PDBQT,
    WT_ASA_COMPLEX,
    WT_ASA_A,
    WT_ASA_B,
    return_components=False
)
print(f"WT ΔG = {G_wt:.3f} kcal/mol")

print("\nCalculating L452R mutant interaction energy...")

G_mut = compute_interaction_energy(
    MUT_PDBQT,
    MUT_ASA_COMPLEX,
    MUT_ASA_A,
    MUT_ASA_B,
    return_components=False
)
print(f"L452R ΔG = {G_mut:.3f} kcal/mol")

ddG = G_mut - G_wt

print("\n==============================")
print(f"ΔΔG (L452R) = {ddG:.3f} kcal/mol")
print("==============================")
