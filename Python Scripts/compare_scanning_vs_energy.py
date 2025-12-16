import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import interface_data
from interaction_energy import (
    compute_interaction_energy,
    compute_interaction_energy_with_ala,
    compute_wt_residue_contribution
)


PDBQT_FILE = "6m0j_fixed.pdbqt"
ASA_COMPLEX = "6m0j_fixed.asa"
ASA_A = "A.asa"
ASA_E = "B.asa"


interface_residues = interface_data.INTERFACE_LIST
print(f"Interface residues: {len(interface_residues)}")


print("Calculando energia WT global...")
WT_total = compute_interaction_energy(
    PDBQT_FILE, ASA_COMPLEX, ASA_A, ASA_E, return_components=False
)

results = []

print(f"{'Residue':<10} {'ddG (Ala)':<12} {'Direct Energy':<15}")

for chain, res in interface_residues:
    mut_total = compute_interaction_energy_with_ala(
        PDBQT_FILE, chain, res, ASA_COMPLEX, ASA_A, ASA_E, return_components=False
    )
    ddG = mut_total - WT_total 
    direct_energy = compute_wt_residue_contribution(
        PDBQT_FILE, chain, res, ASA_COMPLEX, ASA_A, ASA_E
    )

    print(f"{chain}:{res:<7} {ddG:12.4f} {direct_energy:15.4f}")

    results.append({
        "Label": f"{chain}:{res}",
        "ddG": ddG,
        "Direct_Energy": direct_energy
    })

df = pd.DataFrame(results)


plt.figure(figsize=(8, 8))
plt.scatter(df["Direct_Energy"], df["ddG"], color='blue', alpha=0.7, edgecolors='k')


z = np.polyfit(df["Direct_Energy"], df["ddG"], 1)
p = np.poly1d(z)
plt.plot(df["Direct_Energy"], p(df["Direct_Energy"]), "r--", label=f"Fit: y={z[0]:.2f}x + {z[1]:.2f}")


plt.xlabel("Direct Interaction Energy (kcal/mol)\n(Negative = Attractive)")
plt.ylabel("Alanine Scanning $\Delta\Delta G$ (kcal/mol)\n(Positive = Destabilizing)")
plt.title("Correlation: Direct Energy vs Alanine Scanning")
plt.grid(True, linestyle='--', alpha=0.5)
plt.axhline(0, color='black', linewidth=0.8)
plt.axvline(0, color='black', linewidth=0.8)


corr = df["Direct_Energy"].corr(df["ddG"])
plt.legend([f"Correlation (R): {corr:.4f}", "Trendline"])


for i, row in df.iterrows():
    if row["ddG"] > 1.0 or row["Direct_Energy"] < -2.0: 
        plt.text(row["Direct_Energy"], row["ddG"], row["Label"], fontsize=8)

plt.tight_layout()
plt.savefig("correlation_plot.png", dpi=300)
print(f"\nCorrelation Plot saved to correlation_plot.png")
print(f"Correlation Coefficient: {corr:.4f}")


output_filename = "comparison_plot.png"
plt.savefig(output_filename, dpi=300)
print(f"\nGráfico guardado como: {output_filename}")