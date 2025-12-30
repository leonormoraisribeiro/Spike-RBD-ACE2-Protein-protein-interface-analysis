import csv
import matplotlib.pyplot as plt
import pandas as pd  
import interface_data
from interaction_energy import (
    compute_interaction_energy,
    compute_interaction_energy_with_ala
)

# load data
PDBQT_FILE = "6m0j_fixed.pdbqt"
ASA_COMPLEX = "6m0j_fixed.asa"
ASA_A = "A.asa"
ASA_E = "B.asa"
DDG_THRESHOLD = 1.0  # threshold to consider a residue a "hotspot"
CSV_FILENAME = "alanine_scanning_results.csv"
PLOT_FILENAME = "alanine_scanning_plot.png"

# load interface residues
interface_residues = interface_data.INTERFACE_LIST
print(f"Interface residues: {len(interface_residues)} loaded.")

# calculate WT energy
print("\n1. Calculating WT Energy (Wild Type)")
WT_total, WT_lj, WT_elec, WT_solv = compute_interaction_energy(
    PDBQT_FILE, ASA_COMPLEX, ASA_A, ASA_E, return_components=True, verbose=False
)

print(f"WT Total Energy: {WT_total:.4f} kcal/mol")

# 2. Alanine Scanning
results = []
print("\n 2. Strarting Alanine Scanning")
print(f"{'Residue':<10} {'Total Mut':<12} {'ddG':<12} {'Status'}")

for chain, res in interface_residues:
    # Calculate mutant energy
    mut_total = compute_interaction_energy_with_ala(
        PDBQT_FILE, chain, res, ASA_COMPLEX, ASA_A, ASA_E, return_components=False, verbose=False
    )
    
    ddG = mut_total - WT_total
    
    # Verify Hotspot
    status = "HOTSPOT" if ddG > DDG_THRESHOLD else ""
    print(f"{chain}:{res:<7} {mut_total:12.4f} {ddG:+12.4f} {status}")
    
    # Save results
    results.append({
        "Chain": chain,
        "ResNum": res,
        "Label": f"{chain}:{res}",
        "Energy_Mutant": mut_total,
        "ddG": ddG,
        "Is_Hotspot": ddG > DDG_THRESHOLD
    })

# 3. Save results to CSV
# Create DataFrame
df = pd.DataFrame(results)

# Save to CSV
df.to_csv(CSV_FILENAME, index=False)
print(f"\nCSV saved to: {CSV_FILENAME}")


# 4. Generate Plot
print("\n 3. Generating Plot...")

# Order DataFrame by ddG
df_sorted = df.sort_values("ddG", ascending=False).reset_index(drop=True)

# Define scores colors
colors = ['red' if x > DDG_THRESHOLD else 'skyblue' for x in df_sorted["ddG"]]

plt.figure(figsize=(14, 7))

# Index for x-axis
plt.bar(df_sorted.index, df_sorted["ddG"], color=colors, edgecolor='black', alpha=0.8)
plt.xticks(df_sorted.index, df_sorted["Label"], rotation=90, fontsize=9)

# Reference lines
plt.axhline(y=0, color='black', linewidth=1)
plt.axhline(y=DDG_THRESHOLD, color='gray', linestyle='--', linewidth=1, label=f'Threshold ({DDG_THRESHOLD} kcal/mol)')

plt.ylabel('$\Delta\Delta G$ (kcal/mol) - (Positive = Destabilizing)')
plt.title('Alanine Scanning In-Silico: Contribution per Residue')
plt.legend()

for i, row in df_sorted.iterrows():
    if row["ddG"] > DDG_THRESHOLD:
        plt.text(i, row["ddG"] + 0.1, 
                 row["Label"], 
                 ha='center', va='bottom', fontsize=8, fontweight='bold', color='darkred')

plt.tight_layout()

# Save image
plt.savefig(PLOT_FILENAME, dpi=300)
plt.show()
