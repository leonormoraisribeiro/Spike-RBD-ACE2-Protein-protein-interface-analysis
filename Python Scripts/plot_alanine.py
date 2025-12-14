import pandas as pd
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("alanine_scan_results.csv")

# Sort residues by ΔΔG (highest = most important for stability)
df = df.sort_values("dDG", ascending=False)

# Define threshold for “important” residues
threshold = df["dDG"].mean() + df["dDG"].std()

# Colors: red for important, blue for others
colors = ["red" if x >= threshold else "steelblue" for x in df["dDG"]]

# Plot
plt.figure(figsize=(14, 6))
plt.bar(df["Residue"], df["dDG"], color=colors)
plt.xticks(rotation=90)
plt.ylabel("ΔΔG (kcal/mol)")
plt.title("Alanine Scanning – Residue Contributions to Interface Stability")

# Highlight annotation
for i, row in df.iterrows():
    if row["dDG"] >= threshold:
        plt.text(i, row["dDG"] + 0.5,
                 f"{row['Residue']}",
                 ha="center", va="bottom", fontsize=8)

plt.tight_layout()
plt.savefig("alanine_scan_plot.png", dpi=300)
plt.show()

print("Plot saved as alanine_scan_plot.png")
