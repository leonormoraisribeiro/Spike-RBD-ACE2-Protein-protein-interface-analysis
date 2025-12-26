#!/usr/bin/env bash
set -e

# Usage:
#   ./run_variant.sh E484K
#   ./run_variant.sh L452R

VAR="$1"
if [ -z "$VAR" ]; then
  echo "Usage: $0 <VARIANT_TAG>  (e.g., E484K or L452R)"
  exit 1
fi

COMPLEX="mut_${VAR}_complex.pdb"
A="mut_${VAR}_A.pdb"
B="mut_${VAR}_B.pdb"

# Check required inputs
for f in "$COMPLEX" "$A" "$B" "6m0j_fixed.pdbqt" "6m0j_fixed.asa" "A.asa" "B.asa" "interaction_energy.py"; do
  if [ ! -f "$f" ]; then
    echo "Missing file: $f"
    exit 1
  fi
done

echo "== NACCESS for ${VAR} =="
./naccess "$COMPLEX"
./naccess "$A"
./naccess "$B"

echo "== PDBQT for ${VAR} (OpenBabel) =="
obabel -ipdb "$COMPLEX" -opdbqt -O "mut_${VAR}_complex.pdbqt"

echo "== Energy for ${VAR} =="
python3 - <<'PY'
import sys
from interaction_energy import compute_interaction_energy

VAR = sys.argv[1]

WT_PDBQT = "6m0j_fixed.pdbqt"
WT_ASA_COMPLEX = "6m0j_fixed.asa"
WT_ASA_A = "A.asa"
WT_ASA_B = "B.asa"

MUT_PDBQT = f"mut_{VAR}_complex.pdbqt"
MUT_ASA_COMPLEX = f"mut_{VAR}_complex.asa"
MUT_ASA_A = f"mut_{VAR}_A.asa"
MUT_ASA_B = f"mut_{VAR}_B.asa"

G_wt = compute_interaction_energy(WT_PDBQT, WT_ASA_COMPLEX, WT_ASA_A, WT_ASA_B, return_components=False)
G_mut = compute_interaction_energy(MUT_PDBQT, MUT_ASA_COMPLEX, MUT_ASA_A, MUT_ASA_B, return_components=False)

ddG = G_mut - G_wt

print(f"WT ΔG = {G_wt:.3f} kcal/mol")
print(f"{VAR} ΔG = {G_mut:.3f} kcal/mol")
print(f"ΔΔG ({VAR}) = {ddG:.3f} kcal/mol")
PY "$VAR"
