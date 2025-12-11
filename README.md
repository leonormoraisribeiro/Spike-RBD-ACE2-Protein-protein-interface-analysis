# Spike-RBD-ACE2-Protein-protein-interface-analysis

**Project:** Bio-Physics Exercise 2025-26
**Objective:** To evaluate the relative contribution of interface residues to the interaction energy in the RBD-ACE2 protein-protein complex using the structure **6m0j**.

---

##  Repository Structure

```text
├── data/
│   ├── 6m0j_fixed.pdb       # Cleaned structure (no heteroatoms, H added)
│   ├── 6m0j_fixed.pdbqt     # Structure with partial charges and atom types
│   ├── *.asa                # Atomic solvent accessibility (for Energy calc)
│   └── *.rsa                # Residue solvent accessibility (for Analysis)
├── scripts/
│   ├── structure_setup.py   # Setup pipeline using biobb
│   ├── interface_analysis.py # Interface detection (Distance & ASA)
│   └── energy_calculation.py # (In progress)
├── README.md
└── requirements.txt
```

## Step 1: Structure Preparation & Interface Definition

**Objective:** To prepare the biological system for simulation and identify the amino acid residues forming the protein-protein interface between SARS-CoV-2 RBD and ACE2.

### 1.1 Structure Setup (BioBB)
The raw structure **6m0j** was processed to fix experimental artifacts and prepare it for energy calculations.

* [cite_start]**Tool:** `biobb_structure_checking`[cite: 33, 77].
* **Workflow:**
    1.  [cite_start]**Chain Selection:** Isolated chains **E** (RBD) and **A** (ACE2)[cite: 30].
    2.  [cite_start]**Cleaning:** Removed crystallographic waters and heteroatoms[cite: 32].
    3.  **Fixing:** Added missing side-chain atoms and hydrogens.
    4.  [cite_start]**Charges:** Added partial charges and atom types, generating a **PDBQT** file for electrostatic calculations[cite: 33].
* **Key Outputs:**
    * `6m0j_fixed.pdb` (Clean PDB)
    * `6m0j_fixed.pdbqt` (PDB with charges/atom types)

### 1.2 Solvent Accessibility (NACCESS)
[cite_start]To evaluate solvation effects, we calculated the Accessible Surface Area (ASA)[cite: 82].

* **Calculation:** Performed on the bound complex (`A+E`) and unbound components (`A` and `E` isolated).
* **Outputs:**
    * `.rsa` files: Residue-level data (used for quick visual analysis).
    * [cite_start]`.asa` files: Atomic-level data (used for $\Delta G_{solv}$ calculation)[cite: 38].

### 1.3 Defining the Interface
[cite_start]We implemented two complementary methods to define the interface residues[cite: 14, 23, 49]:

1.  **Distance Criterion (Geometric):**
    * [cite_start]Residues are considered part of the interface if any of their atoms are within a cut-off distance (e.g., 5Å) of an atom in the opposing chain[cite: 23].
    * *Implementation:* Python script using `Bio.PDB.NeighborSearch`.

2.  **Solvation Criterion ($\Delta$ASA):**
    * Residues are defined as interface if they lose solvent accessibility upon complex formation.
    * $\Delta ASA = ASA_{unbound} - ASA_{bound} > \text{cutoff}$
    * *Status:* Implemented using custom Python parser for `.rsa`/`.asa` files.
