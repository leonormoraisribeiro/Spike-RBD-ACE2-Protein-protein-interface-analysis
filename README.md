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

* **Tool:** `biobb_structure_checking`.
* **Workflow:**
    1.  **Chain Selection:** Isolated chains **E** (RBD) and **A** (ACE2).
    2.  **Cleaning:** Removed crystallographic waters and heteroatoms.
    3.  **Fixing:** Added missing side-chain atoms and hydrogens.
    4.  **Charges:** Added partial charges and atom types, generating a **PDBQT** file for electrostatic calculations.
* **Key Outputs:**
    * `6m0j_fixed.pdb` (Clean PDB)
    * `6m0j_fixed.pdbqt` (PDB with charges/atom types)

### 1.2 Solvent Accessibility (NACCESS)
To evaluate solvation effects, we calculated the Accessible Surface Area (ASA).

* **Calculation:** Performed on the bound complex (`A+E`) and unbound components (`A` and `E` isolated).
* **Outputs:**
    * `.rsa` files: Residue-level data (used for quick visual analysis).
    * `.asa` files: Atomic-level data (used for $\Delta G_{solv}$ calculation).

### 1.3 Defining the Interface
We implemented two complementary methods to define the interface residues:

1.  **Distance Criterion (Geometric):**
    * Residues are considered part of the interface if any of their atoms are within a cut-off distance (e.g., 5Å) of an atom in the opposing chain.
    * *Implementation:* Python script using `interface_distance.py`.

2.  **Solvation Criterion (△ASA):**
    * Residues are defined as interface if they lose solvent accessibility upon complex formation.
    * $\Delta ASA = ASA_{unbound} - ASA_{bound} > \text{cutoff}$
    * *Implementation:* Python script using `interface_asa_variation.py`.

3. **Create a Pymol Visualization of the Interface**
    * With the residues calculated with the △ASA, a python script was created to visualize in python the residues selected
    * *Implementation:* Python script using `create_interface_pymol.py`.
