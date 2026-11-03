# Alzheimers-MTDL-sEH-HDAC6
# 🧠 Rational Design & Optimization of CNS-Permeable Dual sEH/HDAC6 Inhibitors (NXI Series) for Alzheimer's Disease

## 📌 Project Overview
This repository contains the end-to-end computational medicinal chemistry pipeline for designing **Multi-Target Directed Ligands (MTDLs)** that simultaneously inhibit **Soluble Epoxide Hydrolase (sEH)** and **Histone Deacetylase 6 (HDAC6)** to target neuroinflammation and tau pathology in Alzheimer's disease.

Starting from the peripheral neuropathic pain scaffold **M9** (*Chen & Hammock et al., J. Med. Chem. 2024*), we performed iterative structure-based optimization across multiple generations to solve critical Blood-Brain Barrier (BBB) permeability, ionization, and toxicity bottlenecks.


---

## 🔬 Key Scientific Milestones (Stage 0 to Stage 5)
* **Stage 0 (Parent M9):** Identified CNS bottlenecks — high polarity ($\text{TPSA} = 119.56\text{ \AA}^2$) and high conformational flexibility (11 rotatable bonds).
* **Stage 1 (NXI-01 / F01):** Reduced molecular weight using a bicyclic Bornyl cage and hydroxylated ether linker ($\text{TPSA} = 107.89\text{ \AA}^2$), but identified physiological ionization ($\text{-COO}^-$ at pH 7.4) and rapid 1 ns explicit-water MD unbinding.
* **Stage 2 (Gen-1 Redesign):** Replaced ether oxygens with pure hydrocarbon linkers; diagnosed Lossen-rearrangement mutagenicity risks in free hydroxamates (`NXI-05`) and immunotoxicity alerts in free thiols (`NXI-07`) via ProTox-3.0.
* **Stage 3 & 4 (Gen-2 Clean Leads — NXI-06-Clean, NXI-07-Clean, NXI-09-Gold):** Engineered non-hydroxamate Zinc-Binding Groups (**2-Pyridylamide** and **$\alpha$-Hydroxy Ketone / Ketol warheads**), achieving strictly optimal CNS polarity ($\text{TPSA} = 78.43 - 83.12\text{ \AA}^2$), balanced dual-target binding ($\Delta \text{Affinity} < 0.1\text{ kcal/mol}$), and Class 5 safety ($\text{LD}_{50} = 5000\text{ mg/kg}$).
* **Stage 5 (Conformational Rigidification — NXI-12 / NXI-13):** Locked flexible linkers into Cinnamoyl and Indane rings, reducing rotatable bonds from 8 down to 5 to minimize entropic penalties.

---

## ⚙️ Computational Pipeline & Tools Used
* **Cheminformatics & CNS Profiling:** Python, RDKit (Ertl TPSA, Crippen cLogP, Wager Pfizer CNS MPO Score, BOILED-Egg, PAINS/Brenk Toxicophore Catalogs)
* **Structure-Based Dual & Reverse Docking:** AutoDock Vina v1.2.7, Meeko, OpenBabel (Targeted against human sEH `PDB: 3OTQ/3I28`, human HDAC6 CD2 `PDB: 5EDU` with $\text{Zn}^{2+}$ coordination, and off-target safety panels including `hERG`, `CYP3A4`, `AChE`, and `MAO-B`)
* **Machine Learning Toxicity Audit:** ProTox-3.0 (`tox.charite.de`)
* **Molecular Dynamics:** OpenMM (Apple Silicon GPU pipeline)

---

## 🚀 How to Run the Pipeline
1. Install dependencies:
   ```bash
   pip install -r requirements.txt
