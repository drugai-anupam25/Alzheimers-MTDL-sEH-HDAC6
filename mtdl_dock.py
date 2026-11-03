import os
import subprocess
import urllib.request
import numpy as np
import pandas as pd
import warnings
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation

warnings.filterwarnings("ignore")

# --- 1. SETUP & CANDIDATES ---
OUTPUT_DIR = "CDRI_Final_Dual_Docking"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, "ligands"), exist_ok=True)
RESULTS_CSV = os.path.join(OUTPUT_DIR, "Master_Dual_Docking_Report.csv")

TARGETS = {
    "sEH": {"pdb": "3OTQ", "ref_ligand": "CIU"},
    "HDAC6": {"pdb": "5EDU", "ref_ligand": "ZN"} # Removed strict spacing for safer parsing
}

CANDIDATES = {
    "M9 (Stage-0 Parent)": "O=C(NC1=CC=C(C(NCCCCCCC(NO)=O)=O)C=C1)N[C@]2(C3)C[C@]4(C)C[C@](C[C@H]3C4)(C)C2",
    "NXI-01 (Stage-1 Orig)": "O=C(NC1=CC=C(OCC(O)C(=O)O)C=C1)N[C@]2(C3)C[C@]4(C)C[C@H]3C4(C)C2",
    "NXI-05-Clean (Stage-3A)": "O=C(CCc1ccc(NC(=O)N[C@]2(C3)C[C@]4(C)C[C@H]3C4(C)C2)cc1)Nc5ccccn5",
    "NXI-06-Clean (Stage-3B)": "CC1CCC(NC(=O)Nc2ccc(CCC(=O)Nc3ccccn3)cc2)CC1",
    "NXI-07-Clean (Stage-3C)": "CC1CCC(NC(=O)Nc2ccc(CCC(=O)CO)cc2)CC1",
    "NXI-08-Clean (Stage-3D)": "CC1CCC(NC(=O)Nc2ccc(CCC(=O)N(C)O)cc2)CC1"
}

# --- 2. ZINC FIX FOR VINA CRASH ---
def fix_zinc_pdbqt(pdbqt_file):
    if not os.path.exists(pdbqt_file): return
    with open(pdbqt_file, 'r') as f:
        lines = f.readlines()
    with open(pdbqt_file, 'w') as f:
        for line in lines:
            if " ZN " in line or " Zn " in line:
                if "HETATM99999" in line:
                    line = line.replace("HETATM99999", "HETATM 9999 ")
                if len(line) >= 77:
                    line = line[:77] + " Zn\n"
            f.write(line)

# --- 3. SMART PROTEIN CLEANER & GRID BOX (FIXED) ---
def clean_and_get_box(pdb_id, ref_name, out_pdb):
    coords = []
    protein_coords = []
    raw_pdb = os.path.join(OUTPUT_DIR, f"{pdb_id}_raw.pdb")
    
    if not os.path.exists(raw_pdb):
        urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pdb_id}.pdb', raw_pdb)
    
    with open(raw_pdb, 'r') as fin, open(out_pdb, 'w') as fout:
        for line in fin:
            if line.startswith("ATOM"):
                fout.write(line)
                try:
                    protein_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except ValueError:
                    pass
            elif line.startswith("HETATM"):
                # Loose matching for metal/ligand to avoid zero-size array crashes
                if "ZN" in line or "Zn" in line:
                    fout.write(line)
                    if ref_name == "ZN":
                        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                elif ref_name in line and "HOH" not in line:
                    coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            elif line.startswith("TER") or line.startswith("END"):
                fout.write(line)
                
    # FAIL-SAFE: If ligand not found, use protein center
    if len(coords) == 0:
        print(f"    [!] Warning: {ref_name} not found. Using whole protein center.")
        coords = protein_coords

    coords = np.array(coords)
    center = np.mean(coords, axis=0)
    sizes = np.max(coords, axis=0) - np.min(coords, axis=0) + 20.0 
    sizes = np.clip(sizes, a_min=22.0, a_max=30.0) 
    return center, sizes

# --- 4. PREPARE RECEPTORS ---
print("🔥 DOWNLOADING, CLEANING & FIXING RECEPTORS 🔥")
target_boxes = {}
for target_name, data in TARGETS.items():
    clean_pdb = os.path.join(OUTPUT_DIR, f"{data['pdb']}_clean.pdb")
    receptor_pdbqt = os.path.join(OUTPUT_DIR, f"{data['pdb']}_receptor.pdbqt")
    
    center, sizes = clean_and_get_box(data['pdb'], data['ref_ligand'], clean_pdb)
    
    if not os.path.exists(receptor_pdbqt):
        os.system(f"obabel {clean_pdb} -O {receptor_pdbqt} -xr -p 7.4 2> /dev/null")
        if target_name == "HDAC6":
            fix_zinc_pdbqt(receptor_pdbqt)
            
    target_boxes[target_name] = {"pdbqt": receptor_pdbqt, "center": center, "sizes": sizes}

# --- 5. BATCH PREPARATION & DOCKING ---
print("\n⚡ STARTING BATCH DUAL-TARGET DOCKING ⚡")
all_results = []

for lig_name, smiles in CANDIDATES.items():
    print(f"\n🚀 Processing: {lig_name}")
    safe_name = lig_name.split()[0]
    ligand_pdbqt = os.path.join(OUTPUT_DIR, "ligands", f"{safe_name}.pdbqt")
    
    ligand_scores = {"Candidate": lig_name, "sEH_Score": None, "HDAC6_Score": None}
    
    # Generate 3D Conformation (FIXED RDKIT BUG)
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42) # Removed ETKDGv3() to prevent C++ signature crash
        AllChem.MMFFOptimizeMolecule(mol)
        prep = MoleculePreparation()
        prep.prepare(mol)
        with open(ligand_pdbqt, 'w') as f:
            f.write(prep.write_pdbqt_string())
    except Exception as e:
        print(f"    [!] Error preparing ligand {safe_name}: {e}")
        all_results.append(ligand_scores)
        continue

    # Dock into both targets
    for target_name, box in target_boxes.items():
        docked_out = os.path.join(OUTPUT_DIR, f"{safe_name}_in_{target_name}.pdbqt")
        log_file = os.path.join(OUTPUT_DIR, f"{safe_name}_in_{target_name}.log")
        
        print(f"    [->] Docking into {target_name}...")
        vina_cmd = (
            f"vina --receptor {box['pdbqt']} --ligand {ligand_pdbqt} "
            f"--center_x {box['center'][0]:.3f} --center_y {box['center'][1]:.3f} --center_z {box['center'][2]:.3f} "
            f"--size_x {box['sizes'][0]:.3f} --size_y {box['sizes'][1]:.3f} --size_z {box['sizes'][2]:.3f} "
            f"--cpu 8 --exhaustiveness 8 --out {docked_out} > {log_file} 2>/dev/null"
        )
        subprocess.run(vina_cmd, shell=True)
        
        # Parse output
        if os.path.exists(log_file):
            with open(log_file, "r") as f:
                for line in f:
                    if line.strip().startswith("1"):
                        ligand_scores[f"{target_name}_Score"] = float(line.split()[1])
                        break
                        
    all_results.append(ligand_scores)

# --- 6. SAVE MASTER CSV REPORT ---
print("\n📊 SAVING MASTER REPORT...")
df = pd.DataFrame(all_results)
df.to_csv(RESULTS_CSV, index=False)
print(f"Mission Complete! Saare 6 molecules ka data '{RESULTS_CSV}' mein save ho gaya hai.")
print(df.to_string(index=False))