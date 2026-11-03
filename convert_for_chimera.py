import os
import urllib.request
from rdkit import Chem
from rdkit.Chem import AllChem

# Output folder for ChimeraX files
out_dir = "ChimeraX_3D_Structures"
os.makedirs(out_dir, exist_ok=True)

# Tere 6 Evolution Molecules (Exact SMILES + 3-letter Residue Names for ChimeraX)
ligands = {
    "00_M9_Parent": {
        "smiles": "CC12CC3CC(C)(C1)CC(C3)(C2)NC(=O)Nc4ccc(C(=O)NCCCCCC(=O)NO)cc4",
        "resname": "M9P"
    },
    "01_NXI_01_Initial_Lead": {
        "smiles": "CC1(C)C2CCC1(C)C(NC(=O)Nc3ccc(OCC(O)C(=O)O)cc3)C2",
        "resname": "NX1"
    },
    "02_NXI_07_Thiol_Tox": {
        "smiles": "FC(F)(F)C1CCC(NC(=O)Nc2ccc(CCCC(=O)NCCS)cc2)CC1",
        "resname": "NX7"
    },
    "03_NXI_07_Clean_Pyridyl": {
        "smiles": "CC1(C)C2CCC1(C)C(NC(=O)Nc3ccc(CCCCC(=O)Nc4ccccn4)cc3)C2",
        "resname": "N7C"
    },
    "04_NXI_09_Gold_Ketol": {
        "smiles": "CC1(C)C2CCC1(C)C(NC(=O)Nc3ccc(CCCCC(=O)CO)cc3)C2",
        "resname": "NX9"
    },
    "05_NXI_13_Rigid_Indane": {
        "smiles": "CC1(C)[C@@H]2CC[C@@]1(C)[C@H](NC(=O)Nc3ccc4CC[C@H](CC(=O)CO)c4c3)C2",
        "resname": "N13"
    }
}

print("="*65)
print("  GENERATING 3D OPTIMIZED .PDB & .SDF FILES FOR CHIMERAX")
print("="*65)

for name, data in ligands.items():
    mol = Chem.MolFromSmiles(data["smiles"])
    if mol is None:
        print(f"❌ Error in SMILES for {name}")
        continue
        
    # Hydrogens add karna zaroori hai 3D geometry aur H-bonds ke liye
    mol_3d = Chem.AddHs(mol)
    
    # ETKDGv3 algorithm se best 3D conformer generate karna
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMolecule(mol_3d, params)
    
    # MMFF94 Force Field se 3D structure ko relax/minimize karna
    AllChem.MMFFOptimizeMolecule(mol_3d, maxIters=1000)
    
    # PDB Residue name set karna taaki ChimeraX mein select karna aasan ho
    for atom in mol_3d.GetAtoms():
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName(data["resname"])
        info.SetResidueNumber(1)
        info.SetChainId("L")
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
        
    pdb_path = os.path.join(out_dir, f"{name}.pdb")
    sdf_path = os.path.join(out_dir, f"{name}.sdf")
    
    Chem.MolToPDBFile(mol_3d, pdb_path)
    
    writer = Chem.SDWriter(sdf_path)
    writer.write(mol_3d)
    writer.close()
    
    print(f"✅ Saved: {pdb_path}  (Residue ID: :{data['resname']})")

# Target Proteins (HDAC6: 5EDU & sEH: 3OTQ) bhi download kar le agar folder mein nahi hain
proteins = {"5EDU_HDAC6.pdb": "https://files.rcsb.org/download/5EDU.pdb",
            "3OTQ_sEH.pdb": "https://files.rcsb.org/download/3OTQ.pdb"}

for p_name, url in proteins.items():
    p_path = os.path.join(out_dir, p_name)
    if not os.path.exists(p_path):
        try:
            print(f"⬇️ Downloading {p_name} from RCSB PDB...")
            urllib.request.urlretrieve(url, p_path)
            print(f"✅ Downloaded: {p_path}")
        except Exception as e:
            print(f"⚠️ Could not auto-download {p_name}: {e}")

# ChimeraX Automation Script (.cxc) generate karna
cxc_script = """# ChimeraX Visualisation Script for NXI Series
set bgColor white
lighting soft
graphics silhouettes true

# Open all generated 3D Ligands in a clean grid
open 00_M9_Parent.pdb
open 01_NXI_01_Initial_Lead.pdb
open 02_NXI_07_Thiol_Tox.pdb
open 03_NXI_07_Clean_Pyridyl.pdb
open 04_NXI_09_Gold_Ketol.pdb
open 05_NXI_13_Rigid_Indane.pdb

show atoms
style ball
size atomRadius 0.4
size stickRadius 0.2
color byhetero
color #1 & C gray
color #2 & C orange
color #3 & C magenta
color #4 & C cyan
color #5 & C gold
color #6 & C limegreen

tile
label #1-6 text "{0.name}" height 1.2 color black yoffset 4
"""

cxc_path = os.path.join(out_dir, "view_all_ligands.cxc")
with open(cxc_path, "w") as f:
    f.write(cxc_script)

print("\n" + "="*65)
print(f"🎯 DONE! Ab ChimeraX khol aur '{cxc_path}' open kar le!")
print("="*65)