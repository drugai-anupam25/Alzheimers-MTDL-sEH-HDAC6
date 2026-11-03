import os
import sys
import urllib.request
from rdkit import Chem

# --- 1. AUTO-DOWNLOAD SA_SCORE FILES DIRECTLY FROM RDKIT GITHUB ---
def setup_sascorer():
    base_url = "https://raw.githubusercontent.com/rdkit/rdkit/master/Contrib/SA_Score/"
    files = {
        "sascorer.py": base_url + "sascorer.py",
        "fpscores.pkl.gz": base_url + "fpscores.pkl.gz"
    }
    
    for filename, url in files.items():
        if not os.path.exists(filename):
            print(f"🔄 Downloading required RDKit file: {filename}...")
            try:
                urllib.request.urlretrieve(url, filename)
            except Exception as e:
                print(f"❌ Download failed for {filename}: {e}")
                sys.exit(1)

setup_sascorer()

try:
    import sascorer
except ImportError:
    print("❌ Error: SA Score module import fail ho gaya.")
    sys.exit(1)

# --- 2. SMILES DICTIONARY ---
SMILES_DICT = {
    "NXI-07 (Baseline)": "CC1CCC(NC(=O)Nc2ccc(CCC(=O)CO)cc2)CC1",
    "NXI-12-Rigid": "O=C(/C=C/c1ccc(NC(=O)N[C@]2(C3)C[C@]4(C)C[C@H]3C4(C)C2)cc1)CO",
    "NXI-13-Rigid": "O=C(Cc1ccc(NC(=O)N[C@]2(C3)C[C@]4(C)C[C@H]3C4(C)C2)cc1)CO",
    "NXI-14A-Direct": "O=C(CO)c1ccc(NC(=O)N[C@]2(C3)C[C@]4(C)C[C@H]3C4(C)C2)cc1",
    "NXI-14B-Indane": "O=C(CO)C1Cc2ccc(NC(=O)N[C@]3(C4)C[C@]5(C)C[C@H]4C5(C)C3)cc2C1"
}

# --- 3. CALCULATE SCORES ---
print("\n🧪 SYNTHETIC ACCESSIBILITY (SA) SCORE CALCULATOR 🧪")
print("-" * 75)
print("Score Guide: 1-3.5 (Easy) | 3.5-5.5 (Medium) | >5.5 (Lab Nightmare)")
print("-" * 75)

for name, smiles in SMILES_DICT.items():
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        sa_score = round(sascorer.calculateScore(mol), 2)
        
        # Visual feedback
        if sa_score <= 3.5:
            difficulty = "🟢 Easy to Synthesize"
        elif sa_score <= 5.5:
            difficulty = "🟡 Moderate (Requires effort)"
        else:
            difficulty = "🔴 Hard (Avoid if possible)"
            
        print(f"🔹 {name: <18} | SA Score: {sa_score: <5} | {difficulty}")
    else:
        print(f"❌ {name: <18} | Error: Invalid SMILES")

print("-" * 75)