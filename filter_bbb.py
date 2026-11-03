import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

INPUT_FILE = "top_leads_new.csv"
OUTPUT_FILE = "new_BBB_Passed_Finalists.csv"

print(f"Loading {INPUT_FILE} for ADMET/BBB screening...")
df = pd.read_csv('new_mutant.csv', sep='\t', header=None, names=['SMILES', 'ID'])

results = []

for index, row in df.iterrows():
    lig_id = row['ID']
    smiles = row['SMILES']
    
    # RDKit se molecule read karna
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        continue
        
    # Pharmacokinetics (PK) Parameters calculate karna
    mw = Descriptors.MolWt(mol)           # Molecular Weight (< 450 target)
    logp = Descriptors.MolLogP(mol)       # Lipophilicity / Tel jaisa nature (2 to 5 target)
    tpsa = Descriptors.TPSA(mol)          # Polar Surface Area (< 90 target for Brain)
    hbd = rdMolDescriptors.CalcNumHBD(mol)# Hydrogen Bond Donors (<= 3 target)
    hba = rdMolDescriptors.CalcNumHBA(mol)# Hydrogen Bond Acceptors
    
    # Strict BBB Filter Logic (Pfizer rules for CNS)
   # MTDL Relaxed BBB Filter Logic (CNS threshold for dual-targets)
    if tpsa <= 120 and hbd <= 5 and mw <= 480 and logp <= 5.5:
        bbb_status = "✅ PASS"
    else:
        bbb_status = "❌ FAIL"
        
    results.append({
        'ID': lig_id,
        'SMILES': smiles,
        'Average_Score': row.get('Average_Score', 'N/A'),
        'MW': round(mw, 2),
        'LogP': round(logp, 2),
        'TPSA': round(tpsa, 2),
        'HBD': hbd,
        'HBA': hba,
        'BBB_Status': bbb_status
    })

# Data ko table mein convert karna
results_df = pd.DataFrame(results)

# Saare results (Pass aur Fail dono) ek jagah save karna taaki check kar sakein
results_df.to_csv("ADMET_Full_Stats.csv", index=False)

# Sirf PASS wale filter karna
passed_df = results_df[results_df['BBB_Status'] == '✅ PASS']
passed_df.to_csv(OUTPUT_FILE, index=False)

print("\n🧠 === BBB FILTER RESULTS === 🧠")
print(f"Total evaluated: {len(results_df)}")
print(f"Successfully passed BBB criteria: {len(passed_df)}\n")

if len(passed_df) > 0:
    print(passed_df[['ID', 'Average_Score', 'MW', 'TPSA', 'LogP']].to_string(index=False))
    print(f"\nBOOM! Tere final filtered candidates '{OUTPUT_FILE}' mein save ho gaye hain.")
else:
    print("Warning: Ek bhi molecule ne strict BBB rules paas nahi kiye! TPSA limit thodi badhani padegi.")