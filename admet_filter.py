import pandas as pd
from rdkit import Chem
from rdkit.Chem import FilterCatalog

INPUT_FILE = "new_advanceNONTOX.csv"
OUTPUT_FILE = "top_leads_new.csv"

print(f"Loading {INPUT_FILE} for ADMET Tox Screening...")
df = pd.read_csv(INPUT_FILE, sep='\t', header=None, names=['SMILES', 'ID'])

# RDKit ke built-in PAINS aur BRENK (Toxicity) filters setup karna
params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
catalog = FilterCatalog.FilterCatalog(params)

# Apna custom Hydroxamate filter (Ames/Lossen risk)
hydroxamate_pattern = Chem.MolFromSmarts('C(=O)NO')

results = []

for index, row in df.iterrows():
    lig_id = row['ID']
    smiles = row['SMILES']
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        continue
        
    # 1. PAINS & BRENK Check (Check for toxic kachra)
    if catalog.HasMatch(mol):
        match = catalog.GetFirstMatch(mol)
        tox_alert = match.GetDescription()
        status = "❌ TOXIC"
    else:
        tox_alert = "Clean"
        status = "✅ PASS"
        
    # 2. Custom Hydroxamate Ames Risk Check
    has_hydroxamate = mol.HasSubstructMatch(hydroxamate_pattern)
    if has_hydroxamate:
        mutagenicity_risk = "High (Hydroxamate)"
    else:
        mutagenicity_risk = "Low (Bioisostere)"
        
    results.append({
        'ID': lig_id,
        'SMILES': smiles,
        'Average_Score': row.get('Average_Score', 'N/A'),
        'Tox_Alert': tox_alert,
        'Mutagenicity_Risk': mutagenicity_risk,
        'ADMET_Status': status
    })

results_df = pd.DataFrame(results)
# Poora data save kar rahe hain taaki reason pata chale kaun kyun fail hua
results_df.to_csv("ADMET_Full_Report.csv", index=False)

# Sirf saaf molecules aage badhenge
passed_df = results_df[results_df['ADMET_Status'] == '✅ PASS']
passed_df.to_csv(OUTPUT_FILE, index=False)

print("\n☠️ === ADMET & TOXICITY FILTER RESULTS === ☠️")
print(f"Total evaluated: {len(results_df)}")
print(f"Passed strict structural toxicity checks: {len(passed_df)}\n")

if len(passed_df) > 0:
    print(passed_df[['ID', 'Tox_Alert', 'Mutagenicity_Risk', 'ADMET_Status']].to_string(index=False))
    print(f"\nAwesome! Tere final safe candidates '{OUTPUT_FILE}' mein save ho gaye hain.")
else:
    print("Warning: Sab mein toxicophores nikal aaye! 'ADMET_Full_Report.csv' check kar.")