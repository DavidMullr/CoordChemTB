import requests
import os

# List of common ligand IDs (3-letter codes)
ligand_ids = [
    'ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'GMP', 'NAD', 'FAD', 'FMN', 'HEM',
    'COA', 'SAM', 'SAH', 'NAG', 'MAN', 'GLC', 'GAL', 'FUC', 'SIA', 'BMA',
    'SO4', 'PO4', 'CL', 'MG', 'CA', 'ZN', 'FE', 'CU', 'MN', 'CO',
    'NO3', 'NH4', 'K', 'NA', 'IOD', 'BR', 'CS', 'SR', 'BA', 'PB',
    'EDO', 'DMS', 'PEG', 'MPD', 'ACT', 'ACE', 'TRS', 'MES', 'HEZ', 'PGE',
    'IMD', 'MLI', 'MLY', 'MLZ', 'MLW', 'MLV', 'MLU', 'MLT', 'MLR', 'MLQ',
    'MLP', 'MLS', 'MLM', 'MLN', 'MLK', 'MLJ', 'MLI', 'MLH', 'MLG', 'MLF',
    'MLE', 'MLD', 'MLC', 'MLB', 'MLA', 'MLZ', 'MLY', 'MLX', 'MLW', 'MLV',
    'MLU', 'MLT', 'MLR', 'MLQ', 'MLP', 'MLS', 'MLN', 'MLM', 'MLL', 'MLK',
    'MLJ', 'MLI', 'MLH', 'MLG', 'MLF', 'MLE', 'MLD', 'MLC', 'MLB', 'MLA'
]

# Directory to save downloaded ligand files
output_dir = "ligands"
os.makedirs(output_dir, exist_ok=True)

def download_ligand_sdf(lig_id):
    url = f"https://files.rcsb.org/ligands/view/{lig_id}_ideal.sdf"
    response = requests.get(url)
    if response.status_code == 200:
        filename = os.path.join(output_dir, f"{lig_id}.sdf")
        with open(filename, "wb") as f:
            f.write(response.content)
        print(f"Downloaded: {lig_id}")
    else:
        print(f"Failed to download: {lig_id}")

# Download ligand structures
for lig_id in ligand_ids:
    download_ligand_sdf(lig_id)
