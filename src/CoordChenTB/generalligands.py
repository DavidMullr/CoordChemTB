# LIGAND SDF FILES DOWNLOADER FOR RCSB.ORG

import requests
import os

# Allowing user to append new ligands by typing their 3-letter RCSB codes
default_ligand_ids = [
    'ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'GMP', 'NAD', 'FAD', 'FMN', 'HEM',
    'COA', 'SAM', 'SAH', 'NAG', 'MAN', 'GLC', 'GAL', 'FUC', 'SIA', 'BMA',
    'SO4', 'PO4', 'CL', 'MG', 'CA', 'ZN', 'FE', 'CU', 'MN', 'CO',
    'NO3', 'NH4', 'K', 'NA', 'IOD', 'BR', 'CS', 'SR', 'BA', 'PB',
    'EDO', 'DMS', 'PEG', 'MPD', 'ACT', 'ACE', 'TRS', 'MES', 'HEZ', 'PGE'
]

print("\nYou can add extra RCSB ligand codes (3-letter), or press ENTER to skip.")
while True:
    new_code = input("Enter ligand code to add (or ENTER to finish): ").strip().upper()
    if not new_code:
        break
    if len(new_code) == 3 and new_code.isalpha():
        if new_code not in default_ligand_ids:
            default_ligand_ids.append(new_code)
            print(f"✅ Added: {new_code}")
        else:
            print("⚠️ Already in the list.")
    else:
        print("❌ Invalid code. Must be a 3-letter alphabetic code.")

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
        print(f"✅ Downloaded: {lig_id}")
    else:
        print(f"❌ Failed to download: {lig_id}")

# Downloading each ligand
for lig_id in default_ligand_ids:
    download_ligand_sdf(lig_id)

#CREATING SDF FILES FOR SPECHTROCHEMICAL SERIES LIGANDS THAT ARE NOT ON RCSB (with their spechtrochemical series info)  

import re
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

# Ligands: name -> (SMILES, field strength)
ligands = {
    "I⁻": ("[I-]", "very weak"),
    "Br⁻": ("[Br-]", "very weak"),
    "S²⁻": ("[S-2]", "very weak"),
    "SCN⁻ (S-bound)": ("S=C=N", "weak"),
    "Cl⁻": ("[Cl-]", "weak"),
    "NO₃⁻": ("O=N(=O)[O-]", "weak"),
    "F⁻": ("[F-]", "weak"),
    "OH⁻": ("[OH-]", "moderate"),
    "H₂O": ("O", "moderate"),
    "C₂O₄²⁻": ("O=C([O-])C(=O)[O-]", "moderate"),
    "NCS⁻ (N-bound)": ("N=C=S", "moderate"),
    "py": ("c1ccncc1", "moderate"),
    "NH₃": ("N", "moderate-strong"),
    "en": ("NCCN", "moderate-strong"),
    "MeCN": ("CC#N", "moderate-strong"),
    "DMF": ("CN(C)C=O", "moderate-strong"),
    "DMSO": ("CS(=O)C", "moderate"),
    "bipy": ("c1ccnc(c1)c2ccccn2", "strong"),
    "phen": ("c1ccc2c(c1)ccnc2c3cccnc3", "strong"),
    "NO₂⁻": ("[O-][N+](=O)", "strong"),
    "PPh₃": ("P(c1ccccc1)(c2ccccc2)(c3ccccc3)", "strong"),
    "CN⁻": ("[C-]#N", "very strong"),
    "CO": ("[C-]#[O+]", "very strong"),
    "H⁻": ("[H-]", "very strong"),
    "CH₃⁻": ("[CH3-]", "very strong"),
    "THF": ("C1CCOC1", "moderate"),
    "EtOH": ("CCO", "moderate"),
    "MeOH": ("CO", "moderate"),
    "Acetate": ("CC(=O)[O-]", "moderate"),
    "Formate": ("C(=O)[O-]", "moderate"),
    "Pyrrolidine": ("C1CCNC1", "strong"),
    "Aniline": ("c1ccccc1N", "moderate"),
    "Imidazole": ("c1cnc[nH]1", "moderate-strong")
}

# Function to clean names to avoid later errors

def clean_filename(name):
    replacements = {
        "⁻": "minus", "⁺": "plus",
        "₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4",
        "₅": "5", "₆": "6", "₇": "7", "₈": "8", "₉": "9",
        "²": "2", "³": "3"
    }
    for uni, ascii_equiv in replacements.items():
        name = name.replace(uni, ascii_equiv)
    name = re.sub(r"[^\w\-]", "_", name)
    return name

# Generating and writing the SDF files

for name, (smiles, field_strength) in ligands.items():
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"❌ Could not parse: {name}")
        continue

  # Adding hydrogens and generating 3D geometry

    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
        print(f"⚠️  3D embedding failed: {name}")
        continue
    AllChem.UFFOptimizeMolecule(mol)

# Setting metadata fields
    mol.SetProp("LigandName", name)
    mol.SetProp("SMILES", smiles)
    mol.SetProp("FieldStrength", field_strength)
    mol.SetProp("MolecularWeight", str(round(Descriptors.MolWt(mol), 2)))
    mol.SetProp("HDonors", str(Descriptors.NumHDonors(mol)))
    mol.SetProp("HAcceptors", str(Descriptors.NumHAcceptors(mol)))
    mol.SetProp("LogP", str(round(Descriptors.MolLogP(mol), 2)))
    mol.SetProp("NumAtoms", str(mol.GetNumAtoms()))

  # Saving with cleaned filename

    sdf_path = os.path.join(output_dir, f"{clean_filename(name)}.sdf")
    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)
    writer.close()
    print(f"✅ SDF saved: {sdf_path}")

# CALCULATING LIGAND FORMAL CHARGE AND PUTING IT IN THE SDF FILE

for filename in os.listdir(output_dir):
    if filename.endswith(".sdf"):
        file_path = os.path.join(output_dir, filename)
        try:
            suppl = Chem.SDMolSupplier(file_path, removeHs=False)
            mol = suppl[0] if suppl and len(suppl) > 0 else None
            if mol:
                net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
                mol.SetProp("NetFormalCharge", str(net_charge))

                writer = Chem.SDWriter(file_path)
                writer.write(mol)
                writer.close()
                print(f"✅ Updated {filename} with NetFormalCharge = {net_charge}")
            else:
                print(f"⚠️ Could not parse {filename}")
        except Exception as e:
            print(f"❌ Error processing {filename}: {e}")

# 
# COMBINING ALL SDF FILES IN 1 FILE FOR EASIER LATER USE

combined_writer = Chem.SDWriter("all_ligands.sdf")

for file in os.listdir(output_dir):
    if file.endswith(".sdf"):
        filepath = os.path.join(output_dir, file)
        try:
            supplier = Chem.SDMolSupplier(filepath)
            mol = supplier[0] if supplier and len(supplier) > 0 else None
            if mol:
                combined_writer.write(mol)
                print(f"✅ Added: {file}")
            else:
                print(f"⚠️ Skipped (unreadable): {file}")
        except Exception as e:
            print(f"❌ Error loading {file}: {e}")

combined_writer.close()
print("✅ Combined SDF file written: all_ligands.sdf")
