# LIGAND SDF FILES DOWNLOADER FROM RCSB & METADATA ENRICHER FOR RCSB.ORG

import os
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

# Hardcoded denticity dictionary
HARDCODED_DENTICITY = {
    "ACE": 1, "ACT": 1, "ADP": 3, "AMP": 3, "ATP": 4, "Acetate": 2, "Aniline": 1, "BA": 1, "BMA": 1,
    "BR": 1, "Br⁻": 1, "CA": 1, "CH₃⁻": 1, "CL": 1, "CN⁻": 1, "CO": 1, "COA": 4, "CS": 1, "CU": 1,
    "Cl⁻": 1, "C₂O₄²⁻": 2, "DMF": 1, "DMS": 1, "DMSO": 1, "EDO": 2, "EtOH": 1, "FAD": 4, "FE": 1,
    "FMN": 3, "FUC": 1, "Formate": 1, "F⁻": 1, "GAL": 1, "GDP": 2, "GLC": 1, "GMP": 2, "GTP": 3,
    "HEM": 4, "HEZ": 1, "H⁻": 1, "H₂O": 1, "IOD": 1, "Imidazole": 1, "I⁻": 1, "K": 1, "MAN": 1,
    "MES": 1, "MG": 1, "MN": 1, "MPD": 1, "MeCN": 1, "MeOH": 1, "NA": 1, "NAD": 3, "NAG": 1,
    "NCS⁻ (N-bound)": 1, "NH4": 1, "NH₃": 1, "NO3": 1, "NO₂⁻": 1, "NO₃⁻": 1, "OH⁻": 1, "PB": 1,
    "PEG": 1, "PGE": 1, "PO4": 1, "PPh₃": 1, "Pyrrolidine": 1, "SAH": 3, "SAM": 3, "SCN⁻ (S-bound)": 1,
    "SIA": 1, "SO4": 1, "SR": 1, "S²⁻": 1, "THF": 1, "TRS": 1, "ZN": 1, "bipy": 2, "en": 2, "phen": 2,
    "py": 1
}

# Describe denticity in words
def describe_denticity(n):
    names = {
        0: "non-dentate", 1: "monodentate", 2: "bidentate", 3: "tridentate",
        4: "tetradentate", 5: "pentadentate", 6: "hexadentate"
    }
    return names.get(n, f"{n}-dentate")

# Estimate denticity and Delta based on donor atoms (WE CHANGED TO HARDCODED DENTICIY SO THIS FUNCTION IS ONLY USED FOR DELTA OCT)
def estimate_denticity_and_delta(mol):
    smarts_rules = [
        ("[NX3;H2,H1;!$(NC=O)]", 23000),
        ("n", 24000),
        ("[O-]", 20000),
        ("C(=O)[O-]", 20000),
        ("[OH]", 20000),
        ("C#N", 30000),
        ("[S-]", 16000),
        ("S(=O)", 18000),
        ("[P]", 28000),
        ("[Cl,Br,I]", 12000),
        ("[C-]#[O+]", 32000),
    ]

    donor_atoms = set()
    matched_deltas = []

    for smarts, delta in smarts_rules:
        pattern = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            donor_atoms.update(match)
            matched_deltas.append(delta)

    denticity = len(donor_atoms)
    dent_desc = describe_denticity(denticity)

    if denticity > 0:
        avg_delta = int(sum(matched_deltas) / len(matched_deltas))
        return denticity, dent_desc, avg_delta, "estimated"
    else:
        return 0, "non-dentate", 20000, "default"

# Ligand input
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

# Download & enrich RCSB ligands
output_dir = "ligands"
os.makedirs(output_dir, exist_ok=True)

for lig_id in default_ligand_ids:
    url = f"https://files.rcsb.org/ligands/view/{lig_id}_ideal.sdf"
    response = requests.get(url)
    if response.status_code != 200:
        print(f"❌ Failed to download: {lig_id}")
        continue

    file_path = os.path.join(output_dir, f"{lig_id}.sdf")
    with open(file_path, "wb") as f:
        f.write(response.content)
    print(f"✅ Downloaded: {lig_id}")

    suppl = Chem.SDMolSupplier(file_path, removeHs=False)
    mol = suppl[0] if suppl and len(suppl) > 0 else None
    if mol is None:
        print(f"⚠️ Skipped unreadable file: {lig_id}")
        continue

    # Enrich molecule
    denticity = HARDCODED_DENTICITY.get(lig_id, 1)
    denticity_desc = describe_denticity(denticity)
    _, _, delta, delta_confidence = estimate_denticity_and_delta(mol)


    mol.SetProp("LigandName", lig_id)
    mol.SetProp("Denticity", str(denticity))
    mol.SetProp("DenticityDescription", denticity_desc)
    mol.SetProp("Delta_cm-1", str(delta))
    mol.SetProp("DeltaConfidence", delta_confidence)

    try:
        smiles = Chem.MolToSmiles(mol)
        mol.SetProp("SMILES", smiles)
        print(f"🔬 SMILES: {smiles}")
    except Exception as e:
        print(f"⚠️ Failed to generate SMILES for {lig_id}: {e}")

    writer = Chem.SDWriter(file_path)
    writer.write(mol)
    writer.close()
    print(f"🧪 Enriched: {lig_id} | Denticity = {denticity_desc}, Δ = {delta} cm⁻¹ ({delta_confidence})")


#CREATING SDF FILES FOR SPECHTROCHEMICAL SERIES LIGANDS THAT ARE NOT ON RCSB (with their spechtrochemical series info)  

import os
import re
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Geometry import Point3D

# Output directory
output_dir = "ligands"
os.makedirs(output_dir, exist_ok=True)

# Describe denticity
def describe_denticity(n):
    names = {
        0: "non-dentate",
        1: "monodentate",
        2: "bidentate",
        3: "tridentate",
        4: "tetradentate",
        5: "pentadentate",
        6: "hexadentate",
    }
    return names.get(n, f"{n}-dentate")

# Estimate denticity based on known donor SMARTS patterns
def estimate_denticity(mol):
    smarts_rules = [
        "[NX3;H2,H1;!$(NC=O)]",
        "n",
        "[O-]",
        "C(=O)[O-]",
        "[OH]",
        "C#N",
        "[S-]",
        "S(=O)",
        "[P]",
        "[Cl,Br,I]",
        "[C-]#[O+]"
    ]

    donor_atoms = set()
    for smarts in smarts_rules:
        pattern = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            donor_atoms.update(match)

    denticity = len(donor_atoms)
    return denticity, describe_denticity(denticity)

# Ligands: name -> (SMILES, field_strength, Delta_cm-1, donor_atom)
ligands = {
    "I⁻": ("[I-]", "very weak", 10000, "I"),
    "Br⁻": ("[Br-]", "very weak", 12000, "Br"),
    "S²⁻": ("[S-2]", "very weak", 11000, "S"),
    "SCN⁻ (S-bound)": ("S=C=N", "weak", 14000, "S"),
    "Cl⁻": ("[Cl-]", "weak", 15000, "Cl"),
    "NO₃⁻": ("O=N(=O)[O-]", "weak", 15000, "O"),
    "F⁻": ("[F-]", "weak", 18000, "F"),
    "OH⁻": ("[OH-]", "moderate", 19000, "O"),
    "H₂O": ("O", "moderate", 20000, "O"),
    "C₂O₄²⁻": ("O=C([O-])C(=O)[O-]", "moderate", 20000, "O,O"),
    "NCS⁻ (N-bound)": ("N=C=S", "moderate", 20000, "N"),
    "py": ("c1ccncc1", "moderate", 22000, "N"),
    "NH₃": ("N", "moderate-strong", 23000, "N"),
    "en": ("NCCN", "moderate-strong", 25000, "N,N"),
    "MeCN": ("CC#N", "moderate-strong", 24000, "N"),
    "DMF": ("CN(C)C=O", "moderate-strong", 22000, "O"),
    "DMSO": ("CS(=O)C", "moderate", 21000, "S"),
    "bipy": ("c1ccnc(c1)c2ccccn2", "strong", 26000, "N,N"),
    "phen": ("c1ccc2c(c1)ccnc2c3cccnc3", "strong", 27000, "N,N"),
    "NO₂⁻": ("[O-][N+](=O)", "strong", 27000, "N"),
    "PPh₃": ("P(c1ccccc1)(c2ccccc2)(c3ccccc3)", "strong", 28000, "P"),
    "CN⁻": ("[C-]#N", "very strong", 30000, "C"),
    "CO": ("[C-]#[O+]", "very strong", 32000, "C"),
    "H⁻": ("[H-]", "very strong", 31000, "H"),
    "CH₃⁻": ("[CH3-]", "very strong", 31000, "C"),
    "THF": ("C1CCOC1", "moderate", 21000, "O"),
    "EtOH": ("CCO", "moderate", 20000, "O"),
    "MeOH": ("CO", "moderate", 20000, "O"),
    "Acetate": ("CC(=O)[O-]", "moderate", 20000, "O,O"),
    "Formate": ("C(=O)[O-]", "moderate", 20000, "O"),
    "Pyrrolidine": ("C1CCNC1", "strong", 25000, "N"),
    "Aniline": ("c1ccccc1N", "moderate", 21000, "N"),
    "Imidazole": ("c1cnc[nH]1", "moderate-strong", 24000, "N")
}

# Clean names for filenames
def clean_filename(name):
    replacements = {
        "⁻": "minus", "⁺": "plus",
        "₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4",
        "₅": "5", "₆": "6", "₇": "7", "₈": "8", "₉": "9",
        "²": "2", "³": "3"
    }
    for uni, ascii_equiv in replacements.items():
        name = name.replace(uni, ascii_equiv)
    return re.sub(r"[^\w\-]", "_", name)

# Generate and write SDF files
for name, (smiles, field_strength, delta_cm1, donor_atom) in ligands.items():
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"❌ Could not parse: {name}")
        continue

    mol = Chem.AddHs(mol)

    if mol.GetNumAtoms() == 1:
        conf = Chem.Conformer(1)
        conf.SetAtomPosition(0, Point3D(4.0, 0.0, 0.0))
        mol.RemoveAllConformers()
        mol.AddConformer(conf)
        print(f"📍 Manually positioned monoatomic ligand: {name}")
    else:
        if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
            print(f"⚠️  3D embedding failed: {name}")
            continue
        AllChem.UFFOptimizeMolecule(mol)

    # Estimate denticity
    denticity = HARDCODED_DENTICITY.get(name, 1)
    denticity_desc = describe_denticity(denticity)


    # Add metadata
    mol.SetProp("LigandName", name)
    mol.SetProp("SMILES", smiles)
    mol.SetProp("FieldStrength", field_strength)
    mol.SetProp("Delta_cm-1", str(delta_cm1))
    mol.SetProp("DonorAtom", donor_atom)
    mol.SetProp("Denticity", str(denticity))
    mol.SetProp("DenticityDescription", denticity_desc)
    mol.SetProp("MolecularWeight", str(round(Descriptors.MolWt(mol), 2)))
    mol.SetProp("HDonors", str(Descriptors.NumHDonors(mol)))
    mol.SetProp("HAcceptors", str(Descriptors.NumHAcceptors(mol)))
    mol.SetProp("LogP", str(round(Descriptors.MolLogP(mol), 2)))
    mol.SetProp("NumAtoms", str(mol.GetNumAtoms()))

    # Write to SDF
    filename = os.path.join(output_dir, f"{clean_filename(name)}.sdf")
    writer = Chem.SDWriter(filename)
    writer.write(mol)
    writer.close()
    print(f"✅ SDF saved: {filename}")


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
