import os
import re
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

# Output folder
output_dir = "ligands"
os.makedirs(output_dir, exist_ok=True)

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

# Helper function to clean names for filenames
def clean_filename(name):
    replacements = {
        "⁻": "minus", "⁺": "plus",
        "₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4",
        "₅": "5", "₆": "6", "₇": "7", "₈": "8", "₉": "9",
        "²": "2", "³": "3"
    }
    for uni, ascii_equiv in replacements.items():
        name = name.replace(uni, ascii_equiv)
    name = re.sub(r"[^\w\-]", "_", name)  # Remove non-safe characters
    return name

# Generate and write SDF files
for name, (smiles, field_strength) in ligands.items():
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"❌ Could not parse: {name}")
        continue

    # Add hydrogens and generate 3D geometry
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
        print(f"⚠️  3D embedding failed: {name}")
        continue
    AllChem.UFFOptimizeMolecule(mol)

    # Set metadata fields
    mol.SetProp("LigandName", name)
    mol.SetProp("SMILES", smiles)
    mol.SetProp("FieldStrength", field_strength)
    mol.SetProp("MolecularWeight", str(round(Descriptors.MolWt(mol), 2)))
    mol.SetProp("HDonors", str(Descriptors.NumHDonors(mol)))
    mol.SetProp("HAcceptors", str(Descriptors.NumHAcceptors(mol)))
    mol.SetProp("LogP", str(round(Descriptors.MolLogP(mol), 2)))
    mol.SetProp("NumAtoms", str(mol.GetNumAtoms()))

    # Save with cleaned filename
    clean_name = clean_filename(name)
    sdf_path = os.path.join(output_dir, f"{clean_name}.sdf")

    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)
    writer.close()
    print(f"✅ SDF saved: {sdf_path}")
