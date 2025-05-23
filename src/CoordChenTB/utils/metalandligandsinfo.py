import os
from rdkit import Chem
from difflib import get_close_matches


METALS = [
    {'symbol': 'Sc',  'atomic_number': 21,  'block': '3d', 'smiles': '[Sc]',  'atomic_weight': 44.956,  'oxidation_states': [3]},
    {'symbol': 'Ti',  'atomic_number': 22,  'block': '3d', 'smiles': '[Ti]',  'atomic_weight': 47.867,  'oxidation_states': [2, 3, 4]},
    {'symbol': 'V',   'atomic_number': 23,  'block': '3d', 'smiles': '[V]',   'atomic_weight': 50.942,  'oxidation_states': [2, 3, 4, 5]},
    {'symbol': 'Cr',  'atomic_number': 24,  'block': '3d', 'smiles': '[Cr]',  'atomic_weight': 51.996,  'oxidation_states': [2, 3, 6]},
    {'symbol': 'Mn',  'atomic_number': 25,  'block': '3d', 'smiles': '[Mn]',  'atomic_weight': 54.938,  'oxidation_states': [2, 4, 7]},
    {'symbol': 'Fe',  'atomic_number': 26,  'block': '3d', 'smiles': '[Fe]',  'atomic_weight': 55.845,  'oxidation_states': [2, 3]},
    {'symbol': 'Co',  'atomic_number': 27,  'block': '3d', 'smiles': '[Co]',  'atomic_weight': 58.933,  'oxidation_states': [2, 3]},
    {'symbol': 'Ni',  'atomic_number': 28,  'block': '3d', 'smiles': '[Ni]',  'atomic_weight': 58.693,  'oxidation_states': [2, 3]},
    {'symbol': 'Cu',  'atomic_number': 29,  'block': '3d', 'smiles': '[Cu]',  'atomic_weight': 63.546,  'oxidation_states': [1, 2]},
    {'symbol': 'Zn',  'atomic_number': 30,  'block': '3d', 'smiles': '[Zn]',  'atomic_weight': 65.38,   'oxidation_states': [2]},

    {'symbol': 'Y',   'atomic_number': 39,  'block': '4d', 'smiles': '[Y]',   'atomic_weight': 88.906,  'oxidation_states': [3]},
    {'symbol': 'Zr',  'atomic_number': 40,  'block': '4d', 'smiles': '[Zr]',  'atomic_weight': 91.224,  'oxidation_states': [2, 3, 4]},
    {'symbol': 'Nb',  'atomic_number': 41,  'block': '4d', 'smiles': '[Nb]',  'atomic_weight': 92.906,  'oxidation_states': [3, 5]},
    {'symbol': 'Mo',  'atomic_number': 42,  'block': '4d', 'smiles': '[Mo]',  'atomic_weight': 95.95,   'oxidation_states': [2, 3, 4, 6]},
    {'symbol': 'Ru',  'atomic_number': 44,  'block': '4d', 'smiles': '[Ru]',  'atomic_weight': 101.07,  'oxidation_states': [2, 3, 4, 6, 8]},
    {'symbol': 'Rh',  'atomic_number': 45,  'block': '4d', 'smiles': '[Rh]',  'atomic_weight': 102.906, 'oxidation_states': [3]},
    {'symbol': 'Pd',  'atomic_number': 46,  'block': '4d', 'smiles': '[Pd]',  'atomic_weight': 106.42,  'oxidation_states': [2, 4]},
    {'symbol': 'Ag',  'atomic_number': 47,  'block': '4d', 'smiles': '[Ag]',  'atomic_weight': 107.868, 'oxidation_states': [1]},
    {'symbol': 'Cd',  'atomic_number': 48,  'block': '4d', 'smiles': '[Cd]',  'atomic_weight': 112.414, 'oxidation_states': [2]},

    {'symbol': 'Hf',  'atomic_number': 72,  'block': '5d', 'smiles': '[Hf]',  'atomic_weight': 178.49,  'oxidation_states': [4]},
    {'symbol': 'W',   'atomic_number': 74,  'block': '5d', 'smiles': '[W]',   'atomic_weight': 183.84,  'oxidation_states': [2, 3, 4, 5, 6]},
    {'symbol': 'Re',  'atomic_number': 75,  'block': '5d', 'smiles': '[Re]',  'atomic_weight': 186.207, 'oxidation_states': [4, 7]},
    {'symbol': 'Os',  'atomic_number': 76,  'block': '5d', 'smiles': '[Os]',  'atomic_weight': 190.23,  'oxidation_states': [2, 3, 4, 6, 8]},
    {'symbol': 'Ir',  'atomic_number': 77,  'block': '5d', 'smiles': '[Ir]',  'atomic_weight': 192.217, 'oxidation_states': [3, 4, 6]},
    {'symbol': 'Pt',  'atomic_number': 78,  'block': '5d', 'smiles': '[Pt]',  'atomic_weight': 195.084, 'oxidation_states': [2, 4]},
    {'symbol': 'Au',  'atomic_number': 79,  'block': '5d', 'smiles': '[Au]',  'atomic_weight': 196.967, 'oxidation_states': [1, 3]},
    {'symbol': 'Hg',  'atomic_number': 80,  'block': '5d', 'smiles': '[Hg]',  'atomic_weight': 200.592, 'oxidation_states': [1, 2]},

    {'symbol': 'La',  'atomic_number': 57,  'block': 'f',  'smiles': '[La]',  'atomic_weight': 138.905, 'oxidation_states': [3]},
    {'symbol': 'Ce',  'atomic_number': 58,  'block': 'f',  'smiles': '[Ce]',  'atomic_weight': 140.116, 'oxidation_states': [3, 4]},
    {'symbol': 'Sm',  'atomic_number': 62,  'block': 'f',  'smiles': '[Sm]',  'atomic_weight': 150.36,  'oxidation_states': [2, 3]},
    {'symbol': 'Eu',  'atomic_number': 63,  'block': 'f',  'smiles': '[Eu]',  'atomic_weight': 151.964, 'oxidation_states': [2, 3]},
    {'symbol': 'Gd',  'atomic_number': 64,  'block': 'f',  'smiles': '[Gd]',  'atomic_weight': 157.25,  'oxidation_states': [3]},
    {'symbol': 'Tb',  'atomic_number': 65,  'block': 'f',  'smiles': '[Tb]',  'atomic_weight': 158.925, 'oxidation_states': [3, 4]},
    {'symbol': 'Dy',  'atomic_number': 66,  'block': 'f',  'smiles': '[Dy]',  'atomic_weight': 162.5,   'oxidation_states': [3]},
    {'symbol': 'Er',  'atomic_number': 68,  'block': 'f',  'smiles': '[Er]',  'atomic_weight': 167.259, 'oxidation_states': [3]},
    {'symbol': 'Yb',  'atomic_number': 70,  'block': 'f',  'smiles': '[Yb]',  'atomic_weight': 173.045, 'oxidation_states': [2, 3]},
    {'symbol': 'Th',  'atomic_number': 90,  'block': 'f',  'smiles': '[Th]',  'atomic_weight': 232.038, 'oxidation_states': [4]},
    {'symbol': 'U',   'atomic_number': 92,  'block': 'f',  'smiles': '[U]',   'atomic_weight': 238.029, 'oxidation_states': [3, 4, 5, 6]},
]



def load_ligand_data(sdf_path):
    ligands = {}
    suppl = Chem.SDMolSupplier(sdf_path)
    for mol in suppl:
        if mol is None:
            continue
        name = mol.GetProp("LigandName") if mol.HasProp("LigandName") else None
        if name:
            ligands[name] = mol
    return ligands


def find_closest_ligand(name, ligands):
    matches = get_close_matches(name, ligands.keys(), n=5, cutoff=0.5)
    if not matches:
        return None
    print(f"\n🔍 Ligand '{name}' not found. Did you mean one of these?")
    for i, m in enumerate(matches, 1):
        print(f"{i}) {m}")
    print("0) None of these")
    try:
        choice = int(input("→ Choose a number [0 to skip]: ").strip())
        if 1 <= choice <= len(matches):
            return matches[choice - 1]
    except Exception:
        pass
    return None


def print_ligand_info(mol, name):
    print(f"\n📘 Ligand: {name}\n{'-'*40}")
    properties_to_show = [
        ("Denticity", "Denticity"),
        ("DenticityDescription", "Denticity Description"),
        ("Delta_cm-1", "Δ₀ (cm⁻¹)"),
        ("DeltaConfidence", "Δ₀ Confidence"),
        ("FieldStrength", "Field Strength"),
        ("LogP", "logP"),
        ("NetFormalCharge", "Charge"),
        ("HDonors", "H-Donors"),
        ("HAcceptors", "H-Acceptors"),
        ("NumAtoms", "Atom Count"),
        ("MolecularWeight", "Mol. Weight"),
        ("SMILES", "SMILES")
    ]

    for prop_key, label in properties_to_show:
        if mol.HasProp(prop_key):
            print(f"{label:20}: {mol.GetProp(prop_key)}")


def print_metal_info(symbol):
    print(f"\n🧲 Metal: {symbol}\n{'-'*40}")
    symbol_upper = symbol.upper()
    for entry in METALS:
        if entry.get("symbol", "").upper() == symbol_upper:
            for key, val in entry.items():
                print(f"{key:20}: {val}")
            return
    print("❌ Metal not found in database.")


if __name__ == "__main__":
    path = os.path.abspath("all_ligands.sdf")
    if not os.path.exists(path):
        print("❌ all_ligands.sdf not found.")
        exit()

    ligands = load_ligand_data(path)
    print("✅ Loaded ligands:", ", ".join(sorted(ligands.keys())[:10]), "...")

    while True:
        name = input("\nEnter ligand name or metal symbol (or 'exit'): ").strip()
        if name.lower() in {"exit", "quit"}:
            break

        if name in ligands:
            print_ligand_info(ligands[name], name)
        else:
            # If not a ligand, try as metal
            found = any(m.get("symbol", "").upper() == name.upper() for m in METALS)
            if found:
                print_metal_info(name)
            else:
                suggestion = find_closest_ligand(name, ligands)
                if suggestion:
                    print_ligand_info(ligands[suggestion], suggestion)
                else:
                    print("❌ Name not recognized as ligand or metal.")
