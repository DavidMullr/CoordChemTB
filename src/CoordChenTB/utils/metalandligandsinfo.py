import os
import sys
from rdkit import Chem
from difflib import get_close_matches

# Make sure imports from CoordChenTB work
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from CoordChenTB.utils.metals_db import METALS


def load_ligand_data(sdf_path):
    """Returns a dict of ligand names to RDKit molecules."""
    ligands = {}
    suppl = Chem.SDMolSupplier(sdf_path)
    for mol in suppl:
        if mol and mol.HasProp("LigandName"):
            name = mol.GetProp("LigandName").strip()
            ligands[name] = mol
    return ligands


def find_closest_ligand(name, ligand_dict):
    """Suggests similar ligand names and lets user choose one."""
    matches = get_close_matches(name, ligand_dict.keys(), n=5, cutoff=0.5)
    if not matches:
        return None

    print(f"\n🔍 Ligand '{name}' not found exactly. Did you mean one of the following?")
    for i, match in enumerate(matches, 1):
        mol = ligand_dict[match]
        delta = mol.GetProp("Delta_cm-1") if mol.HasProp("Delta_cm-1") else "?"
        charge = mol.GetProp("NetFormalCharge") if mol.HasProp("NetFormalCharge") else "?"
        print(f"  {i}) {match} (Δ = {delta} cm⁻¹, Charge = {charge})")
    print("  0) None of these")

    try:
        choice = int(input("→ Choose a number [0 to skip]: ").strip())
        if 1 <= choice <= len(matches):
            return matches[choice - 1]
    except ValueError:
        pass

    return None


def print_ligand_info(mol, name):
    print(f"\n✅ Info for ligand: {name}")
    if mol.HasProp("Delta_cm-1"):
        print(f"Δ (cm⁻¹): {mol.GetProp('Delta_cm-1')}")
    if mol.HasProp("NetFormalCharge"):
        print(f"Charge: {mol.GetProp('NetFormalCharge')}")
    if mol.HasProp("DENTATE"):
        print(f"Denticity: {mol.GetProp('DENTATE')}")
    if mol.HasProp("DeltaConfidence"):
        print(f"Confidence: {mol.GetProp('DeltaConfidence')}")


if __name__ == "__main__":
    # Path to the ligand SDF file at project root
    sdf_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "all_ligands.sdf"))

    if not os.path.isfile(sdf_path):
        print(f"❌ Ligand file not found at {sdf_path}")
        sys.exit(1)

    ligand_data = load_ligand_data(sdf_path)

    if not ligand_data:
        print("❌ No ligands loaded.")
        sys.exit(1)

    print("✅ Loaded ligands:")
    print(", ".join(sorted(ligand_data.keys())))

    while True:
        ligand = input("\nEnter a ligand name (or 'exit'): ").strip()
        if ligand.lower() == 'exit':
            break

        if ligand not in ligand_data:
            suggested = find_closest_ligand(ligand, ligand_data)
            if not suggested:
                print("❌ Ligand not recognized.")
                continue
            ligand = suggested

        mol = ligand_data[ligand]
        print_ligand_info(mol, ligand)
