import os
from rdkit import Chem

input_dir = "ligands"

for filename in os.listdir(input_dir):
    if filename.endswith(".sdf"):
        file_path = os.path.join(input_dir, filename)
        try:
            suppl = Chem.SDMolSupplier(file_path, removeHs=False)
            mol = suppl[0] if suppl and len(suppl) > 0 else None
            if mol:
                # Compute net charge
                net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
                mol.SetProp("NetFormalCharge", str(net_charge))

                # Overwrite original SDF file with updated molecule
                writer = Chem.SDWriter(file_path)
                writer.write(mol)
                writer.close()
                print(f"✅ Updated {filename} with NetFormalCharge = {net_charge}")
            else:
                print(f"⚠️ Could not parse {filename}")
        except Exception as e:
            print(f"❌ Error processing {filename}: {e}")
