from rdkit import Chem
import os

output_dir = "ligands"
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

