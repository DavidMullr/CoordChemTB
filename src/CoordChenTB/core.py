"""
CoordChemTB Core Controller
A CLI for generating ligands, building coordination complexes,
and analyzing ligand field properties.
"""
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from CoordChenTB.utils import complexdrawerfinal as drawer
from CoordChenTB.utils import complexorbitalssplitting as orbital
from CoordChenTB.utils import metalandligandsinfo as info
# Not importing generalligands yet since it runs automatically on import in current version


def analyze_complex():
    print("\n🔬 Ligand Field Analysis")
    orbital.main()


def build_and_draw_complex():
    print("\n🎨 Build Coordination Complex")
    metal = input("Enter metal symbol (e.g., Fe): ").strip()
    geometry = input("Geometry (octahedral/tetrahedral/square_planar): ").strip()

    print("Enter ligands and counts (type 'done' to finish):")
    ligand_counts = {}
    while True:
        name = input("Ligand name: ").strip()
        if name.lower() == 'done':
            break
        try:
            count = int(input(f"How many of '{name}'? "))
            ligand_counts[name] = count
        except ValueError:
            print("❌ Invalid count. Try again.")

    try:
        # Automatically name the output image
        output_filename = f"{metal}_{geometry}_{'_'.join(f'{k}{v}' for k, v in ligand_counts.items())}.png"

        drawer.create_complex_from_ligand_dict(
            metal=metal,
            ligand_counts=ligand_counts,
            geometry=geometry,
            carbon_angles=drawer.CARBON_ANGLE_SELECTION,
            output_file=output_filename
        )
        print(f"✅ Image saved to {output_filename}")

    except Exception as e:
        print(f"❌ Failed to build complex: {e}")



def lookup_ligand():
    print("\n🔍 Ligand Info Lookup")
    sdf_path = os.path.abspath("all_ligands.sdf")
    if not os.path.exists(sdf_path):
        print("❌ all_ligands.sdf not found. Please generate ligands first.")
        return

    ligands = info.load_ligand_data(sdf_path)
    name = input("Enter ligand name: ").strip()
    if name not in ligands:
        suggested = info.find_closest_ligand(name, ligands)
        if not suggested:
            print("❌ Ligand not recognized.")
            return
        name = suggested
    mol = ligands[name]
    info.print_ligand_info(mol, name)


def main():
    while True:
        print("""
CoordChemTB CLI
===============
1) Analyze ligand field / spin state
2) Build & draw coordination complex
3) Lookup ligand info
4) Exit
""")
        choice = input("→ Choose an option: ").strip()
        if choice == '1':
            analyze_complex()
        elif choice == '2':
            build_and_draw_complex()
        elif choice == '3':
            lookup_ligand()
        elif choice == '4':
            print("👋 Goodbye!")
            break
        else:
            print("❌ Invalid choice. Try again.")


if __name__ == "__main__":
    main()
