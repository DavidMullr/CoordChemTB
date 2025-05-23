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

def analyze_complex():
    print("\n🔬 Ligand Field Analysis")
    orbital.main()


def build_and_draw_complex():
    import traceback

    print("\n🎨 Build Coordination Complex")
    metal = input("Enter metal symbol (e.g., Fe): ").strip()
    try:
        oxidation_state = int(input(f"Enter oxidation state for {metal} (default is 2): ").strip() or "2")
    except ValueError:
        oxidation_state = 2

    print("Choose geometry:")
    print("1) Octahedral")
    print("2) Tetrahedral")
    print("3) Square planar")
    geometry_input = input("Enter number [1-3]: ").strip()
    geometry_map = {"1": "octahedral", "2": "tetrahedral", "3": "square_planar"}
    geometry = geometry_map.get(geometry_input)

    if not geometry:
        print("❌ Invalid geometry selection.")
        return

    try:
        bond_length = float(input("Enter bond length (useful for bulky ligands) (default is 1): ").strip() or "1")
    except ValueError:
        print("❌ Invalid bond length. Using default of 1.")
        bond_length = 1.0

    print("Enter ligands and counts (type 'done' to finish):")
    ligand_counts = {}
    lig_map = drawer.load_ligands_from_folder(drawer.LIGAND_FOLDER)

    while True:
        name = input("Ligand name (or 'done' to finish): ").strip()
        if name.lower() == 'done':
            break

        if name not in lig_map:
            suggestion = info.find_closest_ligand(name, lig_map)
            if suggestion:
                print(f"✔ Using closest match: {suggestion}")
                name = suggestion
            else:
                print("❌ Ligand not recognized and no close matches found.")
                continue

        try:
            count = int(input(f"How many of '{name}'? "))
            ligand_counts[name] = count
        except ValueError:
            print("❌ Invalid count. Try again.")



    if not ligand_counts:
        print("❌ No ligands provided. Aborting.")
        return

    try:
        output_filename = f"{metal}_{geometry}_{'_'.join(f'{k}{v}' for k, v in ligand_counts.items())}.png"

        metal_with_charge = f"{metal}+{oxidation_state}" if oxidation_state != 0 else metal

        drawer.create_complex_from_ligand_dict(
            metal=metal_with_charge,
            ligand_counts=ligand_counts,
            geometry=geometry,
            carbon_angles=drawer.CARBON_ANGLE_SELECTION,
            output_file=output_filename,
            bond_length=bond_length
        )

        print(f"✅ Image saved to {output_filename}")

    except Exception as e:
        traceback.print_exc()
        print(f"❌ Failed to build complex: {e}")




def lookup_ligand_or_metal():
    print("\n🔍 Ligand or Metal Info Lookup")
    sdf_path = os.path.abspath("all_ligands.sdf")
    if not os.path.exists(sdf_path):
        print("❌ all_ligands.sdf not found. Please generate ligands first.")
        return

    ligands = info.load_ligand_data(sdf_path)
    name = input("Enter ligand name or metal symbol: ").strip()
    if name in ligands:
        info.print_ligand_info(ligands[name], name)
    else:
        found = any(m.get("symbol", "").upper() == name.upper() for m in METALS)
        if found:
            info.print_metal_info(name)
        else:
            suggestion = info.find_closest_ligand(name, ligands)
            if suggestion:
                info.print_ligand_info(ligands[suggestion], suggestion)
            else:
                print("❌ Name not recognized as ligand or metal.")

####DOESN'T WORK PLEASE LAUNCH DIRECLY FRON INTERFACE2.py

def launch_gui():
    try:
        from CoordChenTB.utils.interfacefinalversion import CoordinationGUI
        app = CoordinationGUI()
        app.mainloop()
    except Exception as e:
        print(f"❌ Failed to launch GUI: {e}")



def main():
    while True:
        print("""
CoordChemTB CLI
===============
1) Analyze ligand field / spin state
2) Build & draw coordination complex
3) Lookup ligand/metal info
4) TO LAUNCH GUI PLEASE launch interfacefinalversion.py directly
5) Exit
""")
        choice = input("→ Choose an option: ").strip()
        if choice == '1':
            analyze_complex()
        elif choice == '2':
            build_and_draw_complex()
        elif choice == '3':
            lookup_ligand_or_metal()
        elif choice == '4':
            launch_gui()
        elif choice == '5':
            print("👋 Goodbye!")
            break
        else:
            print("❌ Invalid choice. Try again.")


if __name__ == "__main__":
    main()
