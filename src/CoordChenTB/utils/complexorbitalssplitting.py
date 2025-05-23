import matplotlib.pyplot as plt
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


def load_ligand_data(sdf_file="all_ligands.sdf"):
    ligand_info = {}
    suppl = Chem.SDMolSupplier(sdf_file)
    for mol in suppl:
        if mol is None:
            continue
        if mol.HasProp("LigandName") and mol.HasProp("Delta_cm-1") and mol.HasProp("NetFormalCharge"):
            try:
                name = mol.GetProp("LigandName").strip()
                delta = float(mol.GetProp("Delta_cm-1"))
                charge = int(float(mol.GetProp("NetFormalCharge")))
                confidence = mol.GetProp("DeltaConfidence").strip().lower() if mol.HasProp("DeltaConfidence") else ""
                ligand_info[name] = {"Delta": delta, "Charge": charge, "Confidence": confidence}
            except Exception:
                continue
    return ligand_info

def find_closest_ligand(name, ligand_data):
    matches = get_close_matches(name, ligand_data.keys(), n=5, cutoff=0.5)
    if not matches:
        return None

    print(f"\n🔍 Ligand '{name}' not found exactly. Did you mean one of the following?")
    for i, match in enumerate(matches, 1):
        print(f"{i}) {match} (Δ = {ligand_data[match]['Delta']} cm⁻¹, Charge = {ligand_data[match]['Charge']})")
    print("0) None of these")

    try:
        choice = int(input("→ Choose a number [0 to skip]: ").strip())
        if 1 <= choice <= len(matches):
            return matches[choice - 1]
    except ValueError:
        pass

    return None

def get_metal_data(symbol):
    return next((m for m in METALS if m["symbol"].lower() == symbol.lower()), None)

def get_d_electron_count(atomic_number, ox_state):
    return atomic_number - 18 - ox_state

def get_pairing_energy(block):
    return {"3d": 17000, "4d": 13000, "5d": 12000, "f": 11000}.get(block, 15000)

def predict_spin_state(delta, pairing_energy):
    return "Low-spin" if delta >= pairing_energy else "High-spin"

def plot_cf(d_elec, delta, spin, geom, geom_block, label=None):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 5))

    if geom == "octahedral":
        levels = {
            "d_xy": (0.5, -0.4 * delta),
            "d_xz": (1.5, -0.4 * delta),
            "d_yz": (2.5, -0.4 * delta),
            "d_z2": (1.0, 0.6 * delta),
            "d_x2-y2": (2.0, 0.6 * delta),
        }
        fill_order = ["d_xy", "d_xz", "d_yz", "d_z2", "d_x2-y2"]
        low_energy_orbitals = ["d_xy", "d_xz", "d_yz"]
        arrow_x = 3.3
        top = 0.6 * delta
        bottom = -0.4 * delta
        ax.annotate("", xy=(arrow_x, 0), xytext=(arrow_x, top), arrowprops=dict(arrowstyle='<->', lw=1.5, color='gray'))
        ax.annotate("", xy=(arrow_x, 0), xytext=(arrow_x, bottom), arrowprops=dict(arrowstyle='<->', lw=1.5, color='gray'))
        ax.text(arrow_x + 0.2, top / 2, "+3/5 Δ0", fontsize=9, va='center')
        ax.text(arrow_x + 0.2, bottom / 2, "–2/5 Δ0", fontsize=9, va='center')

    elif geom == "tetrahedral":
        Dt = 4 / 9 * delta
        levels = {
            "d_z2": (1.0, -0.6 * Dt),
            "d_x2-y2": (2.0, -0.6 * Dt),
            "d_xy": (0.5, 0.4 * Dt),
            "d_xz": (1.5, 0.4 * Dt),
            "d_yz": (2.5, 0.4 * Dt),
        }
        fill_order = ["d_z2", "d_x2-y2", "d_xy", "d_xz", "d_yz"]
        low_energy_orbitals = ["d_z2", "d_x2-y2"]
        arrow_x = 3.3
        top = 0.4 * Dt
        bottom = -0.6 * Dt
        ax.annotate("", xy=(arrow_x, 0), xytext=(arrow_x, top), arrowprops=dict(arrowstyle='<->', lw=1.5, color='gray'))
        ax.annotate("", xy=(arrow_x, 0), xytext=(arrow_x, bottom), arrowprops=dict(arrowstyle='<->', lw=1.5, color='gray'))
        ax.text(arrow_x + 0.2, top / 2, "+2/5 Δₜₐₓ", fontsize=9, va='center')
        ax.text(arrow_x + 0.2, bottom / 2, "–3/5 Δₜₐₓ", fontsize=9, va='center')

    elif geom == "square planar":
        base = delta
        levels = {
            "d_xz": (1.0, -0.5 * base),
            "d_yz": (2.0, -0.5 * base),
            "d_z2": (1.5, 0.0),
            "d_xy": (1.5, 0.5 * base),
            "d_x2-y2": (1.5, 1.0 * base),
        }
        fill_order = ["d_xz", "d_yz", "d_z2", "d_xy", "d_x2-y2"]
        low_energy_orbitals = ["d_xz", "d_yz", "d_z2"]

    else:
        print("❌ Invalid geometry.")
        return

    ax.axhline(0, linestyle='--', color='gray', linewidth=1)
    if d_elec <= len(low_energy_orbitals):
        spin = "Neutral-spin"

    for name, (x, y) in levels.items():
        ax.hlines(y, x - 0.3, x + 0.3, color='black')
        ax.text(x, y + 0.05 * abs(delta), name, ha='center', fontsize=8)

    slots = {orb: [] for orb in levels}
    filled = 0
    if spin in ["High-spin", "Neutral-spin"]:
        for orb in fill_order:
            if filled >= d_elec: break
            slots[orb].append("up"); filled += 1
        if spin == "High-spin":
            for orb in fill_order:
                if filled >= d_elec: break
                slots[orb].append("down"); filled += 1
    elif spin == "Low-spin":
        for orb in fill_order:
            if orb not in low_energy_orbitals: continue
            if filled >= d_elec: break
            slots[orb].append("up"); filled += 1
        for orb in fill_order:
            if orb not in low_energy_orbitals: continue
            if filled >= d_elec: break
            if len(slots[orb]) == 1:
                slots[orb].append("down"); filled += 1
        for orb in fill_order:
            if orb in low_energy_orbitals: continue
            if filled >= d_elec: break
            if len(slots[orb]) == 0:
                slots[orb].append("up"); filled += 1
        for orb in fill_order:
            if orb in low_energy_orbitals: continue
            if filled >= d_elec: break
            if len(slots[orb]) == 1:
                slots[orb].append("down"); filled += 1

    def draw_arrow(x, y, up=True):
        dy = abs(delta) * 0.05 if up else -abs(delta) * 0.05
        ax.arrow(x, y, 0, dy, head_width=0.08, head_length=abs(delta) * 0.02, fc='blue', ec='blue')

    for orb, electrons in slots.items():
        x, y = levels[orb]
        for i, spin_dir in enumerate(electrons):
            dx = (i - 0.5) * 0.15
            draw_arrow(x + dx, y, up=(spin_dir == "up"))

    cfse_cm1 = sum(len(e) * levels[orb][1] for orb, e in slots.items())
    pairing_energy_per_pair = get_pairing_energy(geom_block)
    actual_pairs = sum(1 for e in slots.values() if len(e) == 2)
    expected_pairs = max(0, (d_elec - 5) // 2) if d_elec > 5 else 0
    extra_pairs = max(0, actual_pairs - expected_pairs)
    pairing_penalty = extra_pairs * pairing_energy_per_pair
    total_cfse_cm1 = cfse_cm1 - pairing_penalty
    total_cfse_kjmol = total_cfse_cm1 * 0.01196

    print(f"\n📉 CFSE ({spin}): {cfse_cm1:.0f} cm⁻¹")
    print(f"🔁 Extra pairing penalty: {pairing_penalty:.0f} cm⁻¹ ({extra_pairs} extra pairs)")
    print(f"📉 Total CFSE: {total_cfse_cm1:.0f} cm⁻¹ (≈ {total_cfse_kjmol:.1f} kJ/mol)")

    ax.set_xlim(0, 4)
    ax.set_ylabel("Energy (cm⁻¹)")
    ax.set_xticks([])
    ax.set_title(label or f"{geom.capitalize()} field ({spin})")
    plt.tight_layout()
    plt.show()

def main():
    mode = input("\nMode?\n1) Ligand-based complex analysis\n2) Manual d-electron field visualization\n→ ").strip()
    if mode == "2":
        try:
            d_elec = int(input("Enter d-electron count (1–10): "))
            if not 1 <= d_elec <= 10:
                raise ValueError
        except ValueError:
            print("❌ Invalid number.")
            return
        geom_map = {"1": "octahedral", "2": "tetrahedral", "3": "square planar"}
        choice = input("Geometry? 1) Octahedral  2) Tetrahedral  3) Square planar: ").strip()
        geom = geom_map.get(choice)
        if not geom:
            print("❌ Invalid geometry.")
            return
        delta_in = input("Use default Δ0 = 20000 cm⁻¹? (Y/n): ").strip().lower()
        delta = 20000 if delta_in in ["", "y"] else float(input("Enter Δ0 in cm⁻¹: "))
        geom_block = "3d"
        print("\n🔷 HIGH-SPIN CONFIGURATION:")
        plot_cf(d_elec, delta, "High-spin", geom, geom_block, label=f"High-spin {geom}")
        print("\n🔷 LOW-SPIN CONFIGURATION:")
        plot_cf(d_elec, delta, "Low-spin", geom, geom_block, label=f"Low-spin {geom}")
        return

    ligand_data = load_ligand_data()
    if not ligand_data:
        print("❌ Could not load ligand data from all_ligands.sdf.")
        return

    print("\nAvailable metals:")
    print(", ".join(sorted(m["symbol"] for m in METALS)))
    metal_symbol = input("Enter metal center (e.g., Fe): ").strip()
    metal = get_metal_data(metal_symbol)
    if not metal:
        print(f"❌ Metal not found: {metal_symbol}")
        return

    print(f"✅ Found metal: {metal['symbol']} (block: {metal['block']})")
    print(f"Preferred oxidation states: {metal['oxidation_states']}")

    try:
        oxidation_state = int(input(f"Enter oxidation state of {metal['symbol']} (e.g., 2 or 3): ").strip())
    except ValueError:
        print("❌ Invalid oxidation state.")
        return


    print("\nAvailable ligands:")
    print(", ".join(sorted(ligand_data.keys())))
    ligands_input = input("Enter ligands (comma-separated): ").split(',')
    deltas = []
    ligand_charges = []

    for lig in ligands_input:
        lig = lig.strip()
        matched = find_closest_ligand(lig, ligand_data)
        if matched:
            info = ligand_data[matched]
            note = " (⚠️  warning: default Δ)" if info.get("Confidence") == "default" else ""
            deltas.append(info["Delta"])
            ligand_charges.append(info["Charge"])
            print(f"🔍 Matched '{lig}' → '{matched}' (Δ = {info['Delta']}{note}, Charge = {info['Charge']})")
        else:
            print(f"❌ Ligand not found: '{lig}' — using Δ = 20000, Charge = 0")
            deltas.append(20000)
            ligand_charges.append(0)

    total_ligand_charge = sum(ligand_charges)
    total_charge = oxidation_state + total_ligand_charge
    d_electrons = get_d_electron_count(metal["atomic_number"], oxidation_state)
    if d_electrons > 10 or d_electrons < 0:
        print(f"❌ Invalid d-electron count: {d_electrons}. This is not chemically possible.")
        return
    avg_delta = sum(deltas) / len(deltas)
    pairing_energy = get_pairing_energy(metal["block"])
    spin_state = predict_spin_state(avg_delta, pairing_energy)

    print(f"\n📊 Average Δ0: {avg_delta:.0f} cm⁻¹")
    print(f"Ligand total charge: {total_ligand_charge}")
    print(f"Given oxidation state of {metal['symbol']}: +{oxidation_state}")
    print(f"🔋 Total complex charge: {total_charge:+}")
    print(f"d-electron count: {d_electrons}")
    print(f"Estimated pairing energy: {pairing_energy} cm⁻¹")
    print(f"Predicted spin state: {spin_state}")


    geom_map = {"1": "octahedral", "2": "tetrahedral", "3": "square planar"}
    choice = input("Geometry? 1) Octahedral  2) Tetrahedral  3) Square planar: ").strip()
    geom = geom_map.get(choice)
    if not geom:
        print("❌ Invalid geometry choice.")
        return

    plot_cf(d_electrons, avg_delta, spin_state, geom, metal["block"])

if __name__ == "__main__":
    main()
