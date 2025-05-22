import matplotlib.pyplot as plt
from rdkit import Chem
from difflib import get_close_matches
from .metals_db import METALS


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
        total_charge = int(input("Enter total complex charge (e.g., 3 or 0): ").strip())
    except ValueError:
        print("❌ Invalid complex charge.")
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
    inferred_oxidation_state = total_charge - total_ligand_charge
    d_electrons = get_d_electron_count(metal["atomic_number"], inferred_oxidation_state)
    if d_electrons > 10 or d_electrons < 0:
        print(f"❌ Invalid d-electron count: {d_electrons}. This is not chemically possible.")
        return
    avg_delta = sum(deltas) / len(deltas)
    pairing_energy = get_pairing_energy(metal["block"])
    spin_state = predict_spin_state(avg_delta, pairing_energy)

    print(f"\n📊 Average Δ0: {avg_delta:.0f} cm⁻¹")
    print(f"Ligand total charge: {total_ligand_charge}")
    print(f"Inferred oxidation state of {metal['symbol']}: +{inferred_oxidation_state}")
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
