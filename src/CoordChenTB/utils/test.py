import matplotlib.pyplot as plt
from rdkit import Chem
from difflib import get_close_matches

# === Load ligand data from all_ligands.sdf ===
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
                ligand_info[name] = {"Delta": delta, "Charge": charge}
            except Exception:
                continue
    return ligand_info

# === Periodic table values for neutral atoms ===
electron_configurations = {
    'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25,
    'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
    'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45,
    'Pd': 46, 'Ag': 47, 'Cd': 48,
    'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78,
    'Au': 79, 'Hg': 80
}

# === Estimate d-electron count ===
def get_d_electron_count(metal, oxidation_state):
    total = electron_configurations.get(metal, 0)
    return total - 18 - oxidation_state

# === Estimate pairing energy (based on metal block) ===
def get_pairing_energy(metal):
    # 4d/5d metals typically have lower pairing energies
    if metal in ["Ru", "Rh", "Pd", "Ag", "Cd", "Mo", "Tc", "Re"]:
        return 13000
    elif metal in ["Os", "Ir", "Pt", "Au", "Hg", "W"]:
        return 12000
    elif metal in ["Fe", "Co", "Mn", "Cr", "V", "Ti"]:
        return 17000
    elif metal in ["Ni", "Cu", "Zn", "Sc"]:
        return 15000
    else:
        return 15000  # reasonable default

# === Fuzzy ligand name matching ===
def find_closest_ligand(name, ligand_data):
    if name in ligand_data:
        return name
    matches = get_close_matches(name, ligand_data.keys(), n=1, cutoff=0.6)
    return matches[0] if matches else None

# === Spin state prediction ===
def predict_spin_state(delta, metal, pairing_energy):
    if metal in ["Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
                 "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"]:
        return "Low-spin"
    return "High-spin" if delta < pairing_energy else "Low-spin"

# === d-orbital splitting diagram ===
def plot_d_orbital_splitting(d_electrons, delta, spin_state):
    fig, ax = plt.subplots()
    t2g = -delta / 2
    eg = delta / 2
    ax.hlines(t2g, 0.5, 1.5, colors='blue', label='t2g')
    ax.hlines(eg, 2.5, 3.5, colors='red', label='eg')
    ax.set_ylim(-delta, delta)
    ax.set_xlim(0, 4)
    ax.set_ylabel('Energy (cm⁻¹)')
    ax.set_title(f'd-Orbital Splitting Diagram\n{spin_state} configuration')
    ax.legend()
    plt.show()

# === Main ===
def main():
    ligand_data = load_ligand_data()
    if not ligand_data:
        print("❌ Could not load ligand data from all_ligands.sdf.")
        return

    metal = input("Enter metal center (e.g., Fe): ").strip().capitalize()
    try:
        total_charge = int(input("Enter total complex charge (e.g., 3 or 0): ").strip())
    except ValueError:
        print("❌ Invalid complex charge.")
        return

    ligands_input = input("Enter ligands (comma-separated): ").split(',')
    deltas = []
    ligand_charges = []

    for lig in ligands_input:
        lig = lig.strip()
        matched = find_closest_ligand(lig, ligand_data)
        if matched:
            info = ligand_data[matched]
            deltas.append(info["Delta"])
            ligand_charges.append(info["Charge"])
            print(f"🔍 Matched '{lig}' → '{matched}' (Δ = {info['Delta']}, Charge = {info['Charge']})")
        else:
            print(f"❌ Ligand not found: '{lig}' — using Δ = 20000, Charge = 0")
            deltas.append(20000)
            ligand_charges.append(0)

    if not deltas:
        print("❌ No valid ligands entered.")
        return

    total_ligand_charge = sum(ligand_charges)
    oxidation_state = total_charge - total_ligand_charge
    d_electrons = get_d_electron_count(metal, oxidation_state)
    avg_delta = sum(deltas) / len(deltas)
    pairing_energy = get_pairing_energy(metal)
    spin_state = predict_spin_state(avg_delta, metal, pairing_energy)

    print(f"\n📊 Average Δ₀: {avg_delta:.0f} cm⁻¹")
    print(f"Ligand total charge: {total_ligand_charge}")
    print(f"Inferred oxidation state of {metal}: +{oxidation_state}")
    print(f"d-electron count: {d_electrons}")
    print(f"Estimated pairing energy: {pairing_energy} cm⁻¹")
    print(f"Predicted spin state: {spin_state}")

    plot_d_orbital_splitting(d_electrons, avg_delta, spin_state)

if __name__ == "__main__":
    main()
