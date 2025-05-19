# Dictionary of metals commonly used in coordination compounds
metal_data = {
    'Sc': {'atomic_weight': 44.956,  'oxidation_states': [3]},
    'Ti': {'atomic_weight': 47.867,  'oxidation_states': [2, 3, 4]},
    'V':  {'atomic_weight': 50.942,  'oxidation_states': [2, 3, 4, 5]},
    'Cr': {'atomic_weight': 51.996,  'oxidation_states': [2, 3, 6]},
    'Mn': {'atomic_weight': 54.938,  'oxidation_states': [2, 4, 7]},
    'Fe': {'atomic_weight': 55.845,  'oxidation_states': [2, 3]},
    'Co': {'atomic_weight': 58.933,  'oxidation_states': [2, 3]},
    'Ni': {'atomic_weight': 58.693,  'oxidation_states': [2, 3]},
    'Cu': {'atomic_weight': 63.546,  'oxidation_states': [1, 2]},
    'Zn': {'atomic_weight': 65.38,   'oxidation_states': [2]},
    'Y':  {'atomic_weight': 88.906,  'oxidation_states': [3]},
    'Zr': {'atomic_weight': 91.224,  'oxidation_states': [2, 3, 4]},
    'Nb': {'atomic_weight': 92.906,  'oxidation_states': [3, 5]},
    'Mo': {'atomic_weight': 95.95,   'oxidation_states': [2, 3, 4, 6]},
    'Ru': {'atomic_weight': 101.07,  'oxidation_states': [2, 3, 4, 6, 8]},
    'Rh': {'atomic_weight': 102.906, 'oxidation_states': [3]},
    'Pd': {'atomic_weight': 106.42,  'oxidation_states': [2, 4]},
    'Ag': {'atomic_weight': 107.868, 'oxidation_states': [1]},
    'Cd': {'atomic_weight': 112.414, 'oxidation_states': [2]},
    'Hf': {'atomic_weight': 178.49,  'oxidation_states': [4]},
    'W':  {'atomic_weight': 183.84,  'oxidation_states': [2, 3, 4, 5, 6]},
    'Re': {'atomic_weight': 186.207, 'oxidation_states': [4, 7]},
    'Os': {'atomic_weight': 190.23,  'oxidation_states': [2, 3, 4, 6, 8]},
    'Ir': {'atomic_weight': 192.217, 'oxidation_states': [3, 4, 6]},
    'Pt': {'atomic_weight': 195.084, 'oxidation_states': [2, 4]},
    'Au': {'atomic_weight': 196.967, 'oxidation_states': [1, 3]},
    'Hg': {'atomic_weight': 200.592, 'oxidation_states': [1, 2]},
    # Lanthanides frequently used as coordination centers
    'La': {'atomic_weight': 138.905, 'oxidation_states': [3]},
    'Ce': {'atomic_weight': 140.116, 'oxidation_states': [3, 4]},
    'Sm': {'atomic_weight': 150.36,  'oxidation_states': [2, 3]},
    'Eu': {'atomic_weight': 151.964, 'oxidation_states': [2, 3]},
    'Gd': {'atomic_weight': 157.25,  'oxidation_states': [3]},
    'Tb': {'atomic_weight': 158.925, 'oxidation_states': [3, 4]},
    'Dy': {'atomic_weight': 162.5,   'oxidation_states': [3]},
    'Er': {'atomic_weight': 167.259, 'oxidation_states': [3]},
    'Yb': {'atomic_weight': 173.045, 'oxidation_states': [2, 3]},
    # Actinides often in coordination compounds
    'Th': {'atomic_weight': 232.038, 'oxidation_states': [4]},
    'U':  {'atomic_weight': 238.029, 'oxidation_states': [3, 4, 5, 6]}
}

# Build the METALS list
METALS = []
for symbol, props in metal_data.items():
    METALS.append({
        'symbol': symbol,
        'smiles': f'[{symbol}]',
        'atomic_weight': props['atomic_weight'],
        'oxidation_states': props['oxidation_states']
    })

# Emit the Python module
with open('metals_db.py', 'w') as out:
    out.write('# Auto-generated metals database for coordination chemistry\n')
    out.write('# Contains METALS: list of dicts with symbol, smiles, atomic_weight, oxidation_states\n\n')
    out.write('METALS = [\n')
    for m in METALS:
        out.write(f'    {m!r},\n')
    out.write(']\n')

print(f"Generated metals_db.py with {len(METALS)} coordination-relevant metals.")
