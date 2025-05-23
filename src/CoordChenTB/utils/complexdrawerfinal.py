import os
import io
import math
import importlib.machinery
import pytest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdchem
from rdkit.Geometry import Point2D, Point3D
from rdkit.Chem.Draw import MolToFile
from PIL import Image

# Path to metals database
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from CoordChenTB.utils.metals_db import METALS

# Octahedral 2D positions for donor placement
OCTAHEDRAL_POS = [
    Point2D(0.0, 3.0),    # top (0)
    Point2D(0.0, -3.0),   # bottom (1)
    Point2D(-2.6, 1.5),   # front-left (2)
    Point2D(2.6, 1.5),    # front-right (3)
    Point2D(-2.6, -1.5),  # back-left (4)
    Point2D(2.6, -1.5),   # back-right (5)
]

# Define 2D donor positions for each geometry
GEOMETRY_POSITIONS = {
    'octahedral': [
        Point2D(0.0, 3.0),
        Point2D(0.0, -3.0),
        Point2D(-2.6, 1.5),
        Point2D(2.6, 1.5),
        Point2D(-2.6, -1.5),
        Point2D(2.6, -1.5),
    ],
    'tetrahedral': [
        Point2D(1.5, 1.5),
        Point2D(-1.5, 1.5),
        Point2D(-1.5, -1.5),
        Point2D(1.5, -1.5),
    ],
    'square_planar': [
        Point2D(-2.0, 0.0),
        Point2D(0.0, 2.0),
        Point2D(2.0, 0.0),
        Point2D(0.0, -2.0),
    ],
}

DEFAULT_CARBON_ANGLES_BY_GEOMETRY = {
    'octahedral': [120, 240],
    'tetrahedral': [60, 180],
    'square_planar': [90, 270],
}


# Default angles for explicit carbon placement around each donor atom
DEFAULT_CARBON_ANGLES = [120, 240]

ALLOWED_ANGLES = {30, 90, 150, 210, 270, 330}

def validate_carbon_angles(angles_dict):
    """Validate that all carbon angles are in the allowed set (30, 90, 150, 210, 270, 330 degrees).
    
    Args:
        angles_dict: Dictionary mapping site indices to angle lists
        
    Returns:
        Validated dictionary with all angles rounded to nearest allowed value
        
    Raises:
        ValueError if any angle cannot be reasonably rounded to an allowed value
    """
    if not angles_dict:
        return {}
    
    validated = {}
    for site, angles in angles_dict.items():
        validated_angles = []
        for angle in angles:
            # Find the closest allowed angle
            closest = min(ALLOWED_ANGLES, key=lambda x: abs(x - angle))
            if abs(closest - angle) > 15:  # If more than 15 degrees off, it's probably a mistake
                raise ValueError(f"Angle {angle}° for site {site} is too far from allowed values {ALLOWED_ANGLES}")
            validated_angles.append(closest)
        validated[site] = validated_angles
    return validated

def reset_to_octahedral_positions(mol):
    """Reset all donor atoms to their original OCTAHEDRAL_POS positions."""
    if mol is None:
        raise ValueError("Input molecule is None")
    
    rw = Chem.RWMol(mol)
    conf = rw.GetConformer()
    
    # Find metal center (ensure we have one)
    metal_idx = None
    for atom in rw.GetAtoms():
        if atom.GetSymbol() in METALS:
            metal_idx = atom.GetIdx()
            break
    
    if metal_idx is None:
        raise ValueError("No metal atom found in molecule")
    
    # Set metal to origin
    conf.SetAtomPosition(metal_idx, Point2D(0, 0, 0))
    
    # Find donor atoms (connected to metal)
    donors = []
    for bond in rw.GetBonds():
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a1 == metal_idx:
            donors.append(a2)
        elif a2 == metal_idx:
            donors.append(a1)
    
    # Verify we have exactly 6 donors for octahedral
    if len(donors) != 6:
        raise ValueError(f"Expected 6 donor atoms, found {len(donors)}")
    
    # Create coordinate map for metal and donors
    coord_map = {metal_idx: Point2D(0, 0)}
    for i, donor_idx in enumerate(donors):
        if i < len(OCTAHEDRAL_POS):
            pos = OCTAHEDRAL_POS[i]
            coord_map[donor_idx] = pos
            conf.SetAtomPosition(donor_idx, Point2D(pos.x, pos.y, 0))
    
    # Recompute coordinates while keeping metal and donors fixed
    AllChem.Compute2DCoords(rw, coordMap=coord_map, clearConfs=True)
    
    return rw.GetMol()

def ensure_proper_smiles_drawing(mol):
    """Ensure SMILES are correctly drawn while preserving ligand positions.
    
    Args:
        mol: RDKit molecule object of the coordination complex
        
    Returns:
        RDKit molecule with proper SMILES drawing but preserved ligand positions
    """
    # Create a copy of the molecule to work with
    mol_copy = Chem.Mol(mol)
    
    # Store the original coordinates
    original_conformer = mol_copy.GetConformer()
    original_coords = {}
    for i in range(mol_copy.GetNumAtoms()):
        pos = original_conformer.GetAtomPosition(i)
        original_coords[i] = Point2D(pos.x, pos.y)
    
    # Generate a clean 2D structure from SMILES
    AllChem.Compute2DCoords(mol_copy)
    
    # Now reapply our original coordinates for the metal and donor atoms
    # We'll identify these by their connectivity to the metal
    metal_idx = None
    for atom in mol_copy.GetAtoms():
        if atom.GetSymbol() in METALS:
            metal_idx = atom.GetIdx()
            break
    
    if metal_idx is None:
        return mol_copy  # No metal found, return the cleaned version
    
    # Find all atoms directly bonded to the metal
    donor_indices = []
    for bond in mol_copy.GetBonds():
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a1 == metal_idx:
            donor_indices.append(a2)
        elif a2 == metal_idx:
            donor_indices.append(a1)
    
    # Create a coordinate map for the metal and donors
    coord_map = {}
    for idx in [metal_idx] + donor_indices:
        if idx in original_coords:
            coord_map[idx] = original_coords[idx]
    
    # Recompute coordinates while keeping metal and donors fixed
    AllChem.Compute2DCoords(mol_copy, coordMap=coord_map, clearConfs=True)
    
    return mol_copy
#Load a ligand from SDF, but build from its SMILES if available
def load_ligand(sdf_path):
    suppl = Chem.SDMolSupplier(sdf_path)
    orig = next((m for m in suppl if m), None)
    if orig is None:
        raise FileNotFoundError(f"No valid molecule in {sdf_path}")
    # If the SDF has a SMILES property, use that to build the ligand
    if orig.HasProp("SMILES"):
        smiles = orig.GetProp("SMILES").strip()
        lig = Chem.MolFromSmiles(smiles)
        if lig is None:
            raise ValueError(f"Cannot parse SMILES '{smiles}' from {sdf_path}")
        # Preserve original properties like DENTATE, LigandName, etc.
        for prop in orig.GetPropNames():
            lig.SetProp(prop, orig.GetProp(prop))
        # Generate 2D coordinates from SMILES
        return lig
    # Fallback to original SDF-loaded molecule
    return orig

# Identify donor atom indices
# Define donor atom priority (higher priority first)
DONOR_PRIORITY = ["O", "N", "S", "P", "Cl"]

def get_donor_atom_indices(lig):
    # First check for explicit Denticity property (numeric)
    if lig.HasProp("Denticity"):
        try:
            denticity = int(lig.GetProp("Denticity"))
            if denticity <= 0:
                raise ValueError(f"Invalid Denticity value {denticity} - must be positive")
            
            # Collect all candidate donor atoms with priority
            candidate_atoms = [(atom.GetIdx(), atom.GetSymbol()) 
                             for atom in lig.GetAtoms()
                             if atom.GetSymbol() in DONOR_PRIORITY]
            
            if len(candidate_atoms) < denticity:
                raise ValueError(f"Ligand claims denticity={denticity} but only {len(candidate_atoms)} donor-atoms found in priority set {DONOR_PRIORITY}")
            
            # Sort by priority (O first, then N, etc.)
            candidate_atoms.sort(key=lambda x: DONOR_PRIORITY.index(x[1]))
            
            # Return exactly the first 'denticity' donor indices
            return [idx for idx, sym in candidate_atoms[:denticity]]
        except ValueError as e:
            raise ValueError(f"Invalid Denticity property: {e}")

    # Fallback to DENTATE property (textual or numeric)
    if lig.HasProp("DENTATE"):
        dent_str = lig.GetProp("DENTATE").strip().lower()
        dent_map = {
            "monodentate": 1,
            "bidentate":   2,
            "tridentate":  3,
            "tetradentate":4,
            "pentadentate":5,
            "hexadentate": 6
        }
        if dent_str.isdigit():
            dent_val = int(dent_str)
        else:
            try:
                dent_val = dent_map[dent_str]
            except KeyError:
                raise ValueError(f"Unrecognized DENTATE value '{dent_str}' in ligand {lig.GetProp('LigandName') if lig.HasProp('LigandName') else ''}")
        
        # Collect all candidate donor atoms with priority
        candidate_atoms = [(atom.GetIdx(), atom.GetSymbol())
                         for atom in lig.GetAtoms()
                         if atom.GetSymbol() in DONOR_PRIORITY]
        if len(candidate_atoms) < dent_val:
            raise ValueError(f"Ligand claims denticity={dent_val} but only {len(candidate_atoms)} donor-atoms found in priority set {DONOR_PRIORITY}")
        
        # Sort by priority (O first, then N, etc.)
        candidate_atoms.sort(key=lambda x: DONOR_PRIORITY.index(x[1]))
        
        return [idx for idx, sym in candidate_atoms[:dent_val]]

    # Final fallback - return all candidate donor atoms sorted by priority
    candidate_atoms = [(atom.GetIdx(), atom.GetSymbol())
                      for atom in lig.GetAtoms()
                      if atom.GetSymbol() in DONOR_PRIORITY]
    candidate_atoms.sort(key=lambda x: DONOR_PRIORITY.index(x[1]))
    return [idx for idx, sym in candidate_atoms]
# Compute explicit carbon positioning around a donor atom
def compute_carbon_position(donor_pos, angle_deg, distance=1.5):
    angle_rad = math.radians(angle_deg)
    cx = donor_pos.x + distance * math.cos(angle_rad)
    cy = donor_pos.y + distance * math.sin(angle_rad)
    return Point2D(cx, cy)

# Build coordination complex with optional custom carbon angles per octahedral site
# carbon_angles: dict mapping octahedral position index (0-5) -> list of angles in degrees


def build_coordination_complex(metal_smiles, ligand_mols, geometry='octahedral', carbon_angles=None, chain_delta_angle=60.0, chain_delta_r=0.5):
    from rdkit.Chem import rdMolTransforms

    def get_vector(p1, p2):
        return p2.x - p1.x, p2.y - p1.y

    def angle_between(v1, v2):
        x1, y1 = v1
        x2, y2 = v2
        dot = x1 * x2 + y1 * y2
        det = x1 * y2 - y1 * x2
        return math.atan2(det, dot)

    def transform_ligand_rigid(lig, donor1, donor2, target1, target2):
        lig = Chem.Mol(lig)
        conf = lig.GetConformer()
        p1 = conf.GetAtomPosition(donor1)
        p2 = conf.GetAtomPosition(donor2)
        lig_vec = get_vector(p1, p2)
        tgt_vec = get_vector(target1, target2)

        angle = angle_between(lig_vec, tgt_vec)

        cx = (p1.x + p2.x) / 2
        cy = (p1.y + p2.y) / 2

        tx = (target1.x + target2.x) / 2 - cx
        ty = (target1.y + target2.y) / 2 - cy

        cos_a = math.cos(angle)
        sin_a = math.sin(angle)

        for i in range(lig.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            x, y = pos.x - cx, pos.y - cy
            xr = x * cos_a - y * sin_a + cx + tx
            yr = x * sin_a + y * cos_a + cy + ty
            conf.SetAtomPosition(i, Point3D(xr, yr, 0.0))

        return lig

    def translate_rotate_ligand(lig, donor_idx, target_pos, direction_angle):
        lig = Chem.Mol(lig)
        conf = lig.GetConformer()

        donor_pos = conf.GetAtomPosition(donor_idx)

        # Step 1: Center donor at origin
        dx, dy = -donor_pos.x, -donor_pos.y
        for i in range(lig.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, Point3D(pos.x + dx, pos.y + dy, 0.0))

        # Step 2: Rotate ligand
        cos_theta = math.cos(direction_angle)
        sin_theta = math.sin(direction_angle)
        for i in range(lig.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            x_rot = pos.x * cos_theta - pos.y * sin_theta
            y_rot = pos.x * sin_theta + pos.y * cos_theta
            conf.SetAtomPosition(i, Point3D(x_rot, y_rot, 0.0))

        # Step 3: Translate to target
        for i in range(lig.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, Point3D(pos.x + target_pos.x, pos.y + target_pos.y, 0.0))
        
        return lig
    def place_ligand_conditional_flip(ligand, donor_idx, metal_pos):
        """Flip ligand ONLY if its first heavy neighbor points toward the metal."""
        # 1) Get donor and metal positions
        conf = ligand.GetConformer()
        donor_pos = conf.GetAtomPosition(donor_idx)
        dx, dy = metal_pos.x - donor_pos.x, metal_pos.y - donor_pos.y
        base_angle = math.atan2(dy, dx)  # Angle from donor to metal

        # 2) Place ligand initially
        lig_placed = translate_rotate_ligand(ligand, donor_idx, metal_pos, base_angle)
        conf = lig_placed.GetConformer()

        # 3) Find the first heavy neighbor (not H, not metal)
        donor_atom = lig_placed.GetAtomWithIdx(donor_idx)
        neighbor_idx = None
        for nbr in donor_atom.GetNeighbors():
            if nbr.GetAtomicNum() > 1:  # Exclude hydrogens
                neighbor_idx = nbr.GetIdx()
                break

        if neighbor_idx is None:
            return lig_placed  # No heavy neighbor (e.g., H-only ligand), skip flip check

        # 4) Compare vectors:
        #    v1 = donor -> neighbor (ligand direction)
        #    v2 = donor -> metal
        neighbor_pos = conf.GetAtomPosition(neighbor_idx)
        v_lig = (neighbor_pos.x - donor_pos.x, neighbor_pos.y - donor_pos.y)
        v_met = (metal_pos.x - donor_pos.x, metal_pos.y - donor_pos.y)

        # 5) Flip ONLY if the angle between v_lig and v_met is < 90°
        dot = v_lig[0] * v_met[0] + v_lig[1] * v_met[1]  # Dot product
        if dot > 30:  # cos(θ) > 0 → θ < 90° → pointing toward metal
            lig_placed = translate_rotate_ligand(
                ligand, donor_idx, metal_pos, base_angle + math.pi
            )

        return lig_placed
    

    if geometry not in GEOMETRY_POSITIONS:
        raise ValueError(f"Unsupported geometry '{geometry}'. Choose from {list(GEOMETRY_POSITIONS)}.")

    positions = GEOMETRY_POSITIONS[geometry]
    carbon_angles = validate_carbon_angles(carbon_angles) if carbon_angles else {}

    # Create metal center
    m = Chem.MolFromSmiles(metal_smiles)
    if not m or m.GetNumAtoms() != 1:
        raise ValueError(f"Invalid metal SMILES: {metal_smiles}")

    metal_symbol = m.GetAtomWithIdx(0).GetSymbol()
    metal = Chem.RWMol()
    metal_idx = metal.AddAtom(Chem.Atom(metal_symbol))

    bond_indices = []
    full_coord_map = {metal_idx: Point2D(0.0, 0.0)}
    donor_site_map = {}

    PAIR_POSITIONS_BY_GEOMETRY = {
        'octahedral': [(0, 2), (4, 1), (5, 3)],
        'tetrahedral': [(0, 1), (2, 3)],
        'square_planar': [(0, 1), (2, 3)],
}
    pair_positions = PAIR_POSITIONS_BY_GEOMETRY[geometry]
    used_positions = set()
    bidentate_count = 0
    monodentate_ligs = []
    monodentate_indices = []

    # Handle bidentate ligands
    for i, lig in enumerate(ligand_mols):
        AllChem.Compute2DCoords(lig)
        donors = get_donor_atom_indices(lig)
        if len(donors) > 1:
            if bidentate_count >= len(pair_positions):
                raise ValueError("Too many bidentate ligands for this geometry")
            p0, p1 = pair_positions[bidentate_count]
            lig = transform_ligand_rigid(lig, donors[0], donors[1], positions[p0], positions[p1])
            #lig = place_ligand_conditional_flip(lig, donors[0], positions[p0])
            site_assignments = [p0, p1]
            used_positions.update(site_assignments)
            bidentate_count += 1

            combo = Chem.CombineMols(metal.GetMol(), lig)
            rw = Chem.RWMol(combo)
            offset = metal.GetNumAtoms()
            conf = lig.GetConformer()

            for j, d in enumerate(donors[:2]):
                di = d + offset
                rw.AddBond(metal_idx, di, Chem.BondType.SINGLE)
                bond = rw.GetMol().GetBondBetweenAtoms(metal_idx, di)
                bond_indices.append(bond.GetIdx())
                donor_site_map[di] = site_assignments[j]

            for k in range(lig.GetNumAtoms()):
                p = conf.GetAtomPosition(k)
                full_coord_map[offset + k] = Point2D(p.x, p.y)

            metal = rw
            ligand_mols[i] = lig
        else:
            monodentate_ligs.append((lig, donors))
            monodentate_indices.append(i)

    # Handle monodentate ligands
    available_sites = [i for i in range(len(positions)) if i not in used_positions]

    for i, (lig, donors) in enumerate(monodentate_ligs):
        site = available_sites[i]
        used_positions.add(site)
        angle = math.atan2(positions[site].y, positions[site].x)
        lig = translate_rotate_ligand(lig, donors[0], positions[site], angle)
        lig = place_ligand_conditional_flip(lig, donors[0], positions[site])
        combo = Chem.CombineMols(metal.GetMol(), lig)
        rw = Chem.RWMol(combo)
        offset = metal.GetNumAtoms()
        conf = lig.GetConformer()

        di = donors[0] + offset
        rw.AddBond(metal_idx, di, Chem.BondType.SINGLE)
        bond = rw.GetMol().GetBondBetweenAtoms(metal_idx, di)
        bond_indices.append(bond.GetIdx())
        donor_site_map[di] = site

        for k in range(lig.GetNumAtoms()):
            p = conf.GetAtomPosition(k)
            full_coord_map[offset + k] = Point2D(p.x, p.y)

        metal = rw
        ligand_mols[monodentate_indices[i]] = lig

    # Final coordination map and 2D drawing
    mol = metal.GetMol()
    AllChem.Compute2DCoords(mol, coordMap=full_coord_map, clearConfs=True)

    # Optional: Add wedge/dash for octahedral visualization
    for bidx in bond_indices:
        bond = mol.GetBondWithIdx(bidx)
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        donor = a1 if a1 != metal_idx else a2
        site_idx = donor_site_map.get(donor)
        if geometry == 'octahedral':
            if site_idx in (4, 5):
                bond.SetBondDir(rdchem.BondDir.BEGINWEDGE)
            elif site_idx in (2, 3):
                bond.SetBondDir(rdchem.BondDir.BEGINDASH)
    
    return mol



# Draw molecule to 2D image
def draw_mol_2D(mol, img_size=(500,500), output_path=None):
    drawer = rdMolDraw2D.MolDraw2DCairo(*img_size)
    opts = drawer.drawOptions()
    opts.includeAtomTags = False
    opts.wedgeStrokeWidth = 10
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    img = Image.open(io.BytesIO(drawer.GetDrawingText()))
    if output_path:
        img.save(output_path)
    return img

# Load all ligands from folder
LIGAND_FOLDER = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "ligands"))
def load_ligands_from_folder(path):
    lig_map = {}
    for fn in os.listdir(path):
        if fn.endswith('.sdf'):
            lig_map[os.path.splitext(fn)[0]] = load_ligand(os.path.join(path, fn))
    return lig_map

# Create and draw complex from ligand count dict with optional carbon angles per site
def create_complex_from_ligand_dict(metal, ligand_counts, geometry='octahedral', carbon_angles=None, chain_delta_angle: float = 60.0, chain_delta_r: float = 0.5, output_file=None):
    """
    Create and draw a coordination complex from a ligand dictionary.
    """
    lig_map = load_ligands_from_folder(LIGAND_FOLDER)
    total_sites = sum(len(get_donor_atom_indices(lig_map[n])) * c for n, c in ligand_counts.items())
    expected_sites = len(GEOMETRY_POSITIONS[geometry])
    if total_sites != expected_sites:
        raise ValueError(f"Total donor sites must equal {expected_sites} for {geometry} geometry")

    ligs = []
    for n, c in ligand_counts.items():
        ligs.extend([lig_map[n]] * c)

    mol = build_coordination_complex(
        metal_smiles=f"[{metal}]",  # Oxidation state is now included in the metal string
        ligand_mols=ligs,
        geometry=geometry,
        carbon_angles=carbon_angles,
        chain_delta_angle=chain_delta_angle,
        chain_delta_r=chain_delta_r
    )

    img = draw_mol_2D(mol, output_path=output_file)

    if output_file and os.path.exists(output_file):
        im = Image.open(output_file)
        im.show()

    return img


# --- Angle selection ---
# Define your custom angles here, keyed by site index (0-5). Remove entries to use DEFAULT_CARBON_ANGLES.
CARBON_ANGLE_SELECTION = {
     0: [30, 150],  # top
     1: [210, 330],  # bottom
     2: [90, 210],  # back-left
     3: [330, 90],  # back-right
     4: [150, 270],  # front-left
     5: [270, 30],  # front-right
}


from difflib import get_close_matches

def find_closest_ligand(name, lig_map):
    matches = get_close_matches(name, lig_map.keys(), n=5, cutoff=0.5)
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

if __name__ == "__main__":
    print("🎨 Coordination Complex Builder")
    lig_map = load_ligands_from_folder(LIGAND_FOLDER)
    print(f"✅ Found {len(lig_map)} ligands:", ", ".join(sorted(lig_map.keys())))

    metal = input("Enter metal symbol (e.g., Fe): ").strip().capitalize()
    
    # Find metal in METALS list with error handling
    metal_info = None
    try:
        for m in METALS:
            # Handle case where m might not be a dictionary or lacks 'symbol' key
            if isinstance(m, dict) and 'symbol' in m:
                if m['symbol'].upper() == metal.upper():
                    metal_info = m
                    break
    except Exception as e:
        print(f"⚠️ Error accessing metals database: {e}")

    if not metal_info:
        print(f"⚠️ Warning: Metal {metal} not found in database. Using default +2")
        oxidation_state = 2
    else:
        # Safely get oxidation states with default [2] if not specified
        available_states = metal_info.get('oxidation_states', [2])
        if not isinstance(available_states, list):  # Handle case where oxidation_states isn't a list
            available_states = [2]
            print(f"⚠️ Warning: Invalid oxidation states format for {metal}. Using default [2]")
        
        print(f"Available oxidation states for {metal}: {', '.join(map(str, available_states))}")
        while True:
            try:
                oxidation_state = int(input(f"Enter oxidation state for {metal}: ").strip())
                if oxidation_state in available_states:
                    break
                print(f"❌ {oxidation_state} is not a common oxidation state for {metal}. Try again.")
            except ValueError:
                print("❌ Please enter a valid integer.")
                
    print("\nChoose geometry:")
    print("1) Octahedral")
    print("2) Tetrahedral")
    print("3) Square planar")
    geometry_input = input("Enter number [1-3]: ").strip()
    geometry_map = {"1": "octahedral", "2": "tetrahedral", "3": "square_planar"}
    geometry = geometry_map.get(geometry_input)
    if not geometry:
        print("❌ Invalid selection.")
        exit()

    ligand_counts = {}
    print("➕ Enter ligand names and counts. Type 'done' to finish.")
    while True:
        name = input("Ligand name: ").strip()
        if name.lower() == 'done':
            break
        if name not in lig_map:
            suggested = find_closest_ligand(name, lig_map)
            if not suggested:
                print("❌ Ligand not found.")
                continue
            name = suggested
        try:
            count = int(input(f"How many of '{name}'? "))
            ligand_counts[name] = count
        except ValueError:
            print("❌ Invalid number.")

    if not ligand_counts:
        print("⚠️ No ligands provided. Exiting.")
        exit()

    total_sites = sum(len(get_donor_atom_indices(lig_map[n])) * c for n, c in ligand_counts.items())
    expected_sites = len(GEOMETRY_POSITIONS[geometry])
    print(f"\n🧮 Total donor sites: {total_sites}")
    if total_sites != expected_sites:
        print(f"⚠️ For {geometry} complexes, you need exactly {expected_sites} donor sites.")
        confirm = input("Continue anyway? (y/N): ").strip().lower()
        if confirm != 'y':
            print("❌ Aborted.")
            exit()
    
    out_file = f"{metal}({oxidation_state})_{geometry}_{'_'.join(f'{k}{v}' for k, v in ligand_counts.items())}.png"
    try:
        # Include oxidation state in the metal string
        metal_with_oxidation = f"{metal}+{oxidation_state}" if oxidation_state != 0 else metal
        img = create_complex_from_ligand_dict(
            metal=metal_with_oxidation,  # Now includes oxidation state
            ligand_counts=ligand_counts,
            geometry=geometry,
            carbon_angles=CARBON_ANGLE_SELECTION,
            output_file=out_file
        )
        print(f"✅ Complex image saved to {out_file}")
        img.show()
    except Exception as e:
        print(f"❌ Failed to create complex: {e}")
