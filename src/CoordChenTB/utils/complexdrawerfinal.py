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
from PIL import ImageDraw
from PIL import ImageFont

# Path to metals database
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from CoordChenTB.utils.metals_db import METALS

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
    'square_planar': [
        Point2D(1.5, 1.5),
        Point2D(-1.5, 1.5),
        Point2D(-1.5, -1.5),
        Point2D(1.5, -1.5),
    ],
    'tetrahedral': [
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

bond_length = 1 #To later on change bond length

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

def predict_4coordinate_geometry(metal_smiles, ligand_smiles_list):
    """
    Predicts tetrahedral vs. square planar geometry for a 4-coordinate complex.
    Simplified rules:
    - d⁸ metals (Ni²⁺, Pd²⁺, Pt²⁺) → square planar
    - 4d/5d metals → square planar
    - Everything else → tetrahedral
    """
    # Parse metal
    metal_mol = Chem.MolFromSmiles(metal_smiles)
    if not metal_mol or metal_mol.GetNumAtoms() != 1:
        raise ValueError("Invalid metal SMILES")
    
    metal_symbol = metal_mol.GetAtomWithIdx(0).GetSymbol()
    metal_info = next((m for m in METALS if m.get("symbol") == metal_symbol), None)
    if not metal_info:
        raise ValueError(f"Metal {metal_symbol} not found in database")
    
    # Get metal properties
    oxidation_state = metal_info.get("common_oxidation_state", 2)
    block = metal_info.get("block", "")
    atomic_number = metal_info.get("atomic_number", 0)
    
    # Calculate d-electrons
    if block == '3d':
        d_electrons = (atomic_number - 21) - oxidation_state  # 3d starts at Sc (21)
    elif block in ('4d', '5d'):
        d_electrons = (atomic_number - 39) - oxidation_state  # 4d starts at Y (39)
    else:
        d_electrons = 0
    
    # Decision logic
    if d_electrons == 8:  # Ni²⁺, Pd²⁺, Pt²⁺
        return "square_planar"
    elif block in ('4d', '5d'):  # All other 4d/5d metals
        return "square_planar"
    else:  # Everything else
        return "tetrahedral"
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

# Build coordination complex with optional custom carbon angles per octahedral site
# carbon_angles: dict mapping octahedral position index (0-5) -> list of angles in degrees


def build_coordination_complex(metal_smiles, ligand_mols, geometry='octahedral', carbon_angles=None, chain_delta_angle=60.0, chain_delta_r=0.5, bond_length=1.0):
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


    if geometry not in GEOMETRY_POSITIONS:
        raise ValueError(f"Unsupported geometry '{geometry}'. Choose from {list(GEOMETRY_POSITIONS)}.")

    positions = GEOMETRY_POSITIONS[geometry]
    positions = [Point2D(p.x * bond_length, p.y * bond_length) for p in positions]
    carbon_angles = validate_carbon_angles(carbon_angles) if carbon_angles else {}
    total_donors = sum(len(get_donor_atom_indices(lig)) for lig in ligand_mols)

    m = Chem.MolFromSmiles(metal_smiles)
    if not m or m.GetNumAtoms() != 1:
        raise ValueError(f"Invalid metal SMILES: {metal_smiles}")


    if geometry is None:
        if total_donors == 4:
            ligand_smiles_list = [Chem.MolToSmiles(lig) for lig in ligand_mols]
            geometry = predict_4coordinate_geometry(metal_smiles, ligand_smiles_list)
        elif total_donors == 6:
            geometry = "octahedral"
        else:
            raise ValueError(f"Unsupported number of donor atoms: {total_donors}")
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

    total_donors = sum(len(get_donor_atom_indices(lig)) for lig in ligand_mols)

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
        elif geometry == 'tetrahedral':
            if site_idx == 0:  # position 0 gets dash
                bond.SetBondDir(rdchem.BondDir.BEGINDASH)
            elif site_idx == 3:  # position 3 gets wedge
                bond.SetBondDir(rdchem.BondDir.BEGINWEDGE)
    return mol



# Draw molecule to 2D image
def draw_mol_2D(mol, img_size=(2000,2000), output_path=None, ligand_mols=None, metal_smiles=None):
    drawer = rdMolDraw2D.MolDraw2DCairo(*img_size)
    opts = drawer.drawOptions()
    opts.includeAtomTags = False
    opts.wedgeStrokeWidth = 10
    opts.getBackgroundColour = lambda: (1, 1, 1)  # White background

    opts.bondLineWidth = 2
    opts.scaleBondWidth = False
    opts.scale = 0.8  # Reduced from default 1.0 to zoom out

    total_charge = 0
    if metal_smiles:
        # Parse oxidation state from metal SMILES
        if '+' in metal_smiles:
            charge_part = metal_smiles.split('+')[1].rstrip(']')
            if charge_part.isdigit():
                total_charge = int(charge_part)
            else:
                total_charge = len(charge_part)  # For cases like [Fe++]
        
    if ligand_mols:
        for lig in ligand_mols:
            if lig.HasProp("NetFormalCharge"):
                try:
                    total_charge += int(lig.GetProp("NetFormalCharge"))
                except ValueError:
                    pass

    # Draw the molecule first to get coordinates
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    
    # Get the image and convert to PIL
    img = Image.open(io.BytesIO(drawer.GetDrawingText()))
    draw = ImageDraw.Draw(img)
    
    # Get all atomic positions in image coordinates
    atom_positions = []
    for atom in mol.GetAtoms():
        point = drawer.GetDrawCoords(atom.GetIdx())
        atom_positions.append((point.x, point.y))
    
    if not atom_positions:
        if output_path:
            img.save(output_path)
        return img
    
    # Get bounding box in image coordinates
    min_x = min(p[0] for p in atom_positions)
    max_x = max(p[0] for p in atom_positions)
    min_y = min(p[1] for p in atom_positions)
    max_y = max(p[1] for p in atom_positions)
    
    # Increased bracket parameters for more space
    bracket_padding = 55  # Increased from 15 (space between molecule and brackets)
    bracket_thickness = 3  # Line width
    bracket_extension = 20  # Increased from 20 (length of horizontal lines)
    
    # Left bracket ([)
    left = min_x - bracket_padding
    draw.line([(left, min_y - bracket_padding), 
              (left, max_y + bracket_padding)], 
              fill="black", width=bracket_thickness)
    draw.line([(left, min_y - bracket_padding), 
              (left + bracket_extension, min_y - bracket_padding)], 
              fill="black", width=bracket_thickness)
    draw.line([(left, max_y + bracket_padding), 
              (left + bracket_extension, max_y + bracket_padding)], 
              fill="black", width=bracket_thickness)
    
    # Right bracket (])
    right = max_x + bracket_padding
    draw.line([(right, min_y - bracket_padding), 
              (right, max_y + bracket_padding)], 
              fill="black", width=bracket_thickness)
    draw.line([(right - bracket_extension, min_y - bracket_padding), 
              (right, min_y - bracket_padding)], 
              fill="black", width=bracket_thickness)
    draw.line([(right - bracket_extension, max_y + bracket_padding), 
              (right, max_y + bracket_padding)], 
              fill="black", width=bracket_thickness)
    
    # Draw charge notation (only if non-zero)
    charge_text = None
    if total_charge != 0:
        charge_text = f"{abs(total_charge)}{'+' if total_charge > 0 else '-'}"
        if abs(total_charge) == 1:
            charge_text = charge_text[0]  # Remove the 1 for single charges

    if charge_text:
        # Position at top-right of right bracket
        text_x = right + 20
        text_y = min_y - bracket_padding - 30

        # Calculate font size based on image dimensions
        font_size = int(min(img_size) * 0.08)  # 8% of smaller dimension

        try:
            font = ImageFont.truetype("arial.ttf", font_size)
        except:
            font = ImageFont.load_default()

        # Adjust position to avoid clipping
        temp_img = Image.new('RGB', (1, 1))
        temp_draw = ImageDraw.Draw(temp_img)
        text_bbox = temp_draw.textbbox((0, 0), charge_text, font=font)
        text_width = text_bbox[2] - text_bbox[0]
        text_height = text_bbox[3] - text_bbox[1]

        if text_x + text_width > img_size[0]:
            text_x = img_size[0] - text_width - 30
        if text_y + text_height > img_size[1]:
            text_y = img_size[1] - text_height - 30

        draw.text(
            (text_x, text_y),
            charge_text,
            fill="black",
            font=font,
            stroke_width=max(2, int(font_size * 0.1)),
            stroke_fill="white"
        )

    # Final cropping
    margin = 50
    bbox = (
        max(0, min_x - bracket_padding - bracket_extension - margin),
        max(0, min_y - bracket_padding - margin),
        min(img_size[0], max_x + bracket_padding + bracket_extension + margin + 50),  # Extra space for charge
        min(img_size[1], max_y + bracket_padding + margin)
    )


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
def create_complex_from_ligand_dict(metal, ligand_counts, geometry='octahedral', carbon_angles=None, 
                                  chain_delta_angle: float = 60.0, chain_delta_r: float = 0.5, 
                                  output_file=None, oxidation_state=0, bond_length=1.0):
    lig_map = load_ligands_from_folder(LIGAND_FOLDER)
    total_sites = sum(len(get_donor_atom_indices(lig_map[n])) * c for n, c in ligand_counts.items())
    expected_sites = len(GEOMETRY_POSITIONS[geometry])
    
    if total_sites != expected_sites:
        raise ValueError(f"Total donor sites must equal {expected_sites} for {geometry} geometry")

    ligs = []
    for n, c in ligand_counts.items():
        ligs.extend([lig_map[n]] * c)

    mol = build_coordination_complex(
        metal_smiles=f"[{metal}]",
        ligand_mols=ligs,
        geometry=geometry,
        carbon_angles=carbon_angles,
        chain_delta_angle=chain_delta_angle,
        chain_delta_r=chain_delta_r,
        bond_length=bond_length 
)


    img = draw_mol_2D(
        mol, 
        output_path=output_file,
        ligand_mols=ligs,
        metal_smiles=f"[{metal}+{oxidation_state}]" if oxidation_state != 0 else f"[{metal}]"
    )

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

# Replace the geometry input section in the __main__ block with this:
if __name__ == "__main__":
    print("🎨 Coordination Complex Builder")
    lig_map = load_ligands_from_folder(LIGAND_FOLDER)
    
    metal = input("Enter metal symbol (e.g., Fe): ").strip().capitalize()
    
    # Get metal info (existing code)
    metal_info = None
    try:
        for m in METALS:
            if isinstance(m, dict) and 'symbol' in m and m['symbol'].upper() == metal.upper():
                metal_info = m
                break
    except Exception as e:
        print(f"⚠️ Error accessing metals database: {e}")

    if not metal_info:
        print(f"⚠️ Warning: Metal {metal} not found in database. Using default +2")
        oxidation_state = 2
    else:
        available_states = metal_info.get('oxidation_states', [2])
        if not isinstance(available_states, list):
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

    # --- MODIFIED GEOMETRY SELECTION LOGIC ---
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

    # Calculate total donor sites
    total_sites = sum(len(get_donor_atom_indices(lig_map[n])) * c for n, c in ligand_counts.items())
    
    # Auto-determine geometry for common cases
    geometry = None
    if total_sites == 4:
        metal_with_ox = f"{metal}+{oxidation_state}" if oxidation_state != 0 else metal
        ligand_smiles_list = []
        for name, count in ligand_counts.items():
            lig = lig_map[name]
            ligand_smiles_list.extend([Chem.MolToSmiles(lig)] * count)
        
        try:
            geometry = predict_4coordinate_geometry(f"[{metal_with_ox}]", ligand_smiles_list)
            print(f"🔍 Auto-detected geometry: {geometry} (based on metal/ligands)")
        except Exception as e:
            print(f"⚠️ Could not auto-determine geometry: {e}. Defaulting to tetrahedral.")
            geometry = "tetrahedral"
    elif total_sites == 6:
        geometry = "octahedral"
    else:
        print(f"❌ Unsupported number of donor sites: {total_sites}. Only 4 or 6 supported.")
        exit()


   
    out_file = f"{metal}({oxidation_state})_{geometry}_{'_'.join(f'{k}{v}' for k, v in ligand_counts.items())}.png"
    try:
        metal_with_oxidation = f"{metal}+{oxidation_state}" if oxidation_state != 0 else metal
        try:
            bond_length = float(input("Enter metal-ligand bond length (e.g., 1.0): ").strip())
        except ValueError:
            bond_length = 1.0
        img = create_complex_from_ligand_dict(
            metal=metal_with_oxidation,
            ligand_counts=ligand_counts,
            geometry=geometry,
            carbon_angles=CARBON_ANGLE_SELECTION,
            output_file=out_file,
            bond_length=bond_length 
        )

        print(f"✅ Complex image saved to {out_file}")
        img.show()
    except Exception as e:
        print(f"❌ Failed to create complex: {e}")