import os
import io
import math
import importlib.machinery
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdchem
from rdkit.Geometry import Point2D
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
        AllChem.Compute2DCoords(lig)
        return lig
    # Fallback to original SDF-loaded molecule
    return orig

# Identify donor atom indices
def get_donor_atom_indices(lig):
    donor_elems = {"N", "O", "S", "P", "Cl"}
    # first, collect all atoms whose element is in donor_elems
    candidate_idxs = [atom.GetIdx()
                      for atom in lig.GetAtoms()
                      if atom.GetSymbol() in donor_elems]

    # if the ligand has an explicit DENTATE property, use it
    if lig.HasProp("DENTATE"):
        dent_str = lig.GetProp("DENTATE").strip().lower()
        # map textual denticity to integer
        dent_map = {
            "monodentate": 1,
            "bidentate":   2,
            "tridentate":  3,
            "tetradentate":4,
            "pentadentate":5,
            "hexadentate": 6
        }
        # allow numeric strings too, e.g. "2"
        if dent_str.isdigit():
            dent_val = int(dent_str)
        else:
            try:
                dent_val = dent_map[dent_str]
            except KeyError:
                raise ValueError(f"Unrecognized DENTATE value '{dent_str}' in ligand {lig.GetProp('LigandName') if lig.HasProp('LigandName') else ''}")
        if len(candidate_idxs) < dent_val:
            raise ValueError(f"Ligand claims denticity={dent_val} but only {len(candidate_idxs)} donor‐atoms found in element set {donor_elems}")
        # return exactly the first dent_val donor indices
        return candidate_idxs[:dent_val]

    # otherwise fall back to the old behaviour
    return candidate_idxs

# Compute explicit carbon positioning around a donor atom
def compute_carbon_position(donor_pos, angle_deg, distance=1.5):
    angle_rad = math.radians(angle_deg)
    cx = donor_pos.x + distance * math.cos(angle_rad)
    cy = donor_pos.y + distance * math.sin(angle_rad)
    return Point2D(cx, cy)

# Build coordination complex with optional custom carbon angles per octahedral site
# carbon_angles: dict mapping octahedral position index (0-5) -> list of angles in degrees
def build_coordination_complex(metal_smiles, ligand_mols, carbon_angles=None):
    carbon_angles = validate_carbon_angles(carbon_angles) if carbon_angles else {}

    # Initialize metal center
    rw = Chem.RWMol()
    m = Chem.MolFromSmiles(metal_smiles)
    if not m or m.GetNumAtoms() != 1:
        raise ValueError(f"Invalid metal SMILES: {metal_smiles}")
    metal_symbol = m.GetAtomWithIdx(0).GetSymbol()
    metal_idx = rw.AddAtom(Chem.Atom(metal_symbol))
    mol = rw.GetMol()

    donors_per_lig = []
    bond_indices = []
    # Attach ligands and record bonds
    for lig in ligand_mols:
        combo = Chem.CombineMols(mol, lig)
        rw = Chem.RWMol(combo)
        offset = mol.GetNumAtoms()
        donors = get_donor_atom_indices(lig)
        if not donors:
            raise ValueError(f"No donor atoms in ligand {Chem.MolToSmiles(lig)}")
        current = []
        for d in donors:
            di = d + offset
            rw.AddBond(metal_idx, di, Chem.BondType.SINGLE)
            tmp = rw.GetMol()
            bidx = tmp.GetBondBetweenAtoms(metal_idx, di).GetIdx()
            bond_indices.append(bidx)
            current.append(di)
        donors_per_lig.append(current)
        mol = rw.GetMol()

    # Place donor atoms at octahedral positions
    coord_map = {metal_idx: Point2D(0.0, 0.0)}
    used_positions = set()
    pair_positions = [(0,2), (4,1), (5,3)]
    bidentate_count = 0
    donor_site_map = {}  # maps donor atom idx -> octahedral site index

    # bidentate ligands
    for donors in donors_per_lig:
        if len(donors) > 1:
            p = pair_positions[bidentate_count]
            coord_map[donors[0]] = OCTAHEDRAL_POS[p[0]]
            coord_map[donors[1]] = OCTAHEDRAL_POS[p[1]]
            donor_site_map[donors[0]] = p[0]
            donor_site_map[donors[1]] = p[1]
            used_positions.update(p)
            bidentate_count += 1

    # monodentate ligands
    rem = [i for i in range(6) if i not in used_positions]
    mi = 0
    for donors in donors_per_lig:
        if len(donors) == 1:
            site = rem[mi]
            coord_map[donors[0]] = OCTAHEDRAL_POS[site]
            donor_site_map[donors[0]] = site
            mi += 1

    AllChem.Compute2DCoords(mol, coordMap=coord_map)

    # Explicit carbon placement with custom angles per site
    donor_elems = {"N", "O", "S", "P", "Cl"}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in donor_elems:
            didx = atom.GetIdx()
            dpos = coord_map.get(didx)
            if dpos is None:
                continue
            site_idx = donor_site_map.get(didx)
            # determine angles for this site
            angles = carbon_angles.get(site_idx, DEFAULT_CARBON_ANGLES)
            c_neigh = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetSymbol() == 'C']
            for c_idx, ang in zip(c_neigh, angles):
                coord_map[c_idx] = compute_carbon_position(dpos, ang)

    AllChem.Compute2DCoords(mol, coordMap=coord_map, clearConfs=True)

    # Set wedge/dash styling
    for bidx in bond_indices:
        bond = mol.GetBondWithIdx(bidx)
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        donor = a1 if a1 != metal_idx else a2
        site_idx = donor_site_map.get(donor)
        if site_idx in (4,5):
            bond.SetBondDir(rdchem.BondDir.BEGINWEDGE)
        elif site_idx in (2,3):
            bond.SetBondDir(rdchem.BondDir.BEGINDASH)
   
    mol = ensure_proper_smiles_drawing(mol)
    #mol = reset_to_octahedral_positions(mol)
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
def create_complex_from_ligand_dict(metal_symbol, ligand_counts, carbon_angles=None, output_file=None):
    """
    carbon_angles: dict mapping octahedral site index (0-5) -> [angle1, angle2]
    """
    lig_map = load_ligands_from_folder(LIGAND_FOLDER)
    total_sites = sum(len(get_donor_atom_indices(lig_map[n]))*c for n,c in ligand_counts.items())
    if total_sites != 6:
        raise ValueError("Total donor sites must equal 6 for octahedral geometry")
    ligs = []
    for n,c in ligand_counts.items():
        ligs.extend([lig_map[n]]*c)
    mol = build_coordination_complex(f"[{metal_symbol}+2]", ligs, carbon_angles=carbon_angles)
    return draw_mol_2D(mol, output_path=output_file)

# --- Angle selection ---
# Define your custom angles here, keyed by site index (0-5). Remove entries to use DEFAULT_CARBON_ANGLES.
CARBON_ANGLE_SELECTION = {
     0: [150, 30],  # top
     1: [210, 330],  # bottom
     2: [90, 210],  # back-left
     3: [330, 90],  # back-right
     4: [270, 150],  # front-left
     5: [30, 270],  # front-right
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

if __name__ == '__main__':
    ligands = load_ligands_from_folder(LIGAND_FOLDER)

    if not ligands:
        print("❌ No ligands found in folder.")
        exit()

    print("✅ Available ligands:")
    print(", ".join(sorted(ligands)))

    metal = input("Enter metal symbol (e.g., Fe): ").strip()
    ligand_counts = {}

    print("\nEnter ligands and counts (type 'done' to finish):")
    while True:
        name = input("Ligand name: ").strip()
        if name.lower() == "done":
            break
        if name not in ligands:
            suggested = find_closest_ligand(name, ligands)
            if not suggested:
                print("❌ Ligand not recognized.")
                continue
            name = suggested
        try:
            count = int(input(f"How many '{name}' ligands? "))
            ligand_counts[name] = count
        except ValueError:
            print("❌ Invalid number.")

    if not ligand_counts:
        print("❌ No ligands selected.")
        exit()

    total_sites = sum(len(get_donor_atom_indices(ligands[n])) * c for n, c in ligand_counts.items())
    print(f"\n🧮 Total donor sites: {total_sites}")
    if total_sites != 6:
        proceed = input("⚠️ Expected 6 donor sites. Continue anyway? (y/N): ").strip().lower()
        if proceed != "y":
            print("❌ Aborting.")
            exit()

    try:
        lig_list = [ligands[n] for n in ligand_counts for _ in range(ligand_counts[n])]
        mol = build_coordination_complex(f"[{metal}+2]", lig_list, carbon_angles=CARBON_ANGLE_SELECTION)
        img = draw_mol_2D(mol)
        if input("💾 Save image? (y/N): ").strip().lower() == 'y':
            out = f"{metal}_{'_'.join(f'{k}{v}' for k,v in ligand_counts.items())}.png"
            img.save(out)
            print(f"✅ Image saved as {out}")
        else:
            img.show()
    except Exception as e:
        print(f"❌ Error creating complex: {e}")