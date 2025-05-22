import os
import importlib.machinery
import math
import io
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdchem
from rdkit.Geometry import Point2D, Point3D
from PIL import Image

# Load metals database
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from CoordChenTB.utils.metals_db import METALS

# Octahedral scaffold positions
OCTAHEDRAL_POS = [
    Point2D(0.0, 3.0),    # top (0)
    Point2D(0.0, -3.0),   # bottom (1)
    Point2D(-2.6, 1.5),   # front-left (2)
    Point2D(2.6, 1.5),    # front-right (3)
    Point2D(-2.6, -1.5),  # back-left (4)
    Point2D(2.6, -1.5),   # back-right (5)
]

# Donor elements
DONOR_ELEMENTS = {"N", "O", "S", "P", "Cl"}

# Internal utility functions

def get_donor_atom_indices(lig):
    """Identify donor atom indices by element and DENTATE property."""
    donors = [a.GetIdx() for a in lig.GetAtoms() if a.GetSymbol() in DONOR_ELEMENTS]
    if lig.HasProp("DENTATE"):
        val = lig.GetProp("DENTATE").lower()
        dent = int(val) if val.isdigit() else {"monodentate":1, "bidentate":2, "tridentate":3}.get(val)
        if dent and len(donors) >= dent:
            return donors[:dent]
    return donors


def align_ligand_2D(lig, donor_idxs, targets):
    """Rigidly align a ligand's 2D conformer onto specified donor positions."""
    conf = lig.GetConformer()
    # Translation for monodentate
    if len(donor_idxs) == 1:
        p0 = conf.GetAtomPosition(donor_idxs[0])
        dx, dy = targets[0].x - p0.x, targets[0].y - p0.y
        for atom in lig.GetAtoms():
            p = conf.GetAtomPosition(atom.GetIdx())
            conf.SetAtomPosition(atom.GetIdx(), Point3D(p.x + dx, p.y + dy, 0.0))
        return
    # Rotation for poly-dentate
    p0, p1 = (conf.GetAtomPosition(i) for i in donor_idxs[:2])
    v0 = Point2D(p1.x - p0.x, p1.y - p0.y)
    v1 = Point2D(targets[1].x - targets[0].x, targets[1].y - targets[0].y)
    theta = math.atan2(v1.y, v1.x) - math.atan2(v0.y, v0.x)

    # Pair-specific rotations
    site_pair = tuple(sorted((OCTAHEDRAL_POS.index(targets[0]), OCTAHEDRAL_POS.index(targets[1]))))
    ROT_ANGLES = {
        (0,1): math.pi/2,  (2,3): math.pi/4,
        (4,5): -math.pi/4, (0,2): 0,
    }
    angle = ROT_ANGLES.get(site_pair, 0)
    new_coords = {}

    for atom in lig.GetAtoms():
        idx = atom.GetIdx()
        p = conf.GetAtomPosition(idx)
        x_rel, y_rel = p.x - p0.x, p.y - p0.y
        if idx not in donor_idxs and angle:
            x_rel, y_rel = (x_rel*math.cos(angle)-y_rel*math.sin(angle),
                            x_rel*math.sin(angle)+y_rel*math.cos(angle))
        xr = x_rel*math.cos(theta) - y_rel*math.sin(theta)
        yr = x_rel*math.sin(theta) + y_rel*math.cos(theta)
        new_coords[idx] = (xr + targets[0].x, yr + targets[0].y)
    for idx,(x,y) in new_coords.items():
        conf.SetAtomPosition(idx, Point3D(x, y, 0.0))


def build_coordination_complex(metal_smiles, ligand_mols, carbon_angles=None):
    """Constructs an octahedral coordination complex."""
    # Metal
    m = Chem.MolFromSmiles(metal_smiles)
    if not m or m.GetNumAtoms()!=1:
        raise ValueError(f"Invalid metal SMILES: {metal_smiles}")
    metal = m.GetAtomWithIdx(0).GetSymbol()
    rw = Chem.RWMol()
    metal_idx = rw.AddAtom(Chem.Atom(metal))
    complex_mol = rw.GetMol()
    conf = Chem.Conformer(1)
    conf.SetAtomPosition(0, Point3D(0,0,0))
    complex_mol.AddConformer(conf)

    # Assign sites
    donors_list = [get_donor_atom_indices(l) for l in ligand_mols]
    pair_sites = [(0,2),(4,1),(5,3)]
    site_map, used, bi= [], set(), 0
    for d in donors_list:
        if len(d)>1:
            s=pair_sites[bi]; bi+=1
            site_map.append(s); used.update(s)
        else:
            site_map.append(None)
    rem = [i for i in range(6) if i not in used]
    mi=0
    for i,d in enumerate(donors_list):
        if len(d)==1:
            site_map[i]=(rem[mi],); mi+=1
    
    combo, coord_map = complex_mol, {0: Point2D(0,0)}
    offset=combo.GetNumAtoms()
    donor_to_site={}
    for lig, sites in zip(ligand_mols, site_map):
        l2=Chem.Mol(lig); AllChem.Compute2DCoords(l2)
        idxs = get_donor_atom_indices(l2)
        targets=[OCTAHEDRAL_POS[i] for i in sites]
        align_ligand_2D(l2, idxs, targets)
        for d,s in zip(idxs,sites): donor_to_site[offset+d]=s
        combo = Chem.CombineMols(combo, l2)
        new_conf = Chem.Conformer(combo.GetNumAtoms())
        for i,pos in coord_map.items(): new_conf.SetAtomPosition(i, Point3D(pos.x,pos.y,0))
        c2 = l2.GetConformer()
        for a in l2.GetAtoms():
            ni=a.GetIdx()+offset; p=c2.GetAtomPosition(a.GetIdx())
            new_conf.SetAtomPosition(ni,Point3D(p.x,p.y,0)); coord_map[ni]=Point2D(p.x,p.y)
        combo.RemoveAllConformers(); combo.AddConformer(new_conf)
        rw2=Chem.RWMol(combo)
        for d in idxs: rw2.AddBond(metal_idx, d+offset, Chem.BondType.SINGLE)
        combo=rw2.GetMol(); offset=combo.GetNumAtoms()
    # Wedges/dashes
    for b in combo.GetBonds():
        a1,a2=b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if metal_idx not in (a1,a2): continue
        didx=a1 if a2==metal_idx else a2
        s=donor_to_site.get(didx)
        if s in (4,5): b.SetBondDir(rdchem.BondDir.BEGINWEDGE)
        elif s in (2,3): b.SetBondDir(rdchem.BondDir.BEGINDASH)
    return combo

# Single unified load and create functions

def load_ligands(sdf_folder):
    """Load all ligands from a folder, extracting SMILES and computing 2D coords."""
    lig_map={}
    for fn in os.listdir(sdf_folder):
        if not fn.lower().endswith('.sdf'): continue
        path=os.path.join(sdf_folder, fn)
        suppl=Chem.SDMolSupplier(path)
        mol = next((m for m in suppl if m), None)
        if mol is None:
            continue
        # prefer embedded SMILES
        smi = mol.GetProp('SMILES').strip() if mol.HasProp('SMILES') else Chem.MolToSmiles(mol)
        m2 = Chem.MolFromSmiles(smi)
        if m2 is None:
            continue
        for prop in mol.GetPropNames(): m2.SetProp(prop, mol.GetProp(prop))
        AllChem.Compute2DCoords(m2)
        name=os.path.splitext(fn)[0]
        lig_map[name]=m2
    return lig_map


def create_complex_from_ligands(metal_symbol, lig_input, sdf_folder=None, carbon_angles=None, output_file=None):
    """Create and optionally draw an octahedral complex.
    lig_input: dict{name:count} or list of RDKit Mol
    sdf_folder: required if lig_input is dict
    """
    if isinstance(lig_input, dict):
        if sdf_folder is None:
            raise ValueError("sdf_folder must be provided when passing ligand counts.")
        lig_map = load_ligands(sdf_folder)
        lig_list=[]
        for name,count in lig_input.items():
            if name not in lig_map:
                raise KeyError(f"Ligand '{name}' not found in {sdf_folder}")
            lig_list.extend([lig_map[name]]*count)
    elif isinstance(lig_input, list):
        lig_list = lig_input
    else:
        raise TypeError("lig_input must be a dict or list of RDKit Mol objects.")
    # Check donor count
    total_sites = sum(len(get_donor_atom_indices(l)) for l in lig_list)
    if total_sites != 6:
        raise ValueError("Total donor sites must equal 6 for octahedral geometry")
    # Build
    complex_mol = build_coordination_complex(f"[{metal_symbol}+2]", lig_list, carbon_angles)
    # Draw
    img = rdMolDraw2D.MolDraw2DCairo(500,500)
    opts = img.drawOptions(); opts.includeAtomTags=False; opts.wedgeStrokeWidth=10
    img.DrawMolecule(complex_mol); img.FinishDrawing()
    pil_img = Image.open(io.BytesIO(img.GetDrawingText()))
    if output_file:
        pil_img.save(output_file)
    return pil_img

CARBON_ANGLE_SELECTION = {
     0: [150, 30],  # top
     1: [210, 330],  # bottom
     2: [90, 210],  # back-left
     3: [330, 90],  # back-right
     4: [270, 150],  # front-left
     5: [60, 100],  # front-right
}

if __name__ == "__main__":
    # Use the ligands folder in the repo
    sdf_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "ligands"))

    # Load ligands
    ligands = load_ligands(sdf_folder)
    print(f"Loaded ligands: {list(ligands.keys())}")

from difflib import get_close_matches

def find_closest_ligand(name, lig_map):
    matches = get_close_matches(name, lig_map.keys(), n=5, cutoff=0.5)
    if not matches:
        return None
    print(f"\n🔍 Ligand '{name}' not found exactly. Did you mean one of the following?")
    for i, match in enumerate(matches, 1):
        print(f"{i}) {match}")
    print("0) None of these")

    try:
        choice = int(input("→ Choose a number [0 to skip]: ").strip())
        if 1 <= choice <= len(matches):
            return matches[choice - 1]
    except ValueError:
        pass

    return None


if __name__ == "__main__":
    sdf_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "ligands"))
    ligands = load_ligands(sdf_folder)
    
    if not ligands:
        print("❌ No ligands found.")
        exit()

    print("✅ Available ligands:")
    print(", ".join(ligands.keys()))

    metal = input("Enter metal symbol (e.g., Fe): ").strip()

    ligand_counts = {}
    print("\nEnter ligand names and quantities. Type 'done' to finish.")
    while True:
        name = input("Ligand name: ").strip()
        if name.lower() == 'done':
            break

        if name not in ligands:
            suggested = find_closest_ligand(name, ligands)
            if not suggested:
                print("❌ Ligand not recognized.")
                continue
            name = suggested
        else:
            suggested = name

        try:
            count = int(input(f"How many '{suggested}' ligands? "))
            ligand_counts[suggested] = count
        except ValueError:
            print("❌ Invalid number.")

    if not ligand_counts:
        print("❌ No ligands specified.")
        exit()

# Check total donor site count
total_sites = 0
for name, count in ligand_counts.items():
    mol = ligands[name]
    donor_count = len(get_donor_atom_indices(mol))
    total_sites += donor_count * count

    print(f"\n🧮 Total donor sites: {total_sites}")
    if total_sites != 6:
        print("⚠️ Warning: For an octahedral complex, the total donor sites should be exactly 6.")
        confirm = input("Proceed anyway? (y/N): ").strip().lower()
        if confirm != 'y':
            print("❌ Aborting complex creation.")
            exit()


    print(f"\n🔧 Building complex: {metal} with {ligand_counts}")

    try:
        img = create_complex_from_ligands(metal, ligand_counts, sdf_folder, carbon_angles=CARBON_ANGLE_SELECTION)
        save = input("💾 Save image to file? (y/N): ").strip().lower()
        if save == "y":
            out_path = f"{metal}_{'_'.join([f'{k}{v}' for k,v in ligand_counts.items()])}.png"
            img.save(out_path)
            print(f"✅ Image saved as {out_path}")
        else:
            img.show()
    except Exception as e:
        print(f"❌ Error creating complex: {e}")

