import os
import io
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point2D, Point3D
from PIL import Image

# Ligands as heavy atoms only
LIGAND_SYMBOLS = {
    "NH3": "N",
    "Cl": "Cl",
    "Br": "Br",
    "H2O": "O",
    "CN": "C"
}

# Manual 2D layout coordinates (octahedral)
OCTAHEDRAL_POS = [
    Point2D(0.0, 1.5),    # top
    Point2D(0.0, -1.5),   # bottom
    Point2D(-1.5, 0.0),   # left
    Point2D(1.5, 0.0),    # right
    Point2D(-1.0, -1.0),  # bottom-left
    Point2D(1.0, -1.0)    # bottom-right
]

def draw_coordination_complex_2D(metal="Co", ligands=None, geometry="octahedral", output_file=None):
    if ligands is None:
        raise ValueError("Please provide a list of ligands.")
    if geometry != "octahedral":
        raise NotImplementedError("Only octahedral layout is implemented.")
    if len(ligands) != 6:
        raise ValueError("Octahedral geometry requires exactly 6 ligands.")

    mol = Chem.RWMol()
    pos_dict = {}

    # Add metal atom at center
    metal_idx = mol.AddAtom(Chem.Atom(metal))
    pos_dict[metal_idx] = Point2D(0.0, 0.0)

    # Add ligands as heavy atoms only
    for i, lig in enumerate(ligands):
        symbol = LIGAND_SYMBOLS.get(lig, lig)
        lig_idx = mol.AddAtom(Chem.Atom(symbol))
        pos_dict[lig_idx] = OCTAHEDRAL_POS[i]
        mol.AddBond(metal_idx, lig_idx, Chem.BondType.SINGLE)

    # Create 2D-style conformer using Point3D
    conf = Chem.Conformer(mol.GetNumAtoms())
    for idx, pos in pos_dict.items():
        conf.SetAtomPosition(idx, Point3D(pos.x, pos.y, 0.0))

    mol = mol.GetMol()
    mol.RemoveAllConformers()
    mol.AddConformer(conf)

    # Remove hydrogens before drawing
    mol_no_H = Chem.RemoveHs(mol)

    # Draw image
    drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)
    opts = drawer.drawOptions()
    opts.addAtomIndices = False
    opts.addBondIndices = False
    opts.includeAtomTags = False
    opts.explicitMethyl = True
    opts.addStereoAnnotation = False

    drawer.DrawMolecule(mol_no_H)
    drawer.FinishDrawing()

    img_data = drawer.GetDrawingText()
    img = Image.open(io.BytesIO(img_data))

    if output_file:
        img.save(output_file)
        print(f"✅ Saved: {output_file}")
    img.show()

    return mol

# === Example Usage ===
if __name__ == "__main__":
    draw_coordination_complex_2D(
        metal="Co",
        ligands=["NH3", "NH3", "Cl", "Cl", "NH3", "NH3"],
        geometry="octahedral",
        output_file="Co_complex_clean.png"
    )
