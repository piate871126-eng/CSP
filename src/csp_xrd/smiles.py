from __future__ import annotations

from dataclasses import dataclass

from pymatgen.core import Lattice, Structure


class SmilesConversionError(RuntimeError):
    """Raised when SMILES cannot be converted to a structure."""


@dataclass
class SmilesCrystalOptions:
    min_box: float = 10.0
    padding: float = 4.0


def smiles_to_structure(smiles: str, options: SmilesCrystalOptions | None = None) -> Structure:
    """Convert a SMILES string into a toy molecular crystal structure.

    Notes
    -----
    This is a lightweight approximation for rapid prototyping:
    1) Build a 3D molecule with RDKit.
    2) Place a single molecule in a cubic periodic box.

    It is suitable for generating a synthetic XRD pattern in a web demo,
    not for production-grade crystal prediction.
    """
    options = options or SmilesCrystalOptions()

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception as exc:  # pragma: no cover
        raise SmilesConversionError(
            "RDKit is required for SMILES input. Install dependency: rdkit-pypi"
        ) from exc

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SmilesConversionError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
        raise SmilesConversionError("Failed to generate 3D conformer from SMILES")
    AllChem.UFFOptimizeMolecule(mol, maxIters=300)

    conf = mol.GetConformer()
    species = []
    coords = []
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        species.append(atom.GetSymbol())
        coords.append([p.x, p.y, p.z])

    xs, ys, zs = zip(*coords)
    extent = max(max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))
    box = max(options.min_box, extent + options.padding)

    # Center molecule into periodic box and convert to fractional coordinates.
    centered = []
    for x, y, z in coords:
        centered.append([(x - min(xs)) + (box - (max(xs) - min(xs))) / 2,
                         (y - min(ys)) + (box - (max(ys) - min(ys))) / 2,
                         (z - min(zs)) + (box - (max(zs) - min(zs))) / 2])

    frac = [[x / box, y / box, z / box] for x, y, z in centered]
    return Structure(Lattice.cubic(box), species, frac)
