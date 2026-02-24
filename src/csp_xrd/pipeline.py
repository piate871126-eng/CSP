from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
import json
from typing import Iterable

import numpy as np
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


BRAVAIS_TEMPLATES = {
    "cubic": lambda a, b, c: Lattice.cubic((a * b * c) ** (1 / 3)),
    "tetragonal": lambda a, b, c: Lattice.tetragonal(np.sqrt(a * b), c),
    "orthorhombic": lambda a, b, c: Lattice.orthorhombic(a, b, c),
    "hexagonal": lambda a, b, c: Lattice.hexagonal(np.sqrt(a * b), c),
    "rhombohedral": lambda a, b, c: Lattice.rhombohedral((a * b * c) ** (1 / 3), 75),
    "monoclinic": lambda a, b, c: Lattice.monoclinic(a, b, c, 105),
    "triclinic": lambda a, b, c: Lattice.from_parameters(a, b, c, 82, 95, 102),
}


@dataclass
class CandidateResult:
    template: str
    predicted_crystal_system: str
    energy_score: float
    volume: float
    spacegroup_symbol: str
    spacegroup_number: int


def _short_range_repulsion(structure: Structure) -> float:
    dmat = structure.distance_matrix
    species = structure.species
    n = len(structure)
    score = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            d = dmat[i, j]
            ri = species[i].element.atomic_radius or 1.2
            rj = species[j].element.atomic_radius or 1.2
            cutoff = 0.9 * (ri + rj)
            if d < cutoff:
                score += ((cutoff - d) / cutoff) ** 2 * 100.0
            score += 1.0 / (d**6 + 1e-6)
    return score


def _strain_penalty(structure: Structure) -> float:
    lat = structure.lattice
    angles = np.array([lat.alpha, lat.beta, lat.gamma])
    lengths = np.array([lat.a, lat.b, lat.c])
    angle_term = np.mean(((angles - 90.0) / 90.0) ** 2)
    ratio_term = np.var(lengths / np.mean(lengths))
    return 10.0 * angle_term + 5.0 * ratio_term


def score_structure(structure: Structure) -> float:
    return _short_range_repulsion(structure) + _strain_penalty(structure)


def generate_candidates(structure: Structure, systems: Iterable[str] | None = None) -> list[tuple[str, Structure]]:
    systems = list(systems or BRAVAIS_TEMPLATES.keys())
    a, b, c = structure.lattice.abc
    frac = structure.frac_coords
    species = structure.species
    candidates: list[tuple[str, Structure]] = []
    for system in systems:
        lattice = BRAVAIS_TEMPLATES[system](a, b, c)
        candidate = Structure(lattice, species, frac)
        candidates.append((system, candidate))
    return candidates


def evaluate_candidates(candidates: Iterable[tuple[str, Structure]]) -> list[CandidateResult]:
    results: list[CandidateResult] = []
    for template, candidate in candidates:
        sga = SpacegroupAnalyzer(candidate, symprec=0.1)
        results.append(
            CandidateResult(
                template=template,
                predicted_crystal_system=sga.get_crystal_system(),
                energy_score=float(score_structure(candidate)),
                volume=float(candidate.volume),
                spacegroup_symbol=sga.get_space_group_symbol(),
                spacegroup_number=int(sga.get_space_group_number()),
            )
        )
    return sorted(results, key=lambda x: x.energy_score)


def _simulate_xrd(structure: Structure, output_dir: Path, wavelength: str = "CuKa") -> dict:
    calc = XRDCalculator(wavelength=wavelength)
    pattern = calc.get_pattern(structure, two_theta_range=(5, 90))
    csv_path = output_dir / "xrd_pattern.csv"
    with csv_path.open("w", encoding="utf-8") as f:
        f.write("two_theta,intensity\n")
        for t, i in zip(pattern.x, pattern.y):
            f.write(f"{t:.4f},{i:.4f}\n")

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.vlines(pattern.x, 0, pattern.y, linewidth=1.5)
    ax.set_xlabel("2θ (deg)")
    ax.set_ylabel("Intensity (a.u.)")
    ax.set_title("Simulated XRD")
    fig.tight_layout()
    png_path = output_dir / "xrd_pattern.png"
    fig.savefig(png_path, dpi=160)
    plt.close(fig)

    return {"csv": str(csv_path), "plot": str(png_path), "peak_count": len(pattern.x)}


def run_pipeline(structure_path: str | Path, output_dir: str | Path) -> Path:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    structure = Structure.from_file(str(structure_path))
    scored = evaluate_candidates(generate_candidates(structure))
    best = scored[0]

    best_lattice = BRAVAIS_TEMPLATES[best.template](*structure.lattice.abc)
    best_structure = Structure(best_lattice, structure.species, structure.frac_coords)

    best_cif = output_dir / "predicted_best_structure.cif"
    best_structure.to(filename=str(best_cif))

    xrd_info = _simulate_xrd(best_structure, output_dir)

    result = {
        "input_structure": str(structure_path),
        "candidate_count": len(scored),
        "best": asdict(best),
        "all_candidates": [asdict(s) for s in scored],
        "predicted_structure_cif": str(best_cif),
        "xrd": xrd_info,
    }

    report = output_dir / "report.json"
    report.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return report
