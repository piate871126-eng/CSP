from pathlib import Path
import json

from pymatgen.core import Lattice, Structure

from csp_xrd.pipeline import run_pipeline


def test_pipeline_end_to_end(tmp_path: Path):
    structure = Structure(
        Lattice.cubic(5.64),
        ["Na", "Cl"],
        [[0, 0, 0], [0.5, 0.5, 0.5]],
    )
    input_cif = tmp_path / "nacl.cif"
    structure.to(filename=str(input_cif))

    report_path = run_pipeline(input_cif, tmp_path / "out")
    assert report_path.exists()

    report = json.loads(report_path.read_text())
    assert report["candidate_count"] == 7
    assert Path(report["predicted_structure_cif"]).exists()
    assert Path(report["xrd"]["csv"]).exists()
    assert report["xrd"]["peak_count"] > 0
