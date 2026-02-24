from __future__ import annotations

import argparse
from pathlib import Path

from .pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run crystal structure prediction (CSP) scoring + stable lattice selection + XRD simulation"
    )
    parser.add_argument("--input", required=True, help="Input structure file (e.g., CIF, POSCAR)")
    parser.add_argument("--output-dir", default="outputs", help="Directory to save report and XRD artifacts")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    report = run_pipeline(Path(args.input), Path(args.output_dir))
    print(f"Done. Report: {report}")


if __name__ == "__main__":
    main()
