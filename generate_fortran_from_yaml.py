#!/usr/bin/env python3
"""Generate Fortran code from a YAML model configuration.

Runs the ACGOrchestrator pipeline for a given YAML model file and writes the
generated Fortran source files to an output directory.

Usage:
    python generate_fortran_from_yaml.py [-c YAML_PATH] [-o OUTPUT_DIR]

If no arguments are given, the defaults defined below (YAML_PATH /
OUTPUT_DIR) are used.
"""

import argparse
import warnings
from pathlib import Path
from pprint import pprint

from acg_brns.acg_orchestrator import ACGOrchestrator

ROOT = Path(__file__).resolve().parent

# Default paths, used when no command-line arguments are given.
YAML_PATH = ROOT / "models" / "minimal" / "simple.yaml"
OUTPUT_DIR = ROOT / "generated_fortran" / "simple"

# To hard-code different defaults instead of using the command line, edit the
# values above and comment out the argparse block in main().


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate Fortran code from a YAML model configuration."
    )
    parser.add_argument(
        "-c", "--config",
        type=Path,
        default=YAML_PATH,
        help=f"Path to the YAML model file (default: {YAML_PATH})",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=OUTPUT_DIR,
        help=f"Output directory for generated Fortran files (default: {OUTPUT_DIR})",
    )
    return parser.parse_args()


def generate_fortran(yaml_path: Path, output_dir: Path):
    summary = None
    if yaml_path.exists():
        try:
            orchestrator = ACGOrchestrator(
                yaml_path=str(yaml_path),
                output_dir=str(output_dir),
                verbose=True,
            )
            summary = orchestrator.generate()
            pprint(summary)
        except Exception as exc:
            warnings.warn(f"Code generation failed: {exc}")
    else:
        print(f"⚠ Skipped: YAML file not found: {yaml_path}")
    return summary


def main():
    args = parse_args()
    generate_fortran(args.config, args.output)


if __name__ == "__main__":
    main()
