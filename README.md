# ACG-BRNS

Automatic code generation for BRNS: convert YAML-based biogeochemical models into BRNS-compatible Fortran source files.

## BRNS stands for **Biogeochemical Reaction Network Simulator**.

It is a process-based modeling framework used to simulate coupled transport-reaction systems in sediments and aquatic geochemical environments. In practice, BRNS solves multi-species reaction-transport equations (advection, diffusion/bioturbation, and biogeochemical source/sink terms) and relies on user/model-specific Fortran routines for rates, Jacobians, boundary conditions, and initial conditions.

ACG-BRNS generates these model-specific Fortran files automatically from YAML model definitions.

## Overview

ACG-BRNS is a Python pipeline for BRNS model development.

It provides:
- YAML model parsing and validation
- symbolic formula handling (SymPy)
- Jacobian/residual generation
- Gaussian-elimination-based preprocessing
- Fortran code generation (BRNS-compatible fixed-format files)

The project is intended as a Maple-free workflow for building BRNS model code from reproducible YAML configurations.

## Current Workflow (Recommended)

1. Define or edit a model YAML in [models](models).
2. Generate Fortran files (Notebook or script).
3. Build and run the generated model.
4. Compare or plot results with notebooks.

## Installation

From source:

```bash
cd BRNSPackage
pip install -e .
```

Optional extras:

```bash
# Dev tools
pip install -e ".[dev]"

# Notebook support
pip install -e ".[notebooks]"

# Everything
pip install -e ".[dev,notebooks]"
```

## CLI / Script-Based Usage

### Build and run a single YAML model

Use [build_python.sh](build_python.sh):

```bash
./build_python.sh -c ./models/single_species_example.yaml -i ./path/to/input_files
```

Options:
- `-c` YAML file (required)
- `-i` input directory containing `.inp` files (required)
- `-o` output root directory (optional, default `./build_output`)
- `-n` model/build name (optional)

This script performs:
1. YAML → Fortran generation via `ACGOrchestrator`
2. build directory preparation
3. Fortran compilation (`gfortran`)
4. model execution and result collection (`.dat`)

## Python API Usage

High-level YAML pipeline:

```python
from acg_brns.acg_orchestrator import ACGOrchestrator

orchestrator = ACGOrchestrator(
        yaml_path="models/single_species_example.yaml",
        output_dir="generated_fortran/single_species",
        verbose=True,
)

summary = orchestrator.generate()
print(summary)
```

Low-level generation APIs are still available via `ACGModule` for advanced/custom flows.

## Notebook Guide

Main notebooks in [notebooks](notebooks):

- [generate_fortran_from_yaml.ipynb](notebooks/generate_fortran_from_yaml.ipynb)
    - focused YAML → Fortran generation walkthrough
- [generate_fortran_debug.ipynb](notebooks/generate_fortran_debug.ipynb)
    - advanced/debug-oriented pipeline inspection
- [compare_results.ipynb](notebooks/compare_results.ipynb)
    - compare two result directories (e.g., reference vs Python)
- [plot_results.ipynb](notebooks/plot_results.ipynb)
    - plot results from a single version
- [gaussian_elimination_tutorial.ipynb](notebooks/gaussian_elimination_tutorial.ipynb)
    - preprocessing and elimination concepts

## Typical Output Files

Generated BRNS user files usually include (model-dependent):
- `rates.f`
- `jacobian.f`
- `residual.f`
- `ssrates.f`
- `boundaries.f`
- `initialcond.f`
- `output.f`
- include/config files such as `common_geo.inc`

## Requirements

Core:
- Python >= 3.8
- SymPy >= 1.12
- NumPy >= 1.20

Build/runtime tools (for full BRNS run flow):
- `gfortran`
- BLAS/LAPACK libraries

## Project Structure

- [acg_brns](acg_brns): Python generation pipeline
- [models](models): YAML model definitions
- [generated_fortran](generated_fortran): generated Fortran outputs
- [build](build) / `build_output`: compiled runs and results
- [notebooks](notebooks): interactive workflows (generation/diagnostics/plotting)
- [reference_fortran](reference_fortran): reference model code

## Development

```bash
# Tests
pytest

# Coverage
pytest --cov=acg_brns --cov-report=html

# Formatting
black acg_brns/ tests/

# Linting
flake8 acg_brns/ tests/
```

## License

MIT License. See [LICENSE](LICENSE).

## Literature Hints

- Regnier, P., Jourabchi, P., & Slomp, C.P. (2003). *Reactive-transport modeling as a technique for understanding coupled biogeochemical processes in surface and subsurface environments*. *Netherlands Journal of Geosciences*, 82(1), 5–18.

## Acknowledgments

ACG-BRNS is a Python evolution of the Maple-based BRNS ACG workflow and integrates with the broader thermooptiplan ecosystem.

Funded by ptj, Förderprogramm: Geoforschung und Nachhaltigkeit (GEO:N), Förderkennzeichen: 03G0937B (BMFTR)
