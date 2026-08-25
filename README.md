# BRNS – Biogeochemical Reaction Network Simulator

![GitHub License](https://img.shields.io/github/license/jtecklenburg/BRNS)

BRNS (**Biogeochemical Reaction Network Simulator**) is a flexible modelling framework for simulating coupled, multi-component reaction networks in porous media. It combines an automatic code generator (ACG) with a compiled Fortran simulation core, allowing to define arbitrarily complex kinetic and equilibrium reaction networks without hand-writing solver code for every new problem.

---

## What is BRNS?

BRNS was originally developed to simulate early diagenetic processes in marine sediments, where large numbers of biogeochemical reactions (organic matter degradation, redox cycling of iron, manganese and sulfur, nutrient cycling, pH buffering, etc.) need to be solved simultaneously with transport processes such as diffusion, advection, bioturbation and bioirrigation.

Rather than requiring users to write custom solver code for each reaction network, BRNS uses an **automatic code generator** to translate a user-defined reaction network into ready-to-compile simulation code. This makes BRNS a flexible tool for anyone who needs to couple a reaction network to a transport model, without becoming a numerical-methods specialist first.

## What has BRNS been used for?

Since its introduction, BRNS and its derivatives have been applied to a wide range of subsurface and sediment biogeochemistry problems, including:

- Early diagenesis and redox cycling of carbon, nitrogen, iron, manganese and sulfur in marine sediments [Regnier et al., 2002; Thullner et al., 2005]
- pH dynamics and proton budgets in aquatic sediments [Jourabchi et al., 2005]
- Anaerobic oxidation of methane (AOM) and sulfate-methane transition zones [Blouet et al., 2021]
- Methane hydrate stability and benthic methane escape under permafrost thaw [Sivan et al., 2020]
- Coupling with multidimensional flow and transport codes (e.g. OpenGeoSys, OpenFOAM) for groundwater and pore-scale reactive transport [Centler et al., 2010; Golparvar et al., 2024]
- Redox transformations of trace metals and contaminants (e.g. uranium, arsenic) [e.g. sensitivity analysis of uranium reduction, 2025]

A curated bibliography of BRNS-related publications and applications is being can be found in the project documentation (`docs/`). 

---

## How BRNS works

BRNS follows a three-step workflow:

1. **Define a reaction network.**
   You specify the species, kinetic and equilibrium reactions, and stoichiometry of your problem in a configuration file.

2. **Generate problem-specific Fortran code (ACG).**
   The Automatic Code Generator (ACG) reads your reaction network definition and generates the Fortran source files needed to solve that specific system (mass balances, Jacobian, Newton–Raphson solver routines, etc.).

3. **Compile the generated Fortran code.**
   Compiling the generated files produces a problem-specific simulator. You can use it in two ways:
   - **Stand-alone console program** — run BRNS directly as a batch reaction-transport solver for 1D problems.
   - **Compiled library** — link BRNS against an external transport code (e.g. a groundwater flow and transport model), so BRNS handles the reaction step while the host code handles transport.

### What's new in this version

Earlier versions of BRNS generated Fortran code from a **Maple** worksheet, which required a Maple license and was difficult to script or version-control. **This version replaces the Maple-based ACG with a Python-based code generator.** Reaction networks are now defined in a plain-text **YAML** file instead of a Maple notebook, making the workflow:

- free of proprietary software dependencies for the code-generation step,
- easier to script, version-control and share, and
- more approachable for users without a Maple background.

---

## Installation

BRNS requires two components to be installed, in the following order:

### 1. Install `macrofor`

```bash
pip install "macrofor @ git+https://github.com/jtecklenburg/macrofor.git"
```

### 2. Install BRNS

Clone the repository locally to access critical workflow files:

- Fortran files: Required to compile simulator.
- Bash scripts: Automated build pipelines.
- Jupyter notebooks: Ready-to-use examples.

```bash
git clone https://github.com/jtecklenburg/BRNS.git
cd BRNS
pip install -e .
```

---

## Usage

### 1. Define your reaction network

Create a YAML file describing your species, reactions and transport processes. The following example describes the conversion of species $A$ into $B$ under constant advective flow and diffusion:

$$
A \xrightarrow{k} B,
\qquad r = kA
$$

where

$$
k = 2k_{\mathrm{ref}}.
$$

At the upper boundary, only species $A$ is supplied at a concentration of 1,
while the concentration of $B$ is zero. Initially, the concentrations of both
species are zero throughout the entire model domain.

```yaml
units:
  L: cm
  M: arbitrary
  T: a

grid:
  nnodes: 51
  type: 0

initial_conditions:
  mode: 2

species:
  - name: a
    type: dissolved
    bc_upper_type: 0
    bc_upper_value: 1.0
    bc_lower_type: 1
    bc_lower_value: 0.0
    init_value: 0.0
    transport_D0: 1.0e-9

  - name: b
    type: dissolved
    bc_upper_type: 0
    bc_upper_value: 0.0
    bc_lower_type: 1
    bc_lower_value: 0.0
    init_value: 0.0
    transport_D0: 1.0e-9

parameters:
  kinetics:
    k_ref: 1.0e-3
    k: "2*k_ref"
  physical:
    w0: 1.0e-4
    por0: 0.3
    depthmax: 0.1
    delt: 10.0
    endt: 1.0e4

reactions:
  - id: 1
    rate: "k*a"
    stoichiometry: {a: -1, b: 1}
```

### 2. Run the code generator

#### Option A: script based usage to build and run a single YAML model (Linux only)

Requirements: Python with BRNS, gfortran.

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

#### Option B: Python API Usage

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

Next compile and run BRNS:
1. Create a build directory.
2. Copy static Fortran Files from folder BRNSPackage/FortranFiles into your build directory.
3. Copy the Python generated Fortran Files into your build directory.
4. Compile the Fortran Files in the build directory with your favorite Fortran compiler as stand-alone .console program (or as a library).
5. Add files with initial conditions into the folder with the executable (optional) and run executable.
6. Evaluate results with notebooks/plot_results.ipynb.

---

## Documentation

The documentation can be found under `docs/` or [online](https://jtecklenburg.github.io/BRNS/).

---

## References

The original BRNS concept and code base were introduced and developed in the following foundational publications:

- Aguilera, D.R., Jourabchi, P., Spiteri, C. and Regnier, P.  2005.  A knowledge-based reactive transport approach for the simulation of biogeochemical dynamics in Earth systems. Geochemistry Geophysics Geosystems 6(7), Q07012.
- Regnier, P., O'Kane, J.P., Steefel, C.I. and Vanderborght, P.  2002.  Modeling complex multi-component reactive-transport systems: towards a simulation environment based on the concept of a Knowledge Base. Applied Mathematical Modelling 26(9), 913-927.

A comprehensive bibliography of BRNS-based studies can be found in the [project documentation](https://jtecklenburg.github.io/BRNS/).

---

## Contributing

When you like to contribute to this project, please contact Martin.Thullner@bgr.de.

## Acknowledgments

This project builds upon the work and concepts of several contributors:

- P. Regnier & co-workers: Basic RTM & core concept.
- Florian Centler & Martin Thullner: Coupling and generalization concept.
- Florian Centler: Maple 10+ and DLL versions.
- Jan Tecklenburg: Python version implementation.

The Python based version of BRNS was funded by ptj, Förderprogramm: Geoforschung und Nachhaltigkeit (GEO:N), Förderkennzeichen: 03G0937B (BMFTR)

## License

MIT License. See [LICENSE](LICENSE).