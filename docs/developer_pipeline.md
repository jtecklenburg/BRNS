# Introduction for BRNS developer - Python part 

This page gives a compact overview of how the main Python components cooperate to turn a YAML model into generated Fortran code.

## 1. Orchestration: `acg_orchestrator.py`

`ACGOrchestrator` is the entry point for the full pipeline. It loads and validates the YAML model, resolves formula expressions, and orchestrates each downstream step. In practice it does the following:

- reads the model definition from a YAML file,
- validates required sections such as `species`, `reactions`, and `parameters`,
- evaluates parameter formulas and symbolic expressions,
- builds the mapper input and the reduced reaction system,
- calls the code generator to emit the Fortran files.

This class is the project-level coordinator: it owns the workflow, but it does not implement the numerical details itself.

## 2. YAML-to-structure mapping: `yaml_to_acg_mapper.py`

`YAMLtoACGMapper` translates the YAML schema into the internal arrays and matrix structures expected by the later ACG code generation.

Its main tasks are:

- preserve species ordering from the YAML definition,
- build species indices and variable names for Fortran-compatible indexing,
- assemble the parameter vectors used in `acg3()`-style outputs,
- create the stoichiometric matrix from reaction definitions,
- expose model dimensions such as number of dissolved species, solids, and reactions.

In other words, the mapper is responsible for turning a declarative model description into the structured data needed by the mathematical reduction and code generation layers.

## 3. Algebraic reduction: `gaussian_elimination.py`

`gaussian_elimination.py` implements the Maple-style reaction reduction steps (`p0` to `p10`). This is where the reaction network is transformed into a reduced, solvable form.

The main ideas are:

- create reaction and variable symbols,
- reorder equilibrium and non-equilibrium reactions consistently,
- build the coefficient matrix from the network equations,
- perform Gaussian elimination and pivoting to extract independent stoichiometric relations,
- identify conservation rows and reduced residual expressions,
- compute the Jacobian terms required for the Newton solver.

This module is the numerical heart of the pipeline: it converts the raw stoichiometric model into a reduced system that the generated code can solve efficiently.

## 4. Code generation: `acg.py`

`ACGModule` is the Fortran code generator. It consumes the reduced symbolic data and emits problem-specific Fortran routines, such as residuals, Jacobians, and parameter blocks.

The generator uses SymPy expressions and `macrofor` to create BRNS-compatible source files. The key responsibility is not to solve the model, but to format the reduced symbolic equations into the exact Fortran structures expected by the BRNS runtime.

## How they work together

The pipeline is intentionally layered:

1. `ACGOrchestrator` orchestrates the process.
2. `YAMLtoACGMapper` reshapes the YAML content into model arrays and stoichiometric structures.
3. `gaussian_elimination.py` reduces the reaction network into a minimal symbolic system.
4. `ACGModule` emits the final Fortran code from that reduced representation.

The result is a single workflow: YAML model definition → validated data → structured stoichiometry → reduced algebraic system → generated Fortran solver code.
