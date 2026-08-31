# YAML Model Reference for BRNS

This guide describes the structure of YAML model files for BRNS and explains which blocks are read by the mapper, which values are optional, which default values exist, and how comments should be written.

> Recommendation: place these blocks at the beginning of the YAML file, even if they are not currently consumed by the generator logic.

```yaml
metadata:
  name: "Example reactive transport model"
  description: "Generic BRNS model template for a coupled reaction-transport system"
  author: "Model author"
  version: "1.0"
  date: "2026-08-25"
  reference: "Internal model definition"
  institution: "Research group"
  email: "model.author@example.org"
  keywords:
    - "biogeochemistry"
    - "reaction network"
    - "transport"
    - "example"

# ======================================================================
# UNITS (centralized)
# ======================================================================

units:
  L: cm
  M: g
  T: a
```

These entries provide model metadata and a central unit definition. They are not strictly required for the pipeline, but they are very useful for documentation, traceability, and consistent model maintenance.

---

## 1. General structure of a YAML file

A BRNS YAML file is organized into logically separate blocks. A typical order is:

```yaml
metadata:
  ...

units:
  ...

output:
  ...

initial_conditions:
  ...

species:
  - name: "..."
    type: "dissolved"
    ...

parameters:
  biogeochemical:
    - name: "..."
      value: ...

  stoichiometry:
    x: ...
    y: ...

  physical:
    ...

  physical_flags:
    ...

computed_values:
  ...

grid:
  nnodes: ...
  type: ...

reactions:
  - id: 1
    name: "..."
    rate: "..."
    stoichiometry:
      ...

advanced:
  ...
```

Other block can be added for documentation purpose. These additional blocks are ignored.
---

## 2. Which blocks are read by the mapper?

The mapper primarily works with the blocks required for ACG / Fortran generation.

### 2.1 `species`
- Each species is a model component with a name and physical meaning.
- The order in the YAML is usually treated as a fixed order.
- The mapper builds the variable list and index mapping from this block.

Example:

```yaml
species:
  - name: "o2"
    type: "dissolved"
    transport: true
    bc_upper_type: 0
    bc_upper_value: 132.0e-6
    bc_lower_type: 1
    bc_lower_value: 0.0
    transport_D0: 296.0
    transport_alpha: 0.06
    transport_tortuosity: 0.0
    init_value: 180.0e-6
```

Typical fields and expected data types:

| Field | Type | Meaning / default |
| --- | --- | --- |
| `name` | String | Name of the species |
| `type` | String | Either `"dissolved"` or `"solid"` |
| `transport` | Boolean | Whether the species can be transported in the flow. Optional. Default: `true` |
| `bc_upper_type` | Integer | Type of the upper boundary condition. Code values: `0 = Dirichlet`, `1 = Neumann`, `2 = flux` |
| `bc_upper_value` | Float | Value of the upper boundary condition |
| `bc_lower_type` | Integer | Type of the lower boundary condition. Code values: `0 = Dirichlet`, `1 = Neumann`, `2 = flux` |
| `bc_lower_value` | Float | Value of the lower boundary condition |
| `transport_D0` | Float | Diffusion coefficient |
| `transport_alpha` | Float | Derivative of Diffusion coefficient by temperature. Default = 0 |
| `transport_tortuosity` | Float | Tortuosity. Default = 1 |
| `init_value` | Float | Initial value. This is only evaluated when `initial_conditions.mode = 2` |

### 2.2 `parameters`
The `parameters` section contains model parameters. The mapper reads these entries and builds the ACG parameter lists.

The `parameters` block may contain arbitrary first-level subsections. Each subsection can use either a mapping (parameter_name: value) or a list of objects with name and value fields. Common subsections are biogeochemical, physical, and physical_flags.

Example:

```yaml
parameters:
  biogeochemical:
    - name: "kfox"
      value: 0.221
      description: "Organic matter degradation rate constant"

  physical:
    al: 1.0e-5
    q0: 0.0
    w0: 0.398
    Db0: 27.5
    por0: 0.85
    area0: 1.0
    t_celsius: 10.0
    salin: 35.0
    delt: 10.0
    depthmax: 100.0
    endt: 100000.0

  physical_flags:
    iq: 0
    iw: 0
    iDb: 1
    ipor: 0
    igrid: 1
    iarea: 0
```

Typical fields and expected data types:

| Field | Type | Meaning / default |
| --- | --- | --- |
| `biogeochemical` | List of objects | Biogeochemical parameter definitions. Each item has `name` (String), `value` (Float or Integer), and optional `description` (String). |
| `stoichiometry` | Mapping | Global stoichiometric constants such as `x`, `y`, `z`, `s_dens`; values are Float or Integer. |
| `physical.al` | Float | Cross-sectional area or pore-system area term. |
| `physical.q0` | Float | Vertical advective velocity / seepage velocity. Typical unit: cm/year. `0` means no advection. |
| `physical.w0` | Float | Bioturbation mixing depth or effective mixing velocity. |
| `physical.Db0` | Float | Bioturbation diffusion coefficient. |
| `physical.por0` | Float | Reference porosity. Typical value: `0.3`, i.e. 30% pore space. |
| `physical.area0` | Float | Reference cross-sectional area of the system. Typical value: `1.0` cm². |
| `physical.t_celsius` | Float | Ambient temperature used for temperature corrections. |
| `physical.salin` | Float | Salinity used for salinity-dependent corrections.  |
| `physical.delt` | Float | Time step size for integration.  |
| `physical.depthmax` | Float | Maximum simulated depth / vertical domain height.  |
| `physical.endt` | Float | End of simulation time. |
| `physical_flags.iq` | Integer | Advection switch. `0` = constant advection, `1` = spatially variable advection. Default: `0` if not set. |
| `physical_flags.iw` | Integer | Bioturbation switch. `0` = constant mixing, `1` = spatially variable mixing. Default: `0` if not set. |
| `physical_flags.iDb` | Integer | Bioturbation diffusion switch. `0` = constant diffusion, `1` = spatially variable diffusion profile. Default: `0` if not set. |
| `physical_flags.ipor` | Integer | Porosity switch. `0` = constant porosity, `1` = spatially variable porosity. Default: `0` if not set. |
| `physical_flags.igrid` | Integer | Grid-type switch. `0` = uniform grid, `1` = logarithmic / compressed grid. Default: `0` if not set. |
| `physical_flags.iarea` | Integer | Cross-sectional area switch. `0` = constant area, `1` = spatially variable area. Default: `0` if not set. |

Default behavior:
- Unspecified `physical_flags` values are typically treated as `0` by the mapper or orchestrator.
- `ic` is usually derived automatically from stoichiometric relationships and should not be entered manually.
- If a species has no explicit `transport` value, it defaults to `true`.

### 2.3 `grid`
This block is used when the model is spatially discretized.

Example:

```yaml
grid:
  nnodes: 111
  type: 1
```

Typical fields and expected data types:

| Field | Type | Meaning / default |
| --- | --- | --- |
| `nnodes` | Integer | Number of grid nodes / discretization points |
| `type` | Integer or String | Grid type or discretization mode |

Default behavior:
- If not provided, a model-specific fallback may be used by the pipeline.
- In practice, `nnodes` should be set explicitly for spatial models.

### 2.4 `reactions`
A reaction block is the core of the model.

Example (kinetic reaction):

```yaml
reactions:
  - id: 1
    name: "aerobic_respiration"
    rate: "kch2o * fo2"
    stoichiometry:
      ch2o: -1
      o2: "- SD*(x+2*y)/x"
      no3: "SD*(y/x)"
      co2: "(x+y+2*z)/x * SD"
      hco3: "- (y+2*z)/x * SD"
```

Example (equilibrium reaction):

```yaml
reactions:
  - id: 13
    name: "sulfide_dissociation"
    description: "H2S <-> HS- + H+"
    equilibrium: true
    rate: "F13 - B13"
    rate_components:
      F13: "kfs * hplus * hs"
      B13: "kbs * h2s"
    equilibrium_constraint: "hplus * hs - keqs * h2s"
    stoichiometry:
      h2s: 1
      hs: -1
      hplus: -1
```

Note: In real BRNS models, equilibrium reactions may still carry a `rate` / `rate_components` form as a forward-backward expression, but the decisive algebraic condition is `equilibrium_constraint`. The generator interprets the reaction as an equilibrium constraint, not as a normal kinetic rate law.

Typical fields and expected data types:

| Field | Type | Meaning / default |
| --- | --- | --- |
| `id` | Integer | Unique reaction identifier and sort key. IDs do not need to be consecutive. |
| `name` | String | Reaction name |
| `rate` | String | Rate expression as a symbolic formula |
| `rate_components` | Mapping | Optional helper expressions used by the reaction rate; values are Strings |
| `stoichiometry` | Mapping | Species coefficients for the reaction; values can be Float, Integer, or expression Strings |
| `equilibrium` | Boolean | Whether the reaction is treated as an equilibrium reaction. Default: `false` |
| `equilibrium_constraint` | String or expression | Required when `equilibrium: true` |
| `output` | Boolean | Optional output flag for the reaction; if `true`, the reaction result is written to output |
| `output_filename` | String | Optional file name for reaction output |

### 2.5 `advanced`
This block is partially evaluated but not all content is mandatory.

Example:

```yaml
advanced:
  preprocessing:
    gaussian_elimination: true
    symbolic_simplification: true

  optimization:
    enabled: false

  timestep_parameters:
    tstep: 1.0d-5

```

Typical fields and expected data types:

| Field | Type | Meaning / default |
| --- | --- | --- |
| `advanced.preprocessing.gaussian_elimination` | Boolean | Enables symbolic elimination preprocessing. Optional. |
| `advanced.preprocessing.symbolic_simplification` | Boolean | Enables symbolic simplification. Optional. |
| `advanced.optimization.enabled` | Boolean | Enables optimization steps. Optional. |
| `advanced.timestep_parameters.dtold` | Float | Previous timestep initialization parameter. Default: `1.d-6` |
| `advanced.timestep_parameters.ttol` | Float | Relative tolerance for timestep updates. Default: `5.d-2` |
| `advanced.timestep_parameters.tstep` | Float | Initial timestep size. Default: `1.0d0` |
| `advanced.timestep_parameters.maxconc` | Float | Max concentration threshold. Default: `0.d0` |

### 2.6 `initial_conditions`

Example:

```yaml
initial_conditions:
  mode: 3
  listinput: [1, 2, 3]
  file_in_names: ["o2", "no3", "nh4"]
```

The ACG generator supports three explicit initialization modes, implemented in `acg12`:

| Mode | Meaning in generator | Effect on `species.init_value` |
| --- | --- | --- |
| `1` | Read initial profiles from `initialconc.txt` | `species.init_value` is not used.  |
| `2` | Constant initial conditions everywhere in the domain | Reads `species.init_value`. |
| `3` | Read initial conditions from files | The generator reads one profile file per species with file name species.name with ending `.inp`; example "co2.inp" ; `species.init_value` is not used. |

Typical fields and expected data types:

| Field | Type | Meaning / default |
| --- | --- | --- |
| `mode` | Integer | Global initial-condition mode. Supported generator modes are `1`, `2`, and `3`. |
| `listinput` | List of Integer | Optional explicit species indices used in mode `3`. |
| `file_in_names` | List of String | Optional file names used in mode `3`; if given without extension, the generator appends `.inp`. |

### 2.7 `output`
This block is typically used for generator output configuration.

Example:

```yaml
output:
  directory: "../generated_fortran/multiple_species"
  prefix: "multiple_species"
  timing:
    start: 1.0
    interval: 1.0
```

Typical fields and expected data types:

| Field | Type | Meaning / default |
| --- | --- | --- |
| `directory` | String | Output directory for generated files |
| `prefix` | String | Prefix applied to generated output files |
| `timing.start` | Float | Start time for output events |
| `timing.interval` | Float | Interval between output writes |
| `species[].output` | Boolean | Optional per-species output flag; if enabled, species is written to output |
| `species[].output_filename` | String | Optional custom file name for a species output |
| `reactions[].output` | Boolean | Optional per-reaction output flag; if enabled, reaction output is written |
| `reactions[].output_filename` | String | Optional custom file name for reaction output |

Default behavior:
- This block is optional in many models and may be omitted if the output path is managed elsewhere.

---

## 3. How to add comments

YAML comments begin with `#`.

Example:

```yaml
# This is a comment
species:
  - name: "o2"  # inline comment for a single entry
    type: "dissolved"
```

Rules:
- `#` starts a comment that lasts until the end of the line.
- Comments can stand alone on a line.
- Comments can also appear at the end of a line.