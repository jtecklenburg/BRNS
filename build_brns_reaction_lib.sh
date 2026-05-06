#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat <<'EOF'
Build a Linux BRNS reaction shared library from a YAML model.

Usage:
  ./build_brns_reaction_lib.sh -c MODEL.yaml|MODEL_DIR|MODEL_NAME [-o build_root] [-n model_name] [-p python3] [-f gfortran]

Options:
  -c  YAML model configuration, model directory, or model name under ./models (required)
  -o  Output root directory (default: ./build)
  -n  Build/model name (default: derived from YAML filename)
  -p  Python executable (default: python3)
  -f  Fortran compiler (default: gfortran)
  -h  Show this help

Outputs:
  build/<model>/generated/           YAML-derived Fortran files
  build/<model>/library_src/         staged library sources
  build/<model>/lib/                 compiled .so
  build/<model>/logs/compile.log     compiler log
  build/<model>/brns_library_metadata.json
EOF
}

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FORTRAN_COMMON="$ROOT_DIR/BRNSPackage/FortranFiles"
INVOKEBRNS_SOURCE="$FORTRAN_COMMON/BrnsDll/invokebrns.f"
PYTHON_BIN="python3"
FC="gfortran"
OUTPUT_ROOT="$ROOT_DIR/build"
YAML_CONFIG=""
MODEL_NAME=""

while getopts ":c:o:n:p:f:h" opt; do
  case "$opt" in
    c) YAML_CONFIG="$OPTARG" ;;
    o) OUTPUT_ROOT="$OPTARG" ;;
    n) MODEL_NAME="$OPTARG" ;;
    p) PYTHON_BIN="$OPTARG" ;;
    f) FC="$OPTARG" ;;
    h) show_help; exit 0 ;;
    :) echo "ERROR: Option -$OPTARG requires an argument" >&2; exit 1 ;;
    \?) echo "ERROR: Unknown option -$OPTARG" >&2; exit 1 ;;
  esac
done

if [[ -z "$YAML_CONFIG" ]]; then
  echo "ERROR: YAML config (-c) is required" >&2
  show_help
  exit 1
fi

# Allow shorthand model name from ./models
if [[ ! -e "$YAML_CONFIG" && -d "$ROOT_DIR/models/$YAML_CONFIG" ]]; then
  YAML_CONFIG="$ROOT_DIR/models/$YAML_CONFIG"
fi

# If -c points to a model directory, pick YAML from that directory
if [[ -d "$YAML_CONFIG" ]]; then
  mapfile -t YAML_CANDIDATES < <(find "$YAML_CONFIG" -maxdepth 1 -name "*.yaml" -type f | sort)
  if [[ ${#YAML_CANDIDATES[@]} -eq 0 ]]; then
    echo "ERROR: No YAML file found in model directory: $YAML_CONFIG" >&2
    exit 1
  fi
  if [[ ${#YAML_CANDIDATES[@]} -gt 1 ]]; then
    echo "Warning: Multiple YAML files found in $YAML_CONFIG, using first: ${YAML_CANDIDATES[0]}" >&2
  fi
  YAML_CONFIG="${YAML_CANDIDATES[0]}"
fi

if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "ERROR: Python executable not found: $PYTHON_BIN" >&2
  exit 1
fi

if ! command -v "$FC" >/dev/null 2>&1; then
  echo "ERROR: Fortran compiler not found: $FC" >&2
  exit 1
fi

if [[ ! -f "$YAML_CONFIG" ]]; then
  echo "ERROR: YAML config file not found: $YAML_CONFIG" >&2
  exit 1
fi

YAML_CONFIG="$(cd "$(dirname "$YAML_CONFIG")" && pwd)/$(basename "$YAML_CONFIG")"
OUTPUT_ROOT="$(mkdir -p "$OUTPUT_ROOT" && cd "$OUTPUT_ROOT" && pwd)"

if [[ -z "$MODEL_NAME" ]]; then
  MODEL_NAME="$(basename "$YAML_CONFIG" .yaml | tr '[:upper:]' '[:lower:]' | tr -c '[:alnum:]_' '_')"
fi

BUILD_DIR="$OUTPUT_ROOT/$MODEL_NAME"
GENERATED_DIR="$BUILD_DIR/generated"
STAGE_DIR="$BUILD_DIR/library_src"
LIB_DIR="$BUILD_DIR/lib"
LOG_DIR="$BUILD_DIR/logs"
METADATA_JSON="$BUILD_DIR/brns_library_metadata.json"
LIB_BASENAME="libbrns_reactions_${MODEL_NAME}.so"
LIB_PATH="$LIB_DIR/$LIB_BASENAME"

mkdir -p "$GENERATED_DIR" "$STAGE_DIR" "$LIB_DIR" "$LOG_DIR"
rm -f "$LIB_PATH"

export PYTHONPATH="$ROOT_DIR${PYTHONPATH:+:$PYTHONPATH}"

echo "=========================================="
echo "BRNS Reaction Library Build"
echo "=========================================="
echo "YAML:       $YAML_CONFIG"
echo "Model:      $MODEL_NAME"
echo "Python:     $PYTHON_BIN"
echo "Compiler:   $FC"
echo "Build dir:  $BUILD_DIR"
echo ""

echo "[1/5] Generating Fortran from YAML..."
"$PYTHON_BIN" - "$YAML_CONFIG" "$GENERATED_DIR" "$METADATA_JSON" "$MODEL_NAME" "$LIB_BASENAME" <<'PYEOF'
import json
import sys
from pathlib import Path
import yaml

from acg_brns.acg_orchestrator import ACGOrchestrator

yaml_path = Path(sys.argv[1]).resolve()
generated_dir = Path(sys.argv[2]).resolve()
metadata_json = Path(sys.argv[3]).resolve()
model_name = sys.argv[4]
lib_basename = sys.argv[5]

orchestrator = ACGOrchestrator(str(yaml_path), str(generated_dir), verbose=True)
summary = orchestrator.generate()

config = yaml.safe_load(yaml_path.read_text(encoding="utf-8"))
species_names = list(orchestrator.acg_data.get("variables", []))
parameter_names = list(orchestrator.acg_data.get("bio_name", []))
default_parameter_values = []
for value in orchestrator.acg_data.get("bio_val", []):
    try:
        default_parameter_values.append(float(value))
    except Exception:
        default_parameter_values.append(float(str(value)))

metadata = {
    "model_name": model_name,
    "yaml_path": str(yaml_path),
    "generated_dir": str(generated_dir),
    "library_basename": lib_basename,
    "n_species": int(summary["n_species"]),
    "n_dissolved": int(summary["n_dissolved"]),
    "n_solid": int(summary["n_solid"]),
    "n_reactions": int(summary["n_reactions"]),
    "n_parameters": int(summary["n_parameters"]),
    "species_names": species_names,
    "parameter_names": parameter_names,
    "default_parameter_values": default_parameter_values,
    "reaction_names": [rxn.get("name", f"reaction_{i+1}") for i, rxn in enumerate(config.get("reactions", []))],
    "files_generated": summary["files_generated"],
}
metadata_json.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
print(f"Metadata written: {metadata_json}")
PYEOF

N_SPECIES=$("$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1], encoding="utf-8"))["n_species"])' "$METADATA_JSON")
N_REACTIONS=$("$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1], encoding="utf-8"))["n_reactions"])' "$METADATA_JSON")
N_PARAMETERS=$("$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1], encoding="utf-8"))["n_parameters"])' "$METADATA_JSON")

echo "[2/5] Staging library sources..."
COMMON_FILES=(
  common.inc common_drive.inc common_geo.inc defines.inc
  basic.f biogeo.f boundaries.f drivervalues.f gaussj.f jacobian.f limits.f rates.f residual.f switches.f
  LUBKSB.F LUDCMP.F MPROVE.F newtonsub.f
)

for src in "${COMMON_FILES[@]}"; do
  cp "$FORTRAN_COMMON/$src" "$STAGE_DIR/$src"
done

cp "$GENERATED_DIR"/*.f "$STAGE_DIR/" 2>/dev/null || true
cp "$GENERATED_DIR"/*.F "$STAGE_DIR/" 2>/dev/null || true
cp "$GENERATED_DIR"/*.inc "$STAGE_DIR/" 2>/dev/null || true

sed \
  -e "s|include '../common_geo.inc'|include 'common_geo.inc'|g" \
  -e "s|include '../common.inc'|include 'common.inc'|g" \
  -e "s|include '../common_drive.inc'|include 'common_drive.inc'|g" \
  "$INVOKEBRNS_SOURCE" > "$STAGE_DIR/invokebrns.f"

"$PYTHON_BIN" - "$STAGE_DIR/invokebrns.f" <<'PYEOF'
from pathlib import Path
import sys

path = Path(sys.argv[1])
text = path.read_text(encoding="latin-1")

replacements = {
    "      real*8 concAfterTransport(1)": "      real*8 concAfterTransport(*)",
    "      real*8 concBeforeTransport(1)": "      real*8 concBeforeTransport(*)",
    "      real*8 outputConcentrations(1)": "      real*8 outputConcentrations(*)",
    "      integer fixedConcentrationBoundary(1)": "      integer fixedConcentrationBoundary(*)",
    "      real*8 parameterVector(1)": "      real*8 parameterVector(*)",
    "      open(unit=11,file='ratesAtFinish.dat',access='append',\n     + dispose='DELETE')\n      close(11)": "      open(unit=11,file='ratesAtFinish.dat',status='replace')\n      close(11)",
}

for old, new in replacements.items():
    text = text.replace(old, new)

path.write_text(text, encoding="latin-1")
print(f"Patched for Linux compatibility: {path}")
PYEOF

cat > "$STAGE_DIR/brns_reaction_api.f90" <<EOF
module brns_reaction_api
  use iso_c_binding
  implicit none

  integer(c_int), parameter :: BRNS_N_SPECIES = ${N_SPECIES}
  integer(c_int), parameter :: BRNS_N_REACTIONS = ${N_REACTIONS}
  integer(c_int), parameter :: BRNS_N_PARAMETERS = ${N_PARAMETERS}

contains

  subroutine brns_get_n_species(n_species) bind(C, name="brns_get_n_species")
    integer(c_int), intent(out) :: n_species
    n_species = BRNS_N_SPECIES
  end subroutine brns_get_n_species

  subroutine brns_get_n_reactions(n_reactions) bind(C, name="brns_get_n_reactions")
    integer(c_int), intent(out) :: n_reactions
    n_reactions = BRNS_N_REACTIONS
  end subroutine brns_get_n_reactions

  subroutine brns_get_n_parameters(n_parameters) bind(C, name="brns_get_n_parameters")
    integer(c_int), intent(out) :: n_parameters
    n_parameters = BRNS_N_PARAMETERS
  end subroutine brns_get_n_parameters

  subroutine brns_react_cell(conc_after_transport, conc_before_transport, output_concentrations, &
                             number_of_species, time_step, fixed_boundary, return_code, &
                             pos_x, pos_y, pos_z, porosity, water_saturation, parameter_vector) &
                             bind(C, name="brns_react_cell")
    use iso_c_binding
    implicit none

    real(c_double), intent(in) :: conc_after_transport(*)
    real(c_double), intent(in) :: conc_before_transport(*)
    real(c_double), intent(out) :: output_concentrations(*)
    integer(c_int), value, intent(in) :: number_of_species
    real(c_double), value, intent(in) :: time_step
    integer(c_int), intent(in) :: fixed_boundary(*)
    integer(c_int), intent(out) :: return_code
    real(c_double), value, intent(in) :: pos_x, pos_y, pos_z
    real(c_double), value, intent(in) :: porosity, water_saturation
    real(c_double), intent(in) :: parameter_vector(*)

    integer :: n_species_f
    integer :: return_code_f
    real(c_double) :: dt_f
    real(c_double) :: pos_x_f, pos_y_f, pos_z_f
    real(c_double) :: porosity_f, water_saturation_f
    external :: invokebrns

    n_species_f = number_of_species
    dt_f = time_step
    pos_x_f = pos_x
    pos_y_f = pos_y
    pos_z_f = pos_z
    porosity_f = porosity
    water_saturation_f = water_saturation
    return_code_f = 0

    call invokebrns(conc_after_transport, conc_before_transport, output_concentrations, &
                    n_species_f, dt_f, fixed_boundary, return_code_f, &
                    pos_x_f, pos_y_f, pos_z_f, porosity_f, water_saturation_f, parameter_vector)

    return_code = return_code_f
  end subroutine brns_react_cell

end module brns_reaction_api
EOF

echo "[3/5] Compiling shared library..."
SOURCES=(
  basic.f biogeo.f boundaries.f drivervalues.f gaussj.f jacobian.f limits.f rates.f residual.f switches.f
  LUBKSB.F LUDCMP.F MPROVE.F newtonsub.f
  parameters.f invokebrns.f brns_reaction_api.f90
)

SOURCE_PATHS=()
for src in "${SOURCES[@]}"; do
  SOURCE_PATHS+=("$STAGE_DIR/$src")
done

FFLAGS=(
  -shared -fPIC -cpp -O2 -g -I"$STAGE_DIR"
  -ffixed-line-length-none -std=legacy -fallow-argument-mismatch
)
LDFLAGS=(-llapack -lblas)

set -o pipefail
"$FC" "${FFLAGS[@]}" "${SOURCE_PATHS[@]}" -o "$LIB_PATH" "${LDFLAGS[@]}" 2>&1 | tee "$LOG_DIR/compile.log"
set +o pipefail

echo "[4/5] Finalizing metadata..."
"$PYTHON_BIN" - "$METADATA_JSON" "$LIB_PATH" <<'PYEOF'
import json
import sys
from pathlib import Path
metadata_path = Path(sys.argv[1])
lib_path = Path(sys.argv[2]).resolve()
metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
metadata["library_path"] = str(lib_path)
metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
print(f"Updated metadata: {metadata_path}")
PYEOF

echo "[5/5] Done"
echo "Library:   $LIB_PATH"
echo "Metadata:  $METADATA_JSON"
echo "Compile log: $LOG_DIR/compile.log"
