#!/bin/bash
#
# Simple build and run script for Python-generated BRNS code
# User-configurable: YAML config path, optional input files path, output directory
#
# Usage: ./build_python.sh -c YAML_CONFIG_OR_MODEL_DIR [-i INPUT_DIR] [-o OUTPUT_DIR] [-n NAME]
#   -c YAML_CONFIG_OR_MODEL_DIR  Path to YAML file, model directory, or model name under ./models (required)
#   -i INPUT_DIR                 Directory containing input files (.inp) (default: YAML/model directory)
#   -o OUTPUT_DIR     Output directory for generated code (default: ./build_output)
#   -n NAME           Output name/identifier (default: derived from YAML filename)
#   -h                Show this help message
#
# Examples:
#   ./build_python.sh -c ./models/single_species/single_species_example.yaml
#   ./build_python.sh -c ./models/multiple_species
#   ./build_python.sh -c equilibrium -o ./my_build -n eq_v2
#

set -e

# ==========================================
# Colors
# ==========================================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

# ==========================================
# Default parameters
# ==========================================

YAML_CONFIG=""
INPUT_DIR=""
OUTPUT_DIR="./build_output"
MODEL_NAME=""
FORTRAN_COMMON="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/BRNSPackage/FortranFiles"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FFLAGS="-cpp -I. -O2 -finit-local-zero -fno-unsafe-math-optimizations -fPIC -g -Wall -Wextra -Wuninitialized"
LIBS="-llapack -lblas"

# ==========================================
# Parse command line arguments
# ==========================================

show_help() {
    sed -n '2,20p' "$0" | sed 's/^# //g'
    exit 0
}

if [ $# -eq 0 ]; then
    show_help
fi

while getopts "c:i:o:n:h" opt; do
    case $opt in
        c) YAML_CONFIG="$OPTARG" ;;
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        n) MODEL_NAME="$OPTARG" ;;
        h) show_help ;;
        *)
            echo -e "${RED}ERROR: Unknown option -$OPTARG${NC}"
            echo "Use -h for help"
            exit 1
            ;;
    esac
done

# ==========================================
# Validate required parameters
# ==========================================

if [ -z "$YAML_CONFIG" ]; then
    echo -e "${RED}ERROR: YAML config file (-c) is required!${NC}"
    echo "Use -h for help"
    exit 1
fi

# Allow shorthand model name from ./models
if [ ! -e "$YAML_CONFIG" ] && [ -d "$SCRIPT_DIR/models/$YAML_CONFIG" ]; then
    YAML_CONFIG="$SCRIPT_DIR/models/$YAML_CONFIG"
fi

# If -c points to a model directory, pick YAML from that directory
if [ -d "$YAML_CONFIG" ]; then
    mapfile -t YAML_CANDIDATES < <(find "$YAML_CONFIG" -maxdepth 1 -name "*.yaml" -type f | sort)
    if [ ${#YAML_CANDIDATES[@]} -eq 0 ]; then
        echo -e "${RED}ERROR: No YAML file found in model directory: $YAML_CONFIG${NC}"
        exit 1
    fi
    if [ ${#YAML_CANDIDATES[@]} -gt 1 ]; then
        echo -e "${YELLOW}Warning: Multiple YAML files found in $YAML_CONFIG, using first: ${YAML_CANDIDATES[0]}${NC}"
    fi
    YAML_CONFIG="${YAML_CANDIDATES[0]}"
fi

if [ ! -f "$YAML_CONFIG" ]; then
    echo -e "${RED}ERROR: YAML config file not found: $YAML_CONFIG${NC}"
    exit 1
fi

if [ -z "$INPUT_DIR" ]; then
    INPUT_DIR="$(dirname "$YAML_CONFIG")"
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo -e "${RED}ERROR: Input directory not found: $INPUT_DIR${NC}"
    exit 1
fi

# ==========================================
# Set up paths and names
# ==========================================

YAML_CONFIG="$(cd "$(dirname "$YAML_CONFIG")" && pwd)/$(basename "$YAML_CONFIG")"
INPUT_DIR="$(cd "$(dirname "$INPUT_DIR")" && pwd)/$(basename "$INPUT_DIR")"
OUTPUT_DIR="$(mkdir -p "$OUTPUT_DIR" && cd "$OUTPUT_DIR" && pwd)"

if [ -z "$MODEL_NAME" ]; then
    MODEL_NAME=$(basename "$YAML_CONFIG" .yaml | tr '[:upper:]' '[:lower:]')
fi

BUILD_DIR="$OUTPUT_DIR/$MODEL_NAME"
BUILD_PY="$BUILD_DIR/python"
RESULTS_DIR="$BUILD_DIR/results"
PYTHON_GEN="$BUILD_DIR/generated"

# ==========================================
# Find input files
# ==========================================

find_input_files() {
    local input_files=()
    if [ -d "$INPUT_DIR" ]; then
        while IFS= read -r -d '' file; do
            input_files+=("$file")
        done < <(find "$INPUT_DIR" -maxdepth 1 -name "*.inp" -print0 | sort -z)
    fi
    printf '%s\n' "${input_files[@]}"
}

INPUT_FILES=($(find_input_files))

# ==========================================
# Display configuration
# ==========================================

echo "=========================================="
echo "BRNS Build & Run: Python Code Generator"
echo "=========================================="
echo ""
echo -e "${CYAN}YAML Config:${NC}    $YAML_CONFIG"
echo -e "${CYAN}Model Name:${NC}     $MODEL_NAME"
echo -e "${CYAN}Input Directory:${NC} $INPUT_DIR"
echo -e "${CYAN}Input Files:${NC}    ${#INPUT_FILES[@]} found"
echo -e "${CYAN}Build Directory:${NC} $BUILD_PY"
echo -e "${CYAN}Output Directory:${NC} $RESULTS_DIR"
echo ""

# ==========================================
# Sanity check: input files
# ==========================================

if [ ${#INPUT_FILES[@]} -eq 0 ]; then
    echo -e "${YELLOW}Warning: No input files (.inp) found in $INPUT_DIR${NC}"
    echo ""
fi

# ==========================================
# Check dependencies
# ==========================================

if ! command -v gfortran &>/dev/null; then
    echo -e "${RED}ERROR: gfortran not found!${NC}"
    echo "Please install gfortran (e.g., apt-get install gfortran)"
    exit 1
fi

if ! command -v python3 &>/dev/null; then
    echo -e "${RED}ERROR: python3 not found!${NC}"
    echo "Please install python3"
    exit 1
fi

# ==========================================
# Step 1: Generate Fortran code from YAML
# ==========================================

echo "=========================================="
echo "Step 1: Generate Fortran Code"
echo "=========================================="
echo ""

mkdir -p "$PYTHON_GEN"
cd "$PYTHON_GEN"

echo "Running Python code generator..."
PYTHON_YAML_CONFIG="$YAML_CONFIG"
PYTHON_YAML_DIR="$(dirname "$YAML_CONFIG")"
PYTHON_SCRIPT_DIR="$SCRIPT_DIR"

if command -v cygpath &>/dev/null; then
    PYTHON_YAML_CONFIG="$(cygpath -w "$YAML_CONFIG")"
    PYTHON_YAML_DIR="$(cygpath -w "$PYTHON_YAML_DIR")"
    PYTHON_SCRIPT_DIR="$(cygpath -w "$PYTHON_SCRIPT_DIR")"
fi

python3 -c "
import sys
yaml_config, yaml_dir, script_dir = sys.argv[1:4]
sys.path.insert(0, yaml_dir)
sys.path.insert(0, script_dir)

from acg_brns.acg_orchestrator import ACGOrchestrator

try:
    orchestrator = ACGOrchestrator(yaml_config, '.')
    orchestrator.load_config()
    print(f'✓ Loaded YAML: {orchestrator.config.get(\"model_name\", \"unknown\")}')
    
    orchestrator.evaluate_formulas()
    print('✓ Evaluated formulas')
    
    orchestrator.map_to_acg_structures()
    print('✓ Mapped to ACG structures')
    
    orchestrator.run_preprocessing()
    print('✓ Preprocessing completed')
    
    orchestrator.run_code_generation()
    print('✓ Fortran code generated')
    
except Exception as e:
    print(f'✗ Error: {e}', file=sys.stderr)
    import traceback
    traceback.print_exc()
    sys.exit(1)
" "$PYTHON_YAML_CONFIG" "$PYTHON_YAML_DIR" "$PYTHON_SCRIPT_DIR" || exit 1

echo ""

# ==========================================
# Step 2: Collect files for build
# ==========================================

echo "=========================================="
echo "Step 2: Prepare Build Environment"
echo "=========================================="
echo ""

mkdir -p "$BUILD_PY"
cd "$BUILD_PY"

echo "Copying framework files..."
cp "$FORTRAN_COMMON"/*.f . 2>/dev/null || true
cp "$FORTRAN_COMMON"/*.F . 2>/dev/null || true
cp "$FORTRAN_COMMON"/*.inc . 2>/dev/null || true

echo "Copying generated files..."
cp "$PYTHON_GEN"/*.f . 2>/dev/null || true
cp "$PYTHON_GEN"/*.F . 2>/dev/null || true
cp "$PYTHON_GEN"/*.inc . 2>/dev/null || true

echo "Fixing include paths (relative paths → local paths)..."
for f in *.f *.F; do
    if [ -f "$f" ]; then
        sed -i "s|include[[:space:]]*'../\([^']*\)'|include '\1'|gi" "$f" 2>/dev/null || true
        sed -i 's|include[[:space:]]*"../\([^"]*\)"|include "\1"|gi' "$f" 2>/dev/null || true
    fi
done 2>/dev/null || true

echo "Copying input files (${#INPUT_FILES[@]})..."
for inp_path in "${INPUT_FILES[@]}"; do
    inp=$(basename "$inp_path")
    cp "$inp_path" .
    echo "  ✓ $inp"
done

rm -f printsvnversion_nosvn.f printsvnversion_tmpl.f 2>/dev/null || true

# Compile only the forward-model source set (no optimization files)
CORE_SOURCES=(
    main.f
    basic.f biogeo.f boundaries.f drivervalues.f diagenesis.f
    advdiffcoeff.f gridsetup.f porarea.f molecular.f initialcond.f
    issolid.f jacobian.f limits.f rates.f residual.f ssrates.f steadystate.f switches.f output.f
    notransport.f getdelt.f timestep.f transport.f transcoeff.f transcoeff-MT.f
    gaussj.f LUBKSB.F LUDCMP.F MPROVE.F NEWT.F newtonsub.f TRIDAG.F
    parameters.f printdepth.f printsvnversion.f
)

COMPILE_SOURCES=()
for src in "${CORE_SOURCES[@]}"; do
    if [ -f "$src" ]; then
        COMPILE_SOURCES+=("$src")
    fi
done

if [ ${#COMPILE_SOURCES[@]} -eq 0 ]; then
    echo -e "${RED}ERROR: No Fortran sources selected for compilation!${NC}"
    exit 1
fi

echo ""

# ==========================================
# Step 3: Compile
# ==========================================

echo "=========================================="
echo "Step 3: Compile"
echo "=========================================="
echo ""

EXEC="brns_python"

echo "Compiling $EXEC..."
echo "FFLAGS: $FFLAGS"
echo "Source files: ${#COMPILE_SOURCES[@]}"
echo ""

gfortran $FFLAGS "${COMPILE_SOURCES[@]}" -o "$EXEC" $LIBS 2>&1 | tee compile.log
COMPILE_STATUS=${PIPESTATUS[0]}

if [ $COMPILE_STATUS -eq 0 ] && [ -f "$EXEC" ]; then
    echo ""
    echo -e "${GREEN}✓ Compilation successful${NC}"
    ls -lh "$EXEC"
    echo ""
else
    echo ""
    echo -e "${RED}✗ Compilation failed!${NC}"
    if [ ! -f "$EXEC" ]; then
        echo "  Executable $EXEC was not created"
    fi
    echo ""
    echo "Last 30 lines from compile.log:"
    tail -30 compile.log
    exit 1
fi

# ==========================================
# Step 4: Run
# ==========================================

echo "=========================================="
echo "Step 4: Run"
echo "=========================================="
echo ""

mkdir -p "$RESULTS_DIR"

echo "Starting: ./$EXEC"
echo "Threads: OMP=1 OPENBLAS=1 MKL=1"
echo ""

env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    ./"$EXEC" 2>&1 | tee run.log
RUN_STATUS=${PIPESTATUS[0]}

if [ $RUN_STATUS -eq 0 ]; then
    echo ""
    echo -e "${GREEN}✓ Execution successful${NC}"
    
    cp *.dat "$RESULTS_DIR/" 2>/dev/null || true
    cp *.inp "$RESULTS_DIR/" 2>/dev/null || true
    
    echo ""
    echo "Output files:"
    ls -lh "$RESULTS_DIR"/*.dat 2>/dev/null | awk '{print "  " $9}' || echo "  None found"
    echo ""
    exit 0
else
    echo ""
    echo -e "${RED}✗ Execution failed!${NC}"
    echo "Check run.log for details"
    echo ""
    exit 1
fi
