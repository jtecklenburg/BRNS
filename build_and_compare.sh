#!/bin/bash
#
# Flexible build script for BRNS examples
# Automatically detects available examples from directory structure
#
# Usage: ./build_python.sh [EXAMPLE_NAME] [ACTION]
#   EXAMPLE_NAME: auto-detected or specific (e.g. single_species, multiple_species)
#   ACTION: all (default), build, run, compare
#

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

# ==========================================
# Parameters and configuration
# ==========================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXAMPLE="${1:-}"
ACTION="${2:-all}"

FORTRAN_COMMON="$SCRIPT_DIR/BRNSPackage/FortranFiles"
FFLAGS="-cpp -I. -O2 -finit-local-zero -fno-unsafe-math-optimizations -fPIC -g"
LIBS="-llapack -lblas"
RUN_TIMEOUT_SECONDS="${RUN_TIMEOUT_SECONDS:-0}"

# ==========================================
# Auto-detect available examples
# ==========================================

detect_examples() {
    local examples=()
    
    # Search for generated Python directories
    for dir in "$SCRIPT_DIR/generated_fortran"/*; do
        if [ -d "$dir" ]; then
            local name=$(basename "$dir")
            # Check if at least one .f or .inc file exists
            if ls "$dir"/*.f "$dir"/*.inc >/dev/null 2>&1; then
                examples+=("$name")
            fi
        fi
    done
    
    printf '%s\n' "${examples[@]}"
}

# Determine available examples
AVAILABLE_EXAMPLES=($(detect_examples))

# If no example specified, show list
if [ -z "$EXAMPLE" ]; then
    echo "=========================================="
    echo "BRNS Build System"
    echo "=========================================="
    echo ""
    echo "Available examples:"
    echo ""
    for ex in "${AVAILABLE_EXAMPLES[@]}"; do
        echo "  - $ex"
    done
    echo ""
    echo "Usage: ./build_python.sh [EXAMPLE] [ACTION]"
    echo "  EXAMPLE: $(IFS=', '; echo "${AVAILABLE_EXAMPLES[*]}")"
    echo "  ACTION:  all (default), build, run, compare"
    echo ""
    echo "Examples:"
    echo "  ./build_python.sh single_species"
    echo "  ./build_python.sh multiple_species build"
    echo "  ./build_python.sh single_species run"
    echo ""
    exit 0
fi

# Check if example exists
EXAMPLE_EXISTS=false
for ex in "${AVAILABLE_EXAMPLES[@]}"; do
    if [ "$ex" = "$EXAMPLE" ]; then
        EXAMPLE_EXISTS=true
        break
    fi
done

if [ "$EXAMPLE_EXISTS" = false ]; then
    echo -e "${RED}ERROR: Example '$EXAMPLE' not found!${NC}"
    echo "Available examples: ${AVAILABLE_EXAMPLES[*]}"
    exit 1
fi

# ==========================================
# Directories for selected example
# ==========================================

REFERENCE_GEN="$SCRIPT_DIR/reference_fortran/$EXAMPLE"
PYTHON_GEN="$SCRIPT_DIR/generated_fortran/$EXAMPLE"
MODEL_DIR="$SCRIPT_DIR/models/$EXAMPLE"
BUILD_BASE="$SCRIPT_DIR/build/$EXAMPLE"
BUILD_REF="$BUILD_BASE/reference"
BUILD_PY="$BUILD_BASE/python"
RESULTS_DIR="$BUILD_BASE/results"

# ==========================================
# Auto-detect input files
# ==========================================

detect_input_files() {
    local search_dirs=(
        "$MODEL_DIR"
        "$REFERENCE_GEN"
        "$PYTHON_GEN"
    )
    
    local found_files=()
    
    for dir in "${search_dirs[@]}"; do
        if [ -d "$dir" ]; then
            for inp in "$dir"/*.inp; do
                if [ -f "$inp" ]; then
                    local basename=$(basename "$inp")
                    # Add only if not already in list
                    if [[ ! " ${found_files[@]} " =~ " ${basename} " ]]; then
                        found_files+=("$basename")
                    fi
                fi
            done
        fi
    done
    
    printf '%s\n' "${found_files[@]}"
}

INPUT_FILES=($(detect_input_files))

# ==========================================
# Header
# ==========================================

echo "=========================================="
echo "BRNS Build & Run: $EXAMPLE"
echo "=========================================="
echo ""
echo -e "${CYAN}Example:${NC}        $EXAMPLE"
echo -e "${CYAN}Model dir:${NC}      $MODEL_DIR"
echo -e "${CYAN}Python-Gen:${NC}     $PYTHON_GEN"
echo -e "${CYAN}Reference:${NC}      $REFERENCE_GEN"
echo -e "${CYAN}Input files:${NC}    ${#INPUT_FILES[@]} found"
echo -e "${CYAN}Action:${NC}         $ACTION"
if [ "$RUN_TIMEOUT_SECONDS" -gt 0 ]; then
    echo -e "${CYAN}Timeout:${NC}        ${RUN_TIMEOUT_SECONDS}s"
else
    echo -e "${CYAN}Timeout:${NC}        disabled"
fi
echo ""

# ==========================================
# Build function
# ==========================================

build_version() {
    local VERSION=$1
    local SOURCE=$2
    local TARGET=$3
    local EXEC=$4
    
    echo "=========================================="
    echo "Build: $VERSION ($EXAMPLE)"
    echo "=========================================="
    
    rm -rf "$TARGET"
    mkdir -p "$TARGET"
    cd "$TARGET"
    
    echo "Copying framework..."
    cp "$FORTRAN_COMMON"/*.f . 2>/dev/null || true
    cp "$FORTRAN_COMMON"/*.F . 2>/dev/null || true
    cp "$FORTRAN_COMMON"/*.inc . 2>/dev/null || true

    # Model-specific overrides from reference for Python build only
    # (keeps reference branch unchanged and avoids pulling optimization files into Python compile set)
    if [ "$VERSION" = "Python" ] && [ -d "$REFERENCE_GEN" ]; then
        echo "Copying selected model-specific overrides from reference..."
        REFERENCE_OVERRIDES=(gridsetup.f advdiffcoeff.f porarea.f transcoeff.f transcoeff-MT.f)
        for ref_src in "${REFERENCE_OVERRIDES[@]}"; do
            if [ -f "$REFERENCE_GEN/$ref_src" ]; then
                cp "$REFERENCE_GEN/$ref_src" .
                echo "  ✓ override: $ref_src"
            fi
        done
    fi
    
    echo "Copying generated files..."
    cp "$SOURCE"/*.f . 2>/dev/null || true
    cp "$SOURCE"/*.F . 2>/dev/null || true
    cp "$SOURCE"/*.inc . 2>/dev/null || true
    
    echo "Fixing include paths (relative paths → local paths)..."
    # Replace '../common_geo.inc' with 'common_geo.inc' etc.
    for f in *.f *.F; do
        if [ -f "$f" ]; then
            sed -i "s|include[[:space:]]*'../\([^']*\)'|include '\1'|gi" "$f" 2>/dev/null || true
            sed -i 's|include[[:space:]]*"../\([^"]*\)"|include "\1"|gi' "$f" 2>/dev/null || true
        fi
    done 2>/dev/null || true
    
    echo "Copying input files (${#INPUT_FILES[@]})..."
    for inp in "${INPUT_FILES[@]}"; do
        if [ -f "$MODEL_DIR/$inp" ]; then
            cp "$MODEL_DIR/$inp" .
            echo "  ✓ $inp: $MODEL_DIR/$inp -> $TARGET/$inp"
        else
            # Fallbacks if not present in models
            if [ -f "$REFERENCE_GEN/$inp" ]; then
                cp "$REFERENCE_GEN/$inp" .
                echo "  ⚠ $inp (fallback): $REFERENCE_GEN/$inp -> $TARGET/$inp"
            elif [ -f "$PYTHON_GEN/$inp" ]; then
                cp "$PYTHON_GEN/$inp" .
                echo "  ⚠ $inp (fallback): $PYTHON_GEN/$inp -> $TARGET/$inp"
            else
                echo "  ✗ $inp not found"
            fi
        fi
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
        return 1
    fi
    
    echo ""
    echo "Compiling $EXEC..."
    echo "FFLAGS: $FFLAGS"
    echo "Source files: ${#COMPILE_SOURCES[@]}"
    echo ""
    
    # Compile and save exit code correctly (not overridden by tee)
    gfortran $FFLAGS "${COMPILE_SOURCES[@]}" -o "$EXEC" $LIBS 2>&1 | tee compile.log
    local COMPILE_STATUS=${PIPESTATUS[0]}
    
    # Check both exit code and if executable exists
    if [ $COMPILE_STATUS -eq 0 ] && [ -f "$EXEC" ]; then
        echo ""
        echo -e "${GREEN}✓ $VERSION compiled successfully${NC}"
        ls -lh "$EXEC"
        echo ""
        return 0
    else
        echo ""
        echo -e "${RED}✗ Compilation failed!${NC}"
        if [ ! -f "$EXEC" ]; then
            echo "  Executable $EXEC was not created"
        fi
        echo ""
        echo "Last 30 lines from compile.log:"
        tail -30 compile.log
        return 1
    fi
}

# ==========================================
# Run function
# ==========================================

run_version() {
    local VERSION=$1
    local BUILDDIR=$2
    local EXEC=$3
    local OUTDIR=$4
    local RUN_STATUS=0
    
    echo "=========================================="
    echo "Run: $VERSION ($EXAMPLE)"
    echo "=========================================="
    
    cd "$BUILDDIR"
    
    echo "Copying fresh input files..."
    for inp in "${INPUT_FILES[@]}"; do
        if [ -f "$MODEL_DIR/$inp" ]; then
            cp "$MODEL_DIR/$inp" .
            echo "  ✓ $inp: $MODEL_DIR/$inp -> $BUILDDIR/$inp"
        else
            # Fallbacks if not present in models
            if [ -f "$REFERENCE_GEN/$inp" ]; then
                cp "$REFERENCE_GEN/$inp" .
                echo "  ⚠ $inp (fallback): $REFERENCE_GEN/$inp -> $BUILDDIR/$inp"
            else
                echo "  ✗ $inp not found"
            fi
        fi
    done
    echo ""
    
    echo "Starting: ./$EXEC"
    echo "Threads: OMP=1 OPENBLAS=1 MKL=1"
    echo ""

    if [ "$RUN_TIMEOUT_SECONDS" -gt 0 ]; then
        env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
            timeout "${RUN_TIMEOUT_SECONDS}s" ./"$EXEC" 2>&1 | tee run.log
        RUN_STATUS=${PIPESTATUS[0]}
    else
        env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
            ./"$EXEC" 2>&1 | tee run.log
        RUN_STATUS=${PIPESTATUS[0]}
    fi

    if [ $RUN_STATUS -eq 0 ]; then
        echo ""
        echo -e "${GREEN}✓ $VERSION executed successfully${NC}"
        
        mkdir -p "$OUTDIR"
        cp *.dat "$OUTDIR/" 2>/dev/null || true
        cp *.inp "$OUTDIR/" 2>/dev/null || true
        
        echo ""
        echo "Output files:"
        ls -lh *.dat *.inp 2>/dev/null | head -20 || echo "  None found"
        echo ""
        return 0
    else
        echo ""
        if [ $RUN_STATUS -eq 124 ]; then
            echo -e "${RED}✗ Execution interrupted due to timeout (${RUN_TIMEOUT_SECONDS}s)!${NC}"
            echo "  Hint: Partial .dat files may end differently."
            echo "  Recommendation: RUN_TIMEOUT_SECONDS=0 ./build_python.sh $EXAMPLE run"
        else
            echo -e "${RED}✗ Execution failed!${NC}"
        fi
        return 1
    fi
}

# ==========================================
# Compare function
# ==========================================

compare_results() {
    # Ensure we don't stop on first difference inside comparisons
    local _old_errexit=$-
    set +e
    echo "=========================================="
    echo "Comparing results ($EXAMPLE)"
    echo "=========================================="
    echo ""
    
    cd "$RESULTS_DIR"
    
    local REF_FILES=$(ls reference/*.dat reference/*.inp 2>/dev/null | xargs -n1 basename | sort || echo "")
    local PY_FILES=$(ls python/*.dat python/*.inp 2>/dev/null | xargs -n1 basename | sort || echo "")

    if [ -z "$REF_FILES" ] && [ -z "$PY_FILES" ]; then
        echo -e "${YELLOW}No output files to compare!${NC}"
        return 0
    fi

    # Union of filenames from reference and python outputs
    local ALL_FILES=$(printf "%s\n%s\n" "$REF_FILES" "$PY_FILES" | sort -u)


    local ALL_IDENTICAL=true
    local IDENTICAL=0
    local DIFFERENT=0
    local MISSING=0
    local EXTRA=0

    for file in $ALL_FILES; do
        local REF="reference/$file"
        local PY="python/$file"

        if [ ! -f "$REF" ] && [ -f "$PY" ]; then
            echo -e "${YELLOW}⊕ $file${NC} - Only in Python"
            ALL_IDENTICAL=false
            ((EXTRA++))
            continue
        fi

        if [ -f "$REF" ] && [ ! -f "$PY" ]; then
            echo -e "${RED}✗ $file${NC} - Missing in Python"
            ALL_IDENTICAL=false
            ((MISSING++))
            continue
        fi

        if diff -q "$REF" "$PY" >/dev/null 2>&1; then
            echo -e "${GREEN}✓ $file${NC} - identical"
            ((IDENTICAL++))
        else
            echo -e "${YELLOW}≈ $file${NC} - differences"
            ALL_IDENTICAL=false
            ((DIFFERENT++))

            # Numeric comparison (max abs only)
            if command -v python3 &>/dev/null; then
                python3 <<EOF 2>/dev/null || true
import numpy as np
def load_data(path):
    try:
        return np.loadtxt(path)
    except Exception:
        try:
            return np.loadtxt(path, skiprows=1)
        except Exception:
            return None

ref = load_data('$REF')
py = load_data('$PY')
if ref is not None and py is not None and ref.shape == py.shape:
    diff = np.abs(ref - py)
    rel = diff / (np.abs(ref) + 1e-15)
    print(f"    Max abs: {np.max(diff):.2e}, rel: {np.max(rel):.2e}")
EOF
            fi
        fi
    done
    
    echo ""
    echo "Summary: ✓ $IDENTICAL | ≈ $DIFFERENT | ✗ $MISSING | ⊕ $EXTRA"
    echo ""
    
    if $ALL_IDENTICAL; then
        echo -e "${GREEN}✓ All results identical!${NC}"
        # restore errexit
        [[ $_old_errexit == *e* ]] && set -e
        return 0
    else
        echo -e "${YELLOW}≈ Differences found${NC}"
        # restore errexit
        [[ $_old_errexit == *e* ]] && set -e
        return 1
    fi
}

# ==========================================
# Main program
# ==========================================

if ! command -v gfortran &>/dev/null; then
    echo -e "${RED}ERROR: gfortran not found!${NC}"
    exit 1
fi

SKIP_REF=false
if [ ! -d "$REFERENCE_GEN" ]; then
    echo -e "${YELLOW}Note: No reference available for '$EXAMPLE'${NC}"
    SKIP_REF=true
    echo ""
fi

case "$ACTION" in
    all)
        [ "$SKIP_REF" = false ] && build_version "Reference" "$REFERENCE_GEN" "$BUILD_REF" "brns_reference"
        build_version "Python" "$PYTHON_GEN" "$BUILD_PY" "brns_python" || exit 1

        [ "$SKIP_REF" = false ] && run_version "Reference" "$BUILD_REF" "brns_reference" "$RESULTS_DIR/reference"
        run_version "Python" "$BUILD_PY" "brns_python" "$RESULTS_DIR/python" || exit 1

        if [ "$SKIP_REF" = false ]; then
            compare_results
        else
            echo -e "${YELLOW}Comparison not possible (no reference)${NC}"
        fi
        ;;

    build)
        [ "$SKIP_REF" = false ] && build_version "Reference" "$REFERENCE_GEN" "$BUILD_REF" "brns_reference"
        build_version "Python" "$PYTHON_GEN" "$BUILD_PY" "brns_python" || exit 1
        ;;

    run)
        [ "$SKIP_REF" = false ] && run_version "Reference" "$BUILD_REF" "brns_reference" "$RESULTS_DIR/reference"
        run_version "Python" "$BUILD_PY" "brns_python" "$RESULTS_DIR/python" || exit 1
        ;;

    compare)
        if [ "$SKIP_REF" = false ]; then
            compare_results
        else
            echo -e "${YELLOW}Comparison not possible (no reference)${NC}"
        fi
        ;;

    *)
        echo -e "${RED}Unknown action: $ACTION${NC}"
        echo "Allowed: all, build, run, compare"
        exit 1
        ;;
esac

echo ""
echo -e "${GREEN}✓ Done!${NC}"
