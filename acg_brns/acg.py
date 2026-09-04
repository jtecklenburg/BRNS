"""
Automatic code generation for BRNS.

This module contains the Python implementation of the Maple ACG routines used to
emit problem-specific Fortran source files for BRNS simulations. It mirrors the
Maple naming scheme and keeps the generated code aligned with the legacy BRNS
workflow.

Function list (Maple-compatible names):
- `acg0(...)`: global geometry and parameters (`common_geo.inc`)
- `acg1(...)`: boundary conditions (`boundaries.f`)
- `acg2(...)`: molecular diffusion (`molecular.f`)
- `acg3(...)`: biogeochemical parameters (`biogeo.f`)
- `acg4(...)`: residual vector (`residual.f`)
- `acg5(...)`: Jacobian matrix (`jacobian.f`)
- `acg7(...)`: output (`output.f`)
- `acg8(...)`: physical parameters (`basic.f`)
- `acg12(...)`: initial conditions (`initialcond.f`)
- `acg13(...)`: steady-state rates (`ssrates.f`)
- `acg14(...)`: non-transported species (`notransport.f`)
- `acg15(...)`: reaction rates (`rates.f`)
- `acg16(...)`: solid identification (`issolid.f`)
- `acg17(...)`: terminal electron acceptor cascade (`switches.f` - TEAC)
- `acg17a(...)`: spatial switches (`switches.f`)
- `acg17b(...)`: parameter array (`parameters.f`)
- `acg17c(...)`: variable porosity (`varporosity.f`)
- `acg18(...)`: time-step initialization include (`inittimestep.inc`)

Notes:
- All names follow the Maple original (see `proc0903-M.md`) for direct
  equivalence.
- The function numbering intentionally skips `acg6` and `acg9`-`acg11` to match
  the Maple convention.
- Detailed behavior is documented in the docstrings of the individual methods.
"""

from pathlib import Path
from typing import List, Dict, Union, Set, Any, Tuple
from sympy import Symbol, diff, symbols as sp_symbols
from sympy.core.function import AppliedUndef

# Import Gaussian elimination module
from .gaussian_elimination import run_gaussian_elimination

# Import macrofor with proper error handling
try:
    from macrofor.api import (
        genfor, declaref, equalf, parameterf, commonf,
        subroutinem, callf, if_then_m, if_then_f, elsef, endiff,
        dom, openf, closef, writem, readf, readm,
        formatf, programm, commentf, set_fortran_style
    )
    MACROFOR_AVAILABLE = True
except ImportError:
    MACROFOR_AVAILABLE = False
    print("Warning: macrofor not installed. ACGModule will not work.")

from sympy.printing.fortran import fcode


class ACGModule:
    """
    Main class for Automatic Code Generation (ACG).

    Generates Fortran files for BRNS simulations based on
    symbolic reaction network definitions.

    Attributes:
        output_dir (Path): Output directory for generated Fortran files
    """

    def __init__(self, output_dir: str):
        """
        Initialize the ACG module.
        
        Fixed to F77 format only as required by BRNS.

        Args:
            output_dir: Path to output directory

        Raises:
            ImportError: If macrofor is not installed
        """
        if not MACROFOR_AVAILABLE:
            raise ImportError(
                "macrofor is required for ACGModule. "
                "Install it from the macrofor package."
            )

        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Always use F77 style (fixed format)
        set_fortran_style('f77')
        self.verbose = False  # Default: no verbose output
        print(f"ACG initialized: output={output_dir}, style=f77")

    def _verify_macrofor_available(self) -> None:
        """
        Verify macrofor is available at runtime.
        
        Raises:
            RuntimeError: If macrofor is not installed or not available
        """
        if not MACROFOR_AVAILABLE:
            raise RuntimeError(
                "macrofor is required for ACG code generation but is not installed. "
                "Install it with: pip install macrofor"
            )

    def _python_to_fortran_double(self, value: Union[int, float]) -> str:
        """
        Convert Python number to Fortran double-precision literal.
        Uses Maple-like compact notation (e.g. 0.1D-1 instead of 1.0D-2).

        Args:
            value: Python int or float

        Returns:
            Fortran double-precision string (e.g. '0.1D1', '0.D0')
        """
        # Type checking: ensure value is numeric
        if isinstance(value, str):
            try:
                value = float(value)
            except (ValueError, TypeError):
                raise TypeError(
                    f"Expected numeric value for Fortran conversion, got string '{value}'. "
                    "String values must be evaluated to numeric before code generation."
                )
        elif not isinstance(value, (int, float)):
            raise TypeError(
                f"Expected numeric value (int or float), got {type(value).__name__}: {value}"
            )
        
        if value == 0 or value == 0.0:
            return '0.D0'
        elif value == 1 or value == 1.0:
            return '0.1D1'
        else:
            # Convert to scientific notation
            sci_str = f'{value:.16E}'
            # Parse mantissa and exponent
            if 'E' in sci_str:
                mantissa, exp_str = sci_str.split('E')
                exp_val = int(exp_str)
                mantissa_float = float(mantissa)

                # Maple style: shift decimal point one place left
                # e.g. 1.0e-2 → 0.1e-1
                if mantissa_float != 0:
                    mantissa_float = mantissa_float / 10.0
                    exp_val = exp_val + 1

                # Format mantissa (remove trailing zeros)
                mantissa_str = f'{mantissa_float:.15f}'.rstrip('0').rstrip('.')

                # Build Fortran string (NO + for positive exponent)
                if exp_val == 0:
                    return f'{mantissa_str}D0'
                elif exp_val > 0:
                    return f'{mantissa_str}D{exp_val}'
                else:
                    return f'{mantissa_str}D{exp_val}'
            else:
                # No scientific notation needed
                return f'{value}D0'

    def _sympy_to_fortran(self, expr, precision: str = 'double') -> str:
        """
        Convert SymPy expression to Fortran string.

        Args:
            expr: SymPy expression
            precision: 'double' or 'single'

        Returns:
            Fortran-compatible string representation
        """
        if expr == 0 or expr == 0.0:
            return '0.D0'

        expr, array_ref_map = self._prepare_array_references(expr)

        # Base conversion with SymPy
        fortran_code = fcode(expr, standard=77, source_format='free')

        for placeholder, array_ref in array_ref_map.items():
            fortran_code = fortran_code.replace(placeholder, array_ref)

        # Adjustments for BRNS conventions
        # Flatten to single line and remove free-form continuation markers
        fortran_code = fortran_code.replace('\n', ' ')  # Single line
        fortran_code = fortran_code.replace('&', ' ')   # Remove free-form continuations
        fortran_code = fortran_code.replace('  ', ' ')  # Multiple spaces

        # Fortran77-specific functions (double precision)
        if precision == 'double':
            fortran_code = fortran_code.replace('max(', 'dmax1(')
            fortran_code = fortran_code.replace('min(', 'dmin1(')
            fortran_code = fortran_code.replace('exp(', 'dexp(')
            fortran_code = fortran_code.replace('log(', 'dlog(')
            fortran_code = fortran_code.replace('sqrt(', 'dsqrt(')
            fortran_code = fortran_code.replace('abs(', 'dabs(')

        return fortran_code.strip()

    def _prepare_array_references(self, expr) -> Tuple[Any, Dict[str, str]]:
        """
        Replace SymPy-internal array-like references with printable placeholders.

        SymPy parses strings like ``sp(3,j)`` as an undefined function call,
        which ``fcode`` cannot print in strict Fortran mode. We temporarily
        replace these nodes with plain symbols, let ``fcode`` print the
        expression, and then restore the original Fortran array syntax.

        Returns:
            Tuple of (possibly rewritten expression, placeholder->array ref map)
        """
        if not hasattr(expr, 'atoms'):
            return expr, {}

        replacement_map = {}
        placeholder_map: Dict[str, str] = {}
        placeholder_index = 0

        def register_array_ref(node, ref_text: str) -> None:
            nonlocal placeholder_index
            placeholder_name = f'__arrref_{placeholder_index}__'
            placeholder_index += 1
            replacement_map[node] = Symbol(placeholder_name)
            placeholder_map[placeholder_name] = ref_text

        for atom in sorted(expr.atoms(AppliedUndef), key=str):
            func_name = getattr(atom.func, '__name__', str(atom.func))
            if func_name not in {'sp', 'spold'}:
                continue
            args_text = ','.join(str(arg) for arg in atom.args)
            register_array_ref(atom, f'{func_name}({args_text})')

        for atom in sorted(expr.free_symbols, key=str):
            atom_name = getattr(atom, 'name', str(atom))
            if atom_name.startswith('sp(') or atom_name.startswith('spold('):
                register_array_ref(atom, atom_name)

        if not replacement_map:
            return expr, {}

        return expr.xreplace(replacement_map), placeholder_map

    # ======================================================================
    # ACG0: common_geo.inc - Geometry and global parameters
    # ======================================================================

    def acg0(
        self,
        nsolids: int,
        ndissolved: int,
        nreactions: int,
        nnodes: int,
        bio_names: List[str],
        bio_vals: List[float] = None,
        phys_names: List[str] = None,
        phys_vals: List[float] = None,
        phys_names2: List[str] = None,
        phys_vals2: List[int] = None,
    ) -> None:
        """
        Generate common_geo.inc - Global geometry and parameters.

        **Maple equivalent:** acg0(nsolids, ndissolved, ncompo, nreactions,
                                   bio_name, phys_name, phys_val, phys_name2,
                                   phys_val2, dir_f, nnodes)

        Args:
            nsolids: Number of solid species
            ndissolved: Number of dissolved species
            nreactions: Number of reactions
            nnodes: Number of grid points
            bio_names: Names of biogeochemical parameters (e.g. ['k_decay'])
            bio_vals: Values (optional, documentation only)
            phys_names: Names of physical parameters (real*8)
            phys_vals: Values (optional)
            phys_names2: Names of integer parameters
            phys_vals2: Values (optional)
            
        Raises:
            ValueError: If dimensions are invalid
        """
        # Verify macrofor is available
        self._verify_macrofor_available()
        
        # Validate inputs
        if nsolids < 0 or ndissolved < 0 or nreactions <= 0 or nnodes <= 0:
            raise ValueError(
                f"Invalid dimensions: nsolids={nsolids}, ndissolved={ndissolved}, "
                f"nreactions={nreactions}, nnodes={nnodes}"
            )
        if bio_names is None:
            bio_names = []
        
        ncomp = nsolids + ndissolved
        nx = 2 * nnodes - 1

        if nx < 1:
            nx = 1

        # Parameter list
        param_list = [
            f'nsolid={nsolids}',
            f'ndiss={ndissolved}',
            f'ncomp={ncomp}',
            f'nreac={nreactions}',
            f'nx={nx}',
        ]

        # Build code list
        code = [declaref('implicit real*8', ['a-h', 'o-z']), parameterf(param_list)]

        # Biogeochemical parameters
        if bio_names:
            code.extend(
                [declaref('real*8', bio_names), commonf('kinetics', bio_names)]
            )

        # Physical parameters (real*8)
        if phys_names:
            code.extend(
                [declaref('real*8', phys_names), commonf('physics', phys_names)]
            )

        # Integer parameters
        if phys_names2:
            code.extend(
                [declaref('integer', phys_names2), commonf('physics2', phys_names2)]
            )

        # Generate file
        output_file = self.output_dir / 'common_geo.inc'
        genfor(str(output_file), code)
        print(f"✓ Generated: {output_file.name}")

    # ======================================================================
    # ACG1: boundaries.f - Boundary conditions
    # ======================================================================

    def acg1(
        self,
        ncompo: int,
        type_up: List[int],
        bnddata_up: List[float],
        type_down: List[int],
        bnddata_down: List[float],
    ) -> None:
        """
        Generate boundaries.f - Boundary conditions.

        **Maple equivalent:** acg1(type_up, bnddata_up, type_down, bnddata_down, dir_f)

        Args:
            ncompo: Number of components
            type_up: Boundary condition types top (0=Dirichlet, 1=Neumann)
            bnddata_up: Boundary values top
            type_down: Boundary condition types bottom
            bnddata_down: Boundary values bottom
        """
        code = [
            declaref('include', ["'common_geo.inc'"]),
            declaref('include', ["'common.inc'"]),
            equalf('j', '1'),
        ]

        # Upper boundary - first all spb values, then all ibc values (Maple style)
        for i in range(1, ncompo + 1):
            code.append(
                equalf(f'spb({i},1)', self._python_to_fortran_double(bnddata_up[i - 1]))
            )
        for i in range(1, ncompo + 1):
            code.append(equalf(f'ibc({i},1)', str(type_up[i - 1])))

        code.append(equalf('j', 'nx'))

        # Lower boundary - first all spb values, then all ibc values (Maple style)
        for i in range(1, ncompo + 1):
            code.append(
                equalf(f'spb({i},2)', self._python_to_fortran_double(bnddata_down[i - 1]))
            )
        for i in range(1, ncompo + 1):
            code.append(equalf(f'ibc({i},2)', str(type_down[i - 1])))

        # Generate subroutine
        output_file = self.output_dir / 'boundaries.f'
        genfor(str(output_file), [subroutinem('boundaries', [], code)])
        print(f"✓ Generated: {output_file.name}")

    # ======================================================================
    # ACG2: molecular.f - Molecular diffusion
    # ======================================================================

    def acg2(
        self, ncompo: int, diffdata: List[float], alphadata: List[float]
    ) -> None:
        """
        Generate molecular.f - Molecular diffusion and tortuosity.

        **Maple equivalent:** acg2(ncompo, diffdata, alphadata, dir_f)

        Args:
            ncompo: Number of components
            diffdata: Diffusion coefficients (dsol_0) - must have length ncompo
            alphadata: Tortuosity factors (f_T) - must have length ncompo
            
        Note: Maple writes ALL ncompo entries (including solids with dsol_0=0).
              Python must do the same to match Maple output exactly.
        """
        # Validate input lengths
        if len(diffdata) != ncompo:
            raise ValueError(f"diffdata length {len(diffdata)} != ncompo {ncompo}")
        if len(alphadata) != ncompo:
            raise ValueError(f"alphadata length {len(alphadata)} != ncompo {ncompo}")
        
        code = [
            declaref('include', ["'common_geo.inc'"]),
            declaref('include', ["'common.inc'"]),
        ]

        # Diffusion coefficients - write ALL ncompo entries (like Maple)
        for i in range(1, ncompo + 1):
            code.append(
                equalf(f'dsol_0({i})', self._python_to_fortran_double(diffdata[i - 1]))
            )

        # Tortuosity - write ALL ncompo entries (like Maple)
        for i in range(1, ncompo + 1):
            code.append(
                equalf(f'f_T({i})', self._python_to_fortran_double(alphadata[i - 1]))
            )

        output_file = self.output_dir / 'molecular.f'
        genfor(str(output_file), [subroutinem('molecular', [], code)])
        print(f"✓ Generated: {output_file.name}")

    # ======================================================================
    # ACG3: biogeo.f - Biogeochemical parameters
    # ======================================================================

    def acg3(self, bio_names: List[str], bio_vals: List[float]) -> None:
        """
        Generate biogeo.f - Sets biogeochemical parameter values.

        **Maple equivalent:** acg3(nparam, bio_name, bio_val, dir_f)

        Args:
            bio_names: Parameter names
            bio_vals: Parameter values
        """
        code = [
            declaref('include', ["'common_geo.inc'"]),
            declaref('include', ["'common.inc'"]),
        ]

        for name, val in zip(bio_names, bio_vals):
            code.append(equalf(name, self._python_to_fortran_double(val)))

        output_file = self.output_dir / 'biogeo.f'
        genfor(str(output_file), [subroutinem('biogeo', [], code)])
        print(f"✓ Generated: {output_file.name}")

    # ======================================================================
    # ACG4: residual.f - Residual vector
    # ======================================================================

    def acg4(self, ncompo: int, func_expressions: List) -> None:
        """
        Generate residual.f - Residual vector with implicit time discretization.

        **Maple equivalent:** acg4(ncompo, func, dir_f) with func from p9() modified

        The residual contains both chemical reactions and time discretization
        (implicit Euler):

        funcs(i) = R_i(C) - C_i^(n+1)/delt + C_i^n/delt

        where:
        - R_i(C) = chemical net rate for species i
        - C_i^(n+1) = sp(i,j) = current concentration
        - C_i^n = spold(i,j) = concentration at previous time step

        Args:
            ncompo: Number of components
            func_expressions: List of chemical net rates R_i(C)
        """
        code = [
            declaref('include', ["'common_geo.inc'"]),
            declaref('include', ["'common.inc'"]),
            declaref('dimension', ['funcs(ncomp)']),
            callf('switches', ['j']),
        ]

        for i in range(1, ncompo + 1):
            if hasattr(func_expressions[i - 1], 'free_symbols'):
                # SymPy expression
                func_str = self._sympy_to_fortran(func_expressions[i - 1])
            else:
                # String
                func_str = str(func_expressions[i - 1])

            # Write expression as-is (p9 already applied time discretization)
            code.append(equalf(f'funcs({i})', func_str))

        output_file = self.output_dir / 'residual.f'
        genfor(str(output_file), [subroutinem('residual', ['funcs', 'j'], code)])
        print(f"✓ Generated: {output_file.name}")

    # ======================================================================
    # ACG5: jacobian.f - Jacobian matrix (with SymPy computation)
    # ======================================================================

    def acg5(self, ncompo: int, jacobian_matrix: List[List], 
             conservation_rows: Set[int] = None) -> None:
        """
        Generate jacobian.f - Jacobian matrix for implicit time integration.

        **Maple equivalent:** acg5(ncompo, dir_f)

        IMPORTANT: The input jacobian_matrix should ALREADY contain time discretization
        terms (-1/delt on diagonal) from p9! This function just writes the matrix
        to Fortran, without adding any additional terms.

        The Jacobian from p10 already contains:
        pd(i,j) = ∂R_i/∂C_j                    for i ≠ j
        pd(i,i) = ∂R_i/∂C_i - 1/delt          for reactive rows
        pd(i,j) = conservation_basis[k,j]     for conservation rows

        Maple equivalent: acg5(pd, dir_f) - just writes the Jacobian to file

        Args:
            ncompo: Number of components
            jacobian_matrix: ncompo x ncompo matrix (SymPy or strings)
                            ALREADY contains time discretization from p9!
            conservation_rows: Set of row indices (0-based) that are conservation laws
                             (Currently not used - Jacobian already correct from p10)
        """
        if conservation_rows is None:
            conservation_rows = set()
        
        code = [
            declaref('include', ["'common_geo.inc'"]),
            declaref('include', ["'common.inc'"]),
            declaref('dimension', ['pd(ncomp,ncomp)']),
            callf('switches', ['j']),
        ]

        # Matrix entries - write directly from jacobian_matrix
        # NO additional -1/delt terms needed (already in the matrix from p9!)
        for i in range(1, ncompo + 1):
            for j in range(1, ncompo + 1):
                entry = jacobian_matrix[i - 1][j - 1]
                if hasattr(entry, 'free_symbols'):
                    entry_str = self._sympy_to_fortran(entry)
                else:
                    entry_str = str(entry)

                # Write entry as-is (already contains all necessary terms)
                code.append(equalf(f'pd({i},{j})', entry_str))

        output_file = self.output_dir / 'jacobian.f'
        genfor(str(output_file), [subroutinem('jacobian', ['pd', 'j'], code)])
        print(f"✓ Generated: {output_file.name}")

    # ======================================================================
    # ACG_OPTIMIZE: System optimization via Gaussian elimination (p4-p10)
    # ======================================================================

    def generate_optimized_residual_jacobian(
        self,
        ncompo: int,
        equations: List,
        reaction_rates: List,
        neqrxns: int = 0,
        equilibrium_eqns: List = None,
        inert_components: Set[int] = None,
        verbose: bool = False
    ) -> Dict:
        """
        Generate optimized residual and Jacobian via Gaussian elimination.
        
        This implements the full p4-p10 pipeline from Maple ACG:
        - p4: Build coefficient matrix from equations
        - p5: Gaussian elimination with pivoting
        - p6: Extract stoichiometric matrices
        - p7: Build net rate expressions
        - p8: Handle equilibrium reactions and inert components
        - p9: Add implicit Euler time discretization
        - p10: Compute Jacobian matrix
        
        The optimization reduces the system to minimal form by eliminating
        redundant equations and identifying active reaction pathways.
        
        This is computationally efficient for large systems, especially
        when coupled to 2D/3D transport (evaluated at ~10,000+ grid points).
        
        Maple equivalent: Full ACG pipeline (p4-p10)
        
        Args:
            ncompo: Number of components
            equations: List of symbolic equilibrium equations
            reaction_rates: List of reaction rate symbols [r1, r2, ...]
            neqrxns: Number of equilibrium reactions
            equilibrium_eqns: Expressions for equilibrium reactions
            inert_components: Set of component indices that don't react
            verbose: Print debugging information
        
        Returns:
            Dict with optimized system:
            - 'func': Reduced residual equations
            - 'jacobian': Optimized Jacobian matrix
            - 'matrixM': Stoichiometric matrix
            - 'rightM': Reaction coefficient matrix
            - 'pivot_rows': Rows used in elimination
            - 'optimization_info': Optimization statistics
        """
        # Create variable symbols
        variables = sp_symbols(f'C0:{ncompo}')
        if not isinstance(variables, tuple):
            variables = (variables,)
        variables = list(variables)
        
        variables_old = sp_symbols(f'C0_old:{ncompo}')
        if not isinstance(variables_old, tuple):
            variables_old = (variables_old,)
        variables_old = list(variables_old)
        
        delt = Symbol('delt', positive=True, real=True)
        
        if self.verbose or verbose:
            print("\n" + "="*70)
            print("GAUSSIAN ELIMINATION OPTIMIZATION (p4-p10)")
            print("="*70)
            print(f"System: {ncompo} components, {len(reaction_rates)} reactions")
            print(f"Equilibrium reactions: {neqrxns}")
            print(f"Inert components: {inert_components if inert_components else 'none'}")
        
        # Run Gaussian elimination pipeline
        result = run_gaussian_elimination(
            equations=equations,
            variables=variables,
            reactions=reaction_rates,
            variables_old=variables_old,
            delt=delt,
            neqrxns=neqrxns,
            equilibrium_eqns=equilibrium_eqns,
            inert_components=inert_components if inert_components else set(),
            verbose=(self.verbose or verbose)
        )
        
        # Compute optimization statistics
        original_size = ncompo ** 2  # Original Jacobian size
        # Count non-zero Jacobian entries (after optimization)
        jacobian = result['jacobian']
        nonzero_entries = sum(1 for i in range(ncompo) for j in range(ncompo) 
                             if jacobian[i, j] != 0)
        
        optimization_info = {
            'original_jacobian_size': original_size,
            'nonzero_entries': nonzero_entries,
            'sparsity_ratio': 1 - (nonzero_entries / original_size) if original_size > 0 else 0,
            'pivot_rows_count': len(result['pivot_rows']),
            'inert_components': len(inert_components) if inert_components else 0,
        }
        
        result['optimization_info'] = optimization_info
        
        if self.verbose or verbose:
            print(f"\nOptimization Results:")
            print(f"  Original Jacobian size: {original_size}")
            print(f"  Non-zero entries: {nonzero_entries}")
            print(f"  Sparsity ratio: {optimization_info['sparsity_ratio']:.1%}")
            print(f"  Pivot rows: {result['pivot_rows']}")
            print("="*70 + "\n")
        
        return result

    # ======================================================================
    # ACG16: issolid.f - Identifies solid species
    # ======================================================================

    def acg16(self, nsolids: int = 0, listsolids: List[int] = None) -> None:
        """
        Generate issolid.f - Identifies solid vs. dissolved species.

        **Maple equivalent:** acg16(nsolids, listsolids, dir_f)

        Args:
            nsolids: Number of solid species (default: 0)
            listsolids: Indices of solid species (e.g. [1, 3]) (default: None/[])
        """
        if listsolids is None:
            listsolids = []

        code = [
            declaref('include', ["'common_geo.inc'"]),
            declaref('include', ["'common.inc'"]),
        ]

        # IF-THEN constructs for each solid species (only if present)
        if listsolids:
            for solid_idx in listsolids:
                condition = f'k.eq.{solid_idx}'
                code.append(if_then_m(condition, [equalf('isolid', '1')]))
        # Otherwise the subroutine remains empty (only includes)

        output_file = self.output_dir / 'issolid.f'
        genfor(str(output_file), [subroutinem('issolid', ['k', 'isolid'], code)])
        print(f"✓ Generated: {output_file.name}")

    # ======================================================================
    # ACG7: output.f - Data output (Maple equivalent: acg7)
    # ======================================================================

    def acg7(
        self,
        noutput: int,
        nroutput: int,
        listoutput: List[int],
        listroutput: List[int],
        file_names: List[str],
        file_rnames: List[str],
        time_iniout: float,
        time_intvout: float,
    ) -> None:
        """
        Generate output.f - Writes concentration and rate data.

        **Maple equivalent:** acg7(ncompo, noutput, nroutput, listoutput,
                                   listroutput, file_names, file_rnames,
                                   time_iniout, time_intvout, dir_f)

        Args:
            noutput: Number of concentration outputs
            nroutput: Number of rate outputs
            listoutput: Species indices for concentrations
            listroutput: Reaction indices for rates
            file_names: File names for concentrations
            file_rnames: File names for rates
            time_iniout: Initial output time
            time_intvout: Output interval
        """
        code = [
            declaref('include', ["'common_geo.inc'"]),
            declaref('include', ["'common.inc'"]),
            declaref('real*8', ['time']),
        ]

        # Initialization at t=0, j=1: Create files
        init_block = []

        # Order as in reference: .dat files, then .inp files
        # Concentration files (.dat)
        for fname in file_names:
            init_block.append(openf(11, f'{fname}.dat', 'replace'))
            init_block.append(closef(11))

        # Rate files (.dat)
        for fname in file_rnames:
            init_block.append(openf(11, f'{fname}.dat', 'replace'))
            init_block.append(closef(11))

        # Concentration files (.inp)
        for fname in file_names:
            init_block.append(openf(11, f'{fname}.inp', 'replace'))
            init_block.append(closef(11))

        # Initialize output variables (as integer literals)
        init_block.extend(
            [
                equalf('v_out', str(int(time_iniout))),
                equalf('v_int', str(int(time_intvout))),
            ]
        )

        code.append(if_then_m('nt.eq.1.and.j.eq.1', init_block))

        # Output block: When time reached (condition as in reference)
        output_block = []

        # Write concentrations (.dat)
        for i, fname in enumerate(file_names):
            sp_idx = listoutput[i]
            output_block.extend(
                [
                    openf(11, f'{fname}.dat', 'old', 'append'),
                    writem(11, ['1x,e14.7,2x,f12.4'], [f'sp({sp_idx},j)', 'depth']),
                    closef(11),
                ]
            )

        # Write rates
        for i, fname in enumerate(file_rnames):
            r_idx = listroutput[i]
            output_block.extend(
                [
                    openf(11, f'{fname}.dat', 'old', 'append'),
                    writem(11, ['1x,e14.7,2x,f12.4'], [f'r({r_idx},j)', 'depth']),
                    closef(11),
                ]
            )

        # Increment v_out at j=nx (nested IF block, no spaces around +)
        output_block.append(if_then_m('j.eq.nx', [equalf('v_out', 'v_out+v_int')]))

        code.append(
            if_then_m('time.le.v_out.and.v_out.lt.time+delt', output_block)
        )

        # Final output at time=endt (.inp)
        final_block = []
        for i, fname in enumerate(file_names):
            sp_idx = listoutput[i]
            final_block.extend(
                [
                    openf(11, f'{fname}.inp', 'old', 'append'),
                    writem(11, ['1x,e14.7,2x,f12.4'], [f'sp({sp_idx},j)', 'depth']),
                    closef(11),
                ]
            )

        code.append(if_then_m('time.eq.endt', final_block))

        output_file = self.output_dir / 'output.f'
        genfor(
            str(output_file),
            [subroutinem('out', ['j', 'nt', 'time', 'depth', 'v_out', 'v_int'], code)],
        )
        print(f"✓ Generated: {output_file.name}")

    # ======================================================================
    # ACG8: basic.f - Physical base parameters (Maple equivalent: acg8)
    # ======================================================================

    def acg8(
        self,
        phys_names: List[str],
        phys_vals: List[float],
        phys_names2: List[str] = None,
        phys_vals2: List[int] = None,
    ) -> None:
        """
        Generate basic.f - Sets physical parameter values.

        **Maple equivalent:** acg8(nparphys, phys_name, phys_val,
                                   nparphys2, phys_name2, phys_val2, dir_f)

        Args:
            phys_names: Names of physical parameters (real*8)
            phys_vals: Values
            phys_names2: Names of integer parameters
            phys_vals2: Integer values
        """
        code = [
            declaref('include', ["'common_geo.inc'"]),
            declaref('include', ["'common.inc'"]),
        ]

        # Real parameters
        for name, val in zip(phys_names, phys_vals):
            code.append(equalf(name, self._python_to_fortran_double(val)))

        # Integer parameters
        if phys_names2 and phys_vals2:
            for name, val in zip(phys_names2, phys_vals2):
                code.append(equalf(name, str(val)))

        output_file = self.output_dir / 'basic.f'
        genfor(str(output_file), [subroutinem('basic', [], code)])
        print(f"✓ Generated: {output_file.name}")

    # ======================================================================
    # ACG12: initialcond.f - Initial conditions (Maple equivalent: acg12)
    # ======================================================================

    def acg12(
        self,
        ncompo: int,
        ic_mode: int,
        iniconc: List[float] = None,
        listinput: List[int] = None,
        file_in_names: List[str] = None,
    ) -> None:
        """
        Generate initialcond.f - Initial conditions.

        **Maple equivalent:** acg12(vic, iniconc, ncompo, listinput,
                                    file_in_names, dir_f)

        Args:
            ncompo: Number of components
            ic_mode: Mode (1=file, 2=constant, 3=profiles)
            iniconc: Initial concentrations (for ic_mode=2)
            listinput: Species indices (for ic_mode=3)
            file_in_names: File names (for ic_mode=3)
        """
        code = [
            declaref('include', ["'common_geo.inc'"]),
            declaref('include', ["'common.inc'"]),
            declaref('include', ["'common_drive.inc'"]),
            declaref('real*8', ['spi(ncomp)']),
        ]

        # ic=1: Read from initialconc.txt
        ic1_block = [
            "open(unit=85, file='initialconc.txt', status='old')",
            commentf(''),
            dom('j', 1, 'nx', 'read(85,*) de, (sp(i,j), i=1, ncomp)', 2),
            commentf(''),
            'close(85)',
        ]

        code.append(if_then_m('ic.eq.1', ic1_block))

        # ic=2: Constant values
        if iniconc:
            ic2_block = []
            for i in range(1, ncompo + 1):
                ic2_block.append(
                    equalf(f'spi({i})', self._python_to_fortran_double(iniconc[i - 1]))
                )

            # Nested DO loops
            ic2_block.append(commentf(''))
            ic2_block.append(
                dom(
                    'i',
                    1,
                    'ncomp',
                    [
                        commentf(''),
                        dom('j', 1, 'nx', equalf('sp(i,j)', 'spi(i)'), 2),
                        commentf(''),
                    ],
                )
            )
            ic2_block.append(commentf(''))

            code.append(if_then_m('ic.eq.2', ic2_block))

        # ic=3: Read one profile file per species (Maple-compatible)
        if listinput and file_in_names:
            ic3_block = []

            for i, fname in enumerate(file_in_names):
                unit_num = 101 + i
                file_name = fname if str(fname).lower().endswith('.inp') else f'{fname}.inp'
                ic3_block.append(
                    f"open(unit={unit_num}, file='{file_name}', status='old')"
                )

            ic3_block.append(commentf(''))
            ic3_block.append('do 1003 j=1, nx, 2')

            for i, sp_idx in enumerate(listinput):
                unit_num = 101 + i
                fmt_label = 2000 + i
                ic3_block.append(f'read({unit_num},{fmt_label}) sp({sp_idx},j), depth')
                ic3_block.append(f'{fmt_label} format(1x,e14.7,2x,f8.4)')

            ic3_block.append('1003 continue')
            ic3_block.append(commentf(''))

            for i in range(len(file_in_names)):
                unit_num = 101 + i
                ic3_block.append(f'close({unit_num})')

            code.append(if_then_m('ic.eq.3', ic3_block))

        output_file = self.output_dir / 'initialcond.f'
        genfor(str(output_file), [subroutinem('initialcond', [], code)])
        print(f"✓ Generated: {output_file.name}")

    # ======================================================================
    # ACG13: ssrates.f - Steady-state rates (Maple equivalent: acg13)
    # ======================================================================

    def acg13(
        self,
        ncompo: int,
        rate_expressions: List,
        species_symbols: List[Symbol],
        subst_dict: Dict,
    ) -> None:
        """
        Generate ssrates.f - Steady-state reaction rates and derivatives.

        **Maple equivalent:** acg13(ncompo)

        This function generates the ssrates subroutine that computes for each species
        the net rate (rat) and its derivative (drdc).

        IMPORTANT: rate_expressions must be the ORIGINAL SymPy expressions!
        The function first computes derivatives (in symbolic space),
        then substitutes both expressions to Fortran array notation.

        Workflow:
        1. Compute derivative: ∂(rate)/∂(species) with original symbols
        2. Substitute both: diss_a → sp(1,j), k_deg stays k_deg
        3. Convert to Fortran code

        Args:
            ncompo: Number of components
            rate_expressions: List of SymPy expressions (ORIGINAL, not substituted!)
            species_symbols: List of SymPy symbols (for derivatives)
            subst_dict: Substitution dictionary (Symbol → Fortran array notation)
        """
        code = [
            declaref('include', ["'common_geo.inc'"]),
            declaref('include', ["'common.inc'"]),
            callf('switches', ['j']),
        ]

        # For each species: rate and derivative
        for i in range(1, ncompo + 1):
            rate_expr = rate_expressions[i - 1]  # Original SymPy
            species_sym = species_symbols[i - 1]

            # 1. FIRST compute derivative (before substitution!)
            deriv_expr = diff(rate_expr, species_sym)

            # 2. THEN substitute both expressions
            rate_substituted = rate_expr.subs(subst_dict)
            deriv_substituted = deriv_expr.subs(subst_dict)

            # 3. Convert to Fortran code
            rate_fortran = self._sympy_to_fortran(rate_substituted)
            deriv_fortran = self._sympy_to_fortran(deriv_substituted)

            # IF block for this species
            code.append(
                if_then_m(
                    f'isp.eq.{i}',
                    [equalf('rat', rate_fortran), equalf('drdc', deriv_fortran)],
                )
            )

        output_file = self.output_dir / 'ssrates.f'
        genfor(str(output_file), [subroutinem('ssrates', ['rat', 'drdc', 'isp', 'j'], code)])
        print(f"✓ Generated: {output_file.name}")

    # ======================================================================
    # ACG14: notransport.f - Marks non-transported species (Maple equivalent: acg14)
    # ======================================================================

    def acg14(self, listnotransp: List[int] = None) -> None:
        """
        Generate notransport.f - Marks species without transport.

        **Maple equivalent:** acg14(nrnotransp, listnotransp, dir_f)

        Args:
            listnotransp: Indices of non-transported species (default: None/[])
        """
        if listnotransp is None:
            listnotransp = []

        code = [
            declaref('include', ["'common_geo.inc'"]),
            declaref('include', ["'common.inc'"]),
            declaref('integer', ['k', 'itransp']),
        ]

        # IF-THEN for each non-transported species (only if present)
        if listnotransp:
            for idx in listnotransp:
                code.append(if_then_m(f'k.eq.{idx}', [equalf('itransp', '1')]))
        # Otherwise the subroutine remains empty (only includes and declarations)

        output_file = self.output_dir / 'notransport.f'
        genfor(str(output_file), [subroutinem('notransport', ['k', 'itransp'], code)])
        print(f"✓ Generated: {output_file.name}")

    # ======================================================================
    # ACG15: rates.f - Reaction rates (Maple equivalent: acg15)
    # ======================================================================

    def acg15(self, nreactions: int, rate_expressions: List) -> None:
        """
        Generate rates.f - Computes reaction rates.

        **Maple equivalent:** acg15(nreactions)

        Args:
            nreactions: Number of reactions
            rate_expressions: List of SymPy expressions or strings
        """
        code = [
            declaref('include', ["'common_geo.inc'"]),
            declaref('include', ["'common.inc'"]),
            callf('switches', ['j']),
        ]

        # Rate assignments
        for i in range(1, nreactions + 1):
            rate_expr = rate_expressions[i - 1]

            if hasattr(rate_expr, 'free_symbols'):
                rate_fortran = self._sympy_to_fortran(rate_expr)
            else:
                rate_fortran = str(rate_expr)

            code.append(equalf(f'r({i},j)', rate_fortran))

        output_file = self.output_dir / 'rates.f'
        genfor(str(output_file), [subroutinem('rates', ['j'], code)])
        print(f"✓ Generated: {output_file.name}")

    # ======================================================================
    # ACG17B: parameters.f - Parameter list (Maple equivalent: acg17b)
    # ======================================================================

    def acg17b(self, parameter_names: List[str]) -> None:
        """
        Generate parameters.f - Sets parameters from array.

        **Maple equivalent:** acg17b(nparameters, dir_f)

        Args:
            parameter_names: Parameter names
        """
        code = [
            declaref('include', ["'common_geo.inc'"]),
            declaref('real*8', ['paravalues(*)']),
        ]

        for i, pname in enumerate(parameter_names, 1):
            code.append(equalf(pname, f'paravalues({i})'))

        output_file = self.output_dir / 'parameters.f'
        genfor(str(output_file), [subroutinem('parameters', ['paravalues'], code)])
        print(f"✓ Generated: {output_file.name}")

    # ======================================================================
    # ACG17A: switches.f - Spatial switches (Maple equivalent: acg17a)
    # ======================================================================

    def acg17a(
        self,
        switch_conditions: Dict[str, str],
        switch_names: List[str] = None,
    ) -> None:
        """
        Generate switches.f - Spatial switches based on conditions.

        **Maple equivalent:** acg17a(nswitches, dir_f)

        Args:
            switch_conditions: Dict of {switch_name: condition_string}
                              e.g. {'ho2': 'sp(1,j) > kmo2'}
            switch_names: List of switch names (if None, from dict)
        """
        if switch_names is None:
            switch_names = list(switch_conditions.keys())

        code = [
            declaref('include', ["'common_geo.inc'"]),
            declaref('include', ["'common.inc'"]),
            declaref('integer', ['j']),
            declaref('real*8', ['x_pos']),
            declaref('real*8', ['y_pos']),
            declaref('real*8', ['z_pos']),
            equalf('x_pos', 'x(j)'),
        ]

        # For each switch
        for switch_name in switch_names:
            if switch_name in switch_conditions:
                condition = switch_conditions[switch_name]
                # Use nested if-then-else structure
                # Since if_then_m only supports then-block, we need to build manually
                if_body = [equalf(switch_name, '1.0')]
                else_body = [equalf(switch_name, '0.0')]

                # Create complete if-then-else structure as nested list
                code.extend(
                    [if_then_f(condition), *if_body, elsef(), *else_body, endiff()]
                )

        output_file = self.output_dir / 'switches.f'
        genfor(str(output_file), [subroutinem('switches', ['j'], code)])
        print(f"✓ Generated: {output_file.name}")

    # ======================================================================
    # ACG17: switches.f - Terminal electron acceptor cascade
    # ======================================================================

    def acg17(self, variables: List[Union[str, Symbol]]) -> None:
        """
        Generate switches.f - Terminal electron acceptor cascade.

        **Maple equivalent:** acg17(variables, dir_f)

        The cascade order is:
        O2 → NO3 → MnO2 → FeOH3 → SO4 → CH4

        Args:
            variables: List of species variables (SymPy symbols or strings)
        """

        def _var_name(var: Union[str, Symbol]) -> str:
            return var.name if isinstance(var, Symbol) else str(var)

        def _assignments(pairs: List[tuple]) -> List:
            return [equalf(name, value) for name, value in pairs]

        # Build index map for variable names (1-based for sp(i,j))
        var_names = [_var_name(v) for v in variables]
        index_map = {name: idx + 1 for idx, name in enumerate(var_names)}

        # Assignment blocks (Maple logic)
        assign_o2 = _assignments(
            [
                ('ho2', '1.0'),
                ('hno3', '0.0'),
                ('hmn4', '0.0'),
                ('hfe3', '0.0'),
                ('hso4', '0.0'),
            ]
        )
        assign_no3 = _assignments(
            [
                ('ho2', '0.0'),
                ('hno3', '1.0'),
                ('hno3f2', '1.0'),
                ('hmn4', '0.0'),
                ('hfe3', '0.0'),
                ('hso4', '0.0'),
            ]
        )
        assign_mno2 = _assignments(
            [
                ('ho2', '0.0'),
                ('hno3', '1.0'),
                ('hno3f2', '0.0'),
                ('hmn4', '1.0'),
                ('hmn4f2', '1.0'),
                ('hfe3', '0.0'),
                ('hso4', '0.0'),
            ]
        )
        assign_feoh3 = _assignments(
            [
                ('ho2', '0.0'),
                ('hno3', '1.0'),
                ('hno3f2', '0.0'),
                ('hmn4', '1.0'),
                ('hmn4f2', '0.0'),
                ('hfe3', '1.0'),
                ('hfe3f2', '1.0'),
                ('hso4', '0.0'),
            ]
        )
        assign_so4 = _assignments(
            [
                ('ho2', '0.0'),
                ('hno3', '1.0'),
                ('hno3f2', '0.0'),
                ('hmn4', '1.0'),
                ('hmn4f2', '0.0'),
                ('hfe3', '1.0'),
                ('hfe3f2', '0.0'),
                ('hso4', '1.0'),
                ('hso4f2', '1.0'),
            ]
        )
        assign_ch4 = _assignments(
            [
                ('ho2', '0.0'),
                ('hno3', '1.0'),
                ('hno3f2', '0.0'),
                ('hmn4', '1.0'),
                ('hmn4f2', '0.0'),
                ('hfe3', '1.0'),
                ('hfe3f2', '0.0'),
                ('hso4', '1.0'),
                ('hso4f2', '0.0'),
            ]
        )

        # Build conditional chain for acceptors present in variables
        acceptor_specs = [
            ('o2', 'kmo2', assign_o2),
            ('no3', 'kmno3', assign_no3),
            ('mno2', 'kmno2', assign_mno2),
            ('feoh3', 'kmfeoh3', assign_feoh3),
            ('so4', 'kmso4', assign_so4),
        ]

        conditions = []
        for name, km, assigns in acceptor_specs:
            if name in index_map:
                cond = f'sp({index_map[name]},j) > {km}'
                conditions.append((cond, assigns))

        default_assigns = assign_ch4 if 'ch4' in index_map else None

        code = [
            declaref('include', ["'common_geo.inc'"]),
            declaref('include', ["'common.inc'"]),
            declaref('integer', ['j']),
        ]

        if not conditions:
            if default_assigns:
                code.extend(default_assigns)
            else:
                code.append(commentf('No terminal acceptors found in variables.'))
        else:
            # Build nested if-then-else chain
            first_cond, first_assigns = conditions[0]
            code.append(if_then_f(first_cond))
            code.extend(first_assigns)

            for cond, assigns in conditions[1:]:
                code.append(elsef())
                code.append(if_then_f(cond))
                code.extend(assigns)

            if default_assigns:
                code.append(elsef())
                code.extend(default_assigns)

            # Close all nested IF blocks
            for _ in range(len(conditions)):
                code.append(endiff())

        output_file = self.output_dir / 'switches.f'
        genfor(str(output_file), [subroutinem('switches', ['j'], code)])
        print(f"✓ Generated: {output_file.name}")

    # ======================================================================
    # ACG17c: varporosity.f - Variable porosity update
    # ======================================================================

    def acg17c(
        self,
        varporosity: int,
        porosity_expr: Union[str, Symbol] = None,
        subst_dict: Dict = None,
    ) -> None:
        """
        Generate varporosity.f - Update porosity (space/time-dependent).

        **Maple equivalent:** acg17c(varporosity, dir_f)

        Args:
            varporosity: If > 0, use porosity_expr; else por0 remains unchanged.
            porosity_expr: SymPy expression or string for porosity definition.
            subst_dict: Optional substitution dict for Fortran array notation.
        """
        code = [
            declaref('include', ["'common_geo.inc'"]),
            declaref('include', ["'common.inc'"]),
        ]

        if varporosity > 0 and porosity_expr is not None:
            expr = porosity_expr
            if subst_dict and hasattr(expr, 'subs'):
                expr = expr.subs(subst_dict)

            if hasattr(expr, 'free_symbols'):
                porosity_str = self._sympy_to_fortran(expr)
            else:
                porosity_str = str(expr)

            code.append(equalf('por0', porosity_str))
        else:
            code.append(commentf('varporosity <= 0: por0 unchanged'))

        output_file = self.output_dir / 'varporosity.f'
        genfor(str(output_file), [subroutinem('updateporosity', ['mydummy'], code)])
        print(f"✓ Generated: {output_file.name}")

    # ======================================================================
    # ACG18: inittimestep.inc - Time step parameter initialization include
    # ======================================================================

    def acg18(
        self,
        dtold: Union[str, int, float] = '1.d-6',
        ttol: Union[str, int, float] = '5.d-2',
        tstep: Union[str, int, float] = '1.0d0',
        maxconc: Union[str, int, float] = '0.d0',
    ) -> None:
        """
        Generate inittimestep.inc - BLOCK DATA initialization for /tsparam/.

        Writes all timestep parameters to a dedicated include file:
        dtold, ttol, tstep, maxconc.

        Defaults match the current BRNS values:
        dtold/1.d-6/, ttol/5.d-2/, tstep/1.0d0/, maxconc/0.d0/

        Args:
            dtold: Previous timestep default value
            ttol: Timestep tolerance value
            tstep: Initial/base timestep value
            maxconc: Initial maximum concentration value
        """

        def _format_fortran_data_value(value: Union[str, int, float]) -> str:
            if isinstance(value, str):
                value_str = value.strip()
                if not value_str:
                    raise ValueError("Fortran DATA value string must not be empty")
                return value_str

            if not isinstance(value, (int, float)):
                raise TypeError(
                    "Expected str, int or float for timestep parameter, "
                    f"got {type(value).__name__}: {value}"
                )

            if float(value) == 0.0:
                return '0.d0'

            sci = f"{float(value):.16e}"
            mantissa, exponent = sci.split('e')
            mantissa = mantissa.rstrip('0').rstrip('.')
            if '.' not in mantissa:
                mantissa = f"{mantissa}.0"
            exponent_int = int(exponent)
            return f"{mantissa}d{exponent_int}"

        dtold_s = _format_fortran_data_value(dtold)
        ttol_s = _format_fortran_data_value(ttol)
        tstep_s = _format_fortran_data_value(tstep)
        maxconc_s = _format_fortran_data_value(maxconc)

        lines = [
            "c    $Id: inittimestep.inc 16 2026-07-08 10:56:47Z tecklenburg $",
            "      BLOCK DATA InitTimeStep",
            "        common/tsparam/dtold,ttol,tstep,maxconc",
            "        real*8 dtold,ttol,tstep,maxconc",
            (
                "        DATA "
                f"dtold/{dtold_s}/, ttol/{ttol_s}/, "
                f"tstep/{tstep_s}/, maxconc/{maxconc_s}/"
            ),
            "      END",
            "",
        ]

        output_file = self.output_dir / 'inittimestep.inc'
        output_file.write_text("\n".join(lines), encoding='utf-8')
        print(f"✓ Generated: {output_file.name}")


