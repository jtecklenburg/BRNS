"""
Gaussian Elimination Module for ACG
Python equivalent to Maple ACG functions p4() through p10()

This module reduces reaction network equations to minimal form by:
- Building the stoichiometric coefficient matrix (p4-p6)
- Computing the net rate equations (p7)
- Handling equilibrium reactions (p8)
- Adding implicit time discretization (p9)
- Computing the Jacobian matrix (p10)

Mathematical Framework:
1. p4(): Build augmented matrix [coefficients | identity] and perform Gaussian elimination
2. p5(): Row reduction to get reduced row echelon form
3. p6(): Extract stoichiometric matrices (rightM, matrixM)
4. p7(): Build residual expression from stoichiometric matrix
5. p8(): Handle equilibrium reactions and inert components
6. p9(): Add implicit Euler time discretization terms
7. p10(): Compute Jacobian matrix via symbolic differentiation

Author: Port from Maple ACG module
Date: 2025-02-10
"""

from typing import List, Dict, Tuple, Set, Union, Any
from sympy import (
    Matrix, Symbol, simplify, symbols, Rational,
    zeros, eye, diag, ImmutableMatrix, S, diff
)
from sympy.matrices import Matrix as MatrixType
import copy
import logging


def p0_initialize_old_variables(ncompo: int, variables: List[Union[Symbol, str]]) -> List[Symbol]:
    """
    p0: Initialize old-time variables (Maple p0).

    Args:
        ncompo: Number of components
        variables: List of variable symbols or names

    Returns:
        List of old variable Symbols (e.g., o2_old, no3_old, ...)
    """
    if ncompo <= 0:
        raise ValueError("ncompo must be positive")

    old_vars = []
    for i in range(ncompo):
        var = variables[i]
        name = str(var)
        old_vars.append(Symbol(f"{name}_old"))
    return old_vars


def p1_initialize_reaction_lists(nreactions: int,
                                 ncompo: int,
                                 variables: List[Union[Symbol, str]],
                                 variables_old: List[Union[Symbol, str]]) -> Dict[str, Any]:
    """
    p1: Initialize rate lists and substitutions (Maple p1).

    Returns:
        Dict with ratelist, v1, v2, substi1, substi2
    """
    if nreactions <= 0:
        raise ValueError("nreactions must be positive")
    if ncompo <= 0:
        raise ValueError("ncompo must be positive")

    ratelist = [Symbol(f"r{i}") for i in range(1, nreactions + 1)]
    v1 = [Symbol(f"sp({i},j)") for i in range(1, ncompo + 1)]
    v2 = [Symbol(f"spold({i},j)") for i in range(1, ncompo + 1)]

    substi1 = {Symbol(str(variables[i])): v1[i] for i in range(ncompo)}
    substi2 = {Symbol(str(variables_old[i])): v2[i] for i in range(ncompo)}

    return {
        "ratelist": ratelist,
        "v1": v1,
        "v2": v2,
        "substi1": substi1,
        "substi2": substi2,
    }


def p2_create_equation_names(variables: List[Union[Symbol, str]]) -> List[Symbol]:
    """
    p2: Create equation names (Maple p2).

    Returns:
        List of equation Symbols (e.g., do2dt, dno3dt, ...)
    """
    eqns = []
    for var in variables:
        name = str(var)
        eqns.append(Symbol(f"d{name}dt"))
    return eqns


def p3_reorder_reactions(eqrxn_ids: List[int], ratelist: List[Symbol]) -> List[Symbol]:
    """
    p3: Reorder reaction list so equilibrium reactions come first (Maple p3).

    CRITICAL: Maple's p3 uses set arithmetic and then convert(set, list) which
    produces a LEXICOGRAPHICALLY SORTED list for the non-equilibrium reactions.
    This means r1 < r10 < r11 < ... < r19 < r2 < r5 < r6 < ... (string sort).

    Args:
        eqrxn_ids: List of equilibrium reaction indices (1-based if >0, else 0-based)
        ratelist: List of reaction Symbols

    Returns:
        Reordered ratelist: equilibrium reactions first, then remaining sorted lexicographically
    """
    if not eqrxn_ids:
        return list(ratelist)

    zero_based = any(idx == 0 for idx in eqrxn_ids)
    eq_reactions = []
    used = set()
    for idx in eqrxn_ids:
        i = idx if zero_based else idx - 1
        if 0 <= i < len(ratelist):
            eq_reactions.append(ratelist[i])
            used.add(i)

    # Non-equilibrium reactions: sort LEXICOGRAPHICALLY (Maple set-to-list behavior)
    # Maple's convert(set, list) produces lexicographic order: r1 < r10 < r11 < ... < r2 < r5...
    non_eq_reactions = [r for i, r in enumerate(ratelist) if i not in used]
    non_eq_reactions_sorted = sorted(non_eq_reactions, key=lambda r: str(r))

    return eq_reactions + non_eq_reactions_sorted


class GaussianElimination:
    """
    Gaussian elimination with pivoting for stoichiometric matrix extraction.
    
    Implements the mathematical operations from Maple ACG p4-p9 functions.
    """
    
    def __init__(self, equations: List, variables: List[Symbol], 
                 reactions: List[Symbol], ncompo: int):
        """
        Initialize the Gaussian elimination processor.
        
        Args:
            equations: List of symbolic equations (equilibrium expressions)
            variables: List of SymPy symbols for active variables [C1, C2, ...]
            reactions: List of reaction rate symbols [r1, r2, ...]
            ncompo: Total number of components
            
        Raises:
            TypeError: If equations, variables, or reactions are not list-like
            ValueError: If dimensions don't match or values are invalid
        """
        # Validate input types
        if not isinstance(equations, (list, tuple)):
            raise TypeError(f"equations must be list or tuple, got {type(equations)}")
        if not isinstance(variables, (list, tuple)):
            raise TypeError(f"variables must be list or tuple, got {type(variables)}")
        if not isinstance(reactions, (list, tuple)):
            raise TypeError(f"reactions must be list or tuple, got {type(reactions)}")
        
        # Validate ncompo
        if not isinstance(ncompo, int) or ncompo <= 0:
            raise ValueError(f"ncompo must be positive integer, got {ncompo}")
        
        # Validate dimensions
        if len(equations) != ncompo:
            raise ValueError(f"equations length {len(equations)} != ncompo {ncompo}")
        if len(variables) != ncompo:
            raise ValueError(f"variables length {len(variables)} != ncompo {ncompo}")
        if len(reactions) == 0:
            raise ValueError("reactions list cannot be empty")
        
        self.equations = equations
        self.variables = variables
        self.reactions = reactions
        self.ncompo = ncompo
        self.nreactions = len(reactions)
        
        # State variables
        self.coeff_matrix = None          # p (from p4)
        self.augmented_matrix = None      # q (from p4)
        self.matrixM = None               # stoichiometric coefficients (from p6)
        self.matrixM_before_p8 = None     # stoichiometric matrix BEFORE p8 zeros conservation rows
        self.rightM = None                # reaction coefficients (from p6)
        self.func = None                  # net rates [f1, f2, ...] (from p7-p9)
        self.func_before_p8 = None        # func BEFORE p8 modifies it (for Jacobian of conservation rows)
        self.pd = None                    # jacobian matrix (from p10)
        self.pivot_rows = []              # rows used in elimination
        self.pivot_row_by_col = {}        # map reaction column -> pivot row
        self.conservation_rows = set()    # rows that are conservation laws
        self.conservation_basis = []      # nullspace basis vectors
        self.verbose = False
        self.logger = logging.getLogger(self.__class__.__name__)
        
    def p4_genmatrix(self) -> Matrix:
        """
        p4: Generate coefficient (stoichiometric) matrix from equations.
        
        CRITICAL: Extract coefficients of REACTIONS (r1, r2, ...), not variables!
        
        The equations are linear combinations of reaction rates:
            net_rate_i = sum_j (Q[i,j] * r_j)
        We extract the STOICHIOMETRIC MATRIX Q by taking coefficients w.r.t. the
        reaction list provided to the class.
        
        This matrix is constant (independent of concentrations) and is the
        correct input for Maple-style Gaussian elimination.
        
        Returns:
            Augmented matrix q = [coefficients | identity]
        """
        if self.verbose:
            print("p4: Building coefficient matrix from equations...")
        
        # Extract coefficients of each REACTION from each equation
        # This gives us the stoichiometric coefficient matrix Q
        coeff_list = []
        for eq in self.equations:
            row = []
            expanded_eq = eq.expand()
            for reaction in self.reactions:
                # Get coefficient of this reaction in this equation
                coeff = expanded_eq.coeff(reaction)
                if coeff is None or coeff == 0:
                    coeff = S(0)
                row.append(coeff)
            coeff_list.append(row)
        
        self.coeff_matrix = Matrix(coeff_list)
        if self.verbose:
            print(f"  Coefficient matrix: {self.coeff_matrix.shape}")
            
        # Detect conservation laws from coefficient matrix
        self._detect_conservation_laws_from_coefficients()
        
        # Augment with identity matrix for Gaussian elimination
        identity = eye(self.ncompo)
        q = self.coeff_matrix.row_join(identity)
        self.augmented_matrix = q
        
        if self.verbose:
            print(f"  Augmented matrix shape: {q.shape}")
        
        return q
    
    def p5_gaussian_elimination(self) -> Matrix:
        """
        p5: Gaussian elimination with ONLY ROW PIVOTING (Maple-style).
        
        Implements Maple's linalg[pivot]() behavior EXACTLY:
        - For each reaction column j:
          - Find first unused row i with non-zero entry q[i,j]
          - Call pivot(q, i, j) which:
            1. Does NOT normalize row i (pivot row stays unchanged)
            2. Eliminate column j in ALL other rows k ≠ i:
               q[k,:] = q[k,:] - (q[k,j] / q[i,j]) * q[i,:]
        
        CRITICAL: 
        - NO normalization of the pivot row (Maple linalg[pivot] does not normalize)
        - Eliminates from ALL rows (not just below), i.e. full Gaussian elimination
        - This preserves the original stoichiometric coefficients in pivot rows
        
        Returns: Modified augmented matrix with row operations applied
        """
        if self.verbose:
            print("p5: Gaussian elimination (Maple-style: NO normalization, full elimination)...")
        
        q = self.augmented_matrix.copy()
        self.pivot_rows = []
        self.pivot_row_by_col = {}
        
        # Convert to list of lists for manual manipulation (avoids SymPy overhead)
        q_list = [[q[i, j] for j in range(q.cols)] for i in range(q.rows)]
        
        # For each reaction column (Maple iterates j1 from 1 to nreactions)
        for col in range(self.nreactions):
            # Find first unused row with non-zero entry in this column
            # (Maple: for i1 from 1 to ncompo if q[i1,j1] <> 0 and i1 not in pivotrows)
            pivot_row = None
            for row in range(self.ncompo):
                if row not in self.pivot_rows and q_list[row][col] != 0:
                    pivot_row = row
                    break
            
            if pivot_row is None:
                if self.verbose:
                    print(f"  Column {col}: no pivot found (dependent column)")
                continue
            
            # Record pivot row
            self.pivot_rows.append(pivot_row)
            self.pivot_row_by_col[col] = pivot_row
            
            if self.verbose:
                print(f"  Column {col}: using row {pivot_row} as pivot")
            
            # MAPLE PIVOT: NO normalization — pivot row stays unchanged
            pivot_val = q_list[pivot_row][col]
            
            # Eliminate column in ALL other rows (both above and below pivot)
            # For all rows k ≠ pivot_row:
            #   q[k,:] = q[k,:] - (q[k,col] / pivot_val) * q[pivot_row,:]
            for row in range(self.ncompo):
                if row != pivot_row:
                    factor = q_list[row][col]
                    if factor != 0:
                        for j in range(q.cols):
                            q_list[row][j] = q_list[row][j] - (factor / pivot_val) * q_list[pivot_row][j]
        
        q = Matrix(q_list)
        self.augmented_matrix = q
        
        if self.verbose:
            print(f"  Pivot rows: {self.pivot_rows}")
            print(f"  Matrix transformed (NO column permutation!)")
        return q
    
    def p6_extract_matrices(self) -> Tuple[Matrix, Matrix]:
        """
        p6: Extract stoichiometric matrices from augmented matrix.
        
        After Gaussian elimination, the augmented matrix [A|I] becomes [I'|M].
        After RREF, the first nreactions columns form the identity for independent rows.
        The remaining ncompo columns form matrixM (stoichiometric matrix).
        
        rightM = I' (identity-like, shows which reactions are independent)
        matrixM = M (stoichiometric matrix, shows how components depend on reactions)
        
        Returns:
            (rightM, matrixM) = stoichiometric coefficient matrices
        """
        if self.verbose:
            print("p6: Extracting stoichiometric matrices...")
        
        q = self.augmented_matrix
        
        # Split augmented matrix: first nreactions cols = rightM, next ncompo cols = matrixM
        # The augmented matrix was created as [coefficients (ncompo×ncompo) | identity (ncompo×ncompo)]
        # So: cols 0:nreactions are part of the reduced coefficients
        #     cols nreactions:nreactions+ncompo are the stoichiometric matrix M (from pivoting the identity)
        self.rightM = q[:, :self.nreactions]
        self.matrixM = q[:, self.nreactions:self.nreactions+self.ncompo]
        
        if self.verbose:
            print(f"  rightM shape: {self.rightM.shape}")
            print(f"  matrixM shape: {self.matrixM.shape}")

        # Detect conservation rows from rightM (dependent rows after elimination)
        self._detect_conservation_rows_from_rightM()
        
        return self.rightM, self.matrixM
    
    def _detect_conservation_laws_from_coefficients(self):
        """
        Detect conservation laws from the original coefficient matrix A.
        
        Conservation laws correspond to the nullspace of A^T, which are
        linear combinations of components that don't change in any reaction.
        """
        if self.verbose:
            print("Detecting conservation laws from coefficient matrix...")
        
        # Compute nullspace of A^T (original coefficient matrix transpose)
        A_T = self.coeff_matrix.T
        ns_A_T = A_T.nullspace()
        
        if ns_A_T:
            self.conservation_basis = ns_A_T
            if self.verbose:
                print(f"  Found {len(ns_A_T)} conservation law(s)")
                for i, v in enumerate(ns_A_T):
                    print(f"    Conservation law {i+1}: {v.T}")

    def _detect_conservation_rows_from_rightM(self) -> None:
        """
        Identify conservation rows as dependent rows after elimination.

        A row is considered a conservation row if its rightM row is all zeros
        (i.e., the component equation is linearly dependent on others).
        """
        self.conservation_rows = set()
        for i in range(self.ncompo):
            is_zero_row = True
            for j in range(self.nreactions):
                if self.rightM[i, j] != 0:
                    is_zero_row = False
                    break
            if is_zero_row:
                self.conservation_rows.add(i)

        if self.verbose:
            print(f"  Conservation law rows: {self.conservation_rows}")
    
    def _is_conservation_row(self, row_idx: int) -> bool:
        """
        Check if a row index corresponds to a conservation law.
        
        Args:
            row_idx: Component row index
        
        Returns:
            True if the row is part of a conservation relationship
        """
        return row_idx in self.conservation_rows
    
    def p7_build_residuals(self) -> List:
        """
        p7: Build net rate expressions from rightM and reactions (Maple behavior).
        
        Maple computes:
            func[i] = sum_j rightM[i,j] * r_j
        
        This ensures the reduced system is expressed in terms of reaction rates,
        which can be symbolic (r1, r2, ...) or full expressions in C.
        Conservation rows are handled later in p8.
        
        Returns:
            List of symbolic expressions for func[1..ncompo]
        """
        if self.verbose:
            print("p7: Building net rate expressions from rightM * reactions...")

        # Maple-compatible behavior:
        # p7 uses rightM * rates1 after assignments r_i := rate_i.
        # That means the matrix is multiplied by ACTUAL rate expressions,
        # not by bare marker symbols.
        reaction_terms = getattr(self, 'reaction_values', self.reactions)

        # Build func = rightM * reaction_terms
        func_list = []
        for i in range(self.ncompo):
            expr = S(0)
            for j in range(self.nreactions):
                expr += self.rightM[i, j] * reaction_terms[j]
            func_list.append(simplify(expr))

        self.func = func_list
        self.func_after_p7 = list(func_list)

        if self.verbose:
            print(f"  Generated {len(self.func)} net rate expressions")
        return self.func
    
    def p8_handle_equilibrium(self, 
                              neqrxns: int = 0,
                              eq_columns: List[int] = None,
                              eq_target_rows: List[int] = None,
                              equilibrium_eqns: List = None,
                              inert_components: Set[int] = None,
                              variables_old: List[Symbol] = None,
                              delt: Symbol = None) -> List:
        """
        p8: Handle equilibrium reactions and conservation equations.
        
        Part 1 - Equilibrium reactions:
        - Replace rate expression with equilibrium equation
        - Zero out corresponding row in matrixM
        
        Part 2 - Conservation equations (inert components):
        - CRITICAL (THE p8 TRICK): Compute func[i] = matrixM[i,:] * (C - C_old) FIRST
        - THEN zero out matrixM row
        - This embeds the matrixM coefficients in func, allowing p10 to recover them
          as constant Jacobian entries!
        
        Args:
            neqrxns: Number of equilibrium reactions
            equilibrium_eqns: Expressions for equilibrium reactions
            inert_components: Set of inert component indices
            variables_old: Old concentration symbols (for conservation computation)
            delt: Time step symbol (for conservation computation)
        
        Returns:
            Modified func list
        """
        if self.verbose:
            print(f"p8: Handling equilibrium (neqrxns={neqrxns}) and conservation...")
        
        if equilibrium_eqns is None:
            equilibrium_eqns = []
        if eq_columns is None:
            eq_columns = list(range(neqrxns))
        if eq_target_rows is None:
            eq_target_rows = []
        if inert_components is None:
            inert_components = set()
        
        # Part 1: Handle equilibrium reactions
        # Maple semantics (see acg_optjoined_inpOK.md):
        # for j1 = 1..neqrxns, for all rows i with rightM[i,j1] <> 0:
        #   matrixM[i,:] := 0; func[i] := equilibriumeqns[j1]
        #
        # In Maple, p3 already moved equilibrium reactions to the first neqrxns
        # positions, so j runs over leading columns.
        for j in range(neqrxns):
            col = j
            if eq_columns is not None and j < len(eq_columns):
                col = eq_columns[j]

            if not (0 <= col < self.nreactions):
                continue

            for i in range(self.ncompo):
                if self.rightM[i, col] != 0:
                    if j < len(equilibrium_eqns):
                        self.func[i] = equilibrium_eqns[j]
                    for k in range(self.ncompo):
                        self.matrixM[i, k] = S(0)
                    if self.verbose:
                        print(f"  Component {i} uses equilibrium column {col}")
        
        # Part 2: Handle conservation equations (rows where rightM[i,:] = 0)
        # THE p8 TRICK: Store matrixM coefficients in func BEFORE zeroing matrixM
        for i in range(self.ncompo):
            # Check if this is a conservation row (no reactions affect it)
            is_conservation = True
            for j in range(self.nreactions):
                if self.rightM[i, j] != 0:
                    is_conservation = False
                    break
            
            if is_conservation or (i in inert_components):
                if self.verbose:
                    print(f"  Component {i} is conservation/inert")
                
                # THE CRITICAL FIX (p8 TRICK): 
                # Compute: func[i] = matrixM[i,:] @ (C - C_old)
                # This embeds the matrixM coefficients in func
                # When p10 computes jacobian, it will get: pd[i,j] = d(func[i])/d(C[j]) = matrixM[i,j]
                # This is a CONSTANT (independent of concentration)!
                
                if variables_old is not None:
                    conservation_sum = S(0)
                    for j in range(self.ncompo):
                        # func[i] = matrixM[i,:] @ (C - C_old)
                        conservation_sum += self.matrixM[i, j] * (self.variables[j] - variables_old[j])
                    self.func[i] = conservation_sum
                else:
                    # Fallback if variables_old not provided (shouldn't happen)
                    conservation_sum = S(0)
                    for j in range(self.ncompo):
                        conservation_sum += self.matrixM[i, j] * self.variables[j]
                    self.func[i] = conservation_sum
                
                # AFTER computing func[i] with matrixM embedded, THEN zero out matrixM row
                # This is the clever design: matrixM is "moved into" func, then cleared
                for j in range(self.ncompo):
                    self.matrixM[i, j] = S(0)
        
        if self.verbose:
            print(f"  Processed equilibrium and conservation rows")
        return self.func
    
    def p9_add_time_discretization(self,
                                   variables_old: List[Symbol],
                                   delt: Symbol) -> List:
        """
        p9: Add implicit Euler time discretization.
        
        Uses matrixM to compute time discretization terms:
        func[i] = func[i] - (matrixM @ C_new)[i]/delt + (matrixM @ C_old)[i]/delt
        
        matrixM encodes component coupling in the transformed system.
        For reactive rows, matrixM may have non-trivial entries that couple components.
        For conservation rows, matrixM entries are set to zero after p8.
        
        Args:
            variables_old: List of old concentration symbols [C1_old, C2_old, ...]
            delt: Time step symbol
        
        Returns:
            Modified func list with time discretization terms
        """
        if self.verbose:
            print("p9: Adding implicit Euler time discretization...")
        
        # Create concentration vectors
        c_vector = Matrix(self.variables)
        c_vector_old = Matrix(variables_old)
        
        # Compute time discretization terms using matrixM
        fd_new = self.matrixM * c_vector
        fd_old = self.matrixM * c_vector_old
        
        # Add time discretization to each component
        for i in range(self.ncompo):
            time_term = -(fd_new[i] / delt) + (fd_old[i] / delt)
            self.func[i] = self.func[i] + time_term
            self.func[i] = simplify(self.func[i])
        
        if self.verbose:
            print(f"  Added time discretization to {len(self.func)} equations")
        return self.func
    
    def p10_compute_jacobian(self) -> Matrix:
        """
        p10: Compute Jacobian matrix (symbolic differentiation).
        
        CRITICAL FIX for conservation rows:
        - For conservation rows: Use the conservation_basis vector directly
        - For reactive rows: pd[i,j] = ∂func[i]/∂C[j] (normal differentiation)
        
        The conservation_basis vectors give us the original conservation laws
        before Gaussian elimination scrambled the indices.
        
        For each conservation row i:
        1. Find the basis vector that has component i involved
        2. Use that basis vector's coefficients as the Jacobian row
        
        Returns:
            Jacobian matrix (ncompo x ncompo)
            
        Raises:
            RuntimeError: If func not initialized or required matrices missing
        """
        # Validate state before proceeding
        if self.func is None:
            raise RuntimeError(
                "func not initialized. "
                "Call p9_add_time_discretization() before p10_compute_jacobian()."
            )
        if self.rightM is None or self.matrixM is None:
            raise RuntimeError(
                "rightM or matrixM not initialized. "
                "Call p6_extract_matrices() first."
            )
        
        if self.verbose:
            print("p10: Computing Jacobian matrix...")
        
        jacobian_list = []
        for i in range(self.ncompo):
            row = []
            for j in range(self.ncompo):
                # Maple p10: pd := jacobian(func, v1)
                deriv = self.func[i].diff(self.variables[j])
                deriv = simplify(deriv)
                row.append(deriv)
            jacobian_list.append(row)
        
        self.pd = Matrix(jacobian_list)
        
        if self.verbose:
            print(f"  Jacobian shape: {self.pd.shape}")
        return self.pd
    
    def get_reduced_system(self) -> Dict:
        """
        Get the complete reduced system from p4-p10.
        
        Returns:
            Dict containing:
            - 'func': Residual equations
            - 'jacobian': Jacobian matrix
            - 'matrixM': Stoichiometric matrix (AFTER p8 modifications)
            - 'matrixM_before_p8': Stoichiometric matrix (BEFORE p8, for conservation rows)
            - 'rightM': Reaction coefficients
            - 'pivot_rows': Rows used in elimination
            - 'conservation_rows': Rows that are conservation laws
            - 'conservation_basis': Nullspace basis vectors
            - 'func_after_p7': Residual equations right after p7 (for debugging)
            - 'coeff_matrix': Coefficient matrix from p4 (for debugging)
        """
        result = {
            'func': self.func,
            'jacobian': self.pd,
            'matrixM': self.matrixM,
            'matrixM_before_p8': self.matrixM_before_p8,  # Direct access, not getattr
            'rightM': self.rightM,
            'pivot_rows': self.pivot_rows,
            'conservation_rows': self.conservation_rows,
            'conservation_basis': self.conservation_basis,
            'func_after_p7': self.func_after_p7,  # DEBUG
            'coeff_matrix': self.coeff_matrix,  # DEBUG
        }
        return result


def p4_genmatrix(ge: GaussianElimination) -> Matrix:
    """p4: Wrapper to generate coefficient matrix (Maple p4)."""
    return ge.p4_genmatrix()


def p5_gaussian_elimination(ge: GaussianElimination) -> Matrix:
    """p5: Wrapper to perform Gaussian elimination (Maple p5)."""
    return ge.p5_gaussian_elimination()


def p6_extract_matrices(ge: GaussianElimination) -> Tuple[Matrix, Matrix]:
    """p6: Wrapper to extract stoichiometric matrices (Maple p6)."""
    return ge.p6_extract_matrices()


def p7_build_residuals(ge: GaussianElimination) -> List:
    """p7: Wrapper to build residual expressions (Maple p7)."""
    return ge.p7_build_residuals()


def p8_handle_equilibrium(ge: GaussianElimination,
                          neqrxns: int = 0,
                          eq_columns: List[int] = None,
                          eq_target_rows: List[int] = None,
                          equilibrium_eqns: List = None,
                          inert_components: Set[int] = None,
                          variables_old: List[Symbol] = None,
                          delt: Symbol = None) -> List:
    """p8: Wrapper to handle equilibrium and conservation (Maple p8)."""
    return ge.p8_handle_equilibrium(
        neqrxns=neqrxns,
        eq_columns=eq_columns,
        eq_target_rows=eq_target_rows,
        equilibrium_eqns=equilibrium_eqns,
        inert_components=inert_components,
        variables_old=variables_old,
        delt=delt,
    )


def p9_add_time_discretization(ge: GaussianElimination,
                               variables_old: List[Symbol],
                               delt: Symbol) -> List:
    """p9: Wrapper to add implicit Euler time discretization (Maple p9)."""
    return ge.p9_add_time_discretization(variables_old=variables_old, delt=delt)


def p10_compute_jacobian(ge: GaussianElimination) -> Matrix:
    """p10: Wrapper to compute Jacobian (Maple p10)."""
    return ge.p10_compute_jacobian()


def run_gaussian_elimination(equations: List,
                             variables: List[Symbol],
                             reactions: List[Symbol],
                             variables_old: List[Symbol],
                             delt: Symbol,
                             reaction_values: List = None,
                             eqrxn_ids: List[int] = None,
                             eq_target_rows: List[int] = None,
                             neqrxns: int = 0,
                             equilibrium_eqns: List = None,
                             inert_components: Set[int] = None,
                             verbose: bool = False) -> Dict:
    """
    Execute complete Gaussian elimination pipeline (p4-p10).
    
    This is the high-level interface that runs all steps.
    
    Args:
        equations: Equilibrium conditions/equations
        variables: Current concentration symbols
        reactions: Reaction marker symbols (e.g., r1, r2, ...)
        reaction_values: Optional reaction expressions to evaluate in p7.
                         If None, uses `reactions` directly.
        variables_old: Previous time step concentration symbols
        delt: Time step symbol
        neqrxns: Number of equilibrium reactions
        equilibrium_eqns: Equilibrium expressions
        inert_components: Set of inert component indices
        verbose: Print debugging information
    
    Returns:
        Dict with reduced system information
    """
    ncompo = len(variables)

    # Preprocessing helpers (p0-p3) available for debugging/compatibility
    _ = p0_initialize_old_variables(ncompo, variables)
    _ = p1_initialize_reaction_lists(len(reactions), ncompo, variables, variables_old)
    _ = p2_create_equation_names(variables)
    base_reaction_values = reaction_values if reaction_values is not None else reactions

    # Maple p3: reorder reactions so equilibrium reactions come first.
    # CRITICAL: Use the SAME permutation for both reaction symbols and values!
    # Compute permutation from reaction symbols, apply to both.
    reordered_reactions = p3_reorder_reactions(eqrxn_ids or [], reactions)
    # Derive permutation: for each reordered reaction symbol, find its index in original list
    perm_indices = [reactions.index(r) for r in reordered_reactions]
    reordered_reaction_values = [base_reaction_values[i] for i in perm_indices]
    
    ge = GaussianElimination(equations, variables, reordered_reactions, ncompo)
    ge.reaction_values = reordered_reaction_values
    ge.verbose = verbose
    
    # Run all steps (p4-p10)
    p4_genmatrix(ge)
    p5_gaussian_elimination(ge)
    p6_extract_matrices(ge)
    p7_build_residuals(ge)
    
    # CRITICAL: Save matrixM BEFORE p8 modifies it!
    # For conservation rows, p8 will embed matrixM in func, then zero the matrix
    # But the Jacobian of conservation equations IS matrixM[i,j]!
    # We need to preserve the ORIGINAL matrixM values before p8 changes them
    ge.matrixM_before_p8 = Matrix(ge.matrixM)  # Make a deep copy before p8 modifies it!
    
    # p8 now receives variables_old and delt so it can properly implement the p8 TRICK
    p8_handle_equilibrium(
        ge,
        neqrxns,
        None,
        None,
        equilibrium_eqns,
        inert_components,
        variables_old,
        delt,
    )
    p9_add_time_discretization(ge, variables_old, delt)
    p10_compute_jacobian(ge)
    
    return ge.get_reduced_system()
