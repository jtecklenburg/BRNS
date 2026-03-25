"""
Helper functions for ACG-BRNS.

Provides utility functions for SymPy-to-Fortran substitutions.
"""

from typing import List, Dict
from sympy import Symbol


def create_substitution_dict(
    species_symbols: List[Symbol],
    ncompo: int,
    j_var: str = "j"
) -> Dict[Symbol, Symbol]:
    """
    Create substitution dictionary for SymPy → Fortran conversion.
    
    Converts symbolic species to Fortran array notation.
    Example: Symbol('C') → sp(1,j)
    
    Args:
        species_symbols: List of SymPy symbols representing species
        ncompo: Number of components
        j_var: Fortran index variable name (default: 'j')
        
    Returns:
        Dictionary mapping SymPy symbols to Fortran array symbols
        
    Raises:
        ValueError: If species_symbols length != ncompo
        ValueError: If j_var is not a valid identifier
        TypeError: If species_symbols contains non-Symbol objects
        
    Example:
        >>> from sympy import symbols
        >>> diss_a, diss_b = symbols('diss_a diss_b')
        >>> subst = create_substitution_dict([diss_a, diss_b], 2)
        >>> print(subst[diss_a])
        sp(1,j)
    """
    # Validation
    if len(species_symbols) != ncompo:
        raise ValueError(
            f"species_symbols length ({len(species_symbols)}) "
            f"!= ncompo ({ncompo})"
        )
    
    if not j_var or not j_var.isidentifier():
        raise ValueError(
            f"j_var must be valid Python identifier, got '{j_var}'"
        )
    
    # Build substitution dictionary
    subst_dict = {}
    for i, sym in enumerate(species_symbols, 1):
        if not isinstance(sym, Symbol):
            raise TypeError(
                f"species_symbols[{i-1}] must be SymPy Symbol, "
                f"got {type(sym)}"
            )
        subst_dict[sym] = Symbol(f"sp({i},{j_var})")
    
    return subst_dict


def substitute_for_fortran(expr, subst_dict: Dict):
    """
    Apply substitutions to SymPy expression for Fortran generation.
    
    Replaces symbolic species with Fortran array notation.
    
    Args:
        expr: SymPy expression to substitute
        subst_dict: Substitution dictionary from create_substitution_dict()
        
    Returns:
        Substituted SymPy expression ready for Fortran code generation
        (or original expression if subst_dict is empty)
        
    Raises:
        TypeError: If expr is None
        
    Example:
        >>> from sympy import symbols
        >>> diss_a = symbols('diss_a')
        >>> k_deg = symbols('k_deg')
        >>> rate = k_deg * diss_a
        >>> subst = create_substitution_dict([diss_a], 1)
        >>> result = substitute_for_fortran(rate, subst)
        >>> print(result)
        k_deg*sp(1,j)
    """
    if expr is None:
        raise TypeError("expr cannot be None")
    
    # Empty substitution dict is OK - returns original expression
    if not subst_dict:
        return expr
    
    return expr.subs(subst_dict)
