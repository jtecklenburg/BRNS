"""
ACG-BRNS - Automatic Code Generation for BRNS

Generates Fortran code for Biogeochemical Reaction Network Simulator from symbolic models.

This package provides:
- ACGModule: Main code generation interface
- Helpers for SymPy-to-Fortran conversions
- Gaussian elimination pipeline (p4-p10)
- Support for reaction networks with conservation laws

Quick Start:
    >>> from acg_brns import ACGModule
    >>> acg = ACGModule('output/')
    >>> acg.acg0(
    ...     nsolids=2, 
    ...     ndissolved=3, 
    ...     nreactions=4,
    ...     nnodes=10,
    ...     bio_names=['k_decay', 'k_oxid'],
    ... )

Documentation:
    - CODE_REVIEW_FIXES.md: Code quality improvements
    - UNIT_TESTS_ENHANCEMENT_SUMMARY.md: Test coverage
    - TEST_QUICK_REFERENCE.md: Running tests

Requirements:
    - sympy >= 1.12
    - macrofor >= 0.1.0 (for Fortran code generation)
"""

from .version import __version__, __version_info__
from .acg import ACGModule
from .helpers import create_substitution_dict, substitute_for_fortran

__all__ = [
    '__version__',
    '__version_info__',
    'ACGModule',
    'create_substitution_dict',
    'substitute_for_fortran',
]
