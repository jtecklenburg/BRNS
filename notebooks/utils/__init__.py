"""
Utils Package für BRNS Python Notebooks

Dieses Package enthält Hilfsmodule für die Code-Generierung und
Reaktionsnetzwerk-Verwaltung.
"""

from .acg_module import ACGModule, create_substitution_dict, substitute_for_fortran

__all__ = [
    'ACGModule',
    'create_substitution_dict',
    'substitute_for_fortran'
]

__version__ = '0.1.0'
