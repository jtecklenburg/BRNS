"""Tests for ACG-BRNS helper functions."""

import pytest
from sympy import symbols, Symbol
from acg_brns.helpers import create_substitution_dict, substitute_for_fortran


class TestCreateSubstitutionDict:
    """Test create_substitution_dict function."""

    def test_single_species(self):
        """Test with single species."""
        diss_a = symbols('diss_a')
        subst = create_substitution_dict([diss_a], 1)
        
        assert diss_a in subst
        assert str(subst[diss_a]) == 'sp(1,j)'

    def test_multiple_species(self):
        """Test with multiple species."""
        diss_a, diss_b, diss_c = symbols('diss_a diss_b diss_c')
        species = [diss_a, diss_b, diss_c]
        subst = create_substitution_dict(species, 3)
        
        assert str(subst[diss_a]) == 'sp(1,j)'
        assert str(subst[diss_b]) == 'sp(2,j)'
        assert str(subst[diss_c]) == 'sp(3,j)'

    def test_custom_index_var(self):
        """Test with custom index variable."""
        diss_a = symbols('diss_a')
        subst = create_substitution_dict([diss_a], 1, j_var='i')
        
        assert str(subst[diss_a]) == 'sp(1,i)'

    def test_empty_list(self):
        """Test with empty species list."""
        subst = create_substitution_dict([], 0)
        assert subst == {}


class TestSubstituteForFortran:
    """Test substitute_for_fortran function."""

    def test_simple_substitution(self):
        """Test simple substitution."""
        diss_a = symbols('diss_a')
        k_deg = symbols('k_deg')
        
        expr = k_deg * diss_a
        subst = create_substitution_dict([diss_a], 1)
        
        result = substitute_for_fortran(expr, subst)
        
        # Check that diss_a was substituted
        assert diss_a not in result.free_symbols
        # Check that k_deg remains
        assert k_deg in result.free_symbols

    def test_complex_expression(self):
        """Test complex expression substitution."""
        A, B = symbols('A B')
        k1, k2 = symbols('k1 k2')
        
        expr = k1 * A - k2 * B
        subst = create_substitution_dict([A, B], 2)
        
        result = substitute_for_fortran(expr, subst)
        
        # Original species should be gone
        assert A not in result.free_symbols
        assert B not in result.free_symbols
        # Parameters should remain
        assert k1 in result.free_symbols
        assert k2 in result.free_symbols

    def test_no_substitution_needed(self):
        """Test expression with no species."""
        k_deg = symbols('k_deg')
        expr = k_deg * 2
        
        result = substitute_for_fortran(expr, {})
        assert result == expr


class TestIntegrationHelpers:
    """Integration tests for helper functions."""

    def test_workflow_with_helpers(self):
        """Test complete workflow using helpers."""
        # Define symbolic model
        diss_a, diss_b = symbols('diss_a diss_b', real=True, positive=True)
        k1, k2 = symbols('k1 k2', real=True, positive=True)
        
        # Reaction: A -> B
        rate_a = -k1 * diss_a
        rate_b = k1 * diss_a - k2 * diss_b
        
        # Create substitution
        species = [diss_a, diss_b]
        subst = create_substitution_dict(species, 2)
        
        # Apply substitution
        rate_a_fortran = substitute_for_fortran(rate_a, subst)
        rate_b_fortran = substitute_for_fortran(rate_b, subst)
        
        # Check results
        assert str(subst[diss_a]) == 'sp(1,j)'
        assert str(subst[diss_b]) == 'sp(2,j)'
        assert k1 in rate_a_fortran.free_symbols
        assert k1 in rate_b_fortran.free_symbols
        assert k2 in rate_b_fortran.free_symbols
        assert diss_a not in rate_a_fortran.free_symbols
        assert diss_b not in rate_b_fortran.free_symbols
