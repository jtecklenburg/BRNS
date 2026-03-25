"""
Unit Tests for FormulaEvaluator

Tests formula evaluation, dependency resolution, and scientific notation handling.

Author: Jan Tecklenburg
Date: 2026-03-03
"""

import pytest
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from acg_brns.formula_evaluator import (
    FormulaEvaluator,
    DependencyResolver,
    FormulaEvaluationError,
    CircularDependencyError,
    UndefinedVariableError
)


class TestDependencyResolver:
    """Test dependency extraction and topological sorting"""
    
    def test_extract_dependencies_simple(self):
        """Test extracting dependencies from a simple formula"""
        resolver = DependencyResolver()
        known = {'x', 'y', 'z'}
        
        deps = resolver.extract_dependencies("x + y", known)
        assert deps == {'x', 'y'}
        
        deps = resolver.extract_dependencies("2 * z", known)
        assert deps == {'z'}
    
    def test_extract_dependencies_complex(self):
        """Test extracting dependencies from complex formulas"""
        resolver = DependencyResolver()
        known = {'mlogkp1', 'S', 'T', 'por0', 's_dens'}
        
        # Scientific notation formula
        deps = resolver.extract_dependencies(
            "-13.7201 + 0.031334*T + 3235.67/T + 1.3e-5*S*T - 0.1032*sqrt(S)",
            known
        )
        assert deps == {'T', 'S'}
        
        # Formula with division
        deps = resolver.extract_dependencies(
            "s_dens*(1-por0)/por0*1000",
            known
        )
        assert deps == {'s_dens', 'por0'}
        
        # Power formula
        deps = resolver.extract_dependencies(
            "10^(-mlogkp1)",
            known
        )
        assert deps == {'mlogkp1'}
    
    def test_extract_dependencies_no_deps(self):
        """Test formula with no dependencies (numeric literal)"""
        resolver = DependencyResolver()
        known = {'x', 'y'}
        
        deps = resolver.extract_dependencies("42.5", known)
        assert deps == set()
        
        deps = resolver.extract_dependencies(123, known)
        assert deps == set()
    
    def test_topological_sort_linear(self):
        """Test topological sort with linear dependency chain"""
        resolver = DependencyResolver()
        
        dep_graph = {
            'a': set(),
            'b': {'a'},
            'c': {'b'},
            'd': {'c'}
        }
        
        result = resolver.topological_sort(dep_graph)
        
        # Check that dependencies come before dependents
        assert result.index('a') < result.index('b')
        assert result.index('b') < result.index('c')
        assert result.index('c') < result.index('d')
    
    def test_topological_sort_diamond(self):
        """Test topological sort with diamond-shaped dependencies"""
        resolver = DependencyResolver()
        
        dep_graph = {
            'a': set(),
            'b': {'a'},
            'c': {'a'},
            'd': {'b', 'c'}
        }
        
        result = resolver.topological_sort(dep_graph)
        
        # a must come first
        assert result.index('a') < result.index('b')
        assert result.index('a') < result.index('c')
        # b and c must come before d
        assert result.index('b') < result.index('d')
        assert result.index('c') < result.index('d')
    
    def test_topological_sort_circular_dependency(self):
        """Test that circular dependencies are detected"""
        resolver = DependencyResolver()
        
        dep_graph = {
            'a': {'b'},
            'b': {'c'},
            'c': {'a'}
        }
        
        with pytest.raises(CircularDependencyError):
            resolver.topological_sort(dep_graph)


class TestFormulaEvaluation:
    """Test formula evaluation with various expression types"""
    
    def test_simple_arithmetic(self):
        """Test simple arithmetic expressions"""
        config = {
            'parameters': {
                'physical': {
                    'x': 10.0,
                    'y': 5.0
                }
            }
        }
        
        evaluator = FormulaEvaluator(config)
        evaluator.load_base_parameters()
        
        result = evaluator.evaluate_formula("x + y", evaluator.evaluated)
        assert result == 15.0
        
        result = evaluator.evaluate_formula("x * y", evaluator.evaluated)
        assert result == 50.0
        
        result = evaluator.evaluate_formula("x / y", evaluator.evaluated)
        assert result == 2.0
    
    def test_power_notation(self):
        """Test power notation with ^ operator"""
        config = {
            'parameters': {
                'physical': {
                    'base': 10.0,
                    'exp': 3.0
                }
            }
        }
        
        evaluator = FormulaEvaluator(config)
        evaluator.load_base_parameters()
        
        result = evaluator.evaluate_formula("base^exp", evaluator.evaluated)
        assert result == 1000.0
        
        result = evaluator.evaluate_formula("10^(-2)", evaluator.evaluated)
        assert result == 0.01
    
    def test_scientific_notation_decimal_exponent(self):
        """Test scientific notation with decimal exponents (e.g., 1.0e-8.5)"""
        config = {'parameters': {}}
        
        evaluator = FormulaEvaluator(config)
        
        # Test 1.0e-8.5 = 1.0 * 10^(-8.5)
        result = evaluator.evaluate_formula("1.0e-8.5", {})
        expected = 1.0 * (10 ** -8.5)
        assert abs(result - expected) < 1e-15
        
        # Test 2.5e-2.2
        result = evaluator.evaluate_formula("2.5e-2.2", {})
        expected = 2.5 * (10 ** -2.2)
        assert abs(result - expected) < 1e-15
        
        # Test positive decimal exponent
        result = evaluator.evaluate_formula("1.5e2.3", {})
        expected = 1.5 * (10 ** 2.3)
        assert abs(result - expected) < 1e-10
    
    def test_math_functions(self):
        """Test mathematical functions (sqrt, log, exp)"""
        config = {'parameters': {}}
        
        evaluator = FormulaEvaluator(config)
        
        result = evaluator.evaluate_formula("sqrt(25)", {})
        assert result == 5.0
        
        result = evaluator.evaluate_formula("log10(100)", {})
        assert abs(result - 2.0) < 1e-10

        result = evaluator.evaluate_formula("log10(10 * (4 + 6))", {})
        assert abs(result - 2.0) < 1e-10
        
        result = evaluator.evaluate_formula("exp(0)", {})
        assert result == 1.0
    
    def test_complex_formula_canfield_style(self):
        """Test complex formula similar to Canfield model"""
        config = {
            'parameters': {
                'physical': {
                    't_celsius': 10.3,
                    'salin': 35.0
                }
            }
        }
        
        evaluator = FormulaEvaluator(config)
        evaluator.load_base_parameters()
        
        # Temperature conversion should create T
        assert 'T' in evaluator.evaluated
        assert abs(evaluator.evaluated['T'] - 283.45) < 1e-10
        
        # Evaluate mlogkp1 formula from Canfield
        formula = "-13.7201 + 0.031334*T + 3235.67/T + 1.3e-5*S*T - 0.1032*sqrt(S)"
        result = evaluator.evaluate_formula(formula, evaluator.evaluated)
        
        # Result should be positive (this is -log(K), so K is small)
        assert result > 0
        assert result < 10  # Reasonable range for pK values
    
    def test_undefined_variable_detection(self):
        """Test that undefined variables are detected"""
        config = {'parameters': {}}
        
        evaluator = FormulaEvaluator(config)
        evaluator.load_base_parameters()
        
        with pytest.raises(FormulaEvaluationError) as excinfo:
            evaluator.evaluate_formula("undefined_var * 2", {})
        
        assert "undefined variables" in str(excinfo.value).lower()


class TestFormulaEvaluatorIntegration:
    """Integration tests for complete evaluation pipeline"""
    
    def test_simple_parameter_evaluation(self):
        """Test evaluation of simple parameter dependencies"""
        config = {
            'parameters': {
                'physical': {
                    'por0': 0.85,
                    's_dens': 2.5
                },
                'biogeochemical': [
                    {'name': 'x', 'value': 200.0},
                    {'name': 'y', 'value': 21.0},
                    {'name': 'SD', 'value': 's_dens*(1-por0)/por0*1000'}
                ]
            }
        }
        
        evaluator = FormulaEvaluator(config)
        result = evaluator.evaluate_all()
        
        # Check base parameters
        assert result['por0'] == 0.85
        assert result['s_dens'] == 2.5
        assert result['x'] == 200.0
        
        # Check computed parameter
        expected_SD = 2.5 * (1 - 0.85) / 0.85 * 1000
        assert abs(result['SD'] - expected_SD) < 1e-6
    
    def test_chain_of_dependencies(self):
        """Test evaluation with chain of dependencies"""
        config = {
            'parameters': {
                'physical': {
                    't_celsius': 10.3,
                    'salin': 35.0
                },
                'biogeochemical': [
                    {'name': 'keq1', 'value': '10^(-mlogkp1)'},
                    {'name': 'kbcarb', 'value': 'keq1'}
                ]
            },
            'computed_values': {
                'carbonate_system': {
                    'mlogkp1': '-13.7201 + 0.031334*T + 3235.67/T + 1.3e-5*S*T - 0.1032*sqrt(S)',
                    'k1': '10^(-mlogkp1)'
                }
            }
        }
        
        evaluator = FormulaEvaluator(config)
        result = evaluator.evaluate_all()
        
        # Check temperature conversion
        assert abs(result['T'] - 283.45) < 1e-10
        
        # Check mlogkp1 was evaluated
        assert 'mlogkp1' in result
        assert result['mlogkp1'] > 0
        
        # Check k1 and keq1 are equal (both use same formula)
        assert abs(result['k1'] - result['keq1']) < 1e-15
        
        # Check kbcarb equals keq1
        assert result['kbcarb'] == result['keq1']
    
    def test_scientific_notation_in_yaml(self):
        """Test that scientific notation in YAML is properly handled"""
        config = {
            'parameters': {
                'biogeochemical': [
                    {'name': 'K_mnco3', 'value': '1.0e-8.5'},
                    {'name': 'K_feco3', 'value': '1.0e-8.4'},
                    {'name': 'K_fes', 'value': '1.0e-2.2'}
                ]
            }
        }
        
        evaluator = FormulaEvaluator(config)
        result = evaluator.evaluate_all()
        
        # Check values are correctly evaluated
        assert abs(result['K_mnco3'] - (10 ** -8.5)) < 1e-15
        assert abs(result['K_feco3'] - (10 ** -8.4)) < 1e-15
        assert abs(result['K_fes'] - (10 ** -2.2)) < 1e-15


class TestStoichiometryEvaluation:
    """Test evaluation of stoichiometry coefficients"""
    
    def test_simple_stoichiometry(self):
        """Test evaluation of simple stoichiometry expressions"""
        config = {
            'parameters': {
                'physical': {
                    'por0': 0.85
                },
                'biogeochemical': [
                    {'name': 'SD', 'value': '441.18'},
                    {'name': 'x', 'value': 200.0},
                    {'name': 'y', 'value': 21.0},
                    {'name': 'z', 'value': 0.0}
                ]
            }
        }
        
        evaluator = FormulaEvaluator(config)
        evaluator.evaluate_all()
        
        # Test stoichiometry from reaction 1
        stoich = {
            'ch2o': -1,
            'o2': '- SD*(x+2*y)/x',
            'nh4': '(y/x) * SD'
        }
        
        result = evaluator.evaluate_stoichiometry(stoich)
        
        assert result['ch2o'] == -1.0
        
        # o2: - SD*(x+2*y)/x = - 441.18 * (200+42)/200
        expected_o2 = -441.18 * (200 + 2*21) / 200
        assert abs(result['o2'] - expected_o2) < 0.01
        
        # nh4: (y/x) * SD = (21/200) * 441.18
        expected_nh4 = (21 / 200) * 441.18
        assert abs(result['nh4'] - expected_nh4) < 0.01
    
    def test_stoichiometry_with_division(self):
        """Test stoichiometry with division operators"""
        config = {
            'parameters': {
                'biogeochemical': [
                    {'name': 'SD', 'value': 441.18}
                ]
            }
        }
        
        evaluator = FormulaEvaluator(config)
        evaluator.evaluate_all()
        
        stoich = {
            'mno2': '1/SD',
            'o2': '-1/2'
        }
        
        result = evaluator.evaluate_stoichiometry(stoich)
        
        assert abs(result['mno2'] - (1 / 441.18)) < 1e-6
        assert result['o2'] == -0.5


class TestSpeciesValues:
    """Test evaluation of species initial and boundary values"""
    
    def test_numeric_species_values(self):
        """Test species with numeric values"""
        config = {
            'species': [
                {
                    'name': 'o2',
                    'bc_upper_value': 132.0e-6,
                    'bc_lower_value': 0.0,
                    'init_value': 180.0e-6
                }
            ]
        }
        
        evaluator = FormulaEvaluator(config)
        result = evaluator.get_species_values(config['species'])
        
        assert result['o2']['bc_upper_value'] == 132.0e-6
        assert result['o2']['bc_lower_value'] == 0.0
        assert result['o2']['init_value'] == 180.0e-6
    
    def test_computed_species_values(self):
        """Test species with computed values"""
        config = {
            'parameters': {
                'physical': {'t_celsius': 10.3, 'salin': 35.0}
            },
            'computed_values': {
                'carbonate_system': {
                    'mlogkp1': '-13.7201 + 0.031334*T + 3235.67/T',
                    'k1': '10^(-mlogkp1)'
                }
            },
            'species': [
                {
                    'name': 'test_species',
                    'bc_upper_value': 'k1',
                    'bc_lower_value': 0.0,
                    'init_value': 'k1'
                }
            ]
        }
        
        evaluator = FormulaEvaluator(config)
        result = evaluator.get_species_values(config['species'])
        
        # k1 should be evaluated
        assert isinstance(result['test_species']['bc_upper_value'], float)
        assert result['test_species']['bc_upper_value'] > 0
        assert result['test_species']['bc_upper_value'] < 1e-5
        
        # Should be same for init_value
        assert result['test_species']['init_value'] == result['test_species']['bc_upper_value']
    
    def test_unevaluated_species_values(self):
        """Test that unevaluated species values remain as strings"""
        config = {
            'species': [
                {
                    'name': 'hco3',
                    'bc_upper_value': 'C_hco3',
                    'bc_lower_value': 0.0,
                    'init_value': 'C_hco3'
                }
            ]
        }
        
        evaluator = FormulaEvaluator(config)
        result = evaluator.get_species_values(config['species'])
        
        # C_hco3 cannot be evaluated without alkalinity
        assert result['hco3']['bc_upper_value'] == 'C_hco3'
        assert result['hco3']['init_value'] == 'C_hco3'
        assert result['hco3']['bc_lower_value'] == 0.0


class TestErrorHandling:
    """Test error handling and warning messages"""
    
    def test_circular_dependency_error(self):
        """Test that circular dependencies raise an error"""
        config = {
            'parameters': {
                'biogeochemical': [
                    {'name': 'a', 'value': 'b'},
                    {'name': 'b', 'value': 'c'},
                    {'name': 'c', 'value': 'a'}
                ]
            }
        }
        
        evaluator = FormulaEvaluator(config)
        
        with pytest.raises(CircularDependencyError):
            evaluator.evaluate_all()
    
    def test_unevaluable_formula_warning(self, capsys):
        """Verbose evaluators should emit warnings for unevaluable optional formulas."""
        config = {
            'parameters': {
                'physical': {'t_celsius': 10.3}
            },
            'computed_values': {
                'carbonate_system': {
                    'C_hco3': 'CT * alpha1'  # Cannot evaluate without CT and alpha1
                }
            }
        }

        evaluator = FormulaEvaluator(config, verbose=True)
        evaluator.evaluate_all()

        captured = capsys.readouterr()
        assert 'Warning' in captured.out
        assert 'C_hco3' in captured.out
        assert 'undefined variables' in captured.out.lower()

    def test_unevaluable_formula_silent_by_default(self, capsys):
        """Default evaluator mode should not print warnings to stdout."""
        config = {
            'parameters': {
                'physical': {'t_celsius': 10.3}
            },
            'computed_values': {
                'carbonate_system': {
                    'C_hco3': 'CT * alpha1'
                }
            }
        }

        evaluator = FormulaEvaluator(config)
        evaluator.evaluate_all()

        captured = capsys.readouterr()
        assert captured.out == ''

    def test_required_biogeochemical_parameter_must_evaluate(self):
        """Required biogeochemical parameters must fail fast if dependencies stay undefined."""
        config = {
            'parameters': {
                'biogeochemical': [
                    {'name': 'k_required', 'value': 'C_hco3'}
                ]
            },
            'computed_values': {
                'carbonate_system': {
                    'C_hco3': 'CT * alpha1'
                }
            }
        }

        evaluator = FormulaEvaluator(config)

        with pytest.raises(FormulaEvaluationError) as excinfo:
            evaluator.evaluate_all()
        assert "k_required" in str(excinfo.value)
        assert "could not be evaluated" in str(excinfo.value)
    
    def test_invalid_formula_syntax(self):
        """Test that invalid formula syntax raises an error"""
        config = {'parameters': {}}
        
        evaluator = FormulaEvaluator(config)
        
        with pytest.raises(FormulaEvaluationError):
            evaluator.evaluate_formula("invalid syntax )(", {})


class TestRealWorldCanfieldModel:
    """Integration tests with realistic Canfield model data"""
    
    def test_canfield_parameters(self):
        """Test evaluation of key Canfield model parameters"""
        config = {
            'parameters': {
                'physical': {
                    't_celsius': 10.3,
                    'salin': 35.0,
                    'por0': 0.85
                },
                'biogeochemical': [
                    {'name': 's_dens', 'value': 2.5},
                    {'name': 'x', 'value': 200.0},
                    {'name': 'y', 'value': 21.0},
                    {'name': 'SD', 'value': 's_dens*(1-por0)/por0*1000'},
                    {'name': 'keq1', 'value': '10^(-mlogkp1)'},
                    {'name': 'keq2', 'value': '10^(-mlogkp2)'},
                    {'name': 'kbcarb', 'value': 'keq1'},
                    {'name': 'K_mnco3', 'value': '1.0e-8.5'},
                    {'name': 'K_fes', 'value': '1.0e-2.2'}
                ]
            },
            'computed_values': {
                'carbonate_system': {
                    'mlogkp1': '-13.7201 + 0.031334*T + 3235.67/T + 1.3e-5*S*T - 0.1032*sqrt(S)',
                    'mlogkp2': '5371.9645 + 1.671221*T + 0.22913*S + 18.3802*log10(S) - 128375.28/T - 2194.3055*log10(T) - 8.0944e-4*S*T - 5617.11*log10(S)/T + 2.136*S/T',
                    'k1': '10^(-mlogkp1)',
                    'k2': '10^(-mlogkp2)'
                },
                'sulfide_system': {
                    'mlogkps': '2.527 + 1359.96/T - 0.206*S^(1/3)'
                }
            }
        }
        
        evaluator = FormulaEvaluator(config)
        result = evaluator.evaluate_all()
        
        # Check temperature conversion
        assert abs(result['T'] - 283.45) < 1e-10
        
        # Check SD calculation
        expected_SD = 2.5 * (1 - 0.85) / 0.85 * 1000
        assert abs(result['SD'] - expected_SD) < 1e-6
        
        # Check equilibrium constants are positive and reasonable
        assert result['mlogkp1'] > 0
        assert result['mlogkp2'] > 0
        assert 0 < result['keq1'] < 1e-5  # Small equilibrium constant
        assert 0 < result['keq2'] < 1e-9  # Even smaller
        
        # Check kbcarb equals keq1
        assert result['kbcarb'] == result['keq1']
        
        # Check scientific notation values
        assert abs(result['K_mnco3'] - (10 ** -8.5)) < 1e-15
        assert abs(result['K_fes'] - (10 ** -2.2)) < 1e-15
        
        # Check sulfide system
        assert 'mlogkps' in result
        assert result['mlogkps'] > 0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
