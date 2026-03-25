"""
Unit tests for YAML to ACG Mapper

Tests the conversion of YAML configuration data into ACG-compatible structures.

Author: Jan Tecklenburg
Date: 2026-03-03
"""

import unittest
import numpy as np
import yaml
import tempfile
import os
from acg_brns.yaml_to_acg_mapper import YAMLtoACGMapper
from acg_brns.formula_evaluator import FormulaEvaluator


class TestBioArrayBuilding(unittest.TestCase):
    """Test bio_name and bio_val array construction."""
    
    def setUp(self):
        """Set up test config with biogeochemical parameters."""
        self.config = {
            'parameters': {
                'biogeochemical': [
                    {'name': 'kfox', 'value': 0.221},
                    {'name': 'kmo2', 'value': 8.0e-6},
                    {'name': 'ho2', 'value': 0.0},
                ]
            },
            'species': []
        }
        self.evaluated = {
            'kfox': 0.221,
            'kmo2': 8.0e-6,
            'ho2': 0.0
        }
        self.mapper = YAMLtoACGMapper(self.config, self.evaluated)
    
    def test_bio_arrays_extraction(self):
        """Test extraction of bio parameter names and values."""
        bio_name, bio_val = self.mapper.build_bio_arrays()
        
        self.assertEqual(len(bio_name), 3)
        self.assertEqual(len(bio_val), 3)
        self.assertEqual(bio_name, ['kfox', 'kmo2', 'ho2'])
        self.assertEqual(bio_val, [0.221, 8.0e-6, 0.0])
    
    def test_bio_arrays_types(self):
        """Test that bio_val contains floats."""
        bio_name, bio_val = self.mapper.build_bio_arrays()
        
        for val in bio_val:
            self.assertIsInstance(val, float)
    
    def test_bio_arrays_empty(self):
        """Test handling of empty biogeochemical parameters."""
        config = {'species': []}
        mapper = YAMLtoACGMapper(config, {})
        bio_name, bio_val = mapper.build_bio_arrays()
        
        self.assertEqual(len(bio_name), 0)
        self.assertEqual(len(bio_val), 0)
    
    def test_bio_arrays_with_unevaluated(self):
        """Formula-based parameters must not silently fall back to 0.0."""
        config = {
            'parameters': {
                'biogeochemical': [
                    {'name': 'param1', 'value': 1.5},
                    {'name': 'param2', 'value': 'unevaluable'},
                ]
            },
            'species': []
        }
        evaluated = {'param1': 1.5}  # param2 not evaluated
        mapper = YAMLtoACGMapper(config, evaluated)

        with self.assertRaises(ValueError) as ctx:
            mapper.build_bio_arrays()
        self.assertIn("param2", str(ctx.exception))


class TestVariablesList(unittest.TestCase):
    """Test species ordering in variables list."""
    
    def setUp(self):
        """Set up test config with dissolved and solid species."""
        self.config = {
            'species': [
                {'name': 'o2', 'type': 'dissolved'},
                {'name': 'ch2o', 'type': 'solid'},
                {'name': 'no3', 'type': 'dissolved'},
                {'name': 'mno2', 'type': 'solid'},
                {'name': 'hco3', 'type': 'dissolved'},
            ]
        }
        self.mapper = YAMLtoACGMapper(self.config, {})
    
    def test_variables_ordering(self):
        """Test that species order matches YAML definition order.

        The mapper preserves the YAML species order (model-specific ordering).
        For the Canfield model solid species are interleaved among dissolved ones
        (e.g. ch2o at index 7, mno2 at 11), so no dissolved-first reordering
        is performed.
        """
        variables = self.mapper.build_variables_list()

        self.assertEqual(len(variables), 5)
        # Order must match YAML definition order exactly
        self.assertEqual(variables, ['o2', 'ch2o', 'no3', 'mno2', 'hco3'])
    
    def test_variables_completeness(self):
        """Test that all species are in variables list."""
        variables = self.mapper.build_variables_list()
        
        all_species = ['o2', 'ch2o', 'no3', 'mno2', 'hco3']
        for species in all_species:
            self.assertIn(species, variables)
    
    def test_variables_no_duplicates(self):
        """Test that variables list has no duplicates."""
        variables = self.mapper.build_variables_list()
        self.assertEqual(len(variables), len(set(variables)))
    
    def test_variables_only_dissolved(self):
        """Test handling of config with only dissolved species."""
        config = {
            'species': [
                {'name': 'o2', 'type': 'dissolved'},
                {'name': 'no3', 'type': 'dissolved'},
            ]
        }
        mapper = YAMLtoACGMapper(config, {})
        variables = mapper.build_variables_list()
        
        self.assertEqual(variables, ['o2', 'no3'])
    
    def test_variables_only_solid(self):
        """Test handling of config with only solid species."""
        config = {
            'species': [
                {'name': 'ch2o', 'type': 'solid'},
                {'name': 'mno2', 'type': 'solid'},
            ]
        }
        mapper = YAMLtoACGMapper(config, {})
        variables = mapper.build_variables_list()
        
        self.assertEqual(variables, ['ch2o', 'mno2'])


class TestSpeciesIndexMap(unittest.TestCase):
    """Test species name to Fortran index mapping."""
    
    def setUp(self):
        """Set up test config."""
        self.config = {
            'species': [
                {'name': 'o2', 'type': 'dissolved'},
                {'name': 'no3', 'type': 'dissolved'},
                {'name': 'ch2o', 'type': 'solid'},
            ]
        }
        self.mapper = YAMLtoACGMapper(self.config, {})
    
    def test_index_map_one_indexed(self):
        """Test that indices start at 1 (Fortran convention)."""
        species_map = self.mapper.build_species_index_map()
        
        self.assertEqual(species_map['o2'], 1)
        self.assertEqual(species_map['no3'], 2)
        self.assertEqual(species_map['ch2o'], 3)
    
    def test_index_map_completeness(self):
        """Test that all species are in the map."""
        species_map = self.mapper.build_species_index_map()
        
        self.assertEqual(len(species_map), 3)
        self.assertIn('o2', species_map)
        self.assertIn('no3', species_map)
        self.assertIn('ch2o', species_map)
    
    def test_index_map_consistent_with_variables(self):
        """Test that index map is consistent with variables list."""
        variables = self.mapper.build_variables_list()
        species_map = self.mapper.build_species_index_map()
        
        for idx, name in enumerate(variables):
            self.assertEqual(species_map[name], idx + 1)
    
    def test_get_solid_indices(self):
        """Test extraction of solid species indices."""
        solid_indices = self.mapper.get_solid_indices()
        
        self.assertEqual(len(solid_indices), 1)
        self.assertEqual(solid_indices[0], 3)  # ch2o is 3rd species


class TestStoichiometryMatrix(unittest.TestCase):
    """Test stoichiometry matrix construction."""
    
    def setUp(self):
        """Set up test config with reactions and species."""
        self.config = {
            'species': [
                {'name': 'o2', 'type': 'dissolved'},
                {'name': 'no3', 'type': 'dissolved'},
                {'name': 'ch2o', 'type': 'solid'},
            ],
            'reactions': [
                {
                    'id': 1,
                    'name': 'aerobic_respiration',
                    'stoichiometry': {
                        'o2': -1.0,
                        'ch2o': -1.0,
                        'no3': 0.0,
                    }
                },
                {
                    'id': 2,
                    'name': 'denitrification',
                    'stoichiometry': {
                        'no3': -0.8,
                        'ch2o': -1.0,
                        'o2': 0.0,
                    }
                }
            ]
        }
        self.mapper = YAMLtoACGMapper(self.config, {})
    
    def test_stoich_matrix_dimensions(self):
        """Test that matrix has correct dimensions."""
        stoich = self.mapper.build_stoichiometry_matrix()
        
        self.assertEqual(stoich.shape, (2, 3))  # 2 reactions, 3 species
    
    def test_stoich_matrix_values(self):
        """Test that matrix has correct stoichiometric coefficients."""
        stoich = self.mapper.build_stoichiometry_matrix()
        
        # Reaction 1: o2=-1, no3=0, ch2o=-1
        self.assertEqual(stoich[0, 0], -1.0)  # o2
        self.assertEqual(stoich[0, 1], 0.0)   # no3
        self.assertEqual(stoich[0, 2], -1.0)  # ch2o
        
        # Reaction 2: o2=0, no3=-0.8, ch2o=-1
        self.assertEqual(stoich[1, 0], 0.0)   # o2
        self.assertEqual(stoich[1, 1], -0.8)  # no3
        self.assertEqual(stoich[1, 2], -1.0)  # ch2o
    
    def test_stoich_matrix_with_evaluator(self):
        """Test matrix building with formula evaluator."""
        # Add formula-based stoichiometry
        config = {
            'species': [
                {'name': 'o2', 'type': 'dissolved'},
                {'name': 'ch2o', 'type': 'solid'},
            ],
            'reactions': [
                {
                    'id': 1,
                    'name': 'test_reaction',
                    'stoichiometry': {
                        'o2': '-1',
                        'ch2o': '-2/3',
                    }
                }
            ],
            'computed_values': {}
        }
        
        evaluator = FormulaEvaluator(config)
        mapper = YAMLtoACGMapper(config, {})
        stoich = mapper.build_stoichiometry_matrix(evaluator)
        
        self.assertEqual(stoich.shape, (1, 2))
        self.assertAlmostEqual(stoich[0, 0], -1.0)
        self.assertAlmostEqual(stoich[0, 1], -2.0/3.0, places=10)


class TestSpeciesSubstitution(unittest.TestCase):
    """Test species name substitution in expressions."""
    
    def setUp(self):
        """Set up test config."""
        self.config = {
            'species': [
                {'name': 'o2', 'type': 'dissolved'},
                {'name': 'no3', 'type': 'dissolved'},
                {'name': 'ch2o', 'type': 'solid'},
            ]
        }
        self.mapper = YAMLtoACGMapper(self.config, {})
    
    def test_simple_substitution(self):
        """Test substitution of single species."""
        expr = 'o2/kmo2'
        result = self.mapper.substitute_species_in_expression(expr)
        
        self.assertEqual(result, 'sp(1,j)/kmo2')
    
    def test_multiple_substitutions(self):
        """Test substitution of multiple species."""
        expr = 'o2*no3/ch2o'
        result = self.mapper.substitute_species_in_expression(expr)
        
        self.assertEqual(result, 'sp(1,j)*sp(2,j)/sp(3,j)')
    
    def test_word_boundary_matching(self):
        """Test that word boundaries are respected."""
        # Add species with overlapping names
        config = {
            'species': [
                {'name': 'o2', 'type': 'dissolved'},
                {'name': 'kmo2', 'type': 'dissolved'},  # Contains 'o2'
            ]
        }
        mapper = YAMLtoACGMapper(config, {})
        
        expr = 'o2/kmo2'
        result = mapper.substitute_species_in_expression(expr)
        
        # Should substitute both correctly, not create 'sp(1,j)/kmsp(1,j)'
        self.assertEqual(result, 'sp(1,j)/sp(2,j)')
    
    def test_overlapping_species_names(self):
        """Test handling of species with overlapping names."""
        config = {
            'species': [
                {'name': 'h2s', 'type': 'dissolved'},
                {'name': 'hs', 'type': 'dissolved'},  # 'hs' is contained in 'h2s'
            ]
        }
        mapper = YAMLtoACGMapper(config, {})
        
        expr = 'h2s + hs'
        result = mapper.substitute_species_in_expression(expr)
        
        # Both should be substituted correctly
        self.assertEqual(result, 'sp(1,j) + sp(2,j)')
    
    def test_parameters_not_substituted(self):
        """Test that parameter names are not substituted."""
        expr = 'o2/(kmo2 + o2)'
        result = self.mapper.substitute_species_in_expression(expr)
        
        # Only 'o2' should be substituted, not 'kmo2'
        self.assertEqual(result, 'sp(1,j)/(kmo2 + sp(1,j))')


class TestRateComponentExpansion(unittest.TestCase):
    """Test expansion of rate components."""
    
    def setUp(self):
        """Set up test mapper."""
        config = {'species': []}
        self.mapper = YAMLtoACGMapper(config, {})
    
    def test_simple_expansion(self):
        """Test expansion of simple component."""
        rate = 'kch2o * fo2'
        components = {
            'kch2o': 'kfox * ch2o',
            'fo2': '1.0'
        }
        
        result = self.mapper.expand_rate_components(rate, components)
        
        self.assertIn('(kfox * ch2o)', result)
        self.assertIn('(1.0)', result)
    
    def test_nested_expansion(self):
        """Test expansion of nested components."""
        rate = 'kch2o * fo2'  # Direct rate, not referencing components
        components = {
            'kch2o': 'kfox * ch2o',
            'fo2': 'ho2 + (1.0 - ho2)*(o2/kmo2)'
        }
        
        result = self.mapper.expand_rate_components(rate, components)
        
        # Should contain expanded kch2o and fo2
        self.assertIn('(kfox * ch2o)', result)
        self.assertIn('ho2', result)
        self.assertIn('kmo2', result)
    
    def test_no_components(self):
        """Test handling of rate without components."""
        rate = 'kfox * o2'
        components = {}
        
        result = self.mapper.expand_rate_components(rate, components)
        
        self.assertEqual(result, 'kfox * o2')
    
    def test_numeric_rate(self):
        """Test handling of numeric rate value."""
        rate = 0
        components = {}
        
        result = self.mapper.expand_rate_components(rate, components)
        
        self.assertEqual(result, '0')


class TestRateExpressionBuilding(unittest.TestCase):
    """Test complete rate expression building."""
    
    def setUp(self):
        """Set up test config with reactions."""
        self.config = {
            'species': [
                {'name': 'o2', 'type': 'dissolved'},
                {'name': 'co2', 'type': 'dissolved'},
                {'name': 'hco3', 'type': 'dissolved'},
                {'name': 'ch2o', 'type': 'solid'},
            ],
            'reactions': [
                {
                    'id': 1,
                    'name': 'aerobic_respiration',
                    'rate': 'kch2o * fo2',
                    'rate_components': {
                        'kch2o': 'kfox * ch2o',
                        'fo2': 'ho2 + (1.0 - ho2)*(o2/kmo2)'
                    },
                    'equilibrium': False
                },
                {
                    'id': 2,
                    'name': 'carbonate_eq',
                    'rate': 0,
                    'equilibrium': True,
                    'equilibrium_constraint': 'hco3 - k1 * co2'
                }
            ]
        }
        self.mapper = YAMLtoACGMapper(self.config, {})
    
    def test_rate_expressions_count(self):
        """Test that all reactions are processed."""
        reactions = self.mapper.build_rate_expressions()
        
        self.assertEqual(len(reactions), 2)
    
    def test_kinetic_reaction(self):
        """Test processing of kinetic reaction."""
        reactions = self.mapper.build_rate_expressions()
        rxn = reactions[0]
        
        self.assertEqual(rxn['id'], 1)
        self.assertEqual(rxn['name'], 'aerobic_respiration')
        self.assertEqual(rxn['rate_yaml'], 'kch2o * fo2')
        self.assertFalse(rxn['equilibrium'])
        
        # Check expanded rate
        self.assertIn('kfox', rxn['rate_expanded'])
        self.assertIn('ch2o', rxn['rate_expanded'])
        self.assertNotIn('fo2', rxn['rate_expanded'])
        
        # Check Fortran rate with sp(i,j)
        self.assertIn('sp(', rxn['rate_fortran'])
        self.assertIn('sp(1,j)', rxn['rate_fortran'])  # o2
        self.assertIn('sp(4,j)', rxn['rate_fortran'])  # ch2o (now 4th species)
        self.assertNotIn('fo2', rxn['rate_fortran'])
    
    def test_equilibrium_reaction(self):
        """Test processing of equilibrium reaction."""
        reactions = self.mapper.build_rate_expressions()
        rxn = reactions[1]
        
        self.assertEqual(rxn['id'], 2)
        self.assertTrue(rxn['equilibrium'])
        
        # Check constraint with sp(i,j)
        # hco3 is species 3, co2 is species 2
        self.assertIn('sp(3,j)', rxn['equilibrium_constraint'])  # hco3
        self.assertIn('sp(2,j)', rxn['equilibrium_constraint'])  # co2

    def test_reactions_are_sorted_by_id(self):
        """Reaction order should follow Maple IDs, not YAML declaration order."""
        unsorted_config = {
            'species': self.config['species'],
            'reactions': [
                {
                    'id': 6,
                    'name': 'later_reaction',
                    'rate': 'foo',
                    'equilibrium': False,
                },
                {
                    'id': 2,
                    'name': 'earlier_reaction',
                    'rate': 'bar',
                    'equilibrium': False,
                },
            ],
        }
        mapper = YAMLtoACGMapper(unsorted_config, {})

        reactions = mapper.build_rate_expressions()

        self.assertEqual([reaction['id'] for reaction in reactions], [2, 6])


class TestGlobalRateComponentExpansion(unittest.TestCase):
    """Test Maple-like global helper expansion across reactions."""

    def test_cross_reaction_components_are_expanded(self):
        config = {
            'species': [
                {'name': 'o2', 'type': 'dissolved'},
                {'name': 'no3', 'type': 'dissolved'},
                {'name': 'ch2o', 'type': 'solid'},
            ],
            'parameters': {
                'biogeochemical': [
                    {'name': 'kfox', 'value': 1.0},
                    {'name': 'kmo2', 'value': 1.0},
                    {'name': 'hno3', 'value': 1.0},
                    {'name': 'hno3f2', 'value': 1.0},
                    {'name': 'kmno3', 'value': 1.0},
                ]
            },
            'reactions': [
                {
                    'id': 1,
                    'name': 'r1',
                    'rate': 'kch2o * fo2',
                    'rate_components': {
                        'kch2o': 'kfox * ch2o',
                        'fo2': 'o2/kmo2',
                    },
                },
                {
                    'id': 2,
                    'name': 'r2',
                    'rate': 'kch2o * fno3',
                    'rate_components': {
                        'fno3': 'hno3*(1-fo2)*(hno3f2 + (1.0 - hno3f2)*(no3/kmno3))',
                    },
                },
            ],
        }

        mapper = YAMLtoACGMapper(config, {})
        reactions = mapper.build_rate_expressions()

        self.assertIn('kfox', reactions[1]['rate_expanded'])
        self.assertIn('ch2o', reactions[1]['rate_expanded'])
        self.assertIn('fo2', reactions[0]['rate_yaml'])   # fo2 is in r1's rate, not r2's
        self.assertNotIn('kch2o', reactions[1]['rate_expanded'])

    def test_biogeochemical_parameters_are_not_inlined(self):
        config = {
            'species': [
                {'name': 'mn2', 'type': 'dissolved'},
                {'name': 'co3', 'type': 'dissolved'},
                {'name': 'mnco3', 'type': 'solid'},
            ],
            'parameters': {
                'biogeochemical': [
                    {'name': 'sw17', 'value': 1},
                    {'name': 'k17_1', 'value': 1.0},
                    {'name': 'k17_2', 'value': 1.0},
                    {'name': 'K_mnco3', 'value': 1.0},
                ]
            },
            'reactions': [
                {
                    'id': 17,
                    'name': 'mnco3_precipitation',
                    'rate': '(k17_1*sw17 + k17_2*mnco3*(1-sw17))*(omega_mn-1)',
                    'rate_components': {
                        'omega_mn': 'mn2*co3/K_mnco3',
                        'sw17': '1',
                    },
                }
            ],
        }

        mapper = YAMLtoACGMapper(config, {})
        reactions = mapper.build_rate_expressions()

        self.assertIn('sw17', reactions[0]['rate_expanded'])
        self.assertNotIn('(1)', reactions[0]['rate_expanded'])


class TestIntegrationWithCanfield(unittest.TestCase):
    """Integration tests with full Canfield model."""
    
    @classmethod
    def setUpClass(cls):
        """Load Canfield YAML and create mapper."""
        yaml_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'models',
            'canfield_refactored.yaml'
        )
        
        with open(yaml_path, 'r') as f:
            cls.config = yaml.safe_load(f)
        
        # Evaluate formulas
        cls.evaluator = FormulaEvaluator(cls.config)
        cls.evaluated = cls.evaluator.evaluate_all()
        
        # Create mapper
        cls.mapper = YAMLtoACGMapper(cls.config, cls.evaluated)
    
    def test_canfield_species_count(self):
        """Test that Canfield model has correct species count."""
        variables = self.mapper.build_variables_list()
        
        # 12 dissolved + 6 solid = 18 total
        self.assertEqual(len(variables), 18)
        
        dissolved = [s for s in self.mapper.species_list if s.get('type') == 'dissolved']
        solid = [s for s in self.mapper.species_list if s.get('type') == 'solid']
        
        self.assertEqual(len(dissolved), 12)
        self.assertEqual(len(solid), 6)
    
    def test_canfield_reaction_count(self):
        """Test that Canfield model has correct reaction count."""
        reactions = self.mapper.build_rate_expressions()
        
        # 16 kinetic + 3 equilibrium = 19 total
        self.assertEqual(len(reactions), 19)
        
        equilibrium_rxns = [r for r in reactions if r.get('equilibrium')]
        kinetic_rxns = [r for r in reactions if not r.get('equilibrium')]
        
        self.assertEqual(len(equilibrium_rxns), 3)
        self.assertEqual(len(kinetic_rxns), 16)
    
    def test_canfield_stoichiometry_matrix(self):
        """Test Canfield stoichiometry matrix dimensions."""
        stoich = self.mapper.build_stoichiometry_matrix(self.evaluator)
        
        self.assertEqual(stoich.shape, (19, 18))  # 19 reactions, 18 species
        
        # Check that matrix has reasonable number of non-zero entries
        non_zero = np.count_nonzero(stoich)
        self.assertGreater(non_zero, 50)  # Should have many non-zero coefficients
        self.assertLess(non_zero, 200)    # But not all entries
    
    def test_canfield_bio_arrays(self):
        """Test Canfield biogeochemical parameters."""
        bio_name, bio_val = self.mapper.build_bio_arrays()
        
        # Should have many biogeochemical parameters
        self.assertGreater(len(bio_name), 30)
        
        # Check that all are floats
        for val in bio_val:
            self.assertIsInstance(val, float)
        
        # Check for specific expected parameters
        self.assertIn('kfox', bio_name)
        self.assertIn('kmo2', bio_name)
        self.assertIn('kmno3', bio_name)
    
    def test_canfield_rate_expressions(self):
        """Test that all Canfield reactions have valid rate expressions."""
        reactions = self.mapper.build_rate_expressions()
        
        for rxn in reactions:
            # All reactions should have Fortran rate
            self.assertIn('rate_fortran', rxn)
            
            # Kinetic reactions should have non-empty rate
            if not rxn.get('equilibrium'):
                self.assertTrue(len(rxn['rate_fortran']) > 0)
            
            # Equilibrium reactions should have constraint
            if rxn.get('equilibrium'):
                self.assertIn('equilibrium_constraint', rxn)
                self.assertTrue(len(rxn['equilibrium_constraint']) > 0)


class TestBoundaryConditions(unittest.TestCase):
    """Test boundary condition extraction."""
    
    def setUp(self):
        """Set up test config."""
        self.config = {
            'species': [
                {
                    'name': 'o2',
                    'type': 'dissolved',
                    'bc_upper_type': 0,
                    'bc_upper_value': 0.200,
                    'bc_lower_type': 1,
                    'bc_lower_value': 0.0
                },
                {
                    'name': 'ch2o',
                    'type': 'solid',
                    'bc_upper_type': 2,
                    'bc_upper_value': 100.0,
                    'bc_lower_type': 1,
                    'bc_lower_value': 0.0
                }
            ]
        }
        self.mapper = YAMLtoACGMapper(self.config, {})
    
    def test_boundary_conditions_extraction(self):
        """Test extraction of boundary conditions."""
        bc_data = self.mapper.get_boundary_conditions()
        
        self.assertEqual(len(bc_data), 2)
        
        # Check o2
        self.assertIn('o2', bc_data)
        self.assertEqual(bc_data['o2']['bc_upper_type'], 0)
        self.assertEqual(bc_data['o2']['bc_upper_value'], 0.200)
        self.assertEqual(bc_data['o2']['bc_lower_type'], 1)
        self.assertEqual(bc_data['o2']['bc_lower_value'], 0.0)
        
        # Check ch2o
        self.assertIn('ch2o', bc_data)
        self.assertEqual(bc_data['ch2o']['bc_upper_type'], 2)
        self.assertEqual(bc_data['ch2o']['bc_upper_value'], 100.0)


class TestInitialConditions(unittest.TestCase):
    """Test initial condition extraction."""
    
    def setUp(self):
        """Set up test config."""
        self.config = {
            'species': [
                {
                    'name': 'o2',
                    'type': 'dissolved',
                    'init_value': 0.200,
                    'init_mode': 3
                },
                {
                    'name': 'ch2o',
                    'type': 'solid',
                    'init_value': 1000.0,
                    'init_mode': 3
                }
            ]
        }
        self.mapper = YAMLtoACGMapper(self.config, {})
    
    def test_initial_conditions_extraction(self):
        """Test extraction of initial conditions."""
        init_data = self.mapper.get_initial_conditions()
        
        self.assertEqual(len(init_data), 2)
        
        # Check o2
        self.assertIn('o2', init_data)
        self.assertEqual(init_data['o2']['init_value'], 0.200)
        self.assertEqual(init_data['o2']['init_mode'], 3)
        
        # Check ch2o
        self.assertIn('ch2o', init_data)
        self.assertEqual(init_data['ch2o']['init_value'], 1000.0)
        self.assertEqual(init_data['ch2o']['init_mode'], 3)


if __name__ == '__main__':
    unittest.main()
