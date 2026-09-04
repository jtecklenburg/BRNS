"""
Unit tests for ACG Orchestrator

Tests the complete YAML→Fortran pipeline orchestration.

Author: Jan Tecklenburg
Date: 2026-03-03
"""

import unittest
import tempfile
import shutil
from pathlib import Path
import yaml
from unittest.mock import patch

from sympy import Symbol, eye, zeros

from acg_brns.acg_orchestrator import ACGOrchestrator, ACGOrchestrationError, ValidationIssue


class TestConfigLoading(unittest.TestCase):
    """Test YAML configuration loading and validation."""
    
    def setUp(self):
        """Create temporary directory for test files."""
        self.temp_dir = tempfile.mkdtemp()
        self.output_dir = Path(self.temp_dir) / 'output'
    
    def tearDown(self):
        """Clean up temporary directory."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_load_valid_config(self):
        """Test loading valid YAML configuration."""
        # Create minimal valid YAML
        yaml_path = Path(self.temp_dir) / 'test.yaml'
        config = {
            'species': [
                {'name': 'o2', 'type': 'dissolved'}
            ],
            'reactions': [
                {'id': 1, 'name': 'test', 'stoichiometry': {'o2': -1}}
            ],
            'parameters': {
                'biogeochemical': []
            }
        }
        with open(yaml_path, 'w') as f:
            yaml.dump(config, f)
        
        # Load config
        orchestrator = ACGOrchestrator(str(yaml_path), str(self.output_dir))
        loaded = orchestrator.load_config()
        
        self.assertIsNotNone(loaded)
        self.assertIn('species', loaded)
        self.assertIn('reactions', loaded)
    
    def test_missing_yaml_file(self):
        """Test error handling for missing YAML file."""
        orchestrator = ACGOrchestrator('nonexistent.yaml', str(self.output_dir))
        
        with self.assertRaises(ACGOrchestrationError):
            orchestrator.load_config()
    
    def test_invalid_yaml_format(self):
        """Test error handling for invalid YAML."""
        yaml_path = Path(self.temp_dir) / 'invalid.yaml'
        with open(yaml_path, 'w') as f:
            f.write("invalid: yaml: content: [")
        
        orchestrator = ACGOrchestrator(str(yaml_path), str(self.output_dir))
        
        with self.assertRaises(ACGOrchestrationError):
            orchestrator.load_config()

    def test_empty_yaml_file(self):
        """Empty YAML should fail with a clear validation error."""
        yaml_path = Path(self.temp_dir) / 'empty.yaml'
        yaml_path.write_text('', encoding='utf-8')

        orchestrator = ACGOrchestrator(str(yaml_path), str(self.output_dir))

        with self.assertRaises(ACGOrchestrationError) as ctx:
            orchestrator.load_config()
        self.assertIn('file is empty', str(ctx.exception))

    def test_top_level_yaml_must_be_mapping(self):
        """Top-level YAML sequences/scalars should fail cleanly."""
        yaml_path = Path(self.temp_dir) / 'list_root.yaml'
        yaml_path.write_text('- just\n- a\n- list\n', encoding='utf-8')

        orchestrator = ACGOrchestrator(str(yaml_path), str(self.output_dir))

        with self.assertRaises(ACGOrchestrationError) as ctx:
            orchestrator.load_config()
        self.assertIn('top-level document must be a mapping/object', str(ctx.exception))

    def test_parameters_must_be_mapping(self):
        """Top-level parameters section must be a mapping."""
        yaml_path = Path(self.temp_dir) / 'bad_parameters.yaml'
        config = {
            'species': [{'name': 'o2', 'type': 'dissolved'}],
            'reactions': [{'id': 1, 'name': 'test', 'stoichiometry': {'o2': -1}}],
            'parameters': []
        }
        with open(yaml_path, 'w') as f:
            yaml.dump(config, f)

        orchestrator = ACGOrchestrator(str(yaml_path), str(self.output_dir))

        with self.assertRaises(ACGOrchestrationError) as ctx:
            orchestrator.load_config()
        self.assertIn("'parameters' must be a mapping", str(ctx.exception))

    def test_nested_parameter_sections_must_have_expected_types(self):
        """Malformed nested parameter sections should yield validation errors, not crashes."""
        yaml_path = Path(self.temp_dir) / 'bad_parameter_sections.yaml'
        config = {
            'species': [{'name': 'o2', 'type': 'dissolved'}],
            'reactions': [{'id': 1, 'name': 'test', 'stoichiometry': {'o2': -1}}],
            'parameters': {
                'biogeochemical': {},
                'physical': [],
                'stoichiometry': 'oops',
            }
        }
        with open(yaml_path, 'w') as f:
            yaml.dump(config, f)

        orchestrator = ACGOrchestrator(str(yaml_path), str(self.output_dir))

        with self.assertRaises(ACGOrchestrationError) as ctx:
            orchestrator.load_config()
        msg = str(ctx.exception)
        self.assertIn("parameters.stoichiometry", msg)
        self.assertIn("must be either a mapping or a list", msg)

    def test_parameter_formula_cycle_across_subsections_is_detected(self):
        """Cycles spanning multiple parameters subsections should be reported."""
        yaml_path = Path(self.temp_dir) / 'cycle_across_subsections.yaml'
        config = {
            'species': [{'name': 'o2', 'type': 'dissolved'}],
            'reactions': [{'id': 1, 'name': 'test', 'stoichiometry': {'o2': -1}}],
            'parameters': {
                'physical': {
                    'a': 'b + 1'
                },
                'custom_section': [
                    {'name': 'b', 'value': 'a + 1'}
                ]
            }
        }
        with open(yaml_path, 'w') as f:
            yaml.dump(config, f)

        orchestrator = ACGOrchestrator(str(yaml_path), str(self.output_dir))

        with self.assertRaises(ACGOrchestrationError) as ctx:
            orchestrator.load_config()

        msg = str(ctx.exception)
        self.assertIn("Circular dependency detected in parameter formulas", msg)

    def test_optional_sections_must_have_expected_types(self):
        """Malformed optional top-level sections should yield validation errors."""
        yaml_path = Path(self.temp_dir) / 'bad_optional_sections.yaml'
        config = {
            'species': [{'name': 'o2', 'type': 'dissolved'}],
            'reactions': [{'id': 1, 'name': 'test', 'stoichiometry': {'o2': -1}}],
            'parameters': {'biogeochemical': []},
            'grid': [],
            'time': 'oops',
            'output': [],
            'initial_conditions': 1,
            'advanced': False,
            'transport': 3,
        }
        with open(yaml_path, 'w') as f:
            yaml.dump(config, f)

        orchestrator = ACGOrchestrator(str(yaml_path), str(self.output_dir))

        with self.assertRaises(ACGOrchestrationError) as ctx:
            orchestrator.load_config()
        msg = str(ctx.exception)
        self.assertIn("'grid' must be a mapping", msg)
        self.assertIn("'time' must be a mapping", msg)
        self.assertIn("'output' must be a mapping", msg)
        self.assertIn("'initial_conditions' must be a mapping", msg)
        self.assertIn("'advanced' must be a mapping", msg)
        self.assertIn("'transport' must be a mapping", msg)
    
    def test_missing_required_section(self):
        """Test validation fails for missing required sections."""
        yaml_path = Path(self.temp_dir) / 'incomplete.yaml'
        config = {
            'species': [{'name': 'o2', 'type': 'dissolved'}]
            # Missing 'reactions' and 'parameters'
        }
        with open(yaml_path, 'w') as f:
            yaml.dump(config, f)
        
        orchestrator = ACGOrchestrator(str(yaml_path), str(self.output_dir))
        
        with self.assertRaises(ACGOrchestrationError):
            orchestrator.load_config()

    def test_species_transport_defaults_to_true(self):
        """species.transport is optional and defaults to True."""
        yaml_path = Path(self.temp_dir) / 'species_transport_default.yaml'
        config = {
            'species': [
                {'name': 'o2', 'type': 'dissolved'},
                {'name': 'ch2o', 'type': 'solid'}
            ],
            'reactions': [{'id': 1, 'name': 'test', 'stoichiometry': {'o2': -1}}],
            'parameters': {'biogeochemical': []}
        }
        with open(yaml_path, 'w') as f:
            yaml.dump(config, f)

        orchestrator = ACGOrchestrator(str(yaml_path), str(self.output_dir))
        orchestrator.load_config()
        orchestrator.evaluate_formulas()
        orchestrator.map_to_acg_structures()

        self.assertEqual(orchestrator._build_no_transport_list(), [])

    def test_build_no_transport_list_from_species_transport(self):
        """No-transport list should be derived from species.transport flags."""
        yaml_path = Path(self.temp_dir) / 'species_transport_list.yaml'
        config = {
            'species': [
                {'name': 'o2', 'type': 'dissolved', 'transport': True},
                {'name': 'ch2o', 'type': 'solid', 'transport': False},
                {'name': 'fe2', 'type': 'dissolved', 'transport': False},
            ],
            'reactions': [{'id': 1, 'name': 'test', 'stoichiometry': {'o2': -1}}],
            'parameters': {'biogeochemical': []}
        }
        with open(yaml_path, 'w') as f:
            yaml.dump(config, f)

        orchestrator = ACGOrchestrator(str(yaml_path), str(self.output_dir))
        orchestrator.load_config()
        orchestrator.evaluate_formulas()
        orchestrator.map_to_acg_structures()

        self.assertEqual(orchestrator._build_no_transport_list(), [2, 3])
    
    def test_invalid_species_format(self):
        """Test validation fails for invalid species format."""
        yaml_path = Path(self.temp_dir) / 'badspecies.yaml'
        config = {
            'species': [
                {'name': 'o2'}  # Missing 'type'
            ],
            'reactions': [],
            'parameters': {}
        }
        with open(yaml_path, 'w') as f:
            yaml.dump(config, f)
        
        orchestrator = ACGOrchestrator(str(yaml_path), str(self.output_dir))
        
        with self.assertRaises(ACGOrchestrationError):
            orchestrator.load_config()

    def test_required_biogeochemical_formula_must_evaluate(self):
        """Required biogeochemical formulas must fail during formula evaluation, not map to 0.0."""
        yaml_path = Path(self.temp_dir) / 'bad_biogeochemical_formula.yaml'
        config = {
            'species': [{'name': 'o2', 'type': 'dissolved'}],
            'reactions': [{'id': 1, 'name': 'test', 'stoichiometry': {'o2': -1}}],
            'computed_values': {
                'carbonate_system': {
                    'C_hco3': 'CT * alpha1'
                }
            },
            'parameters': {
                'biogeochemical': [{'name': 'k_required', 'value': 'C_hco3'}]
            }
        }
        with open(yaml_path, 'w') as f:
            yaml.dump(config, f)

        orchestrator = ACGOrchestrator(str(yaml_path), str(self.output_dir))
        orchestrator.load_config()

        with self.assertRaises(ACGOrchestrationError) as ctx:
            orchestrator.evaluate_formulas()
        self.assertIn('Formula evaluation failed', str(ctx.exception))
        self.assertIn('k_required', str(ctx.exception))


class TestPhaseExecution(unittest.TestCase):
    """Test individual pipeline phases."""
    
    def setUp(self):
        """Set up test with the equilibrium model."""
        self.yaml_path = 'models/equilibrium/equilibrium.yaml'
        self.temp_dir = tempfile.mkdtemp()
        self.output_dir = Path(self.temp_dir) / 'output'
        
        # Check if the equilibrium model YAML exists
        if not Path(self.yaml_path).exists():
            self.skipTest(f"The equilibrium model YAML not found: {self.yaml_path}")
        
        self.orchestrator = ACGOrchestrator(
            self.yaml_path,
            str(self.output_dir),
            verbose=False
        )
    
    def tearDown(self):
        """Clean up temporary directory."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_load_config_phase(self):
        """Test Phase 1: Config loading."""
        config = self.orchestrator.load_config()
        
        self.assertIsNotNone(config)
        self.assertIn('species', config)
        self.assertIn('reactions', config)
        self.assertGreater(len(config['species']), 0)
        self.assertGreater(len(config['reactions']), 0)
    
    def test_evaluate_formulas_phase(self):
        """Test Phase 2: Formula evaluation."""
        self.orchestrator.load_config()
        params = self.orchestrator.evaluate_formulas()
        
        self.assertIsNotNone(params)
        self.assertIsInstance(params, dict)
        self.assertGreater(len(params), 0)
    
    def test_map_to_acg_phase(self):
        """Test Phase 3: Mapping to ACG structures."""
        self.orchestrator.load_config()
        self.orchestrator.evaluate_formulas()
        acg_data = self.orchestrator.map_to_acg_structures()
        
        self.assertIsNotNone(acg_data)
        self.assertIn('variables', acg_data)
        self.assertIn('bio_name', acg_data)
        self.assertIn('bio_val', acg_data)
        self.assertIn('stoich_matrix', acg_data)
        self.assertIn('reactions', acg_data)
        
        # Check dimensions
        self.assertEqual(len(acg_data['bio_name']), len(acg_data['bio_val']))
        self.assertEqual(acg_data['ncompo'], len(acg_data['variables']))
        self.assertEqual(acg_data['nreactions'], len(acg_data['reactions']))
    
    def test_preprocessing_phase(self):
        """Test Phase 4: Pre-processing (p0-p10)."""
        self.orchestrator.load_config()
        self.orchestrator.evaluate_formulas()
        self.orchestrator.map_to_acg_structures()
        reduced_system = self.orchestrator.run_preprocessing()
        
        self.assertIsNotNone(reduced_system)
        self.assertIn('func', reduced_system)
        self.assertIn('jacobian', reduced_system)
        self.assertIn('matrixM', reduced_system)
        self.assertIn('rightM', reduced_system)
        
        # Check dimensions
        ncompo = self.orchestrator.acg_data['ncompo']
        self.assertEqual(len(reduced_system['func']), ncompo)
        self.assertEqual(reduced_system['jacobian'].shape, (ncompo, ncompo))


class TestParameterExtraction(unittest.TestCase):
    """Test extraction of parameters from YAML."""
    
    def setUp(self):
        """Set up test with the equilibrium model."""
        self.yaml_path = 'models/equilibrium/equilibrium.yaml'
        self.temp_dir = tempfile.mkdtemp()
        self.output_dir = Path(self.temp_dir) / 'output'
        
        if not Path(self.yaml_path).exists():
            self.skipTest(f"Equilibrium model YAML not found: {self.yaml_path}")
        
        self.orchestrator = ACGOrchestrator(
            self.yaml_path,
            str(self.output_dir),
            verbose=False
        )
        
        # Run first 3 phases
        self.orchestrator.load_config()
        self.orchestrator.evaluate_formulas()
        self.orchestrator.map_to_acg_structures()
    
    def tearDown(self):
        """Clean up."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_extract_boundary_conditions(self):
        """Test extraction of boundary conditions."""
        params = self.orchestrator._extract_generation_parameters()
        
        self.assertIn('type_up', params)
        self.assertIn('bnddata_up', params)
        self.assertIn('type_down', params)
        self.assertIn('bnddata_down', params)
        
        ncompo = self.orchestrator.acg_data['ncompo']
        self.assertEqual(len(params['type_up']), ncompo)
        self.assertEqual(len(params['bnddata_up']), ncompo)
        self.assertEqual(len(params['type_down']), ncompo)
        self.assertEqual(len(params['bnddata_down']), ncompo)
    
    def test_extract_transport_parameters(self):
        """Test extraction of transport parameters."""
        params = self.orchestrator._extract_generation_parameters()
        
        self.assertIn('diffdata', params)
        self.assertIn('alphadata', params)
        
        ncompo = self.orchestrator.acg_data['ncompo']
        self.assertEqual(len(params['diffdata']), ncompo)
        self.assertEqual(len(params['alphadata']), ncompo)
    
    def test_extract_initial_conditions(self):
        """Test extraction of initial conditions."""
        params = self.orchestrator._extract_generation_parameters()
        
        self.assertIn('vic', params)
        self.assertIn('iniconc', params)
        
        ncompo = self.orchestrator.acg_data['ncompo']
        self.assertEqual(len(params['iniconc']), ncompo)


class TestIntegrationEquilibriumModel(unittest.TestCase):
    """Integration test with the full equilibrium model."""
    
    def setUp(self):
        """Set up test."""
        self.yaml_path = 'models/equilibrium/equilibrium.yaml'
        self.temp_dir = tempfile.mkdtemp()
        self.output_dir = Path(self.temp_dir) / 'output'
        
        if not Path(self.yaml_path).exists():
            self.skipTest(f"Equilibrium model YAML not found: {self.yaml_path}")
    
    def tearDown(self):
        """Clean up."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_structure_counts_are_derived_from_model_data(self):
        """Counts must be derived from species and reactions, not stored in YAML."""
        orchestrator = ACGOrchestrator(
            self.yaml_path,
            str(self.output_dir),
            verbose=False
        )

        orchestrator.load_config()
        orchestrator.evaluate_formulas()
        acg_data = orchestrator.map_to_acg_structures()

        self.assertNotIn('structure', orchestrator.config)
        self.assertEqual(acg_data['nsolids'], 6)
        self.assertEqual(acg_data['ndissolved'], 12)
        self.assertEqual(acg_data['nreactions'], 19)
        self.assertEqual(acg_data['neqrxns'], 3)

    def test_complete_pipeline_phases_1_to_4(self):
        """Test complete pipeline through phase 4 (pre-processing)."""
        orchestrator = ACGOrchestrator(
            self.yaml_path,
            str(self.output_dir),
            verbose=False
        )
        
        # Run phases 1-4
        orchestrator.load_config()
        orchestrator.evaluate_formulas()
        orchestrator.map_to_acg_structures()
        orchestrator.run_preprocessing()
        
        # Verify results
        self.assertEqual(orchestrator.acg_data['ncompo'], 18)
        self.assertEqual(orchestrator.acg_data['ndissolved'], 12)
        self.assertEqual(orchestrator.acg_data['nsolids'], 6)
        self.assertEqual(orchestrator.acg_data['nreactions'], 19)
        self.assertEqual(orchestrator.acg_data['neqrxns'], 3)
        
        # Verify reduced system
        self.assertIsNotNone(orchestrator.reduced_system)
        self.assertEqual(len(orchestrator.reduced_system['func']), 18)
        self.assertEqual(orchestrator.reduced_system['jacobian'].shape, (18, 18))


class TestErrorHandling(unittest.TestCase):
    """Test error handling in orchestrator."""
    
    def setUp(self):
        """Create temporary directory."""
        self.temp_dir = tempfile.mkdtemp()
        self.output_dir = Path(self.temp_dir) / 'output'
    
    def tearDown(self):
        """Clean up."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_generate_without_yaml(self):
        """Test that generate() fails gracefully with missing YAML."""
        orchestrator = ACGOrchestrator(
            'nonexistent.yaml',
            str(self.output_dir)
        )
        
        with self.assertRaises(ACGOrchestrationError):
            orchestrator.generate()

    def test_generate_preserves_orchestration_errors(self):
        """generate() should not wrap existing ACGOrchestrationErrors again."""
        yaml_path = Path(self.temp_dir) / 'invalid.yaml'
        config = {
            'species': [{'name': 'o2'}],
            'reactions': [],
            'parameters': {'biogeochemical': []}
        }
        with open(yaml_path, 'w') as f:
            yaml.dump(config, f)

        orchestrator = ACGOrchestrator(str(yaml_path), str(self.output_dir))

        with self.assertRaises(ACGOrchestrationError) as ctx:
            orchestrator.generate()
        msg = str(ctx.exception)
        self.assertIn("'type' is missing.", msg)
        self.assertNotIn('Pipeline failed:', msg)
    
    def test_map_before_load(self):
        """Test that mapping fails if config not loaded."""
        yaml_path = Path(self.temp_dir) / 'test.yaml'
        config = {
            'species': [{'name': 'o2', 'type': 'dissolved'}],
            'reactions': [{'id': 1, 'name': 'test', 'stoichiometry': {'o2': -1}}],
            'parameters': {'biogeochemical': []}
        }
        with open(yaml_path, 'w') as f:
            yaml.dump(config, f)
        
        orchestrator = ACGOrchestrator(str(yaml_path), str(self.output_dir))
        
        # Try to map without loading
        with self.assertRaises(Exception):
            orchestrator.map_to_acg_structures()


class TestEquilibriumConstraintPreprocessing(unittest.TestCase):
    """Tests for equilibrium constraint handling in preprocessing."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.output_dir = Path(self.temp_dir) / 'output'

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def _write_config(self, config, name='test.yaml'):
        yaml_path = Path(self.temp_dir) / name
        with open(yaml_path, 'w') as handle:
            yaml.dump(config, handle)
        return yaml_path

    def _run_with_capture(self, config):
        yaml_path = self._write_config(config)
        orchestrator = ACGOrchestrator(str(yaml_path), str(self.output_dir), verbose=False)
        orchestrator.load_config()
        orchestrator.evaluate_formulas()
        orchestrator.map_to_acg_structures()

        captured = {}

        def fake_run_gaussian_elimination(**kwargs):
            captured['equilibrium_eqns'] = kwargs.get('equilibrium_eqns')
            ncompo = len(kwargs['variables'])
            return {
                'func': [Symbol('f')] * ncompo,
                'jacobian': eye(ncompo),
                'matrixM': eye(ncompo),
                'rightM': zeros(ncompo, 1),
            }

        with patch('acg_brns.acg_orchestrator.run_gaussian_elimination', side_effect=fake_run_gaussian_elimination):
            orchestrator.run_preprocessing()

        return captured['equilibrium_eqns']

    def test_preprocessing_uses_raw_yaml_equilibrium_constraint(self):
        """Preprocessing must differentiate species-name constraints, not sp(i,j) strings."""
        config = {
            'species': [
                {'name': 'co2', 'type': 'dissolved'},
                {'name': 'hco3', 'type': 'dissolved'},
            ],
            'reactions': [
                {
                    'id': 2,
                    'name': 'carbonate_eq',
                    'equilibrium': True,
                    'rate': 0,
                    'equilibrium_constraint': 'hco3 - keq1 * co2',
                    'stoichiometry': {'co2': -1, 'hco3': 1},
                }
            ],
            'parameters': {
                'physical': {},
                'biogeochemical': [
                    {'name': 'keq1', 'value': 1.0},
                ],
            },
        }

        equilibrium_eqns = self._run_with_capture(config)

        self.assertIsNotNone(equilibrium_eqns)
        self.assertEqual(len(equilibrium_eqns), 1)

        expr = equilibrium_eqns[0]
        self.assertIn(Symbol('hco3'), expr.free_symbols)
        self.assertIn(Symbol('co2'), expr.free_symbols)
        self.assertNotIn('sp(', str(expr))

    def test_preprocessing_accepts_integer_d_notation_in_equilibrium_constraint(self):
        """Fortran integer d-notation like 1d0 should parse in preprocessing."""
        config = {
            'species': [
                {'name': 'co2', 'type': 'dissolved'},
                {'name': 'hco3', 'type': 'dissolved'},
            ],
            'reactions': [
                {
                    'id': 2,
                    'name': 'carbonate_eq',
                    'equilibrium': True,
                    'rate': 0,
                    'equilibrium_constraint': 'hco3 - 1d0 * co2',
                    'stoichiometry': {'co2': -1, 'hco3': 1},
                }
            ],
            'parameters': {
                'physical': {},
                'biogeochemical': [],
            },
        }

        equilibrium_eqns = self._run_with_capture(config)

        self.assertIsNotNone(equilibrium_eqns)
        self.assertEqual(len(equilibrium_eqns), 1)

        from sympy import nsimplify
        expr = equilibrium_eqns[0]
        expected = Symbol('hco3') - Symbol('co2')
        # 1d0 is parsed as float 1.0; nsimplify converts to rational so
        # -1.0*co2 + hco3 compares equal to hco3 - co2.
        self.assertEqual(nsimplify(expr), expected)


class TestFailFastParsing(unittest.TestCase):
    """Regression tests for fail-fast parsing in later pipeline phases."""

    def test_run_preprocessing_fails_on_bad_mapped_rate_expression(self):
        orchestrator = ACGOrchestrator.__new__(ACGOrchestrator)
        orchestrator.verbose = False
        orchestrator.config = {
            'reactions': [{'id': 1, 'stoichiometry': {'o2': -1}}],
            'parameters': {},
        }
        orchestrator.evaluated_params = {}
        orchestrator.mapper = None
        orchestrator.acg_data = {
            'variables': ['o2'],
            'bio_name': [],
            'reactions': [{'id': 1, 'rate_expanded': 'bad((', 'equilibrium': False}],
            'reaction_ids': [1],
            'ncompo': 1,
            'nreactions': 1,
        }

        with self.assertRaises(ACGOrchestrationError) as ctx:
            orchestrator.run_preprocessing()
        self.assertIn('Failed to parse mapped rate expression', str(ctx.exception))


class TestSwitchConditionNormalization(unittest.TestCase):
    """Tests for advanced.switches condition normalization."""

    def setUp(self):
        self.orchestrator = ACGOrchestrator.__new__(ACGOrchestrator)

    def test_maple_style_condition_is_preserved(self):
        cond = self.orchestrator._normalize_switch_condition('30.0 - x(j)')
        self.assertEqual(cond, '30.0 - x(j)')

    def test_lt_operator_is_normalized(self):
        cond = self.orchestrator._normalize_switch_condition('x(j) < 30.0')
        self.assertEqual(cond, 'x(j).lt.30.0')

    def test_all_comparison_operators_are_normalized(self):
        self.assertEqual(
            self.orchestrator._normalize_switch_condition('x(j) <= 30.0'),
            'x(j).le.30.0',
        )
        self.assertEqual(
            self.orchestrator._normalize_switch_condition('x(j) > 30.0'),
            'x(j).gt.30.0',
        )
        self.assertEqual(
            self.orchestrator._normalize_switch_condition('x(j) >= 30.0'),
            'x(j).ge.30.0',
        )
        self.assertEqual(
            self.orchestrator._normalize_switch_condition('x(j) == 30.0'),
            'x(j).eq.30.0',
        )
        self.assertEqual(
            self.orchestrator._normalize_switch_condition('x(j) != 30.0'),
            'x(j).ne.30.0',
        )

    def test_logical_words_are_normalized(self):
        cond = self.orchestrator._normalize_switch_condition(
            'x(j) < 30.0 and x(j) > 1.0e-2'
        )
        self.assertEqual(cond, 'x(j).lt.30.0 .and. x(j).gt.1.0e-2')

    def test_empty_condition_raises_error(self):
        with self.assertRaises(ACGOrchestrationError):
            self.orchestrator._normalize_switch_condition('   ')

    def test_build_net_rates_fails_on_bad_stoichiometric_coefficient(self):
        orchestrator = ACGOrchestrator.__new__(ACGOrchestrator)
        orchestrator.verbose = False
        orchestrator.config = {
            'reactions': [{'id': 1, 'stoichiometry': {'o2': 'bad(('}}],
            'parameters': {},
        }
        orchestrator.evaluated_params = {}
        orchestrator.mapper = None
        orchestrator.acg_data = {
            'variables': ['o2'],
            'bio_name': [],
            'reactions': [{'id': 1, 'rate_expanded': '1.0'}],
            'ncompo': 1,
        }

        with self.assertRaises(ACGOrchestrationError) as ctx:
            orchestrator._build_net_reaction_rates([Symbol('o2')])
        self.assertIn('Failed to parse net-rate stoichiometric coefficient', str(ctx.exception))


class TestValidationMessages(unittest.TestCase):
    """Tests for the multi-phase YAML validation with user-friendly messages."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.output_dir = Path(self.temp_dir) / 'output'

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def _make_yaml(self, config, name='test.yaml'):
        """Write a config dict to a temp YAML file and return its path."""
        p = Path(self.temp_dir) / name
        with open(p, 'w') as fh:
            yaml.dump(config, fh)
        return p

    def _minimal_config(self):
        """Return a minimal, fully-valid config for incremental tests."""
        return {
            'species': [
                {'name': 'o2', 'type': 'dissolved'},
                {'name': 'ch2o', 'type': 'solid'},
            ],
            'reactions': [
                {'id': 1, 'name': 'aerobic', 'stoichiometry': {'o2': -1, 'ch2o': -1}},
            ],
            'parameters': {
                'biogeochemical': [{'name': 'kox', 'value': 0.1}],
                'physical': {'por0': 0.8},
                'stoichiometry': {},
            },
        }

    def _load(self, config):
        """Write config and call load_config(), return the orchestrator."""
        p = self._make_yaml(config)
        orch = ACGOrchestrator(str(p), str(self.output_dir), verbose=False)
        orch.load_config()
        return orch

    def _assert_error_contains(self, config, text):
        """Assert that load_config raises ACGOrchestrationError containing `text`."""
        with self.assertRaises(ACGOrchestrationError) as ctx:
            self._load(config)
        self.assertIn(text, str(ctx.exception))

    # ------------------------------------------------------------------
    # I — Infrastructure
    # ------------------------------------------------------------------

    def test_validation_issue_str_with_hint(self):
        """ValidationIssue.__str__ includes severity, path, message and hint."""
        issue = ValidationIssue('ERROR', 'species[0].type', "Unknown type 'blob'.", 'Allowed values: dissolved, solid.')
        s = str(issue)
        self.assertIn('ERROR', s)
        self.assertIn('species[0].type', s)
        self.assertIn("blob", s)
        self.assertIn('dissolved', s)

    def test_validation_issue_str_without_hint(self):
        """ValidationIssue.__str__ without hint has no arrow line."""
        issue = ValidationIssue('WARNING', 'parameters.physical.por0', 'Value > 1.')
        s = str(issue)
        self.assertNotIn('→', s)

    def test_valid_config_raises_no_error(self):
        """Minimal valid config must not raise."""
        self._load(self._minimal_config())  # no exception expected

    def test_error_report_counts_errors_and_warnings(self):
        """Report header must state number of errors."""
        cfg = self._minimal_config()
        cfg['species'][0]['type'] = 'gaseous'  # wrong type → ERROR
        with self.assertRaises(ACGOrchestrationError) as ctx:
            self._load(cfg)
        msg = str(ctx.exception)
        self.assertIn('errors', msg)

    # ------------------------------------------------------------------
    # S — Schema phase
    # ------------------------------------------------------------------

    def test_s5_invalid_species_type(self):
        """S-5: unknown species type produces ERROR."""
        cfg = self._minimal_config()
        cfg['species'][0]['type'] = 'gaseous'
        self._assert_error_contains(cfg, "gaseous")

    def test_s7_equilibrium_missing_constraint(self):
        """S-7: equilibrium reaction without equilibrium_constraint → ERROR."""
        cfg = self._minimal_config()
        cfg['reactions'].append({
            'id': 2, 'name': 'carb_eq', 'equilibrium': True,
            'stoichiometry': {'o2': 1}
        })
        self._assert_error_contains(cfg, 'equilibrium_constraint')

    def test_s7_equilibrium_with_constraint_ok(self):
        """S-7: equilibrium reaction WITH constraint must not raise."""
        cfg = self._minimal_config()
        cfg['reactions'].append({
            'id': 2, 'name': 'carb_eq', 'equilibrium': True,
            'equilibrium_constraint': 'o2 - kox',
            'stoichiometry': {'o2': 1},
            'rate': 0,
        })
        self._load(cfg)  # no exception

    def test_schema_missing_species_name(self):
        """Species without 'name' → ERROR."""
        cfg = self._minimal_config()
        del cfg['species'][0]['name']
        self._assert_error_contains(cfg, "'name' is missing")

    def test_schema_species_name_must_be_string(self):
        """Species name must be a string."""
        cfg = self._minimal_config()
        cfg['species'][0]['name'] = 123
        self._assert_error_contains(cfg, "'name' must be a string")

    def test_schema_missing_species_type(self):
        """Species without 'type' → ERROR."""
        cfg = self._minimal_config()
        del cfg['species'][0]['type']
        self._assert_error_contains(cfg, "'type' is missing")

    def test_schema_reaction_missing_id(self):
        """Reaction without 'id' → ERROR."""
        cfg = self._minimal_config()
        del cfg['reactions'][0]['id']
        self._assert_error_contains(cfg, "'id' is missing")

    def test_schema_reaction_id_must_be_integer(self):
        """Reaction id must be an integer."""
        cfg = self._minimal_config()
        cfg['reactions'][0]['id'] = 'one'
        self._assert_error_contains(cfg, "'id' must be an integer")

    def test_schema_reaction_missing_name_is_allowed(self):
        """Reaction name is optional and falls back to a generated label."""
        cfg = self._minimal_config()
        del cfg['reactions'][0]['name']
        orch = self._load(cfg)
        acg_data = orch.map_to_acg_structures()
        self.assertEqual(acg_data['reactions'][0]['name'], 'reaction_1')

    def test_schema_reaction_name_must_be_string(self):
        """Reaction name must be a string."""
        cfg = self._minimal_config()
        cfg['reactions'][0]['name'] = 5
        self._assert_error_contains(cfg, "'name' must be a string")

    def test_schema_reaction_missing_stoichiometry(self):
        """Reaction without 'stoichiometry' → ERROR."""
        cfg = self._minimal_config()
        del cfg['reactions'][0]['stoichiometry']
        self._assert_error_contains(cfg, "'stoichiometry' is missing")

    def test_schema_reaction_stoichiometry_must_be_mapping(self):
        """Reaction stoichiometry must be a mapping."""
        cfg = self._minimal_config()
        cfg['reactions'][0]['stoichiometry'] = []
        self._assert_error_contains(cfg, "'stoichiometry' must be a mapping")

    # ------------------------------------------------------------------
    # R — Reference phase
    # ------------------------------------------------------------------

    def test_r1_duplicate_species_name(self):
        """R-1: duplicate species name → ERROR."""
        cfg = self._minimal_config()
        cfg['species'].append({'name': 'o2', 'type': 'dissolved'})  # duplicate
        self._assert_error_contains(cfg, "not unique")

    def test_r2_duplicate_reaction_id(self):
        """R-2: duplicate reaction ID → ERROR."""
        cfg = self._minimal_config()
        cfg['reactions'].append({'id': 1, 'name': 'dup', 'stoichiometry': {'o2': 1}})
        self._assert_error_contains(cfg, "not unique")

    def test_reaction_ids_may_be_non_contiguous(self):
        """Unique IDs order reactions but do not define Fortran rate indices."""
        cfg = self._minimal_config()
        cfg['reactions'][0]['id'] = 20
        cfg['reactions'].append({
            'id': 10, 'name': 'earlier', 'output': True,
            'stoichiometry': {'o2': 1},
        })

        orchestrator = self._load(cfg)
        orchestrator.evaluate_formulas()
        orchestrator.map_to_acg_structures()
        params = orchestrator._extract_generation_parameters()

        self.assertEqual(orchestrator.acg_data['reaction_ids'], [10, 20])
        self.assertEqual(params['listroutput'], [1])

    def test_r3_unknown_species_in_stoichiometry(self):
        """R-3: unknown species name in stoichiometry → ERROR with species name."""
        cfg = self._minimal_config()
        cfg['reactions'][0]['stoichiometry']['feii'] = -1  # typo
        self._assert_error_contains(cfg, 'feii')

    def test_r3_known_species_in_stoichiometry_ok(self):
        """R-3: known species in stoichiometry must not raise."""
        cfg = self._minimal_config()
        # ch2o is already a known species — no error expected
        self._load(cfg)

    def test_r4_unknown_symbol_in_equilibrium_constraint(self):
        """R-4: unknown symbol in equilibrium_constraint → ERROR."""
        cfg = self._minimal_config()
        cfg['reactions'].append({
            'id': 3, 'name': 'eq', 'equilibrium': True,
            'equilibrium_constraint': 'o2 - phantom_param',
            'stoichiometry': {'o2': 1}, 'rate': 0,
        })
        self._assert_error_contains(cfg, 'phantom_param')

    def test_r4_known_symbols_in_equilibrium_constraint_ok(self):
        """R-4: all symbols known in equilibrium_constraint must not raise."""
        cfg = self._minimal_config()
        cfg['reactions'].append({
            'id': 3, 'name': 'eq', 'equilibrium': True,
            'equilibrium_constraint': 'o2 - kox',  # kox is in biogeochemical
            'stoichiometry': {'o2': 1}, 'rate': 0,
        })
        self._load(cfg)

    def test_r4_d_notation_in_equilibrium_constraint_ok(self):
        """R-4: Fortran d-notation like '1d0' must not trigger unknown-symbol error."""
        cfg = self._minimal_config()
        cfg['reactions'].append({
            'id': 3, 'name': 'eq', 'equilibrium': True,
            'equilibrium_constraint': 'o2 - 1d0 * kox',
            'stoichiometry': {'o2': 1}, 'rate': 0,
        })
        self._load(cfg)  # no exception

    def test_r5_unknown_symbol_in_rate(self):
        """R-5: unknown symbol in rate string → ERROR."""
        cfg = self._minimal_config()
        cfg['reactions'][0]['rate'] = 'ghost_rate * o2'
        self._assert_error_contains(cfg, 'ghost_rate')

    def test_r5_known_symbol_in_rate_ok(self):
        """R-5: known parameter in rate string must not raise."""
        cfg = self._minimal_config()
        cfg['reactions'][0]['rate'] = 'kox * o2'
        self._load(cfg)

    def test_r5_whitelisted_math_function_in_rate_ok(self):
        """R-5: whitelisted math function names in rate must not raise symbol errors."""
        cfg = self._minimal_config()
        cfg['reactions'][0]['rate'] = 'exp(-kox) * o2'
        self._load(cfg)

    def test_r5_rate_components_unknown_symbol(self):
        """R-5: unknown symbol in rate_components → ERROR."""
        cfg = self._minimal_config()
        cfg['reactions'][0]['rate_components'] = {'helper': 'unknown_sym * o2'}
        self._assert_error_contains(cfg, 'unknown_sym')

    def test_r5_rate_components_local_helper_ok(self):
        """R-5: rate_component names may reference each other without error."""
        cfg = self._minimal_config()
        cfg['reactions'][0]['rate_components'] = {
            'fo2': 'o2 / kox',
            'rate_expr': 'fo2 * kox',
        }
        self._load(cfg)  # fo2 is local helper, no exception

    def test_r5_rate_uses_rate_component_names_ok(self):
        """R-5: top-level 'rate' may reference names defined in rate_components."""
        cfg = self._minimal_config()
        # equilibrium model-style: rate references helpers defined in rate_components
        cfg['reactions'][0]['rate'] = 'kch2o_helper * fo2_helper'
        cfg['reactions'][0]['rate_components'] = {
            'kch2o_helper': 'kox * ch2o',
            'fo2_helper': 'o2 / kox',
        }
        self._load(cfg)  # both helpers are local, no exception

    def test_extract_tokens_d_notation(self):
        """_extract_tokens must not treat Fortran d-notation exponent as a symbol."""
        tokens = ACGOrchestrator._extract_tokens('hco3 - 1d0 * co2')
        self.assertIn('hco3', tokens)
        self.assertIn('co2', tokens)
        self.assertNotIn('E0', tokens)
        self.assertNotIn('d0', tokens)

    def test_extract_tokens_sci_notation(self):
        """_extract_tokens must not treat Python-style scientific notation as symbols."""
        tokens = ACGOrchestrator._extract_tokens('rate * 1.5e-3 + keq1')
        self.assertIn('rate', tokens)
        self.assertIn('keq1', tokens)
        self.assertNotIn('e', tokens)
        self.assertNotIn('E', tokens)

    def test_r6_duplicate_bio_param_name(self):
        """R-6: duplicate biogeochemical parameter name → ERROR."""
        cfg = self._minimal_config()
        cfg['parameters']['biogeochemical'].append({'name': 'kox', 'value': 0.2})
        self._assert_error_contains(cfg, "not unique")

    def test_r7_parameter_name_conflicts_with_species_or_fortran_keywords(self):
        """R-7: parameter names must not override model symbols or reserved Fortran identifiers."""
        cfg = self._minimal_config()
        cfg['parameters']['biogeochemical'].append({'name': 'o2', 'value': 0.2})
        self._assert_error_contains(cfg, "not unique")

        cfg = self._minimal_config()
        cfg['parameters']['biogeochemical'] = [{'name': 'integer', 'value': 2.0}]
        self._assert_error_contains(cfg, "reserved")

    def test_multiple_errors_reported_together(self):
        """Multiple errors are collected and shown in a single exception."""
        cfg = self._minimal_config()
        cfg['species'][0]['type'] = 'gaseous'                        # S-5
        cfg['reactions'][0]['stoichiometry']['phantom'] = -1          # R-3
        with self.assertRaises(ACGOrchestrationError) as ctx:
            self._load(cfg)
        msg = str(ctx.exception)
        self.assertIn('gaseous', msg)          # first error present
        self.assertIn('phantom', msg)          # second error present


class TestSchemaS8(unittest.TestCase):
    """S-8: scalar/formula fields must not be list or dict."""

    def _make_yaml(self, cfg: dict) -> str:
        import yaml
        return yaml.dump(cfg, allow_unicode=True)

    def _minimal_config(self):
        return {
            'species': [{
                'name': 'o2', 'type': 'dissolved',
                'bc_upper_type': 0, 'bc_upper_value': 0.21,
                'bc_lower_type': 1, 'bc_lower_value': 0.0,
                'init_value': 0.0, 'transport_D0': 1e-9, 'transport_alpha': 0.0,
            }],
            'reactions': [],
            'parameters': {
                'biogeochemical': [],
                'stoichiometric': [],
                'physical': {'por0': 0.8},
            },
        }

    def _load(self, cfg):
        import yaml
        from io import StringIO
        orchestrator = ACGOrchestrator.__new__(ACGOrchestrator)
        orchestrator.verbose = False
        orchestrator.config = cfg
        orchestrator._validation_warnings = []
        from acg_brns.acg_orchestrator import ValidationIssue
        return orchestrator._validate_schema()

    def test_s8_list_in_bc_upper_value_raises_error(self):
        cfg = self._minimal_config()
        cfg['species'][0]['bc_upper_value'] = [0.1, 0.2]
        issues = self._load(cfg)
        errors = [i for i in issues if i.severity == 'ERROR' and 'bc_upper_value' in i.location]
        self.assertTrue(len(errors) >= 1, f"Expected S-8 error, got: {issues}")
        self.assertIn('list', errors[0].message)

    def test_s8_dict_in_init_value_raises_error(self):
        cfg = self._minimal_config()
        cfg['species'][0]['init_value'] = {'a': 1}
        issues = self._load(cfg)
        errors = [i for i in issues if i.severity == 'ERROR' and 'init_value' in i.location]
        self.assertTrue(len(errors) >= 1, f"Expected S-8 error, got: {issues}")
        self.assertIn('dict', errors[0].message)

    def test_s8_number_ok(self):
        cfg = self._minimal_config()
        cfg['species'][0]['bc_upper_value'] = 0.21
        issues = self._load(cfg)
        errors = [i for i in issues if i.severity == 'ERROR' and 'bc_upper_value' in i.location]
        self.assertEqual(errors, [])

    def test_s8_formula_string_ok(self):
        cfg = self._minimal_config()
        cfg['species'][0]['bc_upper_value'] = "C_hco3"
        issues = self._load(cfg)
        errors = [i for i in issues if i.severity == 'ERROR' and 'bc_upper_value' in i.location]
        self.assertEqual(errors, [])


class TestReferencesR7(unittest.TestCase):
    """R-7: bc/init formula string references must resolve to known symbols."""

    def _make_yaml(self, cfg: dict) -> str:
        import yaml
        return yaml.dump(cfg, allow_unicode=True)

    def _minimal_config(self):
        return {
            'species': [{
                'name': 'o2', 'type': 'dissolved',
                'bc_upper_type': 0, 'bc_upper_value': 0.21,
                'bc_lower_type': 1, 'bc_lower_value': 0.0,
                'init_value': 0.0, 'transport_D0': 1e-9, 'transport_alpha': 0.0,
            }],
            'reactions': [],
            'parameters': {
                'biogeochemical': [{'name': 'kox', 'value': 0.1}],
                'stoichiometric': [],
                'physical': {'por0': 0.8},
            },
        }

    def _load(self, cfg):
        orchestrator = ACGOrchestrator.__new__(ACGOrchestrator)
        orchestrator.verbose = False
        orchestrator.config = cfg
        orchestrator._validation_warnings = []
        return orchestrator._validate_references()

    def test_r7_unknown_symbol_in_bc_upper_value(self):
        cfg = self._minimal_config()
        cfg['species'][0]['bc_upper_value'] = "totally_unknown_var"
        issues = self._load(cfg)
        errors = [i for i in issues if i.severity == 'ERROR' and 'bc_upper_value' in i.location]
        self.assertTrue(len(errors) >= 1, f"Expected R-7 error, got: {issues}")
        self.assertIn('totally_unknown_var', errors[0].message)

    def test_r7_known_species_in_bc_value_ok(self):
        cfg = self._minimal_config()
        cfg['species'][0]['bc_upper_value'] = "kox"  # defined in biogeochemical params
        issues = self._load(cfg)
        errors = [i for i in issues if i.severity == 'ERROR' and 'bc_upper_value' in i.location]
        self.assertEqual(errors, [], f"Unexpected R-7 error: {errors}")

    def test_r7_computed_values_symbol_ok(self):
        cfg = self._minimal_config()
        cfg['computed_values'] = {'carbonate_system': {'C_hco3': 'some_formula'}}
        cfg['species'][0]['bc_upper_value'] = "C_hco3"
        issues = self._load(cfg)
        errors = [i for i in issues if i.severity == 'ERROR' and 'bc_upper_value' in i.location]
        self.assertEqual(errors, [], f"Unexpected R-7 error: {errors}")

    def test_r7_unknown_symbol_in_init_value(self):
        cfg = self._minimal_config()
        cfg['species'][0]['init_value'] = "ghost_var"
        issues = self._load(cfg)
        errors = [i for i in issues if i.severity == 'ERROR' and 'init_value' in i.location]
        self.assertTrue(len(errors) >= 1, f"Expected R-7 error for init_value, got: {issues}")
        self.assertIn('ghost_var', errors[0].message)

    def test_r5_unknown_function_in_rate_emits_warning(self):
        """Unknown function calls in reference phase should warn instead of error."""
        cfg = self._minimal_config()
        cfg['reactions'] = [{
            'id': 1,
            'name': 'rx1',
            'stoichiometry': {'o2': -1},
            'rate': 'mystery_fn(o2)'
        }]
        issues = self._load(cfg)
        warnings = [
            i for i in issues
            if i.severity == 'WARNING' and i.path == 'reactions[0].rate'
        ]
        errors = [i for i in issues if i.severity == 'ERROR' and i.path == 'reactions[0].rate']
        self.assertTrue(warnings, f"Expected warning for unknown function, got: {issues}")
        self.assertIn("Unknown function 'mystery_fn'", warnings[0].message)
        self.assertEqual(errors, [], f"Unknown function should not be an ERROR in reference phase: {issues}")


class TestFormulaValidation(unittest.TestCase):
    """F-1, F-3, F-4: syntactic formula checks via sympify."""

    def _minimal_config(self):
        return {
            'species': [{
                'name': 'o2', 'type': 'dissolved',
                'bc_upper_type': 0, 'bc_upper_value': 0.21,
                'bc_lower_type': 1, 'bc_lower_value': 0.0,
                'init_value': 0.0, 'transport_D0': 1e-9, 'transport_alpha': 0.0,
            }],
            'reactions': [],
            'parameters': {
                'biogeochemical': [],
                'stoichiometric': [],
                'physical': {'por0': 0.8},
            },
        }

    def _run(self, cfg):
        orchestrator = ACGOrchestrator.__new__(ACGOrchestrator)
        orchestrator.verbose = False
        orchestrator.config = cfg
        orchestrator._validation_warnings = []
        return orchestrator._validate_formulas()

    def test_f1_unparsable_bio_formula(self):
        cfg = self._minimal_config()
        cfg['parameters']['biogeochemical'] = [{'name': 'kbad', 'value': '@@invalid@@'}]
        issues = self._run(cfg)
        errors = [i for i in issues if i.severity == 'ERROR']
        self.assertTrue(len(errors) >= 1, f"Expected F-1 error, got: {issues}")
        self.assertIn('kbad', errors[0].message)

    def test_f1_parsable_bio_formula_ok(self):
        cfg = self._minimal_config()
        cfg['parameters']['biogeochemical'] = [{'name': 'k1', 'value': '10^(-3)'}]
        issues = self._run(cfg)
        self.assertEqual(issues, [], f"Unexpected errors: {issues}")

    def test_f1_d_notation_ok(self):
        cfg = self._minimal_config()
        cfg['parameters']['biogeochemical'] = [{'name': 'k1', 'value': '1.5d-3'}]
        issues = self._run(cfg)
        self.assertEqual(issues, [], f"D-notation should parse fine: {issues}")

    def test_f2_direct_cycle_detected(self):
        cfg = self._minimal_config()
        cfg['parameters']['biogeochemical'] = [
            {'name': 'k1', 'value': 'k2 + 1'},
            {'name': 'k2', 'value': 'k1 + 2'},
        ]
        issues = self._run(cfg)
        cycle_errors = [
            i for i in issues
            if i.severity == 'ERROR' and (
                'Circular dependency' in i.message
            )
        ]
        self.assertTrue(len(cycle_errors) >= 1, f"Expected F-2 cycle error, got: {issues}")
        self.assertEqual(len(cycle_errors), 1, f"Expected one consolidated cycle error, got: {issues}")

    def test_f2_indirect_cycle_detected(self):
        cfg = self._minimal_config()
        cfg['parameters']['biogeochemical'] = [
            {'name': 'k1', 'value': 'k2 + 1'},
            {'name': 'k2', 'value': 'k3 + 2'},
            {'name': 'k3', 'value': 'k1 + 3'},
        ]
        issues = self._run(cfg)
        cycle_errors = [i for i in issues if i.severity == 'ERROR' and 'Circular dependency' in i.message]
        self.assertTrue(len(cycle_errors) >= 1, f"Expected indirect F-2 cycle error, got: {issues}")

    def test_f2_self_cycle_detected(self):
        cfg = self._minimal_config()
        cfg['parameters']['biogeochemical'] = [
            {'name': 'k1', 'value': 'k1 + 1'},
        ]
        issues = self._run(cfg)
        cycle_errors = [i for i in issues if i.severity == 'ERROR' and 'Circular dependency' in i.message]
        self.assertTrue(len(cycle_errors) >= 1, f"Expected self-cycle F-2 error, got: {issues}")

    def test_f2_acyclic_chain_ok(self):
        cfg = self._minimal_config()
        cfg['parameters']['biogeochemical'] = [
            {'name': 'k1', 'value': '2'},
            {'name': 'k2', 'value': 'k1 + 1'},
            {'name': 'k3', 'value': 'k2 + 1'},
        ]
        issues = self._run(cfg)
        cycle_errors = [i for i in issues if i.severity == 'ERROR' and 'Circular dependency' in i.message]
        self.assertEqual(cycle_errors, [], f"No F-2 errors expected for acyclic formulas: {issues}")

    def test_f3_unparsable_equilibrium_constraint(self):
        cfg = self._minimal_config()
        cfg['reactions'] = [{'id': 1, 'type': 'aerobic',
                              'stoichiometry': {}, 'rate_components': {},
                              'equilibrium_constraint': '@@bad@@'}]
        issues = self._run(cfg)
        errors = [i for i in issues if 'equilibrium_constraint' in i.location]
        self.assertTrue(len(errors) >= 1, f"Expected F-3 error, got: {issues}")

    def test_f4_unparsable_rate(self):
        cfg = self._minimal_config()
        cfg['reactions'] = [{'id': 1, 'type': 'aerobic',
                              'stoichiometry': {}, 'rate_components': {},
                              'rate': '@@bad@@'}]
        issues = self._run(cfg)
        errors = [i for i in issues if '.rate' in i.location and 'equilibrium' not in i.location]
        self.assertTrue(len(errors) >= 1, f"Expected F-4 rate error, got: {issues}")

    def test_f4_unparsable_rate_component(self):
        cfg = self._minimal_config()
        cfg['reactions'] = [{'id': 1, 'type': 'aerobic',
                              'stoichiometry': {}, 'rate_components': {'fO2': '@@bad@@'}}]
        issues = self._run(cfg)
        errors = [i for i in issues if 'rate_components' in i.location]
        self.assertTrue(len(errors) >= 1, f"Expected F-4 rate_component error, got: {issues}")
        self.assertIn('fO2', errors[0].message)

    def test_f_valid_rate_with_log10_caret(self):
        cfg = self._minimal_config()
        cfg['reactions'] = [{'id': 1, 'type': 'aerobic',
                              'stoichiometry': {}, 'rate_components': {},
                              'rate': 'kox * log10(o2) + sqrt(2)'}]
        issues = self._run(cfg)
        self.assertEqual(issues, [], f"Should parse without errors: {issues}")

    def test_f_valid_nested_log10_parentheses(self):
        cfg = self._minimal_config()
        cfg['reactions'] = [{'id': 1, 'type': 'aerobic',
                              'stoichiometry': {}, 'rate_components': {},
                              'rate': 'log10(kox * (o2 + 1))'}]
        issues = self._run(cfg)
        self.assertEqual(issues, [], f"Nested log10 expression should parse cleanly: {issues}")


class TestPlausibilityValidation(unittest.TestCase):
    """N-1, N-2, N-3: plausibility warnings."""

    def _minimal_config(self):
        return {
            'species': [{
                'name': 'o2', 'type': 'dissolved',
                'bc_upper_type': 0, 'bc_upper_value': 0.21,
                'bc_lower_type': 1, 'bc_lower_value': 0.0,
                'init_value': 0.0, 'transport_D0': 1e-9, 'transport_alpha': 0.0,
            }],
            'reactions': [],
            'parameters': {
                'biogeochemical': [],
                'stoichiometric': [],
                'physical': {'por0': 0.8},
            },
        }

    def _run(self, cfg):
        orchestrator = ACGOrchestrator.__new__(ACGOrchestrator)
        orchestrator.verbose = False
        orchestrator.config = cfg
        orchestrator._validation_warnings = []
        return orchestrator._validate_plausibility()

    def test_n1_por0_too_high_warning(self):
        cfg = self._minimal_config()
        cfg['parameters']['physical']['por0'] = 1.2
        issues = self._run(cfg)
        warnings = [i for i in issues if i.severity == 'WARNING' and 'por0' in i.location]
        self.assertTrue(len(warnings) >= 1, f"Expected N-1 warning for por0=1.2: {issues}")
        self.assertIn('1.2', warnings[0].message)

    def test_n1_por0_zero_warning(self):
        cfg = self._minimal_config()
        cfg['parameters']['physical']['por0'] = 0.0
        issues = self._run(cfg)
        warnings = [i for i in issues if i.severity == 'WARNING' and 'por0' in i.location]
        self.assertTrue(len(warnings) >= 1, f"Expected N-1 warning for por0=0: {issues}")

    def test_n1_por0_valid_no_warning(self):
        cfg = self._minimal_config()
        cfg['parameters']['physical']['por0'] = 0.8
        issues = self._run(cfg)
        warnings = [i for i in issues if 'por0' in i.location]
        self.assertEqual(warnings, [], f"No warning expected for valid por0: {warnings}")

    def test_n2_negative_D0_warning(self):
        cfg = self._minimal_config()
        cfg['species'][0]['transport_D0'] = -1e-9
        issues = self._run(cfg)
        warnings = [i for i in issues if i.severity == 'WARNING' and 'transport_D0' in i.location]
        self.assertTrue(len(warnings) >= 1, f"Expected N-2 warning for negative D0: {issues}")
        self.assertIn('-', warnings[0].message)

    def test_n2_negative_alpha_warning(self):
        cfg = self._minimal_config()
        cfg['species'][0]['transport_alpha'] = -0.5
        issues = self._run(cfg)
        warnings = [i for i in issues if i.severity == 'WARNING' and 'transport_alpha' in i.location]
        self.assertTrue(len(warnings) >= 1, f"Expected N-2 warning for negative alpha: {issues}")

    def test_n2_positive_transport_ok(self):
        cfg = self._minimal_config()
        cfg['species'][0]['transport_D0'] = 1e-9
        cfg['species'][0]['transport_alpha'] = 0.0
        issues = self._run(cfg)
        warnings = [i for i in issues if 'transport' in i.location]
        self.assertEqual(warnings, [], f"No transport warnings expected: {warnings}")

    def test_n3_invalid_bc_type_warning(self):
        cfg = self._minimal_config()
        cfg['species'][0]['bc_upper_type'] = 99
        issues = self._run(cfg)
        warnings = [i for i in issues if i.severity == 'WARNING' and 'bc_upper_type' in i.location]
        self.assertTrue(len(warnings) >= 1, f"Expected N-3 warning for bc_type=99: {issues}")
        self.assertIn('99', warnings[0].message)

    def test_n3_valid_bc_type_ok(self):
        cfg = self._minimal_config()
        cfg['species'][0]['bc_upper_type'] = 1
        issues = self._run(cfg)
        warnings = [i for i in issues if 'bc_upper_type' in i.location]
        self.assertEqual(warnings, [], f"No warning expected for valid bc_type=1: {warnings}")

    def test_n1_equilibrium_model_por0_no_warning(self):
        """Equilibrium model has por0=0.85 — must produce no por0 warning."""
        cfg = self._minimal_config()
        cfg['parameters']['physical']['por0'] = 0.85
        issues = self._run(cfg)
        warnings = [i for i in issues if 'por0' in i.location]
        self.assertEqual(warnings, [], f"por0=0.85 is valid, no warning: {warnings}")


if __name__ == '__main__':
    unittest.main()
