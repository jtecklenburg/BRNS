"""
ACG Orchestrator - Complete YAML to Fortran Pipeline

Orchestrates the full workflow from YAML configuration to generated Fortran code:
1. Load and validate YAML configuration
2. Evaluate formulas (FormulaEvaluator)
3. Map to ACG structures (YAMLtoACGMapper)
4. Pre-process stoichiometry (p0-p10 from gaussian_elimination)
5. Generate Fortran code (acg0-acg17b from acg)

Author: Jan Tecklenburg
Date: 2026-03-03
"""

from pathlib import Path
from typing import Dict, List, Any, Optional, Set, Tuple
from dataclasses import dataclass, field
import yaml
import os
import re

from sympy import Symbol, symbols, sympify, nsimplify

from acg_brns.formula_evaluator import (
    FormulaEvaluator,
    FormulaEvaluationError,
    CircularDependencyError,
    UndefinedVariableError,
    get_math_function_names,
    get_math_function_namespace,
)
from acg_brns.yaml_to_acg_mapper import YAMLtoACGMapper
from acg_brns.acg import ACGModule
from acg_brns.gaussian_elimination import (
    p0_initialize_old_variables,
    p1_initialize_reaction_lists,
    p2_create_equation_names,
    p3_reorder_reactions,
    run_gaussian_elimination
)


@dataclass
class ValidationIssue:
    """A single validation finding (error or warning)."""
    severity: str          # 'ERROR' or 'WARNING'
    path: str              # e.g. "reactions[2].stoichiometry.feii"
    message: str
    hint: str = ''

    @property
    def location(self) -> str:
        return self.path

    def __str__(self) -> str:
        base = f"{self.severity} [{self.path}]: {self.message}"
        return f"{base}\n  → {self.hint}" if self.hint else base


class ACGOrchestrationError(Exception):
    """Raised when orchestration fails."""
    pass


class ACGOrchestrator:
    """
    Orchestrates complete YAML→Fortran code generation pipeline.
    
    This class manages the entire workflow following proc0903-M.md:
    - Configuration loading and validation
    - Formula evaluation
    - Structure mapping
    - Pre-processing (p0-p10)
    - Code generation (acg0-acg17b)
    """
    
    def __init__(self, yaml_path: str, output_dir: str, verbose: bool = False):
        """
        Initialize orchestrator.
        
        Args:
            yaml_path: Path to YAML configuration file
            output_dir: Directory where Fortran files will be generated
            verbose: Enable verbose output
        """
        self.yaml_path = Path(yaml_path)
        self.output_dir = Path(output_dir)
        self.verbose = verbose
        
        # State variables
        self.config: Optional[Dict] = None
        self.evaluated_params: Optional[Dict] = None
        self.evaluator: Optional[FormulaEvaluator] = None
        self.mapper: Optional[YAMLtoACGMapper] = None
        self.acg: Optional[ACGModule] = None
        self.acg_data: Optional[Dict] = None
        self.reduced_system: Optional[Dict] = None
        self._validation_warnings: List[ValidationIssue] = []  # set by _validate_config()
        
        # Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        if self.verbose:
            print(f"ACGOrchestrator initialized:")
            print(f"  YAML: {self.yaml_path}")
            print(f"  Output: {self.output_dir}")
    
    def load_config(self) -> Dict:
        """
        Load and validate YAML configuration.
        
        Returns:
            Loaded configuration dictionary
            
        Raises:
            ACGOrchestrationError: If YAML is invalid or missing required sections
        """
        if self.verbose:
            print("\n=== Phase 1: Loading Configuration ===")
        
        if not self.yaml_path.exists():
            raise ACGOrchestrationError(f"YAML file not found: {self.yaml_path}")
        
        try:
            with open(self.yaml_path, 'r') as f:
                self.config = yaml.safe_load(f)
        except yaml.YAMLError as e:
            raise ACGOrchestrationError(f"Invalid YAML: {e}")

        if self.config is None:
            raise ACGOrchestrationError(
                "Invalid YAML: file is empty. Expected a top-level mapping with "
                "'species', 'reactions', and 'parameters'."
            )

        if not isinstance(self.config, dict):
            raise ACGOrchestrationError(
                "Invalid YAML: top-level document must be a mapping/object. "
                "Expected keys like 'species', 'reactions', and 'parameters'."
            )
        
        # Validate required sections
        self._validate_config()
        
        if self.verbose:
            print(f"✓ Loaded configuration from {self.yaml_path.name}")
            print(f"  Species: {len(self.config.get('species', []))}")
            print(f"  Reactions: {len(self.config.get('reactions', []))}")
        
        return self.config
    
    def _validate_config(self):
        """
        Validate YAML configuration in multiple phases.

        Collects all issues first, then raises a single ACGOrchestrationError
        with a human-readable report if any ERRORs were found.

        Raises:
            ACGOrchestrationError: If validation finds one or more ERRORs
        """
        issues: List[ValidationIssue] = []
        issues.extend(self._validate_schema())
        issues.extend(self._validate_references())
        issues.extend(self._validate_formulas())
        issues.extend(self._validate_plausibility())

        errors = [i for i in issues if i.severity == 'ERROR']
        warnings = [i for i in issues if i.severity == 'WARNING']
        self._validation_warnings = warnings  # expose for callers / tests

        if warnings and self.verbose:
            for w in warnings:
                print(f"  ⚠ {w}")

        if errors:
            parts = [f"YAML validation failed "
                     f"({len(errors)} errors, {len(warnings)} warnings):"]
            parts += [str(i) for i in issues]
            raise ACGOrchestrationError('\n\n'.join(parts))

    # ------------------------------------------------------------------
    # Phase 1 — Schema: required fields and allowed values
    # ------------------------------------------------------------------

    def _validate_schema(self) -> List[ValidationIssue]:
        """Check required sections, field presence, and allowed values."""
        issues: List[ValidationIssue] = []
        cfg = self.config
        params = cfg.get('parameters')

        # S-1..S-3: required top-level sections
        for section in ('species', 'reactions', 'parameters'):
            if section not in cfg:
                issues.append(ValidationIssue(
                    'ERROR', section,
                    f"Required section '{section}' is missing.",
                    f"Please add section '{section}' to the YAML file."
                ))

        # abort schema checks if fundamental sections missing
        if any(i.path in ('species', 'reactions', 'parameters') for i in issues):
            return issues

        if not isinstance(params, dict):
            issues.append(ValidationIssue(
                'ERROR', 'parameters',
                "'parameters' must be a mapping.",
                "Example: parameters: {biogeochemical: [], physical: {}, stoichiometry: {}}"
            ))
            params = {}

        if isinstance(params, dict):
            for section_name, section_data in params.items():
                section_path = f"parameters.{section_name}"
                if isinstance(section_data, dict):
                    continue
                if isinstance(section_data, list):
                    for i, item in enumerate(section_data):
                        item_path = f"{section_path}[{i}]"
                        if not isinstance(item, dict):
                            issues.append(ValidationIssue(
                                'ERROR', item_path,
                                "List-style parameter entries must be mappings.",
                                "Example: - {name: kox, value: 0.1}"
                            ))
                            continue
                        if 'name' not in item or not isinstance(item.get('name'), str):
                            issues.append(ValidationIssue(
                                'ERROR', f"{item_path}.name",
                                "'name' is required and must be a string.",
                                "Example: name: kox"
                            ))
                        if 'value' not in item:
                            issues.append(ValidationIssue(
                                'ERROR', f"{item_path}.value",
                                "'value' is required for list-style parameter entries.",
                                "Example: value: 0.1"
                            ))
                    continue

                issues.append(ValidationIssue(
                    'ERROR', section_path,
                    "Each parameters subsection must be either a mapping or a list.",
                    "Examples: physical: {por0: 0.8} or biogeochemical: [{name: kox, value: 0.1}]"
                ))

        optional_mapping_sections = {
            'grid': "Example: grid: {nnodes: 100, depth_max: 1.0}",
            'time': "Example: time: {total: 10.0, step: 0.1}",
            'output': "Example: output: {timing: {start: 1.0, interval: 1.0}}",
            'initial_conditions': "Example: initial_conditions: {mode: 3}",
            'advanced': "Example: advanced: {listnotransp: [], switches: []}",
            'transport': "Example: transport: {listnotransp: []}",
        }
        for section_name, hint in optional_mapping_sections.items():
            section_value = cfg.get(section_name)
            if section_value is not None and not isinstance(section_value, dict):
                issues.append(ValidationIssue(
                    'ERROR', section_name,
                    f"'{section_name}' must be a mapping.",
                    hint
                ))

        # species
        if not isinstance(cfg.get('species'), list):
            issues.append(ValidationIssue(
                'ERROR', 'species',
                "'species' must be a list.",
                "Example:\nspecies:\n  - name: o2\n    type: dissolved"
            ))
        else:
            allowed_types = {'dissolved', 'solid'}
            for i, sp in enumerate(cfg['species']):
                path = f"species[{i}]"
                if not isinstance(sp, dict):
                    issues.append(ValidationIssue('ERROR', path, "Entry must be a mapping.", ""))
                    continue
                if 'name' not in sp:
                    issues.append(ValidationIssue(
                        'ERROR', f"{path}.name",
                        "'name' is missing.",
                        "Each species requires a name field, e.g. 'name: o2'."
                    ))
                elif not isinstance(sp['name'], str):
                    issues.append(ValidationIssue(
                        'ERROR', f"{path}.name",
                        "'name' must be a string.",
                        "Example: name: o2"
                    ))
                if 'type' not in sp:
                    issues.append(ValidationIssue(
                        'ERROR', f"{path}.type",
                        "'type' is missing.",
                        f"Allowed values: {sorted(allowed_types)}."
                    ))
                elif not isinstance(sp['type'], str):
                    issues.append(ValidationIssue(
                        'ERROR', f"{path}.type",
                        "'type' must be a string.",
                        f"Allowed values: {sorted(allowed_types)}."
                    ))
                elif sp['type'] not in allowed_types:
                    issues.append(ValidationIssue(
                        'ERROR', f"{path}.type",
                        f"Unknown species type '{sp['type']}'.",
                        f"Allowed values: {sorted(allowed_types)}."
                    ))
                # S-8: scalar/formula fields must not be list or dict
                _scalar_fields = (
                    'bc_upper_value', 'bc_lower_value', 'init_value',
                    'transport_D0', 'transport_alpha', 'transport_tortuosity',
                )
                for fname in _scalar_fields:
                    val = sp.get(fname)
                    if isinstance(val, (list, dict)):
                        issues.append(ValidationIssue(
                            'ERROR', f"{path}.{fname}",
                            f"Field '{fname}' must be a number or formula string,"
                            f" not {type(val).__name__}.",
                            "Lists and mappings are not allowed here."
                        ))

        # reactions
        if not isinstance(cfg.get('reactions'), list):
            issues.append(ValidationIssue(
                'ERROR', 'reactions',
                "'reactions' must be a list.",
                ""
            ))
        else:
            for i, rxn in enumerate(cfg['reactions']):
                path = f"reactions[{i}]"
                if not isinstance(rxn, dict):
                    issues.append(ValidationIssue('ERROR', path, "Entry must be a mapping.", ""))
                    continue
                if 'id' not in rxn:
                    issues.append(ValidationIssue(
                        'ERROR', f"{path}.id",
                        "'id' is missing.",
                        "Each reaction requires a unique integer ID."
                    ))
                elif not isinstance(rxn['id'], int) or isinstance(rxn['id'], bool):
                    issues.append(ValidationIssue(
                        'ERROR', f"{path}.id",
                        "'id' must be an integer.",
                        "Example: id: 1"
                    ))
                if 'name' in rxn and not isinstance(rxn['name'], str):
                    issues.append(ValidationIssue(
                        'ERROR', f"{path}.name",
                        "'name' must be a string.",
                        "Example: name: aerobic_respiration"
                    ))
                if 'stoichiometry' not in rxn:
                    issues.append(ValidationIssue(
                        'ERROR', f"{path}.stoichiometry",
                        "'stoichiometry' is missing.",
                        "Provide stoichiometry as a mapping, e.g. 'stoichiometry: {o2: -1}'."
                    ))
                elif not isinstance(rxn['stoichiometry'], dict):
                    issues.append(ValidationIssue(
                        'ERROR', f"{path}.stoichiometry",
                        "'stoichiometry' must be a mapping.",
                        "Provide stoichiometry as a mapping, e.g. 'stoichiometry: {o2: -1}'."
                    ))
                # S-7: equilibrium reactions need a constraint
                if rxn.get('equilibrium') is True and 'equilibrium_constraint' not in rxn:
                    rid = rxn.get('id', i)
                    issues.append(ValidationIssue(
                        'ERROR', f"{path}.equilibrium_constraint",
                        f"Reaction {rid}: 'equilibrium: true' is set but 'equilibrium_constraint' is missing.",
                        "Please add 'equilibrium_constraint' (e.g. 'hco3 - keq1 * co2')."
                    ))

        return issues

    # ------------------------------------------------------------------
    # Phase 2 — References: unique names, known cross-references
    # ------------------------------------------------------------------

    @staticmethod
    def _extract_tokens(expr_str: str) -> set:
        """
        Extract identifier tokens from an expression string.

        Removes standalone numeric literals (including Fortran d-notation like
        1d0 and Python sci-notation like 1.5e-3) WITHOUT touching digits that
        are embedded inside identifiers (e.g. co2, hco3, fo2, keq1, F3, B3).

        The lookbehind ``(?<![A-Za-z0-9_])`` ensures only a digit sequence that is
        NOT preceded by a letter or underscore is treated as the start of a
        standalone number, so ``fo2`` and ``hco3`` are preserved intact while
        ``1d0`` and ``2.5e-4`` are removed entirely.
        """
        s = str(expr_str)
        # Remove standalone Fortran d-notation numbers (e.g. 1d0, 1.5d-3).
        # Lookbehind includes digits so that e.g. the '7' in 'sw17' (preceded
        # by '1') is NOT treated as the start of a standalone number.
        s = re.sub(r'(?<![A-Za-z0-9_])\d+\.?\d*[dD][+-]?\d+(?![A-Za-z0-9_])', '', s)
        # Remove standalone Python/C sci-notation numbers (e.g. 1.5e-3, 1E0).
        s = re.sub(r'(?<![A-Za-z0-9_])\d+\.?\d*[Ee][+-]?\d+(?![A-Za-z0-9_])', '', s)
        # Remove remaining standalone integer / float literals.
        s = re.sub(r'(?<![A-Za-z0-9_])\d+\.?\d*(?![A-Za-z0-9_])', '', s)
        return set(re.findall(r'[A-Za-z_][A-Za-z0-9_]*', s))

    @staticmethod
    def _extract_function_calls(expr_str: str) -> set:
        """Extract function-like identifiers used as name(...)."""
        return set(re.findall(r'\b([A-Za-z_][A-Za-z0-9_]*)\s*\(', str(expr_str)))

    @staticmethod
    def _iter_parameter_entries(params_cfg: Dict[str, Any]):
        """
        Iterate parameter entries across all `parameters.*` subsections.

        Yields tuples of:
            (section_name, parameter_name, parameter_value, value_path)
        """
        if not isinstance(params_cfg, dict):
            return

        for section_name, section_data in params_cfg.items():
            if isinstance(section_data, dict):
                for pname, pvalue in section_data.items():
                    if isinstance(pname, str):
                        yield section_name, pname, pvalue, f"parameters.{section_name}.{pname}"
            elif isinstance(section_data, list):
                for i, item in enumerate(section_data):
                    if not isinstance(item, dict):
                        continue
                    pname = item.get('name')
                    if isinstance(pname, str) and 'value' in item:
                        yield section_name, pname, item.get('value'), f"parameters.{section_name}[{i}].value"

    @staticmethod
    def _collect_parameter_names(params_cfg: Dict[str, Any]) -> Set[str]:
        """
        Collect parameter names across all `parameters.*` subsections.

        Supported per-subsection syntaxes:
        - Mapping style: key/value pairs
        - List style: entries with `name` and `value`
        """
        names: Set[str] = set()
        if not isinstance(params_cfg, dict):
            return names

        for section_data in params_cfg.values():
            if isinstance(section_data, dict):
                names |= {k for k in section_data.keys() if isinstance(k, str)}
            elif isinstance(section_data, list):
                for item in section_data:
                    if not isinstance(item, dict):
                        continue
                    pname = item.get('name')
                    if isinstance(pname, str) and 'value' in item:
                        names.add(pname)

        return names

    def _validate_references(self) -> List[ValidationIssue]:
        """Check uniqueness of names/IDs and forward references."""
        issues: List[ValidationIssue] = []
        cfg = self.config
        params = cfg.get('parameters', {})
        if not isinstance(params, dict):
            params = {}

        # Collect known species names
        species_list = cfg.get('species', [])
        if not isinstance(species_list, list):
            return issues  # already caught in schema phase

        species_names: List[str] = [
            sp['name'] for sp in species_list
            if isinstance(sp, dict) and isinstance(sp.get('name'), str)
        ]

        # R-1: species name uniqueness
        seen_species: set = set()
        for i, name in enumerate(species_names):
            if name in seen_species:
                issues.append(ValidationIssue(
                    'ERROR', f"species[{i}].name",
                    f"Species name '{name}' is not unique.",
                    "Use each name only once."
                ))
            seen_species.add(name)
        known_species = set(species_names)

        known_math_functions = get_math_function_names()

        # Build complete known-symbol set (species + all parameters subsections)
        known_params = self._collect_parameter_names(params)
        known_symbols = known_species | known_params

        # Add computed_values keys (all nested subsections) to known_symbols so that
        # R-7 accepts formula references like 'C_hco3' defined in computed_values.
        cv = cfg.get('computed_values', {})
        if isinstance(cv, dict):
            for _subsect in cv.values():
                if isinstance(_subsect, dict):
                    known_symbols |= {k for k in _subsect if isinstance(k, str)}

        # R-6: parameter name uniqueness across all parameters subsections
        seen_param_names: Set[str] = set()
        for _, pname, _, value_path in self._iter_parameter_entries(params):
            if pname in seen_param_names:
                issues.append(ValidationIssue(
                    'ERROR', value_path,
                    f"Parameter name '{pname}' is not unique across parameters subsections.",
                    "Use each parameter name only once in parameters.*."
                ))
            seen_param_names.add(pname)

        reactions = cfg.get('reactions', [])
        if not isinstance(reactions, list):
            return issues

        # Collect ALL rate_component keys across all reactions into a global
        # helper pool. In models like Canfield, helpers defined in one reaction
        # (e.g. kch2o, fo2) are referenced in other reactions' rate strings.
        all_rate_comp_names: set = set()
        for rxn in reactions:
            if isinstance(rxn, dict):
                rc = rxn.get('rate_components', {})
                if isinstance(rc, dict):
                    all_rate_comp_names |= set(rc.keys())
        known_symbols = known_symbols | all_rate_comp_names

        # R-2: reaction ID uniqueness
        seen_ids: set = set()
        for i, rxn in enumerate(reactions):
            if not isinstance(rxn, dict):
                continue
            rid = rxn.get('id')
            if rid is not None:
                if rid in seen_ids:
                    issues.append(ValidationIssue(
                        'ERROR', f"reactions[{i}].id",
                        f"Reaction ID {rid} is not unique.",
                        "Each reaction requires a unique ID."
                    ))
                seen_ids.add(rid)

        for i, rxn in enumerate(reactions):
            if not isinstance(rxn, dict):
                continue
            rid = rxn.get('id', i)
            path = f"reactions[{i}]"

            # R-3: stoichiometry references known species
            stoich = rxn.get('stoichiometry', {})
            if isinstance(stoich, dict):
                for sp_name in stoich:
                    if sp_name not in known_species:
                        issues.append(ValidationIssue(
                            'ERROR', f"{path}.stoichiometry.{sp_name}",
                            f"Unknown species '{sp_name}' in stoichiometry of reaction {rid}.",
                            f"Defined species: {sorted(known_species)}."
                        ))

            # R-4: equilibrium_constraint references known species
            ec = rxn.get('equilibrium_constraint')
            if ec is not None:
                call_names = self._extract_function_calls(ec)
                for call_name in call_names:
                    if call_name not in known_math_functions:
                        issues.append(ValidationIssue(
                            'WARNING', f"{path}.equilibrium_constraint",
                            f"Unknown function '{call_name}' in equilibrium_constraint of reaction {rid}.",
                            "Function is not in supported math whitelist and may fail during formula parsing/evaluation."
                        ))

                for token in self._extract_tokens(str(ec)):
                    if token in known_math_functions or token in call_names:
                        continue
                    if token not in known_symbols:
                        issues.append(ValidationIssue(
                            'ERROR', f"{path}.equilibrium_constraint",
                            f"Unknown symbol '{token}' in equilibrium_constraint of reaction {rid}.",
                            f"Known symbols: species + biogeochemical + stoichiometry + physical."
                        ))

            # R-5: rate / rate_components reference known symbols.
            # All rate_component keys (from all reactions) are in known_symbols
            # via the global all_rate_comp_names pool collected above.
            rate_components = rxn.get('rate_components', {})
            local_known = known_symbols

            rate_val = rxn.get('rate')
            if rate_val is not None and isinstance(rate_val, str):
                call_names = self._extract_function_calls(rate_val)
                for call_name in call_names:
                    if call_name not in known_math_functions:
                        issues.append(ValidationIssue(
                            'WARNING', f"{path}.rate",
                            f"Unknown function '{call_name}' in rate of reaction {rid}.",
                            "Function is not in supported math whitelist and may fail during formula parsing/evaluation."
                        ))
                for token in self._extract_tokens(rate_val):
                    if token in known_math_functions or token in call_names:
                        continue
                    if token not in local_known:
                        issues.append(ValidationIssue(
                            'ERROR', f"{path}.rate",
                            f"Unknown symbol '{token}' in rate of reaction {rid}.",
                            f"Symbol must be a species, bio-, stoich-, or physical-parameter."
                        ))

            if isinstance(rate_components, dict):
                for comp_name, comp_val in rate_components.items():
                    if comp_val is not None and isinstance(comp_val, str):
                        call_names = self._extract_function_calls(comp_val)
                        for call_name in call_names:
                            if call_name not in known_math_functions:
                                issues.append(ValidationIssue(
                                    'WARNING',
                                    f"{path}.rate_components.{comp_name}",
                                    f"Unknown function '{call_name}' in rate_components.{comp_name} of reaction {rid}.",
                                    "Function is not in supported math whitelist and may fail during formula parsing/evaluation."
                                ))
                        for token in self._extract_tokens(comp_val):
                            if token in known_math_functions or token in call_names:
                                continue
                            if token not in local_known:
                                issues.append(ValidationIssue(
                                    'ERROR',
                                    f"{path}.rate_components.{comp_name}",
                                    f"Unknown symbol '{token}' in rate_components.{comp_name} of reaction {rid}.",
                                    f"Symbol must be a species, bio-, stoich-, or physical-parameter."
                                ))

        # R-7: bc_upper_value / bc_lower_value / init_value formula string references
        for i, sp in enumerate(species_list):
            if not isinstance(sp, dict):
                continue
            sp_name = sp.get('name', f"?[{i}]")
            for fname in ('bc_upper_value', 'bc_lower_value', 'init_value'):
                val = sp.get(fname)
                if isinstance(val, str):
                    call_names = self._extract_function_calls(val)
                    for call_name in call_names:
                        if call_name not in known_math_functions:
                            issues.append(ValidationIssue(
                                'WARNING', f"species[{i}].{fname}",
                                f"Unknown function '{call_name}' in {fname} of species '{sp_name}'.",
                                "Function is not in supported math whitelist and may fail during formula parsing/evaluation."
                            ))
                    for token in self._extract_tokens(val):
                        if token in known_math_functions or token in call_names:
                            continue
                        if token not in known_symbols:
                            issues.append(ValidationIssue(
                                'ERROR', f"species[{i}].{fname}",
                                f"Unknown symbol '{token}' in {fname} of species '{sp_name}'.",
                                "Known symbols: species + parameter + computed_values."
                            ))

        return issues

    def _validate_formulas(self) -> List[ValidationIssue]:
        """Phase 3: Syntactic validation of all formula strings via sympify."""
        issues: List[ValidationIssue] = []
        cfg = self.config
        params = cfg.get('parameters', {})
        if not isinstance(params, dict):
            params = {}

        _math_ns = get_math_function_namespace()

        def _preprocess(s: str) -> str:
            s = re.sub(
                r'(?<![A-Za-z0-9_])(\d+\.?\d*)[dD]([+-]?\d+)(?![A-Za-z0-9_])',
                r'\1e\2', s)
            s = s.replace('^', '**')
            return s

        def _try_parse(expr_str: str):
            try:
                sympify(_preprocess(str(expr_str)), locals=_math_ns)
                return None
            except Exception as exc:
                return str(exc)[:200]

        # F-1: formula syntax in all parameters subsections
        param_formula_by_name: Dict[str, str] = {}
        param_formula_path_by_name: Dict[str, str] = {}
        all_param_names = self._collect_parameter_names(params)

        for _, pname, pval, value_path in self._iter_parameter_entries(params):
            if isinstance(pval, str):
                err = _try_parse(pval)
                if err:
                    issues.append(ValidationIssue(
                        'ERROR', value_path,
                        f"Formula for parameter '{pname}' is not syntactically parseable.",
                        err
                    ))

                if pname not in param_formula_by_name:
                    param_formula_by_name[pname] = pval
                    param_formula_path_by_name[pname] = value_path

        # F-2: detect circular dependencies among all parameter formulas
        deps: Dict[str, Set[str]] = {name: set() for name in param_formula_by_name.keys()}
        for pname, expr in param_formula_by_name.items():
            tokens = self._extract_tokens(expr)
            deps[pname] = {tok for tok in tokens if tok in all_param_names}

        visited: Dict[str, int] = {}  # 0=unseen, 1=visiting, 2=done
        stack: List[str] = []
        reported_cycles: Set[str] = set()

        def _dfs(node: str):
            visited[node] = 1
            stack.append(node)
            for nxt in deps.get(node, set()):
                state = visited.get(nxt, 0)
                if state == 0:
                    _dfs(nxt)
                elif state == 1:
                    # back-edge: cycle found in current stack
                    try:
                        idx = stack.index(nxt)
                        cycle = stack[idx:] + [nxt]
                    except ValueError:
                        cycle = [nxt, node, nxt]
                    cycle_key = '->'.join(cycle)
                    if cycle_key not in reported_cycles:
                        reported_cycles.add(cycle_key)
                        pretty = ' -> '.join(cycle)
                        issues.append(ValidationIssue(
                            'ERROR', param_formula_path_by_name.get(nxt, 'parameters'),
                            f"Circular dependency detected in parameter formulas: {pretty}.",
                            "Check dependency order and remove cyclic references."
                        ))
            stack.pop()
            visited[node] = 2

        for node in deps.keys():
            if visited.get(node, 0) == 0:
                _dfs(node)

        # F-3 & F-4: per-reaction formula strings
        for i, rxn in enumerate(cfg.get('reactions', [])):
            if not isinstance(rxn, dict):
                continue
            rid = rxn.get('id', i)
            ec = rxn.get('equilibrium_constraint')
            if isinstance(ec, str):
                err = _try_parse(ec)
                if err:
                    issues.append(ValidationIssue(
                        'ERROR', f"reactions[{i}].equilibrium_constraint",
                        f"equilibrium_constraint of reaction {rid} is not syntactically parseable.",
                        err
                    ))
            rate = rxn.get('rate')
            if isinstance(rate, str):
                err = _try_parse(rate)
                if err:
                    issues.append(ValidationIssue(
                        'ERROR', f"reactions[{i}].rate",
                        f"rate of reaction {rid} is not syntactically parseable.",
                        err
                    ))
            rc = rxn.get('rate_components', {})
            if isinstance(rc, dict):
                for comp_name, comp_val in rc.items():
                    if isinstance(comp_val, str):
                        err = _try_parse(comp_val)
                        if err:
                            issues.append(ValidationIssue(
                                'ERROR',
                                f"reactions[{i}].rate_components.{comp_name}",
                                f"rate_components.{comp_name} of reaction {rid} is not syntactically parseable.",
                                err
                            ))
        return issues

    def _validate_plausibility(self) -> List[ValidationIssue]:
        """Phase 4: Plausibility checks; emits warnings only."""
        issues: List[ValidationIssue] = []
        cfg = self.config
        params = cfg.get('parameters', {})
        if not isinstance(params, dict):
            params = {}

        # N-1: Porosity por0 ∈ (0, 1)
        physical = params.get('physical', {})
        if isinstance(physical, dict):
            por0 = physical.get('por0')
            if isinstance(por0, (int, float)) and not (0 < por0 < 1):
                issues.append(ValidationIssue(
                    'WARNING', 'parameters.physical.por0',
                    f"Porosity por0 = {por0} is outside (0, 1).",
                    "Porosity must be between 0 and 1."
                ))

        # N-2 & N-3: per-species transport parameters and bc types
        for i, sp in enumerate(cfg.get('species', [])):
            if not isinstance(sp, dict):
                continue
            sp_name = sp.get('name', f"?[{i}]")
            path = f"species[{i}]"
            for fname in ('transport_D0', 'transport_alpha'):
                val = sp.get(fname)
                if isinstance(val, (int, float)) and val < 0:
                    issues.append(ValidationIssue(
                        'WARNING', f"{path}.{fname}",
                        f"Transport parameter '{fname}' of '{sp_name}' is negative ({val}).",
                        f"'{fname}' should be >= 0."
                    ))
            _valid_bc = {0, 1, 2}
            for fname in ('bc_upper_type', 'bc_lower_type'):
                val = sp.get(fname)
                if isinstance(val, int) and val not in _valid_bc:
                    issues.append(ValidationIssue(
                        'WARNING', f"{path}.{fname}",
                        f"Boundary condition type '{fname}' = {val} for '{sp_name}' is unknown.",
                        f"Allowed values: {sorted(_valid_bc)} (0=Dirichlet, 1=Neumann/zero-flux, 2=...)."
                    ))
        return issues

    def evaluate_formulas(self) -> Dict:
        """
        Evaluate all formulas using FormulaEvaluator.

        Returns:
            Dictionary of evaluated parameters
        """
        if self.verbose:
            print("\n=== Phase 2: Evaluating Formulas ===")

        try:
            self.evaluator = FormulaEvaluator(self.config, verbose=self.verbose)
            self.evaluated_params = self.evaluator.evaluate_all()
        except (FormulaEvaluationError, CircularDependencyError, UndefinedVariableError) as exc:
            raise ACGOrchestrationError(f"Formula evaluation failed: {exc}") from exc

        if self.verbose:
            n_evaluated = len([
                v for v in self.evaluated_params.values()
                if isinstance(v, (int, float))
            ])
            print(f"✓ Evaluated {n_evaluated} parameters")

        return self.evaluated_params

    def _get_mapping_section(self, name: str, default: Optional[Dict] = None) -> Dict:
        """Return a top-level config section as a mapping, or a safe default."""
        if default is None:
            default = {}
        if not isinstance(self.config, dict):
            return dict(default)
        value = self.config.get(name, default)
        return value if isinstance(value, dict) else dict(default)

    def _get_parameters_mapping(self) -> Dict:
        """Return the validated parameters mapping, or an empty mapping."""
        return self._get_mapping_section('parameters', {})

    def _normalize_switch_condition(self, condition: str) -> str:
        """
        Normalize a YAML switch condition into Fortran-compatible syntax.

        Supported forms:
        - Maple style numeric condition: "30.0 - x(j)"
        - Comparison operators: "x(j) < 30.0", "x(j) >= 1.0e-2"
        - Fortran operators are passed through: ".lt.", ".ge.", etc.
        - Logical words are normalized: and/or/not -> .and./.or./.not.
        """
        if not isinstance(condition, str) or not condition.strip():
            raise ACGOrchestrationError(
                "Invalid switch condition: expected non-empty string"
            )

        normalized = condition.strip()

        # Normalize common logical operators to Fortran style.
        normalized = re.sub(r'\band\b', '.and.', normalized, flags=re.IGNORECASE)
        normalized = re.sub(r'\bor\b', '.or.', normalized, flags=re.IGNORECASE)
        normalized = re.sub(r'\bnot\b', '.not.', normalized, flags=re.IGNORECASE)

        # Translate comparison operators to Fortran form. Apply longer operators
        # first so e.g. '<=' is not split by the '<' replacement.
        replacements = [
            (r'\s*<=\s*', '.le.'),
            (r'\s*>=\s*', '.ge.'),
            (r'\s*==\s*', '.eq.'),
            (r'\s*!=\s*', '.ne.'),
            (r'\s*<\s*', '.lt.'),
            (r'\s*>\s*', '.gt.'),
        ]
        for pattern, repl in replacements:
            normalized = re.sub(pattern, repl, normalized)

        return normalized

    def _extract_physical_codegen_data(self) -> Tuple[List, List, List, List]:
        """Build physical parameter arrays for acg0/acg8 consistently."""
        params_cfg = self._get_parameters_mapping()
        physical = params_cfg.get('physical', {})
        if not isinstance(physical, dict):
            physical = {}
        phys_names = list(physical.keys())
        phys_vals = []

        # Prefer evaluator outputs so physical formulas like "0.39*0.86" are
        # already numeric by code-generation time.
        evaluated = self.evaluated_params if isinstance(self.evaluated_params, dict) else {}

        for pname, raw_val in physical.items():
            val = evaluated.get(pname, raw_val)

            if isinstance(val, (int, float)):
                phys_vals.append(float(val))
                continue

            if isinstance(val, str):
                # Fast path: plain numeric strings
                try:
                    phys_vals.append(float(val))
                    continue
                except (ValueError, TypeError):
                    pass

                # Fallback: evaluate the expression string in current context
                if self.evaluator is not None:
                    try:
                        evaluated_val = self.evaluator.evaluate_formula(val, evaluated)
                        phys_vals.append(float(evaluated_val))
                        continue
                    except Exception as e:
                        raise ACGOrchestrationError(
                            f"Physical parameter '{pname}' with value '{val}' "
                            f"could not be evaluated to float: {e}"
                        )

            raise ACGOrchestrationError(
                f"Physical parameter '{pname}' has unsupported/non-numeric value '{val}' "
                f"(type {type(val).__name__})."
            )

        phys_names2 = ['iq', 'iw', 'iDb', 'ipor', 'igrid', 'iarea', 'ic']
        phys_flags = params_cfg.get('physical_flags', {})
        if not isinstance(phys_flags, dict):
            phys_flags = {}

        phys_vals2 = []
        for name in phys_names2:
            if name == 'ic':
                phys_vals2.append(int(self._get_initial_mode()))
            else:
                phys_vals2.append(int(phys_flags.get(name, 0)))

        return phys_names, phys_vals, phys_names2, phys_vals2

    def map_to_acg_structures(self) -> Dict:
        """
        Map YAML to ACG-compatible data structures.

        Returns:
            Dictionary with ACG data structures
        """
        if self.verbose:
            print("\n=== Phase 3: Mapping to ACG Structures ===")

        self.mapper = YAMLtoACGMapper(self.config, self.evaluated_params)

        bio_name, bio_val = self.mapper.build_bio_arrays()
        variables = self.mapper.build_variables_list()
        species_map = self.mapper.build_species_index_map()
        stoich_matrix = self.mapper.build_stoichiometry_matrix(self.evaluator)
        reactions = self.mapper.build_rate_expressions()

        dissolved = [s for s in self.mapper.species_list if s.get('type') == 'dissolved']
        solid = [s for s in self.mapper.species_list if s.get('type') == 'solid']
        ordered_reaction_ids = [reaction['id'] for reaction in self.mapper._get_ordered_reactions()]

        self.acg_data = {
            'bio_name': bio_name,
            'bio_val': bio_val,
            'variables': variables,
            'species_map': species_map,
            'stoich_matrix': stoich_matrix,
            'reactions': reactions,
            'reaction_ids': ordered_reaction_ids,
            'ndissolved': len(dissolved),
            'nsolids': len(solid),
            'ncompo': len(variables),
            'nreactions': len(reactions)
        }

        if self.verbose:
            print("✓ Mapped structures:")
            print(f"  Variables: {self.acg_data['ncompo']}")
            print(f"    Dissolved: {self.acg_data['ndissolved']}")
            print(f"    Solid: {self.acg_data['nsolids']}")
            print(f"  Reactions: {self.acg_data['nreactions']}")
            print(f"  Bio-parameters: {len(bio_name)}")

        return self.acg_data

    def run_preprocessing(self) -> Dict:
        """
        Execute pre-processing functions p0-p10.

        These functions prepare the stoichiometric system for code generation
        by performing Gaussian elimination and computing residuals/Jacobian.

        Returns:
            Reduced system from Gaussian elimination
        """
        if self.verbose:
            print("\n=== Phase 4: Pre-processing (p0-p10) ===")

        ncompo = self.acg_data['ncompo']
        nreactions = self.acg_data['nreactions']

        var_symbols = [Symbol(v) for v in self.acg_data['variables']]
        var_symbols_old = [Symbol(f"{v}_old") for v in self.acg_data['variables']]
        reaction_symbols = [Symbol(f"r{i+1}") for i in range(nreactions)]

        species_index = {name: idx for idx, name in enumerate(self.acg_data['variables'])}
        ordered_reactions = (
            self.mapper._get_ordered_reactions()
            if self.mapper is not None
            else sorted(self.config.get('reactions', []), key=lambda r: r.get('id', 10**9))
        )

        local_symbols = {str(v): Symbol(str(v)) for v in self.acg_data.get('variables', [])}
        for name in self.acg_data.get('bio_name', []):
            local_symbols[str(name)] = Symbol(str(name))

        reaction_expr_by_id = {}
        for rxn in self.acg_data.get('reactions', []):
            rid = rxn.get('id')
            rate_text = str(rxn.get('rate_expanded') or rxn.get('rate_yaml') or '0')
            rate_text = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', rate_text)
            try:
                reaction_expr_by_id[rid] = sympify(rate_text, locals=local_symbols)
            except Exception as exc:
                raise ACGOrchestrationError(
                    f"Failed to parse mapped rate expression for reaction {rid}: {rate_text}. {exc}"
                ) from exc

        reaction_values = []
        for reaction in ordered_reactions:
            if len(reaction_values) >= nreactions:
                break
            rid = reaction.get('id')
            reaction_values.append(reaction_expr_by_id.get(rid, sympify('0')))

        numeric_subs = {}
        params_cfg = self._get_parameters_mapping()
        stoich_cfg = params_cfg.get('stoichiometry', {})
        if not isinstance(stoich_cfg, dict):
            stoich_cfg = {}
        for name in ('x', 'y', 'z'):
            value = stoich_cfg.get(name, None)
            if value is None:
                value = self.evaluated_params.get(name, None) if self.evaluated_params else None
            if isinstance(value, (int, float)):
                numeric_subs[Symbol(name)] = nsimplify(value, rational=True)

        equations = [sympify('0') for _ in range(ncompo)]
        for j, reaction in enumerate(ordered_reactions):
            if j >= nreactions:
                break
            stoich_dict = reaction.get('stoichiometry', {}) or {}
            for sp_name, coeff in stoich_dict.items():
                i = species_index.get(sp_name)
                if i is None:
                    continue
                coeff_text = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', str(coeff))
                try:
                    c = sympify(coeff_text, locals=local_symbols)
                except Exception as exc:
                    raise ACGOrchestrationError(
                        f"Failed to parse stoichiometric coefficient '{coeff}' for species '{sp_name}' "
                        f"in reaction {reaction.get('id', j)}. {exc}"
                    ) from exc
                equations[i] += c * reaction_symbols[j]

        if numeric_subs:
            equations = [sympify(eq).subs(numeric_subs) for eq in equations]

        eq_reactions = [r for r in self.acg_data['reactions'] if r.get('equilibrium')]
        neqrxns = len(eq_reactions)
        eqrxn_ids = []
        ordered_reaction_ids = self.acg_data.get('reaction_ids', [])
        for eq_rxn in eq_reactions:
            rid = eq_rxn.get('id')
            if rid in ordered_reaction_ids:
                eqrxn_ids.append(ordered_reaction_ids.index(rid) + 1)

        equilibrium_eqns = []
        if neqrxns > 0:
            sym_map = {str(v): Symbol(str(v)) for v in self.acg_data['variables']}
            for pname in self.acg_data.get('bio_name', []):
                sym_map[str(pname)] = Symbol(str(pname))
            for pname in self._collect_parameter_names(params_cfg):
                sym_map[str(pname)] = Symbol(str(pname))
            if isinstance(self.evaluated_params, dict):
                for pname in self.evaluated_params.keys():
                    sym_map.setdefault(str(pname), Symbol(str(pname)))

            raw_constraints_by_id = {}
            for raw_rxn in (self.config.get('reactions', []) if self.config else []):
                rid = raw_rxn.get('id')
                if rid is not None and 'equilibrium_constraint' in raw_rxn:
                    raw_constraints_by_id[rid] = str(raw_rxn['equilibrium_constraint'])

            for eq_rxn in eq_reactions:
                rid = eq_rxn.get('id')
                raw = raw_constraints_by_id.get(rid)
                if raw is not None:
                    constraint_str = raw
                elif 'equilibrium_constraint' in eq_rxn:
                    constraint_str = str(eq_rxn['equilibrium_constraint'])
                else:
                    continue
                constraint_str = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', constraint_str)
                try:
                    if '=' in constraint_str:
                        lhs, rhs = constraint_str.split('=', 1)
                        expr = sympify(lhs, locals=sym_map) - sympify(rhs, locals=sym_map)
                    else:
                        expr = sympify(constraint_str, locals=sym_map)
                    equilibrium_eqns.append(expr)
                except Exception as exc:
                    raise ACGOrchestrationError(
                        f"Failed to parse equilibrium constraint for reaction {rid}: {constraint_str}. {exc}"
                    ) from exc

        delt = Symbol('delt')

        if self.verbose:
            print("  Running Gaussian elimination...")

        self.reduced_system = run_gaussian_elimination(
            equations=equations,
            variables=var_symbols,
            reactions=reaction_symbols,
            reaction_values=reaction_values,
            variables_old=var_symbols_old,
            delt=delt,
            eqrxn_ids=eqrxn_ids,
            neqrxns=neqrxns,
            equilibrium_eqns=equilibrium_eqns if equilibrium_eqns else None,
            inert_components=None,
            verbose=self.verbose
        )

        substi_sp = {}
        for i, v in enumerate(var_symbols, 1):
            substi_sp[v] = Symbol(f"sp({i},j)")
        for i, v_old in enumerate(var_symbols_old, 1):
            substi_sp[v_old] = Symbol(f"spold({i},j)")

        self.reduced_system['func'] = [
            sympify(f).subs(substi_sp) if hasattr(f, 'subs') else f
            for f in self.reduced_system.get('func', [])
        ]

        jac = self.reduced_system.get('jacobian', None)
        if jac is not None and hasattr(jac, 'applyfunc'):
            self.reduced_system['jacobian'] = jac.applyfunc(
                lambda x: sympify(x).subs(substi_sp) if hasattr(x, 'subs') else x
            )

        if self.verbose:
            print("✓ Pre-processing complete")
            print(f"  Residuals: {len(self.reduced_system['func'])}")
            print(f"  Jacobian: {self.reduced_system['jacobian'].shape}")

        return self.reduced_system
    
    def run_code_generation(self) -> Dict[str, Any]:
        """
        Execute ACG code generation functions (acg0-acg18).
        
        Generates all Fortran subroutines following proc0903-M.md sequence.
        
        Returns:
            Summary of generated files
        """
        if self.verbose:
            print("\n=== Phase 5: Code Generation (acg0-acg17b) ===")
        
        # Initialize ACG module
        self.acg = ACGModule(str(self.output_dir))
        
        # Extract parameters for code generation
        params = self._extract_generation_parameters()
        
        ncompo = self.acg_data['ncompo']
        nreactions = self.acg_data['nreactions']
        ndissolved = self.acg_data['ndissolved']
        nsolids = self.acg_data['nsolids']
        nnodes = params.get('nnodes', 100)
        
        # Model structure (acg0)
        if self.verbose:
            print("  Generating model structure...")
        phys_names, phys_vals, phys_names2, phys_vals2 = self._extract_physical_codegen_data()
        self.acg.acg0(nsolids, ndissolved, nreactions, nnodes,
                     self.acg_data['bio_name'], self.acg_data['bio_val'],
                     phys_names, phys_vals, phys_names2, phys_vals2)
        
        # Boundary conditions (acg1)
        if self.verbose:
            print("  Generating boundary conditions...")
        self.acg.acg1(
            ncompo,
            params['type_up'],
            params['bnddata_up'],
            params['type_down'],
            params['bnddata_down']
        )
        
        # Diffusion (acg2)
        if self.verbose:
            print("  Generating diffusion...")
        self.acg.acg2(ncompo, params['diffdata'], params['alphadata'])
        
        # Biogeochemical parameters (acg3)
        if self.verbose:
            print("  Generating biogeochemical parameters...")
        self.acg.acg3(self.acg_data['bio_name'], self.acg_data['bio_val'])
        
        # Residual function (acg4)
        if self.verbose:
            print("  Generating residual function...")
        self.acg.acg4(ncompo, self.reduced_system['func'])
        
        # Jacobian (acg5)
        if self.verbose:
            print("  Generating Jacobian...")
        jacobian_list = self.reduced_system['jacobian'].tolist()
        self.acg.acg5(ncompo, jacobian_list, self.acg_data['variables'])
        
        # Output variables (acg7 - Maple: acg7)
        if self.verbose:
            print("  Generating output specifications...")
        noutput = params.get('noutput', 0)
        nroutput = params.get('nroutput', 0)
        output_indices = params.get('listoutput', [])
        routput_indices = params.get('listroutput', [])
        file_names = params.get('file_names', [])
        file_rnames = params.get('file_rnames', [])
        time_iniout = params.get('time_iniout', 1.0)
        time_intvout = params.get('time_intvout', 1.0)
        self.acg.acg7(noutput, nroutput, output_indices, routput_indices,
                     file_names, file_rnames, time_iniout, time_intvout)
        
        # Physical parameters (acg8 - Maple: acg8)
        if self.verbose:
            print("  Generating physical parameters...")
        phys_names, phys_vals, phys_names2, phys_vals2 = self._extract_physical_codegen_data()
        self.acg.acg8(phys_names, phys_vals, phys_names2, phys_vals2)
        
        # Initial conditions (acg12 - Maple: acg12)
        if self.verbose:
            print("  Generating initial conditions...")
        listinput = params.get('listinput', [])
        file_in_names = params.get('file_in_names', [])
        self.acg.acg12(ncompo, params['vic'], params['iniconc'], listinput, file_in_names)
        
        # Net reaction rates (acg13 - Maple: acg13)
        if self.verbose:
            print("  Generating net reaction rates...")
        species_symbols = [Symbol(v) for v in self.acg_data['variables']]
        net_rates = self._build_net_reaction_rates(species_symbols)
        # Maple-compatible substitution: species -> sp(i,j)
        subst_dict = {
            Symbol(v): Symbol(f'sp({i+1},j)')
            for i, v in enumerate(self.acg_data['variables'])
        }
        self.acg.acg13(ncompo, net_rates, species_symbols, subst_dict)
        
        # No transport list (acg14 - Maple: acg14)
        if self.verbose:
            print("  Generating no-transport list...")
        listnotransp = self._build_no_transport_list()
        self.acg.acg14(listnotransp=listnotransp)
        
        # Reaction rates (acg15 - Maple: acg15)
        if self.verbose:
            print("  Generating reaction rates...")
        rate_exprs = [r['rate_fortran'] for r in self.acg_data['reactions']]
        self.acg.acg15(nreactions, rate_exprs)
        
        # Solid species (acg16 - Maple: acg16)
        if self.verbose:
            print("  Generating solid species identification...")
        solid_indices = self.mapper.get_solid_indices()
        self.acg.acg16(self.acg_data['nsolids'], solid_indices)
        
        # Terminal electron acceptor cascade (acg17 - Maple: acg17, optional)
        # In Maple: if (iswitch=1) then acg17(variables, dir_f) fi;
        use_teac = params.get('iswitch', 0) == 1
        if use_teac:
            if self.verbose:
                print("  Generating terminal electron acceptor cascade...")
            self.acg.acg17(self.acg_data['variables'])
        
        # Spatial switches (acg17a - Maple: acg17a)
        # Do not overwrite a previously generated TEAC `switches.f` with an
        # empty stub. Only generate custom switches when explicit conditions are
        # provided and no TEAC switches file was requested.
        advanced_cfg = self._get_mapping_section('advanced', {})
        switch_entries = advanced_cfg.get('switches', [])
        if not isinstance(switch_entries, list):
            switch_entries = []
        switch_conditions = {}
        switch_names = []
        for entry in switch_entries:
            if not isinstance(entry, dict):
                continue
            name = entry.get('name')
            condition = entry.get('condition')
            if name and condition:
                condition = self._normalize_switch_condition(condition)
                switch_names.append(name)
                switch_conditions[name] = condition

        if not use_teac:
            if self.verbose:
                print("  Generating spatial switches (stub)..." if not switch_conditions else "  Generating spatial switches...")
            self.acg.acg17a(switch_conditions, switch_names)
        
        # Parameter switches (acg17b - Maple: acg17b)
        if self.verbose:
            print("  Generating parameter array...")
        parameter_names = self.acg_data.get('bio_name', [])
        self.acg.acg17b(parameter_names=parameter_names)

        # Time step initialization include (acg18)
        if self.verbose:
            print("  Generating timestep initialization include...")
        self.acg.acg18(**params['timestep_parameters'])
        
        if self.verbose:
            print(f"✓ Code generation complete")
            print(f"  Output directory: {self.output_dir}")
        
        return {
            'output_dir': str(self.output_dir),
            'files_generated': self._list_generated_files()
        }
    
    def _extract_generation_parameters(self) -> Dict[str, Any]:
        """
        Extract parameters from YAML for ACG code generation.
        
        Returns:
            Dictionary with all parameters needed for acg functions
        """
        params = {}
        params_cfg = self._get_parameters_mapping()
        
        # Boundary conditions
        bc_data = self.mapper.get_boundary_conditions(self.evaluator)
        params['type_up'] = [bc['bc_upper_type'] for bc in bc_data.values()]
        params['bnddata_up'] = [bc['bc_upper_value'] for bc in bc_data.values()]
        params['type_down'] = [bc['bc_lower_type'] for bc in bc_data.values()]
        params['bnddata_down'] = [bc['bc_lower_value'] for bc in bc_data.values()]
        
        # Transport parameters
        transport = self.mapper.get_transport_parameters()
        params['diffdata'] = [t['D0'] for t in transport.values()]
        params['alphadata'] = [t['alpha'] for t in transport.values()]
        
        # Initial conditions
        init_data = self.mapper.get_initial_conditions(self.evaluator)
        global_ic_mode = self._get_initial_mode()
        if init_data:
            # Maple-compatible behavior: use one global mode (vic/ic),
            # not per-species init_mode.
            params['vic'] = global_ic_mode
            params['iniconc'] = [ic['init_value'] for ic in init_data.values()]
        else:
            params['vic'] = global_ic_mode
            params['iniconc'] = [0.0] * self.acg_data['ncompo']
        
        # Grid parameters
        grid_cfg = self._get_mapping_section('grid', {})
        params['nnodes'] = grid_cfg.get('nnodes', 100)
        params['depth_max'] = grid_cfg.get('depth_max', 1.0)
        
        # Time stepping
        time_cfg = self._get_mapping_section('time', {})
        params['tot_time'] = time_cfg.get('total', 1.0)
        params['delt'] = time_cfg.get('step', 0.1)

        # Timestep initialization parameters (advanced.timestep_parameters)
        params['timestep_parameters'] = self.mapper.get_timestep_parameters()
        
        # Physical/environment aliases (prefer evaluated generic parameters)
        physical = params_cfg.get('physical', {})
        if not isinstance(physical, dict):
            physical = {}
        eval_params = self.evaluated_params if isinstance(self.evaluated_params, dict) else {}

        t_raw = eval_params.get('T_C', eval_params.get('t_celsius', physical.get('T_C', physical.get('temperature', 20.0))))
        s_raw = eval_params.get('S', eval_params.get('salin', physical.get('S', physical.get('salinity', 35.0))))

        try:
            params['T_C'] = float(t_raw)
        except (TypeError, ValueError):
            params['T_C'] = 20.0

        try:
            params['S'] = float(s_raw)
        except (TypeError, ValueError):
            params['S'] = 35.0
        
        # Output configuration - derived from per-species/reaction output flags
        output_cfg = self._get_mapping_section('output', {})
        timing = output_cfg.get('timing', {})
        if not isinstance(timing, dict):
            timing = {}
        params['time_iniout'] = timing.get('start', 1.0)
        params['time_intvout'] = timing.get('interval', 1.0)

        # Species outputs: collect species with output: true (preserve definition order)
        output_sp_indices = []
        output_file_names = []
        for sp in self.mapper.species_list:
            if sp.get('output', False):
                output_sp_indices.append(self.mapper.species_map[sp['name']])
                output_file_names.append(sp.get('output_filename', sp['name']))
        params['noutput'] = len(output_sp_indices)
        params['listoutput'] = output_sp_indices
        params['file_names'] = output_file_names

        # Reaction outputs: collect reactions with output: true, sorted by id (Maple order)
        output_r_indices = []
        output_rfile_names = []
        for rxn in sorted(self.config.get('reactions', []), key=lambda r: r['id']):
            if rxn.get('output', False):
                output_r_indices.append(rxn['id'])
                output_rfile_names.append(rxn.get('output_filename', f"rate{rxn['id']}"))
        params['nroutput'] = len(output_r_indices)
        params['listroutput'] = output_r_indices
        params['file_rnames'] = output_rfile_names
        
        # Input file handling for global initial-condition mode ic=3
        init_cfg = self._get_mapping_section('initial_conditions', {})
        if global_ic_mode == 3:
            explicit_listinput = init_cfg.get('listinput')
            explicit_files = init_cfg.get('file_in_names')

            if explicit_listinput is not None and explicit_files is not None:
                params['listinput'] = explicit_listinput
                params['file_in_names'] = explicit_files
            else:
                params['listinput'] = list(range(1, self.acg_data['ncompo'] + 1))
                params['file_in_names'] = [species['name'] for species in self.mapper.species_list]
        else:
            params['listinput'] = []
            params['file_in_names'] = []
        
        return params

    def _get_initial_mode(self) -> int:
        """
        Return global initial-condition mode (Maple vic/ic semantics).

        Preferred source: top-level initial_conditions.mode
        Backward-compatible fallback: parameters.physical_flags.ic
        """
        if self.config is None:
            return 3

        init_cfg = self._get_mapping_section('initial_conditions', {})
        mode = init_cfg.get('mode', None)
        if mode is None:
            if isinstance(self.evaluated_params, dict):
                mode = self.evaluated_params.get('ic', None)
        if mode is None:
            params_cfg = self._get_parameters_mapping()
            phys_flags = params_cfg.get('physical_flags', {})
            if not isinstance(phys_flags, dict):
                phys_flags = {}
            mode = phys_flags.get('ic', 3)

        try:
            return int(mode)
        except (TypeError, ValueError):
            return 3
    
    def _build_net_reaction_rates(self, species_symbols: List[Symbol]) -> List:
        """
        Build Maple-compatible net reaction rates for each species.

        Maple `acg13` uses the original stoichiometric system (dX/dt equations),
        not the transformed residuals from p8/p9 (which include old-state and delt).
        Therefore we reconstruct:

            dC_i/dt = sum_j stoich[j, i] * rate_j

        where `rate_j` are the (expanded) kinetic rate expressions.
        
        Returns:
            List of SymPy expressions for net rates
        """
        ncompo = self.acg_data.get('ncompo', 0)
        reactions = self.acg_data.get('reactions', [])

        if ncompo == 0 or not reactions:
            return self.reduced_system.get('func', [])

        # Symbols available in kinetic expressions: species + parameters
        local_symbols = {str(s): s for s in species_symbols}
        for name in self.acg_data.get('bio_name', []):
            local_symbols[str(name)] = Symbol(str(name))

        # Parse each reaction rate to SymPy
        rate_exprs = []
        for rxn in reactions:
            rate_text = str(rxn.get('rate_expanded') or rxn.get('rate_yaml') or '0')
            # SymPy expects E-notation instead of Fortran D-notation
            rate_text = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', rate_text)
            try:
                rate_exprs.append(sympify(rate_text, locals=local_symbols))
            except Exception as exc:
                raise ACGOrchestrationError(
                    f"Failed to parse net-rate expression for reaction {rxn.get('id')}: {rate_text}. {exc}"
                ) from exc

        # Build SYMBOLIC stoichiometric matrix from raw reaction stoichiometry.
        # Important: do NOT use evaluated numeric stoich_matrix here; Maple acg13
        # keeps symbolic factors (e.g., SD, x, y) in ssrates.f.
        ordered_reactions = (
            self.mapper._get_ordered_reactions()
            if self.mapper is not None
            else sorted(self.config.get('reactions', []), key=lambda r: r.get('id', 10**9))
        )

        species_index = {name: idx for idx, name in enumerate(self.acg_data.get('variables', []))}
        stoich_sym = [[sympify('0') for _ in range(ncompo)] for _ in range(len(ordered_reactions))]

        for r_idx, reaction in enumerate(ordered_reactions):
            stoich_dict = reaction.get('stoichiometry', {}) or {}
            for sp_name, coeff in stoich_dict.items():
                s_idx = species_index.get(sp_name)
                if s_idx is None:
                    continue
                coeff_text = str(coeff)
                coeff_text = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', coeff_text)
                try:
                    stoich_sym[r_idx][s_idx] = sympify(coeff_text, locals=local_symbols)
                except Exception as exc:
                    raise ACGOrchestrationError(
                        f"Failed to parse net-rate stoichiometric coefficient '{coeff}' for species '{sp_name}' "
                        f"in reaction {reaction.get('id', r_idx)}. {exc}"
                    ) from exc

        # dC_i/dt = sum_j stoich[j, i] * rate_j
        # Maple reference resolves stoichiometric constants x,y,z numerically
        # in ssrates.f, while keeping SD symbolic.
        # Use exact Integer/Rational for same reason as run_preprocessing().
        numeric_subs = {}
        params_cfg = self._get_parameters_mapping()
        stoich_cfg = params_cfg.get('stoichiometry', {}) if self.config else {}
        if not isinstance(stoich_cfg, dict):
            stoich_cfg = {}
        for name in ('x', 'y', 'z'):
            value = stoich_cfg.get(name, None)
            if value is None:
                value = self.evaluated_params.get(name, None) if self.evaluated_params else None
            if isinstance(value, (int, float)):
                numeric_subs[Symbol(name)] = nsimplify(value, rational=True)

        net_rates = []
        nreactions = min(len(rate_exprs), len(stoich_sym))
        for i in range(ncompo):
            expr = sympify('0')
            for j in range(nreactions):
                expr += stoich_sym[j][i] * rate_exprs[j]
            if numeric_subs:
                expr = expr.subs(numeric_subs)
            net_rates.append(expr)

        return net_rates
    
    def _build_no_transport_list(self) -> List[str]:
        """
        Build list of species without transport.

        Maple-compatible behavior:
        - Uses explicit `listnotransp` only
        - Does NOT implicitly add all solids or D0==0 species
          (those are handled via `issolid.f` / diffusion settings)
        
        Returns:
            List of species indices or names for acg14()
        """
        # Preferred location: advanced.listnotransp (Maple naming retained)
        # Optional fallback: transport.listnotransp
        advanced_cfg = self._get_mapping_section('advanced', {})
        raw_list = advanced_cfg.get('listnotransp', None)
        if raw_list is None:
            transport_cfg = self._get_mapping_section('transport', {})
            raw_list = transport_cfg.get('listnotransp', [])

        if not isinstance(raw_list, list):
            return []

        # Normalize entries to species names or integer indices accepted by acg14
        species_map = self.acg_data.get('species_map', {})
        normalized = []
        for item in raw_list:
            if isinstance(item, int):
                normalized.append(item)
            elif isinstance(item, str):
                # Allow either species name ("ch2o") or numeric string ("7")
                if item.isdigit():
                    normalized.append(int(item))
                elif item in species_map:
                    normalized.append(item)

        # De-duplicate while preserving order
        deduped = []
        seen = set()
        for item in normalized:
            if item not in seen:
                deduped.append(item)
                seen.add(item)

        return deduped
    
    def _list_generated_files(self) -> List[str]:
        """
        List all Fortran files generated in output directory.
        
        Returns:
            List of generated file names
        """
        if not self.output_dir.exists():
            return []
        
        fortran_files = []
        for file in self.output_dir.iterdir():
            if file.suffix in ['.f', '.f90']:
                fortran_files.append(file.name)
        
        return sorted(fortran_files)
    
    def generate(self) -> Dict[str, Any]:
        """
        Execute complete YAML→Fortran pipeline.
        
        Main entry point that runs all phases in sequence.
        
        Returns:
            Summary dictionary with generation statistics
            
        Raises:
            ACGOrchestrationError: If any phase fails
        """
        try:
            # Phase 1: Load configuration
            self.load_config()
            
            # Phase 2: Evaluate formulas
            self.evaluate_formulas()
            
            # Phase 3: Map to ACG structures
            self.map_to_acg_structures()
            
            # Phase 4: Pre-processing
            self.run_preprocessing()
            
            # Phase 5: Code generation
            generation_result = self.run_code_generation()
            
            # Build summary
            summary = {
                'success': True,
                'yaml_file': str(self.yaml_path),
                'output_dir': str(self.output_dir),
                'n_species': self.acg_data['ncompo'],
                'n_dissolved': self.acg_data['ndissolved'],
                'n_solid': self.acg_data['nsolids'],
                'n_reactions': self.acg_data['nreactions'],
                'n_parameters': len(self.acg_data['bio_name']),
                'files_generated': generation_result['files_generated']
            }
            
            if self.verbose:
                print("\n" + "="*60)
                print("PIPELINE COMPLETE")
                print("="*60)
                print(f"Species: {summary['n_species']} "
                      f"({summary['n_dissolved']} dissolved + {summary['n_solid']} solid)")
                print(f"Reactions: {summary['n_reactions']}")
                print(f"Parameters: {summary['n_parameters']}")
                print(f"Files: {len(summary['files_generated'])}")
                print(f"Output: {summary['output_dir']}")
            
            return summary
            
        except ACGOrchestrationError:
            raise
        except Exception as e:
            if self.verbose:
                print(f"\n✗ Pipeline failed: {e}")
            raise ACGOrchestrationError(f"Pipeline failed: {e}") from e


def main():
    """Example usage with Canfield model."""
    import sys
    
    # Default paths
    yaml_path = 'models/equilibrium/equilibrium.yaml'
    output_dir = 'build/fortran/equilibrium'
    
    # Parse command line arguments
    if len(sys.argv) > 1:
        yaml_path = sys.argv[1]
    if len(sys.argv) > 2:
        output_dir = sys.argv[2]
    
    # Create orchestrator
    orchestrator = ACGOrchestrator(
        yaml_path=yaml_path,
        output_dir=output_dir,
        verbose=True
    )
    
    # Generate Fortran code
    try:
        summary = orchestrator.generate()
        
        print("\n" + "="*60)
        print("SUCCESS")
        print("="*60)
        print(f"Generated {len(summary['files_generated'])} Fortran files")
        print(f"Location: {summary['output_dir']}")
        
    except ACGOrchestrationError as e:
        print(f"\nERROR: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
