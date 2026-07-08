"""
Formula Evaluator for YAML Model Configurations

This module evaluates formulas in YAML configuration files to concrete numeric values,
replicating Maple's behavior of pre-calculating all formula values before passing them
to ACG functions.

Author: Jan Tecklenburg
Date: 2026-03-03
"""

from typing import Dict, Set, List, Any, Union
import sympy as sp
import re
from collections import defaultdict, deque


MATH_FUNCTION_NAMESPACE = {
    'sqrt': sp.sqrt,
    'log': sp.log,
    'log10': lambda x: sp.log(x, 10),
    'ln': sp.log,
    'exp': sp.exp,
    'sin': sp.sin,
    'cos': sp.cos,
    'tan': sp.tan,
    'asin': sp.asin,
    'acos': sp.acos,
    'atan': sp.atan,
    'sinh': sp.sinh,
    'cosh': sp.cosh,
    'tanh': sp.tanh,
    'abs': sp.Abs,
    'ceil': sp.ceiling,
    'floor': sp.floor,
    'min': sp.Min,
    'max': sp.Max,
}

MATH_FUNCTION_NAMES = frozenset(MATH_FUNCTION_NAMESPACE.keys())


def get_math_function_names() -> Set[str]:
    """Return the supported formula-function names."""
    return set(MATH_FUNCTION_NAMES)


def get_math_function_namespace() -> Dict[str, Any]:
    """Return SymPy locals for supported formula functions."""
    return dict(MATH_FUNCTION_NAMESPACE)


class FormulaEvaluationError(Exception):
    """Raised when formula evaluation fails"""
    pass


class CircularDependencyError(Exception):
    """Raised when circular dependencies are detected"""
    pass


class UndefinedVariableError(Exception):
    """Raised when a variable is referenced but not defined"""
    pass


class DependencyResolver:
    """
    Analyzes formula strings to extract variable dependencies and
    performs topological sorting.
    """
    
    @staticmethod
    def extract_dependencies(formula: str, known_symbols: Set[str]) -> Set[str]:
        """
        Parse formula to find referenced variables.
        
        Args:
            formula: String expression (e.g., "10^(-mlogkp1)")
            known_symbols: Set of all defined variable names
            
        Returns:
            Set of variable names used in formula
        """
        if not isinstance(formula, str):
            # Numeric literal, no dependencies
            return set()
        
        # Find all potential variable names (alphanumeric + underscore)
        potential_vars = re.findall(r'\b[a-zA-Z_][a-zA-Z0-9_]*\b', formula)
        
        # Filter to only known symbols (exclude math functions)
        math_functions = MATH_FUNCTION_NAMES
        
        dependencies = set()
        for var in potential_vars:
            if var in known_symbols and var not in math_functions:
                dependencies.add(var)
        
        return dependencies
    
    @staticmethod
    def topological_sort(dep_graph: Dict[str, Set[str]]) -> List[str]:
        """
        Kahn's algorithm for topological sort.
        
        Args:
            dep_graph: {variable: set_of_dependencies}
            
        Returns:
            List of variables in evaluation order (dependencies first)
            
        Raises:
            CircularDependencyError: If circular dependencies detected
        """
        # Build in-degree map and adjacency list
        in_degree = defaultdict(int)
        adj_list = defaultdict(list)
        all_nodes = set(dep_graph.keys())
        
        # Add all dependencies to node set
        for deps in dep_graph.values():
            all_nodes.update(deps)
        
        # Initialize in-degrees
        for node in all_nodes:
            in_degree[node] = 0
        
        # Build adjacency list and count in-degrees
        for node, deps in dep_graph.items():
            for dep in deps:
                adj_list[dep].append(node)
                in_degree[node] += 1
        
        # Find all nodes with in-degree 0
        queue = deque([node for node in all_nodes if in_degree[node] == 0])
        result = []
        
        while queue:
            node = queue.popleft()
            result.append(node)
            
            # Reduce in-degree for neighbors
            for neighbor in adj_list[node]:
                in_degree[neighbor] -= 1
                if in_degree[neighbor] == 0:
                    queue.append(neighbor)
        
        # Check for cycles
        if len(result) != len(all_nodes):
            # Find remaining nodes (part of cycle)
            remaining = [node for node in all_nodes if in_degree[node] > 0]
            raise CircularDependencyError(
                f"Circular dependency detected involving: {remaining}"
            )
        
        return result


class FormulaEvaluator:
    """
    Evaluates YAML formulas to concrete numeric values.
    Replicates Maple's pre-calculation behavior.
    """
    
    def __init__(self, yaml_config: dict, verbose: bool = False):
        """
        Args:
            yaml_config: Parsed YAML configuration dictionary
            verbose: Emit evaluation warnings to stdout
        """
        self.config = yaml_config
        self.verbose = verbose
        self.symbol_table = {}  # {name: value_or_formula_string}
        self.evaluated = {}     # {name: concrete_numeric_value}
        self.dep_resolver = DependencyResolver()
        self.required_parameters = set()
        self.warnings = []

    def _warn(self, message: str):
        """Record a warning and optionally emit it to stdout."""
        self.warnings.append(message)
        if self.verbose:
            print(f"⚠️  Warning: {message}")

    def _iter_parameter_entries(self):
        """
        Iterate all parameter entries independent of subsection naming.

        Supported syntaxes in any `parameters.<subsection>`:
        - Mapping style: `name: value`
        - List style: `- name: ...; value: ...`
        """
        params_cfg = self.config.get('parameters', {})
        if not isinstance(params_cfg, dict):
            return

        for section_name, section_data in params_cfg.items():
            if isinstance(section_data, dict):
                for pname, pvalue in section_data.items():
                    if isinstance(pname, str):
                        yield section_name, pname, pvalue
            elif isinstance(section_data, list):
                for item in section_data:
                    if not isinstance(item, dict):
                        continue
                    pname = item.get('name')
                    if isinstance(pname, str) and 'value' in item:
                        yield section_name, pname, item.get('value')

    def _sync_derived_aliases(self):
        """Backfill common aliases without overriding explicit user parameters."""
        # Temperature aliases
        if 't_celsius' in self.evaluated and 'T_C' not in self.symbol_table:
            self.symbol_table['T_C'] = self.evaluated['t_celsius']
            self.evaluated['T_C'] = self.evaluated['t_celsius']
        if 'T_C' in self.evaluated and 't_celsius' not in self.symbol_table:
            self.symbol_table['t_celsius'] = self.evaluated['T_C']
            self.evaluated['t_celsius'] = self.evaluated['T_C']

        if 'T' not in self.symbol_table:
            if 'T_C' in self.evaluated:
                self.symbol_table['T'] = self.evaluated['T_C'] + 273.15
                self.evaluated['T'] = self.evaluated['T_C'] + 273.15
            elif 't_celsius' in self.evaluated:
                self.symbol_table['T'] = self.evaluated['t_celsius'] + 273.15
                self.evaluated['T'] = self.evaluated['t_celsius'] + 273.15

        # Salinity alias
        if 'salin' in self.evaluated and 'S' not in self.symbol_table:
            self.symbol_table['S'] = self.evaluated['salin']
            self.evaluated['S'] = self.evaluated['salin']
        if 'S' in self.evaluated and 'salin' not in self.symbol_table:
            self.symbol_table['salin'] = self.evaluated['S']
            self.evaluated['salin'] = self.evaluated['S']
        
    def load_base_parameters(self):
        """
        Load parameters with concrete numeric values.
        These serve as the base of the dependency graph.
        """
        for _, name, value in self._iter_parameter_entries():
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                self.symbol_table[name] = value
                self.evaluated[name] = float(value)
            elif isinstance(value, str):
                self.symbol_table[name] = value
                self.required_parameters.add(name)

        self._sync_derived_aliases()
    
    def load_computed_formulas(self):
        """
        Load formulas from computed_values section.
        These are typically complex expressions that depend on base parameters.
        """
        if 'computed_values' not in self.config:
            return
        
        computed = self.config['computed_values']
        
        # Iterate through all subsections (carbonate_system, sulfide_system, etc.)
        for section_name, section_data in computed.items():
            if section_name == 'description':
                continue
            
            if isinstance(section_data, dict):
                for key, value in section_data.items():
                    if isinstance(value, (int, float)):
                        # Numeric constant
                        self.symbol_table[key] = value
                        self.evaluated[key] = float(value)
                    elif isinstance(value, str):
                        # Skip function definitions and Maple-specific notation
                        if any(skip in value for skip in ['=', 'D(']):
                            continue
                        self.symbol_table[key] = value
    
    def load_formula_parameters(self):
        """
        Backward-compatible no-op.

        Formula strings from all `parameters.*` sections are now loaded in
        load_base_parameters() and evaluated in one unified pass.
        """
        return
    
    def load_stoichiometry_parameters(self):
        """
        Backward-compatible no-op.

        Stoichiometry parameters are part of unified `parameters.*` loading in
        load_base_parameters().
        """
        return
    
    def build_dependency_graph(self) -> Dict[str, Set[str]]:
        """
        Analyze formulas to build dependency graph.
        
        Returns:
            {variable: {dependencies}}
        """
        dep_graph = {}
        known_symbols = set(self.symbol_table.keys())
        
        for name, value in self.symbol_table.items():
            if isinstance(value, str):
                # It's a formula, extract dependencies
                deps = self.dep_resolver.extract_dependencies(value, known_symbols)
                dep_graph[name] = deps
            else:
                # Already evaluated, no dependencies
                dep_graph[name] = set()
        
        return dep_graph
    
    def evaluate_formula(self, formula: str, context: Dict[str, float]) -> float:
        """
        Evaluate a single formula using sympy.
        
        Args:
            formula: String expression (e.g., "10^(-mlogkp1)" or "1.0e-8.5")
            context: Available variables {name: value}
            
        Returns:
            Concrete numeric value
            
        Raises:
            FormulaEvaluationError: If evaluation fails
        """
        try:
            # Handle scientific notation with decimal exponents (e.g., 1.0e-8.5)
            # Convert to proper mathematical form: 1.0 * 10^(-8.5)
            formula_python = re.sub(
                r'(\d+\.?\d*)e([+-]?\d+\.?\d+)',
                r'(\1 * 10**(\2))',
                formula,
                flags=re.IGNORECASE
            )
            
            # Replace ^ with ** for Python/sympy
            formula_python = formula_python.replace('^', '**')
            
            # Create sympy symbols for variables in context
            symbols = {name: sp.Symbol(name) for name in context.keys()}
            
            # Add supported math functions
            local_dict = {
                **get_math_function_namespace(),
                **symbols,
            }
            
            # Parse expression
            expr = sp.sympify(formula_python, locals=local_dict)
            
            # Substitute values
            result = expr.subs(context)
            
            # Check if result still contains symbols (undefined variables)
            if result.free_symbols:
                undefined = [str(s) for s in result.free_symbols]
                raise FormulaEvaluationError(
                    f"Cannot evaluate '{formula}': undefined variables {undefined}"
                )
            
            # Evaluate to float
            return float(result.evalf())
            
        except Exception as e:
            raise FormulaEvaluationError(
                f"Failed to evaluate formula '{formula}': {e}"
            )
    
    def evaluate_all(self) -> Dict[str, float]:
        """
        Main evaluation pipeline:
        1. Load all parameters and formulas
        2. Build dependency graph
        3. Topological sort
        4. Evaluate in order
        
        Returns:
            {name: concrete_value} for all variables
            
        Raises:
            CircularDependencyError: If circular dependencies detected
            UndefinedVariableError: If undefined variable referenced
        """
        # Step 1: Load all data
        self.load_base_parameters()
        self.load_computed_formulas()
        self.load_formula_parameters()
        self.load_stoichiometry_parameters()
        
        # Step 2: Build dependency graph
        dep_graph = self.build_dependency_graph()
        
        # Step 3: Topological sort
        evaluation_order = self.dep_resolver.topological_sort(dep_graph)
        
        # Step 4: Evaluate in order
        skipped_params = []
        for name in evaluation_order:
            if name in self.evaluated:
                # Already evaluated (base parameter)
                continue
            
            if name not in self.symbol_table:
                self._warn(f"Variable '{name}' referenced but not defined in YAML")
                skipped_params.append((name, "not defined"))
                continue
            
            value = self.symbol_table[name]
            
            if isinstance(value, str):
                try:
                    # Evaluate formula
                    result = self.evaluate_formula(value, self.evaluated)
                    self.evaluated[name] = result
                    self._sync_derived_aliases()
                except FormulaEvaluationError as e:
                    detail = str(e).split(':', 1)[-1].strip()
                    if name in self.required_parameters:
                        raise FormulaEvaluationError(
                            f"Required biogeochemical parameter '{name}' could not be evaluated: {detail}"
                        ) from e
                    # Skip formulas that can't be evaluated (likely reference formulas)
                    # These are often documentation-only formulas like pH_from_alkalinity
                    self._warn(f"Cannot evaluate '{name}': {detail}")
                    skipped_params.append((name, "evaluation failed"))
                    continue
            else:
                # Should already be in self.evaluated, but handle just in case
                self.evaluated[name] = float(value)
        
        # Summary of skipped parameters
        if skipped_params and self.verbose:
            print(f"\n📊 Summary: {len(skipped_params)} parameter(s) could not be evaluated:")
            for param_name, reason in skipped_params:
                print(f"   - {param_name}: {reason}")
        
        return self.evaluated.copy()
    
    def get_species_values(self, species_list: List[dict]) -> Dict[str, Dict[str, float]]:
        """
        Evaluate species initial and boundary condition values.
        
        Args:
            species_list: List of species dictionaries from YAML
            
        Returns:
            {species_name: {
                'bc_upper_value': float,
                'bc_lower_value': float,
                'init_value': float
            }}
        """
        if not self.evaluated:
            self.evaluate_all()
        
        species_values = {}
        unevaluated_references = []
        
        for species in species_list:
            name = species['name']
            values = {}
            
            # Boundary conditions
            for bc_key in ['bc_upper_value', 'bc_lower_value']:
                if bc_key in species:
                    bc_val = species[bc_key]
                    if isinstance(bc_val, str):
                        # Reference to computed value
                        if bc_val in self.evaluated:
                            values[bc_key] = self.evaluated[bc_val]
                        else:
                            # Try to evaluate as formula
                            try:
                                values[bc_key] = self.evaluate_formula(bc_val, self.evaluated)
                            except FormulaEvaluationError:
                                # Cannot evaluate, leave as string for later runtime calculation
                                values[bc_key] = bc_val
                                unevaluated_references.append(f"{name}.{bc_key} = '{bc_val}'")
                    else:
                        values[bc_key] = float(bc_val)
            
            # Initial value
            if 'init_value' in species:
                init_val = species['init_value']
                if isinstance(init_val, str):
                    if init_val in self.evaluated:
                        values['init_value'] = self.evaluated[init_val]
                    else:
                        try:
                            values['init_value'] = self.evaluate_formula(init_val, self.evaluated)
                        except FormulaEvaluationError:
                            # Cannot evaluate, leave as string
                            values['init_value'] = init_val
                            unevaluated_references.append(f"{name}.init_value = '{init_val}'")
                else:
                    values['init_value'] = float(init_val)
            
            species_values[name] = values
        
        # Warn about unevaluated references
        if unevaluated_references:
            for ref in unevaluated_references:
                self.warnings.append(f"Species value could not be pre-computed: {ref}")

        if unevaluated_references and self.verbose:
            print(f"\n⚠️  Warning: {len(unevaluated_references)} species value(s) could not be pre-computed:")
            for ref in unevaluated_references:
                print(f"   - {ref}")
            print("   These will need to be computed at runtime.")
        
        return species_values
    
    def evaluate_stoichiometry(self, stoich_dict: Dict[str, Union[str, float]]) -> Dict[str, float]:
        """
        Evaluate stoichiometry coefficients.
        
        Args:
            stoich_dict: {species: coefficient_or_formula}
            
        Returns:
            {species: evaluated_coefficient}
        """
        if not self.evaluated:
            self.evaluate_all()
        
        evaluated_stoich = {}
        
        for species, coeff in stoich_dict.items():
            if isinstance(coeff, str):
                evaluated_stoich[species] = self.evaluate_formula(coeff, self.evaluated)
            else:
                evaluated_stoich[species] = float(coeff)
        
        return evaluated_stoich


def main():
    """
    Example usage with Canfield model.
    """
    import yaml
    import os
    
    # Load YAML configuration
    script_dir = os.path.dirname(os.path.abspath(__file__))
    yaml_path = os.path.join(script_dir, '..', 'models', 'equilibrium', 'equilibrium.yaml')
    with open(yaml_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Create evaluator
    evaluator = FormulaEvaluator(config)
    
    # Evaluate all formulas
    print("Evaluating formulas...")
    evaluated = evaluator.evaluate_all()
    
    # Print some key results
    print("\nKey Parameters:")
    for key in ['T', 'mlogkp1', 'k1', 'keq1', 'kbcarb', 'SD']:
        if key in evaluated:
            print(f"  {key} = {evaluated[key]:.6e}")
    
    # Evaluate species values
    print("\nSpecies Boundary/Initial Values:")
    species_values = evaluator.get_species_values(config['species'])
    for species in ['hco3', 'co3', 'co2', 'hplus']:
        if species in species_values:
            vals = species_values[species]
            print(f"  {species}:")
            for key, val in vals.items():
                if isinstance(val, str):
                    print(f"    {key} = '{val}' (not evaluated)")
                else:
                    print(f"    {key} = {val:.6e}")
    
    # Evaluate a reaction's stoichiometry
    print("\nReaction 1 Stoichiometry:")
    reaction = config['reactions'][0]
    stoich = evaluator.evaluate_stoichiometry(reaction['stoichiometry'])
    for species, coeff in stoich.items():
        print(f"  {species}: {coeff:.6e}")


if __name__ == '__main__':
    main()
