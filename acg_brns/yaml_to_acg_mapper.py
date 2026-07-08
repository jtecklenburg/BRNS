"""
YAML to ACG Mapper

Converts YAML configuration data into ACG-compatible data structures
for Fortran code generation.

Author: Jan Tecklenburg
Date: 2026-03-03
"""

from typing import Dict, List, Tuple, Any, Optional
import numpy as np
import re


class YAMLtoACGMapper:
    """
    Maps YAML configuration to ACG data structures.
    
    Takes evaluated parameters from FormulaEvaluator and YAML config,
    produces ACG-compatible arrays and data structures (bio_name, bio_val,
    variables list, stoichiometry matrix, rate expressions).
    """
    
    def __init__(self, yaml_config: dict, evaluated_params: dict):
        """
        Initialize mapper with YAML config and evaluated parameters.
        
        Args:
            yaml_config: Parsed YAML configuration dictionary
            evaluated_params: Results from FormulaEvaluator.evaluate_all()
        """
        self.config = yaml_config
        self.params = evaluated_params
        self.species_list = []
        self.species_map = {}
        self.variables = []
        
        # Initialize species data
        self._initialize_species()
    
    def _initialize_species(self):
        """Initialize species list and index mapping."""
        if 'species' in self.config:
            species = self.config['species']
            
            # Use species in the order they appear in YAML
            # (YAML must define species in correct Maple/Fortran order)
            self.species_list = species
            
            # Build variables list (just names)
            self.variables = [s['name'] for s in self.species_list]
            
            # Build 1-indexed species map (Fortran convention)
            self.species_map = {
                name: idx + 1 
                for idx, name in enumerate(self.variables)
            }

    def _get_ordered_reactions(self) -> List[Dict[str, Any]]:
        """
        Return reactions in canonical Maple-compatible order.

        Maple's `rate1..rateN` numbering is tied to reaction IDs, not the
        incidental order in which reactions appear in YAML. The whole Python
        pipeline must therefore use the same stable ordering.
        """
        reactions = self.config.get('reactions', [])
        if not reactions:
            return []

        if all('id' in reaction for reaction in reactions):
            return sorted(reactions, key=lambda reaction: reaction['id'])

        return list(reactions)

    def _iter_parameter_entries(self):
        """
        Iterate all parameter entries across first-level `parameters.*` sections.

        Supported per-section syntaxes:
        - mapping style: `name: value`
        - list style: `- name: ...; value: ...`

        Yields:
            Tuple[str, Any]: (parameter_name, raw_value)
        """
        params_cfg = self.config.get('parameters', {})
        if not isinstance(params_cfg, dict):
            return

        for section_data in params_cfg.values():
            if isinstance(section_data, dict):
                for name, value in section_data.items():
                    if isinstance(name, str):
                        yield name, value
            elif isinstance(section_data, list):
                for item in section_data:
                    if not isinstance(item, dict):
                        continue
                    name = item.get('name')
                    if isinstance(name, str) and 'value' in item:
                        yield name, item.get('value')

    def _get_parameter_names(self) -> set:
        """
        Return names of declared parameters across all first-level sections.

        Subsection labels under `parameters` are treated as user-level
        structure only; parameter semantics are defined by the entries.
        """
        return {name for name, _ in self._iter_parameter_entries()}

    def _build_global_rate_components(self) -> Dict[str, str]:
        """
        Collect reusable rate components across all reactions.

        YAML models often define helper expressions such as `kch2o`, `fo2`,
        `fno3`, `fmn4`, `ffe3`, `fso4` locally inside individual reactions,
        while other reactions reference them later. Maple resolves those helper
        symbols globally before generating Fortran, so we mirror that behavior.

        Components that are also declared model parameters (for example
        `sw17`, `sw18`, `sw19`) must NOT be inlined, otherwise the generated
        Fortran hardcodes them and diverges from Maple/reference behavior.
        """
        parameter_names = self._get_parameter_names()
        components: Dict[str, str] = {}

        for reaction in self._get_ordered_reactions():
            for component_name, component_expr in reaction.get('rate_components', {}).items():
                if component_name in parameter_names:
                    continue
                components[component_name] = str(component_expr)

        return components
    
    def build_bio_arrays(self) -> Tuple[List[str], List[float]]:
        """
        Build bio_name and bio_val arrays for acg3() (bioparams).
        
        Returns:
            Tuple of (bio_name, bio_val):
                - bio_name: List of parameter names
                - bio_val: List of parameter values (floats)
        """
        bio_name = []
        bio_val = []
        
        if 'parameters' not in self.config:
            return bio_name, bio_val
        
        for name, val in self._iter_parameter_entries():
            bio_name.append(name)

            if name in self.params:
                bio_val.append(float(self.params[name]))
            elif isinstance(val, (int, float)):
                # Allow explicit numeric literals from YAML only.
                # Formula-based parameters must be evaluated before mapping.
                bio_val.append(float(val))
            else:
                raise ValueError(
                    f"Parameter '{name}' has no evaluated numeric value. "
                    f"Original value: {val!r}"
                )
        
        return bio_name, bio_val
    
    def build_variables_list(self) -> List[str]:
        """
        Build ordered species list (variables).

        The mapper preserves YAML definition order exactly. This order
        determines Fortran array indexing sp(i,j), and some models rely on
        interleaved dissolved/solid species positions.
        
        Returns:
            List of species names in YAML/Fortran array order
        """
        return self.variables.copy()
    
    def build_species_index_map(self) -> Dict[str, int]:
        """
        Build species name -> Fortran index mapping.
        
        Returns:
            Dictionary mapping species name to 1-indexed Fortran array index
        """
        return self.species_map.copy()
    
    def get_solid_indices(self) -> List[int]:
        """
        Get list of solid species indices (1-indexed).
        
        Returns:
            List of Fortran indices for solid species
        """
        solid_indices = []
        for species in self.species_list:
            if species.get('type') == 'solid':
                idx = self.species_map[species['name']]
                solid_indices.append(idx)
        return solid_indices
    
    def build_stoichiometry_matrix(self, evaluator=None) -> np.ndarray:
        """
        Build stoichiometry coefficient matrix.
        
        Matrix has shape (n_reactions, n_species) where:
        - Rows represent reactions
        - Columns represent species
        - Values are stoichiometric coefficients
        
        Args:
            evaluator: Optional FormulaEvaluator instance for evaluating
                      stoichiometry expressions. If None, assumes all
                      coefficients are numeric.
        
        Returns:
            NumPy array of shape (n_reactions, n_species)
        """
        if 'reactions' not in self.config:
            return np.array([])
        
        reactions = self._get_ordered_reactions()
        n_reactions = len(reactions)
        n_species = len(self.variables)
        
        stoich_matrix = np.zeros((n_reactions, n_species))
        
        for r_idx, reaction in enumerate(reactions):
            if 'stoichiometry' not in reaction:
                continue
            
            stoich_dict = reaction['stoichiometry']
            
            # Evaluate stoichiometry coefficients
            if evaluator:
                evaluated_stoich = evaluator.evaluate_stoichiometry(stoich_dict)
            else:
                # Convert to float if possible
                evaluated_stoich = {}
                for species_name, coeff in stoich_dict.items():
                    if isinstance(coeff, (int, float)):
                        evaluated_stoich[species_name] = float(coeff)
                    else:
                        # Cannot evaluate, skip
                        evaluated_stoich[species_name] = 0.0
            
            # Fill matrix
            for species_name, coeff in evaluated_stoich.items():
                if species_name in self.species_map:
                    s_idx = self.species_map[species_name] - 1  # 0-indexed for numpy
                    stoich_matrix[r_idx, s_idx] = coeff
        
        return stoich_matrix
    
    def substitute_species_in_expression(self, expr: str) -> str:
        """
        Replace species names with sp(i,j) in expression.
        
        Uses word boundary matching to avoid partial replacements.
        Processes species in reverse length order to handle overlapping names.
        
        Args:
            expr: Expression string (e.g., "o2/kmo2")
            
        Returns:
            Converted expression (e.g., "sp(1,j)/kmo2")
        """
        result = expr
        
        # Sort species by length (longest first) to handle overlapping names
        # e.g., "h2s" before "hs"
        species_names = sorted(self.species_map.keys(), key=len, reverse=True)
        
        for species_name in species_names:
            # Match whole words only (word boundaries)
            pattern = r'\b' + re.escape(species_name) + r'\b'
            replacement = f'sp({self.species_map[species_name]},j)'
            result = re.sub(pattern, replacement, result)
        
        return result
    
    def expand_rate_components(self, rate: str, components: Dict[str, str]) -> str:
        """
        Expand rate_components into rate expression.
        
        Handles nested components (components referencing other components).
        Uses topological sort to expand in correct order.
        
        Args:
            rate: Rate expression with component references
            components: Dictionary of component_name -> component_expression
            
        Returns:
            Fully expanded rate expression
        """
        # Convert rate to string if needed
        rate = str(rate)
        
        if not components:
            return rate
        
        # Build dependency graph for components
        comp_deps = {}
        for comp_name, comp_expr in components.items():
            comp_expr_str = str(comp_expr)  # Convert to string
            deps = set()
            for other_comp in components.keys():
                pattern = r'\b' + re.escape(other_comp) + r'\b'
                if other_comp != comp_name and re.search(pattern, comp_expr_str):
                    deps.add(other_comp)
            comp_deps[comp_name] = deps
        
        # Topological sort
        sorted_comps = self._topological_sort_components(comp_deps)

        # Fallback for pathological graphs: ensure every component appears at
        # least once in processing order.
        for comp_name in components.keys():
            if comp_name not in sorted_comps:
                sorted_comps.append(comp_name)
        
        # Expand components in dependency order
        expanded = {}
        for comp_name, comp_expr in components.items():
            expanded[comp_name] = str(comp_expr)
        
        # Expand helper expressions transitively until fixed-point
        max_passes = max(1, len(sorted_comps) * 2)
        for _ in range(max_passes):
            changed = False
            for comp_name in sorted_comps:
                comp_expr = expanded[comp_name]
                for other_comp in sorted_comps:
                    if other_comp != comp_name:
                        pattern = r'\b' + re.escape(other_comp) + r'\b'
                        new_expr = re.sub(pattern, f'({expanded[other_comp]})', comp_expr)
                        if new_expr != comp_expr:
                            changed = True
                            comp_expr = new_expr
                expanded[comp_name] = comp_expr
            if not changed:
                break
        
        # Substitute expanded components into rate
        result = rate
        for _ in range(max_passes):
            changed = False
            for comp_name in sorted_comps:
                pattern = r'\b' + re.escape(comp_name) + r'\b'
                new_result = re.sub(pattern, f'({expanded[comp_name]})', result)
                if new_result != result:
                    changed = True
                    result = new_result
            if not changed:
                break
        
        return result
    
    def _topological_sort_components(self, deps: Dict[str, set]) -> List[str]:
        """
        Topological sort of rate components.
        
        Args:
            deps: Dictionary of component -> set of dependencies
            
        Returns:
            List of components in dependency order
        """
        from collections import deque
        
        # Calculate in-degrees (# of dependencies for each node)
        in_degree = {node: 0 for node in deps.keys()}
        for node, node_deps in deps.items():
            in_degree[node] = len([dep for dep in node_deps if dep in in_degree])
        
        # Find nodes with no dependencies
        queue = deque([node for node, degree in in_degree.items() if degree == 0])
        result = []
        
        while queue:
            node = queue.popleft()
            result.append(node)
            
            # Reduce in-degree for dependents
            for other_node in deps.keys():
                if node in deps[other_node]:
                    in_degree[other_node] -= 1
                    if in_degree[other_node] == 0:
                        queue.append(other_node)
        
        return result
    
    def build_rate_expressions(self) -> List[Dict[str, Any]]:
        """
        Build rate expressions with species substitution.
        
        For each reaction:
        1. Expand rate_components
        2. Substitute species names with sp(i,j)
        3. Keep parameter names (they're in bio_val)
        
        Returns:
            List of reaction dictionaries with:
                - id: Reaction ID
                - name: Reaction name
                - rate_yaml: Original YAML rate
                - rate_expanded: Expanded with components
                - rate_fortran: Fortran-ready with sp(i,j)
                - equilibrium: True if equilibrium reaction
                - equilibrium_constraint: Constraint if equilibrium
        """
        if 'reactions' not in self.config:
            return []
        
        reactions = self._get_ordered_reactions()
        global_components = self._build_global_rate_components()
        parameter_names = self._get_parameter_names()
        result = []
        
        for reaction in reactions:
            rxn_data = {
                'id': reaction.get('id'),
                'name': reaction.get('name'),
                'rate_yaml': reaction.get('rate', '')
            }
            
            # Get rate and components
            rate = reaction.get('rate', '')
            components = dict(global_components)
            components.update(
                {
                    name: str(expr)
                    for name, expr in reaction.get('rate_components', {}).items()
                    if name not in parameter_names
                }
            )
            
            # Expand rate components
            rate_expanded = self.expand_rate_components(rate, components)
            rxn_data['rate_expanded'] = rate_expanded
            
            # Substitute species names
            rate_fortran = self.substitute_species_in_expression(rate_expanded)
            rxn_data['rate_fortran'] = rate_fortran
            
            # Handle equilibrium reactions
            if reaction.get('equilibrium', False):
                rxn_data['equilibrium'] = True
                constraint = reaction.get('equilibrium_constraint', '')
                constraint_fortran = self.substitute_species_in_expression(constraint)
                rxn_data['equilibrium_constraint'] = constraint_fortran
            else:
                rxn_data['equilibrium'] = False
            
            result.append(rxn_data)
        
        return result
    
    def get_transport_parameters(self) -> Dict[str, Dict[str, float]]:
        """
        Extract transport parameters for each species.
        
        Returns:
            Dictionary mapping species name to transport params:
                - D0: Molecular diffusion coefficient
                - alpha: Temperature dependence
                - tortuosity: Tortuosity factor
        """
        transport_params = {}
        
        for species in self.species_list:
            name = species['name']
            transport_params[name] = {
                'D0': species.get('transport_D0', 0.0),
                'alpha': species.get('transport_alpha', 0.0),
                'tortuosity': species.get('transport_tortuosity', 0.0)
            }
        
        return transport_params
    
    def get_boundary_conditions(self, evaluator=None) -> Dict[str, Dict[str, Any]]:
        """
        Extract boundary conditions for each species.
        
        Args:
            evaluator: Optional FormulaEvaluator for evaluating BC values
        
        Returns:
            Dictionary mapping species name to boundary conditions:
                - bc_upper_type: Upper BC type (0=Dirichlet, 1=Neumann, 2=flux)
                - bc_upper_value: Upper BC value
                - bc_lower_type: Lower BC type
                - bc_lower_value: Lower BC value
        """
        bc_data = {}
        
        for species in self.species_list:
            name = species['name']
            
            upper_val = species.get('bc_upper_value', 0.0)
            lower_val = species.get('bc_lower_value', 0.0)
            
            # Evaluate if needed
            if evaluator and isinstance(upper_val, str):
                species_vals = evaluator.get_species_values([species])
                upper_val = species_vals[name].get('bc_upper_value', upper_val)
                lower_val = species_vals[name].get('bc_lower_value', lower_val)
            
            bc_data[name] = {
                'bc_upper_type': species.get('bc_upper_type', 0),
                'bc_upper_value': upper_val,
                'bc_lower_type': species.get('bc_lower_type', 1),
                'bc_lower_value': lower_val
            }
        
        return bc_data
    
    def get_initial_conditions(self, evaluator=None) -> Dict[str, Dict[str, Any]]:
        """
        Extract initial conditions for each species.
        
        Args:
            evaluator: Optional FormulaEvaluator for evaluating init values
        
        Returns:
            Dictionary mapping species name to initial conditions:
                - init_value: Initial concentration/amount
                - init_mode: Mode (3=constant, other values for profiles)
        """
        init_data = {}
        global_mode = self.config.get('initial_conditions', {}).get('mode', 3)
        
        for species in self.species_list:
            name = species['name']
            
            init_val = species.get('init_value', 0.0)
            
            # Evaluate if needed
            if evaluator and isinstance(init_val, str):
                species_vals = evaluator.get_species_values([species])
                init_val = species_vals[name].get('init_value', init_val)
            
            init_data[name] = {
                'init_value': init_val,
                # Global mode (Maple vic/ic semantics), with species-level
                # fallback for backward compatibility.
                'init_mode': species.get('init_mode', global_mode)
            }
        
        return init_data

    def get_timestep_parameters(self) -> Dict[str, Any]:
        """
        Extract timestep initialization parameters from advanced YAML block.

        Reads optional `advanced.timestep_parameters` with keys:
        - dtold
        - ttol
        - tstep
        - maxconc

        If a key is missing, default BRNS values are used.

        Returns:
            Dictionary with keys dtold, ttol, tstep, maxconc.
            Values may be numeric or strings (e.g. Fortran literals).
        """
        defaults: Dict[str, Any] = {
            'dtold': '1.d-6',
            'ttol': '5.d-2',
            'tstep': '1.0d0',
            'maxconc': '0.d0',
        }

        advanced_cfg = self.config.get('advanced', {})
        if not isinstance(advanced_cfg, dict):
            return defaults.copy()

        ts_cfg = advanced_cfg.get('timestep_parameters', {})
        if not isinstance(ts_cfg, dict):
            return defaults.copy()

        result = defaults.copy()
        for key in defaults.keys():
            if key in ts_cfg and ts_cfg[key] is not None:
                result[key] = ts_cfg[key]

        return result


def main():
    """Example usage with Canfield model."""
    import yaml
    import os
    import sys
    
    # Add parent directory for imports
    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
    from acg_brns.formula_evaluator import FormulaEvaluator
    
    # Load YAML
    yaml_path = os.path.join(os.path.dirname(__file__), '..', 'models', 'equilibrium', 'equilibrium.yaml')
    with open(yaml_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Evaluate formulas
    print("Evaluating formulas...")
    evaluator = FormulaEvaluator(config)
    evaluated = evaluator.evaluate_all()
    
    # Create mapper
    print("\nMapping to ACG structures...")
    mapper = YAMLtoACGMapper(config, evaluated)
    
    # Build bio arrays
    bio_name, bio_val = mapper.build_bio_arrays()
    print(f"\nBio-parameters: {len(bio_name)}")
    print(f"  Sample: {bio_name[:5]}")
    print(f"  Values: {bio_val[:5]}")
    
    # Build variables
    variables = mapper.build_variables_list()
    print(f"\nSpecies (variables): {len(variables)}")
    print(f"  Dissolved: {variables[:12]}")
    print(f"  Solid: {variables[12:]}")
    
    # Species map
    species_map = mapper.build_species_index_map()
    print(f"\nSpecies index map:")
    for name, idx in list(species_map.items())[:5]:
        print(f"  {name}: sp({idx},j)")
    
    # Stoichiometry matrix
    stoich = mapper.build_stoichiometry_matrix(evaluator)
    print(f"\nStoichiometry matrix: {stoich.shape}")
    print(f"  Non-zero entries: {np.count_nonzero(stoich)}")
    
    # Rate expressions
    reactions = mapper.build_rate_expressions()
    print(f"\nReactions: {len(reactions)}")
    print(f"\nSample reaction (Reaction 1 - Aerobic respiration):")
    rxn = reactions[0]
    print(f"  Name: {rxn['name']}")
    print(f"  YAML rate: {rxn['rate_yaml']}")
    print(f"  Fortran rate: {rxn['rate_fortran'][:100]}...")
    
    print(f"\nEquilibrium reactions:")
    for rxn in reactions:
        if rxn.get('equilibrium'):
            print(f"  - {rxn['name']} (Reaction {rxn['id']})")


if __name__ == '__main__':
    main()
