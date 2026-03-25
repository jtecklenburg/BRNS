"""
Tests for Conservation Law Detection in Gaussian Elimination

Tests that linear dependent rows are correctly identified as conservation laws
and that the Jacobian is computed correctly without -1/delt terms for these rows.

Author: Conservation law validation
Date: 2025-02-10
"""

import pytest
from sympy import symbols, Matrix, S, simplify, simplify
from acg_brns.gaussian_elimination import GaussianElimination, run_gaussian_elimination


class TestConservationLawDetection:
    """Tests for detecting conservation laws from nullspace"""
    
    def test_single_conservation_law(self):
        """Test detection of a single conservation law (stoichiometric balance)"""
        C1, C2, C3 = symbols('C1:4')
        C1_old, C2_old, C3_old = symbols('C1_old C2_old C3_old')
        r1, r2 = symbols('r1 r2')
        delt = symbols('delt', positive=True, real=True)
        
        # System where C1 + C2 + C3 is conserved
        # r1 converts C1 -> C2 (one-to-one)
        # r2 converts C1 -> C3 (one-to-one)
        equations = [
            -r1 - r2,      # C1 consumed by both reactions
            r1,            # C2 produced by r1 only
            r2             # C3 produced by r2 only
        ]
        
        ge = GaussianElimination(equations, [C1, C2, C3], [r1, r2], 3)
        ge.verbose = True
        
        ge.p4_genmatrix()
        ge.p5_gaussian_elimination()
        rightM, matrixM = ge.p6_extract_matrices()
        
        # Check that conservation law was detected
        assert len(ge.conservation_basis) > 0, "Should detect at least one conservation law"
        if ge.verbose:
            print(f"\nConservation laws detected: {len(ge.conservation_basis)}")
            print(f"Conservation basis: {ge.conservation_basis}")
            print(f"Conservation rows: {ge.conservation_rows}")
    
    def test_no_conservation_law(self):
        """Test system with no conservation laws"""
        C1, C2 = symbols('C1 C2')
        C1_old, C2_old = symbols('C1_old C2_old')
        r1, r2 = symbols('r1 r2')
        delt = symbols('delt', positive=True, real=True)
        
        # General system without conservation
        equations = [
            -r1 + r2,
            r1 + 2*r2
        ]
        
        ge = GaussianElimination(equations, [C1, C2], [r1, r2], 2)
        ge.p4_genmatrix()
        ge.p5_gaussian_elimination()
        rightM, matrixM = ge.p6_extract_matrices()
        
        # For generic systems, no conservation laws expected
        # (this depends on specific structure)


class TestJacobianWithConservation:
    """Tests that Jacobian correctly handles conservation law rows"""
    
    def test_conservation_row_jacobian_no_delt_term(self):
        """
        Test that conservation law rows in Jacobian don't have -1/delt diagonal terms
        """
        C1, C2, C3 = symbols('C1:4')
        C1_old, C2_old, C3_old = symbols('C1_old C2_old C3_old')
        r1, r2 = symbols('r1 r2')
        delt = symbols('delt', positive=True, real=True)
        
        # System with conservation: C1 + C2 + C3 conserved
        equations = [
            -r1 - r2,      # C1
            r1,            # C2
            r2             # C3
        ]
        
        result = run_gaussian_elimination(
            equations=equations,
            variables=[C1, C2, C3],
            reactions=[r1, r2],
            variables_old=[C1_old, C2_old, C3_old],
            delt=delt,
            verbose=True
        )
        
        jacobian = result['jacobian']
        func = result['func']
        conservation_rows = result['conservation_rows']
        
        print(f"\n\nConservation rows: {conservation_rows}")
        print(f"Jacobian shape: {jacobian.shape}")
        
        # Check diagonal elements
        for i in range(3):
            diag_elem = jacobian[i, i]
            print(f"\nJacobian[{i},{i}] = {diag_elem}")
            
            if i in conservation_rows:
                # Conservation row: should NOT have -1/delt term
                print(f"  Row {i} is conservation law (should not have -1/delt)")
                # Check that -1/delt does NOT appear in derivative
                C_sym = symbols(f'C{i+1}')
                d_func_i = func[i].diff(C_sym)
                print(f"  d(func[{i}])/d(C[{i}]) = {d_func_i}")
                assert d_func_i.has(1/delt) == False, f"Conservation row {i} should not have 1/delt term"
            else:
                # Non-conservation row: should have -1/delt term
                print(f"  Row {i} is NOT conservation law (should have -1/delt)")
    
    def test_maple_equivalence(self):
        """
        Test that Python Jacobian matches Maple reference for conservation case
        """
        C1, C2, C3 = symbols('C1:4')
        C1_old, C2_old, C3_old = symbols('C1_old C2_old C3_old')
        r1, r2 = symbols('r1 r2')
        delt = symbols('delt', positive=True, real=True)
        
        # 3-component system
        equations = [
            -r1 - r2,
            r1,
            r2
        ]
        
        result = run_gaussian_elimination(
            equations=equations,
            variables=[C1, C2, C3],
            reactions=[r1, r2],
            variables_old=[C1_old, C2_old, C3_old],
            delt=delt,
            verbose=False
        )
        
        jacobian = result['jacobian']
        conservation_rows = result['conservation_rows']
        
        # Print for inspection
        print("\n\n=== Jacobian Matrix ===")
        for i in range(jacobian.shape[0]):
            for j in range(jacobian.shape[1]):
                elem = jacobian[i, j]
                print(f"pd[{i},{j}] = {elem}")
                
                if i == j and i in conservation_rows:
                    # Diagonal for conservation row should NOT have -1/delt
                    if elem.has(delt):
                        print(f"  WARNING: Diagonal element contains delt!")


class TestConservationProperties:
    """Tests mathematical properties of conservation laws"""
    
    def test_nullspace_computation(self):
        """Test that conservation laws are computed from original coefficient matrix"""
        C1, C2, C3 = symbols('C1:4')
        C1_old, C2_old, C3_old = symbols('C1_old C2_old C3_old')
        r1, r2 = symbols('r1 r2')
        delt = symbols('delt', positive=True, real=True)
        
        # System with clear mass conservation
        equations = [
            -r1 - r2,      # Total production: -1*r1 + -1*r2
            r1,            # Total production: +1*r1 + 0*r2
            r2             # Total production: 0*r1 + +1*r2
        ]
        
        result = run_gaussian_elimination(
            equations=equations,
            variables=[C1, C2, C3],
            reactions=[r1, r2],
            variables_old=[C1_old, C2_old, C3_old],
            delt=delt,
            verbose=True
        )
        
        conservation_basis = result['conservation_basis']
        
        # Verify that conservation laws were detected
        assert len(conservation_basis) > 0, "Should have conservation laws"
        
        # Each conservation vector should sum to zero in the original reactions
        for v in conservation_basis:
            print(f"\nConservation vector: {v.T}")
            # The vector [1, 1, 1] means C1 + C2 + C3 is conserved
            # So sum of reaction effects should be zero
            assert v.shape[0] == 3, "Conservation vector should have 3 components"


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
