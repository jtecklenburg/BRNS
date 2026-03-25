"""Tests for ACG-BRNS core functionality."""

import pytest
from pathlib import Path
from sympy import symbols, Symbol
from acg_brns import ACGModule


@pytest.fixture
def temp_output_dir(tmp_path):
    """Provide a temporary output directory."""
    return tmp_path / "fortran_output"


class TestACGModuleInit:
    """Test ACGModule initialization."""

    def test_init_f77(self, temp_output_dir):
        """Test initialization (F77 only)."""
        acg = ACGModule(str(temp_output_dir))
        assert acg.output_dir == temp_output_dir
        assert temp_output_dir.exists()

    def test_creates_directory(self, temp_output_dir):
        """Test that output directory is created."""
        assert not temp_output_dir.exists()
        acg = ACGModule(str(temp_output_dir))
        assert temp_output_dir.exists()


class TestACG0:
    """Test ACG0 - common_geo.inc generation."""

    def test_acg0_basic(self, temp_output_dir):
        """Test basic common_geo.inc generation."""
        acg = ACGModule(str(temp_output_dir))
        acg.acg0(
            nsolids=0,
            ndissolved=1,
            nreactions=1,
            nnodes=51,
            bio_names=['k_deg'],
            bio_vals=[0.01]
        )
        
        output_file = temp_output_dir / 'common_geo.inc'
        assert output_file.exists()
        
        content = output_file.read_text()
        assert 'nsolid=0' in content
        assert 'ndiss=1' in content
        assert 'ncomp=1' in content
        assert 'nreac=1' in content
        assert 'k_deg' in content

    def test_acg0_with_physical_params(self, temp_output_dir):
        """Test with physical parameters."""
        acg = ACGModule(str(temp_output_dir))
        acg.acg0(
            nsolids=1,
            ndissolved=2,
            nreactions=3,
            nnodes=101,
            bio_names=['k1', 'k2'],
            phys_names=['phi', 'tort'],
            phys_vals=[0.3, 1.5],
            phys_names2=['iopt'],
            phys_vals2=[1]
        )
        
        output_file = temp_output_dir / 'common_geo.inc'
        content = output_file.read_text()
        assert 'nsolid=1' in content
        assert 'ndiss=2' in content
        assert 'ncomp=3' in content
        assert 'phi' in content
        assert 'iopt' in content

    def test_acg0_allows_empty_bio_names(self, temp_output_dir):
        """Models without biogeochemical parameters should still generate common_geo.inc."""
        acg = ACGModule(str(temp_output_dir))
        acg.acg0(
            nsolids=0,
            ndissolved=1,
            nreactions=1,
            nnodes=51,
            bio_names=[],
            bio_vals=[]
        )

        output_file = temp_output_dir / 'common_geo.inc'
        assert output_file.exists()

        content = output_file.read_text()
        assert 'nsolid=0' in content
        assert 'ndiss=1' in content
        assert 'nreac=1' in content


class TestACG1:
    """Test ACG1 - boundaries.f generation."""

    def test_acg1_basic(self, temp_output_dir):
        """Test basic boundary conditions."""
        acg = ACGModule(str(temp_output_dir))
        acg.acg1(
            ncompo=1,
            type_up=[0],
            bnddata_up=[1.0],
            type_down=[1],
            bnddata_down=[0.0]
        )
        
        output_file = temp_output_dir / 'boundaries.f'
        assert output_file.exists()
        
        content = output_file.read_text()
        assert 'subroutine boundaries' in content
        assert 'spb(1,1)' in content
        assert 'ibc(1,1)' in content


class TestACG2:
    """Test ACG2 - molecular.f generation."""

    def test_acg2_basic(self, temp_output_dir):
        """Test molecular diffusion."""
        acg = ACGModule(str(temp_output_dir))
        acg.acg2(
            ncompo=2,
            diffdata=[1e-9, 2e-9],
            alphadata=[1.0, 1.0]
        )
        
        output_file = temp_output_dir / 'molecular.f'
        assert output_file.exists()
        
        content = output_file.read_text()
        assert 'subroutine molecular' in content
        assert 'dsol_0(1)' in content
        assert 'f_T(1)' in content


class TestACG3:
    """Test ACG3 - biogeo.f generation."""

    def test_acg3_basic(self, temp_output_dir):
        """Test biogeochemical parameters."""
        acg = ACGModule(str(temp_output_dir))
        acg.acg3(
            bio_names=['k_deg', 'k_ads'],
            bio_vals=[0.01, 0.5]
        )
        
        output_file = temp_output_dir / 'biogeo.f'
        assert output_file.exists()
        
        content = output_file.read_text()
        assert 'subroutine biogeo' in content
        assert 'k_deg' in content


class TestACG15:
    """Test ACG15 - rates.f generation."""

    def test_acg15_with_sympy(self, temp_output_dir):
        """Test rates with SymPy expressions."""
        acg = ACGModule(str(temp_output_dir))
        
        # Create symbolic expression
        diss_a, k_deg = symbols('diss_a k_deg', real=True, positive=True)
        rate = k_deg * Symbol('sp(1,j)')
        
        acg.acg15(nreactions=1, rate_expressions=[rate])
        
        output_file = temp_output_dir / 'rates.f'
        assert output_file.exists()
        
        content = output_file.read_text()
        assert 'subroutine rates' in content
        assert 'r(1,j)' in content


class TestACG6:
    """Test ACG6 - issolid.f generation."""

    def test_acg6_no_solids(self, temp_output_dir):
        """Test with no solid species."""
        acg = ACGModule(str(temp_output_dir))
        acg.acg16(nsolids=0)
        
        output_file = temp_output_dir / 'issolid.f'
        assert output_file.exists()

    def test_acg6_with_solids(self, temp_output_dir):
        """Test with solid species."""
        acg = ACGModule(str(temp_output_dir))
        acg.acg16(nsolids=2, listsolids=[1, 3])
        
        output_file = temp_output_dir / 'issolid.f'
        content = output_file.read_text()
        assert 'k.eq.1' in content
        assert 'k.eq.3' in content


class TestACG14:
    """Test ACG14 - notransport.f generation."""

    def test_acg14_empty(self, temp_output_dir):
        """Test with no non-transported species."""
        acg = ACGModule(str(temp_output_dir))
        acg.acg14()
        
        output_file = temp_output_dir / 'notransport.f'
        assert output_file.exists()

    def test_acg14_with_species(self, temp_output_dir):
        """Test with non-transported species."""
        acg = ACGModule(str(temp_output_dir))
        acg.acg14(listnotransp=[2, 4])
        
        output_file = temp_output_dir / 'notransport.f'
        content = output_file.read_text()
        assert 'k.eq.2' in content
        assert 'k.eq.4' in content


class TestPythonToFortranDouble:
    """Test _python_to_fortran_double conversion."""

    def test_zero(self, temp_output_dir):
        """Test zero conversion."""
        acg = ACGModule(str(temp_output_dir))
        assert acg._python_to_fortran_double(0) == '0.D0'
        assert acg._python_to_fortran_double(0.0) == '0.D0'

    def test_one(self, temp_output_dir):
        """Test one conversion."""
        acg = ACGModule(str(temp_output_dir))
        assert acg._python_to_fortran_double(1) == '0.1D1'
        assert acg._python_to_fortran_double(1.0) == '0.1D1'

    def test_small_number(self, temp_output_dir):
        """Test small number conversion."""
        acg = ACGModule(str(temp_output_dir))
        result = acg._python_to_fortran_double(0.01)
        assert 'D' in result
        assert result.startswith('0.')


class TestIntegration:
    """Integration tests with complete workflow."""

    def test_single_species_workflow(self, temp_output_dir):
        """Test complete single-species workflow."""
        from acg_brns import create_substitution_dict
        
        # Initialize
        acg = ACGModule(str(temp_output_dir))
        
        # Symbolic setup
        diss_a = symbols('diss_a', real=True, positive=True)
        k_deg = symbols('k_deg', real=True, positive=True)
        
        # Substitution
        subst_dict = create_substitution_dict([diss_a], 1)
        
        # Rate expression
        rate_orig = -k_deg * diss_a
        rate_fortran = rate_orig.subs(subst_dict)
        
        # Generate files
        acg.acg0(
            nsolids=0, ndissolved=1, nreactions=1, nnodes=51,
            bio_names=['k_deg'], bio_vals=[0.01]
        )
        
        acg.acg1(
            ncompo=1,
            type_up=[0], bnddata_up=[1.0],
            type_down=[1], bnddata_down=[0.0]
        )
        
        acg.acg2(ncompo=1, diffdata=[1e-9], alphadata=[1.0])
        
        acg.acg3(bio_names=['k_deg'], bio_vals=[0.01])
        
        acg.acg15(nreactions=1, rate_expressions=[rate_fortran])
        
        # Check all files exist
        assert (temp_output_dir / 'common_geo.inc').exists()
        assert (temp_output_dir / 'boundaries.f').exists()
        assert (temp_output_dir / 'molecular.f').exists()
        assert (temp_output_dir / 'biogeo.f').exists()
        assert (temp_output_dir / 'rates.f').exists()
