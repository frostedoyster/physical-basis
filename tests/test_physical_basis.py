import pytest
import numpy as np
from physical_basis import PhysicalBasis

def test_compute():
    physical_basis = PhysicalBasis()
    x = np.array([0.1, 0.2, 0.3])
    n = 2
    l = 1
    result = physical_basis.compute(n, l, x)
    assert isinstance(result, np.ndarray)
    assert result.shape == x.shape
    
    # Test edge cases
    with pytest.raises(ValueError):
        physical_basis.compute(n, physical_basis.l_max + 1, x)
    with pytest.raises(ValueError):
        physical_basis.compute(physical_basis.n_max_big, l, x)

def test_compute_derivative():
    physical_basis = PhysicalBasis()
    x = np.array([0.1, 0.2, 0.3])
    n = 2
    l = 1
    result = physical_basis.compute_derivative(n, l, x)
    assert isinstance(result, np.ndarray)
    assert result.shape == x.shape
    
    # Test edge cases
    with pytest.raises(ValueError):
        physical_basis.compute_derivative(n, physical_basis.l_max + 1, x)
    with pytest.raises(ValueError):
        physical_basis.compute_derivative(physical_basis.n_max_big, l, x)

def test_return_eigenvalue():
    physical_basis = PhysicalBasis()
    n = 2
    l = 1
    result = physical_basis.return_eigenvalue(n, l)
    assert result == physical_basis.E_ln[l, n]
    
    # Test edge cases
    with pytest.raises(ValueError):
        physical_basis.return_eigenvalue(n, physical_basis.l_max + 1)
    with pytest.raises(ValueError):
        physical_basis.return_eigenvalue(physical_basis.n_max_big, l)
