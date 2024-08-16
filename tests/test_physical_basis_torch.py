import pytest
import torch
from physical_basis.torch import PhysicalBasis

def test_forward():
    physical_basis = PhysicalBasis()
    x = torch.tensor([0.1, 0.2, 0.3], requires_grad=True)
    n = 2
    l = 1
    result = physical_basis(n, l, x)
    assert isinstance(result, torch.Tensor)
    assert result.shape == x.shape

    dummy_loss = torch.sum(result)
    dummy_loss.backward()
    assert x.grad is not None
    
    # Test edge cases
    with pytest.raises(ValueError):
        physical_basis(n, physical_basis.l_max + 1, x)
    with pytest.raises(ValueError):
        physical_basis(physical_basis.n_max_big, l, x)

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
