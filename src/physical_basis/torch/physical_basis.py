import numpy as np
import torch
import os


def s(n, x):
    return torch.sin(torch.pi*(n+1.0)*x/10.0)

def c(n, x):
    return torch.cos(torch.pi*(n+0.5)*x/10.0)


class PhysicalBasis(torch.nn.Module):
    def __init__(self):
        super(PhysicalBasis, self).__init__()

        dir_path = os.path.dirname(os.path.realpath(__file__))
        dir_path = os.path.join(dir_path, "..")

        self.l_max = 50
        self.n_max_big = 200

        E_ln = np.load(
            os.path.join(
                dir_path,
                "eigenvalues.npy"
            )
        )
        self.register_buffer("E_ln", torch.tensor(E_ln))

        eigenvectors = np.load(        
            os.path.join(
                dir_path,
                "eigenvectors.npy"
            )
        )
        self.register_buffer("eigenvectors", torch.tensor(eigenvectors))

    def forward(self, n, l, x):
        if l > self.l_max:
            raise ValueError(f"l must be less than or equal to {self.l_max}")
        if n >= self.n_max_big:
            raise ValueError(f"n must be less than or equal to {self.n_max_big}")
        ret = torch.zeros_like(x)
        for m in range(self.n_max_big):
            ret += (self.eigenvectors[l][m, n]*c(m, x) if l == 0 else self.eigenvectors[l][m, n]*s(m, x))
        return ret
    
    def return_eigenvalue(self, n, l):
        if l > self.l_max:
            raise ValueError(f"l must be less than or equal to {self.l_max}")
        if n >= self.n_max_big:
            raise ValueError(f"n must be less than or equal to {self.n_max_big}")
        return self.E_ln[l, n].item()
