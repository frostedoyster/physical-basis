import numpy as np
import os


def s(n, x):
    return np.sin(np.pi*(n+1.0)*x/10.0)

def ds(n, x):
    return np.pi*(n+1.0)*np.cos(np.pi*(n+1.0)*x/10.0)/10.0

def c(n, x):
    return np.cos(np.pi*(n+0.5)*x/10.0)

def dc(n, x):
    return -np.pi*(n+0.5)*np.sin(np.pi*(n+0.5)*x/10.0)/10.0


class PhysicalBasis():
    def __init__(self):
        dir_path = os.path.dirname(os.path.realpath(__file__))

        self.l_max = 50
        self.n_max_big = 200

        self.E_ln = np.load(
            os.path.join(
                dir_path,
                "eigenvalues.npy"
            )
        )
        self.eigenvectors = np.load(        
            os.path.join(
                dir_path,
                "eigenvectors.npy"
            )
        )

    def compute(self, n, l, x):
        if l > self.l_max:
            raise ValueError(f"l must be less than or equal to {self.l_max}")
        if n >= self.n_max_big:
            raise ValueError(f"n must be less than or equal to {self.n_max_big}")
        ret = np.zeros_like(x)
        for m in range(self.n_max_big):
            ret += (self.eigenvectors[l][m, n]*c(m, x) if l == 0 else self.eigenvectors[l][m, n]*s(m, x))
        return ret
    
    def compute_derivative(self, n, l, x):
        if l > self.l_max:
            raise ValueError(f"l must be less than or equal to {self.l_max}")
        if n >= self.n_max_big:
            raise ValueError(f"n must be less than or equal to {self.n_max_big}")
        ret = np.zeros_like(x)
        for m in range(self.n_max_big):
            ret += (self.eigenvectors[l][m, n]*dc(m, x) if l == 0 else self.eigenvectors[l][m, n]*ds(m, x))
        return ret
    
    def return_eigenvalue(self, n, l):
        if l > self.l_max:
            raise ValueError(f"l must be less than or equal to {self.l_max}")
        if n >= self.n_max_big:
            raise ValueError(f"n must be less than or equal to {self.n_max_big}")
        return self.E_ln[l, n]
