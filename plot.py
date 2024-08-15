import numpy as np
import os


# All these periodic functions are zeroed for the (unlikely) case where r > 10*r_0
# which is outside the domain where the eigenvalue equation was solved

def s(n, x):
    return np.sin(np.pi*(n+1.0)*x/10.0)

def ds(n, x):
    return np.pi*(n+1.0)*np.cos(np.pi*(n+1.0)*x/10.0)/10.0

def c(n, x):
    return np.cos(np.pi*(n+0.5)*x/10.0)

def dc(n, x):
    return -np.pi*(n+0.5)*np.sin(np.pi*(n+0.5)*x/10.0)/10.0


l_max = 50
n_max_big = 200

dir_path = os.path.dirname(os.path.realpath(__file__))

E_ln = np.load(
    os.path.join(
        dir_path,
        "eigenvalues.npy"
    )
)
eigenvectors = np.load(        
    os.path.join(
        dir_path,
        "eigenvectors.npy"
    )
)

def function_for_splining(n, l, x):
    ret = np.zeros_like(x)
    for m in range(n_max_big):
        ret += (eigenvectors[l][m, n]*c(m, x) if l == 0 else eigenvectors[l][m, n]*s(m, x))
    return ret


import matplotlib.pyplot as plt

for l in range(9):
    r = np.linspace(0, 10, 1000)
    for n in range(5):
        y = function_for_splining(n, l, r)
        if y[20] < 0:
            y = -y
        plt.plot(r, y, label=f"n={n}")
    plt.legend()
    plt.xlabel("r (A)")
    plt.title(f"l={l}")
    plt.savefig(f"physical_{l}_new.pdf")
    plt.clf()
