import numpy as np
import scipy as sp


def s(n, x):
    return np.sin(np.pi*(n+1.0)*x/10.0)

def ds(n, x):
    return np.pi*(n+1.0)*np.cos(np.pi*(n+1.0)*x/10.0)/10.0

def d2s(n, x):
    return -np.pi**2*(n+1.0)**2*np.sin(np.pi*(n+1.0)*x/10.0)/100.0

def c(n, x):
    return np.cos(np.pi*(n+0.5)*x/10.0)

def dc(n, x):
    return -np.pi*(n+0.5)*np.sin(np.pi*(n+0.5)*x/10.0)/10.0

def d2c(n, x):
    return -np.pi**2*(n+0.5)**2*np.cos(np.pi*(n+0.5)*x/10.0)/100.0


l_max = 8
n_max_big = 30
a = 10.0


O = np.zeros((l_max+1, n_max_big, n_max_big))
S = np.zeros((l_max+1, n_max_big, n_max_big))


def function_to_integrate_O_sin(x, l, n1, n2):
    return s(n1, x)*(- x**2 * d2s(n2, x) - (2.0 * x + x**2) * ds(n2, x) + l*(l+1) * s(n2, x)) * np.exp(x)

def function_to_integrate_S_sin(x, n1, n2):
    return s(n1, x) * s(n2, x) * x**2

def function_to_integrate_O_cos(x, l, n1, n2):
    return c(n1, x)*(- x**2 * d2c(n2, x) - (2.0 * x + x**2) * dc(n2, x) + l*(l+1) * c(n2, x)) * np.exp(x)

def function_to_integrate_S_cos(x, n1, n2):
    return c(n1, x) * c(n2, x) * x**2


for l in range(l_max+1):
    print(l)
    for n1 in range(n_max_big):
        for n2 in range(n_max_big):
            fn_O = function_to_integrate_O_cos if l == 0 else function_to_integrate_O_sin
            fn_S = function_to_integrate_S_cos if l == 0 else function_to_integrate_S_sin
            O[l, n1, n2] = sp.integrate.quad(fn_O, 0, 10.0, args=(l, n1, n2), limit=100)[0]
            S[l, n1, n2] = sp.integrate.quad(fn_S, 0, 10.0, args=(n1, n2), limit=100)[0]

E_ln = np.zeros((l_max+1, n_max_big))
eigenvectors = np.zeros((l_max+1, n_max_big, n_max_big))

for l in range(l_max+1):
    eigval, eigvec = sp.linalg.eigh(O[l], S[l])
    E_ln[l] = eigval
    eigenvectors[l] = eigvec

def physical_basis_function(n, l, x):
    ret = np.zeros_like(x)
    for m in range(n_max_big):
        ret += eigenvectors[l][m, n]*(c(m, x) if l == 0 else s(m, x))
    return ret


import matplotlib.pyplot as plt
for l in range(l_max+1):
    x = np.linspace(0, 10, 1000)
    for n in range(5):
        y = physical_basis_function(n, l, x)
        if y[10] < 0:
            y = -y
        plt.plot(x, y, label=f"n={n}")
    plt.legend()
    plt.xlabel("r (A)")
    plt.title(f"l={l}")
    plt.savefig(f"physical_{l}_new.pdf")
    plt.clf()
