import numpy as np
import scipy as sp


l_max = 50
n_max_big = 200

O = np.loadtxt("O.txt").reshape((l_max+1, n_max_big, n_max_big))
S = np.loadtxt("S.txt").reshape((l_max+1, n_max_big, n_max_big))

eigenvalues = np.zeros((l_max+1, n_max_big))
eigenvectors = np.zeros((l_max+1, n_max_big, n_max_big))

for l in range(0, l_max+1):
    eigval, eigvec = sp.linalg.eigh(O[l], S[l])
    eigenvalues[l] = eigval
    eigenvectors[l] = eigvec

np.save("eigenvalues.npy", eigenvalues)
np.save("eigenvectors.npy", eigenvectors)
