This is a repository showing how to create a physically inspired basis for atomistic machine learning.

- if the required accuracy is not high, the `physical_basis.py` script allows to solve the eigenvalue equation for the physical basis in Python.

- to obtain better accuracy, it is desirable to include more basis functions. This is slow in Python due to the many numerical integration steps that need to be executed. Therefore, we include a C++ program to calculate these integrals much faster. The C++ program depends on GSL and OpenMP, and it can be compiled, for example, with `g++ integrate.cc -o integrate -lgsl -fopenmp`. This program takes less than 2 minutes on a laptop, and it will dump the `O.txt` and `S.txt` matrices, which can be used by the `generate_coefficients.py` script to generate the `eigenvalues.npy` and `eigenvectors.npy` that fully define the solution. Then, the radial basis functions can be loaded as shown in `plot.py`.

- if you want to use the precomputed `eigenvalues.npy` and `eigenvectors.npy`, you can simply run `plot.py` as an example of how to load the physically inspired basis and plot it.
