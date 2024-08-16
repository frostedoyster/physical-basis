import numpy as np
import matplotlib.pyplot as plt
from physical_basis import PhysicalBasis

physical_basis = PhysicalBasis()

for l in range(9):
    r = np.linspace(0, 10, 1000)
    for n in range(5):
        y = physical_basis.compute(n, l, r)
        if y[20] < 0:
            y = -y
        plt.plot(r, y, label=f"n={n}")
    plt.legend()
    plt.xlabel("r (A)")
    plt.title(f"l={l}")
    plt.savefig(f"physical_{l}_new.pdf")
    plt.clf()
