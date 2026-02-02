import matplotlib.pyplot as plt
import numpy as np

p = 1
rho = 1
u = 0
T = 0.2


def schema_euler(n_cell,CFL):
    x = np.linspace(0, 1, n_cell+1)
    dx = x[1]-x[0]
    sol = np.array([[0.0, 0.0, 0.0]]*n_cell)  # [V1, V2, V3] for each cell
    t = 0.0

    # Initialisation
    for i in range(n_cell):
        if x[i] < 0.5:
            sol[i] = [1,1,0]
        else:
            sol[i] = [0.125,0.1,0]

    while (t<T):
        dt = CFL * dx/2/np.max()
    return