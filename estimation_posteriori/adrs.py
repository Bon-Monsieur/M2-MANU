import math
import numpy as np
import matplotlib.pyplot as plt

# PHYSICAL PARAMETERS
K = 0.1      # Diffusion coefficient
L = 1.0      # Domain size
Time = 20.   # Integration time
V = 1
lamda = 1
# NUMERICAL PARAMETERS
NX = 2       # Initial number of grid points
NT = 10000   # Max number of time steps
ifre = 1000000  # Plot every ifre time iterations
eps = 0.001    # Relative convergence ratio
niter_refinement = 10  # Number of mesh refinements
errors = []  # Liste pour stocker les erreurs L2 de chaque raffinement

# Create a figure with subplots
fig, axs = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle("Advection-Diffusion-Reaction-Source 1D", fontsize=16)

for iter in range(niter_refinement):
    NX += 3
    dx = L / (NX - 1)  # Grid step (space)
    dt = dx**2 / (V * dx + 2 * K + dx**2)  # Initial time step (CFL condition)

    # Initialization
    x = np.linspace(0.0, 1.0, NX)
    T = np.zeros((NX))
    F = np.zeros((NX))
    Tex = np.zeros((NX))
    Texx = np.zeros((NX))

    # Exact solution and source term
    for j in range(1, NX - 1):
        Tex[j] = np.sin(2 * j * math.pi / NX)
    for j in range(1, NX - 1):
        Texx[j] = (Tex[j + 1] - Tex[j - 1]) / (2 * dx)
        Txx = (Tex[j + 1] - 2 * Tex[j] + Tex[j - 1]) / (dx**2)
        F[j] = V * Texx[j] - K * Txx + lamda * Tex[j]

    # Adjust time step for stability
    dt = dx**2 / (V * dx + 2 * K + abs(np.max(F)) * dx**2)

    # Main time loop
    n = 0
    res = 1
    res0 = 1
    rest = []

    while n < NT and res / res0 > eps:
        n += 1
        res = 0
        RHS = np.zeros((NX))
        for j in range(1, NX - 1):
            Tx = (T[j + 1] - T[j - 1]) / (2 * dx)
            Txx = (T[j - 1] - 2 * T[j] + T[j + 1]) / (dx**2)
            RHS[j] = dt * (-V * Tx + K * Txx - lamda * T[j] + F[j])
            res += abs(RHS[j])
        for j in range(1, NX - 1):
            T[j] += RHS[j]
        if n == 1:
            res0 = res
        rest.append(res)

    # Plot solution for each refinement
    axs[0, 0].plot(x, T, label=f"NX={NX}, t={n * dt:.2f}")

    # Plot residual convergence for each refinement
    axs[1, 0].plot(np.log10(rest / rest[0]), label=f"NX={NX}")

    # Calculate L2 error for this refinement
    err = np.dot(T - Tex, T - Tex)
    error_l2 = np.sqrt(err)
    errors.append(error_l2)  # Stocker l'erreur L2 pour ce raffinement

    # Plot the L2 error for this refinement on the error plot
    axs[1, 1].scatter(iter, error_l2, label=f"NX={NX}", s=50)

# Finalize plots
axs[0, 0].set_title("Solution numérique pour chaque raffinement")
axs[0, 0].set_xlabel("x")
axs[0, 0].set_ylabel("T(x)")
axs[0, 0].legend()

axs[0, 1].plot(x, Tex, label="Solution exacte", color='red')
axs[0, 1].set_title("Solution exacte")
axs[0, 1].set_xlabel("x")
axs[0, 1].set_ylabel("Tex(x)")
axs[0, 1].legend()

axs[1, 0].set_title("Convergence du résidu (log10) pour chaque raffinement")
axs[1, 0].set_xlabel("Itération")
axs[1, 0].set_ylabel("log10(res/res0)")
axs[1, 0].legend()

# Plot all L2 errors on the same graph
axs[1, 1].plot(range(niter_refinement), errors, marker='o', linestyle='--', color='blue', label="Erreur L2")
axs[1, 1].set_title("Erreur L2 pour chaque raffinement")
axs[1, 1].set_xlabel("Itération de raffinement")
axs[1, 1].set_ylabel("Erreur L2")
axs[1, 1].set_xticks(range(niter_refinement))
axs[1, 1].set_xticklabels([str(2 + 3 * i) for i in range(niter_refinement)])
axs[1, 1].grid(True)
axs[1, 1].legend()

plt.tight_layout()
plt.show()
