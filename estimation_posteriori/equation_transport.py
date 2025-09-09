import numpy as np
import matplotlib.pyplot as plt

# Paramètres
V = 1.0        # vitesse (m/s)
L = 1.0        # longueur du domaine (m)
T = 1.0        # temps final (s)
Dt = 0.01      # pas de temps (s)
Dx = 0.01      # pas d'espace (m)
Nx = int(L / Dx)  # nombre de points en espace
Nt = int(T / Dt)  # nombre de pas de temps

# Initialisation
x = np.linspace(0, L, Nx + 1)
u = np.zeros(Nx + 1)  # u(t=0, x) = 0

# Stockage des solutions pour le tracé
solutions = []
solutions.append(u.copy())

# Schéma d'Euler explicite (amont)
for n in range(1, Nt + 1):
    u_new = u.copy()
    for i in range(1, Nx):
        # Schéma décentré amont : u_t + V u_x = 0
        u_new[i] = u[i] - (V * Dt / Dx) * (u[i] - u[i - 1])

    # Conditions aux limites
    u_new[0] = 1.0  # u(t, x=0) = 1
    u_new[-1] = u_new[-2]  # u_x(t, x=L) = 0 (flux nul)

    u = u_new

    # Sauvegarde toutes les 10 itérations
    if n % 10 == 0:
        solutions.append(u.copy())

# Tracé des solutions
plt.figure(figsize=(10, 6))
for i, sol in enumerate(solutions):
    plt.plot(x, sol, label=f"Itération {i * 10}")

plt.xlabel("Position x (m)")
plt.ylabel("u(t, x)")
plt.title("Solution de l'équation de transport (Euler explicite)")
plt.legend()
plt.grid()
plt.show()
