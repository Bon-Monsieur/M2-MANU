import numpy as np
import matplotlib.pyplot as plt

# Paramètres
V = 1.0          # vitesse (m/s)
nu = 1e-2        # coefficient de diffusion (m²/s)
lambda_ = 0.0    # taux de réaction (s⁻¹)
L = 1.0          # longueur du domaine (m)
T = 1.0          # temps final (s)
Dt = 0.1        # pas de temps (s)
Dx = 0.1        # pas d'espace (m)
Nx = int(L / Dx) # nombre de points en espace
Nt = int(T / Dt) # nombre de pas de temps

# Solution exacte
def u_exacte(x):
    return np.exp(-10 * x**2)

# Calcul de f(t,x) pour que u_exacte soit solution
def calcul_f(x):
    u = u_exacte(x)
    u_x = -20 * x * u
    u_xx = (-20 + 400 * x**2) * u
    f = V * u_x - nu * u_xx - lambda_ * u
    return f

# Initialisation
x = np.linspace(0, L, Nx + 1)
u = u_exacte(x)
f = calcul_f(x)

# Stockage des solutions pour le tracé
solutions = []
solutions.append(u.copy())

# Schéma d'Euler explicite
for n in range(Nt):
    u_new = u.copy()
    for i in range(1, Nx):
        # Dérivée première (décentré amont)
        u_x = (u[i] - u[i - 1]) / Dx
        # Dérivée seconde (centrée)
        u_xx = (u[i + 1] - 2 * u[i] + u[i - 1]) / (Dx**2)
        # Mise à jour
        u_new[i] = u[i] + Dt * (-V * u_x + nu * u_xx + lambda_ * u[i] + f[i])

    # Conditions aux limites
    u_new[0] = u_exacte(0)  # u(t, x=0) = 1
    u_new[-1] = u_exacte(L) # u_x(t, x=L) = 0 (flux nul)

    u = u_new

    # Sauvegarde toutes les 10 itérations
    if n % 10 == 0 or True:
        solutions.append(u.copy())

# Tracé des solutions
plt.figure(figsize=(10, 6))
for i, sol in enumerate(solutions):
    plt.plot(x, sol, label=f"Itération {i * 10}")

plt.plot(x, u_exacte(x), 'k--', label="Solution exacte")
plt.xlabel("Position x (m)")
plt.ylabel("u(t, x)")
plt.title("Solution de l'équation de transport-diffusion-réaction (Euler explicite)")
plt.legend()
plt.grid()
plt.show()
