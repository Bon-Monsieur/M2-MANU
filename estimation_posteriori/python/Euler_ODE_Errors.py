import numpy as np
import matplotlib.pyplot as plt

# Paramètres
a = 1.0
u0 = 1.0
T = 60.0  # 1 minute en secondes

# 1. Solution pour Dt=1s
Dt = 1.0
N = int(T / Dt)
t = np.linspace(0, T, N+1)
u = np.zeros(N+1)
u[0] = u0
for n in range(N):
    u[n+1] = u[n] - a * u[n] * Dt
u_exact = u0 * np.exp(-a * t)

# Tracé 1 : solution exacte, numérique et erreur en temps
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(t, u_exact, label="Solution exacte", linestyle='--')
plt.plot(t, u, label="Solution numérique (Dt=1s)", marker='o', markersize=2)
plt.xlabel("Temps t (s)")
plt.ylabel("u(t)")
plt.title("Solution exacte et numérique (Euler explicite)")
plt.legend()
plt.grid()

plt.subplot(1, 2, 2)
plt.plot(t, np.abs(u - u_exact), label="Erreur en temps", color='red')
plt.xlabel("Temps t (s)")
plt.ylabel("Erreur absolue")
plt.title("Erreur en temps (Dt=1s)")
plt.grid()
plt.tight_layout()
plt.show()

# 2. Calcul des erreurs L2 pour 20 pas de temps décroissants
Dts = np.logspace(-3, 0, 20)[::-1]  # de 1s à 0.001s
errors_u = []
errors_du = []

for dt in Dts:
    N = int(T / dt)
    t = np.linspace(0, T, N+1)
    u = np.zeros(N+1)
    u[0] = u0
    for n in range(N):
        u[n+1] = u[n] - a * u[n] * dt
    u_exact = np.exp(-a * t)
    # Erreur L2 sur u
    errors_u.append(np.sqrt(np.sum((u - u_exact)**2) / N))
    # Dérivée numérique (Euler) et exacte
    du_num = np.zeros(N+1)
    du_num[1:] = (u[1:] - u[:-1]) / dt
    du_exact = -a * u_exact
    # Erreur L2 sur la dérivée
    errors_du.append(np.sqrt(np.sum((du_num - du_exact)**2) / N))

# Tracé 2 : Erreur L2 en fonction du pas de temps
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.loglog(Dts, errors_u, label="Erreur L2 sur u", marker='o')
plt.xlabel("Pas de temps Δt (s)")
plt.ylabel("Erreur L2")
plt.title("Erreur L2 sur la fonction u")
plt.grid()

plt.subplot(1, 2, 2)
plt.loglog(Dts, errors_du, label="Erreur L2 sur la dérivée", marker='o', color='orange')
plt.xlabel("Pas de temps Δt (s)")
plt.ylabel("Erreur L2")
plt.title("Erreur L2 sur la dérivée du")
plt.grid()
plt.tight_layout()
plt.show()
