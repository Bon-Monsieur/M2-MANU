"""
adrs_fixed_mesh.py
Version corrigée pour répondre aux consignes de la séance 2.
- Maillage fixe (ou séries de maillages pour calcul d'erreurs) : on teste 5 maillages partant de 3 points
- Condition CFL et schéma explicite en temps (upwind pour l'advection, centré pour la diffusion)
- Implémentation du terme source f obtenu à partir de u_exact(s) = exp(-10*(s-L/2)**2)
- Condition de Neumann u_s(L)=0 au bord droit, condition Dirichlet u(0)=u_exact(0) au bord gauche
- Critère de convergence vers l'état stationnaire basé sur ||u^{n+1}-u^n||_L2 normalisé
- Calcul des normes L2 et H1 (discrètes) après convergence pour chaque maillage
- Tracés : solution numérique vs exacte, historique de convergence temporelle, erreurs L2 et H1 vs h (log-log)

Usage : lancer le script. Nécessite numpy et matplotlib.
"""

import numpy as np
import math
import matplotlib.pyplot as plt

# Physical parameters
v = 1.0        # advection speed
nu = 0.01      # diffusion (nu = K in votre code)
lam = 1.0      # reaction coefficient
L = 1.0

# Numerical/time parameters
n_meshes = 5                 # 5 maillages partant de 3 points
NX_start = 3
eps_time = 1e-6              # critère de convergence temporelle (relatif)
CFL = 0.4                    # nombre de Courant pour la partie advective (stability)
max_timesteps = 200000

# Exact solution (stationary in time) and computed forcing f(s)
def u_exact(s):
    return np.exp(-10.0 * (s - L/2.0)**2)

# Compute forcing f(s) so that u_exact is an exact steady solution of:
# 0 + v u_s - nu u_ss + lam u = f(s)
# => f(s) = v u_s - nu u_ss + lam u

def compute_f_from_uex(s):
    # compute derivatives analytically for the Gaussian
    # u = exp(-10*(s-L/2)^2)
    z = -10.0 * (s - L/2.0)**2
    u = np.exp(z)
    du = u * (-20.0 * (s - L/2.0))          # u'
    d2u = u * (400.0 * (s - L/2.0)**2 - 20.0)  # u''
    f = v * du - nu * d2u + lam * u
    return f

# Discrete norms
def discrete_L2(u, uh, dx):
    # L2 norm approximated by sqrt(sum((uh- u)^2)*dx)
    return np.sqrt(np.sum((uh - u)**2) * dx)

def discrete_H1(u, uh, dx):
    # H1 semi-norm: sqrt(sum((duh-d u)^2)*dx)
    # Approximate derivatives with central differences in interior and one-sided at boundaries
    N = len(uh)
    # compute continuous derivatives of u (analytical) at nodes
    du_anal = -20.0 * (np.linspace(0, L, N) - L/2.0) * u_exact(np.linspace(0, L, N))
    # discrete derivative of uh
    duh = np.zeros_like(uh)
    for i in range(1, N-1):
        duh[i] = (uh[i+1] - uh[i-1]) / (2.0 * dx)
    # boundaries: one-sided
    duh[0] = (uh[1] - uh[0]) / dx
    duh[-1] = (uh[-1] - uh[-2]) / dx
    return np.sqrt(np.sum((duh - du_anal)**2) * dx)

# Build list of NX (mesh sizes) starting from 3 points and roughly halving h each time
NX_list = []
NX = NX_start
for k in range(n_meshes):
    NX_list.append(NX)
    # double number of intervals (approx halving h) but keep at least 3 nodes
    NX = max(3, 2 * (NX - 1) + 1)

errors_L2 = []
errors_H1 = []
hs = []

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()

for idx, NX in enumerate(NX_list):
    # Mesh
    x = np.linspace(0.0, L, NX)
    dx = x[1] - x[0]
    hs.append(dx)

    # Forcing sampled at nodes (f_h)
    f_h = compute_f_from_uex(x)

    # Initial condition: zero or u_exact? start from zero to test convergence to stationary
    u = np.zeros_like(x)
    u_old = u.copy()

    # Boundary conditions:
    # left: Dirichlet u(0) = u_exact(0) (you can change to other BCs)
    # right: Neumann u_s(L) = 0 implemented by ghost point or one-sided difference

    # choose time step from CFL and diffusion constraint
    # advective CFL: dt <= CFL * dx / |v|
    # diffusion constraint (explicit) dt <= dx^2 / (2*nu)
    if abs(v) > 1e-12:
        dt_adv = CFL * dx / abs(v)
    else:
        dt_adv = 1e6
    if nu > 0:
        dt_diff = 0.5 * dx**2 / nu
    else:
        dt_diff = 1e6
    dt = min(dt_adv, dt_diff)

    # Time march until stationarity
    conv_hist = []
    norm_prev = np.linalg.norm(u)
    for nstep in range(1, max_timesteps + 1):
        # compute RHS and update explicitly (forward Euler)
        rhs = np.zeros_like(u)

        # interior nodes
        for j in range(1, NX-1):
            # upwind for advection
            if v >= 0:
                du_dx = (u[j] - u[j-1]) / dx
            else:
                du_dx = (u[j+1] - u[j]) / dx

            # diffusion central
            d2u_dx2 = (u[j-1] - 2.0 * u[j] + u[j+1]) / (dx**2)

            rhs[j] = -v * du_dx + nu * d2u_dx2 - lam * u[j] + f_h[j]

        # Boundary conditions:
        # left Dirichlet
        u[0] = u_exact(x[0])
        # right Neumann u_s(L)=0 -> implement by ghost point: u_{N} = u_{N-2}
        # Then second derivative at j=N-1 uses u_{N} as u_{N-2} so becomes (u_{N-2}-2u_{N-1}+u_{N-2})/dx^2 = 2*(u_{N-2}-u_{N-1})/dx^2
        # For advection at last interior node j=N-2 handle with upwind as usual

        # update interior points (explicit Euler)
        u_new = u.copy()
        for j in range(1, NX-1):
            u_new[j] = u[j] + dt * rhs[j]

        # enforce left Dirichlet again
        u_new[0] = u_exact(x[0])
        # enforce Neumann at right: set ghost value u_N = u_{N-2} and update last node accordingly
        # but easier: impose zero derivative by one-sided update: u_N-1 = u_N-2
        u_new[-1] = u_new[-2]

        # compute convergence metric: normalized L2 of difference
        diff = u_new - u
        norm_diff = np.sqrt(np.sum(diff**2) * dx)
        norm_u = np.sqrt(np.sum(u_new**2) * dx)
        if norm_u < 1e-16:
            rel = norm_diff
        else:
            rel = norm_diff / norm_u
        conv_hist.append(rel)

        u = u_new

        # check convergence to steady state
        if rel < eps_time:
            print(f"Converged mesh NX={NX} in {nstep} steps, dt={dt:.2e}, rel={rel:.2e}")
            break
    else:
        print(f"Warning: not converged after max steps for NX={NX}, last rel={rel:.2e}")

    # compute errors against exact solution
    u_ex = u_exact(x)
    errL2 = discrete_L2(u_ex, u, dx)
    errH1 = discrete_H1(u_ex, u, dx)
    errors_L2.append(errL2)
    errors_H1.append(errH1)

    # Plot numerical vs exact (first subplot)
    axes[0].plot(x, u, marker='o', label=f'uh NX={NX}')
    axes[0].plot(x, u_ex, '--', label=f'uex NX={NX}')

    # Plot convergence history for this mesh (second subplot)
    axes[1].semilogy(conv_hist, label=f'NX={NX}')

# Finalize subplot 0 and 1
axes[0].set_title('Solution numérique vs exacte pour chaque maillage')
axes[0].set_xlabel('s')
axes[0].set_ylabel('u(s)')
axes[0].legend()

axes[1].set_title('Historique de convergence temporelle (||u^{n+1}-u^n||_2 normalisé)')
axes[1].set_xlabel('itérations temporelles')
axes[1].set_ylabel('||u^{n+1}-u^n||_2 (rel)')
axes[1].legend()

# Plot errors L2 and H1 vs h in log-log (subplots 2 and 3 merged)
hs = np.array(hs)
errors_L2 = np.array(errors_L2)
errors_H1 = np.array(errors_H1)

# log-log plot
fig2, ax = plt.subplots(1, 1, figsize=(6,5))
ax.loglog(hs, errors_L2, 'o--', label='L2 error')
ax.loglog(hs, errors_H1, 's--', label='H1 error')
ax.set_xlabel('h')
ax.set_ylabel('error')
ax.set_title('Convergence en norme L2 et H1')
ax.grid(True, which='both', ls=':')
ax.legend()

# estimate slopes (orders) by linear fit in log-log
coef_L2 = np.polyfit(np.log(hs), np.log(errors_L2), 1)
coef_H1 = np.polyfit(np.log(hs), np.log(errors_H1), 1)
print(f"Estimated order L2 ~ h^{coef_L2[0]:.2f}")
print(f"Estimated order H1 ~ h^{coef_H1[0]:.2f}")

plt.tight_layout()
plt.show()
