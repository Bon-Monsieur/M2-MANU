import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# =========================================================
# PHYSICAL PARAMETERS
# =========================================================
K = 0.1        # diffusion coefficient
L = 1.0        # domain length
Time = 20.0
V = 1.0
lamda = 1.0

# =========================================================
# NUMERICAL PARAMETERS
# =========================================================
NX_init = 10
NT      = 10000
eps     = 1e-3
niter_refinement = 10

# =========================================================
# Exact solution and its derivatives
# =========================================================
def u_exact(x):
    return np.exp(-20*(x-0.5)**2)

def u_x(x):
    return -40*(x-0.5)*np.exp(-20*(x-0.5)**2)

def u_xx(x):
    return (-40 + 1600*(x-0.5)**2)*np.exp(-20*(x-0.5)**2)

# H2 norm of exact solution for relative error
x_fine = np.linspace(0,1,4001)
norm_H2 = np.sqrt(
    np.trapz(u_exact(x_fine)**2 + u_x(x_fine)**2 + u_xx(x_fine)**2, x_fine)
)

# =========================================================
# Storage arrays
# =========================================================
dx_tab, L2_T, L2_Tx, L2_interp, L2_Txx, rel_error = [], [], [], [], [], []

# =========================================================
# Figures (3 subplots)
# =========================================================
fig, axes = plt.subplots(3, 1, figsize=(8, 12))
ax_sol, ax_res, ax_err = axes

NX = NX_init
for _ in range(niter_refinement):
    NX += 5
    dx = L / (NX - 1)
    dt = dx**2 / (V*dx + 4*K + dx**2)
    dx_tab.append(dx)

    # --- Exact solution and source term on coarse mesh ---
    x   = np.linspace(0.0, 1.0, NX)
    Tex = u_exact(x)
    F   = np.zeros(NX)
    for j in range(1, NX - 1):
        Tx  = (Tex[j + 1] - Tex[j - 1]) / (2 * dx)
        Txx = (Tex[j + 1] - 2 * Tex[j] + Tex[j - 1]) / (dx ** 2)
        F[j] = V * Tx - K * Txx + lamda * Tex[j]

    # --- Initial numerical solution ---
    T   = np.zeros(NX)
    RHS = np.zeros(NX)
    rest = []
    n, res, res0 = 0, 1.0, 1.0

    # --- Time stepping ---
    while n < NT and res / res0 > eps:
        n += 1
        res = 0.0
        for j in range(1, NX - 1):
            xnu = K + 0.5 * dx * abs(V)
            Tx  = (T[j + 1] - T[j - 1]) / (2 * dx)
            Txx = (T[j - 1] - 2 * T[j] + T[j + 1]) / (dx ** 2)
            RHS[j] = dt * (-V * Tx + xnu * Txx - lamda * T[j] + F[j])
            res += abs(RHS[j])
        T[1:-1] += RHS[1:-1]
        if n == 1:
            res0 = res
        rest.append(res)

    # --- Plot solution and residual ---
    ax_sol.plot(x, T, label=f'dx={dx:.3f}')
    ax_res.plot(np.log10(np.array(rest) / rest[0]), label=f'dx={dx:.3f}')

    # ===== Norms ==================================================
    # 1) L2(T - Tex)
    err_L2_T = np.sqrt(np.sum((T - Tex)**2) * dx)

    # 2) L2(T' - Tex')
    dTdx  = np.zeros_like(T)
    Tex_x = np.zeros_like(T)
    for j in range(1, NX - 1):
        dTdx[j]  = (T[j + 1] - T[j - 1]) / (2 * dx)
        Tex_x[j] = u_x(x[j])
    err_L2_Tx = np.sqrt(np.sum((dTdx - Tex_x) ** 2) * dx)

    # 3) L2(Tex - I_h Tex) with linear interpolation on fine mesh
    x_fine_sub = np.linspace(0, 1, 10 * (NX - 1) + 1)
    Tex_fine   = u_exact(x_fine_sub)
    IhTex      = np.interp(x_fine_sub, x, Tex)
    err_L2_interp = np.sqrt(np.trapz((Tex_fine - IhTex) ** 2, x_fine_sub))

    # 4) L2(Tex'')
    err_L2_Txx = np.sqrt(np.sum(u_xx(x)**2) * dx)

    # Store
    L2_T.append(err_L2_T)
    L2_Tx.append(err_L2_Tx)
    L2_interp.append(err_L2_interp)
    L2_Txx.append(err_L2_Txx)
    rel_error.append(err_L2_T / norm_H2)

# =========================================================
# Table of norms
# =========================================================
df = pd.DataFrame({
    'dx': dx_tab,
    '||T-Tex||_L2': L2_T,
    '||T\'-Tex\'||_L2': L2_Tx,
    '||Tex - IhTex||_L2': L2_interp,
    '||Tex\'\'||_L2': L2_Txx
})
print("\n=== Normes discrètes ===")
print(df.to_string(index=False))

# =========================================================
# Fit log(error) = log(C) + k log(h)
# =========================================================
lnh = np.log(dx_tab)
lne = np.log(rel_error)
p   = np.polyfit(lnh, lne, 1)
k   = p[0]
C   = math.exp(p[1])
print(f"\n--- Estimation (C,k) ---\nC = {C:.4e}, k = {k:.4f}")

# =========================================================
# Plot of discrete norms
# =========================================================
ax_err.loglog(dx_tab, L2_T,     'o-', label=r'$\|T-T_{ex}\|_{L^2}$')
ax_err.loglog(dx_tab, L2_Tx,    's-', label=r'$\|T^\prime-T_{ex}^\prime\|_{L^2}$')
ax_err.loglog(dx_tab, L2_interp,'d-', label=r'$\|T_{ex}-I_h T_{ex}\|_{L^2}$')
ax_err.loglog(dx_tab, L2_Txx,   'x-', label=r'$\|T_{ex}^{\prime\prime}\|_{L^2}$')

# Titles and labels
ax_sol.set_title('Solution T(x) pour différents raffinements')
ax_sol.set_xlabel('x'); ax_sol.set_ylabel('T'); ax_sol.legend()

ax_res.set_title('Évolution du résidu log10(res/res0)')
ax_res.set_xlabel('Itérations'); ax_res.set_ylabel('log10(res/res0)')
ax_res.legend()

ax_err.set_title('Normes discrètes en fonction de dx (log-log)')
ax_err.set_xlabel('dx'); ax_err.set_ylabel('Norme L2')
ax_err.legend()

plt.tight_layout()
plt.show()

# =========================================================
# Superposition de l'erreur relative et des lois C h^k, C h^{k+1}
# =========================================================
h_plot = np.linspace(min(dx_tab), max(dx_tab), 200)
plt.figure(figsize=(7,5))
plt.loglog(dx_tab, rel_error, 'o', label=r'$\|T-T_{ex}\|_{0,2}/\|T_{ex}\|_{2,2}$')
plt.loglog(h_plot, C * h_plot**k,  '--', label=fr'$C\,h^{{{k:.2f}}}$')
plt.loglog(h_plot, C * h_plot**(k+1), ':', label=fr'$C\,h^{{{k:.2f}+1}}$')
plt.gca().invert_xaxis()
plt.xlabel('Pas de maillage h')
plt.ylabel('Erreur relative')
plt.title('Identification de C et k')
plt.grid(True, which='both', ls=':')
plt.legend()
plt.show()
