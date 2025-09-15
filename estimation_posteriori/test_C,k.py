import numpy as np
import matplotlib.pyplot as plt
from math import log
from numpy import trapz

# -------- Exact solution and its derivatives -------------
def u_exact(x):
    return np.exp(-20*(x-0.5)**2)

def u_x(x):
    return -40*(x-0.5)*np.exp(-20*(x-0.5)**2)

def u_xx(x):
    return (-40 + 1600*(x-0.5)**2)*np.exp(-20*(x-0.5)**2)

# ---- compute ||u||_{H2} = sqrt( ∫ (u^2 + u_x^2 + u_xx^2) )
x_fine = np.linspace(0,1,2001)
norm_H2 = np.sqrt(trapz(u_exact(x_fine)**2 +
                        u_x(x_fine)**2 +
                        u_xx(x_fine)**2, x_fine))

# -------- Solve ADRS and compute L2 errors --------------
def solve_and_error(NX):
    K, V, lam = 0.1, 1.0, 1.0
    dx = 1.0/(NX-1)
    dt = dx**2/(V*dx + 4*K + dx**2)
    x  = np.linspace(0,1,NX)
    Tex = u_exact(x)

    # source term
    F = np.zeros(NX)
    for j in range(1,NX-1):
        Tx  = (Tex[j+1]-Tex[j-1])/(2*dx)
        Txx = (Tex[j+1]-2*Tex[j]+Tex[j-1])/(dx**2)
        F[j] = V*Tx - K*Txx + lam*Tex[j]

    # time stepping explicit
    T   = np.zeros(NX)
    RHS = np.zeros(NX)
    NT, eps = 10000, 1e-3
    n, res0 = 0, 1.0
    while n < NT:
        n += 1
        res = 0.0
        for j in range(1,NX-1):
            xnu = K + 0.5*dx*abs(V)
            Tx  = (T[j+1]-T[j-1])/(2*dx)
            Txx = (T[j-1]-2*T[j]+T[j+1])/(dx**2)
            RHS[j] = dt * (-V*Tx + xnu*Txx - lam*T[j] + F[j])
            res += abs(RHS[j])
        T[1:-1] += RHS[1:-1]
        if n==1: res0=res
        if res/res0 < eps: break

    # L2 error of solution
    err_L2 = np.sqrt(np.sum((T-Tex)**2)*dx)
    return dx, err_L2

# --- refinement loop
dxs, rel_errors = [], []
NX = 10
for _ in range(10):
    NX += 5
    dx, e = solve_and_error(NX)
    dxs.append(dx)
    rel_errors.append(e / norm_H2)   # relative error

dxs = np.array(dxs)
rel_errors = np.array(rel_errors)

# --------- Fit log(error) = log C + k log h ---------------
lnh = np.log(dxs)
lne = np.log(rel_errors)
p = np.polyfit(lnh, lne, 1)
k  = p[0]             # pente
C  = np.exp(p[1])     # ordonnée
print(f"Identifié :  C = {C:.4e} ,  k = {k:.4f}")

# --------- Courbes à superposer --------------------------
h_plot = np.linspace(dxs.min(), dxs.max(), 200)
fit_k   = C * h_plot**k
fit_k1  = C * h_plot**(k+1)

# --------- Affichage unique ------------------------------
plt.figure(figsize=(7,5))
plt.loglog(dxs, rel_errors, 'o', label=r"$\|u-u_h\|_{0,2}/\|u\|_{2,2}$ (données)")
plt.loglog(h_plot, fit_k,  '--', label=fr"$C\,h^{k:.2f}$")
plt.loglog(h_plot, fit_k1, ':',  label=fr"$C\,h^{{{k:.2f}+1}}$")
plt.gca().invert_xaxis()
plt.xlabel("Pas de maillage h")
plt.ylabel(r"Erreur relative")
plt.title(r"Identification $(C,k)$ et superposition $C h^k$ / $C h^{k+1}$")
plt.grid(True, which='both', ls=':')
plt.legend()
plt.show()
