# Running numerical comparison between centered and upwind schemes
import numpy as np
import math
import matplotlib.pyplot as plt

# Problem parameters (stationary exact solution u_ex)
L = 1.0
v = 1.0
nu = 0.01
lam = 1.0

# Exact solution (Gaussian) and source
def exact_and_source(x):
    Tex = np.exp(-10*(x - L/2.0)**2)
    Tex_s = -20*(x - L/2.0)*Tex
    Tex_ss = (-20 + 400*(x - L/2.0)**2)*Tex
    F = v*Tex_s - nu*Tex_ss + lam*Tex
    return Tex, F

# Time stepping solver (explicit) for given scheme
def run_solver(NX, scheme='centered', CFL=0.4, eps=1e-6, NT=20000):
    dx = L/(NX-1)
    x = np.linspace(0, L, NX)
    Tex, F = exact_and_source(x)
    T = np.zeros_like(x)  # initial condition
    # Dirichlet BC enforced to exact
    T[0] = Tex[0]; T[-1] = Tex[-1]
    # time step chosen from advective and diffusive constraints
    dt_adv = CFL * dx / (abs(v) + 1e-16)
    dt_diff = 0.5 * dx**2 / nu
    dt = min(dt_adv, dt_diff)
    # time loop until residual small or NT
    res0 = None
    rest = []
    for n in range(1, NT+1):
        RHS = np.zeros_like(T)
        res = 0.0
        for j in range(1, NX-1):
            # advection derivative
            if scheme == 'centered':
                Tx = (T[j+1] - T[j-1])/(2*dx)
            elif scheme == 'upwind':
                if v >= 0:
                    Tx = (T[j] - T[j-1])/dx
                else:
                    Tx = (T[j+1] - T[j])/dx
            else:
                raise ValueError("scheme must be 'centered' or 'upwind'")
            # diffusion centered
            Txx = (T[j-1] - 2*T[j] + T[j+1])/(dx**2)
            RHS[j] = dt * (-v*Tx + nu*Txx - lam*T[j] + F[j])
            res += abs(RHS[j])
        # update
        T[1:-1] += RHS[1:-1]
        # reapply Dirichlet BC
        T[0] = Tex[0]; T[-1] = Tex[-1]
        rest.append(res)
        if n == 1:
            res0 = res if res>0 else 1.0
        if res/res0 < eps:
            break

    # --- Normes d'erreur ---
    # L2 norm
    err_L2 = np.sqrt(np.sum((T - Tex)**2)*dx)

    # <<< AJOUT H1 : dérivée exacte et dérivée numérique
    dTex_dx = -20*(x - L/2.0)*Tex
    dT_dx   = np.zeros_like(x)
    dT_dx[1:-1] = (T[2:] - T[:-2])/(2*dx)           # centrée
    err_H1 = np.sqrt(np.sum((T - Tex)**2 + (dT_dx - dTex_dx)**2)*dx)
    # -----------------------

    return {'T':T, 'Tex':Tex, 'x':x,
            'err_L2':err_L2,
            'err_H1':err_H1,     # <<< AJOUT H1
            'rest':rest, 'n':n, 'dt':dt}

# Run experiment over refinements
niter_refinement = 10
NX0 = 5  # start number of points (>2)
NXs = []
dxs = []
errs_center_L2 = []
errs_upwind_L2 = []
errs_center_H1 = []   # <<< AJOUT H1
errs_upwind_H1 = []   # <<< AJOUT H1
results_center = []
results_upwind = []

for k in range(niter_refinement):
    NX = NX0 + 3*k
    NXs.append(NX)
    dx = L/(NX-1)
    dxs.append(dx)

    res_c = run_solver(NX, scheme='centered', CFL=0.2)
    res_u = run_solver(NX, scheme='upwind', CFL=0.2)

    errs_center_L2.append(res_c['err_L2'])
    errs_upwind_L2.append(res_u['err_L2'])
    errs_center_H1.append(res_c['err_H1'])   # <<< AJOUT H1
    errs_upwind_H1.append(res_u['err_H1'])   # <<< AJOUT H1

    results_center.append(res_c)
    results_upwind.append(res_u)

# --- Plot L2 & H1 convergence on same figure ---
plt.figure(figsize=(8,6))
plt.loglog(dxs, errs_center_L2, 'o-', label='Centered  L2')
plt.loglog(dxs, errs_upwind_L2, 's-', label='Upwind    L2')
plt.loglog(dxs, errs_center_H1, 'o--', label='Centered  H1')   # <<< AJOUT H1
plt.loglog(dxs, errs_upwind_H1, 's--', label='Upwind    H1')   # <<< AJOUT H1
plt.gca().invert_xaxis()
plt.xlabel('dx')
plt.ylabel('Error norm')
plt.title('Convergence: L2 and H1 errors vs dx (log-log)')
plt.grid(True, which='both', ls=':')
plt.legend()

# --- Example solution comparison for finest mesh ---
fin = -1
res_c = results_center[fin]; res_u = results_upwind[fin]
x = res_c['x']
plt.figure(figsize=(8,4))
plt.plot(x, res_c['Tex'], 'k-', label='Exact')
plt.plot(x, res_c['T'], 'C0--', marker='o', label=f'Centered NX={NXs[fin]}')
plt.plot(x, res_u['T'], 'C1--', marker='s', label=f'Upwind NX={NXs[fin]}')
plt.xlabel('x'); plt.ylabel('T')
plt.title('Solution: exact vs centered vs upwind (finest mesh)')
plt.legend(); plt.grid(True)

# --- Table of errors and estimated orders (L2 only for brevity) ---
print("NX   dx        L2_center   H1_center   L2_upwind   H1_upwind")
for i,NX in enumerate(NXs):
    print(f"{NX:2d}  {dxs[i]:.4e}  {errs_center_L2[i]:.4e}  {errs_center_H1[i]:.4e}  "
          f"{errs_upwind_L2[i]:.4e}  {errs_upwind_H1[i]:.4e}")

plt.show()
