import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import minimize

# ==========================================================
# ===============  ÉQUATION ADRS AVEC MAILLAGE ==============
# ==========================================================
def ADRS(NX, xcontrol, Target, L=1.0, adapt=True):
    """
    Résout une équation ADRS (Advection-Diffusion-Réaction-Source) sur [0,L].
    Si adapt=True, utilise un maillage raffiné autour des zones actives.
    """
    # Paramètres physiques
    K = 0.1
    V = 1.0
    lamda = 1.0
    Time = 20.

    # Paramètres numériques
    NT = 1000
    eps = 1e-4

    # -------------------
    # Maillage adaptatif
    # -------------------
    if adapt:
        # Maillage plus dense autour des zones de source (L/(i+1))
        base = np.linspace(0, L, NX)
        weight = np.ones_like(base)
        for ic in range(len(xcontrol)):
            weight += np.exp(-100 * (base - L/(ic+1))**2)
        x = np.interp(np.linspace(0, 1, NX), np.cumsum(weight)/np.sum(weight), base)
    else:
        x = np.linspace(0, L, NX)

    dx = np.diff(x)
    dx_min = np.min(dx)
    dt = 0.5 * dx_min**2 / (V * dx_min + 2 * K + abs(np.max(xcontrol)) * dx_min**2)

    # -------------------
    # Initialisation
    # -------------------
    T = np.zeros(NX)
    F = np.zeros(NX)

    # Calcul du second membre (terme source)
    for j in range(1, NX-1):
        for ic in range(len(xcontrol)):
            F[j] += xcontrol[ic] * np.exp(-100 * (x[j] - L/(ic+1))**2)

    # -------------------
    # Boucle temporelle
    # -------------------
    res, res0, n = 1, 1, 0
    RHS = np.zeros(NX)
    while n < NT and res > eps * res0:
        n += 1
        res = 0
        for j in range(1, NX-1):
            dxm, dxp = x[j] - x[j-1], x[j+1] - x[j]
            Txx = (T[j-1] - 2*T[j] + T[j+1]) / (0.5*(dxm+dxp))**2
            Tx = (T[j+1] - T[j-1]) / (dxm + dxp)
            RHS[j] = dt * (-V*Tx + K*Txx - lamda*T[j] + F[j])
            res += abs(RHS[j])
        if n == 1:
            res0 = res
        T[1:-1] += RHS[1:-1]
    cost = np.trapz((T - Target[:NX])**2, x)  # intégrale approchée
    return cost, T, x


# ==========================================================
# ======== CALCULS DE Aij et Bi AVEC INTERPOLATION =========
# ==========================================================
def compute_A_B(sol_list, T0_list, Target_list, x_list, nq=5):
    """
    Calcule A_ij et B_i avec interpolation sur un maillage commun + quadrature Gauss-Legendre.
    nq : nombre de points de Gauss (doit être suffisant pour exactitude jusqu'à degré 2p)
    """
    n = len(sol_list)
    # Maillage commun (fin)
    xmin = min(x[0] for x in x_list)
    xmax = max(x[-1] for x in x_list)
    N_ref = max(len(x)*4 for x in x_list)
    x_ref = np.linspace(xmin, xmax, N_ref)

    # Interpolation sur maillage commun
    T_ref = []
    for T, x in zip(sol_list, x_list):
        f = interp1d(x, T, kind='cubic', fill_value="extrapolate")
        T_ref.append(f(x_ref))
    fT0 = interp1d(x_list[0], T0_list[0], kind='cubic', fill_value="extrapolate")
    fTarget = interp1d(x_list[0], Target_list[0], kind='cubic', fill_value="extrapolate")
    T0_ref = fT0(x_ref)
    Target_ref = fTarget(x_ref)
    R_ref = Target_ref - T0_ref

    # Quadrature de Gauss-Legendre
    gauss_x, gauss_w = np.polynomial.legendre.leggauss(nq)
    gauss_x = 0.5 * (gauss_x + 1)  # [0,1]

    def integrate(f_vals):
        """Intégration par Gauss-Legendre"""
        total = 0
        for i in range(len(x_ref)-1):
            xL, xR = x_ref[i], x_ref[i+1]
            xg = xL + (xR - xL) * gauss_x
            fg = np.interp(xg, x_ref, f_vals)
            total += np.dot(fg, gauss_w) * (xR - xL) * 0.5
        return total

    # Calcul de A_ij et B_i
    A = np.zeros((n, n))
    B = np.zeros(n)
    for i in range(n):
        for j in range(i, n):
            Aij = integrate(T_ref[i] * T_ref[j])
            A[i, j] = A[j, i] = Aij
        B[i] = integrate(R_ref * T_ref[i])

    return A, B


# ==========================================================
# ================== PROGRAMME PRINCIPAL ===================
# ==========================================================
if __name__ == "__main__":

    nbc = 4        # nombre de contrôles
    NX = 30
    nb_iter_refine = 1
    best_cost = 1e10

    for iref in range(nb_iter_refine):

        NX += 5
        Target = np.zeros(NX)
        xcible = np.arange(nbc) + 1
        cost_ref, Target, x_target = ADRS(NX, xcible, Target)

        xcontrol = np.zeros(nbc)
        cost, T0, x_T0 = ADRS(NX, xcontrol, Target)

        sol_list, x_list = [], []
        for ic in range(nbc):
            xic = np.zeros(nbc)
            xic[ic] = 1
            cost, Tic, xic_mesh = ADRS(NX, xic, Target)
            sol_list.append(Tic)
            x_list.append(xic_mesh)

        # Calcul Aij, Bi avec interpolation + quadrature
        A, B = compute_A_B(sol_list, [T0], [Target], x_list, nq=6)

        # Résolution du système linéaire
        xopt = np.linalg.solve(A, B)
        print("xopt =", xopt)
        cost_opt, Topt, xopt_mesh = ADRS(NX, xopt, Target)
        print("cost_opt =", cost_opt)

        plt.figure()
        plt.plot(xopt_mesh, Topt, label="Optim")
        plt.plot(x_target, Target, "--", label="Target")
        plt.xlabel("x")
        plt.ylabel("T(x)")
        plt.legend()
        plt.title("Solution optimisée vs Target (maillage adaptatif)")
        plt.show()
