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
import numpy as np
from scipy.interpolate import interp1d

def compute_A_B(T_list, T0_list, Target_list, x_list, nq=6, x_target_ref=None):
    """
    Compute the matrix A and vector B for the least-squares optimal control problem.

    Parameters
    ----------
    T_list : list of np.ndarray
        Solutions on possibly adapted meshes
    T0_list : list of np.ndarray
        Initial solutions (aligned with T_list)
    Target_list : list of np.ndarray
        Target solution (can be on a different mesh)
    x_list : list of np.ndarray
        Meshes corresponding to each solution in T_list
    nq : int
        Number of quadrature points for numerical integration
    x_target_ref : np.ndarray or None
        Optional mesh on which Target_list[0] is defined

    Returns
    -------
    A : np.ndarray
        The matrix of the linear system
    B : np.ndarray
        The right-hand side vector
    """

    # --- Define a reference mesh covering all solution meshes ---
    x_ref = np.linspace(min([x.min() for x in x_list]),
                        max([x.max() for x in x_list]),
                        400)  # you can increase points if needed

    # --- Interpolate T0 onto reference mesh ---
    fT0 = interp1d(x_list[0], T0_list[0], kind='cubic', fill_value="extrapolate")
    T0_ref = fT0(x_ref)

    # --- Interpolate Target onto reference mesh ---
    if x_target_ref is None:
        # assume Target_list[0] is defined on [0,1], scale accordingly
        fTarget = interp1d(np.linspace(0, 1, len(Target_list[0])), Target_list[0],
                           kind='cubic', fill_value="extrapolate")
        Target_ref = fTarget(x_ref / x_ref.max())
    else:
        # interpolate using provided mesh
        fTarget = interp1d(x_target_ref, Target_list[0], kind='cubic', fill_value="extrapolate")
        Target_ref = fTarget(x_ref)

    # --- Interpolate all solutions T_list onto x_ref ---
    T_ref_list = []
    for i, T in enumerate(T_list):
        fT = interp1d(x_list[i], T, kind='cubic', fill_value="extrapolate")
        T_ref_list.append(fT(x_ref))

    # --- Numerical integration weights using Gauss-Legendre quadrature ---
    # Map quadrature points from [-1,1] to [x_i, x_{i+1}]
    xi_q, w_q = np.polynomial.legendre.leggauss(nq)
    dx = np.diff(x_ref)
    w_all = []
    x_all = []
    for i, dx_i in enumerate(dx):
        x_i = 0.5 * dx_i * xi_q + 0.5 * (x_ref[i+1] + x_ref[i])
        w_i = 0.5 * dx_i * w_q
        x_all.append(x_i)
        w_all.append(w_i)
    x_all = np.concatenate(x_all)
    w_all = np.concatenate(w_all)

    # --- Interpolate T_ref_list onto quadrature points ---
    T_quad = [np.interp(x_all, x_ref, T_ref) for T_ref in T_ref_list]
    T0_quad = np.interp(x_all, x_ref, T0_ref)
    Target_quad = np.interp(x_all, x_ref, Target_ref)

    # --- Assemble A and B ---
    n_controls = len(T_list)
    A = np.zeros((n_controls, n_controls))
    B = np.zeros(n_controls)

    for i in range(n_controls):
        B[i] = np.sum(w_all * (Target_quad - T0_quad) * T_quad[i])
        for j in range(n_controls):
            A[i, j] = np.sum(w_all * T_quad[i] * T_quad[j])

    return A, B



# ==========================================================
# ================== PROGRAMME PRINCIPAL ===================
# ==========================================================
import numpy as np
import matplotlib.pyplot as plt

# (Les définitions de ADRS et compute_A_B restent inchangées)
# =======================================================================
# ======================== PROGRAMME PRINCIPAL ==========================
# =======================================================================
if __name__ == "__main__":

    nbc = 4        # nombre de contrôles
    NX = 3000
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
        print("xopt (adaptatif) =", xopt)

        # ---- Solution optimisée sur maillage adaptatif ----
        cost_opt, Topt_adapt, x_adapt = ADRS(NX, xopt, Target)
        print("cost_opt (adaptatif) =", cost_opt)

        # <<< === AJOUT : Calcul de la solution sur maillage fixe === >>>
        NX_fine = 5000  # maillage fixe très fin
        cost_ref_fixe, Target_fixe, x_fixe = ADRS(NX_fine, xcible, np.zeros(NX_fine), adapt=False)
        cost_opt_fixe, Topt_fixe, x_fixe = ADRS(NX_fine, xopt, Target_fixe, adapt=False)

        # <<< === AJOUT : Comparaison graphique === >>>
        plt.figure(figsize=(8,5))
        plt.plot(x_adapt, Topt_adapt, label="Optimisé (maillage adaptatif)")
        plt.plot(x_fixe, Topt_fixe, "--", label="Optimisé (maillage fixe)")
        plt.plot(x_target, Target, "k:", label="Target")
        plt.xlabel("x")
        plt.ylabel("T(x)")
        plt.legend()
        plt.title("Comparaison : contrôle optimal sur maillage adaptatif vs fixe")
        plt.grid(True)
        plt.show()

        # <<< === AJOUT : Comparaison quantitative === >>>
        # interpolation sur le même maillage pour évaluer l'erreur
        Topt_fixe_interp = np.interp(x_adapt, x_fixe, Topt_fixe)
        erreur_L2 = np.sqrt(np.trapz((Topt_fixe_interp - Topt_adapt)**2, x_adapt))
        print(f"Erreur L2 entre maillage adaptatif et fixe = {erreur_L2:.3e}")
