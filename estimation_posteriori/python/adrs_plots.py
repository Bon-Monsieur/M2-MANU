import math
import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------------------------
# 1. Fonction donnant la solution exacte et le terme source
# --------------------------------------------------------------------
def fex(NX, dx, time):
    F = np.zeros(NX)
    Tex = np.zeros(NX)
    Text = np.zeros(NX)
    Texx = np.zeros(NX)
    for j in range(1, NX-1):
        v = (np.exp(-1000*((j-NX/3)/NX)**2)
             + np.exp(-10*np.exp(-1000*((j-NX/3)/NX)**2))) \
            * np.sin(5*j*math.pi/NX)
        Tex[j]  = np.sin(4*math.pi*time) * v
        Text[j] = 4*math.pi*np.cos(4*math.pi*time) * v
    for j in range(1, NX-1):
        Texx[j] = (Tex[j+1]-Tex[j-1])/(2*dx)
        Txx     = (Tex[j+1]-2*Tex[j]+Tex[j-1])/(dx**2)
        F[j]    = V*Texx[j] - K*Txx + lamda*Tex[j] + Text[j]
    return F, Tex, Texx

# --------------------------------------------------------------------
# 2. Paramètres physiques
# --------------------------------------------------------------------
K      = 0.1      # diffusion
L      = 1.0      # domaine [0,1]
V      = 1.0      # advection
lamda  = 1.0      # réaction
Time   = 1.0      # temps final
ifre   = 100      # fréquence d’affichage
eps    = 1e-3
irk_max = 4

# --------------------------------------------------------------------
# 3. Boucle de raffinement : Erreur L2 aux temps 0.5 et 1.0
# --------------------------------------------------------------------
NX_init = 5
niter_refinement = 20
NX_tab, Err_tab1, Err_tab2 = [], [], []
error = np.zeros(niter_refinement)

for iter_ref in range(niter_refinement):
    NX = NX_init + 3*(iter_ref+1)
    NX_tab.append(NX)
    dx = L/(NX-1)
    dt = dx**2/(V*dx + K + dx**2)

    x = np.linspace(0.0, 1.0, NX)
    T = np.zeros(NX)
    time = 0.0
    n = 0
    rest = []
    time_tab = []

    while time < Time:
        n += 1
        F, Tex, Texx = fex(NX, dx, time)
        dt = dx**2/(V*dx + 2*K + abs(np.max(F))*dx**2)
        time += dt
        time_tab.append(time)
        T0 = T.copy()

        alpha = np.array([1/(irk_max - k) for k in range(irk_max)])
        for irk in range(irk_max):
            for j in range(1, NX-1):
                xnu = K + 0.5*dx*abs(V)
                Tx  = (T[j+1]-T[j-1])/(2*dx)
                Txx = (T[j-1]-2*T[j]+T[j+1])/(dx**2)
                RHS = dt*(-V*Tx + xnu*Txx - lamda*T[j] + F[j])
                T[j] = T0[j] + RHS*alpha[irk]

        res = np.sum(np.abs(T - T0))
        rest.append(res)

        # erreurs L2
        err   = np.dot(T-Tex, T-Tex)*dx
        errh1 = 0.0
        for j in range(1, NX-1):
            errh1 += dx*(Texx[j]-(T[j+1]-T[j-1])/(2*dx))**2
        error[iter_ref] = np.sqrt(err)/NX

        if abs(time - 0.3) < dt*0.5:
            Err_tab1.append(error[iter_ref])
        if abs(time - 0.6) < dt*0.5:
            Err_tab2.append(error[iter_ref])

    if iter_ref == niter_refinement - 1:
        # Sauvegarde de la figure 2 : résidu (sans légende)
        plt.figure()
        plt.plot(time_tab, rest, color='blue')
        plt.title("Évolution du résidu au cours du temps")
        plt.xlabel("Temps t")
        plt.ylabel("Résidu")
        plt.grid(True)
        plt.savefig("fig2_Residus.png", dpi=300)

# --------------------------------------------------------------------
# 4. Figure 1 : T(x,t) aux deux instants 0.5 et 1.0 (seulement 2 courbes)
# --------------------------------------------------------------------
NX_plot = 101
dx = L/(NX_plot-1)
x  = np.linspace(0, 1, NX_plot)

def solution_num_exact(time_target):
    """Renvoie la solution numérique à un temps donné."""
    T = np.zeros(NX_plot)
    time = 0.0
    dt = dx**2/(V*dx + K + dx**2)
    while time < time_target:
        F, Tex, Texx = fex(NX_plot, dx, time)
        dt = dx**2/(V*dx + 2*K + abs(np.max(F))*dx**2)
        time += dt
        T0 = T.copy()
        alpha = np.array([1/(irk_max - k) for k in range(irk_max)])
        for irk in range(irk_max):
            for j in range(1, NX_plot-1):
                xnu = K + 0.5*dx*abs(V)
                Tx  = (T[j+1]-T[j-1])/(2*dx)
                Txx = (T[j-1]-2*T[j]+T[j+1])/(dx**2)
                RHS = dt*(-V*Tx + xnu*Txx - lamda*T[j] + F[j])
                T[j] = T0[j] + RHS*alpha[irk]
    F, Tex, Texx = fex(NX_plot, dx, time_target)
    return x, T, Tex

x, T05, Tex05 = solution_num_exact(0.3)
x, T1 , Tex1  = solution_num_exact(0.6)

plt.figure()
plt.plot(x, T05, 'b-', label='Numérique t=0.3')
plt.plot(x, Tex05, 'b--', label='Exacte t=0.3')
plt.plot(x, T1 , 'r-', label='Numérique t=0.6')
plt.plot(x, Tex1 , 'r--', label='Exacte t=0.6')
plt.title("Évolution spatiale de T(x,t) – t=0.3 et t=0.6")
plt.xlabel("x"); plt.ylabel("T")
plt.legend()
plt.grid(True)
plt.savefig("fig1_Temporel_2courbes.png", dpi=300)

# --------------------------------------------------------------------
# 5. Figure 3 : Erreur L2 vs maillage
# --------------------------------------------------------------------
plt.figure()
plt.plot(Err_tab1, NX_tab, 'o-', label="t = 0.3 s")
plt.plot(Err_tab2, NX_tab, 's-', label="t = 0.6 s")
plt.title("Erreur L² en fonction du nombre de points de grille")
plt.xlabel("Erreur L²")
plt.ylabel("Nombre de points NX")
plt.legend()
plt.grid(True)
plt.savefig("fig3_ErreurMaillage.png", dpi=300)

# --------------------------------------------------------------------
# 6. Figure 4 : Erreur au centre pour différents ordres de Runge–Kutta
# --------------------------------------------------------------------
def run_rk_order(rk_order):
    alpha = np.array([1.0/(rk_order - k) for k in range(rk_order)])
    T = np.zeros(NX_plot)
    time = 0.0
    dt = dx**2/(V*dx + K + dx**2)
    idx_mid = np.argmin(np.abs(x - 0.5))
    time_tab, err_mid = [], []
    while time < Time:
        F, Tex, Texx = fex(NX_plot, dx, time)
        dt = dx**2/(V*dx + 2*K + abs(np.max(F))*dx**2)
        time += dt
        time_tab.append(time)
        T0 = T.copy()
        for irk in range(rk_order):
            for j in range(1, NX_plot-1):
                xnu = K + 0.5*dx*abs(V)
                Tx  = (T[j+1]-T[j-1])/(2*dx)
                Txx = (T[j-1]-2*T[j]+T[j+1])/(dx**2)
                RHS = dt*(-V*Tx + xnu*Txx - lamda*T[j] + F[j])
                T[j] = T0[j] + RHS*alpha[irk]
        err_mid.append(abs(T[idx_mid] - Tex[idx_mid]))
    return time_tab, err_mid

plt.figure()
for order in [1,2,3,4]:
    ttab, emid = run_rk_order(order)
    plt.plot(ttab, emid, label=f"RK{order}")
plt.title("Erreur au centre du domaine pour différents ordres de Runge–Kutta")
plt.xlabel("Temps t")
plt.ylabel("Erreur absolue en x=0.5")
plt.legend()
plt.grid(True)
plt.savefig("fig4_ErreurCentre.png", dpi=300)

print("Figures générées :")
print("  fig1_Temporel_2courbes.png")
print("  fig2_Residus.png")
print("  fig3_ErreurMaillage.png")
print("  fig4_ErreurCentre.png")
plt.show()