"""
integration_riemann_lebesgues_v2.py

Comparaison Riemann vs Lebesgue pour

    f(x) = a*x**2 + b*x + c * np.sin(4*pi*x) + 10*exp(-100*(x-0.5)**2)

sur l'intervalle [Left, Right] = [0,1]
avec a=0.5, b=10, c=3.

Objectifs implémentés :
 1) tracé de la fonction
 2) intégration Riemann (maillage uniforme en x)
 3) intégration Lebesgue (maillage uniforme en y)
 4) vérification que l'intégrale ~ 6.94
 5) recherche automatique du plus petit échantillonnage (N pour Riemann, n_levels pour Lebesgue)
    donnant une précision 1e-3
 6) démonstration de stratégies d'adaptation en y (quantiles, adaptatif)

Usage :
    python integration_riemann_lebesgues_v2.py

Sorties : impressions sur la console et figure montrant la fonction et les niveaux y utilisés.

"""

import numpy as np
import matplotlib.pyplot as plt

# ----------------------- définition du problème -----------------------
Left, Right = 0.0, 1.0
a = 0.5
b = 10.0
c = 3.0

def f(x):
    x = np.asarray(x)
    return a*x**2 + b*x + c*np.sin(4*np.pi*x) + 10.0*np.exp(-100.0*(x-0.5)**2)

# ----------------------- Riemann (uniform in x) -----------------------

def riemann_uniform(Left, Right, N):
    x = np.linspace(Left, Right, N)
    y = f(x)
    I = np.trapezoid(y, x)
    return I, x, y

# ----------------------- Lebesgue-style (uniform in y) -----------------------

def measure_above_threshold_on_grid(x_grid, f_grid, t):
    """Mesure (longueur) de l'ensemble {x | f(x) > t} approximée sur (x_grid,f_grid)
    avec interpolation linéaire entre points.
    """
    above = f_grid > t
    if not np.any(above):
        return 0.0
    total = 0.0
    N = len(x_grid)
    i = 0
    while i < N:
        if not above[i]:
            i += 1
            continue
        j = i
        while j+1 < N and above[j+1]:
            j += 1
        xL = x_grid[i]
        xR = x_grid[j]
        # ajuster les extrémités par interpolation si nécessaire
        if i > 0 and not above[i-1]:
            x0, x1 = x_grid[i-1], x_grid[i]
            y0, y1 = f_grid[i-1], f_grid[i]
            if y1 != y0:
                alpha = (t - y0) / (y1 - y0)
                xL = x0 + alpha*(x1 - x0)
        if j+1 < N and not above[j+1]:
            x0, x1 = x_grid[j], x_grid[j+1]
            y0, y1 = f_grid[j], f_grid[j+1]
            if y1 != y0:
                alpha = (t - y0) / (y1 - y0)
                xR = x0 + alpha*(x1 - x0)
        total += max(0.0, xR - xL)
        i = j + 1
    return total


def lebesgue_uniform_in_y(x_grid, f_grid, n_levels):
    """Approche Lebesgue : échantillonnage uniforme des niveaux y entre ymin et ymax
    Retourne l'intégrale approximée et les tableaux t_levels, mpos, mneg.
    """
    ymin = np.min(f_grid)
    ymax = np.max(f_grid)
    t = np.linspace(ymin, ymax, n_levels)
    # On sépare parties positive et négative : on intégrera f+ et f- puis soustraira
    fpos = np.maximum(f_grid, 0.0)
    fneg = np.maximum(-f_grid, 0.0)
    mpos = np.array([measure_above_threshold_on_grid(x_grid, fpos, ti) for ti in t])
    mneg = np.array([measure_above_threshold_on_grid(x_grid, fneg, ti) for ti in t])
    # intégration en t (trapèze)
    Ipos = np.trapezoid(mpos, t)
    Ineg = np.trapezoid(mneg, t)
    return Ipos - Ineg, t, mpos, mneg

# ----------------------- stratégies adaptatives en y -----------------------

def quantile_t_levels_from_samples(f_grid, n_levels):
    """Niveaux t choisis comme quantiles empiriques de f_grid (valeurs réelles).
    On retourne n_levels niveaux uniformément répartis en quantile entre ymin et ymax.
    """
    vals = f_grid.copy()
    q = np.linspace(0.0, 1.0, n_levels)
    t = np.quantile(vals, q)
    # garantir croissant
    t = np.unique(t)
    return t

def lebesgue_with_custom_t(x_grid, f_grid, t_levels):
    fpos = np.maximum(f_grid, 0.0)
    fneg = np.maximum(-f_grid, 0.0)
    mpos = np.array([measure_above_threshold_on_grid(x_grid, fpos, ti) for ti in t_levels])
    mneg = np.array([measure_above_threshold_on_grid(x_grid, fneg, ti) for ti in t_levels])
    # si t_levels non-uniform en espacement, trapèze s'applique directement
    Ipos = np.trapezoid(mpos, t_levels)
    Ineg = np.trapezoid(mneg, t_levels)
    return Ipos - Ineg, t_levels, mpos, mneg

# ----------------------- recherche du plus petit échantillonnage -----------------------

def find_min_N_riemann(tol=1e-3, Nmax=200000):
    # référence de haute précision
    x_ref = np.linspace(Left, Right, 200001)
    I_ref = np.trapezoid(f(x_ref), x_ref)
    N = 10
    while N <= Nmax:
        I, _, _ = riemann_uniform(Left, Right, N)
        if abs(I - I_ref) < tol:
            return N, I, I_ref
        # augmenter N (progression multiplicative)
        N = int(np.ceil(N * 1.3))
    return None, I, I_ref


def find_min_nlevels_lebesgue(tol=1e-3, nlevels_max=2000):
    # pour la mesure, on travaille sur une grille x dense (pour construire m(t) précisément)
    x_dense = np.linspace(Left, Right, 20001)
    f_dense = f(x_dense)
    # référence
    I_ref = np.trapezoid(f_dense, x_dense)
    n = 10
    while n <= nlevels_max:
        # uniform in y
        I, t, mpos, mneg = lebesgue_uniform_in_y(x_dense, f_dense, n)
        if abs(I - I_ref) < tol:
            return n, I, I_ref, t
        n = int(np.ceil(n * 1.3))
    return None, I, I_ref, t

# ----------------------- main démonstration -----------------------

def main():
    # tracer la fonction
    x_plot = np.linspace(Left, Right, 2001)
    y_plot = f(x_plot)

    plt.figure(figsize=(8,4))
    plt.plot(x_plot, y_plot, label='f(x)')
    plt.title(r'f(x) = a x^2 + b x + c sin(4 pi x) + 10 exp(-100(x-0.5)^2)')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Calcul référence numérique pour vérification (très fin)
    x_ref = np.linspace(Left, Right, 200001)
    I_ref = np.trapezoid(f(x_ref), x_ref)
    print(f"Référence (trapezoïde sur {len(x_ref)} points) : I_ref = {I_ref:.10f}")

    # Vérifier la valeur attendue ~ 6.94 (tolérance large)
    print(f"La valeur cherchée (attendue) = 6.94. Erreur réf - 6.94 = {I_ref - 6.94:.6e}")

    # trouver N minimal pour Riemann
    tol = 1e-3
    Nmin, I_Nmin, Iref = find_min_N_riemann(tol=tol)
    if Nmin is not None:
        print(f"Riemann: N minimal pour erreur < {tol} => N = {Nmin}, I = {I_Nmin:.10f}, erreur = {abs(I_Nmin - Iref):.3e}")
    else:
        print(f"Riemann: N minimal non trouvé jusqu'à la limite. Dernière erreur = {abs(I_Nmin - Iref):.3e}")

    # trouver n_levels minimal pour Lebesgue (uniform en y)
    nmin, I_nmin, _, t_used = find_min_nlevels_lebesgue(tol=tol)
    if nmin is not None:
        print(f"Lebesgue (uniform en y): n_levels minimal pour erreur < {tol} => n = {nmin}, I = {I_nmin:.10f}, erreur = {abs(I_nmin - Iref):.3e}")
    else:
        print(f"Lebesgue: n_levels minimal non trouvé jusqu'à la limite. Dernière erreur = {abs(I_nmin - Iref):.3e}")

    # Exemples d'échantillonnage adaptatif en y : quantiles
    x_dense = np.linspace(Left, Right, 20001)
    f_dense = f(x_dense)
    n_levels = 200
    t_quant = quantile_t_levels_from_samples(f_dense, n_levels)
    I_quant, t_q, mpos_q, mneg_q = lebesgue_with_custom_t(x_dense, f_dense, t_quant)
    print(f"Lebesgue (quantile, n={len(t_q)}): I = {I_quant:.10f}, erreur = {abs(I_quant - Iref):.3e}")

    # Affichage comparatif : niveaux y utilisés
    ymin, ymax = np.min(f_dense), np.max(f_dense)
    t_uniform = np.linspace(ymin, ymax, 50)

    plt.figure(figsize=(8,4))
    plt.plot(x_plot, y_plot, label='f(x)')
    # montrer quelques niveaux uniformes
    for tt in t_uniform[::5]:
        plt.axhline(tt, alpha=0.15)
    # montrer niveaux quantiles comme lignes horizontales plus visibles
    for tt in t_q[::max(1,len(t_q)//30)]:
        plt.axhline(tt, alpha=0.5, linestyle='--')
    plt.title('Niveaux y : lignes pleines = uniform en y (ex.), tiretées = quantile')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('fig_f_function.png', dpi=200)
    plt.show()

    # Conseils pour adapter le pas en y (réponse synthétique)
    print('\nConseils pour obtenir un pas "uniforme en y = f(x)" :')
    print(' 1) Échantillonner f sur une grille dense en x pour obtenir un jeu de valeurs y_i = f(x_i).')
    print(' 2) Construire les niveaux t_k comme quantiles empiriques des y_i (séparer f+ et f- si nécessaire).')
    print('    - cela positionne plus de niveaux là où la fonction prend fréquemment ces valeurs (approx. uniforme en distribution de y).')
    print('3) Ou utiliser une stratégie adaptative qui raffine les t_k là où la mesure m({|f|>t}) varie rapidement (critère d\'erreur).')
    print(' 4) En pratique, pour une précision donnée, quantiles donne un bon compromis robustesse/simplicité; adaptatif peut réduire le nombre de niveaux nécessaire.')

if __name__ == '__main__':
    main()
