"""
integration_pi_methods.py

But: Approximating the integral

    f(x) = 1/sqrt(1-x^2)

over (-1,1), which is exactly equal to pi.

We compare:
 1) Riemann (trapezoidal rule, uniform in x)
 2) Lebesgue-style (uniform in y levels)
 3) Monte Carlo (random sampling of x in (-1,1))

We plot the error vs number of function evaluations.

Usage:
    python integration_pi_methods.py
"""

import numpy as np
import matplotlib.pyplot as plt

# exact integral
def f(x):
    return 1.0/np.sqrt(1.0 - x**2)

I_exact = np.pi

# Riemann (trapezoid rule on uniform grid)
def riemann_uniform(N):
    eps = 1e-12
    x = np.linspace(-1.0 + eps, 1.0 - eps, N)
    y = f(x)
    I = np.trapezoid(y, x)
    return I

# Lebesgue-style integration
# integral f = \int_0^\infty m({f>t}) dt
# For monotone decreasing measure in t.

def measure_above_threshold(x_grid, f_grid, t):
    above = f_grid > t
    if not np.any(above):
        return 0.0
    # measure = length of set {x: f(x) > t}
    # Since function symmetric, we can just count measure of intervals
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
        # adjust boundaries
        if i > 0 and not above[i-1]:
            x0, x1 = x_grid[i-1], x_grid[i]
            y0, y1 = f_grid[i-1], f_grid[i]
            if y1 != y0:
                alpha = (t-y0)/(y1-y0)
                xL = x0 + alpha*(x1-x0)
        if j+1 < N and not above[j+1]:
            x0, x1 = x_grid[j], x_grid[j+1]
            y0, y1 = f_grid[j], f_grid[j+1]
            if y1 != y0:
                alpha = (t-y0)/(y1-y0)
                xR = x0 + alpha*(x1-x0)
        total += (xR-xL)
        i = j+1
    return total

def lebesgue_uniform_in_y(Nlevels, Nx=2000):
    eps = 1e-12
    x_grid = np.linspace(-1.0 + eps, 1.0 - eps, Nx)
    f_grid = f(x_grid)
    tmin, tmax = np.min(f_grid), np.max(f_grid)
    t_levels = np.linspace(tmin, tmax, Nlevels)
    mvals = np.array([measure_above_threshold(x_grid, f_grid, t) for t in t_levels])
    I = np.trapezoid(mvals, t_levels)
    return I

# Monte Carlo integration (uniform sampling in x)
def monte_carlo(N):
    eps = 1e-12
    x = np.random.uniform(-1.0 + eps, 1.0 - eps, N)
    y = f(x)
    I = (2.0/N) * np.sum(y)
    return I

# ---------------- main ----------------
def main():
    Nvals = np.unique(np.logspace(1, 4, 20, dtype=int))
    errors_riemann = []
    errors_lebesgue = []
    errors_mc = []

    for N in Nvals:
        I_r = riemann_uniform(N)
        err_r = abs(I_r - I_exact)
        errors_riemann.append(err_r)

        I_l = lebesgue_uniform_in_y(N)
        err_l = abs(I_l - I_exact)
        errors_lebesgue.append(err_l)

        # Monte Carlo: average several trials for stability
        mc_trials = 5
        mc_vals = [monte_carlo(N) for _ in range(mc_trials)]
        I_m = np.mean(mc_vals)
        err_m = abs(I_m - I_exact)
        errors_mc.append(err_m)

    plt.figure(figsize=(7,5))
    plt.loglog(Nvals, errors_riemann, 'o-', label='Riemann (uniform x)')
    plt.loglog(Nvals, errors_lebesgue, 's-', label='Lebesgue (uniform y)')
    plt.loglog(Nvals, errors_mc, '^-', label='Monte Carlo (mean of 5 runs)')
    plt.axhline(1e-3, color='k', linestyle='--', alpha=0.5)
    plt.xlabel('Nombre d\'évaluations de f(x)')
    plt.ylabel('Erreur absolue')
    plt.title('Intégrale de 1/sqrt(1-x^2) sur (-1,1) ≡ π')
    plt.legend()
    plt.grid(True, which='both', ls=':')
    plt.tight_layout()
    plt.savefig('fig_pi_errors.png', dpi=200)
    plt.show()

if __name__ == '__main__':
    main()