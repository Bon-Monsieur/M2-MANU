import numpy as np
import matplotlib.pyplot as plt

# --- Définition du signal ---
def f(x):
    return np.sin(2 * np.pi * x) + 0.5 * np.sin(4 * np.pi * x)

# --- Intégration numérique ---
def integrale(f, a, b, n):
    x = np.linspace(a, b, n)
    y = f(x)
    return (b - a) / n * np.sum(y)

# --- Transformée de Fourier directe ---
def fourier_1d(f, nb_xi):
    res = np.zeros(len(nb_xi), dtype=complex)
    for k in range(len(nb_xi)):
        res[k] = integrale(lambda x: f(x) * np.exp(-2j * np.pi * nb_xi[k] * x), -1, 1, 4000)
    return res

# --- Reconstruction à partir des coefficients ---
def fourier_reconstruct_from_coeff(coeffs, nb_xi, x, L=2):
    f_recons = np.zeros_like(x, dtype=complex)
    for k in range(len(nb_xi)):
        f_recons += coeffs[k] * np.exp(2j * np.pi * nb_xi[k] * x)
    return (f_recons / L).real

# --- Paramètres ---
nb_xi = np.arange(-5, 6)
coeffs = fourier_1d(f, nb_xi)

# --- Reconstruction du signal ---
x = np.linspace(-1, 1, 500)
f_recons = fourier_reconstruct_from_coeff(coeffs, nb_xi, x, L=2)

# --- Tracés ---
fig, axs = plt.subplots(3, 1, figsize=(8, 8))

# 1️⃣ Signal original
axs[0].plot(x, f(x), label="f(x) original")
axs[0].set_title("Signal initial")
axs[0].legend()

# 2️⃣ Spectre des coefficients
axs[1].stem(nb_xi, np.abs(coeffs))
axs[1].set_title("Spectre de Fourier |f̂(ξ)|")
axs[1].set_xlabel("ξ")
axs[1].set_ylabel("|f̂(ξ)|")

# 3️⃣ Reconstruction
axs[2].plot(x, f(x), 'k', label="f(x) original", alpha=0.6)
axs[2].plot(x, f_recons, 'r--', label="f reconstruit")
axs[2].set_title("Comparaison signal original / reconstruit")
axs[2].legend()

plt.tight_layout()
plt.show()
