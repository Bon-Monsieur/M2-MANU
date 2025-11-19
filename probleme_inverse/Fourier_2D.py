import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import math

# Charger l'image
chemin_image = r"C:\Users\Raphael\Desktop\FAC\M2-MANU\probleme_inverse\img\baboon256.png"
img = Image.open(chemin_image)

# Convertir en niveaux de gris
img_gris = img.convert('L')

# Convertir en matrice numpy
matrice_gris = np.array(img_gris)
M = matrice_gris.shape[0]  # Nombre de lignes
N = matrice_gris.shape[1]  # Nombre de colonnes


# Afficher les informations
print(f"Dimensions de l'image: {matrice_gris.shape}")
print(f"\nMatrice des niveaux de gris:\n{matrice_gris}")


#np.savetxt("matrice_gris.txt", matrice_gris, fmt='%d') # Pour sauvegarder la matrice dans un fichier

def transformee_fourier_2D_manuelle(s):
    """
    Calcule S(u,v) = Σ Σ s(m,n) * e^(-2iπ(u*m/M + v*n/N))
    """
    M, N = s.shape
    S = np.zeros((M, N), dtype=complex)
    
    print("Calcul de la transformée de Fourier (peut prendre du temps)...")
    
    for u in range(M):
        for v in range(N):
            somme = 0
            for m in range(M):
                for n in range(N):
                    # Formule: e^(-2iπ(u*m/M + v*n/N))
                    exponent = -2j * np.pi * (u * m / M + v * n / N)
                    somme += s[m, n] * np.exp(exponent)
            S[u, v] = somme
        
        if (u + 1) % 10 == 0:
            print(f"Progression: {u+1}/{M}")
    
    return S

def transformee_fourier_2D_vectorisee(s):
    """
    Version vectorisée de la transformée de Fourier 2D
    Calcule ligne par ligne pour économiser la mémoire
    """
    M, N = s.shape
    S = np.zeros((M, N), dtype=complex)
    
    print("Calcul de la transformée de Fourier vectorisée...")
    
    # Créer les grilles pour m et n
    m = np.arange(M).reshape(M, 1)
    n = np.arange(N).reshape(1, N)
    
    for u in range(M):
        for v in range(N):
            # Calculer l'exponentielle pour tous les m,n
            exponent = -2j * np.pi * (u * m / M + v * n / N)
            # Multiplier par s et sommer
            S[u, v] = np.sum(s * np.exp(exponent))
        
        if (u + 1) % 20 == 0:
            print(f"Progression: {u+1}/{M}")
    
    return S

def transformee_fourier_2D_numpy(s):
    """
    Utilise l'algorithme FFT de NumPy (Fast Fourier Transform)
    """
    return np.fft.fft2(s)

def transformee_fourier_inverse_2D_numpy(S):
    """
    Utilise l'algorithme FFT inverse de NumPy
    """
    return np.fft.ifft2(S)



import numpy as np

def erreur_pourcentage(img1, img2):
    """
    Calcule le pourcentage d'erreur relative moyenne entre deux images.
    
    Paramètres
    ----------
    img1, img2 : ndarray
        Images à comparer (même taille).
    
    Retour
    ------
    err : float
        Pourcentage d'erreur moyenne.
    """
    img1 = np.array(img1, dtype=float)
    img2 = np.array(img2, dtype=float)

    diff = np.abs(img1 - img2)
    denom = np.maximum(np.abs(img1), 1e-12)  # éviter division par zéro

    erreur = np.mean(diff / denom) * 100
    return erreur



def filtreGaussien(P):
    epsilon = 0.05
    sigma = P*1.0/math.sqrt(-2*math.log(epsilon))
    h = np.zeros((2*P+1,2*P+1))
    som = 0
    for m in range(-P,P+1):
        for n in range(-P,P+1):
            h[m+P][n+P] = math.exp(-(n*n+m*m)/(2*sigma*sigma))
            som += h[m+P][n+P]
    h = h/som
    return h

# Calculer la transformée de Fourier 2D
# S = transformée_fourier_2D_manuelle(matrice_gris)
# S = transformee_fourier_2D_vectorisee(matrice_gris)
S = transformee_fourier_2D_numpy(matrice_gris)
np.savetxt("fourier_S.txt", S, fmt='%d')
S_recons = transformee_fourier_inverse_2D_numpy(S)

# Calculer le spectre de magnitude (pour visualisation)
magnitude_spectrum = np.abs(S)

magnitude_spectrum_centree = np.fft.fftshift(magnitude_spectrum) # Centrer le spectre (mettre les basses fréquences au centre)

magnitude_log = np.log(magnitude_spectrum_centree) # Échelle logarithmique pour mieux voir

test = transformee_fourier_inverse_2D_numpy(S * magnitude_spectrum) # fou
# ===== AFFICHAGE =====
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Image originale
axes[0].imshow(matrice_gris, cmap='gray')
axes[0].set_title('Image originale S(m,n)')
axes[0].axis('off')

# Spectre de magnitude
axes[1].imshow(magnitude_spectrum_centree, cmap='gray')
axes[1].set_title('Transformée de Fourier de S(u,v)')
axes[1].axis('off')

# Spectre en échelle log
axes[2].imshow(test.real, cmap='gray')
axes[2].set_title('Transformée de Fourier inverse')
axes[2].axis('off')

print(f"Erreur relative moyenne: {erreur_pourcentage(matrice_gris, test.real):.6f} %")
plt.tight_layout()
plt.show()
