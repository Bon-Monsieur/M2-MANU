import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# PHYSICAL PARAMETERS
K = 0.1      # Diffusion coefficient
Lx, Ly = 1.0, 1.0  # Domain size in x and y
Time = 20.   # Integration time
Vx, Vy = 1.0, 1.0  # Advection velocities
lamda = 1.0  # Reaction rate

# NUMERICAL PARAMETERS
NX, NY = 10, 10  # Initial number of grid points in x and y
NT = 10000      # Max number of time steps
eps = 0.001     # Relative convergence ratio
ifre = 100      # Plot every ifre time iterations

# Grid step (space)
dx = Lx / (NX - 1)
dy = Ly / (NY - 1)

# Time step (CFL condition)
dt = 0.5 * min(dx, dy)**2 / (2 * K + Vx * dx + Vy * dy)

# Initialization
x = np.linspace(0.0, Lx, NX)
y = np.linspace(0.0, Ly, NY)
X, Y = np.meshgrid(x, y)

# Solution exacte (2D sine wave)
Tex = np.sin(2 * np.pi * X) * np.sin(2 * np.pi * Y)

# Source term F (calculated to match the exact solution)
F = np.zeros((NY, NX))
for i in range(1, NX - 1):
    for j in range(1, NY - 1):
        Tex_xx = (Tex[j, i+1] - 2 * Tex[j, i] + Tex[j, i-1]) / (dx**2)
        Tex_yy = (Tex[j+1, i] - 2 * Tex[j, i] + Tex[j-1, i]) / (dy**2)
        Tex_x = (Tex[j, i+1] - Tex[j, i-1]) / (2 * dx)
        Tex_y = (Tex[j+1, i] - Tex[j-1, i]) / (2 * dy)
        F[j, i] = Vx * Tex_x + Vy * Tex_y - K * (Tex_xx + Tex_yy) + lamda * Tex[j, i]

# Initial condition (zero everywhere)
T = np.zeros((NY, NX))

# Main time loop
n = 0
res = 1
res0 = 1
rest = []

while n < NT and res / res0 > eps:
    n += 1
    res = 0
    T_new = T.copy()

    # Update interior points
    for i in range(1, NX - 1):
        for j in range(1, NY - 1):
            Tx = (T[j, i+1] - T[j, i-1]) / (2 * dx)
            Ty = (T[j+1, i] - T[j-1, i]) / (2 * dy)
            Txx = (T[j, i-1] - 2 * T[j, i] + T[j, i+1]) / (dx**2)
            Tyy = (T[j-1, i] - 2 * T[j, i] + T[j+1, i]) / (dy**2)
            T_new[j, i] = T[j, i] + dt * (-Vx * Tx - Vy * Ty + K * (Txx + Tyy) - lamda * T[j, i] + F[j, i])
            res += abs(T_new[j, i] - T[j, i])

    T = T_new
    if n == 1:
        res0 = res
    rest.append(res)

    # Plot every ifre steps or at convergence
    if n % ifre == 0 or (res / res0) < eps:
        print(f"Iteration {n}, Residual: {res:.4e}")

# Calculate L2 error
error_l2 = np.sqrt(np.sum((T - Tex)**2) / (NX * NY))
print(f"L2 Error: {error_l2:.4e}")

# Plotting
fig = plt.figure(figsize=(14, 6))

# Exact solution
ax1 = fig.add_subplot(121, projection='3d')
ax1.plot_surface(X, Y, Tex, cmap='viridis')
ax1.set_title("Solution exacte")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_zlabel("u(x, y)")

# Numerical solution
ax2 = fig.add_subplot(122, projection='3d')
ax2.plot_surface(X, Y, T, cmap='viridis')
ax2.set_title("Solution numérique")
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_zlabel("u(x, y)")

plt.tight_layout()
plt.show()
