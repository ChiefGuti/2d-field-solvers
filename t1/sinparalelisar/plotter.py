import numpy as np
import matplotlib.pyplot as plt

L, h = 5.0, 0.02  # box side [m], grid step [m] — must match placas.c

data = np.genfromtxt('output.txt')
size = int(np.sqrt(data.shape[0]))
potential = data.reshape(size, size)

x = np.linspace(0, L, size)
xi, yi = np.meshgrid(x, x)

# np.gradient on (row, col) gives (-dV/dy, -dV/dx) = (Ey, Ex)
Ey, Ex = np.gradient(-potential, h)
Emag = np.hypot(Ex, Ey)

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 6.2), layout='constrained')

im = ax.imshow(potential, extent=[0, L, L, 0], cmap='RdBu_r', vmin=-50, vmax=50)
ax.contour(xi, yi, potential, levels=np.arange(-45, 50, 5),
           colors='k', linewidths=0.4, alpha=0.5)
ax.streamplot(xi, yi, Ex, Ey, color='k', density=1.4, linewidth=0.6,
              arrowsize=0.8)
fig.colorbar(im, ax=ax, label='V  [volt]')
ax.set_title('Potential + field lines')

im2 = ax2.imshow(Emag, extent=[0, L, L, 0], cmap='inferno', vmin=0, vmax=250)
fig.colorbar(im2, ax=ax2, label='|E|  [V/m]')
ax2.set_title('Field magnitude — edge enhancement at the plate tips')

for a in (ax, ax2):
    a.set_xlim(0, L)
    a.set_ylim(L, 0)
    a.set_xlabel('x [m]')
    a.set_ylabel('y [m]')

fig.suptitle('Parallel-plate capacitor, 2D Laplace (Jacobi, 250x250, h=0.02 m)')
fig.savefig('V.png', dpi=130)
