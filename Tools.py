import matplotlib.pyplot as plt
import numpy as np
from OrbitalMechanics.OrbitalMotionSimulator import planetary_data as pd

d2r = np.pi/180.0


def plot_n_orbits(rs, labels, cb=pd.earth, show_plot=False, save_plot=False, title='Many Orbits', k=1):
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Plot trajectories
    j = 0
    for r in rs:
        ax.plot(r[:, 0], r[:, 1], r[:, 2], label=labels[j])
        ax.scatter3D(r[0, 0], r[0, 1], r[0, 2])
        j += 1


    # Plot Central Body:
    _u, _v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    _x = cb['radius'] * np.cos(_u) * np.sin(_v) * k
    _y = cb['radius'] * np.sin(_u) * np.sin(_v) * k
    _z = cb['radius'] * np.cos(_v) * k
    ax.plot_surface(_x, _y, _z, cmap='Blues')

    # Plot the x,y,z axis:
    x, y, z = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    u, v, w = [[2 * cb['radius'], 0, 0], [0, 2 * cb['radius'], 0], [0, 0, 2 * cb['radius']]]

    ax.quiver(x, y, z, u, v, w, color='k')

    max_val = np.max(np.abs(rs))
    ax.set_xlim([-max_val, max_val])
    ax.set_ylim([-max_val, max_val])
    ax.set_zlim([-max_val, max_val])

    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')

    # ax.set_aspect('equal')

    ax.set_title(title)
    plt.legend()

    if show_plot:
        plt.show()
    if save_plot:
        plt.savefig(title + 'png', dpi=300)
