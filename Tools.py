import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, cos, sin, tan, atan
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


def kep2car(kep, deg=False, mu=pd.earth['mu']):
    if deg:
        a, e, i, om, w, teta = kep
        i *= d2r
        om *= d2r
        w *= d2r
        teta *= d2r
    else:
        a, e, i, om, w, teta = kep

    E = ecc_anomaly([teta, e], 'tae')
    r_norm = a * (1 - e ** 2) / (1 + e * cos(teta))

    # Calculate r and v in the perifocal frame:
    r_perif = r_norm * np.array([cos(teta), sin(teta), 0])
    v_perif = sqrt(mu * a) / r_norm * np.array([-sin(E), cos(E) * sqrt(1 - e ** 2), 0])

    # Rotation matrix from perifocal to ECI:
    perif2eci = np.transpose(eci2perif(om, w, i))

    # Calculate r and v vectors in the ECI frame:
    r = np.dot(perif2eci, r_perif)
    v = np.dot(perif2eci, v_perif)

    return r, v


def ecc_anomaly(arr, method, tol=1e-8):
    if method == 'newton':
        # Newton's method for iteratively finding E
        Me, e = arr
        if Me < np.pi/2:
            E0 = Me + e/2
        else:
            E0 = Me - e

        for n in range(200):  # Arbitrary max number of steps
            ratio = (E0-e*sin(E0)-Me)/(1-e*cos(E0))
            if abs(ratio) < tol:
                if n == 0:
                    return E0
                else:
                    return E1
            else:
                E1 = E0-ratio
                E0 = E1

        # Did not converge
        return False
    elif method == 'tae':
        teta, e = arr
        return 2*atan(sqrt((1-e)/(1+e))*tan(teta/2))
    else:
        print('Invalid method for eccentric anomaly')


def eci2perif(om, w, i):
    row0 = [-sin(om)*cos(i)*sin(w)+cos(om)*cos(w),
            cos(om)*cos(i)*sin(w)+sin(om)*cos(w),
            sin(i)*sin(w)]

    row1 = [-sin(om)*cos(i)*cos(w)-cos(om)*sin(w),
            cos(om)*cos(i)*cos(w)-sin(om)*sin(w),
            sin(i)*cos(w)]

    row2 = [sin(om)*sin(i), -cos(om)*sin(i), cos(i)]

    return np.array([row0, row1, row2])
