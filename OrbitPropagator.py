import numpy as np
from numpy import ceil
from numpy.linalg import norm
import matplotlib.pyplot as plt
from scipy.integrate import ode

import planetary_data as pd
import Tools
import SpiceTools


def perturbations():
    return {
        'J2': False,
        'aero': False,
        'cd': 0,
        'A': 0,
        'CR': 0,
        'A_srp': 0,
        'moon_gravity': False,
        'sun_gravity': False,
        'srp': False,
        'thrust': 0,
        'isp': 0,
        'n_bodies': []
    }


class OrbitPropagator:
    def __init__(self, state0, tspan, dt, kep=False, deg=True, cb=pd.earth, perts=perturbations()):
        if kep:
            self.r0, self.v0 = Tools.kep2car(state0, deg=deg, mu=cb['mu'])
        else:
            self.r0 = state0[:3]
            self.v0 = state0[3:]

        self.tspan = tspan
        self.dt = dt
        self.cb = cb

        # Total number of steps
        self.n_steps = int(ceil(self.tspan / self.dt)) + 1

        # Initialize arrays
        self.ys = np.zeros((self.n_steps, 6))
        self.ts = np.zeros((self.n_steps, 1))
        self.rs = np.zeros((self.n_steps, 3))
        self.vs = np.zeros((self.n_steps, 3))
        self.kep = np.zeros((self.n_steps, 6))

        # Initial conditions
        self.y0 = np.concatenate((self.r0, self.v0))
        self.ys[0] = np.array(self.y0)
        self.steps = 1

        # Initialize the solver
        self.solver = ode(self.dynamics)
        self.solver.set_integrator('lsoda')
        self.solver.set_initial_value(self.y0, 0)

        # Define perturbations dictionary
        self.perts = perts

        self.propagate_orbit()

    def propagate_orbit(self):
        # Propagate the orbit
        while self.solver.successful() and self.steps < self.n_steps:
            self.solver.integrate(self.solver.t + self.dt)
            self.ts[self.steps] = self.solver.t
            self.ys[self.steps] = self.solver.y
            self.steps += 1

        self.rs = self.ys[:, :3]
        self.vs = self.ys[:, 3:]

    def dynamics(self, t, y):
        # Unpack the state
        rx, ry, rz, vx, vy, vz = y
        r = np.array([rx, ry, rz])
        norm_r = norm(r)

        # Two-Body acceleration
        a = -r * self.cb['mu'] / norm_r ** 3

        # J2 perturbation:
        if self.perts['J2']:
            z2 = r[2] ** 2
            r2 = norm_r ** 2
            tx = r[0] / norm_r * (5 * z2 / r2 - 1)
            ty = r[1] / norm_r * (5 * z2 / r2 - 1)
            tz = r[2] / norm_r * (5 * z2 / r2 - 3)

            a_j2 = 1.5 * self.cb['J2'] * self.cb['mu'] * self.cb['radius'] ** 2 / norm_r ** 4 * np.array([tx, ty, tz])

            a += a_j2

        return [vx, vy, vz, a[0], a[1], a[2]]

    def calculate_kep(self, deg=True):
        print('Calculating keplerian elements...')

        for n in range(self.steps):
            self.kep[n, :] = SpiceTools.rv2kep(self.ys[n, :6], mu=self.cb['mu'], deg=deg)

    def plot(self, show_plot=False, save_plot=False, title='Test Title', k=1):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111, projection='3d')

        # Plot trajectory
        ax.plot(self.rs[:, 0], self.rs[:, 1], self.rs[:, 2], 'r', label='Trajectory')
        ax.scatter3D(self.rs[0, 0], self.rs[0, 1], self.rs[0, 2], 'wo', label='Initial conditions')

        # Plot Central Body:
        _u, _v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
        _x = self.cb['radius'] * np.cos(_u) * np.sin(_v) * k
        _y = self.cb['radius'] * np.sin(_u) * np.sin(_v) * k
        _z = self.cb['radius'] * np.cos(_v) * k
        ax.plot_surface(_x, _y, _z, cmap='Blues')

        # Plot the x,y,z axis:
        x, y, z = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        u, v, w = [[2 * self.cb['radius'], 0, 0], [0, 2 * self.cb['radius'], 0], [0, 0, 2 * self.cb['radius']]]

        ax.quiver(x, y, z, u, v, w, color='k')

        max_val = np.max(np.abs(self.rs))
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
            plt.savefig(title+'png', dpi=300)

    def plot_kep(self, hours=False, days=False, show_plot=False, save_plot=False,
                 title='Kep', figsize=(12, 8)):
        print('Plotting keplerian elements...')

        # create figure and axis instances:
        fig, axs = plt.subplots(nrows=2, ncols=3, figsize=figsize)

        # figure title:
        fig.suptitle(title, fontsize=20)

        # x-axis:
        if hours:
            ts = self.ts / 3600
            xlabel = 'Time Elapsed (hours)'
        elif days:
            ts = self.ts / 3600 / 24
            xlabel = 'Time Elapsed (days)'
        else:
            ts = self.ts
            xlabel = 'Time Elapsed (s)'

        # plot true anomaly:
        axs[0, 0].plot(ts, self.kep[:, 5])
        axs[0, 0].set_title('True anomaly vs Time')
        axs[0, 0].grid(True)
        axs[0, 0].set_ylabel('Angle (deg)')

        # plot semi-major axis:
        axs[1, 0].plot(ts, self.kep[:, 0])
        axs[1, 0].set_title('Semi-Major Axis vs Time')
        axs[1, 0].grid(True)
        axs[1, 0].set_ylabel('Semi-Major Axis (km)')
        axs[1, 0].set_xlabel(xlabel)

        # plot eccentricity:
        axs[0, 1].plot(ts, self.kep[:, 1])
        axs[0, 1].set_title('Eccentricity vs Time')
        axs[0, 1].grid(True)
        axs[0, 1].set_xlabel(xlabel)

        # plot inclination:
        axs[1, 1].plot(ts, self.kep[:, 2])
        axs[1, 1].set_title('Inclination vs Time')
        axs[1, 1].grid(True)
        axs[1, 1].set_ylabel('Angle (deg)')
        axs[1, 1].set_xlabel(xlabel)

        # plot RAAN:
        axs[0, 2].plot(ts, self.kep[:, 3])
        axs[0, 2].set_title('RAAN vs Time')
        axs[0, 2].grid(True)
        axs[0, 2].set_ylabel('Angle (deg)')
        axs[0, 2].set_xlabel(xlabel)

        # plot argument of perigee:
        axs[1, 2].plot(ts, self.kep[:, 4])
        axs[1, 2].set_title('Argument of Perigee vs Time')
        axs[1, 2].grid(True)
        axs[1, 2].set_ylabel('Angle (deg)')
        axs[1, 2].set_xlabel(xlabel)

        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title + '.png', dpi=300)
