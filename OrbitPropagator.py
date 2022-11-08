import numpy as np
from numpy import ceil
from numpy.linalg import norm
import matplotlib.pyplot as plt
from scipy.integrate import ode

from OrbitalMechanics.OrbitalMotionSimulator import planetary_data as pd
from OrbitalMechanics.OrbitalMotionSimulator import Tools


class OrbitPropagator:
    def __init__(self, state0, tspan, dt, kep=False, deg=True, cb=pd.earth):
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

        # Initial conditions
        self.y0 = np.concatenate((self.r0, self.v0))
        self.ys[0] = np.array(self.y0)
        self.steps = 1

        # Initialize the solver
        self.solver = ode(self.dynamics)
        self.solver.set_integrator('lsoda')
        self.solver.set_initial_value(self.y0, 0)

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
        ax, ay, az = -r * self.cb['mu'] / norm_r ** 3

        return [vx, vy, vz, ax, ay, az]

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
