import numpy as np
from numpy import ceil
from numpy.linalg import norm
import matplotlib.pyplot as plt
from scipy.integrate import ode
import spiceypy as spice

from Tools.Tools import calc_atmospheric_density, kep2car, esc_v, state2period
from Tools import SpiceTools, planetary_data as pd

g0 = 9.81
km2AU = 149598073


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
    def __init__(self, state0, tspan, dt, et0=None, kep=False, deg=True, cb=pd.earth,
                 perts=perturbations(), mass=0, stop={}, propagator='lsoda',
                 frame='J2000', date0='2022-12-02'):
        if kep:
            self.r0, self.v0 = kep2car(state0, deg=deg, mu=cb['mu'])
        else:
            self.r0 = state0[:3]
            self.v0 = state0[3:]

        if type(tspan) == str:
            self.tspan = float(tspan) * state2period(state0, coes=kep, mu=self.cb['mu'])
        else:
            self.tspan = tspan

        if et0 is not None:
            self.et0 = et0
        else:
            self.et0 = spice.str2et(date0)

        self.dt = dt
        self.cb = cb
        self.mass = mass
        self.stop = stop
        self.propagator = propagator
        self.frame = frame
        self.date0 = date0

        self.ets = np.arange(self.et0, self.et0 + self.tspan + self.dt, self.dt)

        self.coes_calculated = False
        self.latlons_calculated = False
        self.altitudes_calculated = False
        self.ra_rp_calculated = False
        self.eclipses_calculated = False

        # Total number of steps
        self.n_steps = int(ceil(self.tspan / self.dt)) + 1

        # Initialize arrays
        self.ts = np.zeros((self.n_steps, 1))
        self.rs = np.zeros((self.n_steps, 3))
        self.vs = np.zeros((self.n_steps, 3))
        self.kep = np.zeros((self.n_steps, 6))
        self.alts = np.zeros((self.n_steps, 1))
        self.lat_longs = np.zeros((self.n_steps, 3))

        # Initial conditions
        if perts['thrust']:
            self.ys = np.zeros((self.n_steps, 7))
            self.y0 = np.concatenate((self.r0, self.v0, [self.mass]))
        else:
            self.ys = np.zeros((self.n_steps, 6))
            self.y0 = np.concatenate((self.r0, self.v0))
        self.ys[0] = np.array(self.y0)
        self.alts[0] = norm(self.r0) - self.cb['radius']
        self.steps = 1
        # self.ets = np.zeros((self.n_steps, 1))

        # Initialize the solver
        self.solver = ode(self.dynamics)
        self.solver.set_integrator(self.propagator)
        self.solver.set_initial_value(self.y0, 0)

        # List of Spice files
        self.spice_files_loaded = []

        # Define perturbations dictionary
        self.perts = perts

        # Define dictionary to map internal methods:
        self.stop_condition_map = {'max_alt': self.check_max_alt, 'min_alt': self.check_min_alt,
                                   'escape_velocity': self.check_escape_velocity,
                                   'enter_SOI': self.check_enter_SOI, 'exit_SOI': self.check_exit_SOI}
        # Create stop conditions function list with de-orbit always checked:
        self.stop_condition_functions = [self.check_deorbit]
        # Fill in the rest of stop conditions:
        for key in self.stop.keys():
            self.stop_condition_functions.append(self.stop_condition_map[key])

        # load leap seconds kernel:
        leap_seconds = 'Data/spice/lsk/naif0012.tls'

        spice.furnsh(leap_seconds)

        # add to list of loaded spice files:
        self.spice_files_loaded.append(leap_seconds)

        # convert start date to seconds after J2000:
        self.start_time = spice.utc2et(self.date0)

        # create time span array in seconds after J2000:
        self.spice_tspan = np.linspace(self.start_time, self.start_time + self.tspan, self.n_steps)

        # check if loading in spice data:
        if self.perts['n_bodies'] or self.perts['srp']:

            # if srp, get states of the Sun:
            if self.perts['srp']:
                # load spice file for given body:
                spice.furnsh(self.cb['spice_file'])

                # add to spice file list:
                self.spice_files_loaded.append(self.cb['spice_file'])

                # calculate central body's state throughout the entire propagation wrt the Sun:
                self.cb['states'] = SpiceTools.get_ephemeris(self.cb['name'], self.spice_tspan, self.frame, 'SUN')

        # load kernels for each body:
        for body in self.perts['n_bodies']:
            # if spice hasn't already been loaded:
            if body['spice_file'] not in self.spice_files_loaded:
                spice.furnsh(body['spice_file'])

                # add spice file to the list:
                self.spice_files_loaded.append(body['spice_file'])

            # calculate body states wrt central body:
            body['states'] = SpiceTools.get_ephemeris(body['name'], self.spice_tspan, self.frame, self.cb['spice_name'])

        # ----- Propagate the orbit ----- #
        self.propagate_orbit()

    def propagate_orbit(self):
        # Propagate the orbit
        while self.solver.successful() and self.steps < self.n_steps:
            self.solver.integrate(self.solver.t + self.dt)
            self.ts[self.steps] = self.solver.t
            self.ys[self.steps] = self.solver.y
            self.alts[self.steps] = norm(self.solver.y[:3]) - self.cb['radius']

            if self.check_stop_condition():
                self.steps += 1
            else:
                break

        self.ys = self.ys[:self.steps]
        self.rs = self.ys[:, :3]
        self.vs = self.ys[:, 3:]
        self.alts = self.alts[:self.steps]
        self.ets = self.ets[:self.steps]
        self.states = self.ys

    def dynamics(self, t, y):
        if self.perts['thrust']:
            # Unpack the state
            rx, ry, rz, vx, vy, vz, mass = y
        else:
            # Unpack the state
            rx, ry, rz, vx, vy, vz = y

        r = np.array([rx, ry, rz])
        v = np.array([vx, vy, vz])
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

        # Aerodynamic drag:
        if self.perts['aero']:
            # Calculate altitude and air density:
            z = norm_r - self.cb['radius']
            rho = calc_atmospheric_density(z)

            # Calculate motion of the spacecraft relative to rotating atmosphere:
            v_rel = v - np.cross(self.cb['atm_rot_vector'], r)

            # Compute drag force:
            drag = -v_rel * 0.5 * rho * norm(v_rel) * self.perts['cd'] * self.perts['A'] / self.mass

            a += drag

        # N-bodies perturbation:
        for body in self.perts['n_bodies']:
            # Vector pointing from central body to n-body:
            r_cb2nb = body['states'][self.steps, :3]

            # Vector pointing from satellite to body:
            r_sat2body = r_cb2nb - r

            # nth body acceleration vector:
            a += body['mu'] * (
                 r_sat2body / norm(r_sat2body) ** 3 - r_cb2nb / norm(r_cb2nb) ** 3)

        # Solar radiation pressure:
        if self.perts['srp']:
            # Vector pointing from the Sun to the spacecraft:
            r_sun2sc = self.cb['states'][self.steps, :3] + r

            # srp acceleration:
            a += (1 + self.perts['CR']) * pd.sun['G1'] * self.perts['A_srp'] / self.mass / norm(
                 r_sun2sc) ** 3 * r_sun2sc

        # Thrust action:
        if self.perts['thrust']:
            # Thrust vector
            a += self.perts['thrust_direction'] * (np.array(v) / norm(v)) \
                     * self.perts['thrust'] / self.mass / 1000  # km/s^2

            # Derivative of total mass
            dmdt = -self.perts['thrust'] / self.perts['isp'] / g0

            return [vx, vy, vz, a[0], a[1], a[2], dmdt]

        return [vx, vy, vz, a[0], a[1], a[2]]

    def check_deorbit(self):
        if self.alts[self.steps] < self.cb['deorbit_altitude']:
            print('Spacecraft de-orbited after %.1f seconds' % self.ts[self.steps])
            return False
        else:
            return True

    def check_max_alt(self):
        if self.alts[self.steps] > self.stop['max_alt']:
            print('Spacecraft reached maximum altitude after %.1f seconds' % self.ts[self.steps])
            return False
        else:
            return True

    def check_min_alt(self):
        if self.alts[self.steps] < self.stop['min_alt']:
            print('Spacecraft reached minimum altitude after %.1f seconds' % self.ts[self.steps])
            return False
        else:
            return True

    def check_escape_velocity(self):
        # if escape velocity is less than current velocity norm:
        if esc_v(np.linalg.norm(self.ys[self.steps, :3]),
                   mu=self.cb['mu']) < np.linalg.norm(self.ys[self.steps, 3:6]):
            print('Spacecraft reached escape velocity after %.1f seconds' % self.ts[self.steps])
            return False
        return True

    def check_enter_SOI(self):
        body = self.stop['enter_SOI']
        r_cb2body = spice.spkgps(body['SPICE_ID'], self.ets[self.steps],
                                 self.frame,  self.cb['SPICE_ID'])[0]
        r_sc2body = r_cb2body - self.ys[self.steps, :3]

        if norm(r_sc2body) < body['SOI']:
            print(f"Spacecraft enters SOI of {body['name']}")
            return False
        return True

    def check_exit_SOI(self):
        if norm(self.ys[self.steps, :3]) > self.cb['SOI']:
            print(f"Spacecraft exits SOI of {self.cb['name']}")
            return False
        return True

    def check_stop_condition(self):
        for stop_condition in self.stop_condition_functions:
            if not stop_condition():
                return False
        return True

    def calculate_kep(self, deg=True):
        print('Calculating keplerian elements...')

        for n in range(self.steps):
            self.kep[n, :] = SpiceTools.rv2kep(self.ys[n, :6], mu=self.cb['mu'], deg=deg)

    def plot(self, cmap='Blues', AU=False, show_plot=False, save_plot=False, title='Test title', k=1):
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Plot Trajectory:
        max_val = 0
        if AU:
            self.rs /= km2AU
            self.cb['radius'] /= km2AU

        ax.plot(self.rs[:, 0], self.rs[:, 1], self.rs[:, 2], 'r', label='Trajectory', zorder=3)
        ax.scatter3D(self.rs[0, 0], self.rs[0, 1], self.rs[0, 2], 'wo', label='Initial Position')

        max__ = np.max(self.rs)
        if max__ > max_val:
            max_val = max__

        # Plot Central Body:
        _u, _v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
        _x = self.cb['radius'] * np.cos(_u) * np.sin(_v) * k
        _y = self.cb['radius'] * np.sin(_u) * np.sin(_v) * k
        _z = self.cb['radius'] * np.cos(_v) * k
        ax.plot_surface(_x, _y, _z, cmap=cmap, zorder=0)

        # Plot the x,y,z axis:
        x, y, z = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        u, v, w = [[2 * self.cb['radius'], 0, 0], [0, 2 * self.cb['radius'], 0], [0, 0, 2 * self.cb['radius']]]

        ax.quiver(x, y, z, u, v, w, color='k')

        # max_val = np.max(np.abs(self.rs))
        ax.set_xlim([-max_val, max_val])
        ax.set_ylim([-max_val, max_val])
        ax.set_zlim([-max_val, max_val])

        # Set axis label:
        if AU:
            ax.set_xlabel('X (AU)')
            ax.set_ylabel('Y (AU)')
            ax.set_zlabel('Z (AU)')
        else:
            ax.set_xlabel('X (km)')
            ax.set_ylabel('Y (km)')
            ax.set_zlabel('Z (km)')

        ax.set_title(title)

        plt.legend()
        if show_plot:
            plt.show()

        if save_plot:
            plt.savefig(title + '.png', dpi=300)

    def plot_altitudes(self, show_plot=False, save_plot=False, hours=False, days=False,
                       title='Radial distance over time', figsize=(12, 8), dpi=300):
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

        plt.figure(figsize=figsize)
        plt.plot(ts, self.alts, 'w')
        plt.grid(True)
        plt.xlabel(xlabel)
        plt.ylabel('Altitude (km)')
        plt.title(title)
        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title + '.png', dpi=dpi)

    def calculate_apo_peri(self):
        self.apoapses = self.kep[:, 0] * (1 + self.kep[:, 1])
        self.periapses = self.kep[:, 0] * (1 - self.kep[:, 1])

    def plot_apo_peri(self, show_plot=False, hours=False, days=False, save_plot=False,
                      title='Apoapse and Periapse over Time', dpi=300):
        plt.figure(figsize=(16, 8))

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

        plt.plot(ts, self.apoapses, 'b', label='Apoapse')
        plt.plot(ts, self.periapses, 'r', label='Periapse')

        plt.xlabel(xlabel)
        plt.ylabel('Altitude (km)')

        plt.grid(True)
        plt.title(title)
        plt.legend()

        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title + '.png', dpi=dpi)

    def plot_kep(self, hours=False, days=False, show_plot=False, save_plot=False,
                 title='Kep', figsize=(12, 8), dpi=300):
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
            plt.savefig(title + '.png', dpi=dpi)

    def write_traj(self):
        pass

    def calc_lat_lon(self):
        self.lat_longs, self.rs_ecef, _ = SpiceTools.inert2latlong(self.rs, self.spice_tspan, self.frame)
