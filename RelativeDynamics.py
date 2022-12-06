"""
Relative dynamics propagator class definition
"""

# 3rd party libraries
from scipy.integrate import solve_ivp
import numpy as np
from math import sqrt, cos, sin, ceil
import matplotlib.pyplot as plt

from Tools import planetary_data as pd

COLORS = [
    'm', 'deeppink', 'chartreuse', 'w', 'springgreen', 'peachpuff',
    'white', 'lightpink', 'royalblue', 'lime', 'aqua' ] * 100


def _args():
    return {
        'tspan': None,
        'dt': 10,
        'propagator': 'LSODA',
        'atol': 1e-9,
        'rtol': 1e-9,
        'dense_output': False,
        'model': 'CW',
        'mu': pd.earth['mu'],
        'R': 6878,  # h=500km by default
        'a': 6878,  # circular orbit by default
        'e': 0,
    }


def plot_relative_trajectory(rs, args, vectors=[]):
    _args = {
        'figsize': (10, 8),
        'labels': [''] * len(rs),
        'colors': COLORS[:],
        'traj_lws': 3,
        'dist_unit': 'km',
        'axes_mag': 0.8,
        'axes_custom': None,
        'title': 'Relative Trajectories',
        'legend': True,
        'axes_no_fill': True,
        'hide_axes': False,
        'azimuth': False,
        'elevation': False,
        'show': False,
        'filename': False,
        'dpi': 300,
        'vector_colors': [''] * len(vectors),
        'vector_labels': [''] * len(vectors),
        'vector_texts': False
    }
    for key in args.keys():
        _args[key] = args[key]

    fig = plt.figure(figsize=_args['figsize'])
    ax = fig.add_subplot(111, projection='3d')

    max_val = 0
    n = 0
    cs = ['c', 'b', 'r', 'k', 'g', 'm', 'y', 'w', 'r-.']

    for r in rs:
        ax.plot(r[:, 0], r[:, 1], r[:, 2],
                color=_args['colors'][n], label=_args['labels'][n],
                zorder=10, linewidth=_args['traj_lws'])
        ax.plot([r[0, 0]], [r[0, 1]], [r[0, 2]], 'o',
                color=_args['colors'][n])
        ax.scatter3D(0, 0, 0, c='k', marker='*', label='Target')

        max_val = max([r.max(), max_val])
        n += 1

    for vector in vectors:
        ax.quiver(0, 0, 0,
                  vector['r'][0], vector['r'][1], vector['r'][2],
                  color=vector['color'], label=vector['label'])

        if _args['vector_texts']:
            vector['r'] *= _args['vector_text_scale']
            ax.text(vector['r'][0], vector['r'][1], vector['r'][2],
                    vector['label'],
                    color=vector['color'])
    xlabel = 'X (km)'
    ylabel = 'Y (km)'
    zlabel = 'Z (km)'
    if _args['axes_custom'] is not None:
        max_val = _args['axes_custom']
    else:
        max_val *= _args['axes_mag']

    ax.set_xlim([-max_val, max_val])
    ax.set_ylim([-max_val, max_val])
    ax.set_zlim([-max_val, max_val])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.set_box_aspect([1, 1, 1])
    ax.set_aspect('auto')

    if _args['azimuth'] is not False:
        ax.view_init(elev=_args['elevation'],
                     azim=_args['azimuth'])

    if _args['axes_no_fill']:
        ax.w_xaxis.pane.fill = False
        ax.w_yaxis.pane.fill = False
        ax.w_zaxis.pane.fill = False

    if _args['hide_axes']:
        ax.set_axis_off()

    if _args['legend']:
        plt.legend()

    if _args['filename']:
        plt.savefig(_args['filename'], dpi=_args['dpi'])
        print('Saved', _args['filename'])

    if _args['show']:
        plt.show()

    plt.close()


class RelDynPropagtor:
    def __init__(self, args=_args()):

        self.args = args
        for key in args.keys():
            self.args[key] = args[key]

        self.tspan = self.args['tspan']
        self.dt = self.args['dt']
        self.ets = np.arange(0, 0 + self.tspan + self.dt, self.dt)

        # Total number of steps
        self.n_steps = int(ceil(self.tspan / self.dt)) + 1

        if self.args['model'] == 'CW':
            self.states = np.zeros((self.n_steps, 6))
        elif self.args['model'] == 'tillerson':
            self.states = np.zeros((self.n_steps, 7))
        elif self.args['model'] == 'nonlinear':
            self.states = np.zeros((self.n_steps, 11))
        elif self.args['model'] == 'J2 nonlinear':
            self.states = np.zeros((self.n_steps, 11))
        else:
            raise Exception('Invalid model selected.')

        self.steps = 1

    def diffy_q(self, state):

        # Retrieve relative position and velocity
        r, v = state[0:3], state[3:6]

        if self.args['model'] == 'CW':
            if len(state) != 6:
                raise Exception(f'The state vector must have 6 elements! Its length is {len(state)}')

            n = sqrt(self.args['mu']/self.args['R'])

            dx, dy, dz = v[0], v[1], v[2]
            ddx = 2*n*v[1] + 3*r[0]*n**2
            ddy = -2*n*v[0]
            ddz = -r[2]*n**2

            ddstate = np.array([dx, dy, dz, ddx, ddy, ddz])

            return ddstate

        elif self.args['model'] == 'tillerson':
            if len(state) != 7:
                raise Exception(f'The state vector must have 7 elements! Its length is {len(state)}')

            # Update true anomaly
            teta = state[7]

            n = sqrt(self.args['mu']/self.args['R'])

            dteta = (n*(1+self.args['e']*cos(teta))**2) / (1-self.args['e']**2) ** (3/2)
            ddteta = -2*(n**2)*self.args['e']*sin(teta)*((1+self.args['e']*cos(teta))**3) / ((1-self.args['e']**2)**3)

            # Update relative motion
            dx, dy, dz = v[0], v[1], v[2]
            ddx = 2*dteta*v[1] + r[0]*dteta**2 + ddteta*r[1] + 2*r[0]*n**2*((1+self.args['e']*cos(teta)) / (1-self.args['e']**2))**3
            ddy = -2*dteta*v[0] + r[1]*dteta**2 - ddteta*r[0] - r[1]*n**2*((1+self.args['e']*cos(teta)) / (1-self.args['e']**2))**3
            ddz = -r[2]*n**2*((1+self.args['e']*cos(teta)) / (1-self.args['e']**2))**3

            ddstate = np.array([dx, dy, dz, ddx, ddy, ddz, dteta])

            return ddstate

        elif self.args['model'] == 'nonlinear':
            if len(state) != 11:
                raise Exception(f'The state vector must have 11 elements! Its length is {len(state)}')

            p = self.args['a']*(1-self.args['e']**2)
            h = sqrt(self.args['mu']*p)

            # Unpack state variables
            teta, R, i, vr = state[7], state[8], state[9], state[10]

            # Compute frequently used trig functions
            si, ci = sin(i), cos(i)
            st, ct = sin(teta), cos(teta)
            rj = sqrt((R+r[0])**2 + r[1]**2 + r[2]**2)

            # Define constant parameters
            wz = h/(R**2)
            wx = 0
            alfa_x = 0
            alfa_z = -(2*h*vr/(R**3))
            n2 = self.args['mu']/(rj**3)
            n2j = self.args['mu']/(rj**3)
            zetai = 0
            zeta = 0

            # Integrate radial distance
            R_dot = np.array([vr])
            # Integrate radial velocity
            vr_dot = np.array([-self.args['mu']/(R**2) + h**2/(R**3)])
            # Integrate the true anomaly
            dteta = np.array([h/(R**2)])
            # Integrate orbit inclination
            i_dot = np.arra([0])
            # Integrate the angular momentum
            h_dot = np.array([0])

            # Integrate the relative dynamics
            dx, dy, dz = v[0], v[1], v[2]
            ddx = 2*v[1]*wz - r[0]*(n2j-wz**2) + r[1]*alfa_z - r[2]*wz*wz - (zetai-zeta)*si*st - R*(n2j-n2)
            ddy = -2*v[0]*wz + 2*v[2]*wx - r[0]*alfa_z - r[1]*(n2j-wz**2-wx**2) + r[2]*alfa_x - (zetai-zeta)*si*ct
            ddz = -2*v[1]*wx - r[0]*wx*wz - r[1]*alfa_x - r[2]*(n2j-wx**2) - (zetai-zeta)*ci

            ddstate = np.array([dx, dy, dz, ddx, ddy, ddz, dteta, R_dot, i_dot, vr_dot, h_dot])

            return ddstate

        elif self.args['model'] == 'J2 nonlinear':
            if len(state) != 11:
                raise Exception(f'The state vector must have 11 elements! Its length is {len(state)}')

            p = self.args['a'] * (1 - self.args['e'] ** 2)
            h = sqrt(self.args['mu'] * p)

            # Unpack state variables
            teta, R, i, vr = state[7], state[8], state[9], state[10]

            # Compute frequently used trig functions
            si, ci = sin(i), cos(i)
            st, ct = sin(teta), cos(teta)
            si2, ci2, st2 = (sin(i))**2, (cos(i))**2, (sin(teta))**2
            rj = sqrt((R + r[0]) ** 2 + r[1] ** 2 + r[2] ** 2)
            rjz = (R+r[0])*si*st + r[1]*si*ct + r[2]*ci

            # J2 perturbation
            J2 = 1.08262668e-3
            R_e = 6378  # km
            k_J2 = 0.5*3*J2*self.args['mu']*R_e**2

            # Define constant parameters
            wz = h / (R ** 2)
            wx = -(k_J2/(h*R**3)) * (sin(2*i)*st)
            alfa_x = -(k_J2*sin(2*i)*ct/(R**5)) + (3*vr*k_J2*sin(2*i)*st/(h*(R**4))) - (8*(k_J2**2)*(si**3)*ci*(st**2)*ct)/((h**2)*(R**6))
            alfa_z = -(2*h*vr/(R**3)) - (k_J2*si2*sin(2*teta)/(R**5))
            n2 = self.args['mu']/(R**3) + k_J2/(R**5) - ((5*k_J2/(R**5))*(si2*st2))
            n2j = self.args['mu']/(rj**3) + k_J2/(rj**5) - 5*k_J2*rjz**2/(rj**7)
            zetai = (2 * k_J2 * rjz / (rj ** 5))
            zeta = (2 * k_J2 / (R ** 4)) * (si * st)

            # Integrate radial distance
            R_dot = np.array([vr])
            # Integrate radial velocity
            vr_dot = np.array([-self.args['mu']/(R**2) + h**2/(R**3) - (k_J2/(R**4))*(1 - 3*si2*st2)])
            # Integrate the true anomaly
            dteta = np.array([h/(R**2) + (2*k_J2/(h*R**3))*(ci2*st2)])
            # Integrate orbit inclination
            i_dot = np.array([-(k_J2/(2*h*R**3))*(sin(2*i)*sin(2*teta))])
            # Integrate the angular momentum
            h_dot = np.array([-(k_J2/(R**3))*(si2*sin(2*teta))])

            # Integrate the relative dynamics
            dx, dy, dz = v[0], v[1], v[2]
            ddx = 2 * dy * wz - r[0] * (n2j - wz ** 2) + r[1] * alfa_z - r[2] * wx * wz - (
                    zetai - zeta) * si * st - R * (n2j - n2)
            ddy = -2 * dx * wz + 2 * dz * wx - r[0] * alfa_z - r[1] * (
                    n2j - wz ** 2 - wx ** 2) + r[2] * alfa_x - (zetai - zeta) * si * ct
            ddz = -2 * dy * wx - r[0] * wx * wz - r[1] * alfa_x - r[2] * (n2j - wx ** 2) - (zetai - zeta) * ci

            ddstate = np.array([dx, dy, dz, ddx, ddy, ddz, dteta, R_dot, i_dot, vr_dot, h_dot])

            return ddstate

        else:
            raise Exception('Invalid model selected.')

    def propagate_orbit(self, state0):
        print('Propagating orbit..')

        ode_sol = solve_ivp(
            fun=self.diffy_q,
            t_span=(0, self.tspan),
            y0=state0,
            method=self.args['propagator'],
            atol=self.args['atol'],
            rtol=self.args['rtol'],
            dense_output=self.args['dense_output'])

        self.states = ode_sol.y.T
        self.ets = ode_sol.t
        self.n_steps = self.states.shape[0]

        return self.ets, self.states
