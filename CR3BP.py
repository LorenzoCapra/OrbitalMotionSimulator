'''
CR3BP class definition
'''

# 3rd party libraries
from scipy.integrate import solve_ivp
import numpy as np
from numpy.linalg import norm

from Tools import Tools


def null_args():
    return {
        'propagator': 'LSODA',
        'atol': 1e-9,
        'rtol': 1e-9,
        'dense_output': False,
        'tspan': 1e3,
    }


CR3BP_SYSTEMS = {
    'earth-moon': {
        'mu': 0.012277471,
        'L1': 0.8362925909457339,
        'L2': 1.1561681659055247,
        'L3': -1.0051155116068917
    },
    'sun-jupiter': {'mu': 0.000953875}
}


class CR3BP:

    def __init__(self, system, state0):

        if isinstance(system, str):
            if system not in CR3BP_SYSTEMS.keys():
                raise RuntimeError('Invalid CR3BP system.')
            self.mu = CR3BP_SYSTEMS[system]['mu']
        else:
            self.mu = system
        self.one_mu = 1.0 - self.mu

        self.state = state0

        self.propagate_orbit()

    def diffy_q(self, et, state):
        rx, ry, rz, vx, vy, vz = state

        r13_vec = [rx + self.mu, ry, rz]
        r23_vec = [rx - 1 + self.mu, ry, rz]
        r13_3 = norm(r13_vec) ** 3
        r23_3 = norm(r23_vec) ** 3
        omega_x = rx - self.one_mu * (rx + self.mu) / r13_3 - \
                  self.mu * (rx - 1 + self.mu) / r23_3
        omega_y = ry - self.one_mu * ry / r13_3 - self.mu * ry / r23_3
        omega_z = -self.one_mu * rz / r13_3 - self.mu * rz / r23_3

        state_dot = np.zeros(6)
        state_dot[:3] = [vx, vy, vz]
        state_dot[3] = 2 * vy + omega_x
        state_dot[4] = -2 * vx + omega_y
        state_dot[5] = omega_z
        return state_dot

    def propagate_orbit(self, args={}):
        print('Propagating orbit..')

        _args = null_args()
        for key in args.keys():
            _args[key] = args[key]

        self.ode_sol = solve_ivp(
            fun=self.diffy_q,
            t_span=(0, _args['tspan']),
            y0=self.state,
            method=_args['propagator'],
            atol=_args['atol'],
            rtol=_args['rtol'],
            dense_output=_args['dense_output'])

        self.states = self.ode_sol.y.T
        self.ets = self.ode_sol.t
        self.n_steps = self.states.shape[0]

        return self.ets, self.states

    def plot_2d(self, args={'show': True}):
        _args = {
            'show': True
        }
        for key in args.keys():
            _args[key] = args[key]

        Tools.plot_cr3bp_2d(self.mu, [self.states[:, :3]], _args)

    def plot_3d(self, args={'show': True}):
        _args = {
            'show': True
        }
        for key in args.keys():
            _args[key] = args[key]

        Tools.plot_cr3bp_3d(self.mu, [self.states[:, :3]], _args)
