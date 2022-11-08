import matplotlib.pyplot as plt
from math import sqrt
from OrbitalMechanics.OrbitalMotionSimulator import planetary_data as pd
from OrbitalMechanics.OrbitalMotionSimulator import Tools
from OrbitalMechanics.OrbitalMotionSimulator.OrbitPropagator import OrbitPropagator

plt.style.use('dark_background')

if __name__ == '__main__':
    cb = pd.earth

    # Initial orbit parameters
    r_mag = cb['radius'] + 400
    v_mag = sqrt(cb['mu']/r_mag)
    # Initial conditions
    r0 = [r_mag, 0, 0]
    v0 = [0, v_mag, 0]

    # Initial orbit parameters
    r_mag1 = cb['radius'] + 1000
    v_mag1 = sqrt(cb['mu'] / r_mag) * 1.2
    # Initial conditions
    r1 = [r_mag1, 0, 0]
    v1 = [0, v_mag1, 0.4]

    # 1 day tspan
    tspan = 24*3600.0

    # 100 seconds dt
    dt = 100.0

    op0 = OrbitPropagator(r0, v0, tspan, dt)
    op1 = OrbitPropagator(r1, v1, tspan, dt)
    op0.propagate_orbit()
    op1.propagate_orbit()

    Tools.plot_n_orbits([op0.rs, op1.rs], labels=['0', '1'], show_plot=True)
