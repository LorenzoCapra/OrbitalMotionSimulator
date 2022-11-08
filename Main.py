from math import sqrt
from OrbitalMechanics.OrbitalMotionSimulator import planetary_data as pd
from OrbitalMechanics.OrbitalMotionSimulator import Tools
from OrbitalMechanics.OrbitalMotionSimulator.OrbitPropagator import OrbitPropagator, perturbations

# plt.style.use('dark_background')

if __name__ == '__main__':
    cb = pd.earth
    perts = perturbations()
    perts['J2'] = True

    # ISS
    c0 = [cb['radius']+414.0, 0.0006189, 51.64, 0, 234.0, 105.6]
    # GEO
    c1 = [cb['radius']+35800.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    # Random
    c2 = [cb['radius']+3000, 0.3, 20, 0.0, 15.0, 40.0]

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

    # 20 day tspan
    tspan = 24*3600.0*20

    # 100 seconds dt
    dt = 100.0

    op0 = OrbitPropagator(c0, tspan, dt, kep=True, perts=perts)
    op1 = OrbitPropagator(c1, tspan, dt, kep=True, perts=perts)
    op2 = OrbitPropagator(c2, tspan, dt, kep=True, perts=perts)

    Tools.plot_n_orbits([op0.rs, op1.rs, op2.rs], labels=['ISS', 'GEO', 'Random'], show_plot=True)
