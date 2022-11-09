from math import sqrt
import planetary_data as pd
import Tools
from OrbitPropagator import OrbitPropagator, perturbations

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

    # 2 day tspan
    tspan = 24*3600.0*2

    # 100 seconds dt
    dt = 100.0

    op0 = OrbitPropagator(c0, tspan, dt, kep=True, perts=perts)
    op1 = OrbitPropagator(c1, tspan, dt, kep=True, perts=perts)
    op2 = OrbitPropagator(c2, tspan, dt, kep=True, perts=perts)

    op0.calculate_kep()
    op1.calculate_kep()
    op2.calculate_kep()

    Tools.plot_n_orbits([op0.rs, op1.rs, op2.rs], labels=['ISS', 'GEO', 'Random'], show_plot=True)

    op0.plot_kep(show_plot=True)
    op1.plot_kep(show_plot=True)
    op2.plot_kep(show_plot=True)
