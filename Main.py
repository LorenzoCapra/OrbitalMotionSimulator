from Tools import Tools, planetary_data as pd
from OrbitPropagator import OrbitPropagator, perturbations
# import matplotlib.pyplot as plt
# plt.style.use('dark_background')

if __name__ == '__main__':
    cb = pd.earth
    perts = perturbations()
    # perts['J2'] = True
    perts['aero'] = True
    perts['cd'] = 2.2
    perts['A'] = 0.1*0.1

    # ISS
    c0 = [cb['radius']+414.0, 0.0006189, 51.64, 0, 234.0, 105.6]
    # GEO
    c1 = [cb['radius']+35800.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    # Vespa
    c2 = [cb['radius']+700, 0.009156, 98, 230, 250.0, 18.0]

    # 2 day tspan
    tspan = 24*3600*2

    # 100 seconds dt
    dt = 100.0

    op0 = OrbitPropagator(c0, tspan, dt, kep=True, perts=perts, mass=4e5)
    # op1 = OrbitPropagator(c1, tspan, dt, kep=True, perts=perts, mass=1000)
    # op2 = OrbitPropagator(c2, tspan, dt, kep=True, perts=perts)

    op0.calculate_kep()
    # op1.calculate_kep()
    # op2.calculate_kep()

    Tools.plot_n_orbits([op0.rs], labels=['ISS', 'GEO', 'Random'], show_plot=True)

    op0.plot_kep(show_plot=True, figsize=(12,6), show_avg=True)
    # op1.plot_kep(show_plot=True)
    # op2.plot_kep(show_plot=True)

    op0.plot_altitudes(show_plot=True, figsize=(12,6), days=True, show_avg=True)

    op0.plot_apo_peri(show_plot=True, figsize=(12,6))

    op0.calc_lat_lon()
    args = {'n_sat': 1}
    # Tools.plot_groundtracks(op0.lat_longs)
