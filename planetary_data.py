# Dictionaries containing tha data of different celestial bodies:

import numpy as np


sun = {
    'name': 'Sun',
    'mass': 1.989e30,
    'mu': 1.32712e11,
    'radius': 695700,
    'G1': 1e9,
    'deorbit_altitude': 2*695510,
    'spice_file': 'Data/spice/spk/de432s.bsp'
}

atm = np.array([[63.096, 2.059e-4], [251.189, 5.909e-11], [1000.0, 3.561e-15]])
earth = {
    'name': 'Earth',
    'mass': 5.972e24,
    'mu': 398600,  # km^3/s^2
    'radius': 6378,
    'J2': 1.082635854e-3,
    'zs': atm[:, 0],  # km
    'rhos': atm[:, 1]*10**8,  # kg/km^3
    'atm_rot_vector': np.array([0, 0, 72.9211e-6]),  # rad/s
    'deorbit_altitude': 100,  # km
    'spice_file': 'Data/spice/spk/de432s.bsp'
}

moon = {
    'name': 'Moon',
    'mass': 7.34767309e22,
    'mu': 4.9e3,  # km^3/s^2
    'radius': 1737.1,
    'orbit_T': 29*24*3600 + 12*3600 + 44*60 + 2.8,
    'dist2earth': 384400,
    'spice_file': 'Data/spice/spk/de432s.bsp'
}
# moon['orbit_w'] = 2*np.pi/moon['orbit_T']

point = {
    'name': 'Point',
    'mass': 0,
    'mu': 0,
    'radius': 0,
}
