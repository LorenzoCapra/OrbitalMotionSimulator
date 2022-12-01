# Dictionaries containing tha data of different celestial bodies:

import numpy as np

# Gravitational constant
G_meters = 6.67430e-11       # m**3 / kg / s**2
G        = G_meters * 10**-9 # km**3/ kg / s**2


sun = {
    'name': 'Sun',
    'SPICE_ID': 10,
    'mass': 1.989e30,
    'mu': 1.32712e11,
    'radius': 695700,
    'G1': 1e9,
    'deorbit_altitude': 2*695510,
    'spice_file': 'Data/spice/spk/de432s.bsp',
    'cmap':'gist_heat'
}

venus = {
		'name'            : 'Venus',
		'spice_name'      : 'VENUS BARYCENTER',
		'SPICE_ID'        : 2,
		'mass'            : 4.867e24,
		'mu'              : 3.2485859200000006E+05,
		'radius'          : 6051.8,
		'sma'             : 108.209e6,   # km
		'SOI'             : 617183.2511, # km
		'deorbit_altitude': 100.0,       # km
		'cmap'            : 'Wistia',
		'body_fixed_frame': 'IAU_VENUS',
		'traj_color'      : 'y'
		}

atm = np.array([[63.096, 2.059e-4], [251.189, 5.909e-11], [1000.0, 3.561e-15]])
earth = {
    'name': 'Earth',
    'spice_name'      : 'EARTH',
    'SPICE_ID'        : 399,
    'mass': 5.972e24,
    'mu': 398600,  # km^3/s^2
    'radius': 6378,
    'J2': 1.082635854e-3,
    'zs': atm[:, 0],  # km
    'rhos': atm[:, 1]*10**8,  # kg/km^3
    'atm_rot_vector': np.array([0, 0, 72.9211e-6]),  # rad/s
    'sma'             : 149.596e6, # km
	'SOI'             : 926006.6608, # km
    'deorbit_altitude': 100,  # km
    'spice_file': 'Data/spice/spk/de432s.bsp',
    'cmap'            : 'Blues',
	'body_fixed_frame': 'IAU_EARTH',
	'traj_color'      : 'b'
}

moon = {
    'name': 'Moon',
    'spice_name'      : 'MOON',
	'SPICE_ID'        : 301,
    'mass': 7.34767309e22,
    'mu': 4.9e3,  # km^3/s^2
    'radius': 1737.1,
    'orbit_T': 29*24*3600 + 12*3600 + 44*60 + 2.8,
    'dist2earth': 384400,
    'spice_file': 'Data/spice/spk/de432s.bsp',
    'sma'             : 149.596e6, # km
	'SOI'             : 926006.6608, # km
	'deorbit_altitude': 100.0, # km
	'cmap'            : 'Blues',
	'body_fixed_frame': 'IAU_EARTH',
	'traj_color'      : 'b'
}

mars = {
		'name'            : 'Mars',
		'spice_name'      : 'MARS BARYCENTER',
		'SPICE_ID'        : 4,
		'mass'            : 6.39e23,
		'mu'              : 4.282837362069909E+04,
		'radius'          : 3397.0,
		'sma'             : 227.923e6, # km
		'SOI'             : 0.578e6,   # km
		'deorbit_altitude': 50.0,      # km
		'cmap'            : 'Reds',
		'body_fixed_frame': 'IAU_MARS',
		'traj_color'      : 'r'
		}

jupiter = {
		'name'            : 'Jupiter',
		'spice_name'      : 'JUPITER BARYCENTER',
		'SPICE_ID'        : 5,
		'mass'            : 1.898e27,
		'mu'              : 1.26686e8,
		'radius'          : 71490.0,   # km
		'sma'             : 778.570e6, # km
		'deorbit_altitude': 1000.0,    # km
		'SOI'             : 48.2e6,    # km
		'body_fixed_frame': 'IAU_JUPITER',
		'traj_color'      : 'C3'
}

io = {
		'name'            : 'Io',
		'spice_name'      : 'Io',
		'SPICE_ID'        : 501,
		'mass'            : 1.898e27,
		'mu'              : 5.959916033410404E+03,
		'radius'          : 1821.6,   # km
		'deorbit_altitude': 10.0,    # km
		'traj_color'      : 'C1'
}

europa = {
		'name'            : 'Europa',
		'spice_name'      : 'Europa',
		'SPICE_ID'        : 502,
		'mu'              : 3.202738774922892E+03,
		'radius'          : 1560.8,   # km
		'deorbit_altitude': 10.0,    # km
		'traj_color'      : 'C2'
}

ganymede = {
		'name'            : 'Ganymede',
		'spice_name'      : 'Ganymede',
		'SPICE_ID'        : 503,
		'mu'              : 9.887834453334144E+03,
		'radius'          : 2631.2,   # km
		'deorbit_altitude': 100.0,    # km
		'traj_color'      : 'C3'
}

callisto = {
		'name'            : 'Callisto',
		'spice_name'      : 'Callisto',
		'SPICE_ID'        : 504,
		'mu'              : 7.179289361397270E+03,
		'radius'          : 2410.3,   # km
		'deorbit_altitude': 10.0,    # km
		'traj_color'      : 'C4'
}

saturn = {
	'name'            : 'Saturn',
	'spice_name'      : 'SATURN BARYCENTER',
	'SPICE_ID'        : 6,
	'mass'            : 568.34e24,
	'radius'          : 58232.0,
	'mu'              : 37.931e6,
	'sma'             : 1433.529e6,
	'deorbit_altitude': 1000.0,
	'SOI'             : 54890347.727,
	'traj_color'      : 'C2'
}

point = {
    'name': 'Point',
    'mass': 0,
    'mu': 0,
    'radius': 0,
}

bodies = [
	venus, earth, moon, mars,
	jupiter, io, europa, ganymede, callisto,
	saturn, sun ]

for body in bodies:
	body[ 'diameter' ] = body[ 'radius' ] * 2
