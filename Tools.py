import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import sqrt, cos, sin, tan, atan, acos, pi, log, exp
import datetime
import os

import planetary_data as pd

d2r = np.pi/180.0
r2d = 180.0/np.pi
km2AU = 149598073.

time_handler = {
	'seconds': { 'coeff': 1.0,        'xlabel': 'Time (seconds)' },
	'hours'  : { 'coeff': 3600.0,     'xlabel': 'Time (hours)'   },
	'days'   : { 'coeff': 86400.0,    'xlabel': 'Time (days)'    },
	'years'  : { 'coeff': 31536000.0, 'xlabel': 'Time (years)'   }
}

dist_handler = {
	'km'    : 1.0,
	'ER'    : 1 / 6378.0,
	'JR'    : 1 / 71490.0,
	'AU'    : 6.68459e-9,
	r'$\dfrac{km}{s}$': 1.0
}

COLORS = [
	'm', 'deeppink', 'chartreuse', 'w', 'springgreen', 'peachpuff',
	'white', 'lightpink', 'royalblue', 'lime', 'aqua' ] * 100

COASTLINES_COORDINATES_FILE = os.path.join(
	os.path.dirname( os.path.realpath( __file__ ) ),
	os.path.join( '..', '..', 'Data', 'coastlines.csv' )
	)

EARTH_SURFACE_IMAGE = os.path.join(
	os.path.dirname( os.path.realpath( __file__ ) ),
	os.path.join( '..', '..', 'Data', 'earth_surface.png' )
	)

SURFACE_BODY_MAP = {
	'earth'  : EARTH_SURFACE_IMAGE,
}

CITY_COLORS = [
	'w', 'deeppink', 'chartreuse', 'magenta', 'springgreen', 'peachpuff',
	'white', 'lightpink', 'royalblue', 'lime', 'aqua' ] * 100

WORLD_CITIES_FILE = os.path.join(
	os.path.dirname( os.path.realpath( __file__ ) ),
	os.path.join( '..', '..', 'Data', 'world_cities.csv' )
	)

city_list0 = [
	 'Seattle', 'Pasadena',                 # US
	 'New York', 'San Luis Obispo',
	 'Phoenix', 'Cape Canaveral',
	 'Mexico City', 'Villahermosa',         # Mexico
	 'New Delhi', 'Mumbai',                 # India
	 'Tirunelveli', 'Surat', 'Chennai',
	 'Olney', 'Norwich',                    # England
	 'Ponce',                               # Puerto Rico
	 'Berlin',                              # Germany
	 'Lyon',                                # France
	 'Vienna',                              # Austria
	 'Madrid', 'Sevilla', 'Barcelona',      # Spain
	 'Moscow',                              # Russia
	 'Rome', 'Cortemaggiore',               # Italy
	 'Aalborg',                             # Denmark
	 'Sao Paulo',                           # Brazil
	 'Luxembourg City', 'Esch-sur-Alzette', # Luxembourg
	 'Toronto',                             # Canada
	 'Tokyo',                               # Japan
	 'Istanbul',                            # Turkey
	 'Jihlava',                             # Czech Republic
	 'Warsaw',                              # Poland
	 'Zagreb',                              # Croatia
	 'Sydney', 'Adelaide',                  # Australia
	 'Dubai',                               # UAE
	 'Port Louis',                          # Mauritius
	 'Casablanca',                          # Morocco
	 'Khartoum',                            # Sudan
	 'Tunis',                               # Tunisia
	 'Buenos Aires',                        # Argentina
	 'Cape Town',                           # South Africa
	 'Bucharest',                           # Romania
	 'Bogota',                              # Colombia
	 'Quito',                               # Ecuador
	 'Noordwijk',                           # Netherlands
	 'San Jose',                            # Costa Rica
	 'Stockholm',                           # Sweden
	 'Santiago',                            # Chile
	 'Jakarta',                             # Indonesia
	 'Antwerp',                             # Belgium
	 'Geneva',                              # Switzerland
	 'Manila',                              # Phillipines
	 'Porto', 'Ponta Delgada',              # Portugal
	 'Budapest',                            # Hungary
	 'Panama City',                         # Panama
	 'Cairo',                               # Egypt
	 'Seoul',                               # South Korea
	 'Broom Bridge',                        # Ireland
	 'Lima',                                # Peru
	 'Akure'                                # Nigeria
]


def plot_orbits(rs, args, vectors=[]):
    _args = {
        'figsize': (10, 8),
        'labels': [''] * len(rs),
        'colors': COLORS[:],
        'AU': False,
        'traj_lws': 3,
        'dist_unit': 'km',
        'groundtracks': False,
        'cb_radius': 6378.0,
        'cb_SOI': None,
        'cb_SOI_color': 'c',
        'cb_SOI_alpha': 0.7,
        'cb_axes': True,
        'cb_axes_mag': 2,
        'cb_cmap': 'Blues',
        'cb_axes_color': 'w',
        'axes_mag': 0.8,
        'axes_custom': None,
        'title': 'Trajectories',
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

    for r in rs:
        _r = r.copy() * dist_handler[_args['dist_unit']]

        ax.plot(_r[:, 0], _r[:, 1], _r[:, 2],
                color=_args['colors'][n], label=_args['labels'][n],
                zorder=10, linewidth=_args['traj_lws'])
        ax.plot([_r[0, 0]], [_r[0, 1]], [_r[0, 2]], 'o',
                color=_args['colors'][n])

        if _args['groundtracks']:
            rg = _r / np.linalg.norm(r, axis=1).reshape((r.shape[0], 1))
            rg *= _args['cb_radius']

            ax.plot(rg[:, 0], rg[:, 1], rg[:, 2], cs[n], zorder=10)
            ax.plot([rg[0, 0]], [rg[0, 1]], [rg[0, 2]], cs[n] + 'o', zorder=10)

        max_val = max([_r.max(), max_val])
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

    _args['cb_radius'] *= dist_handler[_args['dist_unit']]
    _u, _v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:20j]
    _x = _args['cb_radius'] * np.cos(_u) * np.sin(_v)
    _y = _args['cb_radius'] * np.sin(_u) * np.sin(_v)
    _z = _args['cb_radius'] * np.cos(_v)
    ax.plot_surface(_x, _y, _z, cmap=_args['cb_cmap'], zorder=1)

    if _args['cb_SOI'] is not None:
        _args['cb_SOI'] *= dist_handler[_args['dist_unit']]
        _x *= _args['cb_SOI'] / _args['cb_radius']
        _y *= _args['cb_SOI'] / _args['cb_radius']
        _z *= _args['cb_SOI'] / _args['cb_radius']
        ax.plot_wireframe(_x, _y, _z,
                          color=_args['cb_SOI_color'],
                          alpha=_args['cb_SOI_alpha'])

    if _args['cb_axes']:
        l = _args['cb_radius'] * _args['cb_axes_mag']
        x, y, z = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        u, v, w = [[l, 0, 0], [0, l, 0], [0, 0, l]]
        ax.quiver(x, y, z, u, v, w, color=_args['cb_axes_color'])

    xlabel = 'X (%s)' % _args['dist_unit']
    ylabel = 'Y (%s)' % _args['dist_unit']
    zlabel = 'Z (%s)' % _args['dist_unit']

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

def plot_n_orbits(rs, labels, cb=pd.earth, show_plot=False, save_plot=False, title='Many Orbits',
                  AU=False, k=1):
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Plot trajectories
    j = 0
    for r in rs:
        if AU:
            r /= km2AU
        ax.plot(r[:, 0], r[:, 1], r[:, 2], label=labels[j])
        ax.scatter3D(r[0, 0], r[0, 1], r[0, 2])
        j += 1

    r_plot = cb['radius']
    if AU:
        r_plot /= km2AU

    # Plot Central Body:
    _u, _v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    _x = r_plot * np.cos(_u) * np.sin(_v) * k
    _y = r_plot * np.sin(_u) * np.sin(_v) * k
    _z = r_plot * np.cos(_v) * k
    ax.plot_surface(_x, _y, _z, cmap='Blues')

    # Plot the x,y,z axis:
    x, y, z = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    u, v, w = [[2 * r_plot, 0, 0], [0, 2 * r_plot, 0], [0, 0, 2 * r_plot]]

    ax.quiver(x, y, z, u, v, w, color='k')

    max_val = np.max(np.abs(rs))
    ax.set_xlim([-max_val, max_val])
    ax.set_ylim([-max_val, max_val])
    ax.set_zlim([-max_val, max_val])

    if AU:
        ax.set_xlabel('X (AU)')
        ax.set_ylabel('Y (AU)')
        ax.set_zlabel('Z (AU)')
    else:
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.set_zlabel('Z (km)')

    # ax.set_aspect('equal')

    ax.set_title(title)
    plt.legend()

    if show_plot:
        plt.show()
    if save_plot:
        plt.savefig(title + 'png', dpi=300)


def kep2car(kep, deg=False, mu=pd.earth['mu']):
    if deg:
        a, e, i, om, w, teta = kep
        i *= d2r
        om *= d2r
        w *= d2r
        teta *= d2r
    else:
        a, e, i, om, w, teta= kep

    E = ecc_anomaly([teta, e], 'tae')
    r_norm = a * (1 - e ** 2) / (1 + e * cos(teta))

    # Calculate r and v in the perifocal frame:
    r_perif = r_norm * np.array([cos(teta), sin(teta), 0])
    v_perif = sqrt(mu * a) / r_norm * np.array([-sin(E), cos(E) * sqrt(1 - e ** 2), 0])

    # Rotation matrix from perifocal to ECI:
    perif2eci = np.transpose(eci2perif(om, w, i))

    # Calculate r and v vectors in the ECI frame:
    r = np.dot(perif2eci, r_perif)
    v = np.dot(perif2eci, v_perif)

    return r, v

# Unused -> look at the function from SpiceTools
def rv2kep(state, mu=pd.earth['mu'], deg=False, print_results=False):
    r = state[:3]
    v = state[3:]

    r_norm = norm(r)
    h = np.cross(r, v)
    h_norm = norm(h)

    # Inclination
    i = acos(h[2]/h_norm)

    # Eccentricity vector
    e = ((norm(v)**2-mu/r_norm)*r - np.dot(r, v)*v)/mu
    e_norm = norm(e)

    # Node line
    N = np.cross([0,0,1], h)
    N_norm = norm(N)

    # RAAN
    om = acos(N[0]/N_norm)
    if N[1] < 0:
        om = 2*np.pi - om  # quadrant check

    # Argument of perigee
    w = acos(np.dot(N, e)/N_norm/e_norm)
    if e[2] < 0:
        w = 2*np.pi - w  # quadrant check

    # True anomaly
    teta = acos(np.dot(e, r)/e_norm/r_norm)
    if np.dot(r, v) < 0:
        teta = 2*np.pi - teta  # quadrant check

    # Semi-major axis
    a = r_norm*(1+e_norm*cos(teta))/(1-e_norm**2)

    if print_results:
        print('a', a)
        print('e', e_norm)
        print('i', i*r2d)
        print('om', om*r2d)
        print('w', w*r2d)
        print('teta', teta*r2d)

    if deg:
        return [a, e_norm, i*r2d, om*r2d, w*r2d, teta*r2d]
    else:
        return [a, e_norm, i, om, w, teta]


def ecc_anomaly(arr, method, tol=1e-8):
    if method == 'newton':
        # Newton's method for iteratively finding E
        Me, e = arr
        if Me < np.pi/2:
            E0 = Me + e/2
        else:
            E0 = Me - e

        for n in range(200):  # Arbitrary max number of steps
            ratio = (E0-e*sin(E0)-Me)/(1-e*cos(E0))
            if abs(ratio) < tol:
                if n == 0:
                    return E0
                else:
                    return E1
            else:
                E1 = E0-ratio
                E0 = E1

        # Did not converge
        return False
    elif method == 'tae':
        teta, e = arr
        return 2*atan(sqrt((1-e)/(1+e))*tan(teta/2))
    else:
        print('Invalid method for eccentric anomaly')


def eci2perif(om, w, i):
    row0 = [-sin(om)*cos(i)*sin(w)+cos(om)*cos(w),
            cos(om)*cos(i)*sin(w)+sin(om)*cos(w),
            sin(i)*sin(w)]

    row1 = [-sin(om)*cos(i)*cos(w)-cos(om)*sin(w),
            cos(om)*cos(i)*cos(w)-sin(om)*sin(w),
            sin(i)*cos(w)]

    row2 = [sin(om)*sin(i), -cos(om)*sin(i), cos(i)]

    return np.array([row0, row1, row2])


def tle2kep(tle_filename, mu=pd.earth['mu']):
    # Read TLE file
    with open(tle_filename, 'r') as f:
        lines = f.readlines()

    # Separate into 3 lines
    line0 = lines[0].strip()
    line1 = lines[1].strip().split()
    line2 = lines[2].strip().split()

    # Epoch (year and day)
    epoch = line1[3]
    year, month, day, hour = compute_epoch(epoch)

    # Collect keplerian elements
    i = float(line2[2]) * d2r
    raan = float(line2[3]) * d2r
    e_string = line2[4]
    e = float('0.'+e_string)
    aop = float(line2[5]) * d2r
    Me = float(line2[6]) * d2r
    mean_motion = float(line2[7])
    T = 1/mean_motion*24*3600
    a = (T**2*mu/4/np.pi**2)**(1/3)

    # Compute eccentric anomaly
    E = ecc_anomaly([Me, e], 'newton')

    # Compute true anomaly
    teta = true_anomaly([E, e])

    return a, e, i, raan, aop, teta, [year, month, day, hour]


def compute_epoch(epoch):
    # Epoch year
    year = int('20'+epoch[:2])

    epoch = epoch[2:].split('.')

    # Day of the year
    day_of_year = int(epoch[0]) - 1

    # Decimal hour of the day
    hour = float('0'+epoch[1]) * 24

    # Get year/month/day/hour
    date = datetime.date(year, 1, 1) + datetime.timedelta(day_of_year)

    # Extract month and day
    month = float(date.month)
    day = float(date.day)

    return year, month, day, hour


def true_anomaly(arr):
    E, e = arr
    return 2*atan(sqrt((1+e)/(1-e)) * tan(E/2))


def tle2cart(tle_filename):
    return kep2car(tle2kep(tle_filename))


# calculate atmospheric density from given altitude:
def calc_atmospheric_density(z):
    rhos, zs = find_rho_z(z)
    if rhos[0] == 0:
        return 0

    Hi = -(zs[1]-zs[0])/log(rhos[1]/rhos[0])

    return rhos[0]*exp(-(z-zs[0])/Hi)


# find endpoints of altitude and density surrounding input altitude:
def find_rho_z(z, zs=pd.earth['zs'], rhos=pd.earth['rhos']):
    if not 1 < z < 1000:
        return [[0, 0], [0, 0]]

    for n in range(len(rhos)-1):
        if zs[n] < z < zs[n+1]:
            return [[rhos[n], rhos[n+1]], [zs[n], zs[n+1]]]

    return [[0, 0], [0, 0]]


# Compute escape velocity at a certain point:
def esc_v(r, mu=pd.earth['mu']):
    return sqrt(2*mu/r)


# Hohmann transfer function:
def hohmann_transfer(r0=0, rf=0, coes0=None, coes1=None, propagate=False, altitude=True, cb=pd.earth,
                     write_output=False):
    # Check if coes are passed in
    if coes0:
        # extract r0 and r1 values
        r0 = coes0[0]
        rf = coes1[0]
    # if passing in altitude (not semi-major axis)
    if altitude:
        # add central body radius to r values
        r0 += cb['radius']
        rf += cb['radius']

    # compute semi-major axis of transfer orbit
    a_transfer = (r0 + rf) / 2.

    # compute velocities of circular orbits
    v_circ_init = sqrt(cb['mu']/r0)
    v_circ_fin = sqrt(cb['mu']/rf)

    # compute initial and final transfer orbit velocities
    v0_transfer = sqrt(cb['mu'] * (2/r0 - 1/a_transfer))
    vf_transfer = sqrt(cb['mu'] * (2/rf - 1/a_transfer))

    # compute transfer time
    t_transfer = pi * sqrt(a_transfer**3 / cb['mu'])

    # compute delta V values
    delta_vs = [abs(v0_transfer - v_circ_init), abs(v_circ_fin - vf_transfer)]

    # propagate the orbits
    if propagate:
        # if coes not passed in
        if not coes0:
            # create coes list
            coes0 = [r0, 0, 0, 0, 0, 0]
            coes1 = [rf, 0, 0, 0, 0, 180.0]

        # compute eccentricity of transfer orbit
        e_transfer = 1 - r0/a_transfer

        # coes for transfer orbit
        coes_transfer = [a_transfer, e_transfer, coes1[2], coes0[3], coes0[4], 0]

        # compute periods of initial and final orbit
        T0 = 2*np.pi*(r0**3 / cb['mu']) ** 0.5
        T1 = 2*np.pi*(rf**3 / cb['mu']) ** 0.5

        return T0, T1, t_transfer, coes_transfer, delta_vs

    return delta_vs, t_transfer


# Stump function 1:
def C2(psi):
    return (1-cos(sqrt(psi)))/psi


# Stump function 2:
def C3(psi):
    return (sqrt(psi) - sin(sqrt(psi))) / (psi*sqrt(psi))


def plot_groundtracks( coords, args ):
	_args = {
		'figsize'    : ( 18, 9 ),
		'markersize' : 1,
		'labels'     : [ '' ] * len( coords ),
		'city_names' : cities_lat_long.city_list0,
		'colors'     : [ 'c', 'r', 'b', 'g', 'w', 'y' ],
		'grid'       : True,
		'title'      : 'Groundtracks',
		'show'       : False,
		'filename'   : False,
		'dpi'        : 300,
		'city_colors': CITY_COLORS[ : ],
		'city_msize' : 3,
		'city_fsize' : 10,
		'legend'     : True,
		'surface_image': True,
		'surface_body' : 'earth',
		'plot_coastlines': False
	}
	for key in args.keys():
		_args[ key ] = args[ key ]

	plt.figure( figsize = _args[ 'figsize' ] )

	if _args[ 'surface_image' ]:
		plt.imshow(
			plt.imread( SURFACE_BODY_MAP[ _args[ 'surface_body' ] ] ),
			extent = [ -180, 180, -90, 90 ] )

	if _args[ 'plot_coastlines' ]:
		coast_coords = np.genfromtxt( COASTLINES_COORDINATES_FILE,
			delimiter = ',' )

		plt.plot( coast_coords[ :, 0 ], coast_coords[ :, 1 ], 'mo',
			markersize = 0.3 )

	for n in range( len( coords ) ):
		plt.plot( [ coords[ n ][ 0, 1 ] ], [ coords[ n ][ 0, 2 ] ], 'o',
			color = _args[ 'colors' ][ n ],
			label = _args[ 'labels' ][ n ] )
		plt.plot( coords[ n ][ 1:, 1 ], coords[ n ][ 1:, 2 ], 'o',
			color = _args[ 'colors' ][ n ],
			markersize = _args[ 'markersize' ] )

	# TODO save this as a .json
	cities = city_dict()
	n      = 0

	for city in _args[ 'city_names' ]:
		coords = cities[ city ]
		plt.plot( [ coords[ 1 ] ], [ coords[ 0 ] ], 'o',
			color      = _args[ 'city_colors' ][ n ],
			markersize = _args[ 'city_msize' ] )

		if n % 2 == 0:
			xytext = ( 0, 2 )
		else:
			xytext = ( 0, -8 )

		plt.annotate( city, [ coords[ 1 ], coords[ 0 ] ],
					  textcoords = 'offset points', xytext = xytext,
					  ha = 'center', color = _args[ 'city_colors' ][ n ],
					  fontsize = _args[ 'city_fsize' ]
					)
		n += 1

	plt.xlim( [ -180, 180 ] )
	plt.ylim( [ -90, 90 ] )
	plt.xticks( range( -180, 200, 20 ) )
	plt.yticks( range( -90, 100, 10 ) )
	plt.xlabel( r'Longitude (degrees $^\circ$)' )
	plt.ylabel( r'Latitude (degrees $^\circ$)' )
	plt.tight_layout()

	if _args[ 'legend' ]:
		plt.legend()

	if _args[ 'grid' ]:
		plt.grid( linestyle = 'dotted' )

	if _args[ 'show' ]:
		plt.show()

	if _args[ 'filename' ]:
		plt.savefig( _args[ 'filename' ], dpi = _args[ 'dpi' ] )



def city_dict():
	with open( WORLD_CITIES_FILE, 'r' ) as f:
		lines = f.readlines()

	header = lines[ 0 ]
	cities = {}

	for line in lines[ 1: ]:
		line = line.split(',')

		# create new dictionary for given city
		try:
			# city name and lat/long coordinates
			cities[ line[ 1 ] ] = [ float( line[ 2 ] ), float( line[ 3 ] ) ]

		except:
			pass

	return cities
