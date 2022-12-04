import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import sqrt, cos, sin, tan, atan, acos, pi, log, exp
import datetime
import os
import spiceypy as spice

from Tools import SpiceTools as st, planetary_data as pd

d2r = np.pi/180.0
r2d = 180.0/np.pi
km2AU = 149598073.
sec2day = 1/24/3600

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
    os.path.join('../Data', 'coastlines.csv')
    )

EARTH_SURFACE_IMAGE = os.path.join(
    os.path.dirname( os.path.realpath( __file__ ) ),
    os.path.join('../Data', 'earth_surface.png')
    )

SURFACE_BODY_MAP = {
    'earth'  : EARTH_SURFACE_IMAGE,
}

CITY_COLORS = [
    'w', 'deeppink', 'chartreuse', 'magenta', 'springgreen', 'peachpuff',
    'white', 'lightpink', 'royalblue', 'lime', 'aqua' ] * 100

WORLD_CITIES_FILE = os.path.join(
    os.path.dirname( os.path.realpath( __file__ ) ),
    os.path.join('../Data', 'world_cities.csv')
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
    cs = ['c', 'b', 'r', 'k', 'g', 'm', 'y', 'w', 'r-.']

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
                  AU=False, cmap='Blues', k=1):
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Plot trajectories
    j = 0
    max_val = 0
    for r in rs:
        if AU:
            r /= km2AU
        ax.plot(r[:, 0], r[:, 1], r[:, 2], label=labels[j])
        ax.scatter3D(r[0, 0], r[0, 1], r[0, 2])
        j += 1

        max__ = np.max(r)
        if max__ > max_val:
            max_val = max__

    r_plot = cb['radius']
    if AU:
        r_plot /= km2AU

    # Plot Central Body:
    _u, _v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    _x = r_plot * np.cos(_u) * np.sin(_v) * k
    _y = r_plot * np.sin(_u) * np.sin(_v) * k
    _z = r_plot * np.cos(_v) * k
    ax.plot_surface(_x, _y, _z, cmap=cmap)

    # Plot the x,y,z axis:
    x, y, z = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    u, v, w = [[2 * r_plot, 0, 0], [0, 2 * r_plot, 0], [0, 0, 2 * r_plot]]

    ax.quiver(x, y, z, u, v, w, color='k')

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


# Solve the Lambert Problem:
def lamberts_universal_variables(r0, r1, deltat, tm=1, mu=pd.earth['mu'], tol=1e-6,
                                 max_steps=200, psi=0, psi_u=4*pi**2, psi_l=4*pi):
    # Calculate square root of the mu parameter:
    sqrt_mu = sqrt(mu)

    # Calculate norms of position vectors:
    r0_norm = norm(r0)
    r1_norm = norm(r1)

    # Calculate gamma parameter:
    gamma = np.dot(r0, r1)/r0_norm/r1_norm

    # Calculate beta parameter:
    # beta = tm*sqrt(1-gamma**2)

    # Calculate A parameter:
    A = tm*sqrt(r0_norm*r1_norm*(1+gamma))

    # if A=0, solution can't be computed:
    if A == 0:
        return np.array([0, 0, 0]), np.array([0, 0, 0])

    # Initial values of c2 and c3 parameters:
    c2 = 0.5
    c3 = 1/6.0

    # Counter and Solved variables:
    step = 0
    solved = False

    # While tolerance not met and not at max step:
    for n in range(max_steps):
        # Calculate B parameter:
        B = r0_norm + r1_norm + A*(psi*c3-1)/sqrt(c2)

        # if A and B parameter out of range:
        if A>0.0 and B<0.0:
            # increase lower psi value:
            psi_l += np.pi

            # recalculate B parameter:
            B *= -1

        # Calculate Universal Variable cubed:
        chi3 = sqrt(B/c2)**3

        # Compute deltat variable:
        deltat_ = (chi3*c3 + A*sqrt(B))/sqrt_mu

        # if difference between deltat variables is within the tolerance:
        if abs(deltat - deltat_)<tol:
            # set Solved variable to True:
            solved = True

            # break out of the for loop:
            break

        # if it's not:
        if deltat_ <= deltat:
            # adjust lower psi value:
            psi_l = psi

        # else: deltat_ > deltat
        else:
            # adjust upper psi value:
            psi_u = psi

        # Update psi, c2, c3 values:
        psi = (psi_u + psi_l)/2.0
        c2 = C2(psi)
        c3 = C3(psi)
        print(psi)

    # Check if maximum number of steps is reached:
    if not solved:
        # algorithm did not converge on a psi value:
        print('Lamberts UV variables did not converge')
        return np.array([0, 0, 0]), np.array([0, 0, 0])

    # Calculate coefficients:
    f = 1 - B/r0_norm
    g = A*sqrt(B/mu)
    gdot = 1 - B/r1_norm

    # Calculate velocity vector:
    v0 = (r1 - f*r0)/g
    v1 = (gdot*r1 - r0)/g

    return v0, v1



# Stump function 1:
def C2(psi):
    return (1-cos(sqrt(psi)))/psi


# Stump function 2:
def C3(psi):
    return (sqrt(psi) - sin(sqrt(psi))) / (psi*sqrt(psi))


def groundtracks(coords, labels=None, city_names=None, cs=['w', 'C3', 'b', 'g', 'C1'],
                 surface_image=False, coastlines=True, title='Groundtracks',
                 show_plot=False, save_plot=False, filename='groundtracks.png', dpi=300):
    plt.figure(figsize=(12,6))

    if surface_image:
        plt.imshow(
            plt.imread(SURFACE_BODY_MAP['earth']), extent=[-180, 180, -90, 90]
        )

    if coastlines:
        coast_coords = np.genfromtxt(COASTLINES_COORDINATES_FILE,
                                     delimiter=',')

        plt.plot(coast_coords[:, 0], coast_coords[:, 1], 'mo',
                 markersize=0.3)

    for n in range(len(coords)):
        if labels is None:
            label = str(n)
        else:
            label = labels[n]

        plt.plot(coords[n][0,0], coords[n][0,1], cs[n]+'o', label=label)
        plt.plot(coords[n][1:,0], coords[n][1:,1], cs[n]+'o', markersize=1)

    cities = city_dict()
    n = 0

    for city in city_names:
        coords_ = cities[city]
        plt.plot([coords_[1]], [coords_[0]], cs[n%len(cs)]+'o', markersize=2)

        if n % 2 == 0:
            xytext = (0,2)
        else:
            xytext = (0,-8)

        plt.annotate(city, [coords_[1], coords_[0]], textcoords='offset points', xytext=xytext,
                     ha='center', color=cs[n%len(cs)], fontsize='small')

        n += 1

    plt.grid(linestyle='dotted')
    plt.xlim([-180, 180])
    plt.ylim([-90, 90])
    plt.xlabel(r'Longitude (degrees $^\circ$)')
    plt.ylabel(r'Latitude (degrees $^\circ$)')
    plt.title(title)
    plt.legend()

    if save_plot:
        plt.savefig(filename, dpi=dpi)

    if show_plot:
        plt.show()

def plot_groundtracks( coords, args ):
    _args = {
        'figsize'    : ( 18, 9 ),
        'markersize' : 1,
        'labels'     : [ '' ] * len( coords ),
        'city_names' : city_list0,
        'colors'     : [ 'c', 'r', 'b', 'g', 'w', 'y' ],
        'grid'       : True,
        'title'      : 'Groundtracks',
        'show'       : True,
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
            color      = _args[ 'city_colors' ][ n%len(CITY_COLORS) ],
            markersize = _args[ 'city_msize' ] )

        if n % 2 == 0:
            xytext = ( 0, 2 )
        else:
            xytext = ( 0, -8 )

        plt.annotate( city, [ coords[ 1 ], coords[ 0 ] ],
                      textcoords = 'offset points', xytext = xytext,
                      ha = 'center', color = _args[ 'city_colors' ][ n%len(CITY_COLORS) ],
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


def interplanetary_porkchop(config):
    _config = {
        'planet0': 'Earth',
        'planet1': 'MARS_BARYCENTER',
        'departure0': '2020-07-01',
        'departure1': '2020-09-01',
        'arrival0': '2020-11-01',
        'arrival1': '2022-01-24',
        'mu': pd.sun['mu'],
        'step': 1/sec2day,
        'frame': 'ECLIPJ2000',
        'observer': 'SOLAR_SYSTEM_BARYCENTER',
        'cutoff_v': 20.0,
        'c3_levels': None,
        'vinf_levels': None,
        'tof_levels': None,
        'dv_levels': None,
        'dv_cmap': 'RdPu_r',
        'figsize': (13,7),
        'lw': 1.5,
        'title': 'Porkchop Plot',
        'title_dv': 'Dv Porkchop Plot',
        'fontsize': 15,
        'show': True,
        'filename': None,
        'filename_dv': None,
        'dpi': 300,
        'load': False
    }

    for key in config.keys():
        _config[key] = config[key]

    cutoff_v = _config['cutoff_v']
    cutoff_c3 = cutoff_v**2

    # Arrays of departure and arrival time
    et_departures = np.arange(
        spice.utc2et(_config['departure0']),
        spice.utc2et(_config['departure1']) + _config['step'],
        _config['step']
    )
    et_arrivals = np.arange(
        spice.utc2et(_config['arrival0']),
        spice.utc2et(_config['arrival1']) + _config['step'],
        _config['step']
    )

    # Number of days in each array and total combinations
    ds = len(et_departures)
    as_ = len(et_arrivals)
    total = ds * as_

    print(f'Departure days: {ds}')
    print(f'Arrival days: {as_}')
    print(f'Total combinations: {total}')

    # Create empty array for C3, v infinity and tof
    C3_shorts = np.zeros((as_, ds))
    C3_longs = np.zeros((as_, ds))
    v_inf_shorts = np.zeros((as_, ds))
    v_inf_longs = np.zeros((as_, ds))
    tofs = np.zeros((as_, ds))

    # Create arrays for indexing meshgrid
    x = np.arange(ds)
    y = np.arange(as_)

    for na in y:
        for nd in x:
            # state of planet0 at departure
            state_depart = st.get_ephemeris(_config['planet0'], [et_departures[nd]], _config['frame'],
                                            _config['observer'])[0]
            # state of planet1 at arrival
            state_arrive = st.get_ephemeris(_config['planet1'], [et_arrivals[na]], _config['frame'],
                                            _config['observer'])[0]
            # calculate flight time
            tof = et_arrivals[na] - et_departures[nd]

            try:
                v_sc_depart_short, v_sc_arrive_short = lamberts_universal_variables(
                    state_depart[:3], state_arrive[:3], tof, tm=1, mu=_config['mu']
                )
            except:
                v_sc_depart_short = np.array([1000, 1000, 1000])  # above the cutoff value so not considered
                v_sc_arrive_short = np.array([1000, 1000, 1000])

            try:
                v_sc_depart_long, v_sc_arrive_long = lamberts_universal_variables(
                    state_depart[:3], state_arrive[:3], tof, tm=-1, mu=_config['mu']
                )
            except:
                v_sc_depart_long = np.array([1000, 1000, 1000])
                v_sc_arrive_long = np.array([1000, 1000, 1000])

            # Compute C3 values departing
            C3_short = norm(v_sc_depart_short - state_depart[3:]) ** 2
            C3_long = norm(v_sc_depart_long - state_depart[3:]) ** 2

            # check for unreasonable values
            if C3_short > cutoff_c3: C3_short = cutoff_c3
            if C3_long > cutoff_c3: C3_long = cutoff_c3

            # Compute v infinity values at arrival
            v_inf_short = norm(v_sc_arrive_short - state_arrive[3:]) ** 2
            v_inf_long = norm(v_sc_arrive_long - state_arrive[3:]) ** 2

            # check for unreasonable values
            if v_inf_short > cutoff_v: v_inf_short = cutoff_v
            if v_inf_long > cutoff_v: v_inf_long = cutoff_v

            # Add values to corresponding array
            C3_shorts[na, nd] = C3_short
            C3_longs[na, nd] = C3_long
            v_inf_shorts[na, nd] = v_inf_short
            v_inf_longs[na, nd] = v_inf_long
            tofs[na, nd] = tof

        print(f'{na}    ,   {nd}')

    tofs /= (3600*24)

    # Total delta v
    dv_shorts = v_inf_shorts + np.sqrt(C3_shorts)
    dv_longs = v_inf_longs + np.sqrt(C3_longs)

    # Create levels arrays
    if _config['c3_levels'] is None:
        _config['c3_levels'] = np.arange(10, 50, 2)
    if _config['vinf_levels'] is None:
        _config['vinf_levels'] = np.arange(0, 15, 1)
    if _config['tof_levels'] is None:
        _config['tof_levels'] = np.arange(100, 500, 20)
    if _config['dv_levels'] is None:
        _config['dv_levels'] = np.arange(3, 20, 0.5)

    lw = _config['lw']
    fig, ax = plt.subplots(figsize=_config['figsize'])
    c0 = ax.contour(C3_shorts, levels=_config['c3_levels'], colors='m', linewidths=lw)
    c1 = ax.contour(C3_longs, levels=_config['c3_levels'], colors='m', linewidths=lw)
    c2 = ax.contour(v_inf_shorts, levels=_config['vinf_levels'], colors='deepskyblue', linewidths=lw)
    c3 = ax.contour(v_inf_longs, levels=_config['vinf_levels'], colors='deepskyblue', linewidths=lw)
    c4 = ax.contour(tofs, levels=_config['tof_levels'], colors='white', linewidths=lw*0.5)

    plt.clabel(c0, fmt='%i')
    plt.clabel(c1, fmt='%i')
    plt.clabel(c2, fmt='%i')
    plt.clabel(c3, fmt='%i')
    plt.clabel(c4, fmt='%i')

    plt.plot([0], [0], 'm')
    plt.plot([0], [0], 'c')
    plt.plot([0], [0], 'w')

    plt.legend(
        ['C3', 'v_inf', 'Tof (days)'], bbox_to_anchor=(1.005, 1.01), fontsize=10
    )
    ax.set_title(_config['title'], fontsize=_config['fontsize']*1.5)
    ax.set_ylabel(f"Arrival (days past {_config['arrival0']})", fontsize=_config['fontsize'])
    ax.set_xlabel(f"Departure (days past {_config['departure0']})", fontsize=_config['fontsize'])

    if _config['show']:
        plt.show()

    if _config['filename']:
        plt.savefig(_config['filename'], dpi=_config['dpi'])

    plt.close()

    fig, ax = plt.subplots(figsize=_config['figsize'])
    c0 = ax.contour(
        dv_shorts, levels=_config['dv_levels'], cmap=_config['dv_cmap'], linewidths=lw
    )
    c1 = ax.contour(
        dv_longs, levels=_config['dv_levels'], cmap=_config['dv_cmap'], linewidths=lw
    )
    c2 = ax.contour(
        tofs, levels=_config['tof_levels'], colors='c', linewidths=lw*0.5
    )

    plt.clabel(c0, fmt='%.1f')
    plt.clabel(c1, fmt='%.1f')
    plt.clabel(c2, fmt='%i')

    ax.set_title(_config['title_dv'], fontsize=_config['fontsize']*1.5)
    ax.set_ylabel(f"Arrival (days past {_config['arrival0']}", fontsize=_config['fontsize'])
    ax.set_xlabel(f"Departure (days past {_config['departure0']}", fontsize=_config['fontsize'])

    if _config['show']:
        plt.show()

    if _config['filename_dv']:
        plt.savefig(_config['filename_dv'], dpi=_config['dpi'])


def plot_states( ets, states, args ):
    _args = {
        'figsize'      : ( 16, 8 ),
        'colors'       : COLORS[ : ],
        'dist_unit'    : 'km',
        'time_unit'    : 'seconds',
        'lw'           : 2.5,
        'r_hlines'     : [],
        'v_hlines'     : [],
        'hline_lstyles': 'dashed',
        'title'        : 'Trajectories',
        'xlim'         : None,
        'r_ylim'       : None,
        'v_ylim'       : None,
        'legend'       : True,
        'show'         : False,
        'filename'     : False,
        'dpi'          : 300,
    }
    for key in args.keys():
        _args[ key ] = args[ key ]

    fig, ( ax0, ax1 ) = plt.subplots( 2, 1,
        figsize = _args[ 'figsize' ] )

    _args[ 'xlabel' ]     = time_handler[ _args[ 'time_unit' ] ][ 'xlabel' ]
    _args[ 'time_coeff' ] = time_handler[ _args[ 'time_unit' ] ][ 'coeff' ]
    ts     = ets[:] - ets[0]
    ts    /= _args[ 'time_coeff' ]
    rnorms = np.linalg.norm( states[ :, :3 ], axis = 1 )
    vnorms = np.linalg.norm( states[ :, 3: ], axis = 1 )

    if _args[ 'xlim' ] is None:
        _args[ 'xlim' ] = [ 0, ts[ -1 ] ]

    if _args[ 'r_ylim' ] is None:
        _args[ 'r_ylim' ] = [ states[ :, :3 ].min(), rnorms.max() ]

    if _args[ 'v_ylim' ] is None:
        _args[ 'v_ylim' ] = [ states[ :, 3: ].min(), vnorms.max() ]

    ''' Positions '''
    ax0.plot( ts, states[ :, 0 ], 'r', label = r'$r_x$',
        linewidth = _args[ 'lw' ] )
    ax0.plot( ts, states[ :, 1 ], 'g', label = r'$r_y$',
        linewidth = _args[ 'lw' ] )
    ax0.plot( ts, states[ :, 2 ], 'b', label = r'$r_z$',
        linewidth = _args[ 'lw' ] )
    ax0.plot( ts, rnorms        , 'm', label = r'$Norms$',
        linewidth = _args[ 'lw' ] )

    ax0.grid( linestyle = 'dotted' )
    ax0.set_xlim( _args[ 'xlim'   ] )
    ax0.set_ylim( _args[ 'r_ylim' ] )
    ax0.set_ylabel( r'Position $(km)$')

    for hline in _args[ 'r_hlines' ]:
        ax0.hlines( hline[ 'val' ], ts[ 0 ], ts[ -1 ],
            color     = hline[ 'color' ],
            linestyle = _args[ 'hline_lstyles' ] )

    ''' Velocities '''
    ax1.plot( ts, states[ :, 3 ], 'r', label = r'$r_x$',
        linewidth = _args[ 'lw' ] )
    ax1.plot( ts, states[ :, 4 ], 'g', label = r'$r_y$',
        linewidth = _args[ 'lw' ] )
    ax1.plot( ts, states[ :, 5 ], 'b', label = r'$r_z$',
        linewidth = _args[ 'lw' ] )
    ax1.plot( ts, vnorms        , 'm', label = r'$Norms$',
        linewidth = _args[ 'lw' ] )

    ax1.grid( linestyle = 'dotted' )
    ax1.set_xlim( _args[ 'xlim'   ] )
    ax1.set_ylim( _args[ 'v_ylim' ] )
    ax1.set_ylabel( r'Velocity $(\dfrac{km}{s})$' )
    ax1.set_xlabel( _args[ 'xlabel' ] )

    for hline in _args[ 'v_hlines' ]:
        ax1.hlines( hline[ 'val' ], ts[ 0 ], ts[ -1 ],
            color     = hline[ 'color' ],
            linestyle = _args[ 'hline_lstyles' ] )

    plt.suptitle( _args[ 'title' ] )
    plt.tight_layout()

    if _args[ 'legend' ]:
        ax0.legend()
        ax1.legend()

    if _args[ 'filename' ]:
        plt.savefig( _args[ 'filename' ], dpi = _args[ 'dpi' ] )
        print( 'Saved', _args[ 'filename' ] )

    if _args[ 'show' ]:
        plt.show()

    plt.close()

def state2period(state, coes=True, mu=pd.earth['mu']):
    if coes:
        return 2*pi*sqrt(state[0]**3 / mu)
    else:
        state_ = st.rv2kep(state)
        return 2 * pi * sqrt(state_[0] ** 3 / mu)

def vecs2angle( v0, v1, deg = True ):
    '''
    Calculate angle between 2 vectors
    '''
    angle = acos( np.dot( v0, v1 ) / norm( v0 ) / norm( v1 ) )
    if deg:
        angle *= r2d
    return angle

def calc_close_approach( turn_angle, v_inf, mu = pd.sun[ 'mu' ] ):
	'''
	Calculate periapsis distance in flyby trajectory
	'''
	return mu * ( 1 / sin( turn_angle ) - 1 ) / v_inf ** 2

def calc_vinfinity( tof, args ):

	r1_planet1 = spice.spkgps( args[ 'planet1_ID' ],
		args[ 'et0' ] + tof, args[ 'frame' ], args[ 'center_ID' ] )[ 0 ]

	v0_sc_depart, v1_sc_arrive = lamberts_universal_variables(
		args[ 'state0_planet0' ][ :3 ], r1_planet1, tof, mu=args[ 'mu' ], tm=args[ 'tm' ] )

	vinf = norm( v0_sc_depart - args[ 'state0_planet0' ][ 3: ] )
	return args[ 'vinf' ] - vinf

def vinfinity_match( planet0, planet1, v0_sc, et0, tof0, args={} ):
	'''
	Given an incoming v-infinity vector to planet0, calculate the
	outgoing v-infinity vector that will arrive at planet1 after
	time of flight (tof) where the incoming and outgoing v-infinity
	vectors at planet0 have equal magnitude
	'''
	_args = {
		'et0'       : et0,
		'planet1_ID': planet1,
		'frame'     : 'ECLIPJ2000',
		'center_ID' : 0,
		'mu'        : pd.sun[ 'mu' ],
		'tm'        : 1,
		'diff_step' : 1e-3,
		'tol'       : 1e-4
	}
	for key in args.keys():
		_args[ key ] = args[ key ]

	_args[ 'state0_planet0' ] = spice.spkgeo( planet0, et0,
		_args[ 'frame' ], _args[ 'center_ID' ] )[ 0 ]

	_args[ 'vinf' ] = norm( v0_sc - _args[ 'state0_planet0' ][ 3: ] )

	tof, steps = newton_root_single_fd(
		calc_vinfinity, tof0, _args )

	r1_planet1 = spice.spkgps( planet1, et0 + tof,
		_args[ 'frame' ], _args[ 'center_ID' ] )[ 0 ]

	v0_sc_depart, v1_sc_arrive = lamberts_universal_variables(
		_args[ 'state0_planet0' ][ :3 ], r1_planet1, tof, mu=_args[ 'mu' ], tm=_args[ 'tm' ] )

	return tof, v0_sc_depart, v1_sc_arrive

def newton_root_single( f, fp, x0, args={} ):
	'''
	Calculate root of single variable function
	using explicit derivative function
	'''
	_args = {
		'tol'        : 1e-10,
		'max_steps'  : 50
	}
	for key in args.keys():
		_args[ key ] = args[ key ]

	delta_x = f( x0, args ) / fp( x0, args )

	for n in range( _args[ 'max_steps' ] ):
		x0     -= delta_x
		delta_x = f( x0, args ) / fp( x0, args )

		if abs( delta_x ) < _args[ 'tol' ]:
			return x0, n

def newton_root_single_fd( f, x0, args={} ):
	'''
	Calculate root of single variable function using
	finite differences (no explicit derivative function)
	'''
	_args = {
		'tol'        : 1e-10,
		'max_steps'  : 200,
		'diff_method': 'central',
		'diff_step'  : 1e-6
	}
	for key in args.keys():
		_args[ key ] = args[ key ]

	delta_x = f( x0, _args ) /\
				fdiff_cs( f, x0, _args[ 'diff_step' ], _args )

	for n in range( _args[ 'max_steps' ] ):
		x0     -= delta_x
		delta_x = f( x0, _args ) /\
				  fdiff_cs( f, x0, _args[ 'diff_step' ], _args )

		if abs( delta_x ) < _args[ 'tol' ]:
			return x0, n

	raise RuntimeError( 'Newton\'s root solver FD single variable did not converge.' )

def fdiff_cs( f, x, dx, args={} ):
	'''
	Calculate central finite difference
	of single variable, scalar valued function
	'''
	return ( f( x + dx, args ) - f( x - dx, args ) ) / ( 2 * dx )

def plot_velocities( ets, vs, args ):
	_args = {
		'figsize'          : ( 10, 7 ),
		'dist_unit'        : 'km',
		'time_unit'        : 'seconds',
		'hlines'           : [],
		'hline_lstyles'    : 'dotted',
		'lw'               : 2,
		'labelsize'        : 15,
		'legend_fontsize'  : 20,
		'legend_framealpha': 0.3,
		'title'            : 'Trajectories',
		'xlim'             : None,
		'ylim'             : None,
		'legend'           : True,
		'show'             : False,
		'filename'         : False,
		'dpi'              : 300,
	}
	for key in args.keys():
		_args[ key ] = args[ key ]

	fig, ax0 = plt.subplots( 1, 1, figsize = _args[ 'figsize' ] )

	_args[ 'xlabel' ] = time_handler[ _args[ 'time_unit' ] ][ 'xlabel' ]
	time_coeff        = time_handler[ _args[ 'time_unit' ] ][ 'coeff'  ]

	_ets   = ets.copy() - ets[ 0 ]
	_ets  /= time_coeff
	vnorms = np.linalg.norm( vs, axis = 1 )

	if _args[ 'xlim' ] is None:
		_args[ 'xlim' ] = [ 0, _ets[ -1 ] ]

	if _args[ 'ylim' ] is None:
		_args[ 'ylim' ] = [ vs.min(), vnorms.max() ]

	ax0.plot( _ets, vs[ :, 0 ], 'r', label = r'$v_x$',
		linewidth = _args[ 'lw' ] )
	ax0.plot( _ets, vs[ :, 1 ], 'g', label = r'$v_y$',
		linewidth = _args[ 'lw' ] )
	ax0.plot( _ets, vs[ :, 2 ], 'b', label = r'$v_z$',
		linewidth = _args[ 'lw' ]  )
	ax0.plot( _ets, vnorms    , 'm', label = r'$Norms$',
		linewidth = _args[ 'lw' ] )

	ax0.grid( linestyle = 'dotted' )
	ax0.set_xlim( _args[ 'xlim'   ] )
	ax0.set_ylim( _args[ 'ylim' ] )
	ax0.set_xlabel( _args[ 'xlabel' ], size = _args[ 'labelsize' ] )
	ax0.set_ylabel( r'Velocity $(\dfrac{km}{s})$',
		size = _args[ 'labelsize' ] )

	for hline in _args[ 'hlines' ]:
		ax0.hlines( hline[ 'val' ], _ets[ 0 ], _ets[ -1 ],
			color     = hline[ 'color' ],
			linewidth = _args[ 'lw' ],
			linestyle = _args[ 'hline_lstyles' ] )

	plt.suptitle( _args[ 'title' ] )
	plt.tight_layout()

	if _args[ 'legend' ]:
		ax0.legend( fontsize = _args[ 'legend_fontsize' ],
			loc = 'upper right', framealpha = _args[ 'legend_framealpha' ] )

	if _args[ 'filename' ]:
		plt.savefig( _args[ 'filename' ], dpi = _args[ 'dpi' ] )
		print( 'Saved', _args[ 'filename' ] )

	if _args[ 'show' ]:
		plt.show()

	plt.close()
