import spiceypy as spice
import numpy as np
import os

import planetary_data as pd


r2d = 180 / np.pi

base_dir = os.path.join(
    os.path.dirname( os.path.realpath( __file__ ) ),
    os.path.join( '..', '..', 'Data', 'spice' )
    )

leapseconds_kernel = os.path.join( base_dir, 'lsk/naif0012.tls' )
de432              = os.path.join( base_dir, 'spk/de432s.bsp'   )
pck00010           = os.path.join( base_dir, 'pck/pck00010.tpc' )


# Function to retrieve ID, name and time coverage of all objects in space:
def get_objects(filename, display=False):
    objects = spice.spkobj(filename)
    ids, names, tcs_sec, tcs_cal = [], [], [], []
    n = 0
    if display:
        print('\nObjects in %s:' % filename)

    for o in objects:
        # ID:
        ids.append(o)

        # time coverage in seconds since J2000:
        tc_sec = spice.wnfetd(spice.spkcov(filename, ids[n]), n)

        # convert time coverage to human readible:
        tc_cal = [spice.timout(f, "YYYY MM DD HR:MN:SC.### (TDB) ::TDB") for f in tc_sec]

        # append time coverages to output lists:
        tcs_sec.append(tc_sec)
        tcs_cal.append(tc_cal)

        # get name of body:
        try:
            # add associated name to list:
            names.append(id2body(o))
        except:
            # called if body name does not exixt:
            names.append('Unknown name')

        # print out to console:
        if display:
            print('id: %i\t\tname: %s\t\ttc: %s --> %s' % (ids[-1], names[-1], tc_cal[0], tc_cal[1]))

    return ids, names, tcs_sec, tcs_cal


# Returns body name given the ID:
def id2body(id_):
    return spice.bodc2n(id_)


# Create time array for given time coverages:
def tc2array(tcs, steps):
    arr = np.zeros((steps, 1))
    arr[:, 0] = np.linspace(tcs[0], tcs[1], steps)
    return arr


# Get ephemeris data from a given time array
def get_ephemeris(target, times, frame, observer):
    return np.array(spice.spkezr(target, times, frame, 'NONE', observer)[0])


def calc_ephemeris( target, ets, frame, observer ):
    '''
    Convenience wrapper for spkezr and spkgeo
    '''

    if type( target ) == str:
        return np.array( spice.spkezr( target, ets, frame, 'NONE', observer )[ 0 ] )

    else:
        n_states = len( ets )
        states   = np.zeros( ( n_states, 6 ) )
        for n in range( n_states ):
            states[ n ] = spice.spkgeo( target, ets[ n ], frame, observer )[ 0 ]
        return states

def write_bsp( ets, states, args={} ):
    '''
    Write or append to a BSP / SPK kernel from a NumPy array
    '''
    _args = {
        'bsp_fn'   : 'traj.bsp',
        'spice_id' : -999,
        'center'   : 399,
        'frame'    : 'J2000',
        'degree'   : 5,
        'verbose'  : True,
        'new'      : True,
        'comments' : '',
    }
    for key in args.keys():
        _args[ key ] = args[ key ]

    if _args[ 'new' ]:
        handle = spice.spkopn( _args[ 'bsp_fn' ],
            'SPK_file', len( _args[ 'comments' ] ) )
        action = 'Wrote'
    else:
        handle = spice.spkopa( _args[ 'bsp_fn' ] )
        action = 'Updated'

    spice.spkw09( handle, _args[ 'spice_id' ], _args[ 'center' ],
        _args[ 'frame'  ], ets[ 0 ], ets[ -1 ], '0',
        _args[ 'degree' ], len( ets ),
        states.tolist(), ets.tolist() )

    spice.spkcls( handle )

    if _args[ 'verbose' ]:
        print( f'{action} { _args[ "bsp_fn" ] }.' )


# Compute keplerian elements from the state vector:
def rv2kep(state, et=0, mu=pd.earth['mu'], deg=True):
    rp, e, i, raan, aop, ma, t0, mu, ta, a, T = spice.oscltx(state, et, mu)

    if deg:
        i *= r2d
        ta *= r2d
        raan *= r2d
        aop *= r2d

    return [a, e, i, raan, aop, ta]

def inert2ecef(rs, tspan, frame='J200', filenames=None):
    steps = rs.shape[0]
    Cs = np.zeros((steps, 3, 3))
    rs_ecef = np.zeros(rs.shape)

    for step in range(steps):
        Cs[step, :, :] = spice.pxform(frame, 'ITRF93', tspan[step])
        rs_ecef[step, :] = np.dot(Cs[step, :, :], rs[step, :])

    if filenames is not None:
        pass

    return rs_ecef, Cs

def inert2latlong(rs, tspan, frame='J200', filenames=None):
    steps = rs.shape[0]
    latlongs = np.zeros((steps, 3))
    rs_ecef, Cs = inert2ecef(rs, tspan, frame, filenames)

    for step in range(rs.shape[0]):
        r_norm, lat, lon = spice.reclat(rs_ecef[step, :])
        latlongs[step, :] = [lat*r2d, lon*r2d, r_norm]

    return latlongs, rs_ecef, Cs
