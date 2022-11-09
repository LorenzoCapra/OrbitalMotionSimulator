import spiceypy as spice
import numpy as np
import planetary_data as pd

r2d = 180 / np.pi


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


# Compute keplerian elements from the state vector:
def rv2kep(state, et=0, mu=pd.earth['mu'], deg=True):
    rp, e, i, raan, aop, ma, t0, mu, ta, a, T = spice.oscltx(state, et, mu)

    if deg:
        i *= r2d
        ta *= r2d
        raan *= r2d
        aop *= r2d

    return [a, e, i, raan, aop, ta]
