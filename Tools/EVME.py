import spiceypy as spice
from Tools import planetary_data as pd
from Tools.ITVIM import ITVIM

def calc_EVME_1963():
	'''
	Earth-Venus-Mars-Earth (EVME) 2 year trajectory
	launching in 1966
	Example comes from Richard Battin's book
	called "Astronautical Guidance"
	'''
	spice.furnsh('Data/spice/lsk/naif0012.tls')
	spice.furnsh('Data/spice/spk/de432s.bsp')

	sequence = [
		{
		'planet': 3,
		'time'  : '1966-02-10',
		'tm'    : -1
		},
		{
		'planet'   : 2,
		'planet_mu': pd.venus[ 'mu' ],
		'time'     : '1966-07-07',
		'tm'       : 1,
		'tol'      : 1e-5
		},
		{
		'planet'   : 4,
		'planet_mu': pd.mars[ 'mu' ],
		'time'     : '1967-01-10',
		'tm'       : -1,
		'tol'      : 1e-5
		},
		{
		'planet'   : 3,
		'planet_mu': pd.earth[ 'mu' ],
		'time'     : '1967-12-18',
		'tol'      : 1e-5
		}
	]
	itvim = ITVIM( { 'sequence': sequence } )
	itvim.print_summary()

	return itvim

if __name__ == '__main__':
	calc_EVME_1963()