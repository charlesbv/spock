import numpy as np
import sys

# This scrypt computes common parameters used in astrodynamics
# It uses constants from the WGS84 model
# Units for distance are km

earth_flattening    = 1/298.257223560; # Earth flattening coefficient (no unit)
earth_radius        = 6378.137; # mean equatorial radius (km)
earth_mu    = 398600.4418; # gravitational parameter (km^3/s^2)
earth_j2    = 1.081874e-3; # J2 zonal harmonic coefficient (no unit)
rad2deg = 180./np.pi
deg2rad = 1 / rad2deg
second2year = 3600 * 24 * 365.25

# PARAMETERS TO SET
sma = 6904 # in km 
ecc = 0.001327
inc = 35. # in degrees 

## Convert degree to radian
inc = inc * deg2rad

if 'period' in sys.argv:
    # ORBITAL PERIOD
    # This assumes no perturbation (spherical Earth in particular)
    # Parameters to set are: semi-major axis
    sma = earth_radius + 500 # in km
    period = 2*np.pi * np.sqrt( sma**3 / earth_mu )
    print 'The orbital period is ' + str(period) + ' seconds (= ' + '{0:.2f}'.format(period/60.) + ' minutes).'

if 'precession' in sys.argv:
    # PRECESSION RATE
    # Compute the precession rate taking only J2 into account. 
    # Parameters to set are: semi-major axis, eccentricity, and inclination

    ## Compute angular rate of precession
    omega_dot = -3/2. * earth_radius**2 * earth_j2 * np.sqrt(earth_mu / sma**7.) * 1./(1-ecc**2)**2 * np.cos(inc) # rad/s
    omega_dot_in_degree = omega_dot * rad2deg # degree/s
    nb_rotations_of_ascending_node_in_a_year = np.abs( omega_dot_in_degree * second2year / 360. )
    print 'The precession rate is ' + str(omega_dot) + ' rad/s (= ' + '{0:.2e}'.format(omega_dot_in_degree) +' degree/s). This corresponds to ' + '{0:.2f}'.format(nb_rotations_of_ascending_node_in_a_year) + ' revolutions of the ascending node during a year.'


if 'argument_perigee' in sys.argv: # from http://www.braeunig.us/space/orbmech.htm

    arg_perigee_dot = 1.03237 * 10**14 * sma**(-7./2) * (4-5*(np.sin(inc))**2)/((1-ecc**2)**2) # degree/day
    print arg_perigee_dot
