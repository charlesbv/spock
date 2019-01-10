# This function converts ECEF coordinates to ECI (J2000) coordinates using the NASA Spice library. It uses the same algorithm as in SpOCK to convert ECEF to ECI.
# INPUTS:
# ECEF position r_ecef (numpy array of 3 elements)
# ECEF velocity v_ecef (numpy array of 3 elements)
# time format YYYY-MM-DDTHH:MM:SS.microseconds or YYYY/MM/DD HH:MM:SS.microseconds (.microseconds is optional) 
# OUTPUTS:
# ECI J2000 position r_eci
# ECI J2000 velocity v_eci
# ASSUMPTIONS:
# - the spice kenerls need to be loaded only once. So first call of ecef2eci should be with load_spice = 1. Then next calls should be with load_spice = 0
# - need to install spiceypy (installation works with pip)
# More info on Spice functions at http://spiceypy.readthedocs.io/en/master/documentation.html#module-spiceypy.spiceypy
# TEST: 
# time = "2017/06/02 00:00:00"
# r_ecef = np.array([5076.4938284307, 4632.9859232012, 520.7687183103])
# v_ecef = np.array([-3.6403679900, 4.4706798164, -4.3262456870])
# ecef2eci(r_ecef, v_ecef, time) should return:
# eci_r = [ 2666.33459874, -6334.89739839,   516.02613502]
# eci_v = [ 5.88632231,  2.12759367, -4.33601995]

import spiceypy as spice
import numpy as np
def ecef2eci(r_ecef, v_ecef, time, load_spice):
    # Load Spice kernels
    if load_spice == 1:
        eop = "/Users/cbv/cspice/data/pck00010.tpc"
        planet_ephem = "/Users/cbv/cspice/data/de432s.bsp"
        earth_binary_pck = "/Users/cbv/cspice/data/earth_000101_180418_180125.bpc"
        leap_sec = "/Users/cbv/cspice/data/naif0012.tls"
        spice.furnsh(eop)
        spice.furnsh(planet_ephem)
        spice.furnsh(earth_binary_pck)
        spice.furnsh(leap_sec)
    
    # Time conversion to seconds past J2000, TDB.
    et = 0
    et = spice.str2et(time)
    # Compute matrix that transforms ECEF to ECI J2000: xform
    xform = spice.sxform("ITRF93", "J2000", et)
    # Apply matrix on ECEF position and velocity: estate
    estate = np.zeros([6])
    estate[0:3] = r_ecef
    estate[3:6] = v_ecef
    jstate = spice.mxvg(xform, estate, 6, 6 )
    # Store the results in eci_r and eci_v
    r_eci = jstate[0:3]
    v_eci = jstate[3:6]

    return r_eci, v_eci

