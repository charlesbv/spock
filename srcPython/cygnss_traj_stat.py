# This script was made to answer Scott's email on Dec 4, 2018.
# The script:
# 1- reads the position and velocity (ECEF r/v) of a given FM between date_start_str and date_stop_str
# 2- converts ECEF r/v to ECI r/v
# 3- converts ECI r/v to osculating orbital elements
# 4- for each orbital element, compute distributions
# INPUTS:
# - cygfm: which FM to look at
# - date_start: start date of the simulation
# - date_stop: stop date of the simulation

# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
cygfm = 3
date_start_str = '2018-10-01' # YYYY-MM-DD
date_stop_str = '2018-10-02' # YYYY-MM-DD
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

import sys
sys.path.append('/Users/cbv/work/spock/srcPython')
from datetime import datetime, timedelta
import os
from cygnss_read_netcdf_and_convert_to_eci import *

cygfm_str = str(cygfm)
netcdf_dir = '/Users/cbv/cygnss/netcdf'
if netcdf_dir[-1] != '/':
    netcdf_dir = netcdf_dir + '/'

date_start = datetime.strptime(date_start_str, '%Y-%m-%d')
date_stop = datetime.strptime(date_stop_str, '%Y-%m-%d')
nb_day = (int)((date_stop - date_start).total_seconds()/3600./24)
for iday in range(nb_day):
    print 'day ' + str(iday) + ' out of ' + str(nb_day) + ' day(s).'
    date = date_start + timedelta(days = iday)
    yy = datetime.strftime(date, '%Y')
    doy = datetime.strftime(date, '%j')
    # Read netcdf file to get the ECEF r/v
    netcdf_path = netcdf_dir + yy + '/' + doy + '/'
    filename = [filename for filename in os.listdir(netcdf_path) if filename.startswith('cyg0' + cygfm_str)][0]
    netcdf_filename = netcdf_path + filename
    date_flight_rounded, lon_cyg, lat_cyg, lon_spec, lat_spec, fom, gps,\
        x_cyg, y_cyg, z_cyg, vx_cyg, vy_cyg, vz_cyg,date_flight_rounded_date,\
        r_eci_cyg, v_eci_cyg = cygnss_read_netcdf_and_convert_to_eci(netcdf_filename)
    # Convert ECEF r/v to ECI r/v
    

