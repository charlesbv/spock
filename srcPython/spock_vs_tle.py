# This script compares the position predicted by SpOCK to the one reported in
# the TLEs, between date_start and date_stop. It linear interpolates the
# position of SpOCK between two time steps at the TLE epochs
# No need to previously downlaod the TLEs, this script does it (unless
# download_tle is set to a value different from 1
# Inputs:
# - date_start: start date of the analysis (YYYY-mm-dd)
# - date_stop: stop date of the analysis (YYYY-mm-dd)
# - download_tle: if 1 then download tle between date_start and date_stop
# - run_spock: if 1 then run SpOCK between between date_start and date_stop  
# ASSUMPTIONS:
# - one TLE per day is downloaded (even if there are more than one per day
#   at space-track.org)
# - date_start will be changed to the epoch of the first TLE rounded
#   to the nearest second
#   For example: if date_start = '2018-04-26' but the TLE epoch is
#   '2018-04-26T04:18:39.173902' then date_start is changed to
#   '2018-04-26T04:18:40' 
# - if more than one TLE in a TLE file, then only consider the first TLE of
#   the file

# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
date_start = '2018-04-26' 
date_stop = '2018-04-29'
download_tle = 0
run_spock = 1
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

import sys
sys.path.append("/Users/cbv/work/spock/srcPythons")
import os
from datetime import datetime, timedelta
from convert_tle_date_to_date import *
from spock_main_input import *
from read_input_file import *
from read_output_file import *
from find_in_read_input_order_variables import *

# Get the TLE epochs and the TLE position
date_start_date = datetime.strptime(date_start, "%Y-%m-%d") \
 + timedelta(days = 1) # add one day because of how cygnss_tle.py works
date_stop_date = datetime.strptime(date_stop, "%Y-%m-%d")
nday = (int)((date_stop_date - date_start_date).total_seconds()/3600./24) + 1
tle_epoch = []
tle_filename = []
r_tle = []
for iday in range(nday):
    date_now = date_start_date + timedelta(days = iday)
    date_now_str = datetime.strftime(date_now, "%Y-%m-%d")
    if download_tle == 1:
        os.system("cygnss_tle.py " + date_now_str)
    tle_filename.append( 'cygnss_' + date_now_str + '.txt' )
    tle_file = open(tle_filename[-1])
    read_tle_file = tle_file.readlines()
    tle_epoch_raw = read_tle_file[0].split()[3]
    tle_epoch.append( convert_tle_date_to_date(tle_epoch_raw) )
    ## Run SpOCK for one second after TLE epoch to convert the TLE
    # state elements to ECI position
    spock_tle_filename = 'TLE_' + tle_filename[-1]
    date_start_tle_spock_date = tle_epoch[-1] + timedelta(seconds = 1)
    date_start_tle_spock =  datetime.strftime( date_start_tle_spock_date, \
                          "%Y-%m-%dT%H:%M:%S")
    
    date_stop_tle_spock_date = date_start_tle_spock_date +timedelta(seconds = 1)
    date_stop_tle_spock =  datetime.strftime( date_stop_tle_spock_date, \
                          "%Y-%m-%dT%H:%M:%S")
    # spock_main_input(
    #     spock_tle_filename,
    #     date_start_tle_spock,
    #     date_stop_tle_spock,
    #     60.,
    #     # for SPACECRAFT section
    #     1,
    #     '0',
    #     29.,
    #     'cygnss_geometry_2016_acco09_sp12.txt',
    #     # for ORBIT section
    #     tle_filename[-1],
    #     # for FORCES section
    #     4,
    #     'drag solar_pressure moon_gravity sun_gravity',
    #     'static',
    #     # for OUTPUT section
    #     '~/work/spockOut/scott/out', 
    #     60, 
    #     # for ATTITUDE section
    #     "nadir",
    #     # for GROUNDS_STATIONS section
    #     "0",#"my_ground_stations.txt"
    #     # for SPICE section
    #     '/Users/cbv/cspice/data', 
    #     # for DENSITY_MOD section
    #     1)
    # os.system("mpirun -np 1 spock_grav_read_bin " + spock_tle_filename)
    var_in, var_in_order = read_input_file(spock_tle_filename)
    output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order,\
                            'output_file_path_list')]; 
    output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order,\
                             'output_file_name_list')];
    var_to_read = ["position_tle"]
    isc = 0
    var_out, var_out_order = read_output_file( output_file_path_list[isc] + \
                                               output_file_name_list[isc], var_to_read )
    r_tle_temp = var_out[find_in_read_input_order_variables(var_out_order, 'position_tle')][0]
    r_tle.append(r_tle_temp)
    
    
# Run SpOCK between date_start and date_stop
date_start_raw = date_start
date_stop_raw = date_stop
date_start_date = tle_epoch[0] + timedelta(seconds = 1)
date_start = datetime.strftime(date_start_date, "%Y-%m-%dT%H:%M:%S")
date_stop = date_stop_raw + 'T00:00:00'
main_input_filename = 'spock_' + date_start_raw + '_to_' + date_stop_raw + '.txt'
if run_spock == 1:
    spock_main_input(
        main_input_filename,
        # for TIME section
        date_start,
        date_stop,
        10.,
        # for SPACECRAFT section
        1,
        '0',
        29.,
        'ballistic_coefficient',#'cygnss_geometry_2016_acco09_sp12.txt',
        # for ORBIT section
        tle_filename[0],
        # for FORCES section
        2,
        'drag solar_pressure moon_gravity sun_gravity',
        'dynamic',
        # for OUTPUT section
        '~/work/spockOut/scott/out', 
        10, 
        # for ATTITUDE section
        "nadir",
        # for GROUNDS_STATIONS section
        "0",#"my_ground_stations.txt"
        # for SPICE section
        '/Users/cbv/cspice/data', 
        # for DENSITY_MOD section
        1)
    os.system("mpirun -np 1 spock_dev " + main_input_filename)

# Read the position predicted by SpOCK
var_in, var_in_order = read_input_file( main_input_filename)
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order,\
                            'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order,\
                             'output_file_name_list')];
nb_sc = var_in[find_in_read_input_order_variables(var_in_order, 'nb_sc')];
isc = 0 
var_to_read = ["position"]
var_out, var_out_order = read_output_file( output_file_path_list[isc] + \
                            output_file_name_list[isc], var_to_read )
date_spock = var_out[find_in_read_input_order_variables(var_out_order, \
                                                        'date_datetime')]
r_spock = var_out[find_in_read_input_order_variables(var_out_order, 'position')]
nstep = len(r_spock)

# Interpolate the SpOCK position at the TLE epochs
r_spock_interpo = np.zeros([nday, 3])
iday = 1 # since SpOCK was nitialized with the first tle the position of SpOCK
# at the inizializatino is the same aas the first TLE so ignore the first TLE
istep = 0
error_spock = np.zeros([nday])
while istep < range(nstep-1):
    if ( ( date_spock[istep] <= tle_epoch[iday] ) & \
       ( date_spock[istep + 1] > tle_epoch[iday] ) ):
        #print date_spock[istep], tle_epoch[iday], date_spock[istep + 1]
        dtime = (tle_epoch[iday] - date_spock[istep]).total_seconds()
        interval_time = (date_spock[istep+1] - date_spock[istep]).total_seconds()
        r0 = r_spock[istep, :]
        r1 = r_spock[istep+1, :]
        interval_r = r1-r0
        slope = interval_r/interval_time
        r_spock_interpo[iday] = r0 + slope*dtime
        error_spock[iday] = np.linalg.norm(r_spock_interpo[iday] - \
                                           r_tle[iday])
        iday = iday + 1
        if iday == nday:
            break
    istep = istep + 1





# Do the same with the STK file
stk_filename = '/Users/cbv/Downloads/CYGFM05_41884 J2000 Position Velocity.txt'
stk_file = open(stk_filename)
read_stk_file = stk_file.readlines()
nheader = 7
nsteps_stk = len(read_stk_file) - nheader
r_stk = np.zeros([nsteps_stk, 3])
date_stk = []
for itime in range(nsteps_stk):
    r_stk[itime, 0] = read_stk_file[itime + nheader].split()[4]
    r_stk[itime, 1] = read_stk_file[itime + nheader].split()[5]
    r_stk[itime, 2] = read_stk_file[itime + nheader].split()[6]
    date_stk_raw = ' '.join(read_stk_file[itime + nheader].split()[:4])
    date_stk.append(datetime.strptime(date_stk_raw,\
                            '%d %b %Y %H:%M:%S.%f')) # 29 Apr 2018 00:00:00.000

# Interpolate the STK position at the TLE epochs
r_stk_interpo = np.zeros([nday, 3])
iday = 1 # since STK was nitialized with the first tle the position of Stk
# at the inizializatino is the same aas the first TLE so ignore the first TLE
istep = 0
error_stk = np.zeros([nday])
while istep < range(nstep-1):
    if ( ( date_stk[istep] <= tle_epoch[iday] ) & \
       ( date_stk[istep + 1] > tle_epoch[iday] ) ):
        #print date_stk[istep], tle_epoch[iday], date_stk[istep + 1]
        dtime = (tle_epoch[iday] - date_stk[istep]).total_seconds()
        interval_time = (date_stk[istep+1] - date_stk[istep]).total_seconds()
        r0 = r_stk[istep, :]
        r1 = r_stk[istep+1, :]
        interval_r = r1-r0
        slope = interval_r/interval_time
        r_stk_interpo[iday] = r0 + slope*dtime
        error_stk[iday] = np.linalg.norm(r_stk_interpo[iday] - \
                                           r_tle[iday])
        iday = iday + 1
        if iday == nday:
            break
    istep = istep + 1

