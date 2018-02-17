# This script:
# - reads in a TLE file from spacetrack.org with all TLEs of MinXSS (TLEs are downloaded using download_tle_spacetrack_run_automatically.py)
# - runs SpOCK for N_min_propagation minutes for each TLE. So 20 TLEs in this TLE file, SpOCK will be run 20 times, the initial epoch for each run will be the TLE epcoh, the final epoch will be initial epoch +N minutes
from read_input_file import *
from read_output_file import *
from spock_main_input import *
from orbit_average import *
from convert_tle_date_to_date import *
import matplotlib.ticker as mtick
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
import pickle
import sys
import fileinput
import time
import numpy as np
from matplotlib import pyplot as plt
import os
import subprocess
from datetime import datetime, timedelta
from out_minxss_sma_average import *

comment = 1

# import cygnss_spock_tle_run_automatically; reload(cygnss_spock_tle_run_automatically); from cygnss_spock_tle_run_automatically  import *
mpi_path = 'mpirun-openmpi-gcc49' # path of mpi run for the machine running this script. BE CAREFUL: if you change it here, you also need to change it in spock_cygnss_spec_set_spock_simu_parameters.py
N_min_propagation = 200 # Number of minutes to propagate each past TLEs in cygnss_spock_tle_run_automatically.py  !!!!!!!! go back to 200


######### ALGORITHM #########
earth_mu    = 398600.4418; # gravitational parameter (km^3/s^2)
earth_radius        = 6378.137; # mean equatorial radius (km)

# Download all MinXSS TLEs from the launch (actually, from 2016-07-16) until 2017-05-06 (when it reentered the atmosphere)
date_start = "2016-07-16"
date_end = "2016-07-19"
tle_filename_all_sc_all_time = date_start + "--" + date_end + ".txt"
if comment != 1:
    link = "https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/" + date_start + "--" + date_end + "/NORAD_CAT_ID/41474/format/tle"
    os.system("wget  --post-data='identity=cbv@umich.edu&password=cygnssisawesome' --cookies=on --keep-session-cookies --save-cookies=cookies.txt 'https://www.space-track.org/ajaxauth/login' -olog")
    os.system("wget --limit-rate=100K --keep-session-cookies --load-cookies=cookies.txt " + link + " -O " + tle_filename_all_sc_all_time)

# Read these TLEs to propagate MinXSS for N_min_propagation
tle_file = open(tle_filename_all_sc_all_time, "r")
read_tle_file = tle_file.readlines()


date_tle = []
sc_name = []

n = len(read_tle_file) # number of lines in the TLE file!!!!!!!!!!!go back to len(read_tle_file)
iline = 0
main_input_filename_list  = [] # list of all the main input file used to run SpOCK
while iline < n: # go over each TLE in the file
    print iline, n
    date_tle_per_sc = []
    sc_name.append( read_tle_file[iline+1].split()[1] ) 
    sc_name_temp = sc_name[-1]
    alt_apogee_per_sc = []
    alt_perigee_per_sc = []
    sma_per_sc = []
    radius_apogee_per_sc = []
    radius_perigee_per_sc = []
    main_input_filename_list_per_sc = []

    while ( ( sc_name_temp == sc_name[-1] ) & ( iline < n ) ): # the file includes blocks of TLE for each sc. Here go over a block for a given sc
        ############
        # Get the date of the current TLE, used to initalize SpOCK
        dt = 5 # dt that will be used to run SpOCK !!!!!!!! put back to 5
        dt_output = 60
        date_tle_per_sc.append( read_tle_file[iline].split()[3] ) # TLE epoch -> initial epoch of SpOCK
        date_tle_start_temp = convert_tle_date_to_date(date_tle_per_sc[-1])
        date_tle_start_tle = datetime.strftime(date_tle_start_temp, "%Y-%m-%dT%H:%M:%S.%f")
        date_tle_start_spock = datetime.strftime(date_tle_start_temp + timedelta(seconds = dt), "%Y-%m-%dT%H:%M:%S") # make sur ethe initial epoch starts at least one time step after TLE epoch
        ############
        # # Create a TLE file for SpOCK with only the current TLE
        tle_filename_one_sc_one_time = "tle_" + sc_name[-1] + "_" + date_tle_start_tle.replace(".","_") + ".txt"
        tle_file_one_sc_one_time = open( tle_filename_one_sc_one_time, "w")
        print >> tle_file_one_sc_one_time, read_tle_file[iline].replace("\r", ""), read_tle_file[iline+1].replace("\r", "")
        tle_file_one_sc_one_time.close()
        ############        
        # Create main input file for SpOCK with this TLE
        main_input_filename = sc_name[-1] + "_" + date_tle_start_tle.replace(".","_") + ".txt" # !!!!!!!!!remove _temp
        main_input_filename_list_per_sc.append(main_input_filename)
        date_tle_end_temp = date_tle_start_temp + timedelta(seconds = dt) + timedelta(minutes = N_min_propagation) # + timedelta(seconds = dt) because needed to make sure that the initial epoch started at least one time step after TLE epoch
        date_tle_end_spock = datetime.strftime(date_tle_end_temp, "%Y-%m-%dT%H:%M:%S")
        spock_main_input( # need to be in spokc/srcPython to run this script   
            main_input_filename,
        # for TIME section
        date_tle_start_spock,
        date_tle_end_spock,
        dt,
            # for SPACECRAFT section
            1,
            '0',
            3.5112,
            "MINXSS_geo_double_face.txt",
        # for ORBIT section
        tle_filename_one_sc_one_time,
        # for FORCES section
        20, # !!!!!!!!!!! put back 20
        "drag sun_gravity moon_gravity", # !!!!!!!!!!!!! put back to "drag sun_gravity moon_gravity"
            'static',
        # for OUTPUT section
            "out",
        dt_output, 
    # for ATTITUDE section
    "sun_pointed",
    # for GROUND_STATIONS section
            "0",
                # for SPICE section
            "/Users/cbv/cspice/data"
    )
        ############
        # Run SpOCK with this main input file
        if comment != 1:
            os.system(mpi_path + " -np 1 spock_dev_parallel_kalman_9state " + main_input_filename) # !!!!!!!! uncomment

        ############
        iline = iline + 2
        if iline < n:
            sc_name_temp = read_tle_file[iline+1].split()[1] 
    #raise Exception
        print iline, n        
    date_tle.append( date_tle_per_sc )
    main_input_filename_list.append( main_input_filename_list_per_sc )


# Reads the output of minxss_sma_average.py
date_ini = date_start + "T00:00:00"
date_ini = datetime.strptime(date_ini, "%Y-%m-%dT%H:%M:%S") 
ecc = []
sma_average_historical = []
x_axis_average_historical = []
radius_apogee = []
radius_perigee = []
dsmadday = []
alt_perigee = []
alt_apogee = []
date_tle = []
sc_name = []
nb_sc = len(main_input_filename_list)
for isc in range(nb_sc):
    date_tle_per_sc = []
    sc_name.append( main_input_filename_list[isc][0][0:7] )
    sc_name_temp = sc_name[-1]
    alt_apogee_per_sc = []
    alt_perigee_per_sc = []
    sma_average_historical_per_sc = []
    radius_apogee_per_sc = []
    radius_perigee_per_sc = []
    x_axis_average_historical_per_sc = []
    nb_tle_per_sc = len(main_input_filename_list[isc])
    for itle in range(nb_tle_per_sc): # !!!!!!!!!! go back to range(nb_tle_per_sc)
        ## Reads the output of the N_min_propagation propagations performed in cygnss_spock_tle_run_automatically.py
        main_input_filename = main_input_filename_list[isc][itle]
        var_in, var_in_order = read_input_file( main_input_filename)
        date_start = var_in[find_in_read_input_order_variables(var_in_order, 'date_start')]; 
        date_tle_per_sc.append( date_start )
        output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
        output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
        ### Read SpOCK's outputs of the N_min_propagation propagation
        isc_current = 0 # here only one sc in main input file
        var_to_read = ["altitude", "sma", "radius", "latitude"]
        var_out, var_out_order = read_output_file( output_file_path_list[isc_current] + output_file_name_list[isc_current], var_to_read )
        date_temp = var_out[find_in_read_input_order_variables(var_out_order, 'date')]
        altitude = var_out[find_in_read_input_order_variables(var_out_order, 'altitude')]
        latitude_temp = var_out[find_in_read_input_order_variables(var_out_order, 'latitude')]
        sma_temp = var_out[find_in_read_input_order_variables(var_out_order, 'sma')]
        sma_orbit_averaged, time_averaged, index_time_averaged = orbit_average(sma_temp, latitude_temp, date_temp ) # !!!!!!!! uncomment

        date_average_start_orbit_list = np.array(time_averaged)[:,0]  # take the date at the start of the bin
        date_average_end_orbit_list = np.array(time_averaged)[:,2]  # take the date at the end of the bin
        nb_orbit_for_this_sc_for_this_tle_run = len(time_averaged)
        for iorbit in range(nb_orbit_for_this_sc_for_this_tle_run):
            date_average_start_orbit = date_average_start_orbit_list[iorbit]
            date_average_start_orbit = datetime.strptime( date_average_start_orbit, "%Y/%m/%d %H:%M:%S" )
            nb_seconds_between_start_orbit_and_date_start = ( date_average_start_orbit - date_ini ).total_seconds()
            date_average_end_orbit = date_average_end_orbit_list[iorbit]
            date_average_end_orbit = datetime.strptime( date_average_end_orbit, "%Y/%m/%d %H:%M:%S" )
            nb_seconds_between_end_orbit_and_date_start = ( date_average_end_orbit - date_ini ).total_seconds()

            x_axis_average_historical_per_sc.append( nb_seconds_between_start_orbit_and_date_start )
            x_axis_average_historical_per_sc.append( nb_seconds_between_end_orbit_and_date_start )
            sma_average_historical_per_sc.append( sma_orbit_averaged[iorbit] ) 


        radius_temp = var_out[find_in_read_input_order_variables(var_out_order, 'radius')]
        alt_perigee_per_sc.append( np.min(altitude) )
        alt_apogee_per_sc.append( np.max(altitude) )
        radius_apogee_per_sc.append( np.max( radius_temp ))
        radius_perigee_per_sc.append( np.min( radius_temp ) )
    date_tle.append( date_tle_per_sc ) # note that this is actually the date of the TLE + one time step (see cygnss_spock_tle_run_automatically.py)
    alt_apogee.append( alt_apogee_per_sc )
    alt_perigee.append( alt_perigee_per_sc )
    sma_average_historical.append(sma_average_historical_per_sc)
    x_axis_average_historical.append(x_axis_average_historical_per_sc)
    radius_apogee.append(radius_apogee_per_sc)
    radius_perigee.append(radius_perigee_per_sc)
# end of reads the output of cygnss_spock_tle_run_automatically.py

