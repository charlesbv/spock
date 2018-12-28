# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at

#   http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.
# This script:
# - reads in a TLE file from spacetrack.org with all TLEs of all sc from a start date to an end date (TLEs are downloaded using download_tle_spacetrack_run_automatically.py)
# - for each sc, runs SpOCK for N_min_propagation minutes for each TLE. So if sc1 has 20 TLEs in this file, SpOCK will be run 20 times, the initial epoch for each run will be the TLE epcoh, the final epoch will be initial epoch +N minutes
# This script is tun automatically from cygnss_spock_past_tle_and_future_prediction.py so don't make modifications to it unless you know what you're doing!
from read_input_file import *
from read_output_file import *
from spock_main_input import *
from orbit_average import *
from convert_tle_date_to_date import *
from norad_id_to_cygnss_name import *
import matplotlib.ticker as mtick
import matplotlib.colors as colors
import matplotlib.cm as cmx
from get_prop_dir import *
import matplotlib.gridspec as gridspec
from read_input_file import *
import pickle
import sys
import fileinput
import time
import numpy as np
from matplotlib import pyplot as plt
import os
import subprocess
from get_name_mission import *
from find_in_read_input_order_variables import *
from datetime import datetime, timedelta

# import cygnss_spock_tle_run_automatically; reload(cygnss_spock_tle_run_automatically); from cygnss_spock_tle_run_automatically  import * 
def cygnss_spock_tle_run_automatically(mpi_path, run_dir, N_min_propagation):
    ######### ALGORITHM #########
    earth_mu    = 398600.4418; # gravitational parameter (km^3/s^2)
    earth_radius        = 6378.137; # mean equatorial radius (km)

    # Download all CYGNSS TLEs from the launch (actually, from Dec 20 2016) until today (actually, until the latest TLEs available at space-track.org)
    os.system("python download_tle_spacetrack_run_automatically.py " + run_dir)  # !!!!!!!!! uncomment
    tle_filename_all_sc_all_time = "../" + run_dir + "/input/tle/" + "2016-12-20--" + datetime.strftime( datetime.today(), "%Y-%m-%d") + ".txt" # "./cygnss/2016-12-20--" + datetime.strftime( datetime.today(), "%Y-%m-%d") + ".txt" # name of the file including all TLEs of all CYGNSS from the launch until today

    # Read these TLEs to propagate the CYNGSS sc for N_min_propagation
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
        sc_name.append( norad_id_to_cygnss_name( read_tle_file[iline+1].split()[1] ) )
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
            tle_filename_one_sc_one_time = "CYG" + sc_name[-1] + "_" + date_tle_start_tle.replace(".","_") + ".txt"
            tle_file_one_sc_one_time = open("../" + run_dir + "/input/tle/" + tle_filename_one_sc_one_time, "w")
            print >> tle_file_one_sc_one_time, read_tle_file[iline].replace("\r", ""), read_tle_file[iline+1].replace("\r", "")
            tle_file_one_sc_one_time.close()
            ############        
            # Create main input file for SpOCK with this TLE
            main_input_filename = "CYG" + sc_name[-1] + "_" + date_tle_start_tle.replace(".","_") + ".txt" # !!!!!!!!!remove _temp
            main_input_filename_list_per_sc.append(main_input_filename)
            date_tle_end_temp = date_tle_start_temp + timedelta(seconds = dt) + timedelta(minutes = N_min_propagation) # + timedelta(seconds = dt) because needed to make sure that the initial epoch started at least one time step after TLE epoch
            date_tle_end_spock = datetime.strftime(date_tle_end_temp, "%Y-%m-%dT%H:%M:%S")
            spock_main_input( # need to be in spokc/srcPython to run this script   
                run_dir, main_input_filename,
            # for TIME section
            date_tle_start_spock,
            date_tle_end_spock,
            dt,
                # for SPACECRAFT section
                1,
                '0',
                29,
                "cygnss_geometry_2016.txt",
            # for ORBIT section
            tle_filename_one_sc_one_time,
            # for FORCES section
            10, # !!!!!!!!!!! put back 10
            "drag sun_gravity moon_gravity", # !!!!!!!!!!!!! put back to "drag sun_gravity moon_gravity"
                'static',
            # for OUTPUT section
                "out",
            dt_output, 
        # for ATTITUDE section
        "nadir",
        # for GROUND_STATIONS section
        "0"

        )
            ############
            # Run SpOCK with this main input file
            os.chdir("../" + run_dir)
            os.system(mpi_path + " -np 1 spock " + main_input_filename) # !!!!!!!! uncomment
            os.chdir("../srcPython")

            ############
            iline = iline + 2
            if iline < n:
                sc_name_temp = norad_id_to_cygnss_name( read_tle_file[iline+1].split()[1] )
        #raise Exception
            print iline, n        
        date_tle.append( date_tle_per_sc )
        main_input_filename_list.append( main_input_filename_list_per_sc )

    return main_input_filename_list
