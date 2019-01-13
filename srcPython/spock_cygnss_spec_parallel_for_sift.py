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

# This script predicts the positions of the CYGNSS satellites. It can also predict the positions (latitude and longitude) and the gain of the specular points of the CYGNSS constellation. The inputs are:
# - start date of the simulation (1st argument)
# - end dates of the simulation (2nd argument)
# - write spec as a 3rd argument if you want to predict the positions of the specular points
# Assumptions:
## - the results are put in the SpOCK run directory run.cygnss
## - start_date must have the format "YYYY-MM-DD" or "YYYY-MM-DDTHH:MM:SS". Same for end_date
## - the user run SpOCK with mpicc. The path to mpirun is set up in path_mpirun 
## - start_date must be older than the most recent TLEs (CYGNSS and GPS)
# Notes
## - if start_date (or end_date) has the format "YYYY-MM-DD" then SpOCK will be run starting (or finishing) at midnight of this day.
## - to calculate the atmospheric drag, we use SWPC data (observations and predictions if run in the future) of F10.7 and Ap. However, if the end date is more than 43 days in the future from today, SWPC predictions are not valid anymore so we use a constant F10.7 and Ap during the entire propagation.
## - a log file is generated in ../run.cygnss/output/log/log_spock_spec_start_date_end_date.txt


import os
import sys
from spock_main_input import *
from datetime import datetime, timedelta
from sys import platform
from nb_usable_proc import *
from cygnss_tle_for_sift import *
from gps_tle_for_sift import *
def spock_cygnss_spec_parallel_for_sift(start_date, end_date, spec_or_not):
#     start_date = sys.argv[1]
#     end_date = sys.argv[2]
    if ( 'T' in start_date ) == False:
        start_date = start_date + "T00:00:00"
    if ( 'T' in end_date ) == False:
        end_date = end_date + "T00:00:00"
        
    if spec_or_not == 'spec':
        predict_spec = 1
    else:
        predict_spec = 0

    # Open log file 
    log_filename = "log_spock_spec_start_" + start_date.replace(":", "_") + "_end_" + end_date.replace(":", "_")
    log_file = open(log_filename, "w+")

    # Figure out number of processors available for mpirun
    nb_proc_avail = str(nb_usable_proc(path_mpirun))

    # Download CYGNSS TLEs (the ones generated right 'before' (older) than start_date)
    print >> log_file, "Downloading CYGNSS TLEs from space-track.org..." 
    ## The start_date to download can't have the hour/min/sec in it
    start_date_no_time = start_date[0:10]
    print  >> log_file,"cygnss_tle " + start_date_no_time 
    cygnss_tle_for_sift(  start_date_no_time, '' )
    #os.system("python cygnss_tle.py " + start_date_no_time )
    ## Check that TLEs were found for this period
    tle_cygnss_file = open( 'cygnss_' + start_date_no_time + '.txt')
    if len(tle_cygnss_file.readlines()) < 16: # if not all CYGNSS TLEs have been published yet, then download the latest TLEs 
        print  >> log_file,"***! CYGNSS TLEs for this period are not available at space-track.org. SpOCK will download the latest ones and starts the propagation at that time !***"
        cygnss_tle_for_sift(start_date_no_time, 'latest_tle')
        #os.system("python cygnss_tle.py " + start_date_no_time + " latest_tle")
    tle_cygnss_file.close()

    if predict_spec == 1:
        # Download GPS TLEs (the ones generated right 'before' (older) than start_date)
        print >> log_file, "Downloading GPS TLEs from space-track.org..."
        gps_tle_for_sift(start_date_no_time, '')
        #os.system("python gps_tle.py " + start_date_no_time )
        ## Check that TLEs were found for this period
        tle_gps_file = open("gps_" + start_date_no_time + '.txt')
        if len(tle_gps_file.readlines()) < 1: # if GPS TLEs have not been published yet, then download the latest TLEs 
            print  >> log_file,"***! GPS TLEs for this period are not available at space-track.org. SpOCK will download the latest ones and starts the propagation at that time !***"
            gps_tle_for_sift(start_date_no_time, 'latest_tle')
            #os.system("python gps_tle.py " + start_date_no_time + " latest_tle")
        tle_gps_file.close()

    # Write the SpOCK main input file 
    ## Set up parameters
    main_input_filename = 'spock_spec_start_' + start_date.replace(":", "_") + '_end_' + end_date.replace(":", "_") + '.txt'
    # for TIME section
    date_start = start_date 
    date_end = end_date 
    dt = 10 # !!!!!!!change back to 10
    # for SPACECRAFT section
    nb_sc = 8
    if predict_spec == 1:
        gps_tle = 'gps_' + start_date_no_time + '.txt'
    else:
        gps_tle = '0'
    mass = 28
    geometry_filename = "cygnss_geometry_2016.txt"

    # for ORBIT section
    tle_filename = 'cygnss_' + start_date_no_time + '.txt'
    # for FORCES section
    gravity_order = 4 #!!!!!!!!change back to 4
    forces = 'drag sun_gravity moon_gravity'
    if datetime.strptime(end_date, "%Y-%m-%dT%H:%M:%S") > ( datetime.today() + timedelta(days = 43) ):
        density_mode = 'static' # if the end date is more than 43 days in the future from today, then use a constant F10.7 and Ap
    else:
        density_mode = 'dynamic' # !!!!!!!!!!! should be dynamic # if the end date is less than 43 days in the future from today, then use the F10.7 and Ap from SWPC (this is because SWPC predicts F10.7 and Ap up to 45 days in the future (we take a margin of 2 days)). So in this mode, we use SWPC data (observations and predictions if run in the future) of F10.7 and Ap

    # for OUTPUT section
    name_output = 'spock/out'
    dt_output = 10

    # #for SPICE section put in makeall.sh
    # spice_path = "/Users/cbv/Downloads/cspice/data"

    ## Create main input file
    spock_main_input(
         main_input_filename,
        # for TIME section
        date_start,
        date_end,
        dt,
        # for SPACECRAFT section
        nb_sc,
        gps_tle,
        mass,
        geometry_filename, 
        # for ORBIT section
        tle_filename,
        # for FORCES section
        gravity_order,
        forces,
        density_mode,
        # for OUTPUT section
        name_output,
        dt_output,
        # for ATTITUDE section
        "nadir",
        # for GROUNDS_STATIONS section
        "my_ground_stations.txt",
         # for SPICE section
         spice_path,
         # for DENSITY_MOD section
    1
        )

    # Run SpOCK 
    ## CYGNSS and GPS positions (no GPS position is predicted if the user did not call 'spec' at a 3rd argument of this script)
    if predict_spec == 1:
        print  >> log_file,"Predicting CYGNSS and GPS positions..."
    else:
        print  >> log_file,"Predicting CYGNSS positions..."
    print  >> log_file, path_mpirun + " -np " + nb_proc_avail + " spock " + main_input_filename
    if ( ( 'win' in platform ) & (('darwin' in platform ) == False )):
        os.system(path_mpirun + " -np " + nb_proc_avail + " spock.exe " + main_input_filename) # cant redirect otuput. really you can't.
    else:
        os.system(path_mpirun + " -np " + nb_proc_avail + " spock " + main_input_filename)
        #used to be (11-30-17): os.system(path_mpirun + " -np " + nb_proc_avail + " spock_dev_parallel_kalman_9state " + main_input_filename + " >> " + log_filename)

    if predict_spec == 1:
        ## Specular points positions -> binary files
        print >> log_file, "Predicting specular points positions (output in binary files)..."
        print >> log_file, path_mpirun + " -np " + nb_proc_avail + " spec " + main_input_filename + " -lon=0 -rot=0 -min"
        if ( ( 'win' in platform ) & (('darwin' in platform ) == False )):
            os.system(path_mpirun + " -np " + nb_proc_avail + " spec.exe " + main_input_filename + " -lon=0 -rot=0 -min")
        else:
            os.system(path_mpirun + " -np " + nb_proc_avail + " spec " + main_input_filename + " -lon=0 -rot=0 -min")
            #used to be (11-30-17): os.system(path_mpirun + " -np " + nb_proc_avail + " spec_dev_parallel_kalman_save " + main_input_filename + " -lon=0 -rot=0 -min >> " + log_filename)
        ## Specular points positions -> txt files
        print >> log_file, "Converting binary files into txt files and interpolating CYGNSS and GPS positions every second..."
        print >> log_file, path_mpirun + " -np " + nb_proc_avail + " storm " + main_input_filename + " 0"
        if ( ( 'win' in platform ) & (('darwin' in platform ) == False )):
            os.system(path_mpirun + " -np " + nb_proc_avail + " storm.exe " + main_input_filename + " 0")
        else:
            os.system(path_mpirun + " -np " + nb_proc_avail + " storm " + main_input_filename + " 0")
            #used to be (11-30-17): os.system(path_mpirun + " -np " + nb_proc_avail + " storm_dev_parallel_kalman_save " + main_input_filename + " 0 >> " +  log_filename)
    log_file.close()
    return
