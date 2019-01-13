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

# This script predicts the positions of a satellite given by its NORAD ID. The inputs are:
# - start date of the simulation (1st argument)
# - end dates of the simulation (2nd argument)
# - norad id (3rd argument)
# Assumptions:
## - start_date must have the format "YYYY-MM-DD" or "YYYY-MM-DDTHH:MM:SS". Same for end_date
## - the user run SpOCK with mpicc. The path to mpirun is set up in path_mpirun 
## - start_date must be older than the most recent TLEs
# Notes
## - if start_date (or end_date) has the format "YYYY-MM-DD" then SpOCK will be run starting (or finishing) at midnight of this day.
## - to calculate the atmospheric drag, we use SWPC data (observations and predictions if run in the future) of F10.7 and Ap. However, if the end date is more than 43 days in the future from today, SWPC predictions are not valid anymore so we use a constant F10.7 and Ap during the entire propagation.
## - a log file is generated in log_spock_spec_start_date_end_date_norad_id.txt

# # #for SPICE section and path to mpirun put in makeall.sh !!!!!!! TO COMMENT 
# path_mpirun = '/usr/local/bin/mpirun'
# spice_path = '/Users/cbv/cspice/data'

import os
import sys
from spock_main_input import *
from datetime import datetime, timedelta

start_date = sys.argv[1]
end_date = sys.argv[2]
if ( 'T' in start_date ) == False:
    start_date = start_date + "T00:00:00"
if ( 'T' in end_date ) == False:
    end_date = end_date + "T00:00:00"

norad_id = str(sys.argv[3])
    
# Open log file 
log_filename = "log_spock_spec_start_" + start_date.replace(":", "_") + "_end_" + end_date.replace(":", "_") + "_" + norad_id

# Download the TLE (the ones generated right 'before' (older) than start_date)
print "Downloading the TLE for " + norad_id + " from space-track.org..."
## The start_date to download can't have the hour/min/sec in it
start_date_no_time = start_date[0:10]
print "download_tle.py " + norad_id + " " + start_date_no_time 
os.system("download_tle.py " + norad_id + " " +start_date_no_time  )# !!!!!!!! remove 'python '
## Check that TLEs were found for this period
tle_file = open( norad_id + '_' + start_date_no_time + '.txt')
if len(tle_file.readlines()) < 2: # if the TLE has not been published yet, then download the latest TLE
    print "***! The TLE for this period is not available at space-track.org. SpOCK will download the latest one and starts the propagation at that time !*** > " + log_filename
    os.system("download_tle.py " + norad_id + " " + start_date_no_time + " latest_tle") # !!!!!!!! remove 'python '
tle_file.close()

# Write the SpOCK main input file 
## Set up parameters
main_input_filename = 'spock_spec_start_' + start_date.replace(":", "_") + '_end_' + end_date.replace(":", "_") + '.txt'
# for TIME section
date_start = start_date 
date_end = end_date 
dt = 10 # !!!!!!!change back to 10
# for SPACECRAFT section
nb_sc = 1
gps_tle = '0'
mass = 28
geometry_filename = "cygnss_geometry_2016.txt"

# for ORBIT section
tle_filename = norad_id + '_' + start_date_no_time + '.txt'
# for FORCES section
gravity_order = 4 #!!!!!!!!change back to 4
forces = 'drag sun_gravity moon_gravity' # !!!!!!!!! change back to 'drag sun_gravity moon_gravity'
if datetime.strptime(end_date, "%Y-%m-%dT%H:%M:%S") > ( datetime.today() + timedelta(days = 43) ):
    density_mode = 'static' # if the end date is more than 43 days in the future from today, then use a constant F10.7 and Ap
else:
    density_mode = 'dynamic' # !!!!!!!!!!! should be dynamic # if the end date is less than 43 days in the future from today, then use the F10.7 and Ap from SWPC (this is because SWPC predicts F10.7 and Ap up to 45 days in the future (we take a margin of 2 days)). So in this mode, we use SWPC data (observations and predictions if run in the future) of F10.7 and Ap

# for OUTPUT section
name_output = 'out'
dt_output = dt


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
    "0",#"my_ground_stations.txt"
     # for SPICE section
     spice_path
    )

# Run SpOCK to predict the positions of norad_id
print "Predicting the positions of " + norad_id +  " ..."
print path_mpirun + " -np 4 spock_dev_parallel " + main_input_filename + " >> " + log_filename
os.system(path_mpirun + " -np 4 spock_dev_parallel " + main_input_filename + " >> " + log_filename)
