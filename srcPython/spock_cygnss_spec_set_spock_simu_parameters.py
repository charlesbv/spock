# This script predicts the positions of the CYGNSS satellites. It can also predict the positions (latitude and longitude) and the gain of the specular points of the CYGNSS constellation. The inputs are:
# - start date of the simulation (1st argument)
# - end dates of the simulation (2nd argument)
# - write spec as a 3rd argument if you want to predict the positions of the specular points
# - dt=x to set the integration time step for the SpOCK's simulation, where x is the number of seconds. This argument is optional (if not set, dt is set to a default value). NO SPACE between 'dt' and '=', and between '=' and 'x'
# - order=x to set the order of the gravity model for the SpOCK's simulation, where x is the order (if not set, the order is set to a default value). NO SPACE between 'dt' and '=', and between '=' and 'x'
# - nb_sc=x to set the number of sc to propagate for the SpOCK's simulation, where x is the number of sc (if not set, the nb_sc is set to a default value). NO SPACE between 'dt' and '=', and between '=' and 'x'
# - geo=x to set the geometry filename to propagate for the SpOCK's simulation, where x is the name of the geometry file (if not set, the geometry filename is set to a default value). NO SPACE between 'geo' and '=', and between '=' and 'x'
# - att=x to set the attitude to propagate for the SpOCK's simulation, where x is the name of the attitude parameter in section ATTITUDE (if not set, the attitude is set to a default value). NO SPACE between 'att' and '=', and between '=' and 'x'
# - out_name=x to set the output name of the SpOCK's simulation, where x is the name of the output directory in section OUTPUT (if not set, out_name is set to a default value). NO SPACE between 'out_name' and '=', and between '=' and 'x'
# - in_name=x to set the main input file name of the SpOCK's simulation, where x is the name of the main input file name (if not set, in_name is set to a default value). NO SPACE between 'in_name' and '=', and between '=' and 'x'
# Assumptions:
## - the results are put in the SpOCK run directory run_dir
## - start_date must have the format "YYYY-MM-DD" or "YYYY-MM-DDTHH:MM:SS". Same for end_date
## - the user run SpOCK with mpicc. The path to mpirun is set up in path_mpirun 
## - start_date must be older than the most recent TLEs (CYGNSS and GPS)
# Notes
## - if start_date (or end_date) has the format "YYYY-MM-DD" then SpOCK will be run starting (or finishing) at midnight of this day.
## - to calculate the atmospheric drag, we use SWPC data (observations and predictions if run in the future) of F10.7 and Ap. However, if the end date is more than 43 days in the future from today, SWPC predictions are not valid anymore so we use a constant F10.7 and Ap during the entire propagation.
## - a log file is generated in ../" + run_dir + "/output/log/log_" + main_input_filename where main_input_filename is the name of the main input file for the SpOCK's simulation

# PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT
path_mpirun = 'mpirun-openmpi-gcc49'
run_dir = 'run.cygnss'
# end of PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT

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

if 'spec' in sys.argv:
    predict_spec = 1
else:
    predict_spec = 0

# Open log file 
if any(item.startswith('in_name=') for item in sys.argv): # if in_name is set as an argument
    main_input_filename = sys.argv[[item.startswith('in_name=') for item in sys.argv].index(True)].replace("in_name=", "")
else:  # if not, value by default   
    main_input_filename = 'spock_spec_start_' + start_date.replace(":", "_") + '_end_' + end_date.replace(":", "_") + '.txt'

log_filename = "../" + run_dir + "/output/log/log_" + main_input_filename

# Download CYGNSS TLEs (the ones generated right 'before' (older) than start_date)
print "Downloading CYGNSS TLEs from space-track.org..."
## The start_date to download can't have the hour/min/sec in it
start_date_no_time = start_date[0:10]
os.system("python cygnss_tle.py " + start_date_no_time + " " + run_dir)
## Check that TLEs were found for this period
tle_cygnss_file = open("../" + run_dir + "/input/tle/cygnss_" + start_date_no_time + '.txt')
if len(tle_cygnss_file.readlines()) < 16: # if not all CYGNSS TLEs have been published yet, then download the latest TLEs 
    print "***! CYGNSS TLEs for this period are not available at space-track.org. SpOCK will download the latest ones and starts the propagation at that time !*** > " + log_filename
    os.system("python cygnss_tle.py " + start_date_no_time + " " + run_dir + " latest_tle")
tle_cygnss_file.close()

if predict_spec == 1:
    # Download GPS TLEs (the ones generated right 'before' (older) than start_date)
    print "Downloading GPS TLEs from space-track.org..."
    os.system("python gps_tle.py " + start_date_no_time + " " + run_dir)
    ## Check that TLEs were found for this period
    tle_gps_file = open("../" + run_dir + "/input/tle/constellation_gps_tle/gps_" + start_date_no_time + '.txt')
    if len(tle_gps_file.readlines()) < 1: # if GPS TLEs have not been published yet, then download the latest TLEs 
        print "***! GPS TLEs for this period are not available at space-track.org. SpOCK will download the latest ones and starts the propagation at that time !*** >> " + log_filename
        os.system("python gps_tle.py " + start_date_no_time + " " + run_dir + " latest_tle")
    tle_gps_file.close()

# Write the SpOCK main input file 
## Set up parameters
#run_dir = 'run.cygnss'
# for TIME section
date_start = start_date 
date_end = end_date 
if any(item.startswith('dt=') for item in sys.argv): # if dt is set as an argument
    dt = (int)(sys.argv[[item.startswith('dt=') for item in sys.argv].index(True)].replace("dt=", ""))
else: # if not, value by default
    dt = 10
# for SPACECRAFT section
if any(item.startswith('nb_sc=') for item in sys.argv): # if nb_sc is set as an argument
    nb_sc = (int)(sys.argv[[item.startswith('nb_sc=') for item in sys.argv].index(True)].replace("nb_sc=", ""))
else:  # if not, value by default   
    nb_sc = 8

if predict_spec == 1:
    gps_tle = 'gps_' + start_date_no_time + '.txt'
else:
    gps_tle = '0'
mass = 28

if any(item.startswith('geo=') for item in sys.argv): # if geo is set as an argument
    geometry_filename = sys.argv[[item.startswith('geo=') for item in sys.argv].index(True)].replace("geo=", "")
else:
    geometry_filename = "cygnss_geometry_2016.txt"
# for ORBIT section
tle_filename = 'cygnss_' + start_date_no_time + '.txt'
# for FORCES section

if any(item.startswith('order=') for item in sys.argv): # if order is set as an argument
    gravity_order = (int)(sys.argv[[item.startswith('order=') for item in sys.argv].index(True)].replace("order=", ""))
else:  # if not, value by default   
    gravity_order = 10

if any(item.startswith('att=') for item in sys.argv): # if attitude is set as an argument
    attitude = sys.argv[[item.startswith('att=') for item in sys.argv].index(True)].replace("att=", "")
else:  # if not, value by default   
    attitude = 'nadir'

if any(item.startswith('out_name=') for item in sys.argv): # if out_name is set as an argument
    name_output = sys.argv[[item.startswith('out_name=') for item in sys.argv].index(True)].replace("out_name=", "")
else:  # if not, value by default   
    name_output = 'out'

forces = 'drag sun_gravity moon_gravity'
if datetime.strptime(end_date, "%Y-%m-%dT%H:%M:%S") > ( datetime.today() + timedelta(days = 43) ):
    density_mode = 'static' # if the end date is more than 43 days in the future from today, then use a constant F10.7 and Ap
else:
    density_mode = 'dynamic' # if the end date is less than 43 days in the future from today, then use the F10.7 and Ap from SWPC (this is because SWPC predicts F10.7 and Ap up to 45 days in the future (we take a margin of 2 days)). So in this mode, we use SWPC data (observations and predictions if run in the future) of F10.7 and Ap

# for OUTPUT section
dt_output = 60
## Create main input file
spock_main_input(
    run_dir, main_input_filename,
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
    attitude,
        # for GROUND_STATIONS section
        "0"

    )

# Run SpOCK 
os.chdir("../" + run_dir + "/")
## CYGNSS and GPS positions (no GPS position is predicted if the user did not call 'spec' at a 3rd argument of this script)
if predict_spec == 1:
    print "Predicting CYGNSS and GPS positions..."
else:
    print "Predicting CYGNSS positions..."
os.system(path_mpirun + " -np 1 spock " + main_input_filename + " >> " + log_filename)

if predict_spec == 1:
    ## Specular points positions -> binary files
    print "Predicting specular points positions (output in binary files)..."
    os.system(path_mpirun + " -np 8 find_specular_points " + main_input_filename + " -lon=0 -rot=0 -min >> " + log_filename)
    ## Specular points positions -> txt files
    print "Converting binary files into txt files and interpolating CYGNSS and GPS positions every second..."
    os.system(path_mpirun + " -np 8 run_storm " + main_input_filename + " 0 >> " +  log_filename)
os.chdir("../srcPython")
