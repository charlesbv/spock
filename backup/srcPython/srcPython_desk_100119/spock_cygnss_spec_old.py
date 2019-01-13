# This script predicts the positions (latitude and longitude) and the gain of the specular points of the CYGNSS constellation. The inputs are the start and end dates of the simulation, which are entered as argument of this script. The results are put in the SpOCK run directory run.cygnss
# Assumptions:
## - start_date must have the format "YYYY-MM-DD"
## - the user run SpOCK with mpicc. The path to mpirun is set up in path_mpirun 
## - start_date must be older than the most recent TLEs (CYGNSS and GPS)
# Notes
## - to calculate the atmospheric drag, we use SWPC data (observations and predictions if run in the future) of F10.7 and Ap. However, if the end date is more than 43 days in the future from today, SWPC predictions are not valid anymore so we use a constant F10.7 and Ap during the entire propagation.
## - a log file is generated in ../run.cygnss/output/log/log_spock_spec_start_date_end_date.txt

# PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT
path_mpirun = '/data/cygnss/tools/spock/mpi_installation/bin/mpirun'
# end of PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT

import os
import sys
from spock_main_input import *
from datetime import datetime, timedelta

start_date = sys.argv[1]
end_date = sys.argv[2]

# Open log file 
log_filename = "../run.cygnss/output/log/log_spock_spec_" + start_date + "_" + end_date 

# Download CYGNSS TLEs (the ones generated right 'before' (older) than start_date)
print "Downloading CYGNSS TLEs from space-track.org..."
os.system("python cygnss_tle.py " + start_date + " run.cygnss")
## Check that TLEs were found for this period
tle_cygnss_file = open("../run.cygnss/input/tle/cygnss_" + start_date + '.txt')
if len(tle_cygnss_file.readlines()) < 16: # if not all CYGNSS TLEs have been published yet, then download the latest TLEs 
    print "***! CYGNSS TLEs for this period are not available at space-track.org. SpOCK will download the latest ones and starts the propagation at that time !*** > " + log_filename
    os.system("python cygnss_tle.py " + start_date + " run.cygnss latest_tle")
tle_cygnss_file.close()

# Download GPS TLEs (the ones generated right 'before' (older) than start_date)
print "Downloading GPS TLEs from space-track.org..."
os.system("python gps_tle.py " + start_date + " run.cygnss")
## Check that TLEs were found for this period
tle_gps_file = open("../run.cygnss/input/tle/constellation_gps_tle/gps_" + start_date + '.txt')
if len(tle_gps_file.readlines()) < 1: # if GPS TLEs have not been published yet, then download the latest TLEs 
    print "***! GPS TLEs for this period are not available at space-track.org. SpOCK will download the latest ones and starts the propagation at that time !*** >> " + log_filename
    os.system("python gps_tle.py " + start_date + " run.cygnss latest_tle")
tle_gps_file.close()

# Write the SpOCK main input file 
## Set up parameters
run_dir = 'run.cygnss'; main_input_filename = 'spock_spec_' + start_date + '_' + end_date + '.txt'
# for TIME section
date_start = start_date + 'T00:00:00'
date_end = end_date + 'T00:00:00'
dt = 10
# for SPACECRAFT section
nb_sc = 8
gps_tle = 'gps_' + start_date + '.txt'
mass = 28
# for ORBIT section
tle_filename = 'cygnss_' + start_date + '.txt'
# for FORCES section
gravity_order = 4
forces = 'drag sun_gravity moon_gravity'
if datetime.strptime(end_date + "T00:00:00", "%Y-%m-%dT%H:%M:%S") > ( datetime.today() + timedelta(days = 43) ):
    density_mode = 'static' # if the end date is more than 43 days in the future from today, then use a constant F10.7 and Ap
else:
    density_mode = 'dynamic' # if the end date is less than 43 days in the future from today, then use the F10.7 and Ap from SWPC (this is because SWPC predicts F10.7 and Ap up to 45 days in the future (we take a margin of 2 days)). So in this mode, we use SWPC data (observations and predictions if run in the future) of F10.7 and Ap

# for OUTPUT section
name_output = 'out'
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
    # for ORBIT section
    tle_filename,
    # for FORCES section
    gravity_order,
    forces,
    density_mode,
    # for OUTPUT section
    name_output,
    dt_output
    )

# Run SpOCK 
os.chdir("../run.cygnss/")
## CYGNSS and GPS positions
print "Predicting CYGNSS and GPS positions..."
os.system(path_mpirun + " -np 1 spock " + main_input_filename + " >> " + log_filename)

## Specular points positions -> binary files
print "Predicting specular points positions (output in binary files)..."
os.system(path_mpirun + " -np 8 find_specular_points " + main_input_filename + " -lon=0 -rot=0 -min >> " + log_filename)
## Specular points positions -> txt files
print "Converting binary files into txt files and interpolating CYGNSS and GPS positions every second..."
os.system(path_mpirun + " -np 8 run_storm " + main_input_filename + " 0 >> " +  log_filename)
os.chdir("../srcPython")
