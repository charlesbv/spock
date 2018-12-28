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
# This script runs 365 days of CYGNSS/GPS/specular points. For each day, the latest GPS TLE is calculated. The CYGNSS positions are initialized in state (ECEF position/velocity) as the last position of the previous day (so day n at midnight is initialized looking at the file of day n-1 at the last time step, which corresponds to day n midnight). 
# The only 2 things you need to do before running this script are:
# - check that the rm statement in this code deletes the correct folder and not anything else... (the correct folder to delete is the folder of the propagation of the previous day)
# - initialize the orbits in the input file called name_input_propagator for the first day

from read_input_file import *
import fileinput
import os
from os import listdir
from os.path import isfile, join
from datetime import datetime, timedelta
import numpy as np
raise Exception
date_start  = '2014-07-01T00:00'

name_input_propagator = '/home/cbv/Propagator/run_cygnss/input/main_input/cygnss_gps_one_year.txt'
os.system("cp /home/cbv/Propagator/run_cygnss/input/main_input/cygnss_gps_one_year_copy.txt /home/cbv/Propagator/run_cygnss/input/main_input/cygnss_gps_one_year.txt")

# MAKE A LIST OF THE EPOCH START OF EACH DAILY RUN
date_start = datetime.strptime(date_start, "%Y-%m-%dT%H:%M")
day_list = [date_start + timedelta(days=x) for x in np.arange(0, 367)]

# FOR EACH GPS, FIND THE TLE THAT IS THE MOST RECENT FOR EACH OF THESE DAYS OF RUN
directory_gps_tle = '/home/cbv/Propagator/run_cygnss/input/tle/constellation_gps_tle/gps_one_year/'
all_gps_files = [f for f in listdir(directory_gps_tle) if isfile(join(directory_gps_tle, f))]
ini_day = 0
for iday in range(len(day_list)-1):
    print 'day ' + str(iday+1) + ' out of ' + str(len(day_list)) +  ' days'
    tle_file_name = '/home/cbv/Propagator/run_cygnss/input/tle/constellation_gps_tle/GPS_TLE_cygnss.txt'
    tle_file = open(tle_file_name, "w")
    for gps_filename_temp in all_gps_files:
        gps_filename = directory_gps_tle + gps_filename_temp
        gps_file = open(gps_filename, "r")
        read_gps_file = gps_file.readlines()
        skip_header = 0
        nb_lines_header = 0
        while (skip_header < 2):
            if read_gps_file[nb_lines_header].split()[0] == '=====================================================================':
                skip_header = skip_header + 1
            nb_lines_header = nb_lines_header + 1
        nb_tle = (len(read_gps_file) - nb_lines_header)/2

        itle = 0
        tle_epoch_read = read_gps_file[nb_lines_header].split()[3]
        tle_epoch =  datetime.strptime(tle_epoch_read.split('.')[0], "%y%j") + timedelta(days = np.float("0."+tle_epoch_read.split('.')[1]))
        if ( ( day_list[iday] - tle_epoch ).days*24*3600 + ( day_list[iday] - tle_epoch ).seconds +  ( day_list[iday] - tle_epoch ).microseconds / 10**6. > 0 ):
            while ( (itle < nb_tle*2) & ( ( day_list[iday] - tle_epoch ).days*24*3600 + ( day_list[iday] - tle_epoch ).seconds +  ( day_list[iday] - tle_epoch ).microseconds / 10**6. > 0 ) ):
                if read_gps_file[nb_lines_header+itle].split()[0] == '1':
                    tle_epoch_read = read_gps_file[nb_lines_header+itle].split()[3]
                    tle_epoch =  datetime.strptime(tle_epoch_read.split('.')[0], "%y%j") + timedelta(days = np.float("0."+tle_epoch_read.split('.')[1]))
                itle = itle + 1
            right_tle = itle - 3
            tle_epoch_read = read_gps_file[nb_lines_header+right_tle].split()[3]
            tle_epoch =  datetime.strptime(tle_epoch_read.split('.')[0], "%y%j") + timedelta(days = np.float("0."+tle_epoch_read.split('.')[1]))

# WRITE THE GPS TLE FILE FOR THE PROPAGATION OF THIS DAY
            print >> tle_file, 'GPS ' + read_gps_file[nb_lines_header+right_tle+1].split()[1]
            print >> tle_file, read_gps_file[nb_lines_header+right_tle].replace("\n", "").replace('\r','')
            print >> tle_file, read_gps_file[nb_lines_header+right_tle+1].replace("\n", "").replace('\r','')

    tle_file.close()


# EDIT THE PROPAGATOR INPUT FILE: EPOCH START/END, NAME OUTPUT FILE, ORBIT INITIALIZATION
    input_variables, order_input_variables = read_input_file(name_input_propagator)
    name_output_file_prop = 'day' +str(iday+1).zfill(3)+'_cygnss'

    line_input_count = -1
    for line in fileinput.input(name_input_propagator, inplace = True):
        line_input_count = line_input_count + 1
        if ( ( line_input_count != 1) & ( line_input_count != 2 )  & ( line_input_count != 30 )):
            print "%s" % (line),
        if ( line_input_count == 1):
            print "%s\n" % datetime.strftime(day_list[iday], "%Y-%m-%dT%H:%M:%S"),
        if ( line_input_count == 2):
            print "%s\n" % datetime.strftime(day_list[iday+1], "%Y-%m-%dT%H:%M:%S"),
        if ( line_input_count == 30):
            print "%s\n" % name_output_file_prop,
# CYGNSS ARE INITIALIZED USING THE ECEF POSITION/VELOCITY OF THE LAST TIME STEP OF THE PREVIOUS DAY
    if iday > ini_day:
# READ FROM THE PREVIOUS DAY THE ECEF POSITION/VELOCITY OF THE LAST TIME STEP
        ini_orbit_sat = []
        for i in range(8):
            sat_dir = input_variables[6][i]
            ecef_sat = open(sat_dir + "ECEF_" +input_variables[7][i], "r")
            read_ecef_sat = ecef_sat.readlines()
            ini_orbit_sat_temp = [] 
            for iecef in range(6):
                ini_orbit_sat_temp.append(read_ecef_sat[-1].split()[2+iecef])
            ini_orbit_sat_one_sat = '(' + ini_orbit_sat_temp[0] +'; '+ ini_orbit_sat_temp[1] +'; '+ ini_orbit_sat_temp[2] +')' +' (' + ini_orbit_sat_temp[3] +'; '+ ini_orbit_sat_temp[4] +'; '+ ini_orbit_sat_temp[5] +')'
            ini_orbit_sat.append(ini_orbit_sat_one_sat)
# INITALIZE THE NEW DAY WITH THIS ECEF POSITION/VELOCITY
        line_input_count = -1
        for line in fileinput.input(name_input_propagator, inplace = True):
            line_input_count = line_input_count + 1
            if ( ( line_input_count != 13) & ( line_input_count != 14) & ( line_input_count != 15 )  & ( line_input_count != 16 )  & ( line_input_count != 17 )  & ( line_input_count != 18 )  & ( line_input_count != 19 )  & ( line_input_count != 20 ) & ( line_input_count != 21 ) ):  
                print "%s" % (line),
            if ( line_input_count == 13):
                print "%s\n" % "state_ecef",
            if ( line_input_count == 14):
                print "%s\n" % ini_orbit_sat[0],
            if ( line_input_count == 15):
                print "%s\n" % ini_orbit_sat[1],
            if ( line_input_count == 16):
                print "%s\n" % ini_orbit_sat[2],
            if ( line_input_count == 17):
                print "%s\n" % ini_orbit_sat[3],
            if ( line_input_count == 18):
                print "%s\n" % ini_orbit_sat[4],
            if ( line_input_count == 19):
                print "%s\n" % ini_orbit_sat[5],
            if ( line_input_count == 20):
                print "%s\n" % ini_orbit_sat[6],
            if ( line_input_count == 21):
                print "%s\n" % ini_orbit_sat[7],

# Re;mve the previous day to free space
    if (iday > ini_day):
        run_dir = ""
        for j in range(len(input_variables[6][i].split('/'))-2): # input_variables is still from the previous day
            if (j > 0):
                run_dir = run_dir + "/" + input_variables[6][i].split('/')[j] 
        print 'Deleted ' + run_dir
        
        os.system("rm -Rf " + run_dir)
    
    os.chdir("../run_cygnss")
    os.system("/usr/local/mpi/bin/mpirun -np 1 run_moat "+name_input_propagator.split('/')[-1])
    os.system("/usr/local/mpi/bin/mpirun -np 8 find_specular_points " + name_input_propagator.split('/')[-1] + " -lon=0 -rot=0 -min")
    os.system("/usr/local/bin/mpirun -np 8 run_storm_temp " + name_input_propagator.split('/')[-1] + " 0")
    os.chdir("../srcPython")

    os.system("mkdir /raid3/Armada/Charles/python/CYGNSS/simu/cygnss_gps_one_year_091516/itrf93/specular_points/day"+str((iday+1)).zfill(3))
    os.system("mkdir /raid3/Armada/Charles/python/CYGNSS/simu/cygnss_gps_one_year_091516/itrf93/cygnss_gps/day"+str((iday+1)).zfill(3))
    input_variables, order_input_variables = read_input_file(name_input_propagator)
    for i in range(8):
        spec_file  = input_variables[6][i] + "specular_" + input_variables[7][i]
        os.system("cp " + spec_file + " /raid3/Armada/Charles/python/CYGNSS/simu/cygnss_gps_one_year_091516/itrf93/specular_points/day"+str((iday+1)).zfill(3))
    run_dir = ""
    for j in range(len(input_variables[6][i].split('/'))-2):
        if (j > 0):
            run_dir = run_dir + "/" + input_variables[6][i].split('/')[j]
    constellation_gps_file  = run_dir + "/CONSTELLATION_GPS_for_run_"+ input_variables[7][0] 
    constellation_cygnss_file  = run_dir + "/CONSTELLATION_CYGNSS_for_run_"+ input_variables[7][0] 
    os.system("cp " + constellation_gps_file + " /raid3/Armada/Charles/python/CYGNSS/simu/cygnss_gps_one_year_091516/itrf93/cygnss_gps/day"+str((iday+1)).zfill(3))
    os.system("cp " + constellation_cygnss_file + " /raid3/Armada/Charles/python/CYGNSS/simu/cygnss_gps_one_year_091516/itrf93/cygnss_gps/day"+str((iday+1)).zfill(3))
