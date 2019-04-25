# This script prepares a sat-bop simulation:
# Inputs:
# - start_time (%Y-%m-%dT%H:%M:%S) 
# - end_time
# ASSUMPTIONS:
# - before running this script, you need to manually get the sp3 file online (ftp://cddis.nasa.gov/gnss/products/) and place it in the directory './sp3'. MAKE SURE THIS IS THE ONLY FILE IN THE DIRECTORY.

import sys
import os
sys.path.append('/Users/cbv/work/spock/srcPython')
from cygnss_beacon_write_input_satbop import *
from cygnss_beacon_package import *
from pathlib import Path
from shutil import copyfile, move
from distutils.dir_util import copy_tree


# PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT
start_time = '2019-04-21T00:00:00'
end_time = '2019-04-21T23:59:59'
# end of PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT

current_dir = str(Path().absolute()) + '/'
#os.chdir(current_dir) if script crases afet chanding dir, can go back to initial dir

# Create the package and run structure for these dates
satbop_package_dir, satbop_simu_dir, satbop_fm_dir = cygnss_beacon_package(start_time, end_time)

# For each FM, set up the sat-bop simulation and run it
os.system('gps_tle_beacon.py ' + start_time[0:10])
os.system('cygnss_tle_beacon.py ' + start_time[0:10])
os.system('rm log_gps_tle_' + start_time[0:10] + '.txt')
os.system('rm gps_' + start_time[0:10] + '_beacon_without_prn.txt')
for cygfm in range(1, 9):
    os.chdir(satbop_fm_dir[cygfm-1])
    # set up sat-bop simulation
    cygnss_tle_filename, gps_tle_filename, output_dir = cygnss_beacon_write_input_satbop(start_time, end_time, cygfm)
    if cygfm == 8:        
        move(current_dir + gps_tle_filename, gps_tle_filename)
    else:
        copyfile(current_dir + gps_tle_filename, gps_tle_filename)
    move(current_dir + cygnss_tle_filename, cygnss_tle_filename)
    copy_tree(current_dir + 'simu/satbop_original_adapted', '.') # sat-bop exe and files for the exe
    # run sat-bop simulation
    #os.system('sat-bop.exe')
    # copy sat-bop input files to the package sat-bop directory 
    copyfile(cygnss_tle_filename, current_dir + satbop_package_dir + cygnss_tle_filename)
    copyfile('input-params.txt', current_dir + satbop_package_dir + 'cygnss_' + start_time[0:10] + '_fm0' + str(cygfm)+ '_config.txt')
    if cygfm == 1:
        copyfile(gps_tle_filename, current_dir + satbop_package_dir + gps_tle_filename)
    # copy sat-bop summary.txt tile to the package sat-bop directory
    try:
        copyfile('summary.txt', current_dir + satbop_package_dir + 'cygnss_' + start_time[0:10] + '_fm0' + str(cygfm)+ '_log.txt')
    except IOError:
        print 'The sat-bop output summary.txt file could not be found.\nThe program will stop.'
        sys.exit()
    os.chdir(current_dir)

    
# Determine which pass and whive PRN to select -> with SpOCK
# copy the corresponding 


# sp3_filename = []
# for file in os.listdir("sp3/"): # figure out the name of the sp3 file. There should be only one file
#     if file.endswith(".sp3"):
#         sp3_filename.append(file)
# if len(sp3_filename) != 1:
#     print 'There are more than one sp3 file in the ./sp3 directory.\nThere should only be the one for the current simulation.\nThe program will stop.'; sys.exit()
# sp3_filename = sp3_filename[0]
    
# copyfile(current_dir + 'sp3/' + sp3_filename, sp3_filename)
#        copyfile(sp3_filename, current_dir + satbop_package_dir + sp3_filename)
