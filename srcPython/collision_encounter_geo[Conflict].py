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
# This cretes input files for collision of different encounter geometry. It first propagates backward SpOCK with different orbital elements representing different encouter geometries. It then take the last time step and initialize the collision runs with this last step. It actually creates 10 forward simu for each backward simu: each of these 10 simu is using the option swpc_mod wth a different quantile
# Assumptions:
# - see section #PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT


import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/kalman/spock_development_new_structure_kalman_dev/srcPython")
import matplotlib.gridspec as gridspec
import numpy as np
from struct import *
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
from datetime import datetime, timedelta
from collections import *
import os
from read_input_file import *
from read_output_file import *
from cygnss_read_spock_spec import *
from datetime import datetime, timedelta
from spock_main_input import *
from spock_collision import *


#PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
date_start_str = '2017-12-01T00:00:00' # when to sstart the collision runs (forward (so this is really the start date))
tca = 36. # number of hours after date_start the collision occurs
alt_arr = [[300, 300],[350,350],[400, 400]] # altitude of each sc for each run
inc_arr = [[90,0], [90,30],[90, 45], [90,60]] # inclination of each sc for each run
arg_per_arr = [[0,0]] # argument of perifee  of each sc for each run
raan_arr = [[0,0]] # RAAN  of each sc for each run
true_ano_arr = [[0,0]] # true anomaly of each sc for each run
ecc_arr = [[0,0]] # eccentricity of each sc for each run
f107_arr = [100] # f107 for each run
nb_ens = 32600 # number of ensemble for the collision (collision is over nb_ens*nb_ens)
min_dist_ca = 10000. # distance under which a close approach is flagged. in m
min_dist_coll = 1.8 # distance under which a collision is recorded. in m
path_to_mpirun = "mpirun" # where mpirun is in the system
spice_path_backward = "/Users/cbv/cspice/data" # where spice is in the system to run the backward simu
spice_path_forward = "/home1/cbussy/installation_spock/spock/cspice/data" # where spice is in the system to run the forward simu
send_folder = 'send_quantile' # name of foler that will include the forward simulations to run on Pleiades. if this folder doesn't exist then this scripts creates it
#end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT


if os.path.isdir(send_folder) == False:
    os.system("mkdir " + send_folder)
    "The folder " + send_folder + " has been created."

quantile_array = [1,2,3,4,5,6,7,8,9] # same quantile for ap and f10.7
nb_quantile = len(quantile_array)
quantile_filename_array = []
for iquantile in range(nb_quantile):
    quantile_filename_array.append("quartile_f107_" + str(quantile_array[iquantile]) + "_quartile_ap_" + str(quantile_array[iquantile]) + '.txt')
    
date_start = datetime.strptime(date_start_str, "%Y-%m-%dT%H:%M:%S")
date_end = date_start + timedelta(hours = tca)
nb_alt = len(alt_arr)
nb_inc = len(inc_arr)
nb_arg_per = len(arg_per_arr)
nb_raan = len(raan_arr)
nb_true_ano = len(true_ano_arr)
nb_ecc = len(ecc_arr)
nb_f107 = len(f107_arr)

nb_run = nb_alt*nb_inc*nb_arg_per*nb_raan*nb_true_ano*nb_ecc*nb_f107*nb_quantile

# Create backward and forward simulations  (1. create backward; 2. run it; 3. create forward)
date_end_str = datetime.strftime(date_end, "%Y-%m-%dT%H:%M:%S")
dt_backward = -1. # need one second dt for the backward otherwise when forward then don't get collision (for forard ok to have 10s)
dt_forward = 10.
nb_sc  = 2
mass = 100
geometry_filename = 'one_plate_old.txt'

irun = -1
irun_alt = -1

filename_list_forward_simu = "encounter_geo_runs.txt"
file_list_forward_simu = open(filename_list_forward_simu, "w")
while irun_alt < nb_alt-1:
    irun_alt = irun_alt + 1
    alt = alt_arr[irun_alt]
    irun_inc = -1
    while irun_inc < nb_inc-1:
        irun_inc = irun_inc + 1
        inc = inc_arr[irun_inc]
        irun_arg_per = -1
        while irun_arg_per < nb_arg_per-1:
            irun_arg_per = irun_arg_per + 1
            arg_per = arg_per_arr[irun_arg_per]
            irun_raan = -1
            while irun_raan < nb_raan-1:
                irun_raan = irun_raan + 1
                raan = raan_arr[irun_raan]
                irun_true_ano = -1
                while irun_true_ano < nb_true_ano-1:
                    irun_true_ano = irun_true_ano + 1
                    true_ano = true_ano_arr[irun_true_ano]
                    irun_ecc = -1
                    while irun_ecc < nb_ecc-1:
                        irun_ecc = irun_ecc + 1
                        ecc = ecc_arr[irun_ecc]
                        irun_f107 = -1
                        while irun_f107 < nb_f107-1:
                            irun_f107 = irun_f107 + 1
                            f107 = f107_arr[irun_f107]
                            density_mode_backward = ['static', str(f107), str(f107), '12'] 

                            # Create backward main input file
                            orbit_type = ['oe', str(alt[0])  + ' ' + str(inc[0])  + ' ' + str(arg_per[0])  + ' ' + str(raan[0])  + ' ' + str(true_ano[0])  + ' ' + str(ecc[0]), str(alt[1])  + ' ' + str(inc[1])  + ' ' + str(arg_per[1])  + ' ' + str(raan[1])  + ' ' + str(true_ano[1])  + ' ' + str(ecc[1])]
                            main_input_filename = 'alt' + str(alt[0]) + '-' + str(alt[1]) + '_inc' + str(inc[0])+ '-' + str(inc[1]) + '_arg_per' + str(arg_per[0])+ '-' + str(arg_per[1]) + '_raan' + str(raan[0])+ '-' + str(raan[1]) + '_true_ano' + str(true_ano[0])+ '-' + str(true_ano[1]) + '_ecc' + str(ecc[0])+ '-' + str(ecc[1]) + '_f107' + str(f107) + '_backward.txt'
                            spock_main_input(
                                main_input_filename,
                                # for TIME section
                                date_end_str,
                                date_start_str,
                                dt_backward,
                                # for SPACECRAFT section
                                2,
                                '0',
                                mass,
                                geometry_filename, 
                                # for ORBIT section
                                orbit_type,
                                # for FORCES section
                                2,
                                'drag',
                                density_mode_backward,
                                # for OUTPUT section
                                'out',
                                86400, # the last time step is the only we care about here and is always printed
                                # for ATTITUDE section
                                "nadir",
                                # for GROUNDS_STATIONS section
                                "0",#"my_ground_stations.txt"
                                # for SPICE section
                                spice_path_backward,
                                # for DENSITY_MOD section
                                [1, "\n#OUTPUT_ENSEMBLES", "eci_r, collision"]
                                )
                            # run backward SpOCK simulation
#                             if ( os.path.isfile(main_input_filename) == False ): # don't run the simu if you already ran it
#                                 os.system(path_to_mpirun + " -np 2 spock " + main_input_filename)

                            # read position and velocity at last time step (ie at date_start since the run was backward)
                            #print 'reading results of ' + main_input_filename 
                            ## read main input to get output_file_path_list and output_file_path_path
                            var_in, var_in_order = read_input_file(main_input_filename)
                            output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
                            output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
                            ## read output file to get last ECI r/v
                            var_to_read = ["position", "velocity"]
                            r = np.zeros([2, 3]); v = np.zeros([2, 3])
                            for isc in range(2):
                                var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
                                r[isc,:] = var_out[find_in_read_input_order_variables(var_out_order, 'position')][-1]
                                v[isc, :] = var_out[find_in_read_input_order_variables(var_out_order, 'velocity')][-1]
                                #print 'r for ', str(isc) , r[isc,:] 
                                #print 'v for ', str(isc) , v[isc,:] 
                    
                            # Create main input file for collision run (forward). Same as backward but with the initial r/v equal to the r/v of backward (and in collision mode in order to calculate the probability of collision)
                            ## Create collision file (that's where the initial r/v is written)
                            collision_filename = main_input_filename.replace("backward", "collision")
                            spock_collision(collision_filename, r*1000., v*1000., nb_ens, min_dist_ca, min_dist_coll) # r and v in m and m/s
                            ## Create main input file
                            irun_quantile = -1
                            while irun_quantile < nb_quantile-1:
                                irun = irun + 1
                                print 'irun ' + str(irun) + ' out of ' + str(nb_run-1) + ' runs'
                                irun_quantile = irun_quantile + 1

                                density_mode_forward = ['swpc_mod', quantile_filename_array[irun_quantile]]# note that although backward was static here we use option swpc_mod so you could think f10.7 and Ap vary with time, which could prevent a collision. But I changed the f10.7 and Ap files so that they show constant values, equal to 100 and 12 all the time. Like this, the quantil 50% is the same as static 100/12
                                orbit_type = ['collision', collision_filename]
                                date_end_forward = date_end + timedelta(hours = 4) # add a few hours because otherwise SpOCK won't detect the TCA (bug...)
                                date_end_forward_str = datetime.strftime(date_end_forward, "%Y-%m-%dT%H:%M:%S")
                                main_input_filename = 'alt' + str(alt[0]) + '-' + str(alt[1]) + '_inc' + str(inc[0])+ '-' + str(inc[1]) + '_arg_per' + str(arg_per[0])+ '-' + str(arg_per[1]) + '_raan' + str(raan[0])+ '-' + str(raan[1]) + '_true_ano' + str(true_ano[0])+ '-' + str(true_ano[1]) + '_ecc' + str(ecc[0])+ '-' + str(ecc[1]) + '_f107' + str(f107) + '_quantile' + str(quantile_array[irun_quantile]) + '_forward.txt'
                                spock_main_input(
                                    main_input_filename,
                                    # for TIME section
                                    date_start_str,
                                    date_end_forward_str,
                                    dt_forward,
                                    # for SPACECRAFT section
                                    2,
                                    '0',
                                    mass,
                                    geometry_filename, 
                                    # for ORBIT section
                                    orbit_type,
                                    # for FORCES section
                                    2,
                                    'drag',
                                    density_mode_forward,
                                    # for OUTPUT section
                                    'out',
                                    86400, # the last time step is the only we care about here and is always printed
                                    # for ATTITUDE section
                                    "nadir",
                                    # for GROUNDS_STATIONS section
                                    "0",#"my_ground_stations.txt"
                                    # for SPICE section
                                    spice_path_forward,
                                    # for DENSITY_MOD section
                                    [1, "\n#OUTPUT_ENSEMBLES", "eci_r, collision"]
                                    )

                                # write in a file the filename of the forward simulation. This file will be used on Pleiades to run these forward simulations
                                print >> file_list_forward_simu, main_input_filename

                                # pack the main input and collision files into a .tgz folder that will be sent to Pleiades
                                os.system("cp " + main_input_filename + " " + collision_filename + " " + send_folder)

                            #print ''

file_list_forward_simu.close()
os.system("cp " + filename_list_forward_simu + " " + send_folder)
os.system("tar -zcvf send_pleiades_all/" + send_folder + ".tgz " + send_folder)
