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

# This script creates SpOCK main input files and put them in send_folder. the ony difference between two files is the name of the f107/ap deviation file (swpc_mod). It also write a file with all then anmes for these fmain input files. This file will be sued on Pliades to run all the simulations automatically.
# the f107.ap deviairton files were created in deviation_swpc_f107_ap.py
# it was created in Feb 2018 for paper CA2 and revision ot paper CA1 (to compute PDFs and CDFs of Pc, as suggested by one of the reviewers).
# cbv used it again starting July 2018 for paper CA3 (CYGNSS collision avoidance). For the exact same  version nefore July 2018: spock_main_input_deviation_f107_ap_before_july2018.py

import sys
#sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/srcPython")
sys.path.append("/home1/cbussy/Code/spock/srcPython")
import numpy as np
from spock_main_input import *
import os


# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
send_folder = 'FM07_nov01_30000ens' # name of foler that will include the forward simulations to run on Pleiades. if this folder doesn't exist then this scripts creates it
filename_list_runs = "runlist_FM07_nov01_30000ens.txt"
start_date = '2017-11-01T05:12:00.000000'
end_date = '2017-11-05T06:00:00.000000'
dt = '20.'
vcm_filename1 = 'data/41890_27434_20171105_043841/41890_20171101_051106.txt'
vcm_filename2 = 'data/41890_27434_20171105_043841/27434_20171101_043101.txt'
nb_ensemble = 30000
dist_flag_coll = 2200. # in m. if distance between the 2 main sc is below this number, then the distance between all ensembels are computed
dist_min_coll = 20. # hard body radius. if 2 sc are less than this distance apart, it's a collision
f107_ap_pred_filename = "f107_ap_pred_2017-11-01.txt" # 45d predictions by SWPC 
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT


f107_bin_range = np.array([-15.        , -10.71428571,  -6.42857143,  -2.14285714,         2.14285714,   6.42857143,  10.71428571,  15.        ]) # has to be the same as in deviation_swpc_f107_ap.py
ap_bin_range = f107_bin_range




if os.path.isdir(send_folder) == False:
    os.system("mkdir " + send_folder)
    "The folder " + send_folder + " has been created."

file_list_runs = open(filename_list_runs, "w" )

spice_path_local = "/Users/cbv/cspice/data" # /home1/cbussy/cspice/data
spice_path_big = "/raid4/cbv/cspice/data"
spice_path_pleiades = "/home1/cbussy/cspice/data"


f107_bin_center = (f107_bin_range[:-1] + np.roll(f107_bin_range, -1)[:-1])/2
ap_bin_center = (ap_bin_range[:-1] + np.roll(ap_bin_range, -1)[:-1])/2

nb_f107 = len(f107_bin_center)
nb_ap = len(ap_bin_center)
nb_comb = nb_f107 * nb_ap

coll_param = str(dist_flag_coll) + ' ' + str(dist_min_coll) + ' ' + str(nb_ensemble)
for if107 in range(nb_f107):
    for iap in range(nb_ap):
        main_input_filename = send_folder + '_f107_' + format(f107_bin_center[if107], ".0f") + '_ap_' + format(ap_bin_center[iap], ".0f") + '.txt'
        deviation_f107_ap_filename = 'deviation_f107_' + format(f107_bin_center[if107], ".0f") + '_ap_' + format(ap_bin_center[iap], ".0f") + '.txt'
        density_mode = ['swpc_mod', deviation_f107_ap_filename + ' ' + f107_ap_pred_filename]
        spock_main_input(
            main_input_filename,
            # for TIME section
            start_date,
            end_date,
            dt,
            # for SPACECRAFT section
            2,
            '0',
            28,
            'cygnss_geometry_2016_acco08.txt',
            # for ORBIT section
            ['collision_vcm', vcm_filename1,  vcm_filename2, coll_param],
            # for FORCES section
            '36',
            'drag solar_pressure moon_gravity sun_gravity',
            density_mode,
            # for OUTPUT section
            'out',
            20, 
            # for ATTITUDE section
            "nadir",
            # for GROUNDS_STATIONS section
            "0",#"my_ground_stations.txt"
            # for SPICE section
            spice_path_pleiades,
            # for DENSITY_MOD section
            [1]#, "\n#OUTPUT_ENSEMBLES", "eci_r, collision"]
            )

        print >> file_list_runs, main_input_filename
        os.system("cp " + main_input_filename + " " + send_folder)

file_list_runs.close()
os.system("cp " +  filename_list_runs + " " + send_folder)
#os.system("tar -zcvf " +  send_folder + ".tgz " + send_folder)

