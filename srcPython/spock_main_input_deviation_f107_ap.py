# This script creates SpOCK main input files and put them in send_folder. the ony difference between two files is the name of the f107/ap deviation file (swpc_mod). It also write a file with all then anmes for these fmain input files. This file will be sued on Pliades to run all the simulations automatically.
# the f107.ap deviairton files were created in deviation_swpc_f107_ap.py

import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/kalman/spock_development_new_structure_kalman_dev/srcPython")
import numpy as np
from spock_main_input import *
import os

f107_bin_range = np.array([-15.        , -10.71428571,  -6.42857143,  -2.14285714,         2.14285714,   6.42857143,  10.71428571,  15.        ]) # has to be the same as in deviation_swpc_f107_ap.py
ap_bin_range = f107_bin_range


send_folder = 'sendpaper4_no_coll' # name of foler that will include the forward simulations to run on Pleiades. if this folder doesn't exist then this scripts creates it
filename_list_runs = "run_paper4.txt"


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


dt = '10.'
for if107 in range(nb_f107):
    for iap in range(nb_ap):
        main_input_filename = 'f107_' + format(f107_bin_center[if107], ".0f") + '_ap_' + format(ap_bin_center[iap], ".0f") + '.txt'
        deviation_f107_ap_filename = 'deviation_f107_' + format(f107_bin_center[if107], ".0f") + '_ap_' + format(ap_bin_center[iap], ".0f") + '.txt'
        density_mode = ['swpc_mod', deviation_f107_ap_filename + ' swpc_predictions_f107_ap_paper4.txt']
        spock_main_input(
            main_input_filename,
            # for TIME section
            '2016-11-26T00:00:00',
            '2016-11-29T00:00:00',
            dt,
            # for SPACECRAFT section
            2,
            '0',
            100,
            'one_plate.txt',
            # for ORBIT section
            ['collision', 'gravity_order2_inclination_but_cov_same_as_case6_cov_min_dist_coll_1_3m_for_lighter_sc_50000ens.txt'],
            # for FORCES section
            2,
            'drag',
            density_mode,
            # for OUTPUT section
            'out',
            86400, # the last time step is the only we care about here and is always printed
            # for ATTITUDE section
            "nadir",
            # for GROUNDS_STATIONS section
            "0",#"my_ground_stations.txt"
            # for SPICE section
            spice_path_pleiades,
            # for DENSITY_MOD section
            [1, "\n#OUTPUT_ENSEMBLES", "eci_r, collision"]
            )

        print >> file_list_runs, main_input_filename
        os.system("cp " + main_input_filename + " " + send_folder)

file_list_runs.close()
os.system("cp " +  filename_list_runs + " " + send_folder)
os.system("tar -zcvf " +  send_folder + ".tgz " + send_folder)

