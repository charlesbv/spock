
# This cretes input files for collision of different encounter geometry. It first propagates backward SpOCK with different orbital elements representing different encouter geometries. It then take the last time step and initialize the collision runs with this last step. It actually creates 10 forward simu for each backward simu: each of these 10 simu is using the option swpc_mod wth a different quantile
# Assumptions:
# - see section #PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT


import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/kalman/spock_development_new_structure_kalman_dev/srcPython")
import numpy as np
from struct import *
from datetime import datetime, timedelta
import os
from read_input_file import *
from read_output_file import *
from cygnss_read_spock_spec import *
from datetime import datetime, timedelta
from spock_main_input import *
from spock_collision import *
from spock_swpc_f107_ap_file import *

#PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
date_start_str = '2017-12-01T00:00:00' # when to sstart the collision runs (forward (so this is really the start date)). SO the f107 and ap are supposed to be predicted beyond this date (by this it means that SpOCK considers this date as the start date for predictions (so it can apply spwc_mod on the nominal f107 and ap. This also means that ap and f107 are read from the swpc prediction file. This file is threfore createed using spock_swpc_f107_ap_file wheere ap and f107 are input from f107_arr and ap_arr))
tca = 36. # number of hours after date_start the collision occurs
alt_arr = [[300,300],[350,350],[400,400],[450,450],[500,500]] # altitude of each sc for each run
# alt_arr = [[360,360],[380,380],
#            [400,400],[420,420],[440,440]] # altitude of each sc for each run
inc_arr = [[90,45]] # inclination of each sc for each run
arg_per_arr = [[0,0]] # argument of perifee  of each sc for each run
raan_arr = [[0,0]] # RAAN  of each sc for each run
true_ano_arr = [[0,0]] # true anomaly of each sc for each run
ecc_arr = [[0,0]] # eccentricity of each sc for each run
f107_arr = [100] # f107 for each run
ap_arr = [12] # ap for each run
# f107_arr = [ 80,  90, 100, 110, 120] # f107 for each run
# ap_arr = [ 10, 15, 20, 25, 30 ] # ap for each run
bc_arr = [0.001, 0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128] # for examples of bc, see below ('examples of bc'). to input these bc in SpOCK, a cosntant area of 1 m2 and cd of 2.0 are chosen, and only the mass is varied: m = 2.0 * 1.0 / bc_arr[i]
# examples of bc:
# cygnss bc nadir: 2.2 * 0.1 / 28 = 0.008
# cygnss bc high drag (actually assume pitch at 90 deg): 2.2 * 0.8 / 28 = 0.06
# cubesat 1U bc: 2.2 *  0.01 / 1 = 0.02
# iss: 2.2 * 500 / 420000. = 0.002 (assume 500 m^2 frontal area (70 m (length) by 110 (width) by 20 (height)))
# mir: 2.2 * 50 / 130000. = 0.0008 (assume 50 m^2 frontal area (19 m (length) by 31 m (width)  by 28 m (height)))
# one_plate_old: 2.0 * 1 / 100 = 0.02 = 1U cubesat
# end of examples of bc:
nb_ens = 50000 # number of ensemble for the collision (collision is over nb_ens*nb_ens)
min_dist_ca = 1000. # distance under which a close approach is flagged. in m
min_dist_coll = 1. # distance under which a collision is recorded. in m
path_to_mpirun = "mpirun" # where mpirun is in the system
spice_path_backward = "/Users/cbv/cspice/data" # where spice is in the system to run the backward simu
spice_path_forward = "/home1/cbussy/cspice/data" # where spice is in the system to run the forward simu
send_folder = 'sendbc' # name of foler that will include the forward simulations to run on Pleiades. if this folder doesn't exist then this scripts creates it
filename_list_forward_simu = "encounterbc.txt"
#end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT


cd = 2.0 # don't change unless you know what you're doing
area = 1.0 # in m^2. don't change unless you know what you're doing 
if len(ap_arr) != len(f107_arr):
    print "***! F10.7 and Ap arrays have to be the same length. The program will stop. !***"; exit(0);

if os.path.isdir(send_folder) == False:
    os.system("mkdir " + send_folder)
    "The folder " + send_folder + " has been created."

quantile_array = [1,2,3,4,5,6,7,8,9] # same quantile for ap and f10.7
nb_quantile = len(quantile_array)
quantile_filename_array = []
for iquantile in range(nb_quantile):
    quantile_filename_array.append("quartile_f107_" + str(quantile_array[iquantile]) + "_quartile_ap_" + str(quantile_array[iquantile]) + '.txt')
    
date_start = datetime.strptime(date_start_str, "%Y-%m-%dT%H:%M:%S")
date_start_str_for_f107_ap = date_start_str[2:4] + '-' + date_start_str.split('-')[1] + '-' + date_start_str.split('-')[2].split('T')[0]
date_end = date_start + timedelta(hours = tca)
nb_alt = len(alt_arr)
nb_inc = len(inc_arr)
nb_arg_per = len(arg_per_arr)
nb_raan = len(raan_arr)
nb_true_ano = len(true_ano_arr)
nb_ecc = len(ecc_arr)
nb_f107 = len(f107_arr)
nb_ap = len(ap_arr)
nb_bc = len(bc_arr)

nb_run = nb_alt*nb_inc*nb_arg_per*nb_raan*nb_true_ano*nb_ecc*nb_bc*nb_f107*nb_quantile

# Create backward and forward simulations  (1. create backward; 2. run it; 3. create forward)
date_end_str = datetime.strftime(date_end, "%Y-%m-%dT%H:%M:%S")
dt_backward = -1. # need one second dt for the backward otherwise when forward then don't get collision (for forard ok to have 10s)
dt_forward = 10.
nb_sc  = 2
#mass = 100 # has been commented when cbv started to ballistic coefficientx
geometry_filename = 'one_plate_old.txt'

irun = -1
irun_alt = -1

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
                        irun_bc = -1
                        while irun_bc < nb_bc-1:
                            irun_bc = irun_bc + 1
                            bc = bc_arr[irun_bc]
                            irun_f107 = -1
                            irun_ap = -1
                            while irun_f107 < nb_f107-1: # assume ap and f107 vary the same
                                irun_f107 = irun_f107 + 1
                                f107 = f107_arr[irun_f107]
                                irun_ap = irun_ap + 1
                                ap = ap_arr[irun_ap]
                                density_mode_backward = ['static', str(f107), '100', str(ap)]  # keep f107a at 100 

                                # Create backward main input file
                                mass = cd * area / bc
                                orbit_type = ['oe', str(alt[0])  + ' ' + str(inc[0])  + ' ' + str(arg_per[0])  + ' ' + str(raan[0])  + ' ' + str(true_ano[0])  + ' ' + str(ecc[0]), str(alt[1])  + ' ' + str(inc[1])  + ' ' + str(arg_per[1])  + ' ' + str(raan[1])  + ' ' + str(true_ano[1])  + ' ' + str(ecc[1])]
                                main_input_filename = 'alt' + str(alt[0]) + '-' + str(alt[1]) + '_inc' + str(inc[0])+ '-' + str(inc[1]) + '_arg_per' + str(arg_per[0])+ '-' + str(arg_per[1]) + '_raan' + str(raan[0])+ '-' + str(raan[1]) + '_true_ano' + str(true_ano[0])+ '-' + str(true_ano[1]) + '_ecc' + str(ecc[0])+ '-' + str(ecc[1]) + '_bc' + str(bc) +  '_f107' + str(f107) + '_ap' + str(ap) + '_backward.txt'
    #                             if ( os.path.isfile(main_input_filename) == True ): 
    #                                 run_it = 0
    #                             else:# don't run the simu if you already ran it

    #                             run_it = 1
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

                                #os.system(path_to_mpirun + " -np 2 spock_coll " + main_input_filename)

                                ## read position and velocity at last time step (ie at date_start since the run was backward)
                                print 'reading results of ' + main_input_filename 
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
    #                             ## Create main input file
                                irun_quantile = -1
                                while irun_quantile < nb_quantile-1:
                                    irun = irun + 1
                                    print 'irun ' + str(irun) + ' out of ' + str(nb_run-1) + ' runs'
                                    irun_quantile = irun_quantile + 1
                                    # create the swpc_predictions_f107_ap.txt prediction file with the ap anf f107
                                    f107_ap_pred_filename = 'swpc_predictions_f107' + str(f107) + '_ap' + str(ap) + '.txt'
                                    spock_swpc_f107_ap_file(f107_ap_pred_filename, f107, ap, date_start_str_for_f107_ap)

                                    density_mode_forward = ['swpc_mod', quantile_filename_array[irun_quantile] + ' ' + f107_ap_pred_filename]# note that although backward was static here we use option swpc_mod so you could think f10.7 and Ap vary with time, which could prevent a collision. But I changed the f10.7 and Ap files so that they show constant values, equal to 100 and 12 all the time. Like this, the quantil 50% is the same as static 100/12
                                    orbit_type = ['collision', collision_filename]
                                    date_end_forward = date_end + timedelta(hours = 4) # add a few hours because otherwise SpOCK won't detect the TCA (bug...)
                                    date_end_forward_str = datetime.strftime(date_end_forward, "%Y-%m-%dT%H:%M:%S")
                                    main_input_filename = 'alt' + str(alt[0]) + '-' + str(alt[1]) + '_inc' + str(inc[0])+ '-' + str(inc[1]) + '_arg_per' + str(arg_per[0])+ '-' + str(arg_per[1]) + '_raan' + str(raan[0])+ '-' + str(raan[1]) + '_true_ano' + str(true_ano[0])+ '-' + str(true_ano[1]) + '_ecc' + str(ecc[0])+ '-' + str(ecc[1]) + '_bc' + str(bc) + '_f107' + str(f107) + '_ap' + str(ap) + '_quantile' + str(quantile_array[irun_quantile]) + '.txt'
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

                                    # pack the main input and collision files into a .tgz folder that will be sent to Pleiades
                                    os.system("cp " + main_input_filename + " " + collision_filename + " " + f107_ap_pred_filename + " " + send_folder)
                                # write in a file the filename of the forward simulation. This file will be used on Pleiades to run these forward simulations. Note that we write only the root of the filenmae of all qunatiles
#                                     # !!!!!!!!!! only if quamtile is 1, 5, or 9 (this is teporary)
#                                     if ( ( irun_quantile == 1 ) | ( irun_quantile == 5 ) | ( irun_quantile == 9 ) ):
                                    print >> file_list_forward_simu, main_input_filename
                                







                            #print ''
#                                raise Exception
file_list_forward_simu.close()
os.system("cp " + filename_list_forward_simu + " " + send_folder)
os.system("tar -zcvf send_pleiades_all/" + send_folder + ".tgz " + send_folder)
