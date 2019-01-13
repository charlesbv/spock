# This script runs SpOCK with ensembles on the initial state (r, v ECI). The initial mean r,v is taken from GPS measurements. It's a similiar script ot ~/Google Drive/Work/PhD/Research/Code/cygnss/eclipse/ensemble_ini_state/spock_odtk_ensemble_dev.py but this one here uses GPS measuremnts while the other used Kyle Nave ODTK states.
# ASSUMPTIONS
# first run cygnss_convert_swri_att_to_spock.py to convert the GPS measurements into files readible by SpOCK (attitude)
#- see section "PARAMETERS TO SET UP BEFORE RUNNIG THIS SCRIPT"
#- run SpOCK with a 1s time step 
#- rho_mod_arr must be such that the coeff 1 is included



# PARAMETERS TO SET UP BEFORE RUNNIG THIS SCRIPT
isbig = 1 # if runnign script from Big
interval = 1.5 # interval of time to compare the two trajectories (data and SpOCK). In hours
nb_ensemble_ini_state = 50
sigma_x = 0.0296194067316 #0.0592388134632#1 # in m # standard deviation
sigma_y = 0.0296194067316#0.0592388134632#1 # in m
sigma_z = 0.0296194067316#0.0592388134632#1 # in m
sigma_vx = 0.00296194067316#0.00592388134632#0.1 # in m/s
sigma_vy = 0.00296194067316#0.00592388134632#0.1 # in m/s
sigma_vz = 0.00296194067316#0.00592388134632#0.1 # in m/s
rho_mod_min = 0.7# min rho_mod
rho_mod_max = 1. # max rho_mod
drho_mod = 0.3 # step in rho_mod -> rho_mod varies from rho_mod_min to rho_mod_max by values of drho_mod
#spock_FM5_20171216_eng_adcs_query-13528.txt
# end of PARAMETERS TO SET UP BEFORE RUNNIG THIS SCRIPT

import numpy as np
import sys
if isbig == 1:
    sys.path.append("/home/cbv/code/spock/srcPython")
    path_mpirun = '/usr/local/bin/mpirun'
    spice_path = '/raid4/cbv/cspice/data'
    nb_proc = 12
else:
    sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/srcPython")
    path_mpirun = 'mpirun'
    spice_path = '/Users/cbv/cspice/data'
    nb_proc = 4
import os
from read_input_file import *
from read_output_file import *
from spock_main_input import *
from orbit_average import *
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
if isbig != 1:
    from convert_cygnss_obs_ecef_to_eci import *





# Read r/v of observations
obs_rv_filename = 'HD_data/spock_FM5_20171216_eng_pvt_query-13527_1800tomorrow.txt'
obs_att_filename = 'HD_data/spock_FM5_20171216_eng_adcs_query-13528_1800tomorrow.txt'

# Convert ECEF file to ECI file
if isbig == 1:
    obs_rv_filename_eci = obs_rv_filename.replace('.txt','_eci.txt')
else:
    obs_rv_filename_eci = convert_cygnss_obs_ecef_to_eci(obs_rv_filename)

# Read observation ECI r/v 
obs_rv_file = open(obs_rv_filename_eci)
read_obs_rv_file = obs_rv_file.readlines()
nb_header = 0
while (read_obs_rv_file[nb_header].split()[0] != '#START'):
    nb_header = nb_header + 1
nb_header = nb_header + 1
nb_obs = len(read_obs_rv_file) - nb_header
date_obs = []
date_obs_str = []
r_obs = np.zeros([nb_obs, 3])
v_obs = np.zeros([nb_obs, 3])

for iobs in range(nb_obs):
    date_obs_str.append( read_obs_rv_file[iobs + nb_header].split()[0] )
    date_obs.append( datetime.strptime(date_obs_str[-1], "%Y-%m-%dT%H:%M:%S" ) )
    r_obs[iobs, 0] = np.float( read_obs_rv_file[iobs + nb_header].split()[1] ) 
    r_obs[iobs, 1] = np.float( read_obs_rv_file[iobs + nb_header].split()[2] ) 
    r_obs[iobs, 2] = np.float( read_obs_rv_file[iobs + nb_header].split()[3] ) 
    v_obs[iobs, 0] = np.float( read_obs_rv_file[iobs + nb_header].split()[4] ) 
    v_obs[iobs, 1] = np.float( read_obs_rv_file[iobs + nb_header].split()[5] ) 
    v_obs[iobs, 2] = np.float( read_obs_rv_file[iobs + nb_header].split()[6] ) 


# Run SpOCK: initial r/v is given by observations + ensemble with std given by x_sigma, y_sigma, etc
date_obs_start_str = date_obs_str[0]
date_obs_start= datetime.strptime(date_obs_start_str, "%Y-%m-%dT%H:%M:%S")
date_obs_end_str = date_obs_str[-1]
date_obs_end= datetime.strptime(date_obs_end_str, "%Y-%m-%dT%H:%M:%S")
interval_sec = interval * 3600.
nb_interval = (int) ( ( date_obs_end - date_obs_start ).total_seconds()/ ( interval_sec ) )

nb_seconds_since_start = []
distance = []
distance_ref = []

date_start = date_obs_start
date_end = date_start + timedelta(seconds = interval_sec)
date_end_str = datetime.strftime(date_end, "%Y-%m-%dT%H:%M:%S")
date_start_str = datetime.strftime(date_start, "%Y-%m-%dT%H:%M:%S")

index_obs_interval_start = 0 # !!!!!!!!!!! to change
rho_mod_arr = np.arange(rho_mod_min, rho_mod_max+drho_mod, drho_mod)
nb_rho = len(rho_mod_arr)

irho_equal_1 = np.where(rho_mod_arr == 1)[0][0]
if len(np.where(rho_mod_arr == 1)[0]) == 0:
    print "***! The array of density factors must include 1. The program will stop. !***"; raise Exception


## Create SpOCK main input file: same epoch and initial r/v
dt  = 1
dt_output = 1
gravity_order = 50 # !!!!!!!!!! put 50
rho_mod = 1
it = -1

which_ens_min_dist = []
min_mean_dist_itime_iens = []
min_mean_dist_itime_iens.append(1e10)
sigma_x_list = []
sigma_y_list = []
sigma_z_list = []
sigma_vx_list = []
sigma_vy_list = []
sigma_vz_list = []
sigma_x_list_acceptable = []
sigma_y_list_acceptable = []
sigma_z_list_acceptable = []
sigma_vx_list_acceptable = []
sigma_vy_list_acceptable = []
sigma_vz_list_acceptable = []

while min_mean_dist_itime_iens[-1]*1000 > 15: # iterate to minimize the distance to observation
    it = it + 1
    main_input_filename = date_start_str.replace(":","_") + '_' + date_end_str.replace(":","_")+ '_iteration' + str(it) + '.txt'

    if it == 0: # first iteration is initizliae with the r0/v0 of the observations
        r0 = format(r_obs[index_obs_interval_start, 0]*1000, '.14e')
        r1 = format(r_obs[index_obs_interval_start, 1]*1000, '.14e')
        r2 = format(r_obs[index_obs_interval_start, 2]*1000, '.14e')
        v0 = format(v_obs[index_obs_interval_start, 0]*1000, '.14e')
        v1 = format(v_obs[index_obs_interval_start, 1]*1000, '.14e')
        v2 = format(v_obs[index_obs_interval_start, 2]*1000, '.14e')

        # !!!!!!! remove block below
        # r0 = '-2.54076684997000e+06'
        # r1 = '-5.06267228702000e+06'
        # r2 = '-3.95089226458000e+06'
        # v0 = '6.76839628300000e+03'
        # v1 = '-3.44599994200000e+03'
        # v2 = '4.16930640000000e+01'
        r0 = '-2.54076675527000e+06'
        r1 = '-5.06267223919000e+06'
        r2 = '-3.95089235289000e+06'
        v0 = '6.76839873800000e+03'
        v1 = '-3.44599893200000e+03'
        v2 = '4.16982230000000e+01'

        # !!!!!!! end of remove block below
    else: # following iterations are initizliaze with r0/v0 of the optimized trajectory from the preivous iteration
    # the standard deviation on r0/v0 is divided by iter_std_fac but only ig the new distance to observations is smaller than at the previous iteration

        if min_mean_dist_itime_iens[-1] < min_mean_dist_itime_iens[-2]*0.85:
            if it == 1:
                iter_std_fac = 2.
            else:
                iter_std_fac = min_mean_dist_itime_iens[-2] / min_mean_dist_itime_iens[-1]
            sigma_x = np.sqrt(sigma_x)/iter_std_fac # we taje the square after so need to take the square root here
            sigma_y = np.sqrt(sigma_y)/iter_std_fac
            sigma_z = np.sqrt(sigma_z)/iter_std_fac
            sigma_vx = np.sqrt(sigma_vx)/iter_std_fac
            sigma_vy = np.sqrt(sigma_vy)/iter_std_fac
            sigma_vz = np.sqrt(sigma_vz)/iter_std_fac
            r0 = format(r_spock_ok[0, (int)(which_ens_min_dist[it-1]), 0]*1000, '.14e')
            r1 = format(r_spock_ok[0, (int)(which_ens_min_dist[it-1]), 1]*1000, '.14e')
            r2 = format(r_spock_ok[0, (int)(which_ens_min_dist[it-1]), 2]*1000, '.14e')
            v0 = format(v_spock_ok[0, (int)(which_ens_min_dist[it-1]), 0]*1000, '.14e')
            v1 = format(v_spock_ok[0, (int)(which_ens_min_dist[it-1]), 1]*1000, '.14e')
            v2 = format(v_spock_ok[0, (int)(which_ens_min_dist[it-1]), 2]*1000, '.14e')

        else:  # if the decrease in distance hasn't been good enough, then don't stay too close to the new solution because otherwise you're going to get stuck there. So go back to the previous sigma. Also go back to the r0/v0 that gave a smaller distance
            print 'previous iteration not acceptable'
            sigma_x = sigma_x_list_acceptable[-1] 
            sigma_y = sigma_y_list_acceptable[-1]
            sigma_z = sigma_z_list_acceptable[-1]
            sigma_vx = sigma_vx_list_acceptable[-1]
            sigma_vy = sigma_vy_list_acceptable[-1]
            sigma_vz = sigma_vz_list_acceptable[-1]
            r0 = previous_r0
            r1 = previous_r1
            r2 = previous_r2
            v0 = previous_v0
            v1 = previous_v1
            v2 = previous_v2


        # as the distance to observation decreases, increase order of gravity model. 
       #this is because the error on r0/v0 has been reduced so now the error on gravity model has as much effects as the error on r0/v0
        if min_mean_dist_itime_iens[-1]*1000 < 60:
            gravity_order = 50

    print '\nNOW'
    print r0,r1,r2
    print v0,v1,v2
    print sigma_x, sigma_vx
    # SpOCK inital state uncetainty file
    # In SpOCK collision file, the diagnoal terms of the covaraiance matrix represent the variance, which is the sqaure of the strand deviation
    sigma_x = sigma_x**2
    sigma_y = sigma_y**2
    sigma_z = sigma_z**2
    sigma_vx = sigma_vx**2
    sigma_vy = sigma_vy**2
    sigma_vz = sigma_vz**2

    sigma_x_list.append( np.sqrt( sigma_x ) )
    sigma_y_list.append( np.sqrt( sigma_y ) )
    sigma_z_list.append( np.sqrt( sigma_z ) )

    sigma_vx_list.append( np.sqrt( sigma_vx ) )
    sigma_vy_list.append( np.sqrt( sigma_vy ) )
    sigma_vz_list.append( np.sqrt( sigma_vz ) )


    filename_ini_state = main_input_filename.replace('.txt', '_ini_state.txt')
    file_ini_state = open(filename_ini_state, "w+")
    print >> file_ini_state, "#STATE_ECI"
    print >> file_ini_state, '(' + r0 + '; ' + r1 + '; ' + r2 + ') (' + v0 + '; ' + v1 + '; ' + v2 + ')' 
    print >> file_ini_state, '(1000000000; 1000000000; 1000000000; 5000; 0; 0)'  # don't care about second satellite but hacve to put one because it's the collision mode
    print >> file_ini_state, "\n#COVARIANCE"
    print >> file_ini_state, '((' + str(sigma_x) + ';0 ; 0; 0; 0; 0);'
    print >> file_ini_state, '(0; ' + str(sigma_y) + ';0 ; 0; 0; 0);'
    print >> file_ini_state, '(0; 0; ' + str(sigma_z) + ';0 ; 0; 0);'
    print >> file_ini_state, '(0; 0; 0; ' + str(sigma_vx) + ';0 ; 0);'
    print >> file_ini_state, '(0; 0; 0; 0; ' + str(sigma_vy) + ';0 );'
    print >> file_ini_state, '(0; 0; 0; 0; 0; ' + str(sigma_vz) + '))'
    # don't care about second sc (so put same as first sc)
    print >> file_ini_state, '((' + str(sigma_x) + ';0 ; 0; 0; 0; 0);'
    print >> file_ini_state, '(0; ' + str(sigma_y) + ';0 ; 0; 0; 0);'
    print >> file_ini_state, '(0; 0; ' + str(sigma_z) + ';0 ; 0; 0);'
    print >> file_ini_state, '(0; 0; 0; ' + str(sigma_vx) + ';0 ; 0);'
    print >> file_ini_state, '(0; 0; 0; 0; ' + str(sigma_vy) + ';0 );'
    print >> file_ini_state, '(0; 0; 0; 0; 0; ' + str(sigma_vz) + '))'

    print >> file_ini_state, "\n#NB_ENSEMBLES_COLLISION\n" + str(nb_ensemble_ini_state) 
    print >> file_ini_state, "\n#MIN_DISTANCE_CLOSE_APPROACH\n10000\n\n#MIN_DISTANCE_COLLISION\n1.3"

    file_ini_state.close()
    # SpOCK main input file
    spock_main_input( # need to be in spokc/srcPython to run this script   
        main_input_filename,
        # for TIME section
            date_start_str,
        date_end_str,
        dt,
        # for SPACECRAFT section
                1,
        '0',
        29,
        "cygnss_geometry_2016_acco08.txt", 
        # for ORBIT section
            ['collision', filename_ini_state ],
        # for FORCES section
            gravity_order, # !!!!!!!!!!! put back 20
        "drag solar_pressure sun_gravity moon_gravity", # !!!!!!!!!!!!! put back to "drag sun_gravity moon_gravity"
        "dynamic",
        # for OUTPUT section
                "out",
        dt_output, 
        # for ATTITUDE section
        obs_att_filename,
        # for GROUND_STATIONS section
                "0",
        # for SPICE section
                spice_path,
        # FOR #DENSITY_MOD section
                rho_mod
    )

    # add in SpOCK main input file the section #OUTPUT_ENSEMBLES
    file_spock = open(main_input_filename, "a")
    print >> file_spock, "#OUTPUT_ENSEMBLES\neci_r, eci_v"
    file_spock.close()

    #os.system("rm -Rf " + main_input_filename.replace(".txt", "")) #!!!!!!!! remove this line
    ## Run SpOCK
    #os.system(path_mpirun + ' -np ' + str(nb_proc) + ' spock ' + main_input_filename)

    ## concatenate proc files
    print "Concatenating processor files"
    #os.system(path_mpirun + ' -np ' + str(nb_proc) + ' python new_mpi_concatenate_proc.py ' + main_input_filename )

    # Read the position and velocity predicted by SpOCK
    isc = 0
    var_in, var_in_order = read_input_file(main_input_filename)

    output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
    output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
    var_to_read = ["position", "velocity"]
    var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
    date_spock = np.array(var_out[find_in_read_input_order_variables(var_out_order, 'date')])
    date_datetime_spock = np.array(var_out[find_in_read_input_order_variables(var_out_order, 'date_datetime')])
    date_datetime_round_sec_spock = np.array(var_out[find_in_read_input_order_variables(var_out_order, 'date_datetime_round_sec')])
    n_spock = len(date_spock)
    ensembles_to_output = var_in[find_in_read_input_order_variables(var_in_order, 'ensembles_to_output')];

    if 'eci_r' in ensembles_to_output: #Read eci position of ensembles
        filename_ens_eci_x = output_file_path_list[isc] + 'ensemble/ensemble_x_eci_' + output_file_name_list[isc]
        file_eci_x = open(filename_ens_eci_x)
        read_file_eci_x = file_eci_x.readlines()

        nb_header_ens_eci_x = 0
        while  read_file_eci_x[nb_header_ens_eci_x].split()[0] != '#START':
            nb_header_ens_eci_x = nb_header_ens_eci_x + 1
        nb_header_ens_eci_x = nb_header_ens_eci_x + 1
        nb_ensemble_ini_state_corrected = len(read_file_eci_x[nb_header_ens_eci_x].split())-2
        r_spock = np.zeros([n_spock, nb_ensemble_ini_state_corrected, 3])
        for itime_ens in range(n_spock):
            for iens in range(nb_ensemble_ini_state_corrected):
                r_spock[itime_ens, iens, 0] = read_file_eci_x[nb_header_ens_eci_x+itime_ens].split()[2+iens]
        file_eci_x.close()

        filename_ens_eci_y = output_file_path_list[isc] + 'ensemble/ensemble_y_eci_' + output_file_name_list[isc]
        file_eci_y = open(filename_ens_eci_y)
        read_file_eci_y = file_eci_y.readlines()

        nb_header_ens_eci_y = 0
        while  read_file_eci_y[nb_header_ens_eci_y].split()[0] != '#START':
            nb_header_ens_eci_y = nb_header_ens_eci_y + 1
        nb_header_ens_eci_y = nb_header_ens_eci_y + 1
        nb_ensemble_ini_state_corrected = len(read_file_eci_y[nb_header_ens_eci_y].split())-2

        for itime_ens in range(n_spock):
            for iens in range(nb_ensemble_ini_state_corrected):
                r_spock[itime_ens, iens, 1] = read_file_eci_y[nb_header_ens_eci_y+itime_ens].split()[2+iens]
        file_eci_y.close()

        filename_ens_eci_z = output_file_path_list[isc] + 'ensemble/ensemble_z_eci_' + output_file_name_list[isc]
        file_eci_z = open(filename_ens_eci_z)
        read_file_eci_z = file_eci_z.readlines()

        nb_header_ens_eci_z = 0
        while  read_file_eci_z[nb_header_ens_eci_z].split()[0] != '#START':
            nb_header_ens_eci_z = nb_header_ens_eci_z + 1
        nb_header_ens_eci_z = nb_header_ens_eci_z + 1
        nb_ensemble_ini_state_corrected = len(read_file_eci_z[nb_header_ens_eci_z].split())-2

        r_spock_ref = var_out[find_in_read_input_order_variables(var_out_order, 'position')]
        v_spock_ref = var_out[find_in_read_input_order_variables(var_out_order, 'velocity')]

        for itime_ens in range(n_spock):
            for iens in range(nb_ensemble_ini_state_corrected):
                r_spock[itime_ens, iens, 2] = read_file_eci_z[nb_header_ens_eci_z+itime_ens].split()[2+iens]
        file_eci_z.close()




    if 'eci_v' in ensembles_to_output: #Read eci position of ensembles
        filename_ens_eci_vx = output_file_path_list[isc] + 'ensemble/ensemble_vx_eci_' + output_file_name_list[isc]
        file_eci_vx = open(filename_ens_eci_vx)
        read_file_eci_vx = file_eci_vx.readlines()

        nb_header_ens_eci_vx = 0
        while  read_file_eci_vx[nb_header_ens_eci_vx].split()[0] != '#START':
            nb_header_ens_eci_vx = nb_header_ens_eci_vx + 1
        nb_header_ens_eci_vx = nb_header_ens_eci_vx + 1
        nb_ensemble_ini_state_corrected = len(read_file_eci_vx[nb_header_ens_eci_vx].split())-2
        v_spock = np.zeros([n_spock, nb_ensemble_ini_state_corrected, 3])
        for itime_ens in range(n_spock):
            for iens in range(nb_ensemble_ini_state_corrected):
                v_spock[itime_ens, iens, 0] = read_file_eci_vx[nb_header_ens_eci_vx+itime_ens].split()[2+iens]
        file_eci_vx.close()

        filename_ens_eci_vy = output_file_path_list[isc] + 'ensemble/ensemble_vy_eci_' + output_file_name_list[isc]
        file_eci_vy = open(filename_ens_eci_vy)
        read_file_eci_vy = file_eci_vy.readlines()

        nb_header_ens_eci_vy = 0
        while  read_file_eci_vy[nb_header_ens_eci_vy].split()[0] != '#START':
            nb_header_ens_eci_vy = nb_header_ens_eci_vy + 1
        nb_header_ens_eci_vy = nb_header_ens_eci_vy + 1
        nb_ensemble_ini_state_corrected = len(read_file_eci_vy[nb_header_ens_eci_vy].split())-2

        for itime_ens in range(n_spock):
            for iens in range(nb_ensemble_ini_state_corrected):
                v_spock[itime_ens, iens, 1] = read_file_eci_vy[nb_header_ens_eci_vy+itime_ens].split()[2+iens]
        file_eci_vy.close()

        filename_ens_eci_vz = output_file_path_list[isc] + 'ensemble/ensemble_vz_eci_' + output_file_name_list[isc]
        file_eci_vz = open(filename_ens_eci_vz)
        read_file_eci_vz = file_eci_vz.readlines()

        nb_header_ens_eci_vz = 0
        while  read_file_eci_vz[nb_header_ens_eci_vz].split()[0] != '#START':
            nb_header_ens_eci_vz = nb_header_ens_eci_vz + 1
        nb_header_ens_eci_vz = nb_header_ens_eci_vz + 1
        nb_ensemble_ini_state_corrected = len(read_file_eci_vz[nb_header_ens_eci_vz].split())-2


        for itime_ens in range(n_spock):
            for iens in range(nb_ensemble_ini_state_corrected):
                v_spock[itime_ens, iens, 2] = read_file_eci_vz[nb_header_ens_eci_vz+itime_ens].split()[2+iens]
        file_eci_vz.close()


    # Compare SpOCK and data
    # Assumption: SpOCK was run with a 1s time step to avoid having to do interpolation here: the steps in SpOCK falls at the same time as the steps in data 
    ## Select the time where date_spock = date_obs 

    if it == 0:
        index_spock_same_date_as_obs = []
        iobs = 0
        while iobs < nb_obs:
            if date_obs[iobs] >= date_datetime_round_sec_spock[-1]:
                break
            else:
                index_spock_same_date_as_obs.append(np.where(date_datetime_round_sec_spock == date_obs[iobs])[0][0])
                nb_seconds_since_start.append( ( date_obs[iobs] - date_obs[0] ).total_seconds() )
            iobs = iobs + 1

    n = iobs #!!!!!!!!!! j-index_interval[iinter]
    date_spock_ok = date_spock[index_spock_same_date_as_obs]

    r_spock_ref_ok = np.zeros([n, 3])
    r_spock_ref_ok[:, 0] = r_spock_ref[index_spock_same_date_as_obs, 0]
    r_spock_ref_ok[:, 1] = r_spock_ref[index_spock_same_date_as_obs, 1]
    r_spock_ref_ok[:, 2] = r_spock_ref[index_spock_same_date_as_obs, 2]
    v_spock_ref_ok = np.zeros([n, 3])
    v_spock_ref_ok[:, 0] = v_spock_ref[index_spock_same_date_as_obs, 0]
    v_spock_ref_ok[:, 1] = v_spock_ref[index_spock_same_date_as_obs, 1]
    v_spock_ref_ok[:, 2] = v_spock_ref[index_spock_same_date_as_obs, 2]

    r_spock_ok = np.zeros([n, nb_ensemble_ini_state_corrected, 3])
    r_spock_ok[:, :, 0] = r_spock[index_spock_same_date_as_obs, :, 0]
    r_spock_ok[:, :, 1] = r_spock[index_spock_same_date_as_obs,:, 1]
    r_spock_ok[:,:, 2] = r_spock[index_spock_same_date_as_obs,:, 2]
    v_spock_ok = np.zeros([n, nb_ensemble_ini_state_corrected, 3])
    v_spock_ok[:, :, 0] = v_spock[index_spock_same_date_as_obs, :, 0]
    v_spock_ok[:, :, 1] = v_spock[index_spock_same_date_as_obs, :, 1]
    v_spock_ok[:, :, 2] = v_spock[index_spock_same_date_as_obs, :, 2]


    nb_steps = v_spock_ok.shape[0]
    if it == 0:
        index_obs = 0 # !!!!!!!!!!!!!!! index_interval[itime]
        for i in range(nb_steps):
            distance_ref.append( np.linalg.norm(r_obs[index_obs, :] - r_spock_ref_ok[i, :]) )
            index_obs = index_obs + 1


    distance_sub = []
    for iens in range(nb_ensemble_ini_state_corrected):
        index_obs = 0 # !!!!!!!!! index_interval[itime]
        distance_sub_ens = []
        for i in range(n):
            distance_sub_ens.append( np.linalg.norm(r_obs[index_obs, :] - r_spock_ok[i, iens, :]) )
            index_obs = index_obs + 1                  
        distance_sub.append( distance_sub_ens )


    distance.append( distance_sub )


    #nb_seconds_since_start = np.array(nb_seconds_since_start)
    mean_dist_itime_iens = np.zeros([nb_ensemble_ini_state_corrected]) # mean of the distance for
    # a given internval, aand a given ensemble. We first want to find the min of this variable over all ensembles


    for iens in range(nb_ensemble_ini_state_corrected):
        dist_itime_iens = np.array(distance[-1][iens])
        mean_dist_itime_iens[ iens] = np.mean(dist_itime_iens)
    which_ens_min_dist.append( np.where( mean_dist_itime_iens ==  np.nanmin(mean_dist_itime_iens) )[0][0] )
    min_mean_dist_itime_iens.append( np.nanmin(mean_dist_itime_iens) )

    main_input_filename_root = main_input_filename
    print it, min_mean_dist_itime_iens[-1]*1000.
    print format(r_spock_ok[0, which_ens_min_dist[-1], 0]*1000, '.14e'),\
        format(r_spock_ok[0, which_ens_min_dist[-1], 1]*1000, '.14e'),\
        format(r_spock_ok[0, which_ens_min_dist[-1], 2]*1000, '.14e')
    print format(v_spock_ok[0, which_ens_min_dist[-1], 0]*1000, '.14e'),\
        format(v_spock_ok[0, which_ens_min_dist[-1], 1]*1000, '.14e'),\
        format(v_spock_ok[0, which_ens_min_dist[-1], 2]*1000, '.14e')
    print  np.sqrt(sigma_x), np.sqrt( sigma_vx )
    if min_mean_dist_itime_iens[-1] < min_mean_dist_itime_iens[-2]*0.85:
        print 'acceptable'
        sigma_x_list_acceptable.append( np.sqrt( sigma_x ) )
        sigma_y_list_acceptable.append( np.sqrt( sigma_y ) )
        sigma_z_list_acceptable.append( np.sqrt( sigma_z ) )

        sigma_vx_list_acceptable.append( np.sqrt( sigma_vx ) )
        sigma_vy_list_acceptable.append( np.sqrt( sigma_vy ) )
        sigma_vz_list_acceptable.append( np.sqrt( sigma_vz ) )

        previous_r0 = r0
        previous_r1 = r1
        previous_r2 = r2
        previous_v0 = v0
        previous_v1 = v1
        previous_v2 = v2



nb_it = it + 1

# We now know the inital state that minimizes the distance to the observations.
# Initialize SpOCK with this ooptimized r0/v0 and propagate it for each interval using different values of the coeffieicnt on rho
# For each interval, record the rho coefficient that minimizes the distance to the observation
## The 4 lines below are redundant with before but are here for when  ic opy paste. 
date_start = date_obs_start
date_end = date_start + timedelta(seconds = interval_sec)
date_end_str = datetime.strftime(date_end, "%Y-%m-%dT%H:%M:%S")
date_start_str = datetime.strftime(date_start, "%Y-%m-%dT%H:%M:%S")

distance_rho = []
nb_seconds_since_start_rho = []
date_datetime_round_sec_spock_rho_ok = []
last_r0 = np.zeros([nb_rho])
last_r1 = np.zeros([nb_rho])
last_r2 = np.zeros([nb_rho])
last_v0 = np.zeros([nb_rho])
last_v1 = np.zeros([nb_rho])
last_v2 = np.zeros([nb_rho])
for iinter in range(nb_interval):
    nb_seconds_since_start_rho_inter = []
    distance_rho_interval = []
    for irho in range(nb_rho):
        print 'Interval', iinter, nb_interval-1, '|', 'rho', irho, nb_rho-1
        if iinter == 0:
            r0 = format(r_spock_ok[0, (int)(which_ens_min_dist[-1]), 0], '.14e')
            r1 = format(r_spock_ok[0, (int)(which_ens_min_dist[-1]), 1], '.14e')
            r2 = format(r_spock_ok[0, (int)(which_ens_min_dist[-1]), 2], '.14e')
            v0 = format(v_spock_ok[0, (int)(which_ens_min_dist[-1]), 0], '.14e')
            v1 = format(v_spock_ok[0, (int)(which_ens_min_dist[-1]), 1], '.14e')
            v2 = format(v_spock_ok[0, (int)(which_ens_min_dist[-1]), 2], '.14e')
        else:
            r0 = format(last_r0[irho], '.14e')
            r1 = format(last_r1[irho], '.14e')
            r2 = format(last_r2[irho], '.14e')
            v0 = format(last_v0[irho], '.14e')
            v1 = format(last_v1[irho], '.14e')
            v2 = format(last_v2[irho], '.14e')

        print r0, r1, r2
        print v0, v1, v2

        rho_mod = rho_mod_arr[irho]
        main_input_filename = date_start_str.replace(":","_") + '_' + date_end_str.replace(":","_")+ '_rhomod_' + format(rho_mod, ".2f").replace(".", "_") + '.txt'
        spock_main_input( # need to be in spokc/srcPython to run this script   
            main_input_filename,
            # for TIME section
               date_start_str, # first interval: same as during the r/v optimization. subsequent intervals: last date of previous interval
            date_end_str,
            dt,
            # for SPACECRAFT section
                    1,
            '0',
            29,
            "cygnss_geometry_2016_acco08.txt", #"cygnss_geometry_2016_smaller_solar_radiation_coeff.txt", #"cygnss_geometry_2016.txt",#"cygnss_geometry_2016_acco09.txt",
            # for ORBIT section
                ['state_eci','(' + r0 + '; ' + r1 + '; ' + r2 + ') (' + v0 + '; ' + v1 + '; ' + v2 + ')' ],
            # for FORCES section
              gravity_order, # !!!!!!!!!!! put back 20
            "drag solar_pressure sun_gravity moon_gravity", # !!!!!!!!!!!!! put back to "drag sun_gravity moon_gravity"
            'dynamic',
            # for OUTPUT section
                    "out",
            dt_output, 
            # for ATTITUDE section
            obs_att_filename,
            # for GROUND_STATIONS section
                    "0",
            # for SPICE section
                    spice_path,
            # FOR #DENSITY_MOD section
                    rho_mod
        )

        ## Run SpOCK
        #if iinter == nb_interval - 1:# !!!!!!! remove this if condition
        #os.system(path_mpirun + ' -np 1 spock ' + main_input_filename)
        ## save position and velocity
        #os.system("python state_dev.py ./ " + main_input_filename + " save position velocity")


        # Read the position and velocity predicted by SpOCK
        isc = 0
        var_in, var_in_order = read_input_file(main_input_filename)

        output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
        output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
        var_to_read = ["position", "velocity"]
        var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
        date_spock = np.array(var_out[find_in_read_input_order_variables(var_out_order, 'date')])
        date_datetime_round_sec_spock = np.array(var_out[find_in_read_input_order_variables(var_out_order, 'date_datetime_round_sec')])
        r_spock = var_out[find_in_read_input_order_variables(var_out_order, 'position')]
        v_spock = var_out[find_in_read_input_order_variables(var_out_order, 'velocity')]
        n_spock = len(date_spock)

        if irho == 0:
            index_spock_same_date_as_obs_rho = []
            if iinter == 0: # for the next interval, start date_obs[iobs] at the last observation of the previous interval
                iobs = 0
            while iobs < nb_obs:
                if date_obs[iobs] > date_datetime_round_sec_spock[-1]:
                    break
                else:
                    if len(index_spock_same_date_as_obs_rho) == 0:
                        first_obs = iobs
                    index_spock_same_date_as_obs_rho.append(np.where(date_datetime_round_sec_spock == date_obs[iobs])[0][0])
                    nb_seconds_since_start_rho_inter.append( ( date_obs[iobs] - date_obs[0] ).total_seconds() )
                iobs = iobs + 1
            nb_seconds_since_start_rho.append(nb_seconds_since_start_rho_inter)

            n = len(index_spock_same_date_as_obs_rho) #!!!!!!!!!! j-index_interval[iinter]
            date_datetime_round_sec_spock_rho_ok.append(date_datetime_round_sec_spock[index_spock_same_date_as_obs_rho])
            

        # Compare SpOCK and data
        r_spock_ok_rho = np.zeros([n, 3])
        r_spock_ok_rho[:, 0] = r_spock[index_spock_same_date_as_obs_rho, 0]
        r_spock_ok_rho[:, 1] = r_spock[index_spock_same_date_as_obs_rho, 1]
        r_spock_ok_rho[:, 2] = r_spock[index_spock_same_date_as_obs_rho, 2]
        v_spock_ok_rho = np.zeros([n, 3])
        v_spock_ok_rho[:, 0] = v_spock[index_spock_same_date_as_obs_rho, 0]
        v_spock_ok_rho[:, 1] = v_spock[index_spock_same_date_as_obs_rho, 1]
        v_spock_ok_rho[:, 2] = v_spock[index_spock_same_date_as_obs_rho, 2]

#        if rho_mod_arr[irho] == 1:
        last_r0[irho] = r_spock_ok_rho[-1, 0]
        last_r1[irho] = r_spock_ok_rho[-1, 1]
        last_r2[irho] = r_spock_ok_rho[-1, 2]
        last_v0[irho] = v_spock_ok_rho[-1, 0]
        last_v1[irho] = v_spock_ok_rho[-1, 1]
        last_v2[irho] = v_spock_ok_rho[-1, 2]
        print format(last_r0[irho],".14e"),format(last_r1[irho],".14e"),format(last_r2[irho],".14e")
        print format(last_v0[irho],".14e"),format(last_v1[irho],".14e"),format(last_v2[irho],".14e")


        if irho == (nb_rho - 1):
            date_start = date_datetime_round_sec_spock_rho_ok[-1][-1] 
            date_start_str = datetime.strftime(date_start, "%Y-%m-%dT%H:%M:%S")
            date_end_str = datetime.strftime(date_start + timedelta(seconds = interval_sec), "%Y-%m-%dT%H:%M:%S")


        distance_rho_sub = []
        index_obs = first_obs # !!!!!!!!!index_interval[itime]
        for i in range(n):
            distance_rho_sub.append( np.linalg.norm(r_obs[index_obs, :] - r_spock_ok_rho[i, :]) )
            index_obs = index_obs + 1
        distance_rho_interval.append( distance_rho_sub )
    distance_rho.append( distance_rho_interval )

    if iinter == 0:
        #nb_seconds_since_start = np.array(nb_seconds_since_start)
        mean_dist_itime_irho = np.zeros([nb_interval, nb_rho]) # min of the distance for a given internval and a given rho_mod. This is what we want to mnimize using the optimum rho_mod
        which_rho_min_dist = np.zeros([nb_interval])
        min_mean_dist_itime_irho =  np.zeros([nb_interval])
    for irho in range(nb_rho):
        dist_itime_irho = np.array(distance_rho[-1][irho])
        mean_dist_itime_irho[iinter, irho] = np.mean(dist_itime_irho)
    which_rho_min_dist[iinter] = rho_mod_arr[np.where( mean_dist_itime_irho[iinter, :] ==  np.nanmin(mean_dist_itime_irho[iinter, :]) )[0][0]]
    min_mean_dist_itime_irho[iinter] = np.nanmin(mean_dist_itime_irho[iinter, :])




distance_rho_concantenate = []
nb_seconds_since_start_rho_concatenate = []
for irho in range(nb_rho):
    distance_rho_concantenate_irho = []
    for iinter in range(nb_interval):
        distance_rho_concantenate_irho = distance_rho_concantenate_irho + distance_rho[iinter][irho]
        if irho == 0:
            nb_seconds_since_start_rho_concatenate = nb_seconds_since_start_rho_concatenate + nb_seconds_since_start_rho[iinter]
    distance_rho_concantenate.append(distance_rho_concantenate_irho)



################### FIGURES ###################
height_fig = 11
ratio_fig_size = 4./3
fontsize_plot = 20


#######
fig_title = ''#'Distance between SpOCK and data for different density coefficient' #'Distance with respect to rho = 0.7'#'Distance between SpOCK and data for different density coefficient'
y_label = 'Distance (m)'
x_label = 'Real time'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

x_axis = nb_seconds_since_start_rho_concatenate
if nb_rho > 1:
    alpha_arr = np.arange(0.2,1+0.2/nb_rho,(1-0.2)/(nb_rho-1))
else: 
    alpha_arr = [1]
for irho in range(nb_rho):
    if alpha_arr[irho] >1:
        alpha_arr[irho] = 1

    ax.plot(x_axis, np.array(distance_rho_concantenate[irho])*1000, linewidth = 2, color = 'b', alpha = alpha_arr[irho])
    #ax.plot(x_axis, (np.array(distance_rho[iinter][irho]) - np.array(distance_rho[iinter][0]))*1000, linewidth = 2, color = 'b', alpha = alpha_arr[irho])
    #ax.text(x_axis[-1], (np.array(distance_rho[iinter][irho]) - np.array(distance_rho[iinter][0]))[-1]*1000, format(rho_mod_arr[irho], ".1f"), horizontalalignment = 'left', fontsize = fontsize_plot, weight = 'bold', color = 'b', alpha = alpha_arr[irho], verticalalignment = 'center')
    ax.text(x_axis[-1], distance_rho_concantenate[irho][-1]*1000, format(rho_mod_arr[irho], ".1f"), horizontalalignment = 'left', fontsize = fontsize_plot, weight = 'bold', color = 'b', alpha = alpha_arr[irho], verticalalignment = 'center')
    #ax.text(x_axis[-1], distance_rho[iinter][irho][-1], str(rho_mod_arr[irho]), horizontalalignment = 'left', fontsize = fontsize_plot, weight = 'bold', color = 'b', alpha = alpha_arr[irho], verticalalignment = 'center')
    print irho
# x axis label is in real time
nb_seconds_in_simu = nb_seconds_since_start_rho_concatenate[-1] - nb_seconds_since_start_rho_concatenate[0]
start_xaxis_label = nb_seconds_since_start_rho_concatenate[0]
date_ref = date_datetime_round_sec_spock_rho_ok[0][0]
nb_ticks_xlabel = 10
dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
date_list_str = []
date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
for i in range(len(xticks)):
    if dt_xlabel > nb_ticks_xlabel*24*3600:
        date_list_str.append( str(date_list[i])[5:10] )
    else:
        date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
        ax.xaxis.set_ticks(xticks)
        ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
        ax.margins(0,0); ax.set_xlim([min(xticks), max(xticks)])
#        ax.set_xlim([ax.get_xlim()[0], most_recent_tle_among_all_sc])

legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="", fontsize = fontsize_plot)
#legend.get_title().set_fontsize(str(fontsize_plot))


fig_save_name = 'rho_distance_ens_to_observations_' + main_input_filename_root.replace("txt", "pdf")
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


raise Exception



# Distance between SpOCK ensembles and ODTK
height_fig = 11
ratio_fig_size = 4./3
fontsize_plot = 20

fig_title = 'Distance between SpOCK ensembles and ODTK'
y_label = 'Distance (m)'
x_label = 'Real time'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.96,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold                                                                                                                                       
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                                                                                          
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold                                                                                                                                       
x_axis = nb_seconds_since_start


for it in range(nb_it):
    if it == 0:
        for iens in range(nb_ensemble_ini_state_corrected):
            if iens == 0:
                ax.plot(x_axis, np.array(distance[it][iens])*1000, linewidth = 2, color = 'b', alpha = 0.15, label= 'SpOCK ensemble')
            else:
                ax.plot(x_axis, np.array(distance[it][iens])*1000, linewidth = 2, color = 'b', alpha = 0.15)
            # distance of SpOCK reference sc to ODT
        ax.plot(x_axis, np.array(distance_ref)*1000, linewidth = 4, color = 'r', label = 'SpOCK from raw observations')

    # min mean distance
    ax.plot(x_axis, np.array(distance[it][(int)(which_ens_min_dist[it])])*1000, linewidth = 5, color = 'b', label = 'Iteration ' + str(it))




# x axis label is in real time
nb_seconds_in_simu = nb_seconds_since_start[-1] - nb_seconds_since_start[0]
start_xaxis_label = nb_seconds_since_start[0]
date_ref = date_obs[0]
nb_ticks_xlabel = 10
dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
date_list_str = []
date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
for i in range(len(xticks)):
    if dt_xlabel > nb_ticks_xlabel*24*3600:
        date_list_str.append( str(date_list[i])[5:10] )
    else:
        date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
        ax.xaxis.set_ticks(xticks)
        ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
        ax.margins(0,0); ax.set_xlim([min(xticks), max(xticks)])
#        ax.set_xlim([ax.get_xlim()[0], most_recent_tle_among_all_sc])

legend = ax.legend(loc='upper left', numpoints = 1,  title="", fontsize = fontsize_plot)
#legend.get_title().set_fontsize(str(fontsize_plot))


fig_save_name = 'rv_distance_ens_to_observations_' + main_input_filename_root.replace("txt", "pdf")
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  










# Distribution x eci
bin_width = sigma_x/5. # in m
fig_title = 'X ECI distribution at initial time (bin size ' + str(bin_width) + ' m, ' + str(nb_ensemble_ini_state_corrected) + ' ensembles)'
y_label = '# ensembles in bin'
x_label = 'X (km)'

fig_x = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig_x.suptitle(fig_title, y = 0.958,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax_x = fig_x.add_subplot(gs[0, 0])

ax_x.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax_x.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax_x.spines.itervalues()] # change the width of the frame of the figure
ax_x.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
bin_width = bin_width / 1000. # m to km

index_obs_start = 0#index_interval[itime_start]
index_in_spock_ok = 0
index_in_obs = index_obs_start + index_in_spock_ok


bins_arr = np.arange(min(r_spock_ok[index_in_spock_ok, :, 0]), max(r_spock_ok[index_in_spock_ok, :, 0]) + bin_width, bin_width)
n, bins, patches = ax_x.hist(r_spock_ok[index_in_spock_ok, :, 0], bins_arr,  histtype='stepfilled', alpha = 1, color = 'cornflowerblue',label = 'SpOCK ensembles\nstd dev: ' + format(np.std(r_spock_ok[index_in_spock_ok, :, 0])*1000, ".2f") + ' m') 
# Add Observations position
ax_x.plot([r_obs[index_in_obs, 0], r_obs[index_in_obs, 0]],[0,np.nanmax(n)], linewidth = 6, color = 'b', label = 'Observations', linestyle = 'dotted')

legend = ax_x.legend(loc='top right', numpoints = 1,  title="", fontsize = fontsize_plot)

fig_save_name = 'x_eci_' + main_input_filename.replace(".txt", ".pdf")
fig_x.savefig(fig_save_name, facecolor=fig_x.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# Distribution y eci
bin_width = sigma_y/5. # in m
fig_title = 'Y ECI distribution at initial time (bin size ' + str(bin_width) + ' m, ' + str(nb_ensemble_ini_state_corrected) + ' ensembles)'
y_label = '# ensembles in bin'
x_label = 'Y (km)'

fig_y = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig_y.suptitle(fig_title, y = 0.958,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax_y = fig_y.add_subplot(gs[0, 0])

ax_y.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax_y.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax_y.spines.itervalues()] # change the width of the frame of the figure
ax_y.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
bin_width = bin_width / 1000. # m to km

index_obs_start = 0 # index_interval[itime_start]
index_in_spock_ok = 0
index_in_obs = index_obs_start + index_in_spock_ok


bins_arr = np.arange(min(r_spock_ok[index_in_spock_ok, :, 1]), max(r_spock_ok[index_in_spock_ok, :, 1]) + bin_width, bin_width)
n, bins, patches = ax_y.hist(r_spock_ok[index_in_spock_ok, :, 1], bins_arr,  histtype='stepfilled', alpha = 1, color = 'cornflowerblue',label = 'SpOCK ensembles\nstd dev: ' + format(np.std(r_spock_ok[index_in_spock_ok, :, 1])*1000, ".2f") + ' m') 
# Add Observations position
ax_y.plot([r_obs[index_in_obs, 1], r_obs[index_in_obs, 1]],[0,np.nanmax(n)], linewidth = 6, color = 'b', label = 'Observations', linestyle = 'dotted')

legend = ax_y.legend(loc='top right', numpoints = 1,  title="", fontsize = fontsize_plot)

fig_save_name = 'y_eci_' + main_input_filename.replace(".txt", ".pdf")
fig_y.savefig(fig_save_name, facecolor=fig_y.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# Distribution z eci
bin_width = sigma_z/5. # in m
fig_title = 'Z ECI distribution at initial time (bin size ' + str(bin_width) + ' m, ' + str(nb_ensemble_ini_state_corrected) + ' ensembles)'
y_label = '# ensembles in bin'
x_label = 'Z (km)'

fig_z = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig_z.suptitle(fig_title, y = 0.958,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax_z = fig_z.add_subplot(gs[0, 0])

ax_z.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax_z.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax_z.spines.itervalues()] # change the width of the frame of the figure
ax_z.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
bin_width = bin_width / 1000. # m to km

index_obs_start = 0 # index_interval[itime_start]
index_in_spock_ok = 0
index_in_obs = index_obs_start + index_in_spock_ok


bins_arr = np.arange(min(r_spock_ok[index_in_spock_ok, :, 2]), max(r_spock_ok[index_in_spock_ok, :, 2]) + bin_width, bin_width)
n, bins, patches = ax_z.hist(r_spock_ok[index_in_spock_ok, :, 2], bins_arr,  histtype='stepfilled', alpha = 1, color = 'cornflowerblue',label = 'SpOCK ensembles\nstd dev: ' + format(np.std(r_spock_ok[index_in_spock_ok, :, 2])*1000, ".2f") + ' m') 
# Add Observations position
ax_z.plot([r_obs[index_in_obs, 2], r_obs[index_in_obs, 2]],[0,np.nanmax(n)], linewidth = 6, color = 'b', label = 'Observations', linestyle = 'dotted')

legend = ax_z.legend(loc='top right', numpoints = 1,  title="", fontsize = fontsize_plot)

fig_save_name = 'z_eci_' + main_input_filename.replace(".txt", ".pdf")
fig_z.savefig(fig_save_name, facecolor=fig_z.get_facecolor(), edgecolor='none', bbox_inches='tight')  



# Distribution vx eci
bin_width = sigma_vx/5. # in m
fig_vtitle = 'Vx ECI distribution at initial time (bin size ' + str(bin_width) + ' m, ' + str(nb_ensemble_ini_state_corrected) + ' ensembles)'
y_label = '# ensembles in bin'
x_label = 'X (km)'

fig_vx = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig_vx.suptitle(fig_vtitle, y = 0.958,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax_vx = fig_vx.add_subplot(gs[0, 0])

ax_vx.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax_vx.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax_vx.spines.itervalues()] # change the width of the frame of the figure
ax_vx.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
bin_width = bin_width / 1000. # m to km

index_obs_start = 0 # index_interval[itime_start]
index_in_spock_ok = 0
index_in_obs = index_obs_start + index_in_spock_ok


bins_arr = np.arange(min(v_spock_ok[index_in_spock_ok, :, 0]), max(v_spock_ok[index_in_spock_ok, :, 0]) + bin_width, bin_width)
n, bins, patches = ax_vx.hist(v_spock_ok[index_in_spock_ok, :, 0], bins_arr,  histtype='stepfilled', alpha = 1, color = 'cornflowerblue',label = 'SpOCK ensembles\nstd dev: ' + format(np.std(v_spock_ok[index_in_spock_ok, :, 0])*1000, ".2f") + ' m/s') 
# Add Observations position
ax_vx.plot([v_obs[index_in_obs, 0], v_obs[index_in_obs, 0]],[0,np.nanmax(n)], linewidth = 6, color = 'b', label = 'Observations', linestyle = 'dotted')

legend = ax_vx.legend(loc='top right', numpoints = 1,  title="", fontsize = fontsize_plot)

fig_vsave_name = 'vx_eci_' + main_input_filename.replace(".txt", ".pdf")
fig_vx.savefig(fig_vsave_name, facecolor=fig_vx.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# Distribution vy eci
bin_width = sigma_vy/5. # in m
fig_vtitle = 'Vy ECI distribution at initial time (bin size ' + str(bin_width) + ' m, ' + str(nb_ensemble_ini_state_corrected) + ' ensembles)'
y_label = '# ensembles in bin'
x_label = 'Y (km)'

fig_vy = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig_vy.suptitle(fig_vtitle, y = 0.958,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax_vy = fig_vy.add_subplot(gs[0, 0])

ax_vy.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax_vy.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax_vy.spines.itervalues()] # change the width of the frame of the figure
ax_vy.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
bin_width = bin_width / 1000. # m to km

index_obs_start = 0 # index_interval[itime_start]
index_in_spock_ok = 0
index_in_obs = index_obs_start + index_in_spock_ok


bins_arr = np.arange(min(v_spock_ok[index_in_spock_ok, :, 1]), max(v_spock_ok[index_in_spock_ok, :, 1]) + bin_width, bin_width)
n, bins, patches = ax_vy.hist(v_spock_ok[index_in_spock_ok, :, 1], bins_arr,  histtype='stepfilled', alpha = 1, color = 'cornflowerblue',label = 'SpOCK ensembles\nstd dev: ' + format(np.std(v_spock_ok[index_in_spock_ok, :, 1])*1000, ".2f") + ' m/s') 
# Add Observations position
ax_vy.plot([v_obs[index_in_obs, 1], v_obs[index_in_obs, 1]],[0,np.nanmax(n)], linewidth = 6, color = 'b', label = 'Observations', linestyle = 'dotted')

legend = ax_vy.legend(loc='top right', numpoints = 1,  title="", fontsize = fontsize_plot)

fig_vsave_name = 'vy_eci_' + main_input_filename.replace(".txt", ".pdf")
fig_vy.savefig(fig_vsave_name, facecolor=fig_vy.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# Distribution vz eci
bin_width = sigma_vz/5. # in m
fig_vtitle = 'Vz ECI distribution at initial time (bin size ' + str(bin_width) + ' m, ' + str(nb_ensemble_ini_state_corrected) + ' ensembles)'
y_label = '# ensembles in bin'
x_label = 'Z (km)'

fig_vz = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig_vz.suptitle(fig_vtitle, y = 0.958,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax_vz = fig_vz.add_subplot(gs[0, 0])

ax_vz.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax_vz.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax_vz.spines.itervalues()] # change the width of the frame of the figure
ax_vz.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
bin_width = bin_width / 1000. # m to km

index_obs_start = 0 # index_interval[itime_start]
index_in_spock_ok = 0
index_in_obs = index_obs_start + index_in_spock_ok


bins_arr = np.arange(min(v_spock_ok[index_in_spock_ok, :, 2]), max(v_spock_ok[index_in_spock_ok, :, 2]) + bin_width, bin_width)
n, bins, patches = ax_vz.hist(v_spock_ok[index_in_spock_ok, :, 2], bins_arr,  histtype='stepfilled', alpha = 1, color = 'cornflowerblue',label = 'SpOCK ensembles\nstd dev: ' + format(np.std(v_spock_ok[index_in_spock_ok, :, 2])*1000, ".2f") + ' m/s') 
# Add Observations position
ax_vz.plot([v_obs[index_in_obs, 2], v_obs[index_in_obs, 2]],[0,np.nanmax(n)], linewidth = 6, color = 'b', label = 'Observations', linestyle = 'dotted')

legend = ax_vz.legend(loc='top right', numpoints = 1,  title="", fontsize = fontsize_plot)

fig_vsave_name = 'vz_eci_' + main_input_filename.replace(".txt", ".pdf")
fig_vz.savefig(fig_vsave_name, facecolor=fig_vz.get_facecolor(), edgecolor='none', bbox_inches='tight')  




