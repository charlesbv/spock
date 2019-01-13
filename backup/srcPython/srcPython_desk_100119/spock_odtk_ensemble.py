# This script runs SpOCK with ensembles on the initial state (r, v ECI). The initial mean r,v is taken from the ODTK simulations by Kyle Nave.
# ASSUMPTIONS
#- see section "PARAMETERS TO SET UP BEFORE RUNNIG THIS SCRIPT"
#- run SpOCK with a 1s time step 

import numpy as np
import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/kalman/spock_development_new_structure_kalman_dev/srcPython")
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


#PARAMETERS TO SET UP BEFORE RUNNIG THIS SCRIPT
cygfm = 5 # which CYGFM to look at
path_mpirun = '/opt/local/bin/mpirun-openmpi-gcc49' #'/usr/local/bin/mpirun'# '/opt/local/bin/mpirun-openmpi-gcc49'
interval = 1.5 # interval of time to compare the two trajectories (data and SpOCK). In hours
rho_mod_min = 0.9 # min rho_mod
rho_mod_max = 0.9 # max rho_mod
drho_mod = 0.1 # step in rho_mod -> rho_mod varies from rho_mod_min to rho_mod_max by values of drho_mod
nb_ensemble_ini_state = 40
sigma_x = 1 # in m
sigma_y = 1 # in m
sigma_z = 1 # in m
sigma_vx = 0. # in m/s
sigma_vy = 0.0 # in m/s
sigma_vz = 0.0 # in m/s

# end of PARAMETERS TO SET UP BEFORE RUNNIG THIS SCRIPT

# ALGORITHM
cygfm_to_ccsds = ['F7','F9','2B','2C','2F','36','37','49']
interval_sec = interval * 3600.



# Read data ephemerides
cygccsds = cygfm_to_ccsds[cygfm - 1]
filename = cygccsds + "_20170823_134632_STKdefpred_v001.e"
filename = '../cyg_data/' + filename 
file = open(filename)
read_file = file.readlines()
iline = 0
found_epoch = 0
while (found_epoch == 0):
    if ('YYDDD' in read_file[iline]):
        epochyyddd = np.float(read_file[iline].split(":")[1])
        found_epoch = 1
    iline = iline + 1

found_hdr = 0
while (found_hdr == 0):
    if (len(read_file[iline].split()) > 0):
        if read_file[iline].split()[0] == "EphemerisTimePosVel":
            found_hdr = 1
    iline = iline + 1
nb_header = iline + 1

n = len(read_file) - nb_header
date_data = []
date_data_raw = []
r_data = []
v_data = []

epochyyddd_no_decimal = (int)(epochyyddd)
epochyyddd_decimal_only = epochyyddd - epochyyddd_no_decimal
epochyyddd_date = datetime.strptime( str(epochyyddd_no_decimal), "%y%j" ) + timedelta( hours = epochyyddd_decimal_only*24 )



iline = 0
new_date_start = []
nb_seconds_between_epochyyddd_and_new_date_start = []
new_r_data = []
new_v_data = []
index_interval = []
date_stop = "2017-08-21T14:08:00"
date_stop = datetime.strptime(date_stop, "%Y-%m-%dT%H:%M:%S")
while (iline < n):
    if iline > 0:
        if datetime.strptime(date_data[-1], "%Y/%m/%d %H:%M:%S") > date_stop:
            break
    if len(read_file[iline+nb_header].split()) == 0:
        break
    date_data_temp = np.float( read_file[iline+nb_header].split()[0] )
    if iline == 0:
        date_data.append( datetime.strftime( epochyyddd_date + timedelta( seconds = date_data_temp ) , "%Y/%m/%d %H:%M:%S") )
        date_data_raw.append( read_file[iline+nb_header].split()[0] )
        r_data.append( [np.float(read_file[iline+nb_header].split()[1]), np.float(read_file[iline+nb_header].split()[2]), np.float(read_file[iline+nb_header].split()[3])] )
        v_data.append( [np.float(read_file[iline+nb_header].split()[4]), np.float(read_file[iline+nb_header].split()[5]), np.float(read_file[iline+nb_header].split()[6])] )
        index_interval.append(iline)
    else:
        index_interval.append(iline-1)
    new_date_start.append( date_data[-1] )
    new_r_data.append( r_data[-1] )
    new_v_data.append( v_data[-1] )
    nb_seconds_between_epochyyddd_and_first_new_date_start = date_data_temp
    time_elapsed = 0
    while time_elapsed < interval_sec:
        if len(read_file[iline+nb_header].split()) == 0:
            break
        date_data_temp = np.float( read_file[iline+nb_header].split()[0] )
        if iline > 0:
            date_data.append( datetime.strftime( epochyyddd_date + timedelta( seconds = date_data_temp ) , "%Y/%m/%d %H:%M:%S") )
            date_data_raw.append( read_file[iline+nb_header].split()[0]  )
            r_data.append( [np.float(read_file[iline+nb_header].split()[1]), np.float(read_file[iline+nb_header].split()[2]), np.float(read_file[iline+nb_header].split()[3])] )
            v_data.append( [np.float(read_file[iline+nb_header].split()[4]), np.float(read_file[iline+nb_header].split()[5]), np.float(read_file[iline+nb_header].split()[6])] )

        time_elapsed = date_data_temp - nb_seconds_between_epochyyddd_and_first_new_date_start
        iline = iline + 1

new_r_data = np.array(new_r_data)/1000.
new_v_data = np.array(new_v_data)/1000.

r_data = np.array(r_data)/1000.
v_data = np.array(v_data)/1000.

nb_interval = len(new_date_start)
rho_mod_arr = np.arange(rho_mod_min, rho_mod_max+drho_mod, drho_mod)
nb_rho = len(rho_mod_arr)
distance = []
distance_ref = []
nb_seconds_since_start = []
date_start_simu = "2017-08-21T11:30:00"# run SpOCK only starting at this date
date_start_simu = datetime.strptime(date_start_simu, "%Y-%m-%dT%H:%M:%S")
itime_count = -1
for itime in range(nb_interval-1):
    if datetime.strptime(new_date_start[itime],  "%Y/%m/%d %H:%M:%S") >= date_start_simu:
        print itime, nb_interval-1
        
        itime_count = itime_count + 1
        if itime_count == 0:
            itime_start = itime
        distance_interval = []
        distance_ref_interval = []
        # for each new interval, run SpOCK with different values of rho_mod: the initial date is the start date of the interval, initial r/v is the r/v at this date
        ## Create SpOCK main input file: same epoch and initial r/v
        date_start = new_date_start[itime].replace("/", "-").replace(" ", "T")+'.000'
        date_end = new_date_start[itime+1].replace("/", "-").replace(" ", "T")+'.000'

        nb_seconds_between_epochyyddd_and_new_date_start.append( ( datetime.strptime( new_date_start[itime], "%Y/%m/%d %H:%M:%S") - epochyyddd_date ).total_seconds() ) 
        dt  = 1
        dt_output = 1
        gravity_order = 50

        for irho in range(nb_rho):
            rho_mod = rho_mod_arr[irho]
            main_input_filename = 'FM0' + str(cygfm) + '_' + date_start.replace(":","_") + '_' + date_end.replace(":","_") + '_rhomod_' + str(rho_mod).replace(".", "") + '_sigmas_x_' + str(sigma_x)  + '_y_' + str(sigma_y)  + '_z_' + str(sigma_z)  + '_vx_' + str(sigma_vx)  + '_vy_' + str(sigma_vy) +  '_vz_' + str(sigma_vz)  + '_nbens_' + str(nb_ensemble_ini_state) + '.txt'

            r0 = format(new_r_data[itime, 0]*1000, '.14e')
            r1 = format(new_r_data[itime, 1]*1000, '.14e')
            r2 = format(new_r_data[itime, 2]*1000, '.14e')
            v0 = format(new_v_data[itime, 0]*1000, '.14e')
            v1 = format(new_v_data[itime, 1]*1000, '.14e')
            v2 = format(new_v_data[itime, 2]*1000, '.14e')

            # SpOCK inital state uncetainty file
            filename_ini_state = 'FM0' + str(cygfm) + '_' + date_start.replace(":","_") + '_ini_state.txt' 
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
                    date_start,
                date_end,
                dt,
                # for SPACECRAFT section
                        1,
                '0',
                28,
                "../cygnss_geometry_2016_acco08.txt", #"cygnss_geometry_2016_smaller_solar_radiation_coeff.txt", #"cygnss_geometry_2016.txt",#"cygnss_geometry_2016_acco09.txt",
                # for ORBIT section
                    ['collision', filename_ini_state ],
                # for FORCES section
                    gravity_order, # !!!!!!!!!!! put back 20
                "drag sun_gravity moon_gravity", # !!!!!!!!!!!!! put back to "drag sun_gravity moon_gravity"
                ['../2017Q2Q3_DSD_converted_for_spock.txt','../2017Q2Q3_DGD_converted_for_spock.txt'],#['f107_for_81daverage_20170621_to_20170822.txt', 'ap_20170730_to_20170822.txt'],
                # for OUTPUT section
                        "~/eclipse_ensemble_ini_state/out",
                dt_output, 
                # for ATTITUDE section
                "nadir",
                # for GROUND_STATIONS section
                        "0",
                # for SPICE section
                        "/Users/cbv/cspice/data",
                # FOR #DENSITY_MOD section
                        1
            )
            # add in SpOCK main input file the section #OUTPUT_ENSEMBLES
            file_spock = open(main_input_filename, "a")
            print >> file_spock, "#OUTPUT_ENSEMBLES\neci_r"
            file_spock.close()


            ## Run SpOCK
            os.system(path_mpirun + ' -np 8 spock_dev_parallel_kalman_9state ' + main_input_filename)

            ## concatenate proc files
            print "Concatenating processor files"
            os.system(path_mpirun + ' -np 8 python new_mpi_concatenate_proc.py ' + main_input_filename )

            # Read the position and velocity predicted by SpOCK
            isc = 0
            var_in, var_in_order = read_input_file(main_input_filename)

            output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
            output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
            var_to_read = ["position", "velocity"]
            var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
            date_spock = np.array(var_out[find_in_read_input_order_variables(var_out_order, 'date')])
            n_spock = len(date_spock)
            if ( ( itime_count == 0 ) & ( irho == 0 ) ):
                ensembles_to_output = var_in[find_in_read_input_order_variables(var_in_order, 'ensembles_to_output')];

            if 'eci_r' in ensembles_to_output: #Read eci position of ensembles
                filename_ens_eci_x = output_file_path_list[isc] + 'ensemble/ensemble_x_eci_' + output_file_name_list[isc]
                file_eci_x = open(filename_ens_eci_x)
                read_file_eci_x = file_eci_x.readlines()
                if ( ( itime_count == 0 ) & ( irho == 0 ) ):
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
                if ( ( itime_count == 0 ) & ( irho == 0 ) ):
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
                if ( ( itime_count == 0 ) & ( irho == 0 ) ):
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



            # Compare SpOCK and data
            # Assumption: SpOCK was run with a 1s time step to avoid having to do interpolation here: the steps in SpOCK falls at the same time as the steps in data 
            ## Select the time where date_spock = date_data 
            if irho == 0:
                index_spock_same_date_as_data = []
                j = index_interval[itime]
                nb_seconds_since_start_itime = []
                for i in range(n_spock):
                    if date_spock[i][:-7] == date_data[j]:
                        index_spock_same_date_as_data.append(i)
                        nb_seconds_since_start_itime.append( ( datetime.strptime(date_data[j],"%Y/%m/%d %H:%M:%S") - datetime.strptime(date_data[0],"%Y/%m/%d %H:%M:%S") ).total_seconds() )
                        j = j + 1
                n = j-index_interval[itime]
                nb_seconds_since_start.append( nb_seconds_since_start_itime )
                date_spock_ok = date_spock[index_spock_same_date_as_data]

            r_spock_ref_ok = np.zeros([n, 3])
            r_spock_ref_ok[:, 0] = r_spock_ref[index_spock_same_date_as_data, 0]
            r_spock_ref_ok[:, 1] = r_spock_ref[index_spock_same_date_as_data, 1]
            r_spock_ref_ok[:, 2] = r_spock_ref[index_spock_same_date_as_data, 2]
            v_spock_ref_ok = np.zeros([n, 3])
            v_spock_ref_ok[:, 0] = v_spock_ref[index_spock_same_date_as_data, 0]
            v_spock_ref_ok[:, 1] = v_spock_ref[index_spock_same_date_as_data, 1]
            v_spock_ref_ok[:, 2] = v_spock_ref[index_spock_same_date_as_data, 2]

            r_spock_ok = np.zeros([n, nb_ensemble_ini_state_corrected, 3])
            r_spock_ok[:, :, 0] = r_spock[index_spock_same_date_as_data, :, 0]
            r_spock_ok[:, :, 1] = r_spock[index_spock_same_date_as_data,:, 1]
            r_spock_ok[:,:, 2] = r_spock[index_spock_same_date_as_data,:, 2]
#             v_spock_ok = np.zeros([n, 3])
#             v_spock_ok[:, 0] = v_spock[index_spock_same_date_as_data, 0]
#             v_spock_ok[:, 1] = v_spock[index_spock_same_date_as_data, 1]
#             v_spock_ok[:, 2] = v_spock[index_spock_same_date_as_data, 2]

            distance_ref_sub = []
            index_data = index_interval[itime]
            for i in range(n):
                distance_ref_sub.append( np.linalg.norm(r_data[index_data, :] - r_spock_ref_ok[i, :]) )
                index_data = index_data + 1
            distance_ref_interval.append( distance_ref_sub )


            distance_sub = []
            for iens in range(nb_ensemble_ini_state_corrected):
                index_data = index_interval[itime]
                distance_sub_ens = []
                for i in range(n):
                    distance_sub_ens.append( np.linalg.norm(r_data[index_data, :] - r_spock_ok[i, iens, :]) )
                    index_data = index_data + 1                  
                distance_sub.append( distance_sub_ens )

            distance_interval.append( distance_sub )
        distance.append( distance_interval )
        distance_ref.append( distance_ref_interval )
#nb_seconds_since_start = np.array(nb_seconds_since_start)
mean_dist_itime_irho_iens = np.zeros([nb_interval-1-itime_start, nb_rho, nb_ensemble_ini_state_corrected]) # mean of the distance for a given internval, a given rho_mod, and a given ensemble. We first want to find the min of this variable over all ensembles
#mean_dist_itime_irho = np.zeros([nb_interval-1-itime_start, nb_rho]) # mean of the distance for a given internval and a given rho_mod. This is what we want to mnimize using the optimum rho_mod
which_ens_min_dist = np.zeros([nb_interval-1-itime_start, nb_rho])
min_mean_dist_itime_irho_iens =  np.zeros([nb_interval-1-itime_start, nb_rho])
nan_in_run = np.zeros([nb_interval-1-itime_start]) 
for itime in range(0, nb_interval-1-itime_start):
    for irho in range(nb_rho):
        for iens in range(nb_ensemble_ini_state_corrected):
            dist_itime_irho_iens = np.array(distance[itime][irho][iens])
            mean_dist_itime_irho_iens[itime, irho, iens] = np.mean(dist_itime_irho_iens)
        which_ens_min_dist[itime, irho] = np.where( mean_dist_itime_irho_iens[itime, irho, :] ==  np.min(mean_dist_itime_irho_iens[itime,irho, :]) )[0][0]
        min_mean_dist_itime_irho_iens[itime, irho] = np.min(mean_dist_itime_irho_iens[itime,irho,  :])



################### FIGURES ###################
height_fig = 11
ratio_fig_size = 4./3
fontsize_plot = 20


# Distance between SpOCK ensembles and ODTK
itime = 0
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
x_axis = nb_seconds_since_start[itime]
irho = 0

for iens in range(nb_ensemble_ini_state_corrected):
    if iens == 0:
        ax.plot(x_axis, np.array(distance[itime][irho][iens])*1000, linewidth = 2, color = 'b', alpha = 0.15, label= 'SpOCK ensemble')
    else:
        ax.plot(x_axis, np.array(distance[itime][irho][iens])*1000, linewidth = 2, color = 'b', alpha = 0.15)

# min mean distance
ax.plot(x_axis, np.array(distance[itime][irho][(int)(which_ens_min_dist[itime, irho])])*1000, linewidth = 5, color = 'b', label = 'Closest SpOCK ensemble to ODTK')

# distance of SpOCK reference sc to ODT
ax.plot(x_axis, np.array(distance_ref[itime][irho])*1000, linewidth = 4, color = 'r', label = 'SpOCK reference')


# x axis label is in real time                                                                                                                                                                              
nb_seconds_in_simu = nb_seconds_since_start[itime][-1] - nb_seconds_since_start[itime][0]
start_xaxis_label = nb_seconds_since_start[itime][0]
date_ref = datetime.strptime(date_data[index_interval[itime+itime_start]],"%Y/%m/%d %H:%M:%S")
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



# x axis label is in real time
nb_seconds_in_simu = nb_seconds_since_start[itime][-1] - nb_seconds_since_start[itime][0]
start_xaxis_label = nb_seconds_since_start[itime][0]
date_ref = datetime.strptime(date_data[index_interval[itime+itime_start]],"%Y/%m/%d %H:%M:%S")
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


fig_save_name = 'distance_ens_to_odtk_' + main_input_filename.replace("txt", "pdf")
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

index_data_start = index_interval[itime_start]
index_in_spock_ok = 0
index_in_data = index_data_start + index_in_spock_ok


bins_arr = np.arange(min(r_spock_ok[index_in_spock_ok, :, 0]), max(r_spock_ok[index_in_spock_ok, :, 0]) + bin_width, bin_width)
n, bins, patches = ax_x.hist(r_spock_ok[index_in_spock_ok, :, 0], bins_arr,  histtype='stepfilled', alpha = 1, color = 'cornflowerblue',label = 'SpOCK ensembles\nstd dev: ' + format(np.std(r_spock_ok[index_in_spock_ok, :, 0])*1000, ".2f") + ' m') 
# Add ODTK position
ax_x.plot([r_data[index_in_data, 0], r_data[index_in_data, 0]],[0,np.max(n)], linewidth = 6, color = 'b', label = 'ODTK', linestyle = 'dotted')

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

index_data_start = index_interval[itime_start]
index_in_spock_ok = 0
index_in_data = index_data_start + index_in_spock_ok


bins_arr = np.arange(min(r_spock_ok[index_in_spock_ok, :, 1]), max(r_spock_ok[index_in_spock_ok, :, 1]) + bin_width, bin_width)
n, bins, patches = ax_y.hist(r_spock_ok[index_in_spock_ok, :, 1], bins_arr,  histtype='stepfilled', alpha = 1, color = 'cornflowerblue',label = 'SpOCK ensembles\nstd dev: ' + format(np.std(r_spock_ok[index_in_spock_ok, :, 1])*1000, ".2f") + ' m') 
# Add ODTK position
ax_y.plot([r_data[index_in_data, 1], r_data[index_in_data, 1]],[0,np.max(n)], linewidth = 6, color = 'b', label = 'ODTK', linestyle = 'dotted')

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

index_data_start = index_interval[itime_start]
index_in_spock_ok = 0
index_in_data = index_data_start + index_in_spock_ok


bins_arr = np.arange(min(r_spock_ok[index_in_spock_ok, :, 2]), max(r_spock_ok[index_in_spock_ok, :, 2]) + bin_width, bin_width)
n, bins, patches = ax_z.hist(r_spock_ok[index_in_spock_ok, :, 2], bins_arr,  histtype='stepfilled', alpha = 1, color = 'cornflowerblue',label = 'SpOCK ensembles\nstd dev: ' + format(np.std(r_spock_ok[index_in_spock_ok, :, 2])*1000, ".2f") + ' m') 
# Add ODTK position
ax_z.plot([r_data[index_in_data, 2], r_data[index_in_data, 2]],[0,np.max(n)], linewidth = 6, color = 'b', label = 'ODTK', linestyle = 'dotted')

legend = ax_z.legend(loc='top right', numpoints = 1,  title="", fontsize = fontsize_plot)

fig_save_name = 'z_eci_' + main_input_filename.replace(".txt", ".pdf")
fig_z.savefig(fig_save_name, facecolor=fig_z.get_facecolor(), edgecolor='none', bbox_inches='tight')  


raise Exception
# Distribution vx eci
bin_width = sigma_vx/5. # in m
fig_vtitle = 'Vx ECI distribution at initial time (bin size ' + str(bin_width) + ' m, ' + str(nb_ensemble_ini_state_corrected) + ' ensembles)'
y_label = '# ensembles in bin'
x_label = 'X (km)'

fig_vx = plt.figure(num=None, figsize=(height_fig * ratio_fig_vsize, height_fig), dpi=80, facecolor='w', edgecolor='k')

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

index_data_start = index_interval[itime_start]
index_in_spock_ok = 0
index_in_data = index_data_start + index_in_spock_ok


bins_arr = np.arange(min(v_spock_ok[index_in_spock_ok, :, 0]), max(v_spock_ok[index_in_spock_ok, :, 0]) + bin_width, bin_width)
n, bins, patches = ax_vx.hist(v_spock_ok[index_in_spock_ok, :, 0], bins_arr,  histtype='stepfilled', alpha = 1, color = 'cornflowerblue',label = 'SpOCK ensembles\nstd dev: ' + format(np.std(v_spock_ok[index_in_spock_ok, :, 0])*1000, ".2f") + ' m') 
# Add ODTK position
ax_vx.plot([v_data[index_in_data, 0], v_data[index_in_data, 0]],[0,np.max(n)], linewidth = 6, color = 'b', label = 'ODTK', linestyle = 'dotted')

legend = ax_vx.legend(loc='top right', numpoints = 1,  title="", fontsize = fontsize_plot)

fig_vsave_name = 'vx_eci_' + main_input_filename.replace(".txt", ".pdf")
fig_vx.savefig(fig_vsave_name, facecolor=fig_vx.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# Distribution vy eci
bin_width = sigma_vy/5. # in m
fig_vtitle = 'Vy ECI distribution at initial time (bin size ' + str(bin_width) + ' m, ' + str(nb_ensemble_ini_state_corrected) + ' ensembles)'
y_label = '# ensembles in bin'
x_label = 'Y (km)'

fig_vy = plt.figure(num=None, figsize=(height_fig * ratio_fig_vsize, height_fig), dpi=80, facecolor='w', edgecolor='k')

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

index_data_start = index_interval[itime_start]
index_in_spock_ok = 0
index_in_data = index_data_start + index_in_spock_ok


bins_arr = np.arange(min(v_spock_ok[index_in_spock_ok, :, 1]), max(v_spock_ok[index_in_spock_ok, :, 1]) + bin_width, bin_width)
n, bins, patches = ax_vy.hist(v_spock_ok[index_in_spock_ok, :, 1], bins_arr,  histtype='stepfilled', alpha = 1, color = 'cornflowerblue',label = 'SpOCK ensembles\nstd dev: ' + format(np.std(v_spock_ok[index_in_spock_ok, :, 1])*1000, ".2f") + ' m') 
# Add ODTK position
ax_vy.plot([v_data[index_in_data, 1], v_data[index_in_data, 1]],[0,np.max(n)], linewidth = 6, color = 'b', label = 'ODTK', linestyle = 'dotted')

legend = ax_vy.legend(loc='top right', numpoints = 1,  title="", fontsize = fontsize_plot)

fig_vsave_name = 'vy_eci_' + main_input_filename.replace(".txt", ".pdf")
fig_vy.savefig(fig_vsave_name, facecolor=fig_vy.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# Distribution vz eci
bin_width = sigma_vz/5. # in m
fig_vtitle = 'Vz ECI distribution at initial time (bin size ' + str(bin_width) + ' m, ' + str(nb_ensemble_ini_state_corrected) + ' ensembles)'
y_label = '# ensembles in bin'
x_label = 'Z (km)'

fig_vz = plt.figure(num=None, figsize=(height_fig * ratio_fig_vsize, height_fig), dpi=80, facecolor='w', edgecolor='k')

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

index_data_start = index_interval[itime_start]
index_in_spock_ok = 0
index_in_data = index_data_start + index_in_spock_ok


bins_arr = np.arange(min(v_spock_ok[index_in_spock_ok, :, 2]), max(v_spock_ok[index_in_spock_ok, :, 2]) + bin_width, bin_width)
n, bins, patches = ax_vz.hist(v_spock_ok[index_in_spock_ok, :, 2], bins_arr,  histtype='stepfilled', alpha = 1, color = 'cornflowerblue',label = 'SpOCK ensembles\nstd dev: ' + format(np.std(v_spock_ok[index_in_spock_ok, :, 2])*1000, ".2f") + ' m') 
# Add ODTK position
ax_vz.plot([v_data[index_in_data, 2], v_data[index_in_data, 2]],[0,np.max(n)], linewidth = 6, color = 'b', label = 'ODTK', linestyle = 'dotted')

legend = ax_vz.legend(loc='top right', numpoints = 1,  title="", fontsize = fontsize_plot)

fig_vsave_name = 'vz_eci_' + main_input_filename.replace(".txt", ".pdf")
fig_vz.savefig(fig_vsave_name, facecolor=fig_vz.get_facecolor(), edgecolor='none', bbox_inches='tight')  



