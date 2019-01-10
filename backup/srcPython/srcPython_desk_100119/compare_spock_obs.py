# This script computes the along-track distance between SpOCK input file spock_in and an observation file obs-in. These two files are given as arguments
import sys
spock_in = sys.argv[1]#'try_interval18_0_iinter23_irho0.txt'
obs_in = 'HD_data/122021_eci.txt'#1220_eci.txt'

# PARAMETERS TO SET UP BEFORE RUNNIG THIS SCRIPT
isbig = 0 # if runnign script from Big
ispleiades = 0 # if runnign script from Pleaides

no_prop = 0 # set this variable to 1 to prevent creating SpOCK main input files and propagating them
interval = 18.0 # interval of time to compare the two trajectories (data and SpOCK). In hours
step_move_save = 23.#95./60
step_drho = 0.1 # the rho control will vary by this amount to find the optimum rho over an interval
kplist = [1.] # list of proportional gains for PID
kdlist = [1.] # list of derivative gains for PID
kilist = [0.000] # list of integral gains for PID
plot_or_not = 1
inter_start_algo = -1#2.0
# end of PARAMETERS TO SET UP BEFORE RUNNIG THIS SCRIPT

if isbig == 1 & ispleiades == 1:
    print "***! Choose to run on Pleiades or Big, but not both. The program will stop. !***"; raise Exception



import numpy as np

if isbig == 1:
    sys.path.append("/home/cbv/code/spock/srcPython")
    path_mpirun = '/usr/local/bin/mpirun'
    spice_path = '/raid4/cbv/cspice/data'
    nb_proc = 12
elif ispleiades == 1:
    sys.path.append("/home1/cbussy/Code/spock/srcPython")
    path_mpirun = 'mpiexec'
    spice_path = '/home1/cbussy/cspice/data'
    nb_proc = 0    

else:
    sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/srcPython")
    #sys.path.append("/home/cbv/Code_drive/Code/spock/srcPython")
    path_mpirun = 'mpirun'
    spice_path = '/Users/cbv/cspice/data'
    nb_proc = 4


import pickle
import os

from read_input_file import *
from read_output_file import *
from spock_main_input import *
#if ispleiades != 1:
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
if ((isbig != 1) & (ispleiades !=1)):
    #from convert_cygnss_obs_ecef_to_eci import *
    import ipdb
from eci_to_lvlh import *
#plt.ion()




# Read r/v of observations

# Read observation ECI r/v 
obs_rv_file = open(obs_in)
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
ecc_obs = np.zeros([nb_obs])

earth_mu = 398600.4418 # km^3/s^2
orbit_period = 95 * 60. # in seconds
for iobs in range(nb_obs):
    date_obs_str.append( read_obs_rv_file[iobs + nb_header].split()[0] )
    date_obs.append( datetime.strptime(date_obs_str[-1], "%Y-%m-%dT%H:%M:%S" ) )
    r_obs[iobs, 0] = np.float( read_obs_rv_file[iobs + nb_header].split()[1] ) 
    r_obs[iobs, 1] = np.float( read_obs_rv_file[iobs + nb_header].split()[2] ) 
    r_obs[iobs, 2] = np.float( read_obs_rv_file[iobs + nb_header].split()[3] ) 
    v_obs[iobs, 0] = np.float( read_obs_rv_file[iobs + nb_header].split()[4] ) 
    v_obs[iobs, 1] = np.float( read_obs_rv_file[iobs + nb_header].split()[5] ) 
    v_obs[iobs, 2] = np.float( read_obs_rv_file[iobs + nb_header].split()[6] ) 


    mu_earth = 398600.4418; # km^3/s^2 
    rrss = np.linalg.norm(r_obs[iobs, :])
    vrss = np.linalg.norm(v_obs[iobs, :])
    rdotv = np.dot(r_obs[iobs, :], v_obs[iobs, :])
    tempv1 = v_obs[iobs, :] * rdotv
    
    coeff = vrss*vrss - (mu_earth/rrss);
    tempv2 = r_obs[iobs, :] * coeff
    e_vector = tempv2 - tempv1
    coeff = 1.0/mu_earth;
    e_vector = e_vector * coeff
    ecc_obs[iobs] = np.linalg.norm(e_vector)


# Run SpOCK: initial r/v is given by observations + ensemble with std given by x_sigma, y_sigma, etc
date_obs_start_str = date_obs_str[0]
date_obs_start= datetime.strptime(date_obs_start_str, "%Y-%m-%dT%H:%M:%S")
date_obs_end_str = date_obs_str[-1]
date_obs_end= datetime.strptime(date_obs_end_str, "%Y-%m-%dT%H:%M:%S")
interval_sec = interval * 3600.
nb_interval = 1


print 'nb of intervals:', nb_interval
nb_seconds_since_start = []
distance = []

date_start = date_obs_start
date_end = date_start + timedelta(seconds = interval_sec)
date_end_str = datetime.strftime(date_end, "%Y-%m-%dT%H:%M:%S")
date_start_str = datetime.strftime(date_start, "%Y-%m-%dT%H:%M:%S")

index_obs_interval_start = 0 # !!!!!!!!!!! to change

## SpOCK main input file:
dt  = 1
dt_output = 60 # !!!!!!!!!used to be 1
gravity_order = 20 # !!!!!!!!!! put 50 (or 20)

date_start = date_obs_start#datetime.strptime("2017-12-18T06:00:00", "%Y-%m-%dT%H:%M:%S")#date_obs_start
date_end = date_start + timedelta(seconds = interval_sec)
date_end_str = datetime.strftime(date_end, "%Y-%m-%dT%H:%M:%S")
date_start_str = datetime.strftime(date_start, "%Y-%m-%dT%H:%M:%S")

min_distance_pid = []
ipid_best = [] 
distance_pid = []
distance_lvlh_pid = []
nb_seconds_since_start_pid = []
date_datetime_round_sec_spock_pid_ok = []

nkp = len(kplist); nkd = len(kdlist); nki = len(kilist); 
nk = (int)( 2. / step_drho ) + 1# this is the maximum of iteration to find the optim rho (since rho_control vaies from -1 to 1 with a step step_drho

last_r0_pid = np.zeros([nk])
last_r1_pid = np.zeros([nk])
last_r2_pid = np.zeros([nk])
last_v0_pid = np.zeros([nk])
last_v1_pid = np.zeros([nk])
last_v2_pid = np.zeros([nk])
index_obs_kept = []
date_obs_pid_ok = []
pid_center = 1 # factor to apply to each pid_mod_arr[ipid]
pid_center_list = []
klist = np.zeros([nk,3])
for ikp in range(nkp):
    kp = kplist[ikp]
    for ikd in range(nkd):
        kd = kdlist[ikd]
        for iki in range(nki):
            ki = kilist[iki]
            klist[ikp*nkd*nki + ikd*nki + iki, 0] = kp
            klist[ikp*nkd*nki + ikd*nki + iki, 1] = kd
            klist[ikp*nkd*nki + ikd*nki + iki, 2] = ki


rho_control = np.zeros([nb_interval]) -0.5 # -.5 cause first should starts at -0.5
step_derivative = 3 # how may indices to go bakward to compute slope
curr_input = -1
index_period_spock_all_inter = []
nb_total_index_spock = 0
distance_lvlh_pid_orbit_average_interval_all_inter = []
ecc_orbit_average_interval_all_inter = []
ecc_obs_orbit_average_interval_all_inter = []
ecc_obs_same_spock_all_inter = []
ecc_spock_ok_pid_all_inter = []
#duration_simu = (nb_interval - inter_start_algo) * step_move_save + inter_start_algo * interval # first inter_start_algo are interval hour long, then rest of the intervals are step_move_save hour long
index_period_spock_step_move = np.zeros([nb_interval]) # index in orbit average variables of the end of each interval
index_step_move_save = []
for iinter in range(nb_interval):#!!!!! shoul be nb_interval):
    if iinter < inter_start_algo:
        step_move = interval
    else:
        step_move = step_move_save
    step_move_sec = step_move * 3600.
    variation_slope_previous_rho = 1e30
    last_error_previous_rho = 1e30
    distance_lvlh_pid_orbit_average_interval_sub = []
    ecc_orbit_average_interval_sub = []
    ecc_obs_orbit_average_interval_sub = []
    nb_seconds_since_start_pid_inter = []
    index_obs_kept_inter = []
    distance_pid_interval = []
    distance_lvlh_pid_interval = []
    print ''
    print ''
    print 'NEW INTERVAL', iinter, nb_interval-1,
    if iinter == 0:
        # This calcualted in first aprt of 071318_spock_odtk_ensemble_new_iteration_on_rv
        r0 = '-2.54076587561000e+03'#'-5427.4037595249'#'-2.54076587561000e+03' 
        r1 = '-5.06266991514000e+03'#'4243.4484776488'#'-5.06266991514000e+03' 
        r2 = '-3.95089081204000e+03'#'106.6657082939'#'-3.95089081204000e+03' 
        v0 = '6.76840081300000e+00'#'-3.9026123016'#'6.76840081300000e+00' 
        v1 = '-3.44599707500000e+00'#'-4.8769691960'#'-3.44599707500000e+00'
        v2 = '4.16872280000000e-02'#'-4.3581694849'#'4.16872280000000e-02' 

    else:
        r0 = format(last_r0, '.14e')
        r1 = format(last_r1, '.14e')
        r2 = format(last_r2, '.14e')
        v0 = format(last_v0, '.14e')
        v1 = format(last_v1, '.14e')
        v2 = format(last_v2, '.14e')

    print 'Initial state:', r0, r1, r2, v0, v1, v2
    print ''

    drho = 0 # first iteration eblow: keep the same rho_control as the previous interval
    irho = -1
    error_change_sign = 0
    while ( irho == -1): # just do it once
        irho = irho + 1
        # if iinter == 6: # !!!!!!!!!!! REMOVE
        #     rho_control[iinter-1]  = -0.4

        if iinter >= inter_start_algo:

            rho_control[iinter] = rho_control[iinter-1] + drho
            
        else:# don't apply the PID for the first 3 intervals (36 horus) and take rho_control = 0
            rho_control[iinter] = -0.5 # !!!!!!uncomment
            error_change_sign = 1 #no iteration for this interval
        # !!!!!!! remove block
        print 'irho ' + str(irho) + ' | rho_control ' + str(rho_control[iinter])

        # # rho_control can't be below -1 otherwises the density si negative. 
        # if rho_control[iinter] <= -1:
        #     rho_control[iinter] = -0.99

        main_input_filename = spock_in
        isc = 0
        var_in, var_in_order = read_input_file(main_input_filename)

        output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
        output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
        var_to_read = ["position", "velocity", "eccentricity"]
        var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
        date_spock = np.array(var_out[find_in_read_input_order_variables(var_out_order, 'date')])
        date_datetime_round_sec_spock = np.array(var_out[find_in_read_input_order_variables(var_out_order, 'date_datetime_round_sec')])
        r_spock = var_out[find_in_read_input_order_variables(var_out_order, 'position')]
        v_spock = var_out[find_in_read_input_order_variables(var_out_order, 'velocity')]
        ecc_spock = var_out[find_in_read_input_order_variables(var_out_order, 'eccentricity')]
        n_spock = len(date_spock)

        if irho == 0:
            index_spock_same_date_as_obs_pid = []
            #if iinter == 0: # for the next interval, start date_obs[iobs] at the last observation of the previous interval
            iobs = np.where(np.array(date_obs) == date_start)[0][0]
            #ipdb.set_trace()
            print 'iobs', iobs
            while iobs < nb_obs:
                if date_obs[iobs] > date_datetime_round_sec_spock[-1]:
                    break
                else:
                    if len(index_spock_same_date_as_obs_pid) == 0:
                        first_obs = iobs
                    if len(np.where(date_datetime_round_sec_spock == date_obs[iobs])[0]) != 0:#can be = 0 if an observation is missing at that time
                        index_spock_same_date_as_obs_pid.append(np.where(date_datetime_round_sec_spock == date_obs[iobs])[0][0])
                        nb_seconds_since_start_pid_inter.append( ( date_obs[iobs] - date_obs[0] ).total_seconds() )
                        index_obs_kept_inter.append(iobs)
                        iobs = iobs + 60
                    else: # find next obs
                        while len(np.where(date_datetime_round_sec_spock == date_obs[iobs])[0]) == 0:
                            iobs = iobs + 1
                        index_spock_same_date_as_obs_pid.append(np.where(date_datetime_round_sec_spock == date_obs[iobs])[0][0])
                        nb_seconds_since_start_pid_inter.append( ( date_obs[iobs] - date_obs[0] ).total_seconds() )
                        index_obs_kept_inter.append(iobs)
                        iobs = iobs + 60

                        
            nb_seconds_since_start_pid.append(nb_seconds_since_start_pid_inter)
            index_obs_kept.append(index_obs_kept_inter)
            n = len(index_spock_same_date_as_obs_pid) #!!!!!!!!!! j-index_interval[iinter]
            date_datetime_round_sec_spock_pid_ok.append(date_datetime_round_sec_spock[index_spock_same_date_as_obs_pid])
            date_obs_pid_ok.append(np.array(date_obs)[index_obs_kept[-1]])

        # Compare SpOCK and data
        date_datetime_round_sec_spock_ok_pid = date_datetime_round_sec_spock[index_spock_same_date_as_obs_pid]
        r_spock_ok_pid = np.zeros([n, 3])
        r_spock_ok_pid[:, 0] = r_spock[index_spock_same_date_as_obs_pid, 0]
        r_spock_ok_pid[:, 1] = r_spock[index_spock_same_date_as_obs_pid, 1]
        r_spock_ok_pid[:, 2] = r_spock[index_spock_same_date_as_obs_pid, 2]
        v_spock_ok_pid = np.zeros([n, 3])
        v_spock_ok_pid[:, 0] = v_spock[index_spock_same_date_as_obs_pid, 0]
        v_spock_ok_pid[:, 1] = v_spock[index_spock_same_date_as_obs_pid, 1]
        v_spock_ok_pid[:, 2] = v_spock[index_spock_same_date_as_obs_pid, 2]
        ecc_spock_ok_pid = np.zeros([n])
        ecc_spock_ok_pid = ecc_spock[index_spock_same_date_as_obs_pid]

        #if pid_mod_arr[irho] == 1:
        if irho == 0:
            #ipdb.set_trace()
            index_step_move = -2#np.where(date_datetime_round_sec_spock_ok_pid == (date_start + timedelta(seconds = step_move_sec)))[0][0]
            index_step_move_save.append(index_step_move)
        last_r0_pid[irho] = r_spock_ok_pid[index_step_move, 0]
        last_r1_pid[irho] = r_spock_ok_pid[index_step_move, 1]
        last_r2_pid[irho] = r_spock_ok_pid[index_step_move, 2]
        last_v0_pid[irho] = v_spock_ok_pid[index_step_move, 0]
        last_v1_pid[irho] = v_spock_ok_pid[index_step_move, 1]
        last_v2_pid[irho] = v_spock_ok_pid[index_step_move, 2]
        print 'interval', iinter
        #print 'Final state', format(last_r0_pid[irho],".14e"),format(last_r1_pid[irho],".14e"),format(last_r2_pid[irho],".14e"), format(last_v0_pid[irho],".14e"),format(last_v1_pid[irho],".14e"),format(last_v2_pid[irho],".14e")
            


        distance_lvlh_pid_sub = []
        # distance represents the along-track distance from observation to SpOCK. > 0 means SpOCK is trailing, < 0 means SpOCK is leading
        nb_seconds_since_start_interval = 0
        iperiod_find = 0
        if irho == 0:
            index_period_spock_temp = []
            ecc_obs_same_spock_all_inter.append( ecc_obs[index_obs_kept[-1]] )
        for i in range(n):
            distance_here = r_obs[index_obs_kept[-1]][i, :] - r_spock_ok_pid[i, :]
            distance_lvlh_pid_sub.append( eci_to_lvlh(r_spock_ok_pid[i, :], v_spock_ok_pid[i, :], distance_here)[0] ) #[0]: along-track
            # "if" below used to calculate orbit average error
            if irho == 0:
                nb_seconds_since_start_interval = ( date_datetime_round_sec_spock_pid_ok[-1][i] - date_datetime_round_sec_spock_pid_ok[-1][0] ).total_seconds()
                if nb_seconds_since_start_interval >= orbit_period * iperiod_find:
                    iperiod_find = iperiod_find + 1
                    index_period_spock_temp.append(i)
                    if irho == 0:
                        if ((nb_seconds_since_start_interval >= step_move_sec) & (index_period_spock_step_move[iinter] == 0 )): # don't go twice in here so second condition
                            index_period_spock_step_move[iinter] =  len(index_period_spock_temp) - 1
        
        distance_lvlh_pid_interval.append( distance_lvlh_pid_sub )
        # calculate orbit average error
        if irho == 0:
            nb_period_this_inter = iperiod_find
            index_period_spock = np.zeros([2*(nb_period_this_inter-1)])
            nb_total_index_spock = nb_total_index_spock + len(r_spock_ok_pid)
        distance_lvlh_pid_orbit_average_interval = np.zeros([2*(nb_period_this_inter-1)])
        ecc_orbit_average_interval = np.zeros([2*(nb_period_this_inter-1)])
        ecc_obs_orbit_average_interval = np.zeros([2*(nb_period_this_inter-1)])
        for iperiod in range(nb_period_this_inter-1):
            iorbit_start = index_period_spock_temp[iperiod]
            iorbit_stop = index_period_spock_temp[iperiod+1]
            distance_lvlh_pid_orbit_average_interval[2*iperiod] = np.mean( distance_lvlh_pid_sub[iorbit_start:iorbit_stop] ) 
            distance_lvlh_pid_orbit_average_interval[2*iperiod+1] = np.mean( distance_lvlh_pid_sub[iorbit_start:iorbit_stop] ) 
            ecc_orbit_average_interval[2*iperiod] = np.mean( ecc_spock_ok_pid[iorbit_start:iorbit_stop] ) 
            ecc_orbit_average_interval[2*iperiod+1] = np.mean( ecc_spock_ok_pid[iorbit_start:iorbit_stop] ) 

            ecc_obs_orbit_average_interval[2*iperiod] = np.mean( ecc_obs_same_spock_all_inter[-1][iorbit_start:iorbit_stop] ) 
            ecc_obs_orbit_average_interval[2*iperiod+1] = np.mean( ecc_obs_same_spock_all_inter[-1][iorbit_start:iorbit_stop] ) 

            index_period_spock[2*iperiod] = index_period_spock_temp[iperiod]
            index_period_spock[2*iperiod+1] = index_period_spock_temp[iperiod+1]
        #ipdb.set_trace()    
            #print 'ZZZZZZZZZZZZZZZZZZZZZZZZ', distance_lvlh_pid_orbit_average_interval_sub
        distance_lvlh_pid_orbit_average_interval_sub.append(distance_lvlh_pid_orbit_average_interval)
        ecc_orbit_average_interval_sub.append(ecc_orbit_average_interval)
        ecc_obs_orbit_average_interval_sub.append(ecc_obs_orbit_average_interval)

        #print 'YYYYYYYYYYYYYYYYYYYYYYY', distance_lvlh_pid_orbit_average_interval
        #print 'XXXXXXXXXXXXXXXXXXXXX', distance_lvlh_pid_orbit_average_interval_sub
        # if irho == 1:
        #     raise Exception
        slope_orbit_average_error_inter_start_to_end = distance_lvlh_pid_orbit_average_interval[-1] - distance_lvlh_pid_orbit_average_interval[0]
        last_error = distance_lvlh_pid_orbit_average_interval[-1]
        variation_slope = np.abs(slope_orbit_average_error_inter_start_to_end)

        sign_slope = np.sign(slope_orbit_average_error_inter_start_to_end)
        if irho == 0: # the first iteration tells which direction to go, ie if we need to add or remove density
            sign_step_drho = sign_slope
            index_period_spock_all_inter.append(list(index_period_spock))#+nb_total_index_spock))
        #if ( (sign_slope != sign_step_drho) | (iinter <= 1)): # in this case (first condition), it means that we've added (or removed) too much density and inverted the trend so we found the optimum rho. So we'll move to the next inter. iinter<=2 is here becasue for the first 3 intervals ndistance_lvlh_pid_oro iteratio
        sign_last_error_previous_rho = np.sign(last_error_previous_rho)
        if (irho == 0):
            sign_last_error_previous_rho = np.sign(last_error)

        error_change_sign = 1

        # start and end dates of next interval
        date_start = date_start + timedelta(seconds = step_move_sec)
        date_start_str = datetime.strftime(date_start, "%Y-%m-%dT%H:%M:%S")
        date_end_str = datetime.strftime(date_start + timedelta(seconds = interval_sec), "%Y-%m-%dT%H:%M:%S")



        distance_lvlh_pid_orbit_average_interval_all_inter.append(distance_lvlh_pid_orbit_average_interval_sub)
        ecc_orbit_average_interval_all_inter.append(ecc_orbit_average_interval_sub)
        ecc_obs_orbit_average_interval_all_inter.append(ecc_obs_orbit_average_interval_sub)
        ecc_spock_ok_pid_all_inter.append(ecc_spock_ok_pid)


    distance_lvlh_pid.append( distance_lvlh_pid_interval )
    # if iinter == 3:
    #     raise Exception
    print '-->\n-->'

    #if iinter == 0:
    ################### FIGURES ###################
    height_fig = 11
    ratio_fig_size = 4./3
    fontsize_plot = 20

    ######
    fig_title = ''#'Distance between SpOCK and data for different density coefficient' #'Distance with respect to pid = 0.7'#'Distance between SpOCK and data for different density coefficient'
    y_label = 'Distance (m)'
    x_label = 'Time (hours)' 

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
        


    for iinter_loop in range(iinter+1):

        iii = -1
        #index_period_spock_step_move[iinter_loop] = -0.5 # *2 gives the last index, which is what we want for the first inter_start_algo because they last interval hours
        # DISTANCE
        ax.plot(np.array(nb_seconds_since_start_pid[iinter_loop])[:index_step_move_save[iinter_loop]+1]/3600., np.array(distance_lvlh_pid[iinter_loop][iii])[:index_step_move_save[iinter_loop]+1]*1000, linewidth = 2, color = 'b', alpha = 0.3)
        ax.plot(np.array(nb_seconds_since_start_pid[iinter_loop])[np.array(index_period_spock_all_inter[iinter_loop]).astype(int)]/3600., np.array(distance_lvlh_pid_orbit_average_interval_all_inter[iinter_loop][iii])*1000, linewidth = 4, color = 'magenta') #-2 becaue the last one is to throw away

#         # ECCENTRICITY
#         ax.plot(np.array(nb_seconds_since_start_pid[iinter_loop])[:index_step_move_save[iinter_loop]+1]/3600., ecc_obs_same_spock_all_inter[iinter_loop][:index_step_move_save[iinter_loop]+1], linewidth = 2, color = 'r', alpha = 0.3)
#         ax.plot(np.array(nb_seconds_since_start_pid[iinter_loop])[np.array(index_period_spock_all_inter[iinter_loop]).astype(int)]/3600., np.array(ecc_obs_orbit_average_interval_all_inter[iinter_loop][iii]), linewidth = 4, color = 'r') #-2 becaue the last one is to throw away


#         ax.plot(np.array(nb_seconds_since_start_pid[iinter_loop])[:index_step_move_save[iinter_loop]+1]/3600., ecc_spock_ok_pid_all_inter[iinter_loop][:index_step_move_save[iinter_loop]+1], linewidth = 2, color = 'b', alpha = 0.3)
#         ax.plot(np.array(nb_seconds_since_start_pid[iinter_loop])[np.array(index_period_spock_all_inter[iinter_loop]).astype(int)]/3600., np.array(ecc_orbit_average_interval_all_inter[iinter_loop][iii]), linewidth = 4, color = 'b') #-2 becaue the last one is to throw away




#ecc_spock_ok_pid_all_inter

#         if iinter == 6:
#             ipdb.set_trace()
#             plt.show(); plt.show()
        # ax.plot(np.array(nb_seconds_since_start_pid[-1])[index_period_spock_all_inter[-1]][:-1]/3600, distance_lvlh_pid_orbit_average_interval)
        # ax.plot(np.array(nb_seconds_since_start_pid[1]) / 3600, distance_lvlh_pid_sub)

        #ax.plot([np.array(nb_seconds_since_start_pid[iinter_loop][0])/3600., np.array(nb_seconds_since_start_pid[iinter_loop][0])/3600.], [-200, np.array(distance_lvlh_pid[iinter_loop][iii][0])*1000], linewidth = 2, color = 'r', linestyle = 'dashed')
        #ax.text(np.array(nb_seconds_since_start_pid[iinter_loop][0])/3600., np.array(distance_lvlh_pid[iinter_loop][iii][0])*1000, format(rho_control[iinter_loop], ".2f"), horizontalalignment = 'left', fontsize = fontsize_plot, weight = 'bold', color = 'r', verticalalignment = 'center', label = 'rho_control')
    #ax.plot([0, duration_simu], [0,0], linestyle = 'dashed', linewidth = 1, color = 'black')

    #ax.set_xlim([0, duration_simu]); 
    ax.set_ylim([-600, 200]) 
    ax.margins(0,0)
    #ax.text(duration_simu/2., -200, 'SpOCK in front -> need rho_control < 0', horizontalalignment = 'center', verticalalignment = 'bottom', fontsize = fontsize_plot, weight = 'bold')
    #ax.text(duration_simu/2., 1200, 'SpOCK behind -> need rho_control > 0', horizontalalignment = 'center', verticalalignment = 'top', fontsize = fontsize_plot, weight = 'bold')



    fig_save_name = 'fig/' + spock_in.replace(".txt","_") + str(nb_interval) + ".pdf"
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  



raise Exception

# Optimum rho vs time
height_fig = 11
ratio_fig_size = 4./3
fontsize_plot = 20

######
fig_title = ''
y_label = 'Rho control'
x_label = 'Time (hours)' 

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



time_array = np.zeros([nb_interval])
for iinter in range(nb_interval):
    if iinter <= inter_start_algo:
        time_array[iinter] = iinter*interval
    else:
        time_array[iinter] = time_array[iinter-1] + step_move_save
ax.plot(time_array, rho_control, linewidth = 2, color = 'k')
ax.scatter(time_array, rho_control, linewidth = 2, color = 'k', marker = '.')
ax.set_yticks(np.arange(-1,1+step_drho,step_drho))
ax.set_xlim([0, duration_simu]); ax.set_ylim([-1.1, 1.1])
ax.plot([0, duration_simu], [0, 0], linewidth = 2, color = 'k', linestyle = 'dashed')
ax.text(0,0.07,'MSIS', horizontalalignment = 'left', verticalalignment = 'top', fontsize = fontsize_plot, weight = 'bold' )
ax.margins(0,0)


fig_save_name = 'fig/optim_nbinter' + str(nb_interval) + "_rho_control.pdf"
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

    
raise Exception

distance_lvlh_pid_concantenate = []
nb_seconds_since_start_pid_concatenate = []
for irho in range(nk):
    distance_lvlh_pid_concantenate_irho = []
    for iinter in range(nb_interval):
        distance_lvlh_pid_concantenate_irho = distance_lvlh_pid_concantenate_irho + distance_lvlh_pid[iinter][irho]
        if irho == 0:
            nb_seconds_since_start_pid_concatenate = nb_seconds_since_start_pid_concatenate + nb_seconds_since_start_pid[iinter]
    distance_lvlh_pid_concantenate.append(distance_lvlh_pid_concantenate_irho )



distance_min_concatenate = []
for iinter in range(nb_interval):
    distance_min_concatenate = distance_min_concatenate + distance_lvlh_pid[iinter][ipid_best[iinter]] 
distance_min_concatenate = np.array(distance_min_concatenate)



main_input_filename_root =  date_obs_start_str.replace(":","_") + '_' + date_obs_end_str.replace(":","_") + '_interval' + format(interval, ".1f").replace(".","_")+ '.txt'
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx
# raise Exception
# pickle.dump(nb_seconds_since_start_pid, open('pickle/nb_seconds_since_start_pid'+ '.pickle', 'w'))
# pickle.dump(nk, open('pickle/nk'+ '.pickle', 'w'))
# pickle.dump(distance_lvlh_pid, open('pickle/distance_lvlh_pid'+ '.pickle', 'w'))
# pickle.dump(pid_mod_arr, open('pickle/pid_mod_arr'+ '.pickle', 'w'))
# pickle.dump(date_datetime_round_sec_spock_pid_ok, open('pickle/date_datetime_round_sec_spock_pid_ok'+ '.pickle', 'w'))
# pickle.dump(main_input_filename_root, open('pickle/main_input_filename_root'+ '.pickle', 'w'))
# pickle.dump(distance_lvlh_pid_concantenate, open('pickle/distance_lvlh_pid_concantenate'+ '.pickle', 'w'))
# pickle.dump(nb_seconds_since_start_pid_concatenate, open('pickle/nb_seconds_since_start_pid_concatenate'+ '.pickle', 'w'))

# nb_seconds_since_start_pid = pickle.load(open('pickle/nb_seconds_since_start_pid' + '.pickle'))
# nk = pickle.load(open('pickle/nk' + '.pickle'))
# distance_lvlh_pid = pickle.load(open('pickle/distance_lvlh_pid' + '.pickle'))
# pid_mod_arr = pickle.load(open('pickle/pid_mod_arr' + '.pickle'))
# date_datetime_round_sec_spock_pid_ok = pickle.load(open('pickle/date_datetime_round_sec_spock_pid_ok' + '.pickle'))
# main_input_filename_root = pickle.load(open('pickle/main_input_filename_root' + '.pickle'))
# distance_lvlh_pid_concantenate = pickle.load(open('pickle/distance_lvlh_pid_concantenate' + '.pickle'))
# nb_seconds_since_start_pid_concatenate = pickle.load(open('pickle/nb_seconds_since_start_pid_concatenate' + '.pickle'))

# raise Exception






################### FIGURES ###################
height_fig = 11
ratio_fig_size = 4./3
fontsize_plot = 20

######
fig_title = ''#'Distance between SpOCK and data for different density coefficient' #'Distance with respect to pid = 0.7'#'Distance between SpOCK and data for different density coefficient'
y_label = 'Distance (m)'
x_label = 'Time (hours)' 

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

x_axis = np.array(nb_seconds_since_start_pid_concatenate)
if nk > 1:
    alpha_arr = np.arange(0.2,1+0.2/nk,(1-0.2)/(nk-1))
else: 
    alpha_arr = [1]
min_y = 1.e30
max_y = -1.e30

for iinter in range(nb_interval):
    for ik in range(nk):
        if alpha_arr[ik] >1:
            alpha_arr[ik] = 1
        ax.plot(np.array(nb_seconds_since_start_pid[iinter])/3600., np.array(distance_lvlh_pid[iinter][ik])*1000, linewidth = 2, color = 'b', alpha = alpha_arr[ik])
        ax.plot([np.array(nb_seconds_since_start_pid[iinter][0])/3600., np.array(nb_seconds_since_start_pid[iinter][0])/3600.],[-1e30, np.array(distance_lvlh_pid[iinter][ik][0])*1000], linewidth = 2, linestyle = 'dashed', color = 'red')

        if np.min(np.array(distance_lvlh_pid[iinter][ik])*1000) < min_y:
            min_y = np.min(np.array(distance_lvlh_pid[iinter][ik])*1000)
        if np.max(np.array(distance_lvlh_pid[iinter][ik])*1000) > max_y:
            max_y = np.max(np.array(distance_lvlh_pid[iinter][ik])*1000)        
        print ik

delta_y_text = (max_y - min_y) / 20.
for iinter in range(nb_interval): # can't couple that with previous loop because need to figure out min_y and max_y firs
    for ik in range(nk):
        # rho_control
        ax.text(np.array(nb_seconds_since_start_pid[iinter][0])/3600., np.array(distance_lvlh_pid[iinter][ik][0])*1000, format(rho_control[iinter], ".2f"), horizontalalignment = 'left', fontsize = fontsize_plot, weight = 'bold', color = 'r', alpha = alpha_arr[ik], verticalalignment = 'center', label = 'rho_control')

        # err
        ax.text(np.array(nb_seconds_since_start_pid[iinter][-1])/3600., np.array(distance_lvlh_pid[iinter][ik][-1])*1000 - delta_y_text, format(err[iinter]*1000, ".2f"), horizontalalignment = 'right', fontsize = fontsize_plot, weight = 'bold', color = 'b', alpha = alpha_arr[ik], verticalalignment = 'center', label = 'error')
        # derr/dr
        ax.text(np.array(nb_seconds_since_start_pid[iinter][-1])/3600., np.array(distance_lvlh_pid[iinter][ik][-1])*1000 - 2*delta_y_text/1.2, format(derrdt[iinter]*1000, ".2f"), horizontalalignment = 'right', fontsize = fontsize_plot, weight = 'bold', color = 'magenta', alpha = alpha_arr[ik], verticalalignment = 'center', label = 'derror/dt')
        # interr
        ax.text(np.array(nb_seconds_since_start_pid[iinter][-1])/3600., np.array(distance_lvlh_pid[iinter][ik][-1])*1000 - 3*delta_y_text/1.2, format(interr[iinter]*1000, ".2f"), horizontalalignment = 'right', fontsize = fontsize_plot, weight = 'bold', color = 'grey', alpha = alpha_arr[ik], verticalalignment = 'center', label = 'interror')


ax.plot([np.array(nb_seconds_since_start_pid[0][0])/3600., np.array(nb_seconds_since_start_pid[nb_interval-1][-1])/3600.],[0, 0], linewidth = 2, linestyle = 'dashed', color = 'blue')
ax.plot(x_axis/3600., distance_min_concatenate*1000, linewidth = 2, color = 'r')
ax.text(x_axis[-1]/3600., max_y/2., 'SpOCK\nbehind', horizontalalignment = 'left', fontsize = fontsize_plot, weight = 'bold',verticalalignment = 'center')
ax.text(x_axis[-1]/3600., min_y/2., 'SpOCK\nin front', horizontalalignment = 'left', fontsize = fontsize_plot, weight = 'bold',verticalalignment = 'center')
# # x axis label is in real time
# nb_seconds_in_simu = nb_seconds_since_start_pid_concatenate[-1] - nb_seconds_since_start_pid_concatenate[0]
# start_xaxis_label = nb_seconds_since_start_pid_concatenate[0]
# date_ref = date_datetime_round_sec_spock_pid_ok[0][0]
# nb_ticks_xlabel = 10
# dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
# xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
# date_list_str = []
# date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
# for i in range(len(xticks)):
#     if dt_xlabel > nb_ticks_xlabel*24*3600:
#         date_list_str.append( str(date_list[i])[5:10] )
#     else:
#         date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
#         ax.xaxis.set_ticks(xticks)
#         ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
ax.margins(0,0); ax.set_ylim([min_y, max_y])
#        ax.set_xlim([ax.get_xlim()[0], most_recent_tle_among_all_sc])

legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="", fontsize = fontsize_plot)
fig_save_name = 'fig/distance_optimum_pid_to_obs_' + main_input_filename_root.replace(".txt", ".pdf")
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

os.system("rsync -av " + fig_save_name + " cbv@srbwks2014-0008.engin.umich.edu:")



raise Exception




################### FIGURES ###################
height_fig = 11
ratio_fig_size = 4./3
fontsize_plot = 20

######
fig_title = ''#'Distance_Lvlh between SpOCK and data for different density coefficient' #'Distance_Lvlh with respect to pid = 0.7'#'Distance_Lvlh between SpOCK and data for different density coefficient'
y_label = 'Distance_Lvlh (m)'
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

x_axis = nb_seconds_since_start_pid_concatenate
if nk > 1:
    alpha_arr = np.arange(0.2,1+0.2/nk,(1-0.2)/(nk-1))
else: 
    alpha_arr = [1]

for iinter in range(nb_interval):
    for ik in range(nk):
        if alpha_arr[ik] >1:
            alpha_arr[ik] = 1
        ax.plot(nb_seconds_since_start_pid[iinter], np.array(distance_lvlh_pid[iinter][ik])*1000, linewidth = 2, color = 'b', alpha = alpha_arr[ik])
        #ax.plot(x_axis, (np.array(distance_lvlh_pid[iinter][ik]) - np.array(distance_lvlh_pid[iinter][0]))*1000, linewidth = 2, color = 'b', alpha = alpha_arr[ik])
        #ax.text(x_axis[-1], (np.array(distance_lvlh_pid[iinter][ik]) - np.array(distance_lvlh_pid[iinter][0]))[-1]*1000, format(pid_mod_arr[ik], ".1f"), horizontalalignment = 'left', fontsize = fontsize_plot, weight = 'bold', color = 'b', alpha = alpha_arr[ik], verticalalignment = 'center')
        #ax.text(x_axis[-1], distance_lvlh_min_concatenate[ik][-1]*1000, format(pid_mod_arr[ik], ".1f"), horizontalalignment = 'left', fontsize = fontsize_plot, weight = 'bold', color = 'b', alpha = alpha_arr[ik], verticalalignment = 'center')
        #ax.text(x_axis[-1], distance_lvlh_pid[iinter][ik][-1], str(pid_mod_arr[ik]), horizontalalignment = 'left', fontsize = fontsize_plot, weight = 'bold', color = 'b', alpha = alpha_arr[ik], verticalalignment = 'center')
        print ik
ax.plot(x_axis, distance_lvlh_min_concatenate*1000, linewidth = 2, color = 'r')
# x axis label is in real time
nb_seconds_in_simu = nb_seconds_since_start_pid_concatenate[-1] - nb_seconds_since_start_pid_concatenate[0]
start_xaxis_label = nb_seconds_since_start_pid_concatenate[0]
date_ref = date_datetime_round_sec_spock_pid_ok[0][0]
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


fig_save_name = 'distance_lvlh_optimum_pid_to_obs_' + main_input_filename_root.replace("txt", "pdf")
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  




#######
fig_title = ''#'Distance between SpOCK and data for different density coefficient' #'Distance with respect to pid = 0.7'#'Distance between SpOCK and data for different density coefficient'
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

x_axis = nb_seconds_since_start_pid_concatenate
if nk > 1:
    alpha_arr = np.arange(0.2,1+0.2/nk,(1-0.2)/(nk-1))
else: 
    alpha_arr = [1]
for ik in range(nk):
    if alpha_arr[ik] >1:
        alpha_arr[ik] = 1

    ax.plot(x_axis, np.array(distance_lvlh_pid_concantenate[ik])*1000, linewidth = 2, color = 'b', alpha = alpha_arr[ik])
    #ax.plot(x_axis, (np.array(distance_lvlh_pid[iinter][ik]) - np.array(distance_lvlh_pid[iinter][0]))*1000, linewidth = 2, color = 'b', alpha = alpha_arr[ik])
    #ax.text(x_axis[-1], (np.array(distance_lvlh_pid[iinter][ik]) - np.array(distance_lvlh_pid[iinter][0]))[-1]*1000, format(pid_mod_arr[ik], ".1f"), horizontalalignment = 'left', fontsize = fontsize_plot, weight = 'bold', color = 'b', alpha = alpha_arr[ik], verticalalignment = 'center')
    ax.text(x_axis[-1], distance_lvlh_pid_concantenate[ik][-1]*1000, format(pid_mod_arr[ik], ".1f"), horizontalalignment = 'left', fontsize = fontsize_plot, weight = 'bold', color = 'b', alpha = alpha_arr[ik], verticalalignment = 'center')
    #ax.text(x_axis[-1], distance_lvlh_pid  [iinter][ik][-1], str(pid_mod_arr[ik]), horizontalalignment = 'left', fontsize = fontsize_plot, weight = 'bold', color = 'b', alpha = alpha_arr[ik], verticalalignment = 'center')
    print ik                                           
# x axis label is in real time
nb_seconds_in_simu = nb_seconds_since_start_pid_concatenate[-1] - nb_seconds_since_start_pid_concatenate[0]
start_xaxis_label = nb_seconds_since_start_pid_concatenate[0]
date_ref = date_datetime_round_sec_spock_pid_ok[0][0]
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


fig_save_name = 'pid_distance_ens_to_observations_' + main_input_filename_root.replace("txt", "pdf")
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  





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

fig_save_name = 'x_eci_' + main_input_filename.replace(".txt", "_test.pdf")
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

fig_save_name = 'y_eci_' + main_input_filename.replace(".txt", "_test.pdf")
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

fig_save_name = 'z_eci_' + main_input_filename.replace(".txt", "_test.pdf")
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

fig_vsave_name = 'vx_eci_' + main_input_filename.replace(".txt", "_test.pdf")
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

fig_vsave_name = 'vy_eci_' + main_input_filename.replace(".txt", "_test.pdf")
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

fig_vsave_name = 'vz_eci_' + main_input_filename.replace(".txt", "_test.pdf")
fig_vz.savefig(fig_vsave_name, facecolor=fig_vz.get_facecolor(), edgecolor='none', bbox_inches='tight')  




