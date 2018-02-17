import matplotlib.colors as colors
import matplotlib.cm as cmx
from cadre_read_last_tle import *
from get_prop_dir import *
import matplotlib.gridspec as gridspec
from read_input_file import *
from matplotlib.colors import LogNorm
import pickle
from eci_to_lvlh import *
import sys
import fileinput
import time
import datetime
import numpy as np
from matplotlib import pyplot as plt
import os
import subprocess

show_plots = 1
save_results = 1
root_save_fig_name = 'collision/more_output_but_no_coll_1126_ok_'

# # runs have to be the same epoch end and start and dt
# run_list = ['just_unperturbed_with_tca_storm_dynamic_32600ens_min_dist_coll_1_2m.txt',
#             'just_unperturbed_with_tca_storm_static_32600ens_min_dist_coll_1_2m.txt']

run_list = ['more_output_but_no_coll_1126_ok_quartile_f107_1_quartile_ap_1.txt',
            'more_output_but_no_coll_1126_ok_quartile_f107_2_quartile_ap_2.txt',
            'more_output_but_no_coll_1126_ok_quartile_f107_3_quartile_ap_3.txt',
            'more_output_but_no_coll_1126_ok_quartile_f107_4_quartile_ap_4.txt',
            'more_output_but_no_coll_1126_ok_quartile_f107_5_quartile_ap_5.txt',
            'more_output_but_no_coll_1126_ok_quartile_f107_6_quartile_ap_6.txt',
            'more_output_but_no_coll_1126_ok_quartile_f107_7_quartile_ap_7.txt',
            'more_output_but_no_coll_1126_ok_quartile_f107_8_quartile_ap_8.txt',
            'more_output_but_no_coll_1126_ok_quartile_f107_9_quartile_ap_9.txt']

# run_list = ['devel_1126_ok_quartile_f107_1_quartile_ap_1.txt',
#              '1126_ok_quartile_f107_2_quartile_ap_2.txt',
#              '1126_ok_quartile_f107_3_quartile_ap_3.txt',
#              '1126_ok_quartile_f107_4_quartile_ap_4.txt',
#              '1126_ok_quartile_f107_5_quartile_ap_5.txt'] 
             # '1126_ok_quartile_f107_6_quartile_ap_6.txt',
             # '1126_ok_quartile_f107_7_quartile_ap_7.txt',
             # '1126_ok_quartile_f107_8_quartile_ap_8.txt',
             # '1126_ok_quartile_f107_9_quartile_ap_9.txt']


median_run = 4
nrun = len(run_list)

# GENERATE DISCTINCT COLORS
NCURVES = nrun
np.random.seed(101)
curves = [np.random.random(20) for i in range(NCURVES)]
values = range(NCURVES)
jet = cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

eighty_quartile_along = []
fifty_quartile_along = []
std_quartile_along = []
eighty_quartile_angle = []
fifty_quartile_angle = []
std_quartile_angle = []

rho_ensemble_orbit_average_all_run = []
f107_ensemble_all_run = []
ap_ensemble_all_run = []
path_folder_results = '/raid3/Armada/Charles/python/' 
# set up interactive figures
if show_plots == 1:
    plt.ion()

for irun in range(nrun):
    main_input_file_name = get_prop_dir(1) + 'run_collision/input/main_input/' + run_list[irun]

# read input file
    input_variables, order_input_variables = read_input_file(main_input_file_name)
    date_start = input_variables[0]; date_stop = input_variables[1]; dt = input_variables[2]; nb_steps = input_variables[3]; satellite_to_plot_path = input_variables[6]; satellite_to_plot = input_variables[7]; nb_ensembles_coe = input_variables[12]; nb_ensembles_attitude = input_variables[13]; nb_ensembles_cd = input_variables[11]; ensemble_to_plot_temp = input_variables[14]; 
    nb_ensembles_density = input_variables[17]
    n = nb_steps
    nb_spacecraft = input_variables[4]
    compute_drag = input_variables[19]
    isat = 0


    # Name of the mission of the simulation (choices: 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', 'other')
    isat_temp = 0
    if ( 'cygnss' in satellite_to_plot[isat_temp].lower() ):
        name_mission = 'CYGNSS' 
    elif ( 'cadre' in satellite_to_plot[isat_temp].lower() ):
        name_mission = 'CADRE' 
    elif ( 'aerie' in satellite_to_plot[isat_temp].lower() ):
        name_mission = 'AERIE' 
    elif ( 'scion' in satellite_to_plot[isat_temp].lower() ):
        name_mission = 'SCION' 
    elif ( 'qb50' in satellite_to_plot[isat_temp].lower() ):
        name_mission = 'QB50' 
    else:
        name_mission = 'other' 



    # Ensembles created by the propagator
    ensemble_to_plot = []
    for i in range(len(ensemble_to_plot_temp)):
        if (ensemble_to_plot_temp[i] == 'eci_r'):
            ensemble_to_plot.append('x_eci'); ensemble_to_plot.append('y_eci'); ensemble_to_plot.append('z_eci')
        if (ensemble_to_plot_temp[i] == 'eci_v'):
            ensemble_to_plot.append('vx_eci'); ensemble_to_plot.append('vy_eci'); ensemble_to_plot.append('vz_eci')
        if (ensemble_to_plot_temp[i] == 'geodetic'):
            ensemble_to_plot.append('longitude'); ensemble_to_plot.append('latitude'); ensemble_to_plot.append('altitude')
        if (ensemble_to_plot_temp[i] == 'power'):
            ensemble_to_plot.append('power')
        if (ensemble_to_plot_temp[i] == 'attitude'):
            ensemble_to_plot.append('pitch'); ensemble_to_plot.append('roll'); ensemble_to_plot.append('yaw')
        if (ensemble_to_plot_temp[i] == 'oe'):
            ensemble_to_plot.append('sma'); ensemble_to_plot.append('inclination'); ensemble_to_plot.append('eccentricity'); ensemble_to_plot.append('true_anomaly'); ensemble_to_plot.append('RAAN'); ensemble_to_plot.append('argument_perigee');ensemble_to_plot.append('phase_angle');ensemble_to_plot.append('sma_difference');
        if ((ensemble_to_plot_temp[i] == 'density') & (compute_drag == 1)):
            ensemble_to_plot.append('rho'); ensemble_to_plot.append('f107'); ensemble_to_plot.append('f107a'); ensemble_to_plot.append('ap'); 
        if (ensemble_to_plot_temp[i] == 'collision'):
            ensemble_to_plot.append('tca'); ensemble_to_plot.append('dca');

    pickle_list = []
    pickle_list_name = []
    pickle_list_name_temp = []
    ## Nb of ensembles
    nb_ensembles_array = [nb_ensembles_coe, nb_ensembles_attitude, nb_ensembles_cd, nb_ensembles_density]
    nb_ensembles = np.max(nb_ensembles_array)
    for i in range(len(nb_ensembles_array)):
        if ( (nb_ensembles_array[i] > 0 ) &  (nb_ensembles_array[i] < nb_ensembles) ) :
            nb_ensembles = nb_ensembles_array[i]

    name_subfolder_save = satellite_to_plot[isat][:-5] + "/"
    dir_final_output_ensemble = satellite_to_plot_path[isat] + 'ensemble/'
    ## RESTORE THE RESULTS FROM THE PICKLE
    list_save_pickle = [f for f in os.listdir(path_folder_results + name_mission + '/pickle/' + name_subfolder_save) if ( satellite_to_plot[isat].replace(".txt", "") in f ) | ( '_tca.pickle' in f) | ( '_dca.pickle' in f)]

    nb_pickle = len(list_save_pickle)
    for ipickle in range(nb_pickle):
        save_pickle_name = path_folder_results + name_mission + '/pickle/' + name_subfolder_save + list_save_pickle[ipickle]
        with open(save_pickle_name) as f:
            print "Loading pickle " + save_pickle_name + "..."
            pickle_list, pickle_list_name_temp = pickle.load(f)
            for ielt in range(len(pickle_list_name_temp)):
                if ((pickle_list_name_temp[ielt] in pickle_list_name) == False):
                    pickle_list_name.append(pickle_list_name_temp[ielt])
        if ipickle == 0:
            nProcs = pickle_list[ pickle_list_name_temp.index('nProcs') ]
            nb_ensembles= (int)(nb_ensembles / nProcs) * nProcs
            nb_ensembles_per_proc = (int)(nb_ensembles /  nProcs)
        # if nProcs > np.sqrt( nb_ensembles ):
        #     print "Note: to plot the density, F10.7, F10.7A, and Ap of ensembles as a function of time, it is recommended to run SpOCK with nProcs <= sqrt( nb_ensembles ) (it will still run fine and plot fine)."
        if ( 'x_eci_ensemble' in pickle_list_name_temp ): 
            x_eci_ensemble = pickle_list[ pickle_list_name_temp.index('x_eci_ensemble') ]
        if ( 'y_eci_ensemble' in pickle_list_name_temp ): 
            y_eci_ensemble = pickle_list[ pickle_list_name_temp.index('y_eci_ensemble') ]
        if ( 'z_eci_ensemble' in pickle_list_name_temp ): 
            z_eci_ensemble = pickle_list[ pickle_list_name_temp.index('z_eci_ensemble') ]
        if ( 'pitch_ensemble' in pickle_list_name_temp ): 
            pitch_ensemble = pickle_list[ pickle_list_name_temp.index('pitch_ensemble') ]
        if ( 'roll_ensemble' in pickle_list_name_temp ): 
            roll_ensemble = pickle_list[ pickle_list_name_temp.index('roll_ensemble') ]
        if ( 'yaw_ensemble' in pickle_list_name_temp ): 
            yaw_ensemble = pickle_list[ pickle_list_name_temp.index('yaw_ensemble') ]
        if ( 'x_eci_main_sc' in pickle_list_name_temp ):
            x_eci_main_sc = pickle_list[ pickle_list_name_temp.index('x_eci_main_sc') ]
        if ( 'y_eci_main_sc' in pickle_list_name_temp ):
            y_eci_main_sc = pickle_list[ pickle_list_name_temp.index('y_eci_main_sc') ]
        if ( 'z_eci_main_sc' in pickle_list_name_temp ):
            z_eci_main_sc = pickle_list[ pickle_list_name_temp.index('z_eci_main_sc') ]
        if ( 'vx_eci_main_sc' in pickle_list_name_temp ):
            vx_eci_main_sc = pickle_list[ pickle_list_name_temp.index('vx_eci_main_sc') ]
        if ( 'vy_eci_main_sc' in pickle_list_name_temp ):
            vy_eci_main_sc = pickle_list[ pickle_list_name_temp.index('vy_eci_main_sc') ]
        if ( 'vz_eci_main_sc' in pickle_list_name_temp ):
            vz_eci_main_sc = pickle_list[ pickle_list_name_temp.index('vz_eci_main_sc') ]
        if ( 'algebric_distance_ensemble_main_sc_lvlh_x' in pickle_list_name_temp ):
            algebric_distance_ensemble_main_sc_lvlh_x = pickle_list[ pickle_list_name_temp.index('algebric_distance_ensemble_main_sc_lvlh_x') ]
        if ( 'algebric_distance_ensemble_main_sc_lvlh_y' in pickle_list_name_temp ):
            algebric_distance_ensemble_main_sc_lvlh_y = pickle_list[ pickle_list_name_temp.index('algebric_distance_ensemble_main_sc_lvlh_y') ]
        if ( 'algebric_distance_ensemble_main_sc_lvlh_z' in pickle_list_name_temp ):
            algebric_distance_ensemble_main_sc_lvlh_z = pickle_list[ pickle_list_name_temp.index('algebric_distance_ensemble_main_sc_lvlh_z') ]
        if ( 'sma_ensemble' in pickle_list_name_temp ):
            sma_ensemble = pickle_list[ pickle_list_name_temp.index('sma_ensemble') ]
        if ( 'sma_ensemble_orbit_average' in pickle_list_name_temp ):
            sma_ensemble_orbit_average = pickle_list[ pickle_list_name_temp.index('sma_ensemble_orbit_average') ]    
        if ( 'phase_angle_ensemble_orbit_average' in pickle_list_name_temp ):
            phase_angle_ensemble_orbit_average = pickle_list[ pickle_list_name_temp.index('phase_angle_ensemble_orbit_average') ]    
        if ( 'sma_difference_ensemble_orbit_average' in pickle_list_name_temp ):
            sma_difference_ensemble_orbit_average = pickle_list[ pickle_list_name_temp.index('sma_difference_ensemble_orbit_average') ]    

        if ( 'period_ensemble' in pickle_list_name_temp ):
            period_ensemble = pickle_list[ pickle_list_name_temp.index('period_ensemble') ]
        if ( 'period_ensemble_orbit_average' in pickle_list_name_temp ):
            period_ensemble_orbit_average = pickle_list[ pickle_list_name_temp.index('period_ensemble_orbit_average') ]    
        if ( 'inclination_ensemble' in pickle_list_name_temp ):
            inclination_ensemble = pickle_list[ pickle_list_name_temp.index('inclination_ensemble') ]
        if ( 'eccentricity_ensemble' in pickle_list_name_temp ):
            eccentricity_ensemble = pickle_list[ pickle_list_name_temp.index('eccentricity_ensemble') ]
        if ( 'true_anomaly_ensemble' in pickle_list_name_temp ):
            true_anomaly_ensemble = pickle_list[ pickle_list_name_temp.index('true_anomaly_ensemble') ]
        if ( 'RAAN_ensemble' in pickle_list_name_temp ):
            RAAN_ensemble = pickle_list[ pickle_list_name_temp.index('RAAN_ensemble') ]
        if ( 'argument_perigee_ensemble' in pickle_list_name_temp ):
            argument_perigee_ensemble = pickle_list[ pickle_list_name_temp.index('argument_perigee_ensemble') ]
        if ( 'rho_ensemble' in pickle_list_name_temp ):
            rho_ensemble = pickle_list[ pickle_list_name_temp.index('rho_ensemble') ]    
        if ( 'rho_ensemble_orbit_average' in pickle_list_name_temp ):
            rho_ensemble_orbit_average = pickle_list[ pickle_list_name_temp.index('rho_ensemble_orbit_average') ]    
        if ( 'time_orbit_average' in pickle_list_name_temp ):
            time_orbit_average = pickle_list[ pickle_list_name_temp.index('time_orbit_average') ]    
        if ( 'index_time_orbit_average' in pickle_list_name_temp ):
            index_time_orbit_average = pickle_list[ pickle_list_name_temp.index('index_time_orbit_average') ]    
        if ( 'f107_ensemble' in pickle_list_name_temp ):
            f107_ensemble = pickle_list[ pickle_list_name_temp.index('f107_ensemble') ]
        if ( 'f107a_ensemble' in pickle_list_name_temp ):
            f107a_ensemble = pickle_list[ pickle_list_name_temp.index('f107a_ensemble') ]
        if ( 'ap_ensemble' in pickle_list_name_temp ):
            ap_ensemble = pickle_list[ pickle_list_name_temp.index('ap_ensemble') ]
        if ( 'angle_asc_node_to_sat_ensemble' in pickle_list_name_temp ):
            angle_asc_node_to_sat_ensemble = pickle_list[ pickle_list_name_temp.index('angle_asc_node_to_sat_ensemble') ][isat, :, :]
        if ( 'angular_distance_between_two_clusters' in pickle_list_name_temp ):
            angular_distance_between_two_clusters = pickle_list[ pickle_list_name_temp.index('angular_distance_between_two_clusters') ]
        if ( 'save_index_dt_angular_distance_between_two_clusters' in pickle_list_name_temp ):
            save_index_dt_angular_distance_between_two_clusters = pickle_list[ pickle_list_name_temp.index('save_index_dt_angular_distance_between_two_clusters') ]

        if ( 'angular_spacing_between_ensemble_sat' in pickle_list_name_temp ):
            angular_spacing_between_ensemble_sat = pickle_list[ pickle_list_name_temp.index('angular_spacing_between_ensemble_sat') ]
        if ( 'angular_spacing_between_ensemble_sat_converted_in_a_distance' in pickle_list_name_temp ):
            angular_spacing_between_ensemble_sat_converted_in_a_distance = pickle_list[ pickle_list_name_temp.index('angular_spacing_between_ensemble_sat_converted_in_a_distance') ]

        if ( 'tca_ensemble' in pickle_list_name_temp ):
            tca_ensemble = pickle_list[ pickle_list_name_temp.index('tca_ensemble') ]


        if ipickle == 0:
            nb_steps = pickle_list[1].shape[0] #!!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA

            
    eighty_quartile_along_run = []
    fifty_quartile_along_run = []
    std_quartile_along_run = []
    eighty_quartile_angle_run = []
    fifty_quartile_angle_run = []
    std_quartile_angle_run = []
    
    ########################################################################################################################################################################################################################################################################################################
    # PLOTS ################################################################################################################################################################################################################################################################################################
    ########################################################################################################################################################################################################################################################################################################

    # Parameters of figures
    width_fig = 15.
    height_fig = width_fig * 3 /4
    fontsize_plot = 20 
    when_plot_in_hour = 3*24 #5 + 50/60. # can be overwritten for each plot by uncommenting the same line in each plot section
    index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
    step_std = 1./60 # step in hours to calculate the standard deviation
    hour_time_step_xticks = 24. # time step of ticks when plotting a function as a function of time


    ###

    if index_when_plot > nb_steps:
        when_plot_in_hour = ( nb_steps - 1 ) * dt / 3600
        index_when_plot = (int) (when_plot_in_hour * 3600L / dt)

############################## ALONG-TRACK SEPARATION WITH REFERENCE SPACECRAFT
########################################################################################################################################################################################################################################################################################################
    if "along" in sys.argv:
        if ( 'algebric_distance_ensemble_main_sc_lvlh_x' in pickle_list_name ):
            # Plot the standard deviation of the angular_spacing_between_ensemble_sat_converted_in_a_distance distributions as a function of time

            if ( step_std < dt / 3600. ):
                step_std = dt / 3600.
            step_std_in_index = step_std * 3600. / dt
            nb_steps_adjusted = (int) ( nb_steps * dt / 3600 / step_std) 
            if  nb_steps_adjusted * step_std * 3600 < nb_steps * dt:
                nb_steps_adjusted = nb_steps_adjusted + 1
            if ( nb_steps_adjusted - 1 ) * step_std_in_index >= nb_steps:
                nb_steps_adjusted = nb_steps_adjusted - 1
            for i in range(nb_steps_adjusted):
                index_when_std = i * step_std_in_index
                eighty_quartile_along_run.append([np.subtract(*np.percentile(algebric_distance_ensemble_main_sc_lvlh_x[index_when_std, :], [90, 10]))])
                fifty_quartile_along_run.append([np.subtract(*np.percentile(algebric_distance_ensemble_main_sc_lvlh_x[index_when_std, :], [75, 25]))])
                std_quartile_along_run.append(np.std(algebric_distance_ensemble_main_sc_lvlh_x[index_when_std, :]))
            eighty_quartile_along.append(eighty_quartile_along_run)
            fifty_quartile_along.append(fifty_quartile_along_run)
            std_quartile_along.append(std_quartile_along_run)


    if "angle" in sys.argv:
        if ( 'angle_asc_node_to_sat_ensemble' in pickle_list_name ):
            # Plot the standard deviation of the angular_spacing_between_ensemble_sat_converted_in_a_distance distributions as a function of time

            if ( step_std < dt / 3600. ):
                step_std = dt / 3600.
            step_std_in_index = step_std * 3600. / dt
            nb_steps_adjusted = (int) ( nb_steps * dt / 3600 / step_std) 
            if  nb_steps_adjusted * step_std * 3600 < nb_steps * dt:
                nb_steps_adjusted = nb_steps_adjusted + 1
            if ( nb_steps_adjusted - 1 ) * step_std_in_index >= nb_steps:
                nb_steps_adjusted = nb_steps_adjusted - 1
            for i in range(nb_steps_adjusted):
                index_when_std = i * step_std_in_index
                eighty_quartile_angle_run.append([np.subtract(*np.percentile(angle_asc_node_to_sat_ensemble[index_when_std, :], [90, 10]))])
                fifty_quartile_angle_run.append([np.subtract(*np.percentile(angle_asc_node_to_sat_ensemble[index_when_std, :], [75, 25]))])
                std_quartile_angle_run.append(np.std(angle_asc_node_to_sat_ensemble[index_when_std, :]))
            eighty_quartile_angle.append(eighty_quartile_angle_run)
            fifty_quartile_angle.append(fifty_quartile_angle_run)
            std_quartile_angle.append(std_quartile_angle_run)

#angle_asc_node_to_sat_ensemble

    if "rho" in sys.argv:

        if ( 'rho_ensemble_orbit_average' in pickle_list_name ):
            # Plot the density of each ensemble as a function of time, as well as the median and the inter-quartile range
            if irun == nrun - 1:
                step_std_in_index = 1
                nb_steps_adjusted = (int) ( nb_steps * dt / 3600 / step_std) 
                if  nb_steps_adjusted * step_std * 3600 < nb_steps * dt:
                    nb_steps_adjusted = nb_steps_adjusted + 1
                if ( nb_steps_adjusted - 1 ) * step_std_in_index >= nb_steps:
                    nb_steps_adjusted = nb_steps_adjusted - 1

                nb_steps_orbit_average = len(index_time_orbit_average[:,1])     # !!!!!!!! there might be ore periods for satellites than for others (because they move at different speeds). So nb_steps_orbit_average should actually be the min of the number of periods of all satellites. Since we don't take the min, there might be 0 in rho_ensemble_orbit_average that we should NOT take into account. So basically NEVER TAKE THE STANDARD DEVIATION for this density section of the script (that takes into account outliers)
                index_time_orbit_average_middle = index_time_orbit_average[:,1]
                x_axis_orbit_average = index_time_orbit_average_middle * dt # conversion index to seconds
                x_axis_orbit_average = x_axis_orbit_average / (24*3600.) # conversion seconds to day

            rho_ensemble_orbit_average_all_run.append(rho_ensemble_orbit_average)

    if "f107" in sys.argv:

        if ( 'f107_ensemble' in pickle_list_name ):
            # Plot the density of each ensemble as a function of time, as well as the median and the inter-quartile range

            if ( step_std < dt / 3600. ):
                step_std = dt / 3600.
            step_std_in_index = step_std * 3600. / dt
            nb_steps_adjusted = (int) ( nb_steps * dt / 3600 / step_std) 
            if  nb_steps_adjusted * step_std * 3600 < nb_steps * dt:
                nb_steps_adjusted = nb_steps_adjusted + 1
            if ( nb_steps_adjusted - 1 ) * step_std_in_index >= nb_steps:
                nb_steps_adjusted = nb_steps_adjusted - 1

            f107_ensemble_all_run.append(f107_ensemble[:,0]) # all ensembles have same f107 so take the f107 of the first ensemble

    if "ap" in sys.argv:

        if ( 'ap_ensemble' in pickle_list_name ):
            # Plot the density of each ensemble as a function of time, as well as the median and the inter-quartile range

            if ( step_std < dt / 3600. ):
                step_std = dt / 3600.
            step_std_in_index = step_std * 3600. / dt
            nb_steps_adjusted = (int) ( nb_steps * dt / 3600 / step_std) 
            if  nb_steps_adjusted * step_std * 3600 < nb_steps * dt:
                nb_steps_adjusted = nb_steps_adjusted + 1
            if ( nb_steps_adjusted - 1 ) * step_std_in_index >= nb_steps:
                nb_steps_adjusted = nb_steps_adjusted - 1

            ap_ensemble_all_run.append(ap_ensemble[:,0]) # all ensembles have same ap so take the ap of the first ensembl



if "along" in sys.argv:
    eighty_quartile_along = np.array(eighty_quartile_along)[:,:,0]
    fifty_quartile_along = np.array(fifty_quartile_along)[:,:,0]
    std_quartile_along = np.array(std_quartile_along)


    ratio_fig_size = 4./3
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left=0.1, right=0.82, top = 0.90,bottom = 0.06)
    ax1 = fig.add_subplot(gs[0, 0])

    fifty_or_eighty = 80
    if fifty_or_eighty == 80:
        ax1.set_title('80% width of the along-track distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('80% width (km)', fontsize = fontsize_plot, weight = 'bold')
    if fifty_or_eighty == 50:
        ax1.set_title('50% width of the along-track distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('50% width (km)', fontsize = fontsize_plot, weight = 'bold')

    ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width o


    for irun in range(nrun):
        colorVal = scalarMap.to_rgba(values[irun])
        if fifty_or_eighty == 80:
            ax1.plot(np.arange(0,nb_steps_adjusted), eighty_quartile_along[irun, :] - eighty_quartile_along[median_run, :] , color = colorVal, linewidth = 2, label = run_list[irun].split('f107_')[1].split('_')[0]) 
        if fifty_or_eighty == 50:
            ax1.plot(np.arange(0,nb_steps_adjusted), fifty_quartile_along[irun, :] - fifty_quartile_along[median_run, :] , color = colorVal, linewidth = 2, label = run_list[irun].split('f107_')[1].split('_')[0]) 
#        ax1.plot(np.arange(0,nb_steps_adjusted), std_quartile_along[irun, :] - std_quartile_along[median_run, :] , color = colorVal, linewidth = 2, label = run_list[irun].split('f107_')[1].split('_')[0]) 

    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)


    hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / step_std
    xticks = np.arange(0, nb_steps_adjusted, hour_time_step_xticks_converted_in_index_adjusted)
    date_list_str = []
    nb_hours_simu = nb_steps * dt/ 3600.
    date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
    for i in range(len(xticks)):
        date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] * step_std / 24.)) + ")")
        ax1.xaxis.set_ticks(xticks)
        ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')

    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))#ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol = 3)
    if save_results == 1:
        fig_save_name = 'along_track_iqr_daily'
        fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
        fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
        os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' )






if "angle" in sys.argv:
    eighty_quartile_angle = np.array(eighty_quartile_angle)[:,:,0]
    fifty_quartile_angle = np.array(fifty_quartile_angle)[:,:,0]
    std_quartile_angle = np.array(std_quartile_angle)


    ratio_fig_size = 4./3
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left=0.1, right=0.82, top = 0.90,bottom = 0.06)
    ax1 = fig.add_subplot(gs[0, 0])

    fifty_or_eighty = 80
    if fifty_or_eighty == 80:
        ax1.set_title('80% width of the angular distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('80% width (' + u'\N{DEGREE SIGN}'+ ')', fontsize = fontsize_plot, weight = 'bold')
    if fifty_or_eighty == 50:
        ax1.set_title('50% width of the angular distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('50% width (' + u'\N{DEGREE SIGN}'+ ')', fontsize = fontsize_plot, weight = 'bold')

    ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width o


    for irun in range(nrun):
        colorVal = scalarMap.to_rgba(values[irun])
        # if fifty_or_eighty == 80:
        #     ax1.plot(np.arange(0,nb_steps_adjusted), eighty_quartile_angle[irun, :] - eighty_quartile_angle[median_run, :] , color = colorVal, linewidth = 2, label = run_list[irun].split('f107_')[1].split('_')[0]) 
        # if fifty_or_eighty == 50:
        #     ax1.plot(np.arange(0,nb_steps_adjusted), fifty_quartile_angle[irun, :] - fifty_quartile_angle[median_run, :] , color = colorVal, linewidth = 2, label = run_list[irun].split('f107_')[1].split('_')[0]) 
        ax1.plot(np.arange(0,nb_steps_adjusted), std_quartile_angle[irun, :] - std_quartile_angle[median_run, :] , color = colorVal, linewidth = 2, label = run_list[irun].split('f107_')[1].split('_')[0]) 

    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)


    hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / step_std
    xticks = np.arange(0, nb_steps_adjusted, hour_time_step_xticks_converted_in_index_adjusted)
    date_list_str = []
    nb_hours_simu = nb_steps * dt/ 3600.
    date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
    for i in range(len(xticks)):
        date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] * step_std / 24.)) + ")")
        ax1.xaxis.set_ticks(xticks)
        ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')

    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))#ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol = 3)
    if save_results == 1:
        fig_save_name = 'angle_track_iqr_daily'
        fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
        fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
        os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' )


        # if save_results == 1:
        #     fig_save_name = 'angle_track_iqr_daily'
        #     fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
        #     fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
        #     os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)






############################## DENSITY AS A FUNCTION OF TIME
########################################################################################################################################################################################################################################################################################################
if "rho" in sys.argv:

    if ( 'rho_ensemble_orbit_average' in pickle_list_name ):
        # Plot the density of each ensemble as a function of time, as well as the median and the inter-quartile range
        rho_ensemble_orbit_average_all_run = np.array( rho_ensemble_orbit_average_all_run)[:,:,0]

        step_std_in_index = 1 # every orbit
        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left=0.1, right=0.82, top = 0.90,bottom = 0.06)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('')#Density at satellite position VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('Density (#/m^3)', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        # label_array = ['With storm', 'Without storm']
        # color_array = ['r','b']
        for irun in range(nrun):
            colorVal = scalarMap.to_rgba(values[irun]) #color_array[irun]
            ax1.plot(x_axis_orbit_average, rho_ensemble_orbit_average_all_run[irun, 0:nb_steps:step_std_in_index], color = colorVal, linewidth = 2, label =  run_list[irun].split('f107_')[1].split('_')[0]+ "0%") #label_array[irun])#



        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)


        nb_days_simu = ( date_stop - date_start ).total_seconds()/3600./24.
        nb_hours_simu = nb_steps * dt/ 3600.
        hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / 24. # x_axis_orbit_average is in unit of days
        xticks = np.arange(0, nb_days_simu, hour_time_step_xticks_converted_in_index_adjusted)
        date_list_str = []
        date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
        for i in range(len(xticks)):
            date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] )) + ")")
            ax1.xaxis.set_ticks(xticks)
            ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')


        legend = ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), title="Quartile of\nF10.7/Ap dist.", fontsize = fontsize_plot)
        legend.get_title().set_fontsize(str(fontsize_plot))
#        ax1.legend()

        if save_results == 1:
            fig_save_name = 'rho_vs_time'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' )

        

    ############################## F10.7
    ########################################################################################################################################################################################################################################################################################################
if "f107" in sys.argv:
    if ( 'f107_ensemble' in pickle_list_name ):
        # Plot the f10.7 for each ensemble as a function of time, as well as the median and the inter-quartile range
 
        f107_ensemble_all_run = np.array( f107_ensemble_all_run)


        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left=0.1, right=0.82, top = 0.90,bottom = 0.06)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('F10.7 VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('F10.7', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        for irun in range(nrun):
            colorVal = scalarMap.to_rgba(values[irun])
#            ax1.plot(x_axis_orbit_average, rho_ensemble_orbit_average_all_run[irun, 0:nb_steps:step_std_in_index], color = colorVal, linewidth = 2, label = run_list[irun].split('f107_')[1].split('_')[0])

            ax1.plot(np.arange(0,nb_steps_adjusted), f107_ensemble_all_run[irun,0:nb_steps:step_std_in_index], color = colorVal, linewidth = 2, label = run_list[irun].split('f107_')[1].split('_')[0]+ "0%")
        #ax1.plot(np.arange(0,nb_steps_adjusted), med_daily[:], 'k', linewidth = 3, label = 'Median')
        #ax1.plot(np.arange(0,nb_steps_adjusted), quartile25[:], 'r', linewidth = 3, label = '25 and 75% quartiles')
        #ax1.plot(np.arange(0,nb_steps_adjusted), quartile75[:], 'r', linewidth = 3)

        legend = ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), title="Quartile of\nF10.7 dist.", fontsize = fontsize_plot)
        legend.get_title().set_fontsize(str(fontsize_plot))
        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)


        
        hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / step_std
        xticks = np.arange(0, nb_steps_adjusted, hour_time_step_xticks_converted_in_index_adjusted)
        date_list_str = []
        nb_hours_simu = nb_steps * dt/ 3600.
        date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
        for i in range(len(xticks)):
            date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] * step_std / 24.)) + ")")
        ax1.xaxis.set_ticks(xticks)
        ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')

#        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))#ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol = 3)


        if save_results == 1:
            fig_save_name = 'f107_vs_time'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' )


    ############################## Ap
    ########################################################################################################################################################################################################################################################################################################
if "ap" in sys.argv:
    if ( 'ap_ensemble' in pickle_list_name ):
        # Plot the f10.7 for each ensemble as a function of time, as well as the median and the inter-quartile range
 
        ap_ensemble_all_run = np.array( ap_ensemble_all_run)


        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left=0.1, right=0.82, top = 0.90,bottom = 0.06)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Ap VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('Ap', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        for irun in range(nrun):
            colorVal = scalarMap.to_rgba(values[irun])
#            ax1.plot(x_axis_orbit_average, rho_ensemble_orbit_average_all_run[irun, 0:nb_steps:step_std_in_index], color = colorVal, linewidth = 2, label = run_list[irun].split('ap_')[1].split('_')[0])

            ax1.plot(np.arange(0,nb_steps_adjusted), ap_ensemble_all_run[irun,0:nb_steps:step_std_in_index], color = colorVal, linewidth = 2, label = run_list[irun].split('f107_')[1].split('_')[0]+ "0%")
        #ax1.plot(np.arange(0,nb_steps_adjusted), med_daily[:], 'k', linewidth = 3, label = 'Median')
        #ax1.plot(np.arange(0,nb_steps_adjusted), quartile25[:], 'r', linewidth = 3, label = '25 and 75% quartiles')
        #ax1.plot(np.arange(0,nb_steps_adjusted), quartile75[:], 'r', linewidth = 3)

        legend = ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), title="Quartile of\nAp dist.", fontsize = fontsize_plot)
        legend.get_title().set_fontsize(str(fontsize_plot))
        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)


        
        hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / step_std
        xticks = np.arange(0, nb_steps_adjusted, hour_time_step_xticks_converted_in_index_adjusted)
        date_list_str = []
        nb_hours_simu = nb_steps * dt/ 3600.
        date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
        for i in range(len(xticks)):
            date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] * step_std / 24.)) + ")")
        ax1.xaxis.set_ticks(xticks)
        ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')

#        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))#ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol = 3)

        if save_results == 1:
            fig_save_name = 'ap_vs_time'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' )

if show_plots == 1:
    plt.show(); plt.show()

# Arrays of quartiles to remove from the nominal predicted values in swpc_predictions_f107_ap.txt:
f107_quartile_array  = np.array([[-4. , -3. , -2. , -0.2,  0. ,  1. ,  2. ,  3. ,  4.8],
                             [-6. , -4. , -2.7, -1. ,  0. ,  1. ,  2.7,  4.8,  7. ],
                             [-8. , -5. , -3. , -2. ,  0. ,  1. ,  3. ,  5. ,  8. ]])

ap_quartile_array = np.array([[-7. , -3. , -2. ,  0. ,  1. ,  2. ,  3. ,  6. ,  9.8],
                              [-9. , -4. , -2. ,  0. ,  1. ,  2. ,  3. ,  5. ,  8.9],
                              [-9. , -5. , -2. ,  0. ,  1. ,  1. ,  3. ,  5. ,  8. ]])


# f107 ap at [d+1  d+2  d+3] (values at d+0: 82 15) (equal to nominal values from swpc_predictions_f107_ap.txt plus values from f107_quartile_array and ap_quartile_array)
# 10% [ 79.  77.  73.] [ 3. -1. -1.]
# 20% [ 80.  79.  76.] [ 7.  4.  3.]
# 30% [ 81.   80.3  78. ] [ 8.  6.  6.]
# 40% [ 82.8  82.   79. ] [ 10.   8.   8.]
# 50% [ 83.  83.  81.] [ 11.   9.   9.]
# 60% [ 84.  84.  82.] [ 12.  10.   9.]
# 70% [ 85.   85.7  84. ] [ 13.  11.  11.]
# 80% [ 86.   87.8  86. ] [ 16.  13.  13.]
# 90% [ 87.8  90.   89. ] [ 19.8  16.9  16. ]
