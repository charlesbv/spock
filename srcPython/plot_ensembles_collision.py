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

## NOTE 0: to run this script, you first need to run distance_ensemble_to_main_sc.py (with first_time = 1). This will create the data that plot_ensembles will then use to makes plots
## NOTE 1: to use this script, the only 5 parameters you have to set are:
## - if yes or no you want to save the plots (save_results = 1 if yes, 0 otherwise)
## - if yes or no you want to show the plots (show_plots = 1 if yes, 0 otherwise)
## - which reference satellite (and its ensembles) you want to look at (recall: the number of reference satellites is indicated on the first line of section #SPACECRAFT of SpOCK main input file). Choose this satellite in the variable isat
## - the name of the propagator main input file 
## - the path of the folder where you want to store the results (pickle, image, video): called 'path_folder_results'. In this folder, there must be the following subfolders: 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', and 'other'. In each of these subfolders, there must be the 2 subsubfolders: 'result', and 'pickle'. In the subsubfolder 'result', there must be the 2 subsubsubfolders: 'image', and 'video'.
## NOTE 2: this can be run if only ONE MAIN satellite was run (with ensembles of course)
## NOTE 3: the results (pickle, image, video) are stored in a subfolder of path_folder_results. The subfolder is either 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', or 'other'. This code figures out which folder the simulation corresponds to by reading the name of the output file chosen by the user in the main input file of the propagator (third line of section #SPACECRAFTS): it tries to find 'cygnss', 'cadre', 'aerie', 'scion', or 'qb50' in the name of the output file. If it does not find it, then the results here will be stored in path_folder_results/other/
## NOTE 4: to run this script, and any python script in the propagator, you need to be one subfolder deep from the main folder where the propagator runs are made. So if path_to_propagator/PropSim is the folder where the propagator runs are made, then the python scripts must be run for example in path_to_propagator/PropSim/subfolder_where_python_scipts_are_run


# !!!!!!!!!! SET THE PARAMETER BELOW:
## Save or not the plots
save_results = 0

## Show or not the plots
show_plots = 1

## Which reference satellite to look at (starts at 0 and ends and nb_spacecraft - 1)
isat = 0

## path of the folder where you want to store the results (pickle, image, video)
path_folder_results = '/raid3/Armada/Charles/python/' #get_prop_dir(2) + 'output/python_propagator/'

## main input file (argument in command line)
if len(sys.argv) > 2:
    main_input_file_name = get_prop_dir(1) + sys.argv[1] + '/input/main_input/' + sys.argv[2]  # ex of sys.argv[1]: cadre_tumbling_0panel.txt'
else:
    main_input_file_name = get_prop_dir(1) + 'run/input/main_input/' + sys.argv[1]  # ex of sys.argv[1]: cadre_tumbling_0panel.txt'    


# read input file
input_variables, order_input_variables = read_input_file(main_input_file_name)
date_start = input_variables[0]; date_stop = input_variables[1]; dt = input_variables[2]; nb_steps = input_variables[3]; satellite_to_plot_path = input_variables[6]; satellite_to_plot = input_variables[7]; nb_ensembles_coe = input_variables[12]; nb_ensembles_attitude = input_variables[13]; nb_ensembles_cd = input_variables[11]; ensemble_to_plot_temp = input_variables[14]; 
nb_ensembles_density = input_variables[17]
n = nb_steps
nb_spacecraft = input_variables[4]
compute_drag = input_variables[19]


# set up interactive figures
if show_plots == 1:
    plt.ion()

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
        
if ((compute_drag == 0) & (ensemble_to_plot_temp[i] == 'density')):
    print "Note: in the main input file of SpOCK (" + main_input_file_name + "), you chose to output ensembles for the density. However, the atmospheric drag was not computed. Therefore, the density is not computed either"
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
if ( save_results == 1 ):
    os.system("mkdir " + path_folder_results + name_mission + '/result/image/' + name_subfolder_save )
    os.system('ssh -t srbwks2014-0008.engin.umich.edu "mkdir ' + name_mission + '/' + name_subfolder_save + '"')
#    os.system("mkdir " + path_folder_results + name_mission + '/result/video/' + name_subfolder_save )

    
save_pickle_name = path_folder_results + name_mission + '/pickle/' + name_subfolder_save + satellite_to_plot[isat].replace(".txt",".pickle")
root_save_fig_name = path_folder_results + name_mission + '/result/image/' + name_subfolder_save + satellite_to_plot[isat].replace(".txt","_") 
root_save_video_name = path_folder_results + name_mission + '/result/video/' + name_subfolder_save + satellite_to_plot[isat].replace(".txt","_")



# # #######################################
dir_final_output_ensemble = satellite_to_plot_path[isat] + 'ensemble/'
## RESTORE THE RESULTS FROM THE PICKLE
list_save_pickle = [f for f in os.listdir(path_folder_results + name_mission + '/pickle/' + name_subfolder_save) if satellite_to_plot[isat].replace(".txt", "") in f]
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
        angle_asc_node_to_sat_ensemble = pickle_list[ pickle_list_name_temp.index('angle_asc_node_to_sat_ensemble') ]
    if ( 'angular_distance_between_two_clusters' in pickle_list_name_temp ):
        angular_distance_between_two_clusters = pickle_list[ pickle_list_name_temp.index('angular_distance_between_two_clusters') ]
    if ( 'save_index_dt_angular_distance_between_two_clusters' in pickle_list_name_temp ):
        save_index_dt_angular_distance_between_two_clusters = pickle_list[ pickle_list_name_temp.index('save_index_dt_angular_distance_between_two_clusters') ]

    if ( 'angular_spacing_between_ensemble_sat' in pickle_list_name_temp ):
        angular_spacing_between_ensemble_sat = pickle_list[ pickle_list_name_temp.index('angular_spacing_between_ensemble_sat') ]
    if ( 'angular_spacing_between_ensemble_sat_converted_in_a_distance' in pickle_list_name_temp ):
        angular_spacing_between_ensemble_sat_converted_in_a_distance = pickle_list[ pickle_list_name_temp.index('angular_spacing_between_ensemble_sat_converted_in_a_distance') ]

    if ipickle == 0:
        nb_steps = pickle_list[1].shape[0] #!!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA

########################################################################################################################################################################################################################################################################################################
# PLOTS ################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################

# Parameters of figures
height_fig = 9.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
when_plot_in_hour = 7*24 #5 + 50/60. # can be overwritten for each plot by uncommenting the same line in each plot section
index_when_plot = 0 # (int) (when_plot_in_hour * 3600L / dt)
step_std = 1./60 # step in hours to calculate the standard deviation
hour_time_step_xticks = 24. # time step of ticks when plotting a function as a function of time


############################## ALONG-TRACK (USING ANGULAR SEPARATION)
########################################################################################################################################################################################################################################################################################################
if "angle" in sys.argv:
    if ( 'angle_asc_node_to_sat_ensemble' in pickle_list_name ):
        # Plot the distribution of angular_spacing_between_ensemble_sat_converted_in_a_distance at different times (set by when_plot_in_hour)
        # when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
        

        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(4, 1)
        gs.update(left= 0.12, right=0.90, top = 0.93,bottom = 0.08, hspace = 0.01)

        med = np.median(angle_asc_node_to_sat_ensemble[isat,index_when_plot, :])
        mad = np.median( np.abs( angle_asc_node_to_sat_ensemble[isat,index_when_plot, :] - np.median(angle_asc_node_to_sat_ensemble[isat,index_when_plot, :]) ) )  
        quartile10 = np.percentile(angle_asc_node_to_sat_ensemble[isat,index_when_plot, :], 10) 
        quartile25 = np.percentile(angle_asc_node_to_sat_ensemble[isat,index_when_plot, :], 25) 
        quartile75 = np.percentile(angle_asc_node_to_sat_ensemble[isat,index_when_plot, :], 75) 
        quartile90 = np.percentile(angle_asc_node_to_sat_ensemble[isat,index_when_plot, :], 90) 
        iqr_daily = quartile75 - quartile25
        quartiles_1090_daily = quartile90 - quartile10

        ax2 = fig.add_subplot(gs[:3, 0])
        ax2.set_title('Spacecraft distribution along the orbit (angle from AN to spacecraft) (sat '+ str(isat) + ')', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
        ax2.set_xlabel('', fontsize = fontsize_plot, weight = 'bold')
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
        bin_width = iqr_daily/10.
        hist_data = angle_asc_node_to_sat_ensemble[isat,index_when_plot, :]
        n, bins, patches = ax2.hist(hist_data , bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'r',label = '')  
        plt.vlines(quartile25, min(n), max(n),linewidth = 2); plt.vlines(quartile75, min(n),max(n), linewidth = 2) 
        plt.vlines(quartile10, min(n), max(n), linewidth = 4, linestyle = 'dotted'); plt.vlines(quartile90, min(n), max(n), linewidth = 4, linestyle = 'dotted') 
        ax2.plot(med, min(n), 'k',marker = '.', markersize = 20)
        ax2.get_xaxis().get_major_formatter().set_useOffset(False)
        ax2.set_xlim([med - 10*iqr_daily, med+10*iqr_daily])
        ax2.get_xaxis().set_ticklabels([])
        ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

        ax1 = fig.add_subplot(gs[3, 0])
        ymin = -0.2; ymax = 0.3
        length_vert_iqr_daily = np.abs(ymin) / 3.
        ax1.set_title('', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Angular distance from the AN (' + u'\N{DEGREE SIGN}'+ ')', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

        ax1.hlines(0, med - 10 * iqr_daily, med + 10 * iqr_daily, linewidth = 2, color = 'b');
        ax1.text( med + 10 * iqr_daily, 0 - length_vert_iqr_daily/0.8,'Orbit', horizontalalignment = 'right', fontsize = fontsize_plot, color = 'b' )

        ax1.plot(angle_asc_node_to_sat_ensemble[isat, index_when_plot, :], np.zeros([nb_ensembles]), 'r.', markersize = 10)
        ax1.vlines(quartile25, 0,length_vert_iqr_daily, linewidth = 2); ax1.vlines(quartile75, 0,length_vert_iqr_daily, linewidth = 2) 
        ax1.text((quartile75 + quartile25 ) /2., length_vert_iqr_daily+length_vert_iqr_daily/5, '{0:.2f}'.format(iqr_daily * 110) + ' km (50% of sc)', horizontalalignment = 'center', fontsize = fontsize_plot )
        ax1.annotate ('', (quartile25, length_vert_iqr_daily), (quartile75, length_vert_iqr_daily), arrowprops={'arrowstyle':'<->'})
        ax1.vlines(quartile10, -length_vert_iqr_daily ,0 , linewidth = 4, linestyle = 'dotted'); ax1.vlines(quartile90, -length_vert_iqr_daily ,0, linewidth = 4, linestyle = 'dotted') 
        ax1.annotate ('', (quartile10, -length_vert_iqr_daily), (quartile90, -length_vert_iqr_daily), arrowprops={'arrowstyle':'<->'})
        ax1.text((quartile10 + quartile90 ) /2., -length_vert_iqr_daily-length_vert_iqr_daily/0.8, '{0:.2f}'.format(quartiles_1090_daily * 110) + ' km (80% of sc)', horizontalalignment = 'center', fontsize = fontsize_plot, verticalalignment = 'bottom' )
        ax1.plot(med, 0, 'k',marker = '.', markersize = 15)

        ax1.legend(fontsize = fontsize_plot, loc = 2, scatterpoints=1)
        ax1.get_xaxis().tick_bottom()
        ax1.set_ylim([ymin, ymax])
    #    ax1.get_yaxis().tick_left()
        ax1.yaxis.set_visible(False)
        ax1.spines['left'].set_visible(False); ax1.spines['right'].set_visible(False); ax1.spines['top'].set_visible(False)
        ax1.get_xaxis().get_major_formatter().set_useOffset(False)
        ax1.set_xlim([med - 10*iqr_daily, med+10*iqr_daily])
    #    ax1.set_xlim([min(bins), 0])
        ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

        if save_results == 1:
            fig_save_name = 'sc_along_track_distribution_using_angular_separation_sat_' + str(isat)
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

    ############################## SMA
    ########################################################################################################################################################################################################################################################################################################
if "sma" in sys.argv:
    if ( 'sma_ensemble' in pickle_list_name ):
        # Plot the distribution of the f10.7 at different times (set by when_plot_in_hour)
        # when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
                  

        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

        med = np.median(sma_ensemble[index_when_plot, :])
        mad = np.median( np.abs( sma_ensemble[index_when_plot, :] - np.median(sma_ensemble[index_when_plot, :]) ) )  
        quartile10 = np.percentile(sma_ensemble[index_when_plot, :], 10) 
        quartile25 = np.percentile(sma_ensemble[index_when_plot, :], 25) 
        quartile75 = np.percentile(sma_ensemble[index_when_plot, :], 75) 
        quartile90 = np.percentile(sma_ensemble[index_when_plot, :], 90) 
        std = np.std(sma_ensemble[index_when_plot, :])
        iqr_daily = quartile75 - quartile25
        quartiles_1090_daily = quartile90 - quartile10

        ax2 = fig.add_subplot(gs[0, 0])
        ax2.set_title('Sma distribution ' + 'after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
        bin_width = iqr_daily/10.
        hist_data = sma_ensemble[index_when_plot, :]
        n, bins, patches = ax2.hist(hist_data , bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'r',label = '')  
        plt.vlines(quartile25, 0, max(n),linewidth = 2); plt.vlines(quartile75, 0,max(n), linewidth = 2, label = '25-75% quartiles') 
        plt.vlines(quartile10, 0, max(n), linewidth = 4, linestyle = 'dotted'); plt.vlines(quartile90, 0, max(n), linewidth = 4, linestyle = 'dotted', label = '10-90% deciles') 
        ax2.plot(med, 0, 'k',marker = '.', markersize = 20)
        ax2.get_xaxis().get_major_formatter().set_useOffset(False)
        ax2.set_xlim([med - 10*iqr_daily, med+10*iqr_daily])
        ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
        ax2.text(ax2.get_xlim()[0] + ( ax2.get_xlim()[1] - ax2.get_xlim()[0] )/ 60., ax2.get_ylim()[1] - (ax2.get_ylim()[1] - ax2.get_ylim()[0])/5., 'std = ' + '{0:.2f}'.format(std), fontsize = fontsize_plot)
        ax2.set_xlabel('Sma', fontsize = fontsize_plot, weight = 'bold')

        ax2.legend(fontsize = fontsize_plot, loc = 3, scatterpoints=1)
        if save_results == 1:
            fig_save_name = 'sma_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)
        
if show_plots == 1:
    plt.show(); plt.show()
