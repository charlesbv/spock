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
## NOTE 1: to use this script, the only 3 parameters you have to set are:
## - if yes or no you want to save the plots (save_results = 1 if yes, 0 otherwise)
## - if yes or no you want to show the plots (show_plots = 1 if yes, 0 otherwise)
## - the name of the propagator main input file 
## - the path of the folder where you want to store the results (pickle, image, video): called 'path_folder_results'. In this folder, there must be the following subfolders: 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', and 'other'. In each of these subfolders, there must be the 2 subsubfolders: 'result', and 'pickle'. In the subsubfolder 'result', there must be the 2 subsubsubfolders: 'image', and 'video'.
## NOTE 2: this can be run if only ONE MAIN satellite was run (with ensembles of course)
## NOTE 3: the results (pickle, image, video) are stored in a subfolder of path_folder_results. The subfolder is either 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', or 'other'. This code figures out which folder the simulation corresponds to by reading the name of the output file chosen by the user in the main input file of the propagator (third line of section #SPACECRAFTS): it tries to find 'cygnss', 'cadre', 'aerie', 'scion', or 'qb50' in the name of the output file. If it does not find it, then the results here will be stored in path_folder_results/other/
## NOTE 4: to run this script, and any python script in the propagator, you need to be one subfolder deep from the main folder where the propagator runs are made. So if path_to_propagator/PropSim is the folder where the propagator runs are made, then the python scripts must be run for example in path_to_propagator/PropSim/subfolder_where_python_scipts_are_run


# !!!!!!!!!! SET THE PARAMETER BELOW:
## Save or not the plots
save_results = 1

## Show or not the plots
show_plots = 0

## path of the folder where you want to store the results (pickle, image, video)
path_folder_results = '/raid3/Armada/Charles/python/' #get_prop_dir(2) + 'output/python_propagator/'

## main input file (argument in command line)
if len(sys.argv) > 2:
    main_input_file_name = get_prop_dir(1) + sys.argv[1] + '/input/main_input/' + sys.argv[2]  # ex of sys.argv[1]: cadre_tumbling_0panel.txt'
else:
    main_input_file_name = get_prop_dir(1) + 'run/input/main_input/' + sys.argv[1]  # ex of sys.argv[1]: cadre_tumbling_0panel.txt'    


# read input file
input_variables, order_input_variables = read_input_file(main_input_file_name)
dt = input_variables[2]; nb_steps = input_variables[3]; satellite_to_plot_path = input_variables[6][0]; satellite_to_plot = input_variables[7][0]; nb_ensembles_coe = input_variables[12]; nb_ensembles_attitude = input_variables[13]; nb_ensembles_cd = input_variables[11]; ensemble_to_plot_temp = input_variables[14]; 
nb_ensembles_density = input_variables[17]
n = nb_steps

# set up interactive figures
if show_plots == 1:
    plt.ion()

# Name of the mission of the simulation (choices: 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', 'other')
if ( 'cygnss' in satellite_to_plot.lower() ):
    name_mission = 'CYGNSS' 
elif ( 'cadre' in satellite_to_plot.lower() ):
    name_mission = 'CADRE' 
elif ( 'aerie' in satellite_to_plot.lower() ):
    name_mission = 'AERIE' 
elif ( 'scion' in satellite_to_plot.lower() ):
    name_mission = 'SCION' 
elif ( 'qb50' in satellite_to_plot.lower() ):
    name_mission = 'QB50' 
else:
    name_mission = 'other' 


    
name_subfolder_save = satellite_to_plot[:-5] + "/"
if ( save_results == 1 ):
    os.system("mkdir " + path_folder_results + name_mission + '/result/image/' + name_subfolder_save )
    os.system('ssh -t srbwks2014-0008.engin.umich.edu "mkdir ' + name_mission + '/' + name_subfolder_save + '"')
#    os.system("mkdir " + path_folder_results + name_mission + '/result/video/' + name_subfolder_save )

    
save_pickle_name = path_folder_results + name_mission + '/pickle/' + name_subfolder_save + satellite_to_plot.replace(".txt",".pickle")
root_save_fig_name = path_folder_results + name_mission + '/result/image/' + name_subfolder_save + satellite_to_plot.replace(".txt","_") 
root_save_video_name = path_folder_results + name_mission + '/result/video/' + name_subfolder_save + satellite_to_plot.replace(".txt","_")


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
        ensemble_to_plot.append('sma'); ensemble_to_plot.append('inclination'); ensemble_to_plot.append('eccentricity'); ensemble_to_plot.append('true_anomaly'); ensemble_to_plot.append('RAAN'); ensemble_to_plot.append('argument_perigee');
    if (ensemble_to_plot_temp[i] == 'density'):
        ensemble_to_plot.append('rho'); ensemble_to_plot.append('f107'); ensemble_to_plot.append('f107a'); ensemble_to_plot.append('ap'); 

pickle_list = []
pickle_list_name = []
## Nb of ensembles
nb_spacecraft = 1#int(a_input_file_input[6][0:6])
nb_ensembles_array = [nb_ensembles_coe, nb_ensembles_attitude, nb_ensembles_cd, nb_ensembles_density]
nb_ensembles = np.max(nb_ensembles_array)
for i in range(len(nb_ensembles_array)):
    if ( (nb_ensembles_array[i] > 0 ) &  (nb_ensembles_array[i] < nb_ensembles) ) :
        nb_ensembles = nb_ensembles_array[i]

# # #######################################
dir_final_output_ensemble = satellite_to_plot_path + 'ensemble/'
## RESTORE THE RESULTS FROM THE PICKLE
with open(save_pickle_name) as f:
    pickle_list, pickle_list_name = pickle.load(f)
if ( 'x_eci_ensemble' in pickle_list_name ): 
    x_eci_ensemble = pickle_list[ pickle_list_name.index('x_eci_ensemble') ]
if ( 'y_eci_ensemble' in pickle_list_name ): 
    y_eci_ensemble = pickle_list[ pickle_list_name.index('y_eci_ensemble') ]
if ( 'z_eci_ensemble' in pickle_list_name ): 
    z_eci_ensemble = pickle_list[ pickle_list_name.index('z_eci_ensemble') ]
if ( 'pitch_ensemble' in pickle_list_name ): 
    pitch_ensemble = pickle_list[ pickle_list_name.index('pitch_ensemble') ]
if ( 'roll_ensemble' in pickle_list_name ): 
    roll_ensemble = pickle_list[ pickle_list_name.index('roll_ensemble') ]
if ( 'yaw_ensemble' in pickle_list_name ): 
    yaw_ensemble = pickle_list[ pickle_list_name.index('yaw_ensemble') ]
if ( 'x_eci_main_sc' in pickle_list_name ):
    x_eci_main_sc = pickle_list[ pickle_list_name.index('x_eci_main_sc') ]
if ( 'y_eci_main_sc' in pickle_list_name ):
    y_eci_main_sc = pickle_list[ pickle_list_name.index('y_eci_main_sc') ]
if ( 'z_eci_main_sc' in pickle_list_name ):
    z_eci_main_sc = pickle_list[ pickle_list_name.index('z_eci_main_sc') ]
if ( 'vx_eci_main_sc' in pickle_list_name ):
    vx_eci_main_sc = pickle_list[ pickle_list_name.index('vx_eci_main_sc') ]
if ( 'vy_eci_main_sc' in pickle_list_name ):
    vy_eci_main_sc = pickle_list[ pickle_list_name.index('vy_eci_main_sc') ]
if ( 'vz_eci_main_sc' in pickle_list_name ):
    vz_eci_main_sc = pickle_list[ pickle_list_name.index('vz_eci_main_sc') ]
if ( 'algebric_distance_ensemble_main_sc_lvlh_x' in pickle_list_name ):
    algebric_distance_ensemble_main_sc_lvlh_x = pickle_list[ pickle_list_name.index('algebric_distance_ensemble_main_sc_lvlh_x') ]
if ( 'algebric_distance_ensemble_main_sc_lvlh_y' in pickle_list_name ):
    algebric_distance_ensemble_main_sc_lvlh_y = pickle_list[ pickle_list_name.index('algebric_distance_ensemble_main_sc_lvlh_y') ]
if ( 'algebric_distance_ensemble_main_sc_lvlh_z' in pickle_list_name ):
    algebric_distance_ensemble_main_sc_lvlh_z = pickle_list[ pickle_list_name.index('algebric_distance_ensemble_main_sc_lvlh_z') ]
if ( 'sma_ensemble' in pickle_list_name ):
    sma_ensemble = pickle_list[ pickle_list_name.index('sma_ensemble') ]
if ( 'inclination_ensemble' in pickle_list_name ):
    inclination_ensemble = pickle_list[ pickle_list_name.index('inclination_ensemble') ]
if ( 'eccentricity_ensemble' in pickle_list_name ):
    eccentricity_ensemble = pickle_list[ pickle_list_name.index('eccentricity_ensemble') ]
if ( 'true_anomaly_ensemble' in pickle_list_name ):
    true_anomaly_ensemble = pickle_list[ pickle_list_name.index('true_anomaly_ensemble') ]
if ( 'RAAN_ensemble' in pickle_list_name ):
    RAAN_ensemble = pickle_list[ pickle_list_name.index('RAAN_ensemble') ]
if ( 'argument_perigee_ensemble' in pickle_list_name ):
    argument_perigee_ensemble = pickle_list[ pickle_list_name.index('argument_perigee_ensemble') ]
if ( 'rho_ensemble' in pickle_list_name ):
    rho_ensemble = pickle_list[ pickle_list_name.index('rho_ensemble') ]    
if ( 'rho_ensemble_orbit_average' in pickle_list_name ):
    rho_ensemble_orbit_average = pickle_list[ pickle_list_name.index('rho_ensemble_orbit_average') ]    
if ( 'time_orbit_average' in pickle_list_name ):
    time_orbit_average = pickle_list[ pickle_list_name.index('time_orbit_average') ]    
if ( 'index_time_orbit_average' in pickle_list_name ):
    index_time_orbit_average = pickle_list[ pickle_list_name.index('index_time_orbit_average') ]    
if ( 'f107_ensemble' in pickle_list_name ):
    f107_ensemble = pickle_list[ pickle_list_name.index('f107_ensemble') ]
if ( 'f107a_ensemble' in pickle_list_name ):
    f107a_ensemble = pickle_list[ pickle_list_name.index('f107a_ensemble') ]
if ( 'ap_ensemble' in pickle_list_name ):
    ap_ensemble = pickle_list[ pickle_list_name.index('ap_ensemble') ]
if ( 'angle_asc_node_to_sat_ensemble' in pickle_list_name ):
    angle_asc_node_to_sat_ensemble = pickle_list[ pickle_list_name.index('angle_asc_node_to_sat_ensemble') ]
if ( 'angular_spacing_between_ensemble_sat' in pickle_list_name ):
    angular_spacing_between_ensemble_sat = pickle_list[ pickle_list_name.index('angular_spacing_between_ensemble_sat') ]
if ( 'angular_spacing_between_ensemble_sat_converted_in_a_distance' in pickle_list_name ):
    angular_spacing_between_ensemble_sat_converted_in_a_distance = pickle_list[ pickle_list_name.index('angular_spacing_between_ensemble_sat_converted_in_a_distance') ]


r_eci_ensemble = np.zeros([nb_steps, nb_ensembles])
for i in range(nb_steps):
    for  j in range(nb_ensembles):
        r_eci_ensemble[i,j] = np.sqrt( x_eci_ensemble[i,j]**2 + y_eci_ensemble[i,j]**2 + z_eci_ensemble[i,j]**2 )


########################################################################################################################################################################################################################################################################################################
# PLOTS ################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################

# Parameters of figures
height_fig = 9.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 



############################## ALONG-TRACK SEPARATION WITH REFERENCE SPACECRAFT
########################################################################################################################################################################################################################################################################################################
if "along" in sys.argv:
    if ( 'algebric_distance_ensemble_main_sc_lvlh_x' in pickle_list_name ):
        # Plot the standard deviation of the angular_spacing_between_ensemble_sat_converted_in_a_distance distributions as a function of time
        step_std = 10./60 # step in hours to calculate the standard deviation
        if ( step_std < dt / 3600. ):
            step_std = dt / 3600.
        step_std_in_index = step_std * 3600. / dt
        nb_steps_adjusted = (int) ( nb_steps * dt / 3600 / step_std) 
        if  nb_steps_adjusted * step_std * 3600 < nb_steps * dt:
            nb_steps_adjusted = nb_steps_adjusted + 1
        if ( nb_steps_adjusted - 1 ) * step_std_in_index >= nb_steps:
            nb_steps_adjusted = nb_steps_adjusted - 1
        std_daily = np.zeros([nb_steps_adjusted])
        med_daily = np.zeros([nb_steps_adjusted])
        mad_daily = np.zeros([nb_steps_adjusted])
        iqr_daily = np.zeros([nb_steps_adjusted])
        for i in range(nb_steps_adjusted):
            index_when_std = i * step_std_in_index
            std_daily[i] = np.std(algebric_distance_ensemble_main_sc_lvlh_x[index_when_std, :]) # NOT GOOD TO USE THE STD
            iqr_daily[i] = np.subtract(*np.percentile(algebric_distance_ensemble_main_sc_lvlh_x[index_when_std, :], [75, 25])) 
            mad_daily[i] = np.median( np.abs( algebric_distance_ensemble_main_sc_lvlh_x[index_when_std, :] - np.median(algebric_distance_ensemble_main_sc_lvlh_x[index_when_std, :]) ) ) 

        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.09, right=0.98, top = 0.93,bottom = 0.08)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Inter-quartile range of the along-track distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('IQR (km)', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Time (days)', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.plot(np.arange(0,nb_steps_adjusted), iqr_daily[:], 'k', linewidth = 2) 


        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

        if save_results == 1:
            fig_save_name = 'along_track_iqr_daily'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)
        

########################################################################################################################################################################################################################################################################################################
    if ( 'algebric_distance_ensemble_main_sc_lvlh_x' in pickle_list_name ):
        # Plot the distribution of angular_spacing_between_ensemble_sat_converted_in_a_distance at different times (set by when_plot_in_hour)
        when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
        index_when_plot = (int) (when_plot_in_hour * 3600L / dt)

        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('Spacecraft along-track distribution ' + 'after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(4, 1)
        gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

        med = np.median(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot, :])
        mad = np.median( np.abs( algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot, :] - np.median(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot, :]) ) )  
        quartile10 = np.percentile(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot, :], 10) 
        quartile25 = np.percentile(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot, :], 25) 
        quartile75 = np.percentile(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot, :], 75) 
        quartile90 = np.percentile(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot, :], 90) 
        iqr_daily = quartile75 - quartile25
        quartiles_1090_daily = quartile90 - quartile10

        ax2 = fig.add_subplot(gs[:3, 0])
        ax2.set_title('', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
        ax2.set_xlabel('', fontsize = fontsize_plot, weight = 'bold')
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
        bin_width = iqr_daily/10.
        hist_data = algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot, :]
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
        ax1.set_xlabel('Along-track distance from reference spacecraft (km)', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

        ax1.hlines(0, med - 10 * iqr_daily, med + 10 * iqr_daily, linewidth = 2, color = 'b');
        ax1.text( med + 10 * iqr_daily, 0 - length_vert_iqr_daily/0.8,'Along-track', horizontalalignment = 'right', fontsize = fontsize_plot, color = 'b' )

        ax1.plot(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot, :], np.zeros([nb_ensembles]), 'r.', markersize = 10)
        ax1.vlines(quartile25, 0,length_vert_iqr_daily, linewidth = 2); ax1.vlines(quartile75, 0,length_vert_iqr_daily, linewidth = 2) 
        ax1.text((quartile75 + quartile25 ) /2., length_vert_iqr_daily+length_vert_iqr_daily/5, '{0:.2f}'.format(iqr_daily ) + ' km (50% of sc)', horizontalalignment = 'center', fontsize = fontsize_plot )
        ax1.annotate ('', (quartile25, length_vert_iqr_daily), (quartile75, length_vert_iqr_daily), arrowprops={'arrowstyle':'<->'})
        ax1.vlines(quartile10, -length_vert_iqr_daily ,0 , linewidth = 4, linestyle = 'dotted'); ax1.vlines(quartile90, -length_vert_iqr_daily ,0, linewidth = 4, linestyle = 'dotted') 
        ax1.annotate ('', (quartile10, -length_vert_iqr_daily), (quartile90, -length_vert_iqr_daily), arrowprops={'arrowstyle':'<->'})
        ax1.text((quartile10 + quartile90 ) /2., -length_vert_iqr_daily-length_vert_iqr_daily/0.8, '{0:.2f}'.format(quartiles_1090_daily ) + ' km (80% of sc)', horizontalalignment = 'center', fontsize = fontsize_plot, verticalalignment = 'bottom' )
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
            fig_save_name = 'sc_along_track_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


############################## RADIAL SEPARATION WITH REFERENCE SPACECRAFT
########################################################################################################################################################################################################################################################################################################
if "radial" in sys.argv:
    algebric_distance_ensemble_main_sc_lvlh_z = algebric_distance_ensemble_main_sc_lvlh_z*1000. # !!!!!!!! conversion km to m (for cross track and radial, not for along track)
    if ( 'algebric_distance_ensemble_main_sc_lvlh_z' in pickle_list_name ):
        # Plot the standard deviation of the angular_spacing_between_ensemble_sat_converted_in_a_distance distributions as a function of time
        step_std = 10./60 # step in hours to calculate the standard deviation
        if ( step_std < dt / 3600. ):
            step_std = dt / 3600.
        step_std_in_index = step_std * 3600. / dt
        nb_steps_adjusted = (int) ( nb_steps * dt / 3600 / step_std) 
        if  nb_steps_adjusted * step_std * 3600 < nb_steps * dt:
            nb_steps_adjusted = nb_steps_adjusted + 1
        if ( nb_steps_adjusted - 1 ) * step_std_in_index >= nb_steps:
            nb_steps_adjusted = nb_steps_adjusted - 1
        std_daily = np.zeros([nb_steps_adjusted])
        med_daily = np.zeros([nb_steps_adjusted])
        mad_daily = np.zeros([nb_steps_adjusted])
        iqr_daily = np.zeros([nb_steps_adjusted])
        for i in range(nb_steps_adjusted):
            index_when_std = i * step_std_in_index
            std_daily[i] = np.std(algebric_distance_ensemble_main_sc_lvlh_z[index_when_std, :]) # NOT GOOD TO USE THE STD
            iqr_daily[i] = np.subtract(*np.percentile(algebric_distance_ensemble_main_sc_lvlh_z[index_when_std, :], [75, 25])) 
            mad_daily[i] = np.median( np.abs( algebric_distance_ensemble_main_sc_lvlh_z[index_when_std, :] - np.median(algebric_distance_ensemble_main_sc_lvlh_z[index_when_std, :]) ) ) 

        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.09, right=0.98, top = 0.93,bottom = 0.08)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Inter-quartile range of the radial distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('IQR (m)', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Time (days)', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.plot(np.arange(0,nb_steps_adjusted), iqr_daily[:], 'k', linewidth = 2) 


        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

        if save_results == 1:
            fig_save_name = 'radial_iqr_daily'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)
        

########################################################################################################################################################################################################################################################################################################
    if ( 'algebric_distance_ensemble_main_sc_lvlh_z' in pickle_list_name ):
        # Plot the distribution of angular_spacing_between_ensemble_sat_converted_in_a_distance at different times (set by when_plot_in_hour)
        when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
        index_when_plot = (int) (when_plot_in_hour * 3600L / dt)

        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('Spacecraft radial distribution ' + 'after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(4, 1)
        gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

        med = np.median(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot, :])
        mad = np.median( np.abs( algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot, :] - np.median(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot, :]) ) )  
        quartile10 = np.percentile(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot, :], 10) 
        quartile25 = np.percentile(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot, :], 25) 
        quartile75 = np.percentile(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot, :], 75) 
        quartile90 = np.percentile(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot, :], 90) 
        iqr_daily = quartile75 - quartile25
        quartiles_1090_daily = quartile90 - quartile10

        ax2 = fig.add_subplot(gs[:3, 0])
        ax2.set_title('', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
        ax2.set_xlabel('', fontsize = fontsize_plot, weight = 'bold')
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
        bin_width = iqr_daily/10.
        hist_data = algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot, :]
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
        ax1.set_xlabel('Radial distance from reference spacecraft (m)', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

        ax1.hlines(0, med - 10 * iqr_daily, med + 10 * iqr_daily, linewidth = 2, color = 'b');
        ax1.text( med + 10 * iqr_daily, 0 - length_vert_iqr_daily/0.8,'Radial', horizontalalignment = 'right', fontsize = fontsize_plot, color = 'b' )

        ax1.plot(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot, :], np.zeros([nb_ensembles]), 'r.', markersize = 10)
        ax1.vlines(quartile25, 0,length_vert_iqr_daily, linewidth = 2); ax1.vlines(quartile75, 0,length_vert_iqr_daily, linewidth = 2) 
        ax1.text((quartile75 + quartile25 ) /2., length_vert_iqr_daily+length_vert_iqr_daily/5, '{0:.2f}'.format(iqr_daily ) + ' m (50% of sc)', horizontalalignment = 'center', fontsize = fontsize_plot )
        ax1.annotate ('', (quartile25, length_vert_iqr_daily), (quartile75, length_vert_iqr_daily), arrowprops={'arrowstyle':'<->'})
        ax1.vlines(quartile10, -length_vert_iqr_daily ,0 , linewidth = 4, linestyle = 'dotted'); ax1.vlines(quartile90, -length_vert_iqr_daily ,0, linewidth = 4, linestyle = 'dotted') 
        ax1.annotate ('', (quartile10, -length_vert_iqr_daily), (quartile90, -length_vert_iqr_daily), arrowprops={'arrowstyle':'<->'})
        ax1.text((quartile10 + quartile90 ) /2., -length_vert_iqr_daily-length_vert_iqr_daily/0.8, '{0:.2f}'.format(quartiles_1090_daily ) + ' m (80% of sc)', horizontalalignment = 'center', fontsize = fontsize_plot, verticalalignment = 'bottom' )
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
            fig_save_name = 'sc_radial_track_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)





############################## CROSS-TRACK SEPARATION WITH REFERENCE SPACECRAFT
########################################################################################################################################################################################################################################################################################################
if "cross" in sys.argv:
    algebric_distance_ensemble_main_sc_lvlh_y = algebric_distance_ensemble_main_sc_lvlh_y*1000. # !!!!!!!! conversion km to m (for cross track and radial, not for along track)
    if ( 'algebric_distance_ensemble_main_sc_lvlh_y' in pickle_list_name ):
        # Plot the standard deviation of the angular_spacing_between_ensemble_sat_converted_in_a_distance distributions as a function of time
        step_std = 10./60 # step in hours to calculate the standard deviation
        if ( step_std < dt / 3600. ):
            step_std = dt / 3600.
        step_std_in_index = step_std * 3600. / dt
        nb_steps_adjusted = (int) ( nb_steps * dt / 3600 / step_std) 
        if  nb_steps_adjusted * step_std * 3600 < nb_steps * dt:
            nb_steps_adjusted = nb_steps_adjusted + 1
        if ( nb_steps_adjusted - 1 ) * step_std_in_index >= nb_steps:
            nb_steps_adjusted = nb_steps_adjusted - 1

        std_daily = np.zeros([nb_steps_adjusted])
        med_daily = np.zeros([nb_steps_adjusted])
        mad_daily = np.zeros([nb_steps_adjusted])
        iqr_daily = np.zeros([nb_steps_adjusted])
        for i in range(nb_steps_adjusted):
            index_when_std = i * step_std_in_index
            std_daily[i] = np.std(algebric_distance_ensemble_main_sc_lvlh_y[index_when_std, :]) # NOT GOOD TO USE THE STD
            iqr_daily[i] = np.subtract(*np.percentile(algebric_distance_ensemble_main_sc_lvlh_y[index_when_std, :], [75, 25])) 
            mad_daily[i] = np.median( np.abs( algebric_distance_ensemble_main_sc_lvlh_y[index_when_std, :] - np.median(algebric_distance_ensemble_main_sc_lvlh_y[index_when_std, :]) ) ) 

        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.09, right=0.98, top = 0.93,bottom = 0.08)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Inter-quartile range of the cross-track distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('IQR (m)', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Time (days)', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.plot(np.arange(0,nb_steps_adjusted), iqr_daily[:], 'k', linewidth = 2) 


        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

        if save_results == 1:
            fig_save_name = 'cross_track_iqr_daily'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)
        

########################################################################################################################################################################################################################################################################################################
    if ( 'algebric_distance_ensemble_main_sc_lvlh_y' in pickle_list_name ):
        # Plot the distribution of angular_spacing_between_ensemble_sat_converted_in_a_distance at different times (set by when_plot_in_hour)
        when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
        index_when_plot = (int) (when_plot_in_hour * 3600L / dt)

        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('Spacecraft cross-track distribution ' + 'after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(4, 1)
        gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

        med = np.median(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot, :])
        mad = np.median( np.abs( algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot, :] - np.median(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot, :]) ) )  
        quartile10 = np.percentile(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot, :], 10) 
        quartile25 = np.percentile(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot, :], 25) 
        quartile75 = np.percentile(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot, :], 75) 
        quartile90 = np.percentile(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot, :], 90) 
        iqr_daily = quartile75 - quartile25
        quartiles_1090_daily = quartile90 - quartile10

        ax2 = fig.add_subplot(gs[:3, 0])
        ax2.set_title('', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
        ax2.set_xlabel('', fontsize = fontsize_plot, weight = 'bold')
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
        bin_width = iqr_daily/10.
        hist_data = algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot, :]
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
        ax1.set_xlabel('Cross-track distance from reference spacecraft (m)', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

        ax1.hlines(0, med - 10 * iqr_daily, med + 10 * iqr_daily, linewidth = 2, color = 'b');
        ax1.text( med + 10 * iqr_daily, 0 - length_vert_iqr_daily/0.8,'Cross-track', horizontalalignment = 'right', fontsize = fontsize_plot, color = 'b' )

        ax1.plot(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot, :], np.zeros([nb_ensembles]), 'r.', markersize = 10)
        ax1.vlines(quartile25, 0,length_vert_iqr_daily, linewidth = 2); ax1.vlines(quartile75, 0,length_vert_iqr_daily, linewidth = 2) 
        ax1.text((quartile75 + quartile25 ) /2., length_vert_iqr_daily+length_vert_iqr_daily/5, '{0:.2f}'.format(iqr_daily ) + ' m (50% of sc)', horizontalalignment = 'center', fontsize = fontsize_plot )
        ax1.annotate ('', (quartile25, length_vert_iqr_daily), (quartile75, length_vert_iqr_daily), arrowprops={'arrowstyle':'<->'})
        ax1.vlines(quartile10, -length_vert_iqr_daily ,0 , linewidth = 4, linestyle = 'dotted'); ax1.vlines(quartile90, -length_vert_iqr_daily ,0, linewidth = 4, linestyle = 'dotted') 
        ax1.annotate ('', (quartile10, -length_vert_iqr_daily), (quartile90, -length_vert_iqr_daily), arrowprops={'arrowstyle':'<->'})
        ax1.text((quartile10 + quartile90 ) /2., -length_vert_iqr_daily-length_vert_iqr_daily/0.8, '{0:.2f}'.format(quartiles_1090_daily ) + ' m (80% of sc)', horizontalalignment = 'center', fontsize = fontsize_plot, verticalalignment = 'bottom' )
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
            fig_save_name = 'sc_cross_track_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)




############################## ALONG-TRACK (USING ANGULAR SEPARATION)
########################################################################################################################################################################################################################################################################################################
if "angle" in sys.argv:
    if ( 'angle_asc_node_to_sat_ensemble' in pickle_list_name ):
        # Plot the standard deviation of the angular_spacing_between_ensemble_sat_converted_in_a_distance distributions as a function of time
        step_std = 10./60 # step in hours to calculate the standard deviation
        if ( step_std < dt / 3600. ):
            step_std = dt / 3600.

        step_std_in_index = step_std * 3600. / dt
        nb_steps_adjusted = (int) ( nb_steps * dt / 3600 / step_std) 
        if  nb_steps_adjusted * step_std * 3600 < nb_steps * dt:
            nb_steps_adjusted = nb_steps_adjusted + 1
        if ( nb_steps_adjusted - 1 ) * step_std_in_index >= nb_steps:
            nb_steps_adjusted = nb_steps_adjusted - 1

        std_daily = np.zeros([nb_steps_adjusted])
        med_daily = np.zeros([nb_steps_adjusted])
        mad_daily = np.zeros([nb_steps_adjusted])
        iqr_daily = np.zeros([nb_steps_adjusted])
        for i in range(nb_steps_adjusted):
            index_when_std = i * step_std_in_index
            std_daily[i] = np.std(angle_asc_node_to_sat_ensemble[index_when_std, :]) # NOT GOOD TO USE THE STD
            iqr_daily[i] = np.subtract(*np.percentile(angle_asc_node_to_sat_ensemble[index_when_std, :], [75, 25])) 
            mad_daily[i] = np.median( np.abs( angle_asc_node_to_sat_ensemble[index_when_std, :] - np.median(angle_asc_node_to_sat_ensemble[index_when_std, :]) ) ) 

        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.09, right=0.98, top = 0.93,bottom = 0.08)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Inter-quartile range of the along-track distributions (using angular separation) VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('IQR (m)', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Time (days)', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.plot(np.arange(0,nb_steps_adjusted), iqr_daily[:]*110000, 'k', linewidth = 2)

        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

        if save_results == 1:
            fig_save_name = 'iqr_daily_along_track_using_angular_separation'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)
        

########################################################################################################################################################################################################################################################################################################
    if ( 'angle_asc_node_to_sat_ensemble' in pickle_list_name ):
        # Plot the distribution of angular_spacing_between_ensemble_sat_converted_in_a_distance at different times (set by when_plot_in_hour)
        when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
        index_when_plot = (int) (when_plot_in_hour * 3600L / dt)

        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('Spacecraft distribution along the orbit (using angular separation) ' + 'after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(4, 1)
        gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

        med = np.median(angle_asc_node_to_sat_ensemble[index_when_plot, :])
        mad = np.median( np.abs( angle_asc_node_to_sat_ensemble[index_when_plot, :] - np.median(angle_asc_node_to_sat_ensemble[index_when_plot, :]) ) )  
        quartile10 = np.percentile(angle_asc_node_to_sat_ensemble[index_when_plot, :], 10) 
        quartile25 = np.percentile(angle_asc_node_to_sat_ensemble[index_when_plot, :], 25) 
        quartile75 = np.percentile(angle_asc_node_to_sat_ensemble[index_when_plot, :], 75) 
        quartile90 = np.percentile(angle_asc_node_to_sat_ensemble[index_when_plot, :], 90) 
        iqr_daily = quartile75 - quartile25
        quartiles_1090_daily = quartile90 - quartile10

        ax2 = fig.add_subplot(gs[:3, 0])
        ax2.set_title('', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
        ax2.set_xlabel('', fontsize = fontsize_plot, weight = 'bold')
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
        bin_width = iqr_daily/10.
        hist_data = angle_asc_node_to_sat_ensemble[index_when_plot, :]
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

        ax1.plot(angle_asc_node_to_sat_ensemble[index_when_plot, :], np.zeros([nb_ensembles]), 'r.', markersize = 10)
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
            fig_save_name = 'sc_along_track_distribution_using_angular_separation' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)




    ############################## ORBIT AVERAGE DENSITY
    ########################################################################################################################################################################################################################################################################################################
if "rho" in sys.argv:
    if ( 'rho_ensemble' in pickle_list_name ):
        nb_steps_orbit_average = len(index_time_orbit_average[:,1])     # !!!!!!!! there might be ore periods for satellites than for others (because they move at different speeds). So nb_steps_orbit_average should actually be the min of the number of periods of all satellites. Since we don't take the min, there might be 0 in rho_ensemble_orbit_average that we should NOT take into account. So basically NEVER TAKE THE STANDARD DEVIATION for this density section of the script (that takes into account outliers)
        index_time_orbit_average_middle = index_time_orbit_average[:,1]
        x_axis_orbit_average = index_time_orbit_average_middle * dt # conversion index to seconds
        x_axis_orbit_average = x_axis_orbit_average / (24*3600.) # conversion seconds to day

        # Plot the inter-quartile rang of the density distributions as a function of time
        step_std_in_index = 1 # we want the step to be every orbit (here we are look at orbit average data)
        std_every_orbit = np.zeros([nb_steps_orbit_average]) # DON'T TAKE IT FOR THE DENSITY SECTION
        med_every_orbit = np.zeros([nb_steps_orbit_average])
        mad_every_orbit = np.zeros([nb_steps_orbit_average])
        iqr_every_orbit = np.zeros([nb_steps_orbit_average])
        for i in range(nb_steps_orbit_average):
            index_when_std = i * step_std_in_index
            std_every_orbit[i] = np.std(rho_ensemble_orbit_average[index_when_std, :]) # NOT GOOD TO USE THE STD
            iqr_every_orbit[i] = np.subtract(*np.percentile(rho_ensemble_orbit_average[index_when_std, :], [75, 25])) 
            mad_every_orbit[i] = np.median( np.abs( rho_ensemble_orbit_average[index_when_std, :] - np.median(rho_ensemble_orbit_average[index_when_std, :]) ) ) 


        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.09, right=0.98, top = 0.93,bottom = 0.08)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Inter-quartile range of the density distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('IQR (#/m^3)', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Time (days)', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.plot(x_axis_orbit_average, iqr_every_orbit, 'k', linewidth = 2)

        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

        if save_results == 1:
            fig_save_name = 'rho_iqr_every_orbit'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

    ########################################################################################################################################################################################################################################################################################################
    if ( 'rho_ensemble_orbit_average' in pickle_list_name ):
        # Plot the distribution of the density at different times (set by when_plot_in_hour)
        when_plot_in_hour = 2 * 24. # the few lines below convert this time in orbits 
        time_between_index_and_when_plot_in_hour = np.zeros([len(index_time_orbit_average_middle)])
        icount = -1
        for i in index_time_orbit_average_middle:
            icount = icount + 1
            time_between_index_and_when_plot_in_hour[icount] = np.abs( when_plot_in_hour * 3600. - i * dt )

        index_when_plot = (int)( np.where(time_between_index_and_when_plot_in_hour == min(time_between_index_and_when_plot_in_hour))[0][0] )
        nb_hour_this_orbit = index_time_orbit_average_middle[index_when_plot] * dt/3600./24

        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

        med = np.median(rho_ensemble_orbit_average[index_when_plot, :])
        mad = np.median( np.abs( rho_ensemble_orbit_average[index_when_plot, :] - np.median(rho_ensemble_orbit_average[index_when_plot, :]) ) )  
        quartile10 = np.percentile(rho_ensemble_orbit_average[index_when_plot, :], 10) 
        quartile25 = np.percentile(rho_ensemble_orbit_average[index_when_plot, :], 25) 
        quartile75 = np.percentile(rho_ensemble_orbit_average[index_when_plot, :], 75) 
        quartile90 = np.percentile(rho_ensemble_orbit_average[index_when_plot, :], 90) 
        iqr_every_orbit = quartile75 - quartile25
        quartiles_1090_every_orbit = quartile90 - quartile10
        std = np.std(rho_ensemble_orbit_average[index_when_plot, :])

        ax2 = fig.add_subplot(gs[0, 0])
        ax2.set_title('Density distribution ' + 'after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax2.set_xlabel('Density (#/m^3)', fontsize = fontsize_plot, weight = 'bold')
        ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
        bin_width = iqr_every_orbit/10.
        hist_data = rho_ensemble_orbit_average[index_when_plot, :]
        n, bins, patches = ax2.hist(hist_data , bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'r',label = '')  
        plt.vlines(quartile25, 0, max(n),linewidth = 2); plt.vlines(quartile75, 0,max(n), linewidth = 2, label = '25-75% quartiles') 
        plt.vlines(quartile10, 0, max(n), linewidth = 4, linestyle = 'dotted'); plt.vlines(quartile90, 0, max(n), linewidth = 4, linestyle = 'dotted', label = '10-90% deciles') 
        ax2.plot(med, 0, 'k',marker = '.', markersize = 20)
        ax2.get_xaxis().get_major_formatter().set_useOffset(False)
        ax2.set_xlim([med - 10*iqr_every_orbit, med+10*iqr_every_orbit])
        ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
        ax2.text(ax2.get_xlim()[0] + ( ax2.get_xlim()[1] - ax2.get_xlim()[0] )/ 60., ax2.get_ylim()[1] - (ax2.get_ylim()[1] - ax2.get_ylim()[0])/5., 'IQR = ' + '{0:.2e} #/m^3'.format(iqr_every_orbit), fontsize = fontsize_plot)

        ax2.legend(fontsize = fontsize_plot, loc = 2, scatterpoints=1)
        if save_results == 1:
            fig_save_name = 'rho_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


    if ( 'rho_ensemble_orbit_average' in pickle_list_name ):
        # Plot the density of each ensemble as a function of time, as well as the median and the inter-quartile range
        step_std_in_index = 1 # every orbit
        med_every_orbit = np.zeros([nb_steps_orbit_average])
        mad_every_orbit = np.zeros([nb_steps_orbit_average])
        iqr_every_orbit = np.zeros([nb_steps_orbit_average])
        quartile25 =np.zeros([nb_steps_orbit_average])
        quartile75 =np.zeros([nb_steps_orbit_average])
        for i in range(nb_steps_orbit_average):
            index_when_std = i * step_std_in_index
            iqr_every_orbit[i] = np.subtract(*np.percentile(rho_ensemble_orbit_average[index_when_std, :], [75, 25])) 
            mad_every_orbit[i] = np.median( np.abs( rho_ensemble_orbit_average[index_when_std, :] - np.median(rho_ensemble_orbit_average[index_when_std, :]) ) ) 
            med_every_orbit[i] = np.median(rho_ensemble_orbit_average[index_when_std, :])
            quartile25[i] = np.percentile(rho_ensemble_orbit_average[index_when_std, :], 25) 
            quartile75[i] = np.percentile(rho_ensemble_orbit_average[index_when_std, :], 75) 

        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.09, right=0.98, top = 0.93,bottom = 0.08)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Density at the position of each ensemble member VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('Density (#/m^3)', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Time (days)', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        for iens in range(nb_ensembles):
            if iens == 0:
                ax1.plot(x_axis_orbit_average, rho_ensemble_orbit_average[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 1, alpha = 0.1, label = 'Ensemble')
            else:
                ax1.plot(x_axis_orbit_average, rho_ensemble_orbit_average[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 1, alpha = 0.1)
        ax1.plot(x_axis_orbit_average, med_every_orbit, 'k', linewidth = 3, label = 'Median')
        ax1.plot(x_axis_orbit_average, quartile25, 'r', linewidth = 3, label = '25 and 75% quartiles')
        ax1.plot(x_axis_orbit_average, quartile75, 'r', linewidth = 3)

        ax1.legend(fontsize = fontsize_plot, loc = 2)
        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

        if save_results == 1:
            fig_save_name = 'rho_all_ensemble'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)



    ############################## F10.7
    ########################################################################################################################################################################################################################################################################################################
if "f107" in sys.argv:
    if ( 'f107_ensemble' in pickle_list_name ):
        # Plot the standard deviation of the f10.7 distributions as a function of time
        step_std = 10./60 # step in hours to calculate the standard deviation
        if ( step_std < dt / 3600. ):
            step_std = dt / 3600.
        step_std_in_index = step_std * 3600. / dt
        nb_steps_adjusted = (int) ( nb_steps * dt / 3600 / step_std) 
        if  nb_steps_adjusted * step_std * 3600 < nb_steps * dt:
            nb_steps_adjusted = nb_steps_adjusted + 1
        if ( nb_steps_adjusted - 1 ) * step_std_in_index >= nb_steps:
            nb_steps_adjusted = nb_steps_adjusted - 1

        std_daily = np.zeros([nb_steps_adjusted])
        med_daily = np.zeros([nb_steps_adjusted])
        mad_daily = np.zeros([nb_steps_adjusted])
        iqr_daily = np.zeros([nb_steps_adjusted])
        for i in range(nb_steps_adjusted):
            index_when_std = i * step_std_in_index
            std_daily[i] = np.std(f107_ensemble[index_when_std, :]) # NOT GOOD TO USE THE STD
            iqr_daily[i] = np.subtract(*np.percentile(f107_ensemble[index_when_std, :], [75, 25])) 
            mad_daily[i] = np.median( np.abs( f107_ensemble[index_when_std, :] - np.median(f107_ensemble[index_when_std, :]) ) ) 

        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.09, right=0.98, top = 0.93,bottom = 0.08)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Standard deviation of the f10.7 distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('Standard deviation', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Time (days)', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.plot(np.arange(0,nb_steps_adjusted), std_daily[:], 'k', linewidth = 2)

        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

        if save_results == 1:
            fig_save_name = 'f107_std_daily'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

    ########################################################################################################################################################################################################################################################################################################
    if ( 'f107_ensemble' in pickle_list_name ):
        # Plot the distribution of the f10.7 at different times (set by when_plot_in_hour)
        when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
        index_when_plot = (int) (when_plot_in_hour * 3600L / dt)

        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

        med = np.median(f107_ensemble[index_when_plot, :])
        mad = np.median( np.abs( f107_ensemble[index_when_plot, :] - np.median(f107_ensemble[index_when_plot, :]) ) )  
        quartile10 = np.percentile(f107_ensemble[index_when_plot, :], 10) 
        quartile25 = np.percentile(f107_ensemble[index_when_plot, :], 25) 
        quartile75 = np.percentile(f107_ensemble[index_when_plot, :], 75) 
        quartile90 = np.percentile(f107_ensemble[index_when_plot, :], 90) 
        std = np.std(f107_ensemble[index_when_plot, :])
        iqr_daily = quartile75 - quartile25
        quartiles_1090_daily = quartile90 - quartile10

        ax2 = fig.add_subplot(gs[0, 0])
        ax2.set_title('F10.7 distribution ' + 'after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
        bin_width = iqr_daily/10.
        hist_data = f107_ensemble[index_when_plot, :]
        n, bins, patches = ax2.hist(hist_data , bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'r',label = '')  
        plt.vlines(quartile25, 0, max(n),linewidth = 2); plt.vlines(quartile75, 0,max(n), linewidth = 2, label = '25-75% quartiles') 
        plt.vlines(quartile10, 0, max(n), linewidth = 4, linestyle = 'dotted'); plt.vlines(quartile90, 0, max(n), linewidth = 4, linestyle = 'dotted', label = '10-90% deciles') 
        ax2.plot(med, 0, 'k',marker = '.', markersize = 20)
        ax2.get_xaxis().get_major_formatter().set_useOffset(False)
        ax2.set_xlim([med - 10*iqr_daily, med+10*iqr_daily])
        ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
        ax2.text(ax2.get_xlim()[0] + ( ax2.get_xlim()[1] - ax2.get_xlim()[0] )/ 60., ax2.get_ylim()[1] - (ax2.get_ylim()[1] - ax2.get_ylim()[0])/5., 'std = ' + '{0:.2f}'.format(std), fontsize = fontsize_plot)
        ax2.set_xlabel('F10.7', fontsize = fontsize_plot, weight = 'bold')

        ax2.legend(fontsize = fontsize_plot, loc = 2, scatterpoints=1)
        if save_results == 1:
            fig_save_name = 'f107_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


    if ( 'f107_ensemble' in pickle_list_name ):
        # Plot the f10.7 for each ensemble as a function of time, as well as the median and the inter-quartile range
        step_std = 2/60. # step in hours to calculate the standard deviation
        if ( step_std < dt / 3600. ):
            step_std = dt / 3600.
        step_std_in_index = step_std * 3600. / dt
        nb_steps_adjusted = (int) ( nb_steps * dt / 3600 / step_std) 
        if  nb_steps_adjusted * step_std * 3600 < nb_steps * dt:
            nb_steps_adjusted = nb_steps_adjusted + 1
        if ( nb_steps_adjusted - 1 ) * step_std_in_index >= nb_steps:
            nb_steps_adjusted = nb_steps_adjusted - 1

        med_daily = np.zeros([nb_steps_adjusted])
        mad_daily = np.zeros([nb_steps_adjusted])
        iqr_daily = np.zeros([nb_steps_adjusted])
        quartile25 =np.zeros([nb_steps_adjusted])
        quartile75 =np.zeros([nb_steps_adjusted])
        for i in range(nb_steps_adjusted):
            index_when_std = i * step_std_in_index
            iqr_daily[i] = np.subtract(*np.percentile(f107_ensemble[index_when_std, :], [75, 25])) 
            mad_daily[i] = np.median( np.abs( f107_ensemble[index_when_std, :] - np.median(f107_ensemble[index_when_std, :]) ) ) 
            med_daily[i] = np.median(f107_ensemble[index_when_std, :])
            quartile25[i] = np.percentile(f107_ensemble[index_when_std, :], 25) 
            quartile75[i] = np.percentile(f107_ensemble[index_when_std, :], 75) 

        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.09, right=0.98, top = 0.93,bottom = 0.08)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('F10.7 for each ensemble member VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('F10.7', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Time (days)', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        for iens in range(nb_ensembles):
            if iens == 0:
                ax1.plot(np.arange(0,nb_steps_adjusted), f107_ensemble[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 1, alpha = 0.1, label = 'Ensemble')
            else:
                ax1.plot(np.arange(0,nb_steps_adjusted), f107_ensemble[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 1, alpha = 0.1)
        ax1.plot(np.arange(0,nb_steps_adjusted), med_daily[:], 'k', linewidth = 3, label = 'Median')
        ax1.plot(np.arange(0,nb_steps_adjusted), quartile25[:], 'r', linewidth = 3, label = '25 and 75% quartiles')
        ax1.plot(np.arange(0,nb_steps_adjusted), quartile75[:], 'r', linewidth = 3)

        ax1.legend(fontsize = fontsize_plot, loc = 2)
        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

        if save_results == 1:
            fig_save_name = 'f107_all_ensemble'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)



    ############################## Ap
    ########################################################################################################################################################################################################################################################################################################
if "ap" in sys.argv:
    if ( 'ap_ensemble' in pickle_list_name ):
        # Plot the standard deviation of the Ap distributions as a function of time
        step_std = 10./60 # step in hours to calculate the standard deviation
        if ( step_std < dt / 3600. ):
            step_std = dt / 3600.

        step_std_in_index = step_std * 3600. / dt
        nb_steps_adjusted = (int) ( nb_steps * dt / 3600 / step_std) 
        if  nb_steps_adjusted * step_std * 3600 < nb_steps * dt:
            nb_steps_adjusted = nb_steps_adjusted + 1
        if ( nb_steps_adjusted - 1 ) * step_std_in_index >= nb_steps:
            nb_steps_adjusted = nb_steps_adjusted - 1

        std_daily = np.zeros([nb_steps_adjusted])
        med_daily = np.zeros([nb_steps_adjusted])
        mad_daily = np.zeros([nb_steps_adjusted])
        iqr_daily = np.zeros([nb_steps_adjusted])
        for i in range(nb_steps_adjusted):
            index_when_std = i * step_std_in_index
            std_daily[i] = np.std(ap_ensemble[index_when_std, :]) # NOT GOOD TO USE THE STD
            iqr_daily[i] = np.subtract(*np.percentile(ap_ensemble[index_when_std, :], [75, 25])) 
            mad_daily[i] = np.median( np.abs( ap_ensemble[index_when_std, :] - np.median(ap_ensemble[index_when_std, :]) ) ) 

        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.09, right=0.98, top = 0.93,bottom = 0.08)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Standard deviation of the Ap distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('Standard deviation', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Time (days)', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.plot(np.arange(0,nb_steps_adjusted), std_daily[:], 'k', linewidth = 2)

        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

        if save_results == 1:
            fig_save_name = 'ap_std_daily'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

    ########################################################################################################################################################################################################################################################################################################
    if ( 'ap_ensemble' in pickle_list_name ):
        # Plot the distribution of the Ap at different times (set by when_plot_in_hour)
        when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
        index_when_plot = (int) (when_plot_in_hour * 3600L / dt)

        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

        med = np.median(ap_ensemble[index_when_plot, :])
        mad = np.median( np.abs( ap_ensemble[index_when_plot, :] - np.median(ap_ensemble[index_when_plot, :]) ) )  
        quartile10 = np.percentile(ap_ensemble[index_when_plot, :], 10) 
        quartile25 = np.percentile(ap_ensemble[index_when_plot, :], 25) 
        quartile75 = np.percentile(ap_ensemble[index_when_plot, :], 75) 
        quartile90 = np.percentile(ap_ensemble[index_when_plot, :], 90) 
        std = np.std(ap_ensemble[index_when_plot, :])
        iqr_daily = quartile75 - quartile25
        quartiles_1090_daily = quartile90 - quartile10

        ax2 = fig.add_subplot(gs[0, 0])
        ax2.set_title('Ap distribution ' + 'after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
        bin_width = iqr_daily/10.
        hist_data = ap_ensemble[index_when_plot, :]
        n, bins, patches = ax2.hist(hist_data , bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'r',label = '')  
        plt.vlines(quartile25, 0, max(n),linewidth = 2); plt.vlines(quartile75, 0,max(n), linewidth = 2, label = '25-75% quartiles') 
        plt.vlines(quartile10, 0, max(n), linewidth = 4, linestyle = 'dotted'); plt.vlines(quartile90, 0, max(n), linewidth = 4, linestyle = 'dotted', label = '10-90% deciles') 
        ax2.plot(med, 0, 'k',marker = '.', markersize = 20)
        ax2.get_xaxis().get_major_formatter().set_useOffset(False)
        ax2.set_xlim([med - 10*iqr_daily, med+10*iqr_daily])
        ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
        ax2.text(ax2.get_xlim()[0] + ( ax2.get_xlim()[1] - ax2.get_xlim()[0] )/ 60., ax2.get_ylim()[1] - (ax2.get_ylim()[1] - ax2.get_ylim()[0])/5., 'std = ' + '{0:.2f}'.format(std), fontsize = fontsize_plot)
        ax2.set_xlabel('Ap', fontsize = fontsize_plot, weight = 'bold')

        ax2.legend(fontsize = fontsize_plot, loc = 2, scatterpoints=1)
        if save_results == 1:
            fig_save_name = 'ap_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


    if ( 'ap_ensemble' in pickle_list_name ):
        # Plot the Ap for each ensemble as a function of time, as well as the median and the inter-quartile range
        step_std = 10./60 # step in hours to calculate the standard deviation
        if ( step_std < dt / 3600. ):
            step_std = dt / 3600.

        step_std_in_index = step_std * 3600. / dt
        nb_steps_adjusted = (int) ( nb_steps * dt / 3600 / step_std) 
        if  nb_steps_adjusted * step_std * 3600 < nb_steps * dt:
            nb_steps_adjusted = nb_steps_adjusted + 1
        if ( nb_steps_adjusted - 1 ) * step_std_in_index >= nb_steps:
            nb_steps_adjusted = nb_steps_adjusted - 1

        med_daily = np.zeros([nb_steps_adjusted])
        mad_daily = np.zeros([nb_steps_adjusted])
        iqr_daily = np.zeros([nb_steps_adjusted])
        quartile25 =np.zeros([nb_steps_adjusted])
        quartile75 =np.zeros([nb_steps_adjusted])
        for i in range(nb_steps_adjusted):
            index_when_std = i * step_std_in_index
            iqr_daily[i] = np.subtract(*np.percentile(ap_ensemble[index_when_std, :], [75, 25])) 
            mad_daily[i] = np.median( np.abs( ap_ensemble[index_when_std, :] - np.median(ap_ensemble[index_when_std, :]) ) ) 
            med_daily[i] = np.median(ap_ensemble[index_when_std, :])
            quartile25[i] = np.percentile(ap_ensemble[index_when_std, :], 25) 
            quartile75[i] = np.percentile(ap_ensemble[index_when_std, :], 75) 

        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.09, right=0.98, top = 0.93,bottom = 0.08)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Ap for each ensemble member VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('Ap', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Time (days)', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        for iens in range(nb_ensembles):
            if iens == 0:
                ax1.plot(np.arange(0,nb_steps_adjusted), ap_ensemble[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 1, alpha = 0.1, label = 'Ensemble')
            else:
                ax1.plot(np.arange(0,nb_steps_adjusted), ap_ensemble[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 1, alpha = 0.1)
        ax1.plot(np.arange(0,nb_steps_adjusted), med_daily[:], 'k', linewidth = 3, label = 'Median')
        ax1.plot(np.arange(0,nb_steps_adjusted), quartile25[:], 'r', linewidth = 3, label = '25 and 75% quartiles')
        ax1.plot(np.arange(0,nb_steps_adjusted), quartile75[:], 'r', linewidth = 3)

        ax1.legend(fontsize = fontsize_plot, loc = 2)
        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left() 
        ax1.margins(0,1) # autoscale both axes(fist value is for the x axis, second value for the y axis)

        if save_results == 1:
            fig_save_name = 'ap_all_ensemble'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


    ############################## F10.7A
    ########################################################################################################################################################################################################################################################################################################
if "f107a" in sys.argv:
    if ( 'f107a_ensemble' in pickle_list_name ):
        # Plot the standard deviation of the F10.7A distributions as a function of time
        step_std = 10./60 # step in hours to calculate the standard deviation
        if ( step_std < dt / 3600. ):
            step_std = dt / 3600.
        step_std_in_index = step_std * 3600. / dt
        nb_steps_adjusted = (int) ( nb_steps * dt / 3600 / step_std) 
        if  nb_steps_adjusted * step_std * 3600 < nb_steps * dt:
            nb_steps_adjusted = nb_steps_adjusted + 1
        if ( nb_steps_adjusted - 1 ) * step_std_in_index >= nb_steps:
            nb_steps_adjusted = nb_steps_adjusted - 1

        std_daily = np.zeros([nb_steps_adjusted])
        med_daily = np.zeros([nb_steps_adjusted])
        mad_daily = np.zeros([nb_steps_adjusted])
        iqr_daily = np.zeros([nb_steps_adjusted])
        for i in range(nb_steps_adjusted):
            index_when_std = i * step_std_in_index
            std_daily[i] = np.std(f107a_ensemble[index_when_std, :]) # NOT GOOD TO USE THE STD
            iqr_daily[i] = np.subtract(*np.percentile(f107a_ensemble[index_when_std, :], [75, 25])) 
            mad_daily[i] = np.median( np.abs( f107a_ensemble[index_when_std, :] - np.median(f107a_ensemble[index_when_std, :]) ) ) 

        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.09, right=0.98, top = 0.93,bottom = 0.08)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Standard deviation of the F10.7A distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('Standard deviation', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Time (days)', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.plot(np.arange(0,nb_steps_adjusted), std_daily[:], 'k', linewidth = 2)

        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

        if save_results == 1:
            fig_save_name = 'f107a_std_daily'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

    ########################################################################################################################################################################################################################################################################################################
    if ( 'f107a_ensemble' in pickle_list_name ):
        # Plot the distribution of the F10.7A at different times (set by when_plot_in_hour)
        when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
        index_when_plot = (int) (when_plot_in_hour * 3600L / dt)

        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

        med = np.median(f107a_ensemble[index_when_plot, :])
        mad = np.median( np.abs( f107a_ensemble[index_when_plot, :] - np.median(f107a_ensemble[index_when_plot, :]) ) )  
        quartile10 = np.percentile(f107a_ensemble[index_when_plot, :], 10) 
        quartile25 = np.percentile(f107a_ensemble[index_when_plot, :], 25) 
        quartile75 = np.percentile(f107a_ensemble[index_when_plot, :], 75) 
        quartile90 = np.percentile(f107a_ensemble[index_when_plot, :], 90) 
        std = np.std(f107a_ensemble[index_when_plot, :])
        iqr_daily = quartile75 - quartile25
        quartiles_1090_daily = quartile90 - quartile10

        ax2 = fig.add_subplot(gs[0, 0])
        ax2.set_title('F10.7A distribution ' + 'after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
        bin_width = iqr_daily/10.
        hist_data = f107a_ensemble[index_when_plot, :]
        n, bins, patches = ax2.hist(hist_data , bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'r',label = '')  
        plt.vlines(quartile25, 0, max(n),linewidth = 2); plt.vlines(quartile75, 0,max(n), linewidth = 2, label = '25-75% quartiles') 
        plt.vlines(quartile10, 0, max(n), linewidth = 4, linestyle = 'dotted'); plt.vlines(quartile90, 0, max(n), linewidth = 4, linestyle = 'dotted', label = '10-90% deciles') 
        ax2.plot(med, 0, 'k',marker = '.', markersize = 20)
        ax2.get_xaxis().get_major_formatter().set_useOffset(False)
        ax2.set_xlim([med - 10*iqr_daily, med+10*iqr_daily])
        ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
        ax2.text(ax2.get_xlim()[0] + ( ax2.get_xlim()[1] - ax2.get_xlim()[0] )/ 60., ax2.get_ylim()[1] - (ax2.get_ylim()[1] - ax2.get_ylim()[0])/5., 'std = ' + '{0:.2f}'.format(std), fontsize = fontsize_plot)
        ax2.set_xlabel('F10.7A', fontsize = fontsize_plot, weight = 'bold')

        ax2.legend(fontsize = fontsize_plot, loc = 2, scatterpoints=1)
        if save_results == 1:
            fig_save_name = 'f107a_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


    if ( 'f107a_ensemble' in pickle_list_name ):
        # Plot the F10.7A for each ensemble as a function of time, as well as the median and the inter-quartile range
        step_std = 10./60 # step in hours to calculate the standard deviation
        if ( step_std < dt / 3600. ):
            step_std = dt / 3600.
        step_std_in_index = step_std * 3600. / dt
        nb_steps_adjusted = (int) ( nb_steps * dt / 3600 / step_std) 
        if  nb_steps_adjusted * step_std * 3600 < nb_steps * dt:
            nb_steps_adjusted = nb_steps_adjusted + 1
        if ( nb_steps_adjusted - 1 ) * step_std_in_index >= nb_steps:
            nb_steps_adjusted = nb_steps_adjusted - 1

        med_daily = np.zeros([nb_steps_adjusted])
        mad_daily = np.zeros([nb_steps_adjusted])
        iqr_daily = np.zeros([nb_steps_adjusted])
        quartile25 =np.zeros([nb_steps_adjusted])
        quartile75 =np.zeros([nb_steps_adjusted])
        for i in range(nb_steps_adjusted):
            index_when_std = i * step_std_in_index
            iqr_daily[i] = np.subtract(*np.percentile(f107a_ensemble[index_when_std, :], [75, 25])) 
            mad_daily[i] = np.median( np.abs( f107a_ensemble[index_when_std, :] - np.median(f107a_ensemble[index_when_std, :]) ) ) 
            med_daily[i] = np.median(f107a_ensemble[index_when_std, :])
            quartile25[i] = np.percentile(f107a_ensemble[index_when_std, :], 25) 
            quartile75[i] = np.percentile(f107a_ensemble[index_when_std, :], 75) 

        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.09, right=0.98, top = 0.93,bottom = 0.08)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('F10.7A for each ensemble member VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('F10.7A', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Time (days)', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        for iens in range(nb_ensembles):
            if iens == 0:
                ax1.plot(np.arange(0,nb_steps_adjusted), f107a_ensemble[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 1, alpha = 0.1, label = 'Ensemble')
            else:
                ax1.plot(np.arange(0,nb_steps_adjusted), f107a_ensemble[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 1, alpha = 0.1)
        ax1.plot(np.arange(0,nb_steps_adjusted), med_daily[:], 'k', linewidth = 3, label = 'Median')
        ax1.plot(np.arange(0,nb_steps_adjusted), quartile25[:], 'r', linewidth = 3, label = '25 and 75% quartiles')
        ax1.plot(np.arange(0,nb_steps_adjusted), quartile75[:], 'r', linewidth = 3)

        ax1.legend(fontsize = fontsize_plot, loc = 2)
        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

        if save_results == 1:
            fig_save_name = 'f107a_all_ensemble'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


if show_plots == 1:
    plt.show(); plt.show()
