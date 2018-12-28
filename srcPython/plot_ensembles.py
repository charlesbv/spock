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
from find_in_read_input_order_variables import *
from read_collision_file import *
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
## NOTE 2: nothing
## NOTE 3: the results (pickle, image, video) are stored in a subfolder of path_folder_results. The subfolder is either 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', or 'other'. This code figures out which folder the simulation corresponds to by reading the name of the output file chosen by the user in the main input file of the propagator (third line of section #SPACECRAFTS): it tries to find 'cygnss', 'cadre', 'aerie', 'scion', or 'qb50' in the name of the output file. If it does not find it, then the results here will be stored in path_folder_results/other/
## NOTE 4: to run this script, and any python script in the propagator, you need to be one subfolder deep from the main folder where the propagator runs are made. So if path_to_propagator/PropSim is the folder where the propagator runs are made, then the python scripts must be run for example in path_to_propagator/PropSim/subfolder_where_python_scipts_are_run


# !!!!!!!!!! SET THE PARAMETER BELOW:
## Save or not the plots
save_results = 1

## Show or not the plots
show_plots = 0

## Which reference satellite to look at (starts at 0 and ends and nb_spacecraft - 1)
isat = 0

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
if isat >= nb_spacecraft:
    print "You choose satellite " + str(isat + 1) + " but there is only " + str(nb_spacecraft) + " satellite(s). Therefore, the satellite is set to satellite 1."
    isat = 0

compute_drag = input_variables[19]

## path of the folder where you want to store the results (pickle, image, video)
path_folder_results = '/'.join(satellite_to_plot_path[0].split('/')[:-2]) + '/python_out/'  #get_prop_dir(2) + 'output/python_propagator/'
if (os.path.isdir(path_folder_results) == False):
    os.system("mkdir " + path_folder_results)



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
    if (ensemble_to_plot_temp[i] == 'cd'):
        ensemble_to_plot.append('cd')
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
    if ((ensemble_to_plot_temp[i] == 'cd') & (compute_drag == 1)):
        ensemble_to_plot.append('cd')

        
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
    
    if os.path.isdir(path_folder_results + name_mission + '/result') == False:
        os.system("mkdir " + path_folder_results + name_mission + '/result' )
    if os.path.isdir(path_folder_results + name_mission + '/result/image/' ) == False:
        os.system("mkdir " + path_folder_results + name_mission + '/result/image/' )
    if os.path.isdir(path_folder_results + name_mission + '/result/image/' + name_subfolder_save) == False:
        os.system("mkdir " + path_folder_results + name_mission + '/result/image/' + name_subfolder_save )
        #    os.system('ssh -t srbwks2014-0008.engin.umich.edu "mkdir ' + name_mission + '/' + name_subfolder_save + '"')
        #    os.system("mkdir " + path_folder_results + name_mission + '/result/video/' + name_subfolder_save )

    
save_pickle_name = path_folder_results + name_mission + '/pickle/' + name_subfolder_save + satellite_to_plot[isat].replace(".txt",".pickle")
root_save_fig_name = path_folder_results + name_mission + '/result/image/' + name_subfolder_save + satellite_to_plot[isat].replace(".txt","_") 
root_save_video_name = path_folder_results + name_mission + '/result/video/' + name_subfolder_save + satellite_to_plot[isat].replace(".txt","_")



# # #######################################
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
    if nProcs > np.sqrt( nb_ensembles ):
        print "Note: to plot the density, F10.7, F10.7A, and Ap of ensembles as a function of time, it is recommended to run SpOCK with nProcs <= sqrt( nb_ensembles ) (it will still run fine and plot fine)."
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
    if ( 'cd_ensemble' in pickle_list_name_temp ):
        cd_ensemble = pickle_list[ pickle_list_name_temp.index('cd_ensemble') ]

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

    if ( 'tca_ensemble' in pickle_list_name_temp ):
        tca_ensemble = pickle_list[ pickle_list_name_temp.index('tca_ensemble') ]
    if ( 'dca_ensemble' in pickle_list_name_temp ):
        dca_ensemble = pickle_list[ pickle_list_name_temp.index('dca_ensemble') ]


    if ipickle == 0:
        if (pickle_list_name_temp[1] != "cd_ensemble"):
            nb_steps = pickle_list[1].shape[0] #!!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
        elif len(pickle_list_name_temp) > 1:
            nb_steps = pickle_list[2].shape[0] #!!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
        
########################################################################################################################################################################################################################################################################################################
# PLOTS ################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################

# Parameters of figures
height_fig = 9.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
when_plot_in_hour =0*24# 3*24 #5 + 50/60. # can be overwritten for each plot by uncommenting the same line in each plot section
index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
step_std = 1./60 # step in hours to calculate the standard deviation
hour_time_step_xticks = 24. # time step of ticks when plotting a function as a function of time


###

if index_when_plot > nb_steps:
    when_plot_in_hour = ( nb_steps - 1 ) * dt / 3600
    index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
    
############################## ECI X
########################################################################################################################################################################################################################################################################################################
if "x_eci_ref" in sys.argv:
    # Plot X ECI of reference spacecraft as a function of time
    hour_time_step_xticks = 24. # !!!!!! overwrite previous value
    var_out_list_in = ["position"]
    var_out, var_out_list = read_output_file(satellite_to_plot_path[0] + satellite_to_plot[0], var_out_list_in)
    eci_ref = var_out[1]
    x_eci_ref = eci_ref[:,0]
    nb_steps_ref = len(x_eci_ref) # !!!!!!!!! sometimes, nb_steps_ref is different from nb_steps because the output of the ensembles can be chosen only for a certain amount of time (and not the full propagation) (for example when collision asssessment is made, the output of ensembles is only during the time spanning the closest approach)

    ratio_fig_size = 4./3
    nb_steps_ref_adjusted = (int) ( nb_steps_ref * dt / 3600 / step_std) 
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left= 0.14, right=0.935, top = 0.93,bottom = 0.12)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('X ECI of reference sc VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
    ax1.set_ylabel('X ECI (km)', fontsize = fontsize_plot, weight = 'bold')
    ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.plot(np.arange(0,nb_steps_ref_adjusted), x_eci_ref, 'k', linewidth = 2)

    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)



    hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / step_std
    xticks = np.arange(0, nb_steps_ref_adjusted, hour_time_step_xticks_converted_in_index_adjusted)
    date_list_str = []
    nb_hours_simu = nb_steps_ref * dt/ 3600.
    date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
    for i in range(len(xticks)):
        if hour_time_step_xticks < 12:
            if i == 0:
                date_list_str.append("h+" + str(int(xticks[i] * step_std)))
            else:
                date_list_str.append("+" + str(int(xticks[i] * step_std)))
        else:
            date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] * step_std / 24.)) + ")")
    ax1.xaxis.set_ticks(xticks)
    ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')


    if save_results == 1:
        fig_save_name = 'x_eci_ref_sc_vs_time'
        fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
        fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
#        os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)



                                             
if "x_eci" in sys.argv:
    if ( 'x_eci_ensemble' in pickle_list_name ):
        # Plot the distribution of the ECI X at different times (set by when_plot_in_hour)
        # when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
        

        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

        med = np.median(x_eci_ensemble[index_when_plot, :])
        mad = np.median( np.abs( x_eci_ensemble[index_when_plot, :] - np.median(x_eci_ensemble[index_when_plot, :]) ) )  
        quartile10 = np.percentile(x_eci_ensemble[index_when_plot, :], 10) 
        quartile25 = np.percentile(x_eci_ensemble[index_when_plot, :], 25) 
        quartile75 = np.percentile(x_eci_ensemble[index_when_plot, :], 75) 
        quartile90 = np.percentile(x_eci_ensemble[index_when_plot, :], 90) 
        std = np.std(x_eci_ensemble[index_when_plot, :])
        iqr_daily = quartile75 - quartile25
        quartiles_1090_daily = quartile90 - quartile10

        ax2 = fig.add_subplot(gs[0, 0])
        ax2.set_title('ECI X distribution ' + 'after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
        bin_width = iqr_daily/10.
        hist_data = x_eci_ensemble[index_when_plot, :]
        n, bins, patches = ax2.hist(hist_data , bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'r',label = '')  
        plt.vlines(quartile25, 0, max(n),linewidth = 2); plt.vlines(quartile75, 0,max(n), linewidth = 2, label = '25-75% quartiles') 
        plt.vlines(quartile10, 0, max(n), linewidth = 4, linestyle = 'dotted'); plt.vlines(quartile90, 0, max(n), linewidth = 4, linestyle = 'dotted', label = '10-90% deciles') 
        ax2.plot(med, 0, 'k',marker = '.', markersize = 20)
        ax2.get_xaxis().get_major_formatter().set_useOffset(False)
        ax2.set_xlim([med - 10*iqr_daily, med+10*iqr_daily])
        ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
        ax2.text(ax2.get_xlim()[0] + ( ax2.get_xlim()[1] - ax2.get_xlim()[0] )/ 60., ax2.get_ylim()[1] - (ax2.get_ylim()[1] - ax2.get_ylim()[0])/5., 'std = ' + '{0:.2e} km'.format(std), fontsize = fontsize_plot)
        ax2.set_xlabel('ECI X (km)', fontsize = fontsize_plot, weight = 'bold')

        ax2.legend(fontsize = fontsize_plot, loc = 3, scatterpoints=1)
        if save_results == 1:
            fig_save_name = 'x_eci_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
#            os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

        # Plot the standard deviation of the X_Eci distributions as a function of time
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
            std_daily[i] = np.std(x_eci_ensemble[index_when_std, :]) # NOT GOOD TO USE THE STD
            iqr_daily[i] = np.subtract(*np.percentile(x_eci_ensemble[index_when_std, :], [75, 25])) 
            mad_daily[i] = np.median( np.abs( x_eci_ensemble[index_when_std, :] - np.median(x_eci_ensemble[index_when_std, :]) ) ) 

        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.11, right=0.935, top = 0.93,bottom = 0.12)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Standard deviation of the X ECI distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('Standard deviation (km)', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.plot(np.arange(0,nb_steps_adjusted), std_daily[:], 'k', linewidth = 2)

        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)


        
        hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / step_std
        xticks = np.arange(0, nb_steps_adjusted, hour_time_step_xticks_converted_in_index_adjusted)
        date_list_str = []
        nb_hours_simu = nb_steps * dt/ 3600.
        date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
        for i in range(len(xticks)):
            if hour_time_step_xticks < 12:
                if i == 0:
                    date_list_str.append("h+" + str(int(xticks[i] * step_std)))
                else:
                    date_list_str.append("+" + str(int(xticks[i] * step_std)))
            else:
                date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] * step_std / 24.)) + ")")
        ax1.xaxis.set_ticks(xticks)
        ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')

        
        if save_results == 1:
            fig_save_name = 'x_eci_std_daily'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


############################## ECI Y
########################################################################################################################################################################################################################################################################################################
if "y_eci_ref" in sys.argv:
    # Plot Y ECI of reference spacecraft as a function of time
    hour_time_step_xticks = 24. # !!!!!! overwrite previous value
    if (("x_eci_ref" in sys.argv) == False):
        var_out_list_in = ["position"]
        var_out, var_out_list = read_output_file(satellite_to_plot_path[0] + satellite_to_plot[0], var_out_list_in)
        eci_ref = var_out[1]
        y_eci_ref = eci_ref[:,1]
        nb_steps_ref = len(y_eci_ref) # !!!!!!!!! sometimes, nb_steps_ref is different from nb_steps because the output of the ensembles can be chosen only for a certain amount of time (and not the full propagation) (for example when collision asssessment is made, the output of ensembles is only during the time spanning the closest approach)
    else:
        y_eci_ref = eci_ref[:,1]
    

    ratio_fig_size = 4./3
    nb_steps_ref_adjusted = (int) ( nb_steps_ref * dt / 3600 / step_std) 
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left= 0.14, right=0.935, top = 0.93,bottom = 0.12)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Y ECI of reference sc VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
    ax1.set_ylabel('Y ECI (km)', fontsize = fontsize_plot, weight = 'bold')
    ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.plot(np.arange(0,nb_steps_ref_adjusted), y_eci_ref, 'k', linewidth = 2)

    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)



    hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / step_std
    xticks = np.arange(0, nb_steps_ref_adjusted, hour_time_step_xticks_converted_in_index_adjusted)
    date_list_str = []
    nb_hours_simu = nb_steps_ref * dt/ 3600.
    date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
    for i in range(len(xticks)):
        if hour_time_step_xticks < 12:
            if i == 0:
                date_list_str.append("h+" + str(int(xticks[i] * step_std)))
            else:
                date_list_str.append("+" + str(int(xticks[i] * step_std)))
        else:
            date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] * step_std / 24.)) + ")")
    ax1.xaxis.set_ticks(xticks)
    ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')


    if save_results == 1:
        fig_save_name = 'y_eci_ref_sc_vs_time'
        fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
        fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
        #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)



if "y_eci" in sys.argv:
    if ( 'y_eci_ensemble' in pickle_list_name ):
        # Plot the distribution of the ECI X at different times (set by when_plot_in_hour)
        # when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
        

        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

        med = np.median(y_eci_ensemble[index_when_plot, :])
        mad = np.median( np.abs( y_eci_ensemble[index_when_plot, :] - np.median(y_eci_ensemble[index_when_plot, :]) ) )  
        quartile10 = np.percentile(y_eci_ensemble[index_when_plot, :], 10) 
        quartile25 = np.percentile(y_eci_ensemble[index_when_plot, :], 25) 
        quartile75 = np.percentile(y_eci_ensemble[index_when_plot, :], 75) 
        quartile90 = np.percentile(y_eci_ensemble[index_when_plot, :], 90) 
        std = np.std(y_eci_ensemble[index_when_plot, :])
        iqr_daily = quartile75 - quartile25
        quartiles_1090_daily = quartile90 - quartile10

        ax2 = fig.add_subplot(gs[0, 0])
        ax2.set_title('ECI Y distribution ' + 'after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
        bin_width = iqr_daily/10.
        hist_data = y_eci_ensemble[index_when_plot, :]
        n, bins, patches = ax2.hist(hist_data , bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'r',label = '')  
        plt.vlines(quartile25, 0, max(n),linewidth = 2); plt.vlines(quartile75, 0,max(n), linewidth = 2, label = '25-75% quartiles') 
        plt.vlines(quartile10, 0, max(n), linewidth = 4, linestyle = 'dotted'); plt.vlines(quartile90, 0, max(n), linewidth = 4, linestyle = 'dotted', label = '10-90% deciles') 
        ax2.plot(med, 0, 'k',marker = '.', markersize = 20)
        ax2.get_xaxis().get_major_formatter().set_useOffset(False)
        ax2.set_xlim([med - 10*iqr_daily, med+10*iqr_daily])
        ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
        ax2.text(ax2.get_xlim()[0] + ( ax2.get_xlim()[1] - ax2.get_xlim()[0] )/ 60., ax2.get_ylim()[1] - (ax2.get_ylim()[1] - ax2.get_ylim()[0])/5., 'std = ' + '{0:.2e} km'.format(std), fontsize = fontsize_plot)
        ax2.set_xlabel('ECI Y (km)', fontsize = fontsize_plot, weight = 'bold')

        ax2.legend(fontsize = fontsize_plot, loc = 3, scatterpoints=1)
        if save_results == 1:
            fig_save_name = 'y_eci_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


        # Plot the standard deviation of the Y_Eci distributions as a function of time
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
            std_daily[i] = np.std(y_eci_ensemble[index_when_std, :]) # NOT GOOD TO USE THE STD
            iqr_daily[i] = np.subtract(*np.percentile(y_eci_ensemble[index_when_std, :], [75, 25])) 
            mad_daily[i] = np.median( np.abs( y_eci_ensemble[index_when_std, :] - np.median(y_eci_ensemble[index_when_std, :]) ) ) 

        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.11, right=0.935, top = 0.93,bottom = 0.12)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Standard deviation of the Y ECI distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('Standard deviation (km)', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.plot(np.arange(0,nb_steps_adjusted), std_daily[:], 'k', linewidth = 2)

        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)


        
        hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / step_std
        xticks = np.arange(0, nb_steps_adjusted, hour_time_step_xticks_converted_in_index_adjusted)
        date_list_str = []
        nb_hours_simu = nb_steps * dt/ 3600.
        date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
        for i in range(len(xticks)):
            if hour_time_step_xticks < 12:
                if i == 0:
                    date_list_str.append("h+" + str(int(xticks[i] * step_std)))
                else:
                    date_list_str.append("+" + str(int(xticks[i] * step_std)))
            else:
                date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] * step_std / 24.)) + ")")
        ax1.xaxis.set_ticks(xticks)
        ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')

        
        if save_results == 1:
            fig_save_name = 'y_eci_std_daily'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

############################## ECI Z
########################################################################################################################################################################################################################################################################################################
if "z_eci" in sys.argv:
    if ( 'z_eci_ensemble' in pickle_list_name ):
        # Plot the distribution of the ECI X at different times (set by when_plot_in_hour)
        # when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
        

        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

        med = np.median(z_eci_ensemble[index_when_plot, :])
        mad = np.median( np.abs( z_eci_ensemble[index_when_plot, :] - np.median(z_eci_ensemble[index_when_plot, :]) ) )  
        quartile10 = np.percentile(z_eci_ensemble[index_when_plot, :], 10) 
        quartile25 = np.percentile(z_eci_ensemble[index_when_plot, :], 25) 
        quartile75 = np.percentile(z_eci_ensemble[index_when_plot, :], 75) 
        quartile90 = np.percentile(z_eci_ensemble[index_when_plot, :], 90) 
        std = np.std(z_eci_ensemble[index_when_plot, :])
        iqr_daily = quartile75 - quartile25
        quartiles_1090_daily = quartile90 - quartile10

        ax2 = fig.add_subplot(gs[0, 0])
        ax2.set_title('ECI Z distribution ' + 'after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
        bin_width = iqr_daily/10.
        hist_data = z_eci_ensemble[index_when_plot, :]
        n, bins, patches = ax2.hist(hist_data , bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'r',label = '')  
        plt.vlines(quartile25, 0, max(n),linewidth = 2); plt.vlines(quartile75, 0,max(n), linewidth = 2, label = '25-75% quartiles') 
        plt.vlines(quartile10, 0, max(n), linewidth = 4, linestyle = 'dotted'); plt.vlines(quartile90, 0, max(n), linewidth = 4, linestyle = 'dotted', label = '10-90% deciles') 
        ax2.plot(med, 0, 'k',marker = '.', markersize = 20)
        ax2.get_xaxis().get_major_formatter().set_useOffset(False)
        ax2.set_xlim([med - 10*iqr_daily, med+10*iqr_daily])
        ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
        ax2.text(ax2.get_xlim()[0] + ( ax2.get_xlim()[1] - ax2.get_xlim()[0] )/ 60., ax2.get_ylim()[1] - (ax2.get_ylim()[1] - ax2.get_ylim()[0])/5., 'std = ' + '{0:.2e} km'.format(std), fontsize = fontsize_plot)
        ax2.set_xlabel('ECI Z (km)', fontsize = fontsize_plot, weight = 'bold')

        ax2.legend(fontsize = fontsize_plot, loc = 3, scatterpoints=1)
        if save_results == 1:
            fig_save_name = 'z_eci_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)







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
        gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Inter-quartile range of the along-track distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('IQR (km)', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.plot(np.arange(0,nb_steps_adjusted), iqr_daily[:], 'k', linewidth = 2) 


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

        if save_results == 1:
            fig_save_name = 'along_track_iqr_daily'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)
        

########################################################################################################################################################################################################################################################################################################
    if ( 'algebric_distance_ensemble_main_sc_lvlh_x' in pickle_list_name ):
        # Plot the distribution of angular_spacing_between_ensemble_sat_converted_in_a_distance at different times (set by when_plot_in_hour)
        # when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
        

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
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


############################## RADIAL SEPARATION WITH REFERENCE SPACECRAFT
########################################################################################################################################################################################################################################################################################################
if "radial" in sys.argv:
    algebric_distance_ensemble_main_sc_lvlh_z = algebric_distance_ensemble_main_sc_lvlh_z*1000. # !!!!!!!! conversion km to m (for cross track and radial, not for along track)
    if ( 'algebric_distance_ensemble_main_sc_lvlh_z' in pickle_list_name ):
        # Plot the standard deviation of the angular_spacing_between_ensemble_sat_converted_in_a_distance distributions as a function of time
        
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
        gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Inter-quartile range of the radial distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('IQR (m)', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.plot(np.arange(0,nb_steps_adjusted), iqr_daily[:], 'k', linewidth = 2) 


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


        if save_results == 1:
            fig_save_name = 'radial_iqr_daily'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)
        

########################################################################################################################################################################################################################################################################################################
    if ( 'algebric_distance_ensemble_main_sc_lvlh_z' in pickle_list_name ):
        # Plot the distribution of angular_spacing_between_ensemble_sat_converted_in_a_distance at different times (set by when_plot_in_hour)
        # when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
        

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
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)





############################## CROSS-TRACK SEPARATION WITH REFERENCE SPACECRAFT
########################################################################################################################################################################################################################################################################################################
if "cross" in sys.argv:
    algebric_distance_ensemble_main_sc_lvlh_y = algebric_distance_ensemble_main_sc_lvlh_y*1000. # !!!!!!!! conversion km to m (for cross track and radial, not for along track)
    if ( 'algebric_distance_ensemble_main_sc_lvlh_y' in pickle_list_name ):
        # Plot the standard deviation of the angular_spacing_between_ensemble_sat_converted_in_a_distance distributions as a function of time
        
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
        gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Inter-quartile range of the cross-track distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('IQR (m)', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.plot(np.arange(0,nb_steps_adjusted), iqr_daily[:], 'k', linewidth = 2) 


        
        hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / step_std
        xticks = np.arange(0, nb_steps_adjusted, hour_time_step_xticks_converted_in_index_adjusted)
        date_list_str = []
        nb_hours_simu = nb_steps * dt/ 3600.
        date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
        for i in range(len(xticks)):
            date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] * step_std / 24.)) + ")")
        ax1.xaxis.set_ticks(xticks)
        ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')


        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

        if save_results == 1:
            fig_save_name = 'cross_track_iqr_daily'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)
        

########################################################################################################################################################################################################################################################################################################
    if ( 'algebric_distance_ensemble_main_sc_lvlh_y' in pickle_list_name ):
        # Plot the distribution of angular_spacing_between_ensemble_sat_converted_in_a_distance at different times (set by when_plot_in_hour)
        # when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
        

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
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)




############################## ALONG-TRACK (USING ANGULAR SEPARATION)
########################################################################################################################################################################################################################################################################################################
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

        std_daily = np.zeros([nb_steps_adjusted])
        med_daily = np.zeros([nb_steps_adjusted])
        mad_daily = np.zeros([nb_steps_adjusted])
        iqr_daily = np.zeros([nb_steps_adjusted])
        for i in range(nb_steps_adjusted):
            index_when_std = i * step_std_in_index
            std_daily[i] = np.std(angle_asc_node_to_sat_ensemble[isat,index_when_std, :]) # NOT GOOD TO USE THE STD
            iqr_daily[i] = np.subtract(*np.percentile(angle_asc_node_to_sat_ensemble[isat,index_when_std, :], [75, 25])) 
            mad_daily[i] = np.median( np.abs( angle_asc_node_to_sat_ensemble[isat,index_when_std, :] - np.median(angle_asc_node_to_sat_ensemble[isat,index_when_std, :]) ) ) 

        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Inter-quartile range of the along-track distributions (using angular separation) VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('IQR (m)', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.plot(np.arange(0,nb_steps_adjusted), iqr_daily[:]*110000, 'k', linewidth = 2)

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

        if save_results == 1:
            fig_save_name = 'iqr_daily_along_track_using_angular_separation'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)
        

########################################################################################################################################################################################################################################################################################################
    if ( 'angle_asc_node_to_sat_ensemble' in pickle_list_name ):
        # Plot the distribution of angular_spacing_between_ensemble_sat_converted_in_a_distance at different times (set by when_plot_in_hour)
        # when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
        

        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('Spacecraft distribution along the orbit (using angular separation) ' + 'after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(4, 1)
        gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

        med = np.median(angle_asc_node_to_sat_ensemble[isat,index_when_plot, :])
        mad = np.median( np.abs( angle_asc_node_to_sat_ensemble[isat,index_when_plot, :] - np.median(angle_asc_node_to_sat_ensemble[isat,index_when_plot, :]) ) )  
        quartile10 = np.percentile(angle_asc_node_to_sat_ensemble[isat,index_when_plot, :], 10) 
        quartile25 = np.percentile(angle_asc_node_to_sat_ensemble[isat,index_when_plot, :], 25) 
        quartile75 = np.percentile(angle_asc_node_to_sat_ensemble[isat,index_when_plot, :], 75) 
        quartile90 = np.percentile(angle_asc_node_to_sat_ensemble[isat,index_when_plot, :], 90) 
        iqr_daily = quartile75 - quartile25
        quartiles_1090_daily = quartile90 - quartile10

        ax2 = fig.add_subplot(gs[:3, 0])
        ax2.set_title('', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
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
            fig_save_name = 'sc_along_track_distribution_using_angular_separation' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


# ########################################################################################################################################################################################################################################################################################################
#     if ( 'angle_asc_node_to_sat_ensemble' in pickle_list_name ):
#         # NOT FINISHED Plot the distribution of angular_spacing_between_ensemble_sat_converted_in_a_distance at different times on the same plot
      
#         ratio_fig_size = 4./3   

#         fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
#         fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
#         plt.rc('font', weight='bold') ## make the labels of the ticks in bold
#         gs = gridspec.GridSpec(4, 1)
#         gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

#         time_step_when_plot_in_hour = 48. # in hour
#         nb_steps_when_plot_in_hour = (int)( np.ceil(nb_steps * dt / (time_step_when_plot_in_hour * 3600. )) )

#         ax2 = fig.add_subplot(gs[:3, 0])
#         ax2.set_title('Distribution of the angular distance every ' + str((int) (time_step_when_plot_in_hour)) + ' hours of propagation', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 

#         ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
#         ax2.set_xlabel('', fontsize = fontsize_plot, weight = 'bold')
#         ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
#         [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
#         bin_width = iqr_daily/10.
#         hist_data = angle_asc_node_to_sat_ensemble[isat,index_when_plot, :]
#         n, bins, patches = ax2.hist(hist_data , bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'r',label = '')  
#         plt.vlines(quartile25, min(n), max(n),linewidth = 2); plt.vlines(quartile75, min(n),max(n), linewidth = 2) 
#         plt.vlines(quartile10, min(n), max(n), linewidth = 4, linestyle = 'dotted'); plt.vlines(quartile90, min(n), max(n), linewidth = 4, linestyle = 'dotted') 
#         ax2.plot(med, min(n), 'k',marker = '.', markersize = 20)
#         ax2.get_xaxis().get_major_formatter().set_useOffset(False)
#         ax2.set_xlim([med - 10*iqr_daily, med+10*iqr_daily])
#         ax2.get_xaxis().set_ticklabels([])
#         ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

#         ax1 = fig.add_subplot(gs[3, 0])
#         ymin = -0.2; ymax = 0.3
#         length_vert_iqr_daily = np.abs(ymin) / 3.
#         ax1.set_title('', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
#         ax1.set_ylabel('', fontsize = fontsize_plot, weight = 'bold')
#         ax1.set_xlabel('Angular distance from the AN (' + u'\N{DEGREE SIGN}'+ ')', fontsize = fontsize_plot, weight = 'bold')
#         ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
#         [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

#         ax1.hlines(0, med - 10 * iqr_daily, med + 10 * iqr_daily, linewidth = 2, color = 'b');
#         ax1.text( med + 10 * iqr_daily, 0 - length_vert_iqr_daily/0.8,'Orbit', horizontalalignment = 'right', fontsize = fontsize_plot, color = 'b' )

#         ax1.plot(angle_asc_node_to_sat_ensemble[isat, index_when_plot, :], np.zeros([nb_ensembles]), 'r.', markersize = 10)
#         ax1.vlines(quartile25, 0,length_vert_iqr_daily, linewidth = 2); ax1.vlines(quartile75, 0,length_vert_iqr_daily, linewidth = 2) 
#         ax1.text((quartile75 + quartile25 ) /2., length_vert_iqr_daily+length_vert_iqr_daily/5, '{0:.2f}'.format(iqr_daily * 110) + ' km (50% of sc)', horizontalalignment = 'center', fontsize = fontsize_plot )
#         ax1.annotate ('', (quartile25, length_vert_iqr_daily), (quartile75, length_vert_iqr_daily), arrowprops={'arrowstyle':'<->'})
#         ax1.vlines(quartile10, -length_vert_iqr_daily ,0 , linewidth = 4, linestyle = 'dotted'); ax1.vlines(quartile90, -length_vert_iqr_daily ,0, linewidth = 4, linestyle = 'dotted') 
#         ax1.annotate ('', (quartile10, -length_vert_iqr_daily), (quartile90, -length_vert_iqr_daily), arrowprops={'arrowstyle':'<->'})
#         ax1.text((quartile10 + quartile90 ) /2., -length_vert_iqr_daily-length_vert_iqr_daily/0.8, '{0:.2f}'.format(quartiles_1090_daily * 110) + ' km (80% of sc)', horizontalalignment = 'center', fontsize = fontsize_plot, verticalalignment = 'bottom' )
#         ax1.plot(med, 0, 'k',marker = '.', markersize = 15)

#         ax1.legend(fontsize = fontsize_plot, loc = 2, scatterpoints=1)
#         ax1.get_xaxis().tick_bottom()
#         ax1.set_ylim([ymin, ymax])
#     #    ax1.get_yaxis().tick_left()
#         ax1.yaxis.set_visible(False)
#         ax1.spines['left'].set_visible(False); ax1.spines['right'].set_visible(False); ax1.spines['top'].set_visible(False)
#         ax1.get_xaxis().get_major_formatter().set_useOffset(False)
#         ax1.set_xlim([med - 10*iqr_daily, med+10*iqr_daily])
#     #    ax1.set_xlim([min(bins), 0])
#         ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

#         if save_results == 1:
#             fig_save_name = ''
#             fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
#             fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
#             #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


    ############################## ORBIT AVERAGE DENSITY
    ########################################################################################################################################################################################################################################################################################################
if "rho" in sys.argv:
    # if ( 'rho_ensemble' in pickle_list_name ):
    #     nb_steps_orbit_average = len(index_time_orbit_average[:,1])     # !!!!!!!! there might be ore periods for satellites than for others (because they move at different speeds). So nb_steps_orbit_average should actually be the min of the number of periods of all satellites. Since we don't take the min, there might be 0 in rho_ensemble_orbit_average that we should NOT take into account. So basically NEVER TAKE THE STANDARD DEVIATION for this density section of the script (that takes into account outliers)
    #     index_time_orbit_average_middle = index_time_orbit_average[:,1]
    #     x_axis_orbit_average = index_time_orbit_average_middle * dt # conversion index to seconds
    #     x_axis_orbit_average = x_axis_orbit_average / (24*3600.) # conversion seconds to day

    #     # Plot the inter-quartile rang of the density distributions as a function of time
    #     step_std_in_index = 1 # we want the step to be every orbit (here we are look at orbit average data)
    #     std_every_orbit = np.zeros([nb_steps_orbit_average]) # DON'T TAKE IT FOR THE DENSITY SECTION
    #     med_every_orbit = np.zeros([nb_steps_orbit_average])
    #     mad_every_orbit = np.zeros([nb_steps_orbit_average])
    #     iqr_every_orbit = np.zeros([nb_steps_orbit_average])
    #     for i in range(nb_steps_orbit_average):
    #         index_when_std = i * step_std_in_index
    #         std_every_orbit[i] = np.std(rho_ensemble_orbit_average[index_when_std, :]) # NOT GOOD TO USE THE STD
    #         iqr_every_orbit[i] = np.subtract(*np.percentile(rho_ensemble_orbit_average[index_when_std, :], [75, 25])) 
    #         mad_every_orbit[i] = np.median( np.abs( rho_ensemble_orbit_average[index_when_std, :] - np.median(rho_ensemble_orbit_average[index_when_std, :]) ) ) 


    #     ## Plot
    #     ratio_fig_size = 4./3
    #     fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    #     fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
    #     plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    #     gs = gridspec.GridSpec(1, 1)
    #     gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
    #     ax1 = fig.add_subplot(gs[0, 0])
    #     ax1.set_title('Inter-quartile range of the density distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
    #     ax1.set_ylabel('IQR (#/m^3)', fontsize = fontsize_plot, weight = 'bold')
    #     ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
    #     ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    #     [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    #     ax1.plot(x_axis_orbit_average, iqr_every_orbit, 'k', linewidth = 2)

    #     ax1.get_xaxis().tick_bottom()
    #     ax1.get_yaxis().tick_left()
    #     ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)


    #     nb_days_simu = ( date_stop - date_start ).total_seconds()/3600./24.
    #     nb_hours_simu = nb_steps * dt/ 3600.
    #     hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / 24. # x_axis_orbit_average is in unit of days
    #     xticks = np.arange(0, nb_days_simu, hour_time_step_xticks_converted_in_index_adjusted)
    #     date_list_str = []
    #     date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
    #     for i in range(len(xticks)):
    #         date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] )) + ")")
    #         ax1.xaxis.set_ticks(xticks)
    #         ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')

        

    #     # nb_days_simu = ( date_stop - date_start ).total_seconds()/3600./24.
    #     # nb_hours_simu = nb_steps * dt/ 3600.
    #     # hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / 24. # x_axis_orbit_average is in unit of days
    #     # xticks = np.arange(0, nb_days_simu+1, hour_time_step_xticks_converted_in_index_adjusted)
    #     # date_list_str = []
    #     # date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
    #     # for i in range(len(xticks)):
    #     #     date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] )) + ")")
    #     # ax1.xaxis.set_ticks(xticks)
    #     # ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')


    #     if save_results == 1:
    #         fig_save_name = 'rho_iqr_every_orbit'
    #         fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
    #         fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    #         #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

    ########################################################################################################################################################################################################################################################################################################
    # if ( 'rho_ensemble_orbit_average' in pickle_list_name ):
    #     # Plot the distribution of the density at different times (set by when_plot_in_hour)
    #     # when_plot_in_hour = 2 * 24. # the few lines below convert this time in orbits 
    #     time_between_index_and_when_plot_in_hour = np.zeros([len(index_time_orbit_average_middle)])
    #     icount = -1
    #     for i in index_time_orbit_average_middle:
    #         icount = icount + 1
    #         time_between_index_and_when_plot_in_hour[icount] = np.abs( when_plot_in_hour * 3600. - i * dt )

    #     index_when_plot = (int)( np.where(time_between_index_and_when_plot_in_hour == min(time_between_index_and_when_plot_in_hour))[0][0] )
    #     nb_hour_this_orbit = index_time_orbit_average_middle[index_when_plot] * dt/3600./24

    #     ratio_fig_size = 4./3
    #     fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    #     fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    #     plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    #     gs = gridspec.GridSpec(1, 1)
    #     gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

    #     med = np.median(rho_ensemble_orbit_average[index_when_plot, :])
    #     mad = np.median( np.abs( rho_ensemble_orbit_average[index_when_plot, :] - np.median(rho_ensemble_orbit_average[index_when_plot, :]) ) )  
    #     quartile10 = np.percentile(rho_ensemble_orbit_average[index_when_plot, :], 10) 
    #     quartile25 = np.percentile(rho_ensemble_orbit_average[index_when_plot, :], 25) 
    #     quartile75 = np.percentile(rho_ensemble_orbit_average[index_when_plot, :], 75) 
    #     quartile90 = np.percentile(rho_ensemble_orbit_average[index_when_plot, :], 90) 
    #     iqr_every_orbit = quartile75 - quartile25
    #     quartiles_1090_every_orbit = quartile90 - quartile10
    #     std = np.std(rho_ensemble_orbit_average[index_when_plot, :])

    #     ax2 = fig.add_subplot(gs[0, 0])
    #     ax2.set_title('Density distribution ' + 'after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
    #     ax2.set_xlabel('Density (#/m^3)', fontsize = fontsize_plot, weight = 'bold')
    #     ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
    #     ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    #     [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    #     bin_width = iqr_every_orbit/10.
    #     hist_data = rho_ensemble_orbit_average[index_when_plot, :]
    #     n, bins, patches = ax2.hist(hist_data , bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'r',label = '')  
    #     plt.vlines(quartile25, 0, max(n),linewidth = 2); plt.vlines(quartile75, 0,max(n), linewidth = 2, label = '25-75% quartiles') 
    #     plt.vlines(quartile10, 0, max(n), linewidth = 4, linestyle = 'dotted'); plt.vlines(quartile90, 0, max(n), linewidth = 4, linestyle = 'dotted', label = '10-90% deciles') 
    #     ax2.plot(med, 0, 'k',marker = '.', markersize = 20)
    #     ax2.get_xaxis().get_major_formatter().set_useOffset(False)
    #     ax2.set_xlim([med - 10*iqr_every_orbit, med+10*iqr_every_orbit])
    #     ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
    #     ax2.text(ax2.get_xlim()[0] + ( ax2.get_xlim()[1] - ax2.get_xlim()[0] )/ 60., ax2.get_ylim()[1] - (ax2.get_ylim()[1] - ax2.get_ylim()[0])/5., 'IQR = ' + '{0:.2e} #/m^3'.format(iqr_every_orbit), fontsize = fontsize_plot)

    #     ax2.legend(fontsize = fontsize_plot, loc = 2, scatterpoints=1)
    #     if save_results == 1:
    #         fig_save_name = 'rho_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
    #         fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
    #         fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    #         #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


    if ( 'rho_ensemble_orbit_average' in pickle_list_name ):
        # Plot the density of each ensemble as a function of time, as well as the median and the inter-quartile range
        nb_steps_orbit_average = len(index_time_orbit_average[:,1])     # !!!!!!!! there might be ore periods for satellites than for others (because they move at different speeds). So nb_steps_orbit_average should actually be the min of the number of periods of all satellites. Since we don't take the min, there might be 0 in rho_ensemble_orbit_average that we should NOT take into account. So basically NEVER TAKE THE STANDARD DEVIATION for this density section of the script (that takes into account outliers)
        index_time_orbit_average_middle = index_time_orbit_average[:,1]
        x_axis_orbit_average = index_time_orbit_average_middle * dt # conversion index to seconds
        x_axis_orbit_average = x_axis_orbit_average / (24*3600.) # conversion seconds to day

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
        gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Density at the position of each ensemble member VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('Density (#/m^3)', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        for iens in np.arange(15,nb_ensembles,(int)(nb_ensembles/5)):#np.arange(0,nb_ensembles_per_proc):
            if iens == 0:
                ax1.plot(x_axis_orbit_average, rho_ensemble_orbit_average[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 2, alpha = 1, label = 'Ensemble')
            else:
                ax1.plot(x_axis_orbit_average, rho_ensemble_orbit_average[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 2, alpha = 1)
            print rho_ensemble_orbit_average[0:nb_steps:step_std_in_index,iens]
        # ax1.plot(x_axis_orbit_average, med_every_orbit, 'k', linewidth = 3, label = 'Median')
        # ax1.plot(x_axis_orbit_average, quartile25, 'r', linewidth = 3, label = '25 and 75% quartiles')
        # ax1.plot(x_axis_orbit_average, quartile75, 'r', linewidth = 3)

        ax1.legend(fontsize = fontsize_plot, loc = 2)
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

        
        
        # nb_days_simu = ( date_stop - date_start ).total_seconds()/3600./24.
        # nb_hours_simu = nb_steps * dt/ 3600.
        # hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / 24. # x_axis_orbit_average is in unit of days
        # xticks = np.arange(0, nb_days_simu, hour_time_step_xticks_converted_in_index_adjusted)
        # date_list_str = []
        # date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
        # for i in range(len(xticks)):
        #     date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] )) + ")")
        # ax1.xaxis.set_ticks(xticks)
        # ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')



        if save_results == 1:
            fig_save_name = 'rho_all_ensemble'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)



    ############################## F10.7
    ########################################################################################################################################################################################################################################################################################################
if "f107" in sys.argv:
#    if ( 'f107_ensemble' in pickle_list_name ):
        # Plot the standard deviation of the f10.7 distributions as a function of time
        
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
        gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Standard deviation of the f10.7 distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('Standard deviation', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.plot(np.arange(0,nb_steps_adjusted), std_daily[:], 'k', linewidth = 2)

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


        if save_results == 1:
            fig_save_name = 'f107_std_daily'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

    ########################################################################################################################################################################################################################################################################################################
#    if ( 'f107_ensemble' in pickle_list_name ):
        # Plot the distribution of the f10.7 at different times (set by when_plot_in_hour)
        # when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
        

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

        ax2.legend(fontsize = fontsize_plot, loc = 3, scatterpoints=1)
        if save_results == 1:
            fig_save_name = 'f107_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


#    if ( 'f107_ensemble' in pickle_list_name ):
        # Plot the f10.7 for each ensemble as a function of time, as well as the median and the inter-quartile range
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
        gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('F10.7 for each ensemble member VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('F10.7', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        for iens in np.arange(15,nb_ensembles,(int)(nb_ensembles/5)): #np.arange(nProcs)*((int)((nb_ensembles_per_proc)/np.float(nProcs)) + nb_ensembles_per_proc ):
            if iens == 0:
                ax1.plot(np.arange(0,nb_steps_adjusted), f107_ensemble[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 2, alpha = 1, label = 'Ensemble')
            else:
                ax1.plot(np.arange(0,nb_steps_adjusted), f107_ensemble[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 2, alpha = 1)
        #ax1.plot(np.arange(0,nb_steps_adjusted), med_daily[:], 'k', linewidth = 3, label = 'Median')
        #ax1.plot(np.arange(0,nb_steps_adjusted), quartile25[:], 'r', linewidth = 3, label = '25 and 75% quartiles')
        #ax1.plot(np.arange(0,nb_steps_adjusted), quartile75[:], 'r', linewidth = 3)

        ax1.legend(fontsize = fontsize_plot, loc = 3)
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


        if save_results == 1:
            fig_save_name = 'f107_all_ensemble'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)



    ############################## Ap
    ########################################################################################################################################################################################################################################################################################################
if "ap" in sys.argv:
    #if ( 'ap_ensemble' in pickle_list_name ):
        # Plot the standard deviation of the Ap distributions as a function of time
        
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
        gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Standard deviation of the Ap distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('Standard deviation', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.plot(np.arange(0,nb_steps_adjusted), std_daily[:], 'k', linewidth = 2)

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


        if save_results == 1:
            fig_save_name = 'ap_std_daily'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

    ########################################################################################################################################################################################################################################################################################################
    #if ( 'ap_ensemble' in pickle_list_name ):
        # Plot the distribution of the Ap at different times (set by when_plot_in_hour)
        # when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
        

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

        ax2.legend(fontsize = fontsize_plot, loc = 3, scatterpoints=1)
        if save_results == 1:
            fig_save_name = 'ap_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


    #if ( 'ap_ensemble' in pickle_list_name ):
        # Plot the Ap for each ensemble as a function of time, as well as the median and the inter-quartile range
        
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
        gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Ap for each ensemble member VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('Ap', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        for iens in np.arange(15,nb_ensembles,(int)(nb_ensembles/5)): #np.arange(nProcs)*((int)((nb_ensembles_per_proc)/np.float(nProcs)) + nb_ensembles_per_proc ):
            if iens == 0:
                ax1.plot(np.arange(0,nb_steps_adjusted), ap_ensemble[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 2, alpha = 1, label = 'Ensemble')
            else:
                ax1.plot(np.arange(0,nb_steps_adjusted), ap_ensemble[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 2, alpha = 1)
        #ax1.plot(np.arange(0,nb_steps_adjusted), med_daily[:], 'k', linewidth = 3, label = 'Median')
        #ax1.plot(np.arange(0,nb_steps_adjusted), quartile25[:], 'r', linewidth = 3, label = '25 and 75% quartiles')
        #ax1.plot(np.arange(0,nb_steps_adjusted), quartile75[:], 'r', linewidth = 3)

        ax1.legend(fontsize = fontsize_plot, loc = 3)
        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left() 
        ax1.margins(0,1) # autoscale both axes(fist value is for the x axis, second value for the y axis)
        ax1.set_ylim([0, 35])

        
        hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / step_std
        xticks = np.arange(0, nb_steps_adjusted, hour_time_step_xticks_converted_in_index_adjusted)
        date_list_str = []
        nb_hours_simu = nb_steps * dt/ 3600.
        date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
        for i in range(len(xticks)):
            date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] * step_std / 24.)) + ")")
        ax1.xaxis.set_ticks(xticks)
        ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')


        if save_results == 1:
            fig_save_name = 'ap_all_ensemble'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


    ############################## F10.7A
    ########################################################################################################################################################################################################################################################################################################
if "f107a" in sys.argv:
    #if ( 'f107a_ensemble' in pickle_list_name ):
        # Plot the standard deviation of the F10.7A distributions as a function of time
        
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
        gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Standard deviation of the F10.7A distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('Standard deviation', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.plot(np.arange(0,nb_steps_adjusted), std_daily[:], 'k', linewidth = 2)

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


        if save_results == 1:
            fig_save_name = 'f107a_std_daily'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

    ########################################################################################################################################################################################################################################################################################################
    #if ( 'f107a_ensemble' in pickle_list_name ):
        # Plot the distribution of the F10.7A at different times (set by when_plot_in_hour)
        # when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
         

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

        ax2.legend(fontsize = fontsize_plot, loc = 3, scatterpoints=1)
        if save_results == 1:
            fig_save_name = 'f107a_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


    #if ( 'f107a_ensemble' in pickle_list_name ):
        # Plot the F10.7A for each ensemble as a function of time, as well as the median and the inter-quartile range
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
        gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)     
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('F10.7A for each ensemble member VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('F10.7A', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        for iens in np.arange(15,nb_ensembles,(int)(nb_ensembles/5)): #np.arange(nProcs)*((int)((nb_ensembles_per_proc)/np.float(nProcs)) + nb_ensembles_per_proc ):
            if iens == 0:
                ax1.plot(np.arange(0,nb_steps_adjusted), f107a_ensemble[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 2, alpha = 1, label = 'Ensemble')
            else:
                ax1.plot(np.arange(0,nb_steps_adjusted), f107a_ensemble[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 2, alpha = 1)
        #ax1.plot(np.arange(0,nb_steps_adjusted), med_daily[:], 'k', linewidth = 3, label = 'Median')
        #ax1.plot(np.arange(0,nb_steps_adjusted), quartile25[:], 'r', linewidth = 3, label = '25 and 75% quartiles')
        #ax1.plot(np.arange(0,nb_steps_adjusted), quartile75[:], 'r', linewidth = 3)

        ax1.legend(fontsize = fontsize_plot, loc = 3 )
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


        if save_results == 1:
            fig_save_name = 'f107a_all_ensemble'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)




    ############################## PERIOD
    ########################################################################################################################################################################################################################################################################################################
if "period" in sys.argv:
#    if ( 'period_ensemble' in pickle_list_name ):
        # Plot the standard deviation of the f10.7 distributions as a function of time
        
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
            std_daily[i] = np.std(period_ensemble[index_when_std, :]) # NOT GOOD TO USE THE STD
            iqr_daily[i] = np.subtract(*np.percentile(period_ensemble[index_when_std, :], [75, 25])) 
            mad_daily[i] = np.median( np.abs( period_ensemble[index_when_std, :] - np.median(period_ensemble[index_when_std, :]) ) ) 

        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Standard deviation of the f10.7 distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('Standard deviation', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.plot(np.arange(0,nb_steps_adjusted), std_daily[:], 'k', linewidth = 2)

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


        if save_results == 1:
            fig_save_name = 'period_std_daily'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

    ########################################################################################################################################################################################################################################################################################################
#    if ( 'period_ensemble' in pickle_list_name ):
        # Plot the distribution of the f10.7 at different times (set by when_plot_in_hour)
        # when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
        

        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

        med = np.median(period_ensemble[index_when_plot, :])
        mad = np.median( np.abs( period_ensemble[index_when_plot, :] - np.median(period_ensemble[index_when_plot, :]) ) )  
        quartile10 = np.percentile(period_ensemble[index_when_plot, :], 10) 
        quartile25 = np.percentile(period_ensemble[index_when_plot, :], 25) 
        quartile75 = np.percentile(period_ensemble[index_when_plot, :], 75) 
        quartile90 = np.percentile(period_ensemble[index_when_plot, :], 90) 
        std = np.std(period_ensemble[index_when_plot, :])
        iqr_daily = quartile75 - quartile25
        quartiles_1090_daily = quartile90 - quartile10

        ax2 = fig.add_subplot(gs[0, 0])
        ax2.set_title('Period distribution ' + 'after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
        bin_width = iqr_daily/10.
        hist_data = period_ensemble[index_when_plot, :]
        n, bins, patches = ax2.hist(hist_data , bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'r',label = '')  
        plt.vlines(quartile25, 0, max(n),linewidth = 2); plt.vlines(quartile75, 0,max(n), linewidth = 2, label = '25-75% quartiles') 
        plt.vlines(quartile10, 0, max(n), linewidth = 4, linestyle = 'dotted'); plt.vlines(quartile90, 0, max(n), linewidth = 4, linestyle = 'dotted', label = '10-90% deciles') 
        ax2.plot(med, 0, 'k',marker = '.', markersize = 20)
        ax2.get_xaxis().get_major_formatter().set_useOffset(False)
        ax2.set_xlim([med - 10*iqr_daily, med+10*iqr_daily])
        ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
        ax2.text(ax2.get_xlim()[0] + ( ax2.get_xlim()[1] - ax2.get_xlim()[0] )/ 60., ax2.get_ylim()[1] - (ax2.get_ylim()[1] - ax2.get_ylim()[0])/5., 'std = ' + '{0:.2f}'.format(std), fontsize = fontsize_plot)
        ax2.set_xlabel('Period', fontsize = fontsize_plot, weight = 'bold')

        ax2.legend(fontsize = fontsize_plot, loc = 3, scatterpoints=1)
        if save_results == 1:
            fig_save_name = 'period_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


#    if ( 'period_ensemble' in pickle_list_name ):
        # Plot the f10.7 for each ensemble as a function of time, as well as the median and the inter-quartile range
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
            iqr_daily[i] = np.subtract(*np.percentile(period_ensemble[index_when_std, :], [75, 25])) 
            mad_daily[i] = np.median( np.abs( period_ensemble[index_when_std, :] - np.median(period_ensemble[index_when_std, :]) ) ) 
            med_daily[i] = np.median(period_ensemble[index_when_std, :])
            quartile25[i] = np.percentile(period_ensemble[index_when_std, :], 25) 
            quartile75[i] = np.percentile(period_ensemble[index_when_std, :], 75) 

        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Period for each ensemble member VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('Period', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        for iens in np.arange(15,nb_ensembles,(int)(nb_ensembles/5)): #np.arange(nProcs)*((int)((nb_ensembles_per_proc)/np.float(nProcs)) + nb_ensembles_per_proc ):
            if iens == 0:
                ax1.plot(np.arange(0,nb_steps_adjusted), period_ensemble[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 2, alpha = 1, label = 'Ensemble')
            else:
                ax1.plot(np.arange(0,nb_steps_adjusted), period_ensemble[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 2, alpha = 1)
        #ax1.plot(np.arange(0,nb_steps_adjusted), med_daily[:], 'k', linewidth = 3, label = 'Median')
        #ax1.plot(np.arange(0,nb_steps_adjusted), quartile25[:], 'r', linewidth = 3, label = '25 and 75% quartiles')
        #ax1.plot(np.arange(0,nb_steps_adjusted), quartile75[:], 'r', linewidth = 3)

        ax1.legend(fontsize = fontsize_plot, loc = 3)
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


        if save_results == 1:
            fig_save_name = 'period_all_ensemble'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)




    ############################## ORBIT AVERAGE PERIOD
    ########################################################################################################################################################################################################################################################################################################
if "period_average" in sys.argv:
    # if ( 'period_ensemble' in pickle_list_name ):
    #     nb_steps_orbit_average = len(index_time_orbit_average[:,1])     # !!!!!!!! there might be ore periods for satellites than for others (because they move at different speeds). So nb_steps_orbit_average should actually be the min of the number of periods of all satellites. Since we don't take the min, there might be 0 in period_ensemble_orbit_average that we should NOT take into account. So basically NEVER TAKE THE STANDARD DEVIATION for this period section of the script (that takes into account outliers)
    #     index_time_orbit_average_middle = index_time_orbit_average[:,1]
    #     x_axis_orbit_average = index_time_orbit_average_middle * dt # conversion index to seconds
    #     x_axis_orbit_average = x_axis_orbit_average / (24*3600.) # conversion seconds to day

    #     # Plot the inter-quartile rang of the period distributions as a function of time
    #     step_std_in_index = 1 # we want the step to be every orbit (here we are look at orbit average data)
    #     std_every_orbit = np.zeros([nb_steps_orbit_average]) # DON'T TAKE IT FOR THE PERIOD SECTION
    #     med_every_orbit = np.zeros([nb_steps_orbit_average])
    #     mad_every_orbit = np.zeros([nb_steps_orbit_average])
    #     iqr_every_orbit = np.zeros([nb_steps_orbit_average])
    #     for i in range(nb_steps_orbit_average):
    #         index_when_std = i * step_std_in_index
    #         std_every_orbit[i] = np.std(period_ensemble_orbit_average[index_when_std, :]) # NOT GOOD TO USE THE STD
    #         iqr_every_orbit[i] = np.subtract(*np.percentile(period_ensemble_orbit_average[index_when_std, :], [75, 25])) 
    #         mad_every_orbit[i] = np.median( np.abs( period_ensemble_orbit_average[index_when_std, :] - np.median(period_ensemble_orbit_average[index_when_std, :]) ) ) 


    #     ## Plot
    #     ratio_fig_size = 4./3
    #     fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    #     fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
    #     plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    #     gs = gridspec.GridSpec(1, 1)
    #     gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
    #     ax1 = fig.add_subplot(gs[0, 0])
    #     ax1.set_title('Inter-quartile range of the period distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
    #     ax1.set_ylabel('IQR (#/m^3)', fontsize = fontsize_plot, weight = 'bold')
    #     ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
    #     ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    #     [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    #     ax1.plot(x_axis_orbit_average, iqr_every_orbit, 'k', linewidth = 2)

    #     ax1.get_xaxis().tick_bottom()
    #     ax1.get_yaxis().tick_left()
    #     ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)


        
    #     # nb_days_simu = ( date_stop - date_start ).total_seconds()/3600./24.
    #     # nb_hours_simu = nb_steps * dt/ 3600.
    #     # hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / 24. # x_axis_orbit_average is in unit of days
    #     # xticks = np.arange(0, nb_days_simu+1, hour_time_step_xticks_converted_in_index_adjusted)
    #     # date_list_str = []
    #     # date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
    #     # for i in range(len(xticks)):
    #     #     date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] )) + ")")
    #     # ax1.xaxis.set_ticks(xticks)
    #     # ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')


    #     if save_results == 1:
    #         fig_save_name = 'period_iqr_every_orbit'
    #         fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
    #         fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    #         #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

    ########################################################################################################################################################################################################################################################################################################
    if ( 'period_ensemble_orbit_average' in pickle_list_name ):
        # Plot the distribution of the period at different times (set by when_plot_in_hour)
        # when_plot_in_hour = 2 * 24. # the few lines below convert this time in orbits 
        index_time_orbit_average_middle = index_time_orbit_average[:,1]
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

        med = np.median(period_ensemble_orbit_average[index_when_plot, :])
        mad = np.median( np.abs( period_ensemble_orbit_average[index_when_plot, :] - np.median(period_ensemble_orbit_average[index_when_plot, :]) ) )  
        quartile10 = np.percentile(period_ensemble_orbit_average[index_when_plot, :], 10) 
        quartile25 = np.percentile(period_ensemble_orbit_average[index_when_plot, :], 25) 
        quartile75 = np.percentile(period_ensemble_orbit_average[index_when_plot, :], 75) 
        quartile90 = np.percentile(period_ensemble_orbit_average[index_when_plot, :], 90) 
        iqr_every_orbit = quartile75 - quartile25
        quartiles_1090_every_orbit = quartile90 - quartile10
        std = np.std(period_ensemble_orbit_average[index_when_plot, :])

        ax2 = fig.add_subplot(gs[0, 0])
        ax2.set_title('Period distribution ' + 'after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation  - CYGNSS' + str(isat + 1), weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax2.set_xlabel('Period (s)', fontsize = fontsize_plot, weight = 'bold')
        ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
        bin_width = iqr_every_orbit/10.
        hist_data = period_ensemble_orbit_average[index_when_plot, :]
        n, bins, patches = ax2.hist(hist_data , bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'm',label = '')  
        plt.vlines(quartile25, 0, max(n),linewidth = 2); plt.vlines(quartile75, 0,max(n), linewidth = 2, label = '25-75% quartiles') 
        plt.vlines(quartile10, 0, max(n), linewidth = 4, linestyle = 'dotted'); plt.vlines(quartile90, 0, max(n), linewidth = 4, linestyle = 'dotted', label = '10-90% deciles') 
        ax2.plot(med, 0, 'k',marker = '.', markersize = 20)
        ax2.get_xaxis().get_major_formatter().set_useOffset(False)
        ax2.set_xlim([med - 10*iqr_every_orbit, med+10*iqr_every_orbit])
        ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
#        ax2.text(ax2.get_xlim()[0] + ( ax2.get_xlim()[1] - ax2.get_xlim()[0] )/ 60., ax2.get_ylim()[1] - (ax2.get_ylim()[1] - ax2.get_ylim()[0])/5., 'IQR = ' + '{0:.2e} #/m^3'.format(iqr_every_orbit), fontsize = fontsize_plot)

        ax2.legend(fontsize = fontsize_plot, loc = 2, scatterpoints=1)
        if save_results == 1:
            fig_save_name = 'period_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment_isat_' + str(isat)
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


    # if ( 'period_ensemble_orbit_average' in pickle_list_name ):
    #     # Plot the period of each ensemble as a function of time, as well as the median and the inter-quartile range
    #     step_std_in_index = 1 # every orbit
    #     med_every_orbit = np.zeros([nb_steps_orbit_average])
    #     mad_every_orbit = np.zeros([nb_steps_orbit_average])
    #     iqr_every_orbit = np.zeros([nb_steps_orbit_average])
    #     quartile25 =np.zeros([nb_steps_orbit_average])
    #     quartile75 =np.zeros([nb_steps_orbit_average])
    #     for i in range(nb_steps_orbit_average):
    #         index_when_std = i * step_std_in_index
    #         iqr_every_orbit[i] = np.subtract(*np.percentile(period_ensemble_orbit_average[index_when_std, :], [75, 25])) 
    #         mad_every_orbit[i] = np.median( np.abs( period_ensemble_orbit_average[index_when_std, :] - np.median(period_ensemble_orbit_average[index_when_std, :]) ) ) 
    #         med_every_orbit[i] = np.median(period_ensemble_orbit_average[index_when_std, :])
    #         quartile25[i] = np.percentile(period_ensemble_orbit_average[index_when_std, :], 25) 
    #         quartile75[i] = np.percentile(period_ensemble_orbit_average[index_when_std, :], 75) 

    #     ## Plot
    #     ratio_fig_size = 4./3
    #     fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    #     fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
    #     plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    #     gs = gridspec.GridSpec(1, 1)
    #     gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
    #     ax1 = fig.add_subplot(gs[0, 0])
    #     ax1.set_title('Period at the position of each ensemble member VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
    #     ax1.set_ylabel('Period (s)', fontsize = fontsize_plot, weight = 'bold')
    #     ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
    #     ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    #     [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    #     for iens in range(nb_ensembles):
    #         if iens == 0:
    #             ax1.plot(x_axis_orbit_average, period_ensemble_orbit_average[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 2, alpha = 1, label = 'Ensemble')
    #         else:
    #             ax1.plot(x_axis_orbit_average, period_ensemble_orbit_average[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 2, alpha = 1)
    #     ax1.plot(x_axis_orbit_average, med_every_orbit, 'k', linewidth = 3, label = 'Median')
    #     ax1.plot(x_axis_orbit_average, quartile25, 'r', linewidth = 3, label = '25 and 75% quartiles')
    #     ax1.plot(x_axis_orbit_average, quartile75, 'r', linewidth = 3)

    #     ax1.legend(fontsize = fontsize_plot, loc = 2)
    #     ax1.get_xaxis().tick_bottom()
    #     ax1.get_yaxis().tick_left()
    #     ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)


        
    #     nb_days_simu = ( date_stop - date_start ).total_seconds()/3600./24.
    #     nb_hours_simu = nb_steps * dt/ 3600.
    #     hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / 24. # x_axis_orbit_average is in unit of days
    #     xticks = np.arange(0, nb_days_simu, hour_time_step_xticks_converted_in_index_adjusted)
    #     date_list_str = []
    #     date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
    #     for i in range(len(xticks)):
    #         date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] )) + ")")
    #     ax1.xaxis.set_ticks(xticks)
    #     ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')



    #     if save_results == 1:
    #         fig_save_name = 'period_all_ensemble'
    #         fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
    #         fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    #         #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)



    ############################## SMA
    ########################################################################################################################################################################################################################################################################################################
if "cd" in sys.argv:
    if ( 'cd_ensemble' in pickle_list_name ):
        # Plot the distribution of the Cd
        index_when_plot = 0
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

        med = np.median(cd_ensemble[index_when_plot, :])
        mad = np.median( np.abs( cd_ensemble[index_when_plot, :] - np.median(cd_ensemble[index_when_plot, :]) ) )  
        quartile10 = np.percentile(cd_ensemble[index_when_plot, :], 10) 
        quartile25 = np.percentile(cd_ensemble[index_when_plot, :], 25) 
        quartile75 = np.percentile(cd_ensemble[index_when_plot, :], 75) 
        quartile90 = np.percentile(cd_ensemble[index_when_plot, :], 90) 
        std = np.std(cd_ensemble[index_when_plot, :])
        iqr_daily = quartile75 - quartile25
        quartiles_1090_daily = quartile90 - quartile10

        ax2 = fig.add_subplot(gs[0, 0])
        ax2.set_title('Cd distribution ' + 'after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
        bin_width = iqr_daily/10.
        hist_data = cd_ensemble[index_when_plot, :]
        n, bins, patches = ax2.hist(hist_data , bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'r',label = '')  
        plt.vlines(quartile25, 0, max(n),linewidth = 2); plt.vlines(quartile75, 0,max(n), linewidth = 2, label = '25-75% quartiles') 
        plt.vlines(quartile10, 0, max(n), linewidth = 4, linestyle = 'dotted'); plt.vlines(quartile90, 0, max(n), linewidth = 4, linestyle = 'dotted', label = '10-90% deciles') 
        ax2.plot(med, 0, 'k',marker = '.', markersize = 20)
        ax2.get_xaxis().get_major_formatter().set_useOffset(False)
        ax2.set_xlim([med - 10*iqr_daily, med+10*iqr_daily])
        ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
        ax2.text(ax2.get_xlim()[0] + ( ax2.get_xlim()[1] - ax2.get_xlim()[0] )/ 60., ax2.get_ylim()[1] - (ax2.get_ylim()[1] - ax2.get_ylim()[0])/5., 'std = ' + '{0:.2f}'.format(std), fontsize = fontsize_plot)
        ax2.set_xlabel('Cd', fontsize = fontsize_plot, weight = 'bold')

        ax2.legend(fontsize = fontsize_plot, loc = 3, scatterpoints=1)
        if save_results == 1:
            fig_save_name = 'cd_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)




    ############################## SMA
    ########################################################################################################################################################################################################################################################################################################
if "sma" in sys.argv:
#    if ( 'sma_ensemble' in pickle_list_name ):
        # Plot the standard deviation of the f10.7 distributions as a function of time
        
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
            std_daily[i] = np.std(sma_ensemble[index_when_std, :]) # NOT GOOD TO USE THE STD
            iqr_daily[i] = np.subtract(*np.percentile(sma_ensemble[index_when_std, :], [75, 25])) 
            mad_daily[i] = np.median( np.abs( sma_ensemble[index_when_std, :] - np.median(sma_ensemble[index_when_std, :]) ) ) 

        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Standard deviation of the f10.7 distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('Standard deviation', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.plot(np.arange(0,nb_steps_adjusted), std_daily[:], 'k', linewidth = 2)

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


        if save_results == 1:
            fig_save_name = 'sma_std_daily'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

    ########################################################################################################################################################################################################################################################################################################
#    if ( 'sma_ensemble' in pickle_list_name ):
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
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


#    if ( 'sma_ensemble' in pickle_list_name ):
        # Plot the f10.7 for each ensemble as a function of time, as well as the median and the inter-quartile range
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
            iqr_daily[i] = np.subtract(*np.percentile(sma_ensemble[index_when_std, :], [75, 25])) 
            mad_daily[i] = np.median( np.abs( sma_ensemble[index_when_std, :] - np.median(sma_ensemble[index_when_std, :]) ) ) 
            med_daily[i] = np.median(sma_ensemble[index_when_std, :])
            quartile25[i] = np.percentile(sma_ensemble[index_when_std, :], 25) 
            quartile75[i] = np.percentile(sma_ensemble[index_when_std, :], 75) 

        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Sma for each ensemble member VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('Sma', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        for iens in np.arange(0, nb_ensembles, 1):#np.arange(15,nb_ensembles,(int)(nb_ensembles/5)):   #np.arange(nProcs)*((int)((nb_ensembles_per_proc)/np.float(nProcs)) + nb_ensembles_per_proc ):
            if iens == 0:
                ax1.plot(np.arange(0,nb_steps_adjusted), sma_ensemble[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 2, alpha = 1, label = 'Ensemble')
            else:
                ax1.plot(np.arange(0,nb_steps_adjusted), sma_ensemble[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 2, alpha = 1)
        #ax1.plot(np.arange(0,nb_steps_adjusted), med_daily[:], 'k', linewidth = 3, label = 'Median')
        #ax1.plot(np.arange(0,nb_steps_adjusted), quartile25[:], 'r', linewidth = 3, label = '25 and 75% quartiles')
        #ax1.plot(np.arange(0,nb_steps_adjusted), quartile75[:], 'r', linewidth = 3)

        ax1.legend(fontsize = fontsize_plot, loc = 3)
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


        if save_results == 1:
            fig_save_name = 'sma_all_ensemble'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)



    ############################## ORBIT AVERAGE SMA
    ########################################################################################################################################################################################################################################################################################################
if "sma_average" in sys.argv:
    if ( 'sma_ensemble' in pickle_list_name ):
        nb_steps_orbit_average = len(index_time_orbit_average[:,1])     # !!!!!!!! there might be ore smas for satellites than for others (because they move at different speeds). So nb_steps_orbit_average should actually be the min of the number of smas of all satellites. Since we don't take the min, there might be 0 in sma_ensemble_orbit_average that we should NOT take into account. So basically NEVER TAKE THE STANDARD DEVIATION for this sma section of the script (that takes into account outliers)
        index_time_orbit_average_middle = index_time_orbit_average[:,1]
        x_axis_orbit_average = index_time_orbit_average_middle * dt # conversion index to seconds
        x_axis_orbit_average = x_axis_orbit_average / (24*3600.) # conversion seconds to day

        # Plot the inter-quartile rang of the sma distributions as a function of time
        step_std_in_index = 1 # we want the step to be every orbit (here we are look at orbit average data)
        std_every_orbit = np.zeros([nb_steps_orbit_average]) # DON'T TAKE IT FOR THE SMA SECTION
        med_every_orbit = np.zeros([nb_steps_orbit_average])
        mad_every_orbit = np.zeros([nb_steps_orbit_average])
        iqr_every_orbit = np.zeros([nb_steps_orbit_average])
        for i in range(nb_steps_orbit_average):
            index_when_std = i * step_std_in_index
            std_every_orbit[i] = np.std(sma_ensemble_orbit_average[index_when_std, :]) # NOT GOOD TO USE THE STD
            iqr_every_orbit[i] = np.subtract(*np.percentile(sma_ensemble_orbit_average[index_when_std, :], [75, 25])) 
            mad_every_orbit[i] = np.median( np.abs( sma_ensemble_orbit_average[index_when_std, :] - np.median(sma_ensemble_orbit_average[index_when_std, :]) ) ) 


        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Inter-quartile range of the sma distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('IQR (#/m^3)', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.plot(x_axis_orbit_average, iqr_every_orbit, 'k', linewidth = 2)

        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)


        
        # nb_days_simu = ( date_stop - date_start ).total_seconds()/3600./24.
        # nb_hours_simu = nb_steps * dt/ 3600.
        # hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / 24. # x_axis_orbit_average is in unit of days
        # xticks = np.arange(0, nb_days_simu+1, hour_time_step_xticks_converted_in_index_adjusted)
        # date_list_str = []
        # date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
        # for i in range(len(xticks)):
        #     date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] )) + ")")
        # ax1.xaxis.set_ticks(xticks)
        # ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')


        if save_results == 1:
            fig_save_name = 'sma_iqr_every_orbit'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

    ########################################################################################################################################################################################################################################################################################################
    if ( 'sma_ensemble_orbit_average' in pickle_list_name ):
        # Plot the distribution of the sma at different times (set by when_plot_in_hour)
        # when_plot_in_hour = 2 * 24. # the few lines below convert this time in orbits 
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

        med = np.median(sma_ensemble_orbit_average[index_when_plot, :])
        mad = np.median( np.abs( sma_ensemble_orbit_average[index_when_plot, :] - np.median(sma_ensemble_orbit_average[index_when_plot, :]) ) )  
        quartile10 = np.percentile(sma_ensemble_orbit_average[index_when_plot, :], 10) 
        quartile25 = np.percentile(sma_ensemble_orbit_average[index_when_plot, :], 25) 
        quartile75 = np.percentile(sma_ensemble_orbit_average[index_when_plot, :], 75) 
        quartile90 = np.percentile(sma_ensemble_orbit_average[index_when_plot, :], 90) 
        iqr_every_orbit = quartile75 - quartile25
        quartiles_1090_every_orbit = quartile90 - quartile10
        std = np.std(sma_ensemble_orbit_average[index_when_plot, :])

        ax2 = fig.add_subplot(gs[0, 0])
        ax2.set_title('Sma distribution ' + 'after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax2.set_xlabel('Sma (#/m^3)', fontsize = fontsize_plot, weight = 'bold')
        ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
        bin_width = iqr_every_orbit/10.
        hist_data = sma_ensemble_orbit_average[index_when_plot, :]
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
            fig_save_name = 'sma_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


    if ( 'sma_ensemble_orbit_average' in pickle_list_name ):
        # Plot the sma of each ensemble as a function of time, as well as the median and the inter-quartile range
        step_std_in_index = 1 # every orbit
        med_every_orbit = np.zeros([nb_steps_orbit_average])
        mad_every_orbit = np.zeros([nb_steps_orbit_average])
        iqr_every_orbit = np.zeros([nb_steps_orbit_average])
        quartile25 =np.zeros([nb_steps_orbit_average])
        quartile75 =np.zeros([nb_steps_orbit_average])
        for i in range(nb_steps_orbit_average):
            index_when_std = i * step_std_in_index
            iqr_every_orbit[i] = np.subtract(*np.percentile(sma_ensemble_orbit_average[index_when_std, :], [75, 25])) 
            mad_every_orbit[i] = np.median( np.abs( sma_ensemble_orbit_average[index_when_std, :] - np.median(sma_ensemble_orbit_average[index_when_std, :]) ) ) 
            med_every_orbit[i] = np.median(sma_ensemble_orbit_average[index_when_std, :])
            quartile25[i] = np.percentile(sma_ensemble_orbit_average[index_when_std, :], 25) 
            quartile75[i] = np.percentile(sma_ensemble_orbit_average[index_when_std, :], 75) 

        ## Plot
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Sma at the position of each ensemble member VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('Sma (s)', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        for iens in range(nb_ensembles):
            if iens == 0:
                ax1.plot(x_axis_orbit_average, sma_ensemble_orbit_average[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 2, alpha = 1, label = 'Ensemble')
            else:
                ax1.plot(x_axis_orbit_average, sma_ensemble_orbit_average[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 2, alpha = 1)
        ax1.plot(x_axis_orbit_average, med_every_orbit, 'k', linewidth = 3, label = 'Median')
        ax1.plot(x_axis_orbit_average, quartile25, 'r', linewidth = 3, label = '25 and 75% quartiles')
        ax1.plot(x_axis_orbit_average, quartile75, 'r', linewidth = 3)

        ax1.legend(fontsize = fontsize_plot, loc = 2)
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



        if save_results == 1:
            fig_save_name = 'sma_all_ensemble'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)




    ############################## ANGULAR SPACING BETWEEN 2 CLUSTER
    ########################################################################################################################################################################################################################################################################################################

if "angle_cluster" in sys.argv:
    if isat == 1: # stuff were saved for isat = 1 in mpi_distance_ensemble_to_main_sc.py (because to compute phase_angle you first need to compute the period of isat 0 and then isat 1)
        if ( 'angular_distance_between_two_clusters' in pickle_list_name ):
            dt_angular_distance_between_two_clusters = 95. # we look at angular distance between the 2 clusters only at certain times otherwise it's a too heavy calculation. here in MINUTES. HAS TO BE THE SAME AS IN MPI_DISTANCE_ENSEMBLE_TO_MAIN_SC.PY
            index_dt_angular_distance_between_two_clusters = (int)( dt_angular_distance_between_two_clusters * 60 / dt )
            nb_steps_angular_distance_between_two_clusters = len(save_index_dt_angular_distance_between_two_clusters)
            iqr_angle_cluster = np.zeros([nb_steps_angular_distance_between_two_clusters])
            median_angle_cluster = np.zeros([nb_steps_angular_distance_between_two_clusters])
            for i in range(nb_steps_angular_distance_between_two_clusters):
                iqr_angle_cluster[i] = np.subtract(*np.percentile(angular_distance_between_two_clusters[i,:], [75, 25]))
                median_angle_cluster[i] = np.median(angular_distance_between_two_clusters[i,:])
            x_axis = np.arange(0,  nb_steps_angular_distance_between_two_clusters * dt_angular_distance_between_two_clusters*60, dt_angular_distance_between_two_clusters*60)

            ## Plot
            ratio_fig_size = 4./3
            fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
            fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
            plt.rc('font', weight='bold') ## make the labels of the ticks in bold
            gs = gridspec.GridSpec(1, 1)
            gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
            ax1 = fig.add_subplot(gs[0, 0])
            ax1.set_title('Inter-quartile range of the angular distance  VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
            ax1.set_ylabel('IQR of the angular distance (' + u'\N{DEGREE SIGN}'+ ')', fontsize = fontsize_plot, weight = 'bold')
            ax1.set_xlabel('Time (days)', fontsize = fontsize_plot, weight = 'bold')
            ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
            [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
            ax1.plot(x_axis/24./3600, iqr_angle_cluster, 'k', linewidth = 2)
            

            ax1.get_xaxis().tick_bottom()
            ax1.get_yaxis().tick_left()
            ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

            # nb_days_simu = ( date_stop - date_start ).total_seconds()/3600./24.
            # nb_hours_simu = nb_steps * dt/ 3600.
            # hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / 24. # x_axis_orbit_average is in unit of days
            # xticks = np.arange(0, nb_days_simu+1, hour_time_step_xticks_converted_in_index_adjusted)
            # date_list_str = []
            # date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
            # for i in range(len(xticks)):
            #     date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] )) + ")")
            # ax1.xaxis.set_ticks(xticks)
            # ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')


            if save_results == 1:
                fig_save_name = 'angle_cluster_iqr_every_orbit'
                fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
                #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


            ## Plot
            ratio_fig_size = 4./3
            fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
            fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
            plt.rc('font', weight='bold') ## make the labels of the ticks in bold
            gs = gridspec.GridSpec(1, 1)
            gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
            ax1 = fig.add_subplot(gs[0, 0])
            ax1.set_title('Median of the angular distance distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
            ax1.set_ylabel('Median of the angular distance (' + u'\N{DEGREE SIGN}'+ ')', fontsize = fontsize_plot, weight = 'bold')
            ax1.set_xlabel('Time (days)', fontsize = fontsize_plot, weight = 'bold')
            ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
            [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
            #ax1.plot(x_axis, iqr_angle_cluster, 'k', linewidth = 2)
            ax1.plot(x_axis/3600./24., median_angle_cluster, 'k', linewidth = 2)

            ax1.get_xaxis().tick_bottom()
            ax1.get_yaxis().tick_left()
            ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

            # nb_days_simu = ( date_stop - date_start ).total_seconds()/3600./24.
            # nb_hours_simu = nb_steps * dt/ 3600.
            # hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / 24. # x_axis_orbit_average is in unit of days
            # xticks = np.arange(0, nb_days_simu+1, hour_time_step_xticks_converted_in_index_adjusted)
            # date_list_str = []
            # date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
            # for i in range(len(xticks)):
            #     date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] )) + ")")
            # ax1.xaxis.set_ticks(xticks)
            # ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')


            if save_results == 1:
                fig_save_name = 'angle_cluster_median_every_orbit'
                fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
                #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

        ########################################################################################################################################################################################################################################################################################################
        # if ( 'angular_distance_between_two_clusters' in pickle_list_name ):
        #     # Plot the distribution of the phase_angle at different times (set by when_plot_in_hour)
        #     # when_plot_in_hour = 2 * 24. # the few lines below convert this time in orbits 
        #     time_between_index_and_when_plot_in_hour = np.zeros([len(index_time_orbit_average_middle)])
        #     icount = -1
        #     for i in index_time_orbit_average_middle:
        #         icount = icount + 1
        #         time_between_index_and_when_plot_in_hour[icount] = np.abs( when_plot_in_hour * 3600. - i * dt )

        #     index_when_plot = (int)( np.where(time_between_index_and_when_plot_in_hour == min(time_between_index_and_when_plot_in_hour))[0][0] )
        #     nb_hour_this_orbit = index_time_orbit_average_middle[index_when_plot] * dt/3600./24

        #     ratio_fig_size = 4./3
        #     fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        #     fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        #     plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        #     gs = gridspec.GridSpec(1, 1)
        #     gs.update(left= 0.10, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

        #     med = np.median(angular_distance_between_two_clusters[index_when_plot, :])
        #     mad = np.median( np.abs( angular_distance_between_two_clusters[index_when_plot, :] - np.median(angular_distance_between_two_clusters[index_when_plot, :]) ) )  
        #     quartile10 = np.percentile(angular_distance_between_two_clusters[index_when_plot, :], 10) 
        #     quartile25 = np.percentile(angular_distance_between_two_clusters[index_when_plot, :], 25) 
        #     quartile75 = np.percentile(angular_distance_between_two_clusters[index_when_plot, :], 75) 
        #     quartile90 = np.percentile(angular_distance_between_two_clusters[index_when_plot, :], 90) 
        #     iqr_every_orbit = quartile75 - quartile25
        #     quartiles_1090_every_orbit = quartile90 - quartile10
        #     std = np.std(angular_distance_between_two_clusters[index_when_plot, :])

        #     ax2 = fig.add_subplot(gs[0, 0])
        #     ax2.set_title('Distribution of the difference in period after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        #     ax2.set_xlabel('Difference in period (s)', fontsize = fontsize_plot, weight = 'bold')
        #     ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles**2)+ ')', fontsize = fontsize_plot, weight = 'bold')
        #     ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        #     [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
        #     bin_width = iqr_every_orbit/10.
        #     hist_data = angular_distance_between_two_clusters[index_when_plot, :]
        #     n, bins, patches = ax2.hist(hist_data , bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'r',label = '')  
        #     plt.vlines(quartile25, 0, max(n),linewidth = 2); plt.vlines(quartile75, 0,max(n), linewidth = 2, label = '25-75% quartiles') 
        #     plt.vlines(quartile10, 0, max(n), linewidth = 4, linestyle = 'dotted'); plt.vlines(quartile90, 0, max(n), linewidth = 4, linestyle = 'dotted', label = '10-90% deciles') 
        #     ax2.plot(med, 0, 'k',marker = '.', markersize = 20)
        #     ax2.get_xaxis().get_major_formatter().set_useOffset(False)
        #     ax2.set_xlim([med - 10*iqr_every_orbit, med+10*iqr_every_orbit])
        #     ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
        #     ax2.text(ax2.get_xlim()[0] + ( ax2.get_xlim()[1] - ax2.get_xlim()[0] )/ 60., ax2.get_ylim()[1] - (ax2.get_ylim()[1] - ax2.get_ylim()[0])/5., 'IQR = ' + '{0:.2e} s'.format(iqr_every_orbit), fontsize = fontsize_plot)

        #     ax2.legend(fontsize = fontsize_plot, loc = 2, scatterpoints=1)
        #     if save_results == 1:
        #         fig_save_name = 'phase_angle_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
        #         fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
        #         fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
        #         #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


            # Plot the distribution of the phase_angle every N hours ON THE SAME PLOT             # when_plot_in_hour = 2 * 24. # the few lines below convert this time in orbits                                                                                                                                                                                                                                                                                              
            ratio_fig_size = 4./3
            fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
            fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
            plt.rc('font', weight='bold') ## make the labels of the ticks in bold                                                                                                                                                                                                                                            
            # when_plot_in_hour = 2 * 24. # the few lines below convert this time in orbits 
            ratio_fig_size = 4./3
            fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
            fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
            plt.rc('font', weight='bold') ## make the labels of the ticks in bold
            gs = gridspec.GridSpec(1, 1)
            gs.update(left= 0.12, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

            time_step_when_plot_in_hour = 48. # in hour
            nb_steps_when_plot_in_hour = (int)( np.ceil(nb_steps * dt / (time_step_when_plot_in_hour * 3600. )) )
            ax2 = fig.add_subplot(gs[0, 0])
            ax2.set_title('Distribution of the angular distance every ' + str((int) (time_step_when_plot_in_hour)) + ' hours of propagation', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
            ax2.set_xlabel('Angular distance (' + u'\N{DEGREE SIGN}'+ ')', fontsize = fontsize_plot, weight = 'bold')
            ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles**2)+ ')', fontsize = fontsize_plot, weight = 'bold')
            ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
            [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure


            colors = ['k','r','m','b','c','y','g']
            
            for istep in range(1,nb_steps_when_plot_in_hour):
                when_plot_in_hour = time_step_when_plot_in_hour * istep
                time_between_index_and_when_plot_in_hour = np.zeros([len(save_index_dt_angular_distance_between_two_clusters)])
                icount = -1
                for i in save_index_dt_angular_distance_between_two_clusters:
                    icount = icount + 1
                    time_between_index_and_when_plot_in_hour[icount] = np.abs( when_plot_in_hour * 3600. - i * dt )

                index_when_plot = (int)( np.where(time_between_index_and_when_plot_in_hour == min(time_between_index_and_when_plot_in_hour))[0][0] )
                nb_hour_this_orbit = save_index_dt_angular_distance_between_two_clusters[index_when_plot] * dt/3600./24

                med = np.median(angular_distance_between_two_clusters[index_when_plot, :])
                mad = np.median( np.abs( angular_distance_between_two_clusters[index_when_plot, :] - np.median(angular_distance_between_two_clusters[index_when_plot, :]) ) )  
                quartile10 = np.percentile(angular_distance_between_two_clusters[index_when_plot, :], 10) 
                quartile25 = np.percentile(angular_distance_between_two_clusters[index_when_plot, :], 25) 
                quartile75 = np.percentile(angular_distance_between_two_clusters[index_when_plot, :], 75) 
                quartile90 = np.percentile(angular_distance_between_two_clusters[index_when_plot, :], 90) 
                iqr_every_orbit = quartile75 - quartile25
                quartiles_1090_every_orbit = quartile90 - quartile10
                std = np.std(angular_distance_between_two_clusters[index_when_plot, :])
                bin_width = iqr_every_orbit/10.
                if bin_width == 0:
                    bin_width = 0.000001
                hist_data = angular_distance_between_two_clusters[index_when_plot, :]
                n, bins, patches = ax2.hist(hist_data , bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7,label = '+ ' + str(when_plot_in_hour/24,) + ' d (IQR = ' + '{0:.1e}'.format(iqr_every_orbit)  + u'\N{DEGREE SIGN}'+ ')', color = colors[istep])  
                ax2.plot(med, 0, color = colors[istep],marker = '.', markersize = 20)
                ax2.get_xaxis().get_major_formatter().set_useOffset(False)
#                ax2.set_xlim([med - 10*iqr_every_orbit, med+10*iqr_every_orbit])
#                ax2.text(ax2.get_xlim()[0] + ( ax2.get_xlim()[1] - ax2.get_xlim()[0] )/ 60., ax2.get_ylim()[1] - (ax2.get_ylim()[1] - ax2.get_ylim()[0])/5., 'IQR = ' + '{0:.2e} s'.format(iqr_every_orbit), fontsize = fontsize_plot)
                
                print min(bins), max(bins)
            ax2.legend(fontsize = fontsize_plot, loc = 2, scatterpoints=1)
            ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
            if save_results == 1:
                fig_save_name = 'angle_cluster_distribution_many_steps_days_after_deployment'
                fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
                #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

    ############################## ORBIT AVERAGE PHASE_ANGLE
    ########################################################################################################################################################################################################################################################################################################

if "phase_angle_average" in sys.argv:
    if isat == 1: # stuff were saved for isat = 1 in mpi_distance_ensemble_to_main_sc.py (because to compute phase_angle you first need to compute the period of isat 0 and then isat 1)
        if ( 'phase_angle_ensemble_orbit_average' in pickle_list_name ):
            nb_steps_orbit_average = len(index_time_orbit_average[:,1])     # !!!!!!!! there might be ore phase_angles for satellites than for others (because they move at different speeds). So nb_steps_orbit_average should actually be the min of the number of phase_angles of all satellites. Since we don't take the min, there might be 0 in phase_angle_ensemble_orbit_average that we should NOT take into account. So basically NEVER TAKE THE STANDARD DEVIATION for this phase_angle section of the script (that takes into account outliers)
            index_time_orbit_average_middle = index_time_orbit_average[:,1]
            x_axis_orbit_average = index_time_orbit_average_middle * dt # conversion index to seconds
            x_axis_orbit_average = x_axis_orbit_average / (24*3600.) # conversion seconds to day

            # Plot the inter-quartile rang of the phase_angle distributions as a function of time
            step_std_in_index = 1 # we want the step to be every orbit (here we are look at orbit average data)
            std_every_orbit = np.zeros([nb_steps_orbit_average]) # DON'T TAKE IT FOR THE PHASE_ANGLE SECTION
            med_every_orbit = np.zeros([nb_steps_orbit_average])
            mad_every_orbit = np.zeros([nb_steps_orbit_average])
            iqr_every_orbit = np.zeros([nb_steps_orbit_average])
            for i in range(nb_steps_orbit_average):
                index_when_std = i * step_std_in_index
                med_every_orbit[i] = np.median(phase_angle_ensemble_orbit_average[index_when_std, :])
                std_every_orbit[i] = np.std(phase_angle_ensemble_orbit_average[index_when_std, :]) # NOT GOOD TO USE THE STD
                iqr_every_orbit[i] = np.subtract(*np.percentile(phase_angle_ensemble_orbit_average[index_when_std, :], [75, 25])) 
                mad_every_orbit[i] = np.median( np.abs( phase_angle_ensemble_orbit_average[index_when_std, :] - np.median(phase_angle_ensemble_orbit_average[index_when_std, :]) ) ) 


            ## Plot
            ratio_fig_size = 4./3
            fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
            fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
            plt.rc('font', weight='bold') ## make the labels of the ticks in bold
            gs = gridspec.GridSpec(1, 1)
            gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
            ax1 = fig.add_subplot(gs[0, 0])
            ax1.set_title('Inter-quartile range of the difference in periods VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
            ax1.set_ylabel('IQR of difference in period (s)', fontsize = fontsize_plot, weight = 'bold')
            ax1.set_xlabel('Time (daysx)', fontsize = fontsize_plot, weight = 'bold')
            ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
            [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
            ax1.plot(x_axis_orbit_average, iqr_every_orbit, 'k', linewidth = 2)

            ax1.get_xaxis().tick_bottom()
            ax1.get_yaxis().tick_left()
            ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)



            # nb_days_simu = ( date_stop - date_start ).total_seconds()/3600./24.
            # nb_hours_simu = nb_steps * dt/ 3600.
            # hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / 24. # x_axis_orbit_average is in unit of days
            # xticks = np.arange(0, nb_days_simu+1, hour_time_step_xticks_converted_in_index_adjusted)
            # date_list_str = []
            # date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
            # for i in range(len(xticks)):
            #     date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] )) + ")")
            # ax1.xaxis.set_ticks(xticks)
            # ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')


            if save_results == 1:
                fig_save_name = 'phase_angle_iqr_every_orbit'
                fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
                #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


        if ( 'phase_angle_ensemble_orbit_average' in pickle_list_name ):
            nb_steps_orbit_average = len(index_time_orbit_average[:,1])     # !!!!!!!! there might be ore phase_angles for satellites than for others (because they move at different speeds). So nb_steps_orbit_average should actually be the min of the number of phase_angles of all satellites. Since we don't take the min, there might be 0 in phase_angle_ensemble_orbit_average that we should NOT take into account. So basically NEVER TAKE THE STANDARD DEVIATION for this phase_angle section of the script (that takes into account outliers)
            index_time_orbit_average_middle = index_time_orbit_average[:,1]
            x_axis_orbit_average = index_time_orbit_average_middle * dt # conversion index to seconds
            x_axis_orbit_average = x_axis_orbit_average / (24*3600.) # conversion seconds to day

            # Plot the inter-quartile rang of the phase_angle distributions as a function of time
            step_std_in_index = 1 # we want the step to be every orbit (here we are look at orbit average data)
            std_every_orbit = np.zeros([nb_steps_orbit_average]) # DON'T TAKE IT FOR THE PHASE_ANGLE SECTION
            med_every_orbit = np.zeros([nb_steps_orbit_average])
            mad_every_orbit = np.zeros([nb_steps_orbit_average])
            iqr_every_orbit = np.zeros([nb_steps_orbit_average])
            for i in range(nb_steps_orbit_average):
                index_when_std = i * step_std_in_index
                med_every_orbit[i] = np.median(phase_angle_ensemble_orbit_average[index_when_std, :])
                std_every_orbit[i] = np.std(phase_angle_ensemble_orbit_average[index_when_std, :]) # NOT GOOD TO USE THE STD
                iqr_every_orbit[i] = np.subtract(*np.percentile(phase_angle_ensemble_orbit_average[index_when_std, :], [75, 25])) 
                mad_every_orbit[i] = np.median( np.abs( phase_angle_ensemble_orbit_average[index_when_std, :] - np.median(phase_angle_ensemble_orbit_average[index_when_std, :]) ) ) 


            ## Plot
            ratio_fig_size = 4./3
            fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
            fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
            plt.rc('font', weight='bold') ## make the labels of the ticks in bold
            gs = gridspec.GridSpec(1, 1)
            gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
            ax1 = fig.add_subplot(gs[0, 0])
            ax1.set_title('Median of the difference in periods VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
            ax1.set_ylabel('Median of difference in period (s)', fontsize = fontsize_plot, weight = 'bold')
            ax1.set_xlabel('Time (days)', fontsize = fontsize_plot, weight = 'bold')
            ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
            [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
            ax1.plot(x_axis_orbit_average, med_every_orbit, 'k', linewidth = 2)

            ax1.get_xaxis().tick_bottom()
            ax1.get_yaxis().tick_left()
            ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)



            # nb_days_simu = ( date_stop - date_start ).total_seconds()/3600./24.
            # nb_hours_simu = nb_steps * dt/ 3600.
            # hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / 24. # x_axis_orbit_average is in unit of days
            # xticks = np.arange(0, nb_days_simu+1, hour_time_step_xticks_converted_in_index_adjusted)
            # date_list_str = []
            # date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
            # for i in range(len(xticks)):
            #     date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] )) + ")")
            # ax1.xaxis.set_ticks(xticks)
            # ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')


            if save_results == 1:
                fig_save_name = 'phase_angle_median_every_orbit'
                fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
                #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

        ########################################################################################################################################################################################################################################################################################################
        if ( 'phase_angle_ensemble_orbit_average' in pickle_list_name ):
            # Plot the distribution of the phase_angle at different times (set by when_plot_in_hour)
            # when_plot_in_hour = 2 * 24. # the few lines below convert this time in orbits 
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
            gs.update(left= 0.10, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

            med = np.median(phase_angle_ensemble_orbit_average[index_when_plot, :])
            mad = np.median( np.abs( phase_angle_ensemble_orbit_average[index_when_plot, :] - np.median(phase_angle_ensemble_orbit_average[index_when_plot, :]) ) )  
            quartile10 = np.percentile(phase_angle_ensemble_orbit_average[index_when_plot, :], 10) 
            quartile25 = np.percentile(phase_angle_ensemble_orbit_average[index_when_plot, :], 25) 
            quartile75 = np.percentile(phase_angle_ensemble_orbit_average[index_when_plot, :], 75) 
            quartile90 = np.percentile(phase_angle_ensemble_orbit_average[index_when_plot, :], 90) 
            iqr_every_orbit = quartile75 - quartile25
            quartiles_1090_every_orbit = quartile90 - quartile10
            std = np.std(phase_angle_ensemble_orbit_average[index_when_plot, :])

            ax2 = fig.add_subplot(gs[0, 0])
            ax2.set_title('Distribution of the difference in period after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
            ax2.set_xlabel('Difference in period (s)', fontsize = fontsize_plot, weight = 'bold')
            ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles**2)+ ')', fontsize = fontsize_plot, weight = 'bold')
            ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
            [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
            bin_width = iqr_every_orbit/10.
            hist_data = phase_angle_ensemble_orbit_average[index_when_plot, :]
            n, bins, patches = ax2.hist(hist_data , bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'r',label = '')  
            plt.vlines(quartile25, 0, max(n),linewidth = 2); plt.vlines(quartile75, 0,max(n), linewidth = 2, label = '25-75% quartiles') 
            plt.vlines(quartile10, 0, max(n), linewidth = 4, linestyle = 'dotted'); plt.vlines(quartile90, 0, max(n), linewidth = 4, linestyle = 'dotted', label = '10-90% deciles') 
            ax2.plot(med, 0, 'k',marker = '.', markersize = 20)
            ax2.get_xaxis().get_major_formatter().set_useOffset(False)
            ax2.set_xlim([med - 10*iqr_every_orbit, med+10*iqr_every_orbit])
            ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
            ax2.text(ax2.get_xlim()[0] + ( ax2.get_xlim()[1] - ax2.get_xlim()[0] )/ 60., ax2.get_ylim()[1] - (ax2.get_ylim()[1] - ax2.get_ylim()[0])/5., 'IQR = ' + '{0:.2e} s'.format(iqr_every_orbit), fontsize = fontsize_plot)

            ax2.legend(fontsize = fontsize_plot, loc = 2, scatterpoints=1)
            if save_results == 1:
                fig_save_name = 'phase_angle_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
                fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
                #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


            # Plot the distribution of the phase_angle every N hours ON THE SAME PLOT             # when_plot_in_hour = 2 * 24. # the few lines below convert this time in orbits                                                                                                                                                                                                                                                                                              
            ratio_fig_size = 4./3
            fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
            fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
            plt.rc('font', weight='bold') ## make the labels of the ticks in bold                                                                                                                                                                                                                                            
            # when_plot_in_hour = 2 * 24. # the few lines below convert this time in orbits 
            ratio_fig_size = 4./3
            fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
            fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
            plt.rc('font', weight='bold') ## make the labels of the ticks in bold
            gs = gridspec.GridSpec(1, 1)
            gs.update(left= 0.12, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

            time_step_when_plot_in_hour = 48. # in hour
            nb_steps_when_plot_in_hour = (int)( np.ceil(nb_steps * dt / (time_step_when_plot_in_hour * 3600. )) )
            ax2 = fig.add_subplot(gs[0, 0])
            ax2.set_title('Distribution of the difference in period every ' + str((int) (time_step_when_plot_in_hour)) + ' hours of propagation', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
            ax2.set_xlabel('Difference in period (s)', fontsize = fontsize_plot, weight = 'bold')
            ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles**2)+ ')', fontsize = fontsize_plot, weight = 'bold')
            ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
            [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
            bin_width = iqr_every_orbit/10.

            colors = ['k','r','m','b','c','y','g']
            
            for istep in range(nb_steps_when_plot_in_hour):
                when_plot_in_hour = time_step_when_plot_in_hour * istep
                time_between_index_and_when_plot_in_hour = np.zeros([len(index_time_orbit_average_middle)])
                icount = -1
                for i in index_time_orbit_average_middle:
                    icount = icount + 1
                    time_between_index_and_when_plot_in_hour[icount] = np.abs( when_plot_in_hour * 3600. - i * dt )

                index_when_plot = (int)( np.where(time_between_index_and_when_plot_in_hour == min(time_between_index_and_when_plot_in_hour))[0][0] )
                nb_hour_this_orbit = index_time_orbit_average_middle[index_when_plot] * dt/3600./24

                med = np.median(phase_angle_ensemble_orbit_average[index_when_plot, :])
                mad = np.median( np.abs( phase_angle_ensemble_orbit_average[index_when_plot, :] - np.median(phase_angle_ensemble_orbit_average[index_when_plot, :]) ) )  
                quartile10 = np.percentile(phase_angle_ensemble_orbit_average[index_when_plot, :], 10) 
                quartile25 = np.percentile(phase_angle_ensemble_orbit_average[index_when_plot, :], 25) 
                quartile75 = np.percentile(phase_angle_ensemble_orbit_average[index_when_plot, :], 75) 
                quartile90 = np.percentile(phase_angle_ensemble_orbit_average[index_when_plot, :], 90) 
                iqr_every_orbit = quartile75 - quartile25
                quartiles_1090_every_orbit = quartile90 - quartile10
                std = np.std(phase_angle_ensemble_orbit_average[index_when_plot, :])

                hist_data = phase_angle_ensemble_orbit_average[index_when_plot, :]
                n, bins, patches = ax2.hist(hist_data , bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7,label = '+ ' + str(when_plot_in_hour/24,) + ' d (IQR = ' + '{0:.1e} s'.format(iqr_every_orbit) +')', color = colors[istep])  
                ax2.plot(med, 0, color = colors[istep],marker = '.', markersize = 20)
                ax2.get_xaxis().get_major_formatter().set_useOffset(False)
#                ax2.set_xlim([med - 10*iqr_every_orbit, med+10*iqr_every_orbit])
#                ax2.text(ax2.get_xlim()[0] + ( ax2.get_xlim()[1] - ax2.get_xlim()[0] )/ 60., ax2.get_ylim()[1] - (ax2.get_ylim()[1] - ax2.get_ylim()[0])/5., 'IQR = ' + '{0:.2e} s'.format(iqr_every_orbit), fontsize = fontsize_plot)
                
                print min(bins), max(bins)
            ax2.legend(fontsize = fontsize_plot, loc = 2, scatterpoints=1)
            ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
            if save_results == 1:
                fig_save_name = 'phase_angle_distribution_many_steps_days_after_deployment'
                fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
                #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


        # if ( 'phase_angle_ensemble_orbit_average' in pickle_list_name ):
        #     # Plot the phase_angle of each ensemble as a function of time, as well as the median and the inter-quartile range
        #     step_std_in_index = 1 # every orbit
        #     med_every_orbit = np.zeros([nb_steps_orbit_average])
        #     mad_every_orbit = np.zeros([nb_steps_orbit_average])
        #     iqr_every_orbit = np.zeros([nb_steps_orbit_average])
        #     quartile25 =np.zeros([nb_steps_orbit_average])
        #     quartile75 =np.zeros([nb_steps_orbit_average])
        #     for i in range(nb_steps_orbit_average):
        #         index_when_std = i * step_std_in_index
        #         iqr_every_orbit[i] = np.subtract(*np.percentile(phase_angle_ensemble_orbit_average[index_when_std, :], [75, 25])) 
        #         mad_every_orbit[i] = np.median( np.abs( phase_angle_ensemble_orbit_average[index_when_std, :] - np.median(phase_angle_ensemble_orbit_average[index_when_std, :]) ) ) 
        #         med_every_orbit[i] = np.median(phase_angle_ensemble_orbit_average[index_when_std, :])
        #         quartile25[i] = np.percentile(phase_angle_ensemble_orbit_average[index_when_std, :], 25) 
        #         quartile75[i] = np.percentile(phase_angle_ensemble_orbit_average[index_when_std, :], 75) 

        #     ## Plot
        #     ratio_fig_size = 4./3
        #     fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        #     fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        #     plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        #     gs = gridspec.GridSpec(1, 1)
        #     gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
        #     ax1 = fig.add_subplot(gs[0, 0])
        #     ax1.set_title('Phase_Angle at the position of each ensemble member VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        #     ax1.set_ylabel('Phase_Angle (s)', fontsize = fontsize_plot, weight = 'bold')
        #     ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
        #     ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        #     [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        #     for iens in range(nb_ensembles):
        #         if iens == 0:
        #             ax1.plot(x_axis_orbit_average, phase_angle_ensemble_orbit_average[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 2, alpha = 1, label = 'Ensemble')
        #         else:
        #             ax1.plot(x_axis_orbit_average, phase_angle_ensemble_orbit_average[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 2, alpha = 1)
        #     ax1.plot(x_axis_orbit_average, med_every_orbit, 'k', linewidth = 3, label = 'Median')
        #     ax1.plot(x_axis_orbit_average, quartile25, 'r', linewidth = 3, label = '25 and 75% quartiles')
        #     ax1.plot(x_axis_orbit_average, quartile75, 'r', linewidth = 3)

        #     ax1.legend(fontsize = fontsize_plot, loc = 2)
        #     ax1.get_xaxis().tick_bottom()
        #     ax1.get_yaxis().tick_left()
        #     ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)



        #     nb_days_simu = ( date_stop - date_start ).total_seconds()/3600./24.
        #     nb_hours_simu = nb_steps * dt/ 3600.
        #     hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / 24. # x_axis_orbit_average is in unit of days
        #     xticks = np.arange(0, nb_days_simu, hour_time_step_xticks_converted_in_index_adjusted)
        #     date_list_str = []
        #     date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
        #     for i in range(len(xticks)):
        #         date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] )) + ")")
        #     ax1.xaxis.set_ticks(xticks)
        #     ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')



        #     if save_results == 1:
        #         fig_save_name = 'phase_angle_all_ensemble'
        #         fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
        #         fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
        #         #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


    else:
        print "You need to set isat = 1 at top of script to plot phase_angle"







    ############################## ORBIT AVERAGE SMA_DIFFERENCE
    ########################################################################################################################################################################################################################################################################################################

if "sma_difference_average" in sys.argv:
    if isat == 1: # stuff were saved for isat = 1 in mpi_distance_ensemble_to_main_sc.py (because to compute sma_difference you first need to compute the period of isat 0 and then isat 1)
        if ( 'sma_difference_ensemble_orbit_average' in pickle_list_name ):
            nb_steps_orbit_average = len(index_time_orbit_average[:,1])     # !!!!!!!! there might be ore sma_differences for satellites than for others (because they move at different speeds). So nb_steps_orbit_average should actually be the min of the number of sma_differences of all satellites. Since we don't take the min, there might be 0 in sma_difference_ensemble_orbit_average that we should NOT take into account. So basically NEVER TAKE THE STANDARD DEVIATION for this sma_difference section of the script (that takes into account outliers)
            index_time_orbit_average_middle = index_time_orbit_average[:,1]
            x_axis_orbit_average = index_time_orbit_average_middle * dt # conversion index to seconds
            x_axis_orbit_average = x_axis_orbit_average / (24*3600.) # conversion seconds to day

            # Plot the inter-quartile rang of the sma_difference distributions as a function of time
            step_std_in_index = 1 # we want the step to be every orbit (here we are look at orbit average data)
            std_every_orbit = np.zeros([nb_steps_orbit_average]) # DON'T TAKE IT FOR THE SMA_DIFFERENCE SECTION
            med_every_orbit = np.zeros([nb_steps_orbit_average])
            mad_every_orbit = np.zeros([nb_steps_orbit_average])
            iqr_every_orbit = np.zeros([nb_steps_orbit_average])
            for i in range(nb_steps_orbit_average):
                index_when_std = i * step_std_in_index
                std_every_orbit[i] = np.std(sma_difference_ensemble_orbit_average[index_when_std, :]) # NOT GOOD TO USE THE STD
                iqr_every_orbit[i] = np.subtract(*np.percentile(sma_difference_ensemble_orbit_average[index_when_std, :], [75, 25])) 
                mad_every_orbit[i] = np.median( np.abs( sma_difference_ensemble_orbit_average[index_when_std, :] - np.median(sma_difference_ensemble_orbit_average[index_when_std, :]) ) ) 


            ## Plot
            ratio_fig_size = 4./3
            fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
            fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
            plt.rc('font', weight='bold') ## make the labels of the ticks in bold
            gs = gridspec.GridSpec(1, 1)
            gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
            ax1 = fig.add_subplot(gs[0, 0])
            ax1.set_title('Inter-quartile range of the sma_difference distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
            ax1.set_ylabel('IQR (#/m^3)', fontsize = fontsize_plot, weight = 'bold')
            ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
            ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
            [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
            ax1.plot(x_axis_orbit_average, iqr_every_orbit, 'k', linewidth = 2)

            ax1.get_xaxis().tick_bottom()
            ax1.get_yaxis().tick_left()
            ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)



            # nb_days_simu = ( date_stop - date_start ).total_seconds()/3600./24.
            # nb_hours_simu = nb_steps * dt/ 3600.
            # hour_time_step_xticks_converted_in_index_adjusted = hour_time_step_xticks / 24. # x_axis_orbit_average is in unit of days
            # xticks = np.arange(0, nb_days_simu+1, hour_time_step_xticks_converted_in_index_adjusted)
            # date_list_str = []
            # date_list = [date_start + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
            # for i in range(len(xticks)):
            #     date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] )) + ")")
            # ax1.xaxis.set_ticks(xticks)
            # ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')


            if save_results == 1:
                fig_save_name = 'sma_difference_iqr_every_orbit'
                fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
                #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

        ########################################################################################################################################################################################################################################################################################################
        if ( 'sma_difference_ensemble_orbit_average' in pickle_list_name ):
            # Plot the distribution of the sma_difference at different times (set by when_plot_in_hour)
            # when_plot_in_hour = 2 * 24. # the few lines below convert this time in orbits 
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

            med = np.median(sma_difference_ensemble_orbit_average[index_when_plot, :])
            mad = np.median( np.abs( sma_difference_ensemble_orbit_average[index_when_plot, :] - np.median(sma_difference_ensemble_orbit_average[index_when_plot, :]) ) )  
            quartile10 = np.percentile(sma_difference_ensemble_orbit_average[index_when_plot, :], 10) 
            quartile25 = np.percentile(sma_difference_ensemble_orbit_average[index_when_plot, :], 25) 
            quartile75 = np.percentile(sma_difference_ensemble_orbit_average[index_when_plot, :], 75) 
            quartile90 = np.percentile(sma_difference_ensemble_orbit_average[index_when_plot, :], 90) 
            iqr_every_orbit = quartile75 - quartile25
            quartiles_1090_every_orbit = quartile90 - quartile10
            std = np.std(sma_difference_ensemble_orbit_average[index_when_plot, :])

            ax2 = fig.add_subplot(gs[0, 0])
            ax2.set_title('Sma_Difference distribution ' + 'after ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days of propagation', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
            ax2.set_xlabel('Sma_Difference (#/m^3)', fontsize = fontsize_plot, weight = 'bold')
            ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
            ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
            [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
            bin_width = iqr_every_orbit/10.
            hist_data = sma_difference_ensemble_orbit_average[index_when_plot, :]
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
                fig_save_name = 'sma_difference_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
                fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
                #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


        if ( 'sma_difference_ensemble_orbit_average' in pickle_list_name ):
            # Plot the sma_difference of each ensemble as a function of time, as well as the median and the inter-quartile range
            step_std_in_index = 1 # every orbit
            med_every_orbit = np.zeros([nb_steps_orbit_average])
            mad_every_orbit = np.zeros([nb_steps_orbit_average])
            iqr_every_orbit = np.zeros([nb_steps_orbit_average])
            quartile25 =np.zeros([nb_steps_orbit_average])
            quartile75 =np.zeros([nb_steps_orbit_average])
            for i in range(nb_steps_orbit_average):
                index_when_std = i * step_std_in_index
                iqr_every_orbit[i] = np.subtract(*np.percentile(sma_difference_ensemble_orbit_average[index_when_std, :], [75, 25])) 
                mad_every_orbit[i] = np.median( np.abs( sma_difference_ensemble_orbit_average[index_when_std, :] - np.median(sma_difference_ensemble_orbit_average[index_when_std, :]) ) ) 
                med_every_orbit[i] = np.median(sma_difference_ensemble_orbit_average[index_when_std, :])
                quartile25[i] = np.percentile(sma_difference_ensemble_orbit_average[index_when_std, :], 25) 
                quartile75[i] = np.percentile(sma_difference_ensemble_orbit_average[index_when_std, :], 75) 

            ## Plot
            ratio_fig_size = 4./3
            fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
            fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
            plt.rc('font', weight='bold') ## make the labels of the ticks in bold
            gs = gridspec.GridSpec(1, 1)
            gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
            ax1 = fig.add_subplot(gs[0, 0])
            ax1.set_title('Sma_Difference at the position of each ensemble member VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
            ax1.set_ylabel('Sma_Difference (s)', fontsize = fontsize_plot, weight = 'bold')
            ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
            ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
            [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
            for iens in range(nb_ensembles):
                if iens == 0:
                    ax1.plot(x_axis_orbit_average, sma_difference_ensemble_orbit_average[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 2, alpha = 1, label = 'Ensemble')
                else:
                    ax1.plot(x_axis_orbit_average, sma_difference_ensemble_orbit_average[0:nb_steps:step_std_in_index,iens], color = 'b', linewidth = 2, alpha = 1 )
            ax1.plot(x_axis_orbit_average, med_every_orbit, 'k', linewidth = 3, label = 'Median')
            ax1.plot(x_axis_orbit_average, quartile25, 'r', linewidth = 3, label = '25 and 75% quartiles')
            ax1.plot(x_axis_orbit_average, quartile75, 'r', linewidth = 3)

            ax1.legend(fontsize = fontsize_plot, loc = 2)
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



            if save_results == 1:
                fig_save_name = 'sma_difference_all_ensemble'
                fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
                #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

    else:
        print "You need to set isat = 1 at top of script to plot sma_difference"


if "tca" in sys.argv:
    if ( 'tca_ensemble' in pickle_list_name ):
        # Plot the distribution of the f10.7 at different times (set by when_plot_in_hour)
        # when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
        

        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.2, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)
        nb_collisions_recorded_in_tca = len(tca_ensemble)
        ax2 = fig.add_subplot(gs[0, 0])

        ax2.set_ylabel('Percentage in bin (%)', fontsize = fontsize_plot, weight = 'bold')
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure

        # tca_ensemble is the number of seconds between the initial epoch and the tca of each ensemble. Below we substract to this nb of seconds the number of seconds between the initial epoch and the tca of the unperturbed orbits
        collision_filename = input_variables[find_in_read_input_order_variables(order_input_variables, 'collision_filename')]
        date_collision, nb_collisions_each_dt, cpc, cpc_final, tca = read_collision_file( collision_filename )
        date_ref = datetime.datetime.strptime(tca[0], "%Y-%m-%dT%H:%M:%S.%f") # !!!!!!!!! assumption: if more than one close unperturbed approach then take the first one
        ax2.set_xlabel('Seconds since ' + str(date_ref)[:-3], fontsize = fontsize_plot, weight = 'bold')
        nb_seconds_from_epoch_to_date_ref = ( date_ref - date_start ).total_seconds()
        hist_data = tca_ensemble  - nb_seconds_from_epoch_to_date_ref
        bin_width = np.std(hist_data)/10
#        ax2.set_title('TCA distribution (bin size: '+ '{0:.2f}'.format(bin_width)+' s)', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax2.set_title('', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        eighty_width = np.subtract(*np.percentile(hist_data, [90, 10]))
        fifty_width = np.subtract(*np.percentile(hist_data, [75, 25]))

        weights = np.ones_like(hist_data)/len(hist_data)
        n, bins, patches = ax2.hist(hist_data,bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width) ,  histtype='stepfilled', alpha = 0.7, color = 'r',label = '', weights = weights * 100)  

#        ax2.text(ax2.get_xlim()[0] + ( ax2.get_xlim()[1] - ax2.get_xlim()[0] )/ 60., ax2.get_ylim()[1] - 0.5*(ax2.get_ylim()[1] - ax2.get_ylim()[0])/10., '80% width = ' + '{0:.2e}'.format(eighty_width) + ' s', fontsize = fontsize_plot)
        ax2.text(ax2.get_xlim()[0] + ( ax2.get_xlim()[1] - ax2.get_xlim()[0] )/ 60., ax2.get_ylim()[1] - 0.5*(ax2.get_ylim()[1] - ax2.get_ylim()[0])/10., '50% width = ' + '{0:.2e}'.format(fifty_width) + ' s', fontsize = fontsize_plot)
        if save_results == 1:
            fig_save_name = 'tca_distribution_'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)



if "dca" in sys.argv:
    if ( 'dca_ensemble' in pickle_list_name ):
        # Plot the distribution of the f10.7 at different times (set by when_plot_in_hour)
        # when_plot_in_hour = 2 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
        

        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.2, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)
        nb_collisions_recorded_in_dca = len(dca_ensemble)
        ax2 = fig.add_subplot(gs[0, 0])

        ax2.set_ylabel('Percentage in bin (%)', fontsize = fontsize_plot, weight = 'bold')
        ax2.set_xlabel('Seconds since epoch', fontsize = fontsize_plot, weight = 'bold')
        ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
        hist_data = dca_ensemble 
        bin_width = np.std(hist_data)/10
        ax2.set_title('DCA distribution (bin size: '+ '{0:.2f}'.format(bin_width)+' s)', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        eighty_width = np.subtract(*np.percentile(hist_data, [90, 10]))
        fifty_width = np.subtract(*np.percentile(hist_data, [75, 25]))

        weights = np.ones_like(hist_data)/len(hist_data)
        n, bins, patches = ax2.hist(hist_data,bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width) ,  histtype='stepfilled', alpha = 0.7, color = 'r',label = '', weights = weights * 100)  

        ax2.text(ax2.get_xlim()[0] + ( ax2.get_xlim()[1] - ax2.get_xlim()[0] )/ 60., ax2.get_ylim()[1] - (ax2.get_ylim()[1] - ax2.get_ylim()[0])/10., '80% width = ' + '{0:.2e}'.format(eighty_width) + ' km', fontsize = fontsize_plot)
        ax2.text(ax2.get_xlim()[0] + ( ax2.get_xlim()[1] - ax2.get_xlim()[0] )/ 60., ax2.get_ylim()[1] - 1.5*(ax2.get_ylim()[1] - ax2.get_ylim()[0])/10., '50% width = ' + '{0:.2e}'.format(fifty_width) + ' km', fontsize = fontsize_plot)
        if save_results == 1:
            fig_save_name = 'dca_distribution_'
            fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)


# 
# iens1 = 100
# iens2 = 200
# dist = np.sqrt( ( x_eci_ensemble[:, iens2] - x_eci_ensemble[:, iens1] )**2 + ( y_eci_ensemble[:, iens2] - y_eci_ensemble[:, iens1] )**2 + ( z_eci_ensemble[:, iens2] - z_eci_ensemble[:, iens1] )**2 ) 
        
if show_plots == 1:
    plt.show(); plt.show()
