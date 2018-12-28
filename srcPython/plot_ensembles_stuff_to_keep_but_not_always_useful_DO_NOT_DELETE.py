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

plt.ion()
plt.isinteractive()

## NOTE 0: to run this script, you first need to run distance_ensemble_to_main_sc.py (with first_time = 1). This will create the data that plot_ensembles will then use to makes plots
## NOTE 1: to use this script, the only 2 parameters you have to set are:
## - the name of the propagator main input file 
## - the path of the folder where you want to store the results (pickle, image, video): called 'path_folder_results'. In this folder, there must be the following subfolders: 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', and 'other'. In each of these subfolders, there must be the 2 subsubfolders: 'result', and 'pickle'. In the subsubfolder 'result', there must be the 2 subsubsubfolders: 'image', and 'video'.
## NOTE 2: this can be run if only ONE MAIN satellite was run (with ensembles of course)
## NOTE 3: the results (pickle, image, video) are stored in a subfolder of path_folder_results. The subfolder is either 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', or 'other'. This code figures out which folder the simulation corresponds to by reading the name of the output file chosen by the user in the main input file of the propagator (third line of section #SPACECRAFTS): it tries to find 'cygnss', 'cadre', 'aerie', 'scion', or 'qb50' in the name of the output file. If it does not find it, then the results here will be stored in path_folder_results/other/
## NOTE 4: to run this script, and any python script in the propagator, you need to be one subfolder deep from the main folder where the propagator runs are made. So if path_to_propagator/PropSim is the folder where the propagator runs are made, then the python scripts must be run for example in path_to_propagator/PropSim/subfolder_where_python_scipts_are_run


# !!!!!!!!!! SET THE PARAMETER BELOW: path of the folder where you want to store the results (pickle, image, video)
path_folder_results = '/raid3/Armada/Charles/python/' #get_prop_dir(2) + 'output/python_propagator/'

if len(sys.argv) > 2:
    main_input_file_name = get_prop_dir(1) + sys.argv[1] + '/input/main_input/' + sys.argv[2]  # ex of sys.argv[1]: cadre_tumbling_0panel.txt'
else:
    main_input_file_name = get_prop_dir(1) + 'run/input/main_input/' + sys.argv[1]  # ex of sys.argv[1]: cadre_tumbling_0panel.txt'    


# read input file
input_variables, order_input_variables = read_input_file(main_input_file_name)
dt = input_variables[2]; nb_steps = input_variables[3]; satellite_to_plot_path = input_variables[6][0]; satellite_to_plot = input_variables[7][0]; nb_ensembles_coe = input_variables[12]; nb_ensembles_attitude = input_variables[13]; nb_ensembles_cd = input_variables[11]; ensemble_to_plot_temp = input_variables[14]; 
nb_ensembles_density = input_variables[17]
n = nb_steps


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


    
save_pickle_name = path_folder_results + name_mission + '/pickle/' + satellite_to_plot.replace(".txt",".pickle")
root_save_fig_name = path_folder_results + name_mission + '/result/image/' + satellite_to_plot.replace(".txt","_") 
root_save_video_name = path_folder_results + name_mission + '/result/video/' + satellite_to_plot.replace(".txt","_")

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

raise Exception
########################################################################################################################################################################################################################################################################################################
# PLOTS ################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################

    
########################################################################################################################################################################################################################################################################################################
if ( 'angle_asc_node_to_sat_ensemble' in pickle_list_name ):
    # Plot the standard deviation of the angular_spacing_between_ensemble_sat_converted_in_a_distance distributions as a function of time
    step_std = 24 # step in hours to calculate the standard deviation
    step_std_in_index = step_std * 3600. / dt
    std_daily = np.zeros([nb_steps])
    med_daily = np.zeros([nb_steps])
    mad_daily = np.zeros([nb_steps])
    iqr_daily = np.zeros([nb_steps])
    for i in range(nb_steps):
        index_when_std = i * step_std_in_index
        std_daily[i] = np.std(angle_asc_node_to_sat_ensemble[index_when_std, :]) # NOT GOOD TO USE THE STD
        med_daily[i] = np.median(angular_spacing_between_ensemble_sat_temp[index_when_std,:]) # in meters
        iqr_daily[i] = np.subtract(*np.percentile(angle_asc_node_to_sat_ensemble[index_when_std, :], [75, 25])) 
        mad_daily[i] = np.median( np.abs( angle_asc_node_to_sat_ensemble[index_when_std, :] - np.median(angle_asc_node_to_sat_ensemble[index_when_std, :]) ) ) 

    ## Plot
    fontsize_plot = 14
    height_fig = 5#12 
    ratio_fig_size = 4./3
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left= 0.09, right=0.98, top = 0.93,bottom = 0.08)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Standard deviation of the along-track distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
    ax1.set_ylabel('Standard Deviation (m)', fontsize = fontsize_plot, weight = 'bold')
    ax1.set_xlabel('Time (days)', fontsize = fontsize_plot, weight = 'bold')
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.plot(np.arange(0,nb_steps), iqr_daily[0:nb_steps]*110000, 'k', linewidth = 2)

    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

    fig_save_name = 'std_daily'
    fig_save_name = root_save_fig_name + fig_save_name + '.png'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


########################################################################################################################################################################################################################################################################################################
if ( 'angle_asc_node_to_sat_ensemble' in pickle_list_name ):
    # Plot the distribution of angular_spacing_between_ensemble_sat_converted_in_a_distance at different times (set by when_plot_in_hour)
    when_plot_in_hour = 7 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
    index_when_plot = (int) (when_plot_in_hour * 3600L / dt)

    fontsize_plot = 14
    height_fig = 12
    ratio_fig_size = 4./3
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('Spacecraft distribution along the orbit ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days after deployment', y = 0.96,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(4, 1)
    gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)
    bin_width = 10 # in m

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
    ax1.text( med + 10 * iqr_daily, 0 - length_vert_iqr_daily/1.5,'Orbit', horizontalalignment = 'right', fontsize = 14, color = 'b' )

    ax1.plot(angle_asc_node_to_sat_ensemble[index_when_plot, :], np.zeros([nb_ensembles]), 'r.', markersize = 10)
    ax1.vlines(quartile25, 0,length_vert_iqr_daily, linewidth = 2); ax1.vlines(quartile75, 0,length_vert_iqr_daily, linewidth = 2) 
    ax1.text((quartile75 + quartile25 ) /2., length_vert_iqr_daily+length_vert_iqr_daily/5, '{0:.2f}'.format(iqr_daily * 110) + ' km (50% of sc)', horizontalalignment = 'center', fontsize = 14 )
    ax1.annotate ('', (quartile25, length_vert_iqr_daily), (quartile75, length_vert_iqr_daily), arrowprops={'arrowstyle':'<->'})
    ax1.vlines(quartile10, -length_vert_iqr_daily ,0 , linewidth = 4, linestyle = 'dotted'); ax1.vlines(quartile90, -length_vert_iqr_daily ,0, linewidth = 4, linestyle = 'dotted') 
    ax1.annotate ('', (quartile10, -length_vert_iqr_daily), (quartile90, -length_vert_iqr_daily), arrowprops={'arrowstyle':'<->'})
    ax1.text((quartile10 + quartile90 ) /2., -length_vert_iqr_daily-length_vert_iqr_daily/1.5, '{0:.2f}'.format(quartiles_1090_daily * 110) + ' km (80% of sc)', horizontalalignment = 'center', fontsize = 14, verticalalignment = 'bottom' )
    ax1.plot(med, 0, 'k',marker = '.', markersize = 15)

    plt.legend(loc = 2, scatterpoints=1)
    ax1.get_xaxis().tick_bottom()
    ax1.set_ylim([ymin, ymax])
#    ax1.get_yaxis().tick_left()
    ax1.yaxis.set_visible(False)
    ax1.spines['left'].set_visible(False); ax1.spines['right'].set_visible(False); ax1.spines['top'].set_visible(False)
    ax1.get_xaxis().get_major_formatter().set_useOffset(False)
    ax1.set_xlim([med - 10*iqr_daily, med+10*iqr_daily])
#    ax1.set_xlim([min(bins), 0])
    ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

    fig_save_name = 'sc_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
    fig_save_name = root_save_fig_name + fig_save_name + '.eps'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission)




raise Exception

########################################################################################################################################################################################################################################################################################################
# Distance between 2 satellites 
sat1 = np.where(angle_asc_node_to_sat_ensemble[90,:] == min(angle_asc_node_to_sat_ensemble[90,:]))[0]
sat2 = np.where(angle_asc_node_to_sat_ensemble[90,:] == max(angle_asc_node_to_sat_ensemble[90,:]))[0]

yaw_vel_sat1 = ( yaw_ensemble[90, sat1] - yaw_ensemble[90-1, sat1] ) / (24. * 3600)
pitch_vel_sat1 = ( pitch_ensemble[90, sat1] - pitch_ensemble[90-1, sat1] ) / (24. * 3600)
roll_vel_sat1 = ( roll_ensemble[90, sat1] - roll_ensemble[90-1, sat1] ) / (24. * 3600)

yaw_vel_sat2 = ( yaw_ensemble[90, sat2] - yaw_ensemble[90-1, sat2] ) / (24. * 3600)
pitch_vel_sat2 = ( pitch_ensemble[90, sat2] - pitch_ensemble[90-1, sat2] ) / (24. * 3600)
roll_vel_sat2 = ( roll_ensemble[90, sat2] - roll_ensemble[90-1, sat2] ) / (24. * 3600)

dist_sat1_sat2 = np.sqrt( ( x_eci_ensemble[:, sat1] - x_eci_ensemble[:, sat2] )**2 +  ( y_eci_ensemble[:, sat1] - y_eci_ensemble[:, sat2] )**2 +  ( z_eci_ensemble[:, sat1] - z_eci_ensemble[:, sat2] )**2 )

print pitch_vel_sat1, roll_vel_sat1, yaw_vel_sat1
print pitch_vel_sat2, roll_vel_sat2, yaw_vel_sat2

# 200/80
# In [86]: print pitch_vel_sat1, roll_vel_sat1, yaw_vel_sat1
# -0.118257294919 -2.12316849306 -0.921393159606
# In [87]: print pitch_vel_sat2, roll_vel_sat2, yaw_vel_sat2
# 0.182033980787 0.677633533565 -2.12392994444


fontsize_plot = 16
height_fig = 12 
ratio_fig_size = 4./3
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig_title = 'Distance between sat1 and sat2 after deployment over 90 days'
fig.suptitle(fig_title, y = 0.953,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

ax = fig.add_subplot(gs[0, 0])
ax.plot(dist_sat1_sat2, linewidth = 2, color = 'k')
ax.set_ylabel('Distance (km)', weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel('Time after deployment (days)', weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
fig_save_name = fig_title.replace(" ","_").lower()
fig_save_name = root_save_fig_name + fig_save_name + '.png'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission)


########################################################################################################################################################################################################################################################################################################
if ( 'angle_asc_node_to_sat_ensemble' in pickle_list_name ):
    # Plot the distribution of angular_spacing_between_ensemble_sat_converted_in_a_distance at different times (set by when_plot_in_hour)
    fontsize_plot = 14
    height_fig = 5 
    ratio_fig_size = 4./3
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left= 0.08, right=0.98, top = 0.93,bottom = 0.08)
    bin_width = 10 # in m

    when_plot_in_hour = 24 * 50. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
    index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Distribution along-track after N hours', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
    ax1.set_ylabel('Distribution (# out of ' + str(nb_ensembles) + ')', fontsize = fontsize_plot, weight = 'bold')
    ax1.set_xlabel('Algebraic LVLH_X distance from reference satellite (km)', fontsize = fontsize_plot, weight = 'bold')
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    n_max, bins, patches = ax1.hist(angular_spacing_between_ensemble_sat_converted_in_a_distance[index_when_plot,:], bins = np.arange(min(angular_spacing_between_ensemble_sat_converted_in_a_distance[index_when_plot,:]), max(angular_spacing_between_ensemble_sat_converted_in_a_distance[index_when_plot,:]) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'k',label = '')  

    plt.legend(loc = 2, scatterpoints=1)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
#    ax1.set_ylim([0, 100])
#    ax1.set_xlim([min(bins), 0])
    ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)


    fig_save_name = 'angular_spacing_between_ensemble_sat_converted_in_a_distance'
    fig_save_name = root_save_fig_name + fig_save_name + '.png'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    




########################################################################################################################################################################################################################################################################################################
if ( ( 'x_eci' in ensemble_to_plot ) & ( 'y_eci' in ensemble_to_plot ) & ( 'z_eci' in ensemble_to_plot ) ): 
    # Plot the algebric lvlh_x distance at a given time
    ## Distance along the x axis
    fontsize_plot = 14
    height_fig = 5 
    ratio_fig_size = 4./3
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left= 0.08, right=0.98, top = 0.93,bottom = 0.08)
    bin_width = 5. # in km

#    when_plot_in_hour = 96. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
    index_when_plot = algebric_distance_ensemble_main_sc_lvlh_x.shape[0] - 1 #(int) (when_plot_in_hour * 3600L / dt)
    median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.)
    tenth_algebric_ditance = np.percentile(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000., 10)
    std = np.std(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.)    
    ax1 = fig.add_subplot(gs[0, 0])
    n, bins, patches = ax1.hist(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000., np.arange(min(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.), max(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'g',label = 'N = ' + '{0:.0f}' .format( ( input_variables[1] - input_variables[0] ).days*24 + ( input_variables[1] - input_variables[0] ).seconds / 3600 ) + 'h - ' + 'median = ' + '{0:.3f}' .format(median_algebric_ditance) + ' km, ' + 'std = ' + '{0:.1e}' .format(std) + ' km') 
    ax1.plot([median_algebric_ditance, median_algebric_ditance], [min(n), max(n)], 'g', linewidth = 2, linestyle = 'dotted')
    ## ADD THE TLE POSITION
    ax1.scatter([algebric_distance_tle_main_sc_lvlh_x/1000.],[1.5], s = 200, c = 'g', marker = '*', edgecolor = 'k', label = 'CADRE')


    plt.legend(loc = 2, scatterpoints=1)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.set_ylim([0, 100])
#    ax1.set_xlim([min(bins), 0])
    ax1.margins(0,1) # autoscale both axes(fist value is for the x axis, second value for the y axis)


    fig_save_name = 'distribution_along_track'
    fig_save_name = root_save_fig_name + fig_save_name + '.png'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    
    # when_plot_in_hour = 12. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
    # index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
    # median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.)
    # tenth_algebric_ditance = np.percentile(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000., 10)
    # std = np.std(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.)    
    # ax1 = fig.add_subplot(gs[0, 0])
    # ax1.set_title('Distribution along-track after N hours', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
    # ax1.set_ylabel('Distribution (# out of ' + str(nb_ensembles) + ')', fontsize = fontsize_plot, weight = 'bold')
    # ax1.set_xlabel('Algebraic LVLH_X distance from reference satellite (km)', fontsize = fontsize_plot, weight = 'bold')
    # ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    # [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    # n_max, bins, patches = ax1.hist(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000., bins = np.arange(min(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.), max(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'k',label = 'N = ' + '{0:.0f}' .format( ( index_when_plot * dt / 3600. )) + 'h - ' + 'median = ' + '{0:.3f}' .format(median_algebric_ditance) + ' km, ' + 'std = ' + '{0:.1e}' .format(std) + ' km')  
    # ax1.plot([median_algebric_ditance, median_algebric_ditance], [min(n_max), max(n_max)], 'k', linewidth = 2, linestyle = 'dotted')

    # when_plot_in_hour = 48. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
    # index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
    # median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.)
    # tenth_algebric_ditance = np.percentile(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000., 10)
    # std = np.std(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.)    
    # ax1 = fig.add_subplot(gs[0, 0])
    # n, bins, patches = ax1.hist(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000., np.arange(min(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.), max(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'b',label = 'N = ' + '{0:.0f}' .format( ( index_when_plot * dt / 3600. )) + 'h - ' + 'median = ' + '{0:.3f}' .format(median_algebric_ditance) + ' km, ' + 'std = ' + '{0:.1e}' .format(std) + ' km') 
    # ax1.plot([median_algebric_ditance, median_algebric_ditance], [min(n), max(n)], 'b', linewidth = 2, linestyle = 'dotted')

    # when_plot_in_hour = 96. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
    # index_when_plot = (int) (when_plot_in_hour * 3600L / dt) 
    # median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.)
    # tenth_algebric_ditance = np.percentile(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000., 10)
    # std = np.std(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.)    
    # ax1 = fig.add_subplot(gs[0, 0])
    # n, bins, patches = ax1.hist(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000., np.arange(min(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.), max(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'r',label = 'N = ' + '{0:.0f}' .format( ( index_when_plot * dt / 3600. )) + 'h - ' + 'median = ' + '{0:.3f}' .format(median_algebric_ditance) + ' km, ' + 'std = ' + '{0:.1e}' .format(std) + ' km') 
    # ax1.plot([median_algebric_ditance, median_algebric_ditance], [min(n), max(n)], 'r', linewidth = 2, linestyle = 'dotted')


########################################################################################################################################################################################################################################################################################################
if ( ( 'x_eci' in ensemble_to_plot ) & ( 'y_eci' in ensemble_to_plot ) & ( 'z_eci' in ensemble_to_plot ) ): 
    # Plot the standard deviation of the along-track distributions as a function of time
    ## Standard deviation of the along-track distribution as a function of time (by step of N hours (called step_std))
    step_std = 24 # step in hours to calculate the standard deviation
    step_std_in_index = step_std * 3600. / dt
    std_daily = np.zeros([nb_steps])
    med_daily_z = np.zeros([nb_steps])
    med_daily_x = np.zeros([nb_steps])
    min_r_ensemble_daily = np.zeros([nb_steps])
    for i in range(nb_steps):
        index_when_std = i * step_std_in_index
        std_daily[i] = np.std(algebric_distance_ensemble_main_sc_lvlh_x[index_when_std,:]) # in meters
        med_daily_z[i] = np.median(algebric_distance_ensemble_main_sc_lvlh_z[index_when_std,:]) # in meters
        med_daily_x[i] = np.median(algebric_distance_ensemble_main_sc_lvlh_x[index_when_std,:]) # in meters
        min_r_ensemble_daily[i] = np.median( r_eci_ensemble[index_when_std, :] )

    ## Plot
    fontsize_plot = 14
    height_fig = 8 
    ratio_fig_size = 4./3
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left= 0.09, right=0.98, top = 0.93,bottom = 0.08)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Standard deviation of the along-track distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
    ax1.set_ylabel('Standard Deviation (m)', fontsize = fontsize_plot, weight = 'bold')
    ax1.set_xlabel('Time (days)', fontsize = fontsize_plot, weight = 'bold')
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.plot(np.arange(0,nb_steps), med_daily_x[0:nb_steps]/1000., 'k', linewidth = 2)

    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

    fig_save_name = 'std_daily'
    fig_save_name = root_save_fig_name + fig_save_name + '.png'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


    ########################################################################################################################################################################################################################################################################################################
if ( ( 'pitch' in ensemble_to_plot ) & ( 'roll' in ensemble_to_plot ) & ( 'yaw' in ensemble_to_plot ) ): 
    # Plot the angular velocities distribution
    fontsize_plot = 14
    when_plot_in_hour = 10. * 24 #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
    index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
    median_algebric_ditance = np.median(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))
    tenth_algebric_ditance = np.percentile(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 10)
    std = np.std(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))    
    fig = plt.figure(num=None, figsize=(13, 9), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('Angular velocity distributions (pitch, roll, yaw)', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 3)
    gs.update(left= 0.05, right=0.99, top = 0.90,bottom = 0.09)

    ## Pitch angular velocity distribution
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Pitch ang. vel. dist.', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
    #Algebric LVLH_X distance after '+str( ( index_when_plot * dt / 3600. ))+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
    ax1.set_ylabel('Distribution (# out of ' + str(nb_ensembles) + ')', fontsize = fontsize_plot, weight = 'bold')
    ax1.set_xlabel('Pitch angular velocity' + u' (\N{DEGREE SIGN}/s)', fontsize = fontsize_plot, weight = 'bold')
    #ax1.xaxis.set_ticklabels([])
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
    n, bins, patches = ax1.hist(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 100,  histtype='stepfilled', alpha = 0.7) 
    ax1.plot([median_algebric_ditance/1000., median_algebric_ditance/1000.], [min(n), max(n)], 'k', linewidth = 2, linestyle = 'dotted')
#    ax1.text(ax1.get_xlim()[1], ax1.get_ylim()[1], 'std = ' + '{0:.2f}' .format(std) + u' (\N{DEGREE SIGN}/s)\n' , horizontalalignment = 'right', verticalalignment = 'top', fontsize = fontsize_plot, weight = 'bold')
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()

    ## Roll ang. vel. distribution
    median_algebric_ditance = np.median(roll_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))
    tenth_algebric_ditance = np.percentile(roll_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 10)
    std = np.std(roll_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))    
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_title('Roll ang. vel. dist.', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
    #Algebric LVLH_X distance after '+str( ( index_when_plot * dt / 3600. ))+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
    ax2.set_xlabel('Roll angular velocity' + u' (\N{DEGREE SIGN}/s)', fontsize = fontsize_plot, weight = 'bold')
    #ax2.xaxis.set_ticklabels([])
    ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
    n, bins, patches = ax2.hist(roll_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 100,  histtype='stepfilled', alpha = 0.7) 
    ax2.plot([median_algebric_ditance/1000., median_algebric_ditance/1000.], [min(n), max(n)], 'k', linewidth = 2, linestyle = 'dotted')
#    ax2.text(ax2.get_xlim()[1], ax2.get_ylim()[1], 'std = ' + '{0:.2f}' .format(std) + u' (\N{DEGREE SIGN}/s)\n' , horizontalalignment = 'right', verticalalignment = 'top', fontsize = fontsize_plot, weight = 'bold')
    ax2.get_xaxis().tick_bottom()
    ax2.get_yaxis().tick_left()

    ## Yaw ang. vel. distribution
    median_algebric_ditance = np.median(yaw_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))
    tenth_algebric_ditance = np.percentile(yaw_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 10)
    std = np.std(yaw_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))    
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.set_title('Yaw ang. vel. dist.', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
    #Algebric LVLH_X distance after '+str( ( index_when_plot * dt / 3600. ))+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
    ax3.set_xlabel('Yaw angular velocity' + u' (\N{DEGREE SIGN}/s)', fontsize = fontsize_plot, weight = 'bold')
    #ax3.xaxis.set_ticklabels([])
    ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
    ax3.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
    n, bins, patches = ax3.hist(yaw_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 100,  histtype='stepfilled', alpha = 0.7) 
    ax3.plot([median_algebric_ditance/1000., median_algebric_ditance/1000.], [min(n), max(n)], 'k', linewidth = 2, linestyle = 'dotted')
#    ax3.text(ax3.get_xlim()[1], ax3.get_ylim()[1], 'std = ' + '{0:.2f}' .format(std) + u' (\N{DEGREE SIGN}/s)\n' , horizontalalignment = 'right', verticalalignment = 'top', fontsize = fontsize_plot, weight = 'bold')
    ax3.get_xaxis().tick_bottom()
    ax3.get_yaxis().tick_left()

    fig_save_name = 'distribution_angular_velocity'
    fig_save_name = root_save_fig_name + fig_save_name + '.eps'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

########################################################################################################################################################################################################################################################################################################
if ( ( 'pitch' in ensemble_to_plot ) & ( 'roll' in ensemble_to_plot ) & ( 'yaw' in ensemble_to_plot ) ): 
    # Examples of attitude: pitch/roll/yaw as a function of time for the first orbit
    ## Pitch vs time
    fontsize_plot = 14
    when_plot_in_hour = 30./3600 #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
    index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
    median_algebric_ditance = np.median(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))
    tenth_algebric_ditance = np.percentile(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 10)
    std = np.std(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))    
    fig = plt.figure(num=None, figsize=(15, 10), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('Examples of the attitude of two satellites', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(3, 1)
    gs.update(left= 0.06, right=0.99, top = 0.93,bottom = 0.08)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Pitch', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
    #Algebric LVLH_X distance after '+str( ( index_when_plot * dt / 3600. ))+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
    ax1.set_ylabel('Pitch' + u' (\N{DEGREE SIGN})', fontsize = fontsize_plot, weight = 'bold')
    #ax1.xaxis.set_ticklabels([])
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
    x_axis = np.arange(0, 97, 0.5)
    sat_nb = 100
    ax1.plot(x_axis, np.mod( pitch_ensemble[0:14, sat_nb], 360 ), 'r', linewidth = 2, label = 's'+str(sat_nb))
    sat_nb = 500
    ax1.plot(x_axis, np.mod( pitch_ensemble[0:14, sat_nb], 360 ), 'b', linewidth = 2, label = 's'+str(sat_nb))
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()

    ## Roll vs time
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_title('Roll', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
    #Algebric LVLH_X distance after '+str( ( index_when_plot * dt / 3600. ))+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
    ax2.set_ylabel('Roll' + u' (\N{DEGREE SIGN})', fontsize = fontsize_plot, weight = 'bold')
    #ax2.xaxis.set_ticklabels([])
    ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
    sat_nb = 100
    ax2.plot(x_axis, np.mod( roll_ensemble[0:14, sat_nb], 360 ), 'r', linewidth = 2, label = 's'+str(sat_nb))
    sat_nb = 500
    ax2.plot(x_axis, np.mod( roll_ensemble[0:14, sat_nb], 360 ), 'b', linewidth = 2, label = 's'+str(sat_nb))
    ax2.get_xaxis().tick_bottom()
    ax2.get_yaxis().tick_left()

    ## Yaw vs time
    ax3 = fig.add_subplot(gs[2, 0])
    ax3.set_title('Yaw', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
    #Algebric LVLH_X distance after '+str( ( index_when_plot * dt / 3600. ))+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
    ax3.set_ylabel('Yaw' + u' (\N{DEGREE SIGN})', fontsize = fontsize_plot, weight = 'bold')
    #ax3.xaxis.set_ticklabels([])
    ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
    ax3.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
    sat_nb = 100
    ax3.plot(x_axis, np.mod( yaw_ensemble[0:14, sat_nb], 360 ), 'r', linewidth = 2, label = 's'+str(sat_nb))
    sat_nb = 500
    ax3.plot(x_axis, np.mod( yaw_ensemble[0:14, sat_nb], 360 ), 'b', linewidth = 2, label = 's'+str(sat_nb))
    ax3.get_xaxis().tick_bottom()
    ax3.get_yaxis().tick_left()
    ax3.set_xlabel('Time in orbit (min)', fontsize = fontsize_plot, weight = 'bold')

    fig_save_name = 'example_attitude_vs_time_ang_velo'
    fig_save_name = root_save_fig_name + fig_save_name + '.png'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

    raise Exception

########################################################################################################################################################################################################################################################################################################

if ( ( 'x_eci' in ensemble_to_plot ) & ( 'y_eci' in ensemble_to_plot ) & ( 'z_eci' in ensemble_to_plot ) ): 
    # Plot the algebric lvlh_z distance at a given time
    ## Distance along the x axis
    fontsize_plot = 14
    height_fig = 10 
    ratio_fig_size = 4./3
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left= 0.06, right=0.98, top = 0.93,bottom = 0.08)
    bin_width = 5/1000. # in km

    when_plot_in_hour = 12. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
    index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
    median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.)
    tenth_algebric_ditance = np.percentile(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000., 10)
    std = np.std(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.)    
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Distribution along-track after N hours', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
    ax1.set_ylabel('Distribution (# out of ' + str(nb_ensembles) + ')', fontsize = fontsize_plot, weight = 'bold')
    ax1.set_xlabel('Algebraic distance from reference satellite (km)', fontsize = fontsize_plot, weight = 'bold')
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    n, bins, patches = ax1.hist(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000., bins = np.arange(min(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.), max(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'k',label = 'N = ' + '{0:.0f}' .format( ( index_when_plot * dt / 3600. )) + 'h - ' + 'median = ' + '{0:.3f}' .format(median_algebric_ditance) + ' km, ' + 'std = ' + '{0:.1e}' .format(std) + ' km')  
    ax1.plot([median_algebric_ditance, median_algebric_ditance], [min(n), max(n)], 'k', linewidth = 2, linestyle = 'dotted')

    when_plot_in_hour = 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
    index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
    median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.)
    tenth_algebric_ditance = np.percentile(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000., 10)
    std = np.std(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.)    
    ax1 = fig.add_subplot(gs[0, 0])
    n, bins, patches = ax1.hist(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000., np.arange(min(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.), max(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'b',label = 'N = ' + '{0:.0f}' .format( ( index_when_plot * dt / 3600. )) + 'h - ' + 'median = ' + '{0:.3f}' .format(median_algebric_ditance) + ' km, ' + 'std = ' + '{0:.1e}' .format(std) + ' km') 
    ax1.plot([median_algebric_ditance, median_algebric_ditance], [min(n), max(n)], 'k', linewidth = 2, linestyle = 'dotted')

    when_plot_in_hour = 48. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
    index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
    median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.)
    tenth_algebric_ditance = np.percentile(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000., 10)
    std = np.std(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.)    
    ax1 = fig.add_subplot(gs[0, 0])
    n, bins, patches = ax1.hist(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000., np.arange(min(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.), max(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'r',label = 'N = ' + '{0:.0f}' .format( ( index_when_plot * dt / 3600. )) + 'h - ' + 'median = ' + '{0:.3f}' .format(median_algebric_ditance) + ' km, ' + 'std = ' + '{0:.1e}' .format(std) + ' km') 
    ax1.plot([median_algebric_ditance, median_algebric_ditance], [min(n), max(n)], 'k', linewidth = 2, linestyle = 'dotted')

    plt.legend(loc = 2)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.set_xlim([min(bins), 0])
    ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

    fig_save_name = 'distribution_radial'
    fig_save_name = root_save_fig_name + fig_save_name + '.png'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


raise Exception
########################################################################################################################################################################################################################################################################################################
##  ## EVERYTHING BELOW IS OLD 
########################################################################################################################################################################################################################################################################################################

raise Exception
########################################################################################################################################################################################################################################################################################################
## EVERYTHING BELOW IS OLD

# 2D Histograms: distribution in the along-track and radial directions
when_plot_in_hour = 12
index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
fig = plt.figure(num=None, figsize=(15, 7), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(111)
ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
#Algebric LVLH_X distance after '+str(when_plot_in_hour)+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
ax1.set_ylabel('Radial (m)', fontsize = fontsize_plot, weight = 'bold')
ax1.set_xlabel('Along-track (m)', fontsize = 18, weight = 'bold')
#ax1.xaxis.set_ticklabels([])
ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
x_hist2d = algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]
y_hist2d = algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]
H, xedges, yedges, img = ax1.hist2d(x_hist2d, y_hist2d, bins = 100) #, range = [[-1000, 1000],[-1000, 1000]]
plt.colorbar(img)
#extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
#im = ax1.imshow(H, cmap=plt.cm.jet, extent=extent)
#fig.colorbar(im, ax=ax1)
plt.show()
#cbar = plt.colorbar()


raise Exception

######################################################################################################################################################################################################################################################################################################## 
## EVERYTHING BELOW IS OLD

#y_pos_distance = 5
#plt.scatter(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:], np.zeros(nb_ensembles)-y_pos_distance, linewidth = 2)    

# # Plot the magnitude of the algebric distance oalong the LVLH_X axis
# ax2 = fig.add_subplot(122)
# ax2.set_title('Absolute distance', weight = 'bold', fontsize = 20,  y = 1.008) 
# n, bins, patches = plt.hist(np.abs(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]), 100,  histtype='stepfilled', alpha = 0.7) # 
# ax2.set_xlabel('Distance (m)', fontsize = 18, weight = 'bold')
# ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
# [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
# ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
# median_absolute_ditance = np.median(np.abs(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]))
# tenth_absolute_ditance = np.percentile(np.abs(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]), 10)
# ax2.plot([median_absolute_ditance, median_absolute_ditance], [min(n), max(n)], 'k', linewidth = 2)

raise Exception
## Distance along the z axis
when_plot_in_hour = 12 #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
index_when_plot = -1 #(int) (when_plot_in_hour * 3600L / dt)
median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:])
tenth_algebric_ditance = np.percentile(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:], 10)
fig = plt.figure(num=None, figsize=(20, 10), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle('Radial distribution', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(121)
ax1.set_title('Algebric distance', weight = 'bold', fontsize = 20,  y = 1.008) 
#Algebric LVLH_Z distance after '+str(when_plot_in_hour)+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
ax1.set_ylabel('Distribution (# out of ' + str(nb_ensembles) + ')', fontsize = 18, weight = 'bold')
ax1.set_xlabel('Distance (km)', fontsize = 18, weight = 'bold')
#ax1.xaxis.set_ticklabels([])
ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
n, bins, patches = ax1.hist(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000., 100,  histtype='stepfilled', alpha = 0.7) 
ax1.plot([median_algebric_ditance/1000., median_algebric_ditance/1000.], [min(n), max(n)], 'k', linewidth = 2)
#y_pos_distance = 5
#plt.scatter(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:], np.zeros(nb_ensembles)-y_pos_distance, linewidth = 2)    

# Plot the magnitude of the algebric distance oalong the LVLH_Z axis
ax2 = fig.add_subplot(122)
ax2.set_title('Absolute distance', weight = 'bold', fontsize = 20,  y = 1.008) 
n, bins, patches = plt.hist(np.abs(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]), 100,  histtype='stepfilled', alpha = 0.7) # 
ax2.set_xlabel('Distance (m)', fontsize = 18, weight = 'bold')
ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
[i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
median_absolute_ditance = np.median(np.abs(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]))
tenth_absolute_ditance = np.percentile(np.abs(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]), 10)
ax2.plot([median_absolute_ditance, median_absolute_ditance], [min(n), max(n)], 'k', linewidth = 2)


## Distance along the y axis
when_plot_in_hour = 12 #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
index_when_plot = -1 #(int) (when_plot_in_hour * 3600L / dt)
median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:])
tenth_algebric_ditance = np.percentile(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:], 10)
fig = plt.figure(num=None, figsize=(20, 10), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle('Cross-track distribution', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(121)
ax1.set_title('Algebric distance', weight = 'bold', fontsize = 20,  y = 1.008) 
#Algebric LVLH_Y distance after '+str(when_plot_in_hour)+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
ax1.set_ylabel('Distribution (# out of ' + str(nb_ensembles) + ')', fontsize = 18, weight = 'bold')
ax1.set_xlabel('Distance (km)', fontsize = 18, weight = 'bold')
#ax1.xaxis.set_ticklabels([])
ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
n, bins, patches = ax1.hist(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:]/1000., 100,  histtype='stepfilled', alpha = 0.7) 
ax1.plot([median_algebric_ditance/1000., median_algebric_ditance/1000.], [min(n), max(n)], 'k', linewidth = 2)
#y_pos_distance = 5
#plt.scatter(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:], np.zeros(nb_ensembles)-y_pos_distance, linewidth = 2)    

# Plot the magnitude of the algebric distance oalong the LVLH_Y axis
ax2 = fig.add_subplot(122)
ax2.set_title('Absolute distance', weight = 'bold', fontsize = 20,  y = 1.008) 
n, bins, patches = plt.hist(np.abs(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:]), 100,  histtype='stepfilled', alpha = 0.7) # 
ax2.set_xlabel('Distance (m)', fontsize = 18, weight = 'bold')
ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
[i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
median_absolute_ditance = np.median(np.abs(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:]))
tenth_absolute_ditance = np.percentile(np.abs(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:]), 10)
ax2.plot([median_absolute_ditance, median_absolute_ditance], [min(n), max(n)], 'k', linewidth = 2)


########################################################################################################################################################################################################################################################################################################
## EVERYTHING BELOW IS OLD

## Distance along the y axis
when_plot_in_hour = 12 #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
index_when_plot = -1 #(int) (when_plot_in_hour * 3600L / dt)
index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:])
fig_algebric_distance_ensemble_main_sc_lvlh_y = plt.figure()
plt.scatter(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:], np.zeros(nb_ensembles)-y_pos_distance, linewidth = 2)    
plt.title('Algebric LVLH_Y distance after '+str(when_plot_in_hour)+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m')


## Histogram
n, bins, patches = plt.hist(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:], 100,  histtype='stepfilled', alpha = 0.7) # 
plt.xlim([min(bins), max(bins)])
plt.ylim([-y_pos_distance, max(n)])
plt.plot([median_algebric_ditance, median_algebric_ditance], [min(n), max(n)], 'k', linewidth = 2)

plt.xlabel('Algebric distance (m)')
print median_algebric_ditance
print np.percentile(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:], 10)


## Distance along the z axis
when_plot_in_hour = 4
z_pos_distance = 5
index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:])
fig_algebric_distance_ensemble_main_sc_lvlh_z = plt.figure()
plt.scatter(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:], np.zeros(nb_ensembles)-y_pos_distance, linewidth = 2)    
plt.title('Algebric LVLH_Z distance after '+str(when_plot_in_hour)+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m')

## Histogram
n, bins, patches = plt.hist(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:], 100,  histtype='stepfilled', alpha = 0.7) # 
plt.xlim([min(bins), max(bins)])
plt.ylim([-y_pos_distance, max(n)])
plt.plot([median_algebric_ditance, median_algebric_ditance], [min(n), max(n)], 'k', linewidth = 2)

plt.xlabel('Algebric distance (m)')
print median_algebric_ditance
print np.percentile(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:], 10)

plt.show()

# raise Exception
# Distribution of the distance between the ensemble sc and the main sc     
median_distance_ensemble_main_sc = np.zeros(n)
tenth_percentile_distance_ensemble_main_sc = np.zeros(n)
twenty_fifth_percentile_distance_ensemble_main_sc = np.zeros(n)
seventy_fifth_percentile_distance_ensemble_main_sc = np.zeros(n)
ninetieth_percentile_distance_ensemble_main_sc = np.zeros(n)
for i in range(n):
    median_distance_ensemble_main_sc[i] = np.median( distance_ensemble_main_sc[i,:] )
    tenth_percentile_distance_ensemble_main_sc[i] = np.percentile( distance_ensemble_main_sc[i,:], 10)
    twenty_fifth_percentile_distance_ensemble_main_sc[i] = np.percentile( distance_ensemble_main_sc[i,:], 25)
    seventy_fifth_percentile_distance_ensemble_main_sc[i] = np.percentile( distance_ensemble_main_sc[i,:], 75)
    ninetieth_percentile_distance_ensemble_main_sc[i] = np.percentile( distance_ensemble_main_sc[i,:], 90)

# Plot the median and the percentiles of the distribution of the distance between the ensemble sc and the main sc        
x_axis = range(n)
for i in range(n):
    x_axis[i] = x_axis[i] / 60.0
fig = plt.figure()
plt.plot(x_axis,median_distance_ensemble_main_sc, linewidth = 2, color = 'k',label = 'Median')    
plt.plot(x_axis,tenth_percentile_distance_ensemble_main_sc, linewidth = 2, color = 'r',label = '10th percentile')    
plt.plot(x_axis,twenty_fifth_percentile_distance_ensemble_main_sc, linewidth = 2, color = 'b',label = '25th percentile')    
plt.plot(x_axis,seventy_fifth_percentile_distance_ensemble_main_sc, linewidth = 2, color = 'b',label = '75th percentile')    
plt.plot(x_axis,ninetieth_percentile_distance_ensemble_main_sc, linewidth = 2, color = 'r',label = '90th percentile')    
#plt.legend()
plt.ylabel('Distance (m)')
plt.xlabel('Time (hours)')
plt.title('Distance between the ensembles and the main spacecraft - Cd = 2.2 +- 0.3 - F10.7 = 150')
plt.show()

