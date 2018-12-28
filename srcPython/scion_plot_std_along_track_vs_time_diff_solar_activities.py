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


run_names = ['SCION_2_deg_s_f107_90_ap_7_90_days_10_s_500_ensembles.txt','SCION_2_deg_s_f107_120_ap_15_90_days_10_s_500_ensembles.txt', 'SCION_2_deg_s_f107_200_ap_80_90_days_10_s_500_ensembles.txt']
#run_names = ['SCION_4_0_deg_s_f107_120_ap_15.txt', 'SCION_5_0_deg_s_f107_120_ap_15.txt', 'SCION_0_8_deg_s_f107_120_ap_15.txt', 'SCION_2_5_deg_s_f107_120_ap_15.txt']

nb_runs = len(run_names)
name_mission = 'SCION' 
for irun in range(nb_runs):
    main_input_file_name = get_prop_dir(1) + "run_scion/input/main_input/" + run_names[irun]  # ex of sys.argv[1]: cadre_tumbling_0panel.txt'

    # read input file
    input_variables, order_input_variables = read_input_file(main_input_file_name)
    dt = input_variables[2]; nb_steps = input_variables[3]; satellite_to_plot_path = input_variables[6][0]; satellite_to_plot = input_variables[7][0]; nb_ensembles_coe = input_variables[12]; nb_ensembles_attitude = input_variables[13]; nb_ensembles_cd = input_variables[11]; ensemble_to_plot_temp = input_variables[14]; 
    n = nb_steps

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

    pickle_list = []
    pickle_list_name = []
    ## Nb of ensembles
    nb_spacecraft = 1#int(a_input_file_input[6][0:6])
    nb_ensembles_array = [nb_ensembles_coe, nb_ensembles_attitude, nb_ensembles_cd]
    nb_ensembles = np.max(nb_ensembles_array)
    for i in range(len(nb_ensembles_array)):
        if ( (nb_ensembles_array[i] > 0 ) &  (nb_ensembles_array[i] < nb_ensembles) ) :
            nb_ensembles = nb_ensembles_array[i]

# # #######################################
    save_pickle_name = '/raid3/Armada/Charles/python/' + name_mission + '/pickle/' + satellite_to_plot.replace(".txt",".pickle")
    root_save_fig_name = '/raid3/Armada/Charles/python/' + name_mission + '/result/image/' + satellite_to_plot.replace(".txt","_") 
    root_save_video_name = '/raid3/Armada/Charles/python/' + name_mission + '/result/video/' + satellite_to_plot.replace(".txt","_")


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
    if ( 'angle_asc_node_to_sat_ensemble' in pickle_list_name ):
        angle_asc_node_to_sat_ensemble = pickle_list[ pickle_list_name.index('angle_asc_node_to_sat_ensemble') ]


# ########################################################################################################################################################################################################################################################################################################
    # Plot the standard deviation of the along-track distributions as a function of time
    ## Standard deviation of the along-track distribution as a function of time (by step of N hours (called step_std))
    step_std = 24 # step in hours to calculate the standard deviation
    step_std_in_index = step_std * 3600. / dt
    std_daily = np.zeros([nb_steps])
    median_daily = np.zeros([nb_steps])
    mad_daily = np.zeros([nb_steps])
    iqr_daily = np.zeros([nb_steps])
    quadriles_1090_daily = np.zeros([nb_steps])
    for istep in range(nb_steps):
        index_when_std = istep * step_std_in_index
        std_daily[istep] = np.std(algebric_distance_ensemble_main_sc_lvlh_x[index_when_std,:]) # in meters
        median_daily[istep] = np.median(algebric_distance_ensemble_main_sc_lvlh_x[index_when_std,:]) # in meters
        mad_daily[istep] = np.median( np.abs( angle_asc_node_to_sat_ensemble[index_when_std, :] - np.median(angle_asc_node_to_sat_ensemble[index_when_std, :]) ) ) * 110
        iqr_daily[istep] = np.subtract(*np.percentile(angle_asc_node_to_sat_ensemble[index_when_std, :], [75, 25])) * 110
        quadriles_1090_daily[istep] = np.subtract(*np.percentile(angle_asc_node_to_sat_ensemble[index_when_std, :], [90, 10])) * 110
#np.subtract(*np.percentile(angle_asc_node_to_sat_ensemble[index_when_std, :], [50+34, 50-34])) #np.std(angle_asc_node_to_sat_ensemble[index_when_std, :]) #np.median(angular_spacing_between_ensemble_sat[index_when_std,:]) # in meters
    if irun == 0:
    ## Plot
        color_array = ['b','k','r']#, 'm']
        label_array = ['Quiet','Moderate','Strong']
#        label_array = ['SCION_4_0_deg_s_f107_120_ap_15.txt', 'SCION_5_0_deg_s_f107_120_ap_15.txt', 'SCION_0_8_deg_s_f107_120_ap_15.txt', 'SCION_2_5_deg_s_f107_120_ap_15.txt']
        fontsize_plot = 14
        height_fig = 8 
        ratio_fig_size = 4./3
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left= 0.09, right=0.98, top = 0.93,bottom = 0.08)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Interval width (50% and 80%) of the statistical dispersion of the spacecraft along the orbit as a function time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
        ax1.set_ylabel('Interquartile Range (km)', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_xlabel('Days after deployment', fontsize = fontsize_plot, weight = 'bold')
        ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.plot(np.arange(0,nb_steps), iqr_daily[0:nb_steps],color =  color_array[irun], label = label_array[irun] + ' - 50%', linewidth = 2)
    ax1.plot(np.arange(0,nb_steps), quadriles_1090_daily[0:nb_steps],color =  color_array[irun], label = label_array[irun] + ' - 80%', linewidth = 4, linestyle = "dotted")

    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
    if irun == nb_runs - 1:
        ax1.legend(loc = 2, ncol = 3)

fig_save_name = '/raid3/Armada/Charles/python/SCION/result/image/iqr_daily_diff_solar_activities.png'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission)

