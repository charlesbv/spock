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

## NOTE: this can be run if only one satellite was run (with ensembles of course)

## NOTE 1: to use this script, the only 2 parameters you have to set are:
## - the name of the main input file (first non-commented line below)
## - if yes or no it is the first time that you read the output files of this simulation: set first_time to 1 in that case. Set it to 0 if you already read the output files and saved the results in a pickle
## NOTE 2: this can be run if only ONE MAIN satellite was run (with ensembles of course)
## NOTE 3: the results (pickle, image, video) are stored in a subfolder of /raid3/Armada/Charles/python/. The subfolder is either 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', or 'other'. This code figures out which folder the simulation corresponds to by reading the name of the output file chosen by the user in the main input file of the propagator (third line of section #SPACECRAFTS): it tries to find 'cygnss', 'cadre', 'aerie', 'scion', or 'qb50' in the name of the output file. If it does not find it, then the results here will be stored in /raid3/Armada/Charles/python/other/

main_input_file_name = '/home/cbv/PropSim/input/main_input/' + sys.argv[1]  # ex of sys.argv[1]: cadre_tumbling_0panel.txt'


####### ##### FOR CADRE, WE WANT TO CALCULATE THE POSITION GIVEN BY TLE AND COMPARE IT TO THE MAIN SC
i = -1 # we look at the last position of the main sc
position_main_sc_eci = np.array([x_eci_main_sc[i], y_eci_main_sc[i], z_eci_main_sc[i]])
position_tle_sc_eci = cadre_read_last_tle('/home/cbv/PropSim/input/main_input/cadre_last_tle.txt')
delta_vector_eci = position_main_sc_eci - position_tle_sc_eci
velocity_main_sc_eci = np.array([vx_eci_main_sc[i], vy_eci_main_sc[i], vz_eci_main_sc[i]])
delta_vector_lvlh_of_main_sc = eci_to_lvlh(position_main_sc_eci, velocity_main_sc_eci, delta_vector_eci)
algebric_distance_tle_main_sc_lvlh_x= delta_vector_lvlh_of_main_sc[0] * 1000 # * 1000 to convert the distance from km to m 
algebric_distance_tle_main_sc_lvlh_y= delta_vector_lvlh_of_main_sc[1] * 1000 # * 1000 to convert the distance from km to m 
algebric_distance_tle_main_sc_lvlh_z= delta_vector_lvlh_of_main_sc[2] * 1000 # * 1000 to convert the distance from km to m 

raise Exception
########################################################################################################################################################################################################################################################################################################
# PLOTS ################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################


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
    # ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
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
    med_daily = np.zeros([nb_steps])
    for i in range(nb_steps):
        index_when_std = i * step_std_in_index
        std_daily[i] = np.std(algebric_distance_ensemble_main_sc_lvlh_x[index_when_std,:]) # in meters
        med_daily[i] = np.median(algebric_distance_ensemble_main_sc_lvlh_x[index_when_std,:]) # in meters

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
    ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.plot(np.arange(0,91), std_daily[0:91], 'k', linewidth = 2)

    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

    fig_save_name = 'std_daily'
    fig_save_name = root_save_fig_name + fig_save_name + '.png'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


    ########################################################################################################################################################################################################################################################################################################
if ( ( 'pitch' in ensemble_to_plot ) & ( 'roll' in ensemble_to_plot ) & ( 'yaw' in ensemble_to_plot ) ): 
    # Plot the angular velocities distribution
    ## Pitch angular velocity distribution
    fonsize_plot = 14
    when_plot_in_hour = 30./3600 #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
    index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
    median_algebric_ditance = np.median(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))
    tenth_algebric_ditance = np.percentile(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 10)
    std = np.std(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))    
    fig = plt.figure(num=None, figsize=(15, 7), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('Angular velocity distributions (pitch, roll, yaw)', fontsize = (int)(fonsize_plot*1.1), weight = 'bold')
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 3)
    gs.update(left= 0.05, right=0.99, top = 0.90,bottom = 0.09)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Pitch ang. vel. dist.', weight = 'bold', fontsize = fonsize_plot,  y = 1.008) 
    #Algebric LVLH_X distance after '+str( ( index_when_plot * dt / 3600. ))+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
    ax1.set_ylabel('Distribution (# out of ' + str(nb_ensembles) + ')', fontsize = fonsize_plot, weight = 'bold')
    ax1.set_xlabel('Pitch angular velocity' + u' (\N{DEGREE SIGN}/s)', fontsize = fonsize_plot, weight = 'bold')
    #ax1.xaxis.set_ticklabels([])
    ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
    n, bins, patches = ax1.hist(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 100,  histtype='stepfilled', alpha = 0.7) 
    ax1.plot([median_algebric_ditance/1000., median_algebric_ditance/1000.], [min(n), max(n)], 'k', linewidth = 2, linestyle = 'dotted')
    ax1.text(ax1.get_xlim()[1], ax1.get_ylim()[1], 'std = ' + '{0:.2f}' .format(std) + u' (\N{DEGREE SIGN}/s)\n' + 'median = ' + '{0:.2f}' .format(median_algebric_ditance) + u' (\N{DEGREE SIGN}/s)\n' , horizontalalignment = 'right', verticalalignment = 'top', fontsize = fonsize_plot, weight = 'bold')
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()

    ## Roll ang. vel. distribution
    median_algebric_ditance = np.median(roll_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))
    tenth_algebric_ditance = np.percentile(roll_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 10)
    std = np.std(roll_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))    
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_title('Roll ang. vel. dist.', weight = 'bold', fontsize = fonsize_plot,  y = 1.008) 
    #Algebric LVLH_X distance after '+str( ( index_when_plot * dt / 3600. ))+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
    ax2.set_xlabel('Roll angular velocity' + u' (\N{DEGREE SIGN}/s)', fontsize = fonsize_plot, weight = 'bold')
    #ax2.xaxis.set_ticklabels([])
    ax2.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
    n, bins, patches = ax2.hist(roll_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 100,  histtype='stepfilled', alpha = 0.7) 
    ax2.plot([median_algebric_ditance/1000., median_algebric_ditance/1000.], [min(n), max(n)], 'k', linewidth = 2, linestyle = 'dotted')
    ax2.text(ax2.get_xlim()[1], ax2.get_ylim()[1], 'std = ' + '{0:.2f}' .format(std) + u' (\N{DEGREE SIGN}/s)\n' + 'median = ' + '{0:.2f}' .format(median_algebric_ditance) + u' (\N{DEGREE SIGN}/s)\n' , horizontalalignment = 'right', verticalalignment = 'top', fontsize = fonsize_plot, weight = 'bold')
    ax2.get_xaxis().tick_bottom()
    ax2.get_yaxis().tick_left()

    ## Yaw ang. vel. distribution
    median_algebric_ditance = np.median(yaw_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))
    tenth_algebric_ditance = np.percentile(yaw_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 10)
    std = np.std(yaw_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))    
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.set_title('Yaw ang. vel. dist.', weight = 'bold', fontsize = fonsize_plot,  y = 1.008) 
    #Algebric LVLH_X distance after '+str( ( index_when_plot * dt / 3600. ))+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
    ax3.set_xlabel('Yaw angular velocity' + u' (\N{DEGREE SIGN}/s)', fontsize = fonsize_plot, weight = 'bold')
    #ax3.xaxis.set_ticklabels([])
    ax3.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
    ax3.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
    n, bins, patches = ax3.hist(yaw_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 100,  histtype='stepfilled', alpha = 0.7) 
    ax3.plot([median_algebric_ditance/1000., median_algebric_ditance/1000.], [min(n), max(n)],'k',  linewidth = 2, linestyle = 'dotted')
    ax3.text(ax3.get_xlim()[1], ax3.get_ylim()[1], 'std = ' + '{0:.2f}' .format(std) + u' (\N{DEGREE SIGN}/s)\n' + 'median = ' + '{0:.2f}' .format(median_algebric_ditance) + u' (\N{DEGREE SIGN}/s)\n' , horizontalalignment = 'right', verticalalignment = 'top', fontsize = fonsize_plot, weight = 'bold')
    ax3.get_xaxis().tick_bottom()
    ax3.get_yaxis().tick_left()

    fig_save_name = 'distribution_angular_velocity'
    fig_save_name = root_save_fig_name + fig_save_name + '.png'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

########################################################################################################################################################################################################################################################################################################
if ( ( 'pitch' in ensemble_to_plot ) & ( 'roll' in ensemble_to_plot ) & ( 'yaw' in ensemble_to_plot ) ): 
    # Examples of attitude: pitch/roll/yaw as a function of time for the first orbit
    ## Pitch vs time
    fonsize_plot = 14
    when_plot_in_hour = 30./3600 #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
    index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
    median_algebric_ditance = np.median(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))
    tenth_algebric_ditance = np.percentile(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 10)
    std = np.std(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))    
    fig = plt.figure(num=None, figsize=(15, 10), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('Examples of the attitude of two satellites', fontsize = (int)(fonsize_plot*1.1), weight = 'bold')
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(3, 1)
    gs.update(left= 0.06, right=0.99, top = 0.93,bottom = 0.08)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Pitch', weight = 'bold', fontsize = fonsize_plot,  y = 1.008) 
    #Algebric LVLH_X distance after '+str( ( index_when_plot * dt / 3600. ))+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
    ax1.set_ylabel('Pitch' + u' (\N{DEGREE SIGN})', fontsize = fonsize_plot, weight = 'bold')
    #ax1.xaxis.set_ticklabels([])
    ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
    x_axis = np.arange(0, 97, 0.5)
    sat_nb = 100
    ax1.plot(x_axis, np.mod( pitch_ensemble[0:194, sat_nb], 360 ), 'r', linewidth = 2, label = 's'+str(sat_nb))
    sat_nb = 500
    ax1.plot(x_axis, np.mod( pitch_ensemble[0:194, sat_nb], 360 ), 'b', linewidth = 2, label = 's'+str(sat_nb))
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()

    ## Roll vs time
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_title('Roll', weight = 'bold', fontsize = fonsize_plot,  y = 1.008) 
    #Algebric LVLH_X distance after '+str( ( index_when_plot * dt / 3600. ))+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
    ax2.set_ylabel('Roll' + u' (\N{DEGREE SIGN})', fontsize = fonsize_plot, weight = 'bold')
    #ax2.xaxis.set_ticklabels([])
    ax2.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
    sat_nb = 100
    ax2.plot(x_axis, np.mod( roll_ensemble[0:194, sat_nb], 360 ), 'r', linewidth = 2, label = 's'+str(sat_nb))
    sat_nb = 500
    ax2.plot(x_axis, np.mod( roll_ensemble[0:194, sat_nb], 360 ), 'b', linewidth = 2, label = 's'+str(sat_nb))
    ax2.get_xaxis().tick_bottom()
    ax2.get_yaxis().tick_left()

    ## Yaw vs time
    ax3 = fig.add_subplot(gs[2, 0])
    ax3.set_title('Yaw', weight = 'bold', fontsize = fonsize_plot,  y = 1.008) 
    #Algebric LVLH_X distance after '+str( ( index_when_plot * dt / 3600. ))+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
    ax3.set_ylabel('Yaw' + u' (\N{DEGREE SIGN})', fontsize = fonsize_plot, weight = 'bold')
    #ax3.xaxis.set_ticklabels([])
    ax3.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
    [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
    ax3.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
    sat_nb = 100
    ax3.plot(x_axis, np.mod( yaw_ensemble[0:194, sat_nb], 360 ), 'r', linewidth = 2, label = 's'+str(sat_nb))
    sat_nb = 500
    ax3.plot(x_axis, np.mod( yaw_ensemble[0:194, sat_nb], 360 ), 'b', linewidth = 2, label = 's'+str(sat_nb))
    ax3.get_xaxis().tick_bottom()
    ax3.get_yaxis().tick_left()
    ax3.set_xlabel('Time in orbit (min)', fontsize = fonsize_plot, weight = 'bold')

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
    ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
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
ax1.set_ylabel('Radial (m)', fontsize = fonsize_plot, weight = 'bold')
ax1.set_xlabel('Along-track (m)', fontsize = 18, weight = 'bold')
#ax1.xaxis.set_ticklabels([])
ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
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
# ax2.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
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
ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
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
ax2.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
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
ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
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
ax2.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
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
