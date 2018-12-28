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
# This script saves and plots orbital characteristics as predicted by SpOCK:
# - 'radius'
# - 'speed'
# - 'altitude'
# - 'latitude'
# - 'eccentricity'
# - 'sma'
# - 'sma_average' orbit average sma
# - 'radius_perigee'
# - 'radius_apogee'
# - 'argument_perigee'
# - 'raan'
# To run it:
# python state.py arg_run_dir arg_input_filename1 arg_input_filename2 arg_input_filename3 ... arg_input_filenameM plot arg_var_to_plot1 arg_var_to_plot2 ... arg_var_to_plotN save arg_var_to_save1 arg_var_to_save2 ... arg_var_to_saveP
# where:
# - arg_run_dir is the name of the run directory where the SpOCK simulation was made
# - arg_input_filename1 arg_input_filename2 arg_input_filename3 ... arg_input_filenameM is the name of the SpOCK main input file that was used to run the SpOCK simulations. So here put the name of the M simulations that you want to plot the caracteristics of (the name of the main input file)
#  arg_var_to_plot1 arg_var_to_plot2 ... arg_var_to_plotN the name of orbital characteristics you want to plot (see list at beginning of header). DO NO FORGET TO PUT 'plot' right after the last argument arg_input_filenameM and right before the first argument arg_var_to_plot1
#  arg_var_to_save1 arg_var_to_save2 ... arg_var_to_saveP the name of orbital characteristics you want to save (see list at beginning of header). DO NO FORGET TO PUT 'save' right after the last argument arg_input_filenameM and right before the first argument arg_var_to_save1
# Assumptions:
# - see block 'PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT' below
# - all SpOCK's simulations (defined by run_list below) to analyze here have to be in the same run directory (defined by arg_run_dir)
# Notes:
# - you can invert save and plot (so you can type: python state.py arg_run_dir arg_input_filename1 arg_input_filename2 arg_input_filename3 ... arg_input_filenameM save arg_var_to_save1 arg_var_to_save2 ... arg_var_to_saveP plot arg_var_to_plot1 arg_var_to_plot2 ... arg_var_to_plotN)

from orbit_average import *
import matplotlib.colors as colors
import matplotlib.cm as cmx
from cadre_read_last_tle import *
from get_prop_dir import *
import matplotlib.gridspec as gridspec
from read_input_file import *
from matplotlib.colors import LogNorm
from norad_id_to_cygnss_name import *
import pickle
from eci_to_lvlh import *
import sys
import fileinput
import time
from datetime import datetime
import numpy as np
from matplotlib import pyplot as plt
import os
import sys
import subprocess
from get_name_mission import *
from find_in_read_input_order_variables import *

############ PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT ############
## If cygnss is set to 1 then the label of the plot won't be 1, 2, 3, ..., nb_sc but FM01, FM02, ..., FM08 ASSUMING THAT:
# - sc 1 is FM05
# - sc 2 is FM04
# - sc 3 is FM02
# - sc 4 is FM01
# - sc 5 is FM08
# - sc 6 is FM06
# - sc 7 is FM07
# - sc 8 is FM03
# (this is because if CYGNSS is run with a TLE from space-track.org, this is the order in the TLE file (it is 41884, 41885, 41886, ..., 41891)) 
cygnss = 1

## Save or not the plots
save_plots = 1

## Show or not the plots
show_plots = 0

## path of the folder where you want to store the results (pickle, image, video)
### Run directory of SpOCK's simulations (no path, just the name). CALLED AS ARGUMENT OF SCRIPT SO DO NO MODIFY LINE BELOW
run_dir = sys.argv[1] #'run'
### This is the line you can modify for the path. Make sure this path exists before running this script
path_folder_results = get_prop_dir(1) + run_dir + '/output/python_out/' #get_prop_dir(2) + 'output/python_propagator/'

## If the second spacecraft (and third, fourth, etc.) was propagated from the same main input file in SpOCK as the first spacecraft, set same_spock_input_file to 1. Otherwise, set it to 0
same_spock_input_file = 1

## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3

# date start and stop of plot. If set to 0 then the start and end epochs of the propagation are used by default. PLEASE follow format YYYY-MM-DDTHH:MM:SS
date_start_plot = 0 # "2017-03-14T01:00:00"
date_stop_plot = 0 # "2017-03-14T02:00:00"

############ end of PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT ############


# ############ ALGORITHM ############
if show_plots == 1:
    plt.ion()


if 'plot' in sys.argv:
    plot_arg = sys.argv.index('plot') + 1
else:
    plot_arg = len(sys.argv) - 2

if 'save' in sys.argv:
    save_arg = sys.argv.index('save') + 1
else:
    save_arg = len(sys.argv) - 2


min_save_plot_arg = min([save_arg, plot_arg])
max_save_plot_arg = max([save_arg, plot_arg])

list_elements_can_be_plot_or_save = ['radius', 'speed', 'altitude', 'latitude', 'eccentricity', 'sma', 'sma_average', 'radius_perigee', 'radius_apogee', 'argument_perigee', 'raan']


not_in_list = 1
for iarg in sys.argv[plot_arg:]:
    if iarg in list_elements_can_be_plot_or_save:
        not_in_list = 0
for iarg in sys.argv[save_arg:]:
    if iarg in list_elements_can_be_plot_or_save:
        not_in_list = 0
if ( ( len(sys.argv[max_save_plot_arg:]) == 0 ) | ( not_in_list == 1 ) ):
    print "Choose between: "
    for ilist in range(len(list_elements_can_be_plot_or_save)):
        print list_elements_can_be_plot_or_save[ilist]
    raise Exception


## List of SpOCK's simulations (name of the main input file, without path). 
run_list = [i for i in sys.argv[2:min_save_plot_arg-1] #["ex_main_input_basic.txt"#41891_on_2017-01-22_cd_2_4.txt"#41891_on_2017-01-22.txt"#cd_CYGFM03_2017-01-22T22:15:17_128512.txt"#, "CYGFM05_100days_gravity_order_20.txt"#,"CYGFM05_100days_j2.txt", "CYGFM05_100days.txt" #"CYGFM05_100days.txt"
            ]

earth_mu     = 398600.4418; # gravitational parameter (km^3/s^2)
earth_radius = 6378.137; # mean equatorial radius (km)

nb_run = len(run_list)
for irun in range(nb_run):
    # Read input file sc
    input_filename =    get_prop_dir(1) + run_dir + '/input/main_input/' + run_list[irun]
    var_in, var_in_order = read_input_file(input_filename)
    output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
    output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
    dt = var_in[find_in_read_input_order_variables(var_in_order, 'dt')]; 
    nb_steps = var_in[find_in_read_input_order_variables(var_in_order, 'nb_steps')]; 
    nb_sc = var_in[find_in_read_input_order_variables(var_in_order, 'nb_sc')]; 
    # Convert SC # to CYGNSS name if cygnss is set to 1
    if cygnss == 1:
        label_arr = ['FM05', 'FM04', 'FM02', 'FM01', 'FM08', 'FM06', 'FM07', 'FM03'] 
    else:
        label_arr = []
        for isc in range(nb_sc):
            label_arr.append(str(isc+1))

#    nb_sc = 1 # !!!!!!!!!!!!!!!!!!!
    name_mission = get_name_mission(output_file_name_list[0])
    root_save_fig_name = path_folder_results + output_file_name_list[0].replace(".txt","_")

    # For plots, generate disctinct colors
    NCURVES = nb_sc
    np.random.seed(101)
    curves = [np.random.random(20) for i in range(NCURVES)]
    values = range(NCURVES)
    jet = cm = plt.get_cmap('jet') 
    cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

    # Read output file sc
    print sys.argv[plot_arg:]
    if len(sys.argv[plot_arg:]) > 0:
        var_to_read = list( sys.argv[min_save_plot_arg:max_save_plot_arg] )
        var_to_read.extend( x for x in sys.argv[max_save_plot_arg:] )
        raise Exception
        if 'sma_average' in sys.argv[plot_arg:]:
            if ('latitude' in var_to_read ) == False:
                var_to_read.append('latitude')
            if ('sma' in var_to_read ) == False:
                var_to_read.append('sma')
    for isc in range(nb_sc):
        var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
        if isc == 0: # same date for all sc
            date = var_out[find_in_read_input_order_variables(var_out_order, 'date')]
            date_ref = datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S")
            start_xaxis_label = 0
            # test if date_start_plot and date_stop_plot are within the range of inital epoch and end epoch of the propagation
            if date_start_plot != 0:
                if (datetime.strptime(date_start_plot, "%Y-%m-%dT%H:%M:%S") - date_ref).total_seconds() < 0:
                    print "!*** The start date for the plot has to be earlier than the start epoch of the propagation. Therefore, it is set by default to the start epoch of the propagation ***!"
                    date_start_plot = 0
            if date_stop_plot != 0:
                if (datetime.strptime(date_stop_plot, "%Y-%m-%dT%H:%M:%S") - datetime.strptime(date[-1], "%Y/%m/%d %H:%M:%S")).total_seconds() > 0:
                    print "!*** The stop date for the plot has to be older than the stop epoch of the propagation. Therefore, it is set by default to the stop epoch of the propagation ***!"
                    date_stop_plot = 0

            nb_seconds_in_simu = ( nb_steps - 1 ) * dt
            x_axis = np.arange(0, nb_seconds_in_simu + 1, dt) # all output files of one simulation have the same number of steps 
            if date_start_plot != 0:
                date_ref = datetime.strptime(date_start_plot, "%Y-%m-%dT%H:%M:%S")
                delta_date_start = datetime.strptime(date_start_plot, "%Y-%m-%dT%H:%M:%S") - datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S")  
                nb_seconds_between_epoch_start_and_date_start_plot = delta_date_start.days * 24 * 3600 + delta_date_start.seconds + delta_date_start.microseconds/10**6.
                start_xaxis_label = nb_seconds_between_epoch_start_and_date_start_plot

            if date_stop_plot == 0:
                if date_start_plot != 0:
                    nb_steps_between_epoch_start_and_date_start_plot = (int)(nb_seconds_between_epoch_start_and_date_start_plot/dt) + 1
                    nb_seconds_in_simu = ( nb_steps - nb_steps_between_epoch_start_and_date_start_plot ) * dt 
                else:
                    nb_seconds_in_simu = ( nb_steps - 1 ) * dt
            else:
                if date_start_plot == 0:
                    date_start_plot =  date_ref
                    delta_date = datetime.strptime(date_stop_plot, "%Y-%m-%dT%H:%M:%S") - date_start_plot
                else:
                    delta_date = datetime.strptime(date_stop_plot, "%Y-%m-%dT%H:%M:%S") - datetime.strptime(date_start_plot, "%Y-%m-%dT%H:%M:%S")
                nb_seconds_in_simu = delta_date.days * 24 * 3600 + delta_date.seconds + delta_date.microseconds/10**6.
#                 nb_steps_plot = (int)(nb_seconds_in_simu/dt) + 1
#                 nb_seconds_in_simu = nb_steps_plot * dt


        if 'radius' in var_to_read:
            if isc == 0:
                radius = np.zeros([nb_sc, nb_steps]) # all output files of one simulation have the same number of steps
            radius[isc, :] = var_out[find_in_read_input_order_variables(var_out_order, 'radius')]
        if 'speed' in var_to_read:
            if isc == 0:
                speed = np.zeros([nb_sc, nb_steps]) # all output files of one simulation have the same number of steps
            speed[isc, :] = var_out[find_in_read_input_order_variables(var_out_order, 'speed')]
        if 'altitude' in var_to_read:
            if isc == 0:
                altitude = np.zeros([nb_sc, nb_steps]) # all output files of one simulation have the same number of steps
            altitude[isc, :] = var_out[find_in_read_input_order_variables(var_out_order, 'altitude')]
        if 'latitude' in var_to_read:
            if isc == 0:
                latitude = np.zeros([nb_sc, nb_steps]) # all output files of one simulation have the same number of steps
            latitude[isc, :] = var_out[find_in_read_input_order_variables(var_out_order, 'latitude')]
        if 'eccentricity' in var_to_read:
            if isc == 0:
                eccentricity = np.zeros([nb_sc, nb_steps]) # all output files of one simulation have the same number of steps
            eccentricity[isc, :] = var_out[find_in_read_input_order_variables(var_out_order, 'eccentricity')]
        if 'sma' in var_to_read:
            if isc == 0:
                sma = np.zeros([nb_sc, nb_steps]) # all output files of one simulation have the same number of steps
            sma[isc, :] = var_out[find_in_read_input_order_variables(var_out_order, 'sma')]

        if 'sma_average' in var_to_read:
            if isc == 0:
                sma_average = []
                date_average_start_orbit_list = []
                date_average_end_orbit_list = []
                x_axis_average = []
            sma_orbit_averaged, time_averaged, index_time_averaged = orbit_average(sma[isc, :], latitude[isc, :], date )
            sma_average.append( sma_orbit_averaged ) # each sc might not have the same orbital period so the length of the array might not be the same between each sc
            date_average_start_orbit_list.append( np.array(time_averaged)[:,0] ) # take the date at the start of the bin
            date_average_end_orbit_list.append( np.array(time_averaged)[:,2] ) # take the date at the end of the bin
            x_axis_average_per_sc = []
            nb_orbit_for_this_sc = len(time_averaged)
            for iorbit in range(nb_orbit_for_this_sc):
                date_average_start_orbit = date_average_start_orbit_list[-1][iorbit]
                date_average_start_orbit = datetime.strptime( date_average_start_orbit, "%Y/%m/%d %H:%M:%S" )
                nb_seconds_between_start_orbit_and_date_start = ( date_average_start_orbit - datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S") ).total_seconds()
                date_average_end_orbit = date_average_end_orbit_list[-1][iorbit]
                date_average_end_orbit = datetime.strptime( date_average_end_orbit, "%Y/%m/%d %H:%M:%S" )
                nb_seconds_between_end_orbit_and_date_start = ( date_average_end_orbit - datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S") ).total_seconds()

                x_axis_average_per_sc.append( nb_seconds_between_start_orbit_and_date_start )
                x_axis_average_per_sc.append( nb_seconds_between_end_orbit_and_date_start )
            x_axis_average.append( x_axis_average_per_sc )
        if 'radius_perigee' in var_to_read:
            if isc == 0:
                radius_perigee = np.zeros([nb_sc, nb_steps]) # all output files of one simulation have the same number of steps
            radius_perigee[isc, :] = var_out[find_in_read_input_order_variables(var_out_order, 'radius_perigee')]
        if 'radius_apogee' in var_to_read:
            if isc == 0:
                radius_apogee = np.zeros([nb_sc, nb_steps]) # all output files of one simulation have the same number of steps
            radius_apogee[isc, :] = var_out[find_in_read_input_order_variables(var_out_order, 'radius_apogee')]
        if 'argument_perigee' in var_to_read:
            if isc == 0:
                argument_perigee = np.zeros([nb_sc, nb_steps]) # all output files of one simulation have the same number of steps
            argument_perigee[isc, :] = var_out[find_in_read_input_order_variables(var_out_order, 'argument_perigee')]
        if 'raan' in var_to_read:
            if isc == 0:
                raan = np.zeros([nb_sc, nb_steps]) # all output files of one simulation have the same number of steps
            raan[isc, :] = var_out[find_in_read_input_order_variables(var_out_order, 'raan')]

        # Radius
        if 'radius' in var_to_read:        
            if isc == 0:
                # Plot
                fig_title = 'Radius as a function of time'
                y_label = 'Radius (km)'
                x_label = 'Real time'
                fig_radius = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                fig_radius.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                gs = gridspec.GridSpec(1, 1)
                gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                ax_radius = fig_radius.add_subplot(gs[0, 0])

                ax_radius.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                ax_radius.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                [i.set_linewidth(2) for i in ax_radius.spines.itervalues()] # change the width of the frame of the figure
                ax_radius.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold

            colorVal = scalarMap.to_rgba(isc)
            y_axis = radius[isc,:] - earth_radius
            ax_radius.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

            if isc == nb_sc - 1:
                # x axis label is in real time
                ## all output files of one simulation have the same number of steps, and start at the same date
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
                ax_radius.xaxis.set_ticks(xticks)
                ax_radius.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                ax_radius.margins(0,0); ax_radius.set_xlim([min(xticks), max(xticks)])
        #        ax_radius.set_xlim([ax_radius.get_xlim()[0], most_recent_tle_among_all_sc])

                legend = ax_radius.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                legend.get_title().set_fontsize(str(fontsize_plot))

                if save_plots == 1:
                    fig_save_name = 'radius'
                    fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                    fig_radius.savefig(fig_save_name, facecolor=fig_radius.get_facecolor(), edgecolor='none', bbox_inches='tight')  

        # Speed
        if 'speed' in var_to_read:        
            if isc == 0:
                # Plot
                fig_title = 'Speed (ECI) as a function of time'
                y_label = 'Speed (km/s)'
                x_label = 'Real time'
                fig_speed = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                fig_speed.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                gs = gridspec.GridSpec(1, 1)
                gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                ax_speed = fig_speed.add_subplot(gs[0, 0])

                ax_speed.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                ax_speed.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                [i.set_linewidth(2) for i in ax_speed.spines.itervalues()] # change the width of the frame of the figure
                ax_speed.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold

            colorVal = scalarMap.to_rgba(isc)
            y_axis = speed[isc,:]
            ax_speed.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

            if isc == nb_sc - 1:
                # x axis label is in real time
                ## all output files of one simulation have the same number of steps, and start at the same date
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
                ax_speed.xaxis.set_ticks(xticks)
                ax_speed.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                ax_speed.margins(0,0); ax_speed.set_xlim([min(xticks), max(xticks)])
        #        ax_speed.set_xlim([ax_speed.get_xlim()[0], most_recent_tle_among_all_sc])

                legend = ax_speed.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                legend.get_title().set_fontsize(str(fontsize_plot))

                if save_plots == 1:
                    fig_save_name = 'speed'
                    fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                    fig_speed.savefig(fig_save_name, facecolor=fig_speed.get_facecolor(), edgecolor='none', bbox_inches='tight')  


        # Altitude
        if 'altitude' in var_to_read:        
            if isc == 0:
                # Plot
                fig_title = 'Altitude as a function of time'
                y_label = 'Altitude (km)'
                x_label = 'Real time'
                fig_altitude = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                fig_altitude.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                gs = gridspec.GridSpec(1, 1)
                gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                ax_altitude = fig_altitude.add_subplot(gs[0, 0])

                ax_altitude.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                ax_altitude.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                [i.set_linewidth(2) for i in ax_altitude.spines.itervalues()] # change the width of the frame of the figure
                ax_altitude.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold

            colorVal = scalarMap.to_rgba(isc)
            y_axis = altitude[isc,:]
            ax_altitude.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

            if isc == nb_sc - 1:
                # x axis label is in real time
                ## all output files of one simulation have the same number of steps, and start at the same date
                nb_ticks_xlabel = 10
                dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
                date_list_str = []
                date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                for i in range(len(xticks)):
                    if dt_xlabel >= 3*24*3600:
                        date_list_str.append( str(date_list[i])[5:10] )
                    else:
                        date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                ax_altitude.xaxis.set_ticks(xticks)
                ax_altitude.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                ax_altitude.margins(0,0); ax_altitude.set_xlim([min(xticks), max(xticks)])
        #        ax_altitude.set_xlim([ax_altitude.get_xlim()[0], most_recent_tle_among_all_sc])

                legend = ax_altitude.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                legend.get_title().set_fontsize(str(fontsize_plot))

                if save_plots == 1:
                    fig_save_name = 'altitude'
                    fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                    fig_altitude.savefig(fig_save_name, facecolor=fig_altitude.get_facecolor(), edgecolor='none', bbox_inches='tight')  

        # Latitude
        if 'latitude' in var_to_read:        
            if isc == 0:
                # Plot
                fig_title = 'Latitude as a function of time'
                y_label = 'Latitude ' + u'(\N{DEGREE SIGN})'
                x_label = 'Real time'
                fig_latitude = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                fig_latitude.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                gs = gridspec.GridSpec(1, 1)
                gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                ax_latitude = fig_latitude.add_subplot(gs[0, 0])

                ax_latitude.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                ax_latitude.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                [i.set_linewidth(2) for i in ax_latitude.spines.itervalues()] # change the width of the frame of the figure
                ax_latitude.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold

            colorVal = scalarMap.to_rgba(isc)
            y_axis = latitude[isc,:]
            ax_latitude.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

            if isc == nb_sc - 1:
                # x axis label is in real time
                ## all output files of one simulation have the same number of steps, and start at the same date
                nb_ticks_xlabel = 10
                dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
                date_list_str = []
                date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                for i in range(len(xticks)):
                    if dt_xlabel >= 3*24*3600:
                        date_list_str.append( str(date_list[i])[5:10] )
                    else:
                        date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                ax_latitude.xaxis.set_ticks(xticks)
                ax_latitude.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                ax_latitude.margins(0,0); ax_latitude.set_xlim([min(xticks), max(xticks)])
        #        ax_latitude.set_xlim([ax_latitude.get_xlim()[0], most_recent_tle_among_all_sc])

                legend = ax_latitude.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                legend.get_title().set_fontsize(str(fontsize_plot))

                if save_plots == 1:
                    fig_save_name = 'latitude'
                    fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                    fig_latitude.savefig(fig_save_name, facecolor=fig_latitude.get_facecolor(), edgecolor='none', bbox_inches='tight')  

        # Eccentricity
        if 'eccentricity' in var_to_read:        
            if isc == 0:
                # Plot
                fig_title = 'Eccentricity as a function of time'
                y_label = 'Eccentricity'
                x_label = 'Real time'
                fig_eccentricity = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                fig_eccentricity.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                gs = gridspec.GridSpec(1, 1)
                gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                ax_eccentricity = fig_eccentricity.add_subplot(gs[0, 0])

                ax_eccentricity.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                ax_eccentricity.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                [i.set_linewidth(2) for i in ax_eccentricity.spines.itervalues()] # change the width of the frame of the figure
                ax_eccentricity.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold

            colorVal = scalarMap.to_rgba(isc)
            y_axis = eccentricity[isc,:]
            ax_eccentricity.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

            if isc == nb_sc - 1:
                # x axis label is in real time
                ## all output files of one simulation have the same number of steps, and start at the same date
                nb_ticks_xlabel = 10
                dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
                date_list_str = []
                date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                for i in range(len(xticks)):
                    if dt_xlabel >= 3*24*3600:
                        date_list_str.append( str(date_list[i])[5:10] )
                    else:
                        date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                ax_eccentricity.xaxis.set_ticks(xticks)
                ax_eccentricity.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                ax_eccentricity.margins(0,0); ax_eccentricity.set_xlim([min(xticks), max(xticks)])
        #        ax_eccentricity.set_xlim([ax_eccentricity.get_xlim()[0], most_recent_tle_among_all_sc])

                legend = ax_eccentricity.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                legend.get_title().set_fontsize(str(fontsize_plot))

                if save_plots == 1:
                    fig_save_name = 'eccentricity'
                    fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                    fig_eccentricity.savefig(fig_save_name, facecolor=fig_eccentricity.get_facecolor(), edgecolor='none', bbox_inches='tight')  


        # SMA
        if 'sma' in var_to_read:        
            if isc == 0:
                # Plot
                fig_title = 'SMA as a function of time'
                y_label = 'SMA (km)'
                x_label = 'Real time'
                fig_sma = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                fig_sma.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                gs = gridspec.GridSpec(1, 1)
                gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                ax_sma = fig_sma.add_subplot(gs[0, 0])

                ax_sma.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                ax_sma.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                [i.set_linewidth(2) for i in ax_sma.spines.itervalues()] # change the width of the frame of the figure
                ax_sma.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold

            colorVal = scalarMap.to_rgba(isc)
            y_axis = sma[isc,:] 
            ax_sma.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

#             # !!!!!!!!!!!!!! to overplot the sma_average
#             y_axis = np.zeros(len(sma_average[isc])*2)
#             for iorbit in range(len(sma_average[isc])): # y_axis is constant along the orbit: it has the same value at the beginning and at the end of the orbit so it looks like a straight line during the orbit
#                 y_axis[iorbit*2] = sma_average[isc][iorbit]
#                 y_axis[iorbit*2+1] = sma_average[isc][iorbit]
#             ax_sma.plot(x_axis_average[isc], y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])
#             # !!!!!!!!!!!!!!
            
            if isc == nb_sc - 1:
                # x axis label is in real time
                ## all output files of one simulation have the same number of steps, and start at the same date
                nb_ticks_xlabel = 10
                dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel); 
                date_list_str = []
                date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                for i in range(len(xticks)):
                    if dt_xlabel >= 3*24*3600:
                        date_list_str.append( str(date_list[i])[5:10] )
                    else:
                        date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                ax_sma.xaxis.set_ticks(xticks)
                ax_sma.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                ax_sma.margins(0,0); ax_sma.set_xlim([min(xticks), max(xticks)])
        #        ax_sma.set_xlim([ax_sma.get_xlim()[0], most_recent_tle_among_all_sc])

                legend = ax_sma.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                legend.get_title().set_fontsize(str(fontsize_plot))

                if save_plots == 1:
                    fig_save_name = 'sma'
                    fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                    fig_sma.savefig(fig_save_name, facecolor=fig_sma.get_facecolor(), edgecolor='none', bbox_inches='tight')  

        # Orbit average SMA
        if 'sma_average' in var_to_read:        
            if isc == 0:
                # Plot
                fig_title = 'Orbit average SMA as a function of time'
                y_label = 'SMA (km)'
                x_label = 'Real time'
                fig_sma_average = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                fig_sma_average.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                gs = gridspec.GridSpec(1, 1)
                gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                ax_sma_average = fig_sma_average.add_subplot(gs[0, 0])

                ax_sma_average.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                ax_sma_average.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                [i.set_linewidth(2) for i in ax_sma_average.spines.itervalues()] # change the width of the frame of the figure
                ax_sma_average.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold

            colorVal = scalarMap.to_rgba(isc)
            y_axis = np.zeros(len(sma_average[isc])*2)
            for iorbit in range(len(sma_average[isc])): # y_axis is constant along the orbit: it has the same value at the beginning and at the end of the orbit so it looks like a straight line during the orbit
                y_axis[iorbit*2] = sma_average[isc][iorbit]
                y_axis[iorbit*2+1] = sma_average[isc][iorbit]
            ax_sma_average.plot(x_axis_average[isc], y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

            if isc == nb_sc - 1:
                # x axis label is in real time
                ## all output files of one simulation have the same number of steps, and start at the same date
                nb_ticks_xlabel = 10
                dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
                date_list_str = []
                date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                for i in range(len(xticks)):
                    if dt_xlabel >= 3*24*3600:
                        date_list_str.append( str(date_list[i])[5:10] )
                    else:
                        date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                ax_sma_average.xaxis.set_ticks(xticks)
                ax_sma_average.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                ax_sma_average.margins(0,0); ax_sma_average.set_xlim([min(xticks), max(xticks)])
        #        ax_sma_average.set_xlim([ax_sma_average.get_xlim()[0], most_recent_tle_among_all_sc])

                legend = ax_sma_average.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                legend.get_title().set_fontsize(str(fontsize_plot))

                if save_plots == 1:
                    fig_save_name = 'sma_average'
                    fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                    fig_sma_average.savefig(fig_save_name, facecolor=fig_sma_average.get_facecolor(), edgecolor='none', bbox_inches='tight')  






        # Radius of perigee
        if 'radius_perigee' in var_to_read:        
            if isc == 0:
                # Plot
                fig_title = 'Radius of perigee as a function of time'
                y_label = 'Radius of perigee (km)'
                x_label = 'Real time'
                fig_radius_perigee = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                fig_radius_perigee.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                gs = gridspec.GridSpec(1, 1)
                gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                ax_radius_perigee = fig_radius_perigee.add_subplot(gs[0, 0])

                ax_radius_perigee.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                ax_radius_perigee.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                [i.set_linewidth(2) for i in ax_radius_perigee.spines.itervalues()] # change the width of the frame of the figure
                ax_radius_perigee.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold

            colorVal = scalarMap.to_rgba(isc)
            y_axis = radius_perigee[isc,:] - earth_radius
            ax_radius_perigee.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

            if isc == nb_sc - 1:
                # x axis label is in real time
                ## all output files of one simulation have the same number of steps, and start at the same date
                nb_ticks_xlabel = 10
                dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
                date_list_str = []
                date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                for i in range(len(xticks)):
                    if dt_xlabel >= 3*24*3600:
                        date_list_str.append( str(date_list[i])[5:10] )
                    else:
                        date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                ax_radius_perigee.xaxis.set_ticks(xticks)
                ax_radius_perigee.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                ax_radius_perigee.margins(0,0); ax_radius_perigee.set_xlim([min(xticks), max(xticks)])
        #        ax_radius_perigee.set_xlim([ax_radius_perigee.get_xlim()[0], most_recent_tle_among_all_sc])

                legend = ax_radius_perigee.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                legend.get_title().set_fontsize(str(fontsize_plot))

                if save_plots == 1:
                    fig_save_name = 'radius_perigee'
                    fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                    fig_radius_perigee.savefig(fig_save_name, facecolor=fig_radius_perigee.get_facecolor(), edgecolor='none', bbox_inches='tight')  

        # Radius of apogee
        if 'radius_apogee' in var_to_read:        
            if isc == 0:
                # Plot
                fig_title = 'Radius of apogee as a function of time'
                y_label = 'Radius of apogee (km)'
                x_label = 'Real time'
                fig_radius_apogee = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                fig_radius_apogee.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                gs = gridspec.GridSpec(1, 1)
                gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                ax_radius_apogee = fig_radius_apogee.add_subplot(gs[0, 0])

                ax_radius_apogee.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                ax_radius_apogee.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                [i.set_linewidth(2) for i in ax_radius_apogee.spines.itervalues()] # change the width of the frame of the figure
                ax_radius_apogee.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold

            colorVal = scalarMap.to_rgba(isc)
            y_axis = radius_apogee[isc,:] - earth_radius
            ax_radius_apogee.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

            if isc == nb_sc - 1:
                # x axis label is in real time
                ## all output files of one simulation have the same number of steps, and start at the same date
                nb_ticks_xlabel = 10
                dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
                date_list_str = []
                date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                for i in range(len(xticks)):
                    if dt_xlabel >= 3*24*3600:
                        date_list_str.append( str(date_list[i])[5:10] )
                    else:
                        date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )

                ax_radius_apogee.xaxis.set_ticks(xticks)
                ax_radius_apogee.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                ax_radius_apogee.margins(0,0); ax_radius_apogee.set_xlim([min(xticks), max(xticks)])
        #        ax_radius_apogee.set_xlim([ax_radius_apogee.get_xlim()[0], most_recent_tle_among_all_sc])

                legend = ax_radius_apogee.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                legend.get_title().set_fontsize(str(fontsize_plot))

                if save_plots == 1:
                    fig_save_name = 'radius_apogee'
                    fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                    fig_radius_apogee.savefig(fig_save_name, facecolor=fig_radius_apogee.get_facecolor(), edgecolor='none', bbox_inches='tight')  


        # Argument of perigee
        if 'argument_perigee' in var_to_read:        
            if isc == 0:
                # Plot
                fig_title = 'Argument of perigee as a function of time'
                y_label = 'Argument of perigee '+ u'(\N{DEGREE SIGN})'
                x_label = 'Real time'
                fig_argument_perigee = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                fig_argument_perigee.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                gs = gridspec.GridSpec(1, 1)
                gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                ax_argument_perigee = fig_argument_perigee.add_subplot(gs[0, 0])

                ax_argument_perigee.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                ax_argument_perigee.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                [i.set_linewidth(2) for i in ax_argument_perigee.spines.itervalues()] # change the width of the frame of the figure
                ax_argument_perigee.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold

            colorVal = scalarMap.to_rgba(isc)
            y_axis = argument_perigee[isc,:] 
            ax_argument_perigee.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

            if isc == nb_sc - 1:
                # x axis label is in real time
                ## all output files of one simulation have the same number of steps, and start at the same date
                nb_ticks_xlabel = 10
                dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
                date_list_str = []
                date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                for i in range(len(xticks)):
                    if dt_xlabel >= 3*24*3600:
                        date_list_str.append( str(date_list[i])[5:10] )
                    else:
                        date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                ax_argument_perigee.xaxis.set_ticks(xticks)
                ax_argument_perigee.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                ax_argument_perigee.margins(0,0); ax_argument_perigee.set_xlim([min(xticks), max(xticks)])
        #        ax_argument_perigee.set_xlim([ax_argument_perigee.get_xlim()[0], most_recent_tle_among_all_sc])

                legend = ax_argument_perigee.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                legend.get_title().set_fontsize(str(fontsize_plot))

                if save_plots == 1:
                    fig_save_name = 'argument_perigee'
                    fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                    fig_argument_perigee.savefig(fig_save_name, facecolor=fig_argument_perigee.get_facecolor(), edgecolor='none', bbox_inches='tight')  


        # RAAN
        if 'raan' in var_to_read:        
            if isc == 0:
                # Plot
                fig_title = 'Argument of perigee as a function of time'
                y_label = 'Argument of perigee '+ u'(\N{DEGREE SIGN})'
                x_label = 'Real time'
                fig_raan = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                fig_raan.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                gs = gridspec.GridSpec(1, 1)
                gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                ax_raan = fig_raan.add_subplot(gs[0, 0])

                ax_raan.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                ax_raan.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                [i.set_linewidth(2) for i in ax_raan.spines.itervalues()] # change the width of the frame of the figure
                ax_raan.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold

            colorVal = scalarMap.to_rgba(isc)
            y_axis = raan[isc,:] 
            ax_raan.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc] )

            if isc == nb_sc - 1:
                # x axis label is in real time
                ## all output files of one simulation have the same number of steps, and start at the same date
                nb_ticks_xlabel = 10
                dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel) 
                date_list_str = []
                date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                for i in range(len(xticks)):
                    if dt_xlabel >= 3*24*3600:
                        date_list_str.append( str(date_list[i])[5:10] )
                    else:
                        date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                ax_raan.xaxis.set_ticks(xticks)
                ax_raan.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                ax_raan.margins(0,0); ax_raan.set_xlim([min(xticks), max(xticks)])
        #        ax_raan.set_xlim([ax_raan.get_xlim()[0], most_recent_tle_among_all_sc])

                legend = ax_raan.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                legend.get_title().set_fontsize(str(fontsize_plot))

                if save_plots == 1:
                    fig_save_name = 'raan'
                    fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                    fig_raan.savefig(fig_save_name, facecolor=fig_raan.get_facecolor(), edgecolor='none', bbox_inches='tight')  



if show_plots == 1:
    plt.show(); plt.show();
