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
# - 'beta'
# - 'radius'
# - 'speed'
# - 'altitude'
# - 'latitude'
# - 'longitude'
# - 'eccentricity'
# - 'sma'
# - 'sma_average' orbit average sma
# - 'radius_perigee'
# - 'radius_apogee'
# - 'argument_perigee'
# - 'raan'
# - 'given_output' # !!!!! first date of given output is not necesarrily the same as for the other variables  
# - 'power'
# - 'density'
# - 'density_average'
# - 'temperature' # temperature of the atmosphere at the satellite postion (from NRLSMSIS if the user didn't chose "density_file" or "gitm", in which case SC->INTEGRATOR.Ta = 800K)
# - 'cd' # total normalized Cd over all surfaces of the sc
# - 'tot_area_drag'
# - 'position'
# - 'velocity'
# - 'acceleration_lvlh_drag_mag'
# - 'acceleration'
# - 'phase_angle'
# - 'solar_zenith'
# To run it:
# python state.py arg_run_dir arg_input_filename1 arg_input_filename2 arg_input_filename3 ... arg_input_filenameM plot save arg_var_to_plot1 arg_var_to_plot2 ... arg_var_to_plotN 
# where:
# - arg_input_filename1 arg_input_filename2 arg_input_filename3 ... arg_input_filenameM is the name of the SpOCK main input file that was used to run the SpOCK simulations. So here put the name of the M simulations that you want to plot the caracteristics of (the name of the main input file)
#  arg_var_to_plot1 arg_var_to_plot2 ... arg_var_to_plotN the name of orbital characteristics you want to plot or save (see list at beginning of header). DO NO FORGET TO PUT 'plot' or 'save' right after the last argument arg_input_filenameM and right before the first argument arg_var_to_plot1. If you put only 'plot', then the orbital chracteristic won't be saved. If you put only 'save', then they won't be plot. If you put both 'plot' and 'save', then they'll be both plotted and saved. Plots are in pdf format, saves are in pickle in format.
# Assumptions:
# - see block 'PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT' below

import sys
sys.path.append("/Users/cbv/work/spock/srcPython")
import pickle
from orbit_average import *
import matplotlib.colors as colors
import matplotlib.cm as cmx
from cadre_read_last_tle import *
from get_prop_dir import *
import matplotlib.gridspec as gridspec
from read_input_file import *
from read_output_file import *
from matplotlib.colors import LogNorm
from norad_id_to_cygnss_name import *
import pickle
from eci_to_lvlh import *
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
cygnss = 0

## Save or not the plots
save_plots = 1

## Show or not the plots
show_plots = 0

## path of the folder where you want to store the results (pickle, image, video)
### This is the line you can modify for the path. Make sure this path exists before running this script
path_folder_results = './' #get_prop_dir(2) + 'output/python_propagator/'

## If the second spacecraft (and third, fourth, etc.) was propagated from the same main input file in SpOCK as the first spacecraft, set same_spock_input_file to 1. Otherwise, set it to 0
same_spock_input_file = 0

## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25
ratio_fig_size = 4./3


# date start and stop of plot. If set to 0 then the start and end epochs of the propagation are used by default. PLEASE follow format YYYY-MM-DDTHH:MM:SS
date_start_plot = 0#"2019-04-22T00:00:00.000000"
date_stop_plot = 0#"2019-04-26T00:00:00.000000"

############ end of PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT ############


# ############ ALGORITHM ############
if show_plots == 1:
    plt.ion()

plot_arg = 0; want_plot = 0
save_arg = 0; want_save = 0
if 'plot' in sys.argv:
    plot_arg = sys.argv.index('plot') + 1
    want_plot = 1
if 'save' in sys.argv:
    save_arg = sys.argv.index('save') + 1
    want_save = 1
if ( ( 'plot' in sys.argv ) | ( 'save' in sys.argv ) ):
    plot_or_save_arg = max([plot_arg, save_arg])
else:
    plot_or_save_arg = len(sys.argv) - 2


list_elements_can_be_plot_or_save = ['radius', 'speed', 'altitude', 'latitude', 'longitude', 'eccentricity', 'sma', 'sma_average', 'radius_perigee', 'radius_apogee', 'argument_perigee', 'local_time', 'raan', 'beta', 'given_output', 'power', 'density', 'density_average', 'temperature', 'cd', 'tot_area_drag', 'position','velocity','acceleration_lvlh_drag_mag', 'phase_angle', 'solar_zenith', 'acceleration'] # !!!!! first date of given output is not necesarrily the same as for the other variables  


not_in_list = 1
for iarg in sys.argv[plot_or_save_arg:]:
    if iarg in list_elements_can_be_plot_or_save:
        not_in_list = 0
if ( ( len(sys.argv[plot_or_save_arg:]) == 0 ) | ( not_in_list == 1 ) ):
    print "To plot or save, choose between: "
    for ilist in range(len(list_elements_can_be_plot_or_save)):
        print list_elements_can_be_plot_or_save[ilist]
    raise Exception


## List of SpOCK's simulations (name of the main input file, without path). 
if ( ( 'plot' in sys.argv ) & ( 'save' in sys.argv ) ):
    run_list = [i for i in sys.argv[2:plot_or_save_arg-2] #["ex_main_input_basic.txt"#41891_on_2017-01-22_cd_2_4.txt"#41891_on_2017-01-22.txt"#cd_CYGFM03_2017-01-22T22:15:17_128512.txt"#, "CYGFM05_100days_gravity_order_20.txt"#,"CYGFM05_100days_j2.txt", "CYGFM05_100days.txt" #"CYGFM05_100days.txt"
            ]
else:
    run_list = [i for i in sys.argv[2:plot_or_save_arg-1] #["ex_main_input_basic.txt"#41891_on_2017-01-22_cd_2_4.txt"#41891_on_2017-01-22.txt"#cd_CYGFM03_2017-01-22T22:15:17_128512.txt"#, "CYGFM05_100days_gravity_order_20.txt"#,"CYGFM05_100days_j2.txt", "CYGFM05_100days.txt" #"CYGFM05_100days.txt"
            ]

earth_mu     = 398600.4418; # gravitational parameter (km^3/s^2)
earth_radius = 6378.137; # mean equatorial radius (km)

nb_run = len(run_list)
isc_irun = -1
for irun in range(nb_run):
    # Read input file sc
    input_filename =     run_list[irun]
    var_in, var_in_order = read_input_file(input_filename)
    output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
    output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 

    dt = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')]; 
    nb_steps = var_in[find_in_read_input_order_variables(var_in_order, 'nb_steps')]; 
    nb_sc = var_in[find_in_read_input_order_variables(var_in_order, 'nb_sc')]; 

    # Convert SC # to CYGNSS name if cygnss is set to 1
    if cygnss == 1:
        label_arr = ['FM05', 'FM04', 'FM02', 'FM01', 'FM08', 'FM06', 'FM07', 'FM03'] 
    else:
        label_arr = []
        if 'quantile' in run_list[irun]:
            label_arr.append(run_list[irun].split('quantile')[1][0])
        else:
            for isc in range(nb_sc):
                #label_arr.append(str(isc+1))
                label_arr.append(run_list[irun].replace('.txt',''))
            

    name_mission = get_name_mission(output_file_name_list[0])
    root_save_fig_name = path_folder_results + input_filename.replace(".txt","") + "_"

    # For plots, generate disctinct colors
    NCURVES = nb_sc*nb_run
    np.random.seed(101)
    curves = [np.random.random(20) for i in range(NCURVES)]
    values = range(NCURVES)
    jet = cm = plt.get_cmap('jet') 
    cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

    # Read output file sc
    if len(sys.argv[plot_or_save_arg:]) > 0:
        var_to_read = sys.argv[plot_or_save_arg:]
        if 'density_average' in sys.argv[plot_or_save_arg:]:
            if ('latitude' in var_to_read ) == False:
                var_to_read.append('latitude')
            if ('density' in var_to_read ) == False:
                var_to_read.append('density')

        if 'sma_average' in sys.argv[plot_or_save_arg:]:
            if ('latitude' in var_to_read ) == False:
                var_to_read.append('latitude')
            if ('sma' in var_to_read ) == False:
                var_to_read.append('sma')
    isc_count = -1
    for isc in range(nb_sc):
        isc_count = isc_count + 1
        isc_irun = isc_irun + 1
        var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
        if isc_count == 0: # same date for all sc
            date = var_out[find_in_read_input_order_variables(var_out_order, 'date')]
            orb_number = var_out[find_in_read_input_order_variables(var_out_order, 'orb_number')]
            nb_seconds_since_start = var_out[find_in_read_input_order_variables(var_out_order, 'nb_seconds_since_start')]
#            if len(date) < nb_steps: # this can happen if the sc reentered the atmosphere before the end of the run
            nb_steps_new = len(date) # in case the sc reentered the atmosphere before the end of the run
            date_ref = datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S.%f")
            start_xaxis_label = 0
            # test if date_start_plot and date_stop_plot are within the range of inital epoch and end epoch of the propagation
            if isc_irun == 0:
                if date_start_plot != 0:
                    if (datetime.strptime(date_start_plot, "%Y-%m-%dT%H:%M:%S.%f") - date_ref).total_seconds() < 0:
                        print "!*** The start date for the plot has to be earlier than the start epoch of the propagation. Therefore, it is set by default to the start epoch of the propagation ***!"
                        date_start_plot = 0
                if date_stop_plot != 0:
                    if (datetime.strptime(date_stop_plot, "%Y-%m-%dT%H:%M:%S.%f") - datetime.strptime(date[-1], "%Y/%m/%d %H:%M:%S.%f")).total_seconds() > 0:
                        print "!*** The stop date for the plot has to be older than the stop epoch of the propagation. Therefore, it is set by default to the stop epoch of the propagation ***!"
                        date_stop_plot = 0

                nb_seconds_in_simu = ( nb_steps_new - 1 ) * dt
                x_axis = np.arange(0, nb_seconds_in_simu + 1, dt) # all output files of one simulation have the same number of steps 

                if date_start_plot != 0:
                    date_ref = datetime.strptime(date_start_plot, "%Y-%m-%dT%H:%M:%S.%f")
                    delta_date_start = datetime.strptime(date_start_plot, "%Y-%m-%dT%H:%M:%S.%f") - datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S.%f")  
                    nb_seconds_between_epoch_start_and_date_start_plot = delta_date_start.days * 24 * 3600 + delta_date_start.seconds + delta_date_start.microseconds/10**6.
                    start_xaxis_label = nb_seconds_between_epoch_start_and_date_start_plot

                if date_stop_plot == 0:
                    if date_start_plot != 0:
                        nb_steps_between_epoch_start_and_date_start_plot = (int)(nb_seconds_between_epoch_start_and_date_start_plot/dt) + 1
                        nb_seconds_in_simu = ( nb_steps_new - nb_steps_between_epoch_start_and_date_start_plot ) * dt 
                    else:
                        nb_seconds_in_simu = ( nb_steps_new - 1 ) * dt
                else:
                    if date_start_plot == 0:
                        date_start_plot =  date_ref
                        delta_date = datetime.strptime(date_stop_plot, "%Y-%m-%dT%H:%M:%S.%f") - date_start_plot
                    else:
                        delta_date = datetime.strptime(date_stop_plot, "%Y-%m-%dT%H:%M:%S.%f") - datetime.strptime(date_start_plot, "%Y-%m-%dT%H:%M:%S.%f")
                    nb_seconds_in_simu = delta_date.days * 24 * 3600 + delta_date.seconds + delta_date.microseconds/10**6.
    #                 nb_steps_plot = (int)(nb_seconds_in_simu/dt) + 1
    #                 nb_seconds_in_simu = nb_steps_plot * dt


        if 'radius' in var_to_read:
            if isc_count == 0:
                radius = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            radius[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'radius')]
        if 'speed' in var_to_read:
            if isc_count == 0:
                speed = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            speed[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'speed')]
        if 'altitude' in var_to_read:
            if isc_count == 0:
                altitude = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            altitude[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'altitude')]
        if 'acceleration_lvlh_drag_mag' in var_to_read:
            if isc_count == 0:
                acceleration_lvlh_drag_mag = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            acceleration_lvlh_drag_mag[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'acceleration_lvlh_drag_mag')]

        if 'acceleration' in var_to_read:
            if isc_count == 0:
                acceleration = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            acceleration[isc, :nb_steps_new] = np.linalg.norm( var_out[find_in_read_input_order_variables(var_out_order, 'acceleration')], axis = 1)

        if 'phase_angle' in var_to_read:
            if isc_count == 0:
                phase_angle = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            phase_angle[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'phase_angle')]

        if 'latitude' in var_to_read:
            if isc_count == 0:
                latitude = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            latitude[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'latitude')]
        if 'longitude' in var_to_read:
            if isc_count == 0:
                longitude = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            longitude[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'longitude')]

        if 'eccentricity' in var_to_read:
            if isc_count == 0:
                eccentricity = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            eccentricity[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'eccentricity')]
        if 'solar_zenith' in var_to_read:
            if isc_count == 0:
                solar_zenith = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            solar_zenith[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'solar_zenith')]

        if 'sma' in var_to_read:
            if isc_count == 0:
                sma = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
                sma_ave = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            sma[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'sma')]
            sma_ave[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'sma_ave')]

        if 'sma_average' in var_to_read:
            if isc_count == 0:
                sma_average = []
                date_average_start_orbit_list = []
                date_average_end_orbit_list = []
                x_axis_average = []
            sma_orbit_averaged, time_averaged, index_time_averaged = orbit_average(sma[isc, :nb_steps_new], latitude[isc, :nb_steps_new], date )
            sma_average.append( sma_orbit_averaged ) # each sc might not have the same orbital period so the length of the array might not be the same between each sc
            date_average_start_orbit_list.append( np.array(time_averaged)[:,0] ) # take the date at the start of the bin
            date_average_end_orbit_list.append( np.array(time_averaged)[:,2] ) # take the date at the end of the bin
            x_axis_average_per_sc = []
            nb_orbit_for_this_sc = len(time_averaged)
            for iorbit in range(nb_orbit_for_this_sc):
                date_average_start_orbit = date_average_start_orbit_list[-1][iorbit]
                date_average_start_orbit = datetime.strptime( date_average_start_orbit, "%Y/%m/%d %H:%M:%S.%f" )
                nb_seconds_between_start_orbit_and_date_start = ( date_average_start_orbit - datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S.%f") ).total_seconds()
                date_average_end_orbit = date_average_end_orbit_list[-1][iorbit]
                date_average_end_orbit = datetime.strptime( date_average_end_orbit, "%Y/%m/%d %H:%M:%S.%f" )
                nb_seconds_between_end_orbit_and_date_start = ( date_average_end_orbit - datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S.%f") ).total_seconds()

#                 x_axis_average_per_sc.append( nb_seconds_between_start_orbit_and_date_start )
#                 x_axis_average_per_sc.append( nb_seconds_between_end_orbit_and_date_start )
                x_axis_average_per_sc.append((nb_seconds_between_end_orbit_and_date_start + nb_seconds_between_start_orbit_and_date_start)/2.)
            x_axis_average.append( x_axis_average_per_sc )

        if 'radius_perigee' in var_to_read:
            if isc_count == 0:
                radius_perigee = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            radius_perigee[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'radius_perigee')]
        if 'radius_apogee' in var_to_read:
            if isc_count == 0:
                radius_apogee = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            radius_apogee[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'radius_apogee')]
        if 'argument_perigee' in var_to_read:
            if isc_count == 0:
                argument_perigee = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
                argument_perigee_ave = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            argument_perigee[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'argument_perigee')]
            argument_perigee_ave[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'argument_perigee_ave')]
        if 'local_time' in var_to_read:
            if isc_count == 0:
                local_time = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps

            local_time[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'local_time')]

        if 'raan' in var_to_read:
            if isc_count == 0:
                raan = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            raan[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'raan')]

        if 'beta' in var_to_read:
            if isc_count == 0:
                beta = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            beta[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'beta')]

        if 'given_output' in var_to_read: # !!!!! first date of given output is not necesarrily the same as for the other variables  
            if isc_count == 0:
                given_output = []
            given_output.append(  var_out[find_in_read_input_order_variables(var_out_order, 'given_output')] )
        if 'power' in var_to_read:
            if isc_count == 0:
                power = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            power[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'power')]

        if 'density' in var_to_read:
            if isc_count == 0:
                density = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            density[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'density')]
        if 'density_average' in var_to_read:
            if isc_count == 0:
                density_average = []
                date_average_start_orbit_list = []
                date_average_end_orbit_list = []
                x_axis_average = []
            density_orbit_averaged, time_averaged, index_time_averaged = orbit_average(density[isc, :nb_steps_new], latitude[isc, :nb_steps_new], date )
            density_average.append( density_orbit_averaged ) # each sc might not have the same orbital period so the length of the array might not be the same between each sc
            date_average_start_orbit_list.append( np.array(time_averaged)[:,0] ) # take the date at the start of the bin
            date_average_end_orbit_list.append( np.array(time_averaged)[:,2] ) # take the date at the end of the bin
            x_axis_average_per_sc = []
            nb_orbit_for_this_sc = len(time_averaged)
            for iorbit in range(nb_orbit_for_this_sc):
                date_average_start_orbit = date_average_start_orbit_list[-1][iorbit]
                date_average_start_orbit = datetime.strptime( date_average_start_orbit, "%Y/%m/%d %H:%M:%S.%f" )
                nb_seconds_between_start_orbit_and_date_start = ( date_average_start_orbit - datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S.%f") ).total_seconds()
                x_axis_average_per_sc.append( nb_seconds_between_start_orbit_and_date_start )
            x_axis_average.append( x_axis_average_per_sc )

        if 'temperature' in var_to_read:
            if isc_count == 0:
                temperature = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            temperature[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'temperature')]
        if 'cd' in var_to_read:
            if isc_count == 0:
                cd = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            cd[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'cd')]
        if 'tot_area_drag' in var_to_read:
            if isc_count == 0:
                tot_area_drag = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
            tot_area_drag[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'tot_area_drag')]

        if 'position' in var_to_read:
            if isc_count == 0:
                position = np.zeros([nb_sc, nb_steps_new, 3]) # all output files of one simulation have the same number of steps
            position[isc, :nb_steps_new,:] = var_out[find_in_read_input_order_variables(var_out_order, 'position')]
        if 'velocity' in var_to_read:
            if isc_count == 0:
                velocity = np.zeros([nb_sc, nb_steps_new, 3]) # all output files of one simulation have the same number of steps
            velocity[isc, :nb_steps_new,:] = var_out[find_in_read_input_order_variables(var_out_order, 'velocity')]



######## PLOT
        if want_plot:
        # Radius
            if 'radius' in var_to_read:        
                if ( isc_count == 0 ) & (irun == 0) :
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

                colorVal = scalarMap.to_rgba(isc_irun)
                y_axis = radius[isc,:nb_steps_new] - earth_radius
                x_axis_temp = x_axis#phase_angle[isc,:nb_steps_new] #!!!!x_axis
                ax_radius.plot(x_axis_temp, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])
                #print np.mean(y_axis), np.mean(argument_perigee[isc,:nb_steps_new]), input_filename
                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8
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
                if ( isc_count == 0 ) & (irun == 0) :
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

                colorVal = scalarMap.to_rgba(isc_irun)
                y_axis = speed[isc,:nb_steps_new]
                ax_speed.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8
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
                if ( isc_count == 0 ) & (irun == 0) :
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

                colorVal = scalarMap.to_rgba(isc_irun)
                y_axis = altitude[isc,:nb_steps_new]
                ax_altitude.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8
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


            # Acceleration_Lvlh_Drag_Mag
            if 'acceleration_lvlh_drag_mag' in var_to_read:        
                if ( isc_count == 0 ) & (irun == 0) :
                    # Plot
                    fig_title = 'Acceleration_Lvlh_Drag_Mag as a function of time'
                    y_label = 'Acceleration_Lvlh_Drag_Mag (km)'
                    x_label = 'Real time'
                    fig_acceleration_lvlh_drag_mag = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                    fig_acceleration_lvlh_drag_mag.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                    gs = gridspec.GridSpec(1, 1)
                    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                    ax_acceleration_lvlh_drag_mag = fig_acceleration_lvlh_drag_mag.add_subplot(gs[0, 0])

                    ax_acceleration_lvlh_drag_mag.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                    ax_acceleration_lvlh_drag_mag.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                    [i.set_linewidth(2) for i in ax_acceleration_lvlh_drag_mag.spines.itervalues()] # change the width of the frame of the figure
                    ax_acceleration_lvlh_drag_mag.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

                colorVal = scalarMap.to_rgba(isc_irun)
                y_axis = acceleration_lvlh_drag_mag[isc,:nb_steps_new]
                ax_acceleration_lvlh_drag_mag.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8
                    dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                    xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
                    date_list_str = []
                    date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                    for i in range(len(xticks)):
                        if dt_xlabel >= 3*24*3600:
                            date_list_str.append( str(date_list[i])[5:10] )
                        else:
                            date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                    ax_acceleration_lvlh_drag_mag.xaxis.set_ticks(xticks)
                    ax_acceleration_lvlh_drag_mag.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                    ax_acceleration_lvlh_drag_mag.margins(0,0); ax_acceleration_lvlh_drag_mag.set_xlim([min(xticks), max(xticks)])
            #        ax_acceleration_lvlh_drag_mag.set_xlim([ax_acceleration_lvlh_drag_mag.get_xlim()[0], most_recent_tle_among_all_sc])

                    legend = ax_acceleration_lvlh_drag_mag.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                    legend.get_title().set_fontsize(str(fontsize_plot))

                    if save_plots == 1:
                        fig_save_name = 'acceleration_lvlh_drag_mag'
                        fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                        fig_acceleration_lvlh_drag_mag.savefig(fig_save_name, facecolor=fig_acceleration_lvlh_drag_mag.get_facecolor(), edgecolor='none', bbox_inches='tight')  

            # Acceleration
            if 'acceleration' in var_to_read:        
                if ( isc_count == 0 ) & (irun == 0) :
                    # Plot
                    fig_title = 'Acceleration as a function of time'
                    y_label = 'Acceleration (m/s$^2$)'
                    x_label = 'Longitude ' + u'(\N{DEGREE SIGN})' #'Real time'
                    fig_acceleration = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                    fig_acceleration.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                    gs = gridspec.GridSpec(1, 1)
                    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                    ax_acceleration = fig_acceleration.add_subplot(gs[0, 0])

                    ax_acceleration.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                    ax_acceleration.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                    [i.set_linewidth(2) for i in ax_acceleration.spines.itervalues()] # change the width of the frame of the figure
                    ax_acceleration.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

                colorVal = scalarMap.to_rgba(isc_irun)
                y_axis = acceleration[isc,:nb_steps_new] * 1000. # km/s2 to m/s2
                x_axis_temp = longitude[isc,:nb_steps_new]
                ax_acceleration.scatter(x_axis_temp, y_axis, s = 2, color = colorVal, label = label_arr[isc]) # !!!!!! soule be x_axis
                #ax_acceleration.plot(x_axis_temp, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc]) # !!!!!! soule be x_axis

                if isc == nb_sc - 1:
                    # x axis label is in real time
                    # ## all output files of one simulation have the same number of steps, and start at the same date
                    # nb_ticks_xlabel = 8
                    # dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                    # xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
                    # date_list_str = []
                    # date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                    # for i in range(len(xticks)):
                    #     if dt_xlabel >= 3*24*3600:
                    #         date_list_str.append( str(date_list[i])[5:10] )
                    #     else:
                    #         date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                    # ax_acceleration.xaxis.set_ticks(xticks)
                    # ax_acceleration.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                    # ax_acceleration.margins(0,0); ax_acceleration.set_xlim([min(xticks), max(xticks)])
            #        ax_acceleration.set_xlim([ax_acceleration.get_xlim()[0], most_recent_tle_among_all_sc])

                    legend = ax_acceleration.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                    legend.get_title().set_fontsize(str(fontsize_plot))

                    if save_plots == 1:
                        fig_save_name = 'acceleration'
                        fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                        fig_acceleration.savefig(fig_save_name, facecolor=fig_acceleration.get_facecolor(), edgecolor='none', bbox_inches='tight')  



                        
            # Latitude
            if 'latitude' in var_to_read:        
                if ( isc_count == 0 ) & (irun == 0) :
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

                colorVal = scalarMap.to_rgba(isc_irun)
                y_axis = latitude[isc,:nb_steps_new]
                ax_latitude.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8
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

            # # Longitude
            # if 'longitude' in var_to_read:        
            #     if ( isc_count == 0 ) & (irun == 0) :
            #         # Plot
            #         fig_title = 'Longitude as a function of time'
            #         y_label = 'Longitude ' + u'(\N{DEGREE SIGN})'
            #         x_label = 'Real time'
            #         fig_longitude = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

            #         fig_longitude.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
            #         plt.rc('font', weight='bold') ## make the labels of the ticks in bold
            #         gs = gridspec.GridSpec(1, 1)
            #         gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
            #         ax_longitude = fig_longitude.add_subplot(gs[0, 0])

            #         ax_longitude.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
            #         ax_longitude.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

            #         [i.set_linewidth(2) for i in ax_longitude.spines.itervalues()] # change the width of the frame of the figure
            #         ax_longitude.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
            #         plt.rc('font', weight='bold') ## make the labels of the ticks in bold

            #     colorVal = scalarMap.to_rgba(isc_irun)
            #     y_axis = longitude[isc,:nb_steps_new]
            #     ax_longitude.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

            #     if isc == nb_sc - 1:
            #         # x axis label is in real time
            #         ## all output files of one simulation have the same number of steps, and start at the same date
            #         nb_ticks_xlabel = 8
            #         dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
            #         xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
            #         date_list_str = []
            #         date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
            #         for i in range(len(xticks)):
            #             if dt_xlabel >= 3*24*3600:
            #                 date_list_str.append( str(date_list[i])[5:10] )
            #             else:
            #                 date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
            #         ax_longitude.xaxis.set_ticks(xticks)
            #         ax_longitude.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
            #         ax_longitude.margins(0,0); ax_longitude.set_xlim([min(xticks), max(xticks)])
            # #        ax_longitude.set_xlim([ax_longitude.get_xlim()[0], most_recent_tle_among_all_sc])

            #         legend = ax_longitude.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
            #         legend.get_title().set_fontsize(str(fontsize_plot))

            #         if save_plots == 1:
            #             fig_save_name = 'longitude'
            #             fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
            #             fig_longitude.savefig(fig_save_name, facecolor=fig_longitude.get_facecolor(), edgecolor='none', bbox_inches='tight')  

                        

            # Phase_Angle
            if 'phase_angle' in var_to_read:        
                if ( isc_count == 0 ) & (irun == 0) :
                    # Plot
                    fig_title = 'Phase_Angle as a function of time'
                    y_label = 'Phase_Angle ' + u'(\N{DEGREE SIGN})'
                    x_label = 'Real time'
                    fig_phase_angle = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                    fig_phase_angle.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                    gs = gridspec.GridSpec(1, 1)
                    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                    ax_phase_angle = fig_phase_angle.add_subplot(gs[0, 0])

                    ax_phase_angle.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                    ax_phase_angle.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                    [i.set_linewidth(2) for i in ax_phase_angle.spines.itervalues()] # change the width of the frame of the figure
                    ax_phase_angle.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

                colorVal = scalarMap.to_rgba(isc_irun)
                y_axis = phase_angle[isc,:nb_steps_new]
                ax_phase_angle.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8
                    dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                    xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
                    date_list_str = []
                    date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                    for i in range(len(xticks)):
                        if dt_xlabel >= 3*24*3600:
                            date_list_str.append( str(date_list[i])[5:10] )
                        else:
                            date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                    ax_phase_angle.xaxis.set_ticks(xticks)
                    ax_phase_angle.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                    ax_phase_angle.margins(0,0); ax_phase_angle.set_xlim([min(xticks), max(xticks)])
            #        ax_phase_angle.set_xlim([ax_phase_angle.get_xlim()[0], most_recent_tle_among_all_sc])

                    legend = ax_phase_angle.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                    legend.get_title().set_fontsize(str(fontsize_plot))

                    if save_plots == 1:
                        fig_save_name = 'phase_angle'
                        fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                        fig_phase_angle.savefig(fig_save_name, facecolor=fig_phase_angle.get_facecolor(), edgecolor='none', bbox_inches='tight')  


            # Eccentricity
            if 'eccentricity' in var_to_read:        
                if ( isc_count == 0 ) & (irun == 0) :
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

                colorVal = scalarMap.to_rgba(isc_irun)
                y_axis = eccentricity[isc,:nb_steps_new]
                ax_eccentricity.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8
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


            # Solar_Zenith
            if 'solar_zenith' in var_to_read:        
                if ( isc_count == 0 ) & (irun == 0) :
                    # Plot
                    fig_title = 'Solar zenith angle as a function of time'
                    y_label = 'Zenith angle '+ u'(\N{DEGREE SIGN})'
                    x_label = 'Real time'
                    fig_solar_zenith = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                    fig_solar_zenith.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                    gs = gridspec.GridSpec(1, 1)
                    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                    ax_solar_zenith = fig_solar_zenith.add_subplot(gs[0, 0])

                    ax_solar_zenith.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                    ax_solar_zenith.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                    [i.set_linewidth(2) for i in ax_solar_zenith.spines.itervalues()] # change the width of the frame of the figure
                    ax_solar_zenith.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

                colorVal = scalarMap.to_rgba(isc_irun)
                y_axis = solar_zenith[isc,:nb_steps_new]
                ax_solar_zenith.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8
                    dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                    xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
                    date_list_str = []
                    date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                    for i in range(len(xticks)):
                        if dt_xlabel >= 3*24*3600:
                            date_list_str.append( str(date_list[i])[5:10] )
                        else:
                            date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                    ax_solar_zenith.xaxis.set_ticks(xticks)
                    ax_solar_zenith.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                    ax_solar_zenith.margins(0,0); ax_solar_zenith.set_xlim([min(xticks), max(xticks)]); ax_solar_zenith.set_ylim([0,360])
            #        ax_solar_zenith.set_xlim([ax_solar_zenith.get_xlim()[0], most_recent_tle_among_all_sc])

                    legend = ax_solar_zenith.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                    legend.get_title().set_fontsize(str(fontsize_plot))

                    if save_plots == 1:
                        fig_save_name = 'solar_zenith'
                        fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                        fig_solar_zenith.savefig(fig_save_name, facecolor=fig_solar_zenith.get_facecolor(), edgecolor='none', bbox_inches='tight')  


            # SMA
            if 'sma' in var_to_read:        
                if ( isc_count == 0 ) & (irun == 0) :
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

                colorVal = scalarMap.to_rgba(isc_irun)
                y_axis = sma[isc,:nb_steps_new] -earth_radius
                if isc_count == 0:
                    min_y = np.min(y_axis)
                    max_y = np.max(y_axis)
                ax_sma.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])
                ax_sma.plot(x_axis, sma_ave[isc,:nb_steps_new]-earth_radius, linewidth = 2, color = colorVal, linestyle = 'dashed')
                if np.min(y_axis) < min_y:
                    min_y = np.min(y_axis)
                if np.max(y_axis) > max_y:
                    max_y = np.max(y_axis)

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
                    nb_ticks_xlabel = 8
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
                    ax_sma.margins(0,0); ax_sma.set_xlim([min(xticks), max(xticks)]); ax_sma.set_ylim([min_y, max_y])
            #        ax_sma.set_xlim([ax_sma.get_xlim()[0], most_recent_tle_among_all_sc])

                    legend = ax_sma.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                    legend.get_title().set_fontsize(str(fontsize_plot))

                    if save_plots == 1:
                        fig_save_name = 'sma'
                        fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                        fig_sma.savefig(fig_save_name, facecolor=fig_sma.get_facecolor(), edgecolor='none', bbox_inches='tight')  

            # Orbit average SMA
            if 'sma_average' in var_to_read:        
                if ( isc_count == 0 ) & (irun == 0) :
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

                colorVal = scalarMap.to_rgba(isc_irun)
                # y_axis = np.zeros(len(sma_average[isc])*2)
                # for iorbit in range(len(sma_average[isc])): # y_axis is constant along the orbit: it has the same value at the beginning and at the end of the orbit so it looks like a straight line during the orbit
                #     y_axis[iorbit*2] = sma_average[isc][iorbit]
                #     y_axis[iorbit*2+1] = sma_average[isc][iorbit]
                y_axis = np.zeros(len(sma_average[isc]))
                for iorbit in range(len(sma_average[isc])): # y_axis is constant along the orbit: it has the same value at the beginning and at the end of the orbit so it looks like a straight line during the orbit
                    y_axis[iorbit] = sma_average[isc][iorbit] - earth_radius

                ax_sma_average.plot(x_axis_average[isc], y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8
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




            if 'density_average' in var_to_read:        
                if ( isc_count == 0 ) & (irun == 0) :
                    # Plot
                    fig_title = ''#Orbit average density as a function of time'
                    y_label = 'Density (kg/m$^3$)'
                    x_label = 'Real time'
                    fig_density_average = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                    fig_density_average.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                    gs = gridspec.GridSpec(1, 1)
                    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                    ax_density_average = fig_density_average.add_subplot(gs[0, 0])

                    ax_density_average.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                    ax_density_average.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                    [i.set_linewidth(2) for i in ax_density_average.spines.itervalues()] # change the width of the frame of the figure
                    ax_density_average.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

                colorVal = scalarMap.to_rgba(isc_irun)
                y_axis = density_average[isc]
                ax_density_average.plot(x_axis_average[isc], y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8
                    dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                    xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
                    date_list_str = []
                    date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                    for i in range(len(xticks)):
                        if dt_xlabel >= 3*24*3600:
                            date_list_str.append( str(date_list[i])[5:10] )
                        else:
                            date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                    ax_density_average.xaxis.set_ticks(xticks)
                    ax_density_average.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                    ax_density_average.margins(0,0); ax_density_average.set_xlim([min(xticks), max(xticks)])
            #        ax_density_average.set_xlim([ax_density_average.get_xlim()[0], most_recent_tle_among_all_sc])

                    legend = ax_density_average.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Quantile", fontsize = fontsize_plot)
                    legend.get_title().set_fontsize(str(fontsize_plot))

                    if save_plots == 1:
                        fig_save_name = 'density_average'
                        fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                        fig_density_average.savefig(fig_save_name, facecolor=fig_density_average.get_facecolor(), edgecolor='none', bbox_inches='tight')  






            # Radius of perigee
            if 'radius_perigee' in var_to_read:        
                if ( isc_count == 0 ) & (irun == 0) :
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

                colorVal = scalarMap.to_rgba(isc_irun)
                y_axis = radius_perigee[isc,:nb_steps_new] - earth_radius
                ax_radius_perigee.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8
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
                if ( isc_count == 0 ) & (irun == 0) :
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

                colorVal = scalarMap.to_rgba(isc_irun)
                y_axis = radius_apogee[isc,:nb_steps_new] - earth_radius
                ax_radius_apogee.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])

                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8
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
                if ( isc_count == 0 ) & (irun == 0) :
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

                colorVal = scalarMap.to_rgba(isc_irun)
                y_axis = argument_perigee[isc,:nb_steps_new] 
                ax_argument_perigee.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])
                ax_argument_perigee.plot(x_axis, argument_perigee_ave[isc,:nb_steps_new], linewidth = 2, color = colorVal,  linestyle = 'dashed')

                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8
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



            # Argument of perigee
            if 'local_time' in var_to_read:        
                if ( isc_count == 0 ) & (irun == 0) :
                    # Plot
                    fig_title = 'Local time as a function of time'
                    y_label = 'Local time '+ u'(\N{DEGREE SIGN})'
                    x_label = 'Real time'
                    fig_local_time = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                    fig_local_time.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                    gs = gridspec.GridSpec(1, 1)
                    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                    ax_local_time = fig_local_time.add_subplot(gs[0, 0])

                    ax_local_time.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                    ax_local_time.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                    [i.set_linewidth(2) for i in ax_local_time.spines.itervalues()] # change the width of the frame of the figure
                    ax_local_time.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

                colorVal = scalarMap.to_rgba(isc_irun)
                y_axis = local_time[isc,:nb_steps_new] 
                ax_local_time.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc])


                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8
                    dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                    xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
                    date_list_str = []
                    date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                    for i in range(len(xticks)):
                        if dt_xlabel >= 3*24*3600:
                            date_list_str.append( str(date_list[i])[5:10] )
                        else:
                            date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                    ax_local_time.xaxis.set_ticks(xticks)
                    ax_local_time.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                    ax_local_time.margins(0,0); ax_local_time.set_xlim([min(xticks), max(xticks)]);ax_local_time.set_ylim([0,360])
            #        ax_local_time.set_xlim([ax_local_time.get_xlim()[0], most_recent_tle_among_all_sc])

                    legend = ax_local_time.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                    legend.get_title().set_fontsize(str(fontsize_plot))

                    if save_plots == 1:
                        fig_save_name = 'local_time'
                        fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                        fig_local_time.savefig(fig_save_name, facecolor=fig_local_time.get_facecolor(), edgecolor='none', bbox_inches='tight')  


            # RAAN
            if 'raan' in var_to_read:        
                if ( isc_count == 0 ) & (irun == 0) :
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

                colorVal = scalarMap.to_rgba(isc_irun)
                y_axis = raan[isc,:nb_steps_new] 
                ax_raan.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc] )

                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8
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

            # BETA
            if 'beta' in var_to_read:        
                if ( isc_count == 0 ) & (irun == 0) :
                    # Plot
                    fig_title = ''#Beta angle as a function of time'
                    y_label = 'Beta angle '+ u'(\N{DEGREE SIGN})'
                    x_label = 'Real time'
                    fig_beta = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                    fig_beta.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
                    plt.rc('font', weight='normal') ## make the labels of the ticks in normal
                    gs = gridspec.GridSpec(1, 1)
                    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                    ax_beta = fig_beta.add_subplot(gs[0, 0])

                    ax_beta.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
                    ax_beta.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

                    [i.set_linewidth(2) for i in ax_beta.spines.itervalues()] # change the width of the frame of the figure
                    ax_beta.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                    plt.rc('font', weight='normal') ## make the labels of the ticks in normal

                colorVal = scalarMap.to_rgba(isc_irun)
                y_axis = beta[isc,:nb_steps_new] 
                ax_beta.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc] )
                ax_beta.plot([x_axis[0], x_axis[-1]], [42,42], linewidth = 2, color = 'r', linestyle = 'dashed')
                ax_beta.plot([x_axis[0], x_axis[-1]], [-42,-42], linewidth = 2, color = 'r', linestyle = 'dashed')
#                 ax_beta.text(x_axis[0], 42, 'above, roll by -22' + u'\N{DEGREE SIGN}', fontsize = fontsize_plot, color = 'r', verticalalignment = 'bottom', horizontalalignment = 'left' )
#                 ax_beta.text(x_axis[0], -42.2, 'below, roll by +22' + u'\N{DEGREE SIGN}', fontsize = fontsize_plot, color = 'r', verticalalignment = 'top', horizontalalignment = 'left' )
                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8
                    dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                    xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel) 
                    date_list_str = []
                    date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                    for i in range(len(xticks)):
                        if dt_xlabel >= 3*24*3600:
                            date_list_str.append( str(date_list[i])[5:10] )
                        else:
                            date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                    ax_beta.xaxis.set_ticks(xticks)
                    ax_beta.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                    ax_beta.margins(0,0); ax_beta.set_xlim([min(xticks), max(xticks)]); ax_beta.set_ylim([-50,50])
            #        ax_beta.set_xlim([ax_beta.get_xlim()[0], most_recent_tle_among_all_sc])

#                     legend = ax_beta.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
#                     legend.get_title().set_fontsize(str(fontsize_plot))

                    if save_plots == 1:
                        fig_save_name = 'beta'
                        fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                        fig_beta.savefig(fig_save_name, facecolor=fig_beta.get_facecolor(), edgecolor='none', bbox_inches='tight')  




            # POWER
            if 'power' in var_to_read:         # !!!!! first date of given output is not necesarrily the same as for the other variables  
                if ( isc_count == 0 ) & (irun == 0) :
                    # Plot
                    fig_title = ' as a function of time'
                    y_label = ''
                    x_label = 'Real time'
                    fig_power = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                    fig_power.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                    gs = gridspec.GridSpec(1, 1)
                    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                    ax_power = fig_power.add_subplot(gs[0, 0])

                    ax_power.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                    ax_power.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                    [i.set_linewidth(2) for i in ax_power.spines.itervalues()] # change the width of the frame of the figure
                    ax_power.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

                colorVal = scalarMap.to_rgba(isc_irun)
                y_axis = power[isc]
                ax_power.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc] )

                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8
                    dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                    xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel) 
                    date_list_str = []
                    date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                    for i in range(len(xticks)):
                        if dt_xlabel >= 3*24*3600:
                            date_list_str.append( str(date_list[i])[5:10] )
                        else:
                            date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                    ax_power.xaxis.set_ticks(xticks)
                    ax_power.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                    ax_power.margins(0,0); ax_power.set_xlim([min(xticks), max(xticks)])
            #        ax_power.set_xlim([ax_power.get_xlim()[0], most_recent_tle_among_all_sc])

                    legend = ax_power.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                    legend.get_title().set_fontsize(str(fontsize_plot))

                    if save_plots == 1:
                        fig_save_name = 'power'
                        fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                        fig_power.savefig(fig_save_name, facecolor=fig_power.get_facecolor(), edgecolor='none', bbox_inches='tight')  


            # GIVEN_OUTPUT
            if 'given_output' in var_to_read:         # !!!!! first date of given output is not necesarrily the same as for the other variables  
                if ( isc_count == 0 ) & (irun == 0) :
                    # Plot
                    fig_title = ' as a function of time'
                    y_label = ''
                    x_label = 'Real time'
                    fig_given_output = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                    fig_given_output.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                    gs = gridspec.GridSpec(1, 1)
                    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                    ax_given_output = fig_given_output.add_subplot(gs[0, 0])

                    ax_given_output.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                    ax_given_output.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                    [i.set_linewidth(2) for i in ax_given_output.spines.itervalues()] # change the width of the frame of the figure
                    ax_given_output.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

                colorVal = scalarMap.to_rgba(isc_irun)
                given_output_sc = np.array(given_output[isc])
                y_axis = given_output_sc 
                ax_given_output.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc] )

                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8
                    dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                    xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel) 
                    date_list_str = []
                    date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                    for i in range(len(xticks)):
                        if dt_xlabel >= 3*24*3600:
                            date_list_str.append( str(date_list[i])[5:10] )
                        else:
                            date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                    ax_given_output.xaxis.set_ticks(xticks)
                    ax_given_output.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                    ax_given_output.margins(0,0); ax_given_output.set_xlim([min(xticks), max(xticks)])
            #        ax_given_output.set_xlim([ax_given_output.get_xlim()[0], most_recent_tle_among_all_sc])

                    legend = ax_given_output.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                    legend.get_title().set_fontsize(str(fontsize_plot))

                    if save_plots == 1:
                        fig_save_name = 'given_output'
                        fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                        fig_given_output.savefig(fig_save_name, facecolor=fig_given_output.get_facecolor(), edgecolor='none', bbox_inches='tight')  


            # DENSITY
            if 'density' in var_to_read:         # !!!!! first date of given output is not necesarrily the same as for the other variables  
                if ( isc_count == 0 ) & (irun == 0) :
                    # Plot
                    fig_title ='' #'Density as a function of time'
                    y_label = 'Density (kg/m$^3$)'#'rho_control' #'Density (kg/m$^3$)'
                    x_label = 'Solar zenith angle ' + u'(\N{DEGREE SIGN})'#Local time  'Phase angle ' + u'(\N{DEGREE SIGN})'#Real time'
                    fig_density = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                    fig_density.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                    gs = gridspec.GridSpec(1, 1)
                    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                    ax_density = fig_density.add_subplot(gs[0, 0])

                    ax_density.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                    ax_density.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                    [i.set_linewidth(2) for i in ax_density.spines.itervalues()] # change the width of the frame of the figure
                    ax_density.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

                colorVal = scalarMap.to_rgba(isc_irun)

                nb_orb = orb_number[-1][1] # nb of orbits traveled by sc during simu

                for iorb in range(nb_orb):
                    if iorb == 0:
                        iorb_start = 0
                    else:
                        iorb_start_previous = iorb_start
                        iorb_start = orb_number[iorb-1][0] + 1
                    if iorb > 0:
                        iorb_end_previous = iorb_end
                    iorb_end = orb_number[iorb][0]+1
                    y_axis = density[isc, iorb_start:iorb_end ]  #density[isc]#density[isc]*1e9 - 1
                    x_axis_density = solar_zenith[isc, iorb_start:iorb_end]#local_time[isc, iorb_start:iorb_end]#phase_angle[isc, iorb_start:iorb_end]#argument_perigee[isc,:nb_steps_new] #latitude[isc,:nb_steps_new]
#                     if iorb == 0:
#                         ax_density.scatter(x_axis_density, y_axis, s = 5, color = colorVal, label = label_arr[isc] )
#                     else:
                    if iorb > 0: # ignore first orbit because for density control the desnity is equal to msis for first orbit -> not interested in plotting for the current analysis
                        ax_density.scatter(x_axis_density, y_axis, s = 5, color = colorVal, label = label_arr[isc] )
#                     # If looking whwere the perigee is, color where the phase angle is close to perigee or apogee
#                     if iorb >= 1:
#                         where_arg_per = np.where(np.abs(x_axis_density - argument_perigee_ave[isc,iorb_start_previous]) < 2.5)[0] # where the phase angle is les than 2.5 deg off the arg per
#                         ax_density.scatter(x_axis_density[where_arg_per], y_axis[where_arg_per], s = 15, color = 'r' )

#                         where_arg_apo = np.where(np.abs(x_axis_density - (np.mod(argument_perigee_ave[isc,iorb_start_previous]+180, 360))) < 2.5)[0] # where the phase angle is les than 2.5 deg off the arg per
#                         ax_density.scatter(x_axis_density[where_arg_apo], y_axis[where_arg_apo], s = 15, color = 'limegreen' )

#                         print iorb_start, argument_perigee_ave[isc,iorb_start], argument_perigee_ave[isc,iorb_start_previous]

#                 #                ax_density.plot(argument_perigee_ave[isc,:nb_steps_new], np.arange(), s = 5, color = colorVal, label = label_arr[isc] )

                    if ((irun == 0) & (isc_count == 0) & (iorb == 0)):
                        min_y_density = np.min(y_axis)
                        max_y_density = np.max(y_axis)

                    if np.min(y_axis) < min_y_density:
                        min_y_density = np.min(y_axis)
                    if np.max(y_axis) > max_y_density:
                        max_y_density = np.max(y_axis)
                
                if isc == nb_sc - 1:
                    # x axis label is in real time
                    # ## all output files of one simulation have the same number of steps, and start at the same date
                    # nb_ticks_xlabel = 8
                    # dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                    # xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel) 
                    # date_list_str = []
                    # date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                    # for i in range(len(xticks)):
                    #     if dt_xlabel >= 3*24*3600:
                    #         date_list_str.append( str(date_list[i])[5:10] )
                    #     else:
                    #         date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                    # ax_density.xaxis.set_ticks(xticks)
                    # ax_density.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                    # ax_density.set_xlim([min(xticks), max(xticks)])
                    ax_density.margins(0,0); ax_density.set_ylim([min_y_density*(1-0.1*np.sign(min_y_density)), max_y_density*1.1])
                    #        ax_density.set_xlim([ax_density.get_xlim()[0], most_recent_tle_among_all_sc])
                    
                    legend = ax_density.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="", fontsize = fontsize_plot)
                    legend.get_title().set_fontsize(str(fontsize_plot))
                    
                    for handle in legend.legendHandles:
                        handle.set_sizes([100.0])
#                     legend.legendHandles[0]._sizes = [100,100,100]
                    #legend.legendHandles[1]._sizes = [100]

                    if save_plots == 1:
                        fig_save_name = 'density'
                        fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                        fig_density.savefig(fig_save_name, facecolor=fig_density.get_facecolor(), edgecolor='none', bbox_inches='tight')  


           # TEMPERATURE
            if 'temperature' in var_to_read:         # !!!!! first date of given output is not necesarrily the same as for the other variables  
                if ( isc_count == 0 ) & (irun == 0) :
                    # Plot
                    fig_title = 'Temperature as a function of time'
                    y_label = 'Temperature (K)'
                    x_label = 'Real time'
                    fig_temperature = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                    fig_temperature.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                    gs = gridspec.GridSpec(1, 1)
                    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                    ax_temperature = fig_temperature.add_subplot(gs[0, 0])

                    ax_temperature.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                    ax_temperature.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                    [i.set_linewidth(2) for i in ax_temperature.spines.itervalues()] # change the width of the frame of the figure
                    ax_temperature.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

                colorVal = scalarMap.to_rgba(isc_irun)
                y_axis = temperature[isc]
                ax_temperature.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc] )

                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8
                    dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                    xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel) 
                    date_list_str = []
                    date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                    for i in range(len(xticks)):
                        if dt_xlabel >= 3*24*3600:
                            date_list_str.append( str(date_list[i])[5:10] )
                        else:
                            date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                    ax_temperature.xaxis.set_ticks(xticks)
                    ax_temperature.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                    ax_temperature.margins(0,0); ax_temperature.set_xlim([min(xticks), max(xticks)])
            #        ax_temperature.set_xlim([ax_temperature.get_xlim()[0], most_recent_tle_among_all_sc])

                    legend = ax_temperature.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                    legend.get_title().set_fontsize(str(fontsize_plot))

                    if save_plots == 1:
                        fig_save_name = 'temperature'
                        fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                        fig_temperature.savefig(fig_save_name, facecolor=fig_temperature.get_facecolor(), edgecolor='none', bbox_inches='tight')  

           # CD
            if 'cd' in var_to_read:         # !!!!! first date of given output is not necesarrily the same as for the other variables  
                if ( isc_count == 0 ) & (irun == 0) : 
                    # Plot
                    fig_title = 'Total normalized Cd as a function of time'
                    y_label = 'Cd'
                    x_label = 'Real time'
                    fig_cd = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                    fig_cd.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                    gs = gridspec.GridSpec(1, 1)
                    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                    ax_cd = fig_cd.add_subplot(gs[0, 0])

                    ax_cd.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                    ax_cd.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                    [i.set_linewidth(2) for i in ax_cd.spines.itervalues()] # change the width of the frame of the figure
                    ax_cd.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

                colorVal = scalarMap.to_rgba(isc_irun)        
                y_axis = cd[isc]
                ax_cd.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc] )

                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8
                    dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                    xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel) 
                    date_list_str = []
                    date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                    for i in range(len(xticks)):
                        if dt_xlabel >= 3*24*3600:
                            date_list_str.append( str(date_list[i])[5:10] )
                        else:
                            date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                    ax_cd.xaxis.set_ticks(xticks)
                    ax_cd.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                    ax_cd.margins(0,0); ax_cd.set_xlim([min(xticks), max(xticks)])
            #        ax_cd.set_xlim([ax_cd.get_xlim()[0], most_recent_tle_among_all_sc])

                    legend = ax_cd.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                    legend.get_title().set_fontsize(str(fontsize_plot))

                    if save_plots == 1:
                        fig_save_name = 'cd'
                        fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                        fig_cd.savefig(fig_save_name, facecolor=fig_cd.get_facecolor(), edgecolor='none', bbox_inches='tight')  


           # TOT_AREA_DRAG
            if 'tot_area_drag' in var_to_read:         # !!!!! first date of given output is not necesarrily the same as for the other variables  
                if ( isc_count == 0 ) & (irun == 0) : 
                    # Plot
                    fig_title = 'Total drag area as a function of time'
                    y_label = 'A (cm$^2$)'
                    x_label = 'Real time'
                    fig_tot_area_drag = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

                    fig_tot_area_drag.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                    gs = gridspec.GridSpec(1, 1)
                    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
                    ax_tot_area_drag = fig_tot_area_drag.add_subplot(gs[0, 0])

                    ax_tot_area_drag.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                    ax_tot_area_drag.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                    [i.set_linewidth(2) for i in ax_tot_area_drag.spines.itervalues()] # change the width of the frame of the figure
                    ax_tot_area_drag.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

                colorVal = scalarMap.to_rgba(isc_irun)        
                y_axis = tot_area_drag[isc]
                ax_tot_area_drag.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = label_arr[isc] )

                if isc == nb_sc - 1:
                    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
                    nb_ticks_xlabel = 8        
                    dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                    xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel) 
                    date_list_str = []
                    date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
                    for i in range(len(xticks)):
                        if dt_xlabel >= 3*24*3600:
                            date_list_str.append( str(date_list[i])[5:10] )
                        else:
                            date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                    ax_tot_area_drag.xaxis.set_ticks(xticks)
                    ax_tot_area_drag.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                    ax_tot_area_drag.margins(0,0); ax_tot_area_drag.set_xlim([min(xticks), max(xticks)])
            #        ax_tot_area_drag.set_xlim([ax_tot_area_drag.get_xlim()[0], most_recent_tle_among_all_sc])

                    legend = ax_tot_area_drag.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                    legend.get_title().set_fontsize(str(fontsize_plot))

                    if save_plots == 1:
                        fig_save_name = 'tot_area_drag'
                        fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                        fig_tot_area_drag.savefig(fig_save_name, facecolor=fig_tot_area_drag.get_facecolor(), edgecolor='none', bbox_inches='tight')  



######## SAVE
        if want_save:
            if isc == nb_sc - 1:
                pickle.dump( date, open( root_save_fig_name + 'date' + ".pickle", "w" ) )
                pickle.dump( nb_seconds_since_start, open( root_save_fig_name + 'nb_seconds_since_start' + ".pickle", "w" ) )
                # Radius
                if 'radius' in var_to_read:
                    pickle.dump( radius, open( root_save_fig_name + 'radius' + ".pickle", "w" ) )
                # Speed
                if 'speed' in var_to_read:
                    pickle.dump( speed, open( root_save_fig_name + 'speed' + ".pickle", "w" ) )
                # Altitude
                if 'altitude' in var_to_read:
                    pickle.dump( altitude, open( root_save_fig_name + 'altitude' + ".pickle", "w" ) )
                # Acceleration_Lvlh_Drag_Mag
                if 'acceleration_lvlh_drag_mag' in var_to_read:
                    pickle.dump( acceleration_lvlh_drag_mag, open( root_save_fig_name + 'acceleration_lvlh_drag_mag' + ".pickle", "w" ) )

                # Acceleration
                if 'acceleration' in var_to_read:
                    pickle.dump( acceleration, open( root_save_fig_name + 'acceleration' + ".pickle", "w" ) )

                # Latitude
                if 'latitude' in var_to_read:
                    pickle.dump( latitude, open( root_save_fig_name + 'latitude' + ".pickle", "w" ) )
                # Longitude
                if 'longitude' in var_to_read:
                    pickle.dump( longitude, open( root_save_fig_name + 'longitude' + ".pickle", "w" ) )

                # Phase_Angle
                if 'phase_angle' in var_to_read:
                    pickle.dump( phase_angle, open( root_save_fig_name + 'phase_angle' + ".pickle", "w" ) )

                # Eccentricity
                if 'eccentricity' in var_to_read:
                    pickle.dump( eccentricity, open( root_save_fig_name + 'eccentricity' + ".pickle", "w" ) )
                # Solar_Zenith
                if 'solar_zenith' in var_to_read:
                    pickle.dump( solar_zenith, open( root_save_fig_name + 'solar_zenith' + ".pickle", "w" ) )

                # Sma
                if 'sma' in var_to_read:
                    pickle.dump( sma, open( root_save_fig_name + 'sma' + ".pickle", "w" ) )
                # Orbit average SMA
                if 'sma_average' in var_to_read:
                    pickle.dump( sma_average, open( root_save_fig_name + 'sma_average' + ".pickle", "w" ) )
                    pickle.dump( x_axis_average, open( root_save_fig_name + 'x_axis_average' + ".pickle", "w" ) )
                if 'density_average' in var_to_read:
                    pickle.dump( density_average, open( root_save_fig_name + 'density_average' + ".pickle", "w" ) )
                    pickle.dump( x_axis_average, open( root_save_fig_name + 'x_axis_average' + ".pickle", "w" ) )

                # Radius of perigee
                if 'radius_perigee' in var_to_read:
                    pickle.dump( radius_perigee, open( root_save_fig_name + 'radius_perigee' + ".pickle", "w" ) )
                # Radius of apogee
                if 'radius_apogee' in var_to_read:
                    pickle.dump( radius_apogee, open( root_save_fig_name + 'radius_apogee' + ".pickle", "w" ) )
                # Argument of perigee
                if 'argument_perigee' in var_to_read:
                    pickle.dump( argument_perigee, open( root_save_fig_name + 'argument_perigee' + ".pickle", "w" ) )
                # RAAN
                if 'raan' in var_to_read:
                    pickle.dump( raan, open( root_save_fig_name + 'raan' + ".pickle", "w" )  )
                # BETA
                if 'beta' in var_to_read:
                    pickle.dump( beta, open( root_save_fig_name + 'beta' + ".pickle", "w" )  )


                if 'given_output' in var_to_read:
                    pickle.dump( given_output, open( root_save_fig_name + 'given_output' + ".pickle", "w" )  )
                # TEMPERATURE
                if 'temperature' in var_to_read:
                    pickle.dump( temperature, open( root_save_fig_name + 'temperature' + ".pickle", "w" )  )
                # CD
                if 'cd' in var_to_read:
                    pickle.dump( cd, open( root_save_fig_name + 'cd' + ".pickle", "w" )  )
                # TOT_AREA_DRAG
                if 'tot_area_drag' in var_to_read:
                    pickle.dump( tot_area_drag, open( root_save_fig_name + 'tot_area_drag' + ".pickle", "w" )  )

                # Position
                if 'position' in var_to_read:
                    pickle.dump( position, open( root_save_fig_name + 'position' + ".pickle", "w" )  )
                # Velocity
                if 'velocity' in var_to_read:
                    pickle.dump( velocity, open( root_save_fig_name + 'velocity' + ".pickle", "w" )  )


if ( ( show_plots == 1 ) & ( want_plot == 1 ) ):
    plt.show(); plt.show();
