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
# This script:
# 1- runs cygnss_spock_tle_run_automatically.py to determine the orbit average SMA from TLE observations (see header of cygnss_spock_tle.py for more details) from the CYGNSS launch until today
# 2- predicts the future states of the CYGNSS sc from today until today + 6 months (this number of months is a parameter so it can be easily changed)
# ASSUMPTIONS:
# - see section 'parameters to set before running this script' below before running this script (in particular you also need to change parameters in spock_cygnss_spec_set_spock_simu_parameters.py (see comments in 'parameters to set before running this script') of this script)           
import pickle
from read_input_file import *
from read_output_file import *
from spock_main_input import *
from orbit_average import *
from convert_tle_date_to_date import *
from norad_id_to_cygnss_name import *
import matplotlib.ticker as mtick
import matplotlib.colors as colors
import matplotlib.cm as cmx
from get_prop_dir import *
import matplotlib.gridspec as gridspec
from read_input_file import *
import pickle
import sys
import fileinput
import time
import numpy as np
from matplotlib import pyplot as plt
import os
import subprocess
from get_name_mission import *
from find_in_read_input_order_variables import *
from datetime import datetime, timedelta
from cygnss_spock_tle_run_automatically import *
from cygnss_decay_run_automatically import *

######### PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT #########
mpi_path = 'mpirun-openmpi-gcc49' # path of mpi run for the machine running this script. BE CAREFUL: if you change it here, you also need to change it in spock_cygnss_spec_set_spock_simu_parameters.py
run_dir = "run.cygnss" # SpOCK run directory. BE CAREFUL: if you change it here, you also need to change it in spock_cygnss_spec_set_spock_simu_parameters.py
N_min_propagation = 200 # Number of minutes to propagate each past TLEs in cygnss_spock_tle_run_automatically.py  !!!!!!!! go back to 200
######### end of PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT #########

## Parameters for the figure
root_save_fig_name = './cygnss/' # folder where figures are saved
save_plots = 1 # set to 1 to save the results
show_plots = 0 # set to 1 to show_plot

height_fig = 9.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3
nb_ticks_xlabel = 6.

######### ALGORITHM #########
earth_mu    = 398600.4418; # gravitational parameter (km^3/s^2)
earth_radius        = 6378.137; # mean equatorial radius (km)

if show_plots == 1:
    plt.ion()


# Runs cygnss_spock_tle_run_automatically.py to determine the orbit average SMA from TLE observations (see header of cygnss_spock_tle.py for more details) from the CYGNSS launch until today
print "Computing historical SMA average..."
date_ini = "2016-12-20T00:00:00" # if you change that, you also need to change it when downloading the TLEs in cygnss_spock_tle_run_automatically
date_ini = datetime.strptime(date_ini, "%Y-%m-%dT%H:%M:%S") 
main_input_filename_list = cygnss_spock_tle_run_automatically(mpi_path, run_dir, N_min_propagation) 
# pickle.dump( main_input_filename_list , open( "main_input_filename_list.pickle", "w" ) )
# main_input_filename_list = pickle.load( open( "main_input_filename_list.pickle", "r" ) )
# main_input_filename_list is the list of all the main input file used to run SpOCK in cygnss_spock_tle_run_automatically.py
# end of runs cygnss_spock_tle_run_automatically.py to determine the orbit average SMA from TLE observations


# Reads the output of cygnss_spock_tle_run_automatically.py
ecc = []
sma_average_historical = []
x_axis_average_historical = []
radius_apogee = []
radius_perigee = []
dsmadday = []
alt_perigee = []
alt_apogee = []
date_tle = []
sc_name = []
nb_sc = len(main_input_filename_list)
for isc in range(nb_sc):
    date_tle_per_sc = []
    sc_name.append( main_input_filename_list[isc][0][0:7] )
    sc_name_temp = sc_name[-1]
    alt_apogee_per_sc = []
    alt_perigee_per_sc = []
    sma_average_historical_per_sc = []
    radius_apogee_per_sc = []
    radius_perigee_per_sc = []
    x_axis_average_historical_per_sc = []
    nb_tle_per_sc = len(main_input_filename_list[isc])
    for itle in range(nb_tle_per_sc): # !!!!!!!!!! go back to range(nb_tle_per_sc)
        ## Reads the output of the N_min_propagation propagations performed in cygnss_spock_tle_run_automatically.py
        main_input_filename = main_input_filename_list[isc][itle]
        var_in, var_in_order = read_input_file(get_prop_dir(1) + run_dir + "/input/main_input/" + main_input_filename)
        date_start = var_in[find_in_read_input_order_variables(var_in_order, 'date_start')]; 
        date_tle_per_sc.append( date_start )
        output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
        output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
        ### Read SpOCK's outputs of the N_min_propagation propagation
        isc_current = 0 # here only one sc in main input file
        var_to_read = ["altitude", "sma", "radius", "latitude"]
        var_out, var_out_order = read_output_file( output_file_path_list[isc_current] + output_file_name_list[isc_current], var_to_read )
        date_temp = var_out[find_in_read_input_order_variables(var_out_order, 'date')]
        altitude = var_out[find_in_read_input_order_variables(var_out_order, 'altitude')]
        latitude_temp = var_out[find_in_read_input_order_variables(var_out_order, 'latitude')]
        sma_temp = var_out[find_in_read_input_order_variables(var_out_order, 'sma')]
        sma_orbit_averaged, time_averaged, index_time_averaged = orbit_average(sma_temp, latitude_temp, date_temp ) # !!!!!!!! uncomment

        date_average_start_orbit_list = np.array(time_averaged)[:,0]  # take the date at the start of the bin
        date_average_end_orbit_list = np.array(time_averaged)[:,2]  # take the date at the end of the bin
        nb_orbit_for_this_sc_for_this_tle_run = len(time_averaged)
        for iorbit in range(nb_orbit_for_this_sc_for_this_tle_run):
            date_average_start_orbit = date_average_start_orbit_list[iorbit]
            date_average_start_orbit = datetime.strptime( date_average_start_orbit, "%Y/%m/%d %H:%M:%S" )
            nb_seconds_between_start_orbit_and_date_start = ( date_average_start_orbit - date_ini ).total_seconds()
            date_average_end_orbit = date_average_end_orbit_list[iorbit]
            date_average_end_orbit = datetime.strptime( date_average_end_orbit, "%Y/%m/%d %H:%M:%S" )
            nb_seconds_between_end_orbit_and_date_start = ( date_average_end_orbit - date_ini ).total_seconds()

            x_axis_average_historical_per_sc.append( nb_seconds_between_start_orbit_and_date_start )
            x_axis_average_historical_per_sc.append( nb_seconds_between_end_orbit_and_date_start )
            sma_average_historical_per_sc.append( sma_orbit_averaged[iorbit] ) 
            

        radius_temp = var_out[find_in_read_input_order_variables(var_out_order, 'radius')]
        alt_perigee_per_sc.append( np.min(altitude) )
        alt_apogee_per_sc.append( np.max(altitude) )
        radius_apogee_per_sc.append( np.max( radius_temp ))
        radius_perigee_per_sc.append( np.min( radius_temp ) )
    date_tle.append( date_tle_per_sc ) # note that this is actually the date of the TLE + one time step (see cygnss_spock_tle_run_automatically.py)
    alt_apogee.append( alt_apogee_per_sc )
    alt_perigee.append( alt_perigee_per_sc )
    sma_average_historical.append(sma_average_historical_per_sc)
    x_axis_average_historical.append(x_axis_average_historical_per_sc)
    radius_apogee.append(radius_apogee_per_sc)
    radius_perigee.append(radius_perigee_per_sc)
# end of reads the output of cygnss_spock_tle_run_automatically.py


# Predicts the future states of the CYGNSS sc from today until today + 6 months (this number of months is a parameter so it can be easily changed)
print "Predicting the future states of the CYGNSS satellites from today until today + 6 months..."
## Run cygnss_decay_run_automatically.py to generate SpOCK's outputs
main_input_filename_nadir_prediction = cygnss_decay_run_automatically()
main_input_filename_high_drag_prediction = main_input_filename_nadir_prediction.replace("nadir", "high_drag")
## Read SpOCK's outputs
### Run state.py to save SpOCK's outputs as pickles
os.system("python state.py " + run_dir + " " + main_input_filename_nadir_prediction + " save sma_average")
os.system("python state.py " + run_dir + " " + main_input_filename_high_drag_prediction + " save sma_average") # !!!!!!!!! uncomment
### Load pickles created by state.py
#### Nadir
pickle_sma_average_nadir_prediction = get_prop_dir(1) + run_dir + '/output/python_out/' + main_input_filename_nadir_prediction.replace(".txt","_sma_average.pickle")
sma_average_nadir_prediction = pickle.load( open( pickle_sma_average_nadir_prediction ) )
pickle_x_axis_average_nadir_prediction = get_prop_dir(1) + run_dir + '/output/python_out/' + main_input_filename_nadir_prediction.replace(".txt","_x_axis_average.pickle")
x_axis_average_nadir_prediction = pickle.load( open( pickle_x_axis_average_nadir_prediction ) ) # be careful, this is in number of seconds this the start date of the SpOCK's simulation. So it is in number of seconds since today. But we want it to be in number of seconds since date_ini
#### High_Drag 
# !!!!!! Uncomment block below
pickle_sma_average_high_drag_prediction = get_prop_dir(1) + run_dir + '/output/python_out/' + main_input_filename_high_drag_prediction.replace(".txt","_sma_average.pickle")
sma_average_high_drag_prediction = pickle.load( open( pickle_sma_average_high_drag_prediction ) )
pickle_x_axis_average_high_drag_prediction = get_prop_dir(1) + run_dir + '/output/python_out/' + main_input_filename_high_drag_prediction.replace(".txt","_x_axis_average.pickle")
x_axis_average_high_drag_prediction = pickle.load( open( pickle_x_axis_average_high_drag_prediction ) ) # be careful, this is in number of seconds this the start date of the SpOCK's simulation. So it is in number of seconds since today. But we want it to be in number of seconds since date_ini
# !!!!!! End of uncomment block below
# End of predicts the future states of the CYGNSS sc from today until today + 6 months (this number of months is a parameter so it can be easily changed)

# Concatenates historical and predictions
print "Concatenating historical and predictions..."
## Convert x_axis_average_nadir_prediction and x_axis_average_high_drag_prediction so that they correspond to the number of seconds since date_ini
date_start_prediction  = main_input_filename_nadir_prediction.split('start_')[1].split('_end')[0].replace("_", ":")  # this basically today but to make sure we grab it from the main input file name of the SpOCK's predictions (in case the run was really slow so when we gert to this line today has changed...)
date_end_prediction  = main_input_filename_nadir_prediction.split('end_')[1].split('_nadir')[0].replace("_", ":")  # this basically today but to make sure we grab it from the main input file name of the SpOCK's predictions (in case the run was really slow so when we gert to this line today has changed...)
date_start_prediction = datetime.strptime(date_start_prediction, "%Y-%m-%dT%H:%M:%S")
date_end_prediction = datetime.strptime(date_end_prediction, "%Y-%m-%dT%H:%M:%S")
for isc in range(nb_sc):
    nb_tle_times_two = len(x_axis_average_nadir_prediction[isc])
    for itle_times_two in range(nb_tle_times_two):
        x_axis_average_nadir_prediction[isc][itle_times_two] = x_axis_average_nadir_prediction[isc][itle_times_two] + ( date_start_prediction - date_ini ).total_seconds()
# !!!!!!  uncomment block below 
for isc in range(nb_sc):
    nb_tle_times_two = len(x_axis_average_high_drag_prediction[isc])
    for itle_times_two in range(nb_tle_times_two):
        x_axis_average_high_drag_prediction[isc][itle_times_two] = x_axis_average_high_drag_prediction[isc][itle_times_two] + ( date_start_prediction - date_ini ).total_seconds()
# !!!!!! end of  uncomment block below 
## Concatenate historical and predictions
sma_average_historical_and_nadir_prediction = []
x_axis_average_historical_and_nadir_prediction = []
y_axis_average_historical_and_nadir_prediction = []
sma_average_historical_and_high_drag_prediction = []
x_axis_average_historical_and_high_drag_prediction = []
y_axis_average_historical_and_high_drag_prediction = []

for isc in range(nb_sc):
    # Nadir
    y_axis_average_historical_and_nadir_prediction_per_sc = []
    sma_average_historical_and_nadir_prediction.append( sma_average_historical[isc] + sma_average_nadir_prediction[isc] )
    x_axis_average_historical_and_nadir_prediction.append( x_axis_average_historical[isc] + x_axis_average_nadir_prediction[isc] )
    # y axis is sma_average_historical_and_nadir_prediction but each element of sma_average_historical_and_nadir_prediction is repeated so that they are kept constant during an orbit
    nb_orbit = len(sma_average_historical_and_nadir_prediction[-1])
    for iorbit in range(nb_orbit):
        y_axis_average_historical_and_nadir_prediction_per_sc.append( sma_average_historical_and_nadir_prediction[-1][iorbit] )
        y_axis_average_historical_and_nadir_prediction_per_sc.append( sma_average_historical_and_nadir_prediction[-1][iorbit] )
    y_axis_average_historical_and_nadir_prediction.append( y_axis_average_historical_and_nadir_prediction_per_sc )
    # High drag
    y_axis_average_historical_and_high_drag_prediction_per_sc = []
    sma_average_historical_and_high_drag_prediction.append( sma_average_historical[isc] + sma_average_high_drag_prediction[isc] )
    x_axis_average_historical_and_high_drag_prediction.append( x_axis_average_historical[isc] + x_axis_average_high_drag_prediction[isc] )
    # y axis is sma_average_historical_and_high_drag_prediction but each element of sma_average_historical_and_high_drag_prediction is repeated so that they are kept constant during an orbit
    nb_orbit = len(sma_average_historical_and_high_drag_prediction[-1])
    for iorbit in range(nb_orbit):
        y_axis_average_historical_and_high_drag_prediction_per_sc.append( sma_average_historical_and_high_drag_prediction[-1][iorbit] )
        y_axis_average_historical_and_high_drag_prediction_per_sc.append( sma_average_historical_and_high_drag_prediction[-1][iorbit] )
    y_axis_average_historical_and_high_drag_prediction.append( y_axis_average_historical_and_high_drag_prediction_per_sc )

# end of concatenates historical and predictions

# Plot the semi-major axis of the 8 CYGNSS as a function of time
## For plots, generate disctinct colors
NCURVES = nb_sc
np.random.seed(101)
curves = [np.random.random(20) for i in range(NCURVES)]
values = range(NCURVES)
jet = cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

label_arr = ['FM05', 'FM04', 'FM02', 'FM01', 'FM08', 'FM06', 'FM07', 'FM03'] 

for isc in range(nb_sc):
            if isc == 0:
                # Plot

                fig_title = 'Orbit average SMA - ' + r'$\mathbf{R_E}$ as a function of time'
                y_label = 'SMA - ' + r'$\mathbf{R_E}$'
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
            ax_sma_average.plot(x_axis_average_historical_and_nadir_prediction[isc], np.array(y_axis_average_historical_and_nadir_prediction[isc]) - earth_radius, linewidth = 2, color = colorVal, label = label_arr[isc]) # !!!!!! uncomment
            ax_sma_average.plot(x_axis_average_historical_and_high_drag_prediction[isc], np.array(y_axis_average_historical_and_high_drag_prediction[isc]) - earth_radius, linewidth = 2, color = colorVal, label = '') # !!!!!! uncomment

            if isc == nb_sc - 1:
                # x axis label is in real time
                ## all output files of one simulation have the same number of steps, and start at the same date
                nb_ticks_xlabel = 8
                nb_seconds_in_simu =  ( date_end_prediction - date_ini ).total_seconds()
                dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
                start_xaxis_label = 0 # we start on date_ini
                xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
                date_list_str = []
                date_list = [date_ini + timedelta(seconds=x-xticks[0]) for x in xticks]
                for i in range(len(xticks)):
                    if dt_xlabel >= 3*24*3600:
                        date_list_str.append( str(date_list[i])[5:10] )
                    else:
                        date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
                ax_sma_average.xaxis.set_ticks(xticks)
                ax_sma_average.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                ax_sma_average.margins(0,0); ax_sma_average.set_xlim([min(xticks), max(xticks)])
                nb_seconds_from_date_ini_until_today = ( date_start_prediction - date_ini ).total_seconds()
                ax_sma_average.plot([nb_seconds_from_date_ini_until_today, nb_seconds_from_date_ini_until_today], [ax_sma_average.get_ylim()[0], ax_sma_average.get_ylim()[1]], linewidth = 2, linestyle = 'dotted', color = 'k') 
                ax_sma_average.text(nb_seconds_from_date_ini_until_today / 2., ax_sma_average.get_ylim()[0] + ( ax_sma_average.get_ylim()[1] - ax_sma_average.get_ylim()[0])/50. , 'Historical TLEs' , verticalalignment = 'bottom', horizontalalignment = 'center', fontsize = fontsize_plot, weight = 'bold')
                ax_sma_average.text(( nb_seconds_from_date_ini_until_today + ( nb_seconds_in_simu - nb_seconds_from_date_ini_until_today ) / 2.  ), ax_sma_average.get_ylim()[0] + ( ax_sma_average.get_ylim()[1] - ax_sma_average.get_ylim()[0])/50., "SpOCK's predictions" , verticalalignment = 'bottom', horizontalalignment = 'center', fontsize = fontsize_plot, weight = 'bold')
                ax_sma_average.text(nb_seconds_from_date_ini_until_today, ax_sma_average.get_ylim()[0] + ( ax_sma_average.get_ylim()[1] - ax_sma_average.get_ylim()[0])/12. , str(date_start_prediction)[0:10] , verticalalignment = 'bottom', horizontalalignment = 'right', fontsize = fontsize_plot, weight = 'bold', rotation = 90)


                legend = ax_sma_average.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
                legend.get_title().set_fontsize(str(fontsize_plot))

                if save_plots == 1:
                    fig_save_name = 'sma_average_historical_until_' + str(date_start_prediction)[0:10] + '_spock_prediction_until_' + str(date_end_prediction)[0:10]
                    fig_save_name = root_save_fig_name + fig_save_name + '.jpeg'
                    fig_sma_average.savefig(fig_save_name, facecolor=fig_sma_average.get_facecolor(), edgecolor='none', bbox_inches='tight')  
