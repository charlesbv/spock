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
# this script computes the time to decrease the SMA of of a sc by x km (for isntance x = 10) for different initial altitudes, F10.7, and cross section area (actually difference mass). It first runs the simu with SpOCK then post-process. Note that ic reated a SpOCK execeutbale spock_cyg that stops the propagation when the sma has dropped by 20 km (since the sma varies within an orbit (by a few km) and we want only the sma average to drop by 10 km i had to take a margin)
# ASSUMPTIONS
# - see section 'PARAMETERES TO SET BEFORE RUNNING THIS SCRIPT'
# - all runs are from the same initial epoch to the same final epoch, with the same dt
import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/kalman/spock_development_new_structure_kalman_dev/srcPython")
#sys.path.append("/home/cbv/spock_development_new_structure_kalman_dev/srcPython")
import pickle
from orbit_average import *
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
from read_input_file import *
from read_output_file import *
from matplotlib.colors import LogNorm
import pickle
import fileinput
import time
from datetime import datetime
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker
import os
import sys
import subprocess
from find_in_read_input_order_variables import *
from spock_main_input import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
# PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT
alt_min = 400
alt_max = 550
dalt = 20
f107_min = 70
f107_max = 200
df107 = 15
area_arr = [2,5,10] # factor to apply on the area (wctually the invert will be aplied on mass)
run_spock = 0 # if set to other than 1, won't trun SpOCK and won'r create the SpOCK main input files
delta_sma = 10. # number of km to look at sma decrease due to drag (the 'x' mentioned in the header of this script) !!!!  if you cahnge this you need to make sure to change the executbale spock_cyg accordinginly (see explanation in header of this script))
# end of PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT
earth_mu    = 398600.4418; # gravitational parameter (km^3/s^2)

alt_arr = np.arange(alt_min, alt_max + dalt, dalt)
nb_alt = len(alt_arr)
f107_arr = np.arange(f107_min, f107_max + df107, df107)
nb_f107 = len(f107_arr)
area_arr = np.array(area_arr)
nb_area = len(area_arr)

nb_run = nb_alt * nb_f107 * nb_area

mass_nominal = 29.
mass_nominal = np.float(mass_nominal)
spice_path = '/Users/cbv/cspice/data' #'/raid4/cbv/cspice/data'  # '/Users/cbv/cspice/data'
sma = []
period = []
x_axis_average = []
date_lost_delta_sma = []
time_to_lose_delta_sma = np.zeros([nb_alt, nb_f107, nb_area]) + 1*365*24*3600
index_time_to_lose_delta_sma = np.zeros([nb_alt, nb_f107, nb_area]) - 1
irun_count = -1
for ialt in range(nb_alt):
    sma_ialt = []
    period_ialt = []
    x_axis_average_ialt = []
    date_lost_delta_sma_ialt = []
    alt = alt_arr[ialt]
    for if107 in range(nb_f107):
        sma_if107 = []
        period_if107 = []
        x_axis_average_if107 = []
        date_lost_delta_sma_if107 = []
        f107 = f107_arr[if107]
        for iarea in range(nb_area):
            irun_count = irun_count + 1
            print irun_count, nb_run
            mass = mass_nominal / area_arr[iarea]
            main_input_filename = 'alt' + str(alt) + '_f107' + str(f107) + '_area' + str(area_arr[iarea]) + '.txt' 
            if run_spock == 1:
                # CREATE MAIN INPUT FILES FOR SPOCK AND RUN SPOCK
                spock_main_input(
                main_input_filename,
                # for TIME section
                '2017-01-01T00:00:00',
                '2017-12-31T00:00:00',
                60., #dt
                # for SPACECRAFT section
                1,
                '0',
                mass,
                'cygnss_geometry_2016_acco08.txt',
                # for ORBIT section
                ['oe', str(alt) + ' 35 0 0 0 0'],
                # for FORCES section
                4,
                'drag solar_pressure moon_gravity sun_gravity',
                ['static', str(f107), str(f107), '15'],
                # for OUTPUT section
                'out',
                300, # the last time step is the only we care about here and is always printed
                # for ATTITUDE section
                "nadir",
                # for GROUNDS_STATIONS section
                "0",#"my_ground_stations.txt"
                # for SPICE section
                spice_path, 
                # for DENSITY_MOD section
                1
                )

                os.system("mpirun -np 1 spock_cyg " + main_input_filename)
                
            # POST PROCESS
            input_filename = main_input_filename
            var_in, var_in_order = read_input_file(input_filename)
            output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
            output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
            dt = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')]; 
            nb_steps = var_in[find_in_read_input_order_variables(var_in_order, 'nb_steps')]; 
            # read the sma of each sc and compute the orbit average sma
            var_to_read = ["sma", "latitude"]
            var_out, var_out_order = read_output_file( output_file_path_list[0] + output_file_name_list[0], var_to_read )
            date = var_out[find_in_read_input_order_variables(var_out_order, 'date')]
            sma_temp = var_out[find_in_read_input_order_variables(var_out_order, 'sma')]
            latitude = var_out[find_in_read_input_order_variables(var_out_order, 'latitude')]
            sma_orbit_averaged, time_averaged, index_time_averaged = orbit_average(sma_temp, latitude, date )
            x_axis_average_temp = []
            nb_orbit_for_this_sc = len(time_averaged)
            date_average_start_orbit_list = np.array(time_averaged)[:,0]  # take the date at the start of the bin
            already_found_lost_sma = 0
            for iorbit in range(nb_orbit_for_this_sc):
                date_average_start_orbit = date_average_start_orbit_list[iorbit]
                date_average_start_orbit = datetime.strptime( date_average_start_orbit, "%Y/%m/%d %H:%M:%S.%f" )
                nb_seconds_between_start_orbit_and_date_start = ( date_average_start_orbit - datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S.%f") ).total_seconds()
                x_axis_average_temp.append( nb_seconds_between_start_orbit_and_date_start )
                if ( ( already_found_lost_sma != 1 ) & (sma_orbit_averaged[iorbit] < sma_orbit_averaged[0] - delta_sma ) ):
                    time_to_lose_delta_sma[ialt,if107,iarea] = nb_seconds_between_start_orbit_and_date_start
                    date_lost_delta_sma_if107.append(date_average_start_orbit)
                    index_time_to_lose_delta_sma[ialt,if107,iarea] = iorbit
                    already_found_lost_sma = 1
            x_axis_average_temp = np.array(x_axis_average_temp)

            # convert each sma into a period (assuming keplerian elements)
            period_temp = 2*np.pi * np.sqrt( np.array(sma_orbit_averaged)**3 / earth_mu )
            period_if107.append(period_temp)
            sma_if107.append(sma_orbit_averaged)
            x_axis_average_if107.append(x_axis_average_temp)
#             if ( (if107 == 2) & ( iarea == 2 ) ):
#                 raise Exception


        sma_ialt.append(sma_if107)
        date_lost_delta_sma_ialt.append(date_lost_delta_sma_if107)
        period_ialt.append(period_if107)
        x_axis_average_ialt.append(x_axis_average_if107)
    sma.append(sma_ialt)
    date_lost_delta_sma.append(date_lost_delta_sma_ialt)
    period.append(period_ialt)
    x_axis_average.append(x_axis_average_ialt)

sma = np.array(sma)
period = np.array(period)
x_axis_average = np.array(x_axis_average) # x_axis_average[ialt, if107, iarea] is the list of times for this alt, f107, and area. Same for sma and period

## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3
### For plots, generate disctinct colors
NCURVES = nb_run
np.random.seed(101)
curves = [np.random.random(20) for i in range(NCURVES)]
values = range(NCURVES)
jet = cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
### or choose among color_arr
color_arr = ['b', 'r','k','cornflowerblue','g', 'm', 'gold', 'cyan', 'fuchsia', 'lawngreen', 'darkgray', 'green', 'chocolate']

## Plot the contour of time to decrease sma by delta_sma as a function of f107 and intial apogee altitude for a given area factor(3d plot) 
#iarea = 1
for iarea in range(nb_area):
    fig_title = 'Contour of # days to lose ' + str((int)(delta_sma)) + ' km VS alt. and F10.7 (area ' + r'$\times$ ' + str(area_arr[iarea]) + ')' 
    y_label = 'Altitude (km)'
    x_label = 'F10.7'
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle(fig_title, y = 0.96,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
    ax = fig.add_subplot(gs[0, 0])

    ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)


    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

    x = f107_arr
    y = alt_arr

    X, Y = np.meshgrid(x, y)
    Z = time_to_lose_delta_sma[:,:,iarea] / 3600. / 24


    nr, nc = Z.shape

    Z = np.ma.array(Z)


    levels_lin = [1,3, 7, 12, 18,25,35,45,60, 80, 110, 150, 210,280, 360] #np.arange(1,365+10,20)
    levels  = levels_lin
    origin = 'lower'
    CS1 = ax.contourf(X, Y, Z, levels,
                      #[-1, -0.1, 0, 0.1],
                      #alpha=0.5,
#                      locator=ticker.LogLocator(),
                      norm = LogNorm(),
                      cmap = plt.cm.get_cmap("rainbow"),
                      origin=origin)
    ax.scatter(85, 520, marker = '*', s = 200, color = 'k')
    ax.plot([f107_arr[0], f107_arr[-1]], [520, 520], linestyle = 'dashed', linewidth = 1, color = 'k')
    fmt = matplotlib.ticker.FormatStrFormatter("%g")
    cbar = plt.colorbar(CS1, ax = ax, ticks = levels, format = fmt)
    cbar.ax.set_ylabel('# days', fontsize = fontsize_plot, weight = 'bold')
    cbar.ax.tick_params(labelsize=fontsize_plot) 


    fig_save_name = 'contour_time_to_lose_' + str((int)(delta_sma)) + '_vs_altitude_and_f107_area_' + str(area_arr[iarea])
    fig_save_name =  fig_save_name + '.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  



# for a given altitude and f107 plot the sma average as a function of time for all areas
alt = 520
f107 = 85
if len(np.where(alt_arr == alt)[0]) == 0:
    print "***! The altitude chosen is not in the altitude array. Please choose a different altitude. The program will stop !***"; raise Exception
ialt = np.where(alt_arr == alt)[0][0]
if len(np.where(f107_arr == f107)[0]) == 0:
    print "***! The f107 chosen is not in the f107 array. Please choose a different f107. The program will stop !***"; raise Exception
if107 = np.where(f107_arr == f107 )[0][0]

fig_title = 'Orbit average SMA as a function of time (altitude: ' + str(alt) + ' km, F10.7: ' + str(f107) + ')'
y_label = 'SMA (km)'
x_label = 'Day #'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.963,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

max_time_to_lose_delta_sma = np.max(time_to_lose_delta_sma[ialt, if107, :])

for iarea in range(nb_area):
    ax.plot(x_axis_average[ialt, if107, iarea][:(int)(index_time_to_lose_delta_sma[ialt, if107, iarea])]/3600./24, sma[ialt, if107, iarea][:(int)(index_time_to_lose_delta_sma[ialt, if107, iarea])], linewidth = 2, color = color_arr[iarea], label =  str((int)(area_arr[iarea])))
    ax.text(x_axis_average[ialt, if107, iarea][(int)(index_time_to_lose_delta_sma[ialt, if107, iarea])-1]/3600./24, sma[ialt, if107, iarea][(int)(index_time_to_lose_delta_sma[ialt, if107, iarea])-1], str((int)(time_to_lose_delta_sma[ialt, if107, iarea]/3600./24)), fontsize = fontsize_plot, weight = 'bold', color= color_arr[iarea], horizontalalignment = 'left', verticalalignment = 'bottom')
    if iarea == 0:
        # x axis label is in number of days
        nb_ticks_xlabel = 10
        dt_xlabel =  max_time_to_lose_delta_sma / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
        xticks_float = np.arange(0, max_time_to_lose_delta_sma+1, dt_xlabel)/3600./24
        xticks = xticks_float.astype(np.int)
        ax.xaxis.set_ticks(xticks)
        ax.margins(0,0); ax.set_xlim([min(xticks), max(xticks)])

legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Area " + r'$\times$', fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = 'sma_vs_time_alt' + str(alt) + '_f107' + str(f107) 
fig_save_name =  fig_save_name + '.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

