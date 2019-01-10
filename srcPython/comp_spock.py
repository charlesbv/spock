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

# This script looks at phase angle and phase rate between sc.
# ASSUMPTIONS:
# - all run have the same dates (start end, time step)

import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/srcPython")
from read_output_file import *
from read_input_file import *
import numpy as np
from find_in_read_input_order_variables import *
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
from orbit_average import *

plt.ion()
run_list = ['fm1.txt', 'fm3.txt']


nb_run = len(run_list)
var_to_read = ['phase_angle',  'nb_seconds_since_start', 'sma', 'latitude']
earth_mu    = 398600.4418; # gravitational parameter (km^3/s^2)
earth_radius        = 6378.137; # mean equatorial radius (km)
sma_average = []
date_average_start_orbit_list = []
date_average_end_orbit_list = []
x_axis_average = []
mean_motion_average = []
phase_rate_average = []
for irun in range(nb_run):
    isc = 0
    input_filename = run_list[irun]
    var_in, var_in_order = read_input_file(input_filename)
    output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
    output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
    dt = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')]; 
    nb_steps = var_in[find_in_read_input_order_variables(var_in_order, 'nb_steps')] -1; 
    nb_sc = var_in[find_in_read_input_order_variables(var_in_order, 'nb_sc')]; 
    if irun == 0: # see assumptions in header
        phase_angle = np.zeros([nb_run, nb_steps])
        phase_rate = np.zeros([nb_run, nb_steps])
        sma = np.zeros([nb_run, nb_steps])
        latitude = np.zeros([nb_run, nb_steps])

    var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
    if irun == 0:
        date = var_out[find_in_read_input_order_variables(var_out_order, 'date')]
        nb_seconds_since_start = var_out[find_in_read_input_order_variables(var_out_order, 'nb_seconds_since_start')]
    #phase_angle[irun, :] = var_out[find_in_read_input_order_variables(var_out_order, 'phase_angle')]
    sma[irun, :] = var_out[find_in_read_input_order_variables(var_out_order, 'sma')] 
    latitude[irun, :] = var_out[find_in_read_input_order_variables(var_out_order, 'latitude')] 
    # # calculate phase rate from SMA (keplerian equations) because if from phase_angle problems every orbit
    # mean_motion = np.sqrt( earth_mu / sma[irun, :]**3 ) # in s^-1
    # phase_rate[irun, :] = mean_motion * 180. / np.pi # in deg/s
    sma_orbit_averaged, time_averaged, index_time_averaged = orbit_average(sma[irun, :], latitude[irun, :], date )
    sma_average.append( sma_orbit_averaged ) # each sc might not have the same orbital period so the length of the array might not be the same between each sc
    date_average_start_orbit_list.append( np.array(time_averaged)[:,0] ) # take the date at the start of the bin
    date_average_end_orbit_list.append( np.array(time_averaged)[:,2] ) # take the date at the end of the bin
    x_axis_average_per_sc = []
    nb_orbit_for_this_sc = len(time_averaged)
    mean_motion_average_per_sc = []
    phase_rate_average_per_sc = []
    for iorbit in range(nb_orbit_for_this_sc):
        date_average_start_orbit = date_average_start_orbit_list[-1][iorbit]
        date_average_start_orbit = datetime.strptime( date_average_start_orbit, "%Y/%m/%d %H:%M:%S.%f" )
        nb_seconds_between_start_orbit_and_date_start = ( date_average_start_orbit - datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S.%f") ).total_seconds()
        date_average_end_orbit = date_average_end_orbit_list[-1][iorbit]
        date_average_end_orbit = datetime.strptime( date_average_end_orbit, "%Y/%m/%d %H:%M:%S.%f" )
        nb_seconds_between_end_orbit_and_date_start = ( date_average_end_orbit - datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S.%f") ).total_seconds()

        x_axis_average_per_sc.append((nb_seconds_between_end_orbit_and_date_start + nb_seconds_between_start_orbit_and_date_start)/2.)


        mean_motion_average_per_sc.append( np.sqrt( earth_mu / sma_average[-1][iorbit]**3 ) ) # in s^-1
        phase_rate_average_per_sc.append( mean_motion_average_per_sc[-1] * 180. / np.pi ) # in deg/s

    x_axis_average.append( x_axis_average_per_sc )
    mean_motion_average.append(mean_motion_average_per_sc)
    phase_rate_average.append(phase_rate_average_per_sc)



phase_rate_average_arr = np.array(phase_rate_average)

period = 95*60.
iperiod = 0
index_period = []
j = 0
i = 0
nb_samples_per_period = []
while i+j < len(nb_seconds_since_start):

    while nb_seconds_since_start[i + j] - nb_seconds_since_start[j] < period:
        i = i + 1
        if i+j >=  len(nb_seconds_since_start):
            break
    if i+j < len(nb_seconds_since_start):
        index_period.append(i+j)
        nb_samples_per_period.append(i)
    j = i + j
    i = 0
    if i+j >=  len(nb_seconds_since_start):
        break
    

nb_period = len(index_period)

sma1_average = np.zeros([nb_period])
r1ecimag_average = np.zeros([nb_period])
v1ecimag_average = np.zeros([nb_period])
sma2_average = np.zeros([nb_period])
r2ecimag_average = np.zeros([nb_period])
v2ecimag_average = np.zeros([nb_period])

mean_motion1 = np.zeros([nb_period])
phase_rate1 = np.zeros([nb_period])
mean_motion2 = np.zeros([nb_period])
phase_rate2 = np.zeros([nb_period])

for i in range(nb_period-1):
    sma1_average[i] = np.mean(sma[0, index_period[i]: index_period[i+1]])
    sma2_average[i] = np.mean(sma[1, index_period[i]: index_period[i+1]])

    mean_motion1[i] = np.sqrt( earth_mu / sma1_average[i]**3 ) # in s^-1                                                                              
    phase_rate1[i] = mean_motion1[i] * 180. / np.pi # in deg/s                                                                                

    mean_motion2[i] = np.sqrt( earth_mu / sma2_average[i]**3 ) # in s^-1                                                                              
    phase_rate2[i] = mean_motion2[i] * 180. / np.pi # in deg/s                                                                                

    
    #print ( phase_rate1[i] - phase_rate2[i] ) * 3600 * 24
    #raise Exception
# #print (phase_rate_average[1][0] - phase_rate_average[0][0])* 3600. * 24
delta_phase_rate = ( phase_rate1 - phase_rate2 ) * 3600 * 24  
# print ( phase_rate1[0] - phase_rate2[0] ) * 3600 * 24  
# #fig, ax = plt.subplots(); ax.plot(nb_seconds_since_start[index_period][:-1]/3600/24, sma2_average[:-1]); ax.plot(nb_seconds_since_start[index_period][:-1]/3600/24, sma1_average[:-1]);  # ax.set_ylim([-2.7,-2.3]);
# fig, ax = plt.subplots(); ax.plot(nb_seconds_since_start[index_period][:-1]/3600/24, sma1_average[:-1] - sma2_average[:-1]);
# fig, ax = plt.subplots(); ax.plot(nb_seconds_since_start[index_period]/3600/24, delta_phase_rate);  ax.set_ylim([-2.8,-2.3]);
# raise Exception
# PLOTS
## Parameters for the figure
height_fig = 15.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 28 
ratio_fig_size = 4./3


fig_title = ''#Phase angle difference as a function of time'

x_label = 'Time (days)'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal')
plt.rc('font', weight='normal') ## make the labels of the ticks in normal
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.1)


## Phase rate difference with 1st run in run_list
y_label = '$\Delta \dot \Phi$ (' + u'\N{DEGREE SIGN}' + '/day)'
ax_rate = fig.add_subplot(gs[0, 0])
ax_rate.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax_rate.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax_rate.spines.itervalues()] # change the width of the frame of the figure
ax_rate.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='normal') ## make the labels of the ticks in normal 

x_axis = nb_seconds_since_start / 3600. / 24

ax_rate.plot(( phase_rate_average_arr[0, :] - phase_rate_average_arr[1, :] ) * 3600. * 24, linewidth = 2, color = 'b')
#ax_rate.plot( ( phase_rate1 - phase_rate2) * 3600. * 24, linewidth = 2, color = 'b')

ax_rate.margins(0,0)

fig_save_name = 'spock_phasse_rate_diff_fm1_fm3.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


plt.show(); plt.show()
