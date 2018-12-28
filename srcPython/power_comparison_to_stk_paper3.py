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
from degree_to_time import *
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
from read_input_file import *
from read_output_file import *
from get_prop_dir import *
from compute_power import *
from datetime import datetime, timedelta
from orbit_average import *
from dateutil.relativedelta import relativedelta
plt.ion()

input_filename = 'cygnss_power.txt'
# normalization_coef = 10 # the run is with a area of 0.1 m^2, so multiply by 10 to get a normalized power (= power for an area of 1 m^2)
# coeff_real = 0.4*0.3*4*0.28*0.74 # length + width * nb of panels * solar cell eff * packing eff

# OUR PROPAGATOR
## READ PROPAGATOR INPUT FILE
input_filename = get_prop_dir(1) + 'run_paper3/input/main_input/' + input_filename
input_variables, order_input_variables = read_input_file(input_filename)
date_start = input_variables[0]; date_stop = input_variables[1]; dt = input_variables[2]; nb_steps = input_variables[3]; nb_satellites = input_variables[4]; 
output_propagator_path = input_variables[6]; output_propagator_file = input_variables[7]; gps_name = input_variables[5]; 
nb_surfaces = input_variables[8]

## COMPUTE POWER
power = np.zeros([nb_satellites,nb_steps, nb_surfaces])
sun_elevation = np.zeros([nb_satellites,nb_steps])
for isat in range(nb_satellites):
    power[isat, :, :], sun_elevation[isat, :] = compute_power(output_propagator_path[isat] + 'power_'+ output_propagator_file[isat])
total_power = np.zeros([nb_satellites,nb_steps])
for isat in range(nb_satellites):
    for istep in range(nb_steps):
        total_power[isat, istep] = np.sum(power[isat, istep, :] ) 

## READ LATITUDE TO THEN COMPUTE ORBIT-AVERAGE POWER
latitude = np.zeros([nb_satellites, nb_steps]); local_time = np.zeros([nb_satellites, nb_steps])
list_var = ["latitude", "local_time"]
lat_var_nb = 1; local_time_var_nb = 2
for isat in range(nb_satellites):
    out_var, order_var = read_output_file( output_propagator_path[isat] + output_propagator_file[isat] , list_var )
    if isat == 0:
        date = out_var[0]
    latitude[isat, :] = out_var[lat_var_nb]
    local_time[isat, :] = out_var[local_time_var_nb]
date_array = np.array(date)

## COMPUTE ORBIT-AVERAGE POWER 
power_average = []
time_average = []
nb_orbits = np.zeros([nb_satellites])
index_time_average = []
for isat in range(nb_satellites):
    power_average_sub = []; time_average_sub = []; index_time_average_sub = []
    power_average_sub, time_average_sub, index_time_average_sub = orbit_average(total_power[isat, :], latitude[isat, :], date)
    nb_orbits[isat] = len(power_average_sub)
    power_average.append(power_average_sub); time_average.append(time_average_sub); index_time_average.append(index_time_average_sub)
power_average = np.array(power_average)    
index_time_average = np.array(index_time_average) # index_time_average[isat, iorbit, pos]: shows the index in date of orbit iorbit for satellite isat. pos is either 0 (beginning of orbit), 1 (middle of orbit), 2 (end of orbit)
time_average_array = np.array(time_average)


# STK
file_power_stk = open("paper3_stk_results/cygnss_power_stk_power.txt")
read_file_power_stk = file_power_stk.readlines()
nb_surfaces_stk = 5
power_stk = np.zeros([nb_steps, nb_surfaces_stk])
for isurf in range(nb_surfaces_stk):
    for istep in range(nb_steps):
        power_stk[istep, isurf]  = np.float( read_file_power_stk[istep + 12 + ( nb_steps + 6 ) * isurf - 1  ].split()[4] )
total_power_stk = np.zeros([nb_steps])
for istep in range(nb_steps):
    total_power_stk[istep] = np.sum(power_stk[istep, :] ) 

# PLOT
width_fig = 15.
height_fig = width_fig * 3 /4
fontsize_plot = 20 # 9
fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.05, right=0.97, top = 0.90,bottom = 0.06)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(gs[0, 0])
#ax1 = fig.add_subplot(111)
ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
#ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

isat = 0
#CYGNSS: top is 0+1+2; front is 4; back is -1
index_orbit_start = index_time_average[isat, 0, 0]; index_orbit_stop = index_time_average[isat, 1, 0] 
#    index_orbit_start = index_time_average[isat, local_max[isat][1], 0]; index_orbit_stop = index_time_average[isat, local_max[isat][1]+1, 0]
x_axis = np.arange(0, (index_orbit_stop - index_orbit_start) * dt /60., dt/60.)
ax1.plot(x_axis, power[isat, index_orbit_start:index_orbit_stop, 0] + power[isat, index_orbit_start:index_orbit_stop, 1] + power[isat, index_orbit_start:index_orbit_stop, 2],'b', linewidth = 2, label = 'Zenith' )
ax1.plot(x_axis, power[isat, index_orbit_start:index_orbit_stop, 4],'r', linewidth = 2, label = 'Ram' )
ax1.plot(x_axis, power[isat, index_orbit_start:index_orbit_stop, -1],'g', linewidth = 4, label = 'Wake')
ax1.plot(x_axis, total_power[isat, index_orbit_start:index_orbit_stop],'k', linewidth = 2, label = 'Total' )

ax1.plot(x_axis, power_stk[ index_orbit_start:index_orbit_stop, 2] + power_stk[ index_orbit_start:index_orbit_stop, 3] + power_stk[ index_orbit_start:index_orbit_stop, 4],'b', marker = 'o', linestyle = '', linewidth = 2, label = 'STK Zenith' )
ax1.plot(x_axis, power_stk[ index_orbit_start:index_orbit_stop, 0],'r', linewidth = 2,marker = 'o', linestyle = '',  label = 'STK Ram' )
ax1.plot(x_axis, power_stk[ index_orbit_start:index_orbit_stop, 1],'g', linewidth = 4, marker = 'o', linestyle = '', label = 'STK Wake')
ax1.plot(x_axis, total_power_stk[ index_orbit_start:index_orbit_stop], 'k', marker = 'o', linestyle = '', linewidth = 2, label = 'STK Total' )

ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)



ax1.set_title('CYGNSS power over one orbit (' + date[index_orbit_start][0:16] + ' to ' + date[index_orbit_stop][11:16] + ')', weight = 'bold', fontsize = 20,  y = 1.005) 
ax1.legend(ncol = 2, fontsize = fontsize_plot)
ax1.set_xlabel('Time in orbit (min)', fontsize = fontsize_plot, weight = 'bold')
ax1.set_ylabel('Power (W)', fontsize = fontsize_plot, weight = 'bold')

fig_save_name = '/raid3/Armada/Charles/python/CYGNSS/result/image/power_one_orbit_compare_with_stk.eps'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

kind_of_rms_temp = np.zeros([nb_steps])
for istep in np.arange(index_orbit_start,index_orbit_stop):
    kind_of_rms_temp[istep] = ( total_power_stk[istep] - total_power[0,istep] ) **2
kind_of_rms = np.sqrt( np.mean(kind_of_rms_temp) / np.mean(total_power_stk[index_orbit_start:index_orbit_stop]**2) )
print kind_of_rms

