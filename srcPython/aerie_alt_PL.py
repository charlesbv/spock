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

from get_prop_dir import *
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from read_input_file import *
from read_output_file import *
from orbit_average import *

plt.ion()
plt.isinteractive()

# Read the input file of propagator
input_filename = "/home/cbv/PropSim/input/main_input/aerie_other_angle.txt"
input_variables, order_input_variables = read_input_file(input_filename)
date_start = input_variables[0]; date_stop = input_variables[1]; dt = input_variables[2];
nb_steps = input_variables[3]; nb_satellites = input_variables[4]; output_propagator_path = input_variables[6]
output_propagator_file = input_variables[7]
    
# Read the output file
list_output_variables_to_read = ["altitude", "latitude"]
# Reference sat (M1)
output_variables, list_output_variables_read = read_output_file(output_propagator_path[0] + output_propagator_file[0], list_output_variables_to_read)
alt_M1 = output_variables[2]; lat_M1 = output_variables[1]; date = output_variables[0]
alt_M1_orbit_averaged, time_M1_averaged = orbit_average(alt_M1, lat_M1, date)
# Sink sat (L1)
output_variables, list_output_variables_read = read_output_file(output_propagator_path[2] + output_propagator_file[2], list_output_variables_to_read)
alt_L1 = output_variables[2]; lat_L1 = output_variables[1]; date = output_variables[0]
alt_L1_orbit_averaged, time_L1_averaged = orbit_average(alt_L1, lat_L1, date)
# Sink sat (L2)
output_variables, list_output_variables_read = read_output_file(output_propagator_path[7] + output_propagator_file[7], list_output_variables_to_read)
alt_L2 = output_variables[2]; lat_L2 = output_variables[1]; date = output_variables[0]
alt_L2_orbit_averaged, time_L2_averaged = orbit_average(alt_L2, lat_L2, date)

# PLOT orbit-averaged altitude difference
width_fig = 17
height_fig = width_fig*3/4
fontsize_plot = 14 
fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.05, right=0.99, top = 0.7,bottom = 0.2)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(gs[0, 0])
#ax1 = fig.add_subplot(111)
ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
ax1.tick_params(axis='both', which='major',  width = 2, pad = 3)
#ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

nb_orbits = len(alt_M1_orbit_averaged)
x_axis = np.arange(0, 24, 24./nb_orbits)
delta_alt_M1_L1 = np.zeros([nb_orbits])
delta_alt_M1_L2 = np.zeros([nb_orbits])
for iorbit in range(nb_orbits):
    delta_alt_M1_L1[iorbit] = alt_L1_orbit_averaged[iorbit] - alt_M1_orbit_averaged[iorbit]
    delta_alt_M1_L2[iorbit] = alt_L2_orbit_averaged[iorbit] - alt_M1_orbit_averaged[iorbit]

ax1.plot( x_axis, delta_alt_M1_L1, color = 'b', linewidth = 2, label = 'L1')
ax1.plot( x_axis, delta_alt_M1_L2, color = 'r', linewidth = 2, label = 'L2')
ax1.legend()
ax1.set_ylabel('Delta altitude (km)', fontsize = fontsize_plot, weight = 'bold')
ax1.set_xlabel('Time (month)', fontsize = fontsize_plot, weight = 'bold')
ax1.set_title('Loss in altitude for the satellites in the lower plane over two years', weight = 'bold', fontsize = 20,  y = 1.008)
ax1.margins(0,0)
fig_save_name = get_prop_dir(2) + 'output/python_propagator/aerie/altitude_sink.png'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')
os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./")

