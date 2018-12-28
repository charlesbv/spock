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
import pickle
from datetime import datetime, timedelta
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from read_input_file import *
from read_output_file import *
import colorsys
import sys
from get_name_mission import *

plt.ion()
plt.isinteractive()

save_results = 1

# Read the input file of propagator
input_filename = get_prop_dir(1) + sys.argv[1] + "/input/main_input/" + sys.argv[2]
input_variables, order_input_variables = read_input_file(input_filename)
date_start = input_variables[0]; date_stop = input_variables[1]; dt = input_variables[2]; nb_steps = input_variables[3]; nb_satellites = input_variables[4]; output_path_propagator = input_variables[6]; output_file_propagator = input_variables[7]

path_folder_results = '/raid3/Armada/Charles/python/' #get_prop_dir(2) + 'output/python_propagator/'
name_subfolder_save = output_file_propagator[0][:-5] + "/"
name_mission = get_name_mission(output_file_propagator[0])
save_pickle_name = path_folder_results + name_mission + '/pickle/' +  output_file_propagator[0].replace(".txt",".pickle")
root_save_fig_name = path_folder_results + name_mission + '/result/image/' +  output_file_propagator[0].replace(".txt","_") 
root_save_video_name = path_folder_results + name_mission + '/result/video/' +      output_file_propagator[0].replace(".txt","_")



# #Read the output files
position = np.zeros([nb_satellites, nb_steps, 3])
velocity = np.zeros([nb_satellites, nb_steps, 3])
longitude = np.zeros([nb_satellites, nb_steps])
latitude = np.zeros([nb_satellites, nb_steps])
altitude = np.zeros([nb_satellites, nb_steps])
true_ano = np.zeros([nb_satellites, nb_steps])
raan = np.zeros([nb_satellites, nb_steps])
arg_perigee = np.zeros([nb_satellites, nb_steps])
right_asc = np.zeros([nb_satellites, nb_steps])
local_time = np.zeros([nb_satellites, nb_steps])
angle_asc_node_to_sat = np.zeros([nb_satellites, nb_steps])
date = []
list_output_variables_to_read = ["position","velocity","longitude","latitude","altitude","raan","true_anomaly", "arg_perigee", "right_asc", "local_time"]
for i in range(nb_satellites):
    print i,nb_satellites-1
    output_filename = output_path_propagator[i] + output_file_propagator[i]  
    print output_filename
    output_variables, list_output_variables_read = read_output_file(output_filename, list_output_variables_to_read)
    if (i == 0):
        date = output_variables[0]

    position[i] = output_variables[1]
    velocity[i] = output_variables[2]
    longitude[i,:] = output_variables[3]
    latitude[i,:] = output_variables[4]
    altitude[i,:] = output_variables[5]
    true_ano[i,:] = output_variables[6]
    raan[i,:] = output_variables[7]
    arg_perigee[i,:]  = output_variables[8]
    right_asc[i,:]  = output_variables[9]
    local_time[i,:]  = output_variables[10]
    angle_asc_node_to_sat[i,:] = (true_ano[i,:] + arg_perigee[i,:])%360

# Spacing relative to M1 
spacing_relative_M1_minus180_to_180 = np.zeros([nb_satellites, nb_steps])
for i in range(nb_steps):
    for j in range(nb_satellites):
        spacing_relative_M1_minus180_to_180_temp = (angle_asc_node_to_sat[j,i] - angle_asc_node_to_sat[0,i])%360
        if (spacing_relative_M1_minus180_to_180_temp > 180):
            spacing_relative_M1_minus180_to_180[j,i] = spacing_relative_M1_minus180_to_180_temp - 360
        else:
            spacing_relative_M1_minus180_to_180[j,i] = spacing_relative_M1_minus180_to_180_temp


######################################
################################ PLOTS
width_fig = 13
height_fig = 8
fontsize_plot = 14 
fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.05, right=0.99, top = 0.7,bottom = 0.2)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') 
ax1 = fig.add_subplot(gs[0, 0])
ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] 
colors = ['k','r','b','k','m','c','y','g']
linewidth_array = [2,2,2,2,2,2,2,2]
alpha_array = [1,1,1,1,1,1,1,1]
sat_for_angle_distance = np.arange(1,nb_satellites)
label_array = ['C','G','F','D','B','H','E','A']
nb_days = nb_steps * dt / 3600. / 24.
x_axis = np.arange(0, nb_days, nb_days/len(spacing_relative_M1_minus180_to_180[0,:]))

for i in sat_for_angle_distance:
    ax1.plot(x_axis,spacing_relative_M1_minus180_to_180[i,:], color=colors[i], linewidth = linewidth_array[i], alpha = alpha_array[i],markersize = 0.5)
    ax1.plot([0,0],[0,0], color=colors[i], label = label_array[i], linewidth = linewidth_array[i], alpha = alpha_array[i])

ax1.set_xlabel('Time (days)', fontsize = fontsize_plot, weight = 'bold')
ax1.set_ylabel('Angular\ndistance' + u' (\N{DEGREE SIGN})', fontsize = fontsize_plot, weight = 'bold', labelpad=0.0001)
ax1.set_title('Angular spacing to reference satellite (M1)', weight = 'bold', fontsize = 20,  y = 1.04)
ax1.legend(ncol = len(sat_for_angle_distance), bbox_to_anchor=(0.5, 1.02), loc = 10,borderpad = 0.1, frameon = False, fontsize = fontsize_plot)
gs.update(left=0.1, right=0.99, top = 0.9,bottom = 0.07)
ax1.margins(0,0)
ax1.set_ylim([-180,180])

if save_results == 1:
    fig_save_name = 'spacing'
    fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

