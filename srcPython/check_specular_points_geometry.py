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

import os
from read_input_file import *
from matplotlib import pyplot as plt
from get_prop_dir import *
import sys
from datetime import datetime, timedelta
plt.ion()

run_dir = 'run_cygnss'
dt_output_sat_position_in_spec_file = 60. # in seconds, the time step of ouput of the CYGNSS and GPS positions (every 60 seconds usually)

input_filename = sys.argv[1]

# READ THE PROPAGATOR INPUT FILE
propagator_directory = get_prop_dir(1)
input_filename = propagator_directory + run_dir  + '/input/main_input/' + input_filename

input_variables, order_input_variables = read_input_file(input_filename)
date_start = input_variables[0]; date_stop = input_variables[1]; dt = input_variables[2]; nb_steps = input_variables[3]; nb_satellites = input_variables[4]; 
output_path_propagator = input_variables[6]; output_file_propagator = input_variables[7]; gps_name = input_variables[5]; 


# READ THE POSITIONS OF THE SPECULAR POINTS, THE GPS AND THE CYGNSS SATELLITES
nb_spec_pts = 4;
#interpolation_step = 1 # in second, interpolation step of find_specular_points.c (1 s usually)
nb_steps_interpolation = (int)((nb_steps-1) * dt / dt_output_sat_position_in_spec_file) # the specular output of find_specular_points.c stops one time step before the end of the propagation

nb_gps = len(gps_name)
ecef_spec = np.zeros([nb_spec_pts, nb_satellites, nb_steps_interpolation, 3])
ecef_gps = np.zeros([nb_spec_pts, nb_satellites, nb_steps_interpolation, 3])
ecef_sat = np.zeros([nb_spec_pts, nb_satellites, nb_steps_interpolation, 3])
spec_to_gps = np.zeros([3])
spec_to_sat = np.zeros([3])
angle_vertical_to_gps = np.zeros([nb_spec_pts])
angle_vertical_to_sat = np.zeros([nb_spec_pts])
max_angle_vertical_to_gps = np.zeros([nb_satellites, nb_steps_interpolation])
max_angle_vertical_to_sat = np.zeros([nb_satellites, nb_steps_interpolation])
for i in range(nb_satellites):
    file_specular = open(output_path_propagator[i] + "/specular_" + output_file_propagator[i], "r")
    read_file_specular  = file_specular.readlines()
    # Nb of lines in the spec file header
    if (i == 0):
        nb_lines_header_output_file_spec = 0
        while (read_file_specular[nb_lines_header_output_file_spec].split()[0] != "#START"):
            nb_lines_header_output_file_spec = nb_lines_header_output_file_spec + 1
        nb_lines_header_output_file_spec = nb_lines_header_output_file_spec + 1
    ispec_save = 0
    j = -1
    while (ispec_save < len(read_file_specular)-1-nb_spec_pts):
        time_spec_sublist_temp_ini = read_file_specular[ispec_save+nb_lines_header_output_file_spec].split()[0] 
        time_spec_sublist_temp_ini = datetime.strptime(time_spec_sublist_temp_ini, "%Y-%m-%dT%H:%M:%S")
        time_since_start_temp = datetime.strptime(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[0], "%Y-%m-%dT%H:%M:%S") - date_start
        time_since_start = time_since_start_temp.days*24*3600. + time_since_start_temp.seconds
        if ( time_since_start % dt_output_sat_position_in_spec_file == 0 ): # the positions of the GPS and CYGNSS satellites are output only every 60 seconds
            j = j + 1
            ecef_spec[0,i,j,0] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[1])
            ecef_spec[0,i,j,1] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[2])
            ecef_spec[0,i,j,2] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[3])
            ecef_gps[0,i,j,0] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[8])
            ecef_gps[0,i,j,1] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[9])
            ecef_gps[0,i,j,2] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[10])
            ecef_sat[0,i,j,0] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[14])
            ecef_sat[0,i,j,1] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[15])
            ecef_sat[0,i,j,2] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[16])
            for icoord in range(3):
                spec_to_gps[icoord] = ecef_gps[0,i,j,icoord] - ecef_spec[0,i,j,icoord]
                spec_to_sat[icoord] = ecef_sat[0,i,j,icoord] - ecef_spec[0,i,j,icoord]
            angle_vertical_to_gps[0] = np.arccos( np.dot(spec_to_gps, ecef_spec[0,i,j,:]) / ( np.linalg.norm(spec_to_gps) * np.linalg.norm(ecef_spec[0,i,j,:]) ) )
            angle_vertical_to_sat[0] = np.arccos( np.dot(spec_to_sat, ecef_spec[0,i,j,:]) / ( np.linalg.norm(spec_to_sat) * np.linalg.norm(ecef_spec[0,i,j,:]) ) )
            ispec = 1
            while (datetime.strptime(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[0], "%Y-%m-%dT%H:%M:%S")  == time_spec_sublist_temp_ini):
                ecef_spec[ispec,i,j,0] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[1])
                ecef_spec[ispec,i,j,1] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[2])
                ecef_spec[ispec,i,j,2] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[3])
                ecef_gps[ispec,i,j,0] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[8])
                ecef_gps[ispec,i,j,1] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[9])
                ecef_gps[ispec,i,j,2] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[10])
                ecef_sat[ispec,i,j,0] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[14])
                ecef_sat[ispec,i,j,1] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[15])
                ecef_sat[ispec,i,j,2] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[16])
                for icoord in range(3):
                    spec_to_gps[icoord] = ecef_gps[ispec,i,j,icoord] - ecef_spec[ispec,i,j,icoord]
                    spec_to_sat[icoord] = ecef_sat[ispec,i,j,icoord] - ecef_spec[ispec,i,j,icoord]
                angle_vertical_to_gps[ispec] = np.arccos( np.dot(spec_to_gps, ecef_spec[ispec,i,j,:]) / ( np.linalg.norm(spec_to_gps) * np.linalg.norm(ecef_spec[ispec,i,j,:]) ) ) * 180. /np.pi
                angle_vertical_to_sat[ispec] = np.arccos( np.dot(spec_to_sat, ecef_spec[ispec,i,j,:]) / ( np.linalg.norm(spec_to_sat) * np.linalg.norm(ecef_spec[ispec,i,j,:]) ) ) * 180. /np.pi
                ispec = ispec + 1
                if (nb_lines_header_output_file_spec+ispec_save+ispec == len(read_file_specular)-1):
                    break
            max_angle_vertical_to_gps[i,j] = max(angle_vertical_to_gps)
            max_angle_vertical_to_sat[i,j] = max(angle_vertical_to_sat)
            ispec_save = ispec + ispec_save
        else:
            ispec = 1
            while (datetime.strptime(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[0], "%Y-%m-%dT%H:%M:%S")  == time_spec_sublist_temp_ini):
                ispec = ispec + 1
                if (nb_lines_header_output_file_spec+ispec_save+ispec == len(read_file_specular)-1):
                    break
            ispec_save = ispec + ispec_save

raise Exception
# PLOT 
fig = plt.figure(num=None, figsize=(17, 11), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(111)
ax1.set_title('Absolute value of the angle difference between the local vertical and the GPS/CYGNSS directions', weight = 'bold', fontsize = 20,  y = 1.008) 
ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
i = 3
x_axis = np.arange(0,nb_steps_interpolation * dt_output_sat_position_in_spec_file / 3600., dt_output_sat_position_in_spec_file / 3600.)
ax1.plot(x_axis, abs( max_angle_vertical_to_sat[i,:] - max_angle_vertical_to_gps[i,:] ) )
ax1.set_xlabel('Time (hours)', fontsize = 18, weight = 'bold')
ax1.set_ylabel('Absolute value of the angle difference ('+u'\N{DEGREE SIGN}'+')', fontsize = 18, weight = 'bold')

for i in range(nb_satellites):
    print max(abs( max_angle_vertical_to_sat[i,:] - max_angle_vertical_to_gps[i,:] ) )
#fig.savefig('check_spec_angle.png', facecolor=fig.get_facecolor(), edgecolor='none')  




