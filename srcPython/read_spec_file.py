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
import os
from read_input_file import *
from matplotlib import pyplot as plt

input_filename = 'cygnss_gps_spec_check_loc.txt'

# READ THE PROPAGATOR INPUT FILE
propagator_directory = get_prop_dir(2) 
input_filename = propagator_directory + '/input/main_input/' + input_filename

input_variables, order_input_variables = read_input_file(input_filename)
date_start = input_variables[0]; date_stop = input_variables[1]; dt = input_variables[2]; nb_steps = input_variables[3]; nb_satellites = input_variables[4]; 
output_path_propagator = input_variables[6]; output_file_propagator = input_variables[7]; 
gps_name = input_variables[5]; 

# READ THE POSITIONS OF THE SPECULAR POINTS, THE GPS AND THE CYGNSS SATELLITES
nb_spec_pts = 4;
interpolation_step = 1 # in second, interpolation step of find_specular_points.c (1 s usually)
nb_steps_interpolation = (int)((nb_steps-1) * dt / interpolation_step) # the specular output of find_specular_points.c stops one time step before the end of the propagation

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
    count_less_than_4_spec = 0
    time_less_than_4_spec = []
    file_specular = open(output_path_propagator[i] + "specular_" + output_file_propagator[i].split('/')[-1], "r")
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
        j = j + 1
        time_spec_sublist_temp_ini = read_file_specular[ispec_save+nb_lines_header_output_file_spec].split()[0] 
        time_spec_sublist_temp_ini = datetime.strptime(time_spec_sublist_temp_ini, "%Y-%m-%dT%H:%M:%S")
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
            if (nb_lines_header_output_file_spec+ispec_save+ispec == len(read_file_specular)):
                break
        if ispec < 4:
            count_less_than_4_spec = count_less_than_4_spec + 1
            time_less_than_4_spec.append(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[0])
        max_angle_vertical_to_gps[i,j] = max(angle_vertical_to_gps)
        max_angle_vertical_to_sat[i,j] = max(angle_vertical_to_sat)
        ispec_save = ispec + ispec_save







