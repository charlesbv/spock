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

from cadre_read_last_tle import *
from get_prop_dir import *
import matplotlib.gridspec as gridspec
from read_input_file import *
from matplotlib.colors import LogNorm
import pickle
from eci_to_lvlh import *
import sys
import fileinput
import time
import datetime
import numpy as np
from matplotlib import pyplot as plt
import os
import subprocess

plt.ion()
plt.isinteractive()

path_folder_results = '/raid3/Armada/Charles/python/' #get_prop_dir(2) + 'output/python_propagator/'

same_spock_input_file = 1

# Read input file sat1
input_filename_sat1 =    get_prop_dir(1) + sys.argv[1] + '/input/main_input/' + sys.argv[2] 
input_variables, order_input_variables = read_input_file(input_filename_sat1)
sat1_to_plot_path = input_variables[6][0]; sat1_to_plot = input_variables[7][0]
dt_sat1 = input_variables[2]; nb_steps_sat1 = input_variables[3]; 

if ( 'cygnss' in sat1_to_plot.lower() ):
    name_mission = 'CYGNSS' 
elif ( 'cadre' in sat1_to_plot.lower() ):
    name_mission = 'CADRE' 
elif ( 'aerie' in sat1_to_plot.lower() ):
    name_mission = 'AERIE' 
elif ( 'scion' in sat1_to_plot.lower() ):
    name_mission = 'SCION' 
elif ( 'qb50' in sat1_to_plot.lower() ):
    name_mission = 'QB50' 
else:
    name_mission = 'other' 

save_pickle_name = path_folder_results + name_mission + '/pickle/' + sat1_to_plot.replace(".txt",".pickle")
root_save_fig_name = path_folder_results + name_mission + '/result/image/' + sat1_to_plot.replace(".txt","_")
root_save_video_name = path_folder_results + name_mission + '/result/video/' + sat1_to_plot.replace(".txt","_")

# Read output file sat1
r_eci_sat1 = np.zeros([nb_steps_sat1])
var_to_read = ["position"]
var_out, var_out_order = read_output_file( sat1_to_plot_path + sat1_to_plot, var_to_read )
r_eci_sat1 = var_out[1]



# Read input file sat2
if same_spock_input_file == 1:
    sat2_to_plot_path = input_variables[6][1]; sat2_to_plot = input_variables[7][1]
    dt_sat2 = dt_sat1; nb_steps_sat2 = nb_steps_sat1
else:
    input_filename_sat2 =    get_prop_dir(1) + sys.argv[1] + '/input/main_input/' + sys.argv[3] 
    input_variables, order_input_variables = read_input_file(input_filename_sat2)
    sat2_to_plot_path = input_variables[6][0]; sat2_to_plot = input_variables[7][0]
    dt_sat2 = input_variables[2]; nb_steps_sat2 = input_variables[3]; 

save_pickle_name = save_pickle_name + sat2_to_plot.replace(".txt","_")
root_save_fig_name = root_save_fig_name + sat2_to_plot.replace(".txt","_")
root_save_video_name = root_save_video_name + sat2_to_plot.replace(".txt","_")

# Read output file sat2
r_eci_sat2 = np.zeros([nb_steps_sat2])
var_to_read = ["position"]
var_out, var_out_order = read_output_file( sat2_to_plot_path + sat2_to_plot, var_to_read )
r_eci_sat2 = var_out[1]

