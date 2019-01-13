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
from datetime import datetime
import numpy as np
from matplotlib import pyplot as plt
import os
import subprocess
from get_name_mission import *

plt.ion()
plt.isinteractive()

############ PARAMETERS TO SET ############
## list of SpOCK main input files
list_run = ["41884.txt",
            "41885.txt",
            "41886.txt",
            "41887.txt",
            "41888.txt"]
            

## Save or not the plots
save_results = 1

## Show or not the plots
show_plots = 1

## path of the folder where you want to store the results (pickle, image, video)
path_folder_results = '/raid3/Armada/Charles/python/' #get_prop_dir(2) + 'output/python_propagator/'

## If the second spacecraft was propagated from the same main input file in SpOCK as the first spacecraft
same_spock_input_file = 0

############ ALGORITHM ############
nb_sat = len(list_run)
for isat in range(nb_sat):
# Read input file sat1
input_filename_sat1 =    get_prop_dir(1) +  'run.cygnss/input/main_input/' + list_run[isat]
input_variables, order_input_variables = read_input_file(input_filename_sat1)
sat1_to_plot_path = input_variables[6][0]; sat1_to_plot = input_variables[7][0]
dt_sat1 = input_variables[2]; nb_steps_sat1 = input_variables[3]; 

name_mission = "cygnss"

root_save_fig_name = path_folder_results + name_mission + '/result/image/' + sat1_to_plot.replace(".txt","_")

# Read output file sat1
r_eci_sat1 = np.zeros([nb_steps_sat1])

var_to_read = ["position", "altitude"]
var_out, var_out_order = read_output_file( sat1_to_plot_path + sat1_to_plot, var_to_read )
r_eci_sat1 = var_out[1]
date_start = datetime.strptime(var_out[0][0], "%Y/%m/%d %H:%M:%S")

# Read input file sat2
if same_spock_input_file == 1:
    sat2_to_plot_path = input_variables[6][1]; sat2_to_plot = input_variables[7][1]
    dt_sat2 = dt_sat1; nb_steps_sat2 = nb_steps_sat1
else:
    input_filename_sat2 =    get_prop_dir(1) + sys.argv[1] + '/input/main_input/' + sys.argv[3] 
    input_variables, order_input_variables = read_input_file(input_filename_sat2)
    sat2_to_plot_path = input_variables[6][0]; sat2_to_plot = input_variables[7][0]
    dt_sat2 = input_variables[2]; nb_steps_sat2 = input_variables[3]; 

if nb_steps_sat2 != nb_steps_sat1:
    print "***! Error: the number of time steps in both sumations is different. The program will stop. !***"
    raise Exception

# Read output file sat2
r_eci_sat2 = np.zeros([nb_steps_sat2])
var_to_read = ["position"]
var_out, var_out_order = read_output_file( sat2_to_plot_path + sat2_to_plot, var_to_read )
r_eci_sat2 = var_out[1]

# Distance between both sc
dist_between_sat1_and_sat2 = np.zeros([nb_steps_sat2])
for i in range(nb_steps_sat2):
    dist_between_sat1_and_sat2[i] = np.linalg.norm( r_eci_sat2[i,:] - r_eci_sat1[i,:] )
print min(dist_between_sat1_and_sat2)
