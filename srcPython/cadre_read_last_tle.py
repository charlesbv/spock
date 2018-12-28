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
from read_input_file import *
from read_output_file import *
import sys

def cadre_read_last_tle(main_input_file_name):
    #main_input_file_name = '/home/cbv/PropSim/input/main_input/cadre_last_tle.txt' 

    # read input file
    input_variables, order_input_variables = read_input_file(main_input_file_name)
    dt = input_variables[2]; nb_steps = input_variables[3]; satellite_to_plot_path = input_variables[6][0]; satellite_to_plot = input_variables[7][0]; nb_ensembles_coe = input_variables[12]; nb_ensembles_attitude = input_variables[13]; nb_ensembles_cd = input_variables[11]; ensemble_to_plot_temp = input_variables[14]; 
    n = nb_steps

    # read input file
    list_output_variables_to_read = ["position"]
    output_filename = satellite_to_plot_path + satellite_to_plot
    output_variables, list_output_variables_read = read_output_file(output_filename, list_output_variables_to_read)

    eci_last_tle = output_variables[1][0]

    return eci_last_tle
