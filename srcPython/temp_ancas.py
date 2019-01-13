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
from read_output_file import *
from orbit_average import *

plt.ion()
plt.isinteractive()
# !!!!!!!!!! SET THE PARAMETER BELOW: path of the folder where you want to store the results (pickle, image, video)
path_folder_results = '/raid3/Armada/Charles/python/' #get_prop_dir(2) + 'output/python_propagator/'

main_input_file_name = get_prop_dir(1) + sys.argv[1] + '/input/main_input/' + sys.argv[2]  # ex of sys.argv[1]: cadre_tumbling_0panel.txt'


# read input file
input_variables, order_input_variables = read_input_file(main_input_file_name)
dt = input_variables[2]; nb_steps = input_variables[3]; satellite_to_plot_path = input_variables[6][0]; satellite_to_plot = input_variables[7][0]; nb_ensembles_coe = input_variables[12]; nb_ensembles_attitude = input_variables[13]; nb_ensembles_cd = input_variables[11]; ensemble_to_plot_temp = input_variables[14]; 
nb_ensembles_density = input_variables[17]


# read output file of reference satellite
list_var = ["position", "velocity", "acceleration"]
out_var, order_var = read_output_file( satellite_to_plot_path + satellite_to_plot, list_var  )
date_main = out_var[0]
r_eci = out_var[1]
v_eci = out_var[2]
a_eci = out_var[3]


# Analytical solution acceleration
earth_flattening    = 1/298.257223560; # Earth flattening coefficient (no unit)
earth_radius        = 6378.137; # mean equatorial radius (km)
earth_mu    = 398600.4418; # gravitational parameter (km^3/s^2)
earth_j2    = 1.081874e-3; # J2 zonal harmonic coefficient (no unit)
rad2deg = 180./np.pi
deg2rad = 1 / rad2deg
second2year = 3600 * 24 * 365.25

r_eci_mag = np.zeros([nb_steps])
for istep in range(nb_steps):
    r_eci_mag[istep] = np.linalg.norm(r_eci[istep, :])

a_eci_analytical = np.zeros([nb_steps, 3])
a_eci_analytical[:, 0] = -earth_mu * r_eci[:,0] / r_eci_mag**3  - 3 * earth_mu * earth_j2 * earth_radius**2 * r_eci[:,0] * ( 1 - 5 * r_eci[:,2]**2 / r_eci_mag**2 ) / ( 2 * r_eci_mag**5 )

a_eci_analytical[:, 1] = -earth_mu * r_eci[:,1] / r_eci_mag**3  - 3 * earth_mu * earth_j2 * earth_radius**2 * r_eci[:,1] * ( 1 - 5 * r_eci[:,2]**2 / r_eci_mag**2 ) / ( 2 * r_eci_mag**5 )

a_eci_analytical[:, 2] = -earth_mu * r_eci[:,2] / r_eci_mag**3  - 3 * earth_mu * earth_j2 * earth_radius**2 * r_eci[:,2] * ( 1 - 5 * r_eci[:,2]**2 / r_eci_mag**2 ) / ( 2 * r_eci_mag**5 )

