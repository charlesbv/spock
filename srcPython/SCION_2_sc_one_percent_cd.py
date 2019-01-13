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

import sys
from get_prop_dir import *
from read_input_file import *
from read_output_file import *
from matplotlib import pyplot as plt

plt.ion()

if len(sys.argv) > 2:
    main_input_file_name = get_prop_dir(1) + sys.argv[1] + '/input/main_input/' + sys.argv[2]  # ex of sys.argv[1]: cadre_tumbling_0panel.txt'
else:
    main_input_file_name = get_prop_dir(1) + 'run/input/main_input/' + sys.argv[1]  # ex of sys.argv[1]: cadre_tumbling_0panel.txt'    

# read input file
input_variables, order_input_variables = read_input_file(main_input_file_name)
dt = input_variables[2]; nb_steps = input_variables[3]; satellite_to_plot_path = input_variables[6]; satellite_to_plot = input_variables[7];
nb_satellites = input_variables[4]

# read output files
var_out_choose = ["position", "velocity"]
r = np.zeros([nb_satellites, nb_steps, 3])
v = np.zeros([nb_satellites, nb_steps, 3])
for isat in range(nb_satellites):
    var_out, var_out_name = read_output_file(satellite_to_plot_path[isat] + satellite_to_plot[isat], var_out_choose)
    r[isat, :, :] = var_out[1]
    v[isat, :, :] = var_out[2]

# Distances between the 2 sc
dist = np.zeros([nb_steps])
for istep in range(nb_steps):
    dist[istep] = np.linalg.norm(r[1, istep] - r[0, istep])

# Plot
fig, ax = plt.subplots()
ax.plot(dist)
