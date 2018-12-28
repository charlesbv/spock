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
import numpy as np

def compute_power(filename):
    #filename = '/home/cbv/PropSim/output/run_aerie_power_only_tail_and_normalized/aerie_power_only_tail_and_normalized/power_aerie_power_only_tail_and_normalized.txt'
    file = open(filename, "r")
    read_file = file.readlines()
    # SKIP HEADER
    n_header = 0
    while (read_file[n_header].split()[0] != '#START'):
        n_header = n_header + 1
    n_header = n_header + 1
    n = len(read_file) - n_header
    nb_surfaces = len(read_file[n_header].split()) - 3 # -3 to skip the time and the Sun elevation angle
    power = np.zeros([n, nb_surfaces])
    sun_elevation = np.zeros([n])
    for istep in range(1,n):
        for isurf in range(nb_surfaces):
            power[istep, isurf] = read_file[istep + n_header].split()[isurf + 2]
        sun_elevation[istep] = read_file[istep + n_header].split()[nb_surfaces+1 + 2]
    return power, sun_elevation
