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

# This script converts an obsevation file from ECEF to ECI coordinates and write thee output in a file
import numpy as np
import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/srcPython")
import os
from ecef2eci import *


def convert_cygnss_obs_ecef_to_eci(obs_rv_filename):
    # Read r/v of observations
    #obs_rv_filename = 'HD_data/spock_FM5_20171216_eng_pvt_query-13527_shorter.txt'
    obs_rv_filename_eci = obs_rv_filename.replace(".txt", "_eci.txt")
    obs_rv_file = open(obs_rv_filename)
    obs_rv_file_eci = open(obs_rv_filename_eci, "w")
    read_obs_rv_file = obs_rv_file.readlines()
    nb_header = 0
    while (read_obs_rv_file[nb_header].split()[0] != '#START'):
        print >> obs_rv_file_eci, read_obs_rv_file[nb_header].replace("\n","").replace("\r","")
        nb_header = nb_header + 1
    print >> obs_rv_file_eci, read_obs_rv_file[nb_header].replace("\n","").replace("\r","")
    nb_header = nb_header + 1
    nb_obs = len(read_obs_rv_file) - nb_header
    for iobs in range(nb_obs):
        date_obs_str = read_obs_rv_file[iobs + nb_header].split()[0] 
        r_obs_ecef = np.zeros([3])
        v_obs_ecef = np.zeros([3])
        r_obs_ecef[0] = np.float( read_obs_rv_file[iobs + nb_header].split()[1] ) 
        r_obs_ecef[1] = np.float( read_obs_rv_file[iobs + nb_header].split()[2] ) 
        r_obs_ecef[2] = np.float( read_obs_rv_file[iobs + nb_header].split()[3] ) 
        v_obs_ecef[0] = np.float( read_obs_rv_file[iobs + nb_header].split()[4] ) 
        v_obs_ecef[1] = np.float( read_obs_rv_file[iobs + nb_header].split()[5] ) 
        v_obs_ecef[2] = np.float( read_obs_rv_file[iobs + nb_header].split()[6] ) 
    # Convert ECEF to ECI
        if iobs == 0:
            r_obs, v_obs = ecef2eci(r_obs_ecef, v_obs_ecef, date_obs_str, 1)
        else:
            r_obs, v_obs = ecef2eci(r_obs_ecef, v_obs_ecef, date_obs_str, 0)

        print >> obs_rv_file_eci, date_obs_str, format(r_obs[0], ".14e"), format(r_obs[1], ".14e"), format(r_obs[2], ".14e"), format(v_obs[0], ".14e"), format(v_obs[1], ".14e"), format(v_obs[2], ".14e")

    obs_rv_file_eci.close()
    return obs_rv_filename_eci
