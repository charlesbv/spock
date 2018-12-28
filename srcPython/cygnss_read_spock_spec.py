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
# This function reads the specular point position files output by SpOCK. It returns the time/lon/lat/gain of the specular points as well as the GPS corresponding to the specular point. Also the NormPower of the spec, xyz of CYGNSS, xyz of GPS, xyz spec (ECEF, km)
# The arguments of the function are
# - filename_spec: the name of the specular file to read (path included)
# lon, lat, gain, gps are lists. For each list, each element of the list correspond to one time so it has nb_spec elements (where nb_spec is the number of specular points at that time (4 usually, if not always))


import os
from read_input_file import *
from datetime import datetime, timedelta


def cygnss_read_spock_spec(filename_spec):
    #filename_spec = "run.spock/spock/spock_spec_start_2017-03-29T22_05_22_end_2017-03-29T23_54_11/spock_spec_start_2017-03-29T22_05_22_end_2017-03-29T23_54_118/specular_spock_spec_start_2017-03-29T22_05_22_end_2017-03-29T23_54_118.txt"

    # READ THE POSITIONS OF THE SPECULAR POINTS, THE GPS AND THE CYGNSS SATELLITES
    nb_spec_pts = 4;
    interpolation_step = 1 # in second, interpolation step of find_specular_points.c (1 s usually)

    count_less_than_4_spec = 0
    time_less_than_4_spec = []
    nb_spec_this_time = []
    file_specular = open(filename_spec)
    read_file_specular  = file_specular.readlines()
    # Nb of lines in the spec file header
    nb_lines_header_output_file_spec = 0
    while (read_file_specular[nb_lines_header_output_file_spec].split()[0] != "#START"):
        nb_lines_header_output_file_spec = nb_lines_header_output_file_spec + 1
    nb_lines_header_output_file_spec = nb_lines_header_output_file_spec + 1

    ispec_save = 0
    j = -1
    lon = []; lat = []; gain = []; gps = []; date = []
    normpower = []; x_cyg = []; y_cyg = []; z_cyg = []; x_gps = []; y_gps = []; z_gps = []
    x_spec = []; y_spec = []; z_spec = []

    while (ispec_save < len(read_file_specular)-1-nb_spec_pts):
        lon_this_time = []
        lat_this_time = []
        gain_this_time = []
        gps_this_time = []
        normpower_this_time = []
        x_cyg_this_time = []
        y_cyg_this_time = []
        z_cyg_this_time = []
        x_gps_this_time = []
        y_gps_this_time = []
        z_gps_this_time = []
        x_spec_this_time = []
        y_spec_this_time = []
        z_spec_this_time = []


        j = j + 1
        time_spec_sublist_temp_ini = read_file_specular[ispec_save+nb_lines_header_output_file_spec].split()[0]
        time_spec_sublist_temp_ini = datetime.strptime(time_spec_sublist_temp_ini, "%Y-%m-%dT%H:%M:%S")
        date.append(time_spec_sublist_temp_ini)
        lon_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[1]) )
        lat_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[2]) )
        gain_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[3]) )
        gps_this_time.append( read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[4] )
        normpower_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[5]) )
        x_cyg_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[6]) )
        y_cyg_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[7]) )
        z_cyg_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[8]) )
        x_gps_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[9]) )
        y_gps_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[10]) )
        z_gps_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[11]) )
        x_spec_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[12]) )
        y_spec_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[13]) )
        z_spec_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[14]) )

        ispec = 1
        while (datetime.strptime(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[0], "%Y-%m-%dT%H:%M:%S")  == time_spec_sublist_temp_ini):
            lon_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[1]) )
            lat_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[2]) )
            gain_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[3]) )
            gps_this_time.append( read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[4] )
            normpower_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[5]) )
            x_cyg_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[6]) )
            y_cyg_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[7]) )
            z_cyg_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[8]) )
            x_gps_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[9]) )
            y_gps_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[10]) )
            z_gps_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[11]) )
            x_spec_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[12]) )
            y_spec_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[13]) )
            z_spec_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[14]) )


            ispec = ispec + 1
            if (nb_lines_header_output_file_spec+ispec_save+ispec == len(read_file_specular)):
                break
            if (len(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()) == 0): # sometimes the spec file ends with two blank lines
                break
        if ispec < 4:
            count_less_than_4_spec = count_less_than_4_spec + 1
            time_less_than_4_spec.append(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec-1].split()[0])
            nb_spec_this_time.append(ispec)
        ispec_save = ispec + ispec_save
        lon.append(lon_this_time)
        lat.append(lat_this_time)
        gain.append(gain_this_time)
        gps.append(gps_this_time)
        normpower.append(normpower_this_time)
        x_cyg.append(x_cyg_this_time)
        y_cyg.append(y_cyg_this_time)
        z_cyg.append(z_cyg_this_time)
        x_gps.append(x_gps_this_time)
        y_gps.append(y_gps_this_time)
        z_gps.append(z_gps_this_time)
        x_spec.append(x_spec_this_time)
        y_spec.append(y_spec_this_time)
        z_spec.append(z_spec_this_time)

    file_specular.close()
    return date, lon, lat, gain, gps, normpower, x_cyg, y_cyg, z_cyg, x_gps, y_gps, z_gps,  x_spec, y_spec, z_spec



