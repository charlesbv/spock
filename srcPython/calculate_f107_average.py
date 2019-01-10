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

import datetime
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os # to exectute terminal commands from python
from  get_prop_dir import *
# plt.ion()
# plt.isinteractive()

def calculate_f107_average(input_filename):
    # input file of propagator
#    input_filename = "../input_armada.d"
    input_filename = get_prop_dir(2) + "input/density/density_NRLMSIS00e/" + input_filename
    file_input = open(input_filename, "r")
    read_file_input = file_input.readlines()
    date_start = read_file_input[1].split()[0]
    date_start = datetime.datetime(*time.strptime(date_start, "%Y-%m-%dT%H:%M")[0:6])
    date_stop = read_file_input[2].split()[0]
    date_stop = datetime.datetime(*time.strptime(date_stop, "%Y-%m-%dT%H:%M")[0:6])
    doy = date_start.strftime('%j')
    dt = np.float(read_file_input[3].split()[0])
    output_file_propagator = read_file_input[11].split(',')[0]
    density_type = read_file_input[53].split()[0]
    if (density_type == "density_file"):
        raid3 = read_file_input[54].split()[0]
    elif (density_type == "dynamic"):
        f107_file = read_file_input[54].split()[0]
        ap_file = read_file_input[56].split()[0]
    elif (density_type == "static"):
        f107 = read_file_input[54].split()[0]
        f107A = read_file_input[55].split()[0]
        ap = read_file_input[56].split()[0]


    f107_file_open = open("../input/density/density_NRLMSIS00e/" + f107_file, "r")
    f107_file_read = f107_file_open.readlines()
    file_output_name = f107_file.split('.')[0]+"_average."+f107_file.split('.')[1]
    file_output = open("../input/density/density_NRLMSIS00e/"+ file_output_name, "w+")
    print >> file_output, "#BEGINNINGOFHEADER\n#ENDOFHEADER\nYEAR DOY HR    1"
    # SKIP HEADER OF FILE F10.7
    nb_lines_header_f107_file = 0
    if (f107_file_read[0].split()[0] == '#BEGINNINGOFHEADER'):
        i = 0
        while ( (f107_file_read[i].split()[0] != "YEAR") & (f107_file_read[i].split()[0] != "YYYY-MM-DD") ):
            nb_lines_header_f107_file = nb_lines_header_f107_file + 1
            i = i + 1
        nb_lines_header_f107_file = nb_lines_header_f107_file + 1
    if ( f107_file_read[-1].split()[0] == "#ENDOFFILE" ):
        n = len(f107_file_read) - nb_lines_header_f107_file - 1
    else:
        n = len(f107_file_read) - nb_lines_header_f107_file 
    f107A = np.zeros(n)
    f107 = np.zeros(n)
    # READ FILE F10.7 VALUES 
    for i in range(n):
        date_in = f107_file_read[i+nb_lines_header_f107_file].split()[0] + ' ' + f107_file_read[i+nb_lines_header_f107_file].split()[1] + ' '  + f107_file_read[i+nb_lines_header_f107_file].split()[2]
        date_in = datetime.datetime.strptime(date_in, '%Y %j %H')
        f107[i] = float(f107_file_read[i+nb_lines_header_f107_file].split()[3])
        if ( ( date_in.strftime('%y') == date_start.strftime('%y') ) & ( date_in.strftime('%m') == date_start.strftime('%m') ) & ( date_in.strftime('%d') == date_start.strftime('%d') ) & ( date_in.strftime('%H') == date_start.strftime('%H') ) ): # calculate the 3 months average only for the period of time of the propagation
            k = 0
            which_time = []
            while (( date_in.strftime('%y') != date_stop.strftime('%y') ) | ( date_in.strftime('%m') != date_stop.strftime('%m') ) | ( date_in.strftime('%d') != date_stop.strftime('%d') ) | ( date_in.strftime('%H') != date_stop.strftime('%H') )):
                date_in_asinfile = f107_file_read[i+nb_lines_header_f107_file+k].split()[0] + ' ' + f107_file_read[i+nb_lines_header_f107_file+k].split()[1] + ' '  + f107_file_read[i+nb_lines_header_f107_file+k].split()[2]
                date_in = datetime.datetime.strptime(date_in_asinfile, '%Y %j %H')

            # Calculate the 3 months average of f10.7
                f107_array_40d_before = np.zeros(40*24)
                f107_array_40d_after = np.zeros(40*24)
                for j in range(1,40*24L):
                    f107_array_40d_before[j] = float(f107_file_read[i+nb_lines_header_f107_file+k-j].split()[3])
                    f107_array_40d_after[j] = float(f107_file_read[i+nb_lines_header_f107_file+k+j].split()[3])
                f107A[k] = np.mean( np.concatenate([ f107_array_40d_before[1:-1], np.array([float(f107_file_read[i+nb_lines_header_f107_file+k].split()[3])]), f107_array_40d_after[1:-1] ]) )
                which_time.append(i+k)
                print >> file_output, date_in_asinfile, f107A[k]
                k = k + 1
    print >> file_output,"#ENDOFFILE"
    file_output.close()

    if (nb_lines_header_f107_file == 0):
        # Adding the header in the original F10.7 file IF NECESSARY (so if the file has just be downloaded from omniweb)
        with open("../input/density/density_NRLMSIS00e/" + f107_file) as f:
            with open("../input/density/density_NRLMSIS00e/" + f107_file + "copy", "w") as f1:
                print >> f1, "#BEGINNINGOFHEADER\n#ENDOFHEADER\nYEAR DOY HR    1"
                for line in f:
                    f1.write(line) 
                print >> f1, "#ENDOFFILE"
        os.rename("../input/density/density_NRLMSIS00e/" + f107_file + "copy", "../input/density/density_NRLMSIS00e/" + f107_file)

    nb_lines_header_ap_file = 0
    ap_file_open = open("../input/density/density_NRLMSIS00e/" + ap_file, "r")
    ap_file_read = ap_file_open.readlines()
    if (ap_file_read[0].split()[0] == '#BEGINNINGOFHEADER'):
        i = 0
        while ( (ap_file_read[i].split()[0] != "YEAR") & (ap_file_read[i].split()[0] != "YYYY-MM-DD") ):
            nb_lines_header_ap_file = nb_lines_header_ap_file + 1
            i = i + 1
        nb_lines_header_ap_file = nb_lines_header_ap_file + 1
    if (nb_lines_header_ap_file == 0):
    # Adding the header in the original Ap file IF NECESSARY (so if the file has just be downloaded from omniweb)
        with open("../input/density/density_NRLMSIS00e/" + ap_file) as f:
            with open("../input/density/density_NRLMSIS00e/" + ap_file + "copy", "w") as f1:
                print >> f1, "#BEGINNINGOFHEADER\n#ENDOFHEADER\nYEAR DOY HR    1"
                for line in f:
                    f1.write(line) 
                print >> f1, "#ENDOFFILE"
        os.rename("../input/density/density_NRLMSIS00e/" + ap_file + "copy", "../input/density/density_NRLMSIS00e/" + ap_file )


    return

# n_prop = len(which_time)
# f107_during_prop = np.zeros(n_prop)
# f107_during_prop = f107[which_time[0]-40*24L:which_time[n_prop-1]+40*24L]

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.plot(range(40*24L*2+24),f107_during_prop)
# ax.plot(range(40*24L,41*24L),f107A[0:n_prop-1])


