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
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/kalman/spock_development_new_structure_kalman_dev/srcPython")
from find_in_read_input_order_variables import *
from get_prop_dir import *
import fileinput
import time
import datetime
import numpy as np
import pylab as P
import os
import subprocess
from read_input_file import *

## NOTE 1: to use this script, the only parameter you have to set is the name of the main input file (first non-commented line below)


main_input_file_name =  sys.argv[1]

# read input file
input_variables, order_input_variables = read_input_file(main_input_file_name)
dt = input_variables[find_in_read_input_order_variables(order_input_variables, 'dt')];
nb_steps = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_steps')];
satellite_to_plot_path = input_variables[find_in_read_input_order_variables(order_input_variables, 'output_file_path_list')];
satellite_to_plot = input_variables[find_in_read_input_order_variables(order_input_variables, 'output_file_name_list')];
nb_ensembles_coe = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_ensembles_coe')];
nb_ensembles_attitude = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_ensembles_attitude')];
nb_ensembles_cd = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_ensembles_cd')];
ensemble_to_plot_temp = input_variables[find_in_read_input_order_variables(order_input_variables, 'ensembles_to_output')];
nb_ensembles_density = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_ensembles_density')]
nb_spacecraft = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_sc')]
compute_drag = input_variables[find_in_read_input_order_variables(order_input_variables, 'compute_drag')]

# Ensembles created by the propagator
ensemble_to_plot = []
for i in range(len(ensemble_to_plot_temp)):
    if (ensemble_to_plot_temp[i] == 'eci_r'):
        ensemble_to_plot.append('x_eci'); ensemble_to_plot.append('y_eci'); ensemble_to_plot.append('z_eci')
    if (ensemble_to_plot_temp[i] == 'eci_v'):
        ensemble_to_plot.append('vx_eci'); ensemble_to_plot.append('vy_eci'); ensemble_to_plot.append('vz_eci')
    if (ensemble_to_plot_temp[i] == 'geodetic'):
        ensemble_to_plot.append('longitude'); ensemble_to_plot.append('latitude'); ensemble_to_plot.append('altitude')
    if (ensemble_to_plot_temp[i] == 'power'):
        ensemble_to_plot.append('power')
    if (ensemble_to_plot_temp[i] == 'attitude'):
        ensemble_to_plot.append('pitch'); ensemble_to_plot.append('roll'); ensemble_to_plot.append('yaw')
    if (ensemble_to_plot_temp[i] == 'oe'):
        ensemble_to_plot.append('sma'); ensemble_to_plot.append('inclination'); ensemble_to_plot.append('eccentricity'); ensemble_to_plot.append('true_anomaly'); ensemble_to_plot.append('RAAN'); ensemble_to_plot.append('argument_perigee');
    if ((ensemble_to_plot_temp[i] == 'density') & (compute_drag == 1)):
        ensemble_to_plot.append('rho'); ensemble_to_plot.append('f107'); ensemble_to_plot.append('f107a'); ensemble_to_plot.append('ap'); 
    if ((ensemble_to_plot_temp[i] == 'cd') & (compute_drag == 1)):
        ensemble_to_plot.append('cd')


## Nb of ensembles
nb_ensembles_array = [nb_ensembles_coe, nb_ensembles_attitude, nb_ensembles_cd, nb_ensembles_density]
nb_ensembles = np.max(nb_ensembles_array)
for i in range(len(nb_ensembles_array)):
    if ( (nb_ensembles_array[i] > 0 ) &  (nb_ensembles_array[i] < nb_ensembles) ) :
        nb_ensembles = nb_ensembles_array[i]

## number of processors
for isat in range(nb_spacecraft):
    nb_ensemble_to_plot = len(ensemble_to_plot)
    nProc = 0
    for file in os.listdir(satellite_to_plot_path[isat] + "ensemble/"):
        if file[0:5] == "iproc":
            nProc = (int)(file.split('-')[1].split('_')[0]) 
    if ( nProc == 0 ):
        print "The concatenation of the processors files has probably been already done. Please check in: " + satellite_to_plot_path[isat] + "ensemble/. If it is not the case then please erase all files in: " + satellite_to_plot_path[isat] + "ensemble/ and run the propagator again (and blame me for that)."
        raise Exception
    nb_ensemble_per_proc = int(np.floor(nb_ensembles / nProc))
    nb_ensembles_ok = nb_ensemble_per_proc * nProc

    # Create folder for iproc files
    if os.path.isdir(satellite_to_plot_path[isat]+"ensemble/proc_files") == False:
        os.system("mkdir "+satellite_to_plot_path[isat]+"ensemble/proc_files &>/dev/null")       
    else:
        os.system("rm -Rf " + satellite_to_plot_path[isat]+"ensemble/proc_files/* &>/dev/null")

    # Create folder for final output ensemble files
    dir_final_output_ensemble = satellite_to_plot_path[isat]+"ensemble/"

# gather all output files from different processors
    for eee in range(nb_ensemble_to_plot):
        for i in range(nProc):
            name_file_iproc = satellite_to_plot_path[isat] + 'ensemble/iproc_'+str(i+1)+'-'+str(nProc)+'_'+ensemble_to_plot[eee]+'_'+satellite_to_plot[isat]
            file_iproc = open(name_file_iproc,'r') 
            a = file_iproc.readlines()
            nb_steps = len(a) - 11 # !!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
            if (i == 0):
                line_iproc = ["" for x in range(nb_steps)]
                x_alltime_ensemble = np.zeros([nb_steps, nb_ensembles_ok])
                time_iproc = ["" for x in range(nb_steps)]
                
            for line_count in range(nb_steps):
                line_iproc[line_count] = (a[11+line_count]).split(' ')
                if (i == 0):
                    time_iproc[line_count] = line_iproc[line_count][0] + ' ' + line_iproc[line_count][1]
                for ensemble_temp in range(nb_ensemble_per_proc):
                    ensemble = ensemble_temp + nb_ensemble_per_proc * i
                    x_alltime_ensemble[line_count][ensemble] = float( line_iproc[line_count][ensemble_temp+2]  )
                    
            file_iproc.close()
            subprocess.call(['mv', name_file_iproc, satellite_to_plot_path[isat] + 'ensemble/proc_files/'])

        both_time_ensemble = ["" for x in range(nb_steps)]
        for i in range(nb_steps):
            str_ensemble = ' '
            for j in range(nb_ensembles_ok):
                str_ensemble = str_ensemble +  str(x_alltime_ensemble[i][j]) + ' '
            both_time_ensemble[i] = time_iproc[i] + str_ensemble

        filename = 'ensemble_'+ensemble_to_plot[eee]+'_'+satellite_to_plot[isat]
        header_file = ''
        for i in range(11):
            if i == 10:
                a[i] = a[i].replace("\n","")
            header_file = header_file + a[i]
        np.savetxt(dir_final_output_ensemble+filename,both_time_ensemble, fmt='%s',header = header_file,comments = '')


