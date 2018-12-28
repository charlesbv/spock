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
import fileinput
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
import numpy as np
import sys
from read_input_file import *
from read_output_file import *
from get_prop_dir import *

plt.ion()

# File with the results of the difference in positions and velocities 
# file_results = open("./paper3_stk_results/results_error_forces.txt", "w+")
# print >> file_results, "name_run max_diff_pos max_diff_vel (distances are in kilometers)\n"

input_filename_leo_array = ['iss','iss_msis', 'iss_gravity_j2','iss_gravity','iss_third_body', 'iss_solar_pressure']
forces_arr = ['2body','2body_msis', '2body_j2','2body_gravity2020','2body_third_body', '2body_solar_pressure']
altitude_arr = [300,600,1000,5000,10000,20000,35786]

# # # Runs
# for ialt in range(len(altitude_arr)):
#     for iforce in range(len(forces_arr)):
#         print ialt, len(altitude_arr) - 1
#         print iforce, len(forces_arr) - 1
#         ## Write input file
#         os.chdir("../run_paper3/")
#         input_filename_propagator = forces_arr[iforce] + '_alt_' + str(altitude_arr[ialt] )
#         os.system("cp input/main_input/" + input_filename_leo_array[iforce] + ".txt " + "input/main_input/" + input_filename_propagator + ".txt")
#         line_number = 0
#         for line in fileinput.input("input/main_input/" + input_filename_propagator + ".txt", inplace = True):
#             line_number = line_number + 1
#             if ( ( line_number != 25 ) & ( line_number != 15 ) ):
#                 print "%s" % (line),
#             elif ( line_number == 15 ):
#                 string_here = str(altitude_arr[ialt]) + " 0 0 0 0 0.00000000000000000000000001"
#                 print "%s\n" % string_here,
#             elif ( line_number == 25 ):
#                 string_here = input_filename_propagator
#                 print "%s\n" % string_here,
#         ## Run propagator
#         os.system("/usr/local/bin/mpirun -np 1 spock " + input_filename_propagator + ".txt")
# os.chdir("../srcPython")

max_r_diff_mag = np.zeros([len(forces_arr)-1, len(altitude_arr)]) # -1 because we don't care about 2body difference with 2body
max_v_diff_mag = np.zeros([len(forces_arr)-1, len(altitude_arr)])
##### BEGINNING OF LOOP

for ialt in range(len(altitude_arr)):
    for iforce in range(len(forces_arr)):
        # Read position and velocity from our propagator
        input_filename_prop = forces_arr[iforce] + '_alt_' + str(altitude_arr[ialt] ) + '.txt'
        input_filename_prop_complete = get_prop_dir(1) + "run_paper3/input/main_input/"  + input_filename_prop
        input_variables, order_input_variables = read_input_file(input_filename_prop_complete)
        output_propagator_path = input_variables[6][0]; output_propgator_file = input_variables[7][0]
        to_output = ["position", "velocity"]
        out, out_var = read_output_file(output_propagator_path + output_propgator_file, to_output)
        if iforce == 0: # 2 body
            r_eci_2body = out[1]; v_eci_2body = out[2]
        else:

            # Difference between positions and velocities
            if '2020' in input_filename_prop: # we compare 20, 20 order gravity model to j2 gravity model (not to 2 body). ASSUMPTION: in the list with the names of the runs, the 20,20 gravity model run must be right after the j2 run
                r_eci_prop_all_gravity = out[1]; v_eci_prop_all_gravity = out[2]
                r_diff = r_eci_prop - r_eci_prop_all_gravity # reference is not 2 body but j2
                v_diff = v_eci_prop - v_eci_prop_all_gravity # reference is not 2 body but j2

            else:
                r_eci_prop = out[1]; v_eci_prop = out[2]
                r_diff = r_eci_2body - r_eci_prop
                v_diff = v_eci_2body - v_eci_prop

            for i in range(len(out[0])):
                r_diff_mag = np.linalg.norm(r_diff[i]) # in km
                v_diff_mag = np.linalg.norm(v_diff[i]) # in km/s

            max_r_diff_mag[iforce-1, ialt] = np.max(r_diff_mag)
            max_v_diff_mag[iforce-1, ialt] = np.max(v_diff_mag)

            # # Write results in file
            # print >> file_results, input_filename_prop.split('.')[0] + ' ' + '{0:.3f}'.format(np.max(r_diff_mag)) + ' ' + '{0:.3f}'.format(np.ma
               

    ##### END OF LOOP

# # Close file with results
# file_results.close()

# PLOT
width_fig = 15.
height_fig = width_fig * 3 /4
fontsize_plot = 20 # 9
color_arr = ['k','b','r','g','m', 'y']
fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.05, right=0.97, top = 0.90,bottom = 0.06)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(gs[0, 0])
#ax1 = fig.add_subplot(111)
ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
#ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
label_arr = ['Drag', 'J2','Spherical harmonics - 20','Third body', 'Solar pressure']
marker_arr = ['o', '^', 's', '*', 'D']

for iforce in range(len(forces_arr)-1): # go through the iss. meo, and geo orbits
    if marker_arr[iforce] == '*':
        ax1.semilogy(altitude_arr, max_r_diff_mag[iforce,:], 'k', marker = marker_arr[iforce], linestyle = '-', linewidth = 2, label = label_arr[iforce], color = color_arr[iforce], markersize = 15 )
    else:
        ax1.semilogy(altitude_arr, max_r_diff_mag[iforce,:], 'k', marker = marker_arr[iforce], linestyle = '-', linewidth = 2, label = label_arr[iforce], color = color_arr[iforce], markersize = 10 )

    ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)


    ax1.set_title('Error as a function of altitude', weight = 'bold', fontsize = fontsize_plot , y = 1.005) 
    ax1.legend( fontsize = fontsize_plot, loc = 4)
    ax1.set_xlabel('Altitude (km)', fontsize = fontsize_plot, weight = 'bold')
    ax1.set_ylabel('Error (km)', fontsize = fontsize_plot, weight = 'bold')

raise Exception
fig_save_name = './paper3_stk_results/error_forces_semilogy.eps'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
os.system("rsync -av " + fig_save_name + " srbwks2014-0008.engin.umich.edu:./")


fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.05, right=0.97, top = 0.90,bottom = 0.06)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(gs[0, 0])
#ax1 = fig.add_subplot(111)
ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
#ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
label_arr = ['Drag', 'J2','Spherical harmonics - 20','Third body', 'Solar pressure']
marker_arr = ['o', '^', 's', '*', 'D']

for iforce in range(len(forces_arr)-1): # go through the iss. meo, and geo orbits
    if marker_arr[iforce] == '*':
        ax1.semilogy(altitude_arr, max_r_diff_mag[iforce,:], 'k', marker = marker_arr[iforce], linestyle = '-', linewidth = 2, label = label_arr[iforce], color = color_arr[iforce], markersize = 15 )
    else:
        ax1.semilogy(altitude_arr, max_r_diff_mag[iforce,:], 'k', marker = marker_arr[iforce], linestyle = '-', linewidth = 2, label = label_arr[iforce], color = color_arr[iforce], markersize = 10 )
        
    ax1.set_xlim([300, 2000])
    ax1.margins(1,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)


    ax1.set_title('Error as a function of altitude', weight = 'bold', fontsize = fontsize_plot , y = 1.005) 
    ax1.legend( fontsize = fontsize_plot, loc = 4)
    ax1.set_xlabel('Altitude (km)', fontsize = fontsize_plot, weight = 'bold')
    ax1.set_ylabel('Error (km)', fontsize = fontsize_plot, weight = 'bold')


fig_save_name = './paper3_stk_results/error_forces_semilogy_zoom_under_2000km.eps'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
os.system("rsync -av " + fig_save_name + " srbwks2014-0008.engin.umich.edu:./")
