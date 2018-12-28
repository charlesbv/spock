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
# This script reads a ECI output file from SpOCK and adds gaussian noise to each time step
# To run this script:
# python noise.py spock_main_input
# where spock_main_input is the main input file that was sued to initialize SpOCK to output the ECI output file on which we add noise
# ASSUMPTIONS:
# - see section "PARAMETERS TO SET UP BEFORE RUNNING THE SCRIPT" at beginning of this script

# PARAMETERS TO SET UP BEFORE RUNNING THE SCRIPT
sigma_noise_r = [10., 10., 10.] # in mm, std for each component of position
sigma_noise_v = [0.01, 0.01, 0.01] # in m/s, std for each component of velocity
sigma_noise_r0 = [1., 1., 1.] # in m, std for each component of initial position
sigma_noise_v0 = [0.01, 0.01, 0.01] # in m/s, std for each component of initial velocity

show_plot = 1 # set it to 1 to plot the true and noise orbits
# end of PARAMETERS TO SET UP BEFORE RUNNING THE SCRIPT

from datetime import datetime, timedelta
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import sys
from read_input_file import *
from read_output_file import *

sigma_noise_r = np.array(sigma_noise_r)/1000.
sigma_noise_v = np.array(sigma_noise_v)/1000.
sigma_noise_r0 = np.array(sigma_noise_r0)/1000.
sigma_noise_v0 = np.array(sigma_noise_v0)/1000.

# Read SpOCK main input file to figure out the name of the output file
spock_main_input = sys.argv[1]
var_in, var_in_order = read_input_file(spock_main_input)
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 

# Read SpOCK output ECI r/v
isc = 0
eci_out = open(output_file_path_list[isc] + output_file_name_list[isc])
read_eci_out = eci_out.readlines()
nb_header = 0
while read_eci_out[nb_header].split()[0][0:2] != "20":
    nb_header = nb_header + 1
n = len(read_eci_out) - nb_header
true_r = np.zeros([n,3])
true_v = np.zeros([n,3])
date = []
nb_seconds_since_start = []
for itime in range(n):
    date_temp = read_eci_out[itime+nb_header].split()[0] + " " + read_eci_out[itime+nb_header].split()[1] 
    date.append( datetime.strptime( date_temp, "%Y/%m/%d %H:%M:%S" ) )
    nb_seconds_since_start.append( (date[-1] - date[0]).total_seconds() )
    true_r[itime,0] = np.float( read_eci_out[itime+nb_header].split()[2] )
    true_r[itime,1] = np.float( read_eci_out[itime+nb_header].split()[3] )
    true_r[itime,2] = np.float( read_eci_out[itime+nb_header].split()[4] )
    true_v[itime,0] = np.float( read_eci_out[itime+nb_header].split()[5] )
    true_v[itime,1] = np.float( read_eci_out[itime+nb_header].split()[6] )
    true_v[itime,2] = np.float( read_eci_out[itime+nb_header].split()[7] )

# Time axis (in hours if more than 3 hours, otherwise in minutes)
nb_seconds_since_start = np.array(nb_seconds_since_start)
if ( date[-1] - date[0] ).total_seconds() > 3*3600:
    x_unit = "hour"
    x_axis = nb_seconds_since_start/3600.
else:
    x_unit = "min"
    x_axis = nb_seconds_since_start/60.

# Add random gaussian noise of std dev sigma_noise
noise_r = np.zeros([n,3])
noise_v = np.zeros([n,3])
# noise_r[0, :] = true_r[0, :] #!!!!!!!!!!! to remove
# noise_v[0, :] = true_v[0, :] #!!!!!!!!!!! to remove
for itime in range(0, n): # !!!!!! should be range(n)
    if itime == 0:
        noise_r[itime, 0] = np.random.normal(true_r[itime, 0], sigma_noise_r0[0])
        noise_r[itime, 1] = np.random.normal(true_r[itime, 1], sigma_noise_r0[1])
        noise_r[itime, 2] = np.random.normal(true_r[itime, 2], sigma_noise_r0[2])
        noise_v[itime, 0] = np.random.normal(true_v[itime, 0], sigma_noise_v0[0])
        noise_v[itime, 1] = np.random.normal(true_v[itime, 1], sigma_noise_v0[1])
        noise_v[itime, 2] = np.random.normal(true_v[itime, 2], sigma_noise_v0[2])
    else:        
        noise_r[itime, 0] = np.random.normal(true_r[itime, 0], sigma_noise_r[0])
        noise_r[itime, 1] = np.random.normal(true_r[itime, 1], sigma_noise_r[1])
        noise_r[itime, 2] = np.random.normal(true_r[itime, 2], sigma_noise_r[2])
        noise_v[itime, 0] = np.random.normal(true_v[itime, 0], sigma_noise_v[0])
        noise_v[itime, 1] = np.random.normal(true_v[itime, 1], sigma_noise_v[1])
        noise_v[itime, 2] = np.random.normal(true_v[itime, 2], sigma_noise_v[2])

# Write results of noise in a file
noise_out = open(output_file_path_list[isc] + "noise_" + output_file_name_list[isc], "w+" )
print >> noise_out, "This file corresponds to the true ECI r/v modeled by SpOCK from the main input file " + spock_main_input + ", with a random gaussian noise of standard deviation " + '{0:.10f}'.format(sigma_noise_r[0]) + " km on x_true, "  + '{0:.10f}'.format(sigma_noise_r[1]) + " km on y_true, " + '{0:.10f}'.format(sigma_noise_r[2]) + " km on z_true, " + '{0:.10f}'.format(sigma_noise_v[0]) + " km/s on vx_true, "  + '{0:.10f}'.format(sigma_noise_v[1]) + " km/s on vy_true, " + '{0:.10f}'.format(sigma_noise_v[2]) + " km on vz_true."
print >> noise_out, "#START"
for itime in range(n):
    date_temp = datetime.strftime(date[itime], "%Y-%m-%dT%H:%M:%S.%f")
    print >> noise_out, date_temp, '{0:.10f}'.format(noise_r[itime, 0]), '{0:.10f}'.format(noise_r[itime, 1]), '{0:.10f}'.format(noise_r[itime, 2]), '{0:.10f}'.format(noise_v[itime, 0]), '{0:.10f}'.format(noise_v[itime, 1]), '{0:.10f}'.format(noise_v[itime, 2])

noise_out.close()

# Plot true and noise orbit
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3

if show_plot == 1:
    # Position 
    if x_unit == "hour":
        x_label = 'Time (hr)'
    else:
        x_label = 'Time (min)'
    fig_title = ''
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(3, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.1, wspace = 0.35)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    
    ## X
    y_label = 'X (m)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Position - Noise (b), True (r)', weight = 'bold', fontsize = fontsize_plot)
    ax1.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax1.xaxis.set_ticklabels("")
    ax1.scatter(x_axis, (noise_r[:,0] - true_r[:,0] ) * 1000., linewidth = 2,label='Noise', color = 'b')
    ax1.plot(x_axis, (true_r[:,0] - true_r[:,0] ) * 1000., linewidth = 4,label='True', color = 'r')
    #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax1.margins(0,0)

    ## Y
    y_label = 'Y (m)'
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax2.xaxis.set_ticklabels("")
    ax2.scatter(x_axis, (noise_r[:,1] - true_r[:,1] ) * 1000., linewidth = 2,label='Noise', color = 'b')
    ax2.plot(x_axis, (true_r[:,1] - true_r[:,1] ) * 1000., linewidth = 4,label='True', color = 'r')
    ax2.margins(0,0)

    ## X
    y_label = 'Z (m)'
    ax3 = fig.add_subplot(gs[2, 0])
    ax3.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax3.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
    ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    ax3.scatter(x_axis, (noise_r[:,2] - true_r[:,2] ) * 1000., linewidth = 2,label='Noise', color = 'b')
    ax3.plot(x_axis, (true_r[:,2] - true_r[:,2] ) * 1000., linewidth = 4,label='True', color = 'r')
    ax3.margins(0,0)    

    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_true_vs_noise_r.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


    # Velocity
    if x_unit == "hour":
        x_label = 'Time (hr)'
    else:
        x_label = 'Time (min)'
    fig_title = ''
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(3, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.1, wspace = 0.35)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    
    ## X
    y_label = 'Vx (m/s)'
    ax4 = fig.add_subplot(gs[0, 0])
    ax4.set_title('Velocity - Noise (b), True (r)', weight = 'bold', fontsize = fontsize_plot)
    ax4.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax4.spines.itervalues()] # change the width of the frame of the figure
    ax4.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax4.xaxis.set_ticklabels("")
    ax4.scatter(x_axis, (noise_v[:,0] - true_v[:,0] ) * 1000., linewidth = 2,label='Noise', color = 'b')
    ax4.plot(x_axis, (true_v[:,0] - true_v[:,0] ) * 1000., linewidth = 4,label='True', color = 'r')
    #legend = ax4.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax4.margins(0,0)

    ## Y
    y_label = 'Vy (m/s)'
    ax5 = fig.add_subplot(gs[1, 0])
    ax5.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax5.spines.itervalues()] # change the width of the frame of the figure
    ax5.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax5.xaxis.set_ticklabels("")
    ax5.scatter(x_axis, (noise_v[:,1] - true_v[:,1] ) * 1000., linewidth = 2,label='Noise', color = 'b')
    ax5.plot(x_axis, (true_v[:,1] - true_v[:,1] ) * 1000., linewidth = 4,label='True', color = 'r')
    ax5.margins(0,0)

    ## X
    y_label = 'Vz (m/s)'
    ax6 = fig.add_subplot(gs[2, 0])
    ax6.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax6.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax6.spines.itervalues()] # change the width of the frame of the figure
    ax6.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    ax6.scatter(x_axis, (noise_v[:,2] - true_v[:,2] ) * 1000., linewidth = 2,label='Noise', color = 'b')
    ax6.plot(x_axis, ( true_v[:,2] - true_v[:,2] ) * 1000., linewidth = 4,label='True', color = 'r')
    ax6.margins(0,0)

    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_true_vs_noise_v.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
