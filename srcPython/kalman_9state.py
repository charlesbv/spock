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
# This script is copy paste of kalman_6state.py adpated to 9 state. kalman.py is too old to be used (but has 9 elt in the state too).
from datetime import datetime, timedelta
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import sys
from read_input_file import *
from read_output_file import *


# How to run this script:
# python kalman.py arg
# arg can be:
# - 3d then plot 3d orbit
# - r then plot position on three graphs on same page
# - v then plot speed on three graphs on same page
# - ad
# - rho
# - Pdiag
# - dX
# - y
# - drag_sigma

# Assumptions:
# - the time in the measurement file must be the same as in the kalman file

# Read measurement file
## Read SpOCK main input file to figure out the name of the measurement file
input_filename = sys.argv[1]
var_in, var_in_order = read_input_file(input_filename)
filename_kalman_meas = var_in[find_in_read_input_order_variables(var_in_order, 'filename_kalman_meas')];

# "spock/true/true1/noise_true1.txt" #"/Users/cbv/kalman/cbv/spock/netcdf/output/cyg07.ddmi.s20170609-000000-e20170609-235959.l1.power-brcs.sand004.txt"
file_meas = open(filename_kalman_meas)
read_file_meas = file_meas.readlines()
nb_header = 0
while (read_file_meas[nb_header].split()[0] != '#START'):
    nb_header = nb_header + 1
nb_header = nb_header + 1
n_meas = len(read_file_meas) - nb_header
r_meas = np.zeros([n_meas,3])
v_meas = np.zeros([n_meas,3])
date_meas = []
for itime in range(n_meas):
    date_meas.append(read_file_meas[itime+nb_header].split()[0])
    r_meas[itime,0] = read_file_meas[itime+nb_header].split()[1]
    r_meas[itime,1] = read_file_meas[itime+nb_header].split()[2]
    r_meas[itime,2] = read_file_meas[itime+nb_header].split()[3]
    v_meas[itime,0] = read_file_meas[itime+nb_header].split()[4]
    v_meas[itime,1] = read_file_meas[itime+nb_header].split()[5]
    v_meas[itime,2] = read_file_meas[itime+nb_header].split()[6]

file_meas.close()

# Read Kalman file

#### Read SpOCK main input file to figure out stuff to then read the output
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
dt_kalm = var_in[find_in_read_input_order_variables(var_in_order, 'dt')]; 


## Read kalman filter output file
isc = 0 # !!!!!!! all sc
filename_kalm = output_file_path_list[isc] + "kalman_" + output_file_name_list[isc]
file_kalm = open(filename_kalm)
read_file_kalm = file_kalm.readlines()
nb_header = 0
# while (read_file_kalm[nb_header].split()[0] != '#START'):
#     nb_header = nb_header + 1
# nb_header = nb_header + 1
n_kalm = len(read_file_kalm) - nb_header
if n_kalm != n_meas:
    print "The nnumber of measurements is different from the number of Kalman filter output, that's weird."; raise Exception

r_kalm = np.zeros([n_kalm,3]) # ECI
v_kalm = np.zeros([n_kalm,3]) # ECI
a_kalm = np.zeros([n_kalm,3]) # LVLH
Pdiag = np.zeros([n_kalm,7])
dX = np.zeros([n_kalm,7])
y = np.zeros([n_kalm,3])
sy = np.zeros([n_kalm,3])
y_pf = np.zeros([n_kalm,3])
ad_kalm = np.zeros([n_kalm])
rho_kalm = np.zeros([n_kalm])
drag_sigma_kalm = np.zeros([n_kalm])

date_kalm = []
date = []
nb_seconds_since_start = []
for itime in range(n_kalm):
    date_kalm.append(read_file_kalm[itime+nb_header].split()[0])
    date.append(datetime.strptime( date_kalm[-1], "%Y-%m-%dT%H:%M:%S.%f" ))
    nb_seconds_since_start.append( (date[-1] - date[0]).total_seconds() )
    r_kalm[itime,0] = read_file_kalm[itime+nb_header].split()[1]
    r_kalm[itime,1] = read_file_kalm[itime+nb_header].split()[2]
    r_kalm[itime,2] = read_file_kalm[itime+nb_header].split()[3]
    v_kalm[itime,0] = read_file_kalm[itime+nb_header].split()[4]
    v_kalm[itime,1] = read_file_kalm[itime+nb_header].split()[5]
    v_kalm[itime,2] = read_file_kalm[itime+nb_header].split()[6]
    Pdiag[itime,0] = read_file_kalm[itime+nb_header].split()[7]
    Pdiag[itime,1] = read_file_kalm[itime+nb_header].split()[8]
    Pdiag[itime,2] = read_file_kalm[itime+nb_header].split()[9]
    Pdiag[itime,3] = read_file_kalm[itime+nb_header].split()[10]
    Pdiag[itime,4] = read_file_kalm[itime+nb_header].split()[11]
    Pdiag[itime,5] = read_file_kalm[itime+nb_header].split()[12]
    dX[itime,0] = read_file_kalm[itime+nb_header].split()[13]
    dX[itime,1] = read_file_kalm[itime+nb_header].split()[14]
    dX[itime,2] = read_file_kalm[itime+nb_header].split()[15]
    dX[itime,3] = read_file_kalm[itime+nb_header].split()[16]
    dX[itime,4] = read_file_kalm[itime+nb_header].split()[17]
    dX[itime,5] = read_file_kalm[itime+nb_header].split()[18]
    y[itime,0] = read_file_kalm[itime+nb_header].split()[19]
    y[itime,1] = read_file_kalm[itime+nb_header].split()[20]
    y[itime,2] = read_file_kalm[itime+nb_header].split()[21]
    sy[itime,0] = read_file_kalm[itime+nb_header].split()[22]
    sy[itime,1] = read_file_kalm[itime+nb_header].split()[23]
    sy[itime,2] = read_file_kalm[itime+nb_header].split()[24]
    y_pf[itime,0] = read_file_kalm[itime+nb_header].split()[25]
    y_pf[itime,1] = read_file_kalm[itime+nb_header].split()[26]
    y_pf[itime,2] = read_file_kalm[itime+nb_header].split()[27]
    a_kalm[itime,0] = read_file_kalm[itime+nb_header].split()[28]
    a_kalm[itime,1] = read_file_kalm[itime+nb_header].split()[29]
    a_kalm[itime,2] = read_file_kalm[itime+nb_header].split()[30]
    ad_kalm[itime] = read_file_kalm[itime+nb_header].split()[31]
    rho_kalm[itime] = read_file_kalm[itime+nb_header].split()[32] #* 10**9 # * 10**9 to convert from kg/m^3 to kg/km^3
    dX[itime,6] = read_file_kalm[itime+nb_header].split()[33]
    Pdiag[itime,6] = read_file_kalm[itime+nb_header].split()[34]
    drag_sigma_kalm[itime] = read_file_kalm[itime+nb_header].split()[35]

file_kalm.close()

# test if times are the same in measurement and kalman file
if date_kalm != date_meas:
    print "Times in the measurement file are different from times in the Kalman file. This is ok but not to run this script. The program will stop."; raise Exception

# Time axis (in hours if more than 3 hours, otherwise in minutes)
nb_seconds_since_start = np.array(nb_seconds_since_start)
if ( date[-1] - date[0] ).total_seconds() > 3*3600:
    x_unit = "hour"
    x_axis = nb_seconds_since_start/3600.
else:
    x_unit = "min"
    x_axis = nb_seconds_since_start/60.

# What to plot
plot_3d = 0
plot_r = 0
plot_v = 0
plot_Pdiag = 0
plot_dX = 0
plot_y = 0
plot_pred = 0
if '3d' in sys.argv:
    plot_3d = 1
if 'r' in sys.argv:
    plot_r = 1
if 'v' in sys.argv:
    plot_v = 1
if 'Pdiag' in sys.argv:
    plot_Pdiag = 1
if 'dX' in sys.argv:
    plot_dX = 1
if 'y' in sys.argv:
    plot_y = 1
compare_out = ""
for i in range(2,len(sys.argv)):
    if '.txt' in sys.argv[i]:
        compare_out = sys.argv[i]
if 'pred' in sys.argv:
    plot_pred = 1
    
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3

dt_plot = 60. # in seconds
if dt_plot <= dt_kalm:
    index_plot = 1
else:
    index_plot = (int)(dt_plot / dt_kalm)


# Plot 3d orbits
if plot_3d == 1:
    plt.ioff()
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig_title = '3D orbits'
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.1)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

    ax = fig.gca(projection='3d')
    ax.plot(r_kalm[::index_plot,0], r_kalm[::index_plot,1], r_kalm[::index_plot,2], label='Kalman', color = 'r', linewidth = 4)
    ax.scatter(r_meas[::index_plot,0], r_meas[::index_plot,1], r_meas[::index_plot,2], label='Measurement', color = 'b', linewidth = 6)
    x_label = "X (km)"
    y_label = "Y (km)"
    z_label = "Z (km)"
    ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot, labelpad = 35)
    ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot, labelpad = 35)
    ax.set_zlabel(z_label, weight = 'bold', fontsize  = fontsize_plot, labelpad = 35)
    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)

    legend = ax.legend(loc = 'center', bbox_to_anchor=(0.5, 1),  fontsize = fontsize_plot)
    fig.savefig(output_file_name_list[isc].replace(".txt", "") + "_3d.pdf")
    plt.show()

# Plot position
if plot_r == 1:
    ## Magnitude of difference (not difference of magnitude)
    mag_of_r_diff = np.zeros([n_kalm])
    for itime in range(n_kalm):
        mag_of_r_diff[itime] = np.linalg.norm( r_meas[itime,:] - r_kalm[itime,:] )
    fig_title = 'Magnitude of difference position meas - kalm'
    if x_unit == "hour":
        x_label = 'Time (hr)'
    else:
        x_label = 'Time (min)'
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.1, wspace = 0.35)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    
    ## X
    y_label = 'Magnitude of delta r (km)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax1.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax1.plot(x_axis[::index_plot], mag_of_r_diff[::index_plot], linewidth = 2,label='', color = 'k')
    ax1.margins(0,0)


    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_r_mag.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

    ## All components
    fig_title = ''
    if x_unit == "hour":
        x_label = 'Time (hr)'
    else:
        x_label = 'Time (min)'
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(3, 2)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.1, wspace = 0.35)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    
    ## X
    y_label = 'X (km)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Meas (b), Kalm (r)', weight = 'bold', fontsize = fontsize_plot)
    ax1.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax1.xaxis.set_ticklabels("")
    ax1.plot(x_axis[::index_plot], r_meas[::index_plot,0], linewidth = 2,label='Measurement', color = 'b')
    ax1.plot(x_axis[::index_plot], r_kalm[::index_plot,0], linewidth = 2,label='Kalman', color = 'r')
    #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax1.margins(0,0)

    ## dX
    y_label = 'dX (km)'
    ax4 = fig.add_subplot(gs[0, 1])
    ax4.set_title('Meas - Kalm', weight = 'bold', fontsize = fontsize_plot)
    ax4.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax4.spines.itervalues()] # change the width of the frame of the figure
    ax4.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax4.xaxis.set_ticklabels("")
    ax4.plot(x_axis[::index_plot], r_meas[::index_plot,0] - r_kalm[::index_plot,0], linewidth = 2,label='', color = 'k')
    legend = ax4.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax4.margins(0,0)
    
    ## Y
    y_label = 'Y (km)'
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax2.xaxis.set_ticklabels("")
    ax2.plot(x_axis[::index_plot], r_meas[::index_plot,1], linewidth = 2,label='Measurement', color = 'b')
    ax2.plot(x_axis[::index_plot], r_kalm[::index_plot,1], linewidth = 2,label='Kalman', color = 'r')
    ax2.margins(0,0)

    ## dY
    y_label = 'dY (km)'
    ax5 = fig.add_subplot(gs[1, 1])
    ax5.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax5.spines.itervalues()] # change the width of the frame of the figure
    ax5.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax5.xaxis.set_ticklabels("")
    ax5.plot(x_axis[::index_plot], r_meas[::index_plot,1] - r_kalm[::index_plot,1], linewidth = 2,label='', color = 'k')
    legend = ax5.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax5.margins(0,0)


    ## X
    y_label = 'Z (km)'
    ax3 = fig.add_subplot(gs[2, 0])
    ax3.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax3.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
    ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    ax3.plot(x_axis[::index_plot], r_meas[::index_plot,2], linewidth = 2,label='Measurement', color = 'b')
    ax3.plot(x_axis[::index_plot], r_kalm[::index_plot,2], linewidth = 2,label='Kalman', color = 'r')
    ax3.margins(0,0)    

    ## dZ
    y_label = 'dZ (km)'
    ax6 = fig.add_subplot(gs[2, 1])
    ax6.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax6.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax6.spines.itervalues()] # change the width of the frame of the figure
    ax6.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax6.plot(x_axis[::index_plot], r_meas[::index_plot,2] - r_kalm[::index_plot,2], linewidth = 2,label='', color = 'k')
    legend = ax6.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax6.margins(0,0)


    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_r.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

# Plot velocity
if plot_v == 1:
    ## Magnitude of difference (not difference of magnitude)
    mag_of_v_diff = np.zeros([n_kalm])
    for itime in range(n_kalm):
        mag_of_v_diff[itime] = np.linalg.norm( v_meas[itime,:] - v_kalm[itime,:] )
    fig_title = 'Magnitude of difference velocity meas - kalm'
    if x_unit == "hour":
        x_label = 'Time (hr)'
    else:
        x_label = 'Time (min)'
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.1, wspace = 0.35)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    
    ## X
    y_label = 'Magnitude of delta v (km/s)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax1.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax1.plot(x_axis[::index_plot], mag_of_v_diff[::index_plot], linewidth = 2,label='', color = 'k')
    ax1.margins(0,0)


    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_v_mag.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

    ## All components

    fig_title = ''
    if x_unit == "hour":
        x_label = 'Time (hr)'
    else:
        x_label = 'Time (min)'
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(3, 2)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.1, wspace = 0.35)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    
    ## X
    y_label = 'Vx (km/s)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Meas (b), Kalm (r)', weight = 'bold', fontsize = fontsize_plot)
    ax1.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax1.xaxis.set_ticklabels("")
    ax1.plot(x_axis[::index_plot], v_meas[::index_plot,0], linewidth = 2,label='Measurement', color = 'b')
    ax1.plot(x_axis[::index_plot], v_kalm[::index_plot,0], linewidth = 2,label='Kalman', color = 'r')
    #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax1.margins(0,0)

    ## dX
    y_label = 'dVx (km/s)'
    ax4 = fig.add_subplot(gs[0, 1])
    ax4.set_title('Meas - Kalm', weight = 'bold', fontsize = fontsize_plot)
    ax4.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax4.spines.itervalues()] # change the width of the frame of the figure
    ax4.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax4.xaxis.set_ticklabels("")
    ax4.plot(x_axis[::index_plot], v_meas[::index_plot,0] - v_kalm[::index_plot,0], linewidth = 2,label='', color = 'k')
    legend = ax4.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax4.margins(0,0)
    
    ## Y
    y_label = 'Vy (km/s)'
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax2.xaxis.set_ticklabels("")
    ax2.plot(x_axis[::index_plot], v_meas[::index_plot,1], linewidth = 2,label='Measurement', color = 'b')
    ax2.plot(x_axis[::index_plot], v_kalm[::index_plot,1], linewidth = 2,label='Kalman', color = 'r')
    ax2.margins(0,0)

    ## dY
    y_label = 'dVy (km/s)'
    ax5 = fig.add_subplot(gs[1, 1])
    ax5.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax5.spines.itervalues()] # change the width of the frame of the figure
    ax5.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax5.xaxis.set_ticklabels("")
    ax5.plot(x_axis[::index_plot], v_meas[::index_plot,1] - v_kalm[::index_plot,1], linewidth = 2,label='', color = 'k')
    legend = ax5.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax5.margins(0,0)


    ## X
    y_label = 'Vz (km/s)'
    ax3 = fig.add_subplot(gs[2, 0])
    ax3.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax3.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
    ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    ax3.plot(x_axis[::index_plot], v_meas[::index_plot,2], linewidth = 2,label='Measurement', color = 'b')
    ax3.plot(x_axis[::index_plot], v_kalm[::index_plot,2], linewidth = 2,label='Kalman', color = 'r')
    ax3.margins(0,0)    

    ## dZ
    y_label = 'dVz (km/s)'
    ax6 = fig.add_subplot(gs[2, 1])
    ax6.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax6.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax6.spines.itervalues()] # change the width of the frame of the figure
    ax6.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax6.plot(x_axis[::index_plot], v_meas[::index_plot,2] - v_kalm[::index_plot,2], linewidth = 2,label='', color = 'k')
    legend = ax6.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax6.margins(0,0)


    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_v.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
# Plot Pdiag
if plot_Pdiag == 1:
    fig_title = ''
    if x_unit == "hour":
        x_label = 'Time (hr)'
    else:
        x_label = 'Time (min)'

    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(3, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.27)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    
    ## Position uncertainty
    y_label = 'sigma position (km)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Position uncertainty', fontsize = fontsize_plot, weight = 'bold')
    ax1.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax1.xaxis.set_ticklabels("")
    ax1.plot(x_axis[::index_plot], Pdiag[::index_plot,0], linewidth = 2,label='X', color = 'b')
    ax1.plot(x_axis[::index_plot], Pdiag[::index_plot,1], linewidth = 2,label='Y', color = 'r')
    ax1.plot(x_axis[::index_plot], Pdiag[::index_plot,2], linewidth = 2,label='Z', color = 'k')
    legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax1.margins(0,0)

    ## Velocity uncertainty
    y_label = 'sigma velocity (km/s)'
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_title('Velocity uncertainty', fontsize = fontsize_plot, weight = 'bold')
    ax2.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax2.xaxis.set_ticklabels("")
    ax2.plot(x_axis[::index_plot], Pdiag[::index_plot,3], linewidth = 2,label='Vx', color = 'b')
    ax2.plot(x_axis[::index_plot], Pdiag[::index_plot,4], linewidth = 2,label='Vy', color = 'r')
    ax2.plot(x_axis[::index_plot], Pdiag[::index_plot,5], linewidth = 2,label='Vz', color = 'k')
    legend = ax2.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax2.margins(0,0)

    ## ? uncertainty
    y_label = '?'
    ax3 = fig.add_subplot(gs[2, 0])
    ax3.set_title('? uncertainty', fontsize = fontsize_plot, weight = 'bold')
    ax3.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax3.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
    ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax3.plot(x_axis[::index_plot], Pdiag[::index_plot,6], linewidth = 2,label='?', color = 'b')
    ax3.margins(0,0)    

    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_Pdiag.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# Plot dX
if plot_dX == 1:
    fig_title = ''
    if x_unit == "hour":
        x_label = 'Time (hr)'
    else:
        x_label = 'Time (min)'

    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(3, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.27)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    
    ## Position 
    y_label = 'dr (km)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Position', fontsize = fontsize_plot, weight = 'bold')
    ax1.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax1.xaxis.set_ticklabels("")
    ax1.plot(x_axis[::index_plot], dX[::index_plot,0], linewidth = 2,label='dX', color = 'b')
    ax1.plot(x_axis[::index_plot], dX[::index_plot,1], linewidth = 2,label='dY', color = 'r')
    ax1.plot(x_axis[::index_plot], dX[::index_plot,2], linewidth = 2,label='dZ', color = 'k')
    legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax1.margins(0,0)

    ## Velocity 
    y_label = 'dv (km/s)'
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_title('Velocity', fontsize = fontsize_plot, weight = 'bold')
    ax2.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax2.xaxis.set_ticklabels("")
    ax2.plot(x_axis[::index_plot], dX[::index_plot,3], linewidth = 2,label='dVx', color = 'b')
    ax2.plot(x_axis[::index_plot], dX[::index_plot,4], linewidth = 2,label='dVy', color = 'r')
    ax2.plot(x_axis[::index_plot], dX[::index_plot,5], linewidth = 2,label='dVz', color = 'k')
    legend = ax2.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax2.margins(0,0)

    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_dX.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# Plot y, sy, y_pf
if plot_y == 1:
    fig_title = ''
    if x_unit == "hour":
        x_label = 'Time (hr)'
    else:
        x_label = 'Time (min)'

    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(3, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.27)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    
    ## Position uncertainty
    y_label = 'y (km)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('y', fontsize = fontsize_plot, weight = 'bold')
    ax1.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax1.xaxis.set_ticklabels("")
    ax1.plot(x_axis[::index_plot], y[::index_plot,0], linewidth = 2,label='Y[0]', color = 'b')
    ax1.plot(x_axis[::index_plot], y[::index_plot,1], linewidth = 2,label='Y[1]', color = 'r')
    ax1.plot(x_axis[::index_plot], y[::index_plot,2], linewidth = 2,label='Y[2]', color = 'k')

    legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax1.margins(0,0)

    ## Velocity uncertainty
    y_label = 'sy (?)'
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_title('sy', fontsize = fontsize_plot, weight = 'bold')
    ax2.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax2.xaxis.set_ticklabels("")
    ax2.plot(x_axis[::index_plot], sy[::index_plot,0], linewidth = 2,label='sy[0]', color = 'b')
    ax2.plot(x_axis[::index_plot], sy[::index_plot,1], linewidth = 2,label='sy[1]', color = 'r')
    ax2.plot(x_axis[::index_plot], sy[::index_plot,2], linewidth = 2,label='sz[2]', color = 'k')
    legend = ax2.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax2.margins(0,0)

    ## Velocity uncertainty
    y_label = 'y_pf (km)'
    ax3 = fig.add_subplot(gs[2, 0])
    ax3.set_title('y_pf', fontsize = fontsize_plot, weight = 'bold')
    ax3.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax3.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
    ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax3.plot(x_axis[::index_plot], y_pf[::index_plot,0], linewidth = 2,label='y_pf[0]', color = 'b')
    ax3.plot(x_axis[::index_plot], y_pf[::index_plot,1], linewidth = 2,label='y_pf[1]', color = 'r')
    ax3.plot(x_axis[::index_plot], y_pf[::index_plot,2], linewidth = 2,label='y_pf[2]', color = 'k')
    legend = ax3.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax3.margins(0,0)
    
    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_y_sy_ypf.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  




if compare_out != "":
    var_in, var_in_order = read_input_file(compare_out)
    output_file_path_list_comp = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
    output_file_name_list_comp = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')];
    var_to_read = ["position","velocity", "acceleration_lvlh", "density"]
    var_out, var_out_order = read_output_file( output_file_path_list_comp[isc] + output_file_name_list_comp[isc], var_to_read )
    r_comp = var_out[find_in_read_input_order_variables(var_out_order, 'position')]
    v_comp = var_out[find_in_read_input_order_variables(var_out_order, 'velocity')]
    a_comp = var_out[find_in_read_input_order_variables(var_out_order, 'acceleration_lvlh')]
    rho_comp = var_out[find_in_read_input_order_variables(var_out_order, 'density')] * 10**9 # * 10**9 to convert from kg/m^3 to kg/km^3
    n_comp = len(a_comp)
    ad_comp = np.zeros([n_comp])
    for itime in range(n_comp):
        ad_comp[itime] = np.linalg.norm(a_comp[itime,:])

    if x_unit == "hour":
        x_label = 'Time (hr)'
    else:
        x_label = 'Time (min)'
    fig_title = ''
    # fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    # fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    # plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    # gs = gridspec.GridSpec(3, 1)
    # gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.1, wspace = 0.35)
    # plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    
    # ## X
    # y_label = 'X (m)'
    # ax1 = fig.add_subplot(gs[0, 0])
    # ax1.set_title('Position - Kalman (b), True (r)', weight = 'bold', fontsize = fontsize_plot)
    # ax1.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    # [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    # ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    # ax1.xaxis.set_ticklabels("")
    # ax1.scatter(x_axis[::index_plot], ( r_kalm[::index_plot,0] - r_comp[::index_plot,0] ) * 1000, linewidth = 2,label='Kalman', color = 'b')
    # ax1.plot(x_axis[::index_plot], ( r_comp[::index_plot,0] - r_comp[::index_plot,0] ) * 1000, linewidth = 4,label='True', color = 'r')
    # #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
    # ax1.margins(0,0)

    # ## Y
    # y_label = 'Y (m)'
    # ax2 = fig.add_subplot(gs[1, 0])
    # ax2.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    # [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    # ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    # ax2.xaxis.set_ticklabels("")
    # ax2.scatter(x_axis[::index_plot], ( r_kalm[::index_plot,1] - r_comp[::index_plot,1] ) * 1000, linewidth = 2,label='Kalman', color = 'b')
    # ax2.plot(x_axis[::index_plot], ( r_comp[::index_plot,1] - r_comp[::index_plot,1] ) * 1000, linewidth = 4,label='True', color = 'r')
    # ax2.margins(0,0)

    # ## X
    # y_label = 'Z (m)'
    # ax3 = fig.add_subplot(gs[2, 0])
    # ax3.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    # ax3.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
    # [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
    # ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    # ax3.scatter(x_axis[::index_plot], ( r_kalm[::index_plot,2] - r_comp[::index_plot,2] ) * 1000, linewidth = 2,label='Kalman', color = 'b')
    # ax3.plot(x_axis[::index_plot], ( r_comp[::index_plot,2] - r_comp[::index_plot,2] ) * 1000, linewidth = 4,label='True', color = 'r')
    # ax3.margins(0,0)    

    # fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_true_vs_kalm_r.pdf'
    # fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

    
    # if x_unit == "hour":
    #     x_label = 'Time (hr)'
    # else:
    #     x_label = 'Time (min)'
    # fig_title = ''
    # fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    # fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    # plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    # gs = gridspec.GridSpec(3, 1)
    # gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.1, wspace = 0.35)
    # plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    
    # ## X
    # y_label = 'Vx (m/s)'
    # ax1 = fig.add_subplot(gs[0, 0])
    # ax1.set_title('Velocity - Kalman (b), True (r)', weight = 'bold', fontsize = fontsize_plot)
    # ax1.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    # [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    # ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    # ax1.xaxis.set_ticklabels("")
    # ax1.scatter(x_axis[::index_plot], (v_kalm[::index_plot,0] - v_comp[::index_plot,0])*1000., linewidth = 2,label='Kalman', color = 'b')
    # ax1.plot(x_axis[::index_plot], (v_comp[::index_plot,0] - v_comp[::index_plot,0])*1000., linewidth = 4,label='True', color = 'r')
    # #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
    # ax1.margins(0,0)

    # ## Y
    # y_label = 'Vy (m/s)'
    # ax2 = fig.add_subplot(gs[1, 0])
    # ax2.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    # [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    # ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    # ax2.xaxis.set_ticklabels("")
    # ax2.scatter(x_axis[::index_plot], (v_kalm[::index_plot,1] - v_comp[::index_plot,1])*1000., linewidth = 2,label='Kalman', color = 'b')
    # ax2.plot(x_axis[::index_plot], (v_comp[::index_plot,1] - v_comp[::index_plot,1])*1000., linewidth = 4,label='True', color = 'r')
    # ax2.margins(0,0)

    # ## X
    # y_label = 'Vz (m/s)'
    # ax3 = fig.add_subplot(gs[2, 0])
    # ax3.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    # ax3.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
    # [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
    # ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    # ax3.scatter(x_axis[::index_plot], (v_kalm[::index_plot,2] - v_comp[::index_plot,2])*1000., linewidth = 2,label='Kalman', color = 'b')
    # ax3.plot(x_axis[::index_plot], (v_comp[::index_plot,2] - v_comp[::index_plot,2])*1000., linewidth = 4,label='True', color = 'r')
    # ax3.margins(0,0)    

    # fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_true_vs_kalm_v.pdf'
    # fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


    r_error_meas =  np.zeros([n_kalm, 3])
    r_error_meas_mag = np.zeros([n_kalm])
    v_error_meas =  np.zeros([n_kalm, 3])
    v_error_meas_mag = np.zeros([n_kalm])    
    r_error = np.zeros([n_kalm, 3])
    v_error = np.zeros([n_kalm, 3])
    a_error = np.zeros([n_kalm, 3])
    r_error_mag = np.zeros([n_kalm])
    v_error_mag = np.zeros([n_kalm])
    a_error_mag = np.zeros([n_kalm])
    a_comp_mag = np.zeros([n_kalm])
    for itime in range(n_kalm):
        r_error_meas[itime, :] = r_meas[itime,:] - r_comp[itime,:]
        v_error_meas[itime, :] = v_meas[itime,:] - v_comp[itime,:]
        r_error[itime, :] = r_kalm[itime,:] - r_comp[itime,:]
        v_error[itime, :] = v_kalm[itime,:] - v_comp[itime,:]
        a_error[itime, :] = a_kalm[itime,:] - a_comp[itime,:]
        r_error_meas_mag[itime] = np.linalg.norm(r_error_meas[itime, :])
        v_error_meas_mag[itime] = np.linalg.norm(v_error_meas[itime, :])
        r_error_mag[itime] = np.linalg.norm(r_error[itime, :])
        v_error_mag[itime] = np.linalg.norm(v_error[itime, :])
        a_error_mag[itime] = np.linalg.norm(a_error[itime, :])
        a_comp_mag[itime] = np.linalg.norm(a_comp[itime, :])
    r_rmse_meas = np.sqrt( np.mean(r_error_meas_mag**2) )
    v_rmse_meas = np.sqrt( np.mean(v_error_meas_mag**2) )
    r_rmse = np.sqrt( np.mean(r_error_mag**2) )
    v_rmse = np.sqrt( np.mean(v_error_mag**2) )
    a_rmse = np.sqrt( np.mean(a_error_mag**2) )
    a_nrmse = np.sqrt( np.mean(a_error_mag**2) / np.mean(a_comp_mag**2))
    ax_rmse = np.sqrt( np.mean(a_error[::index_plot,0]**2) )
    ay_rmse = np.sqrt( np.mean(a_error[::index_plot,1]**2) )
    az_rmse = np.sqrt( np.mean(a_error[::index_plot,2]**2) )
    ax_nrmse = np.sqrt( np.mean(a_error[::index_plot,0]**2) / np.mean(a_comp[::index_plot,0]**2))
    ay_nrmse = np.sqrt( np.mean(a_error[::index_plot,1]**2) / np.mean(a_comp[::index_plot,1]**2))
    az_nrmse = np.sqrt( np.mean(a_error[::index_plot,2]**2) / np.mean(a_comp[::index_plot,2]**2))

    # Plot rho VS time
    fig_title = ''
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(4, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.3, wspace = 0.35)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    
    ## X
    y_label = 'Mass density (kg/km$^3$)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Mass density as a function of time', weight = 'bold', fontsize = fontsize_plot)
    ax1.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax1.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    #ax1.xaxis.set_ticklabels("")
    ax1.scatter(x_axis[::index_plot], rho_kalm[::index_plot], linewidth = 2,label='Kalman', color = 'b')
    ax1.scatter(x_axis[::index_plot], rho_comp[::index_plot], linewidth = 2,label='True', color = 'r')
    #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax1.set_ylim([min([min(rho_kalm), min(rho_comp)]), max([max(rho_kalm), max(rho_comp)])])
    ax1.margins(0,0)
    print rho_comp
    print rho_kalm
    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_rho.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


    # Plot ad VS time
    fig_title = ''
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(4, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.3, wspace = 0.35)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    
    ## X
    y_label = '$|\mathrm{a_{LVLH}}|$ (m/s$^2$)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Magnitude of LVLH acceleration as a function of time', weight = 'bold', fontsize = fontsize_plot)
    ax1.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax1.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    #ax1.xaxis.set_ticklabels("")
    ax1.scatter(x_axis[::index_plot], ad_kalm[::index_plot]*1000, linewidth = 2,label='Kalman', color = 'b')
    ax1.scatter(x_axis[::index_plot], a_comp[::index_plot,0]*1000, linewidth = 2,label='True', color = 'r')
    #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax1.set_ylim([min([min(ad_kalm), min(ad_comp)]), max([max(ad_kalm), max(ad_comp)])])
    ax1.margins(0,0)

    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_ad.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


    #raise Exception

    # Plot error VS time
    fig_title = ''
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(4, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.3, wspace = 0.35)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    
    ## X
    y_label = 'Error (m)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Magnitude of position difference r_kalm - r_true', weight = 'bold', fontsize = fontsize_plot)
    ax1.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax1.xaxis.set_ticklabels("")
    ax1.scatter(x_axis[::index_plot], (r_error_mag[::index_plot])*1000., linewidth = 2,label='Kalman', color = 'b')
    #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax1.margins(0,0)
    ax1.text( max(x_axis[::index_plot]), max((r_error_mag)*1000.) - ( max((r_error_mag)*1000.) - min((r_error_mag)*1000.) )/15, "RMSE = " + '{:.2f}'.format(r_rmse*1000) + " m", fontsize = fontsize_plot, weight = 'bold', color = 'r', horizontalalignment = 'right')
    ax1.text( max(x_axis[::index_plot]), max((r_error_mag)*1000.) - 3* ( max((r_error_mag)*1000.) - min((r_error_mag)*1000.) )/15, "RMSE meas = " + '{:.2f}'.format(r_rmse_meas*1000) + " m", fontsize = fontsize_plot, weight = 'bold', color = 'r', horizontalalignment = 'right')

    ## Y
    y_label = 'Error (cm/s)'
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_title('Magnitude of velocity difference v_kalm - v_true', weight = 'bold', fontsize = fontsize_plot)
    ax2.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax2.xaxis.set_ticklabels("")
    ax2.scatter(x_axis[::index_plot], (v_error_mag[::index_plot])*100000., linewidth = 2,label='Kalman', color = 'b')
    ax2.margins(0,0)
    ax2.text( max(x_axis[::index_plot]),  max((v_error_mag)*100000.) - ( max((v_error_mag)*100000.) - min((v_error_mag)*100000.) )/15, "RMSE = " + '{:.2f}'.format(v_rmse*100000) + " cm/s", fontsize = fontsize_plot, weight = 'bold', color = 'r', horizontalalignment = 'right', verticalalignment = 'top')
    ax2.text( max(x_axis[::index_plot]),  max((v_error_mag)*100000.) - 3 * ( max((v_error_mag)*100000.) - min((v_error_mag)*100000.) )/15, "RMSE meas = " + '{:.2f}'.format(v_rmse_meas*100000) + " cm/s", fontsize = fontsize_plot, weight = 'bold', color = 'r', horizontalalignment = 'right', verticalalignment = 'top')


        ## Z
    y_label = 'Error (mm/s$^2$)'
    ax3 = fig.add_subplot(gs[2, 0])
    ax3.set_title('Magnitude of LVLH acceleration difference a_kalm - a_true', weight = 'bold', fontsize = fontsize_plot)
    ax3.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax3.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
    ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    # ax3.xaxis.set_ticklabels("")
    ax3.scatter(x_axis[::index_plot], (a_error_mag[::index_plot])*1000000., linewidth = 2,label='Kalman', color = 'b')
    ax3.margins(0,0)
    #print a_error_mag
    ax3.text( max(x_axis[::index_plot]),  max((a_error_mag)*1000000.) - ( max((a_error_mag)*1000000.) - min((a_error_mag)*1000000.) )/15, "RMSE = " + '{:.2f}'.format(a_rmse*1000000) + " mm/s$^2$", fontsize = fontsize_plot, weight = 'bold', color = 'r', horizontalalignment = 'right', verticalalignment = 'top')
    ax3.text( max(x_axis[::index_plot]),  max((a_error_mag)*1000000.) - 3 * ( max((a_error_mag)*1000000.) - min((a_error_mag)*1000000.) )/15, "NRMSE = " + '{:.3f}'.format(a_nrmse), fontsize = fontsize_plot, weight = 'bold', color = 'r', horizontalalignment = 'right', verticalalignment = 'top')
    ax3.text( max(x_axis[::index_plot]),  max((a_error_mag)*1000000.) - 5 * ( max((a_error_mag)*1000000.) - min((a_error_mag)*1000000.) )/15, "mean( norm( a_true ) ) = " + '{:.3f}'.format(np.mean(a_comp_mag)*1000000.) + " mm/s$^2$", fontsize = fontsize_plot, weight = 'bold', color = 'r', horizontalalignment = 'right', verticalalignment = 'top')

    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_error_mag.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    
    # Plot error acceleration components VS time
    fig_title = ''
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(4, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.8, wspace = 0.35)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    
    ## X
    y_label = 'Error (mm/s$^2$)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('$\mathrm{a_{x,kalm} - a_{x,true}}$', weight = 'bold', fontsize = fontsize_plot)
#    ax1.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
#    ax1.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax1.xaxis.set_ticklabels("")
    ax1.scatter(x_axis[::index_plot], (a_error[::index_plot,0])*1000000., linewidth = 2,label='Kalman', color = 'b')
    ax1.set_ylim([min((a_error[::index_plot,0])*1000000.), max((a_error[::index_plot,0])*1000000.)])
    ax1.margins(0,0)
    ax1.annotate("RMSE = " + '{:.0e}'.format(ax_rmse*1000000) + " mm/s$^2$", (0,-0.35), (0,0), xycoords='axes fraction', textcoords='offset points', va='bottom', fontsize = fontsize_plot, weight = 'bold')
    ax1.annotate("NRMSE = " + '{:.1e}'.format(ax_nrmse), (0.5,-0.35), (0, 0), xycoords='axes fraction', textcoords='offset points', va='bottom', fontsize = fontsize_plot, weight = 'bold', ha='center')
    ax1.annotate("$\overline{|\mathrm{a_{x,true}}|}$  = " + '{:.1f}'.format(np.mean(np.abs(a_comp[::index_plot,0]))*1000000.) + " mm/s$^2$", (1,-0.35), (0, 0), xycoords='axes fraction', textcoords='offset points', va='bottom', fontsize = fontsize_plot, weight = 'bold', ha = 'right')

        ## X
    y_label = 'Error (mm/s$^2$)'
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_title('$\mathrm{a_{y,kalm} - a_{y,true}}$', weight = 'bold', fontsize = fontsize_plot)
    ax2.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
#    ax2.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax2.xaxis.set_ticklabels("")
    ax2.scatter(x_axis[::index_plot], (a_error[::index_plot,1])*1000000., linewidth = 2,label='Kalman', color = 'b')
    ax2.set_ylim([min((a_error[::index_plot,1])*1000000.), max((a_error[::index_plot,1])*1000000.)])
    ax2.margins(0,0)
    ax2.annotate("RMSE = " + '{:.0e}'.format(ay_rmse*1000000) + " mm/s$^2$", (0,-0.35), (0,0), xycoords='axes fraction', textcoords='offset points', va='bottom', fontsize = fontsize_plot, weight = 'bold')
    ax2.annotate("NRMSE = " + '{:.1e}'.format(ay_nrmse), (0.5,-0.35), (0, 0), xycoords='axes fraction', textcoords='offset points', va='bottom', fontsize = fontsize_plot, weight = 'bold', ha='center')
    ax2.annotate("$\overline{|\mathrm{a_{y,true}}|}$  = " + '{:.1f}'.format(np.mean(np.abs(a_comp[::index_plot,1]))*1000000.) + " mm/s$^2$", (1,-0.35), (0, 0), xycoords='axes fraction', textcoords='offset points', va='bottom', fontsize = fontsize_plot, weight = 'bold', ha = 'right')


    # ax2.text( min(x_axis[::index_plot]),  max((a_error[::index_plot,1])*1000000.) + 3. * ( max((a_error[::index_plot,1])*1000000.) - min((a_error[::index_plot,1])*1000000.) )/15, "RMSE = " + '{:.2f}'.format(ax_rmse*1000000) + " mm/s$^2$", fontsize = fontsize_plot, weight = 'bold', color = 'r', horizontalalignment = 'left', verticalalignment = 'top')
    # ax2.text( (max(x_axis[::index_plot])+min(x_axis[::index_plot]))/2.,  max((a_error[::index_plot,1])*1000000.) + 3 * ( max((a_error[::index_plot,1])*1000000.) - min((a_error[::index_plot,1])*1000000.) )/15, "NRMSE = " + '{:.3f}'.format(ax_nrmse), fontsize = fontsize_plot, weight = 'bold', color = 'r', horizontalalignment = 'center', verticalalignment = 'top')
    # ax2.text( max(x_axis[::index_plot]),  max((a_error[::index_plot,1])*1000000.) + 3 * ( max((a_error[::index_plot,1])*1000000.) - min((a_error[::index_plot,1])*1000000.) )/15, "ax_true  ~ " + '{:.3f}'.format(np.mean(np.abs(a_comp[::index_plot,1]))*1000000.) + " mm/s$^2$", fontsize = fontsize_plot, weight = 'bold', color = 'r', horizontalalignment = 'right', verticalalignment = 'top')

        ## X
    y_label = 'Error (mm/s$^2$)'
    ax3 = fig.add_subplot(gs[2, 0])
    ax3.set_title('$\mathrm{a_{z,kalm} - a_{z,true}}$', weight = 'bold', fontsize = fontsize_plot)
#    ax3.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax3.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot, labelpad = 30)
    [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
    ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    #ax3.xaxis.set_ticklabels("")
    ax3.scatter(x_axis[::index_plot], (a_error[::index_plot,2])*1000000., linewidth = 2,label='Kalman', color = 'b')
    ax3.set_ylim([min((a_error[::index_plot,2])*1000000.), max((a_error[::index_plot,2])*1000000.)])
    ax3.margins(0,0)
    # ax3.tick_params(direction='out', pad=15)
    ax3.annotate("RMSE = " + '{:.0e}'.format(az_rmse*1000000) + " mm/s$^2$", (0,-0.60), (0,0), xycoords='axes fraction', textcoords='offset points', va='bottom', fontsize = fontsize_plot, weight = 'bold')
    ax3.annotate("NRMSE = " + '{:.1e}'.format(az_nrmse), (0.5,-0.60), (0, 0), xycoords='axes fraction', textcoords='offset points', va='bottom', fontsize = fontsize_plot, weight = 'bold', ha='center')
    ax3.annotate("$\overline{|\mathrm{a_{z,true}}|}$  = " + '{:.1f}'.format(np.mean(np.abs(a_comp[::index_plot,2]))*1000000.) + " mm/s$^2$", (1,-0.60), (0, 0), xycoords='axes fraction', textcoords='offset points', va='bottom', fontsize = fontsize_plot, weight = 'bold', ha = 'right')

    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_error_acc_components.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


#        LVLV ACCELERATION COMPONENTS OF KALMAN
    # Plot error VS time
    fig_title = ''
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(3, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.3, wspace = 0.35)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    
    ## X
    y_label = 'Acc (m/s$^2$)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('LVLH X component of true acceleration', weight = 'bold', fontsize = fontsize_plot)
    ax1.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax1.xaxis.set_ticklabels("")
    ax1.scatter(x_axis[::index_plot][1:], (a_kalm[1::index_plot,0])*1000., linewidth = 2,label='Kalman', color = 'b')
    #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax1.margins(0,0)

    ## Y
    y_label = 'Acc (m/s$^2$)'
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_title('LVLH Y component of true acceleration', weight = 'bold', fontsize = fontsize_plot)
    ax2.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax2.xaxis.set_ticklabels("")
    ax2.scatter(x_axis[::index_plot][1:], (a_kalm[1::index_plot,1])*1000., linewidth = 2,label='Kalman', color = 'b')
    ax2.margins(0,0)


        ## Z
    y_label = 'Acc (m/s$^2$)'
    ax3 = fig.add_subplot(gs[2, 0])
    ax3.set_title('LVLH Z component of true acceleration', weight = 'bold', fontsize = fontsize_plot)
    ax3.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax3.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
    ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    # ax3.xaxis.set_ticklabels("")
    ax3.scatter(x_axis[::index_plot][1:], (a_kalm[1::index_plot,2])*1000., linewidth = 2,label='Kalman', color = 'b')
    ax3.margins(0,0)

    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_acceleration.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

#    LVLV ACCELERATION COMPONENTS OF COMP
    # Plot error VS time
    fig_title = ''
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(3, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.3, wspace = 0.35)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    
    ## X
    y_label = 'Acc (m/s$^2$)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('LVLH X component of true acceleration', weight = 'bold', fontsize = fontsize_plot)
    ax1.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax1.xaxis.set_ticklabels("")
    ax1.scatter(x_axis[::index_plot][1:], (a_comp[1::index_plot,0])*1000., linewidth = 2,label='Kalman', color = 'b')
    #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax1.margins(0,0)

    ## Y
    y_label = 'Acc (m/s$^2$)'
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_title('LVLH Y component of true acceleration', weight = 'bold', fontsize = fontsize_plot)
    ax2.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax2.xaxis.set_ticklabels("")
    ax2.scatter(x_axis[::index_plot][1:], (a_comp[1::index_plot,1])*1000., linewidth = 2,label='Kalman', color = 'b')
    ax2.margins(0,0)


        ## Z
    y_label = 'Acc (m/s$^2$)'
    ax3 = fig.add_subplot(gs[2, 0])
    ax3.set_title('LVLH Z component of true acceleration', weight = 'bold', fontsize = fontsize_plot)
    ax3.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax3.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
    ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    # ax3.xaxis.set_ticklabels("")
    ax3.scatter(x_axis[::index_plot][1:], (a_comp[1::index_plot,  2])*1000., linewidth = 2,label='Kalman', color = 'b')
    ax3.margins(0,0)

    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_acceleration_' + output_file_name_list_comp[isc].replace(".txt", ".pdf")
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


    itime = 0
    mag_delta_r_temp = r_kalm[itime,:] - r_comp[itime,:]
    mag_delta_v_temp = v_kalm[itime,:] - v_comp[itime,:]
    print "Magnitude of delta r at time " + str(itime) + ": ", np.linalg.norm(mag_delta_r_temp) * 1000, "m"
    print "Magnitude of delta v at time " + str(itime) + ": ", np.linalg.norm(mag_delta_v_temp)*1000., "m/s"



    if plot_pred == 1: # this is the plot same as Kalman but using only the propagation, not the Kalman filter (equivalent to K = 0)
        var_to_read = ["position","velocity", "acceleration_lvlh", "density"]
        var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
        r_pred = var_out[find_in_read_input_order_variables(var_out_order, 'position')]
        v_pred = var_out[find_in_read_input_order_variables(var_out_order, 'velocity')]
        a_pred = var_out[find_in_read_input_order_variables(var_out_order, 'acceleration_lvlh')]
        rho_pred = var_out[find_in_read_input_order_variables(var_out_order, 'density')] * 10**9 # * 10**9 to convert from kg/m^3 to kg/km^3

        # Plot rho VS time
        fig_title = ''
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

        fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(4, 1)
        gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.3, wspace = 0.35)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold

        ## X
        y_label = 'Mass density (kg/km$^3$)'
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Mass density as a function of time', weight = 'bold', fontsize = fontsize_plot)
        ax1.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
        ax1.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
        #ax1.xaxis.set_ticklabels("")
        ax1.scatter(x_axis[::index_plot], rho_pred[::index_plot], linewidth = 2,label='Prediction', color = 'b')
        ax1.scatter(x_axis[::index_plot], rho_comp[::index_plot], linewidth = 2,label='True', color = 'r')
        #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
        ax1.set_ylim([min([min(rho_pred), min(rho_comp)]), max([max(rho_pred), max(rho_comp)])])
        ax1.margins(0,0)

        fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_rho_pred.pdf'
        fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


        # Plot ad VS time
        fig_title = ''
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

        fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(4, 1)
        gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.3, wspace = 0.35)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold

        ## X
        y_label = '$|\mathrm{a_{LVLH}}|$ (m/s$^2$)'
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Magnitude of LVLH acceleration as a function of time', weight = 'bold', fontsize = fontsize_plot)
        ax1.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
        ax1.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
        #ax1.xaxis.set_ticklabels("")
        ax1.scatter(x_axis[::index_plot], a_pred[::index_plot,0]*1000, linewidth = 2,label='Predictio', color = 'b')
        ax1.scatter(x_axis[::index_plot], a_comp[::index_plot,0]*1000, linewidth = 2,label='True', color = 'r')
        #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
        ax1.set_ylim([min([min(a_pred[::index_plot,0]*1000), min(a_comp[::index_plot,0]*1000)]), max([max(a_pred[::index_plot,0]*1000), max(a_comp[::index_plot,0]*1000)])])
        ax1.margins(0,0)

        fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_ad_pred.pdf'
        fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


