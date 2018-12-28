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
# This script runs different scenarios with SpOCK corresponding to different values of density_mod: coefficient to apply on the density modeled by NRLSMIS. The goal is to find the optimun density_mod (called rho_mod in this script) to minimize the distance between the baseline ephemerides and the ephemerides modeled by SpOCK with the given rho_mod
# ASSUMPTIONS
#- see section "PARAMETERS TO SET UP BEFORE RUNNIG THIS SCRIPT"
#- run SpOCK with a 1s time step 

import numpy as np
import sys
sys.path.append("../../kalman/spock_development_new_structure_kalman_dev/srcPython")
import os
from read_input_file import *
from read_output_file import *
from spock_main_input import *
from orbit_average import *
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec


#PARAMETERS TO SET UP BEFORE RUNNIG THIS SCRIPT
cygfm = 5 # which CYGFM to look at
path_mpirun = '/usr/local/bin/mpirun' #'/usr/local/bin/mpirun'# '/opt/local/bin/mpirun-openmpi-gcc49'
interval = 1 # interval of time to compare the two trajectories (data and SpOCK). In hours
rho_mod_min = 0.9 # min rho_mod
rho_mod_max = 0.9 # max rho_mod
drho_mod = 0.1 # step in rho_mod -> rho_mod varies from rho_mod_min to rho_mod_max by values of drho_mod
# end of PARAMETERS TO SET UP BEFORE RUNNIG THIS SCRIPT

# ALGORITHM
cygfm_to_ccsds = ['F7','F9','2B','2C','2F','36','37','49']
interval_sec = interval * 3600.



# Read data ephemerides
cygccsds = cygfm_to_ccsds[cygfm - 1]
filename = cygccsds + "_20170823_134632_STKdefpred_v001.e"
filename = 'cyg_data/' + filename 
file = open(filename)
read_file = file.readlines()
iline = 0
found_epoch = 0
while (found_epoch == 0):
    if ('YYDDD' in read_file[iline]):
        epochyyddd = np.float(read_file[iline].split(":")[1])
        found_epoch = 1
    iline = iline + 1

found_hdr = 0
while (found_hdr == 0):
    if (len(read_file[iline].split()) > 0):
        if read_file[iline].split()[0] == "EphemerisTimePosVel":
            found_hdr = 1
    iline = iline + 1
nb_header = iline + 1

n = len(read_file) - nb_header
date_data = []
r_data = []
v_data = []

epochyyddd_no_decimal = (int)(epochyyddd)
epochyyddd_decimal_only = epochyyddd - epochyyddd_no_decimal
epochyyddd_date = datetime.strptime( str(epochyyddd_no_decimal), "%y%j" ) + timedelta( hours = epochyyddd_decimal_only*24 )
date_data_old = []


iline = 0
new_date_start = []
nb_seconds_between_epochyyddd_and_new_date_start = []
new_r_data = []
new_v_data = []
index_interval = []
date_stop = "2017-08-21T13:08:00"
date_stop = datetime.strptime(date_stop, "%Y-%m-%dT%H:%M:%S")
while (iline < n):
    if iline > 0:
        if datetime.strptime(date_data[-1], "%Y/%m/%d %H:%M:%S") > date_stop:
            break
    if len(read_file[iline+nb_header].split()) == 0:
        break
    date_data_temp = np.float( read_file[iline+nb_header].split()[0] )
    if iline == 0:
        date_data.append( datetime.strftime( epochyyddd_date + timedelta( seconds = date_data_temp ) , "%Y/%m/%d %H:%M:%S") )
        r_data.append( [np.float(read_file[iline+nb_header].split()[1]), np.float(read_file[iline+nb_header].split()[2]), np.float(read_file[iline+nb_header].split()[3])] )
        v_data.append( [np.float(read_file[iline+nb_header].split()[4]), np.float(read_file[iline+nb_header].split()[5]), np.float(read_file[iline+nb_header].split()[6])] )
        index_interval.append(iline)
    else:
        index_interval.append(iline-1)
    new_date_start.append( date_data[-1] )
    new_r_data.append( r_data[-1] )
    new_v_data.append( v_data[-1] )
    nb_seconds_between_epochyyddd_and_first_new_date_start = date_data_temp
    time_elapsed = 0
    while time_elapsed < interval_sec:
        if len(read_file[iline+nb_header].split()) == 0:
            break
        date_data_temp = np.float( read_file[iline+nb_header].split()[0] )
        if iline > 0:
            date_data.append( datetime.strftime( epochyyddd_date + timedelta( seconds = date_data_temp ) , "%Y/%m/%d %H:%M:%S") )
            r_data.append( [np.float(read_file[iline+nb_header].split()[1]), np.float(read_file[iline+nb_header].split()[2]), np.float(read_file[iline+nb_header].split()[3])] )
            v_data.append( [np.float(read_file[iline+nb_header].split()[4]), np.float(read_file[iline+nb_header].split()[5]), np.float(read_file[iline+nb_header].split()[6])] )

        time_elapsed = date_data_temp - nb_seconds_between_epochyyddd_and_first_new_date_start
        iline = iline + 1

new_r_data = np.array(new_r_data)/1000.
new_v_data = np.array(new_v_data)/1000.

r_data = np.array(r_data)/1000.
v_data = np.array(v_data)/1000.

nb_interval = len(new_date_start)
rho_mod_arr = np.arange(rho_mod_min, rho_mod_max+drho_mod, drho_mod)
nb_rho = len(rho_mod_arr)
distance = []
nb_seconds_since_start = []
date_start_simu = "2017-08-21T11:55:00"# run SpOCK only starting at this date
date_start_simu = datetime.strptime(date_start_simu, "%Y-%m-%dT%H:%M:%S")
itime_count = -1
for itime in range(nb_interval-1):
    if datetime.strptime(new_date_start[itime],  "%Y/%m/%d %H:%M:%S") >= date_start_simu:
        itime_count = itime_count + 1
        if itime_count == 0:
            itime_start = itime
        print itime, nb_interval-1
        distance_interval = []
        # for each new interval, run SpOCK with different values of rho_mod: the initial date is the start date of the interval, initial r/v is the r/v at this date
        ## Create SpOCK main input file: same epoch and initial r/v
        date_start = new_date_start[itime].replace("/", "-").replace(" ", "T")+'.000'
        date_end = new_date_start[itime+1].replace("/", "-").replace(" ", "T")+'.000'
        nb_seconds_between_epochyyddd_and_new_date_start.append( ( datetime.strptime( new_date_start[itime], "%Y/%m/%d %H:%M:%S") - epochyyddd_date ).total_seconds() ) 
        dt  = 1
        dt_output = 1
        gravity_order = 50

        for irho in range(nb_rho):
            rho_mod = rho_mod_arr[irho]
            main_input_filename = 'FM0' + str(cygfm) + '_' + date_start.replace(":","_") + '_' + date_end.replace(":","_") + '_rhomod_' + str(rho_mod).replace(".", "") +'.txt'
            fac_velo = 1 - (1.000001-1)
            r0 = format(new_r_data[itime, 0], '.14e')
            r1 = format(new_r_data[itime, 1], '.14e')
            r2 = format(new_r_data[itime, 2], '.14e')
            v0 = format(new_v_data[itime, 0]*fac_velo, '.14e')
            v1 = format(new_v_data[itime, 1]*fac_velo, '.14e')
            v2 = format(new_v_data[itime, 2]*fac_velo, '.14e')

            spock_main_input( # need to be in spokc/srcPython to run this script   
                main_input_filename,
                # for TIME section
                    date_start,
                date_end,
                dt,
                # for SPACECRAFT section
                        1,
                '0',
                28,
                "cygnss_geometry_2016_acco08.txt", #"cygnss_geometry_2016_smaller_solar_radiation_coeff.txt", #"cygnss_geometry_2016.txt",#"cygnss_geometry_2016_acco09.txt",
                # for ORBIT section
                    ['state_eci','(' + r0 + '; ' + r1 + '; ' + r2 + ') (' + v0 + '; ' + v1 + '; ' + v2 + ')' ],
                # for FORCES section
                    gravity_order, # !!!!!!!!!!! put back 20
                "drag solar_pressure sun_gravity moon_gravity", # !!!!!!!!!!!!! put back to "drag sun_gravity moon_gravity"
                ['2017Q2Q3_DSD_converted_for_spock.txt','2017Q2Q3_DGD_converted_for_spock.txt'],#['f107_for_81daverage_20170621_to_20170822.txt', 'ap_20170730_to_20170822.txt'],
                # for OUTPUT section
                        "~/eclipse_velo/"+str(fac_velo).replace(".","_"),
                dt_output, 
                # for ATTITUDE section
                "nadir",
                # for GROUND_STATIONS section
                        "0",
                # for SPICE section
                        "/Users/cbv/cspice/data",
                # FOR #DENSITY_MOD section
                        1
            )

            ## Run SpOCK
            os.system(path_mpirun + ' -np 1 spock_dev_parallel_kalman_9state ' + main_input_filename)
            ## save position and velocity
            #os.system("python state_dev.py ./ " + main_input_filename + " save position velocity")


            # Read the position and velocity predicted by SpOCK
            isc = 0
            var_in, var_in_order = read_input_file(main_input_filename)

            output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
            output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
            var_to_read = ["position", "velocity"]
            var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
            date_spock = np.array(var_out[find_in_read_input_order_variables(var_out_order, 'date')])
            r_spock = var_out[find_in_read_input_order_variables(var_out_order, 'position')]
            v_spock = var_out[find_in_read_input_order_variables(var_out_order, 'velocity')]
            n_spock = len(date_spock)

            # Compare SpOCK and data
            # Assumption: SpOCK was run with a 1s time step to avoid having to do interpolation here: the steps in SpOCK falls at the same time as the steps in data 
            ## Select the time where date_spock = date_data 
            if irho == 0:
                index_spock_same_date_as_data = []
                j = index_interval[itime]
                nb_seconds_since_start_itime = []
                for i in range(n_spock):
                    if date_spock[i] == date_data[j]:
                        index_spock_same_date_as_data.append(i)
                        nb_seconds_since_start_itime.append( ( datetime.strptime(date_data[j],"%Y/%m/%d %H:%M:%S") - datetime.strptime(date_data[0],"%Y/%m/%d %H:%M:%S") ).total_seconds() )
                        j = j + 1
                n = j-index_interval[itime]
                nb_seconds_since_start.append( nb_seconds_since_start_itime )
                date_spock_ok = date_spock[index_spock_same_date_as_data]
            r_spock_ok = np.zeros([n, 3])
            r_spock_ok[:, 0] = r_spock[index_spock_same_date_as_data, 0]
            r_spock_ok[:, 1] = r_spock[index_spock_same_date_as_data, 1]
            r_spock_ok[:, 2] = r_spock[index_spock_same_date_as_data, 2]
            v_spock_ok = np.zeros([n, 3])
            v_spock_ok[:, 0] = v_spock[index_spock_same_date_as_data, 0]
            v_spock_ok[:, 1] = v_spock[index_spock_same_date_as_data, 1]
            v_spock_ok[:, 2] = v_spock[index_spock_same_date_as_data, 2]

            distance_sub = []
            index_data = index_interval[itime]
            for i in range(n):
                distance_sub.append( np.linalg.norm(r_data[index_data, :] - r_spock_ok[i, :]) )
                index_data = index_data + 1
            distance_interval.append( distance_sub )
        distance.append( distance_interval )

#nb_seconds_since_start = np.array(nb_seconds_since_start)
mean_dist_itime_irho = np.zeros([nb_interval-1-itime_start, nb_rho]) # min of the distance for a given internval and a given rho_mod. This is what we want to mnimize using the optimum rho_mod
which_rho_min_dist = np.zeros([nb_interval-1-itime_start])
min_mean_dist_itime_irho =  np.zeros([nb_interval-1-itime_start])
nan_in_run = np.zeros([nb_interval-1-itime_start]) 
for itime in range(0, nb_interval-1-itime_start):
    for irho in range(nb_rho):
        dist_itime_irho = np.array(distance[itime][irho])
        mean_dist_itime_irho[itime, irho] = np.mean(dist_itime_irho)
        # Problem: nan if SpOCK crashed (which doesn't seem to happen anymore)... -> ignore these runs by replacing the nan with a big value (1e8 for instance)
        if np.isnan(mean_dist_itime_irho[itime, irho]) == True:
            mean_dist_itime_irho[itime, irho] = 1e8
            nan_in_run[itime] = nan_in_run[itime] + 1
    which_rho_min_dist[itime] = rho_mod_arr[np.where( mean_dist_itime_irho[itime, :] ==  np.min(mean_dist_itime_irho[itime, :]) )[0][0]]
    min_mean_dist_itime_irho[itime] = np.min(mean_dist_itime_irho[itime, :])

# Plot

## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3


# # OPTIMUM RHO
# fig_title = 'Optimum density coefficient k as a function of time'
# y_label = 'Optimum k'
# x_label = 'Real time'

# fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
# fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
# plt.rc('font', weight='bold') ## make the labels of the ticks in bold
# gs = gridspec.GridSpec(1, 1)
# gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
# ax = fig.add_subplot(gs[0, 0])

# ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
# ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

# [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
# ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
# plt.rc('font', weight='bold') ## make the labels of the ticks in bold

# x_axis = nb_seconds_between_epochyyddd_and_new_date_start
# ax.plot(x_axis, which_rho_min_dist, linewidth = 2, color = 'b')
# ax.scatter(x_axis, which_rho_min_dist, linewidth = 2, color = 'b')


# # x axis label is in real time
# nb_seconds_in_simu = nb_seconds_between_epochyyddd_and_new_date_start[-1] - nb_seconds_between_epochyyddd_and_new_date_start[0]
# start_xaxis_label = nb_seconds_between_epochyyddd_and_new_date_start[0]
# date_ref = date_start_simu#epochyyddd_date
# nb_ticks_xlabel = 10
# dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
# xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
# date_list_str = []
# date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
# for i in range(len(xticks)):
#     if dt_xlabel > nb_ticks_xlabel*24*3600:
#         date_list_str.append( str(date_list[i])[5:10] )
#     else:
#         date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
#         ax.xaxis.set_ticks(xticks)
#         ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
#         ax.margins(0,0); ax.set_xlim([min(xticks), max(xticks)])
# #        ax.set_xlim([ax.get_xlim()[0], most_recent_tle_among_all_sc])

# legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="", fontsize = fontsize_plot)
# #legend.get_title().set_fontsize(str(fontsize_plot))


# fig_save_name = 'new_optimum_rho_coeff_over_' + str(interval) + 'hours.pdf'
# fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# # MIN OF MEAN DISTANCE FOR EACH TIME
# fig_title = 'Average distance error over ' + str(interval) + ' hours with optimum density coefficient'
# y_label = 'Distance error (m)'
# x_label = 'Real time'

# fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
# fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
# plt.rc('font', weight='bold') ## make the labels of the ticks in bold
# gs = gridspec.GridSpec(1, 1)
# gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
# ax = fig.add_subplot(gs[0, 0])

# ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
# ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

# [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
# ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
# plt.rc('font', weight='bold') ## make the labels of the ticks in bold

# x_axis = nb_seconds_between_epochyyddd_and_new_date_start
# ax.plot(x_axis, min_mean_dist_itime_irho*1000, linewidth = 2, color = 'b')
# ax.scatter(x_axis, min_mean_dist_itime_irho*1000, linewidth = 2, color = 'b')


# # x axis label is in real time
# nb_seconds_in_simu = nb_seconds_between_epochyyddd_and_new_date_start[-1] - nb_seconds_between_epochyyddd_and_new_date_start[0]
# start_xaxis_label = nb_seconds_between_epochyyddd_and_new_date_start[0]
# date_ref = date_start_simu#epochyyddd_date
# nb_ticks_xlabel = 10
# dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
# xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
# date_list_str = []
# date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
# for i in range(len(xticks)):
#     if dt_xlabel > nb_ticks_xlabel*24*3600:
#         date_list_str.append( str(date_list[i])[5:10] )
#     else:
#         date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
#         ax.xaxis.set_ticks(xticks)
#         ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
#         ax.margins(0,0); ax.set_xlim([min(xticks), max(xticks)])
# #        ax.set_xlim([ax.get_xlim()[0], most_recent_tle_among_all_sc])

# legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="", fontsize = fontsize_plot)
# #legend.get_title().set_fontsize(str(fontsize_plot))


# fig_save_name = 'new_min_mean_dist_over_' + str(interval) + 'hours.pdf'
# fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  







#######
itime = 0
fig_title = 'Distance between SpOCK and data for different density coefficient'
y_label = 'Distance (km)'
x_label = 'Real time'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

x_axis = nb_seconds_since_start[itime]
if nb_rho > 1:
    alpha_arr = np.arange(0.2,1+0.2/nb_rho,(1-0.2)/(nb_rho-1))
else: 
    alpha_arr = [1]
for irho in range(nb_rho):
    if alpha_arr[irho] >1:
        alpha_arr[irho] = 1
    if mean_dist_itime_irho[itime,irho] < 1e8: # otherwise SpOCK crashed for this run...:
        ax.plot(x_axis, distance[itime][irho], linewidth = 2, color = 'b', alpha = alpha_arr[irho])
        ax.text(x_axis[-1], distance[itime][irho][-1], str(rho_mod_arr[irho]), horizontalalignment = 'left', fontsize = fontsize_plot, weight = 'bold', color = 'b', alpha = alpha_arr[irho], verticalalignment = 'center')

# x axis label is in real time
nb_seconds_in_simu = nb_seconds_since_start[itime][-1] - nb_seconds_since_start[itime][0]
start_xaxis_label = nb_seconds_since_start[itime][0]
date_ref = datetime.strptime(date_data[index_interval[itime+itime_start]],"%Y/%m/%d %H:%M:%S")
nb_ticks_xlabel = 10
dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
xticks = np.arange(start_xaxis_label, start_xaxis_label+nb_seconds_in_simu+1, dt_xlabel)
date_list_str = []
date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
for i in range(len(xticks)):
    if dt_xlabel > nb_ticks_xlabel*24*3600:
        date_list_str.append( str(date_list[i])[5:10] )
    else:
        date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
        ax.xaxis.set_ticks(xticks)
        ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
        ax.margins(0,0); ax.set_xlim([min(xticks), max(xticks)])
#        ax.set_xlim([ax.get_xlim()[0], most_recent_tle_among_all_sc])

legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="", fontsize = fontsize_plot)
#legend.get_title().set_fontsize(str(fontsize_plot))


fig_save_name = output_file_name_list[0].replace("txt","pdf")#'new_all_rho_coeff_itime_' + new_date_start[itime+itime_start].replace("/", "-").replace(" ", "T").replace(":","_") + '_TO_' +  new_date_start[itime+1+itime_start].replace("/", "-").replace(" ", "T").replace(":","_") + '_itime_' +  str(itime) + '.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
