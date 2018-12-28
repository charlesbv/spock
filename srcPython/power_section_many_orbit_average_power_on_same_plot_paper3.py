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
import matplotlib
import sys
sys.path.append('../')
from operator import itemgetter
from itertools import groupby
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
from operator import itemgetter
from degree_to_time import *
from read_input_file import *
from read_output_file import *
from get_prop_dir import *
from compute_power import *
from datetime import datetime, timedelta
from orbit_average import *
import sys
plt.ion()

input_dir = 'run_cubesat_1u'
input_filename_arr = ['cubesat_1u_polar_with_wings_0_degree.txt', 'cubesat_1u_polar_with_wings_30_degree.txt', 'cubesat_1u_polar_with_wings_60_degree.txt', 'cubesat_1u_polar_with_wings_90_degree.txt']

loop_i = -1
for input_filename in input_filename_arr:
    loop_i = loop_i + 1
    print input_filename
    # READ PROPAGATOR INPUT FILE
    input_filename = get_prop_dir(1) + input_dir +  '/' + 'input/main_input/' + input_filename
    input_variables, order_input_variables = read_input_file(input_filename)
    date_start = input_variables[0]; date_stop = input_variables[1]; dt = input_variables[2]; nb_steps = input_variables[3]; nb_satellites = input_variables[4]; 
    output_path_propagator = input_variables[6]; output_file_propagator = input_variables[7]; nb_surfaces = input_variables[8]

    # !!!!!!!!!! SET THE PARAMETER BELOW: path of the folder where you want to store the results (pickle, image, video)
    path_folder_results = '/raid3/Armada/Charles/python/' #get_prop_dir(2) + 'output/python_propagator/'

    # Name of the mission of the simulation (choices: 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', 'other')
    if ( 'cygnss' in output_file_propagator[0].lower() ):
        name_mission = 'CYGNSS' 
    elif ( 'cadre' in output_file_propagator[0].lower() ):
        name_mission = 'CADRE' 
    elif ( 'aerie' in output_file_propagator[0].lower() ):
        name_mission = 'AERIE' 
    elif ( 'scion' in output_file_propagator[0].lower() ):
        name_mission = 'SCION' 
    elif ( 'qb50' in output_file_propagator[0].lower() ):
        name_mission = 'QB50' 
    else:
        name_mission = 'other' 

    save_pickle_name = path_folder_results + name_mission + '/pickle/' + output_file_propagator[0].replace(".txt",".pickle")
    root_save_fig_name = path_folder_results + name_mission + '/result/image/' + output_file_propagator[0].replace(".txt","_") 
    root_save_video_name = path_folder_results + name_mission + '/result/video/' + output_file_propagator[0].replace(".txt","_")


    # COMPUTE POWER
    power = np.zeros([nb_satellites,nb_steps, nb_surfaces])
    for isat in range(nb_satellites):
        power[isat, :, :] = compute_power(output_path_propagator[isat] +'power_' + output_file_propagator[isat])[0]
    total_power = np.zeros([nb_satellites,nb_steps])
    for isat in range(nb_satellites):
        for istep in range(nb_steps):
            total_power[isat, istep] = np.sum( power[isat, istep, :])

    # READ LATITUDE TO THEN COMPUTE ORBIT-AVERAGE POWER
    latitude = np.zeros([nb_satellites, nb_steps]); local_time = np.zeros([nb_satellites, nb_steps])
    list_var = ["latitude", "local_time"]
    lat_var_nb = 1; local_time_var_nb = 2
    for isat in range(nb_satellites):
        out_var, order_var = read_output_file( output_path_propagator[isat] + output_file_propagator[isat] , list_var )
        if isat == 0:
            date = out_var[0]
        latitude[isat, :] = out_var[lat_var_nb]
        local_time[isat, :] = out_var[local_time_var_nb]
    date_array = np.array(date)

    # COMPUTE ORBIT-AVERAGE POWER 
    power_average = []
    time_average = []
    nb_orbits = np.zeros([nb_satellites])
    index_time_average = []
    for isat in range(nb_satellites):
        power_average_sub = []; time_average_sub = []; index_time_average_sub = []
        power_average_sub, time_average_sub, index_time_average_sub = orbit_average(total_power[isat, :], latitude[isat, :], date)
        nb_orbits[isat] = len(power_average_sub)
        power_average.append(power_average_sub); time_average.append(time_average_sub); index_time_average.append(index_time_average_sub)
    power_average = np.array(power_average)    
    index_time_average = np.array(index_time_average) # index_time_average[isat, iorbit, pos]: shows the index in date of orbit iorbit for satellite isat. pos is either 0 (beginning of orbit), 1 (middle of orbit), 2 (end of orbit)
    time_average_array = np.array(time_average)

    #  COMPUTE LOCAL_TIME OF THE ASCENDING NODE (LTAN)
    ltan = []
    for isat in range(nb_satellites):
        ltan_sub_temp = degree_to_time(local_time[ isat, index_time_average[isat,:,0]], 0)
        ltan.append( ltan_sub_temp )

    show_local_min_max = 0
    if show_local_min_max == 1:
        # FIND LOCAL MAXIMUM (HELP WITH THE EYES)
        threshold = [1] # array of threshold per satellite
        new_threshold = [1.1] # array of threshold per satellite to remove peaks due to noise
        local_max = []
        for isat in range(nb_satellites):
            local_max_sub = []
            index_where_power_average_above_threshold = np.where( power_average[isat,:] > threshold[isat] )[0]
            list_where_power_average_above_threshold = []
            for k, g in groupby(enumerate(index_where_power_average_above_threshold), lambda (i, x): i-x):
                list_where_power_average_above_threshold.append(map(itemgetter(1), g))
            for  i in range(len(list_where_power_average_above_threshold)):
                power_local = power_average[isat, list_where_power_average_above_threshold[i][0]: list_where_power_average_above_threshold[i][-1]+1] 
                local_max_sub.append( list_where_power_average_above_threshold[i][0] + np.where( power_local == np.max(power_local) )[0][0] )
            # REMOVE PEAKS DUE TO NOISE
            local_max_copy = local_max_sub[:]
            for imax in range(len(local_max_copy)):
                if power_average[isat,local_max_copy[imax]] < new_threshold[isat]:
                    local_max_sub.remove(local_max_copy[imax])
            local_max.append(local_max_sub) # local_max[isat][imax]: shows the index in time_average (or power_average or any orbit-average variable) of the maximum number imax. So it shows the orbit number of this max


        # FIND LOCAL MINIMUM (HELP WITH THE EYES)
        threshold = [0.6] # array of threshold per satellite
        new_threshold = [0.5] # array of threshold per satellite to remove peaks due to noise
        local_min = []
        for isat in range(nb_satellites):
            local_min_sub = []
            index_where_power_average_above_threshold = np.where( power_average[isat,:] < threshold[isat] )[0]
            list_where_power_average_above_threshold = []
            for k, g in groupby(enumerate(index_where_power_average_above_threshold), lambda (i, x): i-x):
                list_where_power_average_above_threshold.append(map(itemgetter(1), g))
            for  i in range(len(list_where_power_average_above_threshold)):
                power_local = power_average[isat, list_where_power_average_above_threshold[i][0]: list_where_power_average_above_threshold[i][-1]+1] 
                local_min_sub.append( list_where_power_average_above_threshold[i][0] + np.where( power_local == np.min(power_local) )[0][0] )
            # REMOVE PEAKS DUE TO NOISE
            local_min_copy = local_min_sub[:]
            for imin in range(len(local_min_copy)):
                if power_average[isat,local_min_copy[imin]] > new_threshold[isat]:
                    local_min_sub.remove(local_min_copy[imin])
            local_min.append(local_min_sub)


    # PLOT
    ## ORBIT-AVERAGE POWER OVER A YEAR
    if loop_i == 0:
        color_arr = ['k','g','r','b']
        label_arr = ['0'+ u'\N{DEGREE SIGN}', '30'+ u'\N{DEGREE SIGN}', '60'+ u'\N{DEGREE SIGN}', '90'+ u'\N{DEGREE SIGN}']
        width_fig = 15
        height_fig = width_fig * 3./4
        fontsize_plot = 20 # 9
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

    isat = 0
    ax1.plot( power_average[isat],'k', linewidth = 2 , color = color_arr[loop_i], label  = label_arr[loop_i])
    ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
    show_local_min_max = 0
    if show_local_min_max == 1:
        for imax in range(len(local_max[isat])):
            ax1.text(local_max[isat][imax], ax1.get_ylim()[1] + (ax1.get_ylim()[1] - ax1.get_ylim()[0])/200., ltan[isat][local_max[isat][imax]][0:5], weight = 'bold', horizontalalignment = 'center', fontsize = fontsize_plot)
        for imin in range(len(local_min[isat])):
            ax1.text(local_min[isat][imin], power_average[isat, local_min[isat][imin]] + (ax1.get_ylim()[1] - ax1.get_ylim()[0])/15, ltan[isat][local_min[isat][imin]][0:5], weight = 'bold', horizontalalignment = 'center', fontsize = fontsize_plot)
            print ltan[isat][local_min[isat][imin]][0:5]

    x_max = index_time_average[isat,-1,1]
    x_min = index_time_average[isat,0,1]
    delta_x_seconds = ( x_max - x_min ) * dt
    hour_time_step_xticks = 24*3*30
    second_time_step_xticks = hour_time_step_xticks * 3600
    xticks = np.arange(0, len(power_average[isat]), len(power_average[isat]) / delta_x_seconds * second_time_step_xticks)
    date_list_str = []
    date_start = datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S")
    date_list = [date_start + timedelta(hours=x) for x in np.arange(0, 732*24, hour_time_step_xticks)]
    for i in range(len(xticks)):
        date_list_str.append( str(date_list[i])[0:10] )
    ax1.xaxis.set_ticks(xticks)
    ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')

    if loop_i == 0:
        ### SHADING FOR SEASONS
        spring = date.index( '2014/03/21 12:00:00' )  * dt / delta_x_seconds
        summer = date.index( '2014/06/21 00:00:00' )  * dt / delta_x_seconds * len(power_average[isat])
        fall = date.index( '2014/09/21 00:00:00' )  * dt / delta_x_seconds * len(power_average[isat])
        winter = date.index( '2014/12/21 00:00:00' )  * dt / delta_x_seconds * len(power_average[isat])
        ax1.axvspan( spring, summer, color = 'g', alpha = 0.2)#, label = 'Spring')
        ax1.axvspan( summer, fall, color = 'y', alpha = 0.2)#, label = 'Summer')
        ax1.axvspan( fall, winter, color = 'r', alpha = 0.2)#, label = 'Fall')
        ax1.axvspan( winter, len(power_average[isat]-1), color = 'b', alpha = 0.2)#, label = 'Winter')
        
        ax1.set_ylim([0,ax1.get_ylim()[1]])
        
        ax1.set_title('Cubesat 1U Polar With Wings Different Angles orbit-average power', weight = 'bold', fontsize = 20,  y = 1.01) 
        ax1.set_xlabel('Real Time', fontsize = fontsize_plot, weight = 'bold')
        ax1.set_ylabel('Power (W)', fontsize = fontsize_plot, weight = 'bold')

ax1.legend(loc = 3, fontsize = fontsize_plot)
fig_save_name = '/raid3/Armada/Charles/python/other/result/image/cubesat_1u_polar_with_wings_different_angles_power_year.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
os.system("rsync -av " + fig_save_name +' srbwks2014-0008.engin.umich.edu:./' + name_mission)
