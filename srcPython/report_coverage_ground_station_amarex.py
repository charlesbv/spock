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

# copy of report_coverage_ground_station_for_sift_parallel_sftp.py on 082719. It assumes that the function coverage_ground_station in SpOCK doesn't output the elevation of the sc wrt to the gs (noted "elevation_wtr_to_ground_station_in_ground_station_refce*RAD2DEG" in SpOCK) but the angle of view of the gs wrt to the sc (noted "90 - (-elevation_wtr_to_ground_station_in_sc_refce*RAD2DEG)" in SpOCK). The angle of view is the angle between the nadir direction and the direction sc to gs. It is 0 if the gs is right below the sc (and the sc is nadir pointing). This script then selects only the times for which the angle of view is smaller than fov, given as the second argument of this script.
#modified by cbv on 2018-01-12 for local version of SIFT
# This script shows the access start and end times when the spacecraft flies over the ground station, in the half cone defined by the minimum elevation angle. A report is created in the subfolder of the satellite that includes this information for all ground stations.
# To run it: 
# run report_coverage_ground_station_parallel.py main_input_filename
# where main_input_filename is the name of the main input file that were used to run SpOCK. 
# You can also plot the results by setting plot_results to 1, and path_folder_plot to the path where to put the plots (one plot per satellite). The name of the plot is the same as the name of the output folder of the satellite ('pdf' format)
# PARAMETERS TO SET
plot_results = 0 # if you want to create a plot of the coverage
save_results = 1 # if you want to save the coverage stats as a pickle
path_folder_plot = '.' # if plot_results is set to 1, otherwise don't worry about it 
# end of PARAMETERS TO SET

import sys
sys.path.append("/Users/cbv/work/spock/srcPython")
from datetime import datetime, timedelta
import numpy as np
from read_input_file import *
from os import listdir
from os.path import isfile, join
import ipdb
import math
import os


if save_results == 1:
    import pickle
    if (os.path.isdir("pickle") == False):
        os.system('mkdir pickle')
deg_sign = u'\xb0'.encode('utf8')

def report_coverage_ground_station_amarex(main_input_file_name, fov):



    # Algorithm

    # main_input_file_name = sys.argv[1]
    # fov = sys.argv[2] # select only the times for which the angle of view is smaller than fov. fov is a numpy array

    # read input file
    input_variables, order_input_variables = read_input_file(main_input_file_name)
    dt = input_variables[find_in_read_input_order_variables(order_input_variables, 'dt')]
    nb_steps = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_steps')]
    satellite_to_plot_path_arr = input_variables[find_in_read_input_order_variables(order_input_variables, 'output_file_path_list')]
    satellite_to_plot_arr = input_variables[find_in_read_input_order_variables(order_input_variables, 'output_file_name_list')]
    filename_station = input_variables[find_in_read_input_order_variables(order_input_variables, 'filename_ground_station')]
    nb_satellites = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_sc')]
    date_start = input_variables[find_in_read_input_order_variables(order_input_variables, 'date_start')]
    date_stop = input_variables[find_in_read_input_order_variables(order_input_variables, 'date_stop')]
    # read coverage ground station files
    overflight_all_sat_all_station = []
    overflight_all_stations_combined_arr_all_sat = []
    for isat in range(nb_satellites):

        overflight_per_sat_all_station = []
        satellite_to_plot_path = satellite_to_plot_path_arr[isat]
        satellite_to_plot = satellite_to_plot_arr[isat]
        all_station_files = [satellite_to_plot_path + 'coverage/' + f for f in listdir(satellite_to_plot_path + 'coverage/') if ( isfile(join(satellite_to_plot_path + 'coverage/', f)) and f.endswith('.txt') and f.startswith("report") == 0 )]
        nb_stations = len(all_station_files)

        filename_station_out = satellite_to_plot_path + 'coverage/report_all_by_' + satellite_to_plot
        file_station_out = open( filename_station_out, "w+" )
        station_lon = np.zeros([nb_stations]); station_lat = np.zeros([nb_stations]); station_alt = np.zeros([nb_stations])
        overflight = []
        overflight_all_stations_combined = []
        aov_max_per_sat = []
        aov_min_per_sat = []
        for istation in range( nb_stations ):
            overflight_per_sat = []
            file_station = open( all_station_files[istation] )
            read_file_station = file_station.readlines()
            # skip header
            skip_header = 0
            while ( read_file_station[skip_header].split()[0] != '#START'):
                if 'longitude' in read_file_station[skip_header]:
                    station_lon[istation] = np.float( read_file_station[skip_header].split(':')[-1].split()[0] )
                    station_lat[istation] = np.float( read_file_station[skip_header].split(':')[-1].split()[1] )
                    station_alt[istation] = np.float( read_file_station[skip_header].split(':')[-1].split()[2] )
                skip_header = skip_header + 1
            skip_header = skip_header + 1
            # read over the file 
            istep = -1
            while  istep + skip_header < len(read_file_station)-1:
                istep = istep + 1
                if istep + skip_header + 1 == len(read_file_station):
                    break
                overflight_per_sat_per_overflight = []
                if ((istation == 10000) & (isat == 1)):
                    ipdb.set_trace()
                current_aov = np.float(read_file_station[istep + skip_header].split()[1]) # angle of view
                aov_start = 0
                if current_aov <= fov[isat]:
                    overflight_per_sat_per_overflight.append( read_file_station[istep + skip_header].split()[0] )
                    aov_start = 1
                current_date = datetime.strptime(read_file_station[istep + skip_header].split()[0], "%Y-%m-%dT%H:%M:%S")
                if ((np.isnan(current_aov) == False) & (aov_start == 1)):
                    aov_max = current_aov;
                    aov_min = current_aov;
                else:
                    aov_max = -1e6
                    aov_min = 1e6
                next_date = datetime.strptime(read_file_station[istep + skip_header + 1].split()[0], "%Y-%m-%dT%H:%M:%S")
                delta_date = ( next_date - current_date ).total_seconds()
                while ( ( istep + skip_header  < len(read_file_station) - 1) & ( delta_date <= dt ) ):
                    istep = istep + 1
                    if (istep + skip_header + 1 == len(read_file_station)):
                        if (aov_start == 1): # the condition aov_start == 1 is to prevent the case where only the last time step of the file has a aov<=fov and no time step before (very rare xase)
                            current_date = datetime.strptime(read_file_station[istep + skip_header].split()[0], "%Y-%m-%dT%H:%M:%S")
                            current_aov = np.float(read_file_station[istep + skip_header].split()[1])
                            if current_aov <= fov[isat]:
                                if ((current_aov >= aov_max) & (np.isnan(current_aov) == False)):
                                    aov_max = current_aov
                                if ((current_aov <= aov_min) & (np.isnan(current_aov) == False)):
                                    aov_min = current_aov
                                overflight_per_sat_per_overflight.append( read_file_station[istep + skip_header].split()[0] )                        
                                overflight_per_sat_per_overflight.append(all_station_files[istation].split('/')[-1].split('_by')[0])
                                overflight_per_sat_per_overflight.append(aov_min)
                                overflight_per_sat_per_overflight.append(aov_max)
                                overflight_all_stations_combined.append( overflight_per_sat_per_overflight )
                                aov_min_per_sat.append(aov_min)
                                aov_max_per_sat.append(aov_max)
                        if ((istation == 10000) & (isat == 1)):                        
                            ipdb.set_trace()

                        break

                    current_date = datetime.strptime(read_file_station[istep + skip_header].split()[0], "%Y-%m-%dT%H:%M:%S")
                    current_aov = np.float(read_file_station[istep + skip_header].split()[1])
                    if current_aov <= fov[isat]:
                        if ((current_aov >= aov_max) & (np.isnan(current_aov) == False)):
                            aov_max = current_aov
                        if ((current_aov <= aov_min) & (np.isnan(current_aov) == False)):
                            aov_min = current_aov
                        if aov_start == 0:
                            if ((istation == 10000) & (isat == 1)):
                                ipdb.set_trace()
                            overflight_per_sat_per_overflight.append( read_file_station[istep + skip_header].split()[0] )                        
                            aov_start = 1
                    elif ((current_aov > fov[isat]) & (aov_start == 1)):
                        if ((istation == 10000) & (isat == 1)):
                            ipdb.set_trace()
                        overflight_per_sat_per_overflight.append( read_file_station[istep + skip_header ].split()[0] )
                        overflight_per_sat_per_overflight.append(all_station_files[istation].split('/')[-1].split('_by')[0])
                        overflight_per_sat_per_overflight.append(aov_min)
                        overflight_per_sat_per_overflight.append(aov_max)
                        overflight_all_stations_combined.append( overflight_per_sat_per_overflight )
                        aov_min_per_sat.append(aov_min)
                        aov_max_per_sat.append(aov_max)
                        aov_start = 0
                        overflight_per_sat_per_overflight = []
                    next_date = datetime.strptime(read_file_station[istep + skip_header + 1].split()[0], "%Y-%m-%dT%H:%M:%S")
                    delta_date = ( next_date - current_date ).total_seconds()
                    if ((delta_date > dt) & (current_aov <= fov[isat])):
                        if ((current_aov >= aov_max) & (np.isnan(current_aov) == False)):
                            aov_max = current_aov
                        if ((current_aov <= aov_min) & (np.isnan(current_aov) == False)):
                            aov_min = current_aov
                        overflight_per_sat_per_overflight.append( read_file_station[istep + skip_header].split()[0] )                        
                        overflight_per_sat_per_overflight.append(all_station_files[istation].split('/')[-1].split('_by')[0])
                        overflight_per_sat_per_overflight.append(aov_min)
                        overflight_per_sat_per_overflight.append(aov_max)
                        overflight_all_stations_combined.append( overflight_per_sat_per_overflight )
                        aov_min_per_sat.append(aov_min)
                        aov_max_per_sat.append(aov_max)
                        if ((istation == 10000) & (isat == 1)):
                            ipdb.set_trace()

                new_coverage = 0
        overflight_all_sat_all_station.append(overflight)

        # if ((isat == 1)):
        #     ipdb.set_trace()

        # Write report:  HERE IT WIL CREATE 1 BLOCK: COVERAGE OF ALL STATIONS AS A FUNCTION OF TIME

        nb_overflights = len(overflight_all_stations_combined)
        ## Order from older to most recent coverage
        overflight_all_stations_combined_arr = np.array(overflight_all_stations_combined)
        overflight_all_stations_combined_sorted_temp = np.array([ datetime.strptime(t, '%Y-%m-%dT%H:%M:%S') for t in overflight_all_stations_combined_arr[:,0] ])
        index_sort = np.argsort(overflight_all_stations_combined_sorted_temp)
        overflight_all_stations_combined_arr_all_sat.append(overflight_all_stations_combined_arr[index_sort, :])
        overflight_all_stations_combined_sorted = overflight_all_stations_combined_sorted_temp[index_sort]
        print >> file_station_out, "#Coverage of " + satellite_to_plot.split('.')[0] + ":\n"

        for ioverflight in range(nb_overflights):
            icov = index_sort[ioverflight]
            duration_overflight = (datetime.strptime(overflight_all_stations_combined_arr[icov, 1], "%Y-%m-%dT%H:%M:%S" ) - datetime.strptime(overflight_all_stations_combined_arr[icov, 0], "%Y-%m-%dT%H:%M:%S" )).total_seconds()
            duration_overflight_day = (int)(duration_overflight / 3600./ 24.)
            duration_overflight_hour = (int)((duration_overflight - duration_overflight_day * 24 * 3600)/3600.)
            duration_overflight_min = (int)((duration_overflight - (duration_overflight_day * 24. * 3600 + duration_overflight_hour * 3600.))/60.)
            duration_overflight_sec = (int)(duration_overflight - (duration_overflight_day * 24. * 3600 + duration_overflight_hour * 3600 + duration_overflight_min * 60))
            print >> file_station_out, '(' + str(ioverflight + 1) + ') ' + overflight_all_stations_combined_arr[icov, 0].replace("T", " ") + ' to ' + overflight_all_stations_combined_arr[icov, 1].replace("T", " ") + ' for ' + str(duration_overflight_day) + 'd ' + str(duration_overflight_hour) + 'h' + str(duration_overflight_min) + "'" + str(duration_overflight_sec) + '" (' + str(duration_overflight) + ' s)' + '- aov ' + format(aov_min_per_sat[icov], ".1f") + '-' + format(aov_max_per_sat[icov], ".1f") + ' - '  + overflight_all_stations_combined_arr[icov, 2]
    #        raise Exception
        file_station_out.close()

        # end of Write report:  HERE IT WIL CREATE 1 BLOCK: COVERAGE OF ALL STATIONS AS A FUNCTION OF TIME

    
    if save_results == 1:
        pickle.dump( overflight_all_stations_combined_arr_all_sat, open( 'pickle/' + main_input_file_name.replace('.txt', '.pickle'), "w" ) )    

    # Below works for script report_coverage_ground_station.py but haven't tried it for report_coverage_ground_station_for_sift.py
    if plot_results == 1:
        # PLOT COVERAGE: 1 figure per satellite
        width_fig = 15.
        height_fig = width_fig * 3 /4
        fontsize_plot = 20 # 9
        color_arr = ['k','b','r','g','m', 'y']
        fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
        gs = gridspec.GridSpec(1, 1)
        gs.update(left=0.12, right=0.97, top = 0.90,bottom = 0.06)
        fig.suptitle('', fontsize = 22)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        ax1 = fig.add_subplot(gs[0, 0])
        #ax1 = fig.add_subplot(111)
        ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
        ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
        #ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

        for isat in range(nb_satellites):
            date_ref = date_start
            yticks = np.arange(0, 1, 1./nb_stations)
            ytick_labels = []
            for istation in range(nb_stations):
                for icoverage in range( len(overflight_all_sat_all_station[isat][istation]) ):
                    start_coverage_since_ref = ( datetime.strptime( overflight_all_sat_all_station[isat][istation][icoverage][1], "%Y-%m-%dT%H:%M:%S" ) - date_ref ).total_seconds()
                    end_coverage_since_ref = ( datetime.strptime( overflight_all_sat_all_station[isat][istation][icoverage][0], "%Y-%m-%dT%H:%M:%S" ) - date_ref ).total_seconds()

                    ax1.plot([start_coverage_since_ref, end_coverage_since_ref], [istation * 1./nb_stations, istation * 1./nb_stations], linewidth = 20, color = color_arr[istation])
        #            print [start_coverage_since_ref, end_coverage_since_ref], [istation * 1./nb_stations, istation * 1./nb_stations]
                ytick_labels.append(overflight_all_sat_all_station[isat][istation][0][2].title())
            ax1.yaxis.set_ticks(yticks)
            ax1.yaxis.set_ticklabels(ytick_labels, fontsize = fontsize_plot)#, rotation='vertical')
            for ytick, color in zip(ax1.get_yticklabels(), color_arr):
                ytick.set_color(color)

            max_xaxis = (date_stop - date_ref).total_seconds()

            hour_time_step_xticks = 4
            hour_time_step_xticks_converted_in_seconds = hour_time_step_xticks * 3600
            xticks = np.arange(0, max_xaxis, hour_time_step_xticks_converted_in_seconds)
            date_list_str = []
            date_list = [date_start + timedelta(seconds=x) for x in np.arange(0, max_xaxis, hour_time_step_xticks_converted_in_seconds)]
            for i in range(len(xticks)):
                date_list_str.append( str(date_list[i])[5:10] + '\n' + str(date_list[i])[11:16])
            ax1.xaxis.set_ticks(xticks) 
            ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')

            ax1.set_xlim([0,max_xaxis])
            ax1.set_ylim([-1./nb_stations,1])
            ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

            ax1.set_title('Coverage of ' + str(nb_stations) + ' ground stations by ' + satellite_to_plot_arr[isat].replace(".txt","") + ' - Access times', weight = 'bold', fontsize = fontsize_plot , y = 1.005) 
            ax1.legend( fontsize = fontsize_plot, loc = 4)
            ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
            ax1.set_ylabel('', fontsize = fontsize_plot, weight = 'bold')

            fig_save_name = path_folder_plot + '/' + satellite_to_plot_arr[isat].replace(".txt",".pdf")
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


    return 0
