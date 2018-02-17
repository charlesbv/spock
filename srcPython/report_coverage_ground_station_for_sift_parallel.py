
#modified by cbv on 2018-01-12 for local version of SIFT
# This script shows the access start and end times when the spacecraft flies over the ground station, in the half cone defined by the minimum elevation angle. A report is created in the subfolder of the satellite that includes this information for all ground stations.
# To run it: 
# run report_coverage_ground_station_parallel.py main_input_filename
# where main_input_filename is the name of the main input file that were used to run SpOCK. 
# You can also plot the results by setting plot_results to 1, and path_folder_plot to the path where to put the plots (one plot per satellite). The name of the plot is the same as the name of the output folder of the satellite ('pdf' format)
import sys
from datetime import datetime, timedelta
import numpy as np
from read_input_file import *
from os import listdir
from os.path import isfile, join


def report_coverage_ground_station_for_sift_parallel(main_input_file_name):

    # PARAMETERS TO SET
    plot_results = 0 # if you want to create a plot of the coverage
    path_folder_plot = '.' # if plot_results is set to 1, otherwise don't worry about it 
    # end of PARAMETERS TO SET

    # Algorithm
    if plot_results == 1:
        from matplotlib import pyplot as plt
        import matplotlib.gridspec as gridspec
        plt.ion()
    deg_sign = u'\xb0'.encode('utf8')


    #main_input_file_name = sys.argv[1]

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
                overflight_per_sat_per_overflight = []
                overflight_per_sat_per_overflight.append( read_file_station[istep + skip_header].split()[0] )
                current_date = datetime.strptime(read_file_station[istep + skip_header].split()[0], "%Y-%m-%dT%H:%M:%S")
                next_date = datetime.strptime(read_file_station[istep + skip_header + 1].split()[0], "%Y-%m-%dT%H:%M:%S")
                delta_date = ( next_date - current_date ).total_seconds()
                while ( ( istep + skip_header  < len(read_file_station) - 1) & ( delta_date <= dt ) ):
                    istep = istep + 1
                    if istep + skip_header + 1 == len(read_file_station):
                        break

                    current_date = datetime.strptime(read_file_station[istep + skip_header].split()[0], "%Y-%m-%dT%H:%M:%S")
                    next_date = datetime.strptime(read_file_station[istep + skip_header + 1].split()[0], "%Y-%m-%dT%H:%M:%S")
                    delta_date = ( next_date - current_date ).total_seconds()

                overflight_per_sat_per_overflight.append( read_file_station[istep + skip_header ].split()[0] )
                overflight_per_sat.append( overflight_per_sat_per_overflight )
                overflight_per_sat_per_overflight.append(all_station_files[istation].split('/')[-1].split('_by')[0])
                overflight_all_stations_combined.append( overflight_per_sat_per_overflight )
                new_coverage = 0
            overflight.append( overflight_per_sat )


        #     # Write report: HERE IT WIL CREATE NB_GROUND_STATIONS BLOCKS: COVERAGE OF EACH GROUND STATIONS
        #     print >> file_station_out, '#Coverage of ' + all_station_files[istation].split('/')[-1].split('_by')[0] + ' (longitude = ' + '{0:.2f}'.format(station_lon[istation]) + deg_sign + ', latitude = ' + '{0:.2f}'.format(station_lat[istation]) + deg_sign  + ', altitude = ' + '{0:.0f}'.format(station_alt[istation]) +' m)' + '\n'
        #     for icoverage in range( len(overflight_per_sat) ):
        #         duration_overflight = datetime.strptime( overflight_per_sat[icoverage][1], "%Y-%m-%dT%H:%M:%S" ) - datetime.strptime( overflight_per_sat[icoverage][0], "%Y-%m-%dT%H:%M:%S" ) 
        #         print >> file_station_out, overflight_per_sat[icoverage][0].replace("T", ' ') + ' to ' + overflight_per_sat[icoverage][1].replace("T", ' ') + ': ' + str(duration_overflight.days) + 'd ' + str( duration_overflight.seconds / 3600 ) + 'h' + str( duration_overflight.seconds / 60 ) + "'" + str( np.mod(duration_overflight.seconds,60) ) + '"'
        #     print >> file_station_out, ''
        #     overflight_per_sat_all_station.append(overflight_per_sat)
        # file_station_out.close()
        # overflight_all_sat_all_station.append(overflight_per_sat_all_station)
        # # end of Write report: HERE IT WIL CREATE NB_GROUND_STATIONS BLOCKS: COVERAGE OF EACH GROUND STATIONS

        # Write report:  HERE IT WIL CREATE 1 BLOCK: COVERAGE OF ALL STATIONS AS A FUNCTION OF TIME
        nb_overflights = len(overflight_all_stations_combined)
        ## Order from older to most recent coverage
        overflight_all_stations_combined_arr = np.array(overflight_all_stations_combined)
        overflight_all_stations_combined_sorted = [ datetime.strptime(t, '%Y-%m-%dT%H:%M:%S') for t in overflight_all_stations_combined_arr[:,0] ]
        overflight_all_stations_combined_sorted.sort()
        overflight_all_stations_combined_sorted = np.array(overflight_all_stations_combined_sorted)
        print >> file_station_out, "#Coverage of " + satellite_to_plot.split('.')[0] + ":\n"
        for ioverflight in range(nb_overflights):
            duration_overflight = datetime.strptime(overflight_all_stations_combined_arr[np.where( overflight_all_stations_combined_arr[:,0] == datetime.strftime(overflight_all_stations_combined_sorted[ioverflight], "%Y-%m-%dT%H:%M:%S" ) )[0][0] ,1] , "%Y-%m-%dT%H:%M:%S" ) - overflight_all_stations_combined_sorted[ioverflight]
            print >> file_station_out, '(' + str(ioverflight + 1) + ') ' + datetime.strftime(overflight_all_stations_combined_sorted[ioverflight], "%Y-%m-%d %H:%M:%S") + ' to ' + overflight_all_stations_combined_arr[np.where( overflight_all_stations_combined_arr[:,0] == datetime.strftime(overflight_all_stations_combined_sorted[ioverflight], "%Y-%m-%dT%H:%M:%S" ) )[0][0] ,1].replace("T", " ") + ' for ' + str(duration_overflight.days) + 'd ' + str( duration_overflight.seconds / 3600 ) + 'h' + str( duration_overflight.seconds / 60 ) + "'" + str( np.mod(duration_overflight.seconds,60) ) + '" (' + str(duration_overflight.days*24*3600 + duration_overflight.seconds) + ' s)' + ' - ' + overflight_all_stations_combined_arr[np.where( overflight_all_stations_combined_arr[:,0] == datetime.strftime(overflight_all_stations_combined_sorted[ioverflight], "%Y-%m-%dT%H:%M:%S" ) )[0][0] ,2].title() 
        file_station_out.close()
        # end of Write report:  HERE IT WIL CREATE 1 BLOCK: COVERAGE OF ALL STATIONS AS A FUNCTION OF TIME

    # Below works for script report_coverage_ground_station.py but haven't tried it for report_coverage_ground_station_for_sift.py
    # if plot_results == 1:
    #     # PLOT COVERAGE: 1 figure per satellite
    #     width_fig = 15.
    #     height_fig = width_fig * 3 /4
    #     fontsize_plot = 20 # 9
    #     color_arr = ['k','b','r','g','m', 'y']
    #     fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
    #     gs = gridspec.GridSpec(1, 1)
    #     gs.update(left=0.12, right=0.97, top = 0.90,bottom = 0.06)
    #     fig.suptitle('', fontsize = 22)
    #     plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    #     ax1 = fig.add_subplot(gs[0, 0])
    #     #ax1 = fig.add_subplot(111)
    #     ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
    #     ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
    #     #ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
    #     [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

    #     for isat in range(nb_satellites):
    #         date_ref = date_start
    #         yticks = np.arange(0, 1, 1./nb_stations)
    #         ytick_labels = []
    #         for istation in range(nb_stations):
    #             for icoverage in range( len(overflight_all_sat_all_station[isat][istation]) ):
    #                 start_coverage_since_ref = ( datetime.strptime( overflight_all_sat_all_station[isat][istation][icoverage][1], "%Y-%m-%dT%H:%M:%S" ) - date_ref ).total_seconds()
    #                 end_coverage_since_ref = ( datetime.strptime( overflight_all_sat_all_station[isat][istation][icoverage][0], "%Y-%m-%dT%H:%M:%S" ) - date_ref ).total_seconds()
    #                 ax1.plot([start_coverage_since_ref, end_coverage_since_ref], [istation * 1./nb_stations, istation * 1./nb_stations], linewidth = 20, color = color_arr[istation])
    #     #            print [start_coverage_since_ref, end_coverage_since_ref], [istation * 1./nb_stations, istation * 1./nb_stations]
    #             ytick_labels.append(overflight_all_sat_all_station[isat][istation][0][2].title())
    #         ax1.yaxis.set_ticks(yticks)
    #         ax1.yaxis.set_ticklabels(ytick_labels, fontsize = fontsize_plot)#, rotation='vertical')
    #         for ytick, color in zip(ax1.get_yticklabels(), color_arr):
    #             ytick.set_color(color)

    #         max_xaxis = (date_stop - date_ref).total_seconds()

    #         hour_time_step_xticks = 4
    #         hour_time_step_xticks_converted_in_seconds = hour_time_step_xticks * 3600
    #         xticks = np.arange(0, max_xaxis, hour_time_step_xticks_converted_in_seconds)
    #         date_list_str = []
    #         date_list = [date_start + timedelta(seconds=x) for x in np.arange(0, max_xaxis, hour_time_step_xticks_converted_in_seconds)]
    #         for i in range(len(xticks)):
    #             date_list_str.append( str(date_list[i])[5:10] + '\n' + str(date_list[i])[11:16])
    #         ax1.xaxis.set_ticks(xticks) 
    #         ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')

    #         ax1.set_xlim([0,max_xaxis])
    #         ax1.set_ylim([-1./nb_stations,1])
    #         ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

    #         ax1.set_title('Coverage of ' + str(nb_stations) + ' ground stations by ' + satellite_to_plot_arr[isat].replace(".txt","") + ' - Access times', weight = 'bold', fontsize = fontsize_plot , y = 1.005) 
    #         ax1.legend( fontsize = fontsize_plot, loc = 4)
    #         ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'bold')
    #         ax1.set_ylabel('', fontsize = fontsize_plot, weight = 'bold')

    #         fig_save_name = path_folder_plot + '/' + satellite_to_plot_arr[isat].replace(".txt",".pdf")
    #         fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

    return 
