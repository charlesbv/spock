from get_prop_dir import *
import matplotlib.gridspec as gridspec
import pickle
from datetime import datetime, timedelta
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from read_input_file import *
from read_output_file import *
import colorsys


plt.ion()
plt.isinteractive()

# ASSUMPTION: only one ground station for now


# TO BE SET BY THE USER BEFORE RUNNING
#input_filename = "/home/cbv/PropSim/input/main_input/coverage_ground_station.txt"
input_filename = "/home/cbv/PropSim/input/main_input/scott.txt"
time_step_interpolation = 1. # in seconds

# READ THE INPUT FILE OF PROPAGATOR
# input_variables, order_input_variables = read_input_file(input_filename)
# date_start = input_variables[0]; date_stop = input_variables[1]; dt = input_variables[2]; nb_steps = input_variables[3]; nb_satellites = input_variables[4];
# output_path_propagator = input_variables[6];  output_file_propagator = input_variables[7]; 

# nb_time_steps_coverage = (int)((nb_steps-1)*dt/time_step_interpolation + 1)
# nb_ground_stations = 1 # !!!!!!!!!!!!!!!

# !!!!!!!!!!!!!!!!! DELETE
nb_satellites = 1
nb_ground_stations = 1
nb_time_steps_coverage = 338031 - 13 + 1
date_start = datetime.strptime("2016-06-02T03:41:54", "%Y-%m-%dT%H:%M:%S")
date_stop = datetime.strptime("2016-06-06T01:35:32", "%Y-%m-%dT%H:%M:%S")
time_start_xaxis = '2016-06-02T03:41:54'
time_stop_xaxis = '2016-06-06T01:35:32'
# !!!!!!!!!!!!!!!!! END OF DELETE
elevation_wtr_to_ground_station_in_sc_refce = np.zeros([nb_satellites, nb_ground_stations, nb_time_steps_coverage])
azimuth_wtr_to_ground_station_in_sc_refce = np.zeros([nb_satellites, nb_ground_stations, nb_time_steps_coverage])
range_wtr_to_ground_station = np.zeros([nb_satellites, nb_ground_stations, nb_time_steps_coverage])
elevation_wtr_to_ground_station_in_ground_station_refce = np.zeros([nb_satellites, nb_ground_stations, nb_time_steps_coverage])
azimuth_wtr_to_ground_station_in_ground_station_refce = np.zeros([nb_satellites, nb_ground_stations, nb_time_steps_coverage])
time_coverage_yes = [] # list of times when a satellite is in sight of a ground station
time_coverage = [] # just the time of the coverage output file (first  column of the file)
for isat in range(nb_satellites):
    time_coverage_yes_sat = []
    # !!!!!!!!!!!!! UNCOMMENT
#    file_coverage  = open(output_path_propagator[0] +'coverage_ground_station_' + output_file_propagator[0])
    # !!!!!!!!!!!!! END OF UNCOMMENT
    # !!!!!!!!!!!!!! DELETE
    file_coverage  = open('/home/cbv/PropSim/output/run_test_gps_coverage_v2/constellation_GPS/coverage_ground_station_GPS_BIIR-13.txt')  
    # !!!!!!!!!!!!!! END OF DELETE
    read_file_coverage = file_coverage.readlines()
    n_header = 0
    while ( read_file_coverage[n_header].split()[0] != "#START" ):
        n_header = n_header + 1
    n_header = n_header + 1
    nb_time_steps_coverage_temp = len(read_file_coverage) - n_header

    # QUICK CHECK THAT THERE IS NO MISSING STEP IN THE INTERPOLATION
    if nb_time_steps_coverage_temp != nb_time_steps_coverage:
        print "There might be a missing step in the interpolation for satellite "+str(isat) +".\n"
        raise Exception

    # READ THE COVERAGE OUTPUT FILE
    for istation in range(nb_ground_stations):
        time_coverage_yes_sat_station = []
        index_yes = -2
        for istep in range(nb_time_steps_coverage):
            if ( isat == 0 ) & (istation == 0): # the time is the same for all satellites and ground stations
                time_coverage.append(read_file_coverage[istep+n_header].split()[0])
            elevation_wtr_to_ground_station_in_sc_refce[isat, istation, istep] = np.float(read_file_coverage[istep+n_header].split()[4])
            azimuth_wtr_to_ground_station_in_sc_refce[isat, istation, istep] = np.float(read_file_coverage[istep+n_header].split()[5])
            elevation_wtr_to_ground_station_in_ground_station_refce[isat, istation, istep] = np.float(read_file_coverage[istep+n_header].split()[6])
            azimuth_wtr_to_ground_station_in_ground_station_refce[isat, istation, istep] = np.float(read_file_coverage[istep+n_header].split()[7])
            range_wtr_to_ground_station[isat, istation, istep] = np.float(read_file_coverage[istep+n_header].split()[8])
            if np.float(elevation_wtr_to_ground_station_in_sc_refce[isat, istation, istep]) > -999: 
                if (istep > index_yes + 1):# the coverage corresponds to a new fly over
                    time_coverage_yes_sat_station_flyover = []
                    time_coverage_yes_sat_station_flyover.append(read_file_coverage[istep+n_header].split()[0])
                index_yes = istep
                if istep != nb_time_steps_coverage - 1:
                    if np.float(read_file_coverage[istep+1+n_header].split()[4]) <= -999:# end of the fly over
                        time_coverage_yes_sat_station_flyover.append(read_file_coverage[istep+n_header].split()[0])
                        time_coverage_yes_sat_station.append(time_coverage_yes_sat_station_flyover)
                else:
                    time_coverage_yes_sat_station_flyover.append(read_file_coverage[istep+n_header].split()[0])
            
        time_coverage_yes_sat.append(time_coverage_yes_sat_station)
    time_coverage_yes.append(time_coverage_yes_sat)

    file_coverage.close()


# PLOTS
## ## RANGE
width_fig = 13
height_fig = 8
fontsize_plot = 14 
fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.05, right=0.99, top = 0.7,bottom = 0.2)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') 
ax1 = fig.add_subplot(gs[0, 0])
ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] 

isat = 0
istation = 0
delta_time = (date_stop - date_start).days*24*3600 + (date_stop - date_start).seconds
x_axis = np.arange(0, delta_time+1, time_step_interpolation )
ax1.plot(x_axis, range_wtr_to_ground_station[isat, istation, :],'.', linewidth = 1, color = 'k', label = '')

ax1.set_xlabel('', fontsize = fontsize_plot, weight = 'bold')
ax1.set_ylabel('Range (km)', fontsize = fontsize_plot, weight = 'bold', labelpad=0.0001)
ax1.set_title('Range', weight = 'bold', fontsize = 20,  y = 1.02)
gs.update(left=0.1, right=0.99, top = 0.9,bottom = 0.07)
ax1.legend()
ax1.set_ylim([0, max(range_wtr_to_ground_station[0,0])])

hour_time_step_xticks = 6
second_time_step_xticks = hour_time_step_xticks * 3600
xticks = np.arange(0, delta_time, second_time_step_xticks )
date_list_str = []
date_list = [date_start + timedelta(seconds=x) for x in xticks]
for i in range(len(xticks)):
    date_list_str.append( str(date_list[i])[5:11] +'\n'+str(date_list[i])[11:16] )
ax1.xaxis.set_ticks(xticks)
ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')

#1time_start_xaxis = '2014-07-10T12:44:37'
delta_time_start_xaxis = datetime.strptime(time_start_xaxis, "%Y-%m-%dT%H:%M:%S") - date_start
delta_time_start_xaxis = delta_time_start_xaxis.days*24*3600 + delta_time_start_xaxis.seconds
#time_stop_xaxis = '2014-07-10T19:33:00'
delta_time_stop_xaxis = datetime.strptime(time_stop_xaxis, "%Y-%m-%dT%H:%M:%S") - date_start
delta_time_stop_xaxis = delta_time_stop_xaxis.days*24*3600 + delta_time_stop_xaxis.seconds
ax1.set_xlim([delta_time_start_xaxis, delta_time_stop_xaxis])

raise Exception

## ## ELEVATION AND AZIMUTH ANGLES IN THE SATELLITE REFERENCE FRAME
width_fig = 13
height_fig = 8
fontsize_plot = 14 
fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.05, right=0.99, top = 0.7,bottom = 0.2)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') 
ax1 = fig.add_subplot(gs[0, 0])
ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] 

isat = 0
istation = 0
delta_time = (date_stop - date_start).days*24*3600 + (date_stop - date_start).seconds
x_axis = np.arange(0, delta_time+1, time_step_interpolation )
ax1.plot(x_axis, elevation_wtr_to_ground_station_in_sc_refce[isat, istation, :],'.', linewidth = 1, color = 'r', label = 'Elevation')
ax1.plot(x_axis, azimuth_wtr_to_ground_station_in_sc_refce[isat, istation, :],'.', linewidth = 1, color = 'b', label = 'Azimuth')

ax1.set_xlabel('Real Time', fontsize = fontsize_plot, weight = 'bold')
ax1.set_ylabel('Angle' + u' (\N{DEGREE SIGN})', fontsize = fontsize_plot, weight = 'bold', labelpad=0.0001)
ax1.set_title('Elevation and azimuth angles in the satellite reference frame', weight = 'bold', fontsize = 20,  y = 1.02)
gs.update(left=0.1, right=0.99, top = 0.9,bottom = 0.07)
ax1.legend()
hour_time_step_xticks = 1
second_time_step_xticks = hour_time_step_xticks * 3600
xticks = np.arange(0, delta_time, second_time_step_xticks )
date_list_str = []
date_list = [date_start + timedelta(seconds=x) for x in xticks]
for i in range(len(xticks)):
    date_list_str.append( str(date_list[i])[11:16] )
ax1.xaxis.set_ticks(xticks)
ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')

ax1.set_ylim([0,255])
#time_start_xaxis = '2016-06-02T03:41:54'
delta_time_start_xaxis = datetime.strptime(time_start_xaxis, "%Y-%m-%dT%H:%M:%S") - date_start
delta_time_start_xaxis = delta_time_start_xaxis.days*24*3600 + delta_time_start_xaxis.seconds
#time_stop_xaxis = '2016-06-06T02:00:04'
delta_time_stop_xaxis = datetime.strptime(time_stop_xaxis, "%Y-%m-%dT%H:%M:%S") - date_start
delta_time_stop_xaxis = delta_time_stop_xaxis.days*24*3600 + delta_time_stop_xaxis.seconds
ax1.set_xlim([delta_time_start_xaxis, delta_time_stop_xaxis])
raise Exception


## ## ELEVATION AND AZIMUTH ANGLES IN THE GOUND STATION REFERENCE FRAME
width_fig = 13
height_fig = 8
fontsize_plot = 14 
fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.05, right=0.99, top = 0.7,bottom = 0.2)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') 
ax1 = fig.add_subplot(gs[0, 0])
ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] 

isat = 0
istation = 0
delta_time = (date_stop - date_start).days*24*3600 + (date_stop - date_start).seconds
x_axis = np.arange(0, delta_time+1, time_step_interpolation )
ax1.plot(x_axis, elevation_wtr_to_ground_station_in_ground_station_refce[isat, istation, :],'.', linewidth = 1, color = 'r', label = 'Elevation')
ax1.plot(x_axis, azimuth_wtr_to_ground_station_in_ground_station_refce[isat, istation, :],'.', linewidth = 1, color = 'b', label = 'Azimuth')

ax1.set_xlabel('Real Time', fontsize = fontsize_plot, weight = 'bold')
ax1.set_ylabel('Angle' + u' (\N{DEGREE SIGN})', fontsize = fontsize_plot, weight = 'bold', labelpad=0.0001)
ax1.set_title('Elevation and azimuth angles in the ground station reference frame', weight = 'bold', fontsize = 20,  y = 1.02)
gs.update(left=0.1, right=0.99, top = 0.9,bottom = 0.07)
ax1.legend()
hour_time_step_xticks = 1
second_time_step_xticks = hour_time_step_xticks * 3600
xticks = np.arange(0, delta_time, second_time_step_xticks )
date_list_str = []
date_list = [date_start + timedelta(seconds=x) for x in xticks]
for i in range(len(xticks)):
    date_list_str.append( str(date_list[i])[11:16] )
ax1.xaxis.set_ticks(xticks)
ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')

ax1.set_ylim([0,255])
#time_start_xaxis = '2016-06-02T03:41:54'
delta_time_start_xaxis = datetime.strptime(time_start_xaxis, "%Y-%m-%dT%H:%M:%S") - date_start
delta_time_start_xaxis = delta_time_start_xaxis.days*24*3600 + delta_time_start_xaxis.seconds
#time_stop_xaxis = '2016-06-06T02:00:04'
delta_time_stop_xaxis = datetime.strptime(time_stop_xaxis, "%Y-%m-%dT%H:%M:%S") - date_start
delta_time_stop_xaxis = delta_time_stop_xaxis.days*24*3600 + delta_time_stop_xaxis.seconds
ax1.set_xlim([delta_time_start_xaxis, delta_time_stop_xaxis])
raise Exception


