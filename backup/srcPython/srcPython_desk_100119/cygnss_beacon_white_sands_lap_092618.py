# This script was made for the beacon test campaing in White Sands, NM, in November 2018. The results were discussesd with Darren McKague.


import sys
#sys.path.append("/Users/cbv/Library/Enthought/Canopy_64bit/User/lib/python2.7/site-packages/")
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/srcPython")
#sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/kalman/spock_development_new_structure_kalman_dev/srcPython")
import matplotlib
matplotlib.use("Agg") # without this line, when running this script from a cronjob we get an error "Unable to access the X Display, is $DISPLAY set properly?"

import matplotlib.gridspec as gridspec
import numpy as np
from struct import *
from matplotlib import pyplot as plt
from cygnss_read_spock_spec_bin import *
from mpl_toolkits.basemap import Basemap, shiftgrid
from datetime import datetime, timedelta
from collections import *
import os
from read_input_file import *
from read_output_file import *
#from cygnss_read_spock_spec import *




input_filename = 'beacon.txt'

var_in, var_in_order = read_input_file(input_filename)
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
dt = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')]; 
nb_steps = var_in[find_in_read_input_order_variables(var_in_order, 'nb_steps')]; 
nb_sc = var_in[find_in_read_input_order_variables(var_in_order, 'nb_sc')]; 
date_start = var_in[find_in_read_input_order_variables(var_in_order, 'date_start')]
date_stop = var_in[find_in_read_input_order_variables(var_in_order, 'date_stop')]

contact_time = []
for isc in range(nb_sc):
    contact_time_sub = []
    gs_filename = output_file_path_list[isc] + 'coverage/report_all_by_' + output_file_name_list[isc]
    gs_file = open(gs_filename)
    read_gs_file = gs_file.readlines()
    nhead = 2
    nb_ct = len(read_gs_file) - nhead
    for ict in range(nb_ct):
        start = datetime.strptime( read_gs_file[ict + nhead].split()[1] + 'T' + read_gs_file[ict + nhead].split()[2], "%Y-%m-%dT%H:%M:%S" )
        end = datetime.strptime( read_gs_file[ict + nhead].split()[4] + 'T' + read_gs_file[ict + nhead].split()[5], "%Y-%m-%dT%H:%M:%S" )
        duration = np.float(read_gs_file[ict + nhead].split()[9].replace("(","") )
        contact_time_sub.append([start, end, duration ])
    gs_file.close()

    contact_time.append(contact_time_sub)

# Determine how long it takes for the beacon station to see all 8 sc
# If 2 or more sc are in sight, only consider the one that's in contact for the longest time
# This is because the beacon won't be able to be steered from on sc to another so quickly
# Actually, if for instance there's 2 sc, at the first overpass, consider the one that is in conact the longest
#for isc in range(nb_sc):


# FIGURES
satColors = ['black', 'blue', 'red', 'mediumorchid', 'dodgerblue', 'magenta', 'darkgreen', 'limegreen'] #['lime', 'blue', 'indigo', 'purple', 'dodgerblue', 'steelblue', 'seagreen', 'limegreen']
nb_sc = 8 # !!!!!!!!!
label_arr = ['FM05', 'FM04', 'FM02', 'FM01', 'FM08', 'FM06', 'FM07', 'FM03']
label_arr_conversion = [3, 2, 7, 1, 0, 5, 6, 4]

# Plot the coverage VS time for all sc
width_fig = 15.
height_fig = width_fig * 3 /4
fontsize_plot = 20 # 9
color_arr = ['k','b','r','g','m', 'y']

iday = 0
nb_day = (int)( np.ceil((date_stop - date_start).total_seconds() / 3600./24) )

date_day_arr = np.array([date_start + timedelta(days=i) for i in np.arange(0, nb_day +1, 1)])
for iday in range(nb_day):
    date_start_day = date_day_arr[iday]
    date_stop_day = date_day_arr[iday+1]
    max_xaxis = (date_stop_day - date_start_day).total_seconds()
    fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
    gs = gridspec.GridSpec(1, 1)
    gs.update(left=0.12, right=0.97, top = 0.90,bottom = 0.06)
    fig.suptitle('', fontsize = 22)
    plt.rc('font', weight='normal') ## make the labels of the ticks in normal
    ax1 = fig.add_subplot(gs[0, 0])
    #ax1 = fig.add_subplot(111)
    ax1.set_title('', weight = 'normal', fontsize = 20,  y = 1.008) 
    ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
    #ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    path_folder_plot = '.'
    date_ref = date_start_day
    yticks = np.arange(0, 1, 1./nb_sc)
    ytick_labels = []

    

    for isc_temp in range(nb_sc):
        isc = label_arr_conversion[isc_temp]
        for icoverage in range( len(contact_time[isc]) ):
            if ((contact_time[isc][icoverage][0] >= date_ref) & (contact_time[isc][icoverage][0] <= date_stop_day)):
                print label_arr[isc], contact_time[isc][icoverage][0], contact_time[isc][icoverage][1]
                start_coverage_since_ref = ( contact_time[isc][icoverage][0] - date_ref ).total_seconds()
                end_coverage_since_ref = ( contact_time[isc][icoverage][1] - date_ref ).total_seconds()
                #ax1.plot([start_coverage_since_ref, end_coverage_since_ref], [isc * 1./nb_sc, isc * 1./nb_sc], linewidth = 20, color = satColors[isc])
                ax1.fill_between([start_coverage_since_ref, end_coverage_since_ref], [isc_temp * 1./nb_sc-1./(6*nb_sc), isc_temp * 1./nb_sc-1./(6*nb_sc)],[ isc_temp * 1./nb_sc+1./(6*nb_sc),isc_temp * 1./nb_sc+1./(6*nb_sc)], color = satColors[isc])

        ytick_labels.append(label_arr[isc])

    ax1.yaxis.set_ticks(yticks)
    ax1.yaxis.set_ticklabels(ytick_labels, fontsize = fontsize_plot)#, rotation='vertical')
    isc_temp = 0
    for ytick, color in zip(ax1.get_yticklabels(), satColors):
        isc = label_arr_conversion[isc_temp]
        ytick.set_color(satColors[isc])
        isc_temp = isc_temp + 1


    hour_time_step_xticks = 3
    hour_time_step_xticks_converted_in_seconds = hour_time_step_xticks * 3600
    xticks = np.arange(0, max_xaxis+hour_time_step_xticks_converted_in_seconds, hour_time_step_xticks_converted_in_seconds)
    date_list_str = []
    date_list = [date_start_day + timedelta(seconds=x) for x in np.arange(0, max_xaxis+hour_time_step_xticks_converted_in_seconds, hour_time_step_xticks_converted_in_seconds)]
    for i in range(len(xticks)):
        date_list_str.append( str(date_list[i])[11:16] + '\n' + datetime.strftime(date_list[i], "%b %-d"))
    ax1.xaxis.set_ticks(xticks) 
    ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')

    ax1.set_xlim([16*3600,23*3600])#ax1.set_xlim([0,max_xaxis])
    ax1.set_ylim([-1./nb_sc,1])
    ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

    ax1.set_title('Contact times - day ' + str(iday + 1), weight = 'normal', fontsize = fontsize_plot , y = 1.00) 
    ax1.legend( fontsize = fontsize_plot, loc = 4)
    ax1.set_xlabel('Real time', fontsize = fontsize_plot, weight = 'normal')
    ax1.set_ylabel('', fontsize = fontsize_plot, weight = 'normal')

    fig_save_name = path_folder_plot + '/' + input_filename.replace(".txt","_") + 'iday' + str(iday+1) + ".pdf"
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    raise Exception



# label_arr = ['FM05', 'FM04', 'FM02', 'FM01', 'FM08', 'FM06', 'FM07', 'FM03']
# label_arr_conversion = [3, 2, 7, 1, 0, 5, 6, 4]

#         for isc_temp in range(nb_sc):
#             isc = label_arr_conversion[isc_temp]

