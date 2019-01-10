# this sctipt makes the f107 ap density and screenshot plots of paper4


import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/kalman/spock_development_new_structure_kalman_dev/srcPython")
import pickle
from orbit_average import *
import matplotlib.colors as colors
import matplotlib.cm as cmx
from cadre_read_last_tle import *
from get_prop_dir import *
import matplotlib.gridspec as gridspec
from read_input_file import *
from read_output_file import *
from matplotlib.colors import LogNorm
from norad_id_to_cygnss_name import *
import pickle
from eci_to_lvlh import *
import fileinput
import time
from datetime import datetime
import numpy as np
from matplotlib import pyplot as plt
import os
import sys
import subprocess
from get_name_mission import *
from find_in_read_input_order_variables import *

to_plot = 'density' #'f107', or 'ap'

if to_plot == 'f107':
    list_run = ['f107_-13_ap_0.txt', 
            'f107_-9_ap_0.txt',
            'f107_-4_ap_0.txt',
            'f107_0_ap_0.txt',
            'f107_4_ap_0.txt',
            'f107_9_ap_0.txt',
            'f107_13_ap_0.txt'
            ]

if to_plot == 'ap':
    list_run =     ['f107_0_ap_-13.txt',
                    'f107_0_ap_-9.txt',
                    'f107_0_ap_-4.txt',
                    'f107_0_ap_0.txt',
                    'f107_0_ap_4.txt',
                    'f107_0_ap_9.txt',
                    'f107_0_ap_13.txt']


if to_plot == 'density':
    list_run =     ['f107_-13_ap_-13.txt',
                   'f107_-13_ap_-9.txt',
                   'f107_-13_ap_-4.txt',
                   'f107_-13_ap_0.txt',
                   'f107_-13_ap_4.txt',
                   'f107_-13_ap_9.txt',
                   'f107_-13_ap_13.txt',
                   'f107_-9_ap_-13.txt',
                   'f107_-9_ap_-9.txt',
                   'f107_-9_ap_-4.txt',
                   'f107_-9_ap_0.txt',
                   'f107_-9_ap_4.txt',
                   'f107_-9_ap_9.txt',
                   'f107_-9_ap_13.txt',
                   'f107_-4_ap_-13.txt',
                   'f107_-4_ap_-9.txt',
                   'f107_-4_ap_-4.txt',
                   'f107_-4_ap_0.txt',
                   'f107_-4_ap_4.txt',
                   'f107_-4_ap_9.txt',
                   'f107_-4_ap_13.txt',
                   'f107_0_ap_-13.txt',
                   'f107_0_ap_-9.txt',
                   'f107_0_ap_-4.txt',
                   'f107_0_ap_0.txt',
                   'f107_0_ap_4.txt',
                   'f107_0_ap_9.txt',
                   'f107_0_ap_13.txt',
                   'f107_4_ap_-13.txt',
                   'f107_4_ap_-9.txt',
                   'f107_4_ap_-4.txt',
                   'f107_4_ap_0.txt',
                   'f107_4_ap_4.txt',
                   'f107_4_ap_9.txt',
                   'f107_4_ap_13.txt',
                   'f107_9_ap_-13.txt',
                   'f107_9_ap_-9.txt',
                   'f107_9_ap_-4.txt',
                   'f107_9_ap_0.txt',
                   'f107_9_ap_4.txt',
                   'f107_9_ap_9.txt',
                   'f107_9_ap_13.txt',
                   'f107_13_ap_-13.txt',
                   'f107_13_ap_-9.txt',
                   'f107_13_ap_-4.txt',
                   'f107_13_ap_0.txt',
                   'f107_13_ap_4.txt',
                   'f107_13_ap_9.txt',
                   'f107_13_ap_13.txt']




nb_run = len(list_run)
f107 = []
ap = []
density = []
for irun in range(nb_run):
    run_name = list_run[irun]
    input_filename =  run_name
    var_in, var_in_order = read_input_file(input_filename)
    output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
    output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')];
    nb_sc = var_in[find_in_read_input_order_variables(var_in_order, 'nb_sc')];
    nb_steps = var_in[find_in_read_input_order_variables(var_in_order, 'nb_steps')];
    dt = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')];

    nb_seconds_in_simu = ( nb_steps - 1 ) * dt
    density_isc = []
    nb_orbit_isc = []
    isc = 0
    var_to_read = ["latitude", "density", "ap", "f107"]
    var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )

    density_here = np.zeros([nb_sc, nb_steps]) # all output files of one simulation have the same number of steps
    latitude = np.zeros([nb_sc, nb_steps]) # all output files of one simulation have the same number of steps
    date = var_out[find_in_read_input_order_variables(var_out_order, 'date')]
    
    density_here[isc, :nb_steps] = var_out[find_in_read_input_order_variables(var_out_order, 'density')]
    latitude[isc, :nb_steps] = var_out[find_in_read_input_order_variables(var_out_order, 'latitude')]
    ap.append( var_out[find_in_read_input_order_variables(var_out_order, 'ap')] )
    f107.append( var_out[find_in_read_input_order_variables(var_out_order, 'f107')])
    if irun == 0:
        date_given_output = var_out[find_in_read_input_order_variables(var_out_order, 'date_given_output')]
        nb_steps_given_output = len(date_given_output)
        nb_seconds_since_start_for_given_output = np.zeros([nb_steps_given_output])
        for istep_given in range(nb_steps_given_output):
            date_given_output_datetime = datetime.strptime(date_given_output[istep_given],"%Y-%m-%dT%H:%M:%S.%f:") # 2016-11-26T00:00:10.000000:
            nb_seconds_since_start_for_given_output[istep_given] = ( date_given_output_datetime - datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S.%f") ).total_seconds() #'2016/11/26 00:00:00.000000'

    density_orbit_averaged, time_averaged, index_time_averaged = orbit_average(density_here[isc, :nb_steps], latitude[isc, :nb_steps], date )
    if (irun == 0): # assume the same times for the orbit for all runs. This is true within a few min at most so is an ok approximation
        x_axis_average = []
        nb_orbit_for_this_sc = len(time_averaged)
        date_average_start_orbit_list = np.array(time_averaged)[:,0]  # take the date at the start of the bin
        for iorbit in range(nb_orbit_for_this_sc):
            date_average_start_orbit = date_average_start_orbit_list[iorbit]
            date_average_start_orbit = datetime.strptime( date_average_start_orbit, "%Y/%m/%d %H:%M:%S.%f" )
            nb_seconds_between_start_orbit_and_date_start = ( date_average_start_orbit - datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S.%f") ).total_seconds()
            x_axis_average.append( nb_seconds_between_start_orbit_and_date_start )
        x_axis_average = np.array(x_axis_average)
    density.append(density_orbit_averaged)

f107 = np.array(f107).astype(np.float)
ap = np.array(ap).astype(np.float)
density = np.array(density)

# PLOTS
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3
### For plots, generate disctinct colors
NCURVES = nb_run
np.random.seed(101)
curves = [np.random.random(20) for i in range(NCURVES)]
values = range(NCURVES)
jet = cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)




if to_plot == 'f107':
    ## f107 vs time
    fig_title = ''
    y_label = 'F10.7'
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

    max_y = 0
    min_y = 100000
    for irun in range(nb_run):
        print irun 
        run_name = list_run[irun]
        colorVal = scalarMap.to_rgba(irun)
        ax.plot(nb_seconds_since_start_for_given_output, f107[irun, :], linewidth = 2, color = colorVal, label = run_name.split('f107_')[1].split('_')[0]) # run_name.split('ap_')[1].split('.')[0]

        if np.max(f107[irun, :]) > max_y:
            max_y = np.max(f107[irun, :])
        if np.min(f107[irun, :]) < min_y:
            min_y = np.min(f107[irun, :])

        if irun == 0:
            # x axis label is in real time
            ## all output files of one simulation have the same number of steps, and start at the same date
            nb_ticks_xlabel = 3
            dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
            xticks = np.arange(0, nb_seconds_in_simu+1, dt_xlabel)
            date_ref = datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S.%f")
            date_list_str = []
            date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
            for i in range(len(xticks)):
                if dt_xlabel >= 3*24*3600:
                    date_list_str.append( str(date_list[i])[5:10] )
                else:
                    date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str((int)(xticks[i]/3600./24)) + ')' )
            ax.xaxis.set_ticks(xticks)
            ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')

    ax.set_xlim([min(xticks), max(xticks)])
    ax.set_ylim([min_y, max_y+1])
    ax.margins(0,0);

    legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Bin value", fontsize = fontsize_plot)
    legend.get_title().set_fontsize(str(fontsize_plot))

    fig_save_name = 'f107' 
    fig_save_name =  fig_save_name + '.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


if to_plot == 'ap':
    ## ap vs time
    fig_title = ''
    y_label = 'Ap'
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

    max_y = 0
    min_y = 100000
    for irun in range(nb_run):
        print irun 
        run_name = list_run[irun]
        colorVal = scalarMap.to_rgba(irun)
        ax.plot(nb_seconds_since_start_for_given_output, ap[irun, :], linewidth = 2, color = colorVal, label =  run_name.split('ap_')[1].split('.')[0])

        if np.max(ap[irun, :]) > max_y:
            max_y = np.max(ap[irun, :])
        if np.min(ap[irun, :]) < min_y:
            min_y = np.min(ap[irun, :])

        if irun == 0:
            # x axis label is in real time
            ## all output files of one simulation have the same number of steps, and start at the same date
            nb_ticks_xlabel = 3
            dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
            xticks = np.arange(0, nb_seconds_in_simu+1, dt_xlabel)
            date_ref = datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S.%f")
            date_list_str = []
            date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
            for i in range(len(xticks)):
                if dt_xlabel >= 3*24*3600:
                    date_list_str.append( str(date_list[i])[5:10] )
                else:
                    date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str((int)(xticks[i]/3600./24)) + ')' )
            ax.xaxis.set_ticks(xticks)
            ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')

    ax.set_xlim([min(xticks), max(xticks)])
    ax.set_ylim([min_y, max_y+1])
    ax.margins(0,0);

    legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Bin value", fontsize = fontsize_plot)
    legend.get_title().set_fontsize(str(fontsize_plot))

    fig_save_name = 'ap' 
    fig_save_name =  fig_save_name + '.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  



if to_plot == 'density':
    ## density vs time
    fig_title = ''
    y_label = 'Density (kg/m$^3$)'
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

    max_y = 0
    min_y = 10000000000000000000
    for irun in range(nb_run):
        print irun 
        run_name = list_run[irun]
        colorVal = scalarMap.to_rgba(irun)
        ax.plot(x_axis_average, density[irun, :], linewidth = 2, color = 'cornflowerblue')#colorVal)#, label = run_name.split('density_')[1].split('_')[0]) # run_name.split('ap_')[1].split('.')[0]

        if np.max(density[irun, :]) > max_y:
            max_y = np.max(density[irun, :])
        if np.min(density[irun, :]) < min_y:
            min_y = np.min(density[irun, :])

        if irun == 0:
            # x axis label is in real time
            ## all output files of one simulation have the same number of steps, and start at the same date
            nb_ticks_xlabel = 3
            dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
            xticks = np.arange(0, nb_seconds_in_simu+1, dt_xlabel)
            date_ref = datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S.%f")
            date_list_str = []
            date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
            for i in range(len(xticks)):
                if dt_xlabel >= 3*24*3600:
                    date_list_str.append( str(date_list[i])[5:10] )
                else:
                    date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str((int)(xticks[i]/3600./24)) + ')' )
            ax.xaxis.set_ticks(xticks)
            ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')

    ax.set_xlim([min(xticks), max(xticks)])
    ax.set_ylim([min_y, max_y])
    ax.margins(0,0);

    #     legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Bin value", fontsize = fontsize_plot)
    #     legend.get_title().set_fontsize(str(fontsize_plot))

    fig_save_name = 'density' 
    fig_save_name =  fig_save_name + '.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
