# This script compares the orbit average density between different simulations and plot the integrated difference of the orbit average density as a function of the quantile of the F10.7/Ap distributions (used in our papers on collision avoidance). This script has been used for the GNC confenrece paper about how density uncertainties affect the Pc differently depending on the the encounter geometry. Dont'ned to set up anything, just run it.

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

oe_analysis = 'ecc' # set to 'inc' or 'ecc'

if oe_analysis == 'inc':
    root_list_run = ['alt350-350_inc90-0_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_bc0.02_f107100_ap12_quantile.txt',
                     'alt350-350_inc90-10_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_bc0.02_f107100_ap12_quantile.txt',
                     'alt350-350_inc90-20_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_bc0.02_f107100_ap12_quantile.txt',
                     'alt350-350_inc90-30_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_bc0.02_f107100_ap12_quantile.txt',
                     'alt350-350_inc90-45_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_bc0.02_f107100_ap12_quantile.txt',
                     'alt350-350_inc90-60_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_bc0.02_f107100_ap12_quantile.txt',
                     'alt350-350_inc90-70_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_bc0.02_f107100_ap12_quantile.txt',
                     'alt350-350_inc90-80_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_bc0.02_f107100_ap12_quantile.txt']
     # root list run = list run without the quqntil number in name. the diffce between two root run has to be the encounter geometry

elif oe_analysis == 'ecc':
    root_list_run = ['alt350-350_inc90-45_arg_per180-180_raan0-0_true_ano180-180_ecc0-0_bc0.02_f107100_ap12_quantile.txt',
                     'alt350-350_inc90-45_arg_per180-180_raan0-0_true_ano180-180_ecc0-1e-06_bc0.02_f107100_ap12_quantile.txt',
                     'alt350-350_inc90-45_arg_per180-180_raan0-0_true_ano180-180_ecc0-1e-05_bc0.02_f107100_ap12_quantile.txt',
                    'alt350-350_inc90-45_arg_per180-180_raan0-0_true_ano180-180_ecc0-0.0001_bc0.02_f107100_ap12_quantile.txt',
                     'alt350-350_inc90-45_arg_per180-180_raan0-0_true_ano180-180_ecc0-0.0005_bc0.02_f107100_ap12_quantile.txt',
                     'alt350-350_inc90-45_arg_per180-180_raan0-0_true_ano180-180_ecc0-0.001_bc0.02_f107100_ap12_quantile.txt',
                     'alt350-350_inc90-45_arg_per180-180_raan0-0_true_ano180-180_ecc0-0.0025_bc0.02_f107100_ap12_quantile.txt',
                     'alt350-350_inc90-45_arg_per180-180_raan0-0_true_ano180-180_ecc0-0.005_bc0.02_f107100_ap12_quantile.txt']
else:
    print "***! oe_analysis needs to be set to 'inc' or 'ecc'. The program will stop. !***"; raise Exception;
nb_run_main = len(root_list_run)
inc = np.zeros([nb_run_main])
angle_coll = np.zeros([nb_run_main])
ecc = np.zeros([nb_run_main])
list_run = []
for irun_main in range(nb_run_main):
    run_name_no_quan = root_list_run[irun_main]
    inc[irun_main] = np.float( root_list_run[irun_main].split('inc')[1].split('_')[0].split('-')[1] ) 
    angle_coll[irun_main] = np.abs(inc[irun_main] - np.float( root_list_run[irun_main].split('inc')[1].split('_')[0].split('-')[0] ))
    ecc[irun_main] = np.float( '-'.join(root_list_run[irun_main].split('ecc')[1].split('_')[0].split('-')[1:])  )
    for irun in range(1,10): 
        run_name = run_name_no_quan.split('quantile')[0] + 'quantile'+ str(irun) + run_name_no_quan.split('quantile')[1]
        list_run.append(run_name)

nb_run = len(list_run)
density = []
density_diff = [] # difference of density between the 2 sc
density_diff_integ = [] # integral over time of density_diff
nb_orbit = []
irun = 0
for irun_main in range(nb_run_main):
    irun = irun_main * 9
    irun_count = 0
    density_irun = []
    density_diff_irun = []
    density_diff_integ_irun = []
    nb_orbit_irun = []
    while irun_count< 9:
        print irun + irun_count
        
        run_name = list_run[irun + irun_count]
        input_filename =     run_name
        var_in, var_in_order = read_input_file(input_filename)
        output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
        output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')];
        nb_sc = var_in[find_in_read_input_order_variables(var_in_order, 'nb_sc')];
        nb_steps = var_in[find_in_read_input_order_variables(var_in_order, 'nb_steps')];
        dt = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')];

        nb_seconds_in_simu = ( nb_steps - 1 ) * dt
        density_isc = []
        nb_orbit_isc = []
        for isc in range(nb_sc):
            var_to_read = ["latitude", "density"]
            var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
            if isc == 0:
                density_here = np.zeros([nb_sc, nb_steps]) # all output files of one simulation have the same number of steps
                latitude = np.zeros([nb_sc, nb_steps]) # all output files of one simulation have the same number of steps
                date = var_out[find_in_read_input_order_variables(var_out_order, 'date')]

            density_here[isc, :nb_steps] = var_out[find_in_read_input_order_variables(var_out_order, 'density')]
            latitude[isc, :nb_steps] = var_out[find_in_read_input_order_variables(var_out_order, 'latitude')]
        
            density_orbit_averaged, time_averaged, index_time_averaged = orbit_average(density_here[isc, :nb_steps], latitude[isc, :nb_steps], date )
            if ((isc == 0) & (irun + irun_count == 0)): # assume the same times for the orbit for all runs and all sc. This is true within a few min at most so is an ok approximation
                x_axis_average = []
                nb_orbit_for_this_sc = len(time_averaged)
                date_average_start_orbit_list = np.array(time_averaged)[:,0]  # take the date at the start of the bin
                for iorbit in range(nb_orbit_for_this_sc):
                    date_average_start_orbit = date_average_start_orbit_list[iorbit]
                    date_average_start_orbit = datetime.strptime( date_average_start_orbit, "%Y/%m/%d %H:%M:%S.%f" )
                    nb_seconds_between_start_orbit_and_date_start = ( date_average_start_orbit - datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S.%f") ).total_seconds()
                    x_axis_average.append( nb_seconds_between_start_orbit_and_date_start )
                x_axis_average = np.array(x_axis_average)
            density_isc.append(density_orbit_averaged)
            nb_orbit_isc.append(len(density_orbit_averaged))
        density_diff_irun.append(np.array(density_isc[1]) - np.array(density_isc[0]))
        delta_x = np.roll(x_axis_average,-1)[:-1] - x_axis_average[:-1]
        density_diff_irun_arr = np.array(density_diff_irun[-1])
        delta_y = np.roll(density_diff_irun_arr,-1)[:-1] - density_diff_irun_arr[:-1]
        density_diff_integ_irun.append(np.sum(delta_x * delta_y / 2. )) # integral = delta_x * delta_y / 2 (triangle methods)
        irun_count = irun_count + 1
        density_irun.append(density_isc)
        nb_orbit_irun.append(nb_orbit_isc)
    density.append(density_irun)
    nb_orbit.append(nb_orbit_irun)
    density_diff.append(density_diff_irun)
    density_diff_integ.append(density_diff_integ_irun)

density = np.array(density) # index 0 is main run, 1 is quantile run, 2 is which sc, 3 is which orbit
density_diff = np.array(density_diff) # index 0 is main run, 1 is quantile run, 2 is which orbit
density_diff_integ = np.array(density_diff_integ) # index 0 is main run, 1 is quantile run
nb_orbit = np.array(nb_orbit)
nb_orbit_min = np.min(nb_orbit)
if np.max(nb_orbit) != nb_orbit_min:
    print "***! The number of orbits is not the same between each plot/sc. The script still runs but it's weird. !***"


# PLOTS
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 28
ratio_fig_size = 8./3
### For plots, generate disctinct colors
NCURVES = nb_run
np.random.seed(101)
curves = [np.random.random(20) for i in range(NCURVES)]
values = range(NCURVES)
jet = cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

## for a given main run, plot the desnsity for all quantiles for both sc
irun_main = 0
fig_title = ''
y_label = 'Density (kg/m$^3$)'
x_label = 'Real time'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 2)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

irun = irun_main * 9
for irun_count in range(9):
    print irun + irun_count
    run_name = list_run[irun + irun_count]
    colorVal = scalarMap.to_rgba(irun_count*nb_run_main)
    ax.plot(x_axis_average, density[irun_main, irun_count, 0,:], linewidth = 2, color = colorVal, label =  str(irun_count+1))
    ax.plot(x_axis_average, density[irun_main, irun_count, 1,:], linewidth = 2, color = colorVal, linestyle = 'dashed')

    if irun_count == 0:
        # x axis label is in real time
        ## all output files of one simulation have the same number of steps, and start at the same date
        nb_ticks_xlabel = 5
        dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
        xticks = np.arange(0, nb_seconds_in_simu+1, dt_xlabel)
        date_ref = datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S.%f")
        date_list_str = []
        date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
        for i in range(len(xticks)):
            if dt_xlabel >= 3*24*3600:
                date_list_str.append( str(date_list[i])[5:10] )
            else:
                date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
        ax.xaxis.set_ticks(xticks)
        ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
        ax.margins(0,0); ax.set_xlim([min(xticks), max(xticks)])


fig_title = ''
y_label = 'Density difference (kg/m$^3$)'
x_label = 'Real time'
ax = fig.add_subplot(gs[0, 1])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

irun = irun_main * 9
for irun_count in range(9):
    print irun + irun_count
    run_name = list_run[irun + irun_count]
    colorVal = scalarMap.to_rgba(irun_count*nb_run_main)
    ax.plot(x_axis_average, density_diff[irun_main, irun_count, :], linewidth = 2, color = colorVal, label =  str(irun_count+1))

    if irun_count == 0:
        # x axis label is in real time
        ## all output files of one simulation have the same number of steps, and start at the same date
        nb_ticks_xlabel = 5
        dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
        xticks = np.arange(0, nb_seconds_in_simu+1, dt_xlabel)
        date_ref = datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S.%f")
        date_list_str = []
        date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
        for i in range(len(xticks)):
            if dt_xlabel >= 3*24*3600:
                date_list_str.append( str(date_list[i])[5:10] )
            else:
                date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
        ax.xaxis.set_ticks(xticks)
        ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
        ax.margins(0,0); ax.set_xlim([min(xticks), max(xticks)])



legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Quantile", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))


# legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Quantile", fontsize = fontsize_plot)
# legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = 'density_both_sc_angle_coll_' + str(angle_coll[irun_main]) + '_ecc_' + str(ecc[irun_main])  
fig_save_name =  fig_save_name + '.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

raise Exception



## for a given main run, plot the desnsity difference between the 2 sc for all quantiles 
irun_main = 0
fig_title = ''
y_label = 'Density difference (kg/m$^3$)'
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

irun = irun_main * 9
for irun_count in range(9):
    print irun + irun_count
    run_name = list_run[irun + irun_count]
    colorVal = scalarMap.to_rgba(irun_count*nb_run_main)
    ax.plot(x_axis_average, density_diff[irun_main, irun_count, :], linewidth = 2, color = colorVal, label =  str(irun_count+1))

    if irun_count == 0:
        # x axis label is in real time
        ## all output files of one simulation have the same number of steps, and start at the same date
        nb_ticks_xlabel = 10
        dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
        xticks = np.arange(0, nb_seconds_in_simu+1, dt_xlabel)
        date_ref = datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S.%f")
        date_list_str = []
        date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
        for i in range(len(xticks)):
            if dt_xlabel >= 3*24*3600:
                date_list_str.append( str(date_list[i])[5:10] )
            else:
                date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
        ax.xaxis.set_ticks(xticks)
        ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
        ax.margins(0,0); ax.set_xlim([min(xticks), max(xticks)])



legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Quantile", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = 'density_difference_angle_coll_' + str(angle_coll[irun_main]) + '_ecc_' + str(ecc[irun_main])  
fig_save_name =  fig_save_name + '.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  






## for all main runs, plot the integrated desnity difference as a function of the quantile
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 28
ratio_fig_size = 8./3
### For plots, generate disctinct colors

fig_title = ''
y_label = 'Integrated density difference (kg.s/m$^3$)'
x_label = 'Quantile'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 2)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01, wspace  = 0.16)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

color_arr = ['b', 'r','cornflowerblue','g', 'm', 'gold', 'cyan', 'fuchsia', 'lawngreen', 'darkgray', 'green', 'chocolate']
quantile_arr = np.arange(1,10)
for irun_main in range(nb_run_main):
    run_name = list_run[irun + irun_count]
    if oe_analysis == 'inc':
        ax.plot(quantile_arr, density_diff_integ[irun_main, :], linewidth = 2, color = color_arr[irun_main], label = str(angle_coll[irun_main]) + u'\N{DEGREE SIGN}' )
    if oe_analysis == 'ecc':
        ax.plot(quantile_arr, density_diff_integ[irun_main, :], linewidth = 2, color = color_arr[irun_main], label = str(ecc[irun_main]) )
    #ax.margins(0,0)
ax.set_xticks(quantile_arr)

if oe_analysis == 'inc':
    legend = ax.legend(loc='center right', bbox_to_anchor=(-0.17, 0.5), numpoints = 1,  title="Abs. inclination\ndifference", fontsize = fontsize_plot)
if oe_analysis == 'ecc':
    legend = ax.legend(loc='center right', bbox_to_anchor=(-0.17, 0.5), numpoints = 1,  title="Eccentricity\nof spacecraft 2", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))

## show slope of integrated density difference as a function of encounter geometry (collision angle, eccentricity). by slope, actuallymean last value of integrrated dens diff minus first value of it. This is ok because the integrated density difference is pretty linear (see plot integrated_density_difference.pdf made right above)
slope_integrated_dens_diff = np.zeros([nb_run_main])
for irun_main in range(nb_run_main):
    slope_integrated_dens_diff[irun_main] = density_diff_integ[irun_main, -1] - density_diff_integ[irun_main, 0]
fig_title = '' 
y_label = 'Max - min integrated density diff. (kg.s/m$^3$)'
if oe_analysis == 'inc':
    x_label = 'Absolute difference in inclination ' + u'(\N{DEGREE SIGN})'
if oe_analysis == 'ecc':
    x_label = 'Eccentricity of spacecraft 2'
ax = fig.add_subplot(gs[0, 1])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

color_arr = ['b','r','g', 'm', 'gold', 'cyan', 'cornflowerblue']
if oe_analysis == 'inc':
    ax.plot(angle_coll, slope_integrated_dens_diff, linewidth = 2)
    ax.scatter(angle_coll, slope_integrated_dens_diff, linewidth = 2)
    ax.set_xticks(angle_coll)
if oe_analysis == 'ecc':
    ax.plot(ecc, slope_integrated_dens_diff, linewidth = 2)
    ax.scatter(ecc, slope_integrated_dens_diff, linewidth = 2)
    ax.set_xticks(ecc)
    #ax.margins(0,0)

if oe_analysis == 'inc':
    fig_save_name = 'integrated_density_difference_difference_inclination'
if oe_analysis == 'ecc':
    fig_save_name = 'integrated_density_difference_eccentricity'
fig_save_name =  fig_save_name + '.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  




## show slope of integrated density difference as a function of encounter geometry (collision angle, eccentricity). by slope, actuallymean last value of integrrated dens diff minus first value of it. This is ok because the integrated density difference is pretty linear (see plot integrated_density_difference.pdf made right above)
slope_integrated_dens_diff = np.zeros([nb_run_main])
for irun_main in range(nb_run_main):
    slope_integrated_dens_diff[irun_main] = density_diff_integ[irun_main, -1] - density_diff_integ[irun_main, 0]
fig_title = '' 
y_label = 'Slope of integrated density difference (kg.s/m$^3$)'
if oe_analysis == 'inc':
    x_label = 'Absolute difference in inclination ' + u'(\N{DEGREE SIGN})'
if oe_analysis == 'ecc':
    x_label = 'Eccentricity of spacecraft 2'
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

color_arr = ['b','r','g', 'm', 'gold', 'cyan', 'cornflowerblue']
if oe_analysis == 'inc':
    ax.plot(angle_coll, slope_integrated_dens_diff, linewidth = 2)
    ax.scatter(angle_coll, slope_integrated_dens_diff, linewidth = 2)
if oe_analysis == 'ecc':
    ax.plot(ecc, slope_integrated_dens_diff, linewidth = 2)
    ax.scatter(ecc, slope_integrated_dens_diff, linewidth = 2)
    #ax.margins(0,0)

if oe_analysis == 'inc':
    fig_save_name = 'slope_integrated_density_difference_inclination'
if oe_analysis == 'ecc':
    fig_save_name = 'slope_integrated_density_difference_eccentricity'
fig_save_name =  fig_save_name + '.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
