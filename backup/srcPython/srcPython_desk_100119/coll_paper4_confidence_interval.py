# This script compares the Pc from 100 idnetical runs to compute an empiral confidence interval
# the values of Pc are read from a file genreated file_results_pc on Big by /raid4/cbv/collision/list_read_collision_file_new.py (although the SpOCK runs are made on Pleiades)


import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/srcPython")
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
from mpl_toolkits.mplot3d import Axes3D



filename_results_pc = 'results_f107_-13_ap_-13_all_run_confidence_interval.txt'
file_results_pc = open(filename_results_pc)
read_file_results_pc = file_results_pc.readlines()
nb_run = len(read_file_results_pc)

pc = [] 
for irun in range(nb_run):
    pc.append( np.float( read_file_results_pc[irun].split()[1] ) )
file_results_pc.close()

pc = np.array(pc)
index_sort_pc = np.argsort(pc)
pc_sort = pc[index_sort_pc] 

pc_median = np.median(pc)
pc_max = np.max(pc)
pc_min = np.min(pc)
relative_error_actual = ( (pc_max - pc_median) / pc_median + (pc_median - pc_min) / pc_median ) / 2.

# Compute the relative error that corresponds to a confidence interval of 95%

relative_error_target = 1
pc_interval_max = pc_median * (1 + relative_error_target)
pc_interval_min = pc_median * (1 - relative_error_target)
ci = len(np.where((pc >= pc_interval_min) & (pc <= pc_interval_max))[0]) * 100. / nb_run

while ci > 96.0:
    pc_interval_max = pc_median * (1 + relative_error_target)
    pc_interval_min = pc_median * (1 - relative_error_target)
    ci = np.ceil(len(np.where((pc >= pc_interval_min) & (pc <= pc_interval_max))[0]) * 100. / nb_run)
    relative_error_target = relative_error_target*0.9
    #print ci, relative_error_target * 100

print ci, relative_error_target*100., relative_error_actual * 100.

print '***! For the paper submitted to Space Weather, the relative error was 2.5% for a CI of 93% and 2.8% for a CI of 98%. Charles interpolated the error to 2.6% for a CI of 95%. So Charles "hacked" in the code to put these values of CI and the relative error right before the plot. ***!'
#!!!!!!! Comment lines below, unless you know what you're doing
ci = 95.
relative_error_target = 0.026
#!!!!!!! end of Comment lines below, unless you know what you're doing



# Plot the distribution
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20
ratio_fig_size = 4./3
width_fig = 25


## densityof probbiliy and CDF
fig_title = ''
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 2)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)


ax = fig.add_subplot(gs[0, 0])
y_label = 'Histogram (%)'
x_label = '$P_C$'

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

x_axis = pc_sort
y_axis =  np.zeros([nb_run]) + 1 #np.arange(1, nb_run +1 ,1)
nbins = 10
#bin_width = (np.max(pc) - np.min(pc))/nbins 
#bins = np.arange(np.min(pc), np.max(pc), bin_width)
hist = np.histogram(pc, bins = nbins)
hist_pc = hist[0] * 100. / nb_run
bins_array = hist[1][:-1]
width_bin = bins_array[1] - bins_array[0]
ax.bar(bins_array, hist_pc, width_bin, align = 'edge', linewidth = 2)

ax.plot([pc_median, pc_median], [0, np.max(hist_pc)], linewidth = 3, color = 'k', linestyle = 'dashed' )
ax.plot([pc_median * (1 + relative_error_target), pc_median * (1 + relative_error_target)], [0, np.max(hist_pc)], linewidth = 3, color = 'r', linestyle = 'dashed' )
ax.plot([pc_median * (1 - relative_error_target), pc_median * (1 - relative_error_target)], [0, np.max(hist_pc)], linewidth = 3, color = 'r', linestyle = 'dashed' )

ax.text(pc_median,  ( np.max(hist_pc) + np.min(hist_pc) ) / 2., '$\mathrm{P_{c, median}}$', rotation = 90, horizontalalignment = 'right', verticalalignment = 'center', fontsize = fontsize_plot, weight = 'bold', color = 'k') 
ax.text(pc_median * (1 - relative_error_target),  ( np.max(hist_pc) + np.min(hist_pc) ) / 2., '$\mathrm{P_{c, median} - error}$', rotation = 90, horizontalalignment = 'right', verticalalignment = 'center', fontsize = fontsize_plot, weight = 'bold', color = 'r') 
ax.text(pc_median * (1 + relative_error_target), ( np.max(hist_pc) + np.min(hist_pc) ) / 2., '$\mathrm{P_{c, median} + error}$', rotation = 90, horizontalalignment = 'right', verticalalignment = 'center', fontsize = fontsize_plot, weight = 'bold', color = 'r') 


ax.text(0.58, 0.98, 'N = ' + str(nb_run) + '\nCI = ' + format(ci, ".0f") + "%\nerror = " + format(relative_error_target*100, ".1f")  + "%", transform = ax.transAxes, horizontalalignment = 'left', verticalalignment = 'top', fontsize = fontsize_plot, weight = 'bold')
ax.margins(0,0)






filename_results_pc = 'results_f107_0_ap_0_all_run_confidence_interval.txt'
file_results_pc = open(filename_results_pc)
read_file_results_pc = file_results_pc.readlines()
nb_run = len(read_file_results_pc)

pc = [] 
for irun in range(nb_run):
    pc.append( np.float( read_file_results_pc[irun].split()[1] ) )
file_results_pc.close()

pc = np.array(pc)
index_sort_pc = np.argsort(pc)
pc_sort = pc[index_sort_pc] 

pc_median = np.median(pc)
pc_max = np.max(pc)
pc_min = np.min(pc)
relative_error_actual = ( (pc_max - pc_median) / pc_median + (pc_median - pc_min) / pc_median ) / 2.

# Compute the relative error that corresponds to a confidence interval of 95%

relative_error_target = 1
pc_interval_max = pc_median * (1 + relative_error_target)
pc_interval_min = pc_median * (1 - relative_error_target)
ci = len(np.where((pc >= pc_interval_min) & (pc <= pc_interval_max))[0]) * 100. / nb_run

while ci > 96.0:
    pc_interval_max = pc_median * (1 + relative_error_target)
    pc_interval_min = pc_median * (1 - relative_error_target)
    ci = np.ceil(len(np.where((pc >= pc_interval_min) & (pc <= pc_interval_max))[0]) * 100. / nb_run)
    relative_error_target = relative_error_target*0.9
    #print ci, relative_error_target * 100

print ci, relative_error_target*100., relative_error_actual * 100.

ax = fig.add_subplot(gs[0, 1])
y_label = 'Histogram (%)'
x_label = '$P_C$'

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

x_axis = pc_sort
y_axis =  np.zeros([nb_run]) + 1 #np.arange(1, nb_run +1 ,1)
nbins = 10
#bin_width = (np.max(pc) - np.min(pc))/nbins 
#bins = np.arange(np.min(pc), np.max(pc), bin_width)
hist = np.histogram(pc, bins = nbins)
hist_pc = hist[0] * 100. / nb_run
bins_array = hist[1][:-1]
width_bin = bins_array[1] - bins_array[0]
ax.bar(bins_array, hist_pc, width_bin, align = 'edge', linewidth = 2)

ax.plot([pc_median, pc_median], [0, np.max(hist_pc)], linewidth = 3, color = 'k', linestyle = 'dashed' )
ax.plot([pc_median * (1 + relative_error_target), pc_median * (1 + relative_error_target)], [0, np.max(hist_pc)], linewidth = 3, color = 'r', linestyle = 'dashed' )
ax.plot([pc_median * (1 - relative_error_target), pc_median * (1 - relative_error_target)], [0, np.max(hist_pc)], linewidth = 3, color = 'r', linestyle = 'dashed' )

ax.text(pc_median,  ( np.max(hist_pc) + np.min(hist_pc) ) / 2., '$\mathrm{P_{c, median}}$', rotation = 90, horizontalalignment = 'right', verticalalignment = 'center', fontsize = fontsize_plot, weight = 'bold', color = 'k') 
ax.text(pc_median * (1 - relative_error_target),  ( np.max(hist_pc) + np.min(hist_pc) ) / 2., '$\mathrm{P_{c, median} - error}$', rotation = 90, horizontalalignment = 'right', verticalalignment = 'center', fontsize = fontsize_plot, weight = 'bold', color = 'r') 
ax.text(pc_median * (1 + relative_error_target),  ( np.max(hist_pc) + np.min(hist_pc) ) / 2., '$\mathrm{P_{c, median} + error}$', rotation = 90, horizontalalignment = 'right', verticalalignment = 'center', fontsize = fontsize_plot, weight = 'bold', color = 'r') 


ax.text(0.57, 0.98, 'N = ' + str(nb_run) + '\nCI = ' + format(ci, ".0f") + "%\nerror = " + format(relative_error_target*100, ".1f")  + "%", transform = ax.transAxes, horizontalalignment = 'left', verticalalignment = 'top', fontsize = fontsize_plot, weight = 'bold')
ax.margins(0,0)


fig_save_name = 'pc_confidence_interval.pdf' 
fig.set_figheight(height_fig)
fig.set_figwidth(width_fig)

fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  




