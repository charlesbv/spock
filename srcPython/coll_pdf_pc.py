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
# this script genreates a PDF of Pc values. THe probability of a Pc to occur is defined by the probability of the set (f107, Ap) to occur. This is claculated in the script new_uncertainty_f107_ap_bis.py
# so first run new_uncertainty_f107_ap_bis.py to get hist_f107 and hist_ap for iforecast = 2 (we only look at hte distrbutions of f107 an ap prediction errors at day + 2
# report them in hist_f107 and hist_ap in this script (could have amde a function but we don't run new_uncertainty_f107_ap_bis many times, just once, since hist_f107 and hist_ap don't change (always computed over the same period: jan 1 2016 to nov 20 2016)
# the values of Pc are read from the file file_results_pc (genreated on Big by /raid4/cbv/collision/list_read_collision_file_new.py (although the SpOCK runs are made on Pleiades)

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
from mpl_toolkits.mplot3d import Axes3D



# The bin ranges are defined by the number of bins and the range in the functino .hist used in new_uncertainty_f107_ap_bis.py to caclaulte hist_f107 and hist_ap. So need to make sure that this f107_bin_range and ap_bin_range are the same as the result of the function .hist (coordinate [1] of the result of cuntion hist). This is also used in the script deviation_swpc_f107_ap.py
cf107_bin_range = np.array([-15.        , -10.71428571,  -6.42857143,  -2.14285714,         2.14285714,   6.42857143,  10.71428571,  15.        ])
ap_bin_range = f107_bin_range

f107_bin_center = (f107_bin_range[:-1] + np.roll(f107_bin_range, -1)[:-1])/2
ap_bin_center = (ap_bin_range[:-1] + np.roll(ap_bin_range, -1)[:-1])/2
nb_f107 = len(f107_bin_center)
nb_ap = len(ap_bin_center)

# same as in deviation_swpc_f107_ap.py
f107_bin_center_same_as_filename = []
ap_bin_center_same_as_filename = []
for if107 in range(nb_f107):
    f107_bin_center_same_as_filename.append( np.float(format(f107_bin_center[if107], ".0f")) ) 
for iap in range(nb_ap):
    ap_bin_center_same_as_filename.append( np.float(format(ap_bin_center[iap], ".0f")) )
f107_bin_center_same_as_filename = np.array(f107_bin_center_same_as_filename) # [-13.,  -9.,  -4.,   0.,   4.,   9.,  13.]
ap_bin_center_same_as_filename = np.array(ap_bin_center_same_as_filename) # [-13.,  -9.,  -4.,   0.,   4.,   9.,  13.]

                      # [-13.,         -9.,          -4.,         0.,         4.,           9.,         13.]
hist_f107 = np.array([ 0.02173913,  0.0621118 ,  0.20496894,  0.39751553,  0.1863354 ,   0.07763975,  0.0310559 ])
hist_ap = np.array([ 0.04037267,  0.04968944,  0.1552795 ,  0.35093168,  0.20496894,  0.0931677 ,  0.02484472])


filename_results_pc = 'results_run_paper4.txt'
file_results_pc = open(filename_results_pc)
read_file_results_pc = file_results_pc.readlines()
nb_run = len(read_file_results_pc)

pc = [] 
f107 = []
ap = []
tca = []
for irun in range(nb_run):
    pc_irun = []
    run_name = read_file_results_pc[irun].split()[0]
    f107.append( np.float(run_name.split('f107_')[1].split('_')[0]) )
    ap.append( np.float(run_name.split('ap_')[1].split('.')[0]) )
    pc_irun.append( np.float( read_file_results_pc[irun].split()[1] ) )
    tca.append( read_file_results_pc[irun].split()[2] )

    where_f107  = np.where(f107_bin_center_same_as_filename == f107[-1])[0][0]
    prob_f107 = hist_f107[where_f107]

    where_ap  = np.where(ap_bin_center_same_as_filename == ap[-1])[0][0]
    prob_ap = hist_ap[where_ap]

    pc_irun.append(prob_f107 * prob_ap)
    pc_irun.append(prob_f107)
    pc_irun.append(f107[-1])
    pc_irun.append(prob_ap)
    pc_irun.append(ap[-1])
    pc.append(pc_irun)
# for if107 in range(nb_f107):
#     for iap in range(nb_ap):
        
file_results_pc.close()

pc = np.array(pc) # pc[irun,:]:  0 pc, 1 probability of (f107,ap) to occur, 2 probablity of f107 to occur, 3 f107, 4 probablity of ap to occur, 5 ap
index_sort_pc = np.argsort(pc[:,0])
pc_sort = pc[index_sort_pc] # pc array sorted by increasing number of pc

index_sort_prob = np.argsort(pc[:,1])
pc_sort_on_prob = pc[index_sort_prob] # pc array sorted by increasing number of pc



nb_point_interpo_in_between_two_pc = 10
prob_interpolated = np.zeros([(nb_run-1) * nb_point_interpo_in_between_two_pc+1]) # +1 for the last run of pc_sort
pc_interpolated = np.zeros([(nb_run-1) * nb_point_interpo_in_between_two_pc+1])
for irun in range(nb_run-1):
    x1 = pc_sort[irun, 0]
    x2 = pc_sort[irun + 1, 0]
    y1 = pc_sort[irun, 1]
    y2 = pc_sort[irun + 1, 1]
    for iin in range(nb_point_interpo_in_between_two_pc):
        delta_y = y2 - y1
        delta_x = x2 - x1
        a = delta_y / delta_x
        b = y1 - a * x1        
        pc_interpolated[irun*nb_point_interpo_in_between_two_pc + iin] = x1 + delta_x / nb_point_interpo_in_between_two_pc * iin
        prob_interpolated[irun*nb_point_interpo_in_between_two_pc + iin] = a * pc_interpolated[irun*nb_point_interpo_in_between_two_pc + iin] + b
    
pc_interpolated[-1] = pc_sort[-1,0]
prob_interpolated[-1] = pc_sort[-1,1]

bin_size_cdf = 1.e-5
max_pc = 1.4e-4#np.max(pc_sort[:,0])
min_pc = 1e-5#np.min(pc_sort[:,0])
nb_bin_cdf = (int) ( ( max_pc - min_pc ) / bin_size_cdf )
cdf_per_bin  = np.zeros([nb_bin_cdf]) # probability of Pc to be in each bin (= for each bin, integral of 'probability VS Pc')
average_pc_per_bin  = np.zeros([nb_bin_cdf]) # probability of Pc to be in each bin (= for each bin, integral of 'probability VS Pc')
bin_cdf_array = np.zeros([nb_bin_cdf])
density_prob = np.zeros([nb_bin_cdf])
density_prob_with_interpo = np.zeros([nb_bin_cdf])
for ibin in range(nb_bin_cdf):
    bin_min = min_pc + ibin * bin_size_cdf
    bin_max = min_pc + ( ibin + 1 ) * bin_size_cdf
    where_pc_in_bin = np.where( ( pc_interpolated >= bin_min ) & ( pc_interpolated < bin_max ) )[0]
    nb_pc_in_bin = len(where_pc_in_bin) # nb of runs with Pc in this bin    
    pc_in_bin = pc_interpolated[where_pc_in_bin] # pc_interpolated is already sorted in ascending order
    prob_in_bin = prob_interpolated[where_pc_in_bin]
    density_prob_with_interpo[ibin] = np.sum(prob_in_bin) #np.sum(pc_in_bin * prob_in_bin)

    where_pc_in_bin = np.where( ( pc_sort[:,0] >= bin_min ) & ( pc_sort[:,0] < bin_max ) )[0]
    prob_in_bin = pc_sort[where_pc_in_bin, 1]
    density_prob[ibin] = np.sum(prob_in_bin) #np.sum(pc_in_bin * prob_in_bin)

    bin_cdf_array[ibin] = bin_min 
    for ipc in range(nb_pc_in_bin-1):
        x1 = pc_interpolated[ipc]
        x2 = pc_interpolated[ipc+1]
        y1 = prob_interpolated[ipc]
        y2 = prob_interpolated[ipc+1]
        area = y1 * (x2 - x1) + (x2 - x1) * ( y2 - y1 ) / 2.
        cdf_per_bin[ibin] =  cdf_per_bin[ibin] + area
        average_pc_per_bin[ibin] = average_pc_per_bin[ibin] + prob_interpolated[ipc]
    average_pc_per_bin[ibin] = average_pc_per_bin[ibin] + prob_interpolated[nb_pc_in_bin-1]
    average_pc_per_bin[ibin] = average_pc_per_bin[ibin] / nb_pc_in_bin
raise Exception
# PLOTS
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3
### For plots, generate disctinct colors
NCURVES = len(hist_f107)
np.random.seed(101)
curves = [np.random.random(20) for i in range(NCURVES)]
values = range(NCURVES)
jet = cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)


## densityof probbiliy
fig_title = 'Probability per $P_C$ bin'
y_label = 'Probability (%)'
x_label = '$P_C$'
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
 
#hist_pdf = ax.hist(pc_sort[:, 0])[0] / nb_run * 100
ax.bar(bin_cdf_array, density_prob*100., bin_size_cdf, 0 , linewidth = 2, align = 'edge')
ax.margins(0,0)
# pc[irun,:]:  0 pc, 1 probability of (f107,ap) to occur, 2 probablity of f107 to occur, 3 f107, 4 probablity of ap to occur, 5 ap

fig_save_name = 'probability_per_pc_bin'  
fig_save_name =  fig_save_name + '.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


## Probability of Pc ocurring as a funcitno of Pc
fig_title = 'Probability of scenario (F10.7/Ap) VS $P_C$ for this scenario'
y_label = 'Probability (%)'
x_label = '$P_C$'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.98,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

#ax.plot(pc_sort[:, 0], pc_sort[:, 1] , linewidth = 2, color = 'k')
ax.scatter(pc_sort[:, 0], pc_sort[:, 1]*100 , linewidth = 2, color = 'b')
for irun in range(nb_run):
    if  ( ( pc_sort[irun, 1]*100 < 0.1 ) | ( pc_sort[irun, 1]*100 > 6 ) ):
        ax.text(pc_sort[irun, 0], pc_sort[irun, 1]*100, str(pc_sort[irun, 3]) + ', ' + str( pc_sort[irun, 5]), 
                weight = 'bold', fontsize = fontsize_plot/1.2, horizontalalignment = 'center', verticalalignment = 'bottom', color =  'r')

for ibin in range(nb_bin_cdf):
    bin_min = min_pc + ibin * bin_size_cdf
    bin_max = min_pc + ( ibin + 1 ) * bin_size_cdf
    where_pc_in_bin = np.where( ( pc_interpolated >= bin_min ) & ( pc_interpolated < bin_max ) )[0]
    if len(where_pc_in_bin) > 0:
        nb_pc_in_bin = len(where_pc_in_bin) # nb of runs with Pc in this bin    
        pc_in_bin = pc_interpolated[where_pc_in_bin] # pc_interpolated is already sorted in ascending order
        prob_in_bin = prob_interpolated[where_pc_in_bin]
        ax.plot([pc_in_bin[0], pc_in_bin[0]], [0, prob_in_bin[0]*100], color = 'k', linewidth = 2, linestyle = 'dotted')
ax.plot(pc_interpolated, prob_interpolated*100, linewidth = 2, color = 'b')
ax.set_xlim([np.min(pc_sort[:, 0]), np.max(pc_sort[:, 0])])

# hist_pdf = ax.hist(pc_sort[:, 0])[0] / nb_run * 100
ax.set_ylim([0, np.max(pc[:,1])*100])
ax.set_xlim([min_pc, max_pc-bin_size_cdf])
ax.margins(0,0)
# pc[irun,:]:  0 pc, 1 probability of (f107,ap) to occur, 2 probablity of f107 to occur, 3 f107, 4 probablity of ap to occur, 5 ap

fig_save_name = 'probability_vs_pc'  
fig_save_name =  fig_save_name + '.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  





## Probability of F10.7 to occur
fig_title = 'Probability of F10.7 scenario to occur'
y_label = 'Probability (%)'
x_label = 'F10.7'
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
 

bar_f107 = ax.bar(f107_bin_center, hist_f107*100., f107_bin_center[1] - f107_bin_center[0], 0 , linewidth = 2)
nb_bar = len(bar_f107)
for ibar in range(nb_bar):
    colorVal = scalarMap.to_rgba(ibar)
    bar_f107[ibar].set_color(colorVal)

ax.margins(0,0)

fig_save_name = 'f107_probability'  
fig_save_name =  fig_save_name + '.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  



## Probability of Ap to occur
fig_title = 'Probability of Ap scenario to occur'
y_label = 'Probability (%)'
x_label = 'Ap'
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
 

bar_ap = ax.bar(ap_bin_center, hist_ap*100., ap_bin_center[1] - ap_bin_center[0], 0 , linewidth = 2)
nb_bar = len(bar_ap)
for ibar in range(nb_bar):
    colorVal = scalarMap.to_rgba(ibar)
    bar_ap[ibar].set_color(colorVal)

ax.margins(0,0)

fig_save_name = 'ap_probability'  
fig_save_name =  fig_save_name + '.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

raise Exception

## Pc  as a function of Probability of Pc occuring
fig_title = ''
x_label = 'Probability'
y_label = '$P_C$'
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

ax.plot(pc_sort_on_prob[:, 1], pc_sort_on_prob[:, 0] , linewidth = 2, color = 'k')
ax.scatter(pc_sort_on_prob[:, 1], pc_sort_on_prob[:, 0] , linewidth = 2, color = 'k')
for irun in range(nb_run):
    ax.text(pc_sort_on_prob[irun, 1], pc_sort_on_prob[irun, 0], str(pc_sort_on_prob[irun, 3]) + ', ' + str( pc_sort_on_prob[irun, 5]), weight = 'bold', fontsize = fontsize_plot/1.2, horizontalalignment = 'center', verticalalignment = 'bottom', color =  'k')
ax.set_xlim([np.min(pc_sort_on_prob[:, 1]), np.max(pc_sort_on_prob[:, 1])])
ax.margins(0,0)
# pc[irun,:]:  0 pc, 1 probability of (f107,ap) to occur, 2 probablity of f107 to occur, 3 f107, 4 probablity of ap to occur, 5 ap

fig_save_name = 'pc_vs_prob'  
fig_save_name =  fig_save_name + '.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# Probability as a function of Ap anf F10.7
fig_title = ''
z_label = 'Probability'
y_label = 'Ap'
x_label = 'F10.7'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0], projection = '3d')

ax.set_zlabel(z_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

X, Y = np.meshgrid(f107_bin_center, ap_bin_center)


#ax.plot(pc_sort_on_prob[:, 3], pc_sort_on_prob[:, 5], pc_sort_on_prob[:, 1] , linewidth = 2, color = 'k')

fig_save_name = 'prob_vs_f107_ap'  
fig_save_name =  fig_save_name + '.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

