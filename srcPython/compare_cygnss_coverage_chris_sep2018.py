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

# This scirpt is a copy of cygnss_coverage_temporal_spatial_relation_equal_cell_area2 on Sep 16 2018. It was made to answer Chris Ruf's analysis questions by email on Sep 13 2018.                                                          

# IMPORTANT: THE FOLDER WHERE I STORED THE SSIMULATIONS AND PICLES WAS TO BIG TO KEEP ON MY LAPTOP SO I MOVED IT TO BIG: /raid4/cbv/cygnss/coverage_temporal_spatial_relation. I did that after sumbitting the paper to JSTAR. Note that I moved this folder only from my aptop but the one on the desktop is still on the desktop.

import matplotlib
matplotlib.use("Agg") # without this line, when running this script from a cronjob we get an error "Unable to access the X Display, is $DISPLAY set properly?"
import os
import sys 
sys.path.append("/home1/cbussy/spock/srcPython")
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
#from mpl_toolkits.basemap import Basemap, shiftgrid # to comment if pleaides
from collections import *
from matplotlib import colors
import matplotlib.ticker as ticker
import pickle
import matplotlib.patches as mpatches

plt.ion()

# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
run_list = ["all_nadir", "all_nadir_no_gpsIIF", "all_nadir_no_fm01", "all_nadir_no_fm01_no_gpsIIF", "all_nadir_no_fm05_no_fm08", "all_nadir_no_fm05_no_fm08_no_gpsIIF"] #  list of cygnss_coverage_temporal_spatial_relation.py runs to compare. !!!!! need to run cygnss_coverage_temporal_spatial_relation_equal_cell_area2.py for each of these runs first
# for cyggnss extended: run_list = ["1_inc_35","2_ortho_inc_35_35","2_ortho_inc_35_65","2_ortho_inc_35_75","2_ortho_inc_35_85", "2_ortho_inc_35_35_4sceach", "2_ortho_inc_35_45_4sceach"] -> pickles in  /Users/cbv/coverage_temporal_spatial_relation/pickle2weeks/ (generated on on Pleiades)
# for inclinatino at 90 deg: run_list = ["1_inc_90", "2_ortho_inc_90_90", "4_ortho_inc_90_90_90_90" , "1_inc_65" , "1_inc_70", "1_inc_75", "1_inc_80", "1_inc_85"] -> pickles in /Users/cbv/coverage_temporal_spatial_relation/pickle2weeks/ (generated on on Pleiades)
max_lat_grid_array =    [5, 10, 15, 20, 25, 30]       # !!!!!! needs to be the same as in cygnss_coverage_temporal_spatial_relation.py
# for cyggnss extended: max_lat_grid_array =  [5, 10, 15, 20, 25, 30]
# for inclinatino at 90 deg: max_lat_grid_array = [60,65,70,75,80,85] except for :
# for the histogram of revisit time: [5]
# 1_inc_60: [60,65]
# 1_inc_65: [60,65,70]
# 1_inc_70: [60,65,70,75]
root_fig = "sep18"
# for cyggnss extended: # for cygnss_extended
#for inclinatino at 90 deg:: 90
#If pleaides
pickle_folder = "/nobackup/cbussy/coverage_temporal_spatial_relation/sep18/pickle" # with or without the last '/' (doesnt matter) the 2 pcikle to try are pickle_1215 and pickle_1217 (dont remember the difference between both, i think it depedns if you lok at low or high lat). There is also pickle_0127 but this one is for cells of 50 km and 200 km wide to compare two planes vs one plane at low latitudes. I actually had to remove pickle_0127 because no more space on computer. It's still on Pleiades though.
plot_folder = "/nobackup/cbussy/coverage_temporal_spatial_relation/sep18/plot"

# # If not pleaides:
# pickle_folder = "/Users/cbv/cygnss/coverage_temporal_spatial_relation/sep18/pickle" # with or without the last '/' (doesnt matter) the 2 pcikle to try are pickle_1215 and pickle_1217 (dont remember the difference between both, i think it depedns if you lok at low or high lat). There is also pickle_0127 but this one is for cells of 50 km and 200 km wide to compare two planes vs one plane at low latitudes. I actually had to remove pickle_0127 because no more space on computer. It's still on Pleiades though.
# plot_folder = "/Users/cbv/cygnss/coverage_temporal_spatial_relation/sep18/plot"

pickle_visu_subfolder = "2d_visu" # with or without the last '/' (doesnt matter)
plot_visu_subfolder = "2d_visu_35"
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
run_list_arr = np.array(run_list)
label_arr = run_list_arr
label_arr[np.where(run_list_arr == "all_nadir")[0]] = "All FMs"
label_arr[np.where(run_list_arr == "all_nadir_no_gpsIIF")[0]] = 'All FMs, no IIF'
label_arr[np.where(run_list_arr == "all_nadir_no_fm01")[0]] = "No FM01"
label_arr[np.where(run_list_arr == "all_nadir_no_fm01_no_gpsIIF")[0]] = "No FM01, no IIF" 
label_arr[np.where(run_list_arr == "all_nadir_no_fm05_no_fm08")[0]] = "No FM05/FM08"
label_arr[np.where(run_list_arr == "all_nadir_no_fm05_no_fm08_no_gpsIIF")[0]] = "No FM05/FM08, no IIF"


if plot_folder[-1] != '/':
    plot_folder = plot_folder + '/'

if os.path.isdir(plot_folder) == False:
    os.system("mkdir " + plot_folder)


if os.path.isdir(plot_folder + plot_visu_subfolder) == False:
    os.system("mkdir " + plot_folder + plot_visu_subfolder)

if plot_visu_subfolder[-1] != '/':
    plot_visu_subfolder = plot_visu_subfolder + '/'

if os.path.isdir(plot_folder + plot_visu_subfolder + 'ani') == False:
    os.system("mkdir " + plot_folder + plot_visu_subfolder + 'ani')

folder_to_create = plot_folder + 'coverage_vs_time/'
if os.path.isdir(folder_to_create) == False:
    os.system("mkdir " + folder_to_create)

folder_to_create = plot_folder + 'time_saved/'
if os.path.isdir(folder_to_create) == False:
    os.system("mkdir " + folder_to_create)

folder_to_create = plot_folder + 'revisit_coverage/'
if os.path.isdir(folder_to_create) == False:
    os.system("mkdir " + folder_to_create)

folder_to_create = plot_folder + 'revisit_time/'
if os.path.isdir(folder_to_create) == False:
    os.system("mkdir " + folder_to_create)


folder_to_create = plot_folder + 'time_to_reach_cov/'
if os.path.isdir(folder_to_create) == False:
    os.system("mkdir " + folder_to_create)





if pickle_folder[-1] != '/':
    pickle_folder = pickle_folder + '/'

if pickle_visu_subfolder[-1] != '/':
    pickle_visu_subfolder = pickle_visu_subfolder + '/'




max_lat_grid_array = np.array(max_lat_grid_array)

width_cell_array =  [25] # in km#width_cell_array = [0.01, 0.05, 0.1, 0.15, 0.20, 0.25] #!!!!!! should be [0.01, 0.05, 0.1, 0.15, 0.20, 0.25] # in degrees. for the histogram of revisit time: [200]
coverage_array = [95] # in pecentage in ASCENDING ORDER
color_array = ['k','cornflowerblue','r','g', 'm', 'gold', 'cyan']

nb_type_cell = len(width_cell_array)
nb_coverage = len(coverage_array)

width_cell_array = np.array(width_cell_array)
coverage_array = np.array(coverage_array)

deg2rad = np.pi/180
    
# Sel<ecting only the specular ponts in the window (also called grid) defined by lon/lat_grid_center and lon/lat_width, count the number of specular points in each cell (defined by width_cell) as a function of time
nb_lat_center = len(max_lat_grid_array)

nb_spec = 4# ecef_spec.shape[0]
time_to_reach_coverage_average_over_lon_all_run = []
nb_cell_covered_average_over_lon_all_run = []
mean_revisit_time_average_over_lon_all_run = []
percentage_cell_at_least_one_revisit_average_over_lon_all_run = []
time_since_last_visit_of_this_cell_list_all_run = []
nb_run = len(run_list)
raise Exception

for irun in range(nb_run):
    print irun, nb_run-1
    nb_cell_covered_average_over_lon_all_run.append(pickle.load(open(pickle_folder + run_list[irun] + "_nb_cell_covered_average_over_lon_gain_not0.pickle"))) 
    percentage_cell_at_least_one_revisit_average_over_lon_all_run.append(pickle.load(open(pickle_folder + run_list[irun] + "_percentage_cell_at_least_one_revisit_average_over_lon_gain_not0.pickle")))


nb_time = nb_cell_covered_average_over_lon_all_run[0].shape[2] # all run have same duration

#raise Exception
# !!!!!!!!! TEMP
time_to_reach_coverage_average_over_lon_again_all_run = []
for irun in range(nb_run):
    time_to_reach_coverage_average_over_lon_again  = np.zeros([nb_lat_center, nb_type_cell, nb_coverage])
    for ilat_center in range(nb_lat_center):
        for itype in range(nb_type_cell):
            for icoverage in range(nb_coverage):
                if len(np.where(nb_cell_covered_average_over_lon_all_run[irun][ilat_center, itype, :] >= coverage_array[icoverage])[0]) == 0: 
                    # this means the coverage goal coverage_array[icoverage] was not reached for this type of cell at this latitude of the grid for this run
                    time_to_reach_coverage_average_over_lon_again[ilat_center, itype, icoverage] = -1 # in this case set the time to reach coverage goal to -1
                else:
                    time_to_reach_coverage_average_over_lon_again[ilat_center, itype, icoverage] = np.where(nb_cell_covered_average_over_lon_all_run[irun][ilat_center, itype, :] >= coverage_array[icoverage])[0][0]
    time_to_reach_coverage_average_over_lon_again_all_run.append(time_to_reach_coverage_average_over_lon_again)
# !!!!!!!!! end of TEMP



# Histogram of revisted time for a given lat and cell type - ok sep18
ilat_center = 0
itype = 0
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25 
ratio_fig_size = 4./3    
fig_title = ''#PDF of revisit times (latitude of grid: ' + str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN}' + ' - cell size: ' + str(width_cell_array[itype]) + ' km)'
y_label = 'PDF'
x_label = 'Revisit time (hours)'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.961,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold


color_array_hist = ['b','r']
max_y = 0
nbins = 50
binsize = 45 # in minute
binsize =binsize / 60.
for irun in range(nb_run):
    time_since_last_visit_of_this_cell_list = pickle.load(open(pickle_folder + run_list[irun] + "_time_since_last_visit_of_this_cell_list_ilat_" + str(max_lat_grid_array[ilat_center]) + "_cell_" + str(width_cell_array[itype]) + "_gain_not0.pickle"))
    mean_revisit_time_average_over_lon_all_run.append(pickle.load(open(pickle_folder + run_list[irun] + "_mean_revisit_time_average_over_lon_gain_not0.pickle"))) 
    nb_lon = len(time_since_last_visit_of_this_cell_list) # nb longitude the grid was moved over
    time_since_last_visit_of_this_cell_all_lon = time_since_last_visit_of_this_cell_list[0]
    for ilon in range(1,nb_lon):
        time_since_last_visit_of_this_cell_all_lon = time_since_last_visit_of_this_cell_all_lon + time_since_last_visit_of_this_cell_list[ilon]

    time_since_last_visit_of_this_cell_all_lon = np.array(time_since_last_visit_of_this_cell_all_lon) / 3600.
    nbins = (int)((np.max(time_since_last_visit_of_this_cell_all_lon) - np.min(time_since_last_visit_of_this_cell_all_lon))/binsize)
    hist_rev = np.histogram(time_since_last_visit_of_this_cell_all_lon, bins = nbins, density = 1)[0]#ax.hist(time_since_last_visit_of_this_cell_all_lon, bins = 20, color = color_array_hist[irun], label = label_arr[irun])
    #hist_rev = hist_rev / np.sum(time_since_last_visit_of_this_cell_all_lon) 
    bins_hist_rev_temp = np.histogram(time_since_last_visit_of_this_cell_all_lon, bins = nbins)[1]
    bins_hist_rev = ( bins_hist_rev_temp[:-1] + np.roll(bins_hist_rev_temp,-1)[:-1] ) /2.
    binsize_actual = bins_hist_rev[1] - bins_hist_rev[0]
    #print np.mean(time_since_last_visit_of_this_cell_all_lon_greater_one_period), np.min(time_since_last_visit_of_this_cell_all_lon_greater_one_period) * 60.
    ax.plot(bins_hist_rev, hist_rev,color = color_array_hist[irun], label = label_arr[irun], linewidth = 2)
    if np.max(hist_rev) > max_y:
        max_y = np.max(hist_rev)

    
#    ax.plot( [np.mean(time_since_last_visit_of_this_cell_all_lon), np.mean(time_since_last_visit_of_this_cell_all_lon)], [0,1], linewidth = 1, linestyle = 'dashed', color =color_array_hist[irun])

plt.yscale('log')
#ax.set_ylim([0,max_y*1.1])
legend = ax.legend(loc='upper right', bbox_to_anchor=(1, 1), numpoints = 1,  title="Constellation", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = plot_folder + 'revisit_time/histogram_revisit_time_' + "ilat_" + str(max_lat_grid_array[ilat_center]) + "_cell_" + str(width_cell_array[itype]) + ".pdf"
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

# Histogram of revisted time for a given lat and cell type but for revisti times > one orbital perriod (assumed of 95 min)
ilat_center = 0
itype = 0
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25 
ratio_fig_size = 4./3    
fig_title = 'Histogram of revisit times > $T_{orb}$ (latitude of grid: ' + str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN}' + ' - cell size: ' + str(width_cell_array[itype]) + ' km)'
y_label = 'PDF'
x_label = 'Revisit time (hours)'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.961,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold


color_array_hist = ['b','r']
max_y = 0
nbins = 50
binsize = 30 # in minute
binsize =binsize / 60.
min_revisit_time = 100.# ignore reivist time smaller than this. in min
min_revisit_time = min_revisit_time / 60
for irun in range(nb_run):
    time_since_last_visit_of_this_cell_list = pickle.load(open(pickle_folder + run_list[irun] + "_time_since_last_visit_of_this_cell_list_ilat_" + str(max_lat_grid_array[ilat_center]) + "_cell_" + str(width_cell_array[itype]) + "_gain_not0.pickle"))
    mean_revisit_time_average_over_lon_all_run.append(pickle.load(open(pickle_folder + run_list[irun] + "_mean_revisit_time_average_over_lon_gain_not0.pickle"))) 
    nb_lon = len(time_since_last_visit_of_this_cell_list) # nb longitude the grid was moved over
    time_since_last_visit_of_this_cell_all_lon = time_since_last_visit_of_this_cell_list[0]
    for ilon in range(1,nb_lon):
        time_since_last_visit_of_this_cell_all_lon = time_since_last_visit_of_this_cell_all_lon + time_since_last_visit_of_this_cell_list[ilon]

    time_since_last_visit_of_this_cell_all_lon = np.array(time_since_last_visit_of_this_cell_all_lon) / 3600.
    time_since_last_visit_of_this_cell_all_lon_greater_one_period = time_since_last_visit_of_this_cell_all_lon[np.where(time_since_last_visit_of_this_cell_all_lon > (min_revisit_time))]
    nbins = (int)((np.max(time_since_last_visit_of_this_cell_all_lon_greater_one_period) - np.min(time_since_last_visit_of_this_cell_all_lon_greater_one_period))/binsize)
    hist_rev = np.histogram(time_since_last_visit_of_this_cell_all_lon_greater_one_period, bins = nbins, density = 1)[0]#ax.hist(time_since_last_visit_of_this_cell_all_lon_greater_one_period, bins = 20, color = color_array_hist[irun], label = label_arr[irun])
    #hist_rev = hist_rev / np.sum(time_since_last_visit_of_this_cell_all_lon_greater_one_period) 
    bins_hist_rev_temp = np.histogram(time_since_last_visit_of_this_cell_all_lon_greater_one_period, bins = nbins)[1]
    bins_hist_rev = ( bins_hist_rev_temp[:-1] + np.roll(bins_hist_rev_temp,-1)[:-1] ) /2.
    binsize_actual = bins_hist_rev[1] - bins_hist_rev[0]
    print np.mean(time_since_last_visit_of_this_cell_all_lon_greater_one_period), np.min(time_since_last_visit_of_this_cell_all_lon_greater_one_period) * 60.
    ax.bar(bins_hist_rev, hist_rev, binsize_actual,color = color_array_hist[irun], label = label_arr[irun])
    if np.max(hist_rev) > max_y:
        max_y = np.max(hist_rev)

    
    ax.plot( [np.mean(time_since_last_visit_of_this_cell_all_lon_greater_one_period), np.mean(time_since_last_visit_of_this_cell_all_lon_greater_one_period)], [0,np.max(hist_rev)], linewidth = 2, linestyle = 'dashed', color = 'k')#color = color_array_hist[irun])

plt.yscale('log')
#ax.set_ylim([0,max_y*1.02])
legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Constellation", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))


fig_save_name = plot_folder + 'revisit_time/histogram_revisit_time_greater_one_period_' + "ilat_" + str(max_lat_grid_array[ilat_center]) + "_cell_" + str(width_cell_array[itype]) + ".pdf"
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# Latitude of grid vs Time to reach coverage goal for different run, for a given coverage goal
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25 
ratio_fig_size = 4./3

icoverage = 0
irun  = 0 # constellation configuration to look at the time to reach coverage

fig_title = ''#Width of cell VS # of hours to reach ' + str(coverage_array[icoverage]) + '% coverage for different latitudes of grid'
x_label = 'Latitude of grid (' + u'\N{DEGREE SIGN})'
y_label = 'Number of hours to reach coverage goal (hours)'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.96,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
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
min_y = 1e30
for irun in range(nb_run):

    where_gt_0 = np.where((time_to_reach_coverage_average_over_lon_again_all_run[irun][ilat_center, :, icoverage]/3600.  > 0)) # -1 if the coverage goal has never been reached for this type of cell 
    toplot = time_to_reach_coverage_average_over_lon_again_all_run[irun][:, 0, icoverage]
    ax.scatter( max_lat_grid_array, toplot/3600. ,linewidth = 2, color = color_array[irun], marker = 'o', s = 100)
    ax.loglog( max_lat_grid_array, toplot/3600. ,linewidth = 3, color = color_array[irun], label = label_arr[irun])
    if np.max(toplot/3600.)> max_y:
        max_y = np.max(toplot/3600. )
    if np.min(toplot/3600.) < min_y:
        min_y = np.min(toplot/3600. )

ax.minorticks_off()
yticks = [0,5,10,20,40,50, 60, 70, 80, 90, 100,120,200, 360]
yticks_label = []
for i in range(len(yticks)):
        yticks_label.append( format( yticks[i], ".0f" ) )
ax.yaxis.set_ticks(yticks)
ax.yaxis.set_ticklabels(yticks_label, fontsize = fontsize_plot)

xlim_save = [ax.get_xlim()[0], ax.get_xlim()[1]]
xticks = max_lat_grid_array
xticks_label = []
for i in range(len(xticks)):
    xticks_label.append( str( max_lat_grid_array[i] ) )
ax.xaxis.set_ticks(xticks)
ax.xaxis.set_ticklabels(xticks_label, fontsize = fontsize_plot)#, rotation='vertical')

ax.set_ylim([min_y*0.9, max_y*1.1])
#ax.set_xlim(xlim_save)
ax.set_xlim([np.min(max_lat_grid_array), np.max(max_lat_grid_array)])
legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))


fig_save_name = plot_folder + 'time_to_reach_cov/time_to_reach_coverage_' + str(coverage_array[icoverage]).replace(".","_") + '_gain_not0.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  



# Width of bin vs Time to reach coverage goal for different latitude center of the window, for one coverage -> compare each run to oneplane
## Parameters for the figure
height_fig = 12.5  # the width is calculated as height_fig * 4/3.
fontsize_plot = 18
ratio_fig_size = 1#4./3

icoverage = 2

#fig_title = 'Width of cell VS decrease in number of hours to reach coverage goal for different\nconstellation configurations compared to the ' +  str(run_list[0]) + ' configuration\nLatitude of grid ' + str(max_lat_grid_array[ilat_center])+ u'\N{DEGREE SIGN}' + '  - Coverage goal ' + str(coverage_array[icoverage]) + '%'
fig_title = ''
y_label = 'Size of cell (km)'
x_label = 'Decrease in # hrs to get ' + str(coverage_array[icoverage]) + '% cov.'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 1.023,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(2, 2)
gs.update(left = 0., right=1.05, top = 0.93,bottom = 0.12, hspace = 0.21, wspace = 0.03)
ilat_count = -1
for ilat_center in [0,-3,-2,-1]:#range(0,nb_lat_center):#nb_lat_center):    
    ilat_count = ilat_count + 1
    if ilat_count == 0:
        ax = fig.add_subplot(gs[0, 0])
        ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    elif ilat_count == 1:
        ax = fig.add_subplot(gs[0, 1])
    elif ilat_count == 2:
        ax = fig.add_subplot(gs[1, 0])
        ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
        ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
    else:
        ax = fig.add_subplot(gs[1, 1])
        ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

    ax.set_title('Latitude of grid: ' + str(max_lat_grid_array[ilat_center])+ u'\N{DEGREE SIGN}', weight = 'bold', fontsize = fontsize_plot)
    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

    max_x = 0
    for irun in range(1, nb_run):# !!!!!!!! range(1,nb_run):
    #for ilat_center in range(nb_lat_center): # !!!!!!!!!!!# a value is negative if the coverage never reached the coverage goal. SO we don't want it on the plot
        where_gt_0 = np.where((time_to_reach_coverage_average_over_lon_again_all_run[irun][ilat_center, :, icoverage]/3600.  > 0) & (time_to_reach_coverage_average_over_lon_again_all_run[0][ilat_center, :, icoverage]/3600.  > 0)) # -1 if the coverage goal has never been reached for this type of cell 
        toplot = (time_to_reach_coverage_average_over_lon_again_all_run[0][ilat_center, where_gt_0[0], icoverage]-time_to_reach_coverage_average_over_lon_again_all_run[irun][ilat_center, where_gt_0[0], icoverage])
        where_time_save_gt_0 = np.where( toplot > 0 )[0] # otherwise it means that the second constellation took more time to reach the coverage goal than the baseline constellation. Plotting this neagative time screws up the log plot so we don't plot it
        ax.scatter(toplot[where_time_save_gt_0]/3600. , width_cell_array[where_gt_0[0]][where_time_save_gt_0], linewidth = 2, color = color_array[irun], marker = 'o', s = 100)
        ax.loglog(toplot[where_time_save_gt_0]/3600. , width_cell_array[where_gt_0[0]][where_time_save_gt_0], linewidth = 3, color = color_array[irun], label = label_arr[irun] + u'\N{DEGREE SIGN}')
        if len(where_time_save_gt_0)>=1:
            if np.max(toplot[where_time_save_gt_0]/3600. ) > max_x:
                max_x = np.max(toplot[where_time_save_gt_0]/3600. )
#        print width_cell_array[where_gt_0[0]]
#        print np.min(toplot[where_time_save_gt_0]/3600), np.max(toplot[where_time_save_gt_0]/3600.)
#         if irun >= nb_run-2:
#             print time_to_reach_coverage_average_over_lon_again_all_run[0][ilat_center, where_gt_0[0], icoverage]/3600. , time_to_reach_coverage_average_over_lon_again_all_run[irun][ilat_center, where_gt_0[0], icoverage]/3600. , toplot/3600.
#             print ""
#    print max_lat_grid_array[ilat_center],toplot/3600. 
#    print time_to_reach_coverage_average_over_lon_again_all_run[0][ilat_center, where_gt_0[0], icoverage]/3600. , time_to_reach_coverage_average_over_lon_again_all_run[irun][ilat_center, where_gt_0[0], icoverage]/3600. 
    ax.minorticks_off()
#    xticks = np.logspace(0,np.log(np.max(nb_time) / 3600), num=10, base=np.exp(1))
#    xticks = np.arange(0,np.max(nb_time) / 3600, 20)
#    xticks = np.arange(0,max_x, 20)
    xticks = [0,5,10,20,40,60,80,100,140, 180]
    xticks_label = []
    for i in range(len(xticks)):
            xticks_label.append( format( xticks[i], ".0f" ) )
    ax.xaxis.set_ticks(xticks)
    ax.xaxis.set_ticklabels(xticks_label, fontsize = fontsize_plot)
    
    ylim_save = [ax.get_ylim()[0], ax.get_ylim()[1]]
    yticks = width_cell_array
    yticks_label = []
    for i in range(len(yticks)):
        yticks_label.append( str( width_cell_array[i] ) )
    ax.yaxis.set_ticks(yticks)
    if ( ilat_count == 0 ) | ( ilat_count == 2 ):
        ax.yaxis.set_ticklabels(yticks_label, fontsize = fontsize_plot)#, rotation='vertical')
    else:
        ax.yaxis.set_ticklabels([])

    ax.set_xlim([0, 190])#max_x*1.1])
    ax.set_ylim(ylim_save)
    ax.margins(0,0)

legend = ax.legend(title="", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = plot_folder + 'time_saved/' + root_fig + '_compare_to_one_plane_diff_latitude_coverage_' + str(coverage_array[icoverage]).replace(".","_") + '_several_latitude.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


raise Exception




# Coverage VS time for a given constellation configuration. Min/max latitude and min/amx grid cell shown
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25 
ratio_fig_size = 4./3

icoverage = 0
fig_title = ''#Coverage VS time for cells 5 km or 25 km wide\nLatitude of grid at 5'+ u'\N{DEGREE SIGN} or 30'+ u'\N{DEGREE SIGN}'
y_label = 'Coverage (%)'
x_label = 'Time (days)'

irun = 0 # constellation to look at

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.99,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

xaxis = np.arange(0, nb_time,1) / 3600. / 24

itype = 0
ilat_center = np.where(max_lat_grid_array == 5)[0][0]
ax.plot(xaxis, nb_cell_covered_average_over_lon_all_run[irun][ilat_center, itype, :], linewidth = 3, color = 'r', label = str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN} | ' +str( width_cell_array[itype]) + ' km', linestyle = 'dashed')
ilat_center = np.where(max_lat_grid_array == 30)[0][0]
ax.plot(xaxis, nb_cell_covered_average_over_lon_all_run[irun][ilat_center, itype, :], linewidth = 3, color = 'r', label = str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN} | ' +str( width_cell_array[itype]) + ' km')

itype = -1
ilat_center = np.where(max_lat_grid_array == 5)[0][0]
ax.plot(xaxis, nb_cell_covered_average_over_lon_all_run[irun][ilat_center, itype, :], linewidth = 3, color = 'b', label = str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN} | ' +str( width_cell_array[itype]) + ' km', linestyle = 'dashed')
ilat_center = np.where(max_lat_grid_array == 30)[0][0]
ax.plot(xaxis, nb_cell_covered_average_over_lon_all_run[irun][ilat_center, itype, :], linewidth = 3, color = 'b', label = str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN} | ' +str( width_cell_array[itype]) + ' km')


ax.plot([0, xaxis[-1]], [coverage_array[icoverage], coverage_array[icoverage]], linewidth = 2, linestyle = 'dashed', color = 'k')
ax.minorticks_off()
xticks = np.arange(0,np.max(nb_time) / 3600. / 24., 1)
xticks_label = []
for i in range(len(xticks)):
    xticks_label.append( format( xticks[i], ".0f" ) )
ax.xaxis.set_ticks(xticks)
ax.xaxis.set_ticklabels(xticks_label, fontsize = fontsize_plot)

ax.set_xlim([0, np.max(nb_time) / 3600. / 24.])


legend = ax.legend(loc='lower right', bbox_to_anchor=(1, 0), numpoints = 1,  title="Lat. grid | Cell size", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = plot_folder + 'coverage_vs_time/' + run_list[irun]+ 'coverage_vs_time_cell_width_min_and_max_latitude_min_and_max_TRY.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  




# Coverage VS time for a given type cell and latitude of grid. All runs plotted
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25 
ratio_fig_size = 4./3

ilat_center = np.where(max_lat_grid_array == 5)[0][0]       
itype = 0
icoverage = 0
fig_title = ''#Coverage VS time for different constellations\nLatitude of grid ' + str(max_lat_grid_array[ilat_center])+ u'\N{DEGREE SIGN}' + '  - Width of cell ' + str(width_cell_array[itype]) + ' km'
y_label = 'Coverage (%)'
x_label = 'Time (days)'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.99,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_title('Latitude of grid: ' + str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN}', weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

xaxis = np.arange(0, nb_time,1) / 3600. / 24
# ax.plot(0,0, linewidth = 3, color = 'b', label = 'Visited')
# ax.plot(0,0, linewidth = 3, color = 'b', linestyle = 'dashed', label = 'Revisited')
for irun in range(nb_run):# !!!!!!!! range(1,nb_run):
    ax.plot(xaxis, nb_cell_covered_average_over_lon_all_run[irun][ilat_center, itype, :], linewidth = 3, color = color_array[irun], label = label_arr[irun])
    ax.plot(xaxis, percentage_cell_at_least_one_revisit_average_over_lon_all_run[irun][ilat_center, itype, :], linewidth = 3, color = color_array[irun], linestyle = 'dashed')
ax.plot([0, xaxis[-1]], [coverage_array[icoverage], coverage_array[icoverage]], linewidth = 2, linestyle = 'dashed')
ax.minorticks_off()
xticks = np.arange(0,np.max(nb_time) / 3600. / 24., 1)
xticks_label = []
for i in range(len(xticks)):
    xticks_label.append( format( xticks[i], ".0f" ) )
ax.xaxis.set_ticks(xticks)
ax.xaxis.set_ticklabels(xticks_label, fontsize = fontsize_plot)

#ax.set_xlim([0, np.max(nb_time) / 3600. / 24.])
ax.set_xlim([0, 7])


# legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="", fontsize = fontsize_plot)#Orbit inc.
# legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = plot_folder + 'coverage_vs_time/' + root_fig + '_compare_coverage_vs_time_cellwidth_' + str(width_cell_array[itype]).replace(".","_") + '_latitude_' + str(max_lat_grid_array[ilat_center]) + '_gain_not0.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# RevisitcCoverage VS time for a given constellation configuration. Min/max latitude and min/amx grid cell shown
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25 
ratio_fig_size = 4./3

icoverage = 2
fig_title = ''#Revisit coverage VS time for cells 5 km or 25 km wide\nLatitude of grid at 5'+ u'\N{DEGREE SIGN} or 30'+ u'\N{DEGREE SIGN}'

irun = 0 # constellation to look at

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.99,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])
y_label = 'Revisit coverage (%)'
x_label = 'Time (days)'

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

xaxis = np.arange(0, nb_time,1) / 3600. / 24

itype = 0
ilat_center = np.where(max_lat_grid_array == 5)[0][0]
ax.plot(xaxis, percentage_cell_at_least_one_revisit_average_over_lon_all_run[irun][ilat_center, itype, :], linewidth = 3, color = 'r', label = str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN} | ' +str( width_cell_array[itype]) + ' km', linestyle = 'dashed')
ilat_center = np.where(max_lat_grid_array == 30)[0][0]
ax.plot(xaxis, percentage_cell_at_least_one_revisit_average_over_lon_all_run[irun][ilat_center, itype, :], linewidth = 3, color = 'r', label = str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN} | ' +str( width_cell_array[itype]) + ' km')

itype = -1
ilat_center = np.where(max_lat_grid_array == 5)[0][0]
ax.plot(xaxis, percentage_cell_at_least_one_revisit_average_over_lon_all_run[irun][ilat_center, itype, :], linewidth = 3, color = 'b', label = str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN} | ' +str( width_cell_array[itype]) + ' km', linestyle = 'dashed')
ilat_center = np.where(max_lat_grid_array == 30)[0][0]
ax.plot(xaxis, percentage_cell_at_least_one_revisit_average_over_lon_all_run[irun][ilat_center, itype, :], linewidth = 3, color = 'b', label = str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN} | ' +str( width_cell_array[itype]) + ' km')


ax.plot([0, xaxis[-1]], [coverage_array[icoverage], coverage_array[icoverage]], linewidth = 2, linestyle = 'dashed', color = 'k')
ax.minorticks_off()
xticks = np.arange(0,np.max(nb_time) / 3600. / 24., 1)
xticks_label = []
for i in range(len(xticks)):
    xticks_label.append( format( xticks[i], ".0f" ) )
ax.xaxis.set_ticks(xticks)
ax.xaxis.set_ticklabels(xticks_label, fontsize = fontsize_plot)

ax.set_xlim([0, np.max(nb_time) / 3600. / 24.])


legend = ax.legend(loc='lower right', bbox_to_anchor=(1, 0), numpoints = 1,  title="Lat. grid | Cell size", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = plot_folder + 'revisit_coverage/' + run_list[irun]+ '_revisit_coverage_vs_time_cell_width_min_and_max_latitude_min_and_max_TRY.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  




# Coverage and revisit coverage VS time for a given constellation configuration. Min/max latitude and min/amx grid cell shown
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25 
ratio_fig_size = 8./3

icoverage = 2
fig_title = ''
y_label = 'Coverage (%)'
x_label = 'Time (days)'

irun = 0 # constellation to look at

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.99,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 2)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)

ax = fig.add_subplot(gs[0, 0])
#ax.set_title('Coverage VS time for cells 5 km or 25 km wide\nLatitude of grid at 5'+ u'\N{DEGREE SIGN} or 30'+ u'\N{DEGREE SIGN}',  weight = 'bold', fontsize  = fontsize_plot)
ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

xaxis = np.arange(0, nb_time,1) / 3600. / 24

itype = 0
ilat_center = np.where(max_lat_grid_array == 5)[0][0]
ax.plot(xaxis, nb_cell_covered_average_over_lon_all_run[irun][ilat_center, itype, :], linewidth = 3, color = 'r', label = str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN} | ' +str( width_cell_array[itype]) + ' km', linestyle = 'dashed')
ilat_center = np.where(max_lat_grid_array == 30)[0][0]
ax.plot(xaxis, nb_cell_covered_average_over_lon_all_run[irun][ilat_center, itype, :], linewidth = 3, color = 'r', label = str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN} | ' +str( width_cell_array[itype]) + ' km')

itype = -1
ilat_center = np.where(max_lat_grid_array == 5)[0][0]
ax.plot(xaxis, nb_cell_covered_average_over_lon_all_run[irun][ilat_center, itype, :], linewidth = 3, color = 'b', label = str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN} | ' +str( width_cell_array[itype]) + ' km', linestyle = 'dashed')
ilat_center = np.where(max_lat_grid_array == 30)[0][0]
ax.plot(xaxis, nb_cell_covered_average_over_lon_all_run[irun][ilat_center, itype, :], linewidth = 3, color = 'b', label = str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN} | ' +str( width_cell_array[itype]) + ' km')


ax.plot([0, xaxis[-1]], [coverage_array[icoverage], coverage_array[icoverage]], linewidth = 2, linestyle = 'dashed', color = 'k')
ax.minorticks_off()
xticks = np.arange(0,np.max(nb_time) / 3600. / 24., 1)
xticks_label = []
for i in range(len(xticks)):
    xticks_label.append( format( xticks[i], ".0f" ) )
ax.xaxis.set_ticks(xticks)
ax.xaxis.set_ticklabels(xticks_label, fontsize = fontsize_plot)

ax.set_xlim([0, np.max(nb_time) / 3600. / 24.])


ax_revisit = fig.add_subplot(gs[0, 1])

y_label = 'Revisit coverage (%)'

ax_revisit.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax_revisit.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
#ax_revisit.set_title('Revisit coverage VS time for cells 5 km or 25 km wide\nLatitude of grid at 5'+ u'\N{DEGREE SIGN} or 30'+ u'\N{DEGREE SIGN}', weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax_revisit.spines.itervalues()] # change the width of the frame of the figure
ax_revisit.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

xaxis = np.arange(0, nb_time,1) / 3600. / 24

itype = 0
ilat_center = np.where(max_lat_grid_array == 5)[0][0]
ax_revisit.plot(xaxis, percentage_cell_at_least_one_revisit_average_over_lon_all_run[irun][ilat_center, itype, :], linewidth = 3, color = 'r', label = str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN} | ' +str( width_cell_array[itype]) + ' km', linestyle = 'dashed')
ilat_center = np.where(max_lat_grid_array == 30)[0][0]
ax_revisit.plot(xaxis, percentage_cell_at_least_one_revisit_average_over_lon_all_run[irun][ilat_center, itype, :], linewidth = 3, color = 'r', label = str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN} | ' +str( width_cell_array[itype]) + ' km')

itype = -1
ilat_center = np.where(max_lat_grid_array == 5)[0][0]
ax_revisit.plot(xaxis, percentage_cell_at_least_one_revisit_average_over_lon_all_run[irun][ilat_center, itype, :], linewidth = 3, color = 'b', label = str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN} | ' +str( width_cell_array[itype]) + ' km', linestyle = 'dashed')
ilat_center = np.where(max_lat_grid_array == 30)[0][0]
ax_revisit.plot(xaxis, percentage_cell_at_least_one_revisit_average_over_lon_all_run[irun][ilat_center, itype, :], linewidth = 3, color = 'b', label = str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN} | ' +str( width_cell_array[itype]) + ' km')


ax_revisit.plot([0, xaxis[-1]], [coverage_array[icoverage], coverage_array[icoverage]], linewidth = 2, linestyle = 'dashed', color = 'k')
ax_revisit.minorticks_off()
xticks = np.arange(0,np.max(nb_time) / 3600. / 24., 1)
xticks_label = []
for i in range(len(xticks)):
    xticks_label.append( format( xticks[i], ".0f" ) )
ax_revisit.xaxis.set_ticks(xticks)
ax_revisit.xaxis.set_ticklabels(xticks_label, fontsize = fontsize_plot)

ax_revisit.set_xlim([0, np.max(nb_time) / 3600. / 24.])


legend = ax_revisit.legend(loc='lower right', bbox_to_anchor=(1, 0), numpoints = 1,  title="Lat. grid | Cell size", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = plot_folder + 'coverage_vs_time/' + run_list[irun]+ 'coverage_and_revisit_vs_time_cell_width_min_and_max_latitude_min_and_max.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  






# Coverage of cell with at lest one revisit VS time for a given type cell and latitude of grid. All runs plotted
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25 
ratio_fig_size = 4./3


ilat_center = np.where(max_lat_grid_array == 30)[0][0]#2                                                                                                                                                                                                
itype = 0
fig_title = 'Revisit coverage VS time for different constellations\nLatitude of grid ' + str(max_lat_grid_array[ilat_center])+ u'\N{DEGREE SIGN}' + '  - Width of cell ' + str(width_cell_array[itype]) + ' km'
y_label = 'Coverage (%)'
x_label = 'Time (days)'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.99,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

xaxis = np.arange(0, nb_time,1) / 3600. / 24
for irun in range(nb_run):# !!!!!!!! range(1,nb_run):
    ax.plot(xaxis, percentage_cell_at_least_one_revisit_average_over_lon_all_run[irun][ilat_center, itype, :], linewidth = 3, color = color_array[irun], label = str(run_list[irun]))
#ax.plot([0, xaxis[-1]], [coverage_array[icoverage], coverage_array[icoverage]], linewidth = 2, linestyle = 'dashed'                                  )
ax.minorticks_off()
xticks = np.arange(0,np.max(nb_time) / 3600. / 24., 1)
xticks_label = []
for i in range(len(xticks)):
    xticks_label.append( format( xticks[i], ".0f" ) )
ax.xaxis.set_ticks(xticks)
ax.xaxis.set_ticklabels(xticks_label, fontsize = fontsize_plot)

ax.set_xlim([0, np.max(nb_time) / 3600. / 24.])


legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="2nd plane inc.", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = plot_folder + 'revisit_coverage/' +  root_fig + '_compare_revisit_coverage_vs_time_cellwidth_' + str(width_cell_array[itype]).replace(".","_") + '_latitude_' + str(max_lat_grid_array[ilat_center]) + '.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  






# Mean revisit time for one constellation configuration VS latitude of grid. All types of grid cells on one plot (diffeernt colors)
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25 
ratio_fig_size = 4./3
fig_title = ''#Mean revisit time VS latitude of grid for different cell sizes'
y_label = 'Mean revisit time (hours)'
x_label = 'Latitude of grid (' + u'\N{DEGREE SIGN})'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.96,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold


irun = 0 # configuration to look at
mean_revisit_time_average_over_lon_all_run.append(pickle.load(open(pickle_folder + run_list[irun] + "_mean_revisit_time_average_over_lon_gain_not0.pickle"))) 
xaxis = max_lat_grid_array
for itype in range(nb_type_cell):
    ax.plot(xaxis, mean_revisit_time_average_over_lon_all_run[irun][:, itype]/3600., linewidth = 3, color = color_array[itype], label =  str(width_cell_array[itype])+ ' km') 


legend = ax.legend(loc='lower right', bbox_to_anchor=(1, 1), numpoints = 1,  title="Resolution of cell", fontsize = fontsize_plot, ncol = 5, columnspacing = 0.7)
legend.get_title().set_fontsize(str(fontsize_plot))

#ax.set_ylim([8,100])

fig_save_name = plot_folder + 'revisit_time/' + run_list[irun] + '_revisit_time_vs_lat_all_grid_cell.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    

# Mean revisit time for a given type cell VS latitude of grid. All runs plotted
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25 
ratio_fig_size = 4./3

for irun in range(nb_run):
    mean_revisit_time_average_over_lon_all_run.append(pickle.load(open(pickle_folder + run_list[irun] + "_mean_revisit_time_average_over_lon_gain_not0.pickle"))) 

#for itype in range(nb_type_cell):#
itype = 0
fig_title = ''#Mean revisit time VS latitude of grid for different constellations\nWidth of cell ' + str(width_cell_array[itype]) + ' km'
y_label = 'Mean revisit time (hours)'
x_label = 'Latitude of grid (' + u'\N{DEGREE SIGN})'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.99,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

xaxis = max_lat_grid_array
for irun in range(nb_run):# !!!!!!!! range(1,nb_run):
    ax.plot(xaxis, mean_revisit_time_average_over_lon_all_run[irun][:, itype]/3600., linewidth = 3, color = color_array[irun], label = label_arr[irun]) 


legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = plot_folder + 'revisit_time/' + root_fig + '_compare_revisit_time_vs_lat_grid_cell_width_' + str(width_cell_array[itype]).replace(".","_") + 'gain_not0.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


#XXXXXXXXXXXXXXXXXXXXXXXXXXXX
# 2D visualization
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25 
ratio_fig_size = 4./3

irun = 1
ilat_center = np.where(max_lat_grid_array == 5)[0][0]#2
ilon_center = 0 # has to be 0 
rlat_width = 500 # latitude width of the grid, in km. !!!!!! need to be the same as in cygnss_coverage_temporal_spatial_relation_equal_cell_area2.py
rlon_width = 500 # longitude width of the grid, in km  !!!!!! need to be the same as in cygnss_coverage_temporal_spatial_relation_equal_cell_area2.py
re = 6371.0 # mean Earth radius !!!!!! need to be the same as in cygnss_coverage_temporal_spatial_relation_equal_cell_area2.py


#itype = 0 # grid bin type. need to be 0 or -1

width_cell = []
nb_cell_lat = []
nb_cell_lon = []
nb_cell_in_grid = []
icell_max_lat_grid = []
icell_min_lat_grid = []
width_cell_degree_lat = []
max_lat_grid = []
min_lat_grid = []
width_cell_degree_lon = []
icell_min_lon_grid = []
icell_max_lon_grid = []
min_lon_grid = []
max_lon_grid = []
for itype in [0,-1]:
    width_cell.append( width_cell_array[itype] )

    nb_cell_lat.append( (int) (rlat_width / width_cell[-1]) ) # nb cell bins in one grid in the latitu=dinal direction )
    nb_cell_lon.append( (int) (rlon_width / width_cell[-1]) ) # nb cell bins in one grid in the longitudinal direction. )
    nb_cell_in_grid.append( nb_cell_lon[-1] * nb_cell_lat[-1] )


    icell_max_lat_grid.append( (int) ( re * max_lat_grid_array[ilat_center]*deg2rad  / width_cell[-1] ) ) # cell# that corresponds to the top of the grid
    icell_min_lat_grid.append( icell_max_lat_grid[-1] - nb_cell_lat[-1] )  # cell# that corresponds to the bottom of the grid
    width_cell_degree_lat.append( width_cell[-1] / re /deg2rad )
    max_lat_grid.append( ( icell_max_lat_grid[-1]  ) * width_cell_degree_lat[-1] )
    min_lat_grid.append( icell_min_lat_grid[-1] * width_cell_degree_lat[-1] )
 
    width_cell_degree_lon.append( width_cell[-1] / ( re * np.cos(max_lat_grid_array[ilat_center]*deg2rad) ) / deg2rad )
    icell_min_lon_grid.append( ilon_center*nb_cell_lon[-1]  )
    icell_max_lon_grid.append( ilon_center*nb_cell_lon[-1] + nb_cell_lon[-1] )
    min_lon_grid.append( icell_min_lon_grid[-1] * width_cell_degree_lon[-1] )
    max_lon_grid.append( ( icell_max_lon_grid[-1]  ) * width_cell_degree_lon[-1] )


# load pickle
itype = 0
cell_array_type_lon_lat_cumul_small_cell = pickle.load(open(pickle_folder + pickle_visu_subfolder + run_list[irun]  + "_cell_array_type_lon_lat_cumul_ilat_" + str((int)(max_lat_grid_array[ilat_center]*deg2rad/deg2rad)) + '_cell_' + str(width_cell_array[itype])   + "_gain_not0.pickle")) 
itype = -1
cell_array_type_lon_lat_cumul_big_cell = pickle.load(open(pickle_folder + pickle_visu_subfolder + run_list[irun]  + "_cell_array_type_lon_lat_cumul_ilat_" + str((int)(max_lat_grid_array[ilat_center]*deg2rad/deg2rad)) + '_cell_' + str(width_cell_array[itype])   + "_gain_not0.pickle")) 

lon_spec = pickle.load(open(pickle_folder + pickle_visu_subfolder +  run_list[irun] + "_save_lon_spec_gain_not0.pickle")) /deg2rad
lat_spec = pickle.load(open(pickle_folder + pickle_visu_subfolder + run_list[irun] + "_save_lat_spec_gain_not0.pickle")) /deg2rad

nb_day_visu_with_spec_plotted = 2# beyond this limit, the path of spec wont be pltted on the animation, just the grid cells (which are recorded until the end of SpOCK run)
#!!!!!! need to be the same as in cygnss_coverage_temporal_spatial_relation_equal_cell_area2.py 
nb_second_visu = nb_day_visu_with_spec_plotted * 24*3600
#ipdb.set_trace()
dt_visu_hour = 1 # time step for the visu in hour #!!!!!! need to be the same as in cygnss_coverage_temporal_spatial_relation_equal_cell_area2.py 
dt_visu = dt_visu_hour * 3600

nb_time  = (int) (cell_array_type_lon_lat_cumul_small_cell.shape[0] * dt_visu)


y_label = 'Latitude '+ u'(\N{DEGREE SIGN})'
x_label = 'Longitude '+ u'(\N{DEGREE SIGN})'
cmap = colors.ListedColormap(['azure'])#, 'lightblue','blue' ])#'lightskyblue'])
#handles_arr = [mpatches.Patch(color='lightblue', label='Visited'), mpatches.Patch(color='blue', label='Revisited')]
handles_arr = [mpatches.Patch(color='lightblue', label='')]
nb_sc = 8

itime_count = -1
for itime in np.arange(3600, 2*24*3600, 3600):#nb_time, dt_visu):
    itime_count = itime_count + 1
    print '\n',itime,  nb_time
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    gs = gridspec.GridSpec(2, 1)
    gs.update(left = 0.11, right=0.87, top = 1,bottom = 0.12, hspace = 0.15)

    time_since_start_day =  timedelta(seconds = itime).days#(int) ( np.floor(itime / 3600. / 24))
    time_since_start_hour, remainder =  divmod(timedelta(seconds = itime).seconds, 3600)#(int) (  np.floor(  (itime / 3600. / 24 - time_since_start_day )*24) )
    time_since_start_min, time_since_start_sec = divmod(remainder, 60)#(int) ( np.floor(( (itime / 3600. / 24 - time_since_start_day )*24 - time_since_start_hour ) * 60) )
    #time_since_start_sec = (int) ( np.round( ( ( (itime / 3600. / 24 - time_since_start_day )*24 - time_since_start_hour ) * 60 - time_since_start_min ) * 60 ) )
    time_since_start = str(time_since_start_day) + 'd ' + str(time_since_start_hour) + 'h' + str(time_since_start_min) + "'"  + str(time_since_start_sec) + "''"
    print (int) (itime / dt_visu)
    fig.suptitle('', y = 0.96,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold


    # on top plot the data for small size cells
    ax_2d_small = fig.add_subplot(gs[0, 0])
    # draw grid
    data = np.zeros([nb_cell_lat[0], nb_cell_lon[0]])

    data[:,:] = cell_array_type_lon_lat_cumul_small_cell[(int) (itime / dt_visu), : ,:]
    percentage_covered = len(np.where(data >= 1)[0])*100./(nb_cell_lon[0]*nb_cell_lat[0])
    percentage_revisit_covered = len(np.where(data >= 2)[0])*100./(nb_cell_lon[0]*nb_cell_lat[0])

    ax_2d_small_title = time_since_start + '\n' + str(width_cell[0])+ " km cell - Coverage: " + format(percentage_covered, ".0f") + "% (revisited: " + format(percentage_revisit_covered, ".0f") + "%)"


    [i.set_linewidth(2) for i in ax_2d_small.spines.itervalues()] # change the width of the frame of the figure
    ax_2d_small.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

    lon_min_grid_converted = min_lon_grid[0] 
    lon_max_grid_converted = max_lon_grid[0] 

    lon_arr = np.zeros([nb_cell_lon[0]+1])
    for ilon in range(nb_cell_lon[0]+1):
        lon_arr[ilon]  = min_lon_grid[0] + ilon * width_cell_degree_lon[0]
    lat_arr = np.zeros([nb_cell_lat[0]+1])
    for ilat in range(nb_cell_lat[0]+1):
        lat_arr[ilat]  = min_lat_grid[0] + ilat * width_cell_degree_lat[0]


    m = Basemap( projection       = 'cyl',
             llcrnrlon        =  np.min(lon_arr), #Lower Left  CoRNeR Longitude
             urcrnrlon        =  np.max(lon_arr), #Upper Right CoRNeR Longitude
             llcrnrlat        = np.min(lat_arr)  , #Lower Left  CoRNeR Latitude
             urcrnrlat        = np.max(lat_arr),   #Upper Right CoRNeR Latitude
             resolution       = 'i'  , # 'i', 'h'
             suppress_ticks   = False,
             ax = ax_2d_small,
             )
    # !!!!!!!!! the way the mesh is made here is bad because it assumes that a cell has the same longitudinal dimension no matter its latittude. This is not true, particularly for high latituedes. Ok if latitudes are below about 40 degrees. But if ther is a mismatch between the cells and the spec positions, this is the reason.
    m.pcolormesh(lon_arr, lat_arr, data, cmap=cmap, linewidth = 0.01,  edgecolors = 'lightgrey')  # !!!!!! be careful: X and Y need to be the number of element of data + 1 (), otherwise one cell will be missing (https://matplotlib.org/basemap/api/basemap_api.html#mpl_toolkits.basemap.Basemap.pcolormesh)

    ax_2d_small.plot(0,0, color = 'lightblue', linewidth=2)

    specular_list = []
    point = namedtuple('point', ['x', 'y'])
    color = namedtuple('color', 'red green blue')
    specular = namedtuple('specular',('name',) +  point._fields  + color._fields + ('point_plot',))

    if itime < nb_second_visu: # the spec are shown only for the firset nb_second_visu seconds of the simulation (this is because we don't want to have too heavy pickles for save_lon_spec and save_lat_spec)
        for isc in range(nb_sc):
            # Add on plot the specular points over one orbit
            for k in range(nb_spec):
                specular_list.append(specular)
                if itime > 3600:
                    specular_list[k+isc*nb_spec].x, specular_list[k+isc*nb_spec].y =  m(lon_spec[k,isc,itime-3600:itime], lat_spec[k,isc,itime-3600:itime])
                else:
                    specular_list[k+isc*nb_spec].x, specular_list[k+isc*nb_spec].y =  m(lon_spec[k,isc,:itime], lat_spec[k,isc,:itime])
                specular_list[k+isc*nb_spec].point_plot = m.scatter(specular_list[k+isc*nb_spec].x, specular_list[k+isc*nb_spec].y, marker='o', s = 10, color = 'b')#50*np.cos(max_lat_grid_array[ilat_center]*deg2rad) * width_cell[0] / np.max(width_cell_array) , color = 'b')


    m.drawcoastlines(linewidth=0.7, color='blue')

    ax_2d_small.set_title(ax_2d_small_title, weight = 'bold', fontsize  = (int)(fontsize_plot*1.1), y = 1.0)
    ax_2d_small.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax_2d_small.xaxis.set_ticklabels('')



    # on bottom plot the data for big size cells
    ax_2d_big = fig.add_subplot(gs[1, 0])
    # draw grid
    data = np.zeros([nb_cell_lat[1], nb_cell_lon[1]])

    data[:,:] = cell_array_type_lon_lat_cumul_big_cell[(int) (itime / dt_visu), : ,:]
    percentage_covered = len(np.where(data >= 1)[0])*100./(nb_cell_lon[1]*nb_cell_lat[1])
    percentage_revisit_covered = len(np.where(data >= 2)[0])*100./(nb_cell_lon[1]*nb_cell_lat[1])

    ax_2d_big_title = str(width_cell[1])+ " km cell - Coverage: " + format(percentage_covered, ".0f") + "% (revisited: " + format(percentage_revisit_covered, ".0f") + "%)"


    [i.set_linewidth(2) for i in ax_2d_big.spines.itervalues()] # change the width of the frame of the figure
    ax_2d_big.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

    lon_min_grid_converted = min_lon_grid[1] 
    lon_max_grid_converted = max_lon_grid[1] 

    lon_arr = np.zeros([nb_cell_lon[1]+1])
    for ilon in range(nb_cell_lon[1]+1):
        lon_arr[ilon]  = min_lon_grid[1] + ilon * width_cell_degree_lon[1]
    lat_arr = np.zeros([nb_cell_lat[1]+1])
    for ilat in range(nb_cell_lat[1]+1):
        lat_arr[ilat]  = min_lat_grid[1] + ilat * width_cell_degree_lat[1]


    m = Basemap( projection       = 'cyl',
             llcrnrlon        =  np.min(lon_arr), #Lower Left  CoRNeR Longitude
             urcrnrlon        =  np.max(lon_arr), #Upper Right CoRNeR Longitude
             llcrnrlat        = np.min(lat_arr)  , #Lower Left  CoRNeR Latitude
             urcrnrlat        = np.max(lat_arr),   #Upper Right CoRNeR Latitude
             resolution       = 'i'  , # 'i', 'h'
             suppress_ticks   = False,
             ax = ax_2d_big,
             )
    # !!!!!!!!! the way the mesh is made here is bad because it assumes that a cell has the same longitudinal dimension no matter its latittude. This is not true, particularly for high latituedes. Ok if latitudes are below about 40 degrees. But if ther is a mismatch between the cells and the spec positions, this is the reason.
    m.pcolormesh(lon_arr, lat_arr, data, cmap=cmap, linewidth = 0.01,  edgecolors = 'lightgrey')  # !!!!!! be careful: X and Y need to be the number of element of data + 1 (), otherwise one cell will be missing (https://matplotlib.org/basemap/api/basemap_api.html#mpl_toolkits.basemap.Basemap.pcolormesh)

    ax_2d_big.plot(0,0, color = 'lightblue', linewidth=2)

    specular_list = []
    point = namedtuple('point', ['x', 'y'])
    color = namedtuple('color', 'red green blue')
    specular = namedtuple('specular',('name',) +  point._fields  + color._fields + ('point_plot',))

    if itime < nb_second_visu: # the spec are shown only for the firset nb_second_visu seconds of the simulation (this is because we don't want to have too heavy pickles for save_lon_spec and save_lat_spec)
        for isc in range(nb_sc):
            # Add on plot the specular points over one orbit
            for k in range(nb_spec):
                specular_list.append(specular)
                if itime > 3600:
                    specular_list[k+isc*nb_spec].x, specular_list[k+isc*nb_spec].y =  m(lon_spec[k,isc,itime-3600:itime], lat_spec[k,isc,itime-3600:itime])
                else:
                    specular_list[k+isc*nb_spec].x, specular_list[k+isc*nb_spec].y =  m(lon_spec[k,isc,:itime], lat_spec[k,isc,:itime])
                specular_list[k+isc*nb_spec].point_plot = m.scatter(specular_list[k+isc*nb_spec].x, specular_list[k+isc*nb_spec].y, marker='o', s = 10, color = 'b')#50*np.cos(max_lat_grid_array[ilat_center]*deg2rad) * width_cell[1] / np.max(width_cell_array) , color = 'b')


    m.drawcoastlines(linewidth=0.7, color='blue')

    ax_2d_big.set_title(ax_2d_big_title, weight = 'bold', fontsize  = (int)(fontsize_plot*1.1), y = 1.0)
    ax_2d_big.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax_2d_big.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot, zorder = 5)

    legend = ax_2d_big.legend( loc='center',  bbox_to_anchor=(0.5, -0.144), fontsize = fontsize_plot, handles=handles_arr, ncol=2, columnspacing = height_fig * ratio_fig_size*1.1,frameon=False)#    legend = ax_2d_big.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="", fontsize = fontsize_plot, handles=handles_arr)




    fig_save_name = plot_folder + plot_visu_subfolder + run_list[irun] + '_lat_' + str(max_lat_grid_array[ilat_center])  + '_2d_coverage_itime_' + str(itime_count) + '.png'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')


os.system('ffmpeg -y -r 2 -i ' + plot_folder + plot_visu_subfolder + run_list[irun]  + '_lat_' + str(max_lat_grid_array[ilat_center]) + '_2d_coverage_itime_%d.png -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -vcodec libx264 -pix_fmt yuv420p ' + plot_folder + plot_visu_subfolder + 'ani/ani_' + run_list[irun]  + '_lat_' + str(max_lat_grid_array[ilat_center]) + '.mp4')



#XXXXXXXXXXXXXXXXXXXXXXXXXXXX
# 2D visualization: 2 constellation configurations on same animations. One type of cells
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25            
ratio_fig_size = 4./3


ilat_center = np.where(max_lat_grid_array == 5)[0][0]#2
ilon_center = 0 # has to be 0 
rlat_width = 500 # latitude width of the grid, in km. !!!!!! need to be the same as in cygnss_coverage_temporal_spatial_relation_equal_cell_area2.py
rlon_width = 500 # longitude width of the grid, in km  !!!!!! need to be the same as in cygnss_coverage_temporal_spatial_relation_equal_cell_area2.py
re = 6371.0 # mean Earth radius !!!!!! need to be the same as in cygnss_coverage_temporal_spatial_relation_equal_cell_area2.py


itype = 1

width_cell = []
nb_cell_lat = []
nb_cell_lon = []
nb_cell_in_grid = []
icell_max_lat_grid = []
icell_min_lat_grid = []
width_cell_degree_lat = []
max_lat_grid = []
min_lat_grid = []
width_cell_degree_lon = []
icell_min_lon_grid = []
icell_max_lon_grid = []
min_lon_grid = []
max_lon_grid = []

width_cell.append( width_cell_array[itype] )

nb_cell_lat.append( (int) (rlat_width / width_cell[-1]) ) # nb cell bins in one grid in the latitu=dinal direction )
nb_cell_lon.append( (int) (rlon_width / width_cell[-1]) ) # nb cell bins in one grid in the longitudinal direction. )
nb_cell_in_grid.append( nb_cell_lon[-1] * nb_cell_lat[-1] )


icell_max_lat_grid.append( (int) ( re * max_lat_grid_array[ilat_center]*deg2rad  / width_cell[-1] ) ) # cell# that corresponds to the top of the grid
icell_min_lat_grid.append( icell_max_lat_grid[-1] - nb_cell_lat[-1] )  # cell# that corresponds to the bottom of the grid
width_cell_degree_lat.append( width_cell[-1] / re /deg2rad )
max_lat_grid.append( ( icell_max_lat_grid[-1]  ) * width_cell_degree_lat[-1] )
min_lat_grid.append( icell_min_lat_grid[-1] * width_cell_degree_lat[-1] )

width_cell_degree_lon.append( width_cell[-1] / ( re * np.cos(max_lat_grid_array[ilat_center]*deg2rad) ) / deg2rad )
icell_min_lon_grid.append( ilon_center*nb_cell_lon[-1]  )
icell_max_lon_grid.append( ilon_center*nb_cell_lon[-1] + nb_cell_lon[-1] )
min_lon_grid.append( icell_min_lon_grid[-1] * width_cell_degree_lon[-1] )
max_lon_grid.append( ( icell_max_lon_grid[-1]  ) * width_cell_degree_lon[-1] )


# load pickle
irun = 0
cell_array_type_lon_lat_cumul_config1 = pickle.load(open(pickle_folder + pickle_visu_subfolder + run_list[irun]  + "_cell_array_type_lon_lat_cumul_ilat_" + str((int)(max_lat_grid_array[ilat_center]*deg2rad/deg2rad)) + '_cell_' + str(width_cell_array[itype])   + "_gain_not0.pickle")) 
lon_spec_config1 = pickle.load(open(pickle_folder + pickle_visu_subfolder +  run_list[irun] + "_save_lon_spec_gain_not0.pickle")) /deg2rad
lat_spec_config1 = pickle.load(open(pickle_folder + pickle_visu_subfolder + run_list[irun] + "_save_lat_spec_gain_not0.pickle")) /deg2rad

irun = 1
cell_array_type_lon_lat_cumul_config2 = pickle.load(open(pickle_folder + pickle_visu_subfolder + run_list[irun]  + "_cell_array_type_lon_lat_cumul_ilat_" + str((int)(max_lat_grid_array[ilat_center]*deg2rad/deg2rad)) + '_cell_' + str(width_cell_array[itype])   + "_gain_not0.pickle")) 
lon_spec_config2 = pickle.load(open(pickle_folder + pickle_visu_subfolder +  run_list[irun] + "_save_lon_spec_gain_not0.pickle")) /deg2rad
lat_spec_config2 = pickle.load(open(pickle_folder + pickle_visu_subfolder + run_list[irun] + "_save_lat_spec_gain_not0.pickle")) /deg2rad

nb_day_visu_with_spec_plotted = 2# beyond this limit, the path of spec wont be pltted on the animation, just the grid cells (which are recorded until the end of SpOCK run)
#!!!!!! need to be the same as in cygnss_coverage_temporal_spatial_relation_equal_cell_area2.py 
nb_second_visu = nb_day_visu_with_spec_plotted * 24*3600
#ipdb.set_trace()
dt_visu_hour = 1 # time step for the visu in hour #!!!!!! need to be the same as in cygnss_coverage_temporal_spatial_relation_equal_cell_area2.py 
dt_visu = dt_visu_hour * 3600

nb_time  = (int) (cell_array_type_lon_lat_cumul_config1.shape[0] * dt_visu)


y_label = 'Latitude '+ u'(\N{DEGREE SIGN})'
x_label = 'Longitude '+ u'(\N{DEGREE SIGN})'
cmap = colors.ListedColormap(['azure'])#, 'lightblue','blue' ])#'lightskyblue'])
#handles_arr = [mpatches.Patch(color='lightblue', label='Visited'), mpatches.Patch(color='blue', label='Revisited')]
handles_arr = [mpatches.Patch(color='lightblue', label='')]
nb_sc = 8

itime_count = -1
for itime in np.arange(3600, 2*24*3600, 3600):#nb_time, dt_visu):
    itime_count = itime_count + 1
    print '\n',itime,  nb_time
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    gs = gridspec.GridSpec(2, 1)
    gs.update(left = 0.11, right=0.87, top = 1,bottom = 0.12, hspace = 0.2)

    time_since_start_day =  timedelta(seconds = itime).days#(int) ( np.floor(itime / 3600. / 24))
    time_since_start_hour, remainder =  divmod(timedelta(seconds = itime).seconds, 3600)#(int) (  np.floor(  (itime / 3600. / 24 - time_since_start_day )*24) )
    time_since_start_min, time_since_start_sec = divmod(remainder, 60)#(int) ( np.floor(( (itime / 3600. / 24 - time_since_start_day )*24 - time_since_start_hour ) * 60) )
    #time_since_start_sec = (int) ( np.round( ( ( (itime / 3600. / 24 - time_since_start_day )*24 - time_since_start_hour ) * 60 - time_since_start_min ) * 60 ) )
    time_since_start = str(time_since_start_day) + 'd ' + str(time_since_start_hour) + 'h' + str(time_since_start_min) + "'"  + str(time_since_start_sec) + "''"
    print (int) (itime / dt_visu)
    fig.suptitle('', y = 0.96,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold


    # on top plot the data for small size cells
    ax_2d_small = fig.add_subplot(gs[0, 0])
    # draw grid
    data = np.zeros([nb_cell_lat[0], nb_cell_lon[0]])

    data[:,:] = cell_array_type_lon_lat_cumul_config1[(int) (itime / dt_visu), : ,:]
    percentage_covered = len(np.where(data >= 1)[0])*100./(nb_cell_lon[0]*nb_cell_lat[0])
    percentage_revisit_covered = len(np.where(data >= 2)[0])*100./(nb_cell_lon[0]*nb_cell_lat[0])

    ax_2d_small_title = time_since_start + '\n' + str(width_cell[0])+ " km cell - Coverage: " + format(percentage_covered, ".0f") + "% (revisited: " + format(percentage_revisit_covered, ".0f") + "%)\n" + label_arr[0]


    [i.set_linewidth(2) for i in ax_2d_small.spines.itervalues()] # change the width of the frame of the figure
    ax_2d_small.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

    lon_min_grid_converted = min_lon_grid[0] 
    lon_max_grid_converted = max_lon_grid[0] 

    lon_arr = np.zeros([nb_cell_lon[0]+1])
    for ilon in range(nb_cell_lon[0]+1):
        lon_arr[ilon]  = min_lon_grid[0] + ilon * width_cell_degree_lon[0]
    lat_arr = np.zeros([nb_cell_lat[0]+1])
    for ilat in range(nb_cell_lat[0]+1):
        lat_arr[ilat]  = min_lat_grid[0] + ilat * width_cell_degree_lat[0]


    m = Basemap( projection       = 'cyl',
             llcrnrlon        =  np.min(lon_arr), #Lower Left  CoRNeR Longitude
             urcrnrlon        =  np.max(lon_arr), #Upper Right CoRNeR Longitude
             llcrnrlat        = np.min(lat_arr)  , #Lower Left  CoRNeR Latitude
             urcrnrlat        = np.max(lat_arr),   #Upper Right CoRNeR Latitude
             resolution       = 'i'  , # 'i', 'h'
             suppress_ticks   = False,
             ax = ax_2d_small,
             )
    # !!!!!!!!! the way the mesh is made here is bad because it assumes that a cell has the same longitudinal dimension no matter its latittude. This is not true, particularly for high latituedes. Ok if latitudes are below about 40 degrees. But if ther is a mismatch between the cells and the spec positions, this is the reason.
    m.pcolormesh(lon_arr, lat_arr, data, cmap=cmap, linewidth = 0.01,  edgecolors = 'lightgrey')  # !!!!!! be careful: X and Y need to be the number of element of data + 1 (), otherwise one cell will be missing (https://matplotlib.org/basemap/api/basemap_api.html#mpl_toolkits.basemap.Basemap.pcolormesh)

    ax_2d_small.plot(0,0, color = 'lightblue', linewidth=2)

    specular_list = []
    point = namedtuple('point', ['x', 'y'])
    color = namedtuple('color', 'red green blue')
    specular = namedtuple('specular',('name',) +  point._fields  + color._fields + ('point_plot',))

    if itime < nb_second_visu: # the spec are shown only for the firset nb_second_visu seconds of the simulation (this is because we don't want to have too heavy pickles for save_lon_spec and save_lat_spec)
        for isc in range(nb_sc):
            # Add on plot the specular points over one orbit
            for k in range(nb_spec):
                specular_list.append(specular)
                if itime > 3600:
                    specular_list[k+isc*nb_spec].x, specular_list[k+isc*nb_spec].y =  m(lon_spec_config1[k,isc,itime-3600:itime], lat_spec_config1[k,isc,itime-3600:itime])
                else:
                    specular_list[k+isc*nb_spec].x, specular_list[k+isc*nb_spec].y =  m(lon_spec_config1[k,isc,:itime], lat_spec_config1[k,isc,:itime])
                specular_list[k+isc*nb_spec].point_plot = m.scatter(specular_list[k+isc*nb_spec].x, specular_list[k+isc*nb_spec].y, marker='o', s = 10, color = 'b')#50*np.cos(max_lat_grid_array[ilat_center]*deg2rad) * width_cell[0] / np.max(width_cell_array) , color = 'b')


    m.drawcoastlines(linewidth=0.7, color='blue')

    ax_2d_small.set_title(ax_2d_small_title, weight = 'bold', fontsize  = (int)(fontsize_plot*1.1), y = 1.0)
    ax_2d_small.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax_2d_small.xaxis.set_ticklabels('')



    # on bottom plot the data for big size cells
    ax_2d_big = fig.add_subplot(gs[1, 0])
    # draw grid
    data = np.zeros([nb_cell_lat[0], nb_cell_lon[0]])

    data[:,:] = cell_array_type_lon_lat_cumul_config2[(int) (itime / dt_visu), : ,:]
    percentage_covered = len(np.where(data >= 1)[0])*100./(nb_cell_lon[0]*nb_cell_lat[0])
    percentage_revisit_covered = len(np.where(data >= 2)[0])*100./(nb_cell_lon[0]*nb_cell_lat[0])

    ax_2d_big_title = str(width_cell[0])+ " km cell - Coverage: " + format(percentage_covered, ".0f") + "% (revisited: " + format(percentage_revisit_covered, ".0f") + "%)\n" + label_arr[1]


    [i.set_linewidth(2) for i in ax_2d_big.spines.itervalues()] # change the width of the frame of the figure
    ax_2d_big.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

    lon_min_grid_converted = min_lon_grid[0] 
    lon_max_grid_converted = max_lon_grid[0] 

    lon_arr = np.zeros([nb_cell_lon[0]+1])
    for ilon in range(nb_cell_lon[0]+1):
        lon_arr[ilon]  = min_lon_grid[0] + ilon * width_cell_degree_lon[0]
    lat_arr = np.zeros([nb_cell_lat[0]+1])
    for ilat in range(nb_cell_lat[0]+1):
        lat_arr[ilat]  = min_lat_grid[0] + ilat * width_cell_degree_lat[0]


    m = Basemap( projection       = 'cyl',
             llcrnrlon        =  np.min(lon_arr), #Lower Left  CoRNeR Longitude
             urcrnrlon        =  np.max(lon_arr), #Upper Right CoRNeR Longitude
             llcrnrlat        = np.min(lat_arr)  , #Lower Left  CoRNeR Latitude
             urcrnrlat        = np.max(lat_arr),   #Upper Right CoRNeR Latitude
             resolution       = 'i'  , # 'i', 'h'
             suppress_ticks   = False,
             ax = ax_2d_big,
             )
    # !!!!!!!!! the way the mesh is made here is bad because it assumes that a cell has the same longitudinal dimension no matter its latittude. This is not true, particularly for high latituedes. Ok if latitudes are below about 40 degrees. But if ther is a mismatch between the cells and the spec positions, this is the reason.
    m.pcolormesh(lon_arr, lat_arr, data, cmap=cmap, linewidth = 0.01,  edgecolors = 'lightgrey')  # !!!!!! be careful: X and Y need to be the number of element of data + 1 (), otherwise one cell will be missing (https://matplotlib.org/basemap/api/basemap_api.html#mpl_toolkits.basemap.Basemap.pcolormesh)

    ax_2d_big.plot(0,0, color = 'lightblue', linewidth=2)

    specular_list = []
    point = namedtuple('point', ['x', 'y'])
    color = namedtuple('color', 'red green blue')
    specular = namedtuple('specular',('name',) +  point._fields  + color._fields + ('point_plot',))

    if itime < nb_second_visu: # the spec are shown only for the firset nb_second_visu seconds of the simulation (this is because we don't want to have too heavy pickles for save_lon_spec and save_lat_spec)
        for isc in range(nb_sc):
            # Add on plot the specular points over one orbit
            for k in range(nb_spec):
                specular_list.append(specular)
                if itime > 3600:
                    specular_list[k+isc*nb_spec].x, specular_list[k+isc*nb_spec].y =  m(lon_spec_config2[k,isc,itime-3600:itime], lat_spec_config2[k,isc,itime-3600:itime])
                else:
                    specular_list[k+isc*nb_spec].x, specular_list[k+isc*nb_spec].y =  m(lon_spec_config2[k,isc,:itime], lat_spec_config2[k,isc,:itime])
                specular_list[k+isc*nb_spec].point_plot = m.scatter(specular_list[k+isc*nb_spec].x, specular_list[k+isc*nb_spec].y, marker='o', s = 10, color = 'b')#50*np.cos(max_lat_grid_array[ilat_center]*deg2rad) * width_cell[0] / np.max(width_cell_array) , color = 'b')


    m.drawcoastlines(linewidth=0.7, color='blue')

    ax_2d_big.set_title(ax_2d_big_title, weight = 'bold', fontsize  = (int)(fontsize_plot*1.1), y = 1.0)
    ax_2d_big.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax_2d_big.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot, zorder = 5)

#    legend = ax_2d_big.legend( loc='center',  bbox_to_anchor=(0.5, -0.144), fontsize = fontsize_plot, handles=handles_arr, ncol=2, columnspacing = height_fig * ratio_fig_size*1.1,frameon=False)#    legend = ax_2d_big.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="", fontsize = fontsize_plot, handles=handles_arr)




    fig_save_name = plot_folder + plot_visu_subfolder + run_list[0] + '_and_' + run_list[1] +  '_lat_' + str(max_lat_grid_array[ilat_center])  + '_2d_coverage_itime_' + str(itime_count) + '.png'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')


os.system('ffmpeg -y -r 2 -i ' + plot_folder + plot_visu_subfolder + run_list[0] + '_and_' + run_list[1]  + '_lat_' + str(max_lat_grid_array[ilat_center]) + '_2d_coverage_itime_%d.png -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -vcodec libx264 -pix_fmt yuv420p ' + plot_folder + plot_visu_subfolder + 'ani/ani_' + run_list[0] + '_and_' + run_list[1]  + '_lat_' + str(max_lat_grid_array[ilat_center]) + '.mp4')





# MESSED UP
#     # on bottom plot the data for big size cells
#     ax_2d_big = fig.add_subplot(gs[1, 0])
#     # draw grid
#     data = np.zeros([nb_cell_lat[1], nb_cell_lon[1]])

#     data[:,:] = cell_array_type_lon_lat_cumul_config2[(int) (itime / dt_visu), : ,:]
#     percentage_covered = len(np.where(data >= 1)[0])*100./(nb_cell_lon[1]*nb_cell_lat[1])
#     percentage_revisit_covered = len(np.where(data >= 2)[0])*100./(nb_cell_lon[1]*nb_cell_lat[1])

#     ax_2d_big_title = str(width_cell[1])+ " km cell - Coverage: " + format(percentage_covered, ".0f") + "% (revisited: " + format(percentage_revisit_covered, ".0f") + "%)"


#     [i.set_linewidth(2) for i in ax_2d_big.spines.itervalues()] # change the width of the frame of the figure
#     ax_2d_big.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
#     plt.rc('font', weight='bold') ## make the labels of the ticks in bold

#     lon_min_grid_converted = min_lon_grid[1] 
#     lon_max_grid_converted = max_lon_grid[1] 
    
#     lon_arr = np.zeros([nb_cell_lon[1]+1])
#     for ilon in range(nb_cell_lon[1]+1):
#         lon_arr[ilon]  = min_lon_grid[1] + ilon * width_cell_degree_lon[1]
#     lat_arr = np.zeros([nb_cell_lat[1]+1])
#     for ilat in range(nb_cell_lat[1]+1):
#         lat_arr[ilat]  = min_lat_grid[1] + ilat * width_cell_degree_lat[1]

#     m = Basemap( projection       = 'cea',
#                  llcrnrlon        =  np.min(lon_arr), #Lower Left  CoRNeR Longitude
#                  urcrnrlon        =  np.max(lon_arr), #Upper Right CoRNeR Longitude
#                  llcrnrlat        = np.min(lat_arr)  , #Lower Left  CoRNeR Latitude
#                  urcrnrlat        = np.max(lat_arr),   #Upper Right CoRNeR Latitude
#                  resolution       = 'c'  , # 'c','i', 'h'
#                  suppress_ticks   = False,
#                  ax = ax_2d_big,
#                  )
#     m.drawparallels(np.arange(-90.,91.,30.))
#     m.drawmeridians(np.arange(-180.,181.,60.))
#     m.drawcoastlines()
# #     m.pcolormesh(lon_arr, lat_arr, data, cmap=cmap, linewidth = 0.01,  edgecolors = 'lightgrey', latlon=True) # !!!!!! be careful: X and Y need to be the number of element of data + 1 (), otherwise one cell will be missing (https://matplotlib.org/basemap/api/basemap_api.html#mpl_toolkits.basemap.Basemap.pcolormesh)
#     min_x_arr, min_y_arr = m(np.min(lon_arr), np.min(lat_arr))
#     max_x_arr, max_y_arr = m(np.max(lon_arr), np.max(lat_arr))
#     dlon_arr = np.zeros([nb_cell_lon[1]+1])
#     for ilon in range(nb_cell_lon[1]+1):
#         dlon_arr[ilon]  =  min_x_arr + ilon * (max_x_arr - min_x_arr)/nb_cell_lon[1]
#     dlat_arr = np.zeros([nb_cell_lat[1]+1])
#     for ilat in range(nb_cell_lat[1]+1):
#         dlat_arr[ilat]  = 0 + ilat * (max_y_arr - min_y_arr)/nb_cell_lat[1]

#     m.pcolormesh(dlon_arr, dlat_arr, data, cmap=cmap, linewidth = 0.01,  edgecolors = 'lightgrey', latlon=False) # !!!!!! be careful: X and Y need to be the number of element of data + 1 (), otherwise one cell will be missing (https://matplotlib.org/basemap/api/basemap_api.html#mpl_toolkits.basemap.Basemap.pcolormesh)

# #    ax_2d_big.plot(0,0, color = 'lightblue', linewidth=2)
# #    legend = ax_2d_big.legend( loc='center',  bbox_to_anchor=(0.5, -0.144), fontsize = fontsize_plot, handles=handles_arr, ncol=2, columnspacing = height_fig * ratio_fig_size*1.1,frameon=False)#    legend = ax_2d_big.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="", fontsize = fontsize_plot, handles=handles_arr)
#     specular_list = []
#     point = namedtuple('point', ['x', 'y'])
#     color = namedtuple('color', 'red green blue')
#     specular = namedtuple('specular',('name',) +  point._fields  + color._fields + ('point_plot',))

#     if itime < nb_second_visu: # the spec are shown only for the firset nb_second_visu seconds of the simulation (this is because we don't want to have too heavy pickles for save_lon_spec and save_lat_spec)
#         for isc in range(nb_sc):
#             # Add on plot the specular points over one orbit
#             for k in range(nb_spec):
#                 specular_list.append(specular)
# #                 if itime > 3600:
# #                     specular_list[k+isc*nb_spec].x, specular_list[k+isc*nb_spec].y =  m(lon_spec[k,isc,itime-3600:itime], lat_spec[k,isc,itime-3600:itime])
# #                 else:
#                 specular_list[k+isc*nb_spec].x, specular_list[k+isc*nb_spec].y =  m(lon_spec[k,isc,:itime], lat_spec[k,isc,:itime])
#                 specular_list[k+isc*nb_spec].point_plot = m.scatter(specular_list[k+isc*nb_spec].x, specular_list[k+isc*nb_spec].y, marker='o', s = 50*np.cos(max_lat_grid_array[ilat_center]*deg2rad) * width_cell[1] / np.max(width_cell_array) , color = 'b')

#     ax_2d_big.set_title(ax_2d_big_title, weight = 'bold', fontsize  = (int)(fontsize_plot*1.1), y = 1.0)
#     ax_2d_big.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
#     ax_2d_big.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot, zorder = 5)

