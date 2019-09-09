# This script generates a SpOCK simulation and call cygnss_lakes for that simulation

# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT                                                                                                                                
date_start_simu = '2019-08-28'
date_stop_simu = '2020-02-28'
# end of  PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

raise Exception                                                

# Algorithm
from datetime import datetime, timedelta
import sys
import os
import ipdb
sys.path.append("/Users/cbv/work/spock/srcPython")
from read_input_file import *
from cygnss_read_spock_spec_bin import *
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from matplotlib import pyplot as plt
import cygnss_lakes; reload(cygnss_lakes); from cygnss_lakes import *

date_start_simu_date = datetime.strptime(date_start_simu, "%Y-%m-%d")
date_stop_simu_date = datetime.strptime(date_stop_simu, "%Y-%m-%d")
nday = (date_stop_simu_date - date_start_simu_date).days
date_arr_simu = np.array([date_start_simu_date + timedelta(days=i) for i in np.arange(0, nday, 1)])
visit_time = []; visit_which_sp = []; visit_which_sc = []
for iday in range(133, nday):#!!!!!
    date_start_date = date_arr_simu[iday]
    date_stop_date = date_start_date + timedelta(days = 1)
    date_start = datetime.strftime(date_start_date, "%Y-%m-%d") + "T00:00:00"
    spock_input_filename = date_start[:10].replace('-','') + '.txt'
    print '\niday', iday, nday-1,  spock_input_filename, str(datetime.now())[0:19]
    # Create the SpOCK main input file
    os.system("cp 082819.txt " + spock_input_filename)
    date_stop = datetime.strftime(date_stop_date, "%Y-%m-%d") + "T00:00:00" 
    with open(spock_input_filename) as f:
        lines = f.readlines()
        lines[1] = date_start + '\n'
        lines[2] = date_stop + '\n' 
    with open(spock_input_filename, "w") as f:
        f.writelines(lines)
    
    # Run SpOCK to predict the specular point positions
    if iday >= 132:
        os.system("mpirun -np 4 spock_dev " + spock_input_filename)
        os.system("mpirun -np 4 spec " + spock_input_filename + " -lon=0 -rot=0 -min")

    # Call cygnss_lakes to compute the revisit times
    visit_time_day, visit_which_sp_day, visit_which_sc_day = cygnss_lakes(spock_input_filename)
    print visit_time_day
    visit_time.append(visit_time_day); visit_which_sp.append(visit_which_sp_day); visit_which_sc.append(visit_which_sc_day)


nday = len(visit_time)
nlake = len(visit_time[0])
min_time_between_two_visit = np.zeros([nlake])
# !!!!!!!! if you change anything about the lakes you also neeed to change it in the funciton cygnss_lake (yes it's bad)
upper_left = [[32.992112, -107.335678], [13.900429, -89.585032], [13.713906, -89.107452], [37.157854, 99.701568]]
lower_right = [[32.890515, -107.241349], [13.831710, -89.515894], [13.629988, -88.986665], [36.548034, 100.745112]]
name = ['Caballo', 'Coatepeque', 'LLopango', 'Qinghai']
    
for ilake in range(nlake):
    min_lat = upper_left[ilake][0]
    min_lon = upper_left[ilake][1]
    min_lon_corr = min_lon
    if min_lon < 0:
        min_lon_corr = 360 + min_lon
    upper_left0_360.append([min_lat, min_lon_corr])
    max_lat = lower_right[ilake][0]
    max_lon = lower_right[ilake][1]
    max_lon_corr = max_lon
    if max_lon < 0:
        max_lon_corr = 360 + max_lon
    lower_right0_360.append([max_lat, max_lon_corr])
    
    max_lat = upper_left0_360[ilake][0]
    min_lon = upper_left0_360[ilake][1]
    min_lat = lower_right0_360[ilake][0]
    max_lon = lower_right0_360[ilake][1]
    max_dim_lake = np.max([max_lon - min_lon, max_lat - min_lat])
    min_time_between_two_visit[ilake] = max_dim_lake*np.sqrt(2)*110/7.5

visit_time_all_day = []
visit_time_all_day_lake = []
for ilake in range(nlake): # first time step
    visit_time_all_day_lake.append(visit_time[0][ilake])
visit_time_all_day.append(visit_time_all_day_lake)
for iday in range(1, nday):
    visit_time_all_day_lake = []
    for ilake in range(nlake):
        # last visit of previous day
        last_visit_prev = visit_time[iday-1][ilake][-1]
        # first visit current day
        first_visit_curr = visit_time[iday][ilake][0] + 86400*iday # !!!!!!!!asummes that each simu is exactly one day long
        # only keep the first visit of current day if it occurs later than 
        if first_visit_curr - last_visit_prev > [ilake]:
            visit_time_all_day_lake.append(visit_time[iday][ilake][0])
        # then add all the subsequent visits of the day
        nvisit_day_lake = len(visit_time[iday][ilake])
        for ivisit in range(1, nvisit_day_lake):
            visit_time_all_day_lake.append(visit_time[iday][ilake][ivisit])
        visit_time_all_day[ilake] = visit_time_all_day[ilake] + visit_time_all_day_lake

# Calculate the revisit times for each lake
revisit_time = []
for ilake in range(nlake):
    revisit_time_lake = []
    for iday in range(1, nday):
        revisit_time_lake.append(visit_time_all_day[ilake][iday] - visit_time_all_day[ilake][iday-1])
    revisit_time.append(revisit_time_lake)
    
        
# Historgram of revisit times for each lake
height_fig = 15.  # the width is calculated as height_fig * 4/3. 
fontsize_plot = 25
ratio_fig_size = 4./3
fig_title = 'Distribution of revisit times for Lake ' + name[ilake] + ' over ' + str((int)(nday)) + ' days'
y_label = '# day'
x_label = 'Revisit time (hours)'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.958,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                                                                                                                  
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.17)
ax = fig.add_subplot(gs[0, 0])
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                                                                                                                       
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold

#!!!!!! for ilake in range(nlake):
ilake = 0 #!!!!!!!!!remove
revisit_time_lake = np.array(revisit_time[ilake])
data = revisit_time_lake / 3600.
range_min = 0
range_max = (int)(np.ceil(np.max(data)))
nbins = 10
hist_along_data = np.histogram(data, range = [range_min, range_max], bins = nbins)
bin_array_temp = hist_along_data[1]
bin_array = ( bin_array_temp[:-1] + np.roll(bin_array_temp,-1)[:-1] ) /2.
binsize_actual = bin_array[1] - bin_array[0]
hist_along = hist_along_data[0] * 100. / len(cov_spock_both_pole_arr)
ax.bar(bin_array, hist_along, binsize_actual)
ax.set_ylim([0, np.max(hist_along)*1.1])
ax.set_xlim([np.min(bin_array_temp), np.max(bin_array_temp)])

fig.set_figheight(height_fig)
fig.set_figwidth(height_fig*ratio_fig_size)
fig_save_name = 'dist_revisit_time_' + name[ilake].lower() + '.pdf'

fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')


