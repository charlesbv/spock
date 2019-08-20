# This script downloads all TLEs for satellites de-1, viking, polar, and image during their lifetime an plot the argument of apogee (as read from the TLEs) as a function of tim

# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
download_tle = 0
polar_info = ['Polar', '23802', '1996-02-24', '2008-04-28', 'blue'] # important to put the dates when the sc was operational (and not just in space after being decomissionned)
image_info = ['IMAGE', '26113', '2000-03-25', '2005-12-18', 'red']
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

xaxis_start_date = '1996-01-01' # need ot have the day to be 01
xaxis_stop_date = '2008-06-01' # need ot have the day to be 01
xaxis_nticks = 8


import sys
sys.path.append('/Users/cbv/work/spock/srcPython')
import os
from datetime import datetime, timedelta
import ipdb
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from convert_tle_date_to_date import *
from matplotlib.lines import Line2D


arr_info = [polar_info, image_info]
nsc = len(arr_info)
arg_apogee = []
epoch_start = []
epoch_stop = []
tle_epoch = []
handles_arr = []
for isc in range(nsc):
    epoch_start.append(datetime.strptime(arr_info[isc][2], "%Y-%m-%d"))
    epoch_stop.append(datetime.strptime(arr_info[isc][3], "%Y-%m-%d"))
    handles_arr.append(Line2D([0], [0], color =  arr_info[isc][-1], lw=4, label = arr_info[isc][0]))

xaxis_start_date_date = datetime.strptime(xaxis_start_date, "%Y-%m-%d")
date_ref = xaxis_start_date_date#np.min(epoch_start)
nb_seconds_since_ref = []
for isc in range(nsc):
    if download_tle == 1:
        os.system('python download_tle.py ' + arr_info[isc][1] + ' ' + arr_info[isc][2] + ' ' + arr_info[isc][3] )
    tle_filename = arr_info[isc][1] + '_' + arr_info[isc][2] + '_' + arr_info[isc][3] + '.txt'
    tle_file = open(tle_filename)
    read_tle_file = tle_file.readlines()
    ntle = len(read_tle_file) / 2
    arg_apogee_sc = []
    tle_epoch_sc = []
    nb_seconds_since_ref_sc = []
    for itle in range(ntle):
        iline = itle * 2
        arg_apogee_temp = np.mod(np.float(read_tle_file[iline + 1][34:42]) + 180., 360)
        arg_apogee_sc.append(arg_apogee_temp)
        tle_epoch_sc_raw = read_tle_file[iline][18:32]
        tle_epoch_sc_temp = convert_tle_date_to_date(tle_epoch_sc_raw)
        tle_epoch_sc.append(tle_epoch_sc_temp)
        nb_seconds_since_ref_sc.append((tle_epoch_sc_temp - date_ref).total_seconds())
    arg_apogee.append(arg_apogee_sc)
    tle_epoch.append(tle_epoch_sc)
    nb_seconds_since_ref.append(nb_seconds_since_ref_sc)

### Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25
ratio_fig_size = 4./3


fig_title = 'Argument of apogee as a function of time'
y_label = 'Argument of apogee ' + u'(\N{DEGREE SIGN})' #'Real time'
x_label = 'Real time'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='normal') ## make the labels of the ticks in bold

for isc in range(nsc):
    ax.scatter(nb_seconds_since_ref[isc], arg_apogee[isc], s = 5, color = arr_info[isc][-1], label = arr_info[isc][0])

    
# x axis label is in real time
xaxis_stop_date_date = datetime.strptime(xaxis_stop_date, "%Y-%m-%d")
delta_xaxis_start_stop = (xaxis_stop_date_date - xaxis_start_date_date).total_seconds() / 3600./ 24 / 30
nb_month_per_tick = (int)(delta_xaxis_start_stop / xaxis_nticks) + 1


xticks = np.zeros([xaxis_nticks+1])
date_list_str = []
xprevious = xaxis_start_date_date
date_list_str.append(datetime.strftime(xaxis_start_date_date, "%b") + '\n' + str(xaxis_start_date_date)[:4])
print nb_month_per_tick
for itick in range(xaxis_nticks-1):
    #ipdb.set_trace()
    xtick_temp = datetime.strftime(xprevious, "%Y-%m")
    xtick_temp_year = (int)(xtick_temp[:4])
    xtick_temp_month = (int)(xtick_temp[5:7])
    if xtick_temp_month + nb_month_per_tick > 12:
        xtick_temp_year_new = xtick_temp_year + (int)((xtick_temp_month + nb_month_per_tick)/12.)
        xtick_temp_month_new = np.mod(xtick_temp_month + nb_month_per_tick, 12)
        if xtick_temp_month_new == 0: # month 24, 36, 48, etc is Dec and the year should not be incremented by 2, 3, 4, etc but only 1, 2, 3, etc
            xtick_temp_month_new = 12
            xtick_temp_year_new = xtick_temp_year_new - 1
    else:
        xtick_temp_month_new = xtick_temp_month + nb_month_per_tick
        xtick_temp_year_new = xtick_temp_year
    xtick_temp_str = str(xtick_temp_year_new) + '-' + str(xtick_temp_month_new).zfill(2) + '-01'
    xtick_temp_date = datetime.strptime(xtick_temp_str, "%Y-%m-%d")
    xticks[itick+1] = (xtick_temp_date - xaxis_start_date_date).total_seconds() # xticks[0] is xaxis_start_date_date
    date_list_str.append(datetime.strftime(xtick_temp_date, "%b") + '\n' + xtick_temp_str[:4])
    print str(xprevious)[0:10], str(xtick_temp_date)[0:10]
    xprevious = xtick_temp_date
xticks[-1] = (xaxis_stop_date_date  - xaxis_start_date_date).total_seconds()
date_list_str.append(datetime.strftime(xaxis_stop_date_date, "%b") + '\n' + str(xaxis_stop_date_date)[:4])
ax.xaxis.set_ticks(xticks)
ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
ax.margins(0,0); ax.set_xlim([min(xticks), max(xticks)])
legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Satellite", fontsize = fontsize_plot,  handles=handles_arr)
legend.get_title().set_fontsize(str(fontsize_plot))


fig_save_name = 'arg_apog.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


