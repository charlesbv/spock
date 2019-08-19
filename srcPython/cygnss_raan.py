# This script downloads all TLEs for CYGNSS and plot the RAAN VS time

# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
download_tle = 0
fm01_info = ['FM01', '41887', '2019-01-01', '2019-08-18', 'mediumorchid'] # important to put the dates when the sc was operational (and not just in space after being decomissionned)
fm02_info = ['FM02', '41886', '2019-01-01', '2019-08-18', 'red'] # important to put the dates when the sc was operational (and not just in space after being decomissionned)
fm03_info = ['FM03', '41891', '2019-01-01', '2019-08-18', 'limegreen'] # important to put the dates when the sc was operational (and not just in space after being decomissionned)
fm04_info = ['FM04', '41885', '2019-01-01', '2019-08-18', 'blue'] # important to put the dates when the sc was operational (and not just in space after being decomissionned)
fm05_info = ['FM05', '41884', '2019-01-01', '2019-08-18', 'black'] # important to put the dates when the sc was operational (and not just in space after being decomissionned)
fm06_info = ['FM06', '41889', '2019-01-01', '2019-08-18', 'magenta'] # important to put the dates when the sc was operational (and not just in space after being decomissionned)
fm07_info = ['FM07', '41890', '2019-01-01', '2019-08-18', 'darkgreen'] # important to put the dates when the sc was operational (and not just in space after being decomissionned)
fm08_info = ['FM08', '41888', '2019-01-01', '2019-08-18', 'dodgerblue'] # important to put the dates when the sc was operational (and not just in space after being decomissionned)
arr_info = [fm01_info, fm02_info, fm03_info, fm04_info, fm05_info, fm06_info, fm07_info, fm08_info]
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT


# satColors = ['black', 'blue', 'red', 'mediumorchid', 'dodgerblue', 'magenta', 'darkgreen', 'limegreen'] #['lime', 'blue', 'indigo', 'purple', 'dodgerblue', 'steelblue', 'seagreen', 'limegreen']
# nb_sc = 8 # !!!!!!!!!
# label_arr = ['FM05', 'FM04', 'FM02', 'FM01', 'FM08', 'FM06', 'FM07', 'FM03']

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

mu_earth = 398600.4418 # km^3/s^2
earth_radius = 6378.137 # mean equatorial radius (km)   
j2 = 1.08262668e-3

nsc = len(arr_info)
actual_raan = []
j2_raan = []
ecc = []
sma = []
j2_prec_rate = []
actual_prec_rate = []
inc = []
ang_velo = []
epoch_start = []
epoch_stop = []
tle_epoch = []
handles_arr = []
for isc in range(nsc):
    epoch_start.append(datetime.strptime(arr_info[isc][2], "%Y-%m-%d"))
    epoch_stop.append(datetime.strptime(arr_info[isc][3], "%Y-%m-%d"))
    handles_arr.append(Line2D([0], [0], color =  arr_info[isc][-1], lw=4, label = arr_info[isc][0]))

date_ref = np.min(epoch_start)
nb_seconds_since_ref = []
for isc in range(nsc):
    if download_tle == 1:
        os.system('python download_tle.py ' + arr_info[isc][1] + ' ' + arr_info[isc][2] + ' ' + arr_info[isc][3] )
    tle_filename = arr_info[isc][1] + '_' + arr_info[isc][2] + '_' + arr_info[isc][3] + '.txt'
    tle_file = open(tle_filename)
    read_tle_file = tle_file.readlines()
    ntle = len(read_tle_file) / 2
    actual_raan_sc = []
    j2_raan_sc = []
    ecc_sc = []
    sma_sc = []
    j2_prec_rate_sc = []
    actual_prec_rate_sc = []
    inc_sc = []
    ang_velo_sc = []
    tle_epoch_sc = []
    nb_seconds_since_ref_sc = []
    for itle in range(ntle):
        iline = itle * 2
        actual_raan_temp = np.float(read_tle_file[iline + 1][17:25])
        actual_raan_sc.append(actual_raan_temp)
        ecc_temp = np.float('0.' + read_tle_file[iline + 1][26:33])
        ecc_sc.append(ecc_temp)
        mean_motion_temp = np.float(read_tle_file[iline + 1][52:63])        
        sma_temp = ( mu_earth / (mean_motion_temp/24/3600*2*np.pi)**2 )**(1./3)
        sma_sc.append(sma_temp)
        period_temp = 2*np.pi*np.sqrt(sma_temp**3/mu_earth)
        ang_velo_temp = 2*np.pi/period_temp
        ang_velo_sc.append(ang_velo_temp)
        inc_temp = np.float(read_tle_file[iline + 1][8:16])
        inc_sc.append(inc_temp)
        tle_epoch_sc_raw = read_tle_file[iline][18:32]
        tle_epoch_sc_temp = convert_tle_date_to_date(tle_epoch_sc_raw)
        tle_epoch_sc.append(tle_epoch_sc_temp)
        nb_seconds_since_ref_sc.append((tle_epoch_sc_temp - date_ref).total_seconds())
        j2_prec_rate_temp = -3./2*earth_radius**2/((sma_temp*(1-ecc_temp**2))**2) * j2 * ang_velo_temp * np.cos(inc_temp * np.pi/180.) * 180./np.pi # deg/s
        if itle > 0:
            j2_prec_rate_sc.append(j2_prec_rate_temp)
            if nb_seconds_since_ref_sc[-1] != nb_seconds_since_ref_sc[-2]:
                actual_prec_rate_temp = (actual_raan_sc[-1] - actual_raan_sc[-2]) / (nb_seconds_since_ref_sc[-1] - nb_seconds_since_ref_sc[-2])
                if np.abs(actual_raan_sc[-1] - actual_raan_sc[-2]) > 100:# this happens because of the 360 jump issue
                    actual_prec_rate_sc.append(actual_prec_rate_sc[-2])
                else:
                    actual_prec_rate_sc.append(actual_prec_rate_temp)
            else:
                actual_prec_rate_sc.append(actual_prec_rate_sc[-1])
            j2_raan_temp = np.mod(j2_prec_rate_sc[0]*(nb_seconds_since_ref_sc[-1] - nb_seconds_since_ref_sc[0]) + actual_raan_sc[0], 360) # assume that the raan varies linearly with the rate given by the j2 precession rate formula
        else:
            j2_raan_temp = actual_raan_sc[0]
        j2_raan_sc.append(j2_raan_temp)

    actual_raan.append(actual_raan_sc)
    j2_raan.append(j2_raan_sc)
    ecc.append(ecc_sc)
    sma.append(sma_sc)
    j2_prec_rate.append(j2_prec_rate_sc)
    actual_prec_rate.append(actual_prec_rate_sc)
    inc.append(inc_sc)
    ang_velo.append(ang_velo_sc)
    tle_epoch.append(tle_epoch_sc)
    nb_seconds_since_ref.append(nb_seconds_since_ref_sc)

final_prec_rate = np.zeros([nsc])
for isc in range(nsc):
    final_prec_rate[isc] = actual_prec_rate[isc][-1]
index_sort_prec_rate = np.argsort(final_prec_rate)

print 'Actual precession rates at tf'
for isc in range(nsc):
    print arr_info[index_sort_prec_rate[isc]][0], final_prec_rate[index_sort_prec_rate[isc]]
### parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25
ratio_fig_size = 4./3

fig_title = 'Right Ascension of the Ascending Node (RAAN) as a function of time'
y_label = 'RAAN ' + u'(\N{DEGREE SIGN})' #'Real time'
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
    ax.plot(nb_seconds_since_ref[isc], actual_raan[isc],  color = arr_info[isc][-1], label = arr_info[isc][0])
#    ax.plot(nb_seconds_since_ref[isc], j2_raan[isc],  color = 'black', label = 'j2')

# x axis label is in real time
nb_ticks_xlabel = 8
nb_seconds_in_simu = (np.max(epoch_stop) - date_ref).total_seconds()
dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
xticks = np.arange(0, nb_seconds_in_simu+1, dt_xlabel)
date_list_str = []
date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
for i in range(len(xticks)):
    if dt_xlabel > nb_ticks_xlabel*24*3600:
        date_list_str.append( str(date_list[i])[5:10] )
    else:
        date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
ax.xaxis.set_ticks(xticks)
ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
ax.margins(0,0); ax.set_xlim([min(xticks), max(xticks)])
legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Satellite", fontsize = fontsize_plot,  handles=handles_arr)
legend.get_title().set_fontsize(str(fontsize_plot))


fig_save_name = 'cyg_actual_raan.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


fig_title = 'Actual precession rate as a function of time'
y_label = 'Precession rate ' + u'(\N{DEGREE SIGN}/s)' #'Real time'
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
    ax.scatter(nb_seconds_since_ref[isc][1:], actual_prec_rate[isc],  s = 5, color = arr_info[isc][-1], label = arr_info[isc][0])

# x axis label is in real time
nb_ticks_xlabel = 8
nb_seconds_in_simu = (np.max(epoch_stop) - date_ref).total_seconds()
dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
xticks = np.arange(0, nb_seconds_in_simu+1, dt_xlabel)
date_list_str = []
date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
for i in range(len(xticks)):
    if dt_xlabel > nb_ticks_xlabel*24*3600:
        date_list_str.append( str(date_list[i])[5:10] )
    else:
        date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
ax.xaxis.set_ticks(xticks)
ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
ax.margins(0,0); ax.set_xlim([min(xticks), max(xticks)]); #ax.set_ylim([np.min(j2_prec_rate[isc]), np.max(j2_prec_rate[isc])])
legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Satellite", fontsize = fontsize_plot,  handles=handles_arr)
legend.get_title().set_fontsize(str(fontsize_plot))


fig_save_name = 'cyg_actual_prec_rate.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


