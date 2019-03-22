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



import os
from os import listdir
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
import matplotlib.cm as cmx
from os.path import isfile, join
from datetime import datetime, timedelta
import numpy as np
from matplotlib  import pyplot as plt
import sys
import ipdb
import matplotlib.gridspec as gridspec
import scipy.stats
compute_dist = 0

def read_obs(format_data, filename):
    # filename = 'omniWebObs/b_v_n_t_19970101_to_20181231.txt'
    # format_data = 'ccmc'
    if format_data == 'ccmc':
        file_obs = open(filename)
        read_file_obs = file_obs.readlines()
        nheader = 0
        while ('YEAR' in read_file_obs[nheader]) == False:
            nheader = nheader + 1
        nheader = nheader + 1
        nline = len(read_file_obs) - nheader
        date_obs = []
        b_obs = []; v_obs = []; n_obs = []; t_obs = []
        for iline in range(nline):
            if read_file_obs[iline+nheader][0] == '<':
                break
            #1997   1  0 999.9 9999. 999.9 9999999.
            date_temp = read_file_obs[iline+nheader].split()[0] + '-' +read_file_obs[iline+nheader].split()[1] + 'T' + read_file_obs[iline+nheader].split()[2] 
            date_obs.append(datetime.strptime(date_temp, "%Y-%jT%H"))
            b_obs.append(np.float(read_file_obs[iline+nheader].split()[3]))
            v_obs.append(np.float(read_file_obs[iline+nheader].split()[4]))
            n_obs.append(np.float(read_file_obs[iline+nheader].split()[5]))
            t_obs.append(np.float(read_file_obs[iline+nheader].split()[6]))
        b_obs = np.array(b_obs); v_obs = np.array(v_obs); n_obs = np.array(n_obs); t_obs = np.array(t_obs)
        date_obs = np.array(date_obs)
        return date_obs, b_obs, v_obs, n_obs, t_obs

def read_pred(format_data, filename):
    # filename = 'ccmcPred/20160220_190100_2.0_ENLIL_time_line.dat'
    # format_data = 'ccmc'
    if format_data == 'ccmc': 
        file_pred = open(filename)
        read_file_pred = file_pred.readlines()
        nheader = 0
        while read_file_pred[nheader][0] != '2':
            nheader = nheader + 1
        nline = len(read_file_pred) - nheader
        b_pred = []; v_pred = []; n_pred = []; t_pred = []
        date_pred = []
        b_to_ave = []; v_to_ave = []; n_to_ave = []; t_to_ave = []
        b_to_ave_all = []; v_to_ave_all = []; n_to_ave_all = []; t_to_ave_all = []
        for iline in range(nline):
            if iline != 0:
                previous_date = date_temp_date
                previous_hour = hour
            date_temp = read_file_pred[iline+nheader].split()[0] + '-' + read_file_pred[iline+nheader].split()[1] + '-' + read_file_pred[iline+nheader].split()[2] + 'T' + read_file_pred[iline+nheader].split()[3] + ':' + read_file_pred[iline+nheader].split()[4]
            date_temp_date = datetime.strptime(date_temp, "%Y-%m-%dT%H:%M")
            hour = (int)(datetime.strftime(date_temp_date, "%H"))
            if iline > 0:
                if hour != previous_hour:
                    date_str = datetime.strftime(previous_date, "%Y-%m-%dT%H") + ":00:00"
                    date_pred.append(datetime.strptime(date_str, "%Y-%m-%dT%H:%M:%S"))
                    b_to_ave_all.append(b_to_ave); v_to_ave_all.append(v_to_ave); n_to_ave_all.append(n_to_ave); t_to_ave_all.append(t_to_ave)
                    b_pred.append(np.mean(b_to_ave)); v_pred.append(np.mean(v_to_ave)); n_pred.append(np.mean(n_to_ave)); t_pred.append(np.mean(t_to_ave))  
                    b_to_ave = []; v_to_ave = []; n_to_ave = []; t_to_ave = []

            b_to_ave.append(np.float(read_file_pred[iline+nheader].split()[5]))
            v_to_ave.append(np.float(read_file_pred[iline+nheader].split()[6]))
            n_to_ave.append(np.float(read_file_pred[iline+nheader].split()[7]))
            t_to_ave.append(np.float(read_file_pred[iline+nheader].split()[8]))

        b_to_ave_all.append(b_to_ave); v_to_ave_all.append(v_to_ave); n_to_ave_all.append(n_to_ave); t_to_ave_all.append(t_to_ave)
        b_pred.append(np.mean(b_to_ave)); v_pred.append(np.mean(v_to_ave)); n_pred.append(np.mean(n_to_ave)); t_pred.append(np.mean(t_to_ave))  
        date_str = datetime.strftime(date_temp_date, "%Y-%m-%dT%H") + ":00:00"
        date_pred.append(datetime.strptime(date_str, "%Y-%m-%dT%H:%M:%S"))
        b_pred = np.array(b_pred); v_pred = np.array(v_pred); n_pred = np.array(n_pred); t_pred = np.array(t_pred)*1000
        date_pred = np.array(date_pred)
        return date_pred, b_pred, v_pred, n_pred, t_pred

filename_pred = 'ccmcPred/20150315_064600_2.0_ENLIL_time_line.dat'#'ccmcPred/20160220_190100_2.0_ENLIL_time_line.dat'
format_data = 'ccmc'
date_pred, b_pred, v_pred, n_pred, t_pred = read_pred(format_data, filename_pred)
date_start = str(date_pred[0])[:10].replace('-', '')
date_stop = str(date_pred[-1])[:10].replace('-', '')
filename_obs = 'omniWebObs/b_v_n_t_' + date_start + '_to_' + date_stop + '.txt'
if (os.path.isfile(filename_obs) == False):
    os.system('wget --post-data "activity=retrieve&res=hour&spacecraft=omni2&start_date=' + date_start + '&end_date=' + date_stop + '&vars=8&vars=24&vars=23&vars=22&scale=Linear&ymin=&ymax=&view=0&charsize=&xstyle=0&ystyle=0&symbol=0&symsize=&linestyle=solid&table=0&imagex=640&imagey=480&color=&back=" https://omniweb.sci.gsfc.nasa.gov/cgi/nx1.cgi -O ' + filename_obs)
date_obs, b_obs, v_obs, n_obs, t_obs = read_obs(format_data, filename_obs)

# Compare observations to predictions
npred = len(date_pred)
error_b = np.zeros([npred]); error_v = np.zeros([npred]); error_n = np.zeros([npred]); error_t = np.zeros([npred])
nb_seconds_since_start = np.zeros([npred])
index_obs = []
for ipred in range(npred):
    # there might be gpas between predicitons (ie two conseucitve predictions might be separatd by more tah an hour
    iobs = np.where(date_obs == date_pred[ipred])[0][0]
    index_obs.append( iobs )
    error_b[ipred] = b_pred[ipred] - b_obs[iobs]
    error_v[ipred] = v_pred[ipred] - v_obs[iobs]
    error_n[ipred] = n_pred[ipred] - n_obs[iobs]
    error_t[ipred] = t_pred[ipred] - t_obs[iobs]
    nb_seconds_since_start[ipred] = (date_pred[ipred] - date_pred[0]).total_seconds()
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3
width_fig = 25
color_arr = ['blue', 'red', 'mediumorchid', 'dodgerblue', 'magenta', 'darkgreen', 'limegreen', 'black']
        
fig_title = 'WSA/ENLIL/Cone predictions VS observations: B, v, n and T (' + str(date_pred[0])[:10] + ' to ' +  str(date_pred[-1])[:10] + ')'
x_label = 'Time (days)'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal')
plt.rc('font', weight='normal') ## make the labels of the ticks in normal
gs = gridspec.GridSpec(2, 2)
gs.update(left = 0., right=1, top = 0.927,bottom = 0., hspace = 0.15, wspace = 0.2)
# IMF magnitude
y_label = 'B (nT)'
ax = fig.add_subplot(gs[0, 0])
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='normal') ## make the labels of the ticks in normal
ax.plot(nb_seconds_since_start / 3600./24, b_pred, linewidth = 2, color = 'blue', label ='Predictions')
ax.plot(nb_seconds_since_start / 3600./24, b_obs[index_obs], linewidth = 2, color = 'red', label ='Observations')
nrms = np.sqrt( np.mean((b_pred-b_obs[index_obs])**2) / np.mean(b_obs[index_obs]**2))
corr = scipy.stats.pearsonr(b_pred, b_obs[index_obs])[0]
ax.text(0.98, 0.98, 'NRMSE = ' + format(nrms, '.2f') + '\nr = ' + format(corr, '.2f'), transform = ax.transAxes, fontsize = fontsize_plot, horizontalalignment = 'right', verticalalignment = 'top')
legend = ax.legend(loc='upper left', bbox_to_anchor=(0, 1), numpoints = 1,  title="", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))
ax.margins(0,0)
# Velocity
y_label = 'v (km/s)'
ax = fig.add_subplot(gs[1, 0])
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='normal') ## make the labels of the ticks in normal
ax.plot(nb_seconds_since_start / 3600./24, v_pred, linewidth = 2, color = 'blue', label ='Predictions')
ax.plot(nb_seconds_since_start / 3600./24, v_obs[index_obs], linewidth = 2, color = 'red', label ='Observations')
nrms = np.sqrt( np.mean((v_pred-v_obs[index_obs])**2) / np.mean(v_obs[index_obs]**2))
corr = scipy.stats.pearsonr(v_pred, v_obs[index_obs])[0]
ax.text(0.98, 0.98, 'NRMSE = ' + format(nrms, '.2f') + '\nr = ' + format(corr, '.2f'), transform = ax.transAxes, fontsize = fontsize_plot, horizontalalignment = 'right', verticalalignment = 'top')

ax.margins(0,0)

# density
y_label = 'n (cm$^{-3}$)'
ax = fig.add_subplot(gs[0, 1])
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='normal') ## make the labels of the ticks in normal
ax.plot(nb_seconds_since_start / 3600./24, n_pred, linewidth = 2, color = 'blue', label ='Predictions')
ax.plot(nb_seconds_since_start / 3600./24, n_obs[index_obs], linewidth = 2, color = 'red', label ='Observations')
nrms = np.sqrt( np.mean((n_pred-n_obs[index_obs])**2) / np.mean(n_obs[index_obs]**2))
corr = scipy.stats.pearsonr(n_pred, n_obs[index_obs])[0]
ax.text(0.98, 0.98, 'NRMSE = ' + format(nrms, '.2f') + '\nr = ' + format(corr, '.2f'), transform = ax.transAxes, fontsize = fontsize_plot, horizontalalignment = 'right', verticalalignment = 'top')

ax.margins(0,0)

# temperature
y_label = 'T (kK)'
ax = fig.add_subplot(gs[1, 1])
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='normal') ## make the labels of the ticks in normal
ax.plot(nb_seconds_since_start / 3600./24, t_pred/1000., linewidth = 2, color = 'blue', label ='Predictions')
ax.plot(nb_seconds_since_start / 3600./24, t_obs[index_obs]/1000., linewidth = 2, color = 'red', label ='Observations')
ax.margins(0,0)
nrms = np.sqrt( np.mean((t_pred-t_obs[index_obs])**2) / np.mean(t_obs[index_obs]**2))
corr = scipy.stats.pearsonr(t_pred, t_obs[index_obs])[0]
ax.text(0.98, 0.98, 'NRMSE = ' + format(nrms, '.2f') + '\nr = ' + format(corr, '.2f'), transform = ax.transAxes, fontsize = fontsize_plot, horizontalalignment = 'right', verticalalignment = 'top')

fig_save_name = 'fig/ccmc_vs_obs_' + date_start + '_to_' + date_stop
fig_save_name =  fig_save_name + '.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
