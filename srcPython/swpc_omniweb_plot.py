# This script plots a variable (var) from SWPC or Omniweb  between date_start and date_stop
# Assumptions:
# - see section PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT


# PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT
source = 'omniweb' # omniweb, swpc
date_start = '2017-09-01T00:00:00' # YYYY-mm-ddTHH:MM:SS
date_stop = '2017-09-09T00:00:00' # YYYY-mm-ddTHH:MM:SS
var_name = ['f107', 'ap'] # list: f107, ap, dst
# end of PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT

def plot_var(fig_title_h, y_label_h, date_date_h, var_h, nb_seconds_since_start_h, fig_save_name_h):
    x_label = 'Real time'
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle(fig_title_h, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
    ax = fig.add_subplot(gs[0, 0])
    ax.set_ylabel(y_label_h, weight = 'normal', fontsize  = fontsize_plot)
    ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    ax.plot(nb_seconds_since_start_h, var_h, linewidth = 2, color = 'k')
    ax.margins(0,0)
    ax.set_ylim([np.min(var_h)*0.9, np.max(var_h)*1.1])
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

from datetime import datetime, timedelta
import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec


date_start_date = datetime.strptime(date_start, '%Y-%m-%dT%H:%M:%S')
date_stop_date = datetime.strptime(date_stop, '%Y-%m-%dT%H:%M:%S')

date_stop_date_plus = date_stop_date + timedelta(days = 1)
date_stop_plus = datetime.strftime(date_stop_date_plus, '%Y-%m-%dT%H:%M:%S')
date_start_omni = date_start[0:10].replace('-', '') # only vary by step sizes of a day
date_stop_omni = date_stop_plus[0:10].replace('-', '')
nvar = len(var_name)

## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25
ratio_fig_size = 4./3

if source == 'omniweb':
    var_omni = []
    for ivar in range(nvar):
        # Download file
        if var_name[ivar] == 'f107':
            var_omni = '50'
            var_name_plot = 'F10.7'
        elif var_name[ivar] == 'ap':
            var_omni = '49'
            var_name_plot = 'Ap'
        filename_out = date_start_omni + '_to_' + date_stop_omni + '_' + source + '_' + var_name[ivar] + '.txt'
        os.system( 'wget --no-check-certificate --post-data "activity=retrieve&res=hour&spacecraft=omni2&start_date=' + date_start_omni + '&end_date=' + date_stop_omni + '&vars=' + var_omni + '&scale=Linear&ymin=&ymax=&view=0&charsize=&xstyle=0&ystyle=0&symbol=0&symsize=&linestyle=solid&table=0&imagex=640&imagey=480&color=&back=" https://omniweb.sci.gsfc.nasa.gov/cgi/nx1.cgi -O ' + filename_out )

        # Read file
        file_out = open(filename_out)
        read_file_out = file_out.readlines()
        nheader = 0
        while read_file_out[nheader][0] != '2':
            nheader = nheader + 1
            if len(read_file_out[nheader]) == 0:
                nheader = nheader +	1
        ntemp = len(read_file_out) - nheader
        i = 0
        var = []
        date = []
        date_date = []
        nb_seconds_since_start = []
        while read_file_out[nheader + i][0] == '2':
            date_temp = read_file_out[nheader + i].split()[0] + '-' + read_file_out[nheader + i].split()[1] + '-' + read_file_out[nheader + i].split()[2]
            date_date.append(datetime.strptime(date_temp, '%Y-%j-%H'))
            date.append(datetime.strftime(date_date[-1], '%Y-%m-%dT%H:%M:%S'))
            var.append(np.float(read_file_out[nheader + i].split()[3]))
            if date_date[-1] == date_start_date:
                istart = i
            if date_date[-1] == date_stop_date:
                istop = i+1
            nb_seconds_since_start.append((date_date[-1] - date_start_date).total_seconds())
            i = i + 1
        file_out.close()
        date_date = np.array(date_date)
        var = np.array(var)
        nb_seconds_since_start = np.array(nb_seconds_since_start)
        
        # Plot
        istart = np.where(date_date == date_start_date)[0][0]
        istop = np.where(date_date == date_stop_date)[0][0]
        
            
        fig_title = var_name_plot + ' as a function of time'
        fig_save_name = filename_out.replace('.txt', '.pdf')
        plot_var(fig_title, var_name_plot, date_date[istart:istop], var[istart:istop], nb_seconds_since_start[istart:istop], fig_save_name)
