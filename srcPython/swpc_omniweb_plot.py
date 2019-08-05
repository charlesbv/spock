# This script plots a variable (var) from SWPC or Omniweb or user's own file between date_start and date_stop
# Assumptions:
# - see section PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT
# - the user's own filename must have the following format (example):
#BEGINNINGOFHEADER
#ENDOFHEADER
#YEAR DOY HR  1
#2017 239  0   9
#2017 239  1   9
#ENDOFFILE
#  -> this format matches the one used in SpOCK whwen inputing our own f107 or ap files
# - If the users input their own file then this file must contain only one varilable (f107 or ap). Explicitely tell this script which variable it is by setting the varilable var_name to ['f107'] or ['ap']
# - the SWPC option is incomplete because it doesn't treat the case when date start and date stop are for two different swpc files (example:if date_start is in 2018 and date_stop in 2019)

# PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT
source = 'omniweb'#'20170901_to_20170910_omniweb_f107_no_storm.txt'#'20170827_to_20170910_omniweb_ap_no_storm.txt' # omniweb, swpc, [filename]
date_start = '2017-09-01T00:00:00' # YYYY-mm-ddTHH:MM:SS
date_stop = '2017-09-09T00:00:00' # YYYY-mm-ddTHH:MM:SS
var_name = ['f107'] # list: f107, ap, dst
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
    #ax.set_ylim([0, 252])
    ax.set_ylim([80, 200])
    #ax.set_ylim([np.min(var_h)*0.9, np.max(var_h)*1.1])

    nb_ticks_xlabel = 8
    nb_seconds_in_simu = nb_seconds_since_start_h[-1] - nb_seconds_since_start_h[0]
    dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
    xticks = np.arange(0, nb_seconds_in_simu+1, dt_xlabel)
    date_list_str = []
    date_list = [date_start_date + timedelta(seconds=x-xticks[0]) for x in xticks]
    for i in range(len(xticks)):
        if dt_xlabel > nb_ticks_xlabel*24*3600:
            date_list_str.append( str(date_list[i])[5:10] )
        else:
            date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
    ax.xaxis.set_ticks(xticks)
    ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

from datetime import datetime, timedelta
import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import ipdb

date_start_date = datetime.strptime(date_start, '%Y-%m-%dT%H:%M:%S')
date_stop_date = datetime.strptime(date_stop, '%Y-%m-%dT%H:%M:%S')

date_stop_date_plus = date_stop_date + timedelta(days = 1)
date_stop_plus = datetime.strftime(date_stop_date_plus, '%Y-%m-%dT%H:%M:%S')
nvar = len(var_name)

## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25
ratio_fig_size = 4./3

if ((source != 'omniweb') & (source != 'swpc')):
    source_ok = 'user_file'
else:
    source_ok = source
date_start_source = date_start[0:10].replace('-', '') # only vary by step sizes of a day
date_stop_source = date_stop_plus[0:10].replace('-', '')

if ((source_ok == 'omniweb') | (source_ok == 'user_file')):
    var_omni = []
    for ivar in range(nvar):
        if var_name[ivar] == 'f107':
            var_omni = '50'
            var_name_plot = 'F10.7'
        elif var_name[ivar] == 'ap':
            var_omni = '49'
            var_name_plot = 'Ap'
        # Download file
        if source == 'omniweb':
            filename_out = date_start_source + '_to_' + date_stop_source + '_' + source + '_' + var_name[ivar] + '.txt'
            os.system( 'wget --no-check-certificate --post-data "activity=retrieve&res=hour&spacecraft=omni2&start_date=' + date_start_source + '&end_date=' + date_stop_source + '&vars=' + var_omni + '&scale=Linear&ymin=&ymax=&view=0&charsize=&xstyle=0&ystyle=0&symbol=0&symsize=&linestyle=solid&table=0&imagex=640&imagey=480&color=&back=" https://omniweb.sci.gsfc.nasa.gov/cgi/nx1.cgi -O ' + filename_out )
        else:
            filename_out = source
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

        if source_ok == 'omniweb':
            fig_title = var_name_plot + ' as a function of time - Omniweb'
        else:
            fig_title = var_name_plot + ' as a function of time - No storm'
        fig_save_name = filename_out.replace('.txt', '.pdf')
        plot_var(fig_title, var_name_plot, date_date[istart:istop], var[istart:istop], nb_seconds_since_start[istart:istop], fig_save_name)
        
        
else: # swpc !!!!!!!! incomplete because doesn't treat the case when date start and date stop are for two different swpc files
    var_swpc = []
    for ivar in range(nvar):
        if var_name[ivar] == 'f107':
            var_swpc = '_DSD'
            var_name_plot = 'F10.7'
        elif var_name[ivar] == 'ap':
            var_swpc = '_DGD'
            var_name_plot = 'Ap'
        current_year = datetime.strftime(datetime.today(), "%Y")
        if date_start[0:4] == current_year:
            # figure out which quarter is date_start
            if date_start_date <= datetime.strptime('2019-03-31', '%Y-%m-%d'):
                quarter = 'Q1'
            elif date_start_date <= datetime.strptime('2019-06-30', '%Y-%m-%d'):
                quarter = 'Q2'
            elif date_start_date <= datetime.strptime('2019-09-30', '%Y-%m-%d'):
                quarter = 'Q3'
            else:
                quarter = 'Q4'
            var_swpc = quarter + var_swpc
        # Download file
        filename_out = date_start_source + '_to_' + date_stop_source + '_' + source + '_' + var_name[ivar] + '.txt'
        #os.system('wget --no-check-certificate ftp://ftp.swpc.noaa.gov/pub/indices/old_indices/' + date_start[0:4] + var_swpc + '.txt -O ' + filename_out )
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
        while (read_file_out[nheader + i][0] == '2'):
            date_temp = read_file_out[nheader + i].split()[0] + '-' + read_file_out[nheader + i].split()[1] + '-' + read_file_out[nheader + i].split()[2]
            for ihour in range(24):
                date_date.append(datetime.strptime(date_temp, '%Y-%m-%d') + timedelta(hours = ihour))
                date.append(datetime.strftime(date_date[-1], '%Y-%m-%dT%H:%M:%S'))
                if var_name[ivar] == 'ap':
                    var.append(np.float(read_file_out[nheader + i][59:62]))
                elif var_name[ivar] == 'f107':
                    var.append(np.float(read_file_out[nheader + i].split()[3]))
                if date_date[-1] == date_start_date:
                    istart = i
                if date_date[-1] == date_stop_date:
                    istop = i+1
                nb_seconds_since_start.append((date_date[-1] - date_start_date).total_seconds())
            i = i + 1
            if i == ntemp -1 :
                break

        file_out.close()
        date_date = np.array(date_date)
        var = np.array(var)
        nb_seconds_since_start = np.array(nb_seconds_since_start)

        # Plot
        istart = np.where(date_date == date_start_date)[0][0]
        istop = np.where(date_date == date_stop_date)[0][0]


        fig_title = var_name_plot + ' as a function of time - SWPC'
        fig_save_name = filename_out.replace('.txt', '.pdf')
        plot_var(fig_title, var_name_plot, date_date[istart:istop], var[istart:istop], nb_seconds_since_start[istart:istop], fig_save_name)
