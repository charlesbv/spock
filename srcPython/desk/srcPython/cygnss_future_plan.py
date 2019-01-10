# This script plots the future plan diff drag maneuvers for CYGNSS based on the txt files created by the Matlab script drag_plots.m wrote by Kyle's Nave team. This has been used for the paper on differrential drag maneuvers submitted to the CYGNSS special issue of JSTAR


import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from matplotlib import pyplot as plt
from datetime import datetime, timedelta
import numpy as np

nb_maneuver_at_a_time_arr = [2,3,4]
n = len(nb_maneuver_at_a_time_arr)
nb_sc = 7
for i in range(n):
    nb_maneuver_at_a_time = nb_maneuver_at_a_time_arr[i]
    filename = str(nb_maneuver_at_a_time) + 'atatime.txt'
    file = open(filename)
    read_file = file.readlines()
    nb_header = 1
    nb_time = len(read_file) - nb_header
    if i == 0:# assumes all files have the same date data
        phase_rate = np.zeros([n, nb_time, nb_sc])
        x_axis = np.zeros([nb_time])
        date = []
        name = read_file[0].split()[1:]    
    for itime in range(nb_time):
        if i == 0:
            date.append( read_file[nb_header + itime].split()[0])
            x_axis[itime] = (datetime.strptime(date[-1], "%d-%b-%Y") - datetime.strptime(date[0], "%d-%b-%Y")).total_seconds()            
        phase_rate[i, itime, :] = read_file[nb_header + itime].split()[1:]
    file.close()

date_ref = datetime.strptime(date[0], "%d-%b-%Y") 
# read validManeuverDates.csv to get the times when no maneuver is performed
filename = 'validManeuverDates_longer_spock.txt'#'validManeuverDates.csv'
file =open(filename)
read_file = file.readlines()
nb_man_ok = len(read_file) - 1
date_man_ok = []
date_beta  = []
date_beta_seconds  = []
for i in range(nb_man_ok):
    date_man_ok.append( datetime.strptime(read_file[1  + i].rstrip(), "%d %b %Y %H:%M:%S.%f") )
    if i > 0:
        if (date_man_ok[-1] - date_man_ok[-2]).total_seconds() > 24 * 3600: # if more than one day then no man is perfomred because of beta anle
            date_beta.append( [date_man_ok[-2], date_man_ok[-1]-timedelta(days = 1)])
            date_beta_seconds.append([(date_man_ok[-2] - date_ref).total_seconds(), (date_man_ok[-1]-timedelta(days = 1) - date_ref).total_seconds() ] )
nb_beta = len(date_beta)
date_beta_seconds = np.array(date_beta_seconds)
## Parameters for the figure
height_fig = 14  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25
width_fig = 11*4./3

NCURVES = 8
np.random.seed(101)
curves = [np.random.random(20) for i in range(NCURVES)]
values = range(NCURVES)
jet = cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
color_arr = ['b', 'r','k','g', 'm', 'gold', 'cyan', 'darkgray',  'lawngreen', 'green', 'chocolate','cornflowerblue','fuchsia']

fig_title = ''
y_label = '$\Delta \dot \Phi$ (' + u'\N{DEGREE SIGN}' + '/day)'
x_label = 'Real time'
fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in normal
gs = gridspec.GridSpec(n, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.1)

for iman in range(n): # 
    ax = fig.add_subplot(gs[iman, 0])
    ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    plt.rc('font', weight='normal') ## make the labels of the ticks in normal

    for isc in range(8): 
        # so this figure the same color code as the hsitory of sma figure
        if isc < 2:
            ax.plot(x_axis,phase_rate[iman,:,isc] , linewidth = 2, color = color_arr[isc], label = name[isc])
        if isc > 2:
            ax.plot(x_axis,phase_rate[iman,:,isc-1] , linewidth = 2, color = color_arr[isc], label = name[isc-1])

    for ibeta in range(nb_beta):
        ax.axvspan(date_beta_seconds[ibeta, 0], date_beta_seconds[ibeta, 1],  color='b', alpha = 0.05)



    # x axis label is in real time
                    ## all output files of one simulation have the same number of steps, and start at the same date
    nb_ticks_xlabel = 8
    nb_seconds_in_simu = x_axis[-1]
    dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
    xticks = np.arange(0, nb_seconds_in_simu+1, dt_xlabel)
    date_list_str = []
    date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
    for i in range(len(xticks)):
#         if dt_xlabel > nb_ticks_xlabel*24*3600:
#             date_list_str.append( str(date_list[i])[5:10] )
#         else:
        date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[0:4] )
        ax.xaxis.set_ticks(xticks)
        if iman == n - 1:
            ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
            ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
        else:
            ax.xaxis.set_ticklabels([])

    ax.set_xlim([min(xticks), max(xticks)-3*30*24*3600])
    ax.set_ylim([-0.1, np.max(phase_rate)*1.1])

    nb_maneuver_at_a_time = nb_maneuver_at_a_time_arr[iman]
    ax.text(0.98,0.94, str(nb_maneuver_at_a_time)  + ' CYGNSS at a time', transform = ax.transAxes, fontsize = fontsize_plot, horizontalalignment = 'right', verticalalignment = 'top')
#    ax.margins(0,0); 

                    #        ax.set_xlim([ax.get_xlim()[0], most_recent_tle_among_all_sc])

#                     legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
#                     legend.get_title().set_fontsize(str(fontsize_plot))



legend = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.35), numpoints = 1,  title="", fontsize = fontsize_plot, ncol=4)
fig.set_figheight(height_fig)

fig_save_name = 'futureplan.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  





# 1 A 
# 5 E
# 8 H
# 2 B
# 4 D
# 6 F
# 7 G
