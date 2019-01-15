
# THis script plots the distance, amplitude, orbit average of runs amde with pid_algo_v2.py. The pickle were saved in pid_algo_v2.py
# inputs: pickle_root_list stores each pickle to load (one per run in pid_algo_v2.py) (the pickles are assumed ot be in ./pickle)
# (pickle_root =  prefix_name + '_' + rho_more in pid_algo_v2.py)
pickle_root_list = ['solpres05_mid']#['dt02_mid']#['noSolarPressure_mid']#['grav50_mid']#['dt02_mid']#['onlyDrag_mid']
#pickle_root_list = ['dec17_pole', 'dec17_mid', 'dec17_highamp_pole', 'dec17_equator']
#pickle_root_list = ['rhonosine_grav50_mid', 'rho0_grav50_solarzenith_mid','egm08_mid','grav80_mid', 'localtime70percent_mid']
#['rhonosine_grav50_mid', 'rho0_grav50_solarzenith_mid']#['egm08_mid']#['grav80_mid']#['rho0_grav50_solarzenith_mid']     #['dt0_1s_solarzenith_mid']#['grav50_solarzenith_mid']
#['localtime70percent_mid']#['localtime_pole', 'localtime_equator', 'localtime70percent_mid']
#['solarzenith_equator', 'solarzenith_pole', 'localtime70percent_mid']# ['localtime70percentAp2_mid']#

toplot = 'raw' # raw, amplitude
color_arr = ['blue', 'red', 'green', 'black', 'magenta']
isbig = 0
ispleiades = 0
import sys
import numpy as np

if isbig == 1:
    sys.path.append("/home/cbv/code/spock/srcPython")
    path_mpirun = '/usr/local/bin/mpirun'
    spice_path = '/raid4/cbv/cspice/data'
    nb_proc = 12
elif ispleiades == 1:
    sys.path.append("/home1/cbussy/Code/spock/srcPython")
    path_mpirun = 'mpiexec'
    spice_path = '/home1/cbussy/cspice/data'
    nb_proc = 0    

else:
    sys.path.append("/Users/cbv/work/spock/srcPython")
    path_mpirun = 'mpirun'
    spice_path = '/Users/cbv/cspice/data'
    nb_proc = 4


import pickle
import os

from read_input_file import *
from read_output_file import *
from spock_main_input import *
#if ispleiades != 1:
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
if ((isbig != 1) & (ispleiades !=1)):
    #from convert_cygnss_obs_ecef_to_eci import *
    import ipdb
from eci_to_lvlh import *
#plt.ion()



#FIGURES ###################

height_fig = 11
ratio_fig_size = 4./3
fontsize_plot = 20



######
fig_title = ''#'Distance between SpOCK and data for different density coefficient conditions' 
x_label = 'Time (hours)' 

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])


ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold


pickle_root_concatenate = ''
nb_pickle = len(pickle_root_list)
for ipickle in range(nb_pickle):
    pickle_root = 'pickle/' + pickle_root_list[ipickle]

    [duration_simu, nb_interval, nb_seconds_since_start_pid_concatenate_arr, distance_lvlh_pid_concantenate_arr, nb_seconds_since_start_pid_average_concatenate_arr, \
                 distance_lvlh_pid_average_concantenate_arr, nb_seconds_since_start_pid_average_mid_concatenate_arr, \
                 distance_lvlh_pid_average_mid_concantenate_arr, distance_lvlh_pid_amplitude_mid_concantenate_arr, ecc_average_mid_concantenate_arr, \
                 ecc_obs_average_mid_concantenate_arr, localtime_spock_ok_pid_concatenate, phase_spock_ok_pid_concatenate_arr, argper_average_mid_concantenate_arr, \
                     index_period_spock_concatenate_arr, argper_spock_ok_pid_concatenate_arr,\
                 ecc_ave_conc,ecc_obs_ave_conc,localtime_per,longitude_per,latitude_per,nb_seconds_ave_conc_arr]= pickle.load(open(pickle_root + ".pickle"))




#     [duration_simu, nb_interval, nb_seconds_since_start_pid_concatenate_arr, distance_lvlh_pid_concantenate_arr, nb_seconds_since_start_pid_average_concatenate_arr, \
#          distance_lvlh_pid_average_concantenate_arr, nb_seconds_since_start_pid_average_mid_concatenate_arr, distance_lvlh_pid_average_mid_concantenate_arr,\
#          distance_lvlh_pid_amplitude_mid_concantenate_arr, ecc_average_mid_concantenate_arr, ecc_obs_average_mid_concantenate_arr] = pickle.load(open(pickle_root + ".pickle"))
    
    label_temp = pickle_root_list[ipickle].replace("localtime_", "")
#     if 'equator' in pickle_root_list[ipickle]:
#         label = 'zenith'#'midnight'#'perigee'
#     elif 'pole' in pickle_root_list[ipickle]:
#         label = 'nadir' #'noon'
#     elif 'highamp_pole' in pickle_root_list[ipickle]:
#         label = 'highamp_apogee'
#     elif 'mid' in pickle_root_list[ipickle]:
#         label = '210 deg local time'
    label = pickle_root_list[ipickle]                                                                                                                                           

    if ipickle == 0:
        nb_interval_previous = nb_interval
#     if nb_interval != nb_interval_previous:
#         print "***! The number of interval of all runs has to be the same. The program will stop. !***"; raise Exception;

    if toplot == 'raw':
        ax.plot(nb_seconds_since_start_pid_concatenate_arr/3600., distance_lvlh_pid_concantenate_arr * 1000., linewidth = 2,color = 'b', alpha = 0.3)
        ax.plot(nb_seconds_since_start_pid_average_concatenate_arr/3600., distance_lvlh_pid_average_concantenate_arr * 1000., linewidth = 2, color = 'magenta', linestyle = 'dashed')
    #     ax.plot(nb_seconds_since_start_pid_average_mid_concatenate_arr/3600., distance_lvlh_pid_average_mid_concantenate_arr * 1000., linewidth = 2, color = 'magenta')
    #     ax.plot(nb_seconds_since_start_pid_average_mid_concatenate_arr/3600., ( distance_lvlh_pid_amplitude_mid_concantenate_arr + distance_lvlh_pid_average_mid_concantenate_arr )* 1000., linewidth = 2, color = 'red')
    #     ax.scatter(nb_seconds_since_start_pid_average_mid_concatenate_arr/3600., ( distance_lvlh_pid_amplitude_mid_concantenate_arr + distance_lvlh_pid_average_mid_concantenate_arr )* 1000., s= 500, marker = '.', color = 'red')
    elif toplot == 'amplitude':
        ax.plot(nb_seconds_since_start_pid_average_mid_concatenate_arr/3600., ( distance_lvlh_pid_amplitude_mid_concantenate_arr )* 1000., linewidth = 2, color = color_arr[ipickle], label = label)

    if ipickle == 0:
        pickle_root_concatenate = pickle_root_list[ipickle]
    else:
        pickle_root_concatenate = pickle_root_concatenate + '_+_' +  pickle_root_list[ipickle]
ax.margins(0,0)

# ax.plot([0, duration_simu], [0,0], linestyle = 'dashed', linewidth = 2, color = 'black')
ax.set_xlim([0, 198.]); #ax.set_ylim([-20, 20]) 
# ax.text(duration_simu/2., -200, 'SpOCK in front -> need rho_control < 0', horizontalalignment = 'center', verticalalignment = 'bottom', fontsize = fontsize_plot, weight = 'bold')
# ax.text(duration_simu/2., 1200, 'SpOCK behind -> need rho_control > 0', horizontalalignment = 'center', verticalalignment = 'top', fontsize = fontsize_plot, weight = 'bold')
ax.margins(0,0)
legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="", fontsize = fontsize_plot)


if toplot == 'amplitude':
    fig_save_name = 'fig/all_amplitude_' + pickle_root_concatenate + '_nbinter' + str(nb_interval) + ".pdf"
    y_label = 'Amplitude ocillations (m)'
if toplot == 'raw':
    fig_save_name = 'fig/all_raw_' + pickle_root_concatenate + '_nbinter' + str(nb_interval) + ".pdf"
    y_label = 'Distance (m)'
    ax.set_ylim([-200, 1200])
    ax.text(0.5,0.98,label.title(),fontsize = fontsize_plot, weight = 'bold', color = 'k', transform = ax.transAxes, horizontalalignment = 'center', verticalalignment = 'top')
ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

raise Exception


fig_title = ''#local time of perigee vs delta eccentricity
x_label = 'Difference in eccentricity' 

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])


ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold


pickle_root_concatenate = ''
nb_pickle = len(pickle_root_list)
ipickle = 0
pickle_root = 'pickle/' + pickle_root_list[ipickle]


[duration_simu, nb_interval, nb_seconds_since_start_pid_concatenate_arr, distance_lvlh_pid_concantenate_arr, nb_seconds_since_start_pid_average_concatenate_arr, \
             distance_lvlh_pid_average_concantenate_arr, nb_seconds_since_start_pid_average_mid_concatenate_arr, \
             distance_lvlh_pid_average_mid_concantenate_arr, distance_lvlh_pid_amplitude_mid_concantenate_arr, ecc_average_mid_concantenate_arr, \
             ecc_obs_average_mid_concantenate_arr, localtime_spock_ok_pid_concatenate, phase_spock_ok_pid_concatenate_arr, argper_average_mid_concantenate_arr, \
     index_period_spock_concatenate_arr, argper_spock_ok_pid_concatenate_arr,\
                 ecc_ave_conc,ecc_obs_ave_conc,localtime_per,longitude_per,latitude_per,nb_seconds_ave_conc_arr]= pickle.load(open(pickle_root + ".pickle"))


if 'equator' in pickle_root_list[ipickle]:
    label = 'midnight'#'zenith'#'midnight'#'perigee'
elif 'pole' in pickle_root_list[ipickle]:
    label = 'noon#''nadir' #'noon'
elif 'highamp_pole' in pickle_root_list[ipickle]:
    label = 'highamp_apogee'
elif 'mid' in pickle_root_list[ipickle]:
    label = 'mid'#'210 deg local time'



if ipickle == 0:
    nb_interval_previous = nb_interval
if nb_interval != nb_interval_previous:
    print "***! The number of interval of all runs has to be the same. The program will stop. !***"; raise Exception;

x_axis = ecc_ave_conc-ecc_obs_ave_conc
y_axis = localtime_per
ax.set_ylim([np.min(y_axis), np.max(y_axis)])
ax.set_xlim([np.min(x_axis), np.max(x_axis)])
ax.set_ylabel('Local time of perigee', weight = 'bold', fontsize  = fontsize_plot)
#fig_save_name = 'fig/delta_eccentricity_vs_amplitude' + pickle_root_concatenate + '_nbinter' + str(nb_interval) + ".pdf"
fig_save_name = 'fig/localtime_per_vs_delta_ecc' + pickle_root_concatenate + '_nbinter' + str(nb_interval) + ".pdf"
y_label = 'Eccentricity' #'Eccentricity'#'Difference in eccentricity'

ax.scatter(x_axis, y_axis , linewidth = 2, color = color_arr[ipickle], label = label)

fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

raise Exception

fig_title = ''#eccentriciy or delta ecc vs  amplitutde oscilaltion 
x_label = 'Amplitude oscillations (m)' 

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])


ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold


pickle_root_concatenate = ''
nb_pickle = len(pickle_root_list)
ipickle = 0
pickle_root = 'pickle/' + pickle_root_list[ipickle]


[duration_simu, nb_interval, nb_seconds_since_start_pid_concatenate_arr, distance_lvlh_pid_concantenate_arr, nb_seconds_since_start_pid_average_concatenate_arr, \
             distance_lvlh_pid_average_concantenate_arr, nb_seconds_since_start_pid_average_mid_concatenate_arr, \
             distance_lvlh_pid_average_mid_concantenate_arr, distance_lvlh_pid_amplitude_mid_concantenate_arr, ecc_average_mid_concantenate_arr, \
             ecc_obs_average_mid_concantenate_arr, localtime_spock_ok_pid_concatenate, phase_spock_ok_pid_concatenate_arr, argper_average_mid_concantenate_arr, \
                 index_period_spock_concatenate_arr, argper_spock_ok_pid_concatenate_arr,\
                 ecc_ave_conc,ecc_obs_ave_conc,localtime_per,longitude_per,latitude_per,nb_seconds_ave_conc_arr]= pickle.load(open(pickle_root + ".pickle"))


if 'equator' in pickle_root_list[ipickle]:
    label = 'midnight'#'zenith'#'midnight'#'perigee'
elif 'pole' in pickle_root_list[ipickle]:
    label = 'noon#''nadir' #'noon'
elif 'highamp_pole' in pickle_root_list[ipickle]:
    label = 'highamp_apogee'
elif 'mid' in pickle_root_list[ipickle]:
    label = 'mid'#'210 deg local time'



if ipickle == 0:
    nb_interval_previous = nb_interval
if nb_interval != nb_interval_previous:
    print "***! The number of interval of all runs has to be the same. The program will stop. !***"; raise Exception;

x_axis = ( distance_lvlh_pid_amplitude_mid_concantenate_arr )* 1000.
y_axis = ecc_average_mid_concantenate_arr - ecc_obs_average_mid_concantenate_arr#ecc_average_mid_concantenate_arr
ax.set_ylim([np.min(y_axis), np.max(y_axis)])
ax.set_ylabel('Difference in eccentricity', weight = 'bold', fontsize  = fontsize_plot)
#fig_save_name = 'fig/delta_eccentricity_vs_amplitude' + pickle_root_concatenate + '_nbinter' + str(nb_interval) + ".pdf"
fig_save_name = 'fig/delta_eccentricity_vs_amplitude' + pickle_root_concatenate + '_nbinter' + str(nb_interval) + ".pdf"
y_label = 'Eccentricity' #'Eccentricity'#'Difference in eccentricity'

ax.scatter(x_axis, y_axis , linewidth = 2, color = color_arr[ipickle], label = label)

fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


raise Exception
fig_title = ''#eccentriciy and  amplitutde oscilaltion vs time
x_label = 'Time (hours)'#'Amplitude oscillations (m)' 

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])


ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold


pickle_root_concatenate = ''
nb_pickle = len(pickle_root_list)
ipickle = 0
pickle_root = 'pickle/' + pickle_root_list[ipickle]


[duration_simu, nb_interval, nb_seconds_since_start_pid_concatenate_arr, distance_lvlh_pid_concantenate_arr, nb_seconds_since_start_pid_average_concatenate_arr, \
             distance_lvlh_pid_average_concantenate_arr, nb_seconds_since_start_pid_average_mid_concatenate_arr, \
             distance_lvlh_pid_average_mid_concantenate_arr, distance_lvlh_pid_amplitude_mid_concantenate_arr, ecc_average_mid_concantenate_arr, \
             ecc_obs_average_mid_concantenate_arr, localtime_spock_ok_pid_concatenate, phase_spock_ok_pid_concatenate_arr, argper_average_mid_concantenate_arr, \
                 index_period_spock_concatenate_arr, argper_spock_ok_pid_concatenate_arr,\
                 ecc_ave_conc,ecc_obs_ave_conc,localtime_per,longitude_per,latitude_per,nb_seconds_ave_conc_arr]= pickle.load(open(pickle_root + ".pickle"))


if 'equator' in pickle_root_list[ipickle]:
    label = 'midnight'#'zenith'#'midnight'#'perigee'
elif 'pole' in pickle_root_list[ipickle]:
    label = 'noon#''nadir' #'noon'
elif 'highamp_pole' in pickle_root_list[ipickle]:
    label = 'highamp_apogee'
elif 'mid' in pickle_root_list[ipickle]:
    label = 'mid'#'210 deg local time'



if ipickle == 0:
    nb_interval_previous = nb_interval
if nb_interval != nb_interval_previous:
    print "***! The number of interval of all runs has to be the same. The program will stop. !***"; raise Exception;

x_axis = nb_seconds_since_start_pid_average_mid_concatenate_arr/3600.#( distance_lvlh_pid_amplitude_mid_concantenate_arr )* 1000.
y_axis = ecc_average_mid_concantenate_arr# - ecc_obs_average_mid_concantenate_arr#ecc_average_mid_concantenate_arr
#y_axis = -y_axis*() + y_axis[0] + ( distance_lvlh_pid_amplitude_mid_concantenate_arr[0] )* 1000.
ax.set_ylim([np.min(y_axis), np.max(y_axis)])
ax.tick_params('y', colors='blue')

ax.set_ylabel('Eccentricity', weight = 'bold', fontsize  = fontsize_plot, color = 'blue')
#ax.plot(nb_seconds_since_start_pid_average_mid_concatenate_arr/3600., ecc_average_mid_concantenate_arr - ecc_obs_average_mid_concantenate_arr, linewidth = 2, color = color_arr[ipickle], label = label)
fig_save_name = 'fig/eccentricity_and_amplitude' + pickle_root_concatenate + '_nbinter' + str(nb_interval) + ".pdf"
y_label = 'Eccentricity' #'Eccentricity'#'Difference in eccentricity'

#ax.scatter(x_axis, y_axis , linewidth = 2, color = color_arr[ipickle], label = label)

ax.plot(x_axis, y_axis , linewidth = 2, color = 'blue')
ax.margins(0,0)
ax2 = ax.twinx()
ax2.set_ylabel('Amplitude of oscillations (m)', weight = 'bold', fontsize  = fontsize_plot, color = 'red')

[i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
ax2.tick_params('y', colors='red')

y_axis_dist = ( distance_lvlh_pid_amplitude_mid_concantenate_arr )* 1000. 
ax2.plot(x_axis, y_axis_dist , linewidth = 2, color = 'red')

# pickle_root_concatenate = pickle_root_list[ipickle]

ax2.set_ylim([np.min(y_axis_dist), np.max(y_axis_dist)])
ax2.set_xlim([np.min(x_axis), np.max(x_axis)])
#ax.margins(0,0)

# ax.plot([0, duration_simu], [0,0], linestyle = 'dashed', linewidth = 2, color = 'black')
# ax.set_xlim([0, duration_simu]); #ax.set_ylim([-20, 20]) 
# ax.text(duration_simu/2., -200, 'SpOCK in front -> need rho_control < 0', horizontalalignment = 'center', verticalalignment = 'bottom', fontsize = fontsize_plot, weight = 'bold')
# ax.text(duration_simu/2., 1200, 'SpOCK behind -> need rho_control > 0', horizontalalignment = 'center', verticalalignment = 'top', fontsize = fontsize_plot, weight = 'bold')
ax2.margins(0,0)
#legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="", fontsize = fontsize_plot)


ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

raise Exception


######
fig_title = ''#Difference in eccentricy SpOCK - observations
x_label = 'Time (hours)' 

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])


ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold


pickle_root_concatenate = ''
nb_pickle = len(pickle_root_list)
for ipickle in range(nb_pickle):
    pickle_root = 'pickle/' + pickle_root_list[ipickle]

    [duration_simu, nb_interval, nb_seconds_since_start_pid_concatenate_arr, distance_lvlh_pid_concantenate_arr, nb_seconds_since_start_pid_average_concatenate_arr, \
                 distance_lvlh_pid_average_concantenate_arr, nb_seconds_since_start_pid_average_mid_concatenate_arr, \
                 distance_lvlh_pid_average_mid_concantenate_arr, distance_lvlh_pid_amplitude_mid_concantenate_arr, ecc_average_mid_concantenate_arr, \
                 ecc_obs_average_mid_concantenate_arr, localtime_spock_ok_pid_concatenate, phase_spock_ok_pid_concatenate_arr, argper_average_mid_concantenate_arr, \
                     index_period_spock_concatenate_arr, argper_spock_ok_pid_concatenate_arr,\
                 ecc_ave_conc,ecc_obs_ave_conc,localtime_per,longitude_per,latitude_per,nb_seconds_ave_conc_arr]= pickle.load(open(pickle_root + ".pickle")) 
    
    label_temp = pickle_root_list[ipickle].replace("localtime_", "")
    if 'equator' in pickle_root_list[ipickle]:
        label = 'zenith'#'midnight'#'perigee'
    elif 'pole' in pickle_root_list[ipickle]:
        label = 'nadir' #'noon' 
    elif 'highamp_pole' in pickle_root_list[ipickle]:
        label = 'highamp_apogee'
    elif 'mid' in pickle_root_list[ipickle]:
        label = '210 deg local time' 




    if ipickle == 0:
        nb_interval_previous = nb_interval
    if nb_interval != nb_interval_previous:
        print "***! The number of interval of all runs has to be the same. The program will stop. !***"; raise Exception;

    #ax.plot(nb_seconds_since_start_pid_average_mid_concatenate_arr/3600., ecc_average_mid_concantenate_arr - ecc_obs_average_mid_concantenate_arr, linewidth = 2, color = color_arr[ipickle], label = label)
    ax.plot(nb_seconds_ave_conc_arr[:-1]/3600., ecc_ave_conc-ecc_obs_ave_conc, linewidth = 2, color = color_arr[ipickle], label = label)

#     ax.plot(nb_seconds_since_start_pid_average_mid_concatenate_arr/3600., ecc_average_mid_concantenate_arr , linewidth = 2, color = color_arr[ipickle], label = label)
#     if ipickle == 0: # observations same for each run
#         ax.plot(nb_seconds_since_start_pid_average_mid_concatenate_arr/3600., ecc_obs_average_mid_concantenate_arr, linewidth = 2, color = 'black', label = 'Observations')

    if ipickle == 0:
        pickle_root_concatenate = pickle_root_list[ipickle]
    else:
        pickle_root_concatenate = pickle_root_concatenate + '_+_' +  pickle_root_list[ipickle]
ax.margins(0,0)

# ax.plot([0, duration_simu], [0,0], linestyle = 'dashed', linewidth = 2, color = 'black')
# ax.set_xlim([0, duration_simu]); #ax.set_ylim([-20, 20]) 
# ax.text(duration_simu/2., -200, 'SpOCK in front -> need rho_control < 0', horizontalalignment = 'center', verticalalignment = 'bottom', fontsize = fontsize_plot, weight = 'bold')
# ax.text(duration_simu/2., 1200, 'SpOCK behind -> need rho_control > 0', horizontalalignment = 'center', verticalalignment = 'top', fontsize = fontsize_plot, weight = 'bold')
ax.margins(0,0)
legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="", fontsize = fontsize_plot)

fig_save_name = 'fig/all_delta_eccentricity_' + pickle_root_concatenate + '_nbinter' + str(nb_interval) + "again.pdf"
y_label = 'Difference in eccentricity' #'Eccentricity'#'Difference in eccentricity'

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  



raise Exception
# FTT
plt.close('all')
fig,ax = plt.subplots()
signal = distance_lvlh_pid_amplitude_mid_concantenate_arr*1000
# dt_sample about 95 minutes. actually not always the same (sometimes differ by 10 minutes or a couple of minutes)
dt_sample = nb_seconds_since_start_pid_average_mid_concatenate_arr[-1] - nb_seconds_since_start_pid_average_mid_concatenate_arr[-2] 
n = len(signal)
sp = np.fft.fft(signal)
freq = np.fft.fftfreq(n, dt_sample)
freq_pos_index = np.where(freq > 0)[0]
freq_pos_ok = freq[freq_pos_index]
sp_pos = sp[freq_pos_index]
sp_pos_abs = abs(sp_pos) # sp_pos is  a complex number, only the norm is relevant
ax.plot(signal)

fig,axfft = plt.subplots()
axfft.plot(1/freq_pos_ok/3600., sp_pos_abs) # sp.real





#plt.close('all')
fig,ax = plt.subplots()
signal = distance_lvlh_pid_amplitude_mid_concantenate_arr*1000
# dt_sample about 95 minutes. actually not always the same (sometimes differ by 10 minutes or a couple of minutes)
dt_sample = 1#nb_seconds_since_start_pid_average_mid_concatenate_arr[-1] - nb_seconds_since_start_pid_average_mid_concatenate_arr[-2] 
n = len(signal)
sp = np.fft.fft(signal)
freq = np.fft.fftfreq(n, dt_sample)
freq_pos_index = np.where(freq > 0)[0]
freq_pos = freq[freq_pos_index]
sp_pos = sp[freq_pos_index]
sp_pos_abs = abs(sp_pos) # sp_pos is  a complex number, only the norm is relevant
ax.plot(signal)

fig,axfft = plt.subplots()
#axfft.plot(1/freq_pos/3600., sp_pos_abs) # sp.real
axfft.plot(nb_seconds_since_start_pid_average_mid_concatenate_arr[(1/freq_pos).astype(int)-1]/3600., sp_pos_abs) # sp.real




raise Exception

plt.close('all')
t = np.arange(0, 6*np.pi, 0.1)
signal = np.sin(t)
sp = np.fft.fft(signal)
freq = np.fft.fftfreq(t.shape[-1])#, 0.1)
plt.plot(signal)

fig,ax = plt.subplots()
ax.plot(freq, abs(sp)) # sp.real



Fs = 150.0;  # sampling rate
Ts = 1.0/Fs; # sampling interval
t = np.arange(0,1,Ts) # time vector

ff = 5;   # frequency of the signal
y = np.sin(2*np.pi*ff*t)

n = len(y) # length of the signal
k = np.arange(n)
T = n/Fs
frq = k/T # two sides frequency range
frq = frq[range(n/2)] # one side frequency range

Y = np.fft.fft(y)/n # fft computing and normalization
Y = Y[range(n/2)]

fig, ax = plt.subplots(2, 1)
ax[0].plot(t,y)
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Amplitude')
ax[1].plot(frq,abs(Y),'r') # plotting the spectrum
ax[1].set_xlabel('Freq (Hz)')
ax[1].set_ylabel('|Y(freq)|')
