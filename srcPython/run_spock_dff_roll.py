# This scri[t runs SpOCK from date_start to date_spock with a roll angle varying from roll_start to roll_stop with a angle step of drokk
islin = 0 # 1 if running from linux desktop
import ipdb
import matplotlib
from datetime import datetime, timedelta
import numpy as np
import os
import sys
if islin == 1:
    sys.path.append("/home/cbv/Code/spock/srcPython")
else:
    sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/srcPython")
from spock_main_input import *
from orbit_average import *
from os import listdir
from read_input_file import *
from read_output_file import *
from cygnss_read_spock_spec import *
from cygnss_read_spock_spec_bin import *
from netCDF4 import Dataset
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from matplotlib import pyplot as plt
#from ecef2eci import *
#from eci_to_lvlh import *
from ecef_to_lvlh import *
import pickle
from cygnss_name_to_norad_id import *
import os.path
import matplotlib.colors as colors
import matplotlib.cm as cmx
from cygnss_convert_rcg import *
plt.ion()

date_start = '2018-08-20T00:00:00'
date_stop = '2018-08-20T23:59:59'
cygfm = 1


roll_start = 0. # in deg
roll_end = -22.
droll = -2.
roll_arr = np.arange(roll_start, roll_end + droll, droll)
dir_run_spock = '/Users/cbv/cygnss/valsift_temp'#'/Users/cbv/cygnss/sift_temp' # '.' # no slash
nb_roll = len(roll_arr)
gain_spock_all_roll = []

# for iroll in range(0, nb_roll):
#     print iroll, nb_roll - 1
#     roll = (int)(roll_arr[iroll])
#     #input_filename = "roll" + str((int)(roll)) + "_" date_start.replace(":","_").replace("-","") + '_' + date_stop.replace(":","_").replace("-","") + ".txt"
#     # os.chdir(dir_run_spock )
#     # os.system("spock_cygnss_spec_parallel_tell_roll.py " + date_start + " " + date_stop + " spec " + str(roll))
#     # os.chdir("/Users/cbv/Google Drive/Work/PhD/Research/Code/cygnss/validate_sift")
#     spock_input_filename = dir_run_spock + '/roll' + str(int(roll))   + '_spock_spec_start_' + date_start.replace(":", "_") + '_end_' + date_stop.replace(":", "_") + '.txt'

#         # Read specular positio computed by SpOCK
#     var_in, var_in_order = read_input_file(spock_input_filename)
#     dt_spock_output = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')]; 
#     output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
#     output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')];
#     gps_name_list_spock = var_in[find_in_read_input_order_variables(var_in_order, 'gps_name')];
#     cygfm_to_spock_nb = [4,3,8,2,1,6,7,5] # ['41884', '41885', '41886', '41887', '41888', '41889', '41890', '41891']
#     isc =  cygfm_to_spock_nb[cygfm-1] - 1
#     spec_spock_filename = dir_run_spock + "/" + output_file_path_list[isc] + "specular_" + output_file_name_list[isc]
#     #spec_spock_filename = spec_spock_filename.replace(".txt", "_6SPs.txt")  
#     print "Reading SpOCK specular files..."

#     data_spec = cygnss_read_spock_spec_bin(spec_spock_filename.replace('.txt','.bin'), gps_name_list_spock, dt_spock_output, 1) 
#     date_spock = data_spec[0]; lon_spock = data_spec[1]; lat_spock = data_spec[2]; gain_spock = data_spec[3]; gps_spock = data_spec[4]; normpower_spock = data_spec[5]; x_cyg_spock = data_spec[6]; y_cyg_spock = data_spec[7]; z_cyg_spock = data_spec[8]; x_gps_spock = data_spec[9]; y_gps_spock = data_spec[10]; z_gps_spock = data_spec[11];  x_spec_spock = data_spec[12]; y_spec_spock = data_spec[13]; z_spec_spock = data_spec[14]; nb_spec_spock = data_spec[15];  el_spec_spock = data_spec[16]; az_spec_spock = data_spec[17]; el_gps_from_cyg_spock = data_spec[18];  el_spec_not_int_spock = data_spec[19]; az_spec_not_int_spock = data_spec[20]

#     gain_spock_all_roll.append(gain_spock)



#pickle.dump(gain_spock_all_roll, open("gain_spock_all_roll.pickle", "w"))
gain_spock_all_roll = pickle.load(open("gain_spock_all_roll.pickle")) # this is from the different roll simulations
gain_spock_comp_small_scale,fom_netcdf_small_scale,fom_netcdf_temp_small_scale = pickle.load(open("gain_spock,fom_netcdf,fom_netcdf_temp_from_newer_validate_sift.pickle"))# this is from the run with newer_validate_sift.py where I ccmopare SpOCK to netcdf with a roll of -22.



nn = gain_spock_comp_small_scale.shape[0]
gain_spock_comp = np.zeros([nn, 4])
for ispec in range(4):
    gain_spock_comp[:, ispec] = cygnss_convert_rcg(gain_spock_comp_small_scale[:, ispec], 'small')


nn = fom_netcdf_temp_small_scale.shape[0]
fom_netcdf_temp = np.zeros([nn, 4])
for ispec in range(4):
    fom_netcdf_temp[:, ispec] = cygnss_convert_rcg(fom_netcdf_temp_small_scale[:, ispec], 'small')


    


# # Now just run FM01 for these different rolls with solar power on and compute power
# os.chdir(dir_run_spock )
# iroll = 0
# power_orbit_averaged_first = []
# for iroll in range(0, nb_roll):
#     print iroll, nb_roll - 1
#     roll = (int)(roll_arr[iroll])
#     spock_input_filename = dir_run_spock + '/fm01_solarPower_roll' + str(int(roll))   + '_spock_spec_start_' + date_start.replace(":", "_") + '_end_' + date_stop.replace(":", "_") + '.txt'
#     main_input_filename = 'fm01_solarPower_roll' + str(int(roll))   + '_spock_spec_start_' + date_start.replace(":", "_") + '_end_' + date_stop.replace(":", "_") + '.txt'
#     # Write the SpOCK main input file 
#     if iroll == 0:
#         # for TIME section
#         date_start = date_start 
#         date_end = date_stop
#         dt = 10 # !!!!!!!change back to 10
#         # for SPACECRAFT section
#         nb_sc = 1
#         gps_tle = '0'
#         mass = 29
#         geometry_filename = "cygnss_geometry_2016_acco09.txt"
#         # for ORBIT section
#         tle_filename = 'fm01_2018-08-20.txt' #!!!!!is that really the one you want?
#         # for FORCES section
#         gravity_order = 8 #!!!!!!!!change back to 4
#         forces = 'drag sun_gravity moon_gravity'
#         density_mode = 'dynamic'
#         # for OUTPUT section
#         name_output = 'spock_out/out' # if spock/out then conflicts in the SpOCK code (to create the output folder) with the executable spock
#         dt_output = 60
#         # #for SPICE section put in makeall.sh
#         spice_path = "/Users/cbv/cspice/data"

#     ## Create main input file
#     spock_main_input(
#         main_input_filename,
#         # for TIME section
#         date_start,
#         date_end,
#         dt,
#         # for SPACECRAFT section
#         nb_sc,
#         gps_tle,
#         mass,
#         geometry_filename, 
#         # for ORBIT section
#         tle_filename,
#         # for FORCES section
#         gravity_order,
#         forces,
#         density_mode,
#         # for OUTPUT section
#         name_output,
#         dt_output,
#         # for ATTITUDE section
#         "(0;"  + str(int(roll)) + ";0) (0;0;0)",#"(0;-22;0) (0;0;0)",#"nadir", #!!!!!
#         # for GROUNDS_STATIONS section
#         "0",#"my_ground_stations.txt"
#          # for SPICE section
#          spice_path,
#          # for DENSITY_MOD section
#     1
#         )

#     # Run SpOCK 
#     ## CYGNSS and GPS positions (no GPS position is predicted if the user did not call 'spec' at a 3rd argument of this script)

#     print "mpirun -np 4 spock_dev " + main_input_filename#  
#     os.system("mpirun -np 4 spock_dev " + main_input_filename )

#         # Read specular positio computed by SpOCK
#     var_in, var_in_order = read_input_file(spock_input_filename)
#     dt_spock_output = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')]; 
#     output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
#     output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')];
#     isc = 0
#     var_to_read = ["power", "latitude"]
#     var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
#     date = var_out[find_in_read_input_order_variables(var_out_order, 'date')]
#     power = var_out[find_in_read_input_order_variables(var_out_order, 'power')]
#     latitude = var_out[find_in_read_input_order_variables(var_out_order, 'latitude')]
#     power_orbit_averaged, time_averaged, index_time_averaged = orbit_average(power, latitude, date )
#     power_orbit_averaged_first.append( power_orbit_averaged[0] ) # first orbit only (doesn't vary much within a day)
#     print power_orbit_averaged_first
# os.chdir("/Users/cbv/Google Drive/Work/PhD/Research/Code/cygnss/validate_sift")
# # Figure
raise Exception


# Compare SpOCK and netcdf distirbutions on scale 0 to 220
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25      
ratio_fig_size = 8./3
fig_title = ''#Accuracy VS RCG
y_label = '% samples'
x_label = 'RCG'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle('SpOCK and Netcdf RCG distributions', y = 1,fontsize = (int)(fontsize_plot*1.), weight = 'normal')
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                      
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)

to_hist_spock = gain_spock_comp
values, base = np.histogram(to_hist_spock,  bins = 16)
bin_array = ( base[:-1] + np.roll(base,-1)[:-1] ) /2.
binsize_actual = bin_array[1] - bin_array[0]
hist_spock = values * 100. /(4* len(to_hist_spock))

to_hist_netcdf = fom_netcdf_temp
values, base = np.histogram(to_hist_netcdf, bins = 16)
bin_array = ( base[:-1] + np.roll(base,-1)[:-1] ) /2.
binsize_actual = bin_array[1] - bin_array[0]
hist_netcdf = values * 100. /(4* len(to_hist_netcdf))
#CDF
ax_cdf = fig.add_subplot(gs[0, 0])
ax_cdf.set_title('Cumulative Distribution Function', weight = 'normal', fontsize  = fontsize_plot)

ax_cdf.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax_cdf.spines.itervalues()] # change the width of the frame of the figure


cdf_spock = np.zeros([16])
cdf_spock[15] = hist_spock[15]
for i in np.arange(14,-1,-1):
    cdf_spock[i] = cdf_spock[i+1] + hist_spock[i] 

ax_cdf.plot(bin_array, cdf_spock, linewidth = 2, color = 'blue', label = 'SpOCK')

cdf_netcdf = np.zeros([16])
cdf_netcdf[15] = hist_netcdf[15]
for i in np.arange(14,-1,-1):
    cdf_netcdf[i] = cdf_netcdf[i+1] + hist_netcdf[i] 


ax_cdf.plot(bin_array, cdf_netcdf, linewidth = 2, color = 'red', label = 'Netcdf')

#ax_cdf.xaxis.set_ticks(np.arange(0,16,1))
ax_cdf.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)

ax_cdf.set_ylim([0, 100])
ax_cdf.margins(0,0)


legend = ax_cdf.legend(loc='lower right', bbox_to_anchor=(1, 0), numpoints = 1,  title="", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))

plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           

fig_save_name = 'cdf_spock_netcdf_scale_sclae_0_220.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')


# Compare SpOCK and netcdf distirbutions
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 9      
ratio_fig_size = 8./3
fig_title = ''#Accuracy VS RCG
y_label = '% samples'
x_label = 'RCG'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle('SpOCK and Netcdf RCG distributions', y = 1,fontsize = (int)(fontsize_plot*1.), weight = 'normal')
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                      
gs = gridspec.GridSpec(1, 3)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])
ax.set_title('PDF SpOCK', weight = 'normal', fontsize  = fontsize_plot)
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
to_hist_spock = gain_spock_comp_small_scale
values, base = np.histogram(to_hist_spock, range = [-0.5, 15.5], bins = 16)
bin_array = ( base[:-1] + np.roll(base,-1)[:-1] ) /2.
binsize_actual = bin_array[1] - bin_array[0]
hist_spock = values * 100. /(4* len(to_hist_spock))
ax.bar(bin_array, hist_spock, binsize_actual, edgecolor = 'black', color = 'blue', alpha = 1, label = 'SpOCK')
ax.xaxis.set_ticks(np.arange(0,16,1))
bin_array_large_scale = cygnss_convert_rcg(bin_array, "small")
ticklabel_string = []
for itick in range(16):
    ticklabel_string.append(str((int)(bin_array[itick])) + '\n' + format(bin_array_large_scale[itick], ".1f"))
ax.xaxis.set_ticklabels(ticklabel_string, fontsize = fontsize_plot) 
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
ax.set_ylim([0, 25])
ax.margins(0,0)
# Netcdf
ax = fig.add_subplot(gs[0, 1])
ax.set_title('PDF netcdf', weight = 'normal', fontsize  = fontsize_plot)
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
to_hist_netcdf = fom_netcdf_temp_small_scale
values, base = np.histogram(to_hist_netcdf, range = [-0.5, 15.5], bins = 16)
bin_array = ( base[:-1] + np.roll(base,-1)[:-1] ) /2.
binsize_actual = bin_array[1] - bin_array[0]
hist_netcdf = values * 100. /(4* len(to_hist_netcdf))
ax.bar(bin_array, hist_netcdf, binsize_actual, edgecolor = 'black', color = 'red', alpha = 1, label = 'Netcdf')
ax.xaxis.set_ticks(np.arange(0,16,1))
bin_array_large_scale = cygnss_convert_rcg(bin_array, "small")
ticklabel_string = []
for itick in range(16):
    ticklabel_string.append(str((int)(bin_array[itick])) + '\n' + format(bin_array_large_scale[itick], ".1f"))
ax.xaxis.set_ticklabels(ticklabel_string, fontsize = fontsize_plot) 
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
ax.set_ylim([0, 25])
ax.margins(0,0)
#CDF
ax_cdf = fig.add_subplot(gs[0, 2])
ax_cdf.set_title('Cumulative Distribution Function', weight = 'normal', fontsize  = fontsize_plot)
ax_cdf.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax_cdf.spines.itervalues()] # change the width of the frame of the figure
cdf_spock = np.zeros([16])
cdf_spock[15] = hist_spock[15]
for i in np.arange(14,-1,-1):
    cdf_spock[i] = cdf_spock[i+1] + hist_spock[i] 

ax_cdf.plot(bin_array, cdf_spock, linewidth = 2, color = 'blue', label = 'SpOCK')

cdf_netcdf = np.zeros([16])
cdf_netcdf[15] = hist_netcdf[15]
for i in np.arange(14,-1,-1):
    cdf_netcdf[i] = cdf_netcdf[i+1] + hist_netcdf[i] 
ax_cdf.plot(bin_array, cdf_netcdf, linewidth = 2, color = 'red', label = 'Netcdf')
ax_cdf.xaxis.set_ticks(np.arange(0,16,1))
bin_array_large_scale = cygnss_convert_rcg(bin_array, "small")
ticklabel_string = []
for itick in range(16):
    ticklabel_string.append(str((int)(bin_array[itick])) + '\n' + format(bin_array_large_scale[itick], ".1f"))
ax_cdf.xaxis.set_ticklabels(ticklabel_string, fontsize = fontsize_plot) 
ax_cdf.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
ax_cdf.set_ylim([0, 100])
ax_cdf.margins(0,0)
legend = ax_cdf.legend(loc='lower right', bbox_to_anchor=(1, 0), numpoints = 1,  title="", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           
fig_save_name = 'pdf_cdf_spock_netcdf_scale.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')






## pLOT distirbutions at different roll
height_fig = 11.  # the width is calculated as height_fig * 4/3.                                                                             
fontsize_plot = 25      
ratio_fig_size = 6./3
fig_title = ''#Accuracy VS RCG
y_label = '% samples'
x_label = 'RCG'

for iroll in range(nb_roll):
#iroll = 0
    roll = (int)(roll_arr[iroll])
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('Spacecraft roll: ' + str(roll) + u'\N{DEGREE SIGN}', y = 1,fontsize = (int)(fontsize_plot*1.), weight = 'normal')
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                      
    gs = gridspec.GridSpec(1, 2)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
    ax = fig.add_subplot(gs[0, 0])
    ax.set_title('Probability Distribution Function', weight = 'normal', fontsize  = fontsize_plot)
    ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure

    to_hist = gain_spock_all_roll[iroll]
    values, base = np.histogram(to_hist, range = [-0.5, 15.5], bins = 16)
    bin_array = ( base[:-1] + np.roll(base,-1)[:-1] ) /2.
    binsize_actual = bin_array[1] - bin_array[0]
    hist = values * 100. /(4* len(to_hist))
    ax.bar(bin_array, hist, binsize_actual, edgecolor = 'black')

    ax.xaxis.set_ticks(np.arange(0,16,1))
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax.set_ylim([0, 25])
    ax.margins(0,0)
    #CDF
    ax_cdf = fig.add_subplot(gs[0, 1])
    ax_cdf.set_title('Cumulative Distribution Function', weight = 'normal', fontsize  = fontsize_plot)

    ax_cdf.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax_cdf.spines.itervalues()] # change the width of the frame of the figure


    cdf = np.zeros([16])
    cdf[0] = hist[0]
    for i in range(1,16):
        cdf[i] = cdf[i-1] + hist[i] 
    ax_cdf.plot(bin_array, cdf, linewidth = 2, color = 'black')

    ax_cdf.xaxis.set_ticks(np.arange(0,16,1))
    ax_cdf.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)

    ax_cdf.set_ylim([0, 100])
    ax_cdf.margins(0,0)

    plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           

    fig_save_name = 'roll' + str(roll) + '_pdf_cdf.pdf'
    fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')
#plt.show()




## pLOT all CDFs on same plot
# For plots, generate disctinct colors
NCURVES = nb_roll
np.random.seed(101)
curves = [np.random.random(20) for i in range(NCURVES)]
values = range(NCURVES)
jet = cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

height_fig = 11.  # the width is calculated as height_fig * 4/3.                                                                             
fontsize_plot = 25      
ratio_fig_size = 6./3
fig_title = ''#Accuracy VS RCG
y_label = '% samples'
x_label = 'RCG'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle('', y = 1,fontsize = (int)(fontsize_plot*1.), weight = 'normal')
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                      
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax_cdf = fig.add_subplot(gs[0, 0])
ax_cdf.set_title('Cumluative Distribution Function', weight = 'normal', fontsize  = fontsize_plot)
ax_cdf.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax_cdf.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax_cdf.spines.itervalues()] # change the width of the frame of the figure

for iroll in range(nb_roll):
#iroll = 0
    roll = (int)(roll_arr[iroll])

    to_hist = gain_spock_all_roll[iroll]
    values, base = np.histogram(to_hist, range = [-0.5, 15.5], bins = 16)
    bin_array = ( base[:-1] + np.roll(base,-1)[:-1] ) /2.
    binsize_actual = bin_array[1] - bin_array[0]
    hist = values * 100. /(4* len(to_hist))
    cdf = np.zeros([16])
    cdf[15] = hist[15]
    for i in np.arange(14,-1,-1):
        cdf[i] = cdf[i+1] + hist[i] 

    colorVal = scalarMap.to_rgba(iroll)
    ax_cdf.plot(bin_array, cdf, linewidth = 2, color = colorVal, label = str(roll))

ax_cdf.xaxis.set_ticks(np.arange(0,16,1))
bin_array_large_scale = cygnss_convert_rcg(bin_array, "small")
ticklabel_string = []
for itick in range(16):
    ticklabel_string.append(str((int)(bin_array[itick])) + '\n' + format(bin_array_large_scale[itick], ".1f"))
ax_cdf.xaxis.set_ticklabels(ticklabel_string, fontsize = fontsize_plot) 

ax_cdf.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
ax_cdf.set_ylim([0, 100])
ax_cdf.margins(0,0)
legend = ax_cdf.legend(loc='upper right', bbox_to_anchor=(1, 1), numpoints = 1,  title="", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           
fig_save_name = 'all_roll_cdf_scale.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')
#plt.show()



# Plot pwoer VS roll angle
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25      
ratio_fig_size = 6./3
fig_title = ''#Accuracy VS RCG
y_label = 'Power (W)'
x_label = 'Roll angle (degrees)'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle('', y = 0.99,fontsize = (int)(fontsize_plot*1.1), weight = 'normal')
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                      
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])
ax.set_title('Orbit average power VS spacecraft roll angle', weight = 'normal', fontsize  = fontsize_plot)
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)

ax.plot(roll_arr, power_orbit_averaged_first, linewidth = 2, color = 'black')
ax.scatter(roll_arr, power_orbit_averaged_first, linewidth = 2, color = 'black', marker = 'o')

ax.margins(0,0)

plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           

fig_save_name = 'power_vs_roll.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')

raise Exception
