from cadre_read_last_tle import *
from get_prop_dir import *
import matplotlib.gridspec as gridspec
from read_input_file import *
from matplotlib.colors import LogNorm
import pickle
from eci_to_lvlh import *
import sys
import fileinput
import time
import datetime
import numpy as np
from matplotlib import pyplot as plt
import os
import subprocess

plt.ion()
plt.isinteractive()

## NOTE 0: this script plots loads the distributions (created by distance_ensemble_to_main_sc.py) corresponding to different solar panels configurations for CADRE (0, 1, 2, 3, or 4 panels). Each of the distributions corresponds to the distances between each ensemble member and a reference satellite. The reference satellite that we chose is CADRE 4 panels nadir pointing.
## NOTE 1: to use this script, the only parameter you have to set are:
## - the name of the main input file (first non-commented line below)
## NOTE 2: this can be run if only ONE MAIN satellite was run (with ensembles of course)
## NOTE 3: the results (pickle, image, video) are stored in a subfolder of /raid3/Armada/Charles/python/CADRE/result/image

# Give the name of the input file for the reference satellite
main_input_file_name = '/home/cbv/PropSim/input/main_input/' + sys.argv[1]  # ex of sys.argv[1]: cadre_tumbling_0panel.txt'

# read input file
input_variables, order_input_variables = read_input_file(main_input_file_name)
dt = input_variables[2]; nb_steps = input_variables[3]; satellite_to_plot_path = input_variables[6][0]; satellite_to_plot = input_variables[7][0]; nb_ensembles_coe = input_variables[12]; nb_ensembles_attitude = input_variables[13]; nb_ensembles_cd = input_variables[11]; ensemble_to_plot_temp = input_variables[14]; 
n = nb_steps

# Name of the mission of the simulation (choices: 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', 'other')
name_mission = 'CADRE' 

# Nb of ensembles
nb_spacecraft = 1#int(a_input_file_input[6][0:6])
nb_ensembles_array = [nb_ensembles_coe, nb_ensembles_attitude, nb_ensembles_cd]
nb_ensembles = np.max(nb_ensembles_array)
for i in range(len(nb_ensembles_array)):
    if ( (nb_ensembles_array[i] > 0 ) &  (nb_ensembles_array[i] < nb_ensembles) ) :
        nb_ensembles = nb_ensembles_array[i]

# Load the pickles
## Pickle for the reference satellite
save_pickle_name_refce_sat = '/raid3/Armada/Charles/python/' + name_mission + '/pickle/' + satellite_to_plot.replace(".txt",".pickle")
with open(save_pickle_name_refce_sat) as f:
    pickle_list, pickle_list_name = pickle.load(f)
    x_eci_main_sc = pickle_list[ pickle_list_name.index('x_eci_main_sc') ]
    y_eci_main_sc = pickle_list[ pickle_list_name.index('y_eci_main_sc') ]
    z_eci_main_sc = pickle_list[ pickle_list_name.index('z_eci_main_sc') ]
    vx_eci_main_sc = pickle_list[ pickle_list_name.index('vx_eci_main_sc') ]
    vy_eci_main_sc = pickle_list[ pickle_list_name.index('vy_eci_main_sc') ]
    vz_eci_main_sc = pickle_list[ pickle_list_name.index('vz_eci_main_sc') ]

## Pickle for the ensemble satellites 
pickle_to_load = ['cadre_tumbling_0panel', 'cadre_tumbling_1panel', 'cadre_tumbling_2panels', 'cadre_tumbling_3panels', 'cadre_tumbling_4panels']
algebric_distance_ensemble_main_sc_lvlh_x_diff_panel_conf = np.zeros([len(pickle_to_load), len(vz_eci_main_sc), nb_ensembles])
for ipickle in range(len(pickle_to_load)):
    save_pickle_name = '/raid3/Armada/Charles/python/CADRE/pickle/' + pickle_to_load[ipickle] + '_cd_4.pickle'
    with open(save_pickle_name) as f:
        pickle_list, pickle_list_name = pickle.load(f)
        x_eci_ensemble = pickle_list[ pickle_list_name.index('x_eci_ensemble') ]
        y_eci_ensemble = pickle_list[ pickle_list_name.index('y_eci_ensemble') ]
        z_eci_ensemble = pickle_list[ pickle_list_name.index('z_eci_ensemble') ]
        for i in range(nb_steps):
            for j in range(nb_ensembles):
                position_main_sc_eci = np.array([x_eci_main_sc[i], y_eci_main_sc[i], z_eci_main_sc[i]])
                position_ensemble_sc_eci = np.array([x_eci_ensemble[i,j], y_eci_ensemble[i,j], z_eci_ensemble[i,j]])
                delta_vector_eci = position_main_sc_eci - position_ensemble_sc_eci
                velocity_main_sc_eci = np.array([vx_eci_main_sc[i], vy_eci_main_sc[i], vz_eci_main_sc[i]])
                delta_vector_lvlh_of_main_sc = eci_to_lvlh(position_main_sc_eci, velocity_main_sc_eci, delta_vector_eci)
                algebric_distance_ensemble_main_sc_lvlh_x_diff_panel_conf[ipickle, i,j]= delta_vector_lvlh_of_main_sc[0] * 1000 # * 1000 to convert the distance from km to m 


# Calculate the position given by CADRE tle and compare it to the reference sc
index_when_plot = algebric_distance_ensemble_main_sc_lvlh_x_diff_panel_conf.shape[1] - 1 # we look at the last position of the main sc
position_main_sc_eci = np.array([x_eci_main_sc[index_when_plot], y_eci_main_sc[index_when_plot], z_eci_main_sc[index_when_plot]])
position_tle_sc_eci = cadre_read_last_tle('/home/cbv/PropSim/input/main_input/cadre_last_tle.txt')
delta_vector_eci = position_main_sc_eci - position_tle_sc_eci
velocity_main_sc_eci = np.array([vx_eci_main_sc[index_when_plot], vy_eci_main_sc[index_when_plot], vz_eci_main_sc[index_when_plot]])
delta_vector_lvlh_of_main_sc = eci_to_lvlh(position_main_sc_eci, velocity_main_sc_eci, delta_vector_eci)
algebric_distance_tle_main_sc_lvlh_x= delta_vector_lvlh_of_main_sc[0] * 1000 # * 1000 to convert the distance from km to m 

raise Exception

########################################################################################################################################################################################################################################################################################################
# PLOT
color_dist = ['k','b','r','g','m'] 
label_dist = ['0','1','2','3','4']
fontsize_plot = 14
height_fig = 10
ratio_fig_size = 4./3
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left= 0.08, right=0.98, top = 0.93,bottom = 0.08)
bin_width = 5. # in km
ax1 = fig.add_subplot(gs[0, 0])
for ipanel in range(len(pickle_to_load)):
    median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_x_diff_panel_conf[ipanel,index_when_plot,:]/1000.)
    n, bins, patches = ax1.hist(algebric_distance_ensemble_main_sc_lvlh_x_diff_panel_conf[ipanel,index_when_plot,:]/1000., np.arange(min(algebric_distance_ensemble_main_sc_lvlh_x_diff_panel_conf[ipanel,index_when_plot,:]/1000.), max(algebric_distance_ensemble_main_sc_lvlh_x_diff_panel_conf[ipanel,index_when_plot,:]/1000.) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = color_dist[ipanel], label = label_dist[ipanel]) 
    ax1.plot([median_algebric_ditance, median_algebric_ditance], [min(n), max(n)], color = color_dist[ipanel], linewidth = 2, linestyle = 'dotted')

## ADD THE TLE POSITION
ax1.scatter([algebric_distance_tle_main_sc_lvlh_x/1000.],[1.5], s = 200, c = 'r', marker = '*', edgecolor = 'k', label = 'CADRE')


plt.legend(loc = 2, scatterpoints=1)
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()
ax1.set_ylim([0, 100])
#    ax1.set_xlim([min(bins), 0])
ax1.margins(0,1) # autoscale both axes(fist value is for the x axis, second value for the y axis)


fig_save_name = '/raid3/Armada/Charles/python/CADRE/result/image/cadre_distribution_along_track_diff_panel_conf.png'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    




