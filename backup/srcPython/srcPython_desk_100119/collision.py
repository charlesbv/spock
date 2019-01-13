import scipy as sp
import pandas as pd
import seaborn as sns
from scipy.stats import norm
sns.set_style('white')
sns.set_context('talk')
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
from sampler_modif import *

## NOTE 0: to run this script, you first need to run distance_ensemble_to_main_sc.py (with first_time = 1). This will create the data that plot_ensembles will then use to makes plots
## NOTE 1: to use this script, the only 3 parameters you have to set are:
## - if yes or no you want to save the plots (save_results = 1 if yes, 0 otherwise)
## - if yes or no you want to show the plots (show_plots = 1 if yes, 0 otherwise)
## - the name of the propagator main input file 
## - the path of the folder where you want to store the results (pickle, image, video): called 'path_folder_results'. In this folder, there must be the following subfolders: 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', and 'other'. In each of these subfolders, there must be the 2 subsubfolders: 'result', and 'pickle'. In the subsubfolder 'result', there must be the 2 subsubsubfolders: 'image', and 'video'.
## NOTE 2: this can be run if only ONE MAIN satellite was run (with ensembles of course)
## NOTE 3: the results (pickle, image, video) are stored in a subfolder of path_folder_results. The subfolder is either 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', or 'other'. This code figures out which folder the simulation corresponds to by reading the name of the output file chosen by the user in the main input file of the propagator (third line of section #SPACECRAFTS): it tries to find 'cygnss', 'cadre', 'aerie', 'scion', or 'qb50' in the name of the output file. If it does not find it, then the results here will be stored in path_folder_results/other/
## NOTE 4: to run this script, and any python script in the propagator, you need to be one subfolder deep from the main folder where the propagator runs are made. So if path_to_propagator/PropSim is the folder where the propagator runs are made, then the python scripts must be run for example in path_to_propagator/PropSim/subfolder_where_python_scipts_are_run


# !!!!!!!!!! SET THE PARAMETER BELOW:
## Save or not the plots
save_results = 0

## Show or not the plots
show_plots = 1

## path of the folder where you want to store the results (pickle, image, video)
path_folder_results = '/raid3/Armada/Charles/python/' #get_prop_dir(2) + 'output/python_propagator/'


same_spock_input_file = 1
######## SPACECRAFT 1

## main input file (argument in command line)
if len(sys.argv) > 2:
    main_input_sc1_file_name = get_prop_dir(1) + sys.argv[1] + '/input/main_input/' + sys.argv[2]  # ex of sys.argv[1]: cadre_tumbling_0panel.txt'
else:
    main_input_sc1_file_name = get_prop_dir(1) + 'run/input/main_input/' + sys.argv[1]  # ex of sys.argv[1]: cadre_tumbling_0panel.txt'    


# read input_sc1 file
input_sc1_variables, order_input_sc1_variables = read_input_file(main_input_sc1_file_name)
date_start_sc1 = input_sc1_variables[0]; date_stop_sc1 = input_sc1_variables[1]; dt_sc1 = input_sc1_variables[2]; nb_steps_sc1 = input_sc1_variables[3]; satellite_to_plot_path_sc1 = input_sc1_variables[6][0]; satellite_to_plot_sc1 = input_sc1_variables[7][0]; nb_ensembles_sc1_coe = input_sc1_variables[12]; nb_ensembles_sc1_attitude = input_sc1_variables[13]; nb_ensembles_sc1_cd = input_sc1_variables[11]; ensemble_sc1_to_plot_temp = input_sc1_variables[14]; 
nb_ensembles_sc1_density = input_sc1_variables[17]
n_sc1 = nb_steps_sc1

# set up interactive figures
if show_plots == 1:
    plt.ion()

# Name of the mission of the simulation (choices: 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', 'other')
if ( 'cygnss' in satellite_to_plot_sc1.lower() ):
    name_mission_sc1 = 'CYGNSS' 
elif ( 'cadre' in satellite_to_plot_sc1.lower() ):
    name_mission_sc1 = 'CADRE' 
elif ( 'aerie' in satellite_to_plot_sc1.lower() ):
    name_mission_sc1 = 'AERIE' 
elif ( 'scion' in satellite_to_plot_sc1.lower() ):
    name_mission_sc1 = 'SCION' 
elif ( 'qb50' in satellite_to_plot_sc1.lower() ):
    name_mission_sc1 = 'QB50' 
else:
    name_mission_sc1 = 'other' 


    
name_subfolder_save = satellite_to_plot_sc1[:-5] + "/"
if ( save_results == 1 ):
    os.system("mkdir " + path_folder_results + name_mission_sc1 + '/result/image/' + name_subfolder_save )
    os.system('ssh -t srbwks2014-0008.engin.umich.edu "mkdir ' + name_mission_sc1 + '/' + name_subfolder_save + '"')
#    os.system("mkdir " + path_folder_results + name_mission_sc1 + '/result/video/' + name_subfolder_save )

    
save_pickle_name = path_folder_results + name_mission_sc1 + '/pickle/' + name_subfolder_save + satellite_to_plot_sc1.replace(".txt",".pickle")
root_save_fig_name = path_folder_results + name_mission_sc1 + '/result/image/' + name_subfolder_save + satellite_to_plot_sc1.replace(".txt","_") 
root_save_video_name = path_folder_results + name_mission_sc1 + '/result/video/' + name_subfolder_save + satellite_to_plot_sc1.replace(".txt","_")


# Ensembles_Sc1 created by the propagator
ensemble_sc1_to_plot = []
for i in range(len(ensemble_sc1_to_plot_temp)):
    if (ensemble_sc1_to_plot_temp[i] == 'eci_r'):
        ensemble_sc1_to_plot.append('x_eci'); ensemble_sc1_to_plot.append('y_eci'); ensemble_sc1_to_plot.append('z_eci')
    if (ensemble_sc1_to_plot_temp[i] == 'eci_v'):
        ensemble_sc1_to_plot.append('vx_eci'); ensemble_sc1_to_plot.append('vy_eci'); ensemble_sc1_to_plot.append('vz_eci')
    if (ensemble_sc1_to_plot_temp[i] == 'geodetic'):
        ensemble_sc1_to_plot.append('longitude'); ensemble_sc1_to_plot.append('latitude'); ensemble_sc1_to_plot.append('altitude')
    if (ensemble_sc1_to_plot_temp[i] == 'power'):
        ensemble_sc1_to_plot.append('power')
    if (ensemble_sc1_to_plot_temp[i] == 'attitude'):
        ensemble_sc1_to_plot.append('pitch'); ensemble_sc1_to_plot.append('roll'); ensemble_sc1_to_plot.append('yaw')
    if (ensemble_sc1_to_plot_temp[i] == 'oe'):
        ensemble_sc1_to_plot.append('sma'); ensemble_sc1_to_plot.append('inclination'); ensemble_sc1_to_plot.append('eccentricity'); ensemble_sc1_to_plot.append('true_anomaly'); ensemble_sc1_to_plot.append('RAAN'); ensemble_sc1_to_plot.append('argument_perigee');
    if (ensemble_sc1_to_plot_temp[i] == 'density'):
        ensemble_sc1_to_plot.append('rho'); ensemble_sc1_to_plot.append('f107'); ensemble_sc1_to_plot.append('f107a'); ensemble_sc1_to_plot.append('ap'); 

pickle_list = []
pickle_list_name = []
## Nb of ensembles_sc1
nb_spacecraft = 1#int(a_input_sc1_file_input_sc1[6][0:6])
nb_ensembles_sc1_array = [nb_ensembles_sc1_coe, nb_ensembles_sc1_attitude, nb_ensembles_sc1_cd, nb_ensembles_sc1_density]
nb_ensembles_sc1 = np.max(nb_ensembles_sc1_array)
for i in range(len(nb_ensembles_sc1_array)):
    if ( (nb_ensembles_sc1_array[i] > 0 ) &  (nb_ensembles_sc1_array[i] < nb_ensembles_sc1) ) :
        nb_ensembles_sc1 = nb_ensembles_sc1_array[i]

# # #######################################
dir_final_output_ensemble_sc1 = satellite_to_plot_path_sc1 + 'ensemble/'
## RESTORE THE RESULTS FROM THE PICKLE
with open(save_pickle_name) as f:
    pickle_list, pickle_list_name = pickle.load(f)
if ( 'x_eci_ensemble' in pickle_list_name ): 
    x_eci_ensemble_sc1 = pickle_list[ pickle_list_name.index('x_eci_ensemble') ]
if ( 'y_eci_ensemble' in pickle_list_name ): 
    y_eci_ensemble_sc1 = pickle_list[ pickle_list_name.index('y_eci_ensemble') ]
if ( 'z_eci_ensemble' in pickle_list_name ): 
    z_eci_ensemble_sc1 = pickle_list[ pickle_list_name.index('z_eci_ensemble') ]
if ( 'pitch_ensemble' in pickle_list_name ): 
    pitch_ensemble_sc1 = pickle_list[ pickle_list_name.index('pitch_ensemble') ]
if ( 'roll_ensemble' in pickle_list_name ): 
    roll_ensemble_sc1 = pickle_list[ pickle_list_name.index('roll_ensemble') ]
if ( 'yaw_ensemble' in pickle_list_name ): 
    yaw_ensemble_sc1 = pickle_list[ pickle_list_name.index('yaw_ensemble') ]
if ( 'x_eci_main_sc' in pickle_list_name ):
    x_eci_main_sc = pickle_list[ pickle_list_name.index('x_eci_main_sc') ]
if ( 'y_eci_main_sc' in pickle_list_name ):
    y_eci_main_sc = pickle_list[ pickle_list_name.index('y_eci_main_sc') ]
if ( 'z_eci_main_sc' in pickle_list_name ):
    z_eci_main_sc = pickle_list[ pickle_list_name.index('z_eci_main_sc') ]
if ( 'vx_eci_main_sc' in pickle_list_name ):
    vx_eci_main_sc = pickle_list[ pickle_list_name.index('vx_eci_main_sc') ]
if ( 'vy_eci_main_sc' in pickle_list_name ):
    vy_eci_main_sc = pickle_list[ pickle_list_name.index('vy_eci_main_sc') ]
if ( 'vz_eci_main_sc' in pickle_list_name ):
    vz_eci_main_sc = pickle_list[ pickle_list_name.index('vz_eci_main_sc') ]
if ( 'algebric_distance_ensemble_main_sc_lvlh_x' in pickle_list_name ):
    algebric_distance_ensemble_sc1_main_sc_lvlh_x = pickle_list[ pickle_list_name.index('algebric_distance_ensemble_main_sc_lvlh_x') ]
if ( 'algebric_distance_ensemble_main_sc_lvlh_y' in pickle_list_name ):
    algebric_distance_ensemble_sc1_main_sc_lvlh_y = pickle_list[ pickle_list_name.index('algebric_distance_ensemble_main_sc_lvlh_y') ]
if ( 'algebric_distance_ensemble_main_sc_lvlh_z' in pickle_list_name ):
    algebric_distance_ensemble_sc1_main_sc_lvlh_z = pickle_list[ pickle_list_name.index('algebric_distance_ensemble_main_sc_lvlh_z') ]
if ( 'sma_ensemble' in pickle_list_name ):
    sma_ensemble_sc1 = pickle_list[ pickle_list_name.index('sma_ensemble') ]
if ( 'inclination_ensemble' in pickle_list_name ):
    inclination_ensemble_sc1 = pickle_list[ pickle_list_name.index('inclination_ensemble') ]
if ( 'eccentricity_ensemble' in pickle_list_name ):
    eccentricity_ensemble_sc1 = pickle_list[ pickle_list_name.index('eccentricity_ensemble') ]
if ( 'true_anomaly_ensemble' in pickle_list_name ):
    true_anomaly_ensemble_sc1 = pickle_list[ pickle_list_name.index('true_anomaly_ensemble') ]
if ( 'RAAN_ensemble' in pickle_list_name ):
    RAAN_ensemble_sc1 = pickle_list[ pickle_list_name.index('RAAN_ensemble') ]
if ( 'argument_perigee_ensemble' in pickle_list_name ):
    argument_perigee_ensemble_sc1 = pickle_list[ pickle_list_name.index('argument_perigee_ensemble') ]
if ( 'rho_ensemble' in pickle_list_name ):
    rho_ensemble_sc1 = pickle_list[ pickle_list_name.index('rho_ensemble') ]    
if ( 'rho_ensemble_orbit_average' in pickle_list_name ):
    rho_ensemble_sc1_orbit_average = pickle_list[ pickle_list_name.index('rho_ensemble_orbit_average') ]    
if ( 'time_orbit_average' in pickle_list_name ):
    time_orbit_average = pickle_list[ pickle_list_name.index('time_orbit_average') ]    
if ( 'index_time_orbit_average' in pickle_list_name ):
    index_time_orbit_average = pickle_list[ pickle_list_name.index('index_time_orbit_average') ]    
if ( 'f107_ensemble' in pickle_list_name ):
    f107_ensemble_sc1 = pickle_list[ pickle_list_name.index('f107_ensemble') ]
if ( 'f107a_ensemble' in pickle_list_name ):
    f107a_ensemble_sc1 = pickle_list[ pickle_list_name.index('f107a_ensemble') ]
if ( 'ap_ensemble' in pickle_list_name ):
    ap_ensemble_sc1 = pickle_list[ pickle_list_name.index('ap_ensemble') ]
if ( 'angle_asc_node_to_sat_ensemble' in pickle_list_name ):
    angle_asc_node_to_sat_ensemble_sc1 = pickle_list[ pickle_list_name.index('angle_asc_node_to_sat_ensemble') ]
if ( 'angular_spacing_between_ensemble_sat' in pickle_list_name ):
    angular_spacing_between_ensemble_sc1_sat = pickle_list[ pickle_list_name.index('angular_spacing_between_ensemble_sat') ]
if ( 'angular_spacing_between_ensemble_sat_converted_in_a_distance' in pickle_list_name ):
    angular_spacing_between_ensemble_sc1_sat_converted_in_a_distance = pickle_list[ pickle_list_name.index('angular_spacing_between_ensemble_sat_converted_in_a_distance') ]





######## SPACECRAFT 2


# read input_sc2 file
if same_spock_input_file == 1:
    date_start_sc2 = input_sc1_variables[0]; date_stop_sc2 = input_sc1_variables[1]; dt_sc2 = input_sc1_variables[2]; nb_steps_sc2 = input_sc1_variables[3]; satellite_to_plot_path_sc2 = input_sc1_variables[6][1]; satellite_to_plot_sc2 = input_sc1_variables[7][1]; nb_ensembles_sc2_coe = input_sc1_variables[12]; nb_ensembles_sc2_attitude = input_sc1_variables[13]; nb_ensembles_sc2_cd = input_sc1_variables[11]; ensemble_sc2_to_plot_temp = input_sc1_variables[14]; 
    nb_ensembles_sc2_density = input_sc1_variables[17]
    n_sc2 = nb_steps_sc2

else:
    main_input_sc2_file_name = get_prop_dir(1) + sys.argv[1] + '/input/main_input/' + sys.argv[3]  # ex of sys.argv[1]: cadre_tumbling_0panel.txt'
    input_sc2_variables, order_input_sc2_variables = read_input_file(main_input_sc2_file_name)
    date_start_sc2 = input_sc2_variables[0]; date_stop_sc2 = input_sc2_variables[1]; dt_sc2 = input_sc2_variables[2]; nb_steps_sc2 = input_sc2_variables[3]; satellite_to_plot_path_sc2 = input_sc2_variables[6][0]; satellite_to_plot_sc2 = input_sc2_variables[7][0]; nb_ensembles_sc2_coe = input_sc2_variables[12]; nb_ensembles_sc2_attitude = input_sc2_variables[13]; nb_ensembles_sc2_cd = input_sc2_variables[11]; ensemble_sc2_to_plot_temp = input_sc2_variables[14]; 
    nb_ensembles_sc2_density = input_sc2_variables[17]
    n_sc2 = nb_steps_sc2

# set up interactive figures
if show_plots == 1:
    plt.ion()

# Name of the mission of the simulation (choices: 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', 'other')
if ( 'cygnss' in satellite_to_plot_sc2.lower() ):
    name_mission_sc2 = 'CYGNSS' 
elif ( 'cadre' in satellite_to_plot_sc2.lower() ):
    name_mission_sc2 = 'CADRE' 
elif ( 'aerie' in satellite_to_plot_sc2.lower() ):
    name_mission_sc2 = 'AERIE' 
elif ( 'scion' in satellite_to_plot_sc2.lower() ):
    name_mission_sc2 = 'SCION' 
elif ( 'qb50' in satellite_to_plot_sc2.lower() ):
    name_mission_sc2 = 'QB50' 
else:
    name_mission_sc2 = 'other' 


    
name_subfolder_save = satellite_to_plot_sc2[:-5] + "/"
if ( save_results == 1 ):
    os.system("mkdir " + path_folder_results + name_mission_sc2 + '/result/image/' + name_subfolder_save )
    os.system('ssh -t srbwks2014-0008.engin.umich.edu "mkdir ' + name_mission_sc2 + '/' + name_subfolder_save + '"')
#    os.system("mkdir " + path_folder_results + name_mission_sc2 + '/result/video/' + name_subfolder_save )

    
save_pickle_name = path_folder_results + name_mission_sc2 + '/pickle/' + name_subfolder_save + satellite_to_plot_sc2.replace(".txt",".pickle")
root_save_fig_name = path_folder_results + name_mission_sc2 + '/result/image/' + name_subfolder_save + satellite_to_plot_sc2.replace(".txt","_") 
root_save_video_name = path_folder_results + name_mission_sc2 + '/result/video/' + name_subfolder_save + satellite_to_plot_sc2.replace(".txt","_")


# Ensembles_Sc2 created by the propagator
ensemble_sc2_to_plot = []
for i in range(len(ensemble_sc2_to_plot_temp)):
    if (ensemble_sc2_to_plot_temp[i] == 'eci_r'):
        ensemble_sc2_to_plot.append('x_eci'); ensemble_sc2_to_plot.append('y_eci'); ensemble_sc2_to_plot.append('z_eci')
    if (ensemble_sc2_to_plot_temp[i] == 'eci_v'):
        ensemble_sc2_to_plot.append('vx_eci'); ensemble_sc2_to_plot.append('vy_eci'); ensemble_sc2_to_plot.append('vz_eci')
    if (ensemble_sc2_to_plot_temp[i] == 'geodetic'):
        ensemble_sc2_to_plot.append('longitude'); ensemble_sc2_to_plot.append('latitude'); ensemble_sc2_to_plot.append('altitude')
    if (ensemble_sc2_to_plot_temp[i] == 'power'):
        ensemble_sc2_to_plot.append('power')
    if (ensemble_sc2_to_plot_temp[i] == 'attitude'):
        ensemble_sc2_to_plot.append('pitch'); ensemble_sc2_to_plot.append('roll'); ensemble_sc2_to_plot.append('yaw')
    if (ensemble_sc2_to_plot_temp[i] == 'oe'):
        ensemble_sc2_to_plot.append('sma'); ensemble_sc2_to_plot.append('inclination'); ensemble_sc2_to_plot.append('eccentricity'); ensemble_sc2_to_plot.append('true_anomaly'); ensemble_sc2_to_plot.append('RAAN'); ensemble_sc2_to_plot.append('argument_perigee');
    if (ensemble_sc2_to_plot_temp[i] == 'density'):
        ensemble_sc2_to_plot.append('rho'); ensemble_sc2_to_plot.append('f107'); ensemble_sc2_to_plot.append('f107a'); ensemble_sc2_to_plot.append('ap'); 

pickle_list = []
pickle_list_name = []
## Nb of ensembles_sc2
nb_spacecraft = 1#int(a_input_sc2_file_input_sc2[6][0:6])
nb_ensembles_sc2_array = [nb_ensembles_sc2_coe, nb_ensembles_sc2_attitude, nb_ensembles_sc2_cd, nb_ensembles_sc2_density]
nb_ensembles_sc2 = np.max(nb_ensembles_sc2_array)
for i in range(len(nb_ensembles_sc2_array)):
    if ( (nb_ensembles_sc2_array[i] > 0 ) &  (nb_ensembles_sc2_array[i] < nb_ensembles_sc2) ) :
        nb_ensembles_sc2 = nb_ensembles_sc2_array[i]

# # #######################################
dir_final_output_ensemble_sc2 = satellite_to_plot_path_sc2 + 'ensemble/'
## RESTORE THE RESULTS FROM THE PICKLE
with open(save_pickle_name) as f:
    pickle_list, pickle_list_name = pickle.load(f)
if ( 'x_eci_ensemble' in pickle_list_name ): 
    x_eci_ensemble_sc2 = pickle_list[ pickle_list_name.index('x_eci_ensemble') ]
if ( 'y_eci_ensemble' in pickle_list_name ): 
    y_eci_ensemble_sc2 = pickle_list[ pickle_list_name.index('y_eci_ensemble') ]
if ( 'z_eci_ensemble' in pickle_list_name ): 
    z_eci_ensemble_sc2 = pickle_list[ pickle_list_name.index('z_eci_ensemble') ]
if ( 'pitch_ensemble' in pickle_list_name ): 
    pitch_ensemble_sc2 = pickle_list[ pickle_list_name.index('pitch_ensemble') ]
if ( 'roll_ensemble' in pickle_list_name ): 
    roll_ensemble_sc2 = pickle_list[ pickle_list_name.index('roll_ensemble') ]
if ( 'yaw_ensemble' in pickle_list_name ): 
    yaw_ensemble_sc2 = pickle_list[ pickle_list_name.index('yaw_ensemble') ]
if ( 'x_eci_main_sc' in pickle_list_name ):
    x_eci_main_sc = pickle_list[ pickle_list_name.index('x_eci_main_sc') ]
if ( 'y_eci_main_sc' in pickle_list_name ):
    y_eci_main_sc = pickle_list[ pickle_list_name.index('y_eci_main_sc') ]
if ( 'z_eci_main_sc' in pickle_list_name ):
    z_eci_main_sc = pickle_list[ pickle_list_name.index('z_eci_main_sc') ]
if ( 'vx_eci_main_sc' in pickle_list_name ):
    vx_eci_main_sc = pickle_list[ pickle_list_name.index('vx_eci_main_sc') ]
if ( 'vy_eci_main_sc' in pickle_list_name ):
    vy_eci_main_sc = pickle_list[ pickle_list_name.index('vy_eci_main_sc') ]
if ( 'vz_eci_main_sc' in pickle_list_name ):
    vz_eci_main_sc = pickle_list[ pickle_list_name.index('vz_eci_main_sc') ]
if ( 'algebric_distance_ensemble_main_sc_lvlh_x' in pickle_list_name ):
    algebric_distance_ensemble_sc2_main_sc_lvlh_x = pickle_list[ pickle_list_name.index('algebric_distance_ensemble_main_sc_lvlh_x') ]
if ( 'algebric_distance_ensemble_main_sc_lvlh_y' in pickle_list_name ):
    algebric_distance_ensemble_sc2_main_sc_lvlh_y = pickle_list[ pickle_list_name.index('algebric_distance_ensemble_main_sc_lvlh_y') ]
if ( 'algebric_distance_ensemble_main_sc_lvlh_z' in pickle_list_name ):
    algebric_distance_ensemble_sc2_main_sc_lvlh_z = pickle_list[ pickle_list_name.index('algebric_distance_ensemble_main_sc_lvlh_z') ]
if ( 'sma_ensemble' in pickle_list_name ):
    sma_ensemble_sc2 = pickle_list[ pickle_list_name.index('sma_ensemble') ]
if ( 'inclination_ensemble' in pickle_list_name ):
    inclination_ensemble_sc2 = pickle_list[ pickle_list_name.index('inclination_ensemble') ]
if ( 'eccentricity_ensemble' in pickle_list_name ):
    eccentricity_ensemble_sc2 = pickle_list[ pickle_list_name.index('eccentricity_ensemble') ]
if ( 'true_anomaly_ensemble' in pickle_list_name ):
    true_anomaly_ensemble_sc2 = pickle_list[ pickle_list_name.index('true_anomaly_ensemble') ]
if ( 'RAAN_ensemble' in pickle_list_name ):
    RAAN_ensemble_sc2 = pickle_list[ pickle_list_name.index('RAAN_ensemble') ]
if ( 'argument_perigee_ensemble' in pickle_list_name ):
    argument_perigee_ensemble_sc2 = pickle_list[ pickle_list_name.index('argument_perigee_ensemble') ]
if ( 'rho_ensemble' in pickle_list_name ):
    rho_ensemble_sc2 = pickle_list[ pickle_list_name.index('rho_ensemble') ]    
if ( 'rho_ensemble_orbit_average' in pickle_list_name ):
    rho_ensemble_sc2_orbit_average = pickle_list[ pickle_list_name.index('rho_ensemble_orbit_average') ]    
if ( 'time_orbit_average' in pickle_list_name ):
    time_orbit_average = pickle_list[ pickle_list_name.index('time_orbit_average') ]    
if ( 'index_time_orbit_average' in pickle_list_name ):
    index_time_orbit_average = pickle_list[ pickle_list_name.index('index_time_orbit_average') ]    
if ( 'f107_ensemble' in pickle_list_name ):
    f107_ensemble_sc2 = pickle_list[ pickle_list_name.index('f107_ensemble') ]
if ( 'f107a_ensemble' in pickle_list_name ):
    f107a_ensemble_sc2 = pickle_list[ pickle_list_name.index('f107a_ensemble') ]
if ( 'ap_ensemble' in pickle_list_name ):
    ap_ensemble_sc2 = pickle_list[ pickle_list_name.index('ap_ensemble') ]
if ( 'angle_asc_node_to_sat_ensemble' in pickle_list_name ):
    angle_asc_node_to_sat_ensemble_sc2 = pickle_list[ pickle_list_name.index('angle_asc_node_to_sat_ensemble') ]
if ( 'angular_spacing_between_ensemble_sat' in pickle_list_name ):
    angular_spacing_between_ensemble_sc2_sat = pickle_list[ pickle_list_name.index('angular_spacing_between_ensemble_sat') ]
if ( 'angular_spacing_between_ensemble_sat_converted_in_a_distance' in pickle_list_name ):
    angular_spacing_between_ensemble_sc2_sat_converted_in_a_distance = pickle_list[ pickle_list_name.index('angular_spacing_between_ensemble_sat_converted_in_a_distance') ]

######### DISTANCE BETWEEN EACH ENSEMBLE OF SC1 AND EACH ENSEMBLE OF SC2
###### ASSUMPTION: SAME TIME STEP AND SAME DATE_START IN BOTH RUNS
if n_sc1 != n_sc2:
    print "***! The number of output time steps of the two propagations must be equal. The program will stop. ***!"; raise Exception
if dt_sc1 != dt_sc2:
    print "***! The output time step of the two propagations must be equal. The program will stop. ***!"; raise Exception
nb_steps_sc1_sc2 = x_eci_ensemble_sc1.shape[0]#(int)( 2000. / dt_sc1) #min([n_sc1, n_sc2])  # just for the first time step (temporary) min([n_sc1, n_sc2])
d_sc1_sc2_list = []
# NUMBER OF SETS SC1/SC2 WITH DISTANCE d_sc1_sc2 LOWER THAN THRESHOLD
threshold_d_sc1_sc2 = 15./ 1000 # in km
nb_collision = 0
# nb_collision_arr = np.zeros([nb_steps_sc1_sc2])
# cdf_collision = np.zeros([nb_steps_sc1_sc2]) # Cumulative proability distribution
# for istep in range(320, 400):#nb_steps_sc1_sc2):
#     d_sc1_sc2_list_sub = []
#     for isc1 in range(nb_ensembles_sc1):
# #        print isc1, nb_ensembles_sc1
#         for isc2 in range(nb_ensembles_sc2):
#             d_sc1_sc2_temp = np.sqrt( ( x_eci_ensemble_sc1[istep, isc1] - x_eci_ensemble_sc2[istep, isc2] )**2 + ( y_eci_ensemble_sc1[istep, isc1] - y_eci_ensemble_sc2[istep, isc2] )**2 + ( z_eci_ensemble_sc1[istep, isc1] - z_eci_ensemble_sc2[istep, isc2] )**2 )
#             d_sc1_sc2_list_sub.append( [ isc1, isc2, d_sc1_sc2_temp ] )
#             if d_sc1_sc2_temp < threshold_d_sc1_sc2:
#                 nb_collision = nb_collision + 1
#                 nb_collision_arr[istep] = nb_collision_arr[istep] + 1
#     d_sc1_sc2_list.append(d_sc1_sc2_list_sub)

#     if istep == 0:
#         cdf_collision[istep] = nb_collision_arr[istep] / ( nb_ensembles_sc1 * nb_ensembles_sc2) * dt_sc1
#     else:
#         cdf_collision[istep] = cdf_collision[istep-1] + nb_collision_arr[istep] / ( nb_ensembles_sc1 * nb_ensembles_sc2) * dt_sc1
#     print istep, nb_steps_sc1_sc2, cdf_collision[istep]
#d_sc1_sc2 = np.array(d_sc1_sc2_list)
step_start = 320
step_stop = 400
nb_collision_arr = np.zeros([step_stop - step_start])
nb_collision = []
d_sc1_sc2 = np.zeros([nb_ensembles_sc1, nb_ensembles_sc2, step_stop - step_start])
for isc1 in range(nb_ensembles_sc1):
    print isc1, nb_ensembles_sc1 - 1
    for isc2 in range(nb_ensembles_sc2):
        istep = step_start
        distance_temp = 1e6
        while ( ( istep < step_stop ) & ( distance_temp > threshold_d_sc1_sc2 ) ):
            distance_temp = np.sqrt( ( x_eci_ensemble_sc1[istep, isc1] - x_eci_ensemble_sc2[istep, isc2] )**2 + ( y_eci_ensemble_sc1[istep, isc1] - y_eci_ensemble_sc2[istep, isc2] )**2 + ( z_eci_ensemble_sc1[istep, isc1] - z_eci_ensemble_sc2[istep, isc2] )**2 )
            d_sc1_sc2[isc1, isc2,istep - step_start] = distance_temp
            if ( distance_temp < threshold_d_sc1_sc2 ):
                nb_collision.append( [istep, isc1, isc2] )
            istep = istep + 1


print "<<<<<<<<< TOTAL COLLISIONS BUT ONLY AT STEPS: " + str(len(nb_collision)) + " >>>>>>>>>"
with open('d_sc1_sc2.pickle', 'w') as f:
    pickle.dump([d_sc1_sc2], f)
with open('nb_collision_arr.pickle', 'w') as f:
    pickle.dump([nb_collision], f)

# with open('nb_collision_arr.pickle') as f:
#     [nb_collision_arr] = pickle.load(f)
            

# Probability of collision
# collision_probability = np.double( nb_collision ) / ( nb_ensembles_sc1 * nb_ensembles_sc2)
# print "The probability of collision computed by SpOCK between the two spacecraft is: " + '{:.2e}'.format(collision_probability) + " (" + str(nb_collision) + " collisions out of " + '{0:.0e}'.format(nb_ensembles_sc1 * nb_ensembles_sc2) +" spacecraft)."


raise Exception
################## SUBSET SIMULATION
######### SAME NOTATIONS AS IN MORCELLI15
p0 = 0.2
D = threshold_d_sc1_sc2
N = nb_ensembles_sc1 * nb_ensembles_sc2
g = np.zeros([nb_steps_sc1_sc2, N, 3])
# !!!!!!!!!! FROM NOW ON ISTEP = 0
dist_only = d_sc1_sc2[0, :, 2]
m = 1
while (True):
    g = D - dist_only
    indices_g_storted = [i[0] for i in sorted(enumerate(g), key=lambda x:x[1])]
    g_sorted = g[ indices_g_storted ]

    one_minus_p0_times_N = (int)( ( 1 - p0 ) * N )
    D_l_plus_one = D - g_sorted[one_minus_p0_times_N]
    print D_l_plus_one, D
    if D_l_plus_one <= D:
        nb_collisions_subset_simu = len(np.where( g >= 0 )[0])
        break;
    else: # USE MCMC TO GENERATE ADDITIONAL SAMPLES
        data = dist_only[indices_g_storted[one_minus_p0_times_N:]]
        np.random.seed(123)
        dist_only = sampler_modif(data, total_sample_target = N, mu_init=D_l_plus_one, proposal_width = ( g_sorted[-1] - g_sorted[one_minus_p0_times_N] ) , accept_criteria = D_l_plus_one);
        dist_only = np.array(dist_only)
        m = m + 1

N_T = N
collision_probability_mcmc = p0**(m-1) * np.double( nb_collisions_subset_simu ) / ( N )
#print "The probability of collision computed by SpOCK between the two spacecraft, using the subset simulation method, is: " + '{:.2e}'.format(collision_probability_mcmc) + "."
file_result = open("result+mcmc.txt", "w")
print >> file_result, m, nb_collisions_subset_simu, collision_probability_mcmc
file_result.close()
raise Exception
# PLOT
plt.ion()
ax = plt.subplot()
#sns.distplot(data, kde=False, ax=ax, color = 'red', label = 'Data')
sns.distplot(result, ax=ax, color = 'b', label = 'MCMC')


fig, ax = plt.subplots()
ax.hist(d_sc1_sc2[0, :, 2], range = (0, D_l_plus_one) )
plt.show(); plt.show()


fig, ax = plt.subplots()
ax.hist(data)
ax.hist(result, color = 'r')
plt.show(); plt.show()


# g[:, :, 2] = D - d_sc1_sc2[:, :, 2]
# g[:, :, 0] = d_sc1_sc2[:, :, 0]
# g[:, :, 1] = d_sc1_sc2[:, :, 1]
# g_sorted =  np.zeros([nb_steps_sc1_sc2, N, 3])
# for istep in range(nb_steps_sc1_sc2):    
#     indices_g_storted = [i[0] for i in sorted(enumerate(g[istep, :, 2]), key=lambda x:x[1])]
#     g_sorted[istep, :, 2] = g[istep, indices_g_storted, 2]
#     g_sorted[istep, :, 0] = g[istep, indices_g_storted, 0]
#     g_sorted[istep, :, 1] = g[istep, indices_g_storted, 1]


# one_minus_p0_times_N = (int)( ( 1 - p0 ) * N )
# D_l_plus_one = D - g_sorted[0, one_minus_p0_times_N, 2]

