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

# this script plot the position of all the ensembles of each main sc in a scatter 3d plot.
## NOTE 0: to run this script, you first need to run new_mpi_distance_ensemble_to_main_sc.py (with first_time = 1). This will create the data that plot_ensembles will then use to makes plots
## !!!!!!! to compute the r/v at TCA, i modified SpOCK on Big so that it outputs r/v of ensembles (not main sc) at time steps: 0, time step with tca in it - dt, time step with tca in it, time step with tca in it + dt. Then i linear interpolate to get the r/v at tca (the accuracy of lin interpolating is good enough for the plots we make here). Here is the modif I made on big in generate_ephemerides.c:
## line ~ 1445 i commented the if and replaced it with:
##            if ( ( CONSTELLATION->spacecraft[ii][1 + iProc * OPTIONS->nb_ensemble_min_per_proc].et >= et_current_tca - OPTIONS->dt ) && ( CONSTELLATION->spacecraft[ii][1 + iProc * OPTIONS->nb_ensemble_min_per_proc].et <= et_current_tca + OPTIONS->dt  ) ){ // !!!!!!!!! REMOVE THIS IF AND UNCOMMENT THE ONE RIGHT ABOVE                                                                                                                                                         #NOTE: the quier works weirdly on Pleaides but well on linux desktop        

import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/srcPython")
#sys.path.append("/home1/cbussy/Code/spock/srcPython")
from find_in_read_input_order_variables import *

from cadre_read_last_tle import *
from get_prop_dir import *
import matplotlib.gridspec as gridspec
from read_input_file import *
from read_output_file import *
from matplotlib.colors import LogNorm
import pickle
from eci_to_lvlh import *
import fileinput
import time
from datetime import datetime, timedelta
import numpy as np
from matplotlib import pyplot as plt
import os
import subprocess
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
from read_collision_file_new import *
from read_collision_file_new_with_dca import *


# !!!!!!!!!! SET THE PARAMETER BELOW:
## Save or not the plots
save_results = 1

## Show or not the plots
show_plots = 0

## if laready laoded the pickle, set to 1
already_loaded = 1

list_run = ["oct31_h06_f107_ap_actual.txt"]#["oct31_h06_f107_0_ap_0_with_ens.txt"]#["oct30_h07_f107_0_ap_0_bc_variance_vcm.txt"]
#['FM07_july23_hbr_2_5_f107_ap_static_less_ens.txt']
#['FM07_july22_f107_ap_actual_less_ensemble.txt']
#'eq6_FM07_nov01_f107_-13_ap_-13.txt','FM07_nov01_f107_-13_ap_-13.txt']
#['fm07_01nov_bc010_no_srp010_eq9.txt', 'fm07_01nov_bc010_no_srp010.txt']
#['coll_perturb_eq9_f107_-13_ap_-13.txt','coll_perturb_eq9_f107_13_ap_13.txt', 'coll_perturb_eq9.txt']#, 'coll_bc23_no_srp.txt', 'coll_bc05_no_srp.txt']

## main input file (argument in command line)
#list_run = ['alt400-400_inc90-0_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_f107100_quantile5_forward.txt' ,'alt400-400_inc90-0_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_f107100_quantile9_forward.txt', 'alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_f107100_quantile5_forward.txt', 'alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_f107100_quantile9_forward.txt']
# list_run = ['alt400-400_inc90-0_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_f107100_quantile1_forward.txt','alt400-400_inc90-0_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_f107100_quantile3_forward.txt','alt400-400_inc90-0_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_f107100_quantile5_forward.txt', 'alt400-400_inc90-0_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_f107100_quantile7_forward.txt','alt400-400_inc90-0_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_f107100_quantile9_forward.txt'] # need to be 5 runs, all same oe/bc/f107. just difference between run is quantile


#list_run = ['alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_f107100_quantile1_forward.txt','alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_f107100_quantile3_forward.txt','alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_f107100_quantile5_forward.txt', 'alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_f107100_quantile7_forward.txt','alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_f107100_quantile9_forward.txt'] # need to be 5 runs, all same oe/bc/f107. just difference between run is quantile


# list_run = ['alt400-400_inc90-80_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_f107100_quantile1_forward.txt','alt400-400_inc90-80_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_f107100_quantile3_forward.txt','alt400-400_inc90-80_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_f107100_quantile5_forward.txt', 'alt400-400_inc90-80_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_f107100_quantile7_forward.txt','alt400-400_inc90-80_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_f107100_quantile9_forward.txt'] # need to be 5 runs, all same oe/bc/f107. just difference between run is quantile



#list_run = ['alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano180-180_ecc0-0_f107100_ap12_quantile1_forward.txt','alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano180-180_ecc0-0_f107100_ap12_quantile3_forward.txt','alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano180-180_ecc0-0_f107100_ap12_quantile5_forward.txt', 'alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano180-180_ecc0-0_f107100_ap12_quantile7_forward.txt','alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano180-180_ecc0-0_f107100_ap12_quantile9_forward.txt'] # need to be 5 runs, all same oe/bc/f107. just difference between run is quantile

# list_run = ['alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano180-180_ecc0-0.0001_f107100_ap12_quantile1_forward.txt','alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano180-180_ecc0-0.0001_f107100_ap12_quantile3_forward.txt','alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano180-180_ecc0-0.0001_f107100_ap12_quantile5_forward.txt', 'alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano180-180_ecc0-0.0001_f107100_ap12_quantile7_forward.txt','alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano180-180_ecc0-0.0001_f107100_ap12_quantile9_forward.txt'] # need to be 5 runs, all same oe/bc/f107. just difference between run is quantile


# list_run = ['alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano180-180_ecc0-0.005_f107100_ap12_quantile1_forward.txt','alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano180-180_ecc0-0.005_f107100_ap12_quantile3_forward.txt','alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano180-180_ecc0-0.005_f107100_ap12_quantile5_forward.txt', 'alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano180-180_ecc0-0.005_f107100_ap12_quantile7_forward.txt','alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano180-180_ecc0-0.005_f107100_ap12_quantile9_forward.txt'] # need to be 5 runs, all same oe/bc/f107. just difference between run is quantile


irun = 0
main_input_file_name = list_run[irun] #sys.argv[1]
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3

color_arr = ['b','r']


fig_title = ''#main_input_file_name.replace(".txt",".pdf") #'Inc ' + main_input_file_name.split('_inc')[1].split('-')[1].split('_')[0] + u'\N{DEGREE SIGN} - ecc ' + main_input_file_name.split('_ecc')[1].split('-')[1].split('_')[0]# + ' - quantile ' + main_input_file_name.split('quantile')[1].split('.')[0].split('_')[0]
    #alt400-400_inc90-45_arg_per0-0_raan0-0_true_ano0-0_ecc0-0_f107100_quantile5_forward.txt
z_label = 'Z'
y_label = 'Y'
x_label = 'X'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(2, 2)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.1, wspace = 0.01)

nb_run = len(list_run)
color_by_nb_sc_per_bin = 0 # set to 1 if you want to color code by nmber of sc per bin
color_arr = ['blue', 'red'] # this is used only if you don't want to color cod ewith respect to the number of sc per bin (color_by_nb_sc_per_bin = 1)
for irun in range(nb_run):
    main_input_file_name = list_run[irun] #sys.argv[1]
    first_time = 1
    # read input file
    input_variables, order_input_variables = read_input_file(main_input_file_name)
    date_start = input_variables[find_in_read_input_order_variables(order_input_variables, 'date_start')];
    dt = input_variables[find_in_read_input_order_variables(order_input_variables, 'dt')]; # need to be time step of itnerpolation because when computing collision if output_only_at_tca is 1 in generate_ephemerides.c then the output file for the eensmebles will be written every dt interpolation (not every dt otuput)
    n = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_steps')];
    satellite_to_plot_path = input_variables[find_in_read_input_order_variables(order_input_variables, 'output_file_path_list')];
    satellite_to_plot = input_variables[find_in_read_input_order_variables(order_input_variables, 'output_file_name_list')];
    nb_ensembles_coe = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_ensembles_coe')];
    nb_ensembles_attitude = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_ensembles_attitude')];
    nb_ensembles_cd = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_ensembles_cd')];
    ensemble_to_plot_temp = input_variables[find_in_read_input_order_variables(order_input_variables, 'ensembles_to_output')];
    nb_ensembles_density = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_ensembles_density')]
    nb_spacecraft = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_sc')]
    compute_drag = input_variables[find_in_read_input_order_variables(order_input_variables, 'compute_drag')]
    collision_filename = '/'.join(satellite_to_plot_path[0].split('/')[:-2]) + "/" + satellite_to_plot_path[0].split('/')[-2][:-1] + "_collision.txt"
    #date, nb_collisions_each_dt, cpc, cpc_final, tca, tca_no_coll  = read_collision_file_new(collision_filename) 
    date, nb_collisions_each_dt, cpc, cpc_final, tca, tca_no_coll, dca, dca_no_coll  = read_collision_file_new_with_dca(collision_filename) # ic reated this script after i modfied SpOCK to include DCA inthe collison output file. owver, i had made simulations before this modid. in those simu, dca is not included so you need to use read_collision_file_new, not read_collision_file_new_with_dca
    if len(tca) > 0:
        tca_date = datetime.strptime(tca[0], "%Y-%m-%dT%H:%M:%S.%f") # !!!!!!! only one close approach
    elif len(tca_no_coll) > 0:
        tca_date = datetime.strptime(tca_no_coll[0], "%Y-%m-%dT%H:%M:%S.%f") # !!!!!!! only one close approach
    else:
        print "***! There was no close approach. The program will stop. ***!"; raise Exception
        

    out_var, order_var = read_output_file( satellite_to_plot_path[0] + satellite_to_plot[0], []  )
    date_main_datetime_round_sec = out_var[find_in_read_input_order_variables(order_var, 'date_datetime_round_sec')]
    date_main_datetime_round_sec = np.array(date_main_datetime_round_sec)

    path_folder_results = '.'#'/home1/cbussy/Code/collision/postPro'
#    path_folder_results = '/Users/cbv/coll' # !!!!!! need t be the same as new_mpi_distance_ensemble_to_main_sc.py

    if (os.path.isdir(path_folder_results) == False):
        os.system("mkdir " + path_folder_results)

    if path_folder_results[-1] != '!':
        path_folder_results = path_folder_results + '/'



    # set up interactive figures
    if show_plots == 1:
        plt.ion()



    # Ensembles created by the propagator
    ensemble_to_plot = []
    for i in range(len(ensemble_to_plot_temp)):
        if (ensemble_to_plot_temp[i] == 'eci_r'):
            ensemble_to_plot.append('x_eci'); ensemble_to_plot.append('y_eci'); ensemble_to_plot.append('z_eci')
        if (ensemble_to_plot_temp[i] == 'eci_v'):
            ensemble_to_plot.append('vx_eci'); ensemble_to_plot.append('vy_eci'); ensemble_to_plot.append('vz_eci')
        if (ensemble_to_plot_temp[i] == 'geodetic'):
            ensemble_to_plot.append('longitude'); ensemble_to_plot.append('latitude'); ensemble_to_plot.append('altitude')
        if (ensemble_to_plot_temp[i] == 'cd'):
            ensemble_to_plot.append('cd')
        if (ensemble_to_plot_temp[i] == 'power'):
            ensemble_to_plot.append('power')
        if (ensemble_to_plot_temp[i] == 'attitude'):
            ensemble_to_plot.append('pitch'); ensemble_to_plot.append('roll'); ensemble_to_plot.append('yaw')
        if (ensemble_to_plot_temp[i] == 'oe'):
            ensemble_to_plot.append('sma'); ensemble_to_plot.append('inclination'); ensemble_to_plot.append('eccentricity'); ensemble_to_plot.append('true_anomaly'); ensemble_to_plot.append('RAAN'); ensemble_to_plot.append('argument_perigee');ensemble_to_plot.append('phase_angle');ensemble_to_plot.append('sma_difference');
        if ((ensemble_to_plot_temp[i] == 'density') & (compute_drag == 1)):
            ensemble_to_plot.append('rho'); ensemble_to_plot.append('f107'); ensemble_to_plot.append('f107a'); ensemble_to_plot.append('ap'); 
        if (ensemble_to_plot_temp[i] == 'collision'):
            ensemble_to_plot.append('tca'); ensemble_to_plot.append('dca');
        if ((ensemble_to_plot_temp[i] == 'cd') & (compute_drag == 1)):
            ensemble_to_plot.append('cd')


        if ((compute_drag == 0) & (ensemble_to_plot_temp[i] == 'density')):
            print "Note: in the main input file of SpOCK (" + main_input_file_name + "), you chose to output ensembles for the density. However, the atmospheric drag was not computed. Therefore, the density is not computed either"
    pickle_list = []
    pickle_list_name = []
    pickle_list_name_temp = []
    ## Nb of ensembles
    nb_ensembles_array = [nb_ensembles_coe, nb_ensembles_attitude, nb_ensembles_cd, nb_ensembles_density]
    nb_ensembles = np.max(nb_ensembles_array)
    for i in range(len(nb_ensembles_array)):
        if ( (nb_ensembles_array[i] > 0 ) &  (nb_ensembles_array[i] < nb_ensembles) ) :
            nb_ensembles = nb_ensembles_array[i]

    x_eci_ensemble= [] 
    y_eci_ensemble= [] 
    z_eci_ensemble= [] 
    x_eci_main_sc= [] 
    y_eci_main_sc= [] 
    z_eci_main_sc= [] 
    vx_eci_main_sc= [] 
    vy_eci_main_sc= [] 
    vz_eci_main_sc= [] 
    algebric_distance_ensemble_main_sc_lvlh_x= [] 
    algebric_distance_ensemble_main_sc_lvlh_y= [] 
    algebric_distance_ensemble_main_sc_lvlh_z= [] 
    tca_ensemble= [] 
    dca_ensemble= [] 
    nb_steps = []

    for isat in range(nb_spacecraft):
        name_subfolder_save = satellite_to_plot[isat][:-5] + "/"
        if ( save_results == 1 ):

            if os.path.isdir(path_folder_results + 'result') == False:
                os.system("mkdir " + path_folder_results + 'result' )
            if os.path.isdir(path_folder_results + 'result/' + name_subfolder_save) == False:
                os.system("mkdir " + path_folder_results  + 'result/' + name_subfolder_save )
                #    os.system('ssh -t srbwks2014-0008.engin.umich.edu "mkdir ' + name_mission + '/' + name_subfolder_save + '"')
                #    os.system("mkdir " + path_folder_results + name_mission + '/result/video/' + name_subfolder_save )


        save_pickle_name = path_folder_results + 'pickle/' + name_subfolder_save + satellite_to_plot[isat].replace(".txt",".pickle")
        root_save_fig_name = path_folder_results  + 'result/' + name_subfolder_save + satellite_to_plot[isat].replace(".txt","_") 

        # # #######################################
        dir_final_output_ensemble = satellite_to_plot_path[isat] + 'ensemble/'
        ## RESTORE THE RESULTS FROM THE PICKLE
        list_save_pickle = [f for f in os.listdir(path_folder_results +  'pickle/' + name_subfolder_save) if ( satellite_to_plot[isat].replace(".txt", "") in f ) | ( '_tca.pickle' in f) | ( '_dca.pickle' in f)]

        nb_pickle = len(list_save_pickle)
        for ipickle in range(nb_pickle):
            save_pickle_name = path_folder_results + 'pickle/' + name_subfolder_save + list_save_pickle[ipickle]
            with open(save_pickle_name) as f:
                print "Loading pickle " + save_pickle_name + "..."
                pickle_list, pickle_list_name_temp = pickle.load(f)
                for ielt in range(len(pickle_list_name_temp)):
                    if ((pickle_list_name_temp[ielt] in pickle_list_name) == False):
                        pickle_list_name.append(pickle_list_name_temp[ielt])
            if ipickle == 0:
                nProcs = pickle_list[ pickle_list_name_temp.index('nProcs') ]
                nb_ensembles= (int)(nb_ensembles / nProcs) * nProcs
                nb_ensembles_per_proc = (int)(nb_ensembles /  nProcs)
            if nProcs > np.sqrt( nb_ensembles ):
                print "Note: to plot the density, F10.7, F10.7A, and Ap of ensembles as a function of time, it is recommended to run SpOCK with nProcs <= sqrt( nb_ensembles ) (it will still run fine and plot fine)."
            if ( 'x_eci_ensemble' in pickle_list_name_temp ): 
                x_eci_ensemble.append(pickle_list[ pickle_list_name_temp.index('x_eci_ensemble') ])
            if ( 'y_eci_ensemble' in pickle_list_name_temp ): 
                y_eci_ensemble.append(pickle_list[ pickle_list_name_temp.index('y_eci_ensemble') ])
            if ( 'z_eci_ensemble' in pickle_list_name_temp ): 
                z_eci_ensemble.append(pickle_list[ pickle_list_name_temp.index('z_eci_ensemble') ])
            if ( 'x_eci_main_sc' in pickle_list_name_temp ):
                x_eci_main_sc.append(pickle_list[ pickle_list_name_temp.index('x_eci_main_sc') ])
            if ( 'y_eci_main_sc' in pickle_list_name_temp ):
                y_eci_main_sc.append(pickle_list[ pickle_list_name_temp.index('y_eci_main_sc') ])
            if ( 'z_eci_main_sc' in pickle_list_name_temp ):
                z_eci_main_sc.append(pickle_list[ pickle_list_name_temp.index('z_eci_main_sc') ])
            if ( 'vx_eci_main_sc' in pickle_list_name_temp ):
                vx_eci_main_sc.append(pickle_list[ pickle_list_name_temp.index('vx_eci_main_sc') ])
            if ( 'vy_eci_main_sc' in pickle_list_name_temp ):
                vy_eci_main_sc.append(pickle_list[ pickle_list_name_temp.index('vy_eci_main_sc') ])
            if ( 'vz_eci_main_sc' in pickle_list_name_temp ):
                vz_eci_main_sc.append(pickle_list[ pickle_list_name_temp.index('vz_eci_main_sc') ])
            if ( 'algebric_distance_ensemble_main_sc_lvlh_x' in pickle_list_name_temp ):
                algebric_distance_ensemble_main_sc_lvlh_x.append(pickle_list[ pickle_list_name_temp.index('algebric_distance_ensemble_main_sc_lvlh_x') ])
            if ( 'algebric_distance_ensemble_main_sc_lvlh_y' in pickle_list_name_temp ):
                algebric_distance_ensemble_main_sc_lvlh_y.append(pickle_list[ pickle_list_name_temp.index('algebric_distance_ensemble_main_sc_lvlh_y') ])
            if ( 'algebric_distance_ensemble_main_sc_lvlh_z' in pickle_list_name_temp ):
                algebric_distance_ensemble_main_sc_lvlh_z.append(pickle_list[ pickle_list_name_temp.index('algebric_distance_ensemble_main_sc_lvlh_z') ])
            if ( 'tca_ensemble' in pickle_list_name_temp ):
                tca_ensemble.append(pickle_list[ pickle_list_name_temp.index('tca_ensemble') ])
            if ( 'dca_ensemble' in pickle_list_name_temp ):
                dca_ensemble.append(pickle_list[ pickle_list_name_temp.index('dca_ensemble') ])


            if ipickle == 0:
                if (pickle_list_name_temp[1] != "cd_ensemble"):
                    nb_steps.append(pickle_list[1].shape[0]) #!!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
                elif len(pickle_list_name_temp) > 1:
                    nb_steps.append(pickle_list[2].shape[0]) #!!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA

    # Read the date for ensembles: different from the amin in case of collision here because to compute the r/v at TCA, i modified SpOCK on Big so that it outputs r/v of ensembles (not main sc) at time steps: 0, time step with tca in it - dt, time step with tca in it, time step with tca in it + dt.
    filename_ensemble_x = satellite_to_plot_path[0] + 'ensemble/ensemble_x_eci_' + satellite_to_plot[0]# ensmeble from sat 1 or 2 have the msae dates 
    file_ensemble_x = open(filename_ensemble_x)
    read_file_ensemble_x = file_ensemble_x.readlines()
    nb_header = 0
    while read_file_ensemble_x[nb_header].split()[0] != '#START':
        nb_header = nb_header + 1
    nb_header = nb_header + 1
    date_ens = []
    n_ens = len(read_file_ensemble_x) - nb_header

    for iline in range(n_ens):
        date_ens_temp = read_file_ensemble_x[iline + nb_header].split()[0] + ' ' + read_file_ensemble_x[iline + nb_header].split()[1]
        date_ens.append(datetime.strptime(date_ens_temp, "%Y/%m/%d %H:%M:%S"))
    date_ens = np.array(date_ens)
    where_before_tca = np.where( date_ens <= tca_date )[0]
    time_step_with_tca = date_ens[where_before_tca[-1]]
    i = where_before_tca[-1] # !!!!!!!!!!! neeeds to be the time step of tca -> i modified SpOCK on Big so that it outputs r/v of ensembles (not main sc) at time steps: 0, tca - dt, tca, tca + dt
    ## path of the folder where you want to store the results (pickle, image, video)
    # looking at main SSC works only if dt = dt_output in SpOCK                                                                                                                               
    if input_variables[find_in_read_input_order_variables(order_input_variables, 'dt')] != input_variables[find_in_read_input_order_variables(order_input_variables, 'dt_output')]:
        print "***! The output time step of SpOCK needs to be the same as its integration time step. The program will top. !***"; raise Exception
    i_main = np.where(date_main_datetime_round_sec == date_ens[i])[0][0] # same for sc 1 and 2                                                                                                
    pos_main_before = []
    pos_main_after = []
    pos_main_before.append(np.array([x_eci_main_sc[0][i_main], y_eci_main_sc[0][i_main], z_eci_main_sc[0][i_main]]))
    pos_main_after.append(np.array([x_eci_main_sc[0][i_main+1], y_eci_main_sc[0][i_main+1], z_eci_main_sc[0][i_main+1]]))
    pos_main_before.append(np.array([x_eci_main_sc[1][i_main], y_eci_main_sc[1][i_main], z_eci_main_sc[1][i_main]]))
    pos_main_after.append(np.array([x_eci_main_sc[1][i_main+1], y_eci_main_sc[1][i_main+1], z_eci_main_sc[1][i_main+1]]))

    vel_main_before = []
    vel_main_after = []
    vel_main_before.append(np.array([vx_eci_main_sc[0][i_main], vy_eci_main_sc[0][i_main], vz_eci_main_sc[0][i_main]]))
    vel_main_after.append(np.array([vx_eci_main_sc[0][i_main+1], vy_eci_main_sc[0][i_main+1], vz_eci_main_sc[0][i_main+1]]))
    vel_main_before.append(np.array([vx_eci_main_sc[1][i_main], vy_eci_main_sc[1][i_main], vz_eci_main_sc[1][i_main]]))
    vel_main_after.append(np.array([vx_eci_main_sc[1][i_main+1], vy_eci_main_sc[1][i_main+1], vz_eci_main_sc[1][i_main+1]]))




    ########################################################################################################################################################################################################################################################################################################
    # PLOTS ################################################################################################################################################################################################################################################################################################
    ########################################################################################################################################################################################################################################################################################################



    # Parameters of figures

    dr = 0.01

    
    if irun == 0:
        ax = fig.add_subplot(gs[0, 0], projection='3d')
    elif irun == 1:
        ax = fig.add_subplot(gs[0, 1], projection='3d')
    elif irun == 2:
        ax = fig.add_subplot(gs[1, 0], projection='3d')
    elif irun == 3:
        ax = fig.add_subplot(gs[1, 1], projection='3d')
    elif irun == 4:
        ax = fig.add_subplot(gs[2, 1], projection='3d')

#    quantile  = main_input_file_name.split('quantile')[1].split('.')[0].split('_')[0]
#    ax.set_title('Quant. ' + quantile + ' (' + tca[0].split('T')[1] + ')', weight = 'bold', fontsize  = fontsize_plot)

    ax.set_zlabel(z_label, weight = 'normal', fontsize  = fontsize_plot,labelpad=-8)
    ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot,labelpad=1)
    ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot,labelpad=1)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])

    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot/1.5, size = 10, width = 2, pad = 7) 
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    max_nb_sc_in_cell = 0

    dt_to_tca = ( tca_date - time_step_with_tca ).total_seconds() # nb of seconds from the time step that has tca and tca
    nens = nb_ensembles #20000
    min_x_plot = 100000
    min_y_plot = 100000
    min_z_plot = 100000
    max_x_plot = -100000
    max_y_plot = -100000
    max_z_plot = -100000
#    print dt_to_tca, tca_date, time_step_with_tca
    pos_main = []
    vel_main = []

    for isat in range(nb_spacecraft):
        if dt_to_tca != 0: # i corresponds to the time step before tca but if tca falls right at this time step then don't interpolate (could still do it but the problem is that the way SpOCK outputs the dates makes it that this step sould be the last one so step i+1 is not otuput
            x = x_eci_ensemble[isat][i, :nens] +( x_eci_ensemble[isat][i+1, :nens] - x_eci_ensemble[isat][i, :nens] ) / dt *  dt_to_tca 
            y = y_eci_ensemble[isat][i, :nens] +  ( y_eci_ensemble[isat][i+1, :nens] - y_eci_ensemble[isat][i, :nens] ) / dt * dt_to_tca
            z = z_eci_ensemble[isat][i, :nens] + ( z_eci_ensemble[isat][i+1, :nens] - z_eci_ensemble[isat][i, :nens] ) / dt * dt_to_tca
            pos_main.append( pos_main_before[isat] + ( pos_main_after[isat] - pos_main_before[isat] ) / dt *  dt_to_tca )
            vel_main.append( vel_main_before[isat] + ( vel_main_after[isat] - vel_main_before[isat] ) / dt *  dt_to_tca )
            
        else:
            x = x_eci_ensemble[isat][i, :nens]
            y = y_eci_ensemble[isat][i, :nens]
            z = z_eci_ensemble[isat][i, :nens]
            #ipdb.set_trace()
        if np.min(x) < min_x_plot:
            min_x_plot = np.min(x)
        if np.min(y) < min_y_plot:
            min_y_plot = np.min(y)
        if np.min(z) < min_z_plot:
            min_z_plot = np.min(z)
        if np.max(x) > max_x_plot:
            max_x_plot = np.max(x)
        if np.max(y) > max_y_plot:
            max_y_plot = np.max(y)
        if np.max(z) > max_z_plot:
            max_z_plot = np.max(z)


        points = np.transpose(np.array( [x,y,z] ))
        cell_dx = dr/8.#0.05 # 0.1
        cell_dy = dr/8.#0.05
        cell_dz = dr/8.#0.05

        nb_cell_x = (int)(np.abs(np.max(x) - np.min(x))/cell_dx) + 1 # nb of cells in x direction
        nb_cell_y = (int)(np.abs(np.max(y) - np.min(y))/cell_dy) + 1
        nb_cell_z = (int)(np.abs(np.max(z) - np.min(z))/cell_dz) + 1
        hist, binedges = np.histogramdd(points, bins = (nb_cell_x, nb_cell_y, nb_cell_z), normed=False) # np adapts the width of the cells so that the cells fall exactly at the edges of the distribution (if yo don't undersetand this sentence then jsut keep in mind that the width of the cells are not exactl equal to cell_dx, cell_dy, cell_dz)

        # color by number of sc by cell -> create a histogram from min(hist) to max(hist)
        nb_points = points.shape[0] #  (equal to len(x_eci_ensemble[isat][i, :]))
        nb_sc_in_cell_where_i_am = np.zeros([nb_points])
        cell_x_arr = binedges[0]
        cell_y_arr = binedges[1]
        cell_z_arr = binedges[2]
        cell_dx_actual = cell_x_arr[1] - cell_x_arr[0]
        cell_dy_actual = cell_y_arr[1] - cell_y_arr[0]
        cell_dz_actual = cell_z_arr[1] - cell_z_arr[0]
        for ipoint in range(nb_points):
            icell_x_where_i_am = (int) ( ( x[ipoint] - np.min(cell_x_arr) ) / cell_dx_actual )
            icell_y_where_i_am = (int) ( ( y[ipoint] - np.min(cell_y_arr) ) / cell_dy_actual )
            icell_z_where_i_am = (int) ( ( z[ipoint] - np.min(cell_z_arr) ) / cell_dz_actual )
            if ( icell_x_where_i_am == nb_cell_x ): # this can happen only for on sc: the one right at the edge since histogramdd builds the cells wo that the edge of the max one fall exactly at the max sc
                icell_x_where_i_am = nb_cell_x - 1
            if ( icell_x_where_i_am == nb_cell_x ): # this can happen only for on sc: the one right at the edge since histogramdd builds the cells wo that the edge of the max one fall exactly at the max sc
                icell_x_where_i_am = nb_cell_x - 1
            if ( icell_y_where_i_am == nb_cell_y ):
                icell_y_where_i_am = nb_cell_y - 1
            if ( icell_z_where_i_am == nb_cell_z ):
                icell_z_where_i_am = nb_cell_z - 1

            nb_sc_in_cell_where_i_am[ipoint] = hist[icell_x_where_i_am, icell_y_where_i_am, icell_z_where_i_am]
#        if ((isat == 0) & (quantile == '5')):
        vmax_colormap = np.max(nb_sc_in_cell_where_i_am)*1.#8005:50
        vmin_colormap = np.max(nb_sc_in_cell_where_i_am)*0#750
            
        if (isat == 0):
            x_center_plot = np.mean(x)
            y_center_plot = np.mean(y)
#        if isat == 1:
            z_center_plot = np.mean(z)
        if (color_by_nb_sc_per_bin == 1):
            sc = ax.scatter(x, y, z, c = nb_sc_in_cell_where_i_am, cmap = plt.cm.jet, vmin = vmin_colormap, vmax = vmax_colormap)
        else:
            sc = ax.scatter(x, y, z, c = color_arr[isat],s=10) #!!!! before so that color shows number of sc per bin: 
        if max_nb_sc_in_cell < np.max(nb_sc_in_cell_where_i_am):
            max_nb_sc_in_cell = np.max(nb_sc_in_cell_where_i_am)
    if (color_by_nb_sc_per_bin == 1):
        fig.colorbar(sc,ax = ax)#, norm = vmax_colormap)

    range_x_plot = max_x_plot - min_x_plot
    range_y_plot = max_y_plot - min_y_plot
    range_z_plot = max_z_plot - min_z_plot
    max_range_plot = np.max([range_x_plot, range_y_plot, range_z_plot])
    ax.set_xlim([min_x_plot, min_x_plot + max_range_plot])
    ax.set_ylim([min_y_plot, min_y_plot + max_range_plot])
    ax.set_zlim([min_z_plot, min_z_plot + max_range_plot])
    for isat in range(nb_spacecraft):
        if isat == 0:
            ax.quiver(pos_main[isat][0], pos_main[isat][1], pos_main[isat][2],
                      vel_main[isat][0]/ np.linalg.norm(vel_main[isat]) * max_range_plot / 2.5,
                      vel_main[isat][1]/ np.linalg.norm(vel_main[isat]) * max_range_plot / 2.5,
                      vel_main[isat][2]/ np.linalg.norm(vel_main[isat]) * max_range_plot / 2.5,
                      arrow_length_ratio=0.2,
                      color = 'k',#color_arr[isat],
                      linewidth = 2)

        else:
            ax.quiver(pos_main[isat][0], pos_main[isat][1], pos_main[isat][2],
                      vel_main[isat][0]/ np.linalg.norm(vel_main[isat]) * max_range_plot / 1.75,
                      vel_main[isat][1]/ np.linalg.norm(vel_main[isat]) * max_range_plot / 1.75,
                      vel_main[isat][2]/ np.linalg.norm(vel_main[isat]) * max_range_plot / 1.75,
                      arrow_length_ratio=0.15,
                      color = 'k',#color_arr[isat],
                      linewidth = 2)

    # print "\n\n"
    # print min_x_plot, max_x_plot
    # print min_y_plot, max_y_plot
    # print min_z_plot, max_z_plot

#     ax.set_xlim([x_center_plot-dr, x_center_plot + dr])
#     ax.set_ylim([y_center_plot-dr, y_center_plot + dr])
#     ax.set_zlim([z_center_plot-dr,z_center_plot + dr])

    #ax.view_init(elev=35., azim=0)
    #ax.view_init(elev=35., azim=0)
    ax.view_init(elev=10., azim=-235)
    ax.view_init(elev=10., azim=-235)
    #     ax.view_init(elev=90., azim=0)
    #     ax.view_init(elev=90., azim=0)

    ax.grid(False)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
fig.set_figheight(height_fig)
fig.set_figwidth(height_fig * ratio_fig_size)

fig_save_name = root_save_fig_name + main_input_file_name.replace(".txt","_collision_inter.pdf")# main_input_file_name.split('_quantile')[0] +  '_3d_try_again.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
#os.system("scp -p "  + fig_save_name + " desk:")


if show_plots == 1:
    plt.show(); plt.show()



    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # for isat in range(nb_spacecraft):
    #     ax.scatter(algebric_distance_ensemble_main_sc_lvlh_x[isat][i, :],algebric_distance_ensemble_main_sc_lvlh_y[isat][i,:],algebric_distance_ensemble_main_sc_lvlh_z[isat][i,:], color = color_arr[isat])
    # # ax.set_xlim([np.mean(algebric_distance_ensemble_main_sc_lvlh_x[isat][i, :])-dr, np.mean(algebric_distance_ensemble_main_sc_lvlh_x[isat][i, :]) + dr])
    # # ax.set_ylim([np.mean(algebric_distance_ensemble_main_sc_lvlh_y[isat][i, :])-dr, np.mean(algebric_distance_ensemble_main_sc_lvlh_y[isat][i, :]) + dr])
    # # ax.set_zlim([np.mean(algebric_distance_ensemble_main_sc_lvlh_z[isat][i, :])-dr, np.mean(algebric_distance_ensemble_main_sc_lvlh_z[isat][i, :]) + dr])





