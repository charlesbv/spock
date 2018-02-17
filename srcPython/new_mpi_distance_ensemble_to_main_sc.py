import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/kalman/spock_development_new_structure_kalman_dev/srcPython")
#sys.path.append("/home/cbv/spock_development_new_structure_kalman_dev/srcPython")
import matplotlib.gridspec as gridspec
from read_input_file import *
from matplotlib.colors import LogNorm
import pickle
from eci_to_lvlh import *
import fileinput
import time
import datetime
import numpy as np
from matplotlib import pyplot as plt
import os
import subprocess
from read_output_file import *
from orbit_average import *



plt.ion()
plt.isinteractive()

# from mpi4py import MPI
# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# size = comm.Get_size()
rank = 0
size = 1
## NOTE 1: to use this script, the only 3 parameters you have to set are:
## - the name of the propagator main input file 
## - if yes or no it is the first time that you read the output files of this simulation: set first_time to 1 if yes. Set it to 0 if you already read the output files and saved the results in a pickle
## - the path of the folder where you want to store the results (pickle, image, video): called 'path_folder_results'. In this folder, there must be the following subfolders: 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', and 'other'. In each of these subfolders, there must be the 2 subsubfolders: 'result', and 'pickle'. In the subsubfolder 'result', there must be the 2 subsubsubfolders: 'image', and 'video'.
## NOTE 2: this can be run if only ONE MAIN satellite was run (with ensembles of course)
## NOTE 3: the results (pickle, image, video) are stored in a subfolder of path_folder_results. The subfolder is either 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', or 'other'. This code figures out which folder the simulation corresponds to by reading the name of the output file chosen by the user in the main input file of the propagator (third line of section #SPACECRAFTS): it tries to find 'cygnss', 'cadre', 'aerie', 'scion', or 'qb50' in the name of the output file. If it does not find it, then the results here will be stored in path_folder_results/other/
## NOTE 4: to run this script, and any python script in the propagator, you need to be one subfolder deep from the main folder where the propagator runs are made. So if path_to_propagator/PropSim is the folder where the propagator runs are made, then the python scripts must be run for example in path_to_propagator/PropSim/subfolder_where_python_scipts_are_run

path_folder_results = '/Users/cbv/coll'

if (os.path.isdir(path_folder_results) == False):
    os.system("mkdir " + path_folder_results)

if path_folder_results[-1] != '!':
    path_folder_results = path_folder_results + '/'


main_input_file_name = sys.argv[1]
first_time = 1
# read input file
input_variables, order_input_variables = read_input_file(main_input_file_name)
dt = input_variables[find_in_read_input_order_variables(order_input_variables, 'dt')];
n_main = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_steps')];
satellite_to_plot_path = input_variables[find_in_read_input_order_variables(order_input_variables, 'output_file_path_list')];
satellite_to_plot = input_variables[find_in_read_input_order_variables(order_input_variables, 'output_file_name_list')];
nb_ensembles_coe = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_ensembles_coe')];
nb_ensembles_attitude = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_ensembles_attitude')];
nb_ensembles_cd = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_ensembles_cd')];
ensemble_to_plot_temp = input_variables[find_in_read_input_order_variables(order_input_variables, 'ensembles_to_output')];
nb_ensembles_density = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_ensembles_density')]
nb_spacecraft = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_sc')]
compute_drag = input_variables[find_in_read_input_order_variables(order_input_variables, 'compute_drag')]


dt_read_output = dt# in seconds
dt_read_output_in_time_step = (int)(dt_read_output / dt)
earth_mu    = 398600.4418; # gravitational parameter (km^3/s^2)
# Name of the mission of the simulation (choices: 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', 'other')
isat = 0




# Ensembles created by the propagator
ensemble_to_plot_from_main_input_file = []
for i in range(len(ensemble_to_plot_temp)):
    if (ensemble_to_plot_temp[i] == 'eci_r'):
        ensemble_to_plot_from_main_input_file.append('x_eci'); ensemble_to_plot_from_main_input_file.append('y_eci'); ensemble_to_plot_from_main_input_file.append('z_eci')
    if (ensemble_to_plot_temp[i] == 'eci_v'):
        ensemble_to_plot_from_main_input_file.append('vx_eci'); ensemble_to_plot_from_main_input_file.append('vy_eci'); ensemble_to_plot_from_main_input_file.append('vz_eci')
    if (ensemble_to_plot_temp[i] == 'geodetic'):
        ensemble_to_plot_from_main_input_file.append('longitude'); ensemble_to_plot_from_main_input_file.append('latitude'); ensemble_to_plot_from_main_input_file.append('altitude')
    if (ensemble_to_plot_temp[i] == 'power'):
        ensemble_to_plot_from_main_input_file.append('power')
    if (ensemble_to_plot_temp[i] == 'attitude'):
        ensemble_to_plot_from_main_input_file.append('pitch'); ensemble_to_plot_from_main_input_file.append('roll'); ensemble_to_plot_from_main_input_file.append('yaw')
    if (ensemble_to_plot_temp[i] == 'oe'):
        ensemble_to_plot_from_main_input_file.append('sma'); ensemble_to_plot_from_main_input_file.append('inclination'); ensemble_to_plot_from_main_input_file.append('eccentricity'); ensemble_to_plot_from_main_input_file.append('true_anomaly'); ensemble_to_plot_from_main_input_file.append('RAAN'); ensemble_to_plot_from_main_input_file.append('argument_perigee');ensemble_to_plot_from_main_input_file.append('phase_angle');ensemble_to_plot_from_main_input_file.append('sma_difference');

    if ((ensemble_to_plot_temp[i] == 'density') & (compute_drag == 1)):
        ensemble_to_plot_from_main_input_file.append('rho'); ensemble_to_plot_from_main_input_file.append('f107'); ensemble_to_plot_from_main_input_file.append('f107a'); ensemble_to_plot_from_main_input_file.append('ap'); 
    if (ensemble_to_plot_temp[i] == 'collision'):
        ensemble_to_plot_from_main_input_file.append('tca'); ensemble_to_plot_from_main_input_file.append('dca');


ensemble_to_plot = []
skip_arg = 2


# !!!!!!! after cbv messes up on dec 15 2017
rank_alt = 0
rank_lat = 1
rank_x_eci = 2
rank_y_eci = rank_x_eci
rank_z_eci = rank_x_eci
rank_sma = rank_lat
rank_inclination = 7
rank_eccentricity = 8
rank_true_anomaly = 5
rank_raan  = 9
rank_argument_perigee = rank_true_anomaly
rank_rho = rank_lat
rank_f107 = 10
rank_f107a = 11
rank_ap = 6
rank_angle_asc_node_to_sat = rank_true_anomaly
rank_tca = 3
rank_dca = rank_tca
for i in range(skip_arg, len(sys.argv)):
    ensemble_to_plot.append(sys.argv[i])
    if ( ( sys.argv[i] in ensemble_to_plot_from_main_input_file ) == False ):
        print "!*** " + sys.argv[i] + " was not an output of SpOCK. !***"
        ensemble_to_plot.remove(sys.argv[i])

# end of  after cbv messes up on dec 15 2017
# # !!!!!! before cbv messes up on dec 15 2017
# rank_alt = 0
# rank_lat = 1
# rank_x_eci = 2
# rank_y_eci = rank_x_eci
# rank_z_eci = rank_x_eci
# rank_sma = rank_lat
# rank_inclination = 7
# rank_eccentricity = 8
# rank_true_anomaly = 5
# rank_raan  = 9
# rank_argument_perigee = rank_true_anomaly
# rank_rho = rank_lat
# rank_f107 = 3
# rank_f107a = 4
# rank_ap = 6
# rank_angle_asc_node_to_sat = rank_true_anomaly
# rank_tca = 10
# rank_dca = 11

# if size < 10:
#     print "!*** You need to run with -np 12 at least. The program will now stop. !***"; raise Exception
# for i in range(skip_arg, len(sys.argv)):
#     ensemble_to_plot.append(sys.argv[i])
#     if ( ( sys.argv[i] in ensemble_to_plot_from_main_input_file ) == False ):
#         print "!*** " + sys.argv[i] + " was not an output of SpOCK. !***"
#         ensemble_to_plot.remove(sys.argv[i])

# # !!!!!! end before cbv messes up on dec 15 2017

rank_alt = 0
rank_lat = 0
rank_x_eci = 0
rank_y_eci = 0
rank_z_eci = 0
rank_sma = 0
rank_inclination = 0
rank_eccentricity = 0
rank_true_anomaly = 0
rank_raan = 0
rank_argument_perigee = 0
rank_rho = 0
rank_f107 = 0
rank_f107a = 0
rank_ap = 0
rank_angle_asc_node_to_sat = 0
rank_tca = 0
rank_dca = 0


# if 'latitude' in sys.argv:
#     rank_lat = 0
# if 'true_anomaly' in sys.argv:
#     rank_true_anomaly = 0
# for i in range(skip_arg, len(sys.argv)):
#     ensemble_to_plot.append(sys.argv[i])
#     if ( ( sys.argv[i] in ensemble_to_plot_from_main_input_file ) == False ):
#         print "!*** " + sys.argv[i] + " was not an output of SpOCK. !***"
#         ensemble_to_plot.remove(sys.argv[i])
#     else:
#         if sys.argv[i] == 'altitude':
#             rank_alt = (i - skip_arg)%size
#         if sys.argv[i] == 'x_eci':
#             rank_x_eci = (i - skip_arg)%size
#         if sys.argv[i] == 'y_eci':
#             rank_y_eci = (i - skip_arg)%size
#         if sys.argv[i] == 'z_eci':
#             rank_z_eci = (i - skip_arg)%size
#         if sys.argv[i] == 'sma':
#             if (('latitude' in sys.argv) == False):
#                 rank_sma =  (i - skip_arg)%size
#             else:
#                 rank_sma = rank_lat
#         if sys.argv[i] == 'inclination':
#             rank_inclination = (i - skip_arg)%size
#         if sys.argv[i] == 'eccentricity':
#             rank_eccentricity = (i - skip_arg)%size
#         if sys.argv[i] == 'true_anomaly':
#             rank_true_anomaly = (i - skip_arg)%size
#         if sys.argv[i] == 'raan':
#             rank_raan = (i - skip_arg)%size
#         if sys.argv[i] == 'argument_perigee':
#             if (('true_anomaly' in sys.argv) == False):
#                 rank_argument_perigee =  (i - skip_arg)%size
#             else:
#                 rank_argument_perigee = rank_true_anomaly
#                 rank_angle_asc_node_to_sat = rank_true_anomaly
#         if sys.argv[i] == 'rho':
#             if (('latitude' in sys.argv) == False):
#                 print "To compute the orbit density average, the latitudes of the ensembles need to be computed. This means that in the section #OUTPUT_ENSEMBLES of the main input file of the propagator, there has to be the word 'geodetic'. The program will stop."; raise Exception
#             else:
#                 rank_rho = rank_lat
#         if sys.argv[i] == 'f107':
#             rank_f107 = (i - skip_arg)%size
#         if sys.argv[i] == 'f107a':
#             rank_f107a = (i - skip_arg)%size
#         if sys.argv[i] == 'ap':
#             rank_ap = (i - skip_arg)%size

## Nb of ensembles
nb_ensembles_array = [nb_ensembles_coe, nb_ensembles_attitude, nb_ensembles_cd, nb_ensembles_density]
nb_ensembles = np.max(nb_ensembles_array)
for i in range(len(nb_ensembles_array)):
    if ( (nb_ensembles_array[i] > 0 ) &  (nb_ensembles_array[i] < nb_ensembles) ) :
        nb_ensembles = nb_ensembles_array[i]

period_ensemble_orbit_average_all_sc = []
sma_ensemble_orbit_average_all_sc = []
angle_asc_node_to_sat_ensemble_orbit_average_temp_all_sat = []
for isat in range(nb_spacecraft):
    pickle_list = []
    pickle_list_name = []

    # read output file of reference satellite
    list_var = ["latitude"]
    out_var, order_var = read_output_file( satellite_to_plot_path[isat] + satellite_to_plot[isat], list_var  )
    latitude_main = out_var[1]
    date_main = out_var[0]

    name_subfolder_save = satellite_to_plot[isat][:-5] + "/"
    if ( os.path.isdir(path_folder_results + 'pickle/') == False ):
        os.system("mkdir " + path_folder_results + 'pickle/')
    if ( os.path.isdir(path_folder_results + 'result/') == False ):
        os.system("mkdir " + path_folder_results + 'result/')

    os.system("mkdir " + path_folder_results + 'pickle/' + name_subfolder_save + " &>/dev/null")
    #    os.system("mkdir " + path_folder_results + name_mission + '/result/video/' + name_subfolder_save )

    save_pickle_name = path_folder_results +  'pickle/' + name_subfolder_save + satellite_to_plot[isat].replace(".txt","_") 
    root_save_fig_name = path_folder_results +  'result/' + name_subfolder_save + satellite_to_plot[isat].replace(".txt","_") 


    # # #######################################
    dir_final_output_ensemble = satellite_to_plot_path[isat] + 'ensemble/'
    if ( first_time == 1 ):
        # open any file to get the number of processors
        if len([f for f in os.listdir(dir_final_output_ensemble) if os.path.isfile(os.path.join(dir_final_output_ensemble, f))]) > 0:
            ex_ensemble_file = [f for f in os.listdir(dir_final_output_ensemble) if os.path.isfile(os.path.join(dir_final_output_ensemble, f))][0]
            nProcs = (int)(open(dir_final_output_ensemble + ex_ensemble_file).readlines()[5].split('-')[1].split()[0])
        elif os.path.exists('/'.join(satellite_to_plot_path[0].split('/')[0:-2]) + '/collision/dca/'):
            if len([f for f in os.listdir('/'.join(satellite_to_plot_path[0].split('/')[0:-2]) + '/collision/dca/')]) > 0:
                nProcs = (int)([f for f in os.listdir('/'.join(satellite_to_plot_path[0].split('/')[0:-2]) + '/collision/dca/')][0].split('-')[1].split('_')[0])
            else:
                'Could not find any processor file to concatenate. The program will stop.'; raise Exception;
        else:
            'Could not find any processor file to concatenate. The program will stop.'; raise Exception;

        # ex_ensemble_file = [f for f in os.listdir(dir_final_output_ensemble) if os.path.isfile(os.path.join(dir_final_output_ensemble, f))][0]
        # nProcs = (int)(open(dir_final_output_ensemble + ex_ensemble_file).readlines()[5].split('-')[1].split()[0])
        pickle_list.append(nProcs); pickle_list_name.append('nProcs')
        nb_ensembles= (int)(nb_ensembles / nProcs) * nProcs
        if rank  == rank_tca:
            if ( 'tca' in ensemble_to_plot ): 
                if isat == 0:
                    output_run_dir = '/'.join(satellite_to_plot_path[0].split('/')[0:-2])
                # Read the ensembles on tca
                    print "Read the ensembles on tca by iProc", rank
                    file_x = open(output_run_dir + "/collision/tca/ensemble_tca_" + output_run_dir.split('/')[-1] + ".txt", "r")
                    read_file_x = file_x.readlines()
                    tca_ensemble = []
                    for i in range(0,len(read_file_x[0].split())):#,len(read_file_x[0].split())/10000): # always nly one line in tca file
                        print np.float(i)/len(read_file_x[0].split()) * 100
                        tca_ensemble.append( float( read_file_x[0].split()[i] ))
                    file_x.close()
                    tca_ensemble = np.array(tca_ensemble)               
                    pickle_list.append(tca_ensemble); pickle_list_name.append('tca_ensemble')
                    save_pickle_name = save_pickle_name + 'tca_'
        if rank  == rank_dca:
            if ( 'dca' in ensemble_to_plot ): 
                if isat == 0:
                    output_run_dir = '/'.join(satellite_to_plot_path[0].split('/')[0:-2])
                # Read the ensembles on dca
                    print "Read the ensembles on dca by iProc", rank
                    file_x = open(output_run_dir + "/collision/dca/ensemble_dca_" + output_run_dir.split('/')[-1] + ".txt", "r")
                    read_file_x = file_x.readlines()
                    dca_ensemble = []
                    for i in range(0,len(read_file_x[0].split()),len(read_file_x[0].split())/10000): # always nly one line in dca file
                        print np.float(i)/len(read_file_x[0].split()) * 100
                        dca_ensemble.append( float( read_file_x[0].split()[i] ))
                    file_x.close()
                    dca_ensemble = np.array(dca_ensemble)               
                    pickle_list.append(dca_ensemble); pickle_list_name.append('dca_ensemble')
                    save_pickle_name = save_pickle_name + 'dca_'

        if rank  == rank_alt:
            if ( 'altitude' in ensemble_to_plot ): 
                # Read the ensembles on altitude
                print "Read the ensembles on altitude by iProc", rank
                file_x = open(dir_final_output_ensemble+"ensemble_altitude_"+satellite_to_plot[isat], "r")
                read_file_x = file_x.readlines()
                nb_lines_header = 11
                n = len(read_file_x) - nb_lines_header # !!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
                altitude_ensemble = np.zeros([n, nb_ensembles])
                for i in range(0,n,dt_read_output_in_time_step):
                    if rank == 0:
                        print("%d" % ( i/float(n) * 100 ) + " %")
                        sys.stdout.write("\033[F") # Cursor up one line
                    for j in range(nb_ensembles):
#                        print read_file_x[i + nb_lines_header], '|||', i,j,rank
                        altitude_ensemble[i,j] = float( read_file_x[i + nb_lines_header].split()[j + 2] )
                file_x.close()
                pickle_list.append(altitude_ensemble); pickle_list_name.append('altitude_ensemble')
                save_pickle_name = save_pickle_name + 'altitude_'
        if rank == rank_lat:
            if ( 'latitude' in ensemble_to_plot ): 
                # Read the ensembles on latitude
                print "Read the ensembles on latitude by iProc", rank
                file_x = open(dir_final_output_ensemble+"ensemble_latitude_"+satellite_to_plot[isat], "r")
                read_file_x = file_x.readlines()
                nb_lines_header = 11
                n = len(read_file_x) - nb_lines_header # !!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
                latitude_ensemble = np.zeros([n, nb_ensembles])
                for i in range(0,n,dt_read_output_in_time_step):
                    if rank == 0:
                        print("%d" % ( i/float(n) * 100 ) + " %")
                        sys.stdout.write("\033[F") # Cursor up one line
                    for j in range(nb_ensembles):
                        latitude_ensemble[i,j] = float( read_file_x[i + nb_lines_header].split()[j + 2] )
                file_x.close()
                pickle_list.append(latitude_ensemble); pickle_list_name.append('latitude_ensemble')
                save_pickle_name = save_pickle_name + 'latitude_'
        if rank == rank_x_eci:
            if ( 'x_eci' in ensemble_to_plot ): 
                # Read the ensembles on x_eci
                print "Read the ensembles on x_eci by iProc", rank
                file_x = open(dir_final_output_ensemble+"ensemble_x_eci_"+satellite_to_plot[isat], "r")
                read_file_x = file_x.readlines()
                nb_lines_header = 11
                n = len(read_file_x) - nb_lines_header # !!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA

                x_eci_ensemble = np.zeros([n, nb_ensembles])
#                print range(0,n,dt_read_output_in_time_step)
                for i in range(0,n,dt_read_output_in_time_step):
                    if rank == 0:
                        print("%d" % ( i/float(n) * 100 ) + " %")            
                        sys.stdout.write("\033[F") # Cursor up one line
                    for j in range(nb_ensembles):
                        x_eci_ensemble[i,j] = float( read_file_x[i + nb_lines_header].split()[j + 2] )
                file_x.close()
                pickle_list.append(x_eci_ensemble); pickle_list_name.append('x_eci_ensemble')
                save_pickle_name = save_pickle_name + 'x_eci_'
        if rank == rank_y_eci:
            if ( 'y_eci' in ensemble_to_plot ): 
                # Read the ensembles on y_eci
                print "Read the ensembles on y_eci by iProc", rank
                file_y = open(dir_final_output_ensemble+"ensemble_y_eci_"+satellite_to_plot[isat], "r")
                read_file_y = file_y.readlines()
                nb_lines_header = 11
                n = len(read_file_y) - nb_lines_header # !!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
                y_eci_ensemble = np.zeros([n, nb_ensembles])
                for i in range(0,n,dt_read_output_in_time_step):
                    if rank == 0:
                        print("%d" % ( i/float(n) * 100 ) + "%")
                        sys.stdout.write("\033[F") # Cursor up one line
                    for j in range(nb_ensembles):
                        y_eci_ensemble[i,j] = float( read_file_y[i + nb_lines_header].split()[j + 2] )
                file_y.close()
                pickle_list.append(y_eci_ensemble); pickle_list_name.append('y_eci_ensemble')
                save_pickle_name = save_pickle_name + 'y_eci_'
        if rank == rank_z_eci:
            if ( 'z_eci' in ensemble_to_plot ): 
                # Read the ensembles on z_eci
                print "Read the ensembles on z_eci by iProc", rank
                file_z = open(dir_final_output_ensemble+"ensemble_z_eci_"+satellite_to_plot[isat], "r")
                read_file_z = file_z.readlines()
                nb_lines_header = 11
                n = len(read_file_z) - nb_lines_header # !!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
                z_eci_ensemble = np.zeros([n, nb_ensembles])
                for i in range(0,n,dt_read_output_in_time_step):
                    if rank == 0:
                        print("%d" % ( i/float(n) * 100 ) + "%")
                        sys.stdout.write("\033[F") # Cursor up one line
                    for j in range(nb_ensembles):
                        z_eci_ensemble[i,j] = float( read_file_z[i + nb_lines_header].split()[j + 2] )
                file_z.close()
                pickle_list.append(z_eci_ensemble); pickle_list_name.append('z_eci_ensemble')
                save_pickle_name = save_pickle_name + 'z_eci_'
        # if ( 'pitch' in ensemble_to_plot ): 
        #     # Read the ensembles on pitch
        #     print "Read the ensembles on pitch"
        #     file_pitch = open(dir_final_output_ensemble+"ensemble_pitch_"+satellite_to_plot[isat], "r")
        #     read_file_pitch = file_pitch.readlines()
        #     nb_lines_header = 11
        #     n = len(read_file_x) - nb_lines_header # !!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
        #     pitch_ensemble = np.zeros([n, nb_ensembles])
        #     for i in range(0,n,dt_read_output_in_time_step):
        #         print("%d" % ( i/float(n) * 100 ) + "%")
        #             sys.stdout.write("\033[F") # Cursor up one line
        #         for j in range(nb_ensembles):
        #             pitch_ensemble[i,j] = float( read_file_pitch[i + nb_lines_header].split()[j + 2] )
        #     file_pitch.close()
        #     pickle_list.append(pitch_ensemble); pickle_list_name.append('pitch_ensemble')
        # if ( 'roll' in ensemble_to_plot ): 
        #     # Read the ensembles on roll
        #     print "Read the ensembles on roll"
        #     file_roll = open(dir_final_output_ensemble+"ensemble_roll_"+satellite_to_plot[isat], "r")
        #     read_file_roll = file_roll.readlines()
        #     nb_lines_header = 11
        #     n = len(read_file_x) - nb_lines_header # !!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
        #     roll_ensemble = np.zeros([n, nb_ensembles])
        #     for i in range(0,n,dt_read_output_in_time_step):
        #         print("%d" % ( i/float(n) * 100 ) + "%")
        #             sys.stdout.write("\033[F") # Cursor up one line
        #         for j in range(nb_ensembles):
        #             roll_ensemble[i,j] = float( read_file_roll[i + nb_lines_header].split()[j + 2] )
        #     file_roll.close()
        #     pickle_list.append(roll_ensemble); pickle_list_name.append('roll_ensemble')
        # if ( 'yaw' in ensemble_to_plot ): 
        #     # Read the ensembles on yaw
        #     print "Read the ensembles on yaw"
        #     file_yaw = open(dir_final_output_ensemble+"ensemble_yaw_"+satellite_to_plot[isat], "r")
        #     read_file_yaw = file_yaw.readlines()
        #     nb_lines_header = 11
        #     n = len(read_file_x) - nb_lines_header # !!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
        #     yaw_ensemble = np.zeros([n, nb_ensembles])
        #     for i in range(0,n,dt_read_output_in_time_step):
        #         print("%d" % ( i/float(n) * 100 ) + "%")
        #             sys.stdout.write("\033[F") # Cursor up one line
        #         for j in range(nb_ensembles):
        #             yaw_ensemble[i,j] = float( read_file_yaw[i + nb_lines_header].split()[j + 2] ) 
        #     file_yaw.close()
        #     pickle_list.append(yaw_ensemble); pickle_list_name.append('yaw_ensemble')
        if rank == rank_sma:
            if ( 'sma' in ensemble_to_plot ): 
                # Read the ensembles on sma
                print "Read the ensembles on sma by iProc", rank
                file_sma = open(dir_final_output_ensemble+"ensemble_sma_"+satellite_to_plot[isat], "r")
                read_file_sma = file_sma.readlines()
                nb_lines_header = 11
                n = len(read_file_sma) - nb_lines_header # !!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
                sma_ensemble = np.zeros([n, nb_ensembles])
                period_ensemble = np.zeros([n, nb_ensembles])
                for i in range(0,n,dt_read_output_in_time_step):
                    if rank == 0:
                        print("%d" % ( i/float(n) * 100 ) + "%")
                        sys.stdout.write("\033[F") # Cursor up one line
                    for j in range(nb_ensembles):
                        sma_ensemble[i,j] = float( read_file_sma[i + nb_lines_header].split()[j + 2] ) 
                        period_ensemble[i,j] = 2*np.pi * np.sqrt( sma_ensemble[i,j]**3 / earth_mu )

                # ORBIT AVERAGE PERIOD
                if ( ( 'latitude' in ensemble_to_plot ) == False):
                    print "To compute the orbit period average, the latitudes of the ensembles need to be computed. This means that in the section #OUTPUT_ENSEMBLES of the main input file of the propagator, there has to be the word 'geodetic'. The program will stop."; raise Exception
                if (dt > 300):
                    print "The output time step of the propagator (" + str(dt) + " s) has to be smaller than 5 minutes to compute the orbit average period. Therefore, the orbit average period is not computed.";
                else:
                    period_ensemble_orbit_average = []
                    for j in range(nb_ensembles):
                        period_ensemble_orbit_average_temp, time_orbit_average, index_time_orbit_average = orbit_average(period_ensemble[:, j], latitude_ensemble[:,j], date_main)
                        period_ensemble_orbit_average.append(period_ensemble_orbit_average_temp)
                    if ( 'phase_angle' in ensemble_to_plot ): # !!!only if 2 satellites
                        period_ensemble_orbit_average_all_sc.append(period_ensemble_orbit_average)
                        if isat == 1:# !!!! works only if 2 satellites
                            period_ensemble_orbit_average_all_sc = np.array(period_ensemble_orbit_average_all_sc)
                    period_ensemble_orbit_average = np.array(period_ensemble_orbit_average).transpose()
                    time_orbit_average = np.array(time_orbit_average)
                    index_time_orbit_average = np.array(index_time_orbit_average)
                    pickle_list.append(period_ensemble_orbit_average); pickle_list_name.append('period_ensemble_orbit_average')
                    pickle_list.append(time_orbit_average); pickle_list_name.append('time_orbit_average')
                    save_pickle_name = save_pickle_name + 'time_orbit_average_'
                    save_pickle_name = save_pickle_name + 'period_orbit_average_'
                    pickle_list.append(index_time_orbit_average); pickle_list_name.append('index_time_orbit_average')
                    save_pickle_name = save_pickle_name + 'index_time_orbit_average_'
                    if ( 'phase_angle' in ensemble_to_plot ): # !!!only if 2 satellites
                        if isat == 1:
                            print "  -> Calculate the phase angle by iProc", rank
                            nb_orbits = period_ensemble_orbit_average_all_sc.shape[2]
                            phase_angle_ensemble_orbit_average = np.zeros([nb_orbits, nb_ensembles**2])
                            for iorbit in range(nb_orbits):
                                for isc1 in range(nb_ensembles):
                                    for isc2 in range(nb_ensembles):
                                        phase_angle_ensemble_orbit_average[iorbit, isc1*nb_ensembles + isc2] = period_ensemble_orbit_average_all_sc[1, isc2, iorbit] - period_ensemble_orbit_average_all_sc[0, isc1, iorbit] # pahse angle is defined as the difference between the period of a spacecraft isc2 and the period of a spacecraft isc1
                            save_pickle_name = save_pickle_name + 'phase_angle_orbit_average_'
                            pickle_list.append(phase_angle_ensemble_orbit_average); pickle_list_name.append('phase_angle_ensemble_orbit_average')

#                        min_nb_period = min([])
                        
                # ORBIT AVERAGE SMA
                if ( ( 'latitude' in ensemble_to_plot ) == False):
                    print "To compute the orbit sma average, the latitudes of the ensembles need to be computed. This means that in the section #OUTPUT_ENSEMBLES of the main input file of the propagator, there has to be the word 'geodetic'. The program will stop."; raise Exception
                if (dt > 300):
                    print "The output time step of the propagator (" + str(dt) + " s) has to be smaller than 5 minutes to compute the orbit average sma. Therefore, the orbit average sma is not computed.";
                else:
                    sma_ensemble_orbit_average = []
                    for j in range(nb_ensembles):
                        sma_ensemble_orbit_average_temp, time_orbit_average, index_time_orbit_average = orbit_average(sma_ensemble[:, j], latitude_ensemble[:,j], date_main)
                        sma_ensemble_orbit_average.append(sma_ensemble_orbit_average_temp)
                    if ( 'sma_difference' in ensemble_to_plot ): # !!!only if 2 satellites
                        sma_ensemble_orbit_average_all_sc.append(sma_ensemble_orbit_average)
                        if isat == 1:# !!!! works only if 2 satellites
                            sma_ensemble_orbit_average_all_sc = np.array(sma_ensemble_orbit_average_all_sc)

                    sma_ensemble_orbit_average = np.array(sma_ensemble_orbit_average).transpose()
                    save_pickle_name = save_pickle_name + 'sma_orbit_average_'
                    time_orbit_average = np.array(time_orbit_average)
                    index_time_orbit_average = np.array(index_time_orbit_average)
                    pickle_list.append(sma_ensemble_orbit_average); pickle_list_name.append('sma_ensemble_orbit_average')
                    if ('time_orbit_average' in pickle_list_name) == False:
                        pickle_list.append(time_orbit_average); pickle_list_name.append('time_orbit_average')
                        save_pickle_name = save_pickle_name + 'time_orbit_average_'
                    if ('index_time_orbit_average' in pickle_list_name) == False:
                        pickle_list.append(index_time_orbit_average); pickle_list_name.append('index_time_orbit_average')
                        save_pickle_name = save_pickle_name + 'index_time_orbit_average_'
 

                    if ( 'sma_difference' in ensemble_to_plot ): # !!!only if 2 satellites
                        if isat == 1:
                            print "  -> Calculate the sma difference by iProc", rank
                            nb_orbits = period_ensemble_orbit_average_all_sc.shape[2]
                            sma_difference_ensemble_orbit_average = np.zeros([nb_orbits, nb_ensembles**2])
                            for iorbit in range(nb_orbits):
                                for isc1 in range(nb_ensembles):
                                    for isc2 in range(nb_ensembles):
                                        sma_difference_ensemble_orbit_average[iorbit, isc1*nb_ensembles + isc2] = sma_ensemble_orbit_average_all_sc[1, isc2, iorbit] - sma_ensemble_orbit_average_all_sc[0, isc1, iorbit] # pahse angle is defined as the difference between the sma of a spacecraft isc2 and the sma of a spacecraft isc1
                            save_pickle_name = save_pickle_name + 'sma_difference_orbit_average_'
                            pickle_list.append(sma_difference_ensemble_orbit_average); pickle_list_name.append('sma_difference_ensemble_orbit_average')

                file_sma.close()
                pickle_list.append(sma_ensemble); pickle_list_name.append('sma_ensemble')
                save_pickle_name = save_pickle_name + 'sma_'
                pickle_list.append(period_ensemble); pickle_list_name.append('period_ensemble')
                save_pickle_name = save_pickle_name + 'period_'

        if rank == rank_inclination:
            if ( 'inclination' in ensemble_to_plot ): 
                # Read the ensembles on inclination
                print "Read the ensembles on inclination by iProc", rank
                file_inclination = open(dir_final_output_ensemble+"ensemble_inclination_"+satellite_to_plot[isat], "r")
                read_file_inclination = file_inclination.readlines()
                nb_lines_header = 11
                n = len(read_file_inclination) - nb_lines_header # !!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
                inclination_ensemble = np.zeros([n, nb_ensembles])
                for i in range(0,n,dt_read_output_in_time_step):
                    if rank == 0:
                        print("%d" % ( i/float(n) * 100 ) + "%")
                        sys.stdout.write("\033[F") # Cursor up one line
                    for j in range(nb_ensembles):
                        inclination_ensemble[i,j] = float( read_file_inclination[i + nb_lines_header].split()[j + 2] ) 
                file_inclination.close()
                pickle_list.append(inclination_ensemble); pickle_list_name.append('inclination_ensemble')
                save_pickle_name = save_pickle_name + 'inclination_'
        if rank == rank_eccentricity:
            if ( 'eccentricity' in ensemble_to_plot ): 
                # Read the ensembles on eccentricity
                print "Read the ensembles on eccentricity by iProc", rank
                file_eccentricity = open(dir_final_output_ensemble+"ensemble_eccentricity_"+satellite_to_plot[isat], "r")
                read_file_eccentricity = file_eccentricity.readlines()
                nb_lines_header = 11
                n = len(read_file_eccentricity) - nb_lines_header # !!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
                eccentricity_ensemble = np.zeros([n, nb_ensembles])
                for i in range(0,n,dt_read_output_in_time_step):
                    if rank == 0:
                        print("%d" % ( i/float(n) * 100 ) + "%")
                        sys.stdout.write("\033[F") # Cursor up one line
                    for j in range(nb_ensembles):
                        eccentricity_ensemble[i,j] = float( read_file_eccentricity[i + nb_lines_header].split()[j + 2] ) 
                file_eccentricity.close()
                pickle_list.append(eccentricity_ensemble); pickle_list_name.append('eccentricity_ensemble')
                save_pickle_name = save_pickle_name + 'eccentricity_'
        if rank == rank_true_anomaly:
            if ( 'true_anomaly' in ensemble_to_plot ): 
                # Read the ensembles on true_anomaly
                print "Read the ensembles on true_anomaly by iProc", rank
                file_true_anomaly = open(dir_final_output_ensemble+"ensemble_true_anomaly_"+satellite_to_plot[isat], "r")
                read_file_true_anomaly = file_true_anomaly.readlines()
                nb_lines_header = 11
                n = len(read_file_true_anomaly) - nb_lines_header # !!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
                true_anomaly_ensemble = np.zeros([n, nb_ensembles])
                for i in range(0,n,dt_read_output_in_time_step):
                    if rank == 0:
                        print("%d" % ( i/float(n) * 100 ) + "%")
                        sys.stdout.write("\033[F") # Cursor up one line
                    for j in range(nb_ensembles):
                        true_anomaly_ensemble[i,j] = float( read_file_true_anomaly[i + nb_lines_header].split()[j + 2] ) 
                file_true_anomaly.close()
                pickle_list.append(true_anomaly_ensemble); pickle_list_name.append('true_anomaly_ensemble')
                save_pickle_name = save_pickle_name + 'true_anomaly_'
        if rank == rank_raan:
            if ( 'RAAN' in ensemble_to_plot ): 
                # Read the ensembles on RAAN
                print "Read the ensembles on RAAN by iProc", rank
                file_RAAN = open(dir_final_output_ensemble+"ensemble_RAAN_"+satellite_to_plot[isat], "r")
                read_file_RAAN = file_RAAN.readlines()
                nb_lines_header = 11
                n = len(read_file_RAAN) - nb_lines_header # !!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
                RAAN_ensemble = np.zeros([n, nb_ensembles])
                for i in range(0,n,dt_read_output_in_time_step):
                    if rank == 0:
                        print("%d" % ( i/float(n) * 100 ) + "%")
                        sys.stdout.write("\033[F") # Cursor up one line
                    for j in range(nb_ensembles):
                        RAAN_ensemble[i,j] = float( read_file_RAAN[i + nb_lines_header].split()[j + 2] ) 
                file_RAAN.close()
                pickle_list.append(RAAN_ensemble); pickle_list_name.append('RAAN_ensemble')
                save_pickle_name = save_pickle_name + 'RAAN_'
        if rank == rank_argument_perigee:
            if ( 'argument_perigee' in ensemble_to_plot ): 
                # Read the ensembles on argument_perigee
                print "Read the ensembles on argument_perigee by iProc", rank
                file_argument_perigee = open(dir_final_output_ensemble+"ensemble_argument_perigee_"+satellite_to_plot[isat], "r")
                read_file_argument_perigee = file_argument_perigee.readlines()
                nb_lines_header = 11
                n = len(read_file_argument_perigee) - nb_lines_header # !!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
                argument_perigee_ensemble = np.zeros([n, nb_ensembles])
                for i in range(0,n,dt_read_output_in_time_step):
                    if rank == 0:
                        print("%d" % ( i/float(n) * 100 ) + "%")
                        sys.stdout.write("\033[F") # Cursor up one line
                    for j in range(nb_ensembles):
                        argument_perigee_ensemble[i,j] = float( read_file_argument_perigee[i + nb_lines_header].split()[j + 2] ) 
                file_argument_perigee.close()
                pickle_list.append(argument_perigee_ensemble); pickle_list_name.append('argument_perigee_ensemble')
                save_pickle_name = save_pickle_name + 'argument_perigee_'
        if rank == rank_rho:
            if ( 'rho' in ensemble_to_plot ): 
                # Density orbit-average (need latitude for this so need 'geodetic' in #OUTPUT_ENSEMBLES of propagator main input file)
                if ( ( 'latitude' in ensemble_to_plot ) == False):
                    print "To compute the orbit density average, the latitudes of the ensembles need to be computed. This means that in the section #OUTPUT_ENSEMBLES of the main input file of the propagator, there has to be the word 'geodetic'. The program will stop."; raise Exception
                if (dt > 300):
                    print "The output time step of the propagator (" + str(dt) + " s) has to be smaller than 5 minutes to compute the orbit average density. Therefore, the orbit average density is not computed.";
                else:
                    # Read the ensembles on rho
                    print "Read the ensembles on rho by iProc", rank
                    file_rho = open(dir_final_output_ensemble+"ensemble_rho_"+satellite_to_plot[isat], "r")
                    read_file_rho = file_rho.readlines()
                    nb_lines_header = 11
                    n = len(read_file_rho) - nb_lines_header # !!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
                    rho_ensemble = np.zeros([n, nb_ensembles])
                    for i in range(0,n,dt_read_output_in_time_step):
                        if rank == 0:
                            print("%d" % ( i/float(n) * 100 ) + "%")
                            sys.stdout.write("\033[F") # Cursor up one line
                        for j in range(nb_ensembles):
                            rho_ensemble[i,j] = float( read_file_rho[i + nb_lines_header].split()[j + 2] )*10**(-9) # conversion from km^-3 to m^-3 
                            if (i == 1):
                                rho_ensemble[0,j] = rho_ensemble[1,j]  #!!!!! this is because the first line of the output files of SpOCK for F10.7, Ap, F10.7A, and density has 0 everywhere (this is due to the fact that the first line is printed before calling the function propagate_spacecraft in which the values of these parameters are set). This is not correct to do that but it's fine because we baiscally assume f10,7, ap, f10,7a and density to be constant for the 2 first time steps

                    file_rho.close()
                    pickle_list.append(rho_ensemble); pickle_list_name.append('rho_ensemble')
                    save_pickle_name = save_pickle_name + 'rho_'
                    rho_ensemble_orbit_average = []
                    for j in range(nb_ensembles):
                        rho_ensemble_orbit_average_temp, time_orbit_average, index_time_orbit_average = orbit_average(rho_ensemble[:, j], latitude_ensemble[:,j], date_main)
                        rho_ensemble_orbit_average.append(rho_ensemble_orbit_average_temp)
                    rho_ensemble_orbit_average = np.array(rho_ensemble_orbit_average).transpose()
                    time_orbit_average = np.array(time_orbit_average)
                    index_time_orbit_average = np.array(index_time_orbit_average)
                    pickle_list.append(rho_ensemble_orbit_average); pickle_list_name.append('rho_ensemble_orbit_average')
                    save_pickle_name = save_pickle_name + 'rho_orbit_average_'
                    if ('time_orbit_average' in pickle_list_name) == False:
                        pickle_list.append(time_orbit_average); pickle_list_name.append('time_orbit_average')
                        save_pickle_name = save_pickle_name + 'time_orbit_average_'
                    if ('index_time_orbit_average' in pickle_list_name) == False:
                        pickle_list.append(index_time_orbit_average); pickle_list_name.append('index_time_orbit_average')
                        save_pickle_name = save_pickle_name + 'index_time_orbit_average_'

 
        if rank == rank_f107:
            if ( 'f107' in ensemble_to_plot ): 
                # Read the ensembles on f107
                print "Read the ensembles on f107 by iProc", rank
                file_f107 = open(dir_final_output_ensemble+"ensemble_f107_"+satellite_to_plot[isat], "r")
                read_file_f107 = file_f107.readlines()
                nb_lines_header = 11
                n = len(read_file_f107) - nb_lines_header # !!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
                f107_ensemble = np.zeros([n, nb_ensembles])
                for i in range(0,n,dt_read_output_in_time_step):
                    if rank == 0:
                        print("%d" % ( i/float(n) * 100 ) + "%")
                        sys.stdout.write("\033[F") # Cursor up one line
                    for j in range(nb_ensembles):
                        f107_ensemble[i,j] = float( read_file_f107[i + nb_lines_header].split()[j + 2] ) 
                        if (i == 1):
                            f107_ensemble[0,j] = f107_ensemble[1,j]  #!!!!! this is because the first line of the output files of SpOCK for F10.7, Ap, F10.7A, and density has 0 everywhere (this is due to the fact that the first line is printed before calling the function propagate_spacecraft in which the values of these parameters are set). This is not correct to do that but it's fine because we baiscally assume f10,7, ap, f10,7a and density to be constant for the 2 first time steps
                file_f107.close()
                pickle_list.append(f107_ensemble); pickle_list_name.append('f107_ensemble')
                save_pickle_name = save_pickle_name + 'f107_'
        if rank == rank_f107a:
            if ( 'f107a' in ensemble_to_plot ): 
                # Read the ensembles on f107a
                print "Read the ensembles on f107a by iProc", rank
                file_f107a = open(dir_final_output_ensemble+"ensemble_f107a_"+satellite_to_plot[isat], "r")
                read_file_f107a = file_f107a.readlines()
                nb_lines_header = 11
                n = len(read_file_f107a) - nb_lines_header # !!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
                f107a_ensemble = np.zeros([n, nb_ensembles])
                for i in range(0,n,dt_read_output_in_time_step):
                    if rank == 0:
                        print("%d" % ( i/float(n) * 100 ) + "%")
                        sys.stdout.write("\033[F") # Cursor up one line
                    for j in range(nb_ensembles):
                        f107a_ensemble[i,j] = float( read_file_f107a[i + nb_lines_header].split()[j + 2] ) 
                        if (i == 1):
                            f107a_ensemble[0,j] = f107a_ensemble[1,j]  #!!!!! this is because the first line of the output files of SpOCK for F10.7, Ap, F10.7A, and density has 0 everywhere (this is due to the fact that the first line is printed before calling the function propagate_spacecraft in which the values of these parameters are set). This is not correct to do that but it's fine because we baiscally assume f10,7, ap, f10,7a and density to be constant for the 2 first time steps

                file_f107a.close()
                pickle_list.append(f107a_ensemble); pickle_list_name.append('f107a_ensemble')
                save_pickle_name = save_pickle_name + 'f107a_'
        if rank == rank_ap:
            if ( 'ap' in ensemble_to_plot ): 
                # Read the ensembles on ap
                print "Read the ensembles on ap by iProc", rank
                file_ap = open(dir_final_output_ensemble+"ensemble_ap_"+satellite_to_plot[isat], "r")
                read_file_ap = file_ap.readlines()
                nb_lines_header = 11     
                n = len(read_file_ap) - nb_lines_header # !!!!!!!!!!!this has been added between when computing collisions the time of output stop at end of time spanning TCA
                ap_ensemble = np.zeros([n, nb_ensembles])
                for i in range(0,n,dt_read_output_in_time_step):
                    if rank == 0:
                        print("%d" % ( i/float(n) * 100 ) + "%")
                        sys.stdout.write("\033[F") # Cursor up one line
                    for j in range(nb_ensembles):
                        ap_ensemble[i,j] = float( read_file_ap[i + nb_lines_header].split()[j + 2] ) 
                        if (i == 1):
                            ap_ensemble[0,j] = ap_ensemble[1,j]  #!!!!! this is because the first line of the output files of SpOCK for F10.7, Ap, F10.7A, and density has 0 everywhere (this is due to the fact that the first line is printed before calling the function propagate_spacecraft in which the values of these parameters are set). This is not correct to do that but it's fine because we baiscally assume f10,7, ap, f10,7a and density to be constant for the 2 first time steps

                file_ap.close()
                pickle_list.append(ap_ensemble); pickle_list_name.append('ap_ensemble')
                save_pickle_name = save_pickle_name + 'ap_'

        if ( ( 'x_eci' in ensemble_to_plot ) & ( 'y_eci' in ensemble_to_plot ) & ( 'z_eci' in ensemble_to_plot ) ): 
            if rank == rank_x_eci:
                # Read the main spacecraft x_eci, y_eci, and z_eci
                print "Read the the main spacecraft positions"
                file_main_sc = open(satellite_to_plot_path[isat] + satellite_to_plot[isat], "r")
                read_file_main_sc = file_main_sc.readlines()
                nb_lines_header = 10
                x_eci_main_sc = np.zeros(n_main); y_eci_main_sc = np.zeros(n_main); z_eci_main_sc = np.zeros(n_main)
                vx_eci_main_sc = np.zeros(n_main); vy_eci_main_sc = np.zeros(n_main); vz_eci_main_sc = np.zeros(n_main)
                for i in range(0,n_main,dt_read_output_in_time_step):
                    x_eci_main_sc[i] = float( read_file_main_sc[i + nb_lines_header].split()[2] )
                    y_eci_main_sc[i] = float( read_file_main_sc[i + nb_lines_header].split()[3] )
                    z_eci_main_sc[i] = float( read_file_main_sc[i + nb_lines_header].split()[4] )
                    vx_eci_main_sc[i] = float( read_file_main_sc[i + nb_lines_header].split()[5] )
                    vy_eci_main_sc[i] = float( read_file_main_sc[i + nb_lines_header].split()[6] )
                    vz_eci_main_sc[i] = float( read_file_main_sc[i + nb_lines_header].split()[7] )
                file_main_sc.close()
                pickle_list.append(x_eci_main_sc); pickle_list.append(y_eci_main_sc); pickle_list.append(z_eci_main_sc);    pickle_list_name.append('x_eci_main_sc'); pickle_list_name.append('y_eci_main_sc'); pickle_list_name.append('z_eci_main_sc');
                pickle_list.append(vx_eci_main_sc); pickle_list.append(vy_eci_main_sc); pickle_list.append(vz_eci_main_sc);    pickle_list_name.append('vx_eci_main_sc'); pickle_list_name.append('vy_eci_main_sc'); pickle_list_name.append('vz_eci_main_sc');

                # Distance between the ensemble sc and the main sc
                print "Distance between the ensemble sc and the main sc"
                distance_ensemble_main_sc = np.zeros([n, nb_ensembles])
                algebric_distance_ensemble_main_sc_lvlh_x = np.zeros([n, nb_ensembles])
                algebric_distance_ensemble_main_sc_lvlh_y = np.zeros([n, nb_ensembles])
                algebric_distance_ensemble_main_sc_lvlh_z = np.zeros([n, nb_ensembles])

                for i in range(0,np.min([n, n_main]),dt_read_output_in_time_step):
                    print("%d" % ( i/float(n) * 100 ) + "%")
                    sys.stdout.write("\033[F") # Cursor up one line
                    for j in range(nb_ensembles):
                        distance_ensemble_main_sc[i,j] = np.sqrt( (x_eci_ensemble[i,j] - x_eci_main_sc[i])**2 + (y_eci_ensemble[i,j] - y_eci_main_sc[i])**2 +  (z_eci_ensemble[i,j] - z_eci_main_sc[i])**2 )#*1000 
                        # * 1000 to convert the distance from km to m
                        position_main_sc_eci = np.array([x_eci_main_sc[i], y_eci_main_sc[i], z_eci_main_sc[i]])
                        position_ensemble_sc_eci = np.array([x_eci_ensemble[i,j], y_eci_ensemble[i,j], z_eci_ensemble[i,j]])
                        delta_vector_eci = position_ensemble_sc_eci - position_main_sc_eci #position_main_sc_eci - position_ensemble_sc_eci
                        velocity_main_sc_eci = np.array([vx_eci_main_sc[i], vy_eci_main_sc[i], vz_eci_main_sc[i]])
                        delta_vector_lvlh_of_main_sc = eci_to_lvlh(position_main_sc_eci, velocity_main_sc_eci, delta_vector_eci)
                        algebric_distance_ensemble_main_sc_lvlh_x[i,j]= delta_vector_lvlh_of_main_sc[0] #* 1000 # * 1000 to convert the distance from km to m 
                        algebric_distance_ensemble_main_sc_lvlh_y[i,j]= delta_vector_lvlh_of_main_sc[1] #* 1000 # * 1000 to convert the distance from km to m 
                        algebric_distance_ensemble_main_sc_lvlh_z[i,j]= delta_vector_lvlh_of_main_sc[2] #* 1000 # * 1000 to convert the distance from km to m 
                pickle_list.append(algebric_distance_ensemble_main_sc_lvlh_x); pickle_list.append(algebric_distance_ensemble_main_sc_lvlh_y); pickle_list.append(algebric_distance_ensemble_main_sc_lvlh_z);    pickle_list_name.append('algebric_distance_ensemble_main_sc_lvlh_x'); pickle_list_name.append('algebric_distance_ensemble_main_sc_lvlh_y'); pickle_list_name.append('algebric_distance_ensemble_main_sc_lvlh_z');

        if rank == rank_angle_asc_node_to_sat:
            if ( ( 'true_anomaly' in ensemble_to_plot ) & ( 'argument_perigee' in ensemble_to_plot ) ): # NOTE: the user needs to output the sma too because we use sma to convert from an angle (degree) to a distance (km). The sma we use is the median of all the sma of each ensemble, at a give time
                if isat == 0:
                    angle_asc_node_to_sat_ensemble = np.zeros([nb_spacecraft, n, nb_ensembles])
                    angle_asc_node_to_sat_cumulative_ensemble = np.zeros([n, nb_ensembles]) # by cumulative, we mean the following. Let's call alpha the angle between the ascending node and the satellite. Say that the satellite starts at an angle alpha of 90 degrees. Each time that the satellite goes over an entire orbit, it comes back back to this initial alpha (90 degrees). So after the first orbit, it'll have alpha equal to 90 degrees. Same after two orbits, three, etc.. To take into account the number of orbits it has made, we add to alpha 360 at each orbit and we call this angle angle_asc_node_to_sat_cumulative_ensemble. For example, after one orbit: angle_asc_node_to_sat_cumulative_ensemble = 90, after 2 orbits: angle_asc_node_to_sat_cumulative_ensemble = 90 + 360, after n orbits angle_asc_node_to_sat_cumulative_ensemble = 90 + 360 * ( n - 1 ).
                angular_spacing_between_ensemble_sat_temp = np.zeros([n, nb_ensembles-1]) # store the spacing between all satellites except between the one that has the max angle_asc_node_to_sat_ensemble and the one that has the min angle_asc_node_to_sat_ensemble (reason: it the min is 200 and the max is 220 (all the sc are clustered in a 20 degree interval) then we want to ignore the spacing of 340 degrees between the sc at 220 and the sc at 200)
                angular_spacing_between_ensemble_sat = np.zeros([n, nb_ensembles-2]) # if for example half of the sc are 10 degrees south of the AN (called subcluster_south) and the other half 10 degree north of the AN (called subcluster_north) then there is a jump of 340 degrees between the most northern sc of subcluster_north and the mosst souther sc of subcluster_south. We don't want to take this one into account. angular_spacing_between_ensemble_sat does not take it into account

                # UNCOMMENT IF WANT TO COMPUTE index_dt_angular_distance_between_two_clusters
                # if isat == 1:# !!!! WORKS ONLY IF ONE SATELLITE FOR NOW  
                #     dt_angular_distance_between_two_clusters = 95. # we look at angular distance between the 2 clusters only at certain times otherwise it's a too heavy calculation. here in MINUTES
                #     index_dt_angular_distance_between_two_clusters = (int)( dt_angular_distance_between_two_clusters * 60 / dt )
                #     angular_distance_between_two_clusters = np.zeros([(int)(nb_steps / index_dt_angular_distance_between_two_clusters), nb_ensembles*nb_ensembles])
                #     save_index_dt_angular_distance_between_two_clusters = []
                    # print "Computing the angular distance between each cluster..."
                # end of UNCOMMENT IF WANT TO COMPUTE index_dt_angular_distance_between_two_clusters

#                angular_spacing_between_ensemble_sat_converted_in_a_distance = np.zeros([n, nb_ensembles-2]) # if for example half of the sc are 10 degrees south of the AN (called subcluster_south) and the other half 10 degree north of the AN (called subcluster_north) then there is a jump of 340 degrees between the most northern sc of subcluster_north and the mosst souther sc of subcluster_south. We don't want to take this one into account. angular_spacing_between_ensemble_sat_converted_in_a_distance does not take it into account
    
                for i in range(0,n,dt_read_output_in_time_step):
                    angle_asc_node_to_sat_ensemble[isat, i,:] = (true_anomaly_ensemble[i,:] + argument_perigee_ensemble[i,:])%360            
                    # angle_asc_node_to_sat_ensemble_orbit_average = []
                # OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD
                # for j in range(nb_ensembles):
                #     angle_asc_node_to_sat_ensemble_orbit_average_temp, time_orbit_average, index_time_orbit_average = orbit_average(angle_asc_node_to_sat_ensemble[isat,:, j], latitude_ensemble[:,j], date_main)
                #     angle_asc_node_to_sat_ensemble_orbit_average.append(angle_asc_node_to_sat_ensemble_orbit_average_temp)
                # angle_asc_node_to_sat_ensemble_orbit_average_temp_all_sat.append(angle_asc_node_to_sat_ensemble_orbit_average)
                # angle_asc_node_to_sat_ensemble_orbit_average_all_sat = np.array(angle_asc_node_to_sat_ensemble_orbit_average_temp_all_sat).transpose()
                # pickle_list.append(angle_asc_node_to_sat_ensemble_orbit_average_all_sat); pickle_list_name.append('angle_asc_node_to_sat_ensemble_orbit_average_all_sat')
                # save_pickle_name = save_pickle_name + 'angle_asc_node_to_sat_ensemble_orbit_average_all_sat_'
                # if ('time_orbit_average' in pickle_list_name) == False        :
                #     time_orbit_average = np.array(time_orbit_average)
                #     pickle_list.append(time_orbit_average); pickle_list_name.append('time_orbit_average')
                #     save_pickle_name = save_pickle_name + 'time_orbit_average_'
                # if ('index_time_orbit_average' in pickle_list_name) == False:
                #     index_time_orbit_average = np.array(index_time_orbit_average)
                #     pickle_list.append(index_time_orbit_average); pickle_list_name.append('index_time_orbit_average')
                #     save_pickle_name = save_pickle_name + 'index_time_orbit_average_'
                # end of OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD
                    
                    # UNCOMMENT IF WANT TO COMPUTE angular_distance_between_two_clusters
#                     if isat == 1: # !!!! WORKS ONLY IF ONE SATELLITE FOR NOW
#                         if ((i%index_dt_angular_distance_between_two_clusters == 0) & (i/index_dt_angular_distance_between_two_clusters < angular_distance_between_two_clusters.shape[0])):
#                             print("%d" % ( i/float(n) * 100 ) + "%")
#                             sys.stdout.write("\033[F") # Cursor up one line  
                            
#                             for isc1 in range(nb_ensembles):
#                                 angular_distance_between_two_clusters[i/index_dt_angular_distance_between_two_clusters, isc1*nb_ensembles:(isc1+1)*nb_ensembles] = angle_asc_node_to_sat_ensemble[isat, i, :] - angle_asc_node_to_sat_ensemble[isat-1, i, isc1] # !!!! WORKS ONLY IF ONE SATELLITE FOR NOW
# #                                print max(angular_distance_between_two_clusters[i/index_dt_angular_distance_between_two_clusters, isc1*nb_ensembles:(isc1+1)*nb_ensembles]), i/index_dt_angular_distance_between_two_clusters, isc1*nb_ensembles, (isc1+1)*nb_ensembles
#                             save_index_dt_angular_distance_between_two_clusters.append(i)
                    # end of UNCOMMENT IF WANT TO COMPUTE angular_distance_between_two_clusters

                    # OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD
                    # angle_asc_node_to_sat_ensemble_ascending_order = np.zeros([nb_ensembles])
                    # angle_asc_node_to_sat_ensemble_ascending_order = np.sort(angle_asc_node_to_sat_ensemble[isat, i,:])
#                    sma_median_ensemble = np.median(sma_ensemble[i, :])
                    # for j in range(nb_ensembles-1):
                    #     angular_spacing_between_ensemble_sat_temp[i, j] = angle_asc_node_to_sat_ensemble_ascending_order[j+1] - angle_asc_node_to_sat_ensemble_ascending_order[j]
                    # angular_spacing_between_ensemble_sat_temp_sorted = np.zeros([nb_ensembles-1])
                    # angular_spacing_between_ensemble_sat_temp_sorted = np.sort(angular_spacing_between_ensemble_sat_temp[i, :])
                    # angular_spacing_between_ensemble_sat[i,:] = angular_spacing_between_ensemble_sat_temp_sorted[0:nb_ensembles-2]
                    # # for j in range(nb_ensembles-2):
                    # #     angular_spacing_between_ensemble_sat_converted_in_a_distance[i,j] = ( angular_spacing_between_ensemble_sat[i, j] * np.pi / 180. ) * sma_median_ensemble * 1000.
                    # end of OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD
                pickle_list.append(angle_asc_node_to_sat_ensemble); pickle_list_name.append('angle_asc_node_to_sat_ensemble')
                save_pickle_name = save_pickle_name + 'angle_asc_node_to_sat_'

                # UNCOMMENT IF WANT TO COMPUTE angular_distance_between_two_clusters
                # if isat == 1: # !!!! WORKS ONLY IF ONE SATELLITE FOR NOW
                #     pickle_list.append(angular_distance_between_two_clusters); pickle_list_name.append('angular_distance_between_two_clusters')
                #     pickle_list.append(save_index_dt_angular_distance_between_two_clusters); pickle_list_name.append('save_index_dt_angular_distance_between_two_clusters')
                    # save_pickle_name = save_pickle_name + 'angular_distance_between_two_clusters_'
                # end of  UNCOMMENT IF WANT TO COMPUTE angular_distance_between_two_clusters

                # pickle_list.append(angular_spacing_between_ensemble_sat); pickle_list_name.append('angular_spacing_between_ensemble_sat')
                # save_pickle_name = save_pickle_name + 'angular_spacing_between_ensemble_sat_'
                # pickle_list.append(angular_spacing_between_ensemble_sat_converted_in_a_distance); pickle_list_name.append('angular_spacing_between_ensemble_sat_converted_in_a_distance')
                # save_pickle_name = save_pickle_name + 'angular_spacing_between_ensemble_sat_converted_in_a_distance_'

        ## SAVE THE RESULTS IN A PICKLE
        save_pickle_name = save_pickle_name[:-1] +  '.pickle'
        if len(pickle_list_name) > 1:
            print "Saving the results in a pickle: " + save_pickle_name[:-1].split('/')[-1], rank
            with open(save_pickle_name, 'w') as f:
                pickle.dump([pickle_list, pickle_list_name], f)
                # print save_pickle_name, rank
                # print 

#     # else:
#     #     ## RESTORE THE RESULTS FROM THE PICKLE
#     #     with open(save_pickle_name) as f:
#     #         pickle_list, pickle_list_name = pickle.loadf)
#     #     if ( 'x_eci_ensemble' in pickle_list_name ): 
#     #         x_eci_ensemble = pickle_list[ pickle_list_name.index('x_eci_ensemble') ]
#     #     if ( 'y_eci_ensemble' in pickle_list_name ): 
#     #         y_eci_ensemble = pickle_list[ pickle_list_name.index('y_eci_ensemble') ]
#     #     if ( 'z_eci_ensemble' in pickle_list_name ): 
#     #         z_eci_ensemble = pickle_list[ pickle_list_name.index('z_eci_ensemble') ]
#     #     if ( 'pitch_ensemble' in pickle_list_name ): 
#     #         pitch_ensemble = pickle_list[ pickle_list_name.index('pitch_ensemble') ]
#     #     if ( 'roll_ensemble' in pickle_list_name ): 
#     #         roll_ensemble = pickle_list[ pickle_list_name.index('roll_ensemble') ]
#     #     if ( 'yaw_ensemble' in pickle_list_name ): 
#     #         yaw_ensemble = pickle_list[ pickle_list_name.index('yaw_ensemble') ]
#     #     if ( 'x_eci_main_sc' in pickle_list_name ):
#     #         x_eci_main_sc = pickle_list[ pickle_list_name.index('x_eci_main_sc') ]
#     #     if ( 'y_eci_main_sc' in pickle_list_name ):
#     #         y_eci_main_sc = pickle_list[ pickle_list_name.index('y_eci_main_sc') ]
#     #     if ( 'z_eci_main_sc' in pickle_list_name ):
#     #         z_eci_main_sc = pickle_list[ pickle_list_name.index('z_eci_main_sc') ]
#     #     if ( 'vx_eci_main_sc' in pickle_list_name ):
#     #         vx_eci_main_sc = pickle_list[ pickle_list_name.index('vx_eci_main_sc') ]
#     #     if ( 'vy_eci_main_sc' in pickle_list_name ):
#     #         vy_eci_main_sc = pickle_list[ pickle_list_name.index('vy_eci_main_sc') ]
#     #     if ( 'vz_eci_main_sc' in pickle_list_name ):
#     #         vz_eci_main_sc = pickle_list[ pickle_list_name.index('vz_eci_main_sc') ]
#     #     if ( 'algebric_distance_ensemble_main_sc_lvlh_x' in pickle_list_name ):
#     #         algebric_distance_ensemble_main_sc_lvlh_x = pickle_list[ pickle_list_name.index('algebric_distance_ensemble_main_sc_lvlh_x') ]
#     #     if ( 'algebric_distance_ensemble_main_sc_lvlh_y' in pickle_list_name ):
#     #         algebric_distance_ensemble_main_sc_lvlh_y = pickle_list[ pickle_list_name.index('algebric_distance_ensemble_main_sc_lvlh_y') ]
#     #     if ( 'algebric_distance_ensemble_main_sc_lvlh_z' in pickle_list_name ):
#     #         algebric_distance_ensemble_main_sc_lvlh_z = pickle_list[ pickle_list_name.index('algebric_distance_ensemble_main_sc_lvlh_z') ]
#     #     if ( 'sma_ensemble' in pickle_list_name ):
#     #         sma_ensemble = pickle_list[ pickle_list_name.index('sma_ensemble') ]
#     #     if ( 'inclination_ensemble' in pickle_list_name ):
#     #         inclination_ensemble = pickle_list[ pickle_list_name.index('inclination_ensemble') ]
#     #     if ( 'eccentricity_ensemble' in pickle_list_name ):
#     #         eccentricity_ensemble = pickle_list[ pickle_list_name.index('eccentricity_ensemble') ]
#     #     if ( 'true_anomaly_ensemble' in pickle_list_name ):
#     #         true_anomaly_ensemble = pickle_list[ pickle_list_name.index('true_anomaly_ensemble') ]
#     #     if ( 'RAAN_ensemble' in pickle_list_name ):
#     #         RAAN_ensemble = pickle_list[ pickle_list_name.index('RAAN_ensemble') ]
#     #     if ( 'argument_perigee_ensemble' in pickle_list_name ):
#     #         argument_perigee_ensemble = pickle_list[ pickle_list_name.index('argument_perigee_ensemble') ]
#     #     if ( 'rho_ensemble' in pickle_list_name ):
#     #         rho_ensemble = pickle_list[ pickle_list_name.index('rho_ensemble') ]    
#     #     if ( 'rho_ensemble_orbit_average' in pickle_list_name ):
#     #         rho_ensemble_orbit_average = pickle_list[ pickle_list_name.index('rho_ensemble_orbit_average') ]    
#     #     if ( 'time_orbit_average' in pickle_list_name ):
#     #         time_orbit_average = pickle_list[ pickle_list_name.index('time_orbit_average') ]    
#     #     if ( 'index_time_orbit_average' in pickle_list_name ):
#     #         index_time_orbit_average = pickle_list[ pickle_list_name.index('index_time_orbit_average') ]    
#     #     if ( 'f107_ensemble' in pickle_list_name ):
#     #         f107_ensemble = pickle_list[ pickle_list_name.index('f107_ensemble') ]
#     #     if ( 'f107a_ensemble' in pickle_list_name ):
#     #         f107a_ensemble = pickle_list[ pickle_list_name.index('f107a_ensemble') ]
#     #     if ( 'ap_ensemble' in pickle_list_name ):
#     #         ap_ensemble = pickle_list[ pickle_list_name.index('ap_ensemble') ]
#     #     if ( 'angle_asc_node_to_sat_ensemble' in pickle_list_name ):
#     #         angle_asc_node_to_sat_ensemble = pickle_list[ pickle_list_name.index('angle_asc_node_to_sat_ensemble') ]
#     #     if ( 'angular_spacing_between_ensemble_sat' in pickle_list_name ):
#     #         angular_spacing_between_ensemble_sat = pickle_list[ pickle_list_name.index('angular_spacing_between_ensemble_sat') ]
#     #     if ( 'angular_spacing_between_ensemble_sat_converted_in_a_distance' in pickle_list_name ):
#     #         angular_spacing_between_ensemble_sat_converted_in_a_distance = pickle_list[ pickle_list_name.index('angular_spacing_between_ensemble_sat_converted_in_a_distance') ]

# raise Exception
# r_eci_ensemble = np.zeros([nb_steps, nb_ensembles])
# for i in range(nb_steps):
#     for  j in range(nb_ensembles):
#         r_eci_ensemble[i,j] = np.sqrt( x_eci_ensemble[i,j]**2 + y_eci_ensemble[i,j]**2 + z_eci_ensemble[i,j]**2 )


# raise Exception
# ####### ##### FOR CADRE, WE WANT TO CALCULATE THE POSITION GIVEN BY TLE AND COMPARE IT TO THE MAIN SC
# i = -1 # we look at the last position of the main sc
# position_main_sc_eci = np.array([x_eci_main_sc[i], y_eci_main_sc[i], z_eci_main_sc[i]])
# position_tle_sc_eci = cadre_read_last_tle(get_prop_dir(2) + 'input/main_input/cadre_last_tle.txt')
# delta_vector_eci = position_main_sc_eci - position_tle_sc_eci
# velocity_main_sc_eci = np.array([vx_eci_main_sc[i], vy_eci_main_sc[i], vz_eci_main_sc[i]])
# delta_vector_lvlh_of_main_sc = eci_to_lvlh(position_main_sc_eci, velocity_main_sc_eci, delta_vector_eci)
# algebric_distance_tle_main_sc_lvlh_x= delta_vector_lvlh_of_main_sc[0] #* 1000 # * 1000 to convert the distance from km to m 
# algebric_distance_tle_main_sc_lvlh_y= delta_vector_lvlh_of_main_sc[1] #* 1000 # * 1000 to convert the distance from km to m 
# algebric_distance_tle_main_sc_lvlh_z= delta_vector_lvlh_of_main_sc[2] #* 1000 # * 1000 to convert the distance from km to m 

# raise Exception
# ########################################################################################################################################################################################################################################################################################################
# # PLOTS ################################################################################################################################################################################################################################################################################################
# ########################################################################################################################################################################################################################################################################################################
    
# ########################################################################################################################################################################################################################################################################################################
# if ( 'angle_asc_node_to_sat_ensemble' in pickle_list_name ):
#     # Plot the standard deviation of the angular_spacing_between_ensemble_sat_converted_in_a_distance distributions as a function of time
#     step_std = 24 # step in hours to calculate the standard deviation
#     step_std_in_index = step_std * 3600. / dt
#     std_daily = np.zeros([nb_steps])
#     med_daily = np.zeros([nb_steps])
#     mad_daily = np.zeros([nb_steps])
#     iqr_daily = np.zeros([nb_steps])
#     for i in range(nb_steps):
#         index_when_std = i * step_std_in_index
#         std_daily[i] = np.std(angle_asc_node_to_sat_ensemble[index_when_std, :]) # NOT GOOD TO USE THE STD
#         med_daily[i] = np.median(angular_spacing_between_ensemble_sat_temp[index_when_std,:]) # in meters
#         iqr_daily[i] = np.subtract(*np.percentile(angle_asc_node_to_sat_ensemble[index_when_std, :], [75, 25])) 
#         mad_daily[i] = np.median( np.abs( angle_asc_node_to_sat_ensemble[index_when_std, :] - np.median(angle_asc_node_to_sat_ensemble[index_when_std, :]) ) ) 

#     ## Plot
#     fontsize_plot = 14
#     height_fig = 5#12 
#     ratio_fig_size = 4./3
#     fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
#     fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
#     plt.rc('font', weight='bold') ## make the labels of the ticks in bold
#     gs = gridspec.GridSpec(1, 1)
#     gs.update(left= 0.09, right=0.98, top = 0.93,bottom = 0.08)
#     ax1 = fig.add_subplot(gs[0, 0])
#     ax1.set_title('Standard deviation of the along-track distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
#     ax1.set_ylabel('Standard Deviation (m)', fontsize = fontsize_plot, weight = 'bold')
#     ax1.set_xlabel('Time (days)', fontsize = fontsize_plot, weight = 'bold')
#     ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
#     [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
#     ax1.plot(np.arange(0,nb_steps), iqr_daily[0:nb_steps]*110000, 'k', linewidth = 2)

#     ax1.get_xaxis().tick_bottom()
#     ax1.get_yaxis().tick_left()
#     ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

#     fig_save_name = 'std_daily'
#     fig_save_name = root_save_fig_name + fig_save_name + '.png'
#     fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# ########################################################################################################################################################################################################################################################################################################
# if ( 'angle_asc_node_to_sat_ensemble' in pickle_list_name ):
#     # Plot the distribution of angular_spacing_between_ensemble_sat_converted_in_a_distance at different times (set by when_plot_in_hour)
#     when_plot_in_hour = 7 * 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
#     index_when_plot = (int) (when_plot_in_hour * 3600L / dt)

#     fontsize_plot = 14
#     height_fig = 12
#     ratio_fig_size = 4./3
#     fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
#     fig.suptitle('Spacecraft distribution along the orbit ' + '{0:.0f}'.format(when_plot_in_hour / 24.) + ' days after deployment', y = 0.96,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
#     plt.rc('font', weight='bold') ## make the labels of the ticks in bold
#     gs = gridspec.GridSpec(4, 1)
#     gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)
#     bin_width = 10 # in m

#     med = np.median(angle_asc_node_to_sat_ensemble[index_when_plot, :])
#     mad = np.median( np.abs( angle_asc_node_to_sat_ensemble[index_when_plot, :] - np.median(angle_asc_node_to_sat_ensemble[index_when_plot, :]) ) )  
#     quartile10 = np.percentile(angle_asc_node_to_sat_ensemble[index_when_plot, :], 10) 
#     quartile25 = np.percentile(angle_asc_node_to_sat_ensemble[index_when_plot, :], 25) 
#     quartile75 = np.percentile(angle_asc_node_to_sat_ensemble[index_when_plot, :], 75) 
#     quartile90 = np.percentile(angle_asc_node_to_sat_ensemble[index_when_plot, :], 90) 
#     iqr_daily = quartile75 - quartile25
#     quartiles_1090_daily = quartile90 - quartile10

#     ax2 = fig.add_subplot(gs[:3, 0])
#     ax2.set_title('', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
#     ax2.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+ ')', fontsize = fontsize_plot, weight = 'bold')
#     ax2.set_xlabel('', fontsize = fontsize_plot, weight = 'bold')
#     ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
#     [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
#     bin_width = iqr_daily/10.
#     hist_data = angle_asc_node_to_sat_ensemble[index_when_plot, :]
#     n, bins, patches = ax2.hist(hist_data , bins = np.arange(min(hist_data), max(hist_data) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'r',label = '')  
#     plt.vlines(quartile25, min(n), max(n),linewidth = 2); plt.vlines(quartile75, min(n),max(n), linewidth = 2) 
#     plt.vlines(quartile10, min(n), max(n), linewidth = 4, linestyle = 'dotted'); plt.vlines(quartile90, min(n), max(n), linewidth = 4, linestyle = 'dotted') 
#     ax2.plot(med, min(n), 'k',marker = '.', markersize = 20)
#     ax2.get_xaxis().get_major_formatter().set_useOffset(False)
#     ax2.set_xlim([med - 10*iqr_daily, med+10*iqr_daily])
#     ax2.get_xaxis().set_ticklabels([])
#     ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

#     ax1 = fig.add_subplot(gs[3, 0])
#     ymin = -0.2; ymax = 0.3
#     length_vert_iqr_daily = np.abs(ymin) / 3.
#     ax1.set_title('', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
#     ax1.set_ylabel('', fontsize = fontsize_plot, weight = 'bold')
#     ax1.set_xlabel('Angular distance from the AN (' + u'\N{DEGREE SIGN}'+ ')', fontsize = fontsize_plot, weight = 'bold')
#     ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
#     [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

#     ax1.hlines(0, med - 10 * iqr_daily, med + 10 * iqr_daily, linewidth = 2, color = 'b');
#     ax1.text( med + 10 * iqr_daily, 0 - length_vert_iqr_daily/1.5,'Orbit', horizontalalignment = 'right', fontsize = 14, color = 'b' )

#     ax1.plot(angle_asc_node_to_sat_ensemble[index_when_plot, :], np.zeros([nb_ensembles]), 'r.', markersize = 10)
#     ax1.vlines(quartile25, 0,length_vert_iqr_daily, linewidth = 2); ax1.vlines(quartile75, 0,length_vert_iqr_daily, linewidth = 2) 
#     ax1.text((quartile75 + quartile25 ) /2., length_vert_iqr_daily+length_vert_iqr_daily/5, '{0:.2f}'.format(iqr_daily * 110) + ' km (50% of sc)', horizontalalignment = 'center', fontsize = 14 )
#     ax1.annotate ('', (quartile25, length_vert_iqr_daily), (quartile75, length_vert_iqr_daily), arrowprops={'arrowstyle':'<->'})
#     ax1.vlines(quartile10, -length_vert_iqr_daily ,0 , linewidth = 4, linestyle = 'dotted'); ax1.vlines(quartile90, -length_vert_iqr_daily ,0, linewidth = 4, linestyle = 'dotted') 
#     ax1.annotate ('', (quartile10, -length_vert_iqr_daily), (quartile90, -length_vert_iqr_daily), arrowprops={'arrowstyle':'<->'})
#     ax1.text((quartile10 + quartile90 ) /2., -length_vert_iqr_daily-length_vert_iqr_daily/1.5, '{0:.2f}'.format(quartiles_1090_daily * 110) + ' km (80% of sc)', horizontalalignment = 'center', fontsize = 14, verticalalignment = 'bottom' )
#     ax1.plot(med, 0, 'k',marker = '.', markersize = 15)

#     plt.legend(loc = 2, scatterpoints=1)
#     ax1.get_xaxis().tick_bottom()
#     ax1.set_ylim([ymin, ymax])
# #    ax1.get_yaxis().tick_left()
#     ax1.yaxis.set_visible(False)
#     ax1.spines['left'].set_visible(False); ax1.spines['right'].set_visible(False); ax1.spines['top'].set_visible(False)
#     ax1.get_xaxis().get_major_formatter().set_useOffset(False)
#     ax1.set_xlim([med - 10*iqr_daily, med+10*iqr_daily])
# #    ax1.set_xlim([min(bins), 0])
#     ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

#     fig_save_name = 'sc_distribution_' + '{0:.0f}'.format(when_plot_in_hour / 24.) + '_days_after_deployment'
#     fig_save_name = root_save_fig_name + fig_save_name + '.eps'
#     fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
#     os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission)






# ########################################################################################################################################################################################################################################################################################################
# # Distance between 2 satellites 
# sat1 = np.where(angle_asc_node_to_sat_ensemble[90,:] == min(angle_asc_node_to_sat_ensemble[90,:]))[0]
# sat2 = np.where(angle_asc_node_to_sat_ensemble[90,:] == max(angle_asc_node_to_sat_ensemble[90,:]))[0]

# yaw_vel_sat1 = ( yaw_ensemble[90, sat1] - yaw_ensemble[90-1, sat1] ) / (24. * 3600)
# pitch_vel_sat1 = ( pitch_ensemble[90, sat1] - pitch_ensemble[90-1, sat1] ) / (24. * 3600)
# roll_vel_sat1 = ( roll_ensemble[90, sat1] - roll_ensemble[90-1, sat1] ) / (24. * 3600)

# yaw_vel_sat2 = ( yaw_ensemble[90, sat2] - yaw_ensemble[90-1, sat2] ) / (24. * 3600)
# pitch_vel_sat2 = ( pitch_ensemble[90, sat2] - pitch_ensemble[90-1, sat2] ) / (24. * 3600)
# roll_vel_sat2 = ( roll_ensemble[90, sat2] - roll_ensemble[90-1, sat2] ) / (24. * 3600)

# dist_sat1_sat2 = np.sqrt( ( x_eci_ensemble[:, sat1] - x_eci_ensemble[:, sat2] )**2 +  ( y_eci_ensemble[:, sat1] - y_eci_ensemble[:, sat2] )**2 +  ( z_eci_ensemble[:, sat1] - z_eci_ensemble[:, sat2] )**2 )

# print pitch_vel_sat1, roll_vel_sat1, yaw_vel_sat1
# print pitch_vel_sat2, roll_vel_sat2, yaw_vel_sat2

# # 200/80
# # In [86]: print pitch_vel_sat1, roll_vel_sat1, yaw_vel_sat1
# # -0.118257294919 -2.12316849306 -0.921393159606
# # In [87]: print pitch_vel_sat2, roll_vel_sat2, yaw_vel_sat2
# # 0.182033980787 0.677633533565 -2.12392994444


# fontsize_plot = 16
# height_fig = 12 
# ratio_fig_size = 4./3
# fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
# fig_title = 'Distance between sat1 and sat2 after deployment over 90 days'
# fig.suptitle(fig_title, y = 0.953,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
# plt.rc('font', weight='bold') ## make the labels of the ticks in bold
# gs = gridspec.GridSpec(1, 1)
# gs.update(left= 0.08, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

# ax = fig.add_subplot(gs[0, 0])
# ax.plot(dist_sat1_sat2, linewidth = 2, color = 'k')
# ax.set_ylabel('Distance (km)', weight = 'bold', fontsize  = fontsize_plot)
# ax.set_xlabel('Time after deployment (days)', weight = 'bold', fontsize  = fontsize_plot)

# [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
# ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
# plt.rc('font', weight='bold') ## make the labels of the ticks in bold
# fig_save_name = fig_title.replace(" ","_").lower()
# fig_save_name = root_save_fig_name + fig_save_name + '.png'
# fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
# os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission)


# ########################################################################################################################################################################################################################################################################################################
# if ( 'angle_asc_node_to_sat_ensemble' in pickle_list_name ):
#     # Plot the distribution of angular_spacing_between_ensemble_sat_converted_in_a_distance at different times (set by when_plot_in_hour)
#     fontsize_plot = 14
#     height_fig = 5 
#     ratio_fig_size = 4./3
#     fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
#     fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
#     plt.rc('font', weight='bold') ## make the labels of the ticks in bold
#     gs = gridspec.GridSpec(1, 1)
#     gs.update(left= 0.08, right=0.98, top = 0.93,bottom = 0.08)
#     bin_width = 10 # in m

#     when_plot_in_hour = 24 * 50. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
#     index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
#     ax1 = fig.add_subplot(gs[0, 0])
#     ax1.set_title('Distribution along-track after N hours', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
#     ax1.set_ylabel('Distribution (# out of ' + str(nb_ensembles) + ')', fontsize = fontsize_plot, weight = 'bold')
#     ax1.set_xlabel('Algebraic LVLH_X distance from reference satellite (km)', fontsize = fontsize_plot, weight = 'bold')
#     ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
#     [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
#     n_max, bins, patches = ax1.hist(angular_spacing_between_ensemble_sat_converted_in_a_distance[index_when_plot,:], bins = np.arange(min(angular_spacing_between_ensemble_sat_converted_in_a_distance[index_when_plot,:]), max(angular_spacing_between_ensemble_sat_converted_in_a_distance[index_when_plot,:]) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'k',label = '')  

#     plt.legend(loc = 2, scatterpoints=1)
#     ax1.get_xaxis().tick_bottom()
#     ax1.get_yaxis().tick_left()
# #    ax1.set_ylim([0, 100])
# #    ax1.set_xlim([min(bins), 0])
#     ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)


#     fig_save_name = 'angular_spacing_between_ensemble_sat_converted_in_a_distance'
#     fig_save_name = root_save_fig_name + fig_save_name + '.png'
#     fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    




# ########################################################################################################################################################################################################################################################################################################
# if ( ( 'x_eci' in ensemble_to_plot ) & ( 'y_eci' in ensemble_to_plot ) & ( 'z_eci' in ensemble_to_plot ) ): 
#     # Plot the algebric lvlh_x distance at a given time
#     ## Distance along the x axis
#     fontsize_plot = 14
#     height_fig = 5 
#     ratio_fig_size = 4./3
#     fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
#     fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
#     plt.rc('font', weight='bold') ## make the labels of the ticks in bold
#     gs = gridspec.GridSpec(1, 1)
#     gs.update(left= 0.08, right=0.98, top = 0.93,bottom = 0.08)
#     bin_width = 5. # in km

# #    when_plot_in_hour = 96. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
#     index_when_plot = algebric_distance_ensemble_main_sc_lvlh_x.shape[0] - 1 #(int) (when_plot_in_hour * 3600L / dt)
#     median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:])#/1000.)
#     tenth_algebric_ditance = np.percentile(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:])#/1000., 10)
#     std = np.std(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:])#/1000.)    
#     ax1 = fig.add_subplot(gs[0, 0])
#     n, bins, patches = ax1.hist(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:], np.arange(min(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]), max(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'g',label = 'N = ' + '{0:.0f}' .format( ( input_variables[1] - input_variables[0] ).days*24 + ( input_variables[1] - input_variables[0] ).seconds / 3600 ) + 'h - ' + 'median = ' + '{0:.3f}' .format(median_algebric_ditance) + ' km, ' + 'std = ' + '{0:.1e}' .format(std) + ' km') 
#     ax1.plot([median_algebric_ditance, median_algebric_ditance], [min(n), max(n)], 'g', linewidth = 2, linestyle = 'dotted')
#     ## ADD THE TLE POSITION
#     ax1.scatter([algebric_distance_tle_main_sc_lvlh_x],[1.5], s = 200, c = 'g', marker = '*', edgecolor = 'k', label = 'CADRE')


#     plt.legend(loc = 2, scatterpoints=1)
#     ax1.get_xaxis().tick_bottom()
#     ax1.get_yaxis().tick_left()
#     ax1.set_ylim([0, 100])
# #    ax1.set_xlim([min(bins), 0])
#     ax1.margins(0,1) # autoscale both axes(fist value is for the x axis, second value for the y axis)


#     fig_save_name = 'distribution_along_track'
#     fig_save_name = root_save_fig_name + fig_save_name + '.png'
#     fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    
#     # when_plot_in_hour = 12. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
#     # index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
#     # median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.)
#     # tenth_algebric_ditance = np.percentile(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000., 10)
#     # std = np.std(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.)    
#     # ax1 = fig.add_subplot(gs[0, 0])
#     # ax1.set_title('Distribution along-track after N hours', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
#     # ax1.set_ylabel('Distribution (# out of ' + str(nb_ensembles) + ')', fontsize = fontsize_plot, weight = 'bold')
#     # ax1.set_xlabel('Algebraic LVLH_X distance from reference satellite (km)', fontsize = fontsize_plot, weight = 'bold')
#     # ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
#     # [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
#     # n_max, bins, patches = ax1.hist(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000., bins = np.arange(min(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.), max(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'k',label = 'N = ' + '{0:.0f}' .format( ( index_when_plot * dt / 3600. )) + 'h - ' + 'median = ' + '{0:.3f}' .format(median_algebric_ditance) + ' km, ' + 'std = ' + '{0:.1e}' .format(std) + ' km')  
#     # ax1.plot([median_algebric_ditance, median_algebric_ditance], [min(n_max), max(n_max)], 'k', linewidth = 2, linestyle = 'dotted')

#     # when_plot_in_hour = 48. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
#     # index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
#     # median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.)
#     # tenth_algebric_ditance = np.percentile(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000., 10)
#     # std = np.std(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.)    
#     # ax1 = fig.add_subplot(gs[0, 0])
#     # n, bins, patches = ax1.hist(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000., np.arange(min(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.), max(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'b',label = 'N = ' + '{0:.0f}' .format( ( index_when_plot * dt / 3600. )) + 'h - ' + 'median = ' + '{0:.3f}' .format(median_algebric_ditance) + ' km, ' + 'std = ' + '{0:.1e}' .format(std) + ' km') 
#     # ax1.plot([median_algebric_ditance, median_algebric_ditance], [min(n), max(n)], 'b', linewidth = 2, linestyle = 'dotted')

#     # when_plot_in_hour = 96. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
#     # index_when_plot = (int) (when_plot_in_hour * 3600L / dt) 
#     # median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.)
#     # tenth_algebric_ditance = np.percentile(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000., 10)
#     # std = np.std(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.)    
#     # ax1 = fig.add_subplot(gs[0, 0])
#     # n, bins, patches = ax1.hist(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000., np.arange(min(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.), max(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]/1000.) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'r',label = 'N = ' + '{0:.0f}' .format( ( index_when_plot * dt / 3600. )) + 'h - ' + 'median = ' + '{0:.3f}' .format(median_algebric_ditance) + ' km, ' + 'std = ' + '{0:.1e}' .format(std) + ' km') 
#     # ax1.plot([median_algebric_ditance, median_algebric_ditance], [min(n), max(n)], 'r', linewidth = 2, linestyle = 'dotted')


# ########################################################################################################################################################################################################################################################################################################
# if ( ( 'x_eci' in ensemble_to_plot ) & ( 'y_eci' in ensemble_to_plot ) & ( 'z_eci' in ensemble_to_plot ) ): 
#     # Plot the standard deviation of the along-track distributions as a function of time
#     ## Standard deviation of the along-track distribution as a function of time (by step of N hours (called step_std))
#     step_std = 24 # step in hours to calculate the standard deviation
#     step_std_in_index = step_std * 3600. / dt
#     std_daily = np.zeros([nb_steps])
#     med_daily_z = np.zeros([nb_steps])
#     med_daily_x = np.zeros([nb_steps])
#     min_r_ensemble_daily = np.zeros([nb_steps])
#     for i in range(nb_steps):
#         index_when_std = i * step_std_in_index
#         std_daily[i] = np.std(algebric_distance_ensemble_main_sc_lvlh_x[index_when_std,:]) # in meters
#         med_daily_z[i] = np.median(algebric_distance_ensemble_main_sc_lvlh_z[index_when_std,:]) # in meters
#         med_daily_x[i] = np.median(algebric_distance_ensemble_main_sc_lvlh_x[index_when_std,:]) # in meters
#         min_r_ensemble_daily[i] = np.median( r_eci_ensemble[index_when_std, :] )

#     ## Plot
#     fontsize_plot = 14
#     height_fig = 8 
#     ratio_fig_size = 4./3
#     fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
#     fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
#     plt.rc('font', weight='bold') ## make the labels of the ticks in bold
#     gs = gridspec.GridSpec(1, 1)
#     gs.update(left= 0.09, right=0.98, top = 0.93,bottom = 0.08)
#     ax1 = fig.add_subplot(gs[0, 0])
#     ax1.set_title('Standard deviation of the along-track distributions VS time', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
#     ax1.set_ylabel('Standard Deviation (m)', fontsize = fontsize_plot, weight = 'bold')
#     ax1.set_xlabel('Time (days)', fontsize = fontsize_plot, weight = 'bold')
#     ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
#     [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
#     ax1.plot(np.arange(0,nb_steps), med_daily_x[0:nb_steps]/1000., 'k', linewidth = 2)

#     ax1.get_xaxis().tick_bottom()
#     ax1.get_yaxis().tick_left()
#     ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

#     fig_save_name = 'std_daily'
#     fig_save_name = root_save_fig_name + fig_save_name + '.png'
#     fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


#     ########################################################################################################################################################################################################################################################################################################
# if ( ( 'pitch' in ensemble_to_plot ) & ( 'roll' in ensemble_to_plot ) & ( 'yaw' in ensemble_to_plot ) ): 
#     # Plot the angular velocities distribution
#     fontsize_plot = 14
#     when_plot_in_hour = 10. * 24 #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
#     index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
#     median_algebric_ditance = np.median(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))
#     tenth_algebric_ditance = np.percentile(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 10)
#     std = np.std(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))    
#     fig = plt.figure(num=None, figsize=(13, 9), dpi=80, facecolor='w', edgecolor='k')
#     fig.suptitle('Angular velocity distributions (pitch, roll, yaw)', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
#     plt.rc('font', weight='bold') ## make the labels of the ticks in bold
#     gs = gridspec.GridSpec(1, 3)
#     gs.update(left= 0.05, right=0.99, top = 0.90,bottom = 0.09)

#     ## Pitch angular velocity distribution
#     ax1 = fig.add_subplot(gs[0, 0])
#     ax1.set_title('Pitch ang. vel. dist.', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
#     #Algebric LVLH_X distance after '+str( ( index_when_plot * dt / 3600. ))+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
#     ax1.set_ylabel('Distribution (# out of ' + str(nb_ensembles) + ')', fontsize = fontsize_plot, weight = 'bold')
#     ax1.set_xlabel('Pitch angular velocity' + u' (\N{DEGREE SIGN}/s)', fontsize = fontsize_plot, weight = 'bold')
#     #ax1.xaxis.set_ticklabels([])
#     ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
#     [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
#     ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
#     n, bins, patches = ax1.hist(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 100,  histtype='stepfilled', alpha = 0.7) 
#     ax1.plot([median_algebric_ditance/1000., median_algebric_ditance/1000.], [min(n), max(n)], 'k', linewidth = 2, linestyle = 'dotted')
# #    ax1.text(ax1.get_xlim()[1], ax1.get_ylim()[1], 'std = ' + '{0:.2f}' .format(std) + u' (\N{DEGREE SIGN}/s)\n' , horizontalalignment = 'right', verticalalignment = 'top', fontsize = fontsize_plot, weight = 'bold')
#     ax1.get_xaxis().tick_bottom()
#     ax1.get_yaxis().tick_left()

#     ## Roll ang. vel. distribution
#     median_algebric_ditance = np.median(roll_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))
#     tenth_algebric_ditance = np.percentile(roll_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 10)
#     std = np.std(roll_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))    
#     ax2 = fig.add_subplot(gs[0, 1])
#     ax2.set_title('Roll ang. vel. dist.', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
#     #Algebric LVLH_X distance after '+str( ( index_when_plot * dt / 3600. ))+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
#     ax2.set_xlabel('Roll angular velocity' + u' (\N{DEGREE SIGN}/s)', fontsize = fontsize_plot, weight = 'bold')
#     #ax2.xaxis.set_ticklabels([])
#     ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
#     [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
#     ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
#     n, bins, patches = ax2.hist(roll_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 100,  histtype='stepfilled', alpha = 0.7) 
#     ax2.plot([median_algebric_ditance/1000., median_algebric_ditance/1000.], [min(n), max(n)], 'k', linewidth = 2, linestyle = 'dotted')
# #    ax2.text(ax2.get_xlim()[1], ax2.get_ylim()[1], 'std = ' + '{0:.2f}' .format(std) + u' (\N{DEGREE SIGN}/s)\n' , horizontalalignment = 'right', verticalalignment = 'top', fontsize = fontsize_plot, weight = 'bold')
#     ax2.get_xaxis().tick_bottom()
#     ax2.get_yaxis().tick_left()

#     ## Yaw ang. vel. distribution
#     median_algebric_ditance = np.median(yaw_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))
#     tenth_algebric_ditance = np.percentile(yaw_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 10)
#     std = np.std(yaw_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))    
#     ax3 = fig.add_subplot(gs[0, 2])
#     ax3.set_title('Yaw ang. vel. dist.', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
#     #Algebric LVLH_X distance after '+str( ( index_when_plot * dt / 3600. ))+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
#     ax3.set_xlabel('Yaw angular velocity' + u' (\N{DEGREE SIGN}/s)', fontsize = fontsize_plot, weight = 'bold')
#     #ax3.xaxis.set_ticklabels([])
#     ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
#     [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
#     ax3.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
#     n, bins, patches = ax3.hist(yaw_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 100,  histtype='stepfilled', alpha = 0.7) 
#     ax3.plot([median_algebric_ditance/1000., median_algebric_ditance/1000.], [min(n), max(n)], 'k', linewidth = 2, linestyle = 'dotted')
# #    ax3.text(ax3.get_xlim()[1], ax3.get_ylim()[1], 'std = ' + '{0:.2f}' .format(std) + u' (\N{DEGREE SIGN}/s)\n' , horizontalalignment = 'right', verticalalignment = 'top', fontsize = fontsize_plot, weight = 'bold')
#     ax3.get_xaxis().tick_bottom()
#     ax3.get_yaxis().tick_left()

#     fig_save_name = 'distribution_angular_velocity'
#     fig_save_name = root_save_fig_name + fig_save_name + '.eps'
#     fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

# ########################################################################################################################################################################################################################################################################################################
# if ( ( 'pitch' in ensemble_to_plot ) & ( 'roll' in ensemble_to_plot ) & ( 'yaw' in ensemble_to_plot ) ): 
#     # Examples of attitude: pitch/roll/yaw as a function of time for the first orbit
#     ## Pitch vs time
#     fontsize_plot = 14
#     when_plot_in_hour = 30./3600 #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
#     index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
#     median_algebric_ditance = np.median(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))
#     tenth_algebric_ditance = np.percentile(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600), 10)
#     std = np.std(pitch_ensemble[index_when_plot,:]/( ( index_when_plot * dt / 3600. )*3600))    
#     fig = plt.figure(num=None, figsize=(15, 10), dpi=80, facecolor='w', edgecolor='k')
#     fig.suptitle('Examples of the attitude of two satellites', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
#     plt.rc('font', weight='bold') ## make the labels of the ticks in bold
#     gs = gridspec.GridSpec(3, 1)
#     gs.update(left= 0.06, right=0.99, top = 0.93,bottom = 0.08)
#     ax1 = fig.add_subplot(gs[0, 0])
#     ax1.set_title('Pitch', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
#     #Algebric LVLH_X distance after '+str( ( index_when_plot * dt / 3600. ))+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
#     ax1.set_ylabel('Pitch' + u' (\N{DEGREE SIGN})', fontsize = fontsize_plot, weight = 'bold')
#     #ax1.xaxis.set_ticklabels([])
#     ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
#     [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
#     ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
#     x_axis = np.arange(0, 97, 0.5)
#     sat_nb = 100
#     ax1.plot(x_axis, np.mod( pitch_ensemble[0:14, sat_nb], 360 ), 'r', linewidth = 2, label = 's'+str(sat_nb))
#     sat_nb = 500
#     ax1.plot(x_axis, np.mod( pitch_ensemble[0:14, sat_nb], 360 ), 'b', linewidth = 2, label = 's'+str(sat_nb))
#     ax1.get_xaxis().tick_bottom()
#     ax1.get_yaxis().tick_left()

#     ## Roll vs time
#     ax2 = fig.add_subplot(gs[1, 0])
#     ax2.set_title('Roll', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
#     #Algebric LVLH_X distance after '+str( ( index_when_plot * dt / 3600. ))+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
#     ax2.set_ylabel('Roll' + u' (\N{DEGREE SIGN})', fontsize = fontsize_plot, weight = 'bold')
#     #ax2.xaxis.set_ticklabels([])
#     ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
#     [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
#     ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
#     sat_nb = 100
#     ax2.plot(x_axis, np.mod( roll_ensemble[0:14, sat_nb], 360 ), 'r', linewidth = 2, label = 's'+str(sat_nb))
#     sat_nb = 500
#     ax2.plot(x_axis, np.mod( roll_ensemble[0:14, sat_nb], 360 ), 'b', linewidth = 2, label = 's'+str(sat_nb))
#     ax2.get_xaxis().tick_bottom()
#     ax2.get_yaxis().tick_left()

#     ## Yaw vs time
#     ax3 = fig.add_subplot(gs[2, 0])
#     ax3.set_title('Yaw', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
#     #Algebric LVLH_X distance after '+str( ( index_when_plot * dt / 3600. ))+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
#     ax3.set_ylabel('Yaw' + u' (\N{DEGREE SIGN})', fontsize = fontsize_plot, weight = 'bold')
#     #ax3.xaxis.set_ticklabels([])
#     ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
#     [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
#     ax3.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
#     sat_nb = 100
#     ax3.plot(x_axis, np.mod( yaw_ensemble[0:14, sat_nb], 360 ), 'r', linewidth = 2, label = 's'+str(sat_nb))
#     sat_nb = 500
#     ax3.plot(x_axis, np.mod( yaw_ensemble[0:14, sat_nb], 360 ), 'b', linewidth = 2, label = 's'+str(sat_nb))
#     ax3.get_xaxis().tick_bottom()
#     ax3.get_yaxis().tick_left()
#     ax3.set_xlabel('Time in orbit (min)', fontsize = fontsize_plot, weight = 'bold')

#     fig_save_name = 'example_attitude_vs_time_ang_velo'
#     fig_save_name = root_save_fig_name + fig_save_name + '.png'
#     fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

#     raise Exception

# ########################################################################################################################################################################################################################################################################################################

# if ( ( 'x_eci' in ensemble_to_plot ) & ( 'y_eci' in ensemble_to_plot ) & ( 'z_eci' in ensemble_to_plot ) ): 
#     # Plot the algebric lvlh_z distance at a given time
#     ## Distance along the x axis
#     fontsize_plot = 14
#     height_fig = 10 
#     ratio_fig_size = 4./3
#     fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
#     fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
#     plt.rc('font', weight='bold') ## make the labels of the ticks in bold
#     gs = gridspec.GridSpec(1, 1)
#     gs.update(left= 0.06, right=0.98, top = 0.93,bottom = 0.08)
#     bin_width = 5/1000. # in km

#     when_plot_in_hour = 12. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
#     index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
#     median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.)
#     tenth_algebric_ditance = np.percentile(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000., 10)
#     std = np.std(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.)    
#     ax1 = fig.add_subplot(gs[0, 0])
#     ax1.set_title('Distribution along-track after N hours', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
#     ax1.set_ylabel('Distribution (# out of ' + str(nb_ensembles) + ')', fontsize = fontsize_plot, weight = 'bold')
#     ax1.set_xlabel('Algebraic distance from reference satellite (km)', fontsize = fontsize_plot, weight = 'bold')
#     ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
#     [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
#     n, bins, patches = ax1.hist(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000., bins = np.arange(min(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.), max(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'k',label = 'N = ' + '{0:.0f}' .format( ( index_when_plot * dt / 3600. )) + 'h - ' + 'median = ' + '{0:.3f}' .format(median_algebric_ditance) + ' km, ' + 'std = ' + '{0:.1e}' .format(std) + ' km')  
#     ax1.plot([median_algebric_ditance, median_algebric_ditance], [min(n), max(n)], 'k', linewidth = 2, linestyle = 'dotted')

#     when_plot_in_hour = 24. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
#     index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
#     median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.)
#     tenth_algebric_ditance = np.percentile(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000., 10)
#     std = np.std(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.)    
#     ax1 = fig.add_subplot(gs[0, 0])
#     n, bins, patches = ax1.hist(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000., np.arange(min(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.), max(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'b',label = 'N = ' + '{0:.0f}' .format( ( index_when_plot * dt / 3600. )) + 'h - ' + 'median = ' + '{0:.3f}' .format(median_algebric_ditance) + ' km, ' + 'std = ' + '{0:.1e}' .format(std) + ' km') 
#     ax1.plot([median_algebric_ditance, median_algebric_ditance], [min(n), max(n)], 'k', linewidth = 2, linestyle = 'dotted')

#     when_plot_in_hour = 48. #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
#     index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
#     median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.)
#     tenth_algebric_ditance = np.percentile(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000., 10)
#     std = np.std(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.)    
#     ax1 = fig.add_subplot(gs[0, 0])
#     n, bins, patches = ax1.hist(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000., np.arange(min(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.), max(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000.) + bin_width, bin_width),  histtype='stepfilled', alpha = 0.7, color = 'r',label = 'N = ' + '{0:.0f}' .format( ( index_when_plot * dt / 3600. )) + 'h - ' + 'median = ' + '{0:.3f}' .format(median_algebric_ditance) + ' km, ' + 'std = ' + '{0:.1e}' .format(std) + ' km') 
#     ax1.plot([median_algebric_ditance, median_algebric_ditance], [min(n), max(n)], 'k', linewidth = 2, linestyle = 'dotted')

#     plt.legend(loc = 2)
#     ax1.get_xaxis().tick_bottom()
#     ax1.get_yaxis().tick_left()
#     ax1.set_xlim([min(bins), 0])
#     ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

#     fig_save_name = 'distribution_radial'
#     fig_save_name = root_save_fig_name + fig_save_name + '.png'
#     fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# raise Exception
# ########################################################################################################################################################################################################################################################################################################
# ##  ## EVERYTHING BELOW IS OLD 
# ########################################################################################################################################################################################################################################################################################################

# raise Exception
# ########################################################################################################################################################################################################################################################################################################
# ## EVERYTHING BELOW IS OLD

# # 2D Histograms: distribution in the along-track and radial directions
# when_plot_in_hour = 12
# index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
# fig = plt.figure(num=None, figsize=(15, 7), dpi=80, facecolor='w', edgecolor='k')
# fig.suptitle('', fontsize = 22)
# plt.rc('font', weight='bold') ## make the labels of the ticks in bold
# ax1 = fig.add_subplot(111)
# ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
# #Algebric LVLH_X distance after '+str(when_plot_in_hour)+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
# ax1.set_ylabel('Radial (m)', fontsize = fontsize_plot, weight = 'bold')
# ax1.set_xlabel('Along-track (m)', fontsize = 18, weight = 'bold')
# #ax1.xaxis.set_ticklabels([])
# ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
# [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
# ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
# x_hist2d = algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]
# y_hist2d = algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]
# H, xedges, yedges, img = ax1.hist2d(x_hist2d, y_hist2d, bins = 100) #, range = [[-1000, 1000],[-1000, 1000]]
# plt.colorbar(img)
# #extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
# #im = ax1.imshow(H, cmap=plt.cm.jet, extent=extent)
# #fig.colorbar(im, ax=ax1)
# plt.show()
# #cbar = plt.colorbar()


# raise Exception

# ######################################################################################################################################################################################################################################################################################################## 
# ## EVERYTHING BELOW IS OLD

# #y_pos_distance = 5
# #plt.scatter(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:], np.zeros(nb_ensembles)-y_pos_distance, linewidth = 2)    

# # # Plot the magnitude of the algebric distance oalong the LVLH_X axis
# # ax2 = fig.add_subplot(122)
# # ax2.set_title('Absolute distance', weight = 'bold', fontsize = 20,  y = 1.008) 
# # n, bins, patches = plt.hist(np.abs(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]), 100,  histtype='stepfilled', alpha = 0.7) # 
# # ax2.set_xlabel('Distance (m)', fontsize = 18, weight = 'bold')
# # ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
# # [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
# # ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
# # median_absolute_ditance = np.median(np.abs(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]))
# # tenth_absolute_ditance = np.percentile(np.abs(algebric_distance_ensemble_main_sc_lvlh_x[index_when_plot,:]), 10)
# # ax2.plot([median_absolute_ditance, median_absolute_ditance], [min(n), max(n)], 'k', linewidth = 2)

# raise Exception
# ## Distance along the z axis
# when_plot_in_hour = 12 #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
# index_when_plot = -1 #(int) (when_plot_in_hour * 3600L / dt)
# median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:])
# tenth_algebric_ditance = np.percentile(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:], 10)
# fig = plt.figure(num=None, figsize=(20, 10), dpi=80, facecolor='w', edgecolor='k')
# fig.suptitle('Radial distribution', fontsize = 22)
# plt.rc('font', weight='bold') ## make the labels of the ticks in bold
# ax1 = fig.add_subplot(121)
# ax1.set_title('Algebric distance', weight = 'bold', fontsize = 20,  y = 1.008) 
# #Algebric LVLH_Z distance after '+str(when_plot_in_hour)+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
# ax1.set_ylabel('Distribution (# out of ' + str(nb_ensembles) + ')', fontsize = 18, weight = 'bold')
# ax1.set_xlabel('Distance (km)', fontsize = 18, weight = 'bold')
# #ax1.xaxis.set_ticklabels([])
# ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
# [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
# ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
# n, bins, patches = ax1.hist(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]/1000., 100,  histtype='stepfilled', alpha = 0.7) 
# ax1.plot([median_algebric_ditance/1000., median_algebric_ditance/1000.], [min(n), max(n)], 'k', linewidth = 2)
# #y_pos_distance = 5
# #plt.scatter(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:], np.zeros(nb_ensembles)-y_pos_distance, linewidth = 2)    

# # Plot the magnitude of the algebric distance oalong the LVLH_Z axis
# ax2 = fig.add_subplot(122)
# ax2.set_title('Absolute distance', weight = 'bold', fontsize = 20,  y = 1.008) 
# n, bins, patches = plt.hist(np.abs(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]), 100,  histtype='stepfilled', alpha = 0.7) # 
# ax2.set_xlabel('Distance (m)', fontsize = 18, weight = 'bold')
# ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
# [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
# ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
# median_absolute_ditance = np.median(np.abs(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]))
# tenth_absolute_ditance = np.percentile(np.abs(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:]), 10)
# ax2.plot([median_absolute_ditance, median_absolute_ditance], [min(n), max(n)], 'k', linewidth = 2)


# ## Distance along the y axis
# when_plot_in_hour = 12 #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
# index_when_plot = -1 #(int) (when_plot_in_hour * 3600L / dt)
# median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:])
# tenth_algebric_ditance = np.percentile(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:], 10)
# fig = plt.figure(num=None, figsize=(20, 10), dpi=80, facecolor='w', edgecolor='k')
# fig.suptitle('Cross-track distribution', fontsize = 22)
# plt.rc('font', weight='bold') ## make the labels of the ticks in bold
# ax1 = fig.add_subplot(121)
# ax1.set_title('Algebric distance', weight = 'bold', fontsize = 20,  y = 1.008) 
# #Algebric LVLH_Y distance after '+str(when_plot_in_hour)+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
# ax1.set_ylabel('Distribution (# out of ' + str(nb_ensembles) + ')', fontsize = 18, weight = 'bold')
# ax1.set_xlabel('Distance (km)', fontsize = 18, weight = 'bold')
# #ax1.xaxis.set_ticklabels([])
# ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
# [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
# ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
# n, bins, patches = ax1.hist(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:]/1000., 100,  histtype='stepfilled', alpha = 0.7) 
# ax1.plot([median_algebric_ditance/1000., median_algebric_ditance/1000.], [min(n), max(n)], 'k', linewidth = 2)
# #y_pos_distance = 5
# #plt.scatter(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:], np.zeros(nb_ensembles)-y_pos_distance, linewidth = 2)    

# # Plot the magnitude of the algebric distance oalong the LVLH_Y axis
# ax2 = fig.add_subplot(122)
# ax2.set_title('Absolute distance', weight = 'bold', fontsize = 20,  y = 1.008) 
# n, bins, patches = plt.hist(np.abs(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:]), 100,  histtype='stepfilled', alpha = 0.7) # 
# ax2.set_xlabel('Distance (m)', fontsize = 18, weight = 'bold')
# ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
# [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
# ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
# median_absolute_ditance = np.median(np.abs(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:]))
# tenth_absolute_ditance = np.percentile(np.abs(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:]), 10)
# ax2.plot([median_absolute_ditance, median_absolute_ditance], [min(n), max(n)], 'k', linewidth = 2)


# ########################################################################################################################################################################################################################################################################################################
# ## EVERYTHING BELOW IS OLD

# ## Distance along the y axis
# when_plot_in_hour = 12 #uncomment line below (get rid off the -1 if you want to use when_plot_in_hour)!!!!!!!!!
# index_when_plot = -1 #(int) (when_plot_in_hour * 3600L / dt)
# index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
# median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:])
# fig_algebric_distance_ensemble_main_sc_lvlh_y = plt.figure()
# plt.scatter(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:], np.zeros(nb_ensembles)-y_pos_distance, linewidth = 2)    
# plt.title('Algebric LVLH_Y distance after '+str(when_plot_in_hour)+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m')


# ## Histogram
# n, bins, patches = plt.hist(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:], 100,  histtype='stepfilled', alpha = 0.7) # 
# plt.xlim([min(bins), max(bins)])
# plt.ylim([-y_pos_distance, max(n)])
# plt.plot([median_algebric_ditance, median_algebric_ditance], [min(n), max(n)], 'k', linewidth = 2)

# plt.xlabel('Algebric distance (m)')
# print median_algebric_ditance
# print np.percentile(algebric_distance_ensemble_main_sc_lvlh_y[index_when_plot,:], 10)


# ## Distance along the z axis
# when_plot_in_hour = 4
# z_pos_distance = 5
# index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
# median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:])
# fig_algebric_distance_ensemble_main_sc_lvlh_z = plt.figure()
# plt.scatter(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:], np.zeros(nb_ensembles)-y_pos_distance, linewidth = 2)    
# plt.title('Algebric LVLH_Z distance after '+str(when_plot_in_hour)+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m')

# ## Histogram
# n, bins, patches = plt.hist(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:], 100,  histtype='stepfilled', alpha = 0.7) # 
# plt.xlim([min(bins), max(bins)])
# plt.ylim([-y_pos_distance, max(n)])
# plt.plot([median_algebric_ditance, median_algebric_ditance], [min(n), max(n)], 'k', linewidth = 2)

# plt.xlabel('Algebric distance (m)')
# print median_algebric_ditance
# print np.percentile(algebric_distance_ensemble_main_sc_lvlh_z[index_when_plot,:], 10)

# plt.show()

# # raise Exception
# # Distribution of the distance between the ensemble sc and the main sc     
# median_distance_ensemble_main_sc = np.zeros(n)
# tenth_percentile_distance_ensemble_main_sc = np.zeros(n)
# twenty_fifth_percentile_distance_ensemble_main_sc = np.zeros(n)
# seventy_fifth_percentile_distance_ensemble_main_sc = np.zeros(n)
# ninetieth_percentile_distance_ensemble_main_sc = np.zeros(n)
# for i in range(0,n,dt_read_output_in_time_step):
#     median_distance_ensemble_main_sc[i] = np.median( distance_ensemble_main_sc[i,:] )
#     tenth_percentile_distance_ensemble_main_sc[i] = np.percentile( distance_ensemble_main_sc[i,:], 10)
#     twenty_fifth_percentile_distance_ensemble_main_sc[i] = np.percentile( distance_ensemble_main_sc[i,:], 25)
#     seventy_fifth_percentile_distance_ensemble_main_sc[i] = np.percentile( distance_ensemble_main_sc[i,:], 75)
#     ninetieth_percentile_distance_ensemble_main_sc[i] = np.percentile( distance_ensemble_main_sc[i,:], 90)

# # Plot the median and the percentiles of the distribution of the distance between the ensemble sc and the main sc        
# x_axis = range(0,n,dt_read_output_in_time_step)
# for i in range(0,n,dt_read_output_in_time_step)        :
#     x_axis[i] = x_axis[i] / 60.0
# fig = plt.figure()
# plt.plot(x_axis,median_distance_ensemble_main_sc, linewidth = 2, color = 'k',label = 'Median')    
# plt.plot(x_axis,tenth_percentile_distance_ensemble_main_sc, linewidth = 2, color = 'r',label = '10th percentile')    
# plt.plot(x_axis,twenty_fifth_percentile_distance_ensemble_main_sc, linewidth = 2, color = 'b',label = '25th percentile')    
# plt.plot(x_axis,seventy_fifth_percentile_distance_ensemble_main_sc, linewidth = 2, color = 'b',label = '75th percentile')    
# plt.plot(x_axis,ninetieth_percentile_distance_ensemble_main_sc, linewidth = 2, color = 'r',label = '90th percentile')    
# #plt.legend()
# plt.ylabel('Distance (m)')
# plt.xlabel('Time (hours)')
# plt.title('Distance between the ensembles and the main spacecraft - Cd = 2.2 +- 0.3 - F10.7 = 150')
# plt.show()

