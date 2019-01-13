import glob
import shutil
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

########### NOTE: running this program will run 20 times the propagator with different GITM directories, and one with the GITM directory that uses the actual IMF. To run it, you don't need anything besides having the input file named input_gitm.d with the density type being "gitm"

# Read file of main sc (the one with a density calculated with the actual IMF)
file_inputname = '../input_gitm.d'
file_input = open(file_inputname,'r')
a_input_file_input = file_input.readlines()
epoch_start = a_input_file_input[1][0:16]
epoch_end = a_input_file_input[2][0:16]
nb_seconds_simu_temp = datetime.datetime(*time.strptime(epoch_end, "%Y-%m-%dT%H:%M")[0:6]) - datetime.datetime(*time.strptime(epoch_start, "%Y-%m-%dT%H:%M")[0:6])
nb_seconds_simu = nb_seconds_simu_temp.days * 24 * 3600 + nb_seconds_simu_temp.seconds
dt = int(float(a_input_file_input[3][0:10]))
nb_time_steps = (int)(nb_seconds_simu / dt) + 1

# eci_r_main_sc = np.zeros([nb_time_steps, 3])
# eci_v_main_sc = np.zeros([nb_time_steps, 3])
# # Run this satellite
# directory_input = "../"
# change_line_output_file = 0
# change_line_gitm_directory = 0
# nb_lines_header = 10
# line_number = 0
# for line in fileinput.input(directory_input+"input_gitm.d", inplace = True):
#     line_number = line_number + 1
#     if ( ( line_number != 9 ) & ( line_number != 55 ) ):
#         print "%s" % (line),
#     elif ( line_number == 55 ):
#         print "%s\n" % "/raid3/Gitm/Runs/20150317/data/",
#     elif ( line_number == 9 ):
#         print "%s, CYGNSS_2gitm.txt\n" % "Ccmc_actual_IMF.txt",

# os.chdir("../")
# os.system("/usr/local/mpi/bin/mpirun -np 1 run_moat")
# directory_output = "./output/run_Ccmc_actual_IMF/Ccmc_actual_IMF/"
# file_to_read = open(directory_output + "Ccmc_actual_IMF.txt", "r")
# read_file_to_read = file_to_read.readlines()
# date = []
# #density_main_sc = np.zeros([nb_time_steps])
# #file_density = open(directory_output + "density_Ccmc_actual_IMF.txt", "r")
# #read_file_density = file_density.readlines()
# for j in range(nb_time_steps):
#     eci_r_main_sc[j, 0] = np.float(read_file_to_read[j+nb_lines_header].split()[2])
#     eci_r_main_sc[j, 1] = np.float(read_file_to_read[j+nb_lines_header].split()[3])
#     eci_r_main_sc[j, 2] = np.float(read_file_to_read[j+nb_lines_header].split()[4])
#     eci_v_main_sc[j, 0] = np.float(read_file_to_read[j+nb_lines_header].split()[5])
#     eci_v_main_sc[j, 1] = np.float(read_file_to_read[j+nb_lines_header].split()[6])
#     eci_v_main_sc[j, 2] = np.float(read_file_to_read[j+nb_lines_header].split()[7])

#     date.append(read_file_to_read[j+nb_lines_header].split()[0] + ' ' + read_file_to_read[j+nb_lines_header].split()[1])
#     # if (j > 0): # first time step is not in the density file
#     #     density_main_sc[j] = np.float(read_file_density[j-1].split()[1]) / 10**9

# # "Ensembles"
nb_gitm_runs = 20
# gitm_directory_array = []
# output_file_array = []
# eci_r_ensemble = np.zeros([nb_gitm_runs, nb_time_steps, 3])
# eci_v_ensemble = np.zeros([nb_gitm_runs, nb_time_steps, 3])
# algebric_distance_ensemble_main_sc_lvlh_x = np.zeros([nb_gitm_runs, nb_time_steps])
# algebric_distance_ensemble_main_sc_lvlh_y = np.zeros([nb_gitm_runs, nb_time_steps])
# algebric_distance_ensemble_main_sc_lvlh_z = np.zeros([nb_gitm_runs, nb_time_steps])
# #density_ensemble = np.zeros([nb_gitm_runs, nb_time_steps])
# directory_input = "./"
# for i in range(1,nb_gitm_runs+1):
#     if (i < 10):
#         gitm_directory_array.append("/raid3/Gitm/Runs/20150317/Ccmc_Run0"+str(i) + "/data/")
#         output_file_array.append("Ccmc_Run0"+str(i)+".txt")
#     else:
#         gitm_directory_array.append("/raid3/Gitm/Runs/20150317/Ccmc_Run"+str(i) + "/data/")
#         output_file_array.append("Ccmc_Run"+str(i)+".txt")

# for i in range(nb_gitm_runs):
#     change_line_output_file = 0
#     change_line_gitm_directory = 0
#     line_number = 0
#     for line in fileinput.input(directory_input+"input_gitm.d", inplace = True):
#         line_number = line_number + 1
#         if ( ( line_number != 9 ) & ( line_number != 55 ) ):
#             print "%s" % (line),
#         elif ( line_number == 55 ):
#             print "%s\n" % gitm_directory_array[i],
#         elif ( line_number == 9 ):
#             print "%s, CYGNSS_2gitm.txt\n" % output_file_array[i],
#     print "i = "+str(i)
#     os.system("/usr/local/mpi/bin/mpirun -np 1 run_moat")
#     directory_output = "./output/run_" + output_file_array[i].split('.')[0] + "/"+ output_file_array[i].split('.')[0]+ "/"
#     file_to_read = open(directory_output + output_file_array[i], "r")
# #    file_density = open(directory_output + "density_" + output_file_array[i], "r")
#     read_file_to_read = file_to_read.readlines()
# #    read_file_density = file_density.readlines()
#     for j in range(nb_time_steps):
#         eci_r_ensemble[i,j, 0] = np.float(read_file_to_read[j+nb_lines_header].split()[2])
#         eci_r_ensemble[i,j, 1] = np.float(read_file_to_read[j+nb_lines_header].split()[3])
#         eci_r_ensemble[i,j, 2] = np.float(read_file_to_read[j+nb_lines_header].split()[4])
#         # eci_v_ensemble[i,j, 0] = np.float(read_file_to_read[j+nb_lines_header].split()[5])
#         # eci_v_ensemble[i,j, 1] = np.float(read_file_to_read[j+nb_lines_header].split()[6])
#         # eci_v_ensemble[i,j, 2] = np.float(read_file_to_read[j+nb_lines_header].split()[7])
        
#         # Distance between the ensemble sc and the main sc in LVLH (reference sc is the main sc)
#         position_main_sc_eci = eci_r_main_sc[j]
#         position_ensemble_sc_eci =  eci_r_ensemble[i,j]
#         delta_vector_eci = position_main_sc_eci - position_ensemble_sc_eci
#         velocity_main_sc_eci = eci_v_main_sc[j]
#         delta_vector_lvlh_of_main_sc = eci_to_lvlh(position_main_sc_eci, velocity_main_sc_eci, delta_vector_eci)
#         algebric_distance_ensemble_main_sc_lvlh_x[i,j]= delta_vector_lvlh_of_main_sc[0] * 1000 # * 1000 to convert the distance from km to m 
#         algebric_distance_ensemble_main_sc_lvlh_y[i,j]= delta_vector_lvlh_of_main_sc[1] * 1000 # * 1000 to convert the distance from km to m 
#         algebric_distance_ensemble_main_sc_lvlh_z[i,j]= delta_vector_lvlh_of_main_sc[2] * 1000 # * 1000 to convert the distance from km to m 
        
#         # if (j > 0): # first time step is not in the density file
#         #     density_ensemble[i,j] = np.float(read_file_density[j-1].split()[1]) / 10**9



# save_pickle_name = "./python_propgator/GITM_20_runs.pickle"
# with open(save_pickle_name, 'w') as f:
#     pickle.dump([eci_r_main_sc, eci_v_main_sc, eci_r_ensemble, algebric_distance_ensemble_main_sc_lvlh_x, algebric_distance_ensemble_main_sc_lvlh_y, algebric_distance_ensemble_main_sc_lvlh_z], f)

save_pickle_name = "GITM_20_runs.pickle"
with open(save_pickle_name) as f:
    eci_r_main_sc, eci_v_main_sc, eci_r_ensemble, algebric_distance_ensemble_main_sc_lvlh_x, algebric_distance_ensemble_main_sc_lvlh_y, algebric_distance_ensemble_main_sc_lvlh_z = pickle.load(f)

# std_distance_ensemble_main_sc_lvlh_x = np.zeros([nb_time_steps])
# std_distance_ensemble_main_sc_lvlh_y = np.zeros([nb_time_steps])
# std_distance_ensemble_main_sc_lvlh_z = np.zeros([nb_time_steps])
# #std_density_ensemble = np.zeros([nb_time_steps])
# # median_density_ensemble = np.zeros([nb_time_steps])
# # mean_density_ensemble = np.zeros([nb_time_steps])
# for j in range(nb_time_steps):
#     std_distance_ensemble_main_sc_lvlh_x[j] = np.std(algebric_distance_ensemble_main_sc_lvlh_x[:,j])
#     std_distance_ensemble_main_sc_lvlh_y[j] = np.std(algebric_distance_ensemble_main_sc_lvlh_y[:,j])
#     std_distance_ensemble_main_sc_lvlh_z[j] = np.std(algebric_distance_ensemble_main_sc_lvlh_z[:,j])
    
#     # if (j > 1):
#     #     std_density_ensemble[j] = np.std(density_ensemble[:,j])
#     #     median_density_ensemble[j] = np.median(density_ensemble[:,j])
#     #     mean_density_ensemble[j] = np.mean(density_ensemble[:,j])

# # #######################################

# # Plot the standard deviation of the along-track distance as a function of time
# fig = plt.figure(num=None, figsize=(18, 10), dpi=80, facecolor='w', edgecolor='k')
# fig.suptitle('Standard deviation along-track and radial', fontsize = 22)
# plt.rc('font', weight='bold') ## make the labels of the ticks in bold
# ax1 = fig.add_subplot(211)
# ax1.set_title('Along-track', weight = 'bold', fontsize = 20,  y = 1.008) 
# ax1.set_ylabel('Standard deviation (m)', fontsize = 18, weight = 'bold')
# #ax1.set_xlabel('Time', fontsize = 18, weight = 'bold')
# ax1.xaxis.set_ticklabels([])
# ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
# [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
# ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
# ax1.plot(std_distance_ensemble_main_sc_lvlh_x, 'k', linewidth = 2)

# # Plot the standard deviation of the radial distance as a function of time
# ax2 = fig.add_subplot(212)
# ax2.set_title('Radial', weight = 'bold', fontsize = 20,  y = 1.008) 
# ax2.set_ylabel('Standard deviation (m)', fontsize = 18, weight = 'bold')
# ax2.set_xlabel('Real time', fontsize = 18, weight = 'bold')
# #ax2.xaxis.set_ticklabels([])
# ax2.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
# [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
# ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
# ax2.plot(std_distance_ensemble_main_sc_lvlh_z, 'k', linewidth = 2)

# # hour_time_step_xticks = 12
# # second_time_step_xticks = hour_time_step_xticks * 3600
# # xticks = np.arange(0, nb_time_steps , second_time_step_xticks / dt)
# # date_list_str = []
# # date_start = datetime.datetime(*time.strptime(date[0], "%Y/%m/%d %H:%M:%S")[0:6])
# # date_list = [date_start + datetime.timedelta(hours=x) for x in np.arange(0, 100, hour_time_step_xticks)]
# # for i in range(len(xticks)):
# #     date_list_str.append( str(date_list[i]) )
# # ax2.xaxis.set_ticks(xticks)
# # ax1.xaxis.set_ticks(xticks)
# # ax2.xaxis.set_ticklabels(date_list_str)



# Plot the algebric lvlh_x distance at a given time
## Distance along the x axis
when_plot_in_hour = 48
index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_x[:,index_when_plot])
tenth_algebric_ditance = np.percentile(algebric_distance_ensemble_main_sc_lvlh_x[:,index_when_plot], 10)
fig = plt.figure(num=None, figsize=(18, 10), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle('Along-track and Radial Distributions after '+str(when_plot_in_hour) +' hours', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(121)
ax1.set_title('Along-track distance', weight = 'bold', fontsize = 20,  y = 1.008) 
#Algebric LVLH_X distance after '+str(when_plot_in_hour)+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
ax1.set_ylabel('Distribution (# out of '+str(nb_gitm_runs)+')', fontsize = 18, weight = 'bold')
ax1.set_xlabel('Distance (m)', fontsize = 18, weight = 'bold')
#ax1.xaxis.set_ticklabels([])
ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
n, bins, patches = ax1.hist(algebric_distance_ensemble_main_sc_lvlh_x[:,index_when_plot], 100,  histtype='stepfilled', alpha = 0.7) 
ax1.plot([median_algebric_ditance, median_algebric_ditance], [min(n), max(n)], 'k-.', linewidth = 2)
ax1.plot([0,0], [0,np.max(n)], 'r-.', linewidth = 2)
ax1.text(max(bins)/2, np.max(n) - np.max(n)/20., "BEHIND\nREFERENCE",  horizontalalignment='center', verticalalignment = 'center', color = 'red')
ax1.text(min(bins)/2, np.max(n) - np.max(n)/20., "IN FRONT\nOF REFERENCE",  horizontalalignment='center', verticalalignment = 'center', color = 'red')

# ## Distance along the z axis
when_plot_in_hour = 48
index_when_plot = (int) (when_plot_in_hour * 3600L / dt)
median_algebric_ditance = np.median(algebric_distance_ensemble_main_sc_lvlh_z[:,index_when_plot])
tenth_algebric_ditance = np.percentile(algebric_distance_ensemble_main_sc_lvlh_z[:,index_when_plot], 10)
ax3 = fig.add_subplot(122)
ax3.set_title('Radial distance', weight = 'bold', fontsize = 20,  y = 1.008) 
#Algebric LVLH_Z distance after '+str(when_plot_in_hour)+ ' hours - median = '+str("%.2f"  % median_algebric_ditance) + ' m - 10th = '+str("%.2f"  % tenth_algebric_ditance)
#ax3.set_ylabel('Distribution (# out of '+str(nb_gitm_runs)+')', fontsize = 18, weight = 'bold')
ax3.set_xlabel('Distance (m)', fontsize = 18, weight = 'bold')
#ax3.xaxis.set_ticklabels([])
ax3.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
[i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
ax3.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
n, bins, patches = ax3.hist(algebric_distance_ensemble_main_sc_lvlh_z[:,index_when_plot], 100,  histtype='stepfilled', alpha = 0.7) 
ax3.plot([median_algebric_ditance, median_algebric_ditance], [min(n), max(n)], 'k-.', linewidth = 2)
ax3.plot([0,0], [0,np.max(n)], 'r-.', linewidth = 2)
ax3.text(max(bins)/2, np.max(n) - np.max(n)/20., "BELOW\nREFERENCE",  horizontalalignment='center', verticalalignment = 'center', color = 'red')
ax3.text(min(bins)/2, np.max(n) - np.max(n)/20., "ABOVE\nREFERENCE",  horizontalalignment='center', verticalalignment = 'center', color = 'red')



# # Plot the median of the density for the ensembles
# fig = plt.figure(num=None, figsize=(16, 9), dpi=80, facecolor='w', edgecolor='k')
# #fig.suptitle('Standard deviation of the density distribution', fontsize = 22)
# plt.rc('font', weight='bold') ## make the labels of the ticks in bold
# ax2 = fig.add_subplot(211)
# ax2.set_title('Median of the density distribution and density with the actual IMF', weight = 'bold', fontsize = 20,  y = 1.008) 
# ax2.set_ylabel('Median (kg/m^3)', fontsize = 18, weight = 'bold')
# #ax2.set_xlabel('Real time', fontsize = 18, weight = 'bold')
# ax2.xaxis.set_ticklabels([])
# ax2.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
# [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
# ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
# ax2.plot(median_density_ensemble, 'b', linewidth = 2, label = 'Median')
# ax2.plot(density_main_sc, 'r', linewidth = 2, label = 'Actual IMF')
# ax2.legend(loc = 1)

# # Plot the standard deviation of the density for the ensembles
# ax1 = fig.add_subplot(212)
# ax1.set_title('Standard deviation of the density distribution', weight = 'bold', fontsize = 20,  y = 1.008) 
# ax1.set_ylabel('Std deviation (kg/m^3)', fontsize = 18, weight = 'bold')
# ax1.set_xlabel('Real time', fontsize = 18, weight = 'bold')
# ax1.xaxis.set_ticklabels([])
# ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
# [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
# ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
# ax1.plot(std_density_ensemble, 'k', linewidth = 2)
# hour_time_step_xticks = 12
# second_time_step_xticks = hour_time_step_xticks * 3600
# xticks = np.arange(0, nb_time_steps , second_time_step_xticks / dt)
# date_list_str = []
# date_start = datetime.datetime(*time.strptime(date[0], "%Y/%m/%d %H:%M:%S")[0:6])
# date_list = [date_start + datetime.timedelta(hours=x) for x in np.arange(0, 100, hour_time_step_xticks)]
# for i in range(len(xticks)):
#     date_list_str.append( str(date_list[i]) )
# ax1.xaxis.set_ticks(xticks)
# ax2.xaxis.set_ticks(xticks)
# ax1.xaxis.set_ticklabels(date_list_str)





#os.chdir("./python_propagator")
