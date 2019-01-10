import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
import numpy as np
from read_input_file import *
from read_output_file import *
from orbit_average import *

main_input_file_name = '/home/cbv/PropSim/input/main_input/aerie_sink_oneyear_cd5.txt'

# read input file
input_variables, order_input_variables = read_input_file(main_input_file_name)
dt = input_variables[2]; nb_steps = input_variables[3]; satellite_to_plot_path = input_variables[6][0]; satellite_to_plot = input_variables[7][0]; nb_ensembles_coe = input_variables[12]; nb_ensembles_attitude = input_variables[13]; nb_ensembles_cd = input_variables[11]; ensemble_to_plot_temp = input_variables[14]; 

# read main output file
list_output_variables_to_read = ["altitude", "latitude"]
output_variables, list_output_variables_read = read_output_file(satellite_to_plot_path + satellite_to_plot, list_output_variables_to_read)
alt_main_cd5 = output_variables[2][0:output_variables[0].index('2017/01/01 00:00:00')+1]; lat_main_cd5 = output_variables[1][0:output_variables[0].index('2017/01/01 00:00:00')+1]; date_main_cd5 = output_variables[0][0:output_variables[0].index('2017/01/01 00:00:00')+1]


main_input_file_name = '/home/cbv/PropSim/input/main_input/aerie_sink_oneyear.txt'

# read input file
input_variables, order_input_variables = read_input_file(main_input_file_name)
dt = input_variables[2]; nb_steps = input_variables[3]; satellite_to_plot_path = input_variables[6][0]; satellite_to_plot = input_variables[7][0]; nb_ensembles_coe = input_variables[12]; nb_ensembles_attitude = input_variables[13]; nb_ensembles_cd = input_variables[11]; ensemble_to_plot_temp = input_variables[14]; 

# read main output file
list_output_variables_to_read = ["altitude", "latitude"]
output_variables, list_output_variables_read = read_output_file(satellite_to_plot_path + satellite_to_plot, list_output_variables_to_read)
alt_main = output_variables[2]; lat_main = output_variables[1]; date_main = output_variables[0]


# Orbit-average altitude and difference orbit-average altitude
alt_main_cd5_orbit_averaged, time_main_cd5_averaged = orbit_average(alt_main_cd5, lat_main_cd5, date_main_cd5)
alt_main_orbit_averaged, time_main_averaged = orbit_average(alt_main, lat_main, date_main)
delta_alt = []

delta_alt_sub = np.zeros([len(alt_main_orbit_averaged)])
for  isteps in range(len(alt_main_orbit_averaged)):
    delta_alt_sub[isteps] =  alt_main_cd5_orbit_averaged[isteps] - alt_main_orbit_averaged[isteps]

# Plot difference orbit-average altitude after n months VS area ratio
nb_months_for_delta_altitude = 12
area_ratio = np.arange(0, nb_ensembles_cd+1, 1)*2
area_ratio[0] = 1
nb_orbits_for_delta_altitude = (int) (nb_months_for_delta_altitude * 30 * 15.2)
delta_alt_after_n_months = delta_alt_sub[nb_orbits_for_delta_altitude]

height_fig = 8
fontsize_plot = 9
fig = plt.figure(num=None, figsize=(15, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.05, right=0.99, top = 0.7,bottom = 0.2)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(gs[0, 0])
#ax1 = fig.add_subplot(111)
ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
ax1.tick_params(axis='both', which='major',  width = 2, pad = 3)
#ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

ax1.plot( area_ratio,delta_alt_after_n_months)

ax1.margins(0,0)

# Plot difference orbit-average altitude VS time for each surface
height_fig = 8
fontsize_plot = 9
fig = plt.figure(num=None, figsize=(15, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.05, right=0.99, top = 0.7,bottom = 0.2)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(gs[0, 0])
#ax1 = fig.add_subplot(111)
ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
ax1.tick_params(axis='both', which='major',  width = 2, pad = 3)
#ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

ax1.plot( delta_alt_sub )

ax1.margins(0,0)
