#!/home/ridley/Enthought/Canopy_64bit/User/bin/python

import matplotlib.gridspec as gridspec
from datetime import datetime, timedelta
import time
import numpy as np
import matplotlib.pyplot as plt
from read_input_file import *
from read_output_file import *

plt.ion()
plt.isinteractive()

# Read the input file of propagator
input_filename = "../run/input/main_input/aerie.txt"
input_variables, order_input_variables = read_input_file(input_filename)
date_start = input_variables[0]; date_stop = input_variables[1]; dt = input_variables[2]; nb_steps = input_variables[3]; nb_satellites = input_variables[4]; output_file_propagator = input_variables[6]

print(output_file_propagator)

#Read the output files
position = np.zeros([nb_satellites, nb_steps, 3])
velocity = np.zeros([nb_satellites, nb_steps, 3])
longitude = np.zeros([nb_satellites, nb_steps])
latitude = np.zeros([nb_satellites, nb_steps])
altitude = np.zeros([nb_satellites, nb_steps])
true_ano = np.zeros([nb_satellites, nb_steps])
raan = np.zeros([nb_satellites, nb_steps])
arg_perigee = np.zeros([nb_satellites, nb_steps])
right_asc = np.zeros([nb_satellites, nb_steps])
local_time = np.zeros([nb_satellites, nb_steps])
angle_asc_node_to_sat = np.zeros([nb_satellites, nb_steps])
date = []
list_output_variables_to_read = ["position","velocity","longitude","latitude","altitude","raan","true_anomaly", "arg_perigee", "right_asc", "local_time"]
print "READ OUTPUT"
for i in range(nb_satellites):
    print 'Read satellite ' +str(i+1) + ' out of ' + str(nb_satellites) + ' satellites.'
    output_filename = output_file_propagator[i]
    print(output_filename)
    output_variables, list_output_variables_read = read_output_file(output_filename, list_output_variables_to_read)
    if (i == 0):
        date = output_variables[0]
    true_ano[i,:] = output_variables[6]
    arg_perigee[i,:]  = output_variables[8]
    angle_asc_node_to_sat[i,:] = (true_ano[i,:] + arg_perigee[i,:])%360

# Spacing relative to s1 (same as in animation bur going from -180 to +180)
spacing_relative_s1_minus180_to_180 = np.zeros([nb_satellites, nb_steps])
for i in range(nb_steps):
    for j in range(nb_satellites):
        spacing_relative_s1_minus180_to_180_temp = (angle_asc_node_to_sat[j,i] - angle_asc_node_to_sat[0,i])%360
        if (spacing_relative_s1_minus180_to_180_temp > 180):
            spacing_relative_s1_minus180_to_180[j,i] = spacing_relative_s1_minus180_to_180_temp - 360
        else:
            spacing_relative_s1_minus180_to_180[j,i] = spacing_relative_s1_minus180_to_180_temp

nSpaces = nb_satellites * (nb_satellites) / 2
relative_spacing = np.zeros([nb_steps, nSpaces])
relative_spacing = np.zeros(nb_steps*nSpaces)
iSpace = 0
for i in range(nb_steps):
    for j in range(nb_satellites):
        for k in range(nb_satellites):
            if (k > j):
                relative_spacing[iSpace] = (angle_asc_node_to_sat[k,i] - angle_asc_node_to_sat[j,i] + 360)%360
                if (relative_spacing[iSpace] > 180.0):
                    relative_spacing[iSpace] = 360.0 - relative_spacing[iSpace]
                iSpace=iSpace+1

plt.hist(relative_spacing/360.0*95.0, 50)
plt.savefig('./hist.png')  


#raise Exception
######################################
################################ PLOTS
######################################
height_fig = 1
fontsize_plot = 9
fig = plt.figure(num=None, figsize=(17, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.05, right=0.99, top = 0.7,bottom = 0.2)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(gs[0, 0])
ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
ax1.tick_params(axis='both', which='major',  width = 2, pad = 3)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

######## SPACING_S1_TO_OTHER_SAT
colors = ['c','r','c','b','g','m','k','y'] 
linewidth_array = [2,2,1,2,2,2,2,1]
alpha_array = [1,1,0.8,1,1,1,1,0.8]
ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) #Angular distance to s1 VS time
ax1.set_ylabel('Angular\ndistance' + u' (\N{DEGREE SIGN})', fontsize = fontsize_plot, weight = 'bold')
ax1.set_xlabel('', fontsize = fontsize_plot, weight = 'bold')
sat_for_angle_distance = [1,2,3,4,5,6,7]
label_array = []
for i in range(nb_satellites):
    label_array.append("s"+str(i+1))

for i in sat_for_angle_distance:
    ax1.plot(spacing_relative_s1_minus180_to_180[i,:],'.', color=colors[i], linewidth = linewidth_array[i], alpha = alpha_array[i],markersize = 0.5)
    ax1.plot([0,0],[0,0], color=colors[i], label = label_array[i], linewidth = linewidth_array[i], alpha = alpha_array[i])

ax1.legend(ncol = len(sat_for_angle_distance), bbox_to_anchor=(0.5, 1.2), loc = 10,borderpad = 0.1, frameon = False, fontsize = fontsize_plot)
hour_time_step_xticks = 24*3*30
second_time_step_xticks = hour_time_step_xticks * 3600
xticks = np.arange(0, nb_steps , second_time_step_xticks / dt)
date_list_str = []
date_start = datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S")
date_list = [date_start + timedelta(hours=x) for x in np.arange(0, 732*24, hour_time_step_xticks)]
for i in range(len(xticks)):
    date_list_str.append( str(date_list[i])[5:10] )
ax1.yaxis.set_ticks(np.arange(-180,181,120))
for tick in ax1.yaxis.set_ticks(np.arange(-180,181,120)):
    tick.label.set_fontsize(9) 
ax1.xaxis.set_ticks(xticks)
ax1.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
fig.savefig('./spacing_relative_s1_minus180_to_180_17by'+str(height_fig)+'.png', facecolor=fig.get_facecolor(), edgecolor='none')  
