from get_spaced_colors import *
import pickle
from datetime import datetime, timedelta
import time
from eci_to_lvlh import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from read_input_file import *
from read_output_file import *
import colorsys
from color_table import *

plt.ion()
plt.isinteractive()

# Read the input file of propagator
input_filename = "../input_armada_10satellites_2alt_lvlhini.d"
input_variables, order_input_variables = read_input_file(input_filename)
date_start = input_variables[0]; date_stop = input_variables[1]; dt = input_variables[2]; nb_steps = input_variables[3]; nb_satellites = input_variables[4]; output_file_propagator = input_variables[6]
    
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
    print i,nb_satellites-1
    output_filename = output_file_propagator[i]
    print output_filename
    output_variables, list_output_variables_read = read_output_file(output_filename, list_output_variables_to_read)
    if (i == 0):
        date = output_variables[0]

    position[i] = output_variables[1]
    velocity[i] = output_variables[2]
    longitude[i,:] = output_variables[3]
    latitude[i,:] = output_variables[4]
    altitude[i,:] = output_variables[5]
    true_ano[i,:] = output_variables[6]
    raan[i,:] = output_variables[7]
    arg_perigee[i,:]  = output_variables[8]
    right_asc[i,:]  = output_variables[9]
    local_time[i,:]  = output_variables[10]
    angle_asc_node_to_sat[i,:] = (true_ano[i,:] + arg_perigee[i,:])%360

spacing_s1_to_other_sat = np.zeros([nb_satellites,nb_steps])
label_array = []
sat_for_angle_distance = [1,2,3,4,5,6,7]
print "SPACING_S1_TO_OTHER_SAT"
for i in range(nb_satellites):
    label_array.append("s1 with s"+str(i+1))
for i in sat_for_angle_distance:
    print i, len(sat_for_angle_distance)-1
    for j in range(nb_steps):
        absolute_angular_dist = np.abs( angle_asc_node_to_sat[i,j] - angle_asc_node_to_sat[0,j] )
        spacing_s1_to_other_sat[i, j] = min( [absolute_angular_dist, 360 - absolute_angular_dist] )


sat_for_spacing = [0,1,2,3,4,5,6,7]
sat_for_spacing_no_low_plane = [0,1,3,4,5,6]
spacing = np.zeros([len(sat_for_spacing), nb_steps])
spacing_no_low_plane = np.zeros([len(sat_for_spacing_no_low_plane), nb_steps])
for j in range(nb_steps):
    print "SPACING", j, nb_steps - 1
    angle_asc_node_to_sat_ascending_order = np.zeros([len(sat_for_spacing)])
    angle_asc_node_to_sat_ascending_order = np.sort(angle_asc_node_to_sat[sat_for_spacing,j])
    for i in range(len(sat_for_spacing)):
        if ( i == (len(sat_for_spacing) -1) ):
            spacing[i, j] = angle_asc_node_to_sat_ascending_order[0] + 360. - angle_asc_node_to_sat_ascending_order[i]
        else:
            spacing[i, j] = angle_asc_node_to_sat_ascending_order[i+1] - angle_asc_node_to_sat_ascending_order[i]
    angle_asc_node_to_sat_ascending_order_no_low_plane = np.zeros([len(sat_for_spacing_no_low_plane)])
    angle_asc_node_to_sat_ascending_order_no_low_plane = np.sort(angle_asc_node_to_sat[sat_for_spacing_no_low_plane,j])
    for i in range(len(sat_for_spacing_no_low_plane)):
        if ( i == (len(sat_for_spacing_no_low_plane) -1) ):
            spacing_no_low_plane[i, j] = angle_asc_node_to_sat_ascending_order_no_low_plane[0] + 360. - angle_asc_node_to_sat_ascending_order_no_low_plane[i]
        else:
            spacing_no_low_plane[i, j] = angle_asc_node_to_sat_ascending_order_no_low_plane[i+1] - angle_asc_node_to_sat_ascending_order_no_low_plane[i]


name_pickle = '/raid3/Armada/Charles/python/date_position_+_velocity_+_longitude_+_latitude_+_altitude_+_true_ano_+_raan_+_arg_perigee_+_right_asc_+_local_time_+_angle_asc_node_to_sat_+_spacing_s1_to_other_sat_+_spacing_+_spacing_no_low_plane_armada_2alt_2years_lvlhini.pickle'
with open(name_pickle, 'w') as f:
    pickle.dump([date, position, velocity, longitude, latitude, altitude, true_ano, raan, arg_perigee, right_asc, local_time, angle_asc_node_to_sat, spacing_s1_to_other_sat, spacing, spacing_no_low_plane], f)
raise Exception

name_pickle = '/raid3/Armada/Charles/python/date_position_+_velocity_+_longitude_+_latitude_+_altitude_+_true_ano_+_raan_+_arg_perigee_+_right_asc_+_local_time_+_angle_asc_node_to_sat_+_spacing_s1_to_other_sat_+_spacing_+_spacing_no_low_plane_armada_2alt_2years_lvlhini.pickle'
with open(name_pickle) as f:
    date, position, velocity, longitude, latitude, altitude, true_ano, raan, arg_perigee, right_asc, local_time, angle_asc_node_to_sat, spacing_s1_to_other_sat, spacing, spacing_no_low_plane = pickle.load(f)

raise Exception
min_spacing = np.zeros([nb_steps])
max_spacing = np.zeros([nb_steps])
min_spacing_time = np.zeros([nb_steps])
max_spacing_time = np.zeros([nb_steps])
# min_spacing_no_low_plane = np.zeros([nb_steps])
# max_spacing_no_low_plane = np.zeros([nb_steps])
# min_spacing_no_low_plane_time = np.zeros([nb_steps])
# max_spacing_no_low_plane_time = np.zeros([nb_steps])
mu = 398600.4418
r_e = 6378.137
h = 500
period = 2*np.pi*np.sqrt((r_e+h)**3./mu) / 60.# in minutes
# for j in range(nb_steps):
#     min_spacing[j] = min(spacing[:,j])
#     max_spacing[j] = max(spacing[:,j])
#     min_spacing_time[j] = min_spacing[j] / 360. * period
#     max_spacing_time[j] = max_spacing[j] / 360. * period
#     # min_spacing_no_low_plane[j] = min(spacing_no_low_plane[:,j])
#     # max_spacing_no_low_plane[j] = max(spacing_no_low_plane[:,j])
#     # min_spacing_no_low_plane_time[j] = min_spacing_no_low_plane[j] / 360. * period
#     # max_spacing_no_low_plane_time[j] = max_spacing_no_low_plane[j] / 360. * period

######################################
################################ PLOTS
######################################
#fig = plt.figure(num=None, figsize=(15, 8), dpi=80, facecolor='w', edgecolor='k')
fig = plt.figure(num=None, figsize=(17, 11), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(111)
ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

######## SPACING_S1_TO_OTHER_SAT
#SPACING_S1_TO_OTHER_SAT
#colors = color_table(nb_satellites)
colors = ['c','r','c','b','g','m','k','y'] 
#maker_array = []
linewidth_array = [2,2,1,2,2,2,2,1]
alpha_array = [1,1,0.8,1,1,1,1,0.8]
ax1.set_title('Angular distance to s1 VS time', weight = 'bold', fontsize = 20,  y = 1.008) 
ax1.set_ylabel('Angular distance' + u' (\N{DEGREE SIGN})', fontsize = 18, weight = 'bold')
ax1.set_xlabel('Real Time', fontsize = 18, weight = 'bold')
sat_for_angle_distance = [1,2,3,4,5,6,7]
label_array = []
for i in range(nb_satellites):
    label_array.append("s1 with s"+str(i+1))

for i in sat_for_angle_distance:
    print i
#    if ((i != 2) & (i != 7)):
    ax1.plot(spacing_s1_to_other_sat[i,:], color=colors[i], label = label_array[i], linewidth = linewidth_array[i], alpha = alpha_array[i])
#    if ((i != 2) & (i != 7)):


ax1.legend()
hour_time_step_xticks = 24*3*30
second_time_step_xticks = hour_time_step_xticks * 3600
xticks = np.arange(0, nb_steps , second_time_step_xticks / dt)
date_list_str = []
date_start = datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S")
date_list = [date_start + timedelta(hours=x) for x in np.arange(0, 732*24, hour_time_step_xticks)]
for i in range(len(xticks)):
    date_list_str.append( str(date_list[i])[5:10] )

ax1.xaxis.set_ticks(xticks)
ax1.xaxis.set_ticklabels(date_list_str)#, rotation='vertical')
ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
fig.savefig('./spacing_s1_to_other_sat.png', facecolor=fig.get_facecolor(), edgecolor='none')  
#########################################################################################################################
#########################################################################################################################
# ########## HISTOGRAM OF SPACINGS
for m in range(1,25): # START THE LOOP BEFORE THE CREATION OF THE FIGURE
    fig = plt.figure(num=None, figsize=(15, 8), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('', fontsize = 22)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    ax1 = fig.add_subplot(111)
    ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
    ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

    print m
    month_number_start = m
    nb_months = 1
    nb_time_step_in_n_months = (int)(nb_months * 30 * 24 * 3600. / dt)
    spacing_n_months = np.zeros([nb_time_step_in_n_months*len(sat_for_spacing)]) 
    for j in range(len(sat_for_spacing)):
        print "month: "+str(m)
        for i in range(nb_time_step_in_n_months):
            spacing_n_months[i+j*nb_time_step_in_n_months] = spacing[j, (month_number_start-1)*nb_time_step_in_n_months+i]
    weights = np.ones_like(spacing_n_months / 360. * period)/len(spacing_n_months)
    n3, bins, patches = ax1.hist(spacing_n_months / 360. * period, histtype='stepfilled', alpha = 0.7, bins =  np.arange(0,21, 1), weights = weights*100) 
    ax1.set_title('Histogram of spacing for month '+str(month_number_start), weight = 'bold', fontsize = 20,  y = 1.008) 
    ax1.set_ylabel('Distribution (%)', fontsize = 18, weight = 'bold')
    ax1.set_xlabel('Spacing time (min)', fontsize = 18, weight = 'bold')
    ax1.margins(0,0)
    fig.savefig('./images/spacing_histograms/month_'+str(m), facecolor=fig.get_facecolor(), edgecolor='none')  

                           


# ######## RAAN
# ax1.set_title('RAAN difference', weight = 'bold', fontsize = 20,  y = 1.008) 
# ax1.set_ylabel('Delta RAAN ('+u'\N{DEGREE SIGN}'+')', fontsize = 18, weight = 'bold')
# ax1.set_xlabel('Real Time', fontsize = 18, weight = 'bold')
# # colors = ['b','b','r','b','b','b','b','r','k','k']
# # label_array = ['s1','s2','s3','s4','s5','s6','s7','s8','s9','s10']
# ax1.plot((raan[9,:] - raan[0,:])%360, color='b', label = '83'+u'\N{DEGREE SIGN}/500 km - 81'+u'\N{DEGREE SIGN}/500 km', linewidth = 2)
# ax1.plot((raan[2,:] - raan[0,:])%360, color='r', label = '81'+u'\N{DEGREE SIGN}/475 km - 81'+u'\N{DEGREE SIGN}/500 km', linewidth = 2)



######## SPACING_S1_TO_OTHER_SAT
#SPACING_S1_TO_OTHER_SAT
#colors = color_table(nb_satellites)
# colors = ['b','r','k','c','g','m','y','k']
# linewidth_array = [2,2,2,2,2,2,2,2]
# ax1.set_title('Angular distance to s1 VS time', weight = 'bold', fontsize = 20,  y = 1.008) 
# ax1.set_ylabel('Angular distance (degrees)', fontsize = 18, weight = 'bold')
# ax1.set_xlabel('Real Time', fontsize = 18, weight = 'bold')
# sat_for_angle_distance = [1,2,3,4,5,6,7]
# label_array = []
# for i in sat_for_angle_distance:
#     label_array.append("s1 with s"+str(i+1))
#     print i
#     if ((i != 2) & (i != 7)):
#         ax1.plot(spacing_s1_to_other_sat[i,:], color=colors[i], label = label_array[i], linewidth = linewidth_array[i])
# #    if ((i != 2) & (i != 7)):

######### MIN AND MAX SPACING
# ax1.set_title('Min and max spacing time', weight = 'bold', fontsize = 20,  y = 1.008) 
# ax1.set_ylabel('Spacing Time (min)', fontsize = 18, weight = 'bold')
# ax1.set_xlabel('Real Time', fontsize = 18, weight = 'bold')
# ax1.plot(min_spacing_time, color='b', label = 'Min', linewidth = 2)
# ax1.plot(max_spacing_time, color='r', label = 'Max', linewidth = 2)


# ######### NUMBER OF CONJUNCTIONS
# sat_for_spacing_lower_plane_only = [2,7]
# conjunction_angle = 1.0
# conjunction_down_up_list = []
# for m in range(1, 25):
#     month_number_start = m
#     conjunction_down_up_list.append(m)
#     nb_months = 1
#     nb_time_step_in_n_months = (int)(nb_months * 30 * 24 * 3600. / dt)
#     conjunction_down_up_temp = np.zeros([len(sat_for_spacing_lower_plane_only), len(sat_for_spacing_no_low_plane), nb_time_step_in_n_months]) # when a sat in the lower plane flies below a sat in the upper plane (within +- 1 degree)
#     conjunction_down_up = []
#     angle_asc_node_to_sat
#     conjunction_down_up_sub1list = []
#     for k_count in range(len(sat_for_spacing_lower_plane_only)):
#         k = sat_for_spacing_lower_plane_only[k_count]
#         conjunction_down_up_sub1list.append(k)
#         conjunction_down_up_sub2list = []
#         for j_count in range(len(sat_for_spacing_no_low_plane)):
#             j = sat_for_spacing_no_low_plane[j_count]
#             print k_count,j_count
#             conjunction_down_up_sub2list.append(j)
#             conjunction_down_up_sub3list = []
#             for i in range(nb_time_step_in_n_months):
#                 if ( ( np.abs( angle_asc_node_to_sat[j,(month_number_start-1)*nb_time_step_in_n_months+i] - angle_asc_node_to_sat[k,(month_number_start-1)*nb_time_step_in_n_months+i] ) < conjunction_angle ) | (  np.abs( angle_asc_node_to_sat[j,(month_number_start-1)*nb_time_step_in_n_months+i] - angle_asc_node_to_sat[k,(month_number_start-1)*nb_time_step_in_n_months+i] ) > (360 - conjunction_angle) ) ):
#                     conjunction_down_up_temp[k_count,j_count,i] = i
#                 if (i > 0):
#                     if ( (conjunction_down_up_temp[k_count,j_count,i-1] == 0) & (conjunction_down_up_temp[k_count,j_count,i] != 0) ):
#                         conjunction_down_up_sub4list = []
#                         conjunction_down_up_sub4list.append([i])
#                     if ( (conjunction_down_up_temp[k_count,j_count,i-1] != 0) & (conjunction_down_up_temp[k_count,j_count,i] == 0) ):
#                         conjunction_down_up_sub4list.append([i-1])
#                         conjunction_down_up_sub3list.append(conjunction_down_up_sub4list)
#             conjunction_down_up_sub2list.append(conjunction_down_up_sub3list)
#         conjunction_down_up_sub1list.append(conjunction_down_up_sub2list)
#     conjunction_down_up_list.append(conjunction_down_up_sub1list)
# nb_conjunction_per_down_per_month_per_lower_per_upper_sat = np.zeros([24, len(sat_for_spacing_lower_plane_only), len(sat_for_spacing_no_low_plane)])
# nb_conjunction_per_down_per_month_per_lower_sat = np.zeros([24, len(sat_for_spacing_lower_plane_only)])
# total_nb_conjunction_per_down_per_lower_sat =  np.zeros([24, len(sat_for_spacing_lower_plane_only)])
# m_count = -1
# for m in range(1, 2*24,2):
#     m_count = m_count + 1
#     d_count = -1
#     for d in range(1,2*len(sat_for_spacing_lower_plane_only),2):
#         d_count = d_count + 1
#         u_count = -1
#         for u in range(1,2*len(sat_for_spacing_no_low_plane),2):
#             u_count = u_count + 1
#             nb_conjunction_per_down_per_month_per_lower_per_upper_sat[m_count, d_count, u_count] = len(conjunction_down_up_list[m][d][u])
#         nb_conjunction_per_down_per_month_per_lower_sat[m_count, d_count] = sum(nb_conjunction_per_down_per_month_per_lower_per_upper_sat[m_count, d_count,:])

# for m in range(1,25):
#     for d in range(len(sat_for_spacing_lower_plane_only)):
#         total_nb_conjunction_per_down_per_lower_sat[m-1, d] = sum(nb_conjunction_per_down_per_month_per_lower_sat[0:m, d])
# ax1.plot(np.arange(1,25,1), total_nb_conjunction_per_down_per_lower_sat[:,0], linewidth = 2, color = 'b', label = 's3 conjunctions')
# ax1.plot(np.arange(1,25,1), total_nb_conjunction_per_down_per_lower_sat[:,1], linewidth = 2, color = 'r', label = 's8 conjunctions')
# ax1.plot(np.arange(1,25,1), total_nb_conjunction_per_down_per_lower_sat[:,1]+total_nb_conjunction_per_down_per_lower_sat[:,0], linewidth = 2, color = 'k', label = 'Total conjunctions')
# ax1.set_title('Conjunctions', weight = 'bold', fontsize = 20,  y = 1.008) 
# ax1.set_ylabel('Total number of conjunctions', fontsize = 18, weight = 'bold')
# ax1.set_xlabel('Time (months)', fontsize = 18, weight = 'bold')


raise Exception

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################


spacing_file_name = open("armada_time_spacing_everyday.txt", "w+")
for j in range(0, nb_steps,(int)( 24*3600/dt )):
    print >> spacing_file_name, j*dt/3600/24., spacing[0, j]/ 360. * period,spacing[1, j]/ 360. * period,spacing[2, j]/ 360. * period,spacing[3, j]/ 360. * period,spacing[4, j]/ 360. * period,spacing[5, j]/ 360. * period,spacing[6, j]/ 360. * period,spacing[7, j]/ 360. * period
spacing_file_name.close()

# HISTOGRAM OF SPACING
day_number = 365
j = (int)(day_number * 24 * 3600 /dt)
fig = plt.figure(num=None, figsize=(15, 7), dpi=80, facecolor='w', edgecolor='k')
plt.rc('font', weight='bold') ## make the labels of the ticks in bold                                                                                                      
ax2 = fig.add_subplot(111)
ax2.set_title('Histogram of time spacing at day ' + str(day_number), weight = 'bold', fontsize = 20,  y = 1.008) 
ax2.set_ylabel('Number of satellites', fontsize = 18, weight = 'bold')
ax2.set_xlabel('Spacing (minutes)', fontsize = 18, weight = 'bold')
#ax2.xaxis.set_ticklabels([])
ax2.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
[i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
ax2.hist(spacing[:, j]/ 360. * period, histtype='stepfilled', alpha = 0.7, bins = np.arange(0, period +1, period/period)) 
#n, bins, patches = ax2.hist(spacing[:, j]/ 360. * period, 100,  histtype='stepfilled', alpha = 0.7, bins = np.arange(0, period +1, period/10.)) 
ax2.margins(0,0)


raise Exception

# Plots
x_axis = np.zeros(nb_steps)
for i in range(nb_steps):
    x_axis[i] = i * dt / (3600.* 24)
fig = plt.figure(num=None, figsize=(15, 7), dpi=80, facecolor='w', edgecolor='k')
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

## SPACING ALL YEAR
ax2 = fig.add_subplot(111)
ax2.set_title('Min and max time spacing', weight = 'bold', fontsize = 20,  y = 1.008) 
ax2.set_ylabel('Spacing Angle (s)', fontsize = 18, weight = 'bold')
ax2.set_xlabel('Time (days)', fontsize = 18, weight = 'bold')
#ax2.xaxis.set_ticklabels([])
ax2.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
[i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
ax2.set_ylim([0, period])
ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
ax2.plot(x_axis, min_spacing_time,'b-', linewidth = 2, label = 'Min')
ax2.plot(x_axis, max_spacing_time,'r-', linewidth = 2, label= 'Max')
ax2.legend()

raise Exception
### EXAMPLE 3 SATELLITES FOR SHORT PERIOD
first_orbit_to_show = 1500
nb_orbits_to_show = 5
#ANGLE
ax1 = fig.add_subplot(212)
ax1.set_title('Ascending node to satellite', weight = 'bold', fontsize = 20,  y = 1.008) 
ax1.set_ylabel('Angle (' + u'\N{DEGREE SIGN}'+')', fontsize = 18, weight = 'bold')
ax1.set_xlabel('Time (days)', fontsize = 18, weight = 'bold')
#ax1.xaxis.set_ticklabels([])
ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
for i in range(nb_satellites):
    ax1.plot(x_axis[first_orbit_to_show * 95: (first_orbit_to_show + nb_orbits_to_show )* 95 ], angle_asc_node_to_sat[i, first_orbit_to_show * 95: (first_orbit_to_show + nb_orbits_to_show )* 95 ],'k-', linewidth = 2)

## SPACING
ax2 = fig.add_subplot(211)
ax2.set_title('Min and max spacing', weight = 'bold', fontsize = 20,  y = 1.008) 
ax2.set_ylabel('Spacing Angle ('+u'\N{DEGREE SIGN}'+')', fontsize = 18, weight = 'bold')
#ax2.set_xlabel('Time (days)', fontsize = 18, weight = 'bold')
ax2.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
ax2.xaxis.set_ticklabels([])
[i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
ax2.set_ylim([50, 250])
ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
ax2.plot(x_axis[first_orbit_to_show * 95: (first_orbit_to_show + nb_orbits_to_show )* 95 ], min_spacing[ first_orbit_to_show * 95: (first_orbit_to_show + nb_orbits_to_show )* 95 ],'b-', linewidth = 2, label ='Min')
ax2.plot(x_axis[first_orbit_to_show * 95: (first_orbit_to_show + nb_orbits_to_show )* 95 ], max_spacing[ first_orbit_to_show * 95: (first_orbit_to_show + nb_orbits_to_show )* 95 ],'r-', linewidth = 2, label = 'Max')
ax2.legend()

