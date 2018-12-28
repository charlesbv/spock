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
from get_prop_dir import *
import matplotlib.gridspec as gridspec
import pickle
from datetime import datetime, timedelta
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from read_input_file import *
from read_output_file import *
import colorsys

plt.ion()
plt.isinteractive()

# Read the input file of propagator
input_filename = "/home/cbv/PropSim/input/main_input/aerie_other_angle.txt"
input_variables, order_input_variables = read_input_file(input_filename)
date_start = input_variables[0]; date_stop = input_variables[1]; dt = input_variables[2]; nb_steps = input_variables[3]; nb_satellites = input_variables[4]; output_file_propagator = input_variables[6]
    
# # #Read the output files
# position = np.zeros([nb_satellites, nb_steps, 3])
# velocity = np.zeros([nb_satellites, nb_steps, 3])
# longitude = np.zeros([nb_satellites, nb_steps])
# latitude = np.zeros([nb_satellites, nb_steps])
# altitude = np.zeros([nb_satellites, nb_steps])
# true_ano = np.zeros([nb_satellites, nb_steps])
# raan = np.zeros([nb_satellites, nb_steps])
# arg_perigee = np.zeros([nb_satellites, nb_steps])
# right_asc = np.zeros([nb_satellites, nb_steps])
# local_time = np.zeros([nb_satellites, nb_steps])
# angle_asc_node_to_sat = np.zeros([nb_satellites, nb_steps])
# date = []
# list_output_variables_to_read = ["position","velocity","longitude","latitude","altitude","raan","true_anomaly", "arg_perigee", "right_asc", "local_time"]
# print "READ OUTPUT"
# for i in range(nb_satellites):
#     print i,nb_satellites-1
#     output_filename = output_file_propagator[i]  + output_file_propagator[i].split('/')[-2] + '.txt'
#     print output_filename
#     output_variables, list_output_variables_read = read_output_file(output_filename, list_output_variables_to_read)
#     if (i == 0):
#         date = output_variables[0]

#     position[i] = output_variables[1]
#     velocity[i] = output_variables[2]
#     longitude[i,:] = output_variables[3]
#     latitude[i,:] = output_variables[4]
#     altitude[i,:] = output_variables[5]
#     true_ano[i,:] = output_variables[6]
#     raan[i,:] = output_variables[7]
#     arg_perigee[i,:]  = output_variables[8]
#     right_asc[i,:]  = output_variables[9]
#     local_time[i,:]  = output_variables[10]
#     angle_asc_node_to_sat[i,:] = (true_ano[i,:] + arg_perigee[i,:])%360

# spacing_M1_to_other_sat = np.zeros([nb_satellites,nb_steps])
# label_array = []
# sat_for_angle_distance = [1,2,3,4,7]
# print "SPACING_M1_TO_OTHER_SAT"
# for i in range(nb_satellites):
#     label_array.append("M1 with s"+str(i+1))
# for i in sat_for_angle_distance:
#     print i, len(sat_for_angle_distance)-1
#     for j in range(nb_steps):
#         absolute_angular_dist = np.abs( angle_asc_node_to_sat[i,j] - angle_asc_node_to_sat[0,j] )
#         spacing_M1_to_other_sat[i, j] = min( [absolute_angular_dist, 360 - absolute_angular_dist] )


# sat_for_spacing = [0,1,2,3,4,7] #[3,1,0,7]
# sat_for_spacing_no_low_plane = [0,1,3,4] #[3,1,0] 
# sat_for_spacing_only_low_plane = [7]#[2, 7]
# sat_for_spacing_only_other_inclination_plane = [6]#[5, 6]
# spacing = np.zeros([len(sat_for_spacing), nb_steps])
# spacing_no_low_plane = np.zeros([len(sat_for_spacing_no_low_plane), nb_steps])
# spacing_only_low_plane = np.zeros([len(sat_for_spacing_no_low_plane), nb_steps, 2])
# spacing_main_plane_and_other_inc_plane = np.zeros([len(sat_for_spacing_no_low_plane + sat_for_spacing_only_other_inclination_plane), nb_steps])
# spacing_only_other_inclination_plane = np.zeros([len(sat_for_spacing_only_other_inclination_plane), nb_steps, 2]) 
# spacing_only_other_inclination_plane_with_pmpl =  np.zeros([len(sat_for_spacing_only_other_inclination_plane), nb_steps, 2])
# spacing_main_plane_lower_plane_and_other_inc_plane =  np.zeros([len(sat_for_spacing_no_low_plane + sat_for_spacing_only_low_plane + sat_for_spacing_only_other_inclination_plane), nb_steps])
# for j in range(nb_steps):
#     print "SPACING", j, nb_steps - 1
#     # BOTH UPPER AND LOWER PLANE
#     angle_asc_node_to_sat_ascending_order = np.zeros([len(sat_for_spacing)])
#     angle_asc_node_to_sat_ascending_order = np.sort(angle_asc_node_to_sat[sat_for_spacing,j])
#     for i in range(len(sat_for_spacing)):
#         if ( i == (len(sat_for_spacing) -1) ):
#             spacing[i, j] = angle_asc_node_to_sat_ascending_order[0] + 360. - angle_asc_node_to_sat_ascending_order[i]
#         else:
#             spacing[i, j] = angle_asc_node_to_sat_ascending_order[i+1] - angle_asc_node_to_sat_ascending_order[i]
#     # UPPER PLANE ONLY
#     angle_asc_node_to_sat_ascending_order_no_low_plane = np.zeros([len(sat_for_spacing_no_low_plane)])
#     angle_asc_node_to_sat_ascending_order_no_low_plane = np.sort(angle_asc_node_to_sat[sat_for_spacing_no_low_plane,j])
#     for i in range(len(sat_for_spacing_no_low_plane)):
#         if ( i == (len(sat_for_spacing_no_low_plane) -1) ):
#             spacing_no_low_plane[i, j] = angle_asc_node_to_sat_ascending_order_no_low_plane[0] + 360. - angle_asc_node_to_sat_ascending_order_no_low_plane[i]
#         else:
#             spacing_no_low_plane[i, j] = angle_asc_node_to_sat_ascending_order_no_low_plane[i+1] - angle_asc_node_to_sat_ascending_order_no_low_plane[i]
#     # LOWER PLANE ONLY: !!!!!! WORKS ONLY IF THERE ARE TWO SATELLITES IN THE LOWER PLANE
#     if j > 0:
#         where_is_first_lower_sat = np.where( angle_asc_node_to_sat_ascending_order == angle_asc_node_to_sat[sat_for_spacing_only_low_plane[0],j] )
#         where_is_first_lower_sat = where_is_first_lower_sat[0]
#         spacing_only_low_plane[0,j,1] = spacing[where_is_first_lower_sat, j]
#         spacing_only_low_plane[0,j,0] = spacing[where_is_first_lower_sat-1, j]
#         # where_is_second_lower_sat = np.where( angle_asc_node_to_sat_ascending_order == angle_asc_node_to_sat[sat_for_spacing_only_low_plane[1],j] )
#         # where_is_second_lower_sat = where_is_second_lower_sat[0]
#         # spacing_only_low_plane[1,j,1] = spacing[where_is_second_lower_sat, j]
#         # spacing_only_low_plane[1,j,0] = spacing[where_is_second_lower_sat-1, j]
#         # if ( ( where_is_second_lower_sat == where_is_first_lower_sat + 1 ) | ( ( where_is_first_lower_sat == len(sat_for_spacing) - 1 ) & ( where_is_second_lower_sat == 0 ) ) ):
#         #     spacing_only_low_plane[1,j,0] = 1000000
#         # if ( (  where_is_second_lower_sat == where_is_first_lower_sat - 1 ) | ( ( where_is_second_lower_sat == len(sat_for_spacing) - 1 ) & ( where_is_first_lower_sat == 0 ) ) ):
#         #     spacing_only_low_plane[1,j,1] = 1000000
        

#     # OTHER INCLINATION PLANE ONLY - WITH PM !!!!!! WORKS ONLY IF THERE ARE TWO SATELLITES IN THE OTHER INCLINATION PLANE
#     angle_asc_node_to_sat_ascending_order = np.zeros([len(sat_for_spacing_no_low_plane + sat_for_spacing_only_other_inclination_plane)])
#     angle_asc_node_to_sat_ascending_order = np.sort( angle_asc_node_to_sat[sat_for_spacing_no_low_plane + sat_for_spacing_only_other_inclination_plane,j] )

#     if j > 0:
#         where_is_first_other_inclination_sat = np.where( angle_asc_node_to_sat_ascending_order == angle_asc_node_to_sat[sat_for_spacing_only_other_inclination_plane[0],j] )
#         where_is_first_other_inclination_sat = where_is_first_other_inclination_sat[0]
#         if len(where_is_first_other_inclination_sat) == 1:
#             for i in range(len(sat_for_spacing_no_low_plane  + sat_for_spacing_only_other_inclination_plane)):
#                 if ( i == (len(sat_for_spacing_no_low_plane + sat_for_spacing_only_other_inclination_plane) -1) ):
#                     spacing_main_plane_and_other_inc_plane[i, j] = angle_asc_node_to_sat_ascending_order[0] + 360. - angle_asc_node_to_sat_ascending_order[i]
#                 else:
#                     spacing_main_plane_and_other_inc_plane[i, j] = angle_asc_node_to_sat_ascending_order[i+1] - angle_asc_node_to_sat_ascending_order[i]

#             spacing_only_other_inclination_plane[0,j,1] = spacing_main_plane_and_other_inc_plane[where_is_first_other_inclination_sat, j]
#             spacing_only_other_inclination_plane[0,j,0] = spacing_main_plane_and_other_inc_plane[where_is_first_other_inclination_sat-1, j]
#             # where_is_second_other_inclination_sat = np.where( angle_asc_node_to_sat_ascending_order == angle_asc_node_to_sat[sat_for_spacing_only_other_inclination_plane[1],j] )
#             # where_is_second_other_inclination_sat = where_is_second_other_inclination_sat[0]
#             # if len(where_is_second_other_inclination_sat) == 1:
#             #     spacing_only_other_inclination_plane[1,j,1] = spacing_main_plane_and_other_inc_plane[where_is_second_other_inclination_sat, j]
#             #     spacing_only_other_inclination_plane[1,j,0] = spacing_main_plane_and_other_inc_plane[where_is_second_other_inclination_sat-1, j]
#             #     if ( ( where_is_second_other_inclination_sat == where_is_first_other_inclination_sat + 1 ) | ( ( where_is_first_other_inclination_sat == len(sat_for_spacing_no_low_plane + sat_for_spacing_only_other_inclination_plane) - 1 ) & ( where_is_second_other_inclination_sat == 0 ) ) ):
#             #         spacing_only_other_inclination_plane[1,j,0] = 1000000
#             #     if ( (  where_is_second_other_inclination_sat == where_is_first_other_inclination_sat - 1 ) | ( ( where_is_second_other_inclination_sat == len(sat_for_spacing_no_low_plane + sat_for_spacing_only_other_inclination_plane) - 1 ) & ( where_is_first_other_inclination_sat == 0 ) ) ):
#             #         spacing_only_other_inclination_plane[1,j,1] = 1000000
        
            

#     # OTHER INCLINATION PLANE ONLY - WITH PM AND PL: !!!!!! WORKS ONLY IF THERE ARE TWO SATELLITES IN THE OTHER INCLINATION PLANE
#     angle_asc_node_to_sat_ascending_order = np.zeros([len(sat_for_spacing_no_low_plane + sat_for_spacing_only_other_inclination_plane + sat_for_spacing_only_low_plane)])
#     angle_asc_node_to_sat_ascending_order = np.sort( angle_asc_node_to_sat[sat_for_spacing_no_low_plane + sat_for_spacing_only_low_plane + sat_for_spacing_only_other_inclination_plane,j] )

#     if j > 0:
#         where_is_first_other_inclination_sat = np.where( angle_asc_node_to_sat_ascending_order == angle_asc_node_to_sat[sat_for_spacing_only_other_inclination_plane[0],j] )
#         where_is_first_other_inclination_sat = where_is_first_other_inclination_sat[0]
#         if len(where_is_first_other_inclination_sat) == 1:
#             for i in range(len(sat_for_spacing_no_low_plane + sat_for_spacing_only_low_plane  + sat_for_spacing_only_other_inclination_plane)):
#                 if ( i == (len(sat_for_spacing_no_low_plane + sat_for_spacing_only_low_plane + sat_for_spacing_only_other_inclination_plane) -1) ):
#                     spacing_main_plane_lower_plane_and_other_inc_plane[i, j] = angle_asc_node_to_sat_ascending_order[0] + 360. - angle_asc_node_to_sat_ascending_order[i]
#                 else:
#                     spacing_main_plane_lower_plane_and_other_inc_plane[i, j] = angle_asc_node_to_sat_ascending_order[i+1] - angle_asc_node_to_sat_ascending_order[i]

#             spacing_only_other_inclination_plane_with_pmpl[0,j,1] = spacing_main_plane_lower_plane_and_other_inc_plane[where_is_first_other_inclination_sat, j]
#             spacing_only_other_inclination_plane_with_pmpl[0,j,0] = spacing_main_plane_lower_plane_and_other_inc_plane[where_is_first_other_inclination_sat-1, j]
#             # where_is_second_other_inclination_sat = np.where( angle_asc_node_to_sat_ascending_order == angle_asc_node_to_sat[sat_for_spacing_only_other_inclination_plane[1],j] )
#             # where_is_second_other_inclination_sat = where_is_second_other_inclination_sat[0]
#             # if len(where_is_second_other_inclination_sat) == 1:
#             #     spacing_only_other_inclination_plane_with_pmpl[1,j,1] = spacing_main_plane_lower_plane_and_other_inc_plane[where_is_second_other_inclination_sat, j]
#             #     spacing_only_other_inclination_plane_with_pmpl[1,j,0] = spacing_main_plane_lower_plane_and_other_inc_plane[where_is_second_other_inclination_sat-1, j]
#             #     if ( ( where_is_second_other_inclination_sat == where_is_first_other_inclination_sat + 1 ) | ( ( where_is_first_other_inclination_sat == len(sat_for_spacing_no_low_plane + sat_for_spacing_only_other_inclination_plane + sat_for_spacing_only_low_plane) - 1 ) & ( where_is_second_other_inclination_sat == 0 ) ) ):
#             #         spacing_only_other_inclination_plane_with_pmpl[1,j,0] = 1000000
#             #     if ( (  where_is_second_other_inclination_sat == where_is_first_other_inclination_sat - 1 ) | ( ( where_is_second_other_inclination_sat == len(sat_for_spacing_no_low_plane + sat_for_spacing_only_other_inclination_plane + sat_for_spacing_only_low_plane) - 1 ) & ( where_is_first_other_inclination_sat == 0 ) ) ):
#             #         spacing_only_other_inclination_plane_with_pmpl[1,j,1] = 1000000
        
            

# # # Spacing relative to M1 (same as in animation but going from -180 to +180)
# # spacing_relative_M1_minus180_to_180 = np.zeros([nb_satellites, nb_steps])
# # for i in range(nb_steps):
# #     for j in range(nb_satellites):
# #         spacing_relative_M1_minus180_to_180_temp = (angle_asc_node_to_sat[j,i] - angle_asc_node_to_sat[0,i])%360
# #         if (spacing_relative_M1_minus180_to_180_temp > 180):
# #             spacing_relative_M1_minus180_to_180[j,i] = spacing_relative_M1_minus180_to_180_temp - 360
# #         else:
# #             spacing_relative_M1_minus180_to_180[j,i] = spacing_relative_M1_minus180_to_180_temp


# # name_pickle = '/raid3/Armada/Charles/python/aerie_8sat_ok.pickle'
# # with open(name_pickle, 'w') as f:
# #     pickle.dump([date, position, velocity, longitude, latitude, altitude, true_ano, raan, arg_perigee, right_asc, local_time, angle_asc_node_to_sat, spacing_M1_to_other_sat, spacing, spacing_no_low_plane, spacing_only_low_plane, spacing_only_other_inclination_plane, spacing_only_other_inclination_plane_with_pmpl, spacing_relative_M1_minus180_to_180], f)
# # raise Exception

#### HERE spacing, spacing_no_low_plane, spacing_only_other_inclination_plane, AND spacing_only_other_inclination_plane_with_pmpl WRE SAVED FOR THE CONFIGURATION 422 (NO LOSS OF ANY SATELLITE)
name_pickle = '/raid3/Armada/Charles/python/aerie_8sat_ok.pickle'
with open(name_pickle) as f:
    date, position, velocity, longitude, latitude, altitude, true_ano, raan, arg_perigee, right_asc, local_time, angle_asc_node_to_sat, spacing_M1_to_other_sat, spacing, spacing_no_low_plane,spacing_only_low_plane, spacing_only_other_inclination_plane,spacing_only_other_inclination_plane_with_pmpl, spacing_relative_M1_minus180_to_180 = pickle.load(f)

raise Exception

# NUMBER OF ORBITS WITH A SPACING BELOW 5, 10, OR 20 MINUTES
sat_for_spacing = [0,1,2,3,4,7] #[3,1,0,7]
sat_for_spacing_no_low_plane = [0,1,3,4] #[3,1,0] 
sat_for_spacing_only_low_plane = [7]#[2, 7]
sat_for_spacing_only_other_inclination_plane = [6]#[5, 6]
list_nb_spacing_below_10_min_only_upper_plane = []
list_nb_spacing_below_10_min_only_lower_plane = []
list_nb_spacing_below_10_min_only_other_inclination_plane = []
list_nb_spacing_below_10_min_only_other_inclination_plane_with_pmpl = []
mu = 398600.4418
r_e = 6378.137
h = 500
period = 2*np.pi*np.sqrt((r_e+h)**3./mu) / 60.# in minutes
for istep in range(1,nb_steps):
    if ( (latitude[0, istep] >=0) & (latitude[0, istep-1] < 0) ): # when the reference satellite (here chosen to be the first sat) crosses the equator northward, look at the spacing between the satellites
        # ONLY UPPER PLANE (= MAIN PLANE)
        list_nb_spacing_below_10_min_only_upper_plane_temp = []
        for isat in range(len(sat_for_spacing_no_low_plane)):
            if spacing_no_low_plane[isat, istep] < 10. / period * 360:
                list_nb_spacing_below_10_min_only_upper_plane_temp.append([istep, sat_for_spacing_no_low_plane[isat]])
        list_nb_spacing_below_10_min_only_upper_plane.append(list_nb_spacing_below_10_min_only_upper_plane_temp)
        
        # LOWER PLANE SPACING WITH UPPER PLANE (= MAIN PLANE)
        list_nb_spacing_below_10_min_only_lower_plane_temp = []
        for isat in range(len(sat_for_spacing_only_low_plane)):
            for k in range(2):
                if spacing_only_low_plane[isat, istep, k] < 10. / period * 360:
                    list_nb_spacing_below_10_min_only_lower_plane_temp.append([istep, sat_for_spacing_only_low_plane[isat], k])
        list_nb_spacing_below_10_min_only_lower_plane.append(list_nb_spacing_below_10_min_only_lower_plane_temp)
        
        # OTHER INCLINATION PLANE WITH MAIN PLANE
        list_nb_spacing_below_10_min_only_other_inclination_plane_temp = []
        for isat in range(len(sat_for_spacing_only_other_inclination_plane)):
            for k in range(2):
                if spacing_only_other_inclination_plane[isat, istep, k] < 10. / period * 360:
                    list_nb_spacing_below_10_min_only_other_inclination_plane_temp.append([istep, sat_for_spacing_only_other_inclination_plane[isat], k])
        list_nb_spacing_below_10_min_only_other_inclination_plane.append(list_nb_spacing_below_10_min_only_other_inclination_plane_temp)

        # OTHER INCLINATION PLANE WITH MAIN PLANE AND LOWER PLANE
        list_nb_spacing_below_10_min_only_other_inclination_plane_with_pmpl_temp = []
        for isat in range(len(sat_for_spacing_only_other_inclination_plane)):
            for k in range(2):
                if spacing_only_other_inclination_plane_with_pmpl[isat, istep, k] < 10. / period * 360:
                    list_nb_spacing_below_10_min_only_other_inclination_plane_with_pmpl_temp.append([istep, sat_for_spacing_only_other_inclination_plane[isat], k])
        list_nb_spacing_below_10_min_only_other_inclination_plane_with_pmpl.append(list_nb_spacing_below_10_min_only_other_inclination_plane_with_pmpl_temp)
        
nb_spacing_below_10_min_only_upper_plane = np.zeros([len(list_nb_spacing_below_10_min_only_upper_plane)])
nb_spacing_below_10_min_only_lower_plane = np.zeros([len(list_nb_spacing_below_10_min_only_upper_plane)])
nb_spacing_below_10_min_only_other_inclination_plane = np.zeros([len(list_nb_spacing_below_10_min_only_upper_plane)])
nb_spacing_below_10_min_only_other_inclination_plane_with_pmpl = np.zeros([len(list_nb_spacing_below_10_min_only_upper_plane)])

for iave in range(len(list_nb_spacing_below_10_min_only_upper_plane)):
    nb_spacing_below_10_min_only_upper_plane[iave] = len(list_nb_spacing_below_10_min_only_upper_plane[iave])
    nb_spacing_below_10_min_only_lower_plane[iave] = len(list_nb_spacing_below_10_min_only_lower_plane[iave])
    nb_spacing_below_10_min_only_other_inclination_plane[iave] = len(list_nb_spacing_below_10_min_only_other_inclination_plane[iave])
    nb_spacing_below_10_min_only_other_inclination_plane_with_pmpl[iave] = len(list_nb_spacing_below_10_min_only_other_inclination_plane_with_pmpl[iave])

nb_day_each_month = [31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30.,31,30,31,31,30,31,30,31]
month_nb_spacing_below_10_min_only_upper_plane = np.zeros([len(nb_day_each_month)])
month_nb_spacing_below_10_min_only_lower_plane = np.zeros([len(nb_day_each_month)])
month_nb_spacing_below_10_min_only_other_inclination_plane = np.zeros([len(nb_day_each_month)])
month_nb_spacing_below_10_min_only_other_inclination_plane_with_pmpl = np.zeros([len(nb_day_each_month)])
nb_day_per_month = np.sum(nb_day_each_month) / len(nb_day_each_month)
nb_orbits_in_one_month =  nb_day_per_month * 24 * 60. /period 
nb_orbits_elapsed = 0
for imonth in range(len(nb_day_each_month)):
    nb_orbits_in_this_month = (int) ( nb_day_each_month[imonth] * 24 * 60. / period )
    month_nb_spacing_below_10_min_only_upper_plane[imonth] = np.sum(nb_spacing_below_10_min_only_upper_plane[ nb_orbits_elapsed : nb_orbits_elapsed + nb_orbits_in_this_month])
    month_nb_spacing_below_10_min_only_lower_plane[imonth] = np.sum(nb_spacing_below_10_min_only_lower_plane[ nb_orbits_elapsed : nb_orbits_elapsed + nb_orbits_in_this_month])
    month_nb_spacing_below_10_min_only_other_inclination_plane[imonth] = np.sum(nb_spacing_below_10_min_only_other_inclination_plane[ nb_orbits_elapsed : nb_orbits_elapsed + nb_orbits_in_this_month])
    month_nb_spacing_below_10_min_only_other_inclination_plane_with_pmpl[imonth] = np.sum(nb_spacing_below_10_min_only_other_inclination_plane_with_pmpl[ nb_orbits_elapsed : nb_orbits_elapsed + nb_orbits_in_this_month])
    nb_orbits_elapsed = nb_orbits_elapsed + nb_orbits_in_this_month


############ ONLY FOR CONFIGURATION 311 (LOSS OF ONE SATELLITE IN EACH OF THE THREE PLANES)
# # #### HERE spacing, spacing_no_low_plane, spacing_only_other_inclination_plane, AND spacing_only_other_inclination_plane_with_pmpl WRE SAVED FOR THE CONFIGURATION 311: ...lose_M1 or ...lose_M2 or ...lose_M3 or ...lose_M4
# # name_pickle = '/raid3/Armada/Charles/python/aerie_311_lose_M1.pickle'
# # with open(name_pickle, 'w') as f:
# #     pickle.dump([ spacing_no_low_plane, spacing_only_low_plane, spacing_only_other_inclination_plane, spacing_only_other_inclination_plane_with_pmpl, month_nb_spacing_below_10_min_only_upper_plane, month_nb_spacing_below_10_min_only_lower_plane, month_nb_spacing_below_10_min_only_other_inclination_plane, month_nb_spacing_below_10_min_only_other_inclination_plane_with_pmpl ], f)
# # raise Exception
month_conf311_nb_spacing_below_10_min_only_upper_plane = np.zeros([4, len(nb_day_each_month)])
month_conf311_nb_spacing_below_10_min_only_lower_plane = np.zeros([4, len(nb_day_each_month)])
month_conf311_nb_spacing_below_10_min_only_other_inclination_plane = np.zeros([4, len(nb_day_each_month)])
month_conf311_nb_spacing_below_10_min_only_other_inclination_plane_with_pmpl = np.zeros([4, len(nb_day_each_month)])
#### HERE spacing, spacing_no_low_plane, spacing_only_other_inclination_plane, AND spacing_only_other_inclination_plane_with_pmpl WRE SAVED FOR THE CONFIGURATION 311: ...lose_M1 or ...lose_M2 or ...lose_M3 or ...lose_M4
for i in range(4):
    name_pickle = '/raid3/Armada/Charles/python/aerie_311_lose_M' + str(i+1) + '.pickle'
    with open(name_pickle) as f:
        spacing_no_low_plane, spacing_only_low_plane, spacing_only_other_inclination_plane, spacing_only_other_inclination_plane_with_pmpl, month_conf311_nb_spacing_below_10_min_only_upper_plane[i, :], month_conf311_nb_spacing_below_10_min_only_lower_plane[i, :], month_conf311_nb_spacing_below_10_min_only_other_inclination_plane[i, :], month_conf311_nb_spacing_below_10_min_only_other_inclination_plane_with_pmpl[i, :] = pickle.load(f)

month_conf311_average_nb_spacing_below_10_min_only_upper_plane = np.zeros([len(nb_day_each_month)])
month_conf311_average_nb_spacing_below_10_min_only_lower_plane = np.zeros([len(nb_day_each_month)])
month_conf311_average_nb_spacing_below_10_min_only_other_inclination_plane = np.zeros([len(nb_day_each_month)])
month_conf311_average_nb_spacing_below_10_min_only_other_inclination_plane_with_pmpl = np.zeros([len(nb_day_each_month)])
for imonth in range(len(nb_day_each_month)):
    month_conf311_average_nb_spacing_below_10_min_only_upper_plane[imonth] = np.mean(month_conf311_nb_spacing_below_10_min_only_upper_plane[:, imonth])
    month_conf311_average_nb_spacing_below_10_min_only_lower_plane[imonth] = np.mean(month_conf311_nb_spacing_below_10_min_only_lower_plane[:, imonth])
    month_conf311_average_nb_spacing_below_10_min_only_other_inclination_plane[imonth] = np.mean(month_conf311_nb_spacing_below_10_min_only_other_inclination_plane[:, imonth])
    month_conf311_average_nb_spacing_below_10_min_only_other_inclination_plane_with_pmpl[imonth] = np.mean(month_conf311_nb_spacing_below_10_min_only_other_inclination_plane_with_pmpl[:, imonth])
############ END OF ONLY FOR CONFIGURATION 311 (LOSS OF ONE SATELLITE IN EACH OF THE THREE PLANES)

######################################
################################ PLOTS
######################################
#fig = plt.figure(num=None, figsize=(15, 8), dpi=80, facecolor='w', edgecolor='k')
width_fig = 3# 17
height_fig = 10#width_fig*3/4
fontsize_plot = 11 # 9
fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.05, right=0.99, top = 0.7,bottom = 0.2)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(gs[0, 0])
#ax1 = fig.add_subplot(111)
ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
#ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

######### SAT BELOW ANOTHER
colors = ['b','b','r','b','b','k','k','r'] 
label_array = ['M1','M2','L1','M3','M4','I1','I2','L2']
marker_array = ['s','^','s','o','D','s','^','^']
index_plot = 1221003
isat_reference = 4
mu = 398600.4418
r_e = 6378.137
h = altitude[isat_reference,index_plot]
period = 2*np.pi*np.sqrt((r_e+h)**3./mu) # in minutes
isat = 2; ax1.scatter( [(angle_asc_node_to_sat[isat,index_plot] - angle_asc_node_to_sat[isat_reference,index_plot])/360 * period], [altitude[isat,index_plot]], marker = marker_array[isat], label = label_array[isat], s = 150, color= colors[isat],linewidth = 4, zorder = 5)
isat = 3; ax1.scatter( [(angle_asc_node_to_sat[isat,index_plot]-angle_asc_node_to_sat[isat_reference,index_plot])/360 * period], [altitude[isat,index_plot]], marker = marker_array[isat], label = label_array[isat], s = 150, color= colors[isat],linewidth = 4, zorder = 5)
isat = 4; ax1.scatter( [(angle_asc_node_to_sat[isat,index_plot] - angle_asc_node_to_sat[isat_reference,index_plot])/360 * period], [altitude[isat,index_plot]], marker = marker_array[isat], label = label_array[isat], s = 150, color= colors[isat],linewidth = 4, zorder = 5)
ax1.set_xlabel('Separation time (s)', fontsize = fontsize_plot, weight = 'bold')
ax1.set_ylabel('Altitude (km)', fontsize = fontsize_plot, weight = 'bold')
#ax1.set_title('Altitude VS separation time - PM and PL', weight = 'bold', fontsize = 20,  y = 1.04) # y = 1.09 for 2' figure height / 1.04 for normal height 
ax1.legend(ncol = 3, bbox_to_anchor=(0.5, 1.03), loc = 10,borderpad = 0.1, frameon = False, fontsize = fontsize_plot, scatterpoints = 1) #  1.10 for 2' figure height / 1.02 for normal height
gs.update(left=0.2, right=0.95, top = 0.95,bottom = 0.05)
fig_save_name = get_prop_dir(1) + 'output/python_propagator/aerie/altitude_vs_separation_time.png' # '/raid3/Armada/Charles/python/spacing_M1_to_other_sat_17by2.png' #
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./")


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
fig.savefig('./spacing_relative_M1_minus180_to_180_17by'+str(height_fig)+'.png', facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
#########################################################################################################################
#########################################################################################################################
######## SPACING_M1_TO_OTHER_SAT
#colors = color_table(nb_satellites)
colors = ['c','r','c','b','g','m','k','y'] 
linewidth_array = [2,2,1,2,2,2,2,1]
alpha_array = [1,1,0.8,1,1,1,1,0.8]
sat_for_angle_distance = [1,2,3,4,7]
label_array = ['M1','M2','L1','M3','M4','I1','I2','L2']
x_axis = np.arange(0, 24, 24./len(spacing_relative_M1_minus180_to_180[0,:]))

for i in sat_for_angle_distance:
    ax1.plot(x_axis,spacing_relative_M1_minus180_to_180[i,:], '.', color=colors[i], linewidth = linewidth_array[i], alpha = alpha_array[i],markersize = 0.5)
    ax1.plot([0,0],[0,0], color=colors[i], label = label_array[i], linewidth = linewidth_array[i], alpha = alpha_array[i])

ax1.set_xlabel('Time (months)', fontsize = fontsize_plot, weight = 'bold')
ax1.set_ylabel('Angular\ndistance' + u' (\N{DEGREE SIGN})', fontsize = fontsize_plot, weight = 'bold', labelpad=0.0001)
ax1.set_title('Angular spacing between M1 and the other satellites in PM and PL', weight = 'bold', fontsize = 20,  y = 1.04) # y = 1.09 for 2' figure height / 1.04 for normal height 
gs.update(left=0.063, right=1-0.06, top = 0.93,bottom = 0.055) # left=0.045, right=1-0.05, top = 0.85,bottom = 0.1 for 2' figure height / left=0.063, right=1-0.06, top = 0.93,bottom = 0.055 for normal height
ax1.legend(ncol = len(sat_for_angle_distance), bbox_to_anchor=(0.5, 1.02), loc = 10,borderpad = 0.1, frameon = False, fontsize = fontsize_plot) #  1.10 for 2' figure height / 1.02 for normal height
ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
y_ticks = np.arange(-180,181,120)
ax1.yaxis.set_ticks(y_ticks)
for tick in ax1.yaxis.set_ticks(y_ticks):
    tick.label.set_fontsize(fontsize_plot) 
ax2 = ax1.twinx()
y_ax2_ticks = y_ticks / 360. * period
ax2.yaxis.set_ticks(y_ax2_ticks)
ax2.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
ax2.set_ylabel('Time\nseparation (min)', fontsize = fontsize_plot, weight = 'bold', labelpad=22, rotation = -90) # labelpad = 22 for 2' figure height
fig_save_name = get_prop_dir(1) + 'output/python_propagator/aerie/spacing_M1_to_other_sat.png' # '/raid3/Armada/Charles/python/spacing_M1_to_other_sat_17by2.png' #
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./")

######## NUMBER OF ORBITS WITH A SPACING BELOW 10 MINUTES
x_axis = np.arange(0, 24)

# ax1.plot(x_axis, month_nb_spacing_below_10_min_only_upper_plane, color='b', linewidth = 2 , label ='PM')
# ax1.plot(x_axis, month_nb_spacing_below_10_min_only_lower_plane, color='r', linewidth = 2 , label ='PL with PM')
# ax1.plot(x_axis, month_nb_spacing_below_10_min_only_other_inclination_plane, color='k', linewidth = 2 , label ='PI with PM')
# ax1.plot(x_axis, month_nb_spacing_below_10_min_only_other_inclination_plane_with_pmpl, color='k', linewidth = 2 , label ='PI with PM & PL', linestyle = 'dotted')

ax1.plot(x_axis, month_conf311_average_nb_spacing_below_10_min_only_upper_plane, color='b', linewidth = 2 , label ='PM')
ax1.plot(x_axis, month_conf311_average_nb_spacing_below_10_min_only_lower_plane, color='r', linewidth = 2 , label ='PL with PM')
ax1.plot(x_axis, month_conf311_average_nb_spacing_below_10_min_only_other_inclination_plane, color='k', linewidth = 2 , label ='PI with PM')
ax1.plot(x_axis, month_conf311_average_nb_spacing_below_10_min_only_other_inclination_plane_with_pmpl, color='k', linewidth = 2 , label ='PI with PM & PL', linestyle = 'dotted')

ax1.margins(0,1) # autoscale both axes(fist value is for the x axis, second value for the y axis)
ax1.set_ylim([0, np.max([np.max(month_nb_spacing_below_10_min_only_upper_plane), np.max(month_nb_spacing_below_10_min_only_lower_plane), np.max(month_nb_spacing_below_10_min_only_other_inclination_plane), np.max(month_nb_spacing_below_10_min_only_other_inclination_plane_with_pmpl)])])
ax1.set_ylabel('# Orbits per month', fontsize = fontsize_plot, weight = 'bold', labelpad =  0.001)
ax1.set_xlabel('Time (months)', fontsize = fontsize_plot, weight = 'bold')
#ax1.set_title('Number of orbits per month with a separation time below 10 minutes - Configuration 4-2-2', weight = 'bold', fontsize = 20,  y = 1.05) 
gs.update(left=0.040, right=0.997, top = 0.85,bottom = 0.1)
ax1.legend(ncol = 4, bbox_to_anchor=(0.5, 1.06), loc = 10,borderpad = 0.1, frameon = False, fontsize = fontsize_plot) #1.02 for normal figure height # 1.06 for 2' figure height
fig_save_name = '/raid3/Armada/Charles/python/spacing_under_10_min_conf422_17by2.png' #get_prop_dir(1) + 'output/python_propagator/aerie/spacing_under_10_min_conf422.png' # 
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./")


# ########## HISTOGRAM OF SPACINGS
mu = 398600.4418
r_e = 6378.137
h = 500
period = 2*np.pi*np.sqrt((r_e+h)**3./mu) / 60.# in minutes
nb_months = 3
sat_for_spacing = [0,1,2,3,4,7]
sat_for_spacing_no_low_plane = [0,1,3,4]
sat_for_spacing_only_low_plane = [2, 7]
for m in range(1,(nb_months +1)): # START THE LOOP BEFORE THE CREATION OF THE FIGURE
    fig = plt.figure(num=None, figsize=(15, 8), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('', fontsize = 22)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    ax1 = fig.add_subplot(111)
    ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
    ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

    print m
    month_number_start = m
    nb_time_step_in_one_month = (int)( 30 * 24 * 3600. / dt)

    # # SPACING BOTH LOWER AND UPPER PLANE
    # spacing_n_months = np.zeros([nb_time_step_in_one_month*len(sat_for_spacing)]) 
    # for j in range(len(sat_for_spacing)):
    #     print "month: "+str(m)
    #     for i in range(nb_time_step_in_one_month):
    #         spacing_n_months[i+j*nb_time_step_in_one_month] = spacing[j, (month_number_start-1)*nb_time_step_in_one_month+i]
    # weights = np.ones_like(spacing_n_months / 360. * period)/len(spacing_n_months)
    # n3, bins, patches = ax1.hist(spacing_n_months / 360. * period, histtype='stepfilled', alpha = 0.7, bins =  np.arange(0,21, 1), weights = weights*100) 

    # SPACING UPPER PLANE ONLY
    spacing_n_months = np.zeros([nb_time_step_in_one_month*len(sat_for_spacing_no_low_plane)]) 
    for j in range(len(sat_for_spacing_no_low_plane)):
        print "month: "+str(m)
        for i in range(nb_time_step_in_one_month):
            spacing_n_months[i+j*nb_time_step_in_one_month] = spacing_no_low_plane[j, (month_number_start-1)*nb_time_step_in_one_month+i]
    weights = np.ones_like(spacing_n_months / 360. * period)/len(spacing_n_months)
    n3, bins, patches = ax1.hist(spacing_n_months / 360. * period, histtype='stepfilled', alpha = 0.7, bins =  np.arange(0,21, 1), weights = weights*100, color = 'b') 

    # SPACING LOWER PLANE ONLY
    spacing_n_months = np.zeros([nb_time_step_in_one_month*len(sat_for_spacing_only_low_plane)*2])
    for p in range(2):
        for k in range(2):
#            print "month: "+str(m)
            for i in range(nb_time_step_in_one_month):
                spacing_n_months[ i + k*nb_time_step_in_one_month + 2*p*nb_time_step_in_one_month] = spacing_only_low_plane[p, (month_number_start-1)*nb_time_step_in_one_month+i, k]
    weights = np.ones_like(spacing_n_months / 360. * period)/len(spacing_n_months)
    n3, bins, patches = ax1.hist(spacing_n_months / 360. * period, histtype='stepfilled', alpha = 0.7, bins =  np.arange(0,21, 1), weights = weights*100, color = 'r') 
    
    
    ax1.set_title('Histogram of spacing for month '+str(month_number_start), weight = 'bold', fontsize = 20,  y = 1.008) 
    ax1.set_ylabel('Distribution (%)', fontsize = 18, weight = 'bold')
    ax1.set_xlabel('Spacing time (min)', fontsize = 18, weight = 'bold')
    ax1.margins(0,0)
    fig.savefig('./images/spacing_histograms/month_'+str(m), facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

# ######### NUMBER OF CONJUNCTIONS
nb_months = 3
sat_for_spacing = [0,1,2,3,4,7]
sat_for_spacing_no_low_plane = [0,1,3,4]
sat_for_spacing_lower_plane_only = [2,7]
conjunction_angle = 1.0
conjunction_down_up_list = []
for m in range(1, 4):
    month_number_start = m
    conjunction_down_up_list.append(m)
    nb_time_step_in_one_month = (int)( 30 * 24 * 3600. / dt)
    conjunction_down_up_temp = np.zeros([len(sat_for_spacing_lower_plane_only), len(sat_for_spacing_no_low_plane), nb_time_step_in_one_month]) # when a sat in the lower plane flies below a sat in the upper plane (within +- 1 degree)
    conjunction_down_up = []
    conjunction_down_up_sub1list = []
    for k_count in range(len(sat_for_spacing_lower_plane_only)):
        k = sat_for_spacing_lower_plane_only[k_count]
        conjunction_down_up_sub1list.append(k)
        conjunction_down_up_sub2list = []
        for j_count in range(len(sat_for_spacing_no_low_plane)):
            j = sat_for_spacing_no_low_plane[j_count]
            print k_count,j_count
            conjunction_down_up_sub2list.append(j)
            conjunction_down_up_sub3list = []
            for i in range(nb_time_step_in_one_month):
                if ( ( np.abs( angle_asc_node_to_sat[j,(month_number_start-1)*nb_time_step_in_one_month+i] - angle_asc_node_to_sat[k,(month_number_start-1)*nb_time_step_in_one_month+i] ) < conjunction_angle ) | (  np.abs( angle_asc_node_to_sat[j,(month_number_start-1)*nb_time_step_in_one_month+i] - angle_asc_node_to_sat[k,(month_number_start-1)*nb_time_step_in_one_month+i] ) > (360 - conjunction_angle) ) ):
                    conjunction_down_up_temp[k_count,j_count,i] = i
                if (i > 0):
                    if ( (conjunction_down_up_temp[k_count,j_count,i-1] == 0) & (conjunction_down_up_temp[k_count,j_count,i] != 0) ):
                        conjunction_down_up_sub4list = []
                        conjunction_down_up_sub4list.append([i])
                    if ( (conjunction_down_up_temp[k_count,j_count,i-1] != 0) & (conjunction_down_up_temp[k_count,j_count,i] == 0) ):
                        conjunction_down_up_sub4list.append([i-1])
                        conjunction_down_up_sub3list.append(conjunction_down_up_sub4list)
            conjunction_down_up_sub2list.append(conjunction_down_up_sub3list)
        conjunction_down_up_sub1list.append(conjunction_down_up_sub2list)
    conjunction_down_up_list.append(conjunction_down_up_sub1list)
nb_conjunction_per_down_per_month_per_lower_per_upper_sat = np.zeros([nb_months, len(sat_for_spacing_lower_plane_only), len(sat_for_spacing_no_low_plane)])
nb_conjunction_per_down_per_month_per_lower_sat = np.zeros([nb_months, len(sat_for_spacing_lower_plane_only)])
total_nb_conjunction_per_down_per_lower_sat =  np.zeros([nb_months, len(sat_for_spacing_lower_plane_only)])
m_count = -1
for m in range(1, 2*nb_months,2):
    m_count = m_count + 1
    d_count = -1
    for d in range(1,2*len(sat_for_spacing_lower_plane_only),2):
        d_count = d_count + 1
        u_count = -1
        for u in range(1,2*len(sat_for_spacing_no_low_plane),2):
            u_count = u_count + 1
            nb_conjunction_per_down_per_month_per_lower_per_upper_sat[m_count, d_count, u_count] = len(conjunction_down_up_list[m][d][u])
        nb_conjunction_per_down_per_month_per_lower_sat[m_count, d_count] = sum(nb_conjunction_per_down_per_month_per_lower_per_upper_sat[m_count, d_count,:])

for m in range(1,(nb_months+1)):
    for d in range(len(sat_for_spacing_lower_plane_only)):
        total_nb_conjunction_per_down_per_lower_sat[m-1, d] = sum(nb_conjunction_per_down_per_month_per_lower_sat[0:m, d])
ax1.plot(np.arange(1,(nb_months+1),1), total_nb_conjunction_per_down_per_lower_sat[:,0], linewidth = 2, color = 'b', label = 'L1 conjunctions')
ax1.plot(np.arange(1,(nb_months+1),1), total_nb_conjunction_per_down_per_lower_sat[:,1], linewidth = 2, color = 'r', label = 'L2 conjunctions')
ax1.plot(np.arange(1,(nb_months+1),1), total_nb_conjunction_per_down_per_lower_sat[:,1]+total_nb_conjunction_per_down_per_lower_sat[:,0], linewidth = 2, color = 'k', label = 'Total conjunctions')
ax1.set_title('Conjunctions', weight = 'bold', fontsize = 20,  y = 1.008) 
ax1.set_ylabel('Total number of conjunctions', fontsize = 18, weight = 'bold')
ax1.set_xlabel('Time (months)', fontsize = 18, weight = 'bold')





# ######## RAAN
# ax1.set_title('RAAN difference', weight = 'bold', fontsize = 20,  y = 1.008) 
# ax1.set_ylabel('Delta RAAN ('+u'\N{DEGREE SIGN}'+')', fontsize = 18, weight = 'bold')
# ax1.set_xlabel('Real Time', fontsize = 18, weight = 'bold')
# # colors = ['b','b','r','b','b','b','b','r','k','k']
# # label_array = ['M1','M2','L1','M3','M4','I1','I2','L2','s9','s10']
# ax1.plot((raan[9,:] - raan[0,:])%360, color='b', label = '83'+u'\N{DEGREE SIGN}/500 km - 81'+u'\N{DEGREE SIGN}/500 km', linewidth = 2)
# ax1.plot((raan[2,:] - raan[0,:])%360, color='r', label = '81'+u'\N{DEGREE SIGN}/475 km - 81'+u'\N{DEGREE SIGN}/500 km', linewidth = 2)



######## SPACING_M1_TO_OTHER_SAT
#SPACING_M1_TO_OTHER_SAT
#colors = color_table(nb_satellites)
# colors = ['b','r','k','c','g','m','y','k']
# linewidth_array = [2,2,2,2,2,2,2,2]
# ax1.set_title('Angular distance to M1 VS time', weight = 'bold', fontsize = 20,  y = 1.008) 
# ax1.set_ylabel('Angular distance (degrees)', fontsize = 18, weight = 'bold')
# ax1.set_xlabel('Real Time', fontsize = 18, weight = 'bold')
# sat_for_angle_distance = [1,2,3,4,5,6,7]
# label_array = []
# for i in sat_for_angle_distance:
#     label_array.append("M1 with s"+str(i+1))
#     print i
#     if ((i != 2) & (i != 7)):
#         ax1.plot(spacing_M1_to_other_sat[i,:], color=colors[i], label = label_array[i], linewidth = linewidth_array[i])
# #    if ((i != 2) & (i != 7)):

######### MIN AND MAX SPACING
# ax1.set_title('Min and max spacing time', weight = 'bold', fontsize = 20,  y = 1.008) 
# ax1.set_ylabel('Spacing Time (min)', fontsize = 18, weight = 'bold')
# ax1.set_xlabel('Real Time', fontsize = 18, weight = 'bold')
# ax1.plot(min_spacing_time, color='b', label = 'Min', linewidth = 2)
# ax1.plot(max_spacing_time, color='r', label = 'Max', linewidth = 2)



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

