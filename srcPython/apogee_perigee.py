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
import numpy as np

#rmag = np.sqrt( x_eci_ensemble[:,j]**2 + y_eci_ensemble[:,j]**2 + z_eci_ensemble[:,j]**2 )
rmag = np.sqrt( pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2 )
index_start_orbit = index_time_orbit_average[:,0]
index_stop_orbit = index_time_orbit_average[:,2]
nb_orbit = len(index_start_orbit)
r_a = np.zeros([nb_orbit])
r_p = np.zeros([nb_orbit])
sma_ok = np.zeros([nb_orbit])
for iorbit in range(nb_orbit):
    r_a[iorbit] = np.max( rmag[index_start_orbit[iorbit]: index_stop_orbit[iorbit] + 1] )
    r_p[iorbit] = np.min( rmag[index_start_orbit[iorbit]: index_stop_orbit[iorbit] + 1] )
    sma_ok[iorbit] = ( r_a[iorbit] + r_p[iorbit] ) / 2.
# def orbit_average(x, y, z, latitude, time_before_average): # x, y, and z in ECI
#     # isat = 0
#     # var_to_average = power[isat, :, 0]
#     # time_before_average = date
#     #latitude = latitude[isat, :] 
#     n = len(var_to_average)
#     var_orbit_averaged = []
#     time_averaged = []
#     index_time_averaged = []
#     istep = 0
#     passed_first_ascending_node = 0
#     while istep < n:
#         var_orbit_averaged_temp = []
#         time_averaged_temp = []
#         index_time_averaged_temp = []
#         while ( ( latitude[istep ] > 0 ) | ( latitude[istep  + 1] <= 0 ) ):
#             # start average at the first ascending node (ignore the steps before crossing the ascending node for the first time)
#             if ( passed_first_ascending_node == 1 ):
#                 var_orbit_averaged_temp.append( var_to_average[ istep ] )
#                 if len(var_orbit_averaged_temp) == 1:
#                     time_averaged_temp.append( time_before_average[istep] )
#                     index_time_averaged_temp.append( istep )
#                     istep_beginning = istep
#             istep = istep + 1
#             if istep + 1 >= n:
#                 break
#         if ( passed_first_ascending_node == 1 ):
#             var_orbit_averaged_temp.append( var_to_average[ istep ] )
#             time_averaged_temp.append( time_before_average[(istep_beginning + istep) / 2] )
#             time_averaged_temp.append( time_before_average[istep] )
#             index_time_averaged_temp.append( (istep_beginning + istep) / 2 )
#             index_time_averaged_temp.append( istep )
#             var_orbit_averaged.append( np.mean( var_orbit_averaged_temp ) )
#             time_averaged.append( time_averaged_temp )
#             index_time_averaged.append( index_time_averaged_temp )
#         passed_first_ascending_node = 1
#         istep = istep + 1
#         if istep + 1 >= n:
#             break

#     # Delete the last element because it likely corresponds to an average over only a part of the orbit
#     del time_averaged[-1]
#     del var_orbit_averaged[-1]
#     del index_time_averaged[-1]

#     return var_orbit_averaged, time_averaged, index_time_averaged # time_averaged[0] start of the orbit over which the average was made
#                                                                  # time_averaged[1] middle of the orbit over which the average was made
#                                                                  # time_averaged[2] end of the orbit over which the average was made
#                                                                  # index_time_averaged index in time_before_average of the orbit over which the average was made
