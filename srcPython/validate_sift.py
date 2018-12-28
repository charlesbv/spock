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
# This script compares the psotion of the specular points prediced by SpOCK VS reported in the l1 netcdf files at cygnss-sftp-1:/data/cygnss/products/l1/2017
# Assumptions:
# - see section PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

from datetime import datetime, timedelta
import numpy as np
import os
import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/kalman/spock_development_new_structure_kalman_dev/srcPython")
#sys.path.append("/home/cbv/spock_development_new_structure_kalman_dev/srcPython")
from os import listdir
from read_input_file import *
from cygnss_read_spock_spec import *
from netCDF4 import Dataset
import numpy.ma as ma

import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from matplotlib import pyplot as plt
import pickle

cygfm = 1 # which CYGNSS to look at
downlowd_netcdf =1 # set this variable to 1 if the entcdf files have not been download yet for the interval of time specified by [date_start_val, date_stop_val]

# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
date_start_val_start = "2017-06-01T00:00:00"
date_start_val_start = datetime.strptime(date_start_val_start, "%Y-%m-%dT%H:%M:%S")
date_start_val_array = np.array([date_start_val_start + timedelta(days=i) for i in np.arange(1,30*5, 7)])
nb_date = len(date_start_val_array)
# for idate in range(nb_date):
#     date_start_val = datetime.strftime(date_start_val_array[idate], "%Y-%m-%dT%H:%M:%S")#"2017-08-20T00:00:00" # date to start the validation
#     date_stop_val = datetime.strftime(date_start_val_array[idate] + timedelta(seconds = 6*24*3600+23*3600+59*60+59), "%Y-%m-%dT%H:%M:%S")#"2017-08-20T06:00:00" # date to stop the validation

#     # date_start_val = "2017-08-20T00:00:00" # date to start the validation
#     # date_stop_val = "2017-08-26T23:59:59" # date to stop the validation
#     # end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

#     # Run SpOCK from the date_start_val to date_stop_val
#     os.chdir("/Users/cbv/cygnss/spock")
#     if idate != 5: #!!!!!!!!!!!!!!!!!!!! REMOVE THIS IF condition (leave the os.system right below though)
#         print "spock_cygnss_spec_parallel.py " + date_start_val + " " + date_stop_val + " spec"
#         os.system("spock_cygnss_spec_parallel.py " + date_start_val + " " + date_stop_val + " spec")
#     os.chdir("/Users/cbv/Google Drive/Work/PhD/Research/Code/cygnss/validate_sift")
#     print "Done running SpOCK"
#     # Read specular positio computed by SpOCK
#     spock_input_filename = "/Users/cbv/cygnss/spock/spock_spec_start_" + date_start_val.replace(":", "_") + "_end_" + date_stop_val.replace(":", "_") + ".txt" # spock_spec_start_2017-08-20T00_00_00_end_2017-08-20T06_00_00.txt
#     var_in, var_in_order = read_input_file(spock_input_filename)
#     output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
#     output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')];
#     cygfm_to_spock_nb = [4,3,8,2,1,6,7,5] # ['41884', '41885', '41886', '41887', '41888', '41889', '41890', '41891']
#     isc =  cygfm_to_spock_nb[cygfm-1] - 1
#     spec_spock_filename = "/Users/cbv/cygnss/spock/" + output_file_path_list[isc] + "specular_" + output_file_name_list[isc]
#     date_spock, lon_spock, lat_spock, gain_spock, gps_spock, normpower_spock, x_cyg_spock, y_cyg_spock, z_cyg_spock, x_gps_spock, y_gps_spock, z_gps_spock,  x_spec_spock, y_spec_spock, z_spec_spock = cygnss_read_spock_spec(spec_spock_filename)

#     # Download netcdf files from sftp-1
#     date_start_val_doy = (int)(datetime.strptime( date_start_val, "%Y-%m-%dT%H:%M:%S" ).strftime('%j'))
#     date_stop_val_doy = (int)(datetime.strptime( date_stop_val, "%Y-%m-%dT%H:%M:%S" ).strftime('%j'))
#     if datetime.strptime( date_stop_val, "%Y-%m-%dT%H:%M:%S" ).strftime('%Y') != '2017':
#         print "***!The analysis has to be for data in 2017. The program will stop. !***"; raise Exception
#     doy_array = np.arange(date_start_val_doy, date_stop_val_doy+1, 1)
#     nb_day = len(doy_array)
#     day_remove = [] # list of days to remove from the analysis
#     day_list = []
#     for iday in range(nb_day):
#         if (os.path.isdir("/Users/cbv/cygnss/netcdf/" + str(doy_array[iday]).zfill(3)) == False):
#             os.system("mkdir /Users/cbv/cygnss/netcdf/" + str(doy_array[iday]).zfill(3))
#         if downlowd_netcdf == 1:
#             os.system("scp -p cygnss-sftp-1.engin.umich.edu:/data/cygnss/products/l1/2017/" + str(doy_array[iday]).zfill(3) + "/cyg0" + str(cygfm) + "* /Users/cbv/cygnss/netcdf/" + str(doy_array[iday]).zfill(3))
#         if len([x for x in os.listdir("/Users/cbv/cygnss/netcdf/" + str(doy_array[iday]).zfill(3)) if x.endswith(".nc") and x.startswith("cyg0" + str(cygfm))]) != 1: # if more than one file for this sc then don't take this day into account in the analsyis OR if no netcdf
#             day_remove.append(iday)
#         else:
#             day_list.append(iday)
            
#     nb_day = len(day_list)

#     ################
#     ################
#     ################
#     # For each day, read the specular point position from the netcdf file
#     x_spec_netcdf = []
#     y_spec_netcdf = []
#     z_spec_netcdf = []
#     x_cyg_netcdf = []
#     y_cyg_netcdf = []
#     z_cyg_netcdf = []
#     vx_cyg_netcdf = []
#     vy_cyg_netcdf = []
#     vz_cyg_netcdf = []
#     gain_netcdf = []
#     date_flight = []
#     date_flight_rounded = []
#     index_in_spock_date_netcdf_same = [] 
#     nb_seconds_since_initial_epoch_spock = []
#     iday_count = -1
#     while iday_count < nb_day-1:
#         iday_count = iday_count + 1
#         iday_here = day_list[iday_count]
#         filename_spec_flight = "/Users/cbv/cygnss/netcdf/" + str(doy_array[iday_here]).zfill(3) + "/" + [x for x in os.listdir("/Users/cbv/cygnss/netcdf/" + str(doy_array[iday_here]).zfill(3)) if x.endswith(".nc") and x.startswith("cyg0" + str(cygfm))][0]
#         fh = Dataset(filename_spec_flight, mode='r')
#         # nc_attrs = fh.ncattrs()
#         # nc_dims = [dim for dim in fh.dimensions]  # list of nc dimensions
#         # nc_vars = [var for var in fh.variables]  # list of nc variables

#         x_spec_netcdf_temp = fh.variables['sp_pos_x'][:] # X component of the specular point position in the ECEF coordinate system, in meters, at ddm_timestamp_utc, as calculated on the ground.
#         y_spec_netcdf_temp = fh.variables['sp_pos_y'][:]
#         z_spec_netcdf_temp = fh.variables['sp_pos_z'][:]
#         x_cyg_netcdf_temp = fh.variables['sc_pos_x'][:]
#         y_cyg_netcdf_temp = fh.variables['sc_pos_y'][:]
#         z_cyg_netcdf_temp= fh.variables['sc_pos_z'][:]
#         gain_netcdf_temp = fh.variables['sp_rx_gain'][:] # The receive antenna gain in the direction of the specular point, in dBi, at ddm_timestamp_utc
#         vx_cyg_netcdf_temp = fh.variables['sc_vel_x'][:]
#         vy_cyg_netcdf_temp = fh.variables['sc_vel_y'][:]
#         vz_cyg_netcdf_temp= fh.variables['sc_vel_z'][:]


#         list_are_masked_array = [] # sometimes the netcdf varaible below are maksed array and soemtimes they are not (depending on which netcdf file)...
#         if type(x_spec_netcdf_temp) == ma.core.MaskedArray:
#             list_are_masked_array.append(x_spec_netcdf_temp)
#         if type(y_spec_netcdf_temp) == ma.core.MaskedArray:
#             list_are_masked_array.append(y_spec_netcdf_temp)
#         if type(z_spec_netcdf_temp) == ma.core.MaskedArray:
#             list_are_masked_array.append(z_spec_netcdf_temp)
#         if type(gain_netcdf_temp) == ma.core.MaskedArray:
#             list_are_masked_array.append(gain_netcdf_temp)
#         nb_mask_array = len(list_are_masked_array)

#         time_flight = fh.variables['ddm_timestamp_utc'][:]
#         time_coverage_start = fh.getncattr(fh.ncattrs()[fh.ncattrs().index('time_coverage_start')])
#         time_coverage_start_datetime = datetime.strptime(time_coverage_start[:-4], "%Y-%m-%dT%H:%M:%S.%f") 
#         fh.close()
#         nb_time_flight_temp = len(x_cyg_netcdf_temp)
#         date_flight_t = []
#         date_flight_rounded_temp = []
#         time_remove_list = []
#         itime = -1
#         date_flight_raw = []
#         while itime < nb_time_flight_temp - 60:
#             itime = itime + 60
#             print itime, nb_time_flight_temp-1, iday_count, nb_day-1
#             time_remove = 0
#             date_flight_temp = datetime.strftime(time_coverage_start_datetime + timedelta(microseconds = round(time_flight[itime]*10**6)), "%Y-%m-%dT%H:%M:%S.%f" )
#             date_flight_raw.append(date_flight_temp)
#             imask_arr = 0
#             while (imask_arr < nb_mask_array):
#                         #if ( ( True in x_spec_netcdf_temp.mask[itime] ) | ( True in y_spec_netcdf_temp.mask[itime] ) | ( True in z_spec_netcdf_temp.mask[itime] ) | ( True in gain_netcdf_temp.mask[itime] ) ):# id one of the 4 spec is masked then ignore this time
#                 if (  True in list_are_masked_array[imask_arr].mask[itime] ):
#                      time_remove = 1
#                 imask_arr = imask_arr + 1     
#             # round to neared second but only if the actual time is less than 100 ms from the nearest second, otherwise ignore this time (don't compare to SpOCK).This is because SpOCK propagates with 1 s time step. So to compare to the netcdf file, we assume that the netcdf is also exactly at each second (so no millisecond). 100 ms wrong is not too bad because the satellite movesby less than 1 km.
#             if ( date_flight_temp.split('.')[1][0] == '9' ): # round to next second
#                 date_flight_temp_date = datetime.strptime(date_flight_temp, "%Y-%m-%dT%H:%M:%S.%f")
#                 date_flight_date = date_flight_temp_date + timedelta(seconds = 1)
#                 date_flight_date_rounded_temp = datetime.strftime(date_flight_date, "%Y-%m-%dT%H:%M:%S.%f").split('.')[0]
#                 date_flight_date_rounded = datetime.strptime(date_flight_date_rounded_temp, "%Y-%m-%dT%H:%M:%S")
#             elif ( date_flight_temp.split('.')[1][0] == '0' ): # round to next second
#                 date_flight_date_rounded = datetime.strptime(date_flight_temp.split('.')[0], "%Y-%m-%dT%H:%M:%S" )
#             else: #if time can't be rounded by less than 100 ms
#                 time_remove = 1
#             if ( (time_remove == 0) & (date_flight_date_rounded in date_spock)): # if this time is not in date_spock
#                 index_in_spock_date_netcdf_same.append(date_spock.index(date_flight_date_rounded))
#             else:
#                 time_remove = 1
#             if ( time_remove == 1 ): # remove time if can't be rounded by ess than 100 ms or if is not in date_spock or if masked
#                 time_remove_list.append(itime)
#             else:
#                 if type(x_spec_netcdf_temp) == ma.core.MaskedArray:
#                     x_spec_netcdf.append(x_spec_netcdf_temp.data[itime]/1000.)
#                 else:
#                     x_spec_netcdf.append(x_spec_netcdf_temp[itime]/1000.)
#                 if type(y_spec_netcdf_temp) == ma.core.MaskedArray:
#                     y_spec_netcdf.append(y_spec_netcdf_temp.data[itime]/1000.)
#                 else:
#                     y_spec_netcdf.append(y_spec_netcdf_temp[itime]/1000.)
#                 if type(z_spec_netcdf_temp) == ma.core.MaskedArray:
#                     z_spec_netcdf.append(z_spec_netcdf_temp.data[itime]/1000.)
#                 else:
#                     z_spec_netcdf.append(z_spec_netcdf_temp[itime]/1000.)
#                 if type(gain_netcdf_temp) == ma.core.MaskedArray:
#                     gain_netcdf.append(gain_netcdf_temp.data[itime])
#                 else:
#                     gain_netcdf.append(gain_netcdf_temp[itime])
#                 x_cyg_netcdf.append(x_cyg_netcdf_temp[itime]/1000.)
#                 y_cyg_netcdf.append(y_cyg_netcdf_temp[itime]/1000.)
#                 z_cyg_netcdf.append(z_cyg_netcdf_temp[itime]/1000.)
#                 date_flight_rounded.append(date_flight_date_rounded)
#                 date_flight.append(date_flight_temp)
#                 nb_seconds_since_initial_epoch_spock.append( ( date_flight_date_rounded - date_spock[0] ).total_seconds() )

#                 vx_cyg_netcdf.append(vx_cyg_netcdf_temp[itime]/1000.)
#                 vy_cyg_netcdf.append(vy_cyg_netcdf_temp[itime]/1000.)
#                 vz_cyg_netcdf.append(vz_cyg_netcdf_temp[itime]/1000.)

#             if date_flight_date_rounded > date_spock[-1]:
#                 iday_count = nb_day + 1
#                 break



#     nb_time_flight_rounded = len(date_flight_rounded)
#     date_netcdf = np.array( date_flight_rounded )
#     date_spock = np.array(date_spock)
#     nb_seconds_since_initial_epoch_spock = np.array(nb_seconds_since_initial_epoch_spock)

#     # Compare SpOCK and netcdf
#     ## Select the dates that are both in SpOCK and netcdf
#     nb_time_spock_netcdf = len(date_flight) # number of time netcdf same as spock
#     x_cyg_spock = np.array(x_cyg_spock); y_cyg_spock = np.array(y_cyg_spock); z_cyg_spock = np.array(z_cyg_spock)
#     x_spec_spock = np.array(x_spec_spock); y_spec_spock = np.array(y_spec_spock); z_spec_spock = np.array(z_spec_spock)
#     gain_spock = np.array(gain_spock)
#     x_cyg_spock_same_time_as_netcdf = x_cyg_spock[index_in_spock_date_netcdf_same]
#     y_cyg_spock_same_time_as_netcdf = y_cyg_spock[index_in_spock_date_netcdf_same]
#     z_cyg_spock_same_time_as_netcdf = z_cyg_spock[index_in_spock_date_netcdf_same]
#     x_spec_spock_same_time_as_netcdf = x_spec_spock[index_in_spock_date_netcdf_same]
#     y_spec_spock_same_time_as_netcdf = y_spec_spock[index_in_spock_date_netcdf_same]
#     z_spec_spock_same_time_as_netcdf = z_spec_spock[index_in_spock_date_netcdf_same]
#     gain_spock_same_time_as_netcdf = gain_spock[index_in_spock_date_netcdf_same]
#     ## Calculate difference in position of sc and spec
#     r_cyg_diff_mag = np.zeros([nb_time_spock_netcdf])
#     r_spec_diff_mag_max_gain = np.zeros([nb_time_spock_netcdf])
#     r_spec_diff_mag_min_diff = np.zeros([nb_time_spock_netcdf])
#     gain_netcdf_that_min_error = np.zeros([nb_time_spock_netcdf])
#     gain_spock_that_min_error = np.zeros([nb_time_spock_netcdf])
#     which_comb_min = []
#     r_spec_diff_mag_min_to_max = []
#     r_spec_diff_mag_min_to_max_array = np.zeros([nb_time_spock_netcdf, 4]) - 1 # will be -1 if spock predicts less than 4 spec (if netcdf has less than 4 points than this time is ignored (see mask filter previously in the code))
#     r_cyg_spock = np.zeros([nb_time_spock_netcdf, 3])
#     r_cyg_netcdf = np.zeros([nb_time_spock_netcdf, 3])
#     for itime in range(nb_time_spock_netcdf):
#         # difference sc position
#         r_cyg_spock[itime,:] = np.array([x_cyg_spock_same_time_as_netcdf[itime][0], y_cyg_spock_same_time_as_netcdf[itime][0], z_cyg_spock_same_time_as_netcdf[itime][0]])
#         r_cyg_netcdf[itime,:] = np.array([x_cyg_netcdf[itime], y_cyg_netcdf[itime], z_cyg_netcdf[itime]])
#         r_cyg_diff = r_cyg_spock[itime,:] - r_cyg_netcdf[itime,:]
#         r_cyg_diff_mag[itime] = np.linalg.norm(r_cyg_diff)

#         # difference spec position
#         nb_spec_spock = len(x_spec_spock_same_time_as_netcdf[itime])
#         nb_spec_netcdf = len(x_spec_netcdf[itime])
#         r_spec_diff_mag_temp = []
#         r_spec_diff_mag_temp_all = []
#         where_max_gain_spock = np.where(gain_spock_same_time_as_netcdf[itime] == np.max(gain_spock_same_time_as_netcdf[itime]))[0][0]
#         where_max_gain_netcdf = np.where(gain_netcdf[itime] == np.max(gain_netcdf[itime]))[0][0]
#         r_spec_spock = np.array([x_spec_spock_same_time_as_netcdf[itime][where_max_gain_spock], y_spec_spock_same_time_as_netcdf[itime][where_max_gain_spock], z_spec_spock_same_time_as_netcdf[itime][where_max_gain_spock]])
#         r_spec_netcdf = np.array([x_spec_netcdf[itime][where_max_gain_netcdf], y_spec_netcdf[itime][where_max_gain_netcdf], z_spec_netcdf[itime][where_max_gain_netcdf]])
#         r_spec_diff = r_spec_spock - r_spec_netcdf
#         r_spec_diff_mag_max_gain[itime] =  np.linalg.norm(r_spec_diff)
#         dist_all_spec_to_all_spec = np.zeros([nb_spec_spock, nb_spec_netcdf])
#         for ispec_spock in range(nb_spec_spock):
#             r_spec_spock = np.array([x_spec_spock_same_time_as_netcdf[itime][ispec_spock], y_spec_spock_same_time_as_netcdf[itime][ispec_spock], z_spec_spock_same_time_as_netcdf[itime][ispec_spock]])
#             r_spec_diff_mag_temp_per_spock = []
#             for ispec_netcdf in range(nb_spec_netcdf):
#                 r_spec_netcdf = np.array([x_spec_netcdf[itime][ispec_netcdf], y_spec_netcdf[itime][ispec_netcdf], z_spec_netcdf[itime][ispec_netcdf]])
#                 r_spec_diff = r_spec_spock - r_spec_netcdf
#                 dist_all_spec_to_all_spec[ispec_spock, ispec_netcdf] = np.linalg.norm(r_spec_diff)
#                 r_spec_diff_mag_temp_per_spock.append( np.linalg.norm(r_spec_diff) )
#                 r_spec_diff_mag_temp_all.append(np.linalg.norm(r_spec_diff))
#             r_spec_diff_mag_temp.append(r_spec_diff_mag_temp_per_spock)
#         r_spec_diff_mag_min_diff[itime] = np.min(r_spec_diff_mag_temp_all)
#         # Take the min distance of the 4*4 distance array new_dist_all_spec_to_all_spec. Save the min distance, remove the correspodning spock and netcdf spec (they are the ones that min the difference) from new_dist_all_spec_to_all_spec. Now new_dist_all_spec_to_all_spec is 3*3. Save the new min, and remove the corresponding spock and netcdf. Now new_dist_all_spec_to_all_spec is 2*2. Proceeed again twice. So r_spec_diff_mag_min_to_max_temp shows 4 differences, from min to max
#         new_dist_all_spec_to_all_spec = dist_all_spec_to_all_spec
#         ispec_min = 0
#         r_spec_diff_mag_min_to_max_temp = []
#         #
#         while ispec_min < nb_spec_spock:
#             comb_min = np.unravel_index(new_dist_all_spec_to_all_spec.argmin(), new_dist_all_spec_to_all_spec.shape)
#             r_spec_diff_mag_min_to_max_temp.append( new_dist_all_spec_to_all_spec[comb_min] )
#             r_spec_diff_mag_min_to_max_array[itime, ispec_min] = new_dist_all_spec_to_all_spec[comb_min]
#             new_dist_all_spec_to_all_spec = np.delete(np.delete(new_dist_all_spec_to_all_spec, comb_min[0],0), comb_min[1],1)
#             ispec_min = ispec_min + 1
#         r_spec_diff_mag_min_to_max.append(r_spec_diff_mag_min_to_max_temp)

#         # Look at which spock spec and which netcdf spec minimize the difference
#         which_spec_spock_min = comb_min[0] # # which spock spec that minimizes the distance to netdf
#         which_spec_netcdf_min = comb_min[1] # which netcdf spec that minimizes the distance to spock
#         which_comb_min.append( [which_spec_spock_min, which_spec_netcdf_min] ) # combination [spock_spec, netcdf_spec] that minimizes the difference in spec position
#         gain_netcdf_that_min_error[itime] = gain_netcdf[itime][which_spec_netcdf_min] # netcdf gain of the spec that minimizes the distance to spock
#         gain_spock_that_min_error[itime] = gain_spock_same_time_as_netcdf[itime][which_spec_spock_min] # spock gain of the spec that minimizes the distance to netcdf





#     pickle.dump(nb_seconds_since_initial_epoch_spock, open("/Users/cbv/cygnss/pickle/" + date_start_val[:10] + "_nb_seconds_since_initial_epoch_spock.pickle", "w"))
#     pickle.dump(r_spec_diff_mag_min_to_max_array, open("/Users/cbv/cygnss/pickle/" + date_start_val[:10] + "_r_spec_diff_mag_min_to_max_array.pickle", "w"))
#     pickle.dump(r_cyg_diff_mag, open("/Users/cbv/cygnss/pickle/" + date_start_val[:10] + "_r_cyg_diff_mag.pickle", "w"))
#     pickle.dump(r_cyg_spock, open("/Users/cbv/cygnss/pickle/" + date_start_val[:10] + "_r_cyg_spock.pickle", "w"))
#     pickle.dump(r_cyg_netcdf, open("/Users/cbv/cygnss/pickle/" + date_start_val[:10] + "_r_cyg_netcdf.pickle", "w"))
#     print date_start_val
# raise Exception

nb_seconds_since_initial_epoch_spock_all = []
r_spec_diff_mag_min_to_max_array_all = []
r_cyg_diff_mag_all = []
r_cyg_spock_all = []
r_cyg_netcdf_all = []
nb_time_spock_netcdf_all = []
for idate in range(nb_date):
    print idate, nb_date-1
    date_start_val = datetime.strftime(date_start_val_array[idate], "%Y-%m-%dT%H:%M:%S")#"2017-08-20T00:00:00" # date to start the validation

    nb_seconds_since_initial_epoch_spock_all.append( pickle.load(open("./pickle/" + date_start_val[:10] + "_nb_seconds_since_initial_epoch_spock.pickle")))
    r_spec_diff_mag_min_to_max_array_all.append( pickle.load(open("./pickle/" + date_start_val[:10] + "_r_spec_diff_mag_min_to_max_array.pickle")))
    r_cyg_diff_mag_all.append( pickle.load(open("./pickle/" + date_start_val[:10] + "_r_cyg_diff_mag.pickle")))
    r_cyg_spock_all.append( pickle.load(open("./pickle/" + date_start_val[:10] + "_r_cyg_spock.pickle")))
    r_cyg_netcdf_all.append( pickle.load(open("./pickle/" + date_start_val[:10] + "_r_cyg_netcdf.pickle")))
    nb_time_spock_netcdf_all.append( len(r_cyg_netcdf_all[-1]))


# nb_seconds_since_initial_epoch_spock = pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val[:10] + "_nb_seconds_since_initial_epoch_spock.pickle"))
# r_spec_diff_mag_min_to_max_array = pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val[:10] + "_r_spec_diff_mag_min_to_max_array.pickle"))
# r_cyg_diff_mag = pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val[:10] + "_r_cyg_diff_mag.pickle"))
# r_cyg_spock = pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val[:10] + "_r_cyg_spock.pickle"))
# r_cyg_netcdf = pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val[:10] + "_r_cyg_netcdf.pickle"))
# nb_time_spock_netcdf = len(r_cyg_netcdf)

height_fig = 11.  # the width is calculated as height_fig * 4/3.                                                                                                                                                                                                                                                                                                                           
fontsize_plot = 20
ratio_fig_size = 4./3

fig_title = 'Difference in specular point position SpOCK vs L1' # L1 is netcdf                                                                                                                                                                                                                                         
y_label = 'Difference (km)'
x_label = 'Prediction time (days)'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold                                                                                                                                                                                                                                                                                                                      
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                                                                                                                                                                                                                                                                         
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold                           

for idate in range(nb_date):
    r_spec_diff_mag_min_to_max_array = r_spec_diff_mag_min_to_max_array_all[idate]
    if r_spec_diff_mag_min_to_max_array[-1,0] < 60:
        nb_seconds_since_initial_epoch_spock = nb_seconds_since_initial_epoch_spock_all[idate]
        x_axis = nb_seconds_since_initial_epoch_spock / 3600. / 24.
        ax.plot(x_axis, r_spec_diff_mag_min_to_max_array[:,0], linewidth = 2, color = 'b', label = '1')
        print idate,  r_spec_diff_mag_min_to_max_array[-1,0]
    #ax.plot(x_axis, r_spec_diff_mag_min_to_max_array[:,1], linewidth = 2, color = 'r', label = '2')
    #ax.plot(x_axis, r_spec_diff_mag_min_to_max_array[:,2], linewidth = 2, color = 'k', label = '3')
    #ax.plot(x_axis, r_spec_diff_mag_min_to_max_array[:,3], linewidth = 2, color = 'g', label = '4')                                                                                                                                                                                                                                                                                           

ax.margins(0,0)
ax.set_ylim([0,70])
#ax.set_xlim([0,5])

# legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
# legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = 'all_dates_difference_spec_position_spock_vs_netcdf.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')



# sc position difference
fig_title = 'Difference in CYGNSS position SpOCK vs L1' # L1 is netcdf                                                                                                                                                                                                                                         
y_label = 'Difference (km)'
x_label = 'Prediction time (days)'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold                                                                                                                                                                                                                                                                                                                      
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                                                                                                                                                                                                                                                                         
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold                                                                                                                                                                                                                                                                                                                     
x_axis = nb_seconds_since_initial_epoch_spock / 3600. / 24
ax.plot(x_axis, r_cyg_diff_mag, linewidth = 2, color = 'b', label = '1')
ax.margins(0,0)
ax.set_ylim([0,80])

# legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
# legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = date_start_val[:10] + '_difference_sat_position_spock_vs_netcdf.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')




# radus difference
r_cyg_spock_mag = np.zeros([nb_time_spock_netcdf])
r_cyg_netcdf_mag = np.zeros([nb_time_spock_netcdf])
for itime in range(nb_time_spock_netcdf):
    r_cyg_spock_mag[itime] = np.linalg.norm(r_cyg_spock[itime, :])
    r_cyg_netcdf_mag[itime] = np.linalg.norm(r_cyg_netcdf[itime, :])


fig_title = 'Difference in CYGNSS radius SpOCK vs L1' # L1 is netcdf
y_label = 'Difference (km)'
x_label = 'Prediction time (days)'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
x_axis = nb_seconds_since_initial_epoch_spock / 3600. / 24
ax.plot(x_axis, r_cyg_spock_mag - r_cyg_netcdf_mag, linewidth = 2, color = 'b', label = '1')
ax.margins(0,0)
#ax.set_ylim([0,80])

# legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
# legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = date_start_val[:10] + '_difference_sat_radius_spock_vs_netcdf.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')



