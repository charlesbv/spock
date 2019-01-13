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

# Important: 
# - the presentation "SIFT_accuracy3.pdf" (SOC meeting on May 15 2018) has the final results. The statistics and plots that were made for this presentation are between blocks "PRESENTATION 051518". These are the (only) statisctics and plots you should comput to validate futures versions of SIFT. The powerpoint presentation "SIFT_accuracy3.pptx" includes the names of the variables in the section notes (below the slides) used to compute these statistics.
# - there's a lot of code that I ended up not using anymore. I didn't delete it, though. The code that I've ended up *NOT* using to assess SIFT accuracy is in blocks labeled "IGNORE". Don't delete or comment these blocks because some variables defined in them might be used in other blocks. It just means that the variables calculated in these blcoks were used before to evaluate the accuracy of SIFT but not anymore
# - the code up to the line "break #######################################################################################################################################" (~line 510) is to read SpOCK and netcdf files and select the steps when  SpOCK and netcdf sample times match
# - the code after this line is to compare SpOCK and netcdf (position, PRN, gain, etc)
# - a long time ago, norm_power and gain were two different numbers in SpOCK. With the new antenna gain pattern file (the one that's on-board), they are now the same. Moreover, the word "FOM" (figure of merit) is used interchangeably with the word gain and the word normpower. THey represent the same variable
# - this code is a mess but it has been extensively tested, so it can be trusted (which doesn't mean there's absolutely no error, though!)

# Assumptions:
# - see section PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
import ipdb
from datetime import datetime, timedelta
import numpy as np
import os
import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/srcPython")
#sys.path.append("/home/cbv/spock_development_new_structure_kalman_dev/srcPython")
from os import listdir
from read_input_file import *
from read_output_file import *
from cygnss_read_spock_spec import *
from cygnss_read_spock_spec_bin import *
from netCDF4 import Dataset
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from matplotlib import pyplot as plt
#from ecef2eci import *
#from eci_to_lvlh import *
from ecef_to_lvlh import *
import pickle
from cygnss_name_to_norad_id import *
import os.path

# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
cygfm = 4 #1 # which CYGNSS to look at
download_netcdf = 0 # set this variable to 1 if the entcdf files have not been download yet for the interval of time specified by [date_start_val, date_stop_val]
date_start_val_start = '2018-08-31T00:00:00'#'2018-03-20T00:00:00'#"2017-06-01T00:00:00"
dir_run_spock = '/Users/cbv/cygnss/valsift_temp'#'/Users/cbv/cygnss/sift_temp' # '.' # no slash
dir_run_spock_out = '/Users/cbv/cygnss/valsift_temp/spock_out' #'spock' # name of output folder, no slash
pitch_max = 2. # filter out when pitch is greater than this value (in magnitude)
roll_max = 2. # filter out when roll is greater than this value (in magnitude)
yaw_max = 2. # filter out when yaw is greater than this value (in magnitude)
load_pickle = 0 # set to 1 if results have been previously saved in a pickle so no computation is made here
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
date_start_val_start = datetime.strptime(date_start_val_start, "%Y-%m-%dT%H:%M:%S")
date_start_val_array =  np.array([date_start_val_start + timedelta(days=i) for i in np.arange(1,2, 1)])
nb_date = len(date_start_val_array)


correct_vs_fom_all_date= []
correct_vs_fom_percent_all_date= []
time_diff_prn_all_date= []
time_same_prn_all_date= []
which_spec_diff_all_date= []
gps_spock_all_date= []
fom_spock_all_date= []
gps_netcdf_all_date= []
fom_netcdf_all_date= []
index_netcdf_all_date= []
fom_netcdf_diff_prn_all_all_date= []
nb_spec_diff_all_date= []
nb_spec_diff_gain_not_0_all_date= []
r_spec_diff_lvlh_same_prn_all_date= []
nb_seconds_since_initial_epoch_spock_all_date= []

diff_fom_same_prn_vs_fom_time_all_date = []
az_spec_netcdf_all_date = []
el_spec_netcdf_all_date = []
order_gain_netcdf_diff_prn_all_date = []
time_diff_prn_gain_not_0_all_date = []
nb_spec_spock_got_matched_vs_fom_all_date = []
nb_spec_spock_vs_fom_all_date = []

for idate in range(0,nb_date):# !!!!!!! should be: nb_date):
    date_start_val = datetime.strftime(date_start_val_array[idate], "%Y-%m-%dT%H:%M:%S")#'2017-06-02T00:15:00' # !!! should be datetime.strftime(date_start_val_array[idate], "%Y-%m-%dT%H:%M:%S")#"2017-08-20T00:00:00" # date to start the validation
    date_stop_val = datetime.strftime(date_start_val_array[idate]+ timedelta(seconds = 3.0*24*3600+23*3600+59*60+59), "%Y-%m-%dT%H:%M:%S") #'2017-06-06T01:02:00'#'2017-06-03T00:45:00'#'2017-06-06T01:02:00'#'2017-06-02T00:56:00'#'2017-06-06T01:02:00'#'2017-06-02T00:43:00'  # !!!! shoul be + timedelta(seconds = 6*24*3600+23*3600+59*60+59), "%Y-%m-%dT%H:%M:%S")#"2017-08-20T06:00:00" # date to stop the validation
    print idate, nb_date-1, date_start_val, date_stop_val

#     date_start_val = '2018-04-13T00:00:00'#'2018-04-13T00:00:00'#'2017-06-03T00:34:00'#'2017-06-02T00:15:00' # !!! should be datetime.strftime(date_start_val_array[idate], "%Y-%m-%dT%H:%M:%S")#"2017-08-20T00:00:00" # date to start the validation
#     date_stop_val = '2018-04-16T23:59:59' #'2017-06-06T01:02:00'  #'2017-06-03T00:50:00' #'2017-06-06T01:02:00'#'2017-06-03T00:45:00'#'2017-06-06T01:02:00'#'2017-06-02T00:56:00'#'2017-06-06T01:02:00'#'2017-06-02T00:43:00'  # !!!! shoul be + timedelta(seconds = 6*24*3600+23*3600+59*60+59), "%Y-%m-%dT%H:%M:%S")#"2017-08-20T06:00:00" # date to stop the validation

    # date_start_val = "2017-08-20T00:00:00" # date to start the validation
    # date_stop_val = "2017-08-26T23:59:59" # date to stop the validation
    # end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
    spock_input_filename = dir_run_spock  + "/spock_spec_start_" + date_start_val.replace(":", "_") + "_end_" + date_stop_val.replace(":", "_") + ".txt" # spock_spec_start_2017-08-20T00_00_00_end_2017-08-20T06_00_00.txt
    if load_pickle != 1:
        # Run SpOCK from the date_start_val to date_stop_val
        if idate < 9: # !!!!!!!! remove this if condition (keep what's below the 'if ' line)
            os.chdir(dir_run_spock )
            #if idate != 5: #!!!!!!!!!!!!!!!!!!!! REMOVE THIS IF condition (leave the os.system right below though)
            print "spock_cygnss_spec_parallel.py " + date_start_val + " " + date_stop_val + " spec"
            os.system("spock_cygnss_spec_parallel.py " + date_start_val + " " + date_stop_val + " spec")
            os.chdir("/Users/cbv/Google Drive/Work/PhD/Research/Code/cygnss/validate_sift")

        print "Done running SpOCK"
        # Read specular positio computed by SpOCK
        var_in, var_in_order = read_input_file(spock_input_filename)
        dt_spock_output = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')]; 
        output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
        output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')];
        gps_name_list_spock = var_in[find_in_read_input_order_variables(var_in_order, 'gps_name')];
        cygfm_to_spock_nb = [4,3,8,2,1,6,7,5] # ['41884', '41885', '41886', '41887', '41888', '41889', '41890', '41891']
        isc =  cygfm_to_spock_nb[cygfm-1] - 1
        spec_spock_filename = dir_run_spock + "/" + output_file_path_list[isc] + "specular_" + output_file_name_list[isc]
        #spec_spock_filename = spec_spock_filename.replace(".txt", "_6SPs.txt")  
        print "Reading SpOCK specular files..."
#         print "cp -p " + spec_spock_filename + " " + spec_spock_filename.replace(".txt", "_6SPs.txt")
#         os.system("cp -p " + spec_spock_filename + " " + spec_spock_filename.replace(".txt", "_6SPs.txt"))
#        date_spock, lon_spock, lat_spock, gain_spock, gps_spock, normpower_spock, x_cyg_spock, y_cyg_spock, z_cyg_spock, x_gps_spock, y_gps_spock, z_gps_spock,  x_spec_spock, y_spec_spock, z_spec_spock, nb_spec_spock,  el_spec_spock, az_spec_spock, el_gps_from_cyg_spock,  el_spec_not_int_spock, az_spec_not_int_spock = cygnss_read_spock_spec(spec_spock_filename)

        data_spec = cygnss_read_spock_spec_bin(spec_spock_filename.replace('.txt','.bin'), gps_name_list_spock, dt_spock_output, 1) 
        date_spock = data_spec[0]; lon_spock = data_spec[1]; lat_spock = data_spec[2]; gain_spock = data_spec[3]; gps_spock = data_spec[4]; normpower_spock = data_spec[5]; x_cyg_spock = data_spec[6]; y_cyg_spock = data_spec[7]; z_cyg_spock = data_spec[8]; x_gps_spock = data_spec[9]; y_gps_spock = data_spec[10]; z_gps_spock = data_spec[11];  x_spec_spock = data_spec[12]; y_spec_spock = data_spec[13]; z_spec_spock = data_spec[14]; nb_spec_spock = data_spec[15];  el_spec_spock = data_spec[16]; az_spec_spock = data_spec[17]; el_gps_from_cyg_spock = data_spec[18];  el_spec_not_int_spock = data_spec[19]; az_spec_not_int_spock = data_spec[20]

        print "Done reading SpOCK specular files..."
        ###############
        ###############
        ##########
        # Rad the ECEF files. Although the satllite position is output in the specular file, the velcoity is not. When I write that, I've already run SpOCK (I'm commenting the spock_cygnss_spec_parallel.py line above). I don't want ot re run SpOCK to get the print the velcoity in the spec files so I read it from the ECEF file.
        var_out, var_out_order = read_output_file(dir_run_spock + "/" + output_file_path_list[isc] + output_file_name_list[isc], ["position_ecef","velocity_ecef","position_tle","velocity_tle"] )
        r_cyg_spock_ecef_file = var_out[find_in_read_input_order_variables(var_out_order, 'position_ecef')]
        v_cyg_spock_ecef_file = var_out[find_in_read_input_order_variables(var_out_order, 'velocity_ecef')]
        r_cyg_spock_eci_tle = var_out[find_in_read_input_order_variables(var_out_order, 'position_tle')]
        v_cyg_spock_eci_tle = var_out[find_in_read_input_order_variables(var_out_order, 'velocity_tle')]
#         date_spock_not_interpolated_temp = var_out[find_in_read_input_order_variables(var_out_order, 'date')] # date_spock is every second because it corresponds to the sampling time of the spec points. 
#         date_spock_not_interpolated = []
#         for ii in range(len(date_spock_not_interpolated_temp)):
#             date_here = date_spock_not_interpolated_temp[ii]
#             date_here_date = datetime.strptime(date_here[:19], "%Y/%m/%d %H:%M:%S") # !!!!!! don't look at microseconds since SpOCK outputs at best eevey second anyway (so microseeconds should be 0)
#             date_spock_not_interpolated.append(date_here_date)
        date_spock_not_interpolated = var_out[find_in_read_input_order_variables(var_out_order, 'date_datetime_round_sec')]  # !!!!! used to be block above but probmen with blcok above is that if 59.9999 then rounded to 59 instead of next minute
        date_spock_not_interpolated = list(date_spock_not_interpolated)
        #date_spock_not_interpolated is the date directly output by SpOCK so it's with touput time step chosen in the main input file

        # Download netcdf files from sftp-1
        date_start_val_doy = (int)((datetime.strptime( date_start_val, "%Y-%m-%dT%H:%M:%S" )).strftime('%j'))# - timedelta(days = 2)).strftime('%j')) # start two day earlier so you can compare 
        # to the position given by the TLE that was used to initialize SpOCK. TLEs are usually output every day but go 2 days back to make sure
        date_stop_val_doy = (int)(datetime.strptime( date_stop_val, "%Y-%m-%dT%H:%M:%S" ).strftime('%j'))
#         if datetime.strptime( date_stop_val, "%Y-%m-%dT%H:%M:%S" ).strftime('%Y') != '2017':
#             print "***!The analysis has to be for data in 2017. The program will stop. !***"; raise Exception
        doy_array = np.arange(date_start_val_doy, date_stop_val_doy+1, 1)
        nb_day = len(doy_array)
        day_remove = [] # list of days to remove from the analysis
        day_list = []
        if  datetime.strptime( date_stop_val, "%Y-%m-%dT%H:%M:%S" ).strftime('%Y') == '2018':
            yy = '2018/'
        else:
            yy = ''
        for iday in range(nb_day):
            if (os.path.isdir("/Users/cbv/cygnss/netcdf/" + yy  + str(doy_array[iday]).zfill(3)) == False):
                os.system("mkdir /Users/cbv/cygnss/netcdf/" + yy  + str(doy_array[iday]).zfill(3))
            if download_netcdf == 1:
                os.system("scp -p cygnss-sftp-1.engin.umich.edu:/data/cygnss/products/l1/" + datetime.strptime( date_stop_val, "%Y-%m-%dT%H:%M:%S" ).strftime('%Y') + " /" + str(doy_array[iday]).zfill(3) + "/cyg0" + str(cygfm) + "* /Users/cbv/cygnss/netcdf/" + yy  + str(doy_array[iday]).zfill(3))
            if len([x for x in os.listdir("/Users/cbv/cygnss/netcdf/" + yy  + str(doy_array[iday]).zfill(3)) if x.endswith(".nc") and x.startswith("cyg0" + str(cygfm))]) != 1: # if more than one file for this sc then don't take this day into account in the analsyis OR if no netcd
                day_remove.append(iday)
            else:
                day_list.append(iday)

        nb_day = len(day_list)

        ################
        ################
        ################
        # For each day, read the specular point position from the netcdf file
        time_attitude_bad = []
        time_gain_0 = []
        x_spec_netcdf = []
        y_spec_netcdf = []
        z_spec_netcdf = []
        lat_spec_netcdf = [] 
        lon_spec_netcdf = [] 

        x_cyg_netcdf = []
        y_cyg_netcdf = []
        z_cyg_netcdf = []

        pitch_cyg_netcdf = []
        roll_cyg_netcdf = []
        yaw_cyg_netcdf = []

        x_gps_netcdf = []
        y_gps_netcdf = []
        z_gps_netcdf = []

        vx_cyg_netcdf = []
        vy_cyg_netcdf = []
        vz_cyg_netcdf = []

        x_cyg_netcdf_dt_output = [] # this is to compare to the r/v output in ECEF_ output SpOCK files
        y_cyg_netcdf_dt_output = []
        z_cyg_netcdf_dt_output = []

        pitch_cyg_netcdf_dt_output = [] 
        roll_cyg_netcdf_dt_output = [] 
        yaw_cyg_netcdf_dt_output = [] 

        x_gps_netcdf_dt_output = [] 
        y_gps_netcdf_dt_output = []
        z_gps_netcdf_dt_output = []

        vx_cyg_netcdf_dt_output = []
        vy_cyg_netcdf_dt_output = []
        vz_cyg_netcdf_dt_output = []
        x_spec_netcdf_dt_output = [] 
        y_spec_netcdf_dt_output = []
        z_spec_netcdf_dt_output = []
        lon_spec_netcdf_dt_output = [] 
        lat_spec_netcdf_dt_output = [] 


        gain_netcdf = []
        az_spec_netcdf = []
        el_spec_netcdf = []
        az_orbit_spec_netcdf = []
        el_orbit_spec_netcdf = []

        fom_netcdf = []
        index_netcdf = []
        gps_netcdf = []
        rx_to_sp_range_netcdf = []
        tx_to_sp_range_netcdf = []
        rcg_netcdf = []
        date_flight = []
        date_flight_rounded = []
        index_in_spock_date_netcdf_same = [] 
        index_in_spock_date_netcdf_same_dt_output = [] 
        date_flight_rounded_dt_output = []
        index_in_spock_not_interpolated_date_netcdf_same = []
        nb_seconds_since_initial_epoch_spock = []
        nb_seconds_since_initial_epoch_spock_dt_output = []
        iday_count = -1
        while iday_count < nb_day-1:# !!!!! should be nb_day-1:
            iday_count = iday_count + 1
            iday_here = day_list[iday_count]
            filename_spec_flight = "/Users/cbv/cygnss/netcdf/" + yy  + np.str(doy_array[iday_here]).zfill(3) + "/" + [x for x in os.listdir("/Users/cbv/cygnss/netcdf/" + yy  + np.str(doy_array[iday_here]).zfill(3)) if x.endswith(".nc") and x.startswith("cyg0" + np.str(cygfm))][0]
            fh = Dataset(filename_spec_flight, mode='r')
            # nc_attrs = fh.ncattrs()
            # nc_dims = [dim for dim in fh.dimensions]  # list of nc dimensions
            # nc_vars = [var for var in fh.variables]  # list of nc variables

            x_spec_netcdf_temp = fh.variables['sp_pos_x'][:] # X component of the specular point position in the ECEF coordinate system, in meters, at ddm_timestamp_utc, as calculated on the ground.
            y_spec_netcdf_temp = fh.variables['sp_pos_y'][:]
            z_spec_netcdf_temp = fh.variables['sp_pos_z'][:]

            lat_spec_netcdf_temp = fh.variables['sp_lat'][:]
            lon_spec_netcdf_temp = fh.variables['sp_lon'][:]

            x_cyg_netcdf_temp = fh.variables['sc_pos_x'][:]
            y_cyg_netcdf_temp = fh.variables['sc_pos_y'][:]
            z_cyg_netcdf_temp= fh.variables['sc_pos_z'][:]

            pitch_cyg_netcdf_temp = fh.variables['sc_pitch'][:] # Spacecraft pitch angle relative to the orbit frame, in radians at ddm_timestamp_utc
            roll_cyg_netcdf_temp = fh.variables['sc_roll'][:]
            yaw_cyg_netcdf_temp = fh.variables['sc_yaw'][:]

            x_gps_netcdf_temp = fh.variables['tx_pos_x'][:]
            y_gps_netcdf_temp = fh.variables['tx_pos_y'][:]
            z_gps_netcdf_temp= fh.variables['tx_pos_z'][:]




            gain_netcdf_temp = fh.variables['sp_rx_gain'][:] # The receive antenna gain in the direction of the specular point, in dBi, at ddm_timestamp_utc
            az_spec_netcdf_temp = fh.variables['sp_az_body'][:] # Let line A be the line that extends from the spacecraft to the specular point, at ddm_timestamp_utc. Let line B be the projection of line A onto the spacecraft body frame XY plane. sp_az_body is the angle between the spacecraft body frame +X axis and line B, in degrees, at ddm_timestamp_utc. See UM Doc. 148-0336, CYGNSS Science Data Processing Coordinate Systems Definitions.
            el_spec_netcdf_temp = fh.variables['sp_theta_body'][:] # The angle between the spacecraft body frame +Z axis and the line extending from the spacecraft to the specular point, in degrees, at ddm_timestamp_utc. See UM Doc. 148-0336, CYGNSS Science Data Processing Coordinate Systems Definitions.

            az_orbit_spec_netcdf_temp = fh.variables['sp_az_orbit'][:] # Let line A be the line that extends from the spacecraft to the specular point at ddm_timestamp_utc. Let line B be the projection of line A onto the orbit frame XY plane. sp_az_orbit is the angle between the orbit frame +X axis (the velocity vector) and line B, in degrees, at ddm_timestamp_utc. See UM Doc. 148-0336, CYGNSS Science Data Processing Coordinate Systems D
            el_orbit_spec_netcdf_temp = fh.variables['sp_theta_orbit'][:] # The angle between the orbit frame +Z axis and the line extending from the spacecraft to the specular point, in degrees, at ddm_timestamp_utc. See UM Doc. 148-0336, CYGNSS Science Data Processing Coordinate Systems Definitions.

            fom_netcdf_temp = fh.variables['prn_fig_of_merit'][:] # The RCG Figure of Merit (FOM) for the DDM. Ranges from 0 through 15.The DDMI selects the four strongest specular points (SP) for DDM production. It ranks the strength of SPs using an antenna RCG map. The map converts the position of the SP in antenna azimuth and declination angles to an RCG FOM. 0 represents the least FOM value. 15 represents the greatest FOM value.

            gps_netcdf_temp = fh.variables['prn_code'][:] # The PRN code of the GPS signal associated with the DDM. Ranges from 0 to 32. 0 = reflectometry channel idle. 1 through 32 = GPS PRN codes.
            vx_cyg_netcdf_temp = fh.variables['sc_vel_x'][:]
            vy_cyg_netcdf_temp = fh.variables['sc_vel_y'][:]
            vz_cyg_netcdf_temp= fh.variables['sc_vel_z'][:]

            rx_to_sp_range_netcdf_temp = np.double(fh.variables['rx_to_sp_range'][:])
            tx_to_sp_range_netcdf_temp = np.double(fh.variables['tx_to_sp_range'][:])




            list_are_masked_array = [] # sometimes the netcdf varaible below are maksed array and soemtimes they are not (depending on which netcdf file)...
            if type(x_spec_netcdf_temp) == ma.core.MaskedArray:
                list_are_masked_array.append(x_spec_netcdf_temp)
            if type(y_spec_netcdf_temp) == ma.core.MaskedArray:
                list_are_masked_array.append(y_spec_netcdf_temp)
            if type(z_spec_netcdf_temp) == ma.core.MaskedArray:
                list_are_masked_array.append(z_spec_netcdf_temp)
            if type(gain_netcdf_temp) == ma.core.MaskedArray:
                list_are_masked_array.append(gain_netcdf_temp)
            if type(az_spec_netcdf_temp) == ma.core.MaskedArray:
                list_are_masked_array.append(az_spec_netcdf_temp)
            if type(el_spec_netcdf_temp) == ma.core.MaskedArray:
                list_are_masked_array.append(el_spec_netcdf_temp)

            if type(az_orbit_spec_netcdf_temp) == ma.core.MaskedArray:
                list_are_masked_array.append(az_orbit_spec_netcdf_temp)
            if type(el_orbit_spec_netcdf_temp) == ma.core.MaskedArray:
                list_are_masked_array.append(el_orbit_spec_netcdf_temp)


            if type(fom_netcdf_temp) == ma.core.MaskedArray:
                list_are_masked_array.append(fom_netcdf_temp)

            if type(rx_to_sp_range_netcdf_temp) == ma.core.MaskedArray:
                list_are_masked_array.append(rx_to_sp_range_netcdf_temp)
            if type(tx_to_sp_range_netcdf_temp) == ma.core.MaskedArray:
                list_are_masked_array.append(tx_to_sp_range_netcdf_temp)

            if type(lat_spec_netcdf_temp) == ma.core.MaskedArray:
                list_are_masked_array.append(lat_spec_netcdf_temp)
            if type(lon_spec_netcdf_temp) == ma.core.MaskedArray:
                list_are_masked_array.append(lon_spec_netcdf_temp)

            nb_mask_array = len(list_are_masked_array)

            time_flight = fh.variables['ddm_timestamp_utc'][:]
            time_coverage_start = fh.getncattr(fh.ncattrs()[fh.ncattrs().index('time_coverage_start')])
            time_coverage_start_datetime = datetime.strptime(time_coverage_start[:-4], "%Y-%m-%dT%H:%M:%S.%f") 
            fh.close()
            nb_time_flight_temp = len(x_cyg_netcdf_temp)
            date_flight_t = []
            date_flight_rounded_temp = []
            time_remove_list = []
            itime = -1
            date_flight_raw = []
            date_flight_raw_every_second = []
            print iday_count, nb_day-1
            while itime < nb_time_flight_temp - 60:
                itime = itime + 1
                #print itime, nb_time_flight_temp-1, iday_count, nb_day-1
                time_remove = 0
                date_flight_temp_date = time_coverage_start_datetime + timedelta(microseconds = round(time_flight[itime]*10**6))

                date_flight_temp = datetime.strftime(date_flight_temp_date, "%Y-%m-%dT%H:%M:%S.%f" )
                # round to neared second but only if the actual time is less than 100 ms from the nearest second, otherwise ignore this time (don't compare to SpOCK)
                # .This is because SpOCK propagates with 1 s time step. So to compare to the netcdf file, we assume that the netcdf is also exactly at each second (so no millisecond). 
                #100 ms wrong is not too bad because the satellite movesby less than 1 km.
                if ( date_flight_temp.split('.')[1][0] == '9' ): # round to next second
                    date_flight_temp_date = datetime.strptime(date_flight_temp, "%Y-%m-%dT%H:%M:%S.%f")
                    date_flight_date = date_flight_temp_date + timedelta(seconds = 1)
                    date_flight_date_rounded_temp = datetime.strftime(date_flight_date, "%Y-%m-%dT%H:%M:%S.%f").split('.')[0]
                    date_flight_date_rounded = datetime.strptime(date_flight_date_rounded_temp, "%Y-%m-%dT%H:%M:%S")
                elif ( date_flight_temp.split('.')[1][0] == '0' ): # round to next second
                    date_flight_date_rounded = datetime.strptime(date_flight_temp.split('.')[0], "%Y-%m-%dT%H:%M:%S" )
                else: #if time can't be rounded by less than 100 ms
                    time_remove = 1

                # only select times in netcdf that are also date of SpOCK output (SpOCK output as specified in main input  file, not output imes of spec (these are very second))
                    #print np.mod( (date_flight_date_rounded - date_spock_not_interpolated[0]).total_seconds(), dt_spock_output ), date_flight_date_rounded, date_flight_temp_date

                date_flight_raw_every_second.append(date_flight_temp)
                while ( np.mod( (date_flight_date_rounded - date_spock_not_interpolated[0]).total_seconds(), dt_spock_output ) != 0 ):
                    itime = itime + 1
                    date_flight_temp_date = time_coverage_start_datetime + timedelta(microseconds = round(time_flight[itime]*10**6))
                    date_flight_temp = datetime.strftime(date_flight_temp_date, "%Y-%m-%dT%H:%M:%S.%f" )
                    date_flight_raw_every_second.append(date_flight_temp)
                    # round to neared second but only if the actual time is less than 100 ms from the nearest second, otherwise ignore this time (don't compare to SpOCK)
                    # .This is because SpOCK propagates with 1 s time step. So to compare to the netcdf file, we assume that the netcdf is also exactly at each second (so no millisecond). 
                    #100 ms wrong is not too bad because the satellite movesby less than 1 km.
                    if ( date_flight_temp.split('.')[1][0] == '9' ): # round to next second
                        date_flight_temp_date = datetime.strptime(date_flight_temp, "%Y-%m-%dT%H:%M:%S.%f")
                        date_flight_date = date_flight_temp_date + timedelta(seconds = 1)
                        date_flight_date_rounded_temp = datetime.strftime(date_flight_date, "%Y-%m-%dT%H:%M:%S.%f").split('.')[0]
                        date_flight_date_rounded = datetime.strptime(date_flight_date_rounded_temp, "%Y-%m-%dT%H:%M:%S")
                    elif ( date_flight_temp.split('.')[1][0] == '0' ): # round to next second
                        date_flight_date_rounded = datetime.strptime(date_flight_temp.split('.')[0], "%Y-%m-%dT%H:%M:%S" )
                    else: #if time can't be rounded by less than 100 ms
                        time_remove = 1


                        #print itime, time_coverage_start_datetime + timedelta(microseconds = round(time_flight[itime]*10**6)), x_cyg_netcdf_temp[itime]/1000.
                #print 'XXXXXXXXXXXXX', date_flight_temp

                if ( ( np.abs( pitch_cyg_netcdf_temp[itime] * 180. /np.pi ) > pitch_max ) | ( np.abs( roll_cyg_netcdf_temp[itime] * 180. /np.pi ) > roll_max ) | ( np.abs( yaw_cyg_netcdf_temp[itime] * 180. /np.pi ) > yaw_max ) ):
                    time_remove = 1
                    time_attitude_bad.append(itime)

#                 # !!!!!!!!! remove block below                                                     
#                 if type(fom_netcdf_temp) == ma.core.MaskedArray:
#                     if (np.min(fom_netcdf_temp.data[itime]) == 0):
#                         time_remove = 1
#                         time_gain_0.append(itime)
#                 else:
#                     if (np.min(fom_netcdf_temp[itime]) == 0):
#                         time_remove = 1
#                         time_gain_0.append(itime)
                    
#                 # !!!!!!!!! end of remove block below                                                     

                date_flight_raw.append(date_flight_temp)
                imask_arr = 0
                while (imask_arr < nb_mask_array):
                            #if ( ( True in x_spec_netcdf_temp.mask[itime] ) | ( True in y_spec_netcdf_temp.mask[itime] ) | ( True in z_spec_netcdf_temp.mask[itime] ) | ( True in gain_netcdf_temp.mask[itime] ) ):# id one of the 4 spec is masked then ignore this time
                    if (  True in list_are_masked_array[imask_arr].mask[itime] ):
                         time_remove = 1
                    imask_arr = imask_arr + 1     
                if ( (time_remove == 0) & (date_flight_date_rounded in date_spock)): # if this time is not in date_spock
                    index_in_spock_date_netcdf_same.append(date_spock.index(date_flight_date_rounded))
                else:
                    time_remove = 1
                if ( (time_remove == 0) & (date_flight_date_rounded in date_spock_not_interpolated)): # if this time is in date_spock_not_interpolated. If it's not, don't remove the time though. the time should be removed if it's not a time when the spec position 
                    #were predicted by spock. here we're looking at the sc position (every output time step chosen in the main input file).
                    index_in_spock_not_interpolated_date_netcdf_same.append(date_spock_not_interpolated.index(date_flight_date_rounded)) # save the index in date_spock_not_interpolated (ECEF_ file)
                    x_cyg_netcdf_dt_output.append(x_cyg_netcdf_temp[itime]/1000.)
                    y_cyg_netcdf_dt_output.append(y_cyg_netcdf_temp[itime]/1000.)
                    z_cyg_netcdf_dt_output.append(z_cyg_netcdf_temp[itime]/1000.)

                    pitch_cyg_netcdf_dt_output.append(pitch_cyg_netcdf_temp[itime] * 180. /np.pi)
                    roll_cyg_netcdf_dt_output.append(roll_cyg_netcdf_temp[itime] * 180. /np.pi)
                    yaw_cyg_netcdf_dt_output.append(yaw_cyg_netcdf_temp[itime] * 180. /np.pi)

                    x_gps_netcdf_dt_output.append(x_gps_netcdf_temp[itime]/1000.)
                    y_gps_netcdf_dt_output.append(y_gps_netcdf_temp[itime]/1000.)
                    z_gps_netcdf_dt_output.append(z_gps_netcdf_temp[itime]/1000.)

                    vx_cyg_netcdf_dt_output.append(vx_cyg_netcdf_temp[itime]/1000.)
                    vy_cyg_netcdf_dt_output.append(vy_cyg_netcdf_temp[itime]/1000.)
                    vz_cyg_netcdf_dt_output.append(vz_cyg_netcdf_temp[itime]/1000.)
                    x_spec_netcdf_dt_output.append(x_spec_netcdf_temp[itime]/1000.)
                    y_spec_netcdf_dt_output.append(y_spec_netcdf_temp[itime]/1000.)
                    z_spec_netcdf_dt_output.append(z_spec_netcdf_temp[itime]/1000.)
                    lon_spec_netcdf_dt_output.append(lon_spec_netcdf_temp[itime:itime+20]) # save 20 time steps here
                    lat_spec_netcdf_dt_output.append(lat_spec_netcdf_temp[itime:itime+20])
                    index_in_spock_date_netcdf_same_dt_output.append(date_spock.index(date_flight_date_rounded)) # save the index in date_spock (spec file)
                    date_flight_rounded_dt_output.append(date_flight_date_rounded)
                    nb_seconds_since_initial_epoch_spock_dt_output.append( ( date_flight_date_rounded - date_spock_not_interpolated[0] ).total_seconds() )

                if ( time_remove == 1 ): # remove time if can't be rounded by ess than 100 ms or if is not in date_spock or if masked or if attitude greater than roll_max, pitch_max, yaw_max
                    time_remove_list.append(itime)
                else:
                    if type(x_spec_netcdf_temp) == ma.core.MaskedArray:
                        x_spec_netcdf.append(x_spec_netcdf_temp.data[itime]/1000.)
                    else:
                        x_spec_netcdf.append(x_spec_netcdf_temp[itime]/1000.)
                    if type(y_spec_netcdf_temp) == ma.core.MaskedArray:
                        y_spec_netcdf.append(y_spec_netcdf_temp.data[itime]/1000.)
                    else:
                        y_spec_netcdf.append(y_spec_netcdf_temp[itime]/1000.)
                    if type(z_spec_netcdf_temp) == ma.core.MaskedArray:
                        z_spec_netcdf.append(z_spec_netcdf_temp.data[itime]/1000.)
                    else:
                        z_spec_netcdf.append(z_spec_netcdf_temp[itime]/1000.)
                    if type(gain_netcdf_temp) == ma.core.MaskedArray:
                        gain_netcdf.append(gain_netcdf_temp.data[itime])
                    else:
                        gain_netcdf.append(gain_netcdf_temp[itime])
                    if type(az_spec_netcdf_temp) == ma.core.MaskedArray:
                        az_spec_netcdf.append(az_spec_netcdf_temp.data[itime])
                    else:
                        az_spec_netcdf.append(az_spec_netcdf_temp[itime])

                    if type(el_spec_netcdf_temp) == ma.core.MaskedArray:
                        el_spec_netcdf.append(el_spec_netcdf_temp.data[itime])
                    else:
                        el_spec_netcdf.append(el_spec_netcdf_temp[itime])

                    if type(az_orbit_spec_netcdf_temp) == ma.core.MaskedArray:
                        az_orbit_spec_netcdf.append(az_orbit_spec_netcdf_temp.data[itime])
                    else:
                        az_orbit_spec_netcdf.append(az_orbit_spec_netcdf_temp[itime])

                    if type(el_orbit_spec_netcdf_temp) == ma.core.MaskedArray:
                        el_orbit_spec_netcdf.append(el_orbit_spec_netcdf_temp.data[itime])
                    else:
                        el_orbit_spec_netcdf.append(el_orbit_spec_netcdf_temp[itime])

                    if type(fom_netcdf_temp) == ma.core.MaskedArray:
                        fom_netcdf.append(fom_netcdf_temp.data[itime])
                    else:
                        fom_netcdf.append(fom_netcdf_temp[itime])

                    if type(gps_netcdf_temp) == ma.core.MaskedArray:
                        gps_netcdf.append(gps_netcdf_temp.data[itime])
                    else:
                        gps_netcdf.append(gps_netcdf_temp[itime])

                    if type(rx_to_sp_range_netcdf_temp) == ma.core.MaskedArray:
                        rx_to_sp_range_netcdf.append(rx_to_sp_range_netcdf_temp.data[itime])
                    else:
                        rx_to_sp_range_netcdf.append(rx_to_sp_range_netcdf_temp[itime])
                    if type(tx_to_sp_range_netcdf_temp) == ma.core.MaskedArray:
                        tx_to_sp_range_netcdf.append(tx_to_sp_range_netcdf_temp.data[itime])
                    else:
                        tx_to_sp_range_netcdf.append(tx_to_sp_range_netcdf_temp[itime])

                    if type(lat_spec_netcdf_temp) == ma.core.MaskedArray:
                        lat_spec_netcdf.append(lat_spec_netcdf_temp.data[itime])
                    else:
                        lat_spec_netcdf.append(lat_spec_netcdf_temp[itime])
                    if type(lon_spec_netcdf_temp) == ma.core.MaskedArray:
                        lon_spec_netcdf.append(lon_spec_netcdf_temp.data[itime])
                    else:
                        lon_spec_netcdf.append(lon_spec_netcdf_temp[itime])
                        
                    index_netcdf.append(itime) # can be used only for the first day of the simu otherwise doesnt make sense
                    rcg_temp_temp = 10**( gain_netcdf[-1] / 10. )
                    rcg_temp = rcg_temp_temp / (((tx_to_sp_range_netcdf[-1])**2) * ((rx_to_sp_range_netcdf[-1])**2)) * 10**27
                    rcg_netcdf.append(rcg_temp)
                    x_cyg_netcdf.append(x_cyg_netcdf_temp[itime]/1000.)
                    y_cyg_netcdf.append(y_cyg_netcdf_temp[itime]/1000.)
                    z_cyg_netcdf.append(z_cyg_netcdf_temp[itime]/1000.)

                    pitch_cyg_netcdf.append(pitch_cyg_netcdf_temp[itime] * 180. /np.pi)
                    roll_cyg_netcdf.append(roll_cyg_netcdf_temp[itime] * 180. /np.pi)
                    yaw_cyg_netcdf.append(yaw_cyg_netcdf_temp[itime] * 180. /np.pi)

                    x_gps_netcdf.append(x_gps_netcdf_temp[itime]/1000.)
                    y_gps_netcdf.append(y_gps_netcdf_temp[itime]/1000.)
                    z_gps_netcdf.append(z_gps_netcdf_temp[itime]/1000.)

                    date_flight_rounded.append(date_flight_date_rounded)
                    date_flight.append(date_flight_temp)
                    nb_seconds_since_initial_epoch_spock.append( ( date_flight_date_rounded - date_spock[0] ).total_seconds() )

                    vx_cyg_netcdf.append(vx_cyg_netcdf_temp[itime]/1000.)
                    vy_cyg_netcdf.append(vy_cyg_netcdf_temp[itime]/1000.)
                    vz_cyg_netcdf.append(vz_cyg_netcdf_temp[itime]/1000.)

                if date_flight_date_rounded > date_spock[-1]:
                    iday_count = nb_day + 1
                    break #######################################################################################################################################


        if len(date_flight_rounded) != 0: # if no data in the netcdf went trhough our different filters (attitude, masks, etc) then go to next date
            nb_time_flight_rounded = len(date_flight_rounded)

            date_netcdf = np.array( date_flight_rounded )
            date_spock = np.array(date_spock)
            nb_seconds_since_initial_epoch_spock = np.array(nb_seconds_since_initial_epoch_spock)
            nb_seconds_since_initial_epoch_spock_dt_output = np.array(nb_seconds_since_initial_epoch_spock_dt_output)

            # Compare SpOCK and netcdf
            ## Select the dates that are both in SpOCK and netcdf
            nb_time_spock_netcdf = len(date_flight) # number of time netcdf same as spock
            x_cyg_spock = np.array(x_cyg_spock); y_cyg_spock = np.array(y_cyg_spock); z_cyg_spock = np.array(z_cyg_spock)
            x_spec_spock = np.array(x_spec_spock); y_spec_spock = np.array(y_spec_spock); z_spec_spock = np.array(z_spec_spock); el_spec_spock = np.array(el_spec_spock); az_spec_spock = np.array(az_spec_spock); el_gps_from_cyg_spock = np.array(el_gps_from_cyg_spock); el_spec_not_int_spock = np.array(el_spec_not_int_spock); az_spec_not_int_spock = np.array(az_spec_not_int_spock)
            lat_spock = np.array(lat_spock); lon_spock = np.array(lon_spock)
            gain_spock = np.array(gain_spock)
            gps_spock = np.array(gps_spock)
            normpower_spock = np.array(normpower_spock)
            x_cyg_spock_same_time_as_netcdf = x_cyg_spock[index_in_spock_date_netcdf_same]
            y_cyg_spock_same_time_as_netcdf = y_cyg_spock[index_in_spock_date_netcdf_same]
            z_cyg_spock_same_time_as_netcdf = z_cyg_spock[index_in_spock_date_netcdf_same]
            x_spec_spock_same_time_as_netcdf = x_spec_spock[index_in_spock_date_netcdf_same]
            y_spec_spock_same_time_as_netcdf = y_spec_spock[index_in_spock_date_netcdf_same]
            z_spec_spock_same_time_as_netcdf = z_spec_spock[index_in_spock_date_netcdf_same]
            el_spec_spock_same_time_as_netcdf = el_spec_spock[index_in_spock_date_netcdf_same]
            el_gps_from_cyg_spock_same_time_as_netcdf = el_gps_from_cyg_spock[index_in_spock_date_netcdf_same]
            az_spec_spock_same_time_as_netcdf = az_spec_spock[index_in_spock_date_netcdf_same]
            el_spec_not_int_spock_same_time_as_netcdf = el_spec_not_int_spock[index_in_spock_date_netcdf_same]
            az_spec_not_int_spock_same_time_as_netcdf = az_spec_not_int_spock[index_in_spock_date_netcdf_same]
            gain_spock_same_time_as_netcdf = gain_spock[index_in_spock_date_netcdf_same]
            gps_spock_same_time_as_netcdf = gps_spock[index_in_spock_date_netcdf_same]
            normpower_spock_same_time_as_netcdf = normpower_spock[index_in_spock_date_netcdf_same]
            lon_spock_same_time_as_netcdf = lon_spock[index_in_spock_date_netcdf_same]
            lat_spock_same_time_as_netcdf = lat_spock[index_in_spock_date_netcdf_same]
            date_spock_same_time_as_netcdf  = date_spock[index_in_spock_date_netcdf_same]
            r_cyg_spock_ecef_file_same_time_as_netcdf = np.zeros([len(index_in_spock_not_interpolated_date_netcdf_same), 3])
            v_cyg_spock_ecef_file_same_time_as_netcdf = np.zeros([len(index_in_spock_not_interpolated_date_netcdf_same), 3])
            r_cyg_spock_ecef_file_same_time_as_netcdf[:, 0] = r_cyg_spock_ecef_file[index_in_spock_not_interpolated_date_netcdf_same, 0]
            r_cyg_spock_ecef_file_same_time_as_netcdf[:, 1] = r_cyg_spock_ecef_file[index_in_spock_not_interpolated_date_netcdf_same, 1]
            r_cyg_spock_ecef_file_same_time_as_netcdf[:, 2] = r_cyg_spock_ecef_file[index_in_spock_not_interpolated_date_netcdf_same, 2]
            v_cyg_spock_ecef_file_same_time_as_netcdf[:, 0] = v_cyg_spock_ecef_file[index_in_spock_not_interpolated_date_netcdf_same, 0]
            v_cyg_spock_ecef_file_same_time_as_netcdf[:, 1] = v_cyg_spock_ecef_file[index_in_spock_not_interpolated_date_netcdf_same, 1]
            v_cyg_spock_ecef_file_same_time_as_netcdf[:, 2] = v_cyg_spock_ecef_file[index_in_spock_not_interpolated_date_netcdf_same, 2]

            date_spock_not_interpolated = np.array(date_spock_not_interpolated)
            date_spock_not_interpolated_same_time_as_netcdf = date_spock_not_interpolated[index_in_spock_not_interpolated_date_netcdf_same]

            
            # IGNORE
            ## Calculate difference in position of sc and spec
            r_cyg_diff_mag = np.zeros([nb_time_spock_netcdf])
            r_spec_diff_mag_max_gain = np.zeros([nb_time_spock_netcdf])
            r_spec_diff_mag_min_diff = np.zeros([nb_time_spock_netcdf])
            gain_netcdf_that_min_error = np.zeros([nb_time_spock_netcdf])
            gain_spock_that_min_error = np.zeros([nb_time_spock_netcdf])
            which_comb_min = []
            r_spec_diff_mag_min_to_max = []
        #     r_spec_diff_mag_min_to_max_array = np.zeros([nb_time_spock_netcdf, 4]) - 1 # will be -1 if spock predicts less than 4 spec (if netcdf has less than 4 points than this time is ignored (see mask filter previously in the code))
        #     r_spec_diff_min_to_max_array_lvlh = np.zeros([nb_time_spock_netcdf, 4, 3]) - 1e6 # will be -1e6 if spock predicts less than 4 spec (if netcdf has less than 4 points than this time is ignored (see mask filter previously in the code))
            r_cyg_spock = np.zeros([nb_time_spock_netcdf, 3])
            r_cyg_netcdf = np.zeros([nb_time_spock_netcdf, 3])

            # Compare sc posiitons between netcdf and output by SpOCK in ECEF files. 
            # We do that here because we look at distance in lvlh so we need velocity in ECEF to convert to LVLH
            # Whenn I did the runs with SpOCK (and I didn't want to run everything again), the velcoity was not output 
            # in the spec files, so we need to look at the ECEF_ files 
            r_cyg_eci_spock = np.zeros([len(index_in_spock_not_interpolated_date_netcdf_same), 3])
            v_cyg_eci_spock = np.zeros([len(index_in_spock_not_interpolated_date_netcdf_same), 3])
            r_cyg_netcdf_dt_output = np.zeros([len(index_in_spock_not_interpolated_date_netcdf_same), 3])
            v_cyg_netcdf_dt_output = np.zeros([len(index_in_spock_not_interpolated_date_netcdf_same), 3])
            r_cyg_eci_netcdf_dt_output = np.zeros([len(index_in_spock_not_interpolated_date_netcdf_same), 3])
            v_cyg_eci_netcdf_dt_output = np.zeros([len(index_in_spock_not_interpolated_date_netcdf_same), 3])
            r_cyg_lvlh_diff = np.zeros([len(index_in_spock_not_interpolated_date_netcdf_same), 3])
            for itime in range(len(index_in_spock_not_interpolated_date_netcdf_same)):
        #         # ECEF to ECI for SpOCK
        #         if itime == 0:
        #             r_cyg_eci_spock[itime, :], v_cyg_eci_spock[itime, :] = ecef2eci(r_cyg_spock_ecef_file_same_time_as_netcdf[itime, :], v_cyg_spock_ecef_file_same_time_as_netcdf[itime, :], datetime.strftime(date_spock_not_interpolated_same_time_as_netcdf[itime], "%Y-%m-%dT%H:%M:%S"), 1)
        #         else:
        #             r_cyg_eci_spock[itime, :], v_cyg_eci_spock[itime, :] = ecef2eci(r_cyg_spock_ecef_file_same_time_as_netcdf[itime, :], v_cyg_spock_ecef_file_same_time_as_netcdf[itime, :], datetime.strftime(date_spock_not_interpolated_same_time_as_netcdf[itime], "%Y-%m-%dT%H:%M:%S"), 0)
                # ECEF to ECI for netcdf
                r_cyg_netcdf_dt_output[itime,:] = np.array([x_cyg_netcdf_dt_output[itime], y_cyg_netcdf_dt_output[itime], z_cyg_netcdf_dt_output[itime]])
                v_cyg_netcdf_dt_output[itime,:] = np.array([vx_cyg_netcdf_dt_output[itime], vy_cyg_netcdf_dt_output[itime], vz_cyg_netcdf_dt_output[itime]])
        #        r_cyg_eci_netcdf_dt_output[itime,:], v_cyg_eci_netcdf_dt_output[itime,:] = ecef2eci( r_cyg_netcdf_dt_output[itime,:], v_cyg_netcdf_dt_output[itime,:], datetime.strftime(date_spock_not_interpolated_same_time_as_netcdf[itime], "%Y-%m-%dT%H:%M:%S"), 0)
                # Distance in ECI from netcdf to SpOCK
                # r_cyg_eci_diff = r_cyg_eci_netcdf_dt_output[itime,:] - r_cyg_eci_spock[itime, :] 
                # Distance in LVLH (refce SpOCK) 
                #r_cyg_lvlh_diff[itime,:] = eci_to_lvlh(r_cyg_eci_spock[itime, :], v_cyg_eci_spock[itime, :], r_cyg_eci_diff)
                r_cyg_ecef_diff = r_cyg_netcdf_dt_output[itime,:] - r_cyg_spock_ecef_file_same_time_as_netcdf[itime, :]
                r_cyg_lvlh_diff[itime,:] = ecef_to_lvlh(r_cyg_spock_ecef_file_same_time_as_netcdf[itime, :], v_cyg_spock_ecef_file_same_time_as_netcdf[itime, :], r_cyg_ecef_diff)
            # end of IGNORE



            # Compare spec positions and sc positions from spec files.
            prn_comb_min_vs_time = []
            comb_min_vs_time = []
            spock_comb_min_vs_time = []
            spock_spec_lowest_to_highest_normpower = []
            netcdf_comb_min_vs_time = []
            netcdf_spec_lowest_to_highest_normpower = []
            which_normpower_rank_is_biggest_dist_spock_spec = np.zeros([nb_time_spock_netcdf])
            which_normpower_rank_is_biggest_dist_netcdf_spec = np.zeros([nb_time_spock_netcdf])

            r_spec_diff_lvlh_same_prn = np.zeros([nb_time_spock_netcdf, 4, 3]) - 1e30
            which_spec_diff = np.zeros([nb_time_spock_netcdf, 4]) - 1
            nb_prn_diff_at_given_time = np.zeros([nb_time_spock_netcdf])
            order_gain_messed_up =  np.zeros([nb_time_spock_netcdf]) - 1
            order_spec_netcdf = np.zeros([nb_time_spock_netcdf, 4])
            order_spec_spock = np.zeros([nb_time_spock_netcdf, nb_spec_spock])
            time_diff_prn = []
            time_diff_prn_gain_not_0 = []
            time_same_prn = []
            correct_vs_fom = np.zeros([16, 2]) # for each value of RCG (also noted FOM), compute percentage of correct PRN (correct = PRn necdf is same as PRN SpOCK). correct_vs_fom[ifom, 0] is the number of correct for this ifom, correct_vs_fom[ifor, 1] is the total number of spec with this fom. So percentrage correct = correct_vs_fom[ifom, 0] / correct_vs_fom[ifor, 1 * 100.
            diff_fom_same_prn_vs_fom_time = []#np.zeros([nb_time_spock_netcdf, 16]) - 1e30 # diffce betwe RCG netcdf and RCG SpOCK for spec with same PRN. Keep track of time and of which netcdf RCG the difference is for
            nb_spec_spock_not_get_matched_vs_fom = np.zeros([16])
            nb_spec_spock_got_matched_vs_fom = np.zeros([16])
            nb_spec_spock_got_matched_vs_fom_remove_ties = np.zeros([16]) # same as nb_spec_spock_got_matched_vs_fom but w/o ties. For instance, if there are 2 SpOCK RCG at 12 and both have their PRN matching two of the netcdf PRN for this time, only increment once nb_spec_spock_got_matched_vs_fom_remove_ties for this RCG, not twice (nb_spec_spock_got_matched_vs_fom is incremented twice)
            nb_spec_spock_vs_fom = np.zeros([16])
            nb_spec_spock_vs_fom_remove_ties = np.zeros([16]) # same as nb_spec_spock_vs_fom but w/o ties. For instance, if there are 2 SpOCK RCG at 12 and only one has its PRN matching one of the netcdf PRN for this time, only increment once nb_spec_spock_vs_fom_remove_ties once for this RCG, not twice (nb_spec_spock_vs_fom is incremented twice)
            list_prn_spock_got_matched_with_prn_netcdf = []
            list_prn_spock  = []
            spec_spock_order_descending_got_matched = np.zeros([nb_time_spock_netcdf, 4]) + 1e6
            spec_spock_order_descending_not_get_matched = []
            nb_time_4_highest_gain_got_matched = np.zeros([nb_time_spock_netcdf]) - 1
            nb_time_3_highest_gain_got_matched = np.zeros([nb_time_spock_netcdf]) - 1
            nb_ties_per_rcg_per_time = np.zeros([nb_time_spock_netcdf, 16]) - 1
            time_with_tie = []
            for itime in range(nb_time_spock_netcdf):#!!!!!!nb_time_spock_netcdf):
                #itime = 5 #### !!!!!! remove
            
                # IGNORE
                # difference sc position
                r_cyg_spock[itime,:] = np.array([x_cyg_spock_same_time_as_netcdf[itime][0], y_cyg_spock_same_time_as_netcdf[itime][0], z_cyg_spock_same_time_as_netcdf[itime][0]])
                r_cyg_netcdf[itime,:] = np.array([x_cyg_netcdf[itime], y_cyg_netcdf[itime], z_cyg_netcdf[itime]])
                r_cyg_diff = r_cyg_spock[itime,:] - r_cyg_netcdf[itime,:]
                r_cyg_diff_mag[itime] = np.linalg.norm(r_cyg_diff)
                # difference spec position
                #nb_spec_spock = len(x_spec_spock_same_time_as_netcdf[itime])
                nb_spec_netcdf = len(x_spec_netcdf[itime])
                r_spec_diff_mag_temp = []
                r_spec_diff_mag_temp_all = []
                where_max_gain_spock = np.where(gain_spock_same_time_as_netcdf[itime] == np.max(gain_spock_same_time_as_netcdf[itime]))[0][0]
                where_max_gain_netcdf = np.where(gain_netcdf[itime] == np.max(gain_netcdf[itime]))[0][0]
                r_spec_spock = np.array([x_spec_spock_same_time_as_netcdf[itime][where_max_gain_spock], y_spec_spock_same_time_as_netcdf[itime][where_max_gain_spock], z_spec_spock_same_time_as_netcdf[itime][where_max_gain_spock]])

                r_spec_netcdf = np.array([x_spec_netcdf[itime][where_max_gain_netcdf], y_spec_netcdf[itime][where_max_gain_netcdf], z_spec_netcdf[itime][where_max_gain_netcdf]])
                r_spec_diff = r_spec_spock - r_spec_netcdf
                r_spec_diff_mag_max_gain[itime] =  np.linalg.norm(r_spec_diff)
                dist_all_spec_to_all_spec = np.zeros([nb_spec_spock, nb_spec_netcdf])
                dist_all_spec_to_all_spec_lvlh = np.zeros([nb_spec_spock, nb_spec_netcdf,3])
                for ispec_spock in range(nb_spec_spock):
                    r_spec_spock = np.array([x_spec_spock_same_time_as_netcdf[itime][ispec_spock], y_spec_spock_same_time_as_netcdf[itime][ispec_spock], z_spec_spock_same_time_as_netcdf[itime][ispec_spock]])
                    r_spec_diff_mag_temp_per_spock = []
                    for ispec_netcdf in range(nb_spec_netcdf):
                        r_spec_netcdf = np.array([x_spec_netcdf[itime][ispec_netcdf], y_spec_netcdf[itime][ispec_netcdf], z_spec_netcdf[itime][ispec_netcdf]])
                        r_spec_diff = r_spec_spock - r_spec_netcdf
                        # Calculate disteance in LVLH frame 
                        # In theory, the LVLH frame hsould be taken in the frame of reference of each SpOCK spec point, but this is basically
                        # the same as the SpOCK satellite reference frame since the direction of motion of the satellite is very similar
                        # to the direction of motion of the specular point (the GPS don't move in a interval of a few seconds)
                        #r_spec_netcdf_eci, bla = ecef2eci( r_spec_netcdf, np.zeros([3]), datetime.strftime(date_flight_rounded[itime], "%Y-%m-%dT%H:%M:%S"), 0 ) # we don't have the velcooity of the spec but we don't care abwout it. 
                        #It's not needded to convert the velocity from ecef to eci (it's called by the funciton ecef2eci but used only to convert the ecef velcoity into eci velocity, so it doesn't  have any effect on the postiion)
                        r_spec_ecef_diff = r_spec_netcdf - r_spec_spock
                        r_spec_lvlh_diff = ecef_to_lvlh(r_cyg_spock_ecef_file_same_time_as_netcdf[itime, :], v_cyg_spock_ecef_file_same_time_as_netcdf[itime, :], r_spec_ecef_diff) # ok to take r_cyg_eci_spock at itime 
                        # because we sampled the ECEF and spec files at the same time (every 60s). 
                        # You can check by comparing the length of r_spec_diff_mag_min_to_max_array and 
                        # the length index_in_spock_date_netcdf_same_dt_output
                        dist_all_spec_to_all_spec_lvlh[ispec_spock, ispec_netcdf, :] = r_spec_lvlh_diff
                        dist_all_spec_to_all_spec[ispec_spock, ispec_netcdf] = np.linalg.norm(r_spec_diff)
                        r_spec_diff_mag_temp_per_spock.append( np.linalg.norm(r_spec_diff) )
                        r_spec_diff_mag_temp_all.append(np.linalg.norm(r_spec_diff))
                    r_spec_diff_mag_temp.append(r_spec_diff_mag_temp_per_spock)
                r_spec_diff_mag_min_diff[itime] = np.min(r_spec_diff_mag_temp_all)
                # end of IGNORE

                # different thing here. look at the dsitance between 2 spec with the same prn (this is what we should have always done)
                ## oder netcdf spec lowest to highest FOM
                index_sort_fom = np.argsort(fom_netcdf[itime])
                order_now = 0
                ispec_order_now = 0
                order_spec_netcdf[itime, index_sort_fom[ispec_order_now]] = 0
                already_diff_prn = 0
                already_diff_prn_gain_not_0 = 0
                for ispec_order in range(1,4):
                    if fom_netcdf[itime][index_sort_fom[ispec_order]] == fom_netcdf[itime][index_sort_fom[ispec_order-1]] :
                        order_spec_netcdf[itime, index_sort_fom[ispec_order]] = order_now
                    else:
                        order_now = order_now + 1
                        order_spec_netcdf[itime, index_sort_fom[ispec_order]] = order_now


#                 index_sort_gain_spock = np.argsort(gain_spock_same_time_as_netcdf[itime])
#                 order_now = 0
#                 ispec_order_now = 0
#                 order_spec_spock[itime, index_sort_gain_spock[ispec_order_now]] = 0
#                 for ispec_order in range(1,nb_spec_spock):
#                     if gain_spock_same_time_as_netcdf[itime][index_sort_gain_spock[ispec_order]] == gain_spock_same_time_as_netcdf[itime][index_sort_gain_spock[ispec_order-1]] :
#                         order_spec_spock[itime, index_sort_gain_spock[ispec_order]] = order_now
#                     else:
#                         order_now = order_now + 1
#                         order_spec_spock[itime, index_sort_gain_spock[ispec_order]] = order_now

                # Sort SpOCK gain in descending order. Doesn't matter if 2 gains are the same
                index_sort_gain_spock = np.argsort(-gain_spock_same_time_as_netcdf[itime])
                order_now = 0
                ispec_order_now = 0
                order_spec_spock[itime, index_sort_gain_spock[ispec_order_now]] = 0
                for ispec_order in range(1,nb_spec_spock):
                    if gain_spock_same_time_as_netcdf[itime][index_sort_gain_spock[ispec_order]] == gain_spock_same_time_as_netcdf[itime][index_sort_gain_spock[ispec_order-1]] :
                        order_spec_spock[itime, index_sort_gain_spock[ispec_order]] = order_now
                        if ispec_order + 1 < nb_spec_spock:
                            if gain_spock_same_time_as_netcdf[itime][index_sort_gain_spock[ispec_order+1]] != gain_spock_same_time_as_netcdf[itime][index_sort_gain_spock[ispec_order]] :
                                order_now = order_now + 1
                        
                    else:
                        order_now = order_now + 1
                        order_spec_spock[itime, index_sort_gain_spock[ispec_order]] = order_now


                ispec_count = 0 
                list_gain_spock_this_time = gain_spock_same_time_as_netcdf[itime]
                list_prn_spock_this_time = gps_spock_same_time_as_netcdf[itime]
                list_prn_spock_got_matched_with_prn_netcdf_this_time = []
                spec_spock_order_descending_not_get_matched_this_time = []

                for ispec in range(4):
                    ifom = fom_netcdf[itime][ispec]
                    correct_vs_fom[ifom, 1] = correct_vs_fom[ifom, 1] + 1
                    prn_netcdf_here = gps_netcdf[itime][ispec]
                    if len(np.where(gps_spock_same_time_as_netcdf[itime] == prn_netcdf_here)[0]) > 0:
                        correct_vs_fom[ifom, 0] = correct_vs_fom[ifom, 0] + 1
                        where_prn_spock = np.where(gps_spock_same_time_as_netcdf[itime] == prn_netcdf_here)[0][0]
                        list_prn_spock_got_matched_with_prn_netcdf_this_time.append(gps_spock_same_time_as_netcdf[itime][where_prn_spock])
                        diff_fom_same_prn_vs_fom_time.append( [itime, ispec, fom_netcdf[itime][ispec], gain_spock_same_time_as_netcdf[itime][where_prn_spock] - fom_netcdf[itime][ispec],  where_prn_spock] )
                        r_spec_spock = np.array([x_spec_spock_same_time_as_netcdf[itime][where_prn_spock], y_spec_spock_same_time_as_netcdf[itime][where_prn_spock], z_spec_spock_same_time_as_netcdf[itime][where_prn_spock]])
                        r_spec_netcdf = np.array([x_spec_netcdf[itime][ispec], y_spec_netcdf[itime][ispec], z_spec_netcdf[itime][ispec]])
                        r_spec_diff = r_spec_spock - r_spec_netcdf
                        r_spec_ecef_diff = r_spec_netcdf - r_spec_spock
                        r_spec_lvlh_diff = ecef_to_lvlh(r_cyg_spock_ecef_file_same_time_as_netcdf[itime, :], v_cyg_spock_ecef_file_same_time_as_netcdf[itime, :], r_spec_ecef_diff) # ok to take r_cyg_eci_spock at itime 
                        r_spec_diff_lvlh_same_prn[itime, ispec, :] = r_spec_lvlh_diff
                    else:
                        which_spec_diff[itime, ispec] = order_spec_netcdf[itime, ispec]
                        if already_diff_prn == 0:
                            time_diff_prn.append(itime)
                            if ifom > 3:
                                time_diff_prn_gain_not_0.append(itime)
                                already_diff_prn_gain_not_0 = 1
                        already_diff_prn = 1
                #nb_spec_spock_not_matched = len(list_prn_spock_this_time) - len(list_prn_spock_got_matched_with_prn_netcdf_this_time) # nb of spock spec that didn't get matched with a netcdf prn for this time
                ispec_pos = 0

                fom_spock_already_seen_this_time = np.zeros([16])
                fom_spock_already_matched_this_time = np.zeros([16])
                for ispec_spock in range(nb_spec_spock):
                    ifom_spock = (int)(gain_spock_same_time_as_netcdf[itime][ispec_spock])
                    # look at number of ties
                    nb_ties_per_rcg_per_time[itime, ifom_spock] = len(np.where(gain_spock_same_time_as_netcdf[itime] == ifom_spock)[0]) - 1 # -1 bcause if len is 1 it means there is one spec with this rcg, which is not a tie. if len is 2 then 2 spec have this rcg, so it should be counted as 1 tie. if len is 3 then 3 spec have the same rcgg so it is counted as 2 ties
                    # end of look at number of ties
                    nb_spec_spock_vs_fom[ifom_spock] = nb_spec_spock_vs_fom[ifom_spock] + 1 
                    if fom_spock_already_seen_this_time[ifom_spock] != 1:
                        nb_spec_spock_vs_fom_remove_ties[ifom_spock] = nb_spec_spock_vs_fom_remove_ties[ifom_spock] + 1 
                        
                    if (list_prn_spock_this_time[ispec_spock] in list_prn_spock_got_matched_with_prn_netcdf_this_time) == True:
                        
                        nb_spec_spock_got_matched_vs_fom[ifom_spock] = nb_spec_spock_got_matched_vs_fom[ifom_spock] + 1
                        spec_spock_order_descending_got_matched[itime, ispec_pos] = order_spec_spock[itime, ispec_spock] 
                        ispec_pos = ispec_pos + 1
                        if fom_spock_already_matched_this_time[ifom_spock] != 1:
                            nb_spec_spock_got_matched_vs_fom_remove_ties[ifom_spock] = nb_spec_spock_got_matched_vs_fom_remove_ties[ifom_spock] + 1
                        fom_spock_already_matched_this_time[ifom_spock] = 1

                    else:
                        nb_spec_spock_not_get_matched_vs_fom[ifom_spock] = nb_spec_spock_not_get_matched_vs_fom[ifom_spock] + 1
                        spec_spock_order_descending_not_get_matched_this_time.append( order_spec_spock[itime, ispec_spock] )

                    fom_spock_already_seen_this_time[ifom_spock] = 1
                if len(np.where(nb_ties_per_rcg_per_time[itime, :] >= 1)[0]): # at least one tie for this time
                    time_with_tie.append(itime)

                if np.max(spec_spock_order_descending_got_matched[itime, :]) <= 3:
                    nb_time_4_highest_gain_got_matched[itime] = 1

                if len(np.where(spec_spock_order_descending_got_matched[itime, :] <= 2)[0]) >= 3:
                    nb_time_3_highest_gain_got_matched[itime] = 1
                
                # block belwo is to ignore times when RCG 0 is once or twice in netcdf RCG
                if ( len(np.where(fom_netcdf[itime] == 0)[0]) >= 1 ):
                    nb_time_4_highest_gain_got_matched[itime] = 1
                    if ( len(np.where(fom_netcdf[itime] == 0)[0]) >= 2  ):
                        nb_time_3_highest_gain_got_matched[itime] = 1


                list_prn_spock.append(list_prn_spock_this_time)
                list_prn_spock_got_matched_with_prn_netcdf.append(list_prn_spock_got_matched_with_prn_netcdf_this_time)
                spec_spock_order_descending_not_get_matched.append(spec_spock_order_descending_not_get_matched_this_time)


            time_diff_prn_arr = np.array(time_diff_prn)
            time_same_prn_arr = np.array(time_same_prn)
            gps_spock = gps_spock_same_time_as_netcdf 
            fom_spock = gain_spock_same_time_as_netcdf
            which_spec_diff_time_diff_prn = which_spec_diff[time_diff_prn, :]
            fom_netcdf_arr = np.array(fom_netcdf)
            fom_netcdf_diff_prn_temp = fom_netcdf_arr[time_diff_prn, :]
            fom_netcdf_diff_prn = []
            fom_netcdf_diff_prn_all = []
            order_gain_netcdf_diff_prn = np.zeros([4])
            nb_spec_diff = np.zeros([len(time_diff_prn)])
            nb_spec_diff_gain_not_0 = 0
            for itime in range(len(time_diff_prn)):
                fom_netcdf_diff_prn_itime = []
                if 0 in which_spec_diff_time_diff_prn[itime, :]:
                    order_gain_netcdf_diff_prn[0] = order_gain_netcdf_diff_prn[0] + 1
                if 1 in which_spec_diff_time_diff_prn[itime, :]:
                    order_gain_netcdf_diff_prn[1] = order_gain_netcdf_diff_prn[1] + 1
                if 2 in which_spec_diff_time_diff_prn[itime, :]:
                    order_gain_netcdf_diff_prn[2] = order_gain_netcdf_diff_prn[2] + 1
                if 3 in which_spec_diff_time_diff_prn[itime, :]:
                    order_gain_netcdf_diff_prn[3] = order_gain_netcdf_diff_prn[3] + 1
                fom_netcdf_diff_prn.append( fom_netcdf_diff_prn_temp[itime, np.where(which_spec_diff_time_diff_prn[itime, :] != -1)[0] ] )
                for  idiff in range(len(np.where(which_spec_diff_time_diff_prn[itime, :] != -1)[0])):
                    fom_netcdf_diff_prn_all.append( fom_netcdf_diff_prn_temp[itime, np.where(which_spec_diff_time_diff_prn[itime, :] != -1)[0][idiff] ] )
                    if fom_netcdf_diff_prn_temp[itime, np.where(which_spec_diff_time_diff_prn[itime, :] != -1)[0][idiff] ] != 0:
                        nb_spec_diff_gain_not_0 = nb_spec_diff_gain_not_0 +1 
                nb_spec_diff[itime] = len(np.where(which_spec_diff_time_diff_prn[itime, :] != -1)[0])

            correct_vs_fom_percent = np.zeros([16])
            for ifom in range(16):
                correct_vs_fom_percent[ifom] = correct_vs_fom[ifom, 0] * 100. / correct_vs_fom[ifom, 1]

            pickle.dump(correct_vs_fom, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_correct_vs_fom_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))
            pickle.dump(correct_vs_fom_percent, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_correct_vs_fom_percent_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))
            pickle.dump(time_diff_prn, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_time_diff_prn_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))
            pickle.dump(time_same_prn, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_time_same_prn_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))
            pickle.dump(which_spec_diff, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_which_spec_diff_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))
            pickle.dump(gps_spock, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_gps_spock_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))
            pickle.dump(fom_spock, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_fom_spock_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))
            pickle.dump(gps_netcdf, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_gps_netcdf_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))
            pickle.dump(fom_netcdf, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_fom_netcdf_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))

            pickle.dump(fom_netcdf_diff_prn_all, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_fom_netcdf_diff_prn_all_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))
            pickle.dump(nb_spec_diff, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_nb_spec_diff_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))
            pickle.dump(nb_spec_diff_gain_not_0, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_nb_spec_diff_gain_not_0_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))
            pickle.dump(r_spec_diff_lvlh_same_prn, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_r_spec_diff_lvlh_same_prn_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))
            pickle.dump(nb_seconds_since_initial_epoch_spock, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_nb_seconds_since_initial_epoch_spock_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))
            pickle.dump(diff_fom_same_prn_vs_fom_time, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_diff_fom_same_prn_vs_fom_time_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))
            pickle.dump(az_spec_netcdf, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_az_spec_netcdf_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))
            pickle.dump(el_spec_netcdf, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_el_spec_netcdf_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))
            pickle.dump(order_gain_netcdf_diff_prn, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_order_gain_netcdf_diff_prn_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))
            pickle.dump(time_diff_prn_gain_not_0, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_time_diff_prn_gain_not_0_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))

            pickle.dump(nb_spec_spock_got_matched_vs_fom, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_nb_spec_spock_got_matched_vs_fom_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))
            pickle.dump(nb_spec_spock_vs_fom, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_nb_spec_spock_vs_fom_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle", "w"))

            correct_vs_fom_all_date.append( correct_vs_fom)
            correct_vs_fom_percent_all_date.append(  correct_vs_fom_percent)
            time_diff_prn_all_date.append(  time_diff_prn)
            time_same_prn_all_date.append(  time_same_prn)
            which_spec_diff_all_date.append(  which_spec_diff)
            gps_spock_all_date.append(  gps_spock)
            fom_spock_all_date.append(  fom_spock)
            gps_netcdf_all_date.append(  gps_netcdf)
            fom_netcdf_all_date.append(  fom_netcdf)
            fom_netcdf_diff_prn_all_all_date.append( fom_netcdf_diff_prn_all   )
            nb_spec_diff_all_date.append(  nb_spec_diff )
            nb_spec_diff_gain_not_0_all_date.append(  nb_spec_diff_gain_not_0 )
            r_spec_diff_lvlh_same_prn_all_date.append(  r_spec_diff_lvlh_same_prn )
            nb_seconds_since_initial_epoch_spock_all_date.append(  nb_seconds_since_initial_epoch_spock )

            diff_fom_same_prn_vs_fom_time_all_date.append( diff_fom_same_prn_vs_fom_time )
            az_spec_netcdf_all_date.append( az_spec_netcdf )
            el_spec_netcdf_all_date.append( el_spec_netcdf )
            order_gain_netcdf_diff_prn_all_date.append( order_gain_netcdf_diff_prn )
            time_diff_prn_gain_not_0_all_date.append( time_diff_prn_gain_not_0 )

            nb_spec_spock_got_matched_vs_fom_all_date.append(nb_spec_spock_got_matched_vs_fom)
            nb_spec_spock_vs_fom_all_date.append(nb_spec_spock_vs_fom)

            print correct_vs_fom_percent
            print nb_spec_spock_got_matched_vs_fom * 100. / nb_spec_spock_vs_fom
#             hist_here = np.histogram(fom_netcdf_diff_prn_all, bins = np.arange(-.5, 16.5, 1))[0]*100./len(fom_netcdf_diff_prn_all)
#             print 'distribution', hist_here, np.max(hist_here[1:])
            print '% correct if gain not 0', 100 -  nb_spec_diff_gain_not_0 *100. / (4*len(fom_netcdf))
            # below should be similar to nb_spec_diff_gain_not_0 *100. / (4*len(fom_netcdf))
            print np.sum(correct_vs_fom[2:,1] / np.sum(correct_vs_fom[2:,1]) * correct_vs_fom_percent[2:]) 
            print '% correct all gains', 100 - len(fom_netcdf_diff_prn_all) * 100. / (4*len(fom_netcdf))



    else:
        nb_spec_spock = 4
        if os.path.exists("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_correct_vs_fom_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"): # if netcdf data didn't go thorugh the differerent filters (attitude, mask,...) then these pickles were not created)
        
            correct_vs_fom = pickle.load(open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_correct_vs_fom_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))
            correct_vs_fom_percent = pickle.load( open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_correct_vs_fom_percent_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))
            time_diff_prn = pickle.load( open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_time_diff_prn_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))
            time_same_prn = pickle.load( open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_time_same_prn_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))
            which_spec_diff = pickle.load(open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_which_spec_diff_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))
            gps_spock = pickle.load(open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_gps_spock_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))
            fom_spock = pickle.load(open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_fom_spock_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))
            gps_netcdf = pickle.load(open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_gps_netcdf_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))
            fom_netcdf = pickle.load(open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_fom_netcdf_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))

            fom_netcdf_diff_prn_all = pickle.load(open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_fom_netcdf_diff_prn_all_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))
            nb_spec_diff = pickle.load( open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_nb_spec_diff_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))
            nb_spec_diff_gain_not_0 = pickle.load( open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_nb_spec_diff_gain_not_0_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))
            r_spec_diff_lvlh_same_prn = pickle.load( open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_r_spec_diff_lvlh_same_prn_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))
            nb_seconds_since_initial_epoch_spock = pickle.load( open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_nb_seconds_since_initial_epoch_spock_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))


            diff_fom_same_prn_vs_fom_time = pickle.load(open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_diff_fom_same_prn_vs_fom_time_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))
            az_spec_netcdf = pickle.load( open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_az_spec_netcdf_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))
            el_spec_netcdf = pickle.load( open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_el_spec_netcdf_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))
            order_gain_netcdf_diff_prn = pickle.load( open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_order_gain_netcdf_diff_prn_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))
            time_diff_prn_gain_not_0 = pickle.load(open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_time_diff_prn_gain_not_0_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))

            nb_spec_spock_got_matched_vs_fom = pickle.load(open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_nb_spec_spock_got_matched_vs_fom_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))
            nb_spec_spock_vs_fom = pickle.load(open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_nb_spec_spock_vs_fom_nSP_" + str(nb_spec_spock) + "FM" + str(cygfm) + ".pickle"))


            correct_vs_fom_all_date.append( correct_vs_fom)
            correct_vs_fom_percent_all_date.append(  correct_vs_fom_percent)
            time_diff_prn_all_date.append(  time_diff_prn)
            time_same_prn_all_date.append(  time_same_prn)
            which_spec_diff_all_date.append(  which_spec_diff)
            gps_spock_all_date.append(  gps_spock)
            fom_spock_all_date.append(  fom_spock)
            gps_netcdf_all_date.append(  gps_netcdf)
            fom_netcdf_all_date.append(  fom_netcdf)
            fom_netcdf_diff_prn_all_all_date.append( fom_netcdf_diff_prn_all   )
            nb_spec_diff_all_date.append(  nb_spec_diff )
            nb_spec_diff_gain_not_0_all_date.append(  nb_spec_diff_gain_not_0 )
            r_spec_diff_lvlh_same_prn_all_date.append(  r_spec_diff_lvlh_same_prn )
            nb_seconds_since_initial_epoch_spock_all_date.append(  nb_seconds_since_initial_epoch_spock )

            diff_fom_same_prn_vs_fom_time_all_date.append( diff_fom_same_prn_vs_fom_time )
            az_spec_netcdf_all_date.append( az_spec_netcdf )
            el_spec_netcdf_all_date.append( el_spec_netcdf )
            order_gain_netcdf_diff_prn_all_date.append( order_gain_netcdf_diff_prn )
            time_diff_prn_gain_not_0_all_date.append( time_diff_prn_gain_not_0 )

            nb_spec_spock_got_matched_vs_fom_all_date.append(nb_spec_spock_got_matched_vs_fom)
            nb_spec_spock_vs_fom_all_date.append(nb_spec_spock_vs_fom)


            print correct_vs_fom_percent
            print nb_spec_spock_got_matched_vs_fom * 100. / nb_spec_spock_vs_fom
            #hist_here = np.histogram(fom_netcdf_diff_prn_all, bins = np.arange(-.5, 16.5, 1))[0]*100./len(fom_netcdf_diff_prn_all)
            #print 'distribution', hist_here, np.max(hist_here[1:])
            print '% correct if gain not 0', 100 -  nb_spec_diff_gain_not_0 *100. / (4*len(fom_netcdf))
            # below should be similar to nb_spec_diff_gain_not_0 *100. / (4*len(fom_netcdf))
            print np.sum(correct_vs_fom[2:,1] / np.sum(correct_vs_fom[2:,1]) * correct_vs_fom_percent[2:]) 
            print '% correct all gains', 100 - len(fom_netcdf_diff_prn_all) * 100. / (4*len(fom_netcdf))






#for idate in range(nb_date):
    
correct_vs_fom_percent_all_date_arr = np.array(correct_vs_fom_percent_all_date)
correct_vs_fom_percent_all_date_average = np.mean(correct_vs_fom_percent_all_date_arr, axis = 0)

nb_spec_spock_got_matched_vs_fom_all_date_arr = np.array(nb_spec_spock_got_matched_vs_fom_all_date)
nb_spec_spock_vs_fom_all_date_arr = np.array(nb_spec_spock_vs_fom_all_date)

percentage_nb_spec_spock_got_matched_vs_fom_all_date_arr = nb_spec_spock_got_matched_vs_fom_all_date_arr / nb_spec_spock_vs_fom_all_date_arr * 100.
average_percentage_nb_spec_spock_got_matched_vs_fom_all_date_arr = np.mean(percentage_nb_spec_spock_got_matched_vs_fom_all_date_arr, axis = 0)

print nb_spec_spock_got_matched_vs_fom * 100. / nb_spec_spock_vs_fom

nb_date_ok  = len(nb_spec_diff_gain_not_0_all_date)

# nb_spec_diff_gain_not_0_all_date_arr = np.array(nb_spec_diff_gain_not_0_all_date)
# percentage_nb_spec_diff_gain_not_0_all_date = np.zeros([nb_date_ok])
# percentage_nb_spec_diff_all_gain_all_date = np.zeros([nb_date_ok])
# nb_spec_per_gain_all = np.zeros([nb_date_ok, 16])
# for idate in range(nb_date_ok):
#     percentage_nb_spec_diff_gain_not_0_all_date[idate] = 100. -  nb_spec_diff_gain_not_0_all_date_arr[idate]  *100. / (4*len(fom_netcdf_all_date[idate]))
#     nb_spec_per_gain_all[idate, :] = correct_vs_fom_all_date[idate][:, 1]
#     percentage_nb_spec_diff_all_gain_all_date[idate] = 100. - len(fom_netcdf_diff_prn_all_all_date[idate]) * 100. / (4*len(fom_netcdf_all_date[idate]))
# nb_spec_per_gain_all_average = np.mean(nb_spec_per_gain_all, axis = 0)




# median and quantile over a day for a particular date
median_r_spec_diff_lvlh_all_date = []
q10_r_spec_diff_lvlh_all_date = []
q25_r_spec_diff_lvlh_all_date = []
q75_r_spec_diff_lvlh_all_date = []
q90_r_spec_diff_lvlh_all_date = []
itime_day_list_all_date = []
for idate in range(nb_date_ok):
    nb_day_simu = (int)(np.ceil((nb_seconds_since_initial_epoch_spock_all_date[idate][-1] - nb_seconds_since_initial_epoch_spock_all_date[idate][0] )/ 3600. / 24.))
    #nb_day_simu = 7 # !!!! remove this line and uncomment line above
    median_r_spec_diff_lvlh = np.zeros([nb_day_simu, 4, 3])
    q10_r_spec_diff_lvlh = np.zeros([nb_day_simu, 4, 3])
    q25_r_spec_diff_lvlh = np.zeros([nb_day_simu, 4, 3])
    q75_r_spec_diff_lvlh = np.zeros([nb_day_simu, 4, 3])
    q90_r_spec_diff_lvlh = np.zeros([nb_day_simu, 4, 3])
    x_axis_day = np.arange(0.5, nb_day_simu+0.5)+nb_seconds_since_initial_epoch_spock_all_date[idate][0] / 3600. / 24. # start at 0.5 and go up to nb_day_simu + 0.5 
    #because say that the median and quantiles are over a day so take it at the middle of the day
    #x_axis_day = np.arange(0.5, nb_day_simu+0.5) # !!!!!!!!! remove this line and uncomment line above

    iday = 0
    itime_previous = 0
    itime_day_list = [] # list of itime for that falls every day
    itime_day_list.append(itime_previous)
    while iday < nb_day_simu:
        itime = 0
        while (int)( ( nb_seconds_since_initial_epoch_spock_all_date[idate][itime+itime_previous] - nb_seconds_since_initial_epoch_spock_all_date[idate][itime_previous] )/ 3600. / 24.) < 1:
            itime = itime + 1
            if itime+itime_previous == r_spec_diff_lvlh_same_prn_all_date[idate].shape[0]:
                break
#         # !!!!!!!!! remove this block and ucmmoent block above
#         while (itime < len(fom_netcdf_all_date[idate]) / nb_day_simu):
#             itime = itime + 1
#             if itime+itime_previous == r_spec_diff_lvlh_same_prn_all_date[idate].shape[0]:
#                 break
#         # !!!!!!!!! end of remove this block and ucmmoent block above

        median_r_spec_diff_lvlh[iday, :, :] = np.median(r_spec_diff_lvlh_same_prn_all_date[idate][itime_previous:itime+itime_previous, :, :], axis = 0)
        q10_r_spec_diff_lvlh[iday, :, :] = np.percentile(r_spec_diff_lvlh_same_prn_all_date[idate][itime_previous:itime+itime_previous, :, :], 10, axis = 0)
        q25_r_spec_diff_lvlh[iday, :, :] = np.percentile(r_spec_diff_lvlh_same_prn_all_date[idate][itime_previous:itime+itime_previous, :, :], 25, axis = 0)
        q75_r_spec_diff_lvlh[iday, :, :] = np.percentile(r_spec_diff_lvlh_same_prn_all_date[idate][itime_previous:itime+itime_previous, :, :], 75, axis = 0)
        q90_r_spec_diff_lvlh[iday, :, :] = np.percentile(r_spec_diff_lvlh_same_prn_all_date[idate][itime_previous:itime+itime_previous, :, :], 90, axis = 0)
        #print itime_previous, itime+itime_previous, median_r_spec_diff_lvlh[iday, 0, 0] 
        iday = iday + 1
        itime_previous = itime + itime_previous
        itime_day_list.append(itime_previous)

    median_r_spec_diff_lvlh_all_date.append(median_r_spec_diff_lvlh)
    q10_r_spec_diff_lvlh_all_date.append(q10_r_spec_diff_lvlh)
    q25_r_spec_diff_lvlh_all_date.append(q25_r_spec_diff_lvlh)
    q75_r_spec_diff_lvlh_all_date.append(q75_r_spec_diff_lvlh)
    q90_r_spec_diff_lvlh_all_date.append(q90_r_spec_diff_lvlh)
    itime_day_list_all_date.append(itime_day_list)



raise Exception




#PRESENTATION 051518
# Accuracy of the selection algorithm
nPrn = 0
nPrnG = 0
nPrnGainNot0 = 0
nPrnGainNot0G = 0
percentage_spec_gain_0 = np.zeros([nb_date_ok])
order_gain_netcdf_diff_prn_all_date_arr = np.array(order_gain_netcdf_diff_prn_all_date)
percentage_exactly_one_spec_wrong = np.zeros([nb_date_ok])
for idate in range(nb_date_ok):
    nPrnG = nPrnG + np.sum(nb_spec_spock_got_matched_vs_fom_all_date[idate])
    nPrn = nPrn + np.sum(nb_spec_spock_vs_fom_all_date[idate])
    nPrnGainNot0 = nPrnGainNot0 + np.sum(nb_spec_spock_vs_fom_all_date[idate][1:])
    nPrnGainNot0G = nPrnGainNot0G + np.sum(nb_spec_spock_got_matched_vs_fom_all_date[idate][1:])
    percentage_spec_gain_0[idate] = len(np.where(fom_spock_all_date[idate] == 0)[0]) * 100. / (fom_spock_all_date[idate].shape[0] * fom_spock_all_date[idate].shape[1] )
    percentage_exactly_one_spec_wrong[idate] = len(np.where(nb_spec_diff_all_date[idate] == 1)[0])*100./len(nb_spec_diff_all_date[idate])

print 'Accuracy of the selection algorithm', nPrnG * 100./nPrn
print 'If SPs with a gain of 0 were excluded, accuracy of selection algorithm', nPrnGainNot0G * 100./ nPrnGainNot0
print 'Percentage of SPs with gain 0', np.mean(percentage_spec_gain_0)
print 'Percentage of time steps with exactly one SP misidentified', np.mean(percentage_exactly_one_spec_wrong)
# right below: if element 0 is 94 it means that 94% of the SPs not correctly identified by SpOCK were the SPs with the lowest gain at their corresponding time step. If elemetn 1 is 5 it 
# means that 5% of the SPs not correctly identified by SpOCK were the SPs with the second to lowest gain at their corresponding time step. etc until element 4 (highest gain)
print 'Percentage of the SPs not correctly identified per RCG order',  np.sum(order_gain_netcdf_diff_prn_all_date_arr, axis = 0) * 100./np.sum(order_gain_netcdf_diff_prn_all_date_arr) 
#end of PRESENTATION 051518








#     # Look at x y z of spock and netcdf of each set of spec point (spock, netcdf)
#     x1_spock = np.zeros([nb_time_spock_netcdf])-1e10; x1_netcdf = np.zeros([nb_time_spock_netcdf])-1e10
#     y1_spock = np.zeros([nb_time_spock_netcdf])-1e10; y1_netcdf = np.zeros([nb_time_spock_netcdf])-1e10
#     z1_spock = np.zeros([nb_time_spock_netcdf])-1e10; z1_netcdf = np.zeros([nb_time_spock_netcdf])-1e10
#     lon1_spock = np.zeros([nb_time_spock_netcdf])-1e10; lon1_netcdf = np.zeros([nb_time_spock_netcdf])-1e10
#     lat1_spock = np.zeros([nb_time_spock_netcdf])-1e10; lat1_netcdf = np.zeros([nb_time_spock_netcdf])-1e10

#     x2_spock = np.zeros([nb_time_spock_netcdf])-1e10; x2_netcdf = np.zeros([nb_time_spock_netcdf])-1e10
#     y2_spock = np.zeros([nb_time_spock_netcdf])-1e10; y2_netcdf = np.zeros([nb_time_spock_netcdf])-1e10
#     z2_spock = np.zeros([nb_time_spock_netcdf])-1e10; z2_netcdf = np.zeros([nb_time_spock_netcdf])-1e10
#     lon2_spock = np.zeros([nb_time_spock_netcdf])-1e10; lon2_netcdf = np.zeros([nb_time_spock_netcdf])-1e10
#     lat2_spock = np.zeros([nb_time_spock_netcdf])-1e10; lat2_netcdf = np.zeros([nb_time_spock_netcdf])-1e10

#     x3_spock = np.zeros([nb_time_spock_netcdf])-1e10; x3_netcdf = np.zeros([nb_time_spock_netcdf])-1e10
#     y3_spock = np.zeros([nb_time_spock_netcdf])-1e10; y3_netcdf = np.zeros([nb_time_spock_netcdf])-1e10
#     z3_spock = np.zeros([nb_time_spock_netcdf])-1e10; z3_netcdf = np.zeros([nb_time_spock_netcdf])-1e10
#     lon3_spock = np.zeros([nb_time_spock_netcdf])-1e10; lon3_netcdf = np.zeros([nb_time_spock_netcdf])-1e10
#     lat3_spock = np.zeros([nb_time_spock_netcdf])-1e10; lat3_netcdf = np.zeros([nb_time_spock_netcdf])-1e10

#     x4_spock = np.zeros([nb_time_spock_netcdf])-1e10; x4_netcdf = np.zeros([nb_time_spock_netcdf])-1e10
#     y4_spock = np.zeros([nb_time_spock_netcdf])-1e10; y4_netcdf = np.zeros([nb_time_spock_netcdf])-1e10
#     z4_spock = np.zeros([nb_time_spock_netcdf])-1e10; z4_netcdf = np.zeros([nb_time_spock_netcdf])-1e10
#     lon4_spock = np.zeros([nb_time_spock_netcdf])-1e10; lon4_netcdf = np.zeros([nb_time_spock_netcdf])-1e10
#     lat4_spock = np.zeros([nb_time_spock_netcdf])-1e10; lat4_netcdf = np.zeros([nb_time_spock_netcdf])-1e10


#     for itime in range(nb_time_spock_netcdf):
#         n_spec_both = len(comb_min_vs_time[itime])
#         x1_spock[itime] = x_spec_spock_same_time_as_netcdf[itime][comb_min_vs_time[itime][0][0]]
#         x1_netcdf[itime] = x_spec_netcdf[itime][comb_min_vs_time[itime][0][1]]
#         y1_spock[itime] = y_spec_spock_same_time_as_netcdf[itime][comb_min_vs_time[itime][0][0]]
#         y1_netcdf[itime] = y_spec_netcdf[itime][comb_min_vs_time[itime][0][1]]
#         z1_spock[itime] = z_spec_spock_same_time_as_netcdf[itime][comb_min_vs_time[itime][0][0]]
#         z1_netcdf[itime] = z_spec_netcdf[itime][comb_min_vs_time[itime][0][1]]
#         lon1_spock[itime] = lon_spock_same_time_as_netcdf[itime][comb_min_vs_time[itime][0][0]]
#         lon1_netcdf[itime] = lon_spec_netcdf[itime][comb_min_vs_time[itime][0][1]]
#         lat1_spock[itime] = lat_spock_same_time_as_netcdf[itime][comb_min_vs_time[itime][0][0]]
#         lat1_netcdf[itime] = lat_spec_netcdf[itime][comb_min_vs_time[itime][0][1]]


#         if n_spec_both > 1:
#             x2_spock[itime] = x_spec_spock_same_time_as_netcdf[itime][comb_min_vs_time[itime][1][0]]
#             x2_netcdf[itime] = x_spec_netcdf[itime][comb_min_vs_time[itime][1][1]]
#             y2_spock[itime] = y_spec_spock_same_time_as_netcdf[itime][comb_min_vs_time[itime][1][0]]
#             y2_netcdf[itime] = y_spec_netcdf[itime][comb_min_vs_time[itime][1][1]]
#             z2_spock[itime] = z_spec_spock_same_time_as_netcdf[itime][comb_min_vs_time[itime][1][0]]
#             z2_netcdf[itime] = z_spec_netcdf[itime][comb_min_vs_time[itime][1][1]]
#             lon2_spock[itime] = lon_spock_same_time_as_netcdf[itime][comb_min_vs_time[itime][1][0]]
#             lon2_netcdf[itime] = lon_spec_netcdf[itime][comb_min_vs_time[itime][1][1]]
#             lat2_spock[itime] = lat_spock_same_time_as_netcdf[itime][comb_min_vs_time[itime][1][0]]
#             lat2_netcdf[itime] = lat_spec_netcdf[itime][comb_min_vs_time[itime][1][1]]


#         if n_spec_both > 2:
#             x3_spock[itime] = x_spec_spock_same_time_as_netcdf[itime][comb_min_vs_time[itime][2][0]]
#             x3_netcdf[itime] = x_spec_netcdf[itime][comb_min_vs_time[itime][2][1]]
#             y3_spock[itime] = y_spec_spock_same_time_as_netcdf[itime][comb_min_vs_time[itime][2][0]]
#             y3_netcdf[itime] = y_spec_netcdf[itime][comb_min_vs_time[itime][2][1]]
#             z3_spock[itime] = z_spec_spock_same_time_as_netcdf[itime][comb_min_vs_time[itime][2][0]]
#             z3_netcdf[itime] = z_spec_netcdf[itime][comb_min_vs_time[itime][2][1]]
#             lon3_spock[itime] = lon_spock_same_time_as_netcdf[itime][comb_min_vs_time[itime][2][0]]
#             lon3_netcdf[itime] = lon_spec_netcdf[itime][comb_min_vs_time[itime][2][1]]
#             lat3_spock[itime] = lat_spock_same_time_as_netcdf[itime][comb_min_vs_time[itime][2][0]]
#             lat3_netcdf[itime] = lat_spec_netcdf[itime][comb_min_vs_time[itime][2][1]]


#         if n_spec_both > 3:
#             x4_spock[itime] = x_spec_spock_same_time_as_netcdf[itime][comb_min_vs_time[itime][3][0]]
#             x4_netcdf[itime] = x_spec_netcdf[itime][comb_min_vs_time[itime][3][1]]
#             y4_spock[itime] = y_spec_spock_same_time_as_netcdf[itime][comb_min_vs_time[itime][3][0]]
#             y4_netcdf[itime] = y_spec_netcdf[itime][comb_min_vs_time[itime][3][1]]
#             z4_spock[itime] = z_spec_spock_same_time_as_netcdf[itime][comb_min_vs_time[itime][3][0]]
#             z4_netcdf[itime] = z_spec_netcdf[itime][comb_min_vs_time[itime][3][1]]
#             lon4_spock[itime] = lon_spock_same_time_as_netcdf[itime][comb_min_vs_time[itime][3][0]]
#             lon4_netcdf[itime] = lon_spec_netcdf[itime][comb_min_vs_time[itime][3][1]]
#             lat4_spock[itime] = lat_spock_same_time_as_netcdf[itime][comb_min_vs_time[itime][3][0]]
#             lat4_netcdf[itime] = lat_spec_netcdf[itime][comb_min_vs_time[itime][3][1]]



 





#     pickle.dump(nb_seconds_since_initial_epoch_spock, open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_nb_seconds_since_initial_epoch_spock_againFM" + str(cygfm) + ".pickle", "w"))
#     pickle.dump(r_spec_diff_mag_min_to_max_array, open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_spec_diff_mag_min_to_max_array_againFM" + str(cygfm) + ".pickle", "w"))
#     pickle.dump(r_cyg_diff_mag, open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_cyg_diff_mag_againFM" + str(cygfm) + ".pickle", "w"))
#     pickle.dump(r_cyg_spock, open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_cyg_spock_againFM" + str(cygfm) + ".pickle", "w"))
#     pickle.dump(r_cyg_netcdf, open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_cyg_netcdf_againFM" + str(cygfm) + ".pickle", "w"))
#     #print date_start_val



# min_x1_spock = np.min(x1_spock[np.where(x1_spock != -1e10)[0]])
# min_y1_spock = np.min(y1_spock[np.where(y1_spock != -1e10)[0]])
# min_z1_spock = np.min(z1_spock[np.where(z1_spock != -1e10)[0]])
# min_x2_spock = np.min(x2_spock[np.where(x2_spock != -1e10)[0]])
# min_y2_spock = np.min(y2_spock[np.where(y2_spock != -1e10)[0]])
# min_z2_spock = np.min(z2_spock[np.where(z2_spock != -1e10)[0]])
# min_x3_spock = np.min(x3_spock[np.where(x3_spock != -1e10)[0]])
# min_y3_spock = np.min(y3_spock[np.where(y3_spock != -1e10)[0]])
# min_z3_spock = np.min(z3_spock[np.where(z3_spock != -1e10)[0]])
# min_x4_spock = np.min(x4_spock[np.where(x4_spock != -1e10)[0]])
# min_y4_spock = np.min(y4_spock[np.where(y4_spock != -1e10)[0]])
# min_z4_spock = np.min(z4_spock[np.where(z4_spock != -1e10)[0]])
# min_x1_netcdf = np.min(x1_netcdf[np.where(x1_netcdf != -1e10)[0]])
# min_y1_netcdf = np.min(y1_netcdf[np.where(y1_netcdf != -1e10)[0]])
# min_z1_netcdf = np.min(z1_netcdf[np.where(z1_netcdf != -1e10)[0]])
# min_x2_netcdf = np.min(x2_netcdf[np.where(x2_netcdf != -1e10)[0]])
# min_y2_netcdf = np.min(y2_netcdf[np.where(y2_netcdf != -1e10)[0]])
# min_z2_netcdf = np.min(z2_netcdf[np.where(z2_netcdf != -1e10)[0]])
# min_x3_netcdf = np.min(x3_netcdf[np.where(x3_netcdf != -1e10)[0]])
# min_y3_netcdf = np.min(y3_netcdf[np.where(y3_netcdf != -1e10)[0]])
# min_z3_netcdf = np.min(z3_netcdf[np.where(z3_netcdf != -1e10)[0]])
# min_x4_netcdf = np.min(x4_netcdf[np.where(x4_netcdf != -1e10)[0]])
# min_y4_netcdf = np.min(y4_netcdf[np.where(y4_netcdf != -1e10)[0]])
# min_z4_netcdf = np.min(z4_netcdf[np.where(z4_netcdf != -1e10)[0]])
# min_x = np.min([min_x1_spock, min_x2_spock, min_x3_spock, min_x4_spock, min_x1_netcdf, min_x2_netcdf, min_x3_netcdf, min_x4_netcdf])
# min_y = np.min([min_y1_spock, min_y2_spock, min_y3_spock, min_y4_spock, min_y1_netcdf, min_y2_netcdf, min_y3_netcdf, min_y4_netcdf])
# min_z = np.min([min_z1_spock, min_z2_spock, min_z3_spock, min_z4_spock, min_z1_netcdf, min_z2_netcdf, min_z3_netcdf, min_z4_netcdf])

# max_x1_spock = np.max(x1_spock[np.where(x1_spock != -1e10)[0]])
# max_y1_spock = np.max(y1_spock[np.where(y1_spock != -1e10)[0]])
# max_z1_spock = np.max(z1_spock[np.where(z1_spock != -1e10)[0]])
# max_x2_spock = np.max(x2_spock[np.where(x2_spock != -1e10)[0]])
# max_y2_spock = np.max(y2_spock[np.where(y2_spock != -1e10)[0]])
# max_z2_spock = np.max(z2_spock[np.where(z2_spock != -1e10)[0]])
# max_x3_spock = np.max(x3_spock[np.where(x3_spock != -1e10)[0]])
# max_y3_spock = np.max(y3_spock[np.where(y3_spock != -1e10)[0]])
# max_z3_spock = np.max(z3_spock[np.where(z3_spock != -1e10)[0]])
# max_x4_spock = np.max(x4_spock[np.where(x4_spock != -1e10)[0]])
# max_y4_spock = np.max(y4_spock[np.where(y4_spock != -1e10)[0]])
# max_z4_spock = np.max(z4_spock[np.where(z4_spock != -1e10)[0]])
# max_x1_netcdf = np.max(x1_netcdf[np.where(x1_netcdf != -1e10)[0]])
# max_y1_netcdf = np.max(y1_netcdf[np.where(y1_netcdf != -1e10)[0]])
# max_z1_netcdf = np.max(z1_netcdf[np.where(z1_netcdf != -1e10)[0]])
# max_x2_netcdf = np.max(x2_netcdf[np.where(x2_netcdf != -1e10)[0]])
# max_y2_netcdf = np.max(y2_netcdf[np.where(y2_netcdf != -1e10)[0]])
# max_z2_netcdf = np.max(z2_netcdf[np.where(z2_netcdf != -1e10)[0]])
# max_x3_netcdf = np.max(x3_netcdf[np.where(x3_netcdf != -1e10)[0]])
# max_y3_netcdf = np.max(y3_netcdf[np.where(y3_netcdf != -1e10)[0]])
# max_z3_netcdf = np.max(z3_netcdf[np.where(z3_netcdf != -1e10)[0]])
# max_x4_netcdf = np.max(x4_netcdf[np.where(x4_netcdf != -1e10)[0]])
# max_y4_netcdf = np.max(y4_netcdf[np.where(y4_netcdf != -1e10)[0]])
# max_z4_netcdf = np.max(z4_netcdf[np.where(z4_netcdf != -1e10)[0]])
# max_x = np.max([max_x1_spock, max_x2_spock, max_x3_spock, max_x4_spock, max_x1_netcdf, max_x2_netcdf, max_x3_netcdf, max_x4_netcdf])
# max_y = np.max([max_y1_spock, max_y2_spock, max_y3_spock, max_y4_spock, max_y1_netcdf, max_y2_netcdf, max_y3_netcdf, max_y4_netcdf])
# max_z = np.max([max_z1_spock, max_z2_spock, max_z3_spock, max_z4_spock, max_z1_netcdf, max_z2_netcdf, max_z3_netcdf, max_z4_netcdf])

time_diff_prn_arr = np.array(time_diff_prn)
which_spec_diff_time_diff_prn = which_spec_diff[time_diff_prn, :]
fom_netcdf_arr = np.array(fom_netcdf)
fom_netcdf_diff_prn_temp = fom_netcdf_arr[time_diff_prn, :]
fom_netcdf_diff_prn = []
fom_netcdf_diff_prn_all = []
order_gain_netcdf_diff_prn = np.zeros([4])
nb_spec_diff = np.zeros([len(time_diff_prn)])
nb_spec_diff_gain_not_0 = 0
for itime in range(len(time_diff_prn)):
    fom_netcdf_diff_prn_itime = []
    if 0 in which_spec_diff_time_diff_prn[itime, :]:
        order_gain_netcdf_diff_prn[0] = order_gain_netcdf_diff_prn[0] + 1
    if 1 in which_spec_diff_time_diff_prn[itime, :]:
        order_gain_netcdf_diff_prn[1] = order_gain_netcdf_diff_prn[1] + 1
    if 2 in which_spec_diff_time_diff_prn[itime, :]:
        order_gain_netcdf_diff_prn[2] = order_gain_netcdf_diff_prn[2] + 1
    if 3 in which_spec_diff_time_diff_prn[itime, :]:
        order_gain_netcdf_diff_prn[3] = order_gain_netcdf_diff_prn[3] + 1
    fom_netcdf_diff_prn.append( fom_netcdf_diff_prn_temp[itime, np.where(which_spec_diff_time_diff_prn[itime, :] != -1)[0] ] )
    for  idiff in range(len(np.where(which_spec_diff_time_diff_prn[itime, :] != -1)[0])):
        fom_netcdf_diff_prn_all.append( fom_netcdf_diff_prn_temp[itime, np.where(which_spec_diff_time_diff_prn[itime, :] != -1)[0][idiff] ] )
        if fom_netcdf_diff_prn_temp[itime, np.where(which_spec_diff_time_diff_prn[itime, :] != -1)[0][idiff] ] != 0:
            nb_spec_diff_gain_not_0 = nb_spec_diff_gain_not_0 +1 
    nb_spec_diff[itime] = len(np.where(which_spec_diff_time_diff_prn[itime, :] != -1)[0])

correct_vs_fom_percent = np.zeros([16])
for ifom in range(16):
    correct_vs_fom_percent[ifom] = correct_vs_fom[ifom, 0] * 100. / correct_vs_fom[ifom, 1]


# compare gain for a given phi and thetea to theoretical gain from antenna gain pattern

filename_gain_list = ['/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/src/ant_1_port_ddmi_v1.agm', '/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/src/ant_1_starboard_ddmi_v1.agm'] # list of antenna gain files

nb_file = len(filename_gain_list)
for ifile in range(nb_file):
    #ifile = 0
    filename_gain = filename_gain_list[ifile]
    file_gain = open(filename_gain, "r")
    read_file_gain = file_gain.readlines()
    if ifile == 0:
        nb_theta = len(read_file_gain)
        nb_phi = len(read_file_gain[0].split(','))
        theta_max = 90. 
        dtheta = (int)( theta_max/nb_theta )
        theta_arr = np.arange(0, theta_max, dtheta)
        phi_max = 360. 
        dphi = (int)( phi_max/nb_phi )
        phi_arr = np.arange(0, phi_max, dphi)
        theo_gain = np.zeros([nb_file, nb_theta, nb_phi])
    for itheta in range(nb_theta):
        for iphi in range(nb_phi):
            theo_gain[ifile, itheta, iphi] = np.float( read_file_gain[itheta].split(',')[iphi] )
max_gain = np.max(theo_gain)

# Out[31]: (31.268536, 42.866177, [44, 0])
theta = 35
phi = 107
itheta = np.where((theta_arr > theta))[0][0] - 1
iphi = np.where((phi_arr > phi))[0][0] - 1
print 'theoretical gain', np.max(theo_gain[:, itheta, iphi]), theo_gain[:, itheta, iphi]
print theo_gain[:, itheta-1:itheta+2, iphi-1:iphi+2]


time_where_diff_fom_average = 0
for idate in range(nb_date_ok):
    diff_fom_same_prn_vs_fom_time_arr = np.array(diff_fom_same_prn_vs_fom_time_all_date[idate]) 
    # diff_fom_same_prn_vs_fom_time.append( [itime, ispec, fom_netcdf[itime][ispec], gain_spock_same_time_as_netcdf[itime][where_prn_spock] - fom_netcdf[itime][ispec] ] )
    nb_spec_correct = diff_fom_same_prn_vs_fom_time_arr.shape[0]
    percentrage_fom_correct_vs_fom = np.zeros([16])
    time_where_diff_fom = []                        
    az_when_gain_diff = []
    el_when_gain_diff = []

    az_spock_when_gain_diff = []
    el_spock_when_gain_diff = []

    el_when_gain_diff_spock_and_diff_theory = []
    az_when_gain_diff_spock_and_diff_theory = []
    fom_netcdf_diff_theory = []
    fom_netcdf_diff_theory_even_corr= []
    fom_netcdf_diff_theory_even_other_corr= []
    fom_netcdf_diff_theory_even_any_neighbor = []
    fom_netcdf_diff_theory_not_in_other_but_in_neighbor = []
    fom_netcdf_diff_theory_near_edges = []
    nb_fom_correct_vs_fom = np.zeros([16])
    nb_fom_incorrect_vs_fom = np.zeros([16])
    nb_fom_tot_vs_fom = np.zeros([16])
    mean_absolute_fom_diff_vs_fom = np.zeros([16])

    padding_theta_up = 1
    padding_theta_down = 0
    padding_phi_up = 3
    padding_phi_down = 2
    for icorrect in range(nb_spec_correct):
        itime = (int)(diff_fom_same_prn_vs_fom_time_arr[icorrect, 0])
        ispec_netcdf = (int)(diff_fom_same_prn_vs_fom_time_arr[icorrect, 1])
        ifom = (int)(diff_fom_same_prn_vs_fom_time_arr[icorrect, 2])
        fom_diff = diff_fom_same_prn_vs_fom_time_arr[icorrect, 3]
        ispec_spock = (int)(diff_fom_same_prn_vs_fom_time_arr[icorrect, 4])
        nb_fom_tot_vs_fom[ifom] = nb_fom_tot_vs_fom[ifom] + 1
        if fom_diff == 0:
            nb_fom_correct_vs_fom[ifom] = nb_fom_correct_vs_fom[ifom] + 1
        else:
            mean_absolute_fom_diff_vs_fom[ifom] = mean_absolute_fom_diff_vs_fom[ifom] + np.abs(fom_diff)
            nb_fom_incorrect_vs_fom[ifom] = nb_fom_incorrect_vs_fom[ifom] + 1
#             az_when_gain_diff.append(az_spec_netcdf[itime][ispec_netcdf])
#             el_when_gain_diff.append(el_spec_netcdf[itime][ispec_netcdf])
            time_where_diff_fom.append([itime, ispec_netcdf])
    #         el_spock_when_gain_diff.append( el_spec_spock_same_time_as_netcdf[itime][ispec_spock] )
    #         az_spock_when_gain_diff.append( az_spec_spock_same_time_as_netcdf[itime][ispec_spock] )
    time_where_diff_fom_average = time_where_diff_fom_average + len(time_where_diff_fom) * 100./nb_spec_correct
time_where_diff_fom_average = time_where_diff_fom_average/nb_date_ok

    
#         # compare to theoretical gain from antenna gain pattern
#         if el_spec_netcdf[itime][ispec_netcdf] < theta_arr[-1]:
#             itheta = (int)(el_spec_netcdf[itime][ispec_netcdf]/5.)   # np.where((theta_arr > el_spec_netcdf[itime][ispec_netcdf]))[0][0] - 1
#             itheta_up = itheta + 1
#         else:
#             itheta = len(theta_arr) - 1
#             itheta_up = 0
#         if el_spec_netcdf[itime][ispec_netcdf] != 0:
#             itheta_down = itheta - 1
#         else:
#             itheta_down = len(theta_arr) - 1


#         if az_spec_netcdf[itime][ispec_netcdf] < phi_arr[-1]:
#             iphi = (int)(az_spec_netcdf[itime][ispec_netcdf]/15)  #np.where((phi_arr > az_spec_netcdf[itime][ispec_netcdf]))[0][0] - 1
#             iphi_up = iphi + 1
#         else:
#             iphi = len(phi_arr) - 1
#             iphi_up = 0
#         if az_spec_netcdf[itime][ispec_netcdf] != 0:
#             iphi_down = iphi - 1
#         else:
#             iphi_down = len(phi_arr) - 1

#         if np.mod(el_spec_netcdf[itime][ispec_netcdf], 5) >= 5 - padding_theta_up:
#             itheta_corr = itheta_up
#         elif np.mod(el_spec_netcdf[itime][ispec_netcdf], 5) <= padding_theta_down:
#             itheta_corr = itheta_down
#         else:
#             itheta_corr = itheta
#         if np.mod(az_spec_netcdf[itime][ispec_netcdf], 15) >= 15 - padding_phi_up:
#             iphi_corr = iphi_up
#         elif np.mod(az_spec_netcdf[itime][ispec_netcdf], 15) <= padding_phi_down:
#             iphi_corr = iphi_down
#         else:
#             iphi_corr = iphi
            
#         #print 'theoretical gain', np.max(theo_gain[:, itheta, iphi]), theo_gain[:, itheta, iphi]
#         if np.max(theo_gain[:, itheta, iphi]) != fom_netcdf[itime][ispec_netcdf]: # if netcdf RCG is different fro wha tit should be in theory
#             fom_netcdf_diff_theory.append([itime,  ispec_netcdf, theo_gain[:, itheta, iphi], fom_netcdf[itime][ispec_netcdf], el_spec_netcdf[itime][ispec_netcdf], az_spec_netcdf[itime][ispec_netcdf]])


#             az_when_gain_diff_spock_and_diff_theory.append(az_spec_netcdf[itime][ispec_netcdf])
#             el_when_gain_diff_spock_and_diff_theory.append(el_spec_netcdf[itime][ispec_netcdf])

#             # look only at the SPs are are near the edges
#             if ( (np.mod(el_spec_netcdf[itime][ispec_netcdf],5) >= 5 - padding_theta_up) | (np.mod(az_spec_netcdf[itime][ispec_netcdf],15) >= 15 - padding_phi_up) | (np.mod(az_spec_netcdf[itime][ispec_netcdf],15) <= padding_phi_down) | (np.mod(el_spec_netcdf[itime][ispec_netcdf],5) <= padding_theta_down) ):
#                 fom_netcdf_diff_theory_near_edges.append([itime,  ispec_netcdf, theo_gain[:, itheta, iphi], fom_netcdf[itime][ispec_netcdf], el_spec_netcdf[itime][ispec_netcdf], az_spec_netcdf[itime][ispec_netcdf]])

#                 if ( ( np.max(theo_gain[:, itheta_corr, iphi_corr]) != fom_netcdf[itime][ispec_netcdf]) & \
#                          ( np.max(theo_gain[:, itheta, iphi_corr]) != fom_netcdf[itime][ispec_netcdf]) & \
#                          ( np.max(theo_gain[:, itheta_corr, iphi]) != fom_netcdf[itime][ispec_netcdf]) ):                     
#                     fom_netcdf_diff_theory_even_corr.append([itime,  ispec_netcdf, theo_gain[:, itheta, iphi], fom_netcdf[itime][ispec_netcdf], el_spec_netcdf[itime][ispec_netcdf], az_spec_netcdf[itime][ispec_netcdf]])

#                     itheta_other_corr = (int)(round(el_spec_netcdf[itime][ispec_netcdf]/5.))
#                     iphi_other_corr = (int)(round(az_spec_netcdf[itime][ispec_netcdf]/15.))

#                     if ( ( np.max(theo_gain[:, itheta_other_corr, iphi_other_corr]) != fom_netcdf[itime][ispec_netcdf]) & \
#                              ( np.max(theo_gain[:, itheta, iphi_other_corr]) != fom_netcdf[itime][ispec_netcdf]) & \
#                              ( np.max(theo_gain[:, itheta_other_corr, iphi]) != fom_netcdf[itime][ispec_netcdf]) ):                     
#                         fom_netcdf_diff_theory_even_other_corr.append([itime,  ispec_netcdf, theo_gain[:, itheta, iphi], fom_netcdf[itime][ispec_netcdf], el_spec_netcdf[itime][ispec_netcdf], az_spec_netcdf[itime][ispec_netcdf]])

#                         if (  ( fom_netcdf[itime][ispec_netcdf] in theo_gain[:, itheta-1:itheta+2, iphi-1:iphi+2] ) == False ):
#                             fom_netcdf_diff_theory_even_any_neighbor.append([itime,  ispec_netcdf, theo_gain[:, itheta, iphi], fom_netcdf[itime][ispec_netcdf], el_spec_netcdf[itime][ispec_netcdf], az_spec_netcdf[itime][ispec_netcdf]])
#                         else:
#                             fom_netcdf_diff_theory_not_in_other_but_in_neighbor.append([itime,  ispec_netcdf, theo_gain[:, itheta, iphi], fom_netcdf[itime][ispec_netcdf], el_spec_netcdf[itime][ispec_netcdf], az_spec_netcdf[itime][ispec_netcdf]])

# #                         if itime == 3248: 
# #                             raise Exception
# #             else:
# #                 az_when_gain_diff_spock_and_diff_theory.append(az_spec_netcdf[itime][ispec_netcdf])
# #                 el_when_gain_diff_spock_and_diff_theory.append(el_spec_netcdf[itime][ispec_netcdf])

print len(fom_netcdf_diff_theory), len(time_where_diff_fom)
#print len(fom_netcdf_diff_theory_near_edges), len( fom_netcdf_diff_theory_even_corr), len(fom_netcdf_diff_theory_even_other_corr),  len(fom_netcdf_diff_theory_even_any_neighbor)





percentrage_fom_correct_vs_fom = nb_fom_correct_vs_fom * 100. / nb_fom_tot_vs_fom
mean_absolute_fom_diff_vs_fom = mean_absolute_fom_diff_vs_fom / nb_fom_incorrect_vs_fom

#a = np.where((np.mod(el_when_gain_diff,5) > 2 ) & ((np.mod(el_when_gain_diff,5) < 3 )))[0]

#ax.scatter(np.arange(0, len(np.mod(np.array(az_when_gain_diff)[a],15))), np.mod(np.array(az_when_gain_diff)[a],15))
#np.where( (np.mod(el_when_gain_diff,5) < 4) & (np.mod(az_when_gain_diff,15) < 13) & (np.mod(az_when_gain_diff,15) > 2) & (np.mod(el_when_gain_diff,5) > 1))[0][0]




fig, ax = plt.subplots()
plt.ion()
ax.scatter(np.mod(az_spock_when_gain_diff, 15), np.mod(el_spock_when_gain_diff,5))
ax.plot([15 - padding_phi_up,15 - padding_phi_up],[padding_theta_down,  5 - padding_theta_up], 'r')
ax.plot([padding_phi_down,15 - padding_phi_up],[ 5 - padding_theta_up,  5 - padding_theta_up], 'r')
ax.plot([15 - padding_phi_up,15 - padding_phi_up],[padding_theta_down,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([15 - padding_phi_up,15 - padding_phi_up],[padding_theta_down,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([15 - padding_phi_up,15 - padding_phi_up],[padding_theta_down,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([padding_phi_down,15 - padding_phi_up],[ 5 - padding_theta_up,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([padding_phi_down,padding_phi_down],[padding_theta_down,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([padding_phi_down,15 - padding_phi_up],[padding_theta_down, padding_theta_down], 'r', linewidth = 3)

ax.set_title('SIFT ' + format(len( np.where( (np.mod(az_spock_when_gain_diff, 15) >= 15 - padding_phi_up) | (np.mod(az_spock_when_gain_diff, 15) <= padding_phi_down) | (np.mod(el_spock_when_gain_diff, 5) >=  5 - padding_theta_up ) | (np.mod(el_spock_when_gain_diff, 5) <=  padding_theta_down )   )[0] ) * 100. / len(az_spock_when_gain_diff), ".0f") + '% near edges')


fig, ax = plt.subplots()
plt.ion()
ax.scatter(np.mod(az_when_gain_diff, 15), np.mod(el_when_gain_diff,5))
ax.plot([15 - padding_phi_up,15 - padding_phi_up],[padding_theta_down,  5 - padding_theta_up], 'r')
ax.plot([padding_phi_down,15 - padding_phi_up],[ 5 - padding_theta_up,  5 - padding_theta_up], 'r')
ax.plot([15 - padding_phi_up,15 - padding_phi_up],[padding_theta_down,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([15 - padding_phi_up,15 - padding_phi_up],[padding_theta_down,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([15 - padding_phi_up,15 - padding_phi_up],[padding_theta_down,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([padding_phi_down,15 - padding_phi_up],[ 5 - padding_theta_up,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([padding_phi_down,padding_phi_down],[padding_theta_down,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([padding_phi_down,15 - padding_phi_up],[padding_theta_down, padding_theta_down], 'r', linewidth = 3)

ax.set_title('Netcdf ' + format(len( np.where( (np.mod(az_when_gain_diff, 15) > 15 - padding_phi_up) | (np.mod(az_when_gain_diff, 15) < padding_phi_down) | (np.mod(el_when_gain_diff, 5) >  5 - padding_theta_up )| (np.mod(el_spock_when_gain_diff, 5) <=  padding_theta_down ))[0] ) * 100. / len(az_when_gain_diff), ".0f") + '% near edges')






fig, ax = plt.subplots()
plt.ion()
ax.scatter(np.mod(az_when_gain_diff_spock_and_diff_theory, 15), np.mod(el_when_gain_diff_spock_and_diff_theory,5))
ax.plot([15 - padding_phi_up,15 - padding_phi_up],[padding_theta_down,  5 - padding_theta_up], 'r')
ax.plot([padding_phi_down,15 - padding_phi_up],[ 5 - padding_theta_up,  5 - padding_theta_up], 'r')
ax.plot([15 - padding_phi_up,15 - padding_phi_up],[padding_theta_down,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([15 - padding_phi_up,15 - padding_phi_up],[padding_theta_down,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([15 - padding_phi_up,15 - padding_phi_up],[padding_theta_down,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([padding_phi_down,15 - padding_phi_up],[ 5 - padding_theta_up,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([padding_phi_down,padding_phi_down],[padding_theta_down,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([padding_phi_down,15 - padding_phi_up],[padding_theta_down, padding_theta_down], 'r', linewidth = 3)

ax.set_title('Theory ' + format(len( np.where( (np.mod(az_when_gain_diff_spock_and_diff_theory, 15) > 15 - padding_phi_up) | (np.mod(az_when_gain_diff_spock_and_diff_theory, 15) < padding_phi_down) | (np.mod(el_when_gain_diff_spock_and_diff_theory, 5) >  5 - padding_theta_up ) | (np.mod(el_spock_when_gain_diff, 5) <=  padding_theta_down ))[0] ) * 100. / len(az_when_gain_diff_spock_and_diff_theory), ".0f") + '% near edges')

#len( np.where( (np.mod(az_when_gain_diff_spock_and_diff_theory, 15) < 15 - padding_phi_up) & (np.mod(az_when_gain_diff_spock_and_diff_theory, 15) > padding_phi_down) & (np.mod(el_when_gain_diff_spock_and_diff_theory, 5) <  5 - padding_theta_up  ))[0] ) * 100. / len(az_when_gain_diff_spock_and_diff_theory)







# Look at ties (remember: fom is the another notation nor rcg so the'=y're used the same)
nb_time_at_least_one_tie = 0
nb_tie_per_rcg = np.zeros([16])
time_at_least_one_tie = []
prn_selected_by_onboard = []
#order_prn_selected_by_onboard = np.zeros([nb_time_spock_netcdf, nb_spec_spock]) - 1 # 2nd coord is the RCG that ties
highest_prn_not_selected_among_ties  = np.zeros([nb_time_spock_netcdf, nb_spec_spock]) - 1
at_least_one_but_not_all_spec_of_the_tie_got_selected = np.zeros([nb_time_spock_netcdf, nb_spec_spock]) - 1# -1 if no tie, 0 if a tie (keep track of fom that ties) but no SpOCK spec of the tie matches the netcdf PRNs,
                                                                                               # 1 if at least one of the SpOCK spec that tie matches ta netcdf PRN but not all ties -> called a "valid tie"
                                                                                               # 2 if all SpOCK spec that tie match the netcdf PRN 
                                                                                               # ex: SpOCK PRNs 2,23,17 tie at RCG = 13. 
                                                                                               # -> if none of the 3 PRN is among the entcdf PRN then 0
                                                                                               # -> if 1 or 2 of the 3 PRN is/are among the netcdf PRN then 1
                                                                                               # -> if all 3 match the netcdf PRN then 2
                                                                                               # -> we look at times when at least one of the PRNis sleected but there is a selection among all the PRN that tie (at least 1 is rejected)
                                                                                               # -> at_least_one_but_not_all_spec_of_the_tie_got_selected needs to be 1
nb_spec_of_tie_got_matched = np.zeros([nb_time_spock_netcdf, nb_spec_spock]) - 1
time_at_least_one_tie_gain_not_0  = []
time_valid_tie = [] # times for which at least one SpOCK PRN of the tie is selected in the netcdf but not all -> at least one SpOCK PRN of the tie is not selected. 
nb_fom_tie_this_time = []
for itime in range(nb_time_spock_netcdf):
    prn_selected_by_onboard_this_time = []
    which_fom_tie = np.where( nb_ties_per_rcg_per_time[itime, :]  >= 1 )[0]
    nb_fom_tie_this_time.append( len(which_fom_tie ) )
    if nb_fom_tie_this_time[-1] >= 1:# at leastone tie for this time
        nb_time_at_least_one_tie = nb_time_at_least_one_tie + 1
        time_at_least_one_tie.append(itime)
        already_valid_tie = 0
        for ifom in range(nb_fom_tie_this_time[-1]):
            nb_spec_of_tie_got_matched[itime, ifom] = 0
            prn_spock_tie_list= []
            ifom_tie = which_fom_tie[ifom]
            nb_tie_per_rcg[ifom_tie] = nb_tie_per_rcg[ifom_tie] + nb_ties_per_rcg_per_time[itime, ifom_tie]
            # find which (if any) of the tie is in netcdf
            where_spec_spock_tie = np.where(fom_spock[itime] == ifom_tie)[0]
            for ispec in range(len(where_spec_spock_tie)):
                prn_spock_here = gps_spock[itime][where_spec_spock_tie[ispec]]
                prn_spock_tie_list.append(prn_spock_here)
            prn_selected_by_onboard_this_time_sorted_higher_to_lower_prn = -np.sort(-np.array(prn_spock_tie_list))
            for ispec in range(len(prn_spock_tie_list)):
                if prn_spock_tie_list[ispec] in gps_netcdf[itime]:
                    nb_spec_of_tie_got_matched[itime, ifom]  = nb_spec_of_tie_got_matched[itime, ifom]  + 1
            if ( nb_spec_of_tie_got_matched[itime, ifom] == 0 ): # none of the SpOCK PRN of the tie matches a PRN in netcdf
                at_least_one_but_not_all_spec_of_the_tie_got_selected[itime, ifom] = 0
            elif (nb_spec_of_tie_got_matched[itime, ifom] == (nb_ties_per_rcg_per_time[itime, ifom_tie] + 1)): # + 1 because a tie corresponds to 2 same PRN . all SpOCK PRN of the tie are among netcdf PRNs
                at_least_one_but_not_all_spec_of_the_tie_got_selected[itime, ifom] = 2
            else:
                at_least_one_but_not_all_spec_of_the_tie_got_selected[itime, ifom] = 1
            if at_least_one_but_not_all_spec_of_the_tie_got_selected[itime, ifom] == 1: # at least ont SpOCK PRN is selected but not all -> at least one SpOCK PRN is rejected. Called a "valid tie"
                if already_valid_tie == 0:
                    time_valid_tie.append(itime)
                # look at if the highest PRN of the tie got selected
                if prn_selected_by_onboard_this_time_sorted_higher_to_lower_prn[0] in gps_netcdf[itime]: # the highest PRN is selected
                    highest_prn_not_selected_among_ties[itime, ifom] = 0
                else: # the highest PRN is not selected
                    highest_prn_not_selected_among_ties[itime, ifom] = 1
                already_valid_tie = 1

nb_time_valid_tie = len(time_valid_tie)

# look at which spec got selected among all spec that tied, and which one did not
valid_tie_got_matched = np.zeros([nb_time_valid_tie, np.max(nb_fom_tie_this_time), nb_spec_spock]) - 1 # 2nd coord: which FOM is the tie; 3rd coord: which SpOCK spec ties at this FOM
valid_tie_not_get_matched = np.zeros([nb_time_valid_tie, np.max(nb_fom_tie_this_time), nb_spec_spock]) - 1 # 2nd coord: which FOM is the tie; 3rd coord: which SpOCK spec ties at this FOM
for itime_count in range(nb_time_valid_tie):
    itime = time_valid_tie[itime_count]
    prn_selected_by_onboard_this_time = []
    which_fom_tie = np.where( nb_ties_per_rcg_per_time[itime, :]  >= 1 )[0]
    for ifom in range(nb_fom_tie_this_time[itime]):
        prn_spock_tie_list= []
        ifom_tie = which_fom_tie[ifom]
        where_spec_spock_tie = np.where(fom_spock[itime] == ifom_tie)[0]
        for ispec in range(len(where_spec_spock_tie)):
            prn_spock_here = gps_spock[itime][where_spec_spock_tie[ispec]]
            prn_spock_tie_list.append(prn_spock_here)
        ispec_count_match = 0
        ispec_count_not_match = 0
        for ispec in range(len(prn_spock_tie_list)):
            which_spock_spec_number = np.where(gps_spock[itime] == prn_spock_tie_list[ispec])[0]
            if prn_spock_tie_list[ispec] in gps_netcdf[itime]: 
                valid_tie_got_matched[itime_count, ifom, ispec_count_match] = which_spock_spec_number
                ispec_count_match = ispec_count_match + 1
            else:
                valid_tie_not_get_matched[itime_count, ifom, ispec_count_not_match] = which_spock_spec_number
                ispec_count_not_match = ispec_count_not_match + 1

# Now look at properties of the spec: azimuth, elevation, PRN, ...
az_valid_tie_got_matched = []
az_valid_tie_not_get_matched = []
el_valid_tie_got_matched = []
el_valid_tie_not_get_matched = []
el_gps_valid_tie_got_matched = []
el_gps_valid_tie_not_get_matched = []

prn_valid_tie_got_matched = []
prn_valid_tie_not_get_matched = []
for itime_count in range(nb_time_valid_tie):
    az_valid_tie_got_matched_this_time = []
    az_valid_tie_not_get_matched_this_time = []
    el_valid_tie_got_matched_this_time = []
    el_valid_tie_not_get_matched_this_time = []
    el_gps_valid_tie_got_matched_this_time = []
    el_gps_valid_tie_not_get_matched_this_time = []

    prn_valid_tie_got_matched_this_time = []
    prn_valid_tie_not_get_matched_this_time = []

    itime = time_valid_tie[itime_count]
    for ifom in range(nb_fom_tie_this_time[itime]):
        az_valid_tie_got_matched_this_time_this_fom = []
        az_valid_tie_not_get_matched_this_time_this_fom = []
        el_valid_tie_got_matched_this_time_this_fom = []
        el_valid_tie_not_get_matched_this_time_this_fom = []
        el_gps_valid_tie_got_matched_this_time_this_fom = []
        el_gps_valid_tie_not_get_matched_this_time_this_fom = []

        prn_valid_tie_got_matched_this_time_this_fom = []
        prn_valid_tie_not_get_matched_this_time_this_fom = []
        nb_match = len( np.where( valid_tie_got_matched[itime_count, ifom, :] != -1 )[0] )
        nb_no_match = len( np.where( valid_tie_not_get_matched[itime_count, ifom, :] != -1 )[0] )
        if ( ( nb_match > 0 ) & ( nb_no_match > 0 ) ): # valie tie for this fom. (at this time there's at least one valie tie but there could be another FOM that is not a valid tie, in which case ignore it
            # spec got matched in netcdf 
            for ispec_count in range(nb_match):
                ispec = (int)(valid_tie_got_matched[itime_count, ifom, ispec_count])
                az_valid_tie_got_matched_this_time_this_fom.append(az_spec_spock_same_time_as_netcdf[itime, ispec])
                el_valid_tie_got_matched_this_time_this_fom.append(el_spec_spock_same_time_as_netcdf[itime, ispec])
                el_gps_valid_tie_got_matched_this_time_this_fom.append(el_gps_from_cyg_spock_same_time_as_netcdf[itime, ispec])
                prn_valid_tie_got_matched_this_time_this_fom.append(gps_spock[itime, ispec])
            # spec did not get matched in netcdf
            for ispec_count in range(nb_no_match):
                ispec = (int)(valid_tie_not_get_matched[itime_count, ifom, ispec_count])
                az_valid_tie_not_get_matched_this_time_this_fom.append(az_spec_spock_same_time_as_netcdf[itime, ispec])
                el_valid_tie_not_get_matched_this_time_this_fom.append(el_spec_spock_same_time_as_netcdf[itime, ispec])
                el_gps_valid_tie_not_get_matched_this_time_this_fom.append(el_gps_from_cyg_spock_same_time_as_netcdf[itime, ispec])
                prn_valid_tie_not_get_matched_this_time_this_fom.append(gps_spock[itime, ispec])

            az_valid_tie_got_matched_this_time.append(az_valid_tie_got_matched_this_time_this_fom)
            az_valid_tie_not_get_matched_this_time.append(az_valid_tie_not_get_matched_this_time_this_fom)
            el_valid_tie_got_matched_this_time.append(el_valid_tie_got_matched_this_time_this_fom)
            el_valid_tie_not_get_matched_this_time.append(el_valid_tie_not_get_matched_this_time_this_fom)
            el_gps_valid_tie_got_matched_this_time.append(el_gps_valid_tie_got_matched_this_time_this_fom)
            el_gps_valid_tie_not_get_matched_this_time.append(el_gps_valid_tie_not_get_matched_this_time_this_fom)

            prn_valid_tie_got_matched_this_time.append(prn_valid_tie_got_matched_this_time_this_fom)
            prn_valid_tie_not_get_matched_this_time.append(prn_valid_tie_not_get_matched_this_time_this_fom)

    az_valid_tie_got_matched.append(az_valid_tie_got_matched_this_time)
    az_valid_tie_not_get_matched.append(az_valid_tie_not_get_matched_this_time)
    el_valid_tie_got_matched.append(el_valid_tie_got_matched_this_time)
    el_valid_tie_not_get_matched.append(el_valid_tie_not_get_matched_this_time)
    el_gps_valid_tie_got_matched.append(el_gps_valid_tie_got_matched_this_time)
    el_gps_valid_tie_not_get_matched.append(el_gps_valid_tie_not_get_matched_this_time)

    prn_valid_tie_got_matched.append(prn_valid_tie_got_matched_this_time)
    prn_valid_tie_not_get_matched.append(prn_valid_tie_not_get_matched_this_time)




# to be read if nb_spec_spock >= 6: look at if, for a simu with 4 SpOCK spec, when a PRN is not iamong netcdf PRn, it is because it's a tie. assumes the pickles were saved from a simu with 4 spec
# pickle.dump(which_spec_diff,open("which_spec_diff_4sp_2018FM" + str(cygfm) + ".pickle", "w"))
# pickle.dump(time_diff_prn,open("time_diff_prn_4sp_2018FM" + str(cygfm) + ".pickle", "w"))
# pickle.dump(fom_spock,open("fom_spock_4sp_2018FM" + str(cygfm) + ".pickle", "w"))
# pickle.dump(gps_spock,open("gps_spock_4sp_2018FM" + str(cygfm) + ".pickle", "w"))
# pickle.dump(az_spec_not_int_spock_same_time_as_netcdf,open("az_spec_not_int_spock_same_time_as_netcdf_4sp_2018_test.pickle", "w"))
# pickle.dump(el_spec_not_int_spock_same_time_as_netcdf,open("el_spec_not_int_spock_same_time_as_netcdf_4sp_2018_test.pickle", "w"))
# pickle.dump(el_gps_from_cyg_spock_same_time_as_netcdf, open("el_gps_from_cyg_spock_same_time_as_netcdf_4sp_2018_test.pickle", "w"))



#PRESENTATION 051518 here is where I look at ties. It works but you might need to adapt the names of the pickles. It requires this section of the code to be run with nb_spec_spock = 6. These pickles are loaded from simulations with nb_spec_spock = 4.
which_spec_diff_4_spec = pickle.load(open("which_spec_diff_4sp_2018_test.pickle"))
time_diff_prn_4_spec = pickle.load(open("time_diff_prn_4sp_2018_test.pickle"))
fom_spock_4_spec = pickle.load(open("fom_spock_4sp_2018_test.pickle"))
gps_spock_4_spec = pickle.load(open("gps_spock_4sp_2018_test.pickle"))
az_spec_not_int_spock_same_time_as_netcdf_4_spec = pickle.load(open("az_spec_not_int_spock_same_time_as_netcdf_4sp_2018_test.pickle"))
el_spec_not_int_spock_same_time_as_netcdf_4_spec = pickle.load(open("el_spec_not_int_spock_same_time_as_netcdf_4sp_2018_test.pickle"))
el_gps_from_cyg_spock_same_time_as_netcdf_4_spec = pickle.load(open("el_gps_from_cyg_spock_same_time_as_netcdf_4sp_2018_test.pickle"))

diff_prn_is_not_tie = 0
diff_prn_is_tie = 0
az_wrong_prn_and_not_tie = []
el_wrong_prn_and_not_tie = []
el_gps_wrong_prn_and_not_tie = []
diff_prn_is_not_tie_small_fom = 0
diff_prn_is_tie_small_fom = 0
az_wrong_prn_and_not_tie_small_fom = []
el_wrong_prn_and_not_tie_small_fom = []
el_gps_wrong_prn_and_not_tie_small_fom = []
fom_wrong_prn_and_not_tie = np.zeros([16])
fom_wrong_prn = np.zeros([16])
fom_wrong_prn_tie = np.zeros([16])
padding_theta_up = 1
padding_theta_down = 1
padding_phi_up = 1
padding_phi_down = 1
prn_wrong = []
prn_wrong_arr = np.zeros([33])
for itime_temp in range(len(time_diff_prn_4_spec)):
    #itime_temp = 10
    itime = time_diff_prn_4_spec[itime_temp]
    for ispec in range(4):
        if ( ( gps_spock_4_spec[itime][ispec] in gps_netcdf[itime] ) == False ): 
            prn_wrong.append(gps_spock_4_spec[itime][ispec])
            prn_wrong_arr[gps_spock_4_spec[itime][ispec]] = prn_wrong_arr[gps_spock_4_spec[itime][ispec]] + 1
            fom = (int)(fom_spock_4_spec[itime][ispec])
            if fom < 4:
                nb_spec_have_this_small_fom = len(np.where( fom_spock[itime] == fom )[0])
                if nb_spec_have_this_small_fom >= 2:
                    diff_prn_is_tie_small_fom = diff_prn_is_tie_small_fom + 1
                else:
                    diff_prn_is_not_tie_small_fom = diff_prn_is_not_tie_small_fom + 1
                    az_wrong_prn_and_not_tie_small_fom.append(np.mod(az_spec_not_int_spock_same_time_as_netcdf_4_spec[itime][ispec], 15))
                    el_wrong_prn_and_not_tie_small_fom.append(np.mod(el_spec_not_int_spock_same_time_as_netcdf_4_spec[itime][ispec], 5))
                    el_gps_wrong_prn_and_not_tie_small_fom.append(el_gps_from_cyg_spock_same_time_as_netcdf_4_spec[itime][ispec])

            else:
                nb_spec_have_this_fom = len(np.where( fom_spock[itime] == fom )[0])
                fom_wrong_prn[fom] = fom_wrong_prn[fom] + 1
                if nb_spec_have_this_fom >= 2:
                    diff_prn_is_tie = diff_prn_is_tie + 1
                    fom_wrong_prn_tie[fom] = fom_wrong_prn_tie[fom] + 1
                else:
#                     if (gps_spock_4_spec[itime][ispec] == 18):
#                         raise Exception
                    diff_prn_is_not_tie = diff_prn_is_not_tie + 1
                    az_wrong_prn_and_not_tie.append(np.mod(az_spec_not_int_spock_same_time_as_netcdf_4_spec[itime][ispec], 15))
                    el_wrong_prn_and_not_tie.append(np.mod(el_spec_not_int_spock_same_time_as_netcdf_4_spec[itime][ispec], 5))
                    el_gps_wrong_prn_and_not_tie.append(el_gps_from_cyg_spock_same_time_as_netcdf_4_spec[itime][ispec])
                    fom_wrong_prn_and_not_tie[fom] = fom_wrong_prn_and_not_tie[fom] + 1
#                     if ((np.mod(az_wrong_prn_and_not_tie[-1], 15) < 15 - padding_phi_up) & (np.mod(az_wrong_prn_and_not_tie[-1], 15) > padding_phi_down) & (np.mod(el_wrong_prn_and_not_tie[-1], 5) <  5 - padding_theta_up ) & (np.mod(el_wrong_prn_and_not_tie[-1], 5) >  padding_theta_down ) ):
#                         raise Exception
#                     if diff_prn_is_not_tie > 300:
#                         raise Exception
#end of PRESENTATION 051518 


print format(len( np.where( (np.mod(az_wrong_prn_and_not_tie, 15) >= 15 - padding_phi_up) | (np.mod(az_wrong_prn_and_not_tie, 15) <= padding_phi_down) | (np.mod(el_wrong_prn_and_not_tie, 5) >=  5 - padding_theta_up ) | (np.mod(el_wrong_prn_and_not_tie, 5) <=  padding_theta_down )   )[0] ) * 100. / len(az_wrong_prn_and_not_tie), ".0f") + '% near edges'
fig, ax = plt.subplots()
plt.ion()
ax.scatter(np.mod(az_wrong_prn_and_not_tie, 15), np.mod(el_wrong_prn_and_not_tie,5))
ax.plot([15 - padding_phi_up,15 - padding_phi_up],[padding_theta_down,  5 - padding_theta_up], 'r')
ax.plot([padding_phi_down,15 - padding_phi_up],[ 5 - padding_theta_up,  5 - padding_theta_up], 'r')
ax.plot([15 - padding_phi_up,15 - padding_phi_up],[padding_theta_down,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([15 - padding_phi_up,15 - padding_phi_up],[padding_theta_down,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([15 - padding_phi_up,15 - padding_phi_up],[padding_theta_down,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([padding_phi_down,15 - padding_phi_up],[ 5 - padding_theta_up,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([padding_phi_down,padding_phi_down],[padding_theta_down,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([padding_phi_down,15 - padding_phi_up],[padding_theta_down, padding_theta_down], 'r', linewidth = 3)

ax.set_title('SIFT ' + format(len( np.where( (np.mod(az_wrong_prn_and_not_tie, 15) >= 15 - padding_phi_up) | (np.mod(az_wrong_prn_and_not_tie, 15) <= padding_phi_down) | (np.mod(el_wrong_prn_and_not_tie, 5) >=  5 - padding_theta_up ) | (np.mod(el_wrong_prn_and_not_tie, 5) <=  padding_theta_down )   )[0] ) * 100. / len(az_wrong_prn_and_not_tie), ".0f") + '% near edges')

fig, ax = plt.subplots()
plt.ion()
ax.scatter(np.mod(az_wrong_prn_and_not_tie_small_fom, 15), np.mod(el_wrong_prn_and_not_tie_small_fom,5))
ax.plot([15 - padding_phi_up,15 - padding_phi_up],[padding_theta_down,  5 - padding_theta_up], 'r')
ax.plot([padding_phi_down,15 - padding_phi_up],[ 5 - padding_theta_up,  5 - padding_theta_up], 'r')
ax.plot([15 - padding_phi_up,15 - padding_phi_up],[padding_theta_down,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([15 - padding_phi_up,15 - padding_phi_up],[padding_theta_down,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([15 - padding_phi_up,15 - padding_phi_up],[padding_theta_down,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([padding_phi_down,15 - padding_phi_up],[ 5 - padding_theta_up,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([padding_phi_down,padding_phi_down],[padding_theta_down,  5 - padding_theta_up], 'r', linewidth = 3)
ax.plot([padding_phi_down,15 - padding_phi_up],[padding_theta_down, padding_theta_down], 'r', linewidth = 3)
plt.show();plt.show();

#ax.set_title('SIFT ' + format(len( np.where( (np.mod(az_wrong_prn_and_not_tie_small_fom, 15) >= 15 - padding_phi_up) | (np.mod(az_wrong_prn_and_not_tie_small_fom, 15) <= padding_phi_down) | (np.mod(el_wrong_prn_and_not_tie_small_fom, 5) >=  5 - padding_theta_up ) | (np.mod(el_wrong_prn_and_not_tie_small_fom, 5) <=  padding_theta_down )   )[0] ) * 100. / len(az_wrong_prn_and_not_tie_small_fom), ".0f") + '% near edges')




height_fig = 11.  # the width is calculated as height_fig * 4/3.                                                                             
fontsize_plot = 25      
ratio_fig_size = 4./3
fig_title = ''#Accuracy VS RCG
y_label = 'Azimuth'
x_label = 'Time'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                      
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                       
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           
for itime_count in range(nb_time_valid_tie):
    itime = time_valid_tie[itime_count]
    ifom_count = 0
    for ifom in range(nb_fom_tie_this_time[itime]):
        nb_match = len( np.where( valid_tie_got_matched[itime_count, ifom, :] != -1 )[0] )
        nb_no_match = len( np.where( valid_tie_not_get_matched[itime_count, ifom, :] != -1 )[0] )
        if ( ( nb_match > 0 ) & ( nb_no_match > 0 ) ): # valie tie for this fom. (at this time there's at least one valie tie but there could be another FOM that is not a valid tie, in which case ignore it
            for imatch in range(nb_match):
                for inomatch in range(nb_no_match):
                    ax.scatter(itime_count, az_valid_tie_got_matched[itime_count][ifom_count][imatch] - az_valid_tie_not_get_matched[itime_count][ifom_count][inomatch], color = 'b')

            ifom_count = ifom_count + 1
ax.margins(0,0)
fig_save_name = date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_") + '_tie_azimuth.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')





height_fig = 11.  # the width is calculated as height_fig * 4/3.                                                                             
fontsize_plot = 25      
ratio_fig_size = 4./3
fig_title = ''#Accuracy VS RCG
y_label = 'Elevation'
x_label = 'Time'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                      
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                       
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           
nb_tie_match_elev_higher_than_no_match = 0
nb_tie_match_elev_lower_than_no_match = 0
nb_tie_match_elev_equal_no_match = 0
for itime_count in range(nb_time_valid_tie):
    itime = time_valid_tie[itime_count]
    ifom_count = 0
    for ifom in range(nb_fom_tie_this_time[itime]):
        nb_match = len( np.where( valid_tie_got_matched[itime_count, ifom, :] != -1 )[0] )
        nb_no_match = len( np.where( valid_tie_not_get_matched[itime_count, ifom, :] != -1 )[0] )
        if ( ( nb_match > 0 ) & ( nb_no_match > 0 ) ): # valie tie for this fom. (at this time there's at least one valie tie but there could be another FOM that is not a valid tie, in which case ignore it
            for imatch in range(nb_match):
                for inomatch in range(nb_no_match):
                    ax.scatter(itime_count, el_valid_tie_got_matched[itime_count][ifom_count][imatch] - el_valid_tie_not_get_matched[itime_count][ifom_count][inomatch], color = 'b')
                    if  el_valid_tie_got_matched[itime_count][ifom_count][imatch] - el_valid_tie_not_get_matched[itime_count][ifom_count][inomatch] < 0:
                        nb_tie_match_elev_higher_than_no_match = nb_tie_match_elev_higher_than_no_match + 1
                    elif el_valid_tie_got_matched[itime_count][ifom_count][imatch] - el_valid_tie_not_get_matched[itime_count][ifom_count][inomatch] > 0:
                        nb_tie_match_elev_lower_than_no_match = nb_tie_match_elev_lower_than_no_match + 1
                    else:
                        nb_tie_match_elev_equal_no_match = nb_tie_match_elev_equal_no_match + 1 
            ifom_count = ifom_count + 1
ax.margins(0,0)
fig_save_name = date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_") + '_tie_elevation.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')
print nb_tie_match_elev_higher_than_no_match, nb_tie_match_elev_lower_than_no_match, nb_tie_match_elev_equal_no_match







height_fig = 11.  # the width is calculated as height_fig * 4/3.                                                                             
fontsize_plot = 25      
ratio_fig_size = 4./3
fig_title = ''#Accuracy VS RCG
y_label = 'Elevation GPS'
x_label = 'Time'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                      
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                       
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           
nb_tie_match_elev_gps_higher_than_no_match = 0
nb_tie_match_elev_gps_lower_than_no_match = 0
nb_tie_match_elev_gps_equal_no_match = 0
for itime_count in range(nb_time_valid_tie):
    itime = time_valid_tie[itime_count]
    ifom_count = 0
    for ifom in range(nb_fom_tie_this_time[itime]):
        nb_match = len( np.where( valid_tie_got_matched[itime_count, ifom, :] != -1 )[0] )
        nb_no_match = len( np.where( valid_tie_not_get_matched[itime_count, ifom, :] != -1 )[0] )
        if ( ( nb_match > 0 ) & ( nb_no_match > 0 ) ): # valie tie for this fom. (at this time there's at least one valie tie but there could be another FOM that is not a valid tie, in which case ignore it
            for imatch in range(nb_match):
                for inomatch in range(nb_no_match):
                    ax.scatter(itime_count, el_gps_valid_tie_got_matched[itime_count][ifom_count][imatch] - el_gps_valid_tie_not_get_matched[itime_count][ifom_count][inomatch], color = 'b')
                    if  el_gps_valid_tie_got_matched[itime_count][ifom_count][imatch] - el_gps_valid_tie_not_get_matched[itime_count][ifom_count][inomatch] < 0:
                        nb_tie_match_elev_gps_higher_than_no_match = nb_tie_match_elev_gps_higher_than_no_match + 1
                    elif el_gps_valid_tie_got_matched[itime_count][ifom_count][imatch] - el_gps_valid_tie_not_get_matched[itime_count][ifom_count][inomatch] > 0:
                        nb_tie_match_elev_gps_lower_than_no_match = nb_tie_match_elev_gps_lower_than_no_match + 1
                    else:
                        nb_tie_match_elev_gps_equal_no_match = nb_tie_match_elev_gps_equal_no_match + 1 

            ifom_count = ifom_count + 1
ax.margins(0,0)
fig_save_name = date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_") + '_tie_elevation_gps.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')
print nb_tie_match_elev_gps_higher_than_no_match, nb_tie_match_elev_gps_lower_than_no_match, nb_tie_match_elev_gps_equal_no_match








height_fig = 11.  # the width is calculated as height_fig * 4/3.                                                                             
fontsize_plot = 25      
ratio_fig_size = 4./3
fig_title = ''#Accuracy VS RCG
y_label = 'PRN'
x_label = 'Time'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                      
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                       
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           
for itime_count in range(nb_time_valid_tie):
    itime = time_valid_tie[itime_count]
    ifom_count = 0
    for ifom in range(nb_fom_tie_this_time[itime]):
        nb_match = len( np.where( valid_tie_got_matched[itime_count, ifom, :] != -1 )[0] )
        nb_no_match = len( np.where( valid_tie_not_get_matched[itime_count, ifom, :] != -1 )[0] )
        if ( ( nb_match > 0 ) & ( nb_no_match > 0 ) ): # valie tie for this fom. (at this time there's at least one valie tie but there could be another FOM that is not a valid tie, in which case ignore it
            for imatch in range(nb_match):
                for inomatch in range(nb_no_match):
                    ax.scatter(itime_count, prn_valid_tie_got_matched[itime_count][ifom_count][imatch] - prn_valid_tie_not_get_matched[itime_count][ifom_count][inomatch], color = 'b')
            ifom_count = ifom_count + 1
ax.margins(0,0)
fig_save_name = date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_") + '_tie_prn.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')



height_fig = 11.  # the width is calculated as height_fig * 4/3.                                                                             
fontsize_plot = 25      
ratio_fig_size = 4./3
fig_title = ''#Accuracy VS RCG
y_label = 'Elevation'
x_label = 'Time'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                      
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                       
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           
nb_tie_match_elev_higher_than_no_match = 0
nb_tie_match_elev_lower_than_no_match = 0
nb_tie_match_elev_equal_no_match = 0
for itime_count in range(nb_time_valid_tie):
    itime = time_valid_tie[itime_count]
    ifom_count = 0
    for ifom in range(nb_fom_tie_this_time[itime]):
        nb_match = len( np.where( valid_tie_got_matched[itime_count, ifom, :] != -1 )[0] )
        nb_no_match = len( np.where( valid_tie_not_get_matched[itime_count, ifom, :] != -1 )[0] )
        if ( ( nb_match > 0 ) & ( nb_no_match > 0 ) ): # valie tie for this fom. (at this time there's at least one valie tie but there could be another FOM that is not a valid tie, in which case ignore it
            for imatch in range(nb_match):
                for inomatch in range(nb_no_match):
                    ax.scatter(itime_count, np.abs(np.mod(az_valid_tie_got_matched[itime_count][ifom_count][imatch],15)-7.5) -  np.abs(np.mod(az_valid_tie_not_get_matched[itime_count][ifom_count][inomatch],15)-7.5), color = 'b')
            ifom_count = ifom_count + 1
ax.margins(0,0)
fig_save_name = date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_") + '_tie_azim_dist_center.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')



# In [242]: el_when_gain_diff_spock_and_diff_theory[540]
# Out[242]: 37.15166

# In [243]: az_when_gain_diff_spock_and_diff_theory[540]
# Out[243]: 126.469894

# In [246]: fom_netcdf_diff_theory[540]
# Out[246]: [869, 3, array([6., 0.]), 7, 37.15166, 126.469894]
# -> should be 6 but is 7, but there's no 7 in proximity at all!

print nb_spec_spock_got_matched_vs_fom * 100. / nb_spec_spock_vs_fom
print nb_spec_spock_got_matched_vs_fom_remove_ties * 100. / nb_spec_spock_vs_fom_remove_ties
raise Exception


#             if np.max(theo_gain[:, itheta_up, iphi]) == fom_netcdf[itime][ispec_netcdf]:
#                 fom_netcdf_same_theory_theta_up.append([itime,  ispec_netcdf, theo_gain[:, itheta, iphi], fom_netcdf[itime][ispec_netcdf]])
#             elif np.max(theo_gain[:, itheta, iphi_up]) == fom_netcdf[itime][ispec_netcdf]:
#                 fom_netcdf_same_theory_phi_up.append([itime,  ispec_netcdf, theo_gain[:, itheta, iphi], fom_netcdf[itime][ispec_netcdf]])
#             elif np.max(theo_gain[:, itheta_down, iphi]) == fom_netcdf[itime][ispec_netcdf]:
#                 fom_netcdf_same_theory_theta_down.append([itime,  ispec_netcdf, theo_gain[:, itheta, iphi], fom_netcdf[itime][ispec_netcdf]])
#             elif np.max(theo_gain[:, itheta, iphi_down]) == fom_netcdf[itime][ispec_netcdf]:
#                 fom_netcdf_same_theory_phi_down.append([itime,  ispec_netcdf, theo_gain[:, itheta, iphi], fom_netcdf[itime][ispec_netcdf]])
#             else:
#                 fom_netcdf_diff_theory_even_up_or_down.append([itime,  ispec_netcdf, theo_gain[:, itheta, iphi], fom_netcdf[itime][ispec_netcdf]])

# In [315]: el_when_gain_diff[9]
# Out[315]: 31.438923

# In [316]: az_when_gain_diff[9]
# Out[316]: 311.3164


# In [31]: el_when_gain_diff[31], az_when_gain_diff[31], time_where_diff_fom[31]
# Out[31]: (31.268536, 42.866177, [44, 0])


itime_here = 58
print gps_netcdf[itime_here]
print fom_netcdf[itime_here]
print fom_spock[itime_here]
print gps_spock[itime_here]








############
############
###########
nb_seconds_since_initial_epoch_spock_all = []
r_spec_diff_mag_min_to_max_array_all = []
r_cyg_diff_mag_all = []
r_cyg_spock_all = []
r_cyg_netcdf_all = []
nb_time_spock_netcdf_all = []
date_start_val_save = date_start_val
#nb_date_ok = 1# !!!!!!!!! remove this line
for idate in range(nb_date_ok):
    #print idate, nb_date_ok-1
    #date_start_val = datetime.strftime(date_start_val_array[idate], "%Y-%m-%dT%H:%M:%S")#"2017-08-20T00:00:00" # date to start the validation
    date_start_val = date_start_val #'2017-06-02T00:15:00' # !!!!!! shoule be line right above and not this line
    nb_seconds_since_initial_epoch_spock_all.append( pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_nb_seconds_since_initial_epoch_spock_again_test.pickle")))
    r_spec_diff_mag_min_to_max_array_all.append( pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_spec_diff_mag_min_to_max_array_again_test.pickle")))
    r_cyg_diff_mag_all.append( pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_cyg_diff_mag_again_test.pickle")))
    r_cyg_spock_all.append( pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_cyg_spock_again_test.pickle")))
    r_cyg_netcdf_all.append( pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_cyg_netcdf_again_test.pickle")))
    nb_time_spock_netcdf_all.append( len(r_cyg_netcdf_all[-1]))


# nb_seconds_since_initial_epoch_spock = pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_nb_seconds_since_initial_epoch_spock_test.pickle"))
# r_spec_diff_mag_min_to_max_array = pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_spec_diff_mag_min_to_max_array_test.pickle"))
# r_cyg_diff_mag = pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_cyg_diff_mag_test.pickle"))
# r_cyg_spock = pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_cyg_spock_test.pickle"))
# r_cyg_netcdf = pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_cyg_netcdf_test.pickle"))
# nb_time_spock_netcdf = len(r_cyg_netcdf)




# if befre you run the script cygnss_spec_netcdf_to_spock_format.py, the following pickle was saved. it includes the dates to ignore to look at the spec position. Indeed those dates were missing dates in the netcdf 
idate = 0
date_ignore_missing_netcdf_data_all_sc = pickle.load(open("date_ignore_missing_netcdf_data_all_sc_" + date_start_ok_str.replace(":","_") + '_to_' + date_stop_ok_str.replace(":","_")	 + "_test.pickle", "r"))
which_sc = (int)(cygnss_name_to_norad_id('FM0' + str(cygfm) ) ) - 41884
date_ignore_missing_netcdf_data = date_ignore_missing_netcdf_data_all_sc[isc]

nb_time_to_remove_because_missing_netcdf_data = 0
for itime in range(nb_time_flight_rounded):
    if date_flight_rounded[itime] in date_ignore_missing_netcdf_data:
        nb_time_to_remove_because_missing_netcdf_data = nb_time_to_remove_because_missing_netcdf_data + 1
    

# print len(np.where(which_normpower_rank_is_biggest_dist_spock_spec == 0)[0]) * 100./ len(which_normpower_rank_is_biggest_dist_spock_spec)
# print len(np.where(which_normpower_rank_is_biggest_dist_spock_spec == 1)[0]) * 100./ len(which_normpower_rank_is_biggest_dist_spock_spec)

# print len(np.where(which_normpower_rank_is_biggest_dist_netcdf_spec == 0)[0]) * 100./ len(which_normpower_rank_is_biggest_dist_netcdf_spec)
# print len(np.where(which_normpower_rank_is_biggest_dist_netcdf_spec == 1)[0]) * 100./ len(which_normpower_rank_is_biggest_dist_netcdf_spec)




prn_spock_netcdf_different_vs_time_not_combin =  np.zeros([nb_time_flight_rounded])
for itime in range(nb_time_flight_rounded):
    nb_spec_here = len(gps_spock_same_time_as_netcdf[itime])
    for ispec in range(nb_spec_here):
        spock_spec_here = (int)(gps_spock_same_time_as_netcdf[itime][ispec].replace('PRN_',''))
        if ((spock_spec_here in gps_netcdf[itime] ) == False):
            prn_spock_netcdf_different_vs_time_not_combin[itime] = prn_spock_netcdf_different_vs_time_not_combin[itime] + 1
        

# prn_spock_netcdf_different_vs_time = np.zeros([nb_time_flight_rounded ,4])
# prn_spock_netcdf_different = np.zeros([ 4])
# for itime in range(nb_time_flight_rounded):
#     for ispec in range(len(prn_comb_min_vs_time[itime])):
#         prn_spock_here = prn_comb_min_vs_time[itime][ispec][0]
#         prn_netcdf_here = prn_comb_min_vs_time[itime][ispec][1]
#         if prn_spock_here != prn_netcdf_here:
#             if ((ispec != 2 ) | (len(prn_comb_min_vs_time[itime]) != 3)):
#                 prn_spock_netcdf_different_vs_time[itime, ispec] = prn_spock_netcdf_different_vs_time[itime, ispec]  + 1
#                 prn_spock_netcdf_different[ispec] = prn_spock_netcdf_different[ispec]  + 1
            

# PRESENTATION 051518
height_fig = 11.  # the width is calculated as height_fig * 4/3.                                                                             
fontsize_plot = 25      
ratio_fig_size = 4./3
fig_title = ''#Accuracy VS RCG
y_label = 'Selection accuracy (%)'
x_label = 'RCG'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                      
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                       
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           
x_here = np.arange(0, 16)
#whereok = np.where(np.isnan(correct_vs_fom_percent) == False)[0] # netcdf never has FOM = 1 for some reason...
#ax.plot(x_here[whereok], correct_vs_fom_percent_all_date_average[whereok], linewidth = 2, color = 'k', label = '3')

whereok = np.where(np.isnan(average_percentage_nb_spec_spock_got_matched_vs_fom_all_date_arr) == False)[0] # netcdf never has FOM = 1 for some reason...
ax.plot(x_here[whereok], average_percentage_nb_spec_spock_got_matched_vs_fom_all_date_arr[whereok], linewidth = 2, color = 'b', label = '3')
ax.plot([x_here[whereok][0], 4], [92,92], linestyle = 'dashed', linewidth = 2, color = 'k')
ax.plot([4, x_here[whereok][-1]], [92,92], linestyle = 'dashed', linewidth = 2, color = 'r')
ax.plot([4, 4], [0,92], linestyle = 'dashed', linewidth = 2, color = 'k')
ax.plot([4, 4], [92,100], linestyle = 'dashed', linewidth = 2, color = 'r')
percentage_spec_gain_ge_4 = np.sum(np.sum(np.array(nb_spec_spock_vs_fom_all_date), axis = 0)[4:])*100./np.sum(np.sum(np.array(nb_spec_spock_vs_fom_all_date), axis = 0))
percentage_spec_gain_le_3 = 100-percentage_spec_gain_ge_4
percentage_spec_gain_0 = np.sum(np.sum(np.array(nb_spec_spock_vs_fom_all_date), axis = 0)[0])*100./np.sum(np.sum(np.array(nb_spec_spock_vs_fom_all_date), axis = 0))
ax.text(19./2, 96, format(percentage_spec_gain_ge_4, ".0f") + '% of SPs RCG ' + '$\geq$' +' 4', color = 'r', fontsize = fontsize_plot, horizontalalignment = 'center', verticalalignment = 'top')
ax.text(2, 46, format(percentage_spec_gain_le_3, ".0f") + '% of SPs\nRCG < 4', color = 'k', fontsize = fontsize_plot, horizontalalignment = 'center', verticalalignment = 'center')

# correct_vs_fom_percent_all_date_average
# whereok = np.where(np.isnan(correct_vs_fom_percent_all_date_average) == False)[0] # netcdf never has FOM = 1 for some reason...
# ax.plot(x_here[whereok], correct_vs_fom_percent_all_date_average[whereok], linewidth = 2, color = 'b', label = '3')

ax.set_ylim([0,102])
ax.yaxis.set_ticks(np.arange(0,110,10))
ax.xaxis.set_ticks(np.arange(0,16))

ax.margins(0,0)
fig_save_name = date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_") + '_accuracy_vs_rcg.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')
# end of PRESENTATION 051518

height_fig = 11.  # the width is calculated as height_fig * 4/3.                                                                             
fontsize_plot = 25      
ratio_fig_size = 4./3
fig_title = ''#Number of SPs per gain
y_label = '# SPs'
x_label = 'RCG'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                      
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                       
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           
ax.plot(np.arange(0,16), nb_spec_per_gain_all_average, linewidth = 2, color = 'k', label = '3')
ax.scatter/Users/cbv/Go(np.arange(0,16), nb_spec_per_gain_all_average, linewidth = 2, color = 'k', label = '3')
#ax.set_ylim([0,100])
ax.margins(0,0)
fig_save_name = date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_") + '_nb_spec_per_rcg.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')




height_fig = 11.  # the width is calculated as height_fig * 4/3.                                                                                                                                                                                                                                                                                                                           
fontsize_plot = 25      
ratio_fig_size = 4./3

fig_title = ''#Difference in specular point position SpOCK vs L1' # L1 is netcdf                                                                                                                                                                                                                                         
y_label = 'Error (km)'
x_label = 'Prediction time (days)'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                                                                                                                                                                                                                                                                      
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                                                                                                                                                                                                                                                                         
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           

count_gt_threshold = 0
for idate in range(nb_date_ok):
    r_spec_diff_mag_min_to_max_array = r_spec_diff_mag_min_to_max_array_all[idate]

    if r_spec_diff_mag_min_to_max_array[-1,0] < 60:
        nb_seconds_since_initial_epoch_spock = nb_seconds_since_initial_epoch_spock_all[idate]
        x_axis = nb_seconds_since_initial_epoch_spock / 3600. / 24.
        ax.plot(x_axis, r_spec_diff_mag_min_to_max_array[:,0], linewidth = 2, color = 'b', label = '1')
        print idate,  r_spec_diff_mag_min_to_max_array[-1,0]
    #ax.plot(x_axis, r_spec_diff_mag_min_to_max_array[:,1], linewidth = 2, color = 'r', label = '2')
    #ax.plot(x_axis, r_spec_diff_mag_min_to_max_array[:,2], linewidth = 2, color = 'k', label = '3')
    #ax.plot(x_axis, r_spec_diff_mag_min_to_max_array[:,3], linewidth = 2, color = 'g', label = '4')  
    else:
        count_gt_threshold = count_gt_threshold + 1

    
ax.margins(0,0)
ax.set_ylim([0,10])
ax.set_xlim([0,3])

# legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
# legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_") + '_all_dates_difference_spec_position_spock_vs_netcdf_new.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')



#PRESENTATION 051518
height_fig = 11.  # the width is calculated as height_fig * 4/3.                                                                             
fontsize_plot = 25      
ratio_fig_size = 4./3
fig_title = ''#Difference in specular point position SpOCK vs L1' # L1 is netcdf                                                                            
x_label = 'Prediction time (days)'
for ispec in range(2,3):

    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                                                                                                                                                                                                                                      
    gs = gridspec.GridSpec(1, 2)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01, wspace = 0.5)
    count_gt_threshold = 0
    #for idate in range(nb_date_ok):
    ax = fig.add_subplot(gs[0, 0])
    y_label = 'Along-track difference (km)'
    ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                                                                                                                                                                                                                                                                         
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           

    for idate in range(nb_date_ok):
        diff_here = r_spec_diff_lvlh_same_prn_all_date[idate][:,ispec, 0]
        where_same_prn = np.where(diff_here > -1e29)[0] # =1e30 when prn SpOCK different from prn netcdf                                                    
        ax.scatter(nb_seconds_since_initial_epoch_spock_all_date[idate][where_same_prn]/ 3600. / 24., r_spec_diff_lvlh_same_prn_all_date[idate][where_same_prn,ispec,0], marker = '+', s = 5, color = 'k', label = 'Along-track', alpha = 0.4)
        # !!!!!!! reomve line below and uncomment line above      
        #ax.scatter(np.arange(0, len(diff_here), 1) / (len(diff_here) - 1) * x_axis_day[-1] ,diff_here, marker = '+', s = 5, color = 'k', label = 'Along-track', alpha = 0.4)
        #ax.plot(x_axis_day, median_r_spec_diff_lvlh_all_date[idate][:,ispec,0], linewidth = 3, color = 'k', label = '')
        #ax.plot(x_axis_day, q10_r_spec_diff_lvlh_all_date[idate][:,ispec,0], linewidth = 3, color = 'r', label = '', linestyle = 'dashed')
        #ax.plot(x_axis_day, q25_r_spec_diff_lvlh_all_date[idate][:,ispec,0], linewidth = 3, color = 'b', label = '', linestyle = 'dashed')
        #ax.plot(x_axis_day, q75_r_spec_diff_lvlh_all_date[idate][:,ispec,0], linewidth = 3, color = 'b', label = '', linestyle = 'dashed')
        #ax.plot(x_axis_day, q90_r_spec_diff_lvlh_all_date[idate][:,ispec,0], linewidth = 3, color = 'r', label = '', linestyle = 'dashed')

    ax.margins(0,0)
    ax.xaxis.set_ticks(np.arange(0,5))
    #ax.set_ylim([np.min(q25_r_spec_diff_lvlh[:,ispec,0]), np.max(q75_r_spec_diff_lvlh[:,ispec,0])])
    #ax.set_ylim([np.min(q10_r_spec_diff_lvlh[:,ispec,0]), np.max(q90_r_spec_diff_lvlh[:,ispec,0])])

    #ax.text(0.02, 0.98, 'spec ' + str(ispec), fontsize = fontsize_plot, transform = ax.transAxes, verticalalignment = 'top' )
    ax = fig.add_subplot(gs[0, 1])
    y_label = 'Cross-track difference (km)'
    ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                                                                                                                                                                                                                                                                         
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           

    for idate in range(nb_date_ok):
        diff_here = r_spec_diff_lvlh_same_prn_all_date[idate][:,ispec, 1]
        where_same_prn = np.where(diff_here > -1e29)[0] # =1e30 when prn SpOCK different from prn netcdf                                                    
        ax.scatter(nb_seconds_since_initial_epoch_spock_all_date[idate][where_same_prn]/ 3600. / 24., r_spec_diff_lvlh_same_prn_all_date[idate][where_same_prn,ispec,1], marker = '+', s = 5, color = 'k', label = 'Along-track', alpha = 0.4)
        # !!!!!!! reomve line below and uncomment line above
        #ax.scatter(np.arange(0, len(r_spec_diff_lvlh_same_prn_all_date[idate][:,ispec,1]), 1) / (len(r_spec_diff_lvlh_same_prn_all_date[idate][:,ispec,0]) - 1) *x_axis_day[-1] ,r_spec_diff_lvlh_same_prn_all_date[idate][:,ispec,1], marker = '+', s = 5, color = 'k', label = 'Along-track', alpha = 0.4)

        #ax.plot(x_axis_day, median_r_spec_diff_lvlh_all_date[idate][:,ispec,1], linewidth = 3, color = 'k', label = '')
        #ax.plot(x_axis_day, q10_r_spec_diff_lvlh_all_date[idate][:,ispec,1], linewidth = 3, color = 'r', label = '', linestyle = 'dashed')
        #ax.plot(x_axis_day, q25_r_spec_diff_lvlh_all_date[idate][:,ispec,1], linewidth = 3, color = 'b', label = '', linestyle = 'dashed')
        #ax.plot(x_axis_day, q75_r_spec_diff_lvlh_all_date[idate][:,ispec,1], linewidth = 3, color = 'b', label = '', linestyle = 'dashed')
        #ax.plot(x_axis_day, q90_r_spec_diff_lvlh_all_date[idate][:,ispec,1], linewidth = 3, color = 'r', label = '', linestyle = 'dashed')


    ax.margins(0,0)
    ax.set_ylim([-20,20])
    ax.xaxis.set_ticks(np.arange(0,5))
#     ax.set_ylim([np.min(q10_r_spec_diff_lvlh[:,ispec,1]), np.max(q90_r_spec_diff_lvlh[:,ispec,1])])
#     # legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
#     # legend.get_title().set_fontsize(str(fontsize_plot))

    fig_save_name = date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_") +  'specular_lvlh_error_ispec_' + str(ispec) + '.pdf'
    fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')
#end of PRESENTATION 051518


fig_title = ''#Error on atellitepositoin in LVLH frame of refecnce (SpOCK referece)
y_label = 'Error (km)'
x_label = 'Prediction time (days)'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                                                                                                                                                                                                                                                                      
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                                                                                                                                                                                                                                                                         
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           

count_gt_threshold = 0
idate = 0

x_axis = nb_seconds_since_initial_epoch_spock_dt_output / 3600. / 24.
ax.plot(x_axis, r_cyg_lvlh_diff[:,0], linewidth = 2, color = 'b', label = 'Along-track')
ax.plot(x_axis, r_cyg_lvlh_diff[:,1], linewidth = 2, color = 'm', label = 'Cross-track')
ax.plot(x_axis, r_cyg_lvlh_diff[:,2], linewidth = 2, color = 'r', label = 'Radial')

    
ax.margins(0,0)
#ax.set_ylim([0,70])
#ax.set_xlim([0,5])

legend = ax.legend(loc='upper left', bbox_to_anchor=(0, 1), numpoints = 1,  title="", fontsize = fontsize_plot)
#legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = 'satellite_lvlh_error.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')










# # median and quantile over a day for a particular date
# nb_day_simu = (int)(np.ceil((nb_seconds_since_initial_epoch_spock[-1] - nb_seconds_since_initial_epoch_spock[0] )/ 3600. / 24.))
# median_r_spec_diff_lvlh = np.zeros([nb_day_simu, 4, 3])
# q10_r_spec_diff_lvlh = np.zeros([nb_day_simu, 4, 3])
# q25_r_spec_diff_lvlh = np.zeros([nb_day_simu, 4, 3])
# q75_r_spec_diff_lvlh = np.zeros([nb_day_simu, 4, 3])
# q90_r_spec_diff_lvlh = np.zeros([nb_day_simu, 4, 3])
# x_axis_day = np.arange(0.5, nb_day_simu+0.5)+nb_seconds_since_initial_epoch_spock[0] / 3600. / 24. # start at 0.5 and go up to nb_day_simu + 0.5 because say that the median and quantiles are over a day so take it at the middle of the day
# r_spec_diff_lvlh_same_prn = r_spec_diff_lvlh_same_prn# r_spec_diff_min_to_max_array_lvlh # !!!!!!!! np.abs(r_spec_diff_min_to_max_array_lvlh)

# iday = 0
# itime_previous = 0
# itime_day_list = [] # list of itime for that falls every day
# itime_day_list.append(itime_previous)
# while iday < nb_day_simu:
#     itime = 0
#     while (int)( ( nb_seconds_since_initial_epoch_spock[itime+itime_previous] - nb_seconds_since_initial_epoch_spock[itime_previous] )/ 3600. / 24.) < 1:
#         itime = itime + 1
#         if itime+itime_previous == r_spec_diff_lvlh_same_prn.shape[0]:
#             break
#     median_r_spec_diff_lvlh[iday, :, :] = np.median(r_spec_diff_lvlh_same_prn[itime_previous:itime+itime_previous, :, :], axis = 0)
#     q10_r_spec_diff_lvlh[iday, :, :] = np.percentile(r_spec_diff_lvlh_same_prn[itime_previous:itime+itime_previous, :, :], 10, axis = 0)
#     q25_r_spec_diff_lvlh[iday, :, :] = np.percentile(r_spec_diff_lvlh_same_prn[itime_previous:itime+itime_previous, :, :], 25, axis = 0)
#     q75_r_spec_diff_lvlh[iday, :, :] = np.percentile(r_spec_diff_lvlh_same_prn[itime_previous:itime+itime_previous, :, :], 75, axis = 0)
#     q90_r_spec_diff_lvlh[iday, :, :] = np.percentile(r_spec_diff_lvlh_same_prn[itime_previous:itime+itime_previous, :, :], 90, axis = 0)
#     #print itime_previous, itime+itime_previous, median_r_spec_diff_lvlh[iday, 0, 0] 
#     iday = iday + 1
#     itime_previous = itime + itime_previous
#     itime_day_list.append(itime_previous)


#PRESENTATION 051518
height_fig = 15.  # the width is calculated as height_fig * 4/3. 
fontsize_plot = 25
ratio_fig_size = 4./3
fig_title = ''#VLLH distribution specular points for a particular ispec at a particular time
y_label = 'Percentage (%)'
itime_day = 4 # which day to look at the distirbution over (distibution oof all spec over this day). from 1 (not 0) to nb_day (not nb_day-1). ex: if itime_day = 2 then it'll show the distribution from day 1 to day 2
coordtype = ['Along-track', 'Cross-track']
#for itime_day in range(nb_day):
for idate in range(1):#!!!!!nb_date_ok):
    for icoord in range(0,2):
        #icoord = 0 # 0 for along-track, 1 for cross-track
        x_label = coordtype[icoord] + ' difference (km)'
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
        plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                                                                                                                  
        gs = gridspec.GridSpec(2, 2)
        gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.17)
        for ispec in range(0, 4):
            diff_here = r_spec_diff_lvlh_same_prn_all_date[idate][itime_day_list_all_date[idate][itime_day - 1]:itime_day_list_all_date[idate][itime_day],ispec, icoord]
            if ispec == 0:
                irow = 0; icol = 0
            if ispec == 1:
                irow = 0; icol = 1
            if ispec == 2:
                irow = 1; icol = 0
            if ispec == 3:
                irow = 1; icol = 1
            ax = fig.add_subplot(gs[irow, icol])
            if ispec > 1:
                ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
            if ( ( ispec == 0 ) | ( ispec == 2 ) ):
                ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
            [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                                                                                                                       
            ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
            plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           
            #where_same_nb_spec_spock_netcdf = np.where(diff_here > -1e6)[0] # 1e6 when spock has 3 or 4 spec and netcdf has 4 or 3 spec
            where_same_prn = np.where(diff_here > -1e29)[0] # =1e30 when prn SpOCK different from prn netcdf
            #print np.min(diff_here[where_same_prn]), np.max(diff_here[where_same_prn]),len(where_same_prn) * 100. / len(diff_here)r_spec_diff_lvlh_same_prn_all_date[idate]
            range_max_temp = q75_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord]#q75_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord]*1.5 #np.max(r_spec_diff_lvlh_same_prn_all_date[idate][where_same_nb_spec_spock_netcdf, ispec, icoord])
            range_min_temp = median_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord] - (q75_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord] - median_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord] ) 
                                                                    #q25_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord]#q25_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord]*0.5 
            #         half_width = (range_max_temp - range_min_temp)/2 #(q75_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord] - q25_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord])/2.
            #         range_min  = np.min(diff_here[where_same_prn])#range_min_temp - half_width*3
            #         range_max  = np.max(diff_here[where_same_prn])#range_max_temp + half_width*3
            half_width = (q75_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord] - q25_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord])/2.
            range_min  = q10_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord]-1#range_min_temp - half_width*3
            range_max  = q90_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord]+1#range_max_temp + half_width*3

            hist_along_data = np.histogram(diff_here, range = [range_min, range_max])
            bin_array_temp = hist_along_data[1]
            bin_array = ( bin_array_temp[:-1] + np.roll(bin_array_temp,-1)[:-1] ) /2.
            binsize_actual = bin_array[1] - bin_array[0]
            hist_along = hist_along_data[0] * 100. / len(diff_here)
            ax.bar(bin_array, hist_along, binsize_actual)
    #         ax.plot([median_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord], median_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord]], [0, 100], linestyle = 'dashed', linewidth = 2, color = 'k')
    #         ax.plot([q25_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord], q25_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord]], [0, 100], linestyle = 'dashed', linewidth = 2, color = 'b')
    #         ax.plot([q75_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord], q75_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord]], [0, 100], linestyle = 'dashed', linewidth = 2, color = 'b')
    #         ax.plot([q10_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord], q10_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord]], [0, 100], linestyle = 'dashed', linewidth = 2, color = 'r')
    #         ax.plot([q90_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord], q90_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord]], [0, 100], linestyle = 'dashed', linewidth = 2, color = 'r')
            ax.set_ylim([np.min(hist_along), np.max(hist_along)])
            ax.set_xlim([np.min(bin_array_temp), np.max(bin_array_temp)])
            ax.text(0.02, 0.98, '50% width:\n' +r'$\pm$' + format(half_width, ".1f") + ' km', fontsize = fontsize_plot, transform = ax.transAxes, verticalalignment = 'top' )
            if ispec == 0:
                ax.text(0.98, 0.98, 'spec ' + str(ispec + 1) + '\nday' + str(itime_day), fontsize = fontsize_plot, transform = ax.transAxes, verticalalignment = 'top', horizontalalignment = 'right' )
            else:
                ax.text(0.98, 0.98, 'spec ' + str(ispec + 1), fontsize = fontsize_plot, transform = ax.transAxes, verticalalignment = 'top', horizontalalignment = 'right' )
        fig.set_figheight(height_fig)
        fig.set_figwidth(height_fig*ratio_fig_size)
        fig_save_name = 'distribution_specular_lvlh_error_ispec_' + coordtype[icoord].replace("-","").lower() + '_day_' + str(itime_day) +'_idate_' + str(idate) +  '_piston.pdf'

        fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')
# end of PRESENTATION 051518

len(np.where(diff_here < q10_r_spec_diff_lvlh_all_date[idate][itime_day - 1, ispec,icoord])[0]) * 100. / len(diff_here)

len(np.where(diff_here < 10.7149526383351)[0]) * 100. / len(diff_here)
disthere= 80
print len(np.where(np.abs(r_spec_diff_min_to_max_array_lvlh[itime_day_list_all_date[idate][itime_day - 1]:itime_day_list_all_date[idate][itime_day],ispec,icoord]) < disthere )[0]) * 100. / len(r_spec_diff_min_to_max_array_lvlh[itime_day_list_all_date[idate][itime_day - 1]:itime_day_list_all_date[idate][itime_day],ispec,icoord])

#####
height_fig = 11.  # the width is calculated as height_fig * 4/3.                                                                                                                                                                                                             
fig_title = 'ECEF specular point position SpOCK and L1' # L1 is netcdf                                                                                                                                                                                                
x_label = 'Prediction time (days)'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                                                                                                                                                                                                                                                                      
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)

ax = fig.add_subplot(gs[0, 0])
y_label = 'X (km)'
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
#ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                                                                                                                                                                                                                                                                         
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           

count_gt_threshold = 0
idate = 0
r_spec_diff_mag_min_to_max_array = r_spec_diff_mag_min_to_max_array_all[idate]

if r_spec_diff_mag_min_to_max_array[-1,0] < 60:
    nb_seconds_since_initial_epoch_spock = nb_seconds_since_initial_epoch_spock_all[idate]
    x_axis = nb_seconds_since_initial_epoch_spock / 3600. /24
    ax.plot(x_axis, x4_spock, linewidth = 2, color = 'b', label = '1')
    ax.scatter(x_axis, x4_netcdf, linewidth = 2, color = 'r', label = '1')

    ax.plot(x_axis, y4_spock, linewidth = 2, color = 'b', label = '1')
    ax.scatter(x_axis, y4_netcdf, linewidth = 2, color = 'r', label = '1')

    ax.plot(x_axis, z4_spock, linewidth = 2, color = 'b', label = '1')
    ax.scatter(x_axis, z4_netcdf, linewidth = 2, color = 'r', label = '1')

    print idate,  r_spec_diff_mag_min_to_max_array[-1,0]
    #ax.plot(x_axis, r_spec_diff_mag_min_to_max_array[:,1], linewidth = 2, color = 'r', label = '2')
    #ax.plot(x_axis, r_spec_diff_mag_min_to_max_array[:,2], linewidth = 2, color = 'k', label = '3')
    #ax.plot(x_axis, r_spec_diff_mag_min_to_max_array[:,3], linewidth = 2, color = 'g', label = '4')  
else:
    count_gt_threshold = count_gt_threshold + 1
ax.set_ylim([min_x,max_x])
ax.set_xlim([1/24., 2.5/24.])
fig_save_name = 'spec_position_component_spock_netcdf.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')


ax = fig.add_subplot(gs[1, 0])
y_label = 'Y (km)'
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
#ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                                                                                                                                                                                                                                                                         
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           

count_gt_threshold = 0
idate = 0
r_spec_diff_mag_min_to_max_array = r_spec_diff_mag_min_to_max_array_all[idate]

if r_spec_diff_mag_min_to_max_array[-1,0] < 60:
    nb_seconds_since_initial_epoch_spock = nb_seconds_since_initial_epoch_spock_all[idate]
    x_axis = nb_seconds_since_initial_epoch_spock / 3600. /24
    ax.plot(x_axis, y1_spock, linewidth = 2, color = 'b', label = '1')
    ax.scatter(x_axis, y1_netcdf, linewidth = 2, color = 'r', label = '1')
    print idate,  r_spec_diff_mag_min_to_max_array[-1,0]
    #ax.plot(x_axis, r_spec_diff_mag_min_to_max_array[:,1], linewidth = 2, color = 'r', label = '2')
    #ax.plot(x_axis, r_spec_diff_mag_min_to_max_array[:,2], linewidth = 2, color = 'k', label = '3')
    #ax.plot(x_axis, r_spec_diff_mag_min_to_max_array[:,3], linewidth = 2, color = 'g', label = '4')  
else:
    count_gt_threshold = count_gt_threshold + 1
ax.set_ylim([min_y,max_y])

ax = fig.add_subplot(gs[2, 0])
y_label = 'Z (km)'
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                                                                                                                                                                                                                                                                         
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           

count_gt_threshold = 0
idate = 0
r_spec_diff_mag_min_to_max_array = r_spec_diff_mag_min_to_max_array_all[idate]

if r_spec_diff_mag_min_to_max_array[-1,0] < 60:
    nb_seconds_since_initial_epoch_spock = nb_seconds_since_initial_epoch_spock_all[idate]
    x_axis = nb_seconds_since_initial_epoch_spock / 3600. /24
    ax.plot(x_axis, z1_spock, linewidth = 2, color = 'b', label = '1')
    ax.scatter(x_axis, z1_netcdf, linewidth = 2, color = 'r', label = '1')
    print idate,  r_spec_diff_mag_min_to_max_array[-1,0]
    #ax.plot(x_axis, r_spec_diff_mag_min_to_max_array[:,1], linewidth = 2, color = 'r', label = '2')
    #ax.plot(x_axis, r_spec_diff_mag_min_to_max_array[:,2], linewidth = 2, color = 'k', label = '3')
    #ax.plot(x_axis, r_spec_diff_mag_min_to_max_array[:,3], linewidth = 2, color = 'g', label = '4')  
else:
    count_gt_threshold = count_gt_threshold + 1
ax.set_ylim([min_z,max_z])


fig_save_name = 'spec_position_component_spock_netcdf.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')


raise Exception


# sc position difference
fig_title = 'Difference in CYGNSS position SpOCK vs L1' # L1 is netcdf                                                                                                                                                                                                                                         
y_label = 'Difference (km)'
x_label = 'Prediction time (days)'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                                                                                                                                                                                                                                                                      
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                                                                                                                                                                                                                                                                         
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                                                                                                                                                                                                                                                                     
x_axis = nb_seconds_since_initial_epoch_spock / 3600. / 24
ax.plot(x_axis, r_cyg_diff_mag, linewidth = 2, color = 'b', label = '1')
ax.margins(0,0)
ax.set_ylim([0,80])

# legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
# legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + '_difference_sat_position_spock_vs_netcdf_new.pdf'
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
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold
x_axis = nb_seconds_since_initial_epoch_spock / 3600. / 24
ax.plot(x_axis, r_cyg_spock_mag - r_cyg_netcdf_mag, linewidth = 2, color = 'b', label = '1')
ax.margins(0,0)
#ax.set_ylim([0,80])

# legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
# legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + '_difference_sat_radius_spock_vs_netcdf_new.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')



# Visu on a 2d map. 



fig_title = '' # L1 is netcdf           
y_label = ''
x_label = ''
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
font = {'family' : 'normal',
        'size'   : fontsize_plot}

plt.rc('font', **font) ## make the labels of the ticks in bold                                                                                                                                                                                                                      
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

itime_netcdf =  858#len(index_in_spock_date_netcdf_same_dt_output)-1 # max is  len(index_in_spock_date_netcdf_same_dt_output)
itime_spock = index_in_spock_date_netcdf_same_dt_output[itime_netcdf] # SpOCK otuput from spec file

dt_map = 20 #can't be greater than 20

min_lon_spock = np.min(np.min(lon_spock[itime_spock:itime_spock+dt_map]))
max_lon_spock = np.max(np.max(lon_spock[itime_spock:itime_spock+dt_map]))
min_lat_spock = np.min(np.min(lat_spock[itime_spock:itime_spock+dt_map]))
max_lat_spock = np.max(np.max(lat_spock[itime_spock:itime_spock+dt_map]))
min_lon_netcdf = np.min(lon_spec_netcdf_dt_output[itime_netcdf])
max_lon_netcdf = np.max(lon_spec_netcdf_dt_output[itime_netcdf])
min_lat_netcdf = np.min(lat_spec_netcdf_dt_output[itime_netcdf])
max_lat_netcdf = np.max(lat_spec_netcdf_dt_output[itime_netcdf])
min_lon = np.min([min_lon_spock, min_lon_netcdf])-1
max_lon = np.max([max_lon_spock, max_lon_netcdf])+1
min_lat = np.min([min_lat_spock, min_lat_netcdf])-1
max_lat = np.max([max_lat_spock, max_lat_netcdf])+1
m = Basemap( projection       = 'cyl',
             llcrnrlon        = min_lon , #Lower Left  CoRNeR Longitude
             urcrnrlon        = max_lon  , #Upper Right CoRNeR Longitude
             llcrnrlat        = min_lat  , #Lower Left  CoRNeR Latitude
             urcrnrlat        = max_lat,   #Upper Right CoRNeR Latitude
             resolution       = 'l'  ,
             suppress_ticks   = False,
             ax = ax,
             )

m.drawcoastlines(linewidth=0.7, color='blue')
lon1_spock
lon1_netcdf


print date_spock[itime_spock], r_cyg_lvlh_diff[itime_netcdf, :]
#error_spec_this_time = r_spec_diff_mag_min_to_max_array[itime_netcdf, :] # theoretically r_spec_diff_mag_min_to_max_array should not be taken at itime_netcdf 
                                                                        # but it's actually correct because we sampled the ECEF and spec files at the same time (every 60s). 
                                                                        # You can check by comparing the length of r_spec_diff_mag_min_to_max_array and 
                                                                        # the length index_in_spock_date_netcdf_same_dt_output

# for ispec in range(len(error_spec_this_time)):
#     if ispec == 0:
#         str_spec = 'Specular distance (km): ' + format(error_spec_this_time[ispec], ".0f")
#     else:
#         str_spec = str_spec + ', ' +  format(error_spec_this_time[ispec], ".0f")
# ax.text(0.02, 0.60, str_spec, transform = ax.transAxes, fontsize = fontsize_plot, horizontalalignment = 'left', verticalalignment = 'top')

error_spec_this_time = r_spec_diff_lvlh_same_prn[itime_netcdf, :, :]
str_spec = 'Specular error km):\nalong cross\n'
for ispec in range(len(error_spec_this_time)):
    for  ic in range(2):
        str_spec = str_spec + ' ' + format(error_spec_this_time[ispec, ic], ".0f")
    str_spec = str_spec + '\n'
ax.text(0.02, 0.77, str_spec, transform = ax.transAxes, fontsize = fontsize_plot, horizontalalignment = 'left', verticalalignment = 'top')


str_sc = 'Satellite error:\n' + \
         ' - along-track: ' + format(r_cyg_lvlh_diff[itime_netcdf, 0], ".1f") + ' km\n' +\
         ' - cross-track: ' + format(r_cyg_lvlh_diff[itime_netcdf, 1], ".1f") + ' km\n' +\
         ' - radial: ' + format(r_cyg_lvlh_diff[itime_netcdf, 2], ".1f") + ' km' 
ax.text(0.02, 0.94, str_sc, transform = ax.transAxes, fontsize = fontsize_plot, horizontalalignment = 'left', verticalalignment = 'top')


str_date = 'Prediction time: ' + format(( date_spock[itime_spock] - date_spock[0] ).total_seconds() / 3600. / 24, ".1f") + ' day'
ax.text(0.02, 0.98, str_date, transform = ax.transAxes, fontsize = fontsize_plot, horizontalalignment = 'left', verticalalignment = 'top')

#print 'SPOCK'
for i in range(dt_map):
    for j in range(len(lon_spock[itime_spock+i])):
        x_spock_h, y_spock_h = m(lon_spock[itime_spock+i][j], lat_spock[itime_spock+i][j]) 
        m.plot(x_spock_h, y_spock_h,  marker='.', markersize=10,color = 'b')[0]
        #print '\n',lon_spock[itime_spock+i][j], lat_spock[itime_spock+i][j]

#print '\nNETCDF'
for i in range(dt_map):    
    for j in range(4):
        x_netcdf_h, y_netcdf_h = m(lon_spec_netcdf_dt_output[itime_netcdf][i][j], lat_spec_netcdf_dt_output[itime_netcdf][i][j]) 
        #print '\n', lon_spec_netcdf_dt_output[itime_netcdf+i][j], lat_spec_netcdf_dt_output[itime_netcdf+i][j]
        m.plot(x_netcdf_h, y_netcdf_h,  marker='.', markersize=10,color = 'r')[0]

fig_save_name = '2d_spock_netcdf_day' + format(( date_spock[itime_spock] - date_spock[0] ).total_seconds() / 3600. / 24, ".1f").replace(".", "_")+ '.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')


#itime_netcdf = 100; itime_spock = index_in_spock_date_netcdf_same_dt_output[itime_netcdf]; print date_spock[itime_spock]
