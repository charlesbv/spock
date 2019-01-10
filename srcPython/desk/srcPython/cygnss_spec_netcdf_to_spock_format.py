# This script reads in a netcdf file (cygnss-sftp-1:/data/cygnss/products/l1/2017) and writes a file with position of CYGNSS so it can be input of the find_specular_points.c code: CYGNSS_CONSTELLATION_  file (usually created by SpOCK). Can't do the same with the GPS bevause the netcdf files only have the 4 GPS recorded (the ones selecte from the 4 highest gains). 
# It also creates attitude files for SpOCK with the same attitude as the one reported in the netcdf files.
# Goal: compare the specular point position output by find_specular_points.c to the one reported in the netcdf files
# Assumptions:
# - see section PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
# - SpOCK is run with a 60 s time step output from spock_cygnss_spec_parallel.py
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
from netCDF4 import Dataset
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from matplotlib import pyplot as plt
from ecef2eci import *
from eci_to_lvlh import *
import pickle
from norad_id_to_cygnss_name import *


# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

downlowd_netcdf = 0 # set this variable to 1 if the entcdf files have not been download yet for the interval of time specified by [date_start, date_stop]
date_start = '2017-07-16T00:00:00'#'2017-06-03T00:34:00'#"2017-06-02T00:00:00"
date_stop = '2017-07-16T23:59:59'#'2017-06-03T00:50:00' #2017-06-06T01:02:00'#'2017-06-03T00:45:00'#"2017-06-02T00:56:00"  #"2017-06-02T00:43:00"
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

first_date = []
last_date = []

date_flight_rounded_dt_output_all_sc = []
x_cyg_netcdf_dt_output_all_sc = []
pitch_cyg_netcdf_dt_output_all_sc = []
roll_cyg_netcdf_dt_output_all_sc = []
yaw_cyg_netcdf_dt_output_all_sc = []
y_cyg_netcdf_dt_output_all_sc = []
z_cyg_netcdf_dt_output_all_sc = []
vx_cyg_netcdf_dt_output_all_sc = []
vy_cyg_netcdf_dt_output_all_sc = []
vz_cyg_netcdf_dt_output_all_sc = []
x_spec_netcdf_dt_output_all_sc = [] 
y_spec_netcdf_dt_output_all_sc = []
z_spec_netcdf_dt_output_all_sc = []
lon_spec_netcdf_dt_output_all_sc = [] 
lat_spec_netcdf_dt_output_all_sc = [] 
date_ignore_missing_netcdf_data_all_sc = []
nb_sc = 8
for isc in range(nb_sc):
    cygfm = (int)(norad_id_to_cygnss_name(str(41884+isc)).split('FM0')[1])
    #cygfm = 1
    # Download netcdf files from sftp-1
    date_start_date = datetime.strptime( date_start, "%Y-%m-%dT%H:%M:%S" )
    date_stop_date = datetime.strptime( date_stop, "%Y-%m-%dT%H:%M:%S" )
    date_start_doy = (int)((datetime.strptime( date_start, "%Y-%m-%dT%H:%M:%S" )).strftime('%j'))
    date_stop_doy = (int)(datetime.strptime( date_stop, "%Y-%m-%dT%H:%M:%S" ).strftime('%j'))
#     if datetime.strptime( date_stop, "%Y-%m-%dT%H:%M:%S" ).strftime('%Y') != '2017':
#         print "***!The analysis has to be for data in 2017. The program will stop. !***"; raise Exception        
    doy_array = np.arange(date_start_doy, date_stop_doy+1, 1) #!!!!! shoule be np.arange(date_start_doy, date_stop_doy+1, 1)
    nb_day = len(doy_array)
    day_remove = [] # list of days to remove from the analysis
    day_list = []
    if  datetime.strptime( date_stop, "%Y-%m-%dT%H:%M:%S" ).strftime('%Y') == '2018':
        yy = '2018/'
    elif  datetime.strptime( date_stop, "%Y-%m-%dT%H:%M:%S" ).strftime('%Y') == '2017':
        yy = '2017/'
    else:
        yy = ''
    for iday in range(nb_day):
        if (os.path.isdir("/Users/cbv/cygnss/netcdf/" + yy  + str(doy_array[iday]).zfill(3)) == False):
            os.system("mkdir /Users/cbv/cygnss/netcdf/" + yy  + str(doy_array[iday]).zfill(3))
        if downlowd_netcdf == 1:
            os.system("scp -p cygnss-sftp-1.engin.umich.edu:/data/cygnss/products/l1/" + datetime.strptime( date_stop, "%Y-%m-%dT%H:%M:%S" ).strftime('%Y') +"/" + str(doy_array[iday]).zfill(3) + "/cyg0" + str(cygfm) + "* /Users/cbv/cygnss/netcdf/" + yy  + str(doy_array[iday]).zfill(3))
        if len([x for x in os.listdir("/Users/cbv/cygnss/netcdf/" + yy  + str(doy_array[iday]).zfill(3)) if x.endswith(".nc") and x.startswith("cyg0" + str(cygfm))]) != 1: # if more than one file for this sc then don't take this day into account in the analsyis OR if no netcd
            day_remove.append(iday)
        else:
            day_list.append(iday)
        
    nb_day = len(day_list)

    ################
    ################
    ################
    # For each day, read the specular point position from the netcdf file
    x_spec_netcdf = []
    y_spec_netcdf = []
    z_spec_netcdf = []
    lat_spec_netcdf = [] 
    lon_spec_netcdf = [] 

    x_cyg_netcdf = []
    pitch_cyg_netcdf = []
    roll_cyg_netcdf = []
    yaw_cyg_netcdf = []
    y_cyg_netcdf = []
    z_cyg_netcdf = []
    vx_cyg_netcdf = []
    vy_cyg_netcdf = []
    vz_cyg_netcdf = []

    x_cyg_netcdf_dt_output_temp = [] # this is to compare to the r/v output in ECEF_ output SpOCK files
    pitch_cyg_netcdf_dt_output_temp = [] 
    roll_cyg_netcdf_dt_output_temp = [] 
    yaw_cyg_netcdf_dt_output_temp = [] 
    y_cyg_netcdf_dt_output_temp = []
    z_cyg_netcdf_dt_output_temp = []
    vx_cyg_netcdf_dt_output_temp = []
    vy_cyg_netcdf_dt_output_temp = []
    vz_cyg_netcdf_dt_output_temp = []
    x_spec_netcdf_dt_output_temp = [] 
    y_spec_netcdf_dt_output_temp = []
    z_spec_netcdf_dt_output_temp = []
    lon_spec_netcdf_dt_output_temp = [] 
    lat_spec_netcdf_dt_output_temp = [] 


    gain_netcdf = []
    date_flight = []
    date_flight_rounded = []
    index_in_spock_date_netcdf_same = [] 
    index_in_spock_date_netcdf_same_dt_output = [] 
    date_flight_rounded_dt_output_temp = []
    index_in_spock_not_interpolated_date_netcdf_same = []
    nb_seconds_since_initial_epoch_spock = []
    nb_seconds_since_initial_epoch_spock_dt_output_temp = []
    iday_count = -1
    while iday_count < nb_day-1:# !!!!! should be nb_day-1:
        iday_count = iday_count + 1
        print cygfm, doy_array[iday_count]
        iday_here = day_list[iday_count]
        filename_spec_flight = "/Users/cbv/cygnss/netcdf/" + yy + np.str(doy_array[iday_here]).zfill(3) + "/" + [x for x in os.listdir("/Users/cbv/cygnss/netcdf/" + yy  + np.str(doy_array[iday_here]).zfill(3)) if x.endswith(".nc") and x.startswith("cyg0" + np.str(cygfm))][0]
        fh = Dataset(filename_spec_flight, mode='r')
        # nc_attrs = fh.ncattrs()
        # nc_dims = [dim for dim in fh.dimensions]  # list of nc dimensions
        # nc_vars = [var for var in fh.variables]  # list of nc variables

        x_spec_netcdf_temp = fh.variables['sp_pos_x'][:] # X component of the specular point position in the ECEF coordinate system, in meters, at ddm_timestamp_utc, as calculated on the ground.
        y_spec_netcdf_temp = fh.variables['sp_pos_y'][:]
        z_spec_netcdf_temp = fh.variables['sp_pos_z'][:]

        lat_spec_netcdf_temp = fh.variables['sp_lat'][:]
        lon_spec_netcdf_temp = fh.variables['sp_lon'][:]

        x_cyg_netcdf_temp = fh.variables['sc_pos_x'][:] # for the 4 GPS: tx_pos_x
        pitch_cyg_netcdf_temp = fh.variables['sc_pitch'][:] 
        roll_cyg_netcdf_temp = fh.variables['sc_roll'][:] 
        yaw_cyg_netcdf_temp = fh.variables['sc_yaw'][:] 
        y_cyg_netcdf_temp = fh.variables['sc_pos_y'][:]
        z_cyg_netcdf_temp= fh.variables['sc_pos_z'][:]
        gain_netcdf_temp = fh.variables['sp_rx_gain'][:] # The receive antenna gain in the direction of the specular point, in dBi, at ddm_timestamp_utc
        vx_cyg_netcdf_temp = fh.variables['sc_vel_x'][:]
        vy_cyg_netcdf_temp = fh.variables['sc_vel_y'][:]
        vz_cyg_netcdf_temp= fh.variables['sc_vel_z'][:]


        list_are_masked_array = [] # sometimes the netcdf varaible below are maksed array and soemtimes they are not (depending on which netcdf file)...
        if type(x_spec_netcdf_temp) == ma.core.MaskedArray:
            list_are_masked_array.append(x_spec_netcdf_temp)
        if type(y_spec_netcdf_temp) == ma.core.MaskedArray:
            list_are_masked_array.append(y_spec_netcdf_temp)
        if type(z_spec_netcdf_temp) == ma.core.MaskedArray:
            list_are_masked_array.append(z_spec_netcdf_temp)
        if type(gain_netcdf_temp) == ma.core.MaskedArray:
            list_are_masked_array.append(gain_netcdf_temp)
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
        # The CYGNSS_CONSTELLATION_  files are usually output 
        # every minute so select the cyngss every minute. If it's 
        # not available in the netcdf files, then keep the position constant. 
        # The post processing will ignore those times, though

        while itime < nb_time_flight_temp - 60:
            itime = itime + 1
            #print itime, nb_time_flight_temp-1, iday_count, nb_day-1
            time_remove = 0
            date_flight_temp_date = time_coverage_start_datetime + timedelta(microseconds = round(time_flight[itime]*10**6))
            date_flight_temp = datetime.strftime(date_flight_temp_date, "%Y-%m-%dT%H:%M:%S.%f" )
            # round to neared second but only if the actual time is less than 100 ms from the nearest second, otherwise ignore this time 
            if ( date_flight_temp.split('.')[1][0] == '9' ): # round to next second
                date_flight_temp_date = datetime.strptime(date_flight_temp, "%Y-%m-%dT%H:%M:%S.%f")
                date_flight_date = date_flight_temp_date + timedelta(seconds = 1)
                date_flight_date_rounded_temp = datetime.strftime(date_flight_date, "%Y-%m-%dT%H:%M:%S.%f").split('.')[0]
                date_flight_date_rounded = datetime.strptime(date_flight_date_rounded_temp, "%Y-%m-%dT%H:%M:%S")
            elif ( date_flight_temp.split('.')[1][0] == '0' ): # round to next second
                date_flight_date_rounded = datetime.strptime(date_flight_temp.split('.')[0], "%Y-%m-%dT%H:%M:%S" )
            else: #if time can't be rounded by less than 100 ms
                time_remove = 1




            date_flight_raw.append(date_flight_temp)
            imask_arr = 0
            while (imask_arr < nb_mask_array):
                        #if ( ( True in x_spec_netcdf_temp.mask[itime] ) | ( True in y_spec_netcdf_temp.mask[itime] ) | ( True in z_spec_netcdf_temp.mask[itime] ) | ( True in gain_netcdf_temp.mask[itime] ) ):# id one of the 4 spec is masked then ignore this time
                if (  True in list_are_masked_array[imask_arr].mask[itime] ):
                     time_remove = 1
                imask_arr = imask_arr + 1     
            if ( (time_remove == 0) & (np.mod( ( date_flight_date_rounded - date_start_date ).total_seconds(), 60. ) == 0)): # selecet every minute
                x_cyg_netcdf_dt_output_temp.append(x_cyg_netcdf_temp[itime]/1000.)
                pitch_cyg_netcdf_dt_output_temp.append(pitch_cyg_netcdf_temp[itime]*180./np.pi)
                roll_cyg_netcdf_dt_output_temp.append(roll_cyg_netcdf_temp[itime]*180./np.pi)
                yaw_cyg_netcdf_dt_output_temp.append(yaw_cyg_netcdf_temp[itime]*180./np.pi)
                y_cyg_netcdf_dt_output_temp.append(y_cyg_netcdf_temp[itime]/1000.)
                z_cyg_netcdf_dt_output_temp.append(z_cyg_netcdf_temp[itime]/1000.)
                vx_cyg_netcdf_dt_output_temp.append(vx_cyg_netcdf_temp[itime]/1000.)
                vy_cyg_netcdf_dt_output_temp.append(vy_cyg_netcdf_temp[itime]/1000.)
                vz_cyg_netcdf_dt_output_temp.append(vz_cyg_netcdf_temp[itime]/1000.)
                x_spec_netcdf_dt_output_temp.append(x_spec_netcdf_temp[itime]/1000.)
                y_spec_netcdf_dt_output_temp.append(y_spec_netcdf_temp[itime]/1000.)
                z_spec_netcdf_dt_output_temp.append(z_spec_netcdf_temp[itime]/1000.)
                lon_spec_netcdf_dt_output_temp.append(lon_spec_netcdf_temp[itime:itime+20]) # save 20 time steps here
                lat_spec_netcdf_dt_output_temp.append(lat_spec_netcdf_temp[itime:itime+20])
                date_flight_rounded_dt_output_temp.append(date_flight_date_rounded)
                if date_flight_rounded_dt_output_temp[-1] > date_stop_date:
                    break


    nb_time_temp = len(date_flight_rounded_dt_output_temp)
    date_flight_rounded_dt_output = []
    x_cyg_netcdf_dt_output = []
    pitch_cyg_netcdf_dt_output = []
    roll_cyg_netcdf_dt_output = []
    yaw_cyg_netcdf_dt_output = []
    y_cyg_netcdf_dt_output = []
    z_cyg_netcdf_dt_output = []
    vx_cyg_netcdf_dt_output = []
    vy_cyg_netcdf_dt_output = []
    vz_cyg_netcdf_dt_output = []
    x_spec_netcdf_dt_output = [] 
    y_spec_netcdf_dt_output = []
    z_spec_netcdf_dt_output = []
    lon_spec_netcdf_dt_output = [] 
    lat_spec_netcdf_dt_output = [] 
    date_ignore_missing_netcdf_data = []
    for itime in range(nb_time_temp-1):
        date_flight_rounded_dt_output.append(date_flight_rounded_dt_output_temp[itime])
        x_cyg_netcdf_dt_output.append(x_cyg_netcdf_dt_output_temp[itime])
        pitch_cyg_netcdf_dt_output.append(pitch_cyg_netcdf_dt_output_temp[itime])
        roll_cyg_netcdf_dt_output.append(roll_cyg_netcdf_dt_output_temp[itime])
        yaw_cyg_netcdf_dt_output.append(yaw_cyg_netcdf_dt_output_temp[itime])
        y_cyg_netcdf_dt_output.append(y_cyg_netcdf_dt_output_temp[itime])
        z_cyg_netcdf_dt_output.append(z_cyg_netcdf_dt_output_temp[itime])
        vx_cyg_netcdf_dt_output.append(vx_cyg_netcdf_dt_output_temp[itime])
        vy_cyg_netcdf_dt_output.append(vy_cyg_netcdf_dt_output_temp[itime])
        vz_cyg_netcdf_dt_output.append(vz_cyg_netcdf_dt_output_temp[itime])
        x_spec_netcdf_dt_output.append(x_spec_netcdf_dt_output_temp[itime]) 
        y_spec_netcdf_dt_output.append(y_spec_netcdf_dt_output_temp[itime])
        z_spec_netcdf_dt_output.append(z_spec_netcdf_dt_output_temp[itime])
        lon_spec_netcdf_dt_output.append(lon_spec_netcdf_dt_output_temp[itime]) 
        lat_spec_netcdf_dt_output.append(lat_spec_netcdf_dt_output_temp[itime]) 

        if ( date_flight_rounded_dt_output_temp[itime+1] - date_flight_rounded_dt_output_temp[itime] ).total_seconds() > 60:
            #print date_flight_rounded_dt_output_temp[itime+1] , date_flight_rounded_dt_output_temp[itime]
            # The CYGNSS_CONSTELLATION_ files are usually output 
            # every minute so select the cyngss psition every minute. If it's 
            # not available in the netcdf files, then keep the position constant. 
            # The post processing will ignore those times, though
            while ( date_flight_rounded_dt_output[-1] < date_flight_rounded_dt_output_temp[itime+1] - timedelta(seconds = 60)):
                date_flight_rounded_dt_output.append( date_flight_rounded_dt_output[-1] + timedelta(seconds = 60) )
                date_ignore_missing_netcdf_data.append(date_flight_rounded_dt_output[-1])
                x_cyg_netcdf_dt_output.append(x_cyg_netcdf_dt_output_temp[itime])
                pitch_cyg_netcdf_dt_output.append(pitch_cyg_netcdf_dt_output_temp[itime])
                roll_cyg_netcdf_dt_output.append(roll_cyg_netcdf_dt_output_temp[itime])
                yaw_cyg_netcdf_dt_output.append(yaw_cyg_netcdf_dt_output_temp[itime])
                y_cyg_netcdf_dt_output.append(y_cyg_netcdf_dt_output_temp[itime])
                z_cyg_netcdf_dt_output.append(z_cyg_netcdf_dt_output_temp[itime])
                vx_cyg_netcdf_dt_output.append(vx_cyg_netcdf_dt_output_temp[itime])
                vy_cyg_netcdf_dt_output.append(vy_cyg_netcdf_dt_output_temp[itime])
                vz_cyg_netcdf_dt_output.append(vz_cyg_netcdf_dt_output_temp[itime])
                x_spec_netcdf_dt_output.append(x_spec_netcdf_dt_output_temp[itime]) 
                y_spec_netcdf_dt_output.append(y_spec_netcdf_dt_output_temp[itime])
                z_spec_netcdf_dt_output.append(z_spec_netcdf_dt_output_temp[itime])
                lon_spec_netcdf_dt_output.append(lon_spec_netcdf_dt_output_temp[itime]) 
                lat_spec_netcdf_dt_output.append(lat_spec_netcdf_dt_output_temp[itime]) 


    nb_time = len(date_flight_rounded_dt_output)
    for itime in range(nb_time-1):
        if ( date_flight_rounded_dt_output[itime+1] - date_flight_rounded_dt_output[itime] ).total_seconds() > 60:
            print '***! ERROR ', date_flight_rounded_dt_output[itime+1] , date_flight_rounded_dt_output[itime] , '!***'


    date_flight_rounded_dt_output_all_sc.append(date_flight_rounded_dt_output)
    date_ignore_missing_netcdf_data_all_sc.append(date_ignore_missing_netcdf_data)
    x_cyg_netcdf_dt_output_all_sc.append(x_cyg_netcdf_dt_output)
    pitch_cyg_netcdf_dt_output_all_sc.append(pitch_cyg_netcdf_dt_output)
    roll_cyg_netcdf_dt_output_all_sc.append(roll_cyg_netcdf_dt_output)
    yaw_cyg_netcdf_dt_output_all_sc.append(yaw_cyg_netcdf_dt_output)
    y_cyg_netcdf_dt_output_all_sc.append(y_cyg_netcdf_dt_output)
    z_cyg_netcdf_dt_output_all_sc.append(z_cyg_netcdf_dt_output)
    vx_cyg_netcdf_dt_output_all_sc.append(vx_cyg_netcdf_dt_output)
    vy_cyg_netcdf_dt_output_all_sc.append(vy_cyg_netcdf_dt_output)
    vz_cyg_netcdf_dt_output_all_sc.append(vz_cyg_netcdf_dt_output)
    x_spec_netcdf_dt_output_all_sc.append(x_spec_netcdf_dt_output) 
    y_spec_netcdf_dt_output_all_sc.append(y_spec_netcdf_dt_output)
    z_spec_netcdf_dt_output_all_sc.append(z_spec_netcdf_dt_output)
    lon_spec_netcdf_dt_output_all_sc.append(lon_spec_netcdf_dt_output) 
    lat_spec_netcdf_dt_output_all_sc.append(lat_spec_netcdf_dt_output) 

    first_date.append(date_flight_rounded_dt_output[0])
    last_date.append(date_flight_rounded_dt_output[-1])

#date_ignore_missing_netcdf_data_all_sc



first_date = np.array(first_date)
last_date = np.array(last_date)

date_start_ok = np.max(first_date)
if date_start_ok < date_start_date:
    date_start_ok = date_start_date
date_stop_ok = np.min(last_date)
date_start_ok_str = datetime.strftime( date_start_ok, "%Y-%m-%dT%H:%M:%S" )
date_stop_ok_str = datetime.strftime( date_stop_ok, "%Y-%m-%dT%H:%M:%S" )

# dump_file_date_ignore_missing_netcdf_data_all_sc = open("date_ignore_missing_netcdf_data_all_sc_" + date_start_ok_str.replace(":","_") + '_to_' + date_stop_ok_str.replace(":","_")     + ".pickle", "w")
# pickle.dump(date_ignore_missing_netcdf_data_all_sc, dump_file_date_ignore_missing_netcdf_data_all_sc)
# dump_file_date_ignore_missing_netcdf_data_all_sc.close()
raise Exception
# Run SpOCK from date_start_ok_str to date_stop_ok_str so that you have the main input file and the 
# GPS_CONSTELLATION_ file (both are necessary inputs for find_specular_points.c)
#os.system("spock_cygnss_spec_parallel.py " + date_start_ok_str + " " + date_stop_ok_str + " spec") #!!!!! dt_output needs to be 60s
main_input_filename = 'spock_spec_start_' + date_start_ok_str.replace(":","_") + '_end_' + \
    date_stop_ok_str.replace(":","_") + '.txt'
filename_cygnss_out = 'spock_out/spock_spec_start_' + date_start_ok_str.replace(":","_") + '_end_' + \
    date_stop_ok_str.replace(":","_")  + '/CONSTELLATION_CYGNSS_for_run_spock_spec_start_' + \
    date_start_ok_str.replace(":","_") + '_end_' + date_stop_ok_str.replace(":","_") + '1.txt' 

# Make a copy of the CONSTELLATION_CYGNSS_ written by SpOCK since it's going ot be overwritten by the netcdf equivalent file
os.system("cp " +  filename_cygnss_out + " " + filename_cygnss_out.replace('.txt', '_by_spock.txt'))
# Write CYGNSS constllation ECEF position/velocity in output file: name same as one created by SpOCK 
# using spock_cygnss_spec_parallel.py date_start_ok_str date_stop_ok_str spec

file_cygnss_out = open(filename_cygnss_out, "w")

for isc in range(nb_sc):
    cygfm = (int)(norad_id_to_cygnss_name(str(41884+isc)).split('FM0')[1])

    print >> file_cygnss_out, '// ----------------------------------------------------------------------------------\n' + \
                                  '//\n' + \
                                  '// Please note this file is auto generated by the Spacecraft Orbital Characterization Kit (SpOCK)\n' +\
                                  '// Trajectory Spacecraft spock_spec_start_' + date_start_ok_str.replace(":","_") + '_end_' + date_stop_ok_str.replace(":","_")  + \
                                  str(isc+1) + ' in a constellation of 8 spacecraft and 31 GPS satellites (iProc 0)\n' +\
                                  '//\n' +  \
                                  '// Trajectory is specified with:\n' + \
                                  '//	TIME  		      r_ecef2cg_ECEF(3)(KM) 		  v_ecef2cg_ECEF(3)(KM/S)\n' +\
                                  '//\n' + \
                                  "// Version control of SpOCK is under Joel Getchius's Mac and Charles Bussy-Virat's Mac\n" + \
                                  '// ----------------------------------------------------------------------------------\n' +\
                                   '#START'

    nb_time = len(date_flight_rounded_dt_output_all_sc[isc])
    itime = 0
    while (itime < nb_time):
        if date_flight_rounded_dt_output_all_sc[isc][itime] >= date_start_ok:
            print >>  file_cygnss_out, str(date_flight_rounded_dt_output_all_sc[isc][itime]).replace("-", "/"), \
                x_cyg_netcdf_dt_output_all_sc[isc][itime], y_cyg_netcdf_dt_output_all_sc[isc][itime], z_cyg_netcdf_dt_output_all_sc[isc][itime], \
                vx_cyg_netcdf_dt_output_all_sc[isc][itime], vy_cyg_netcdf_dt_output_all_sc[isc][itime], vz_cyg_netcdf_dt_output_all_sc[isc][itime]
        itime = itime + 1
        if itime == nb_time:
            break
        if (date_flight_rounded_dt_output_all_sc[isc][itime] > date_stop_ok):
            break
    print >>  file_cygnss_out,''
    print >>  file_cygnss_out,''
file_cygnss_out.close()


# Now compute the specular point locations using the new CYGNSS_CONSTELLATION_ file (written from the r/v of CYGNSS reported in the netcdf files)
# and using the GPS_CONSTELLATION_ from the SpOCK propagation (since the r/v of GPS is reported only for 4 GPS in the netcdf files so they can'be used)
## first make a copy of the spec position calculated by SpOCK previously
for isc in range(nb_sc):
    filename_spec_from_netcdf_cyngss_input = 'spock_out/spock_spec_start_' + date_start_ok_str.replace(":","_") + '_end_' + \
    date_stop_ok_str.replace(":","_")  + '/spock_spec_start_' + date_start_ok_str.replace(":","_") + '_end_' + \
    date_stop_ok_str.replace(":","_") + str(isc+1) + '/specular_spock_spec_start_' + date_start_ok_str.replace(":","_") + '_end_' + \
    date_stop_ok_str.replace(":","_") +str(isc+1) + '.txt'
    filename_spec_spock_copy = filename_spec_from_netcdf_cyngss_input.replace('.txt', '_by_spock.txt')
    
    os.system("cp " +  filename_spec_from_netcdf_cyngss_input + " " + filename_spec_spock_copy)
## Now run the spaculr point code
os.system("mpirun -np 4 spec_dev " + main_input_filename + " -lon=0 -rot=0 -min")
os.system("mpirun -np 4 storm " + main_input_filename )





# ATTTIUDE file

for isc in range(nb_sc):
    cygfm = (int)(norad_id_to_cygnss_name(str(41884+isc)).split('FM0')[1])
    filename_cygnss_out = 'attitude/attitude_FM0' + str(cygfm) + '_spock_spec_start_' + date_start_ok_str.replace(":","_") + '_end_' + \
        date_stop_ok_str.replace(":","_")  + '.txt'

    file_cygnss_out = open(filename_cygnss_out, "w")

    print >> file_cygnss_out, '#BEGINNINGOFHEADER\n#ENDOFHEADER'

    nb_time = len(date_flight_rounded_dt_output_all_sc[isc])
    itime = 0
    while (itime < nb_time):
        if date_flight_rounded_dt_output_all_sc[isc][itime] >= date_start_ok:
            print >>  file_cygnss_out, str(date_flight_rounded_dt_output_all_sc[isc][itime]).replace(" ", "T"), \
                '(', pitch_cyg_netcdf_dt_output_all_sc[isc][itime], ';', roll_cyg_netcdf_dt_output_all_sc[isc][itime], ';' ,yaw_cyg_netcdf_dt_output_all_sc[isc][itime], ')', \
                ' (1; 2; 3)'
        itime = itime + 1
        if itime == nb_time:
            break
        if (date_flight_rounded_dt_output_all_sc[isc][itime] > date_stop_ok):
            break
    print >>  file_cygnss_out,'#ENDOFFILE'
    file_cygnss_out.close()
