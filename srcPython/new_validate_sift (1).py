# This script compares the psotion of the specular points prediced by SpOCK VS reported in the l1 netcdf files at cygnss-sftp-1:/data/cygnss/products/l1/2017
# Assumptions:
# - see section PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

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
#from ecef2eci import *
#from eci_to_lvlh import *
from ecef_to_lvlh import *
import pickle
from cygnss_name_to_norad_id import *

# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
cygfm = 1 # which CYGNSS to look at
download_netcdf = 0 # set this variable to 1 if the entcdf files have not been download yet for the interval of time specified by [date_start_val, date_stop_val]
date_start_val_start = '2017-06-02T00:34:00'#"2017-06-01T00:00:00"
dir_run_spock = '/Users/cbv/cygnss/spock' # '.' # no slash
dir_run_spock_out = 'spock_out' #'spock' # name of output folder, no slash
pitch_max = 2. # filter out when pitch is greater than this value (in magnitude)
roll_max = 2. # filter out when roll is greater than this value (in magnitude)
yaw_max = 2. # filter out when yaw is greater than this value (in magnitude)
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
date_start_val_start = datetime.strptime(date_start_val_start, "%Y-%m-%dT%H:%M:%S")
date_start_val_array =  np.array([date_start_val_start + timedelta(days=i) for i in np.arange(1,30*5, 7)])
nb_date = len(date_start_val_array)

for idate in range(11):# !!!!!!! should be: range(nb_date):
    print idate, nb_date-1
    date_start_val = datetime.strftime(date_start_val_array[idate], "%Y-%m-%dT%H:%M:%S")#'2017-06-02T00:15:00' # !!! should be datetime.strftime(date_start_val_array[idate], "%Y-%m-%dT%H:%M:%S")#"2017-08-20T00:00:00" # date to start the validation
    date_stop_val = datetime.strftime(date_start_val_array[idate]+ timedelta(seconds = 6*24*3600+23*3600+59*60+59), "%Y-%m-%dT%H:%M:%S") #'2017-06-06T01:02:00'#'2017-06-03T00:45:00'#'2017-06-06T01:02:00'#'2017-06-02T00:56:00'#'2017-06-06T01:02:00'#'2017-06-02T00:43:00'  # !!!! shoul be + timedelta(seconds = 6*24*3600+23*3600+59*60+59), "%Y-%m-%dT%H:%M:%S")#"2017-08-20T06:00:00" # date to stop the validation

    # date_start_val = "2017-08-20T00:00:00" # date to start the validation
    # date_stop_val = "2017-08-26T23:59:59" # date to stop the validation
    # end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

    # # Run SpOCK from the date_start_val to date_stop_val
    # os.chdir(dir_run_spock )
    # #if idate != 5: #!!!!!!!!!!!!!!!!!!!! REMOVE THIS IF condition (leave the os.system right below though)
    # print "spock_cygnss_spec_parallel.py " + date_start_val + " " + date_stop_val + " spec"
    # os.system("spock_cygnss_spec_parallel.py " + date_start_val + " " + date_stop_val + " spec")
    # os.chdir("/Users/cbv/Google Drive/Work/PhD/Research/Code/cygnss/validate_sift")
    # print "Done running SpOCK"
    # # Read specular positio computed by SpOCK
    spock_input_filename = dir_run_spock  + "/spock_spec_start_" + date_start_val.replace(":", "_") + "_end_" + date_stop_val.replace(":", "_") + ".txt" # spock_spec_start_2017-08-20T00_00_00_end_2017-08-20T06_00_00.txt
    var_in, var_in_order = read_input_file(spock_input_filename)
    dt_spock_output = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')]; 
    output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
    output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')];
    cygfm_to_spock_nb = [4,3,8,2,1,6,7,5] # ['41884', '41885', '41886', '41887', '41888', '41889', '41890', '41891']
    isc =  cygfm_to_spock_nb[cygfm-1] - 1
    spec_spock_filename = dir_run_spock + "/" + output_file_path_list[isc] + "specular_" + output_file_name_list[isc]
    print "Reading SpOCK specular files..."
    date_spock, lon_spock, lat_spock, gain_spock, gps_spock, normpower_spock, x_cyg_spock, y_cyg_spock, z_cyg_spock, x_gps_spock, y_gps_spock, z_gps_spock,  x_spec_spock, y_spec_spock, z_spec_spock, nb_spec_spock, order_spec_spock = cygnss_read_spock_spec(spec_spock_filename)

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
    date_spock_not_interpolated_temp = var_out[find_in_read_input_order_variables(var_out_order, 'date')] # date_spock is every second because it corresponds to the sampling time of the spec points. 
    date_spock_not_interpolated = []
    for ii in range(len(date_spock_not_interpolated_temp)):
        date_here = date_spock_not_interpolated_temp[ii]
        date_here_date = datetime.strptime(date_here[:19], "%Y/%m/%d %H:%M:%S") # !!!!!! don't look at microseconds since SpOCK outputs at best eevey second anyway (so microseeconds should be 0)
        date_spock_not_interpolated.append(date_here_date)
    #date_spock_not_interpolated is the date directly output by SpOCK so it's with touput time step chosen in the main input file

    # Download netcdf files from sftp-1
    date_start_val_doy = (int)((datetime.strptime( date_start_val, "%Y-%m-%dT%H:%M:%S" )).strftime('%j'))# - timedelta(days = 2)).strftime('%j')) # start two day earlier so you can compare 
    # to the position given by the TLE that was used to initialize SpOCK. TLEs are usually output every day but go 2 days back to make sure
    date_stop_val_doy = (int)(datetime.strptime( date_stop_val, "%Y-%m-%dT%H:%M:%S" ).strftime('%j'))
    if datetime.strptime( date_stop_val, "%Y-%m-%dT%H:%M:%S" ).strftime('%Y') != '2017':
        print "***!The analysis has to be for data in 2017. The program will stop. !***"; raise Exception
    doy_array = np.arange(date_start_val_doy, date_stop_val_doy+1, 1)
    nb_day = len(doy_array)
    day_remove = [] # list of days to remove from the analysis
    day_list = []
    for iday in range(nb_day):
        if (os.path.isdir("/Users/cbv/cygnss/netcdf/" + str(doy_array[iday]).zfill(3)) == False):
            os.system("mkdir /Users/cbv/cygnss/netcdf/" + str(doy_array[iday]).zfill(3))
        if download_netcdf == 1:
            os.system("scp -p cygnss-sftp-1.engin.umich.edu:/data/cygnss/products/l1/2017/" + str(doy_array[iday]).zfill(3) + "/cyg0" + str(cygfm) + "* /Users/cbv/cygnss/netcdf/" + str(doy_array[iday]).zfill(3))
        if len([x for x in os.listdir("/Users/cbv/cygnss/netcdf/" + str(doy_array[iday]).zfill(3)) if x.endswith(".nc") and x.startswith("cyg0" + str(cygfm))]) != 1: # if more than one file for this sc then don't take this day into account in the analsyis OR if no netcd
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
        filename_spec_flight = "/Users/cbv/cygnss/netcdf/" + np.str(doy_array[iday_here]).zfill(3) + "/" + [x for x in os.listdir("/Users/cbv/cygnss/netcdf/" + np.str(doy_array[iday_here]).zfill(3)) if x.endswith(".nc") and x.startswith("cyg0" + np.str(cygfm))][0]
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
        while itime < nb_time_flight_temp - 60:
            itime = itime + 1
            print itime, nb_time_flight_temp-1, iday_count, nb_day-1
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
            # # !!!!!!!!! remove block below
            # if type(fom_netcdf_temp) == ma.core.MaskedArray:
            #     if (np.min(fom_netcdf_temp.data[itime]) == 0):
            #         time_remove = 1
            #         time_gain_0.append(itime)
            # else:
            #     if (np.min(fom_netcdf_temp[itime]) == 0):
            #         time_remove = 1
            #         time_gain_0.append(itime)
            # # !!!!!!!!! end of remove block below
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
                break


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
        x_spec_spock = np.array(x_spec_spock); y_spec_spock = np.array(y_spec_spock); z_spec_spock = np.array(z_spec_spock)
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
        time_diff_prn = []
        correct_vs_fom = np.zeros([16, 2]) # for each value of RCG (also noted FOM), compute percentage of correct PRN (correct = PRn necdf is same as PRN SpOCK). correct_vs_fom[ifom, 0] is the number of correct for this ifom, correct_vs_fom[ifor, 1] is the total number of spec with this fom. So percentrage correct = correct_vs_fom[ifom, 0] / correct_vs_fom[ifor, 1 * 100.
        for itime in range(nb_time_spock_netcdf):
            #itime = 5 #### !!!!!! remove
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
                    #r_spec_spock_eci, bla = ecef2eci( r_spec_spock, np.zeros([3]), datetime.strftime(date_flight_rounded[itime], "%Y-%m-%dT%H:%M:%S"), 0 )
                    #r_spec_eci_diff = r_spec_netcdf_eci - r_spec_spock_eci
                    #r_spec_lvlh_diff = eci_to_lvlh(r_cyg_eci_spock[itime, :], v_cyg_eci_spock[itime, :], r_spec_eci_diff) # ok to take r_cyg_eci_spock at itime 
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
    #         # Take the min distance of the 4*4 distance array new_dist_all_spec_to_all_spec. Save the min distance, remove the correspodning spock and netcdf spec (they are the ones that min the difference) 
    #         # from new_dist_all_spec_to_all_spec. Now new_dist_all_spec_to_all_spec is 3*3. Save the new min, and remove the corresponding spock and netcdf. 
    #         # Now new_dist_all_spec_to_all_spec is 2*2. Proceeed again twice. So r_spec_diff_mag_min_to_max_temp shows 4 differences, from min to max
    #         new_dist_all_spec_to_all_spec = np.copy(dist_all_spec_to_all_spec)
    #         new_dist_all_spec_to_all_spec_lvlh = np.copy(dist_all_spec_to_all_spec_lvlh)
    #         ispec_min = 0
    #         r_spec_diff_mag_min_to_max_temp = []        
    #         comb_min_all = []
    #         prn_comb_min = []
    #         spock_comb_min = [] # each time, record spec 0, 1, 2, 3
    #         netcdf_comb_min = [] # each time, record spec 0, 1, 2, 3
    #         while ispec_min < nb_spec_spock:
    #             comb_min = np.unravel_index(new_dist_all_spec_to_all_spec.argmin(), new_dist_all_spec_to_all_spec.shape)
    #             r_spec_diff_mag_min_to_max_temp.append( new_dist_all_spec_to_all_spec[comb_min] )
    #             r_spec_diff_mag_min_to_max_array[itime, ispec_min] = new_dist_all_spec_to_all_spec[comb_min]
    #             r_spec_diff_min_to_max_array_lvlh[itime, ispec_min, :] = new_dist_all_spec_to_all_spec_lvlh[comb_min]
    #             #new_dist_all_spec_to_all_spec = np.delete(np.delete(new_dist_all_spec_to_all_spec, comb_min[0],0), comb_min[1],1)
    #             spock_spec_here = comb_min[0]
    #             netcdf_spec_here = comb_min[1]
    #             new_dist_all_spec_to_all_spec[spock_spec_here, :] = 1e10
    #             new_dist_all_spec_to_all_spec[:, netcdf_spec_here] = 1e10
    #             new_dist_all_spec_to_all_spec_lvlh[spock_spec_here, :] = 1e10
    #             new_dist_all_spec_to_all_spec_lvlh[:, netcdf_spec_here] = 1e10
    #             prn_comb_min_netcdf = gps_netcdf[itime][netcdf_spec_here]
    #             prn_comb_min_spock = (int)(gps_spock_same_time_as_netcdf[itime][spock_spec_here].replace('PRN_',''))
    #             ispec_min = ispec_min + 1
    #             comb_min_all.append(comb_min)
    #             prn_comb_min.append([prn_comb_min_spock, prn_comb_min_netcdf])
    #             spock_comb_min.append(spock_spec_here)
    #             netcdf_comb_min.append(netcdf_spec_here)
    # #            print new_dist_all_spec_to_all_spec
    #         r_spec_diff_mag_min_to_max.append(r_spec_diff_mag_min_to_max_temp)
    #         comb_min_vs_time.append(comb_min_all)
    #         prn_comb_min_vs_time.append(prn_comb_min)
    #         spock_comb_min_vs_time.append(spock_comb_min)
    #         netcdf_comb_min_vs_time.append(netcdf_comb_min)
    #         spock_spec_lowest_to_highest_normpower.append(np.argsort(normpower_spock_same_time_as_netcdf[itime]))
    #         netcdf_spec_lowest_to_highest_normpower.append(np.argsort(rcg_netcdf[itime]))

    #         which_normpower_rank_is_biggest_dist_spock_spec[itime] = np.where(spock_spec_lowest_to_highest_normpower[itime] == spock_comb_min_vs_time[itime][-1])[0][0] # spock_comb_min_vs_time[itime][-1] is the spock spec that gave the
    #         # max distance with the netcdf spec. if which_normpower_rank_is_biggest_dist_spock_spec is 0, it means that the spock spec with the lowest (because which_normpower_rank_is_biggest_dist_spock_spec is 0) normpower
    #         # is this spec that gives the max distance with the netcdf spec 
    #         which_normpower_rank_is_biggest_dist_netcdf_spec[itime] = np.where(netcdf_spec_lowest_to_highest_normpower[itime] == netcdf_comb_min_vs_time[itime][-1])[0][0]


    #         # Look at which spock spec and which netcdf spec minimize the difference
    #         which_spec_spock_min = comb_min[0] # # which spock spec that minimizes the distance to netdf
    #         which_spec_netcdf_min = comb_min[1] # which netcdf spec that minimizes the distance to spock
    #         which_comb_min.append( [which_spec_spock_min, which_spec_netcdf_min] ) # combination [spock_spec, netcdf_spec] that minimizes the difference in spec position
    #         gain_netcdf_that_min_error[itime] = gain_netcdf[itime][which_spec_netcdf_min] # netcdf gain of the spec that minimizes the distance to spock
    #         gain_spock_that_min_error[itime] = gain_spock_same_time_as_netcdf[itime][which_spec_spock_min] # spock gain of the spec that minimizes the distance to netcdf



            # different thing here. look at the dsitance between 2 spec with the same prn (this is what we should have always done)
            ## oder netcdf spec lowest to highest FOM
            index_sort_fom = np.argsort(fom_netcdf[itime])
            order_now = 0
            ispec_order_now = 0
            order_spec_netcdf[itime, index_sort_fom[ispec_order_now]] = 0
            already_diff_prn = 0
            for ispec_order in range(1,4):
                if fom_netcdf[itime][index_sort_fom[ispec_order]] == fom_netcdf[itime][index_sort_fom[ispec_order-1]] :
                    order_spec_netcdf[itime, index_sort_fom[ispec_order]] = order_now
                else:
                    order_now = order_now + 1
                    order_spec_netcdf[itime, index_sort_fom[ispec_order]] = order_now

            index_sort_normpower_spock = np.argsort(normpower_spock_same_time_as_netcdf[itime])
            ispec_count = 0 
            for ispec in range(4):
                ifom = fom_netcdf[itime][ispec]
                correct_vs_fom[ifom, 1] = correct_vs_fom[ifom, 1] + 1
                prn_netcdf_here = gps_netcdf[itime][ispec]
                if len(np.where(gps_spock_same_time_as_netcdf[itime] == prn_netcdf_here)[0]) > 0:
                    correct_vs_fom[ifom, 0] = correct_vs_fom[ifom, 0] + 1
                    where_prn_spock = np.where(gps_spock_same_time_as_netcdf[itime] == prn_netcdf_here)[0][0]
                    r_spec_spock = np.array([x_spec_spock_same_time_as_netcdf[itime][where_prn_spock], y_spec_spock_same_time_as_netcdf[itime][where_prn_spock], z_spec_spock_same_time_as_netcdf[itime][where_prn_spock]])
                    r_spec_netcdf = np.array([x_spec_netcdf[itime][ispec], y_spec_netcdf[itime][ispec], z_spec_netcdf[itime][ispec]])
                    r_spec_diff = r_spec_spock - r_spec_netcdf
                    # Calculate disteance in LVLH frame 
                    # In theory, the LVLH frame hsould be taken in the frame of reference of each SpOCK spec point, but this is basically
                    # the same as the SpOCK satellite reference frame since the direction of motion of the satellite is very similar
                    # to the direction of motion of the specular point (the GPS don't move in a interval of a few seconds)
                    #r_spec_netcdf_eci, bla = ecef2eci( r_spec_netcdf, np.zeros([3]), datetime.strftime(date_flight_rounded[itime], "%Y-%m-%dT%H:%M:%S"), 0 ) # we don't have the velcooity of the spec but we don't care abwout it. 
                    #It's not needded to convert the velocity from ecef to eci (it's called by the funciton ecef2eci but used only to convert the ecef velcoity into eci velocity, so it doesn't  have any effect on the postiion)
                    #r_spec_spock_eci, bla = ecef2eci( r_spec_spock, np.zeros([3]), datetime.strftime(date_flight_rounded[itime], "%Y-%m-%dT%H:%M:%S"), 0 )
                    #r_spec_eci_diff = r_spec_netcdf_eci - r_spec_spock_eci
                    #r_spec_lvlh_diff = eci_to_lvlh(r_cyg_eci_spock[itime, :], v_cyg_eci_spock[itime, :], r_spec_eci_diff) # ok to take r_cyg_eci_spock at itime 
                    r_spec_ecef_diff = r_spec_netcdf - r_spec_spock
                    r_spec_lvlh_diff = ecef_to_lvlh(r_cyg_spock_ecef_file_same_time_as_netcdf[itime, :], v_cyg_spock_ecef_file_same_time_as_netcdf[itime, :], r_spec_ecef_diff) # ok to take r_cyg_eci_spock at itime 


                    # because we sampled the ECEF and spec files at the same time (every 60s). 
                    # You can check by comparing the length of r_spec_diff_mag_min_to_max_array and 
                    # the length index_in_spock_date_netcdf_same_dt_output
                    r_spec_diff_lvlh_same_prn[itime, ispec, :] = r_spec_lvlh_diff
                else:
                    which_spec_diff[itime, ispec] = order_spec_netcdf[itime, ispec]
                    if already_diff_prn == 0:
                        time_diff_prn.append(itime)
                    already_diff_prn = 1



        time_diff_prn_arr = np.array(time_diff_prn)
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

        pickle.dump(correct_vs_fom, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_correct_vs_fom.pickle", "w"))
        pickle.dump(correct_vs_fom_percent, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_correct_vs_fom_percent.pickle", "w"))
        pickle.dump(time_diff_prn, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_time_diff_prn.pickle", "w"))
        pickle.dump(which_spec_diff, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_which_spec_diff.pickle", "w"))
        pickle.dump(gps_spock, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_gps_spock.pickle", "w"))
        pickle.dump(fom_spock, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_fom_spock.pickle", "w"))
        pickle.dump(gps_netcdf, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_gps_netcdf.pickle", "w"))
        pickle.dump(fom_netcdf, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_fom_netcdf.pickle", "w"))

        pickle.dump(fom_netcdf_diff_prn_all, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_fom_netcdf_diff_prn_all.pickle", "w"))
        pickle.dump(nb_spec_diff_gain_not_0, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_nb_spec_diff_gain_not_0.pickle", "w"))
        pickle.dump(r_spec_diff_lvlh_same_prn, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_r_spec_diff_lvlh_same_prn.pickle", "w"))
        pickle.dump(nb_seconds_since_initial_epoch_spock, open("/Users/cbv/cygnss/pickle/" + spock_input_filename.split('/')[-1].replace(".txt", "")  + "_nb_seconds_since_initial_epoch_spock.pickle", "w"))

raise Exception
        
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



 





#     pickle.dump(nb_seconds_since_initial_epoch_spock, open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_nb_seconds_since_initial_epoch_spock_again.pickle", "w"))
#     pickle.dump(r_spec_diff_mag_min_to_max_array, open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_spec_diff_mag_min_to_max_array_again.pickle", "w"))
#     pickle.dump(r_cyg_diff_mag, open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_cyg_diff_mag_again.pickle", "w"))
#     pickle.dump(r_cyg_spock, open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_cyg_spock_again.pickle", "w"))
#     pickle.dump(r_cyg_netcdf, open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_cyg_netcdf_again.pickle", "w"))
#     #print date_start_val


# median and quantile over a day for a particular date
nb_day_simu = (int)(np.ceil((nb_seconds_since_initial_epoch_spock[-1] - nb_seconds_since_initial_epoch_spock[0] )/ 3600. / 24.))
median_r_spec_diff_lvlh = np.zeros([nb_day_simu, 4, 3])
q10_r_spec_diff_lvlh = np.zeros([nb_day_simu, 4, 3])
q25_r_spec_diff_lvlh = np.zeros([nb_day_simu, 4, 3])
q75_r_spec_diff_lvlh = np.zeros([nb_day_simu, 4, 3])
q90_r_spec_diff_lvlh = np.zeros([nb_day_simu, 4, 3])
x_axis_day = np.arange(0.5, nb_day_simu+0.5)+nb_seconds_since_initial_epoch_spock[0] / 3600. / 24. # start at 0.5 and go up to nb_day_simu + 0.5 because say that the median and quantiles are over a day so take it at the middle of the day
mag_r_spec_diff_min_to_max_array_lvlh = r_spec_diff_lvlh_same_prn# r_spec_diff_min_to_max_array_lvlh # !!!!!!!! np.abs(r_spec_diff_min_to_max_array_lvlh)

iday = 0
itime_previous = 0
itime_day_list = [] # list of itime for that falls every day
itime_day_list.append(itime_previous)
while iday < nb_day_simu:
    itime = 0
    while (int)( ( nb_seconds_since_initial_epoch_spock[itime+itime_previous] - nb_seconds_since_initial_epoch_spock[itime_previous] )/ 3600. / 24.) < 1:
        itime = itime + 1
        if itime+itime_previous == mag_r_spec_diff_min_to_max_array_lvlh.shape[0]:
            break
    median_r_spec_diff_lvlh[iday, :, :] = np.median(mag_r_spec_diff_min_to_max_array_lvlh[itime_previous:itime+itime_previous, :, :], axis = 0)
    q10_r_spec_diff_lvlh[iday, :, :] = np.percentile(mag_r_spec_diff_min_to_max_array_lvlh[itime_previous:itime+itime_previous, :, :], 10, axis = 0)
    q25_r_spec_diff_lvlh[iday, :, :] = np.percentile(mag_r_spec_diff_min_to_max_array_lvlh[itime_previous:itime+itime_previous, :, :], 25, axis = 0)
    q75_r_spec_diff_lvlh[iday, :, :] = np.percentile(mag_r_spec_diff_min_to_max_array_lvlh[itime_previous:itime+itime_previous, :, :], 75, axis = 0)
    q90_r_spec_diff_lvlh[iday, :, :] = np.percentile(mag_r_spec_diff_min_to_max_array_lvlh[itime_previous:itime+itime_previous, :, :], 90, axis = 0)
    #print itime_previous, itime+itime_previous, median_r_spec_diff_lvlh[iday, 0, 0] 
    iday = iday + 1
    itime_previous = itime + itime_previous
    itime_day_list.append(itime_previous)


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
raise Exception
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
nb_date = 1# !!!!!!!!! remove this line
for idate in range(nb_date):
    #print idate, nb_date-1
    #date_start_val = datetime.strftime(date_start_val_array[idate], "%Y-%m-%dT%H:%M:%S")#"2017-08-20T00:00:00" # date to start the validation
    date_start_val = date_start_val #'2017-06-02T00:15:00' # !!!!!! shoule be line right above and not this line
    nb_seconds_since_initial_epoch_spock_all.append( pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_nb_seconds_since_initial_epoch_spock_again.pickle")))
    r_spec_diff_mag_min_to_max_array_all.append( pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_spec_diff_mag_min_to_max_array_again.pickle")))
    r_cyg_diff_mag_all.append( pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_cyg_diff_mag_again.pickle")))
    r_cyg_spock_all.append( pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_cyg_spock_again.pickle")))
    r_cyg_netcdf_all.append( pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_cyg_netcdf_again.pickle")))
    nb_time_spock_netcdf_all.append( len(r_cyg_netcdf_all[-1]))


# nb_seconds_since_initial_epoch_spock = pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_nb_seconds_since_initial_epoch_spock.pickle"))
# r_spec_diff_mag_min_to_max_array = pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_spec_diff_mag_min_to_max_array.pickle"))
# r_cyg_diff_mag = pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_cyg_diff_mag.pickle"))
# r_cyg_spock = pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_cyg_spock.pickle"))
# r_cyg_netcdf = pickle.load(open("/Users/cbv/cygnss/pickle/" + date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_")	 + "_r_cyg_netcdf.pickle"))
# nb_time_spock_netcdf = len(r_cyg_netcdf)




# if befre you run the script cygnss_spec_netcdf_to_spock_format.py, the following pickle was saved. it includes the dates to ignore to look at the spec position. Indeed those dates were missing dates in the netcdf 
idate = 0
date_ignore_missing_netcdf_data_all_sc = pickle.load(open("date_ignore_missing_netcdf_data_all_sc_" + date_start_ok_str.replace(":","_") + '_to_' + date_stop_ok_str.replace(":","_")	 + ".pickle", "r"))
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
    else:
        count_gt_threshold = count_gt_threshold + 1

    
ax.margins(0,0)
ax.set_ylim([0,10])
ax.set_xlim([0,3])

# legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
# legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_") + '_all_dates_difference_spec_position_spock_vs_netcdf_new.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')




fig_title = ''#Difference in specular point position SpOCK vs L1' # L1 is netcdf                                                                            
x_label = 'Prediction time (days)'
for ispec in range(4):

    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                                                                                                                                                                                                                                      
    gs = gridspec.GridSpec(1, 2)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01, wspace = 0.5)
    count_gt_threshold = 0
    #for idate in range(nb_date):
    ax = fig.add_subplot(gs[0, 0])
    y_label = 'Along-track difference (km)'
    ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                                                                                                                                                                                                                                                                         
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           

    ax.scatter(nb_seconds_since_initial_epoch_spock/ 3600. / 24., mag_r_spec_diff_min_to_max_array_lvlh[:,ispec,0], marker = '+', s = 5, color = 'k', label = 'Along-track', alpha = 0.4)
    ax.plot(x_axis_day, median_r_spec_diff_lvlh[:,ispec,0], linewidth = 3, color = 'k', label = '')
    ax.plot(x_axis_day, q10_r_spec_diff_lvlh[:,ispec,0], linewidth = 3, color = 'r', label = '', linestyle = 'dashed')
    ax.plot(x_axis_day, q25_r_spec_diff_lvlh[:,ispec,0], linewidth = 3, color = 'b', label = '', linestyle = 'dashed')
    ax.plot(x_axis_day, q75_r_spec_diff_lvlh[:,ispec,0], linewidth = 3, color = 'b', label = '', linestyle = 'dashed')
    ax.plot(x_axis_day, q90_r_spec_diff_lvlh[:,ispec,0], linewidth = 3, color = 'r', label = '', linestyle = 'dashed')

    ax.margins(0,0)
    #ax.set_ylim([np.min(q25_r_spec_diff_lvlh[:,ispec,0]), np.max(q75_r_spec_diff_lvlh[:,ispec,0])])
    ax.set_ylim([np.min(q10_r_spec_diff_lvlh[:,ispec,0]), np.max(q90_r_spec_diff_lvlh[:,ispec,0])])

    ax.text(0.02, 0.98, 'spec ' + str(ispec), fontsize = fontsize_plot, transform = ax.transAxes, verticalalignment = 'top' )
    ax = fig.add_subplot(gs[0, 1])
    y_label = 'Cross-track difference (km)'
    ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                                                                                                                                                                                                                                                                         
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           

    ax.scatter(nb_seconds_since_initial_epoch_spock/ 3600. / 24., mag_r_spec_diff_min_to_max_array_lvlh[:,ispec,1], marker = '+', s = 5, color = 'k', label = 'Along-track', alpha = 0.4)
    ax.plot(x_axis_day, median_r_spec_diff_lvlh[:,ispec,1], linewidth = 3, color = 'k', label = '')
    ax.plot(x_axis_day, q10_r_spec_diff_lvlh[:,ispec,1], linewidth = 3, color = 'r', label = '', linestyle = 'dashed')
    ax.plot(x_axis_day, q25_r_spec_diff_lvlh[:,ispec,1], linewidth = 3, color = 'b', label = '', linestyle = 'dashed')
    ax.plot(x_axis_day, q75_r_spec_diff_lvlh[:,ispec,1], linewidth = 3, color = 'b', label = '', linestyle = 'dashed')
    ax.plot(x_axis_day, q90_r_spec_diff_lvlh[:,ispec,1], linewidth = 3, color = 'r', label = '', linestyle = 'dashed')


    ax.margins(0,0)
    #ax.set_ylim([np.min(q25_r_spec_diff_lvlh[:,ispec,1]), np.max(q75_r_spec_diff_lvlh[:,ispec,1])])
    ax.set_ylim([np.min(q10_r_spec_diff_lvlh[:,ispec,1]), np.max(q90_r_spec_diff_lvlh[:,ispec,1])])
    # legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot)
    # legend.get_title().set_fontsize(str(fontsize_plot))

    fig_save_name = date_start_val.replace(":","_") + '_to_' + date_stop_val.replace(":","_") +  'specular_lvlh_error_ispec_' + str(ispec) + '.pdf'
    fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')



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










height_fig = 15.  # the width is calculated as height_fig * 4/3.                                                                                                                                                                                                             
fontsize_plot = 25      
ratio_fig_size = 4./3
fig_title = ''#VLLH distribution specular points for a particular ispec at a particular time
y_label = 'Percentage (%)'
itime_day = 2 # which day to look at the distirbution over (distibution oof all spec over this day). from 1 (not 0) to nb_day (not nb_day-1). ex: if itime_day = 2 then it'll show the distribution from day 1 to day 2
coordtype = ['Along-track', 'Cross-track']
for icoord in range(1,2):
    #icoord = 0 # 0 for along-track, 1 for cross-track
    x_label = coordtype[icoord] + ' difference (km)'
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                                                                                                                  
    gs = gridspec.GridSpec(2, 2)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.17)
    for ispec in range(4):
        diff_here = mag_r_spec_diff_min_to_max_array_lvlh[itime_day_list[itime_day - 1]:itime_day_list[itime_day],ispec, icoord]
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
        print np.min(diff_here[where_same_prn]), np.max(diff_here[where_same_prn]),len(where_same_prn) * 100. / len(diff_here)
        range_max_temp = q75_r_spec_diff_lvlh[itime_day - 1, ispec,icoord]#q75_r_spec_diff_lvlh[itime_day - 1, ispec,icoord]*1.5 #np.max(mag_r_spec_diff_min_to_max_array_lvlh[where_same_nb_spec_spock_netcdf, ispec, icoord])
        range_min_temp = median_r_spec_diff_lvlh[itime_day - 1, ispec,icoord] - (q75_r_spec_diff_lvlh[itime_day - 1, ispec,icoord] - median_r_spec_diff_lvlh[itime_day - 1, ispec,icoord] ) #q25_r_spec_diff_lvlh[itime_day - 1, ispec,icoord]#q25_r_spec_diff_lvlh[itime_day - 1, ispec,icoord]*0.5 
        half_width = (range_max_temp - range_min_temp)/2 #(q75_r_spec_diff_lvlh[itime_day - 1, ispec,icoord] - q25_r_spec_diff_lvlh[itime_day - 1, ispec,icoord])/2.
        range_min  = np.min(diff_here[where_same_prn])#range_min_temp - half_width*3
        range_max  = np.max(diff_here[where_same_prn])#range_max_temp + half_width*3
        hist_along_data = np.histogram(diff_here, range = [range_min, range_max])
        bin_array_temp = hist_along_data[1]
        bin_array = ( bin_array_temp[:-1] + np.roll(bin_array_temp,-1)[:-1] ) /2.
        binsize_actual = bin_array[1] - bin_array[0]
        hist_along = hist_along_data[0] * 100. / len(diff_here)
        ax.bar(bin_array, hist_along, binsize_actual)
        ax.plot([median_r_spec_diff_lvlh[itime_day - 1, ispec,icoord], median_r_spec_diff_lvlh[itime_day - 1, ispec,icoord]], [0, 100], linestyle = 'dashed', linewidth = 2, color = 'k')
        ax.plot([q25_r_spec_diff_lvlh[itime_day - 1, ispec,icoord], q25_r_spec_diff_lvlh[itime_day - 1, ispec,icoord]], [0, 100], linestyle = 'dashed', linewidth = 2, color = 'b')
        ax.plot([q75_r_spec_diff_lvlh[itime_day - 1, ispec,icoord], q75_r_spec_diff_lvlh[itime_day - 1, ispec,icoord]], [0, 100], linestyle = 'dashed', linewidth = 2, color = 'b')
        ax.plot([q10_r_spec_diff_lvlh[itime_day - 1, ispec,icoord], q10_r_spec_diff_lvlh[itime_day - 1, ispec,icoord]], [0, 100], linestyle = 'dashed', linewidth = 2, color = 'r')
        ax.plot([q90_r_spec_diff_lvlh[itime_day - 1, ispec,icoord], q90_r_spec_diff_lvlh[itime_day - 1, ispec,icoord]], [0, 100], linestyle = 'dashed', linewidth = 2, color = 'r')
        ax.set_ylim([np.min(hist_along), np.max(hist_along)])
        ax.set_xlim([np.min(bin_array_temp), np.max(bin_array_temp)])
        ax.text(0.02, 0.98, '50% width:\n' +r'$\pm$' + format(half_width, ".1f") + ' km', fontsize = fontsize_plot, transform = ax.transAxes, verticalalignment = 'top' )
        if ispec == 0:
            ax.text(0.98, 0.98, 'spec ' + str(ispec + 1) + '\nday' + str(itime_day), fontsize = fontsize_plot, transform = ax.transAxes, verticalalignment = 'top', horizontalalignment = 'right' )
        else:
            ax.text(0.98, 0.98, 'spec ' + str(ispec + 1), fontsize = fontsize_plot, transform = ax.transAxes, verticalalignment = 'top', horizontalalignment = 'right' )
    fig.set_figheight(height_fig)
    fig.set_figwidth(height_fig*ratio_fig_size)
    fig_save_name = 'distribution_specular_lvlh_error_ispec_' + coordtype[icoord].replace("-","").lower() + '_day_' + str(itime_day) +'.pdf'
    fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')


disthere= 80
print len(np.where(np.abs(r_spec_diff_min_to_max_array_lvlh[itime_day_list[itime_day - 1]:itime_day_list[itime_day],ispec,icoord]) < disthere )[0]) * 100. / len(r_spec_diff_min_to_max_array_lvlh[itime_day_list[itime_day - 1]:itime_day_list[itime_day],ispec,icoord])

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

error_spec_this_time = mag_r_spec_diff_min_to_max_array_lvlh[itime_netcdf, :, :]
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
