# This script compares the psotion of the specular points prediced by SpOCK VS reported in the l1 netcdf files at cygnss-sftp-1:/data/cygnss/products/l1/2017
# Assumptions:
# - see section PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

from datetime import datetime, timedelta
import numpy as np
import os
import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/kalman/spock_development_new_structure_kalman_dev/srcPython")
from os import listdir
from read_input_file import *
from cygnss_read_spock_spec import *
from netCDF4 import Dataset
import numpy.ma as ma



# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
date_start_val = "2017-08-20T00:00:00" # date to start the validation
date_stop_val = "2017-08-20T06:00:00" # date to stop the validation
cygfm = 1 # which CYGNSS to look at
downlowd_netcdf = 0 # set this variable to 1 if the entcdf files have not been download yet for the interval of time specified by [date_start_val, date_stop_val]
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

# Run SpOCK from the date_start_val to date_stop_val
os.chdir("spock")
#os.system("spock_cygnss_spec_parallel.py " + date_start_val + " " + date_stop_val + " spec")
os.chdir("../")

# Read specular positio computed by SpOCK
spock_input_filename = "spock/spock_spec_start_" + date_start_val.replace(":", "_") + "_end_" + date_stop_val.replace(":", "_") + ".txt" # spock_spec_start_2017-08-20T00_00_00_end_2017-08-20T06_00_00.txt
var_in, var_in_order = read_input_file(spock_input_filename)
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')];
cygfm_to_spock_nb = [4,3,8,2,1,6,7,5] # ['41884', '41885', '41886', '41887', '41888', '41889', '41890', '41891']
isc = cygfm_to_spock_nb[cygfm-1] - 1
spec_spock_filename = "spock/" + output_file_path_list[isc] + "specular_" + output_file_name_list[isc]
date_spock, lon_spock, lat_spock, gain_spock, gps_spock, normpower_spock, x_cyg_spock, y_cyg_spock, z_cyg_spock, x_gps_spock, y_gps_spock, z_gps_spock,  x_spec_spock, y_spec_spock, z_spec_spock = cygnss_read_spock_spec(spec_spock_filename)



# Download netcdf files from sftp-1
date_start_val_doy = (int)(datetime.strptime( date_start_val, "%Y-%m-%dT%H:%M:%S" ).strftime('%j'))
date_stop_val_doy = (int)(datetime.strptime( date_stop_val, "%Y-%m-%dT%H:%M:%S" ).strftime('%j'))
if datetime.strptime( date_stop_val, "%Y-%m-%dT%H:%M:%S" ).strftime('%Y') != '2017':
    print "***!The analysis has to be for data in 2017. The program will stop. !***"; raise Exception
doy_array = np.arange(date_start_val_doy, date_stop_val_doy+1, 1)
nb_day = len(doy_array)
day_remove = [] # list of days to remove from the analysis
day_list = []
for iday in range(nb_day):
    if (os.path.isdir("netcdf/" + str(doy_array[iday]).zfill(3)) == False):
        os.system("mkdir netcdf/" + str(doy_array[iday]).zfill(3))
    if downlowd_netcdf == 1:
        os.system("scp -p cygnss-sftp-1.engin.umich.edu:/data/cygnss/products/l1/2017/" + str(doy_array[iday]).zfill(3) + "/cyg0" + str(cygfm) + "* ./netcdf/" + str(doy_array[iday]).zfill(3))
    if len([x for x in os.listdir("./netcdf/" + str(doy_array[iday]).zfill(3)) if x.endswith(".nc") and x.startswith("cyg0" + str(cygfm))]) > 1: # if more than one file for this sc then don't take this day into account in the analsyis
        day_remove.append(iday)
    else:
        day_list.append(iday)

nb_day = len(day_list)

# For each day, read the specular point position from the netcdf file
iday_count = -1

y_spec_netcdf = []
z_spec_netcdf = []
x_cyg_netcdf = []
y_cyg_netcdf = []
z_cyg_netcdf = []
date_flight = []
date_flight_rounded = []
index_in_spock_date_netcdf_same = [] 
for iday in day_list:
    iday_count = iday_count + 1
    filename_spec_flight = "./netcdf/" + str(doy_array[iday]).zfill(3) + "/" + [x for x in os.listdir("./netcdf/" + str(doy_array[iday]).zfill(3)) if x.endswith(".nc") and x.startswith("cyg0" + str(cygfm))][0]
    fh = Dataset(filename_spec_flight, mode='r')
    # nc_attrs = fh.ncattrs()
    # nc_dims = [dim for dim in fh.dimensions]  # list of nc dimensions
    # nc_vars = [var for var in fh.variables]  # list of nc variables

    # X component of the specular point position in the ECEF coordinate system, in meters, at ddm_timestamp_utc, as calculated on the ground.
    x_spec_netcdf_temp = fh.variables['sp_pos_x'][:]
    y_spec_netcdf_temp = fh.variables['sp_pos_y'][:]
    z_spec_netcdf_temp = fh.variables['sp_pos_z'][:]
    x_cyg_netcdf_temp = fh.variables['sc_pos_x'][:]
    y_cyg_netcdf_temp = fh.variables['sc_pos_y'][:]
    z_cyg_netcdf_temp= fh.variables['sc_pos_z'][:]
    # sc_vel_x = fh.variables['sc_vel_x'][:]
    # sc_vel_y = fh.variables['sc_vel_y'][:]
    # sc_vel_z = fh.variables['sc_vel_z'][:]

    time_flight = fh.variables['ddm_timestamp_utc'][:]
    time_coverage_start = fh.getncattr(fh.ncattrs()[fh.ncattrs().index('time_coverage_start')])
    time_coverage_start_datetime = datetime.strptime(time_coverage_start[:-4], "%Y-%m-%dT%H:%M:%S.%f") 
    #fh.close()
    nb_time_flight_temp = len(x_cyg_netcdf_temp)
    date_flight_t = []
    date_flight_rounded_temp = []
    time_remove_list = []
    for itime in range(nb_time_flight_temp):
        time_remove = 0
        date_flight_temp = datetime.strftime(time_coverage_start_datetime + timedelta(microseconds = round(time_flight[itime]*10**6)), "%Y-%m-%dT%H:%M:%S.%f" )
        # round to neared second but only if the actual time is less than 100 ms from the nearest second, otherwise ignore this time (don't compare to SpOCK).This is because SpOCK propagates with 1 s time step. So to compare to the netcdf file, we assume that the netcdf is also exactly at each second (so no millisecond). 100 ms wrong is not too bad because the satellite movesby less than 1 km.
        if ( ( True in x_spec_netcdf_temp.mask[itime] ) || ( True in y_spec_netcdf_temp.mask[itime] ) || ( True in z_spec_netcdf_temp.mask[itime] ) || ( True in x_cyg_netcdf_temp.mask[itime] ) || ( True in x_cyg_netcdf_temp.mask[itime] ) || ( True in x_cyg_netcdf_temp.mask[itime] ) ): # id one of the 4 spec is masked then ignore this time
            time_remove = 1
        if ( date_flight_temp.split('.')[1][0] == '9' ): # round to next second
            date_flight_temp_date = datetime.strptime(date_flight_temp, "%Y-%m-%dT%H:%M:%S.%f")
            date_flight_date = date_flight_temp_date + timedelta(seconds = 1)
            date_flight_date_rounded_temp = datetime.strftime(date_flight_date, "%Y-%m-%dT%H:%M:%S.%f").split('.')[0]
            date_flight_date_rounded = datetime.strptime(date_flight_date_rounded_temp, "%Y-%m-%dT%H:%M:%S")
            # if date_flight_date_rounded in date_spock:
            #     index_in_spock_date_netcdf_same.append(date_spock.index(date_flight_date_rounded))
                # date_flight_rounded_temp.append( date_flight_date_rounded )
                # date_flight_t.append( date_flight_temp )
            # else:# if this time is not in date_spock
            #     time_remove = 1
        elif ( date_flight_temp.split('.')[1][0] == '0' ): # round to next second
            date_flight_date_rounded = datetime.strptime(date_flight_temp.split('.')[0], "%Y-%m-%dT%H:%M:%S" )
            # if date_flight_date_rounded in date_spock:
            #     index_in_spock_date_netcdf_same.append(date_spock.index(date_flight_date_rounded))
                # date_flight_rounded_temp.append( date_flight_date_rounded )
                # date_flight_t.append( date_flight_temp )
            # else: # if this time is not in date_spock
            #     time_remove = 1
        else: #if time can't be rounded by less than 100 ms
            # del x_spec_netcdf_temp[itime]
            # del y_spec_netcdf_temp[itime]
            # del z_spec_netcdf_temp[itime]
            # del x_cyg_netcdf_temp[itime]
            # del y_cyg_netcdf_temp[itime]
            # del z_cyg_netcdf_temp[itime]
            # time_remove_list.append(itime)
            # time_already_removed = 1
            time_remove = 1
        if date_flight_date_rounded in date_spock: # if this time is not in date_spock
            index_in_spock_date_netcdf_same.append(date_spock.index(date_flight_date_rounded))
        else:
            time_remove = 1
        if ( time_remove == 1 ): # remove time if can't be rounded by ess than 100 ms or if is not in date_spock or if masked
            # del x_spec_netcdf_temp[itime]
            # del y_spec_netcdf_temp[itime]
            # del z_spec_netcdf_temp[itime]
            # del x_cyg_netcdf_temp[itime]
            # del y_cyg_netcdf_temp[itime]
            # del z_cyg_netcdf_temp[itime]
            time_remove_list.append(itime)
        else:
            # if iday == 0:
            #     x_spec_netcdf = x_spec_netcdf_temp
            #     y_spec_netcdf = y_spec_netcdf_temp
            #     z_spec_netcdf = z_spec_netcdf_temp
            #     x_cyg_netcdf = x_cyg_netcdf_temp
            #     y_cyg_netcdf = y_cyg_netcdf_temp
            #     z_cyg_netcdf = z_cyg_netcdf_temp
    # else:
    #     x_spec_netcdf = ma.concatenate([x_spec_netcdf, x_spec_netcdf_temp])
    #     y_spec_netcdf = ma.concatenate([y_spec_netcdf, y_spec_netcdf_temp])
    #     z_spec_netcdf = ma.concatenate([z_spec_netcdf, z_spec_netcdf_temp])
    #     x_cyg_netcdf = ma.concatenate([x_cyg_netcdf, x_cyg_netcdf_temp])
    #     y_cyg_netcdf = ma.concatenate([y_cyg_netcdf, y_cyg_netcdf_temp])
    #     z_cyg_netcdf = ma.concatenate([z_cyg_netcdf, z_cyg_netcdf_temp])
    date_flight = date_flight + date_flight_t
    date_flight_rounded = date_flight_rounded + date_flight_rounded_temp
    raise Exception
nb_time_flight_rounded = len(date_flight_rounded)
date_netcdf = date_flight_rounded # change of notation
#True x_spec_netcdf[itime].mask



# Compare SpOCK and netcdf
## Find the dates that are both in SpOCK and netcdf









