# This script:
# - predicts the specular point locations from date_start until date_end using SpOCK
# - from the position files created by SpOCK, create new specular point location files with only the times when specular points are between a certain range of dates (set by date_range) and a certain range of latitudes [min_lat_range, max_lat_range] and longitude [min_lon_range, max_lon_range]
# To run it:
# python cygnss_airfield_campaign_may2017.py
# BEFORE running it, there are a few parameters to set up: see section "PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT" below
# Assumptions:
# - see assumptions from script cygnss_animation.py
# - all dates should be UTC 
# - date_range is a list. Each element corresponds to a range of time (so each element of date_range is a list of two elements: start and end date of range). Each date must be a time HH:MM:SS. Example: date_range = [["10:00:00", "14:00:00"], ["22:00:00", "02:00:00"]]

# PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT
date_start = "2017-04-24T00:00:00" # !!!!!!! UTC
date_end = "2017-04-24T23:59:59" # !!!!!!! UTC

date_range = [["10:00:00", "14:00:00"], ["22:00:00", "02:00:00"]] # !!!!!!! UTC
min_lat_range = 30.
max_lat_range = 40.
max_lon_range = -72.
min_lon_range = -82.


# ALGORITHM
from datetime import datetime, timedelta
import sys
import os
sys.path.append("/usr/local/bin/spock/bin")
from read_input_file import *
from cygnss_read_spock_spec import *

if max_lon_range < 0:
    max_lon_range = max_lon_range + 360
if min_lon_range < 0:
    min_lon_range = min_lon_range + 360

# Predicts the specular point locations from date_start until date_end using SpOCK
## Run SpOCK
start_date_str = date_start
end_date_str =  date_end
if start_date_str == "now":
    start_date_str = datetime.strftime(  datetime.utcnow(), "%Y-%m-%dT%H:%M:00")
    start_date_temp =  datetime.strptime(start_date_str, "%Y-%m-%dT%H:%M:%S")
    end_date_str = datetime.strftime( start_date_temp + timedelta(hours = np.float(sys.argv[2])), "%Y-%m-%dT%H:%M:00")
start_date = datetime.strptime(start_date_str, "%Y-%m-%dT%H:%M:%S")
end_date = datetime.strptime(end_date_str, "%Y-%m-%dT%H:%M:%S")

os.system("spock_cygnss_spec_parallel.py " + start_date_str + " " + end_date_str + " spec")

## Read specular point locations
### Read SpOCK main input file to figure out stuff to then read the output
spock_dir = "."
input_filename = spock_dir + '/spock_spec_start_' + start_date_str.replace(":", "_") + '_end_' + end_date_str.replace(":", "_") + '.txt'
var_in, var_in_order = read_input_file(input_filename)
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
dt = var_in[find_in_read_input_order_variables(var_in_order, 'dt')]; 
nb_steps = var_in[find_in_read_input_order_variables(var_in_order, 'nb_steps')]; 
nb_sc = var_in[find_in_read_input_order_variables(var_in_order, 'nb_sc')]; 
### Read SpOCK output files
nb_sc = 8 # !!!!!!!!!
label_arr = ['FM05', 'FM04', 'FM02', 'FM01', 'FM08', 'FM06', 'FM07', 'FM03']
lon_spec = []; lat_spec = []; gain_spec = []; gps_spec = []; date_spec = []
filename_spec_spock = []
nb_time_this_sc = []
for isc in range(nb_sc):
    which_sc = isc
    cyg = format(isc + 1, "02")
    filename_spec_spock.append( spock_dir + "/" + output_file_path_list[which_sc] + "specular_" + output_file_name_list[which_sc] )
    date_spec_this_sc, lon_spec_this_sc, lat_spec_this_sc, gain_spec_this_sc, gps_spec_this_sc = cygnss_read_spock_spec(filename_spec_spock[-1])
    lon_spec.append(lon_spec_this_sc)
    lat_spec.append(lat_spec_this_sc)
    gain_spec.append(gain_spec_this_sc)
    gps_spec.append(gps_spec_this_sc)
    date_spec.append(date_spec_this_sc) # should be the same for all sc
    nb_time_this_sc.append(len(date_spec_this_sc))

# From the position files created by SpOCK, create new specular point location files with only the times when specular points are between a certain range of dates (set by date_range) and a certain range of latitudes [min_lat_range, max_lat_range]
nb_date_range = len(date_range)
## Make an array of ranges of dates
nb_day =  (int) (np.ceil( ( end_date - start_date ).total_seconds() / 3600. / 24))
range_start_this_day = []
range_start_this_day_str = []
range_end_this_day = []
range_end_this_day_str = []
for iday in range(nb_day):
    for irange in range(nb_date_range):
        range_start_this_day_temp_str = datetime.strftime( start_date + timedelta(days = iday), "%Y-%m-%dT%H:%M:%S" )[0:11] + date_range[irange][0]
        range_start_this_day_temp = datetime.strptime(range_start_this_day_temp_str, "%Y-%m-%dT%H:%M:%S")
        range_start_this_day.append(range_start_this_day_temp)
        range_start_this_day_str.append(range_start_this_day_temp_str)

        range_end_this_day_temp_str = datetime.strftime( start_date + timedelta(days = iday), "%Y-%m-%dT%H:%M:%S" )[0:11] + date_range[irange][1]
        range_end_this_day_temp = datetime.strptime(range_end_this_day_temp_str, "%Y-%m-%dT%H:%M:%S")
        if range_end_this_day_temp < range_start_this_day_temp: # if the end date for a date range is the day after (like in ["22:00:00", "02:00:00"])
            range_end_this_day.append( range_end_this_day_temp + timedelta(days = 1))
            range_end_this_day_str.append(datetime.strftime(range_end_this_day[-1], "%Y-%m-%dT%H:%M:%S"))
        else:
            range_end_this_day.append(range_end_this_day_temp)
            range_end_this_day_str.append(range_end_this_day_temp_str)
            
        # print range_start_this_day[-1], range_end_this_day[-1]
## Make new files with specular point locations only in range of latitude and in range of dates
### Each spec file should have the same number of time but we just make sure here
if min(nb_time_this_sc) != max(nb_time_this_sc):
    print "***! Weird, all spec files don't have the same number of lines. The program will stop.***!"; raise Exception
    
nb_time = nb_time_this_sc[0]
# one file for all CYGNSS and all spec. Each line is one time. It shows "FM01 [lon_spec1, lat_spec1] [lon_spec2, lat_spec2] [lon_spec3, lat_spec3] [lon_spec4, lat_spec4]" for each FM (there might not be 4 spec each time)
filename_out = start_date_str.replace(":","").replace("-","") + "_to_" + end_date_str.replace(":","").replace("-","") + ".txt" 
file_out = open(filename_out, "w")
label_arr = ['FM05', 'FM04', 'FM02', 'FM01', 'FM08', 'FM06', 'FM07', 'FM03']
lon_save = []
for itime in range(nb_time):
    at_leat_one_spec_this_time = 0
    line_out = date_spec[0][itime] + " "
    for isc in range(nb_sc):
        date_current = datetime.strptime(date_spec[isc][itime], "%Y-%m-%dT%H:%M:%S")
        for irange in range(nb_date_range): 
            if ( ( date_current >= range_start_this_day[irange] ) & ( date_current <= range_end_this_day[irange] ) ): # current date is in range of dates
                #print date_current
                nb_spec = len(lat_spec[isc][itime])
                already_put_sc_name = 0
                for ispec in range(nb_spec):
                    if ( ( lat_spec[isc][itime][ispec] >= min_lat_range ) & ( lat_spec[isc][itime][ispec] <= max_lat_range ) & ( lon_spec[isc][itime][ispec] >= min_lon_range ) & ( lon_spec[isc][itime][ispec] <= max_lon_range ) ): # spec location in range of latitudes and longitudes
                        # Two lines below to uncomment if you want to show which CYGNSS the spec correspond to
                        # if already_put_sc_name != 1: 
                        #     line_out = line_out + label_arr[isc] + " "
                        lon_to_print = lon_spec[isc][itime][ispec]
                        if lon_to_print > 180:
                            lon_to_print = lon_to_print - 360
                        #line_out = line_out + "[" + "{0:.2f}".format(lon_to_print) + ", " + "{0:.2f}".format(lat_spec[isc][itime][ispec])  + "] "
                        line_out = line_out + "{0:.2f}".format(lon_to_print) + " " + "{0:.2f}".format(lat_spec[isc][itime][ispec])  + " "
                        already_put_sc_name = 1
                        at_leat_one_spec_this_time = 1
#                         if ( ( lon_spec[isc][itime][ispec] >= min_lon_range ) & ( lon_spec[isc][itime][ispec] <= max_lon_range ) ): 
#                             lon_save.append(lon_spec[isc][itime][ispec] )
#                             print lon_spec[isc][itime][ispec] , "I"
#                         else:
#                             print isc, date_spec[isc][itime], lon_spec[isc][itime][ispec] , "O"
                            
    if at_leat_one_spec_this_time == 1:
        print >> file_out, line_out

file_out.close()














# import pytz
# # Converts local date into UTC
# local_time = pytz.timezone ("America/Detroit")
# ## date_local_start to date_start (which is UTC)
# naive_date_local_start = datetime.strptime(date_local_start, "%Y-%m-%dT%H:%M:%S") # local time
# local_date_local_start = local_time.localize(naive_date_local_start, is_dst=None)
# date_start_str = datetime.strftime( local_date_local_start.astimezone(pytz.utc), "%Y-%m-%d %H:%M:%S" )
# date_start = datetime.strptime(date_start_str, "%Y-%m-%d %H:%M:%S") # utc
# date_start_str = datetime.strftime( date_start, "%Y-%m-%dT%H:%M:%S" ) # utc
# ## date_local_end to date_end (which is UTC)
# naive_date_local_end = datetime.strptime(date_local_end, "%Y-%m-%dT%H:%M:%S") # local time
# local_date_local_end = local_time.localize(naive_date_local_end, is_dst=None)
# date_end_str = datetime.strftime( local_date_local_end.astimezone(pytz.utc), "%Y-%m-%d %H:%M:%S" )
# date_end = datetime.strptime(date_end_str, "%Y-%m-%d %H:%M:%S") # utc
# date_end_str = datetime.strftime( date_end, "%Y-%m-%dT%H:%M:%S" ) # utc
# ## convert date_local_start_range to UTC date
# nb_date_range = len(date_local_start_range)
# date_start_range = []; date_start_range_str = []
# for idate in range(nb_date_range):
#     # take the local time corresponding to the date when tthis script is run
#     date_today = datetime.strftime(datetime.now(), "%Y-%m-%dT%H:%M:%S")[0:10]
#     date_local_start_range_with_day = date_today + "T" + date_local_start_range[idate]
#     naive_date_local_start_range_with_day = datetime.strptime(date_local_start_range_with_day, "%Y-%m-%dT%H:%M:%S") # local time
#     local_date_local_start_range_with_day = local_time.localize(naive_date_local_start_range_with_day, is_dst=None)
#     date_start_range_str_i = datetime.strftime( local_date_local_start_range_with_day.astimezone(pytz.utc), "%Y-%m-%d %H:%M:%S" )
#     date_start_range = datetime.strptime(date_start_range_str_i, "%Y-%m-%d %H:%M:%S")  # utc 
#     date_start_range_str.append( datetime.strftime( date_start_range, "%Y-%m-%dT%H:%M:%S" )[11:] ) # utc
    
