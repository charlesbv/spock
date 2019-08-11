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

# July 25 2019: this script is originally a copy of cygnss_piston_campaign_2018.py
# June 21 2018: this script is originally a copy of cygnss_airfield_campaign_may2017.py.
# This script:
# - predicts the specular point locations from date_start until date_end using SpOCK
# - from the position files created by SpOCK, create new specular point location files with only the times when specular points are between a certain range of dates (set by date_range) and a certain range of latitudes [min_lat_range, max_lat_range] and longitude [min_lon_range, max_lon_range]
# To run it:
# python cygnss_piston_campaign_2018.py
# BEFORE running it, there are a few parameters to set up: see section "PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT" below
# Assumptions:
# - see assumptions from script cygnss_animation.py
# - all dates should be UTC 
# - date_range is a list. Each element corresponds to a range of time (so each element of date_range is a list of two elements: start and end date of range). Each date must be a time HH:MM:SS. Example: date_range = [["10:00:00", "14:00:00"], ["22:00:00", "02:00:00"]]

# PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT
date_start = "2019-08-25T00:00:00" #"2019-07-25T00:00:00" #"2019-08-03T00:00:00"#"2019-08-18T00:00:00"#"2019-08-25T00:00:00" # !!!!!!! UTC
date_end = "2019-09-08T23:59:59" #"2019-07-26T00:00:00" #"2019-08-17T23:59:59"#"2019-08-24T23:59:59"#"2019-09-08T23:59:59" # !!!!!!! UTC

date_range = [["00:00:00", "23:59:59"]] # !!!!!!! UTC
min_lat_range = 36.75
max_lat_range = 38.0
min_lon_range = -106.75
max_lon_range = -105.3


ncell_lat = 15 # nb of cells in the lat direction
ncell_lon = 15

# ALGORITHM
from datetime import datetime, timedelta
import sys
import os
import ipdb
sys.path.append("/Users/cbv/work/spock/srcPython")
from read_input_file import *
from cygnss_read_spock_spec_bin import *
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from matplotlib import pyplot as plt


min_lat_range = np.float(min_lat_range)
max_lat_range = np.float(max_lat_range)
min_lon_range = np.float(min_lon_range)
max_lon_range = np.float(max_lon_range)

if max_lon_range < 0:
    max_lon_range = max_lon_range + 360
if min_lon_range < 0:
    min_lon_range = min_lon_range + 360


dcell_lat = (max_lat_range - min_lat_range) / ncell_lat
dcell_lon = (max_lon_range - min_lon_range) / ncell_lon
# Predicts the specular point locations from date_start until date_end using SpOCK
## Run SpOCK
start_date_str = date_start
end_date_str =  date_end

start_date = datetime.strptime(start_date_str, "%Y-%m-%dT%H:%M:%S")
end_date = datetime.strptime(end_date_str, "%Y-%m-%dT%H:%M:%S")

#os.system("spock_cygnss_spec_parallel.py " + start_date_str + " " + end_date_str + " spec")

#ipdb.set_trace()
## Read specular point locations
### Read SpOCK main input file to figure out stuff to then read the output
spock_dir = "."
input_filename = spock_dir + '/spock_spec_start_' + start_date_str.replace(":", "_") + '_end_' + end_date_str.replace(":", "_") + '_sgp4.txt'
var_in, var_in_order = read_input_file(input_filename)
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
dt = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')]; 
nb_steps = var_in[find_in_read_input_order_variables(var_in_order, 'nb_steps')]; 
nb_sc = var_in[find_in_read_input_order_variables(var_in_order, 'nb_sc')]; 
gps_name = var_in[find_in_read_input_order_variables(var_in_order, 'gps_name')];
### Read SpOCK output files
nb_sc = 8 # !!!!!!!!!
label_arr = ['FM05', 'FM04', 'FM02', 'FM01', 'FM08', 'FM06', 'FM07', 'FM03']
lon_spec = []; lat_spec = []; gain_spec = []; gps_spec = []; date_spec = []
filename_spec_spock = []
nb_time_this_sc = []
print "Reading the SP positions..."
for isc in range(nb_sc):
    print isc, nb_sc-1
    which_sc = isc
    cyg = format(isc + 1, "02")
    filename_spec_spock.append( spock_dir + "/" + output_file_path_list[which_sc] + "specular_" + output_file_name_list[which_sc].replace(".txt",".bin") )

    data_spec = cygnss_read_spock_spec_bin(filename_spec_spock[-1], gps_name, dt, 0) 
    date_spec_this_sc = data_spec[0]; lon_spec_this_sc = data_spec[1]; lat_spec_this_sc = data_spec[2]; gain_spec_this_sc = data_spec[3]; gps_spec_this_sc = data_spec[4]

    lon_spec.append(lon_spec_this_sc)
    lat_spec.append(lat_spec_this_sc)
    gain_spec.append(gain_spec_this_sc)
    gps_spec.append(gps_spec_this_sc)
    date_spec.append(date_spec_this_sc) # should be the same for all sc
    nb_time_this_sc.append(len(date_spec_this_sc))
print "Done reading the SP positions"
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
filename_out = start_date_str.replace(":","").replace("-","") + "_to_" + end_date_str.replace(":","").replace("-","") + "_sgp4.txt" 
file_out = open(filename_out, "w")
label_arr = ['FM05', 'FM04', 'FM02', 'FM01', 'FM08', 'FM06', 'FM07', 'FM03']
lon_save = []
nb_visit_per_cell = np.zeros([ncell_lon, ncell_lat])
for itime in range(nb_time):
    if np.mod(itime, nb_time/10) == 0:
        print itime, nb_time - 1
    at_leat_one_spec_this_time = 0
    line_out = datetime.strftime(date_spec[0][itime], "%Y-%m-%dT%H:%M:%S") + " "
    for isc in range(nb_sc):
        date_current = date_spec[isc][itime]
        # for irange in range(nb_date_range): 
        #     if ( ( date_current >= range_start_this_day[irange] ) & ( date_current <= range_end_this_day[irange] ) ): # current date is in range of dates
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
                line_out = line_out + "{0:.2f}".format(lon_to_print) + " " + "{0:.2f}".format(lat_spec[isc][itime][ispec])  + " " + format(gain_spec[isc][itime][ispec], ".0f") + " "# + ' PRN' + str(gps_spec[isc][itime][ispec])
                already_put_sc_name = 1
                at_leat_one_spec_this_time = 1

                # determine which cell the sp is in
                icell_lat =  (int)((lat_spec[isc][itime][ispec] - min_lat_range) / dcell_lat)
                icell_lon =  (int)((lon_spec[isc][itime][ispec] - min_lon_range) / dcell_lon)
                nb_visit_per_cell[icell_lon, icell_lat] = nb_visit_per_cell[icell_lon, icell_lat] + 1
    if at_leat_one_spec_this_time == 1:
        print >> file_out, line_out

file_out.close()


# PLOT
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3
color_arr = ['b', 'r','cornflowerblue','g', 'm', 'gold', 'cyan', 'fuchsia', 'lawngreen', 'darkgray', 'green', 'chocolate']

fig_title = ''
y_label = 'Latitude '  + u'(\N{DEGREE SIGN})'
x_label = 'Longitude '  + u'(\N{DEGREE SIGN})'


# rectangular coordinates: x-axis azimuth, y-axis elevation
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.99,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)

ax = fig.add_subplot(gs[0, 0])

ax.set_title(date_start + ' to ' + date_end, weight = 'bold', fontsize  = fontsize_plot)
ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold


origin = 'lower'

lon_arr = np.arange(min_lon_range, max_lon_range, dcell_lon)
lat_arr = np.arange(min_lat_range, max_lat_range, dcell_lat)
x = lon_arr
y = lat_arr
X, Y = np.meshgrid(x, y)
Z = np.zeros([ncell_lat, ncell_lon])
for ilat in range(ncell_lat):
    Z[ilat, :] =  nb_visit_per_cell[:,ncell_lat-1-ilat]

nr, nc = Z.shape

Z = np.ma.array(Z)
CS1 = ax.imshow(Z, extent=[min_lon_range, max_lon_range, min_lat_range, max_lat_range], cmap = 'jet', aspect = 'auto')#,
#                vmin = 0, vmax = max_gain, origin='upper')
cbar = plt.colorbar(CS1, ax = ax)
cbar.ax.set_ylabel('# visits', fontsize = fontsize_plot, weight = 'bold')
fig_save_name = filename_out.replace('txt', 'pdf')
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
