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
# This script plots the CYGNSS sub-satellite location predicted by SpOCK on the 2D map with the IR brightness temperatures on the background.
# The input, given as argument of this script, is the date at which the orbits need to be plot. Example:
# python ir_brightness_with_cygnss.py start_date
# where start_date is the date where to start the animation.
# if 'now' instead of a date then the date is set to the current time (but the seconds are set to 0)
# The IR brightness temperature corresponding to this date are downloaded from NCEP at NOAA. Only the first orbit of each CYGNSS is shown
# Assumptions:
# - the format of the date must be YYYY-MM-DDTHH
# - results are put in the run directory of SpOCK simulation 
# - the IR data from NCEP is not stored
# - see section 'paramaters to set up before running this script' in this script
# Notes:
# - if the tart_date is too soon so that no data is available at NCEP yet, then this sciprt determines the earliest data available at NCEP and sets the new start date for this date
# A few notes on the IR brightness temperature data:
# -------------------------------------------------------------------------------------------
# Data at: ftp://ftp.cpc.ncep.noaa.gov/precip/global_full_res_IR
# -------------------------------------------------------------------------------------------
# Documentation at http://www.cpc.ncep.noaa.gov/products/global_precip/html/README.
# Part of the documentation:

# Each file contains 2 records: the 1st for the "on the hour" images (":00") and the 2nd for
# the "on the half hour" images (":30")

# Each record is a 9896 x 3298 Fortran array of IR brightness temperatures that have
# been scaled to fit into 1-byte by subtracting "75" from each datum.  Therefore it
# is necessary for the user to add a value of "75" to each data value when using the
# data.

# The orientation of the data in the array is EASTward from 0.0182E (center of
# gridbox) and the grid increment in the east-west direction is 0.036378335 degrees
# of longitude.  The data proceed from North -> South beginning at 59.982N (center of
# gridbox) and the grid increment is 0.036383683 degrees of latitude in the
# north-sout direction.
# -------------------------------------------------------------------------------------------
# Examples of images at: http://www.cloudsat.cira.colostate.edu/quicklooks
# -------------------------------------------------------------------------------------------
# End of a few notes on the IR brightness temperature data:
import matplotlib
matplotlib.use("Agg") # without this line, when running this script from a cronjob we get an error "Unable to access the X Display, is $DISPLAY set properly?"
import sys
sys.path.append("/usr/local/bin/spock/bin")
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/kalman/spock_development_new_structure_kalman_dev/srcPython")

import matplotlib.gridspec as gridspec
import numpy as np
from struct import *
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
from datetime import datetime, timedelta
from collections import *
from read_input_file import *
from read_output_file import *

import os, errno


# PARAMATERS TO SET UP BEFORE RUNNING THIS SCRIPT
save_results = 1 # set to 1 to save the plot
show_plots = 0 # set to 1 to show the plot (in interactive mode)
# end of PARAMATERS TO SET UP BEFORE RUNNING THIS SCRIPT


if show_plots == 1:
    plt.ion()

# Date to plot
start_date_input = sys.argv[1] # the date is fhte first and only argument when running this script

####################################
# # LOAD THE IR BRIGHTNESS TEMPERATURE
print "Downloading the IR brightness temperature image file..."
if start_date_input == "now":
    start_date_input = datetime.strftime( datetime.now(), "%Y-%m-%dT%H:%M:00")
start_date = (start_date_input.replace("-", "").replace("T", ""))[:-6] # change the format of the string
IR_filename = "merg_" + start_date + "_4km-pixel"
## Download the IR data from NCEP
os.system("wget ftp://ftp.cpc.ncep.noaa.gov/precip/global_full_res_IR/" + IR_filename + ".Z")
## Uncompress the file
if os.path.exists(IR_filename) ==  False: 
    os.system("uncompress " + IR_filename + ".Z")  # if this file has not been downloaded in the past
else:
    os.system("rm -f " + IR_filename + ".Z")

## Read IR data
try:
    IR_file = open(IR_filename, "rb")
except IOError as ioex: # if start_date is too recent (so if there is no data available at NCEP yet for start_date) then downlaod the most recent data available 
    # print 'errno:', ioex.errno
    # print 'err code:', errno.errorcode[ioex.errno]
    #print 'err message:', os.strerror(ioex.errno)
    start_date_input_save = start_date_input
    filename_all_files = "filename_all_files.txt"
    os.system("wget ftp://ftp.cpc.ncep.noaa.gov/precip/global_full_res_IR/ -O " + filename_all_files )
    file_all_files = open(filename_all_files)
    read_file_all_files = file_all_files.readlines()
    found_earliest_file = 0
    iup = 1
    while found_earliest_file == 0:
        if ('_4km-pixel.Z">merg_' in read_file_all_files[len(read_file_all_files) - iup]):
            found_earliest_file = 1
        else:                
            iup = iup + 1
    new_start_date = read_file_all_files[len(read_file_all_files) - iup].split('_4km-pixel.Z">merg_')[0].split('merg_')[-1]
    file_all_files.close()
    os.system("rm -f " + filename_all_files)
    start_date_input = datetime.strftime( datetime.strptime(new_start_date, "%Y%m%d%H"), "%Y-%m-%dT%H:%M:00")   
    start_date = (start_date_input.replace("-", "").replace("T", ""))[:-6] # change the format of the string
    IR_filename = "merg_" + start_date + "_4km-pixel"
    ## Download the IR data from NCEP
    os.system("wget ftp://ftp.cpc.ncep.noaa.gov/precip/global_full_res_IR/" + IR_filename + ".Z")
    ## Uncompress the file
    if os.path.exists(IR_filename) ==  False: # if this file has not been downloaded in the past
        os.system("uncompress " + IR_filename + ".Z")
    else:
        os.system("rm -f " + IR_filename + ".Z")
    IR_file = open(IR_filename, "rb")
    print "***! " + start_date_input_save  + " is too recent: the data for this date is not available at NOAA NCEP yet. The most recent data is on " +  start_date_input + ". Therefore, the start date has been set to " +  start_date_input + " !***\n"

print "Reading the IR brightness temperature image file..."
nb_lon = 9896
nb_lat = 3298

min_lon = 0.0182
dlon = 0.036378335
max_lat = 59.982
dlat = 0.036383683

lon_arr = np.arange(min_lon, min_lon + nb_lon*dlon, dlon)
lat_arr = np.arange(max_lat, max_lat - nb_lat*dlat, -dlat)

val = np.zeros([nb_lat, nb_lon])

for ilat in range(nb_lat):
#    print ilat, nb_lat-1
    for ilon in range(nb_lon):
        byte = IR_file.read(1)
        val[ilat, ilon] = unpack('B', byte)[0] + 75 # data have been scaled to fit into 1-byte by subtracting "75" from each datum so need to add 75 back

IR_file.close()
## Delete IR data
os.system("rm -f " + IR_filename)

####################################
print "Running SpOCK to get the CYGNSS sub-satellite locations over one orbit..."
# Run SpOCK to get the CYGNSS sub-satellite locations over one orbit (the orbit that starts at the date specified as an argument of this script)
start_date_spock = start_date_input 
end_date_spock = datetime.strftime( datetime.strptime(start_date_spock, "%Y-%m-%dT%H:%M:%S") + timedelta(minutes = 180), "%Y-%m-%dT%H:%M:%S") # end_date of SpOCK's simulation is one orbit after the start date called as the first argument of this script
os.system("spock_cygnss_spec_parallel_web.py " + start_date_spock + " " + end_date_spock + " spec")

####################################
# PLOT RESULTS
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3

ax_title = start_date_spock.replace("T", " ") + " to " + end_date_spock.replace("T", " ") 
y_label = 'Latitude '+ u'(\N{DEGREE SIGN})'
x_label = 'Longitude '+ u'(\N{DEGREE SIGN})'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig.suptitle('', y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_title(ax_title, weight = 'bold', fontsize  = (int)(fontsize_plot*1.1), y = 1.03)
ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

## Plot the 2D map (continents) with the IR brightness temperature in background

ax.imshow(val,cmap='Greys', extent=(min(lon_arr), max(lon_arr),min(lat_arr),max(lat_arr)))

m = Basemap( projection       = 'cyl',
	     llcrnrlon        = min(lon_arr) , #Lower Left  CoRNeR Longitude
	     urcrnrlon        = max(lon_arr)  , #Upper Right CoRNeR Longitude
	     llcrnrlat        = min(lat_arr)  , #Lower Left  CoRNeR Latitude
	     urcrnrlat        = max(lat_arr),   #Upper Right CoRNeR Latitude
	     resolution       = 'l'  ,
	     suppress_ticks   = False,
	     ax = ax,
	     )
m.drawcoastlines(linewidth=0.7, color='white')

## Plot the CYGNSS sub-satellite locations over one orbit (the orbit that starts at the date specified as an argument of this script)
### Read the SpOCK position output files
#### Read SpOCK main input file to figure out stuff to then read the output
input_filename = 'spock_spec_start_' + start_date_spock.replace(":", "_") + '_end_' + end_date_spock.replace(":", "_") + '.txt'
input_filename =    input_filename
var_in, var_in_order = read_input_file(input_filename)
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
dt = var_in[find_in_read_input_order_variables(var_in_order, 'dt')]; 
nb_steps = var_in[find_in_read_input_order_variables(var_in_order, 'nb_steps')]; 
nb_sc = var_in[find_in_read_input_order_variables(var_in_order, 'nb_sc')]; 
#### Read SpOCK output files
var_to_read = ["date","latitude", "longitude"]
lat = np.zeros([nb_sc, nb_steps]); lon = np.zeros([nb_sc, nb_steps])
spacecraft_list = []
point = namedtuple('point', ['x', 'y'])
color = namedtuple('color', 'red green blue')

satColors = ['lime', 'blue', 'indigo', 'purple', 'dodgerblue', 'steelblue', 'seagreen', 'limegreen']

nb_spec_pts = 4
interpolation_step = 1 # in second, interpolation step of find_specular_points.c (1 s usually)
nb_steps_interpolation = (int)((nb_steps-1) * dt / interpolation_step) +1
lon_spec, lat_spec = np.zeros([nb_spec_pts, nb_sc, nb_steps_interpolation]), np.zeros([nb_spec_pts, nb_sc, nb_steps_interpolation])
specular_list = []

nb_sc = 8 # !!!!!!!!!
label_arr = ['FM05', 'FM04', 'FM02', 'FM01', 'FM08', 'FM06', 'FM07', 'FM03'] 
for isc in range(nb_sc):
    ### Satellite locations
#     var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
#     date = var_out[find_in_read_input_order_variables(var_out_order, 'date')]
#     lat[isc, :] = var_out[find_in_read_input_order_variables(var_out_order, 'latitude')]
#     lon[isc, :] = var_out[find_in_read_input_order_variables(var_out_order, 'longitude')]

#     spacecraft = namedtuple('spacecraft',('name',) +  point._fields + ('point_plot',) + ('marker_spacecraft',))
#     spacecraft_list.append(spacecraft)
#     #  positions over one orbit
#     spacecraft_list[isc].x, spacecraft_list[isc].y =  m(lon[isc,:], lat[isc,:])
#     spacecraft_list[isc].marker_spacecraft = '.'
#     # point on the plot
#     spacecraft_list[isc].point_plot = m.scatter(spacecraft_list[isc].x, spacecraft_list[isc].y,  marker=spacecraft_list[isc].marker_spacecraft, color = satColors[isc], s = 20)

   ###Add in specular points
    time_spec_sublist = [] 
    file_specular = open(output_file_path_list[isc] + "specular_" + output_file_name_list[isc], "r")
    #file_specular = open(spec_dir + "/coverage/storm/coverage_specular_" + output_file_name_list[isc], "r")
    read_file_specular  = file_specular.readlines()
    # Nb of lines in the spec file header
    if (isc == 0):
        nb_lines_header_output_file_spec = 0
        while (read_file_specular[nb_lines_header_output_file_spec].split()[0] != "#START"):
            nb_lines_header_output_file_spec = nb_lines_header_output_file_spec + 1
        nb_lines_header_output_file_spec = nb_lines_header_output_file_spec + 1
    ispec_save = 0
    j = -1 # cbv
    while (ispec_save < len(read_file_specular)-1-nb_spec_pts):
        j = j + 1
        time_spec_sublist_temp_ini = read_file_specular[ispec_save+nb_lines_header_output_file_spec].split()[0] 
        time_spec_sublist_temp_ini = datetime.strptime(time_spec_sublist_temp_ini, "%Y-%m-%dT%H:%M:%S")
        time_spec_sublist.append(datetime.strftime(time_spec_sublist_temp_ini, "%Y-%m-%dT%H:%M:%S"))

        name_spec_sublist = []
        lon_spec[0,isc,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[1]) #CHANGED
#         if lon_spec[0,isc,j] > 180:
#             lon_spec[0,isc,j] = lon_spec[0,isc,j] - 360.
        lat_spec[0,isc,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[2]) #CHANGED

        ispec = 1
        while (datetime.strptime(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[0], "%Y-%m-%dT%H:%M:%S")  == time_spec_sublist_temp_ini):
            lon_spec[ispec,isc,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[1])
#             if lon_spec[ispec,isc,j] > 180:
#                 lon_spec[ispec,isc,j] = lon_spec[ispec,isc,j] - 360.
            lat_spec[ispec,isc,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[2])
            ispec = ispec + 1
            if (nb_lines_header_output_file_spec+ispec_save+ispec == len(read_file_specular) - 1):
                break
        ispec_save = ispec + ispec_save
    # if up to here we still ahve read the entire spec file
    if ( datetime.strptime(time_spec_sublist[-1], "%Y-%m-%dT%H:%M:%S") != datetime.strptime(read_file_specular[len(read_file_specular)-2].split()[0], "%Y-%m-%dT%H:%M:%S") ):        
        j = j + 1
        first_spec_of_last_time_step = +ispec_save
        time_spec_sublist_temp_ini = read_file_specular[first_spec_of_last_time_step].split()[0]
        time_spec_sublist_temp_ini = datetime.strptime(time_spec_sublist_temp_ini, "%Y-%m-%dT%H:%M:%S")
        time_spec_sublist.append(datetime.strftime(time_spec_sublist_temp_ini, "%Y-%m-%dT%H:%M:%S"))
        lon_spec[0,isc,j] = np.float(read_file_specular[first_spec_of_last_time_step].split()[1])
#         if lon_spec[0,isc,j] > 180:
#             lon_spec[0,isc,j] = lon_spec[0,isc,j] - 360.
        lat_spec[0,isc,j] = np.float(read_file_specular[first_spec_of_last_time_step].split()[2])
        ispec = 1
        while (datetime.strptime(read_file_specular[first_spec_of_last_time_step+ispec].split()[0], "%Y-%m-%dT%H:%M:%S")  == time_spec_sublist_temp_ini):
            lon_spec[ispec,isc,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[1])
#             if lon_spec[ispec,isc,j] > 180:
#                 lon_spec[ispec,isc,j] = lon_spec[ispec,isc,j] - 360.
            lat_spec[ispec,isc,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[2])
            if (first_spec_of_last_time_step+ispec < len(read_file_specular) - 1):
                ispec = ispec + 1
            else: 
                break
# Build the tuples for the visualization of the specular points
    for k in range(nb_spec_pts):
        specular = namedtuple('specular',('name',) +  point._fields  + color._fields + ('point_plot',))
        specular_list.append(specular)
# Add on plot the specular points over one orbit
        specular_list[k+isc*nb_spec_pts].x, specular_list[k+isc*nb_spec_pts].y =  m(lon_spec[k,isc,:], lat_spec[k,isc,:])
        if k == 0:
            specular_list[k+isc*nb_spec_pts].point_plot = m.scatter(specular_list[k+isc*nb_spec_pts].x, specular_list[k+isc*nb_spec_pts].y, marker='o', s = 0.1, color = satColors[isc], label = label_arr[isc])

        else:
            specular_list[k+isc*nb_spec_pts].point_plot = m.scatter(specular_list[k+isc*nb_spec_pts].x, specular_list[k+isc*nb_spec_pts].y, marker='o', s = 0.1, color = satColors[isc])


legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot,handlelength=0, handletextpad=0)
legend.get_title().set_fontsize(str(fontsize_plot))
        
for isc in range(nb_sc):
    legend.get_texts()[isc].set_color(satColors[isc]) # set the label the same color as the plot
    legend.legendHandles[isc].set_visible(False) # hide the line in the label



if save_results == 1:
    path_save_results = "/Users/cbv/cygnss/website/data/ir_brightness_temperature/"
    if nb_sc == 8:
        fig_save_name = start_date_spock.replace(":", "_") + '_CYGNSS_IR_all_CYG.png'
    else:
        fig_save_name = start_date_spock.replace(":", "_") + '_CYGNSS_IR.png'
    fig.savefig(path_save_results + fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

