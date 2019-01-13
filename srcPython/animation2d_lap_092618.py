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

# This script makes an 2D animation of trajectories from e simulation input_filename between date_start and date_end in the region [lon_min, lat_min]/[lon_max, lat_max] with a time step of dt_ani
# It will also add the ground station masks if a ground station file was included in the file input_filename
# How to call the script: python animation2d input_filename date_start date_stop lat_min lat_max lon_min lon_max dt_ani is_cygnss
# date_start and date_stop: YYYY-MM-DDTHH:MM:SS. If -1 then this takes the start date of simu
# lon: -180 to 180 
# is_cygnss: 1 if the simulation inlcudes 8 CYGNSS satellites. 0 otherwise. This is to properly label each in sc in case is_cygnss is 1

import sys
#sys.path.append("/Users/cbv/Library/Enthought/Canopy_64bit/User/lib/python2.7/site-packages/")
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/srcPython")
#sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/kalman/spock_development_new_structure_kalman_dev/srcPython")
import matplotlib
matplotlib.use("Agg") # without this line, when running this script from a cronjob we get an error "Unable to access the X Display, is $DISPLAY set properly?"

import matplotlib.gridspec as gridspec
import numpy as np
from struct import *
from matplotlib import pyplot as plt
from cygnss_read_spock_spec_bin import *
from mpl_toolkits.basemap import Basemap, shiftgrid
from datetime import datetime, timedelta
from collections import *
import os
from read_input_file import *
from read_output_file import *
#from cygnss_read_spock_spec import *


def radius_for_tissot(dist_km):
    return np.rad2deg(dist_km/6378.) # this is calculating using the Haversine formula


if (os.path.isdir("ani") == False):
    os.system("mkdir ani")


input_filename = sys.argv[1] #'beaconshort.txt'
date_start = sys.argv[2] #-1 # if -1 then this takes the start date of simu 
date_stop = sys.argv[3] #-1 # if -1 then this takes the start date of simu 
lat_min = np.float(sys.argv[4]) #-90#-20#-90#
lat_max = np.float(sys.argv[5]) #90#-40#90#
lon_min = np.float(sys.argv[6]) #-180#80#-180# #-180 to 180
lon_max = np.float(sys.argv[7]) #180#90#180#
dt_ani = np.float(sys.argv[8]) #60.
is_cygnss = sys.argv[9]

var_in, var_in_order = read_input_file(input_filename)
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
dt = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')]; 
nb_steps = var_in[find_in_read_input_order_variables(var_in_order, 'nb_steps')]; 
nb_sc = var_in[find_in_read_input_order_variables(var_in_order, 'nb_sc')]; 
date_start_simu = var_in[find_in_read_input_order_variables(var_in_order, 'date_start')]
date_stop_simu = var_in[find_in_read_input_order_variables(var_in_order, 'date_stop')]
filename_ground_station = var_in[find_in_read_input_order_variables(var_in_order, 'filename_ground_station')]

if dt_ani < dt:
    print "***! The time step of the animation can't be smaller than the time step of the animation (" + str(dt) + " seconds). Therefore, it was set to " + str(dt) + " seconds. !***"
    dt_ani  = dt

if str(date_start) == '-1':
    date_start = date_start_simu
else:
    date_start = datetime.strptime(date_start, "%Y-%m-%dT%H:%M:%S")
if str(date_stop) == '-1':
    date_stop = date_stop_simu
else:
    date_stop = datetime.strptime(date_stop, "%Y-%m-%dT%H:%M:%S")

date_start_str = datetime.strftime(date_start,  "%Y-%m-%dT%H:%M:%S")
date_stop_str = datetime.strftime(date_stop,  "%Y-%m-%dT%H:%M:%S")
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3

y_label = 'Latitude '+ u'(\N{DEGREE SIGN})'
x_label = 'Longitude '+ u'(\N{DEGREE SIGN})'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig.suptitle('', y = 0.970,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in normal
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='normal') ## make the labels of the ticks in normal 

## Plot the 2D map (continents) 
m = Basemap( projection       = 'cyl',
	     llcrnrlon        = lon_min , #Lower Left  CoRNeR Longitude
	     urcrnrlon        = lon_max  , #Upper Right CoRNeR Longitude
	     llcrnrlat        = lat_min  , #Lower Left  CoRNeR Latitude
	     urcrnrlat        = lat_max,   #Upper Right CoRNeR Latitude
	     resolution       = 'l'  ,
	     suppress_ticks   = False,
	     ax = ax,
	     )

# color_continents = [65,105,225]
# color_continents = np.array(color_continents) / 256.
# color_water  = [100,149,237]
# color_water = np.array(color_water) / 256.
# m.fillcontinents(color=tuple(color_continents),lake_color=tuple(color_water))
# m.drawmapboundary(fill_color=tuple(color_water))

m.drawcoastlines(linewidth=0.3, color = 'blue')

var_to_read = ["date","latitude", "longitude"]
spacecraft_list = []
point = namedtuple('point', ['x', 'y'])
color = namedtuple('color', 'red green blue')

satColors = ['black', 'blue', 'red', 'mediumorchid', 'dodgerblue', 'magenta', 'darkgreen', 'limegreen'] #['lime', 'blue', 'indigo', 'purple', 'dodgerblue', 'steelblue', 'seagreen', 'limegreen']
nb_sc = 8 # !!!!!!!!!
if is_cygnss == '1':
    label_arr = ['FM05', 'FM04', 'FM02', 'FM01', 'FM08', 'FM06', 'FM07', 'FM03']
    label_arr_conversion = [3, 2, 7, 1, 0, 5, 6, 4]
    

for isc in range(nb_sc):
    var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
    if isc == 0:
        date = var_out[find_in_read_input_order_variables(var_out_order, 'date')]
        nb_steps = len(date)
        lat = np.zeros([nb_sc, nb_steps]); lon = np.zeros([nb_sc, nb_steps])
    lat[isc, :] = var_out[find_in_read_input_order_variables(var_out_order, 'latitude')]
    lon[isc, :] = var_out[find_in_read_input_order_variables(var_out_order, 'longitude')]
    
    spacecraft = namedtuple('spacecraft',('name',) +  point._fields + ('point_plot',) + ('marker_spacecraft',))
    spacecraft_list.append(spacecraft)

    spacecraft_list[isc].marker_spacecraft = '.'

    dt_index_sc = (int)(dt_ani / dt) # step for the index in lon and lat
    nb_steps_ani_sc = nb_steps / dt_index_sc 
for istep in range(1,nb_steps_ani_sc):
    ax_title = date[istep*dt_index_sc][:19] + ' UTC' # cheating with ':00'# str(date_start) + ' to ' + str(date_stop)
    ax.set_title(ax_title, weight = 'normal', fontsize  = (int)(fontsize_plot*1.1), y = 1.00)

    print "Step " + str(istep) + ' out of ' + str(nb_steps_ani_sc-1) 
    #  positions over one orbit
    for isc_temp in range(nb_sc):
        if is_cygnss == '1':
            isc = label_arr_conversion[isc_temp]
        else: 
            isc = isc_temp
        spacecraft_list[isc_temp].x, spacecraft_list[isc_temp].y =  m(lon[isc,istep*dt_index_sc], lat[isc,istep*dt_index_sc])
    # point on the plot                                                        
        if is_cygnss == '1':
            spacecraft_list[isc_temp].point_plot = m.scatter(spacecraft_list[isc_temp].x, spacecraft_list[isc_temp].y,  marker=spacecraft_list[isc_temp].marker_spacecraft, color = satColors[isc], s = 200, zorder = 4, label = label_arr[isc])
        else:
            spacecraft_list[isc_temp].point_plot = m.scatter(spacecraft_list[isc_temp].x, spacecraft_list[isc_temp].y,  marker=spacecraft_list[isc_temp].marker_spacecraft, color = satColors[isc], s = 200, zorder = 4, label = str(isc + 1))

    #date_map = ax.text((lon_min + lon_max)/2.,lat_min + (lat_max - lat_min)/30.,date[istep*dt_index_sc][:16]+ ':00', fontsize = fontsize_plot, weight = 'bold', horizontalalignment = 'center') # cheating with ':00'

    if (istep == 1):
        legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="", fontsize = fontsize_plot,handlelength=0, handletextpad=0)
        #legend.get_title().set_fontsize(str(fontsize_plot))

        for isc_temp in range(len(legend.get_texts())):
            if is_cygnss == '1':
                isc = label_arr_conversion[isc_temp]
            else:
                isc = isc_temp
            legend.get_texts()[isc_temp].set_color(satColors[isc]) # set the label the same color as the plot
            legend.legendHandles[isc_temp].set_visible(False) # hide the line in the label

    m.tissot(gs_lon[idx], gs_lat[idx], radius_for_tissot(2000), 100, linestyle='dashdot', fill=False, edgecolor='grey')

    fig_save_name = 'ani/' + str(istep) + '_' + str(nb_steps-1) + '.png'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    for isc in range(nb_sc):
        spacecraft_list[isc].point_plot.remove()

    # Remove the date from plot so that next one does not write on top of it but replaces it instead
    #date_map.remove()

os.system('ffmpeg -r 10 -i ani/%d_' +  str(nb_steps-1) + '.png -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -vcodec libx264 -pix_fmt yuv420p ' + input_filename.replace(".txt", "") + "_" + date_start_str.replace("-", "_").replace(":", "_") +  "_to_" +  date_stop_str.replace("-", "_").replace(":", "_") + ".mp4")

