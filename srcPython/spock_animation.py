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
# This script creates an animation of the sub satellite locations calculated by SpOCK on a 2D map from start_date until end_date. The propagated satellite is given by its NORAD ID.
# To run this script:
# python spock_animation.py start_date end_date norad_id
# where start_date is the date where to start the animation.
# if 'now' instead of a date for start_date then the start date is set to the current time (but the seconds are set to 0) and instead of an end date indicate the number of hours: python spock_animation.py now 5 41884 will run from now until now + 5 hours the sc 41884
# ASSUMPTIONS:
# - start_date and end_date must have the format YYYY-mm-ddTHH:MM:SS (for example: 2017-03-26T12:30:00)
# - all dates have to be UTC
# - see assumptions of spock_run_tle.py

import sys
#sys.path.append("/Users/cbv/Library/Enthought/Canopy_64bit/User/lib/python2.7/site-packages/")
import matplotlib
matplotlib.use("Agg") # without this line, when running this script from a cronjob we get an error "Unable to access the X Display, is $DISPLAY set properly?"

import matplotlib.gridspec as gridspec
import numpy as np
from struct import *
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
from datetime import datetime, timedelta
from collections import *
import os
from read_input_file import *
from read_output_file import *


# Run SpOCK to get the sub-satellite locations from start_date until end_date of satellite norad_id
start_date_str = sys.argv[1]
end_date_str =  sys.argv[2]
norad_id = sys.argv[3]
if start_date_str == "now":
    start_date_str = datetime.strftime(  datetime.utcnow(), "%Y-%m-%dT%H:%M:00")
    start_date_temp =  datetime.strptime(start_date_str, "%Y-%m-%dT%H:%M:%S")
    end_date_str = datetime.strftime( start_date_temp + timedelta(hours = np.float(sys.argv[2])), "%Y-%m-%dT%H:%M:00")
start_date = datetime.strptime(start_date_str, "%Y-%m-%dT%H:%M:%S")
end_date = datetime.strptime(end_date_str, "%Y-%m-%dT%H:%M:%S")

# Run SpOCK
spock_dir = "."
os.system("spock_run_tle.py " + start_date_str + " " + end_date_str + " " + norad_id) # !!!!!!!! remove 'python '

# Make the plot
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3


ax_title = start_date_str.replace("T", " ") + ' to ' + end_date_str.replace("T", " ")

y_label = 'Latitude '+ u'(\N{DEGREE SIGN})'
x_label = 'Longitude '+ u'(\N{DEGREE SIGN})'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig.suptitle('', y = 0.970,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_title(ax_title, weight = 'bold', fontsize  = (int)(fontsize_plot*1.1), y = 1.02)
ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

## Plot the 2D map (continents) with the IR brightness temperature in background
min_lon = 0; max_lon = 360
min_lat = -60; max_lat = 60
# dlon = 50
# dlat = 20
# lon_arr = np.arange(min_lon, min_lon + nb_lon*dlon, dlon)
# lat_arr = np.arange(max_lat, max_lat - nb_lat*dlat, -dlat)



m = Basemap( projection       = 'cyl',
	     llcrnrlon        = min_lon , #Lower Left  CoRNeR Longitude
	     urcrnrlon        = max_lon  , #Upper Right CoRNeR Longitude
	     llcrnrlat        = min_lat  , #Lower Left  CoRNeR Latitude
	     urcrnrlat        = max_lat,   #Upper Right CoRNeR Latitude
	     resolution       = 'l'  ,
	     suppress_ticks   = False,
	     ax = ax,
	     )

color_continents = [65,105,225]
color_continents = np.array(color_continents) / 256.
color_water  = [100,149,237]
color_water = np.array(color_water) / 256.
m.fillcontinents(color=tuple(color_continents),lake_color=tuple(color_water))
m.drawmapboundary(fill_color=tuple(color_water))

#m.drawcoastlines(linewidth=0.7, color=tuple(color_continents))


## Plot the sub-satellite locations over one orbit (the orbit that starts at the date specified as an argument of this script)
### Read the SpOCK position output files
#### Read SpOCK main input file to figure out stuff to then read the output
input_filename = spock_dir + '/spock_spec_start_' + start_date_str.replace(":", "_") + '_end_' + end_date_str.replace(":", "_") + '.txt'
var_in, var_in_order = read_input_file(input_filename)
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
dt = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')]; 
nb_steps = var_in[find_in_read_input_order_variables(var_in_order, 'nb_steps')]; 
nb_sc = var_in[find_in_read_input_order_variables(var_in_order, 'nb_sc')]; 
#### Read SpOCK output files
var_to_read = ["date","latitude", "longitude"]
lat = np.zeros([nb_sc, nb_steps]); lon = np.zeros([nb_sc, nb_steps])
spacecraft_list = []
    
point = namedtuple('point', ['x', 'y'])
color = namedtuple('color', 'red green blue')

satColors = ['lime']
label_arr = [norad_id]

for isc in range(nb_sc):
    ###Satellite locations
    var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
    date = var_out[find_in_read_input_order_variables(var_out_order, 'date')]
    lat[isc, :] = var_out[find_in_read_input_order_variables(var_out_order, 'latitude')]
    lon[isc, :] = var_out[find_in_read_input_order_variables(var_out_order, 'longitude')]

    spacecraft = namedtuple('spacecraft',('name',) +  point._fields + ('point_plot',) + ('marker_spacecraft',))
    spacecraft_list.append(spacecraft)

    spacecraft_list[isc].marker_spacecraft = '.'

#whatever dt_output is in SpOCK, we plot only every 60 seconds for satellites, and every 10 seconds for specular points. The output of specular points is always 1 s time step. However, you need to make sure that the output for the satllites (secion #OUTPUT of main input file) is smaller or equal to 60 s.


dt_ani_sc = 60. # in s
dt_index_sc = (int)(dt_ani_sc / dt) # step for the index in lon and lat
nb_steps_ani_sc = nb_steps / dt_index_sc 
for istep in range(nb_steps_ani_sc):
    print "Step " + str(istep) + ' out of ' + str(nb_steps_ani_sc-1) 
    #  positions over one orbit
    for isc in range(nb_sc):
        spacecraft_list[isc].x, spacecraft_list[isc].y =  m(lon[isc,:istep*dt_index_sc:dt_index_sc], lat[isc,:istep*dt_index_sc:dt_index_sc])
    # point on the plot
        spacecraft_list[isc].point_plot = m.scatter(spacecraft_list[isc].x, spacecraft_list[isc].y,  marker=spacecraft_list[isc].marker_spacecraft, color = satColors[isc], s = 20, zorder = 4, label = label_arr[isc])

    date_map = ax.text((min_lon + max_lon)/2.,min_lat + (max_lat - min_lat)/30.,date[istep*dt_index_sc], fontsize = fontsize_plot, weight = 'bold', horizontalalignment = 'center')

    if (istep == 0):
        legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot,handlelength=0, handletextpad=0)
        legend.get_title().set_fontsize(str(fontsize_plot))

        for isc in range(len(legend.get_texts())):
            legend.get_texts()[isc].set_color(satColors[isc]) # set the label the same color as the plot
            legend.legendHandles[isc].set_visible(False) # hide the line in the label

    if os.path.isdir("animation_temp") == False:
        os.system("mkdir animation_temp")
    fig_save_name = 'animation_temp/' + (start_date_str + '_' + end_date_str + '_' + norad_id +'_' + str(istep) + '_' + str(nb_steps_ani_sc-1) + '.png').replace(":","_")
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

    # Remove the date from plot so that next one does not write on top of it but replaces it instead
    date_map.remove()

# # Create video from previous plots
nb_steps_spec_or_not = nb_steps_ani_sc

os.system('ffmpeg -r 10 -i ' + fig_save_name.split(norad_id + '_')[0] + norad_id + '_%d_' + str(nb_steps_spec_or_not-1) + '.png' + ' -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -vcodec libx264 -pix_fmt yuv420p ' +fig_save_name.split(norad_id + '_')[0].split('/')[-1] + norad_id + '.mp4')

print "rm -Rf "+ '/'.join(fig_save_name.split('/')[:-1])
os.system("rm -Rf "+ '/'.join(fig_save_name.split('/')[:-1]))
