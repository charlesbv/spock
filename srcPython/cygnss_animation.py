# This script creates an animation of the CYGNSS sub satellite locations calculated by SpOCK on a 2D map from start_date until end_date.
# To run this script:
# python cygnss_animation.py start_date end_date
# where start_date is the date where to start the animation.
# if 'now' instead of a date for start_date then the start date is set to the current time (but the seconds are set to 0) and instead of an end date indicate the number of hours: python cygnss_animation.py now 5 will run from now until now + 5 hours
# Instead of showing the sub-satellite locations, you can show the specular point locations. It will show only the current ones, not keeping with the past ones on the map though (so the plot doesn't get too busy). To do, after start_date write "spec": python cygnss_animation.py start_date spec
# ASSUMPTIONS:
# - start_date and end_date must have the format YYYY-mm-ddTHH:MM:SS (for example: 2017-03-26T12:30:00)
# - all dates have to be UTC
# - see assumptions of spock_cygnss_spec.py

import sys
#sys.path.append("/Users/cbv/Library/Enthought/Canopy_64bit/User/lib/python2.7/site-packages/")
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/kalman/spock_development_new_structure_kalman_dev/srcPython")
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
from cygnss_read_spock_spec import *

# Show specular points instead of sub-satellite locations
show_spec = 0
if len(sys.argv) >=4:
    if (sys.argv[3] == "spec"):
        show_spec = 1
# Run SpOCK to get the CYGNSS sub-satellite locations from start_date until end_date
start_date_str = sys.argv[1]
end_date_str =  sys.argv[2]
if start_date_str == "now":
    start_date_str = datetime.strftime(  datetime.utcnow(), "%Y-%m-%dT%H:%M:00")
    start_date_temp =  datetime.strptime(start_date_str, "%Y-%m-%dT%H:%M:%S")
    end_date_str = datetime.strftime( start_date_temp + timedelta(hours = np.float(sys.argv[2])), "%Y-%m-%dT%H:%M:00")
start_date = datetime.strptime(start_date_str, "%Y-%m-%dT%H:%M:%S")
end_date = datetime.strptime(end_date_str, "%Y-%m-%dT%H:%M:%S")

# Run SpOCK
spock_dir = "."
if show_spec == 0:
    os.system("spock_cygnss_spec_parallel_web.py " + start_date_str + " " + end_date_str)
else:
    os.system("spock_cygnss_spec_parallel_web.py " + start_date_str + " " + end_date_str + " spec")
    
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


## Plot the CYGNSS sub-satellite locations over one orbit (the orbit that starts at the date specified as an argument of this script)
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
if show_spec == 0:
    var_to_read = ["date","latitude", "longitude"]
    lat = np.zeros([nb_sc, nb_steps]); lon = np.zeros([nb_sc, nb_steps])
    spacecraft_list = []
else:
    specular_list = []
    
point = namedtuple('point', ['x', 'y'])
color = namedtuple('color', 'red green blue')

satColors = ['lime', 'blue', 'indigo', 'purple', 'dodgerblue', 'steelblue', 'seagreen', 'limegreen']
nb_sc = 8 # !!!!!!!!!
label_arr = ['FM05', 'FM04', 'FM02', 'FM01', 'FM08', 'FM06', 'FM07', 'FM03']
if show_spec == 1:
    lon_spec = []; lat_spec = []; gain_spec = []; gps_spec = []; date_spec = []
for isc in range(nb_sc):
    ###Satellite locations
    if show_spec == 0:
        var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
        date = var_out[find_in_read_input_order_variables(var_out_order, 'date')]
        lat[isc, :] = var_out[find_in_read_input_order_variables(var_out_order, 'latitude')]
        lon[isc, :] = var_out[find_in_read_input_order_variables(var_out_order, 'longitude')]

        spacecraft = namedtuple('spacecraft',('name',) +  point._fields + ('point_plot',) + ('marker_spacecraft',))
        spacecraft_list.append(spacecraft)

        spacecraft_list[isc].marker_spacecraft = '.'
    else:
        which_sc = isc
        cyg = format(isc + 1, "02")
        filename_spec_spock = spock_dir + "/" + output_file_path_list[which_sc] + "specular_" + output_file_name_list[which_sc]
        date_spec_this_sc, lon_spec_this_sc, lat_spec_this_sc, gain_spec_this_sc, gps_spec_this_sc, normpower_useless,x_cyg_useless, y_cyg_useless, z_cyg_useless, x_gps_useless, y_gps_useless, z_gps_useless, x_spec_useless, y_spec_useless, z_spec_useless = cygnss_read_spock_spec(filename_spec_spock)
        lon_spec.append(lon_spec_this_sc)
        lat_spec.append(lat_spec_this_sc)
        gain_spec.append(gain_spec_this_sc)
        gps_spec.append(gps_spec_this_sc)
        date_spec.append(date_spec_this_sc) # should be the same for all sc

        specular = namedtuple('specular',('name',) +  point._fields + ('point_plot',) + ('marker_specular',))
        specular_list.append(specular)

        specular_list[isc].marker_specular = '.'

        if isc == 0:
            nb_steps_spec = len(date_spec_this_sc)
        if len(date_spec_this_sc) > nb_steps_spec:
            nb_steps_spec = len(date_spec_this_sc)


#whatever dt_output is in SpOCK, we plot only every 60 seconds for satellites, and every 10 seconds for specular points. The output of specular points is always 1 s time step. However, you need to make sure that the output for the satllites (secion #OUTPUT of main input file) is smaller or equal to 60 s.

if show_spec == 0:
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

        date_map = ax.text((min_lon + max_lon)/2.,min_lat + (max_lat - min_lat)/30.,date[istep*dt_index_sc][:16]+ ':00', fontsize = fontsize_plot, weight = 'bold', horizontalalignment = 'center') # cheating with ':00'

        if (istep == 0):
            legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot,handlelength=0, handletextpad=0)
            legend.get_title().set_fontsize(str(fontsize_plot))

            for isc in range(len(legend.get_texts())):
                legend.get_texts()[isc].set_color(satColors[isc]) # set the label the same color as the plot
                legend.legendHandles[isc].set_visible(False) # hide the line in the label

        fig_save_name = 'animation_temp/' + (start_date_str + '_' + end_date_str + '_CYGNSS_' + str(istep) + '_' + str(nb_steps_ani_sc-1) + '.png').replace(":","_")
        fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

        # Remove the date from plot so that next one does not write on top of it but replaces it instead
        date_map.remove()
else:
    dt_ani_spec = 10. # in s
    dt_index_spec = (int) (dt_ani_spec) # this is directly dt_ani_spec because the time step of the output file for the spec is always 1 s
    nb_steps_ani_spec = nb_steps_spec / dt_index_spec

## Read specular point files
    for istep in range(nb_steps_ani_spec):
        print "Step " + str(istep) + ' out of ' + str(nb_steps_ani_spec-1) 
        #  positions over one orbit
        for isc in range(nb_sc):
            specular_list[isc].x, specular_list[isc].y =  m(lon_spec[isc][istep*dt_index_spec], lat_spec[isc][istep*dt_index_spec])
        # point on the plot
            specular_list[isc].point_plot = m.scatter(specular_list[isc].x, specular_list[isc].y,  marker=specular_list[isc].marker_specular, color = satColors[isc], s = 20, zorder = 4, label = label_arr[isc])

        date_map = ax.text((min_lon + max_lon)/2.,min_lat + (max_lat - min_lat)/30.,date_spec[isc][istep*dt_index_spec], fontsize = fontsize_plot, weight = 'bold', horizontalalignment = 'center')

        if (istep == 0):
            legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="SC #", fontsize = fontsize_plot,handlelength=0, handletextpad=0)
            legend.get_title().set_fontsize(str(fontsize_plot))

            for isc in range(len(legend.get_texts())):
                legend.get_texts()[isc].set_color(satColors[isc]) # set the label the same color as the plot
                legend.legendHandles[isc].set_visible(False) # hide the line in the label

        fig_save_name = 'animation_temp/' + (start_date_str + '_' + end_date_str + '_CYGNSS_' + str(istep) + '_' + str(nb_steps_ani_spec-1) + '.png').replace(":","_")
        fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

        # Remove the date from plot so that next one does not write on top of it but replaces it instead
        date_map.remove()
        if show_spec == 1:
            for isc in range(nb_sc):
                specular_list[isc].point_plot.remove()

# # Create video from previous plots
if show_spec == 0:
    nb_steps_spec_or_not = nb_steps_ani_sc
else:
    nb_steps_spec_or_not = nb_steps_ani_spec
os.system('ffmpeg -r 10 -i ' + fig_save_name.split('CYGNSS_')[0] + 'CYGNSS_%d_' + str(nb_steps_spec_or_not-1) + '.png' + ' -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -vcodec libx264 -pix_fmt yuv420p ' + fig_save_name.split('CYGNSS_')[0] + 'CYGNSS.mp4')

if show_spec == 0:
    os.system("mv " + fig_save_name.split('CYGNSS_')[0] + 'CYGNSS.mp4 ../../../data/animation/satellites')
    os.system("mv " + fig_save_name.split('CYGNSS_')[0] + 'CYGNSS_' + str(nb_steps_spec_or_not-1) +'_' + str(nb_steps_spec_or_not-1) + '.png ../../../data/animation/satellites/ani_sc.png')
else:
    os.system("mv " + fig_save_name.split('CYGNSS_')[0] + 'CYGNSS.mp4 ../../../data/animation/specular_points')
    os.system("mv " + fig_save_name.split('CYGNSS_')[0] + 'CYGNSS_' + str(nb_steps_spec_or_not-1) +'_' + str(nb_steps_spec_or_not-1) + '.png ../../../data/animation/specular_points/ani_spec.png')
    
os.system("rm -f "+ '/'.join(fig_save_name.split('/')[:-1]) + "/*png")
