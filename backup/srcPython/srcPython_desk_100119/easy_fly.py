import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import FixedLocator
import h5py
from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
import matplotlib.animation as animation
from netCDF4 import Dataset
import os
from collections import *
import time
import datetime
from itertools import groupby
from operator import itemgetter

import urllib
import gdal
from gdalconst import *

import matplotlib.gridspec as gridspec
from get_ellipse_coords import *
from radius_for_tissot import *

from read_input_file import *


# NOTE: this code is not completely ready to be run if there is more than one storm
# ASSUMPTION: if running storm coverage, then the storm file can not start later than the satellites files (but it can start earlier)

plt.ion() #DO NOT UNCOMMENT THIS LINE AND THE LINE BELOW OTHERWISE THE ANIMATION DOES NOT WORK


input_filename = 'cygnss_example_visu.txt'
step_visualization = 1

zoom = 0
nb_storms = 1
show_clouds = 0
show_storm = 1 
show_graph = 0 #you can't show graph if you don't plot coverage. but you can plot coverage and not show graph
plot_coverage = 0

#########################################################################################
################################################################################## 2D MAP
# Figure
height_fig_with_map = 14
width_fig_with_map = 8
fig_with_map = plt.figure(figsize=(height_fig_with_map, width_fig_with_map))
background_color = (0/555,76./255,153/255.)
fig_with_map.set_facecolor(background_color)
# Axes
gs_with_map = gridspec.GridSpec(1, show_graph+1)
gs_with_map.update(left=0.04, right=0.99, top = 0.98,bottom = 0.04,hspace=0.08,wspace = 0.06)
ax_map = fig_with_map.add_subplot(gs_with_map[0, show_graph])
#width_half = 0.90
#ax_map = plt.axes([0.05,0.05,0.9,0.9])
min_lon = [-180, -90] #0
min_lat = [-90, 10] #-10
max_lon = [180, -50] #60
max_lat = [90,50]#40
step_lon = [60,10]#10
step_lat = [30, 10]#10
# array_lon = [['180W', '120W', '60W', '0', '60E', '120E', '180E'],['90W', '80W', '70W', '60W', '50W']]
# array_lat = [['90S', '60S', '30S', 'EQ', '30N', '60N', '90N'],['10N', '20N', '30N','40N','50N']]
array_lon = [ str(ss) for ss in np.arange(min_lon[zoom], max_lon[zoom],step_lon[zoom]) ]
array_lat = [ str(ss) for ss in np.arange(min_lat[zoom], max_lat[zoom],step_lat[zoom]) ]
ax_map.xaxis.set_major_locator(FixedLocator(np.arange(min_lon[zoom], max_lon[zoom]+1, step_lon[zoom])))
ax_map.yaxis.set_major_locator(FixedLocator(np.arange(min_lat[zoom], max_lat[zoom]+1, step_lat[zoom])))
ax_map.set_xticklabels(array_lon, color = 'w')
ax_map.set_yticklabels(array_lat, color = 'w')
# ax_map.set_xticklabels([])
# ax_map.set_yticklabels([])


m = Basemap( projection       = 'cyl',
	     llcrnrlon        = min_lon[zoom] , #Lower Left  CoRNeR Longitude
	     urcrnrlon        = max_lon[zoom]  , #Upper Right CoRNeR Longitude
	     llcrnrlat        = min_lat[zoom]  , #Lower Left  CoRNeR Latitude
	     urcrnrlat        = max_lat[zoom],   #Upper Right CoRNeR Latitude
	     resolution       = 'l'  ,
	     suppress_ticks   = False,
	     ax = ax_map,
	     )
m.drawcoastlines(linewidth=0.7, color='blue')

#########################################################################################
################################################################################## CLOUDS
if (show_clouds == 1):
    # # Everything between the double lines of '#' are the insert from other code
    # # The newest GOES images are always called 'latest':

    # #filnam_e = 'GoesEast1V_latest.tif'
    # # Now get the most current files from the website:
    # #testfile_e = urllib.URLopener()
    # #testfile_e.retrieve("ftp://satepsanone.nesdis.noaa.gov/GIS/GOESeast/"+filnam_e, filnam_e)

    filnam_e = 'IRCOMP_20151002_2100.tif' #'GoesEast1V2781845.tif'#'COMPOSITE0001.tif'

    # # PLOTTING GOES EAST
    ds = gdal.Open(filnam_e)
    
    #Inquire about the projection (may not be required)
    proj = ds.GetProjection()

    # Use GetGeoTransform() to obtain the geoinformation
    gt = ds.GetGeoTransform()

# Examine the extents of the image based on pixel resolution
    width  = ds.RasterXSize
    height = ds.RasterYSize
    xres   = gt[1]
    yres   = gt[5]

# Get the edge coordinates (note, values are already mid-pixel, so there is no need to shift by xres, yres)
    minx = gt[0]
    miny = gt[3] + width*gt[4] + height*yres 
    maxx = gt[0] + width*xres + height*gt[2]
    maxy = gt[3] 

# Now read in the pixel values of the whole image
    band = ds.GetRasterBand(1)
    image = band.ReadAsArray(0,0,width,height)
# value = image[42,94] <-- For example, the pixel brightness value at 42,94

# Establish a lat lon array for this image
    lonE = np.arange(minx,maxx,xres)  
    latE = np.arange(miny,maxy,-1*yres) # <-- Has to be -1*yres to be positive

#Generate a lon x lat grid for the tif file
    xE, yE = np.meshgrid(lonE, latE)

#Now the image has to be rotated 180 degrees before plotting
    image_r = np.zeros([height,width])
    for i in range(height):
        image_r[i,:] = image[(height-1)-i,:]

    Tf  = m.pcolormesh(xE,yE,image_r,cmap='Greys_r')



#########################################################################################
########################################################## SPACECRAFT AND SPECULAR POINTS
propagator_directory_list = os.getcwd().split('/')[0:-2]
propagator_directory = ""
for i in range(len(propagator_directory_list)):
    if (i > 0):
        propagator_directory = propagator_directory + "/" + propagator_directory_list[i]

input_filename = propagator_directory + '/input/main_input/' + input_filename

input_variables, order_input_variables = read_input_file(input_filename)
date_start = input_variables[0]; date_stop = input_variables[1]; dt = input_variables[2]; nb_steps = input_variables[3]; nb_satellites = input_variables[4]; 
output_file_path_propagator = input_variables[6]; output_filename_propagator = input_variables[7];
gps_name = input_variables[5]; 
nb_storms = input_variables[9]; storm_name = input_variables[10]

#Read the output files
nb_spec_pts = 4; color_spec = ['c', 'r','b', 'k','g','m','y','w']
interpolation_step = 1 # in second, interpolation step of find_specular_points.c (1 s usually)
nb_steps_interpolation = (int)((nb_steps-1) * dt / interpolation_step) +1

lon_sat = np.zeros([nb_satellites, nb_steps_interpolation])
lat_sat = np.zeros([nb_satellites, nb_steps_interpolation])
ecef_sat = np.zeros([nb_satellites, nb_steps_interpolation, 3])
time_sat = []
lon_spec, lat_spec, gain_spec, distance_spec_to_storm, specular_over_storm = np.zeros([nb_spec_pts, nb_satellites, nb_steps_interpolation]), np.zeros([nb_spec_pts, nb_satellites, nb_steps_interpolation]), np.zeros([nb_spec_pts, nb_satellites, nb_steps_interpolation]), np.zeros([nb_spec_pts, nb_satellites, nb_steps_interpolation]), np.zeros([nb_spec_pts, nb_satellites, nb_steps_interpolation])
ecef_spec = np.zeros([nb_spec_pts, nb_satellites, nb_steps_interpolation, 3])
name_spec = []
time_spec = []
date = []
list_output_variables_to_read = ["longitude","latitude"]
point = namedtuple('point', ['x', 'y'])
color = namedtuple('color', 'red green blue')
spacecraft_list = []
specular_list = []
init_index = 3 #27, 5
for i in range(nb_satellites):
    time_spec_sublist = []
    name_spec_between_list_and_sublist = []
        ### SATELLITES ###
    output_filename = output_file_path_propagator[i] + "/interpolated_position_LLA_ECEF_" + output_filename_propagator[i]
    output_file = open(output_filename, "r")
    read_output_file = output_file.readlines()
    nb_lines_header_output_file_sat = 0
    for j in range(nb_steps_interpolation):
        lon_sat[i,j] = np.float(read_output_file[j+nb_lines_header_output_file_sat].split()[1])
        if (lon_sat[i,j] > 180):
            lon_sat[i,j] = lon_sat[i,j] - 360.
        lat_sat[i,j] = np.float(read_output_file[j+nb_lines_header_output_file_sat].split()[2])
        ecef_sat[i,j,0] = np.float(read_output_file[j+nb_lines_header_output_file_sat].split()[4])
        ecef_sat[i,j,1] = np.float(read_output_file[j+nb_lines_header_output_file_sat].split()[5])
        ecef_sat[i,j,2] = np.float(read_output_file[j+nb_lines_header_output_file_sat].split()[6])
        time_sat.append(read_output_file[j+nb_lines_header_output_file_sat].split()[0])

# Build the tuples for the visualization of the satellites
    spacecraft = namedtuple('spacecraft',('name',) +  point._fields + color._fields + ('point_plot',) + ('marker_spacecraft',))
    spacecraft_list.append(spacecraft)
    name_temp = output_filename_propagator[i].replace(".txt","")
    spacecraft_list[i].name = name_temp
# initial position
    spacecraft_list[i].x, spacecraft_list[i].y =  m(lon_sat[i,init_index], lat_sat[i,init_index])
# color of the spacecraft in the visualization
    spacecraft_list[i].red   = 1
    spacecraft_list[i].green = 0
    spacecraft_list[i].blue  = 1
    color_sat = (spacecraft_list[i].red, spacecraft_list[i].green, spacecraft_list[i].blue)
    spacecraft_list[i].marker_spacecraft = 'o'
    # point on the plot
#	spacecraft_list[i].point_plot = m.plot(spacecraft_list[i].x, spacecraft_list[i].y,  marker=spacecraft_list[i].marker_spacecraft, markersize=15,color = (spacecraft_list[i].red, spacecraft_list[i].green, spacecraft_list[i].blue))[0]
    spacecraft_list[i].point_plot = m.plot([],[],  marker=spacecraft_list[i].marker_spacecraft, markersize=15,color = color_sat)[0]
#	spacecraft_list[i].point_plot needs to be a line object for the animation, that's why we use [0] (at the end of the line) because othwerwise it'd be a list of lines

#         ### SPECULAR POINTS ###
    spec_dir = ""
    for j in range(len(output_file_path_propagator[i].split('/'))-2):
        if (j > 0):
            spec_dir = spec_dir + "/" + output_file_path_propagator[i].split('/')[j]
    file_specular = open(spec_dir + "/coverage/storm/coverage_specular_" + output_filename_propagator[i], "r")
    read_file_specular  = file_specular.readlines()
    # Nb of lines in the spec file header
    if (i == 0):
        nb_lines_header_output_file_spec = 0
        while (read_file_specular[nb_lines_header_output_file_spec].split()[0] != "#START"):
            nb_lines_header_output_file_spec = nb_lines_header_output_file_spec + 1
        nb_lines_header_output_file_spec = nb_lines_header_output_file_spec + 1
    ispec_save = 0
    ## the output of find_specular_points.c does not start at the same time as the propagation 
    j = 0
    while ( datetime.strptime(read_file_specular[nb_lines_header_output_file_spec].split()[0], "%Y-%m-%dT%H:%M:%S") != datetime.strptime(time_sat[j], "%Y-%m-%dT%H:%M:%S") ):
        j = j +1
        time_spec_sublist.append('')
        name_spec_between_list_and_sublist.append('')
    j = j-1
    while (ispec_save < len(read_file_specular)-1-nb_spec_pts):
        j = j + 1
        time_spec_sublist_temp_ini = read_file_specular[ispec_save+nb_lines_header_output_file_spec].split()[0] 
        time_spec_sublist_temp_ini = datetime.strptime(time_spec_sublist_temp_ini, "%Y-%m-%dT%H:%M:%S")
        time_spec_sublist.append(datetime.strftime(time_spec_sublist_temp_ini, "%Y-%m-%dT%H:%M:%S"))
        name_spec_sublist = []
        lon_spec[0,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[4])
        if lon_spec[0,i,j] > 180:
            lon_spec[0,i,j] = lon_spec[0,i,j] - 360.
        lat_spec[0,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[5])
        ecef_spec[0,i,j,0] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[1])
        ecef_spec[0,i,j,1] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[2])
        ecef_spec[0,i,j,2] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[3])
        distance_spec_to_storm[0,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[8])
        specular_over_storm[0,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[9])
        name_spec_sublist.append(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[7])
        ispec = 1
        while (datetime.strptime(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[0], "%Y-%m-%dT%H:%M:%S")  == time_spec_sublist_temp_ini):
            lon_spec[ispec,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[4])
            if lon_spec[ispec,i,j] > 180:
                lon_spec[ispec,i,j] = lon_spec[ispec,i,j] - 360.
            lat_spec[ispec,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[5])
            ecef_spec[ispec,i,j,0] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[1])
            ecef_spec[ispec,i,j,1] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[2])
            ecef_spec[ispec,i,j,2] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[3])
            distance_spec_to_storm[ispec,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[8])
            specular_over_storm[ispec,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[9])
            name_spec_sublist.append(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[7])
            ispec = ispec + 1
            if (nb_lines_header_output_file_spec+ispec_save+ispec == len(read_file_specular)):
                break
        ispec_save = ispec + ispec_save
        name_spec_between_list_and_sublist.append(name_spec_sublist)
    # if up to here we still ahve read the entire spec file
    if ( datetime.strptime(time_spec_sublist[-1], "%Y-%m-%dT%H:%M:%S") != datetime.strptime(read_file_specular[len(read_file_specular)-1].split()[0], "%Y-%m-%dT%H:%M:%S") ):
        j = j + 1
        first_spec_of_last_time_step = +ispec_save
        time_spec_sublist_temp_ini = read_file_specular[first_spec_of_last_time_step].split()[0]
        time_spec_sublist_temp_ini = datetime.strptime(time_spec_sublist_temp_ini, "%Y-%m-%dT%H:%M:%S")
        time_spec_sublist.append(datetime.strftime(time_spec_sublist_temp_ini, "%Y-%m-%dT%H:%M:%S"))
        name_spec_sublist = []
        lon_spec[0,i,j] = np.float(read_file_specular[first_spec_of_last_time_step].split()[4])
        if lon_spec[0,i,j] > 180:
            lon_spec[0,i,j] = lon_spec[0,i,j] - 360.
        lat_spec[0,i,j] = np.float(read_file_specular[first_spec_of_last_time_step].split()[5])
        ecef_spec[0,i,j,0] = np.float(read_file_specular[first_spec_of_last_time_step].split()[1])
        ecef_spec[0,i,j,1] = np.float(read_file_specular[first_spec_of_last_time_step].split()[2])
        ecef_spec[0,i,j,2] = np.float(read_file_specular[first_spec_of_last_time_step].split()[3])
        distance_spec_to_storm[0,i,j] = np.float(read_file_specular[first_spec_of_last_time_step].split()[8])
        specular_over_storm[0,i,j] = np.float(read_file_specular[first_spec_of_last_time_step].split()[9])
        name_spec_sublist.append(read_file_specular[first_spec_of_last_time_step].split()[7])
        ispec = 1
        while (datetime.strptime(read_file_specular[first_spec_of_last_time_step+ispec].split()[0], "%Y-%m-%dT%H:%M:%S")  == time_spec_sublist_temp_ini):
            lon_spec[ispec,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[4])
            if lon_spec[ispec,i,j] > 180:
                lon_spec[ispec,i,j] = lon_spec[ispec,i,j] - 360.
            lat_spec[ispec,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[5])
            ecef_spec[ispec,i,j,0] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[1])
            ecef_spec[ispec,i,j,1] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[2])
            ecef_spec[ispec,i,j,2] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[3])
            distance_spec_to_storm[ispec,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[8])
            specular_over_storm[ispec,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[9])
            name_spec_sublist.append(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[7])
            if (first_spec_of_last_time_step+ispec < len(read_file_specular) - 1):
                ispec = ispec + 1
            else: 
                break
        name_spec_between_list_and_sublist.append(name_spec_sublist)

    ## the output of find_specular_points.c does not end at the same time as the propagation 
    j_end = 0
    while ( datetime.strptime(time_spec_sublist[-1-j_end], "%Y-%m-%dT%H:%M:%S") != datetime.strptime(time_sat[-1-j_end], "%Y-%m-%dT%H:%M:%S") ):
        j_end = j_end +1
        time_spec_sublist.append('')
        name_spec_between_list_and_sublist.append('')
    time_spec.append(time_spec_sublist)
# Build the tuples for the visualization of the specular points
    for k in range(nb_spec_pts):
        specular = namedtuple('specular',('name',) +  point._fields  + color._fields + ('point_plot',))
        specular_list.append(specular)
    # name 
    #		specular_list[i+j*nb_gps].name = specular_names[i]
    # initial position                    
        specular_list[k+i*nb_spec_pts].x, specular_list[k+i*nb_spec_pts].y =  m(lon_spec[k,i,init_index], lat_spec[k,i,init_index]) 
	# point on the plot
#        specular_list[k+i*nb_spec_pts].point_plot = m.plot(specular_list[k+i*nb_spec_pts].x, specular_list[k+i*nb_spec_pts].y, marker='o', markersize = 15, color = 'chartreuse', fillstyle = 'none', mew = 2)[0]#fillstyle='none',
        specular_list[k+i*nb_spec_pts].point_plot = m.plot([],[], marker='o', markersize = 15, color = 'chartreuse', fillstyle = 'none', mew = 2)[0]#fillstyle='none',
    name_spec.append(name_spec_between_list_and_sublist)

#########################################################################################
################################################################################## STORMS

if (show_storm == 1):
# # Read the storm interpolated files.  ASSUMPTION: all stroms are run at the same epoch and for the same amount of time

# ##
    spec_dir = ""
    for j in range(len(output_file_path_propagator[i].split('/'))-2):
        if (j > 0):
            spec_dir = spec_dir + "/" + output_file_path_propagator[i].split('/')[j]
    storm_filenames = []
    for i in range(nb_storms):
        storm_filenames.append(spec_dir + "/coverage/storm/"+ storm_name[i].split('.')[0]+"_interpolated."+storm_name[i].split('.')[1])
        
## Loop over the storms
### Read first filename to initialize the size of x, y, z
    file = open(storm_filenames[0],'r')
    a = file.readlines()
    nb_lines_storm = len(a) 
    line = ["" for x in range( nb_lines_storm )]
    time_storm = ["" for x in range( (int)( nb_lines_storm / step_visualization )+1 )]
    time_storm = ["" for x in range( (int)( nb_lines_storm / step_visualization )+1 )]
    lon_storm,lat_storm,radius_uncertainty_storm = np.zeros([(int)(nb_lines_storm/step_visualization)+1,nb_storms]), np.zeros([(int)(nb_lines_storm/step_visualization)+1,nb_storms]), np.zeros([(int)(nb_lines_storm/step_visualization)+1,nb_storms])
    file.close()

### Now read each file
    lon_storm_index_line = 2
    lat_storm_index_line = 1
    radius_storm_index_line = 6
    for i in range(nb_storms):
	file = open(storm_filenames[i],'r')
	a = file.readlines()
	for line_count in range(0,nb_lines_storm,step_visualization ):
		line[line_count] = (a[line_count]).split(' ')
		time_storm[(int)(line_count / step_visualization)] = line[line_count][0]
		lon_storm[(int)(line_count / step_visualization),i] = float( line[line_count][lon_storm_index_line]  )
		if (lon_storm[(int)(line_count / step_visualization),i] > 180.0):
			lon_storm[(int)(line_count / step_visualization),i] = lon_storm[(int)(line_count / step_visualization),i] - 360.0
		lat_storm[(int)(line_count / step_visualization),i] = float( line[line_count][lat_storm_index_line]  )
		radius_uncertainty_storm[(int)(line_count / step_visualization),i] = float( line[line_count][radius_storm_index_line]  )
	file.close()

# Build the tuples for the visualization of the storms
    point = namedtuple('point', ['x', 'y'])
    point_circle = namedtuple('point_circle', ['x_circle', 'y_circle'])
    color = namedtuple('color', 'red green blue')
    storm_list = []
    for i in range(nb_storms):
	storm = namedtuple('storm',('name',) +  point._fields + point_circle._fields + color._fields + ('point_plot',) + ('point_circle_plot',) + ('marker_storm',))
	storm_list.append(storm)
    # name
	storm_list[i].name = storm_name[i].split('.')[0].title()
    # initial position
	storm_list[i].x, storm_list[i].y =  m(lon_storm[0,i], lat_storm[0,i])
    # color of the storm in the visualization
        storm_list[i].red   = 0
        storm_list[i].green = 0
        storm_list[i].blue  = 0
        storm_list[i].marker_storm = 'o'
        ell = get_ellipse_coords(a=radius_for_tissot(radius_uncertainty_storm[0,i]), b=radius_for_tissot(radius_uncertainty_storm[0,i]), x=lon_storm[0,i], y=lat_storm[0,i])
        storm_list[i].point_plot = m.plot([], [], color='red',linewidth=3.5)[0]

# this plots all the circles of uncertainty during the entire animation. Like this, it gives a cone (that is the superimposition of all these circles of uncertainty)

        for j in range(0, (int)(3*24*3600),600):
            if float(a[j].split()[lon_storm_index_line])  > 180:
                m.tissot(float(a[j].split()[lon_storm_index_line]) - 360, float(a[j].split()[lat_storm_index_line]), radius_for_tissot(float(a[j].split()[radius_storm_index_line])), 256, facecolor='b', alpha = 0.002)
            else:
                m.tissot(float(a[j].split()[lon_storm_index_line]), float(a[j].split()[lat_storm_index_line]), radius_for_tissot(float(a[j].split()[radius_for_tissot])), 256, facecolor='b', alpha = 0.002)

        for j in range(0, (int)(1*24*3600),100):
            if float(a[j].split()[lon_storm_index_line])  > 180:
                m.tissot(float(a[j].split()[lon_storm_index_line]) - 360, float(a[j].split()[lat_storm_index_line]), radius_for_tissot(float(a[j].split()[radius_storm_index_line])), 256, facecolor='b', alpha = 0.002)
            else:
                m.tissot(float(a[j].split()[lon_storm_index_line]), float(a[j].split()[lat_storm_index_line]), radius_for_tissot(float(a[j].split()[radius_for_tissot])), 256, facecolor='b', alpha = 0.002)


    delay_storm = ( datetime.strptime( time_storm[0], "%Y-%m-%dT%H:%M:%S" ) - datetime.strptime( time_sat[0], "%Y-%m-%dT%H:%M:%S" ) ).days * 24 * 3600 + ( datetime.strptime( time_storm[0], "%Y-%m-%dT%H:%M:%S" ) - datetime.strptime( time_sat[0], "%Y-%m-%dT%H:%M:%S" ) ).seconds # delay in seconds between the start of the storm file and the start of the satellites files. !!!!!!!! ASSUMPTION: THE STORM FILE CAN NOT START LATER THAN THE SATELLITES FILES
    time_step_storm_file = ( datetime.strptime( time_storm[1], "%Y-%m-%dT%H:%M:%S" ) - datetime.strptime( time_storm[0], "%Y-%m-%dT%H:%M:%S" ) ).days * 24 * 3600 + ( datetime.strptime( time_storm[1], "%Y-%m-%dT%H:%M:%S" ) - datetime.strptime( time_storm[0], "%Y-%m-%dT%H:%M:%S" ) ).seconds
    delay_storm = - delay_storm / (time_step_storm_file) + 1 

#########################################################################################
################################################################################ COVERAGE

if ( (plot_coverage) | (show_graph == 1) ):
# Compute the times when each specular point of each satellites is on each storm
# example: list_times_spec_on_storm[0][2][1] are all the times when spec 1 of satellite 2 is on storm 0. If this spec of this sat flies 3 times in this storm then to have the list of time steps of the second time (for ex), print list_times_spec_on_storm[0][2][1][1]: this gives the indices of the min and max times of this second time... just try it, easier to understand!
    min_prob = 0.66
    max_prob = 1
    nb_times_at_least_one_spec_on_storm = np.zeros([nb_satellites, nb_storms])
    nb_days_run = (int) (np.ceil( ( (date_stop-date_start).days * 24 * 3600 + (date_stop-date_start).seconds ) / 3600.  / 24. )) 
    # nb_times_per_day_spec_on_storm = np.zeros([nb_satellites, nb_storms, nb_days_run, nb_spec_pts])
    # nb_times_per_day_at_least_one_spec_on_storm = np.zeros([nb_satellites, nb_storms, nb_days_run])
    # ex: nb_times_per_day_spec_on_storm[2,0,2,1] gives the number of times spec 1 of sat 2 is in storm 0 during the third day
    list_times_spec_on_storm = [] 
    for istorm in range(nb_storms):# Actually the rest of the code is not adapted for more than one storm (for now!)...
        list_times_spec_on_storm_sublist_storm = []
        for isat in range(nb_satellites):
            list_times_spec_on_storm_sublist_storm_sublist_sat = []
            nb_times_spec_on_storm_save = 0
            for ispec in range(nb_spec_pts):
                list_times_spec_on_storm_sublist_storm_sublist_sat_sublist_spec = []
                where_spec_on_storm = (np.where(specular_over_storm[ispec, isat,:] == 1))[0]
                if len(where_spec_on_storm) > 0:
                    for k, g in groupby(enumerate(where_spec_on_storm), lambda (i, x): i-x):
                        list_times_spec_on_storm_sublist_storm_sublist_sat_sublist_spec_temp = map(itemgetter(1), g)
                        list_times_spec_on_storm_sublist_storm_sublist_sat_sublist_spec.append([ min(list_times_spec_on_storm_sublist_storm_sublist_sat_sublist_spec_temp), max(list_times_spec_on_storm_sublist_storm_sublist_sat_sublist_spec_temp) ])
                list_times_spec_on_storm_sublist_storm_sublist_sat.append(list_times_spec_on_storm_sublist_storm_sublist_sat_sublist_spec)
                # if (len(list_times_spec_on_storm_sublist_storm_sublist_sat_sublist_spec) > nb_times_spec_on_storm_save):
                #     nb_times_spec_on_storm_save = len(list_times_spec_on_storm_sublist_storm_sublist_sat_sublist_spec)
                # # for each times a spec is in storm, look at which day it occurs
                # for j in range(len(list_times_spec_on_storm_sublist_storm_sublist_sat_sublist_spec)):
                #     time_index_spec_enters_in_storm = list_times_spec_on_storm_sublist_storm_sublist_sat_sublist_spec[j][0]
                #     seconds_after_start_run_spec_enters_in_storm = (int)(time_index_spec_enters_in_storm * interpolation_step)
                #     which_day_after_start_run_spec_enters_in_storm = seconds_after_start_run_spec_enters_in_storm / (24*3600)
                #     nb_times_per_day_spec_on_storm[isat, istorm, which_day_after_start_run_spec_enters_in_storm, ispec] = nb_times_per_day_spec_on_storm[isat, istorm, which_day_after_start_run_spec_enters_in_storm, ispec] + 1
            list_times_spec_on_storm_sublist_storm.append(list_times_spec_on_storm_sublist_storm_sublist_sat)
            nb_times_at_least_one_spec_on_storm[isat, istorm] = nb_times_spec_on_storm_save
            # for iday in range(nb_days_run):
            #     nb_times_per_day_at_least_one_spec_on_storm[isat, istorm,iday] = np.max(nb_times_per_day_spec_on_storm[isat, istorm, iday, :])
        list_times_spec_on_storm.append(list_times_spec_on_storm_sublist_storm)


# The goal of this block is to find the first spec that enters the storm (and the last spec that leaves the storm) for each fly over of the sat 
    min_delay_between_two_fly_over = 40 * 60 / interpolation_step # a satellite does not come back over the storm in less than 90 minutes (we take 40 minutes here)        
    list_all_start_spec_on_storm = []
    list_all_stop_spec_on_storm = []
    for istorm in range(nb_storms):
        list_all_start_spec_on_storm_sublist_sat =[]
        list_all_stop_spec_on_storm_sublist_sat = []
        for isat in range(nb_satellites):
            list_all_start_spec_on_storm_sublist_sat_sublist =[]
            list_all_stop_spec_on_storm_sublist_sat_sublist = []
            list_all_start_spec_on_storm_sublist_sat_sublist_temp = []
            list_all_stop_spec_on_storm_sublist_sat_sublist_temp = []
            for ispec in range(len(list_times_spec_on_storm[istorm][isat])):
                if len(list_times_spec_on_storm[istorm][isat][ispec]) > 0:
                    for itimes in range(len(list_times_spec_on_storm[istorm][isat][ispec])):
                        list_all_start_spec_on_storm_sublist_sat_sublist_temp.append(list_times_spec_on_storm[istorm][isat][ispec][itimes][0])
                        list_all_stop_spec_on_storm_sublist_sat_sublist_temp.append(list_times_spec_on_storm[istorm][isat][ispec][itimes][1])

            if len(list_all_stop_spec_on_storm_sublist_sat_sublist_temp)>0:
                j_save = 0
                stop_adding_stops = 0
                list_all_stop_spec_on_storm_sublist_sat_sublist.append(sorted(list_all_stop_spec_on_storm_sublist_sat_sublist_temp, reverse = True)[0])
                while (j_save < len(sorted(list_all_stop_spec_on_storm_sublist_sat_sublist_temp, reverse = True))):
                    j = 0 
                    while (sorted(list_all_stop_spec_on_storm_sublist_sat_sublist_temp, reverse = True)[j+j_save] > sorted(list_all_stop_spec_on_storm_sublist_sat_sublist_temp, reverse = True)[j_save] - min_delay_between_two_fly_over):
                        j = j + 1
                        if j + j_save >= len(sorted(list_all_stop_spec_on_storm_sublist_sat_sublist_temp, reverse = True)):
                            stop_adding_stops = 1
                            break
                    j_save = j + j_save
                    if stop_adding_stops == 0:
                        list_all_stop_spec_on_storm_sublist_sat_sublist.append(sorted(list_all_stop_spec_on_storm_sublist_sat_sublist_temp, reverse = True)[j_save])

                j_save = 0
                stop_adding_starts = 0
                list_all_start_spec_on_storm_sublist_sat_sublist.append(sorted(list_all_start_spec_on_storm_sublist_sat_sublist_temp)[0])
                while (j_save < len(sorted(list_all_start_spec_on_storm_sublist_sat_sublist_temp))):
                    j = 0 
                    while (sorted(list_all_start_spec_on_storm_sublist_sat_sublist_temp)[j+j_save] < sorted(list_all_start_spec_on_storm_sublist_sat_sublist_temp)[j_save] + min_delay_between_two_fly_over):
                        j = j + 1
                        if j + j_save >= len(sorted(list_all_start_spec_on_storm_sublist_sat_sublist_temp)):
                            stop_adding_starts = 1
                            break
                    j_save = j + j_save
                    if stop_adding_starts == 0:
                        list_all_start_spec_on_storm_sublist_sat_sublist.append(sorted(list_all_start_spec_on_storm_sublist_sat_sublist_temp)[j_save])
                list_all_start_spec_on_storm_sublist_sat.append(list_all_start_spec_on_storm_sublist_sat_sublist)
                list_all_stop_spec_on_storm_sublist_sat.append(sorted(list_all_stop_spec_on_storm_sublist_sat_sublist))
            else:
                list_all_start_spec_on_storm_sublist_sat.append([])
                list_all_stop_spec_on_storm_sublist_sat.append([])
        list_all_stop_spec_on_storm.append(list_all_stop_spec_on_storm_sublist_sat)
        list_all_start_spec_on_storm.append(list_all_start_spec_on_storm_sublist_sat)

    nb_times_per_day_at_least_one_spec_on_storm = np.zeros([nb_satellites, nb_storms, nb_days_run])  
    for istorm in range(nb_storms):
        for isat in range(nb_satellites):
            if len(list_all_start_spec_on_storm[istorm][isat]) > 0:
                for j in range(len(list_all_start_spec_on_storm[istorm][isat])):
                    time_index_first_spec_enters_in_storm = list_all_start_spec_on_storm[istorm][isat][j]
                    seconds_after_start_run_first_spec_enters_in_storm = (int)(time_index_first_spec_enters_in_storm * interpolation_step)
                    which_day_after_start_first_spec_enters_storm  =  seconds_after_start_run_first_spec_enters_in_storm / (24*3600)
                    nb_times_per_day_at_least_one_spec_on_storm[isat, istorm, which_day_after_start_first_spec_enters_storm] = nb_times_per_day_at_least_one_spec_on_storm[isat, istorm, which_day_after_start_first_spec_enters_storm] + 1


# Creates N new figures that each shows the coverage for the next day of 4 satellites. So if there are 8 satellites that are run for 3 days, then 2 (8/4) * 3 = 6 figures will be created
    if (plot_coverage == 1):
        istorm = 0 # !!!!!!!!!!!!!!!! PLOT COVERAGE CAN BE DONE ONLY FOR ON STORM (FOR NOW!)
        height_fig_coverage = height_fig_with_map
        width_fig_coverage = width_fig_with_map
        list_of_fig_coverage = []
        list_ax_coverage = []
        nb_sat_per_figure = 4
        nb_sets_of_sat = (int) (np.ceil( nb_satellites / np.float(nb_sat_per_figure) )) # total nb of sat divided by the nb of sat per figure. So this is the nb of figure per day
        nb_fig_coverage = nb_sets_of_sat * nb_days_run
        sat_on_this_fig_array = np.arange(0, nb_satellites, nb_sat_per_figure)
        which_gps_and_color = []
        for iday in range(nb_days_run):
            for iset_sat in range(nb_sets_of_sat): # For a given day, loop over the figures for each set of nb_sat_per_figure satellites
                # Here, you are on a particular figure, coresponding to the day and the set of nb_sat_per_figure satellites 
                if (iset_sat < nb_sets_of_sat -1):
                    sat_on_this_fig = np.arange(sat_on_this_fig_array[iset_sat], sat_on_this_fig_array[iset_sat+1],1)
                else:
                    sat_on_this_fig = np.arange(sat_on_this_fig_array[-1], nb_satellites,1)

                list_of_fig_coverage.append( plt.figure(figsize=(height_fig_coverage, width_fig_coverage)) )
                name_sat_on_this_fig = output_filename_propagator[sat_on_this_fig[0]].split('.')[0]
                for j in range(1, len(sat_on_this_fig)):
                    name_sat_on_this_fig = name_sat_on_this_fig + ', ' + output_filename_propagator[sat_on_this_fig[j]].split('.')[0]
                list_of_fig_coverage[-1].canvas.set_window_title('Day ' + str(iday+1) + ' - ' + name_sat_on_this_fig)
                list_of_fig_coverage[-1].set_facecolor(background_color)
                # Axes
                nb_col = (int)(max(nb_times_per_day_at_least_one_spec_on_storm[sat_on_this_fig, istorm, iday]))
                gs_coverage = gridspec.GridSpec( len(sat_on_this_fig), nb_col )
                                                # one line per satellite, one column per time when at least one spec is in the storm
                gs_coverage.update(left=0.025, right=0.99, top = 0.99,bottom = 0.03,hspace=0.08,wspace = 0.05)
                isat_count = -1
                for isat in sat_on_this_fig: # we are on a given line now
                    isat_count = isat_count + 1
                    for itime in range((int)(nb_times_per_day_at_least_one_spec_on_storm[isat,istorm,iday])): 
                        # we are on a given line and column so we are looking at a time when at lest one spec of the 4 spec of this sat flies in the storm for this day
                        # the plot here starts at list_all_start_spec_on_storm[istorm][isat][itime] and ends at list_all_stop_spec_on_storm[istorm][isat][itime]
                        list_ax_coverage.append( list_of_fig_coverage[-1].add_subplot(gs_coverage[isat_count, itime]) )
                        which_gps_list_temp = []
                        for ispec in range(nb_spec_pts): # make a list of all GPS that have spec in storm for this plot (so between  list_all_start_spec_on_storm[istorm][isat][itime] and list_all_stop_spec_on_storm[istorm][isat][itime] )
                            if len(list_times_spec_on_storm[istorm][isat][ispec]) > 0:
                                for j in range(len(list_times_spec_on_storm[istorm][isat][ispec])):
                                    if ( ( list_times_spec_on_storm[istorm][isat][ispec][j][0] >= list_all_start_spec_on_storm[istorm][isat][itime+(int)(np.sum(nb_times_per_day_at_least_one_spec_on_storm[isat,istorm,0:iday]))] ) & ( list_times_spec_on_storm[istorm][isat][ispec][j][1] <= list_all_stop_spec_on_storm[istorm][isat][itime+(int)(np.sum(nb_times_per_day_at_least_one_spec_on_storm[isat,istorm,0:iday]))] ) ):
                                        which_gps_list_temp.append( name_spec[isat][list_times_spec_on_storm[istorm][isat][ispec][j][0]][ispec] )
                        which_gps = list(set(which_gps_list_temp)) # in case the same GPS appears twice (or more) in the list which_gps_list_temp
                        which_gps_and_color_sublist = []
                        which_gps_and_color_sublist.append([list_all_start_spec_on_storm[istorm][isat][itime+(int)(np.sum(nb_times_per_day_at_least_one_spec_on_storm[isat,istorm,0:iday]))]])
                        for igps in range(len(which_gps)):
                            which_gps_and_color_sublist.append([which_gps[igps], color_spec[igps]])
                        color_of_spec_save = []
                        for ispec in range(nb_spec_pts):
                            if len(list_times_spec_on_storm[istorm][isat][ispec]) > 0:
                                for j in range(len(list_times_spec_on_storm[istorm][isat][ispec])):
                                    if ( ( list_times_spec_on_storm[istorm][isat][ispec][j][0] >= list_all_start_spec_on_storm[istorm][isat][itime+(int)(np.sum(nb_times_per_day_at_least_one_spec_on_storm[isat,istorm,0:iday]))] ) & ( list_times_spec_on_storm[istorm][isat][ispec][j][1] <= list_all_stop_spec_on_storm[istorm][isat][itime+(int)(np.sum(nb_times_per_day_at_least_one_spec_on_storm[isat,istorm,0:iday]))] ) ):
                                        start_prob = list_times_spec_on_storm[istorm][isat][ispec][j][0]
                                        end_prob = list_times_spec_on_storm[istorm][isat][ispec][j][1] + 1
                                        prob_spec = np.exp( np.log(2./3.) * distance_spec_to_storm[ispec, isat,start_prob:end_prob ]**2 / (radius_uncertainty_storm[start_prob + delay_storm : end_prob + delay_storm , istorm]**2) )
                                        color_of_spec = which_gps_and_color_sublist[which_gps.index(name_spec[isat][list_times_spec_on_storm[istorm][isat][ispec][j][0]][ispec])+1][1]
                                        color_of_spec_save.append(color_of_spec)
                                        spec_already_plot = 0
                                        for c in range(len(color_of_spec_save)-1):
                                            if color_of_spec_save[c] == color_of_spec:
                                                spec_already_plot = 1
                                        if spec_already_plot == 0:
                                            list_ax_coverage[-1].scatter(np.arange(start_prob, end_prob), prob_spec, color = color_of_spec, marker = '.', label = which_gps_and_color_sublist[which_gps.index(name_spec[isat][list_times_spec_on_storm[istorm][isat][ispec][j][0]][ispec])+1][0].split('_')[-1])
                                        else:
                                            list_ax_coverage[-1].scatter(np.arange(start_prob, end_prob), prob_spec, color = color_of_spec, marker = '.')

                        which_gps_and_color.append(which_gps_and_color_sublist)
                        list_ax_coverage[-1].margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
                        list_ax_coverage[-1].set_ylim([min_prob,max_prob])
                        list_ax_coverage[-1].xaxis.set_ticklabels([])
                        list_ax_coverage[-1].tick_params(axis='y', colors='w')
                        list_ax_coverage[-1].legend(loc= 2, ncol =  nb_spec_pts / 2,columnspacing = 0.5,scatterpoints = 1, frameon = False, handletextpad = 0.001,borderaxespad = 0.001 )
                        if (itime > 0):
                            list_ax_coverage[-1].yaxis.set_ticklabels([])
                        if len(sat_on_this_fig) == 1 :
                            coeff_height_bottom_text = 2
                        else: 
                            coeff_height_bottom_text = 1
                        list_ax_coverage[-1].text(list_all_start_spec_on_storm[istorm][isat][itime+(int)(np.sum(nb_times_per_day_at_least_one_spec_on_storm[isat,istorm,0:iday]))], min_prob-(max_prob-min_prob)/20/coeff_height_bottom_text, time_sat[list_all_start_spec_on_storm[istorm][isat][itime+(int)(np.sum(nb_times_per_day_at_least_one_spec_on_storm[isat,istorm,0:iday]))]].split('T')[-1], horizontalalignment ='left', weight = 'bold', color = 'w')
                        list_ax_coverage[-1].text(list_all_stop_spec_on_storm[istorm][isat][itime+(int)(np.sum(nb_times_per_day_at_least_one_spec_on_storm[isat,istorm,0:iday]))], min_prob-(max_prob-min_prob)/20/coeff_height_bottom_text, time_sat[list_all_stop_spec_on_storm[istorm][isat][itime+(int)(np.sum(nb_times_per_day_at_least_one_spec_on_storm[isat,istorm,0:iday]))]].split('T')[-1], horizontalalignment ='right', weight = 'bold', color = 'w')

                # plt.show()
                # plt.show()

# If show_storm is 1 then on the left of the map, plot the coverage of the specular points that are currently in the storm

    if show_graph == 1:
        ax_graph_left_of_map = fig_with_map.add_subplot(gs_with_map[0, 0])
        ax_graph_left_of_map.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
        ax_graph_left_of_map.set_ylim([min_prob,max_prob])
        ax_graph_left_of_map.set_xlim([0,1])
        ax_graph_left_of_map.xaxis.set_ticklabels([])
        ax_graph_left_of_map.tick_params(axis='y', colors='w')
        graph_left_of_map = ax_graph_left_of_map.scatter([],[])

#########################################################################################
############################################################################### ANIMATION

# Text for time
pos_time_text_ax_map_x = max_lon[zoom]
pos_time_text_ax_map_y = min_lat[zoom] + ( height_fig_with_map - width_fig_with_map ) / 40.
time_text_ax_map = ax_map.text(pos_time_text_ax_map_x, pos_time_text_ax_map_y,'',horizontalalignment ='right', weight = 'bold',fontsize = 15, color = 'b')

# Text for specular points
pos_spec_text_ax_map_x = min_lon[zoom] 
pos_spec_text_ax_map_y_0 = min_lat[zoom] + ( height_fig_with_map - width_fig_with_map ) / 40.
dy_spec_text = ( height_fig_with_map - width_fig_with_map ) / 5.
pos_spec_text_ax_map_y = np.arange(pos_spec_text_ax_map_y_0,pos_spec_text_ax_map_y_0+nb_spec_pts*dy_spec_text,dy_spec_text)
spec_text_ax_map_list = []
for ispec in range(nb_spec_pts):
    spec_text_ax_map_list.append(ax_map.text(pos_spec_text_ax_map_x, pos_spec_text_ax_map_y[ispec],'',horizontalalignment ='left', weight = 'bold',fontsize = 15, color = color_spec[ispec]))

if (show_storm == 1):
    # Vertical line for time on graph left of map
    if show_graph == 1:
        time_now, = ax_graph_left_of_map.plot([], [], lw=2)
    # Text for time on graph left of map  
        time_text_left_ax_graph_left_of_map = ax_graph_left_of_map.text(0, min_prob+(max_prob-min_prob)/500,'',horizontalalignment ='left', weight = 'bold',fontsize = 15, color = 'b')
        time_text_right_ax_graph_left_of_map = ax_graph_left_of_map.text(1, min_prob+(max_prob-min_prob)/500,'',horizontalalignment ='right', weight = 'bold',fontsize = 15, color = 'b')

# Initialization of the animation
def init_prop():
    # SATELLITES AND SPECULAR POINTS
    tuple_spacecraft_point_to_plot = ()
    tuple_specular_point_to_plot = ()
    for isat in range(nb_satellites):
        spacecraft_list[isat].point_plot.set_data([], []) # set_data must be applied to a line object
        tuple_spacecraft_point_to_plot = tuple_spacecraft_point_to_plot + (spacecraft_list[isat].point_plot,)
        for ispec in range(nb_spec_pts):
            specular_list[ispec+isat*nb_spec_pts].point_plot.set_data([], [])
            tuple_specular_point_to_plot = tuple_specular_point_to_plot + (specular_list[ispec+isat*nb_spec_pts].point_plot,)
    # STORMS
    tuple_storm_point_to_plot = ()
    if (show_storm == 1):
        for istorm in range(nb_storms):
            storm_list[istorm].point_plot.set_data([], []) 
            tuple_storm_point_to_plot = tuple_storm_point_to_plot + (storm_list[istorm].point_plot,) 
        # Vertical line for time on graph left of map
        if show_graph == 1:
            time_now.set_data([], [])
            # Text for time on graph left of map 
            time_text_left_ax_graph_left_of_map.set_text('')
            time_text_right_ax_graph_left_of_map.set_text('')
    # TIME
    time_text_ax_map.set_text('')
    tuple_spec_text_ax_map = ()
    for i in range(nb_spec_pts):
        spec_text_ax_map_list[ispec].set_text('')
        tuple_spec_text_ax_map = tuple_spec_text_ax_map + (spec_text_ax_map_list[ispec],)


    # GRAPH ON LEFT OF MAP
    if show_graph == 1:
        graph_left_of_map.set_offsets(([],[]))
    
    if show_graph == 1 :
        tuple_spacecraft_and_specular_and_storm_and_text = tuple_spacecraft_point_to_plot + tuple_specular_point_to_plot + tuple_storm_point_to_plot + (time_text_ax_map,) + tuple_spec_text_ax_map + (time_now,) + (time_text_right_ax_graph_left_of_map,) + (time_text_left_ax_graph_left_of_map,)
    else:
        tuple_spacecraft_and_specular_and_storm_and_text = tuple_spacecraft_point_to_plot + tuple_specular_point_to_plot + tuple_storm_point_to_plot + (time_text_ax_map,) + tuple_spec_text_ax_map 


    return tuple_spacecraft_and_specular_and_storm_and_text
# init_prop() is just here to tell the animator which objects on the plot to update after each frame. init_prop() needs to return a tuple. Each element of the tuple must be a line object. What init_prop() returns tells the animator which objects on the plot to update after each frame    

raise Exception
# Animation
start_visu = '2016-02-05T00:00:00'
stop_visu = '2016-02-05T00:05:00'
when_start_visu = time_sat.index(start_visu)
when_stop_visu = time_sat.index(stop_visu)
delay_show_prob_before = 0 # in seconds # the probability distributino functions of the spec on the left of the map will showp up before the spec enters 
delay_show_prob_after = 0 # in seconds # the probability distributino functions of the spec on the left of the map will stay after the spec leaves
delay_show_prob_before = delay_show_prob_before / interpolation_step
delay_show_prob_after = delay_show_prob_after / interpolation_step
if (show_storm == 1):
    when_start_visu_in_storm_file = time_storm.index(start_visu)
print time_sat[when_start_visu] + ' -> ' + time_sat[when_stop_visu]
def animate_prop(j): # i is the frame number
    i  = j + when_start_visu
    tuple_spacecraft_point_to_plot = ()
    tuple_specular_point_to_plot = ()
    tuple_storm_point_to_plot = ()
    tuple_spec_text_ax_map = ()
#    global time_first_spec_enter
    for isat in range(nb_satellites):
        # SPACECRAFT
        lon_a, lat_a = lon_sat[isat,i], lat_sat[isat,i]
        spacecraft_list[isat].x,  spacecraft_list[isat].y = m(lon_a, lat_a)
        spacecraft_list[isat].point_plot.set_data(spacecraft_list[isat].x, spacecraft_list[isat].y)
        tuple_spacecraft_point_to_plot = tuple_spacecraft_point_to_plot + (spacecraft_list[isat].point_plot,)        
        # SPECULAR POINTS
        kspec = 0
        # if there is less specular points than nb_spec_pts for for this iteration and satellite, then the old specular point stays on the image. So for loop blow is to make sure that does not happen (there are smarter ways to do that...)
        for ispec in range(nb_spec_pts):
            specular_list[ispec + isat*nb_spec_pts].point_plot.set_data([], []) 
        for ispec in range(len(name_spec[isat][i])):
            lon_spec_a, lat_spec_a = lon_spec[ispec,isat,i], lat_spec[ispec,isat,i]
            specular_list[ispec + isat*nb_spec_pts].x,  specular_list[ispec + isat*nb_spec_pts].y = m(lon_spec_a, lat_spec_a)
            specular_list[ispec+isat*nb_spec_pts].point_plot.set_data(specular_list[ispec+isat*nb_spec_pts].x, specular_list[ispec+isat*nb_spec_pts].y)
            if specular_over_storm[ispec,isat,i] ==  1:
                specular_list[ispec+isat*nb_spec_pts].point_plot.set_color('r')
                spec_text_ax_map_list[ispec].set_text(output_filename_propagator[isat].split('.')[0] + '-' + name_spec[isat][i][ispec].split('_')[-1])
            else:
                specular_list[ispec+isat*nb_spec_pts].point_plot.set_color('chartreuse')
                spec_text_ax_map_list[ispec].set_text('')
            tuple_spec_text_ax_map = tuple_spec_text_ax_map + (spec_text_ax_map_list[ispec],)
            tuple_specular_point_to_plot = tuple_specular_point_to_plot + (specular_list[ispec+isat*nb_spec_pts].point_plot,)

    # STORMS
    if (show_storm == 1):
        for istorm in range(nb_storms):       
            # Hurricane changes color when a spec is in it
            if (len(np.where(specular_over_storm[:, :, i] == 1)[0]) > 0):
                color_storm = 'r'
            else:
                color_storm = 'b'
            storm_list[istorm].point_plot.set_color(color_storm)
            no_spec_from_any_sat = 1
            if show_graph == 1:
                time_text_left_ax_graph_left_of_map.set_text('')
                time_text_right_ax_graph_left_of_map.set_text('')
            for isat2 in range(nb_satellites):
                # On the left of the map, plot the probability distribution functions of the spec in storm
                if ( (i < nb_steps_interpolation - delay_show_prob_before) & (i > delay_show_prob_after) ):
                    list_scatter = []
                    color_spec_tuple = () 
                    if ( (len(np.where(specular_over_storm[:, isat2, i+delay_show_prob_before] == 1)[0]) > 0) | (len(np.where(specular_over_storm[:, isat2, i-delay_show_prob_after] == 1)[0]) > 0) ):                        
                        # determine the first spec that enters the storm (and the last spec that leaves the storm) for each fly over of the sat 
                        for j in range(len(list_all_stop_spec_on_storm[istorm][isat2])):
                            if ((i+delay_show_prob_before >= list_all_start_spec_on_storm[istorm][isat2][j]) & (i-delay_show_prob_after<= list_all_stop_spec_on_storm[istorm][isat2][j])):
                                time_first_spec_enter = list_all_start_spec_on_storm[istorm][isat2][j]
                                time_last_spec_leave = list_all_stop_spec_on_storm[istorm][isat2][j]

                                if show_graph == 1:
                                    # Vertical line for time on graph left of map
                                    time_now.set_data([ ( i - time_first_spec_enter) /np.float(time_last_spec_leave-time_first_spec_enter), ( i - time_first_spec_enter) /np.float(time_last_spec_leave-time_first_spec_enter)   ], [min_prob, max_prob])
                                # Text for time on graph left of map 
                                    time_text_left_ax_graph_left_of_map.set_text(time_sat[time_first_spec_enter][-8:19])
                                    time_text_right_ax_graph_left_of_map.set_text(time_sat[time_last_spec_leave][-8:19])

                        # Find the color of the spec so they are the same as the plot of the coverage
                        for j in range(len(which_gps_and_color)):
                            if ( which_gps_and_color[j][0] == time_first_spec_enter ):
                                list_of_color_for_spec = []
                                list_of_gps_for_spec = []
                                for k in range(1,len(which_gps_and_color[j])):
                                    list_of_gps_for_spec.append(which_gps_and_color[j][k][0])
                                    list_of_color_for_spec.append(which_gps_and_color[j][k][1])
                        # now draw the probability distribution functinos for each spec
                        for ispec in range(nb_spec_pts):
                            for j in range(len(list_times_spec_on_storm[istorm][isat2][ispec])):
                                if ( ( list_times_spec_on_storm[istorm][isat2][ispec][j][0] <= i+delay_show_prob_before ) & ( list_times_spec_on_storm[istorm][isat2][ispec][j][1] >= i-delay_show_prob_after ) ):
                                    start_prob = list_times_spec_on_storm[istorm][isat2][ispec][j][0]
                                    end_prob = list_times_spec_on_storm[istorm][isat2][ispec][j][1] + 1
                                    prob_spec = np.exp( np.log(2./3.) * distance_spec_to_storm[ispec, isat2,start_prob:end_prob ]**2 / (radius_uncertainty_storm[start_prob + delay_storm : end_prob + delay_storm , istorm]**2) )
                                    color_of_spec = list_of_color_for_spec[list_of_gps_for_spec.index(name_spec[isat2][list_times_spec_on_storm[istorm][isat2][ispec][j][0]][ispec])]
                                    color_of_spec_save.append(color_of_spec)
                                    spec_already_plot = 0
                                    for c in range(len(color_of_spec_save)-1):
                                        if color_of_spec_save[c] == color_of_spec:
                                            spec_already_plot = 1
#                                    graph_left_of_map.set_offsets(([0.5], [0.8]))
                                    no_spec_from_any_sat = 0
                                    color_spec_tuple =  color_spec_tuple + len(prob_spec)*(color_of_spec,)
                                    for ilist in range(len(prob_spec)):
                                        list_scatter.append( [ np.arange(0+(start_prob - time_first_spec_enter)/(np.float(time_last_spec_leave-time_first_spec_enter)), 1- (time_last_spec_leave - end_prob)/(np.float(time_last_spec_leave-time_first_spec_enter)),1./(np.float(time_last_spec_leave-time_first_spec_enter)))[ilist], prob_spec[ilist] ] )
                                    if show_graph == 1:
                                        if spec_already_plot == 0:
                                            graph_left_of_map.set_offsets(list_scatter)
                                            graph_left_of_map.set_color(color_spec_tuple)
                                        else:
                                            graph_left_of_map.set_offsets(list_scatter)
                                            graph_left_of_map.set_color(color_spec_tuple)
                    elif (no_spec_from_any_sat == 1):
                        if show_graph == 1:
                            graph_left_of_map.set_offsets(([], []))

            ell = get_ellipse_coords(a=radius_for_tissot(radius_uncertainty_storm[j+when_start_visu_in_storm_file,istorm]), b=radius_for_tissot(radius_uncertainty_storm[j+when_start_visu_in_storm_file,istorm]), x=lon_storm[j+when_start_visu_in_storm_file,istorm], y=lat_storm[j+when_start_visu_in_storm_file,istorm])
 
 
            storm_list[istorm].point_plot.set_data(ell[:,0], ell[:,1])
            tuple_storm_point_to_plot = tuple_storm_point_to_plot + (storm_list[istorm].point_plot,)                        
                                    
    # TIME
    time_text_ax_map.set_text(time_sat[i])

    if show_graph == 1:
        tuple_spacecraft_and_specular_and_storm_and_text = tuple_spacecraft_point_to_plot + tuple_specular_point_to_plot + tuple_storm_point_to_plot + (time_text_ax_map,) + tuple_spec_text_ax_map + (graph_left_of_map,) + (time_now,)  + (time_text_right_ax_graph_left_of_map,) + (time_text_left_ax_graph_left_of_map,)
    else:
        tuple_spacecraft_and_specular_and_storm_and_text = tuple_spacecraft_point_to_plot + tuple_specular_point_to_plot + tuple_storm_point_to_plot + (time_text_ax_map,) + tuple_spec_text_ax_map       

#    fig_with_map.savefig('./test_ani/'+str(j), facecolor=fig_with_map.get_facecolor(), edgecolor='none')

    return tuple_spacecraft_and_specular_and_storm_and_text
# what animate_prop returns tells the animation framework what parts of the plot should be animated.
	

fps_we_want = 20L

anim = animation.FuncAnimation(fig_with_map, animate_prop,init_func=init_prop, frames = when_stop_visu - when_start_visu+1, interval=1000L/fps_we_want, blit= True, repeat = True)
# blit tells the animation to only re-draw the pieces of the plot which have changed

# This block below (2 lines) makes the window fills the entire screen but does not work on every computer...
# mng = plt.get_current_fig_manager()
# mng.resize(*mng.window.maxsize())


# plt.show(block = False)
# plt.show(block = False)


